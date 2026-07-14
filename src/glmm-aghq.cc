#include "glmm-aghq.h"

#include <R_ext/Error.h>
#include <R_ext/Utils.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>

namespace {

const double LOG_2PI = std::log(2.0 * std::acos(-1.0));
const double LOG_DBL_MAX = std::log(std::numeric_limits<double>::max());


struct QuadratureRules {
  const double **nodes;
  const double **log_weights;
  int *orders;
  int n_rules;
};


struct QuadratureDiagnostics {
  int max_order;
  int max_mode_iterations;
  int exact_count;
  double max_change;
  int *order_counts;
};


struct IntegralResult {
  double value;
  double change;
  int order;
  int rule_index;
  int mode_iterations;
  bool mode_converged;
  bool quadrature_converged;
};


struct Cholesky2 {
  double r11;
  double r12;
  double r22;
};


int matrix_nrow(SEXP x, const char *name)
{
  SEXP dim = Rf_getAttrib(x, R_DimSymbol);
  if (TYPEOF(dim) != INTSXP || Rf_length(dim) != 2) {
    Rf_error("'%s' must be a matrix.", name);
  }
  return INTEGER(dim)[0];
}


int matrix_ncol(SEXP x, const char *name)
{
  SEXP dim = Rf_getAttrib(x, R_DimSymbol);
  if (TYPEOF(dim) != INTSXP || Rf_length(dim) != 2) {
    Rf_error("'%s' must be a matrix.", name);
  }
  return INTEGER(dim)[1];
}


void check_real(SEXP x, const char *name)
{
  if (TYPEOF(x) != REALSXP) {
    Rf_error("'%s' must be numeric.", name);
  }
}


void check_integer(SEXP x, const char *name)
{
  if (TYPEOF(x) != INTSXP) {
    Rf_error("'%s' must be integer.", name);
  }
}


double scalar_real(SEXP x, const char *name, bool positive)
{
  check_real(x, name);
  if (Rf_length(x) != 1) {
    Rf_error("'%s' must have length one.", name);
  }

  const double value = REAL(x)[0];
  if (!std::isfinite(value) || (positive && value <= 0.0)) {
    Rf_error("'%s' must be %s.", name,
             positive ? "finite and positive" : "finite");
  }
  return value;
}


int scalar_integer(SEXP x, const char *name, int lower)
{
  check_integer(x, name);
  if (Rf_length(x) != 1 || INTEGER(x)[0] < lower) {
    Rf_error("'%s' must be an integer greater than or equal to %d.",
             name, lower);
  }
  return INTEGER(x)[0];
}


void validate_common(SEXP mu_samples, SEXP tau_within, SEXP weights,
                     int *S, int *K)
{
  check_real(mu_samples, "mu_samples");
  check_real(tau_within, "tau_within");

  *S = matrix_nrow(mu_samples, "mu_samples");
  *K = matrix_ncol(mu_samples, "mu_samples");
  if (*S <= 0 || *K <= 0) {
    Rf_error("'mu_samples' must have positive dimensions.");
  }
  if (matrix_nrow(tau_within, "tau_within") != *S ||
      matrix_ncol(tau_within, "tau_within") != *K) {
    Rf_error("'tau_within' dimensions must match 'mu_samples'.");
  }

  // Diagnostic counts are returned as R integers. Cap the cell count before
  // multiplication so neither those counts nor native indexes can overflow.
  const int max_cells = std::numeric_limits<int>::max();
  if (*S > max_cells / *K) {
    Rf_error("AGHQ supports at most %d matrix cells.", max_cells);
  }
  const R_xlen_t n_cells = static_cast<R_xlen_t>(*S) *
    static_cast<R_xlen_t>(*K);
  if (XLENGTH(mu_samples) != n_cells || XLENGTH(tau_within) != n_cells) {
    Rf_error("AGHQ matrix dimensions do not match their storage lengths.");
  }

  const double *mu_p  = REAL(mu_samples);
  const double *tau_p = REAL(tau_within);
  for (R_xlen_t i = 0; i < n_cells; ++i) {
    if (!std::isfinite(mu_p[i]) || !std::isfinite(tau_p[i])) {
      Rf_error("'mu_samples' and 'tau_within' must contain finite values.");
    }
  }

  if (weights != R_NilValue) {
    check_real(weights, "weights");
    if (Rf_length(weights) != *K) {
      Rf_error("'weights' must have one value per observation.");
    }
    const double *weights_p = REAL(weights);
    for (int k = 0; k < *K; ++k) {
      if (!std::isfinite(weights_p[k]) || weights_p[k] <= 0.0) {
        Rf_error("'weights' must contain finite positive values.");
      }
    }
  }
}


QuadratureRules validate_rules(SEXP nodes, SEXP log_weights, int consecutive)
{
  if (TYPEOF(nodes) != VECSXP || TYPEOF(log_weights) != VECSXP ||
      Rf_length(nodes) != Rf_length(log_weights)) {
    Rf_error("Quadrature nodes and log weights must be lists of equal length.");
  }

  const int n_rules = Rf_length(nodes);
  if (n_rules < consecutive + 1) {
    Rf_error("Too few quadrature rules for the convergence requirement.");
  }

  QuadratureRules out;
  out.nodes = reinterpret_cast<const double **>(R_alloc(
    static_cast<size_t>(n_rules),
    sizeof(const double *)
  ));
  out.log_weights = reinterpret_cast<const double **>(R_alloc(
    static_cast<size_t>(n_rules),
    sizeof(const double *)
  ));
  out.orders = reinterpret_cast<int *>(R_alloc(
    static_cast<size_t>(n_rules),
    sizeof(int)
  ));
  out.n_rules = n_rules;

  int previous_order = 0;
  for (int r = 0; r < n_rules; ++r) {
    SEXP nodes_r       = VECTOR_ELT(nodes, r);
    SEXP log_weights_r = VECTOR_ELT(log_weights, r);
    check_real(nodes_r, "nodes[[r]]");
    check_real(log_weights_r, "log_weights[[r]]");

    const int order = Rf_length(nodes_r);
    if (order < 3 || order % 2 == 0 || order <= previous_order ||
        Rf_length(log_weights_r) != order) {
      Rf_error("Quadrature rules must have matching, increasing odd orders.");
    }

    const double *nodes_p = REAL(nodes_r);
    const double *log_w_p = REAL(log_weights_r);
    double max_log_weight = -std::numeric_limits<double>::infinity();
    for (int j = 0; j < order; ++j) {
      if (!std::isfinite(nodes_p[j]) || std::isnan(log_w_p[j]) ||
          log_w_p[j] == std::numeric_limits<double>::infinity()) {
        Rf_error("Quadrature nodes and nonzero weights must be finite.");
      }
      max_log_weight = std::max(max_log_weight, log_w_p[j]);
    }

    if (!std::isfinite(max_log_weight)) {
      Rf_error("Each quadrature rule must contain a positive weight.");
    }

    double weight_sum = 0.0;
    for (int j = 0; j < order; ++j) {
      weight_sum += std::exp(log_w_p[j] - max_log_weight);
    }
    const double log_weight_sum = max_log_weight + std::log(weight_sum);
    if (std::fabs(log_weight_sum) > 1e-8) {
      Rf_error("Each quadrature rule must have weights summing to one.");
    }

    out.nodes[r]       = nodes_p;
    out.log_weights[r] = log_w_p;
    out.orders[r]      = order;
    previous_order = order;
  }

  return out;
}


double logistic(double x)
{
  if (x >= 0.0) {
    return 1.0 / (1.0 + std::exp(-x));
  }
  const double ex = std::exp(x);
  return ex / (1.0 + ex);
}


double softplus(double x)
{
  if (x >= 0.0) {
    return x + std::log1p(std::exp(-x));
  }
  return std::log1p(std::exp(x));
}


double log_inv_logit(double x)
{
  return -softplus(-x);
}


double log_one_minus_inv_logit(double x)
{
  return -softplus(x);
}


double log_add_exp(double x, double y)
{
  const double maximum = std::max(x, y);
  return maximum + std::log(std::exp(x - maximum) + std::exp(y - maximum));
}


bool cholesky2(double h11, double h12, double h22, Cholesky2 *out)
{
  if (!std::isfinite(h11) || !std::isfinite(h12) || !std::isfinite(h22) ||
      h11 <= 0.0) {
    return false;
  }

  out->r11 = std::sqrt(h11);
  out->r12 = h12 / out->r11;
  const double remainder = h22 - out->r12 * out->r12;
  if (!std::isfinite(remainder) || remainder <= 0.0) {
    return false;
  }
  out->r22 = std::sqrt(remainder);
  return true;
}


void solve_hessian(const Cholesky2 &chol, double g1, double g2,
                   double *step1, double *step2)
{
  const double y1 = g1 / chol.r11;
  const double y2 = (g2 - chol.r12 * y1) / chol.r22;
  *step2 = y2 / chol.r22;
  *step1 = (y1 - chol.r12 * *step2) / chol.r11;
}


void solve_cholesky(const Cholesky2 &chol, double z1, double z2,
                    double *delta1, double *delta2)
{
  *delta2 = z2 / chol.r22;
  *delta1 = (z1 - chol.r12 * *delta2) / chol.r11;
}


struct BinomialProblem {
  int a;
  int c;
  int n1;
  int n2;
  double mu;
  double tau;
  double weight;
  double alpha;
  double beta;
  double log_coefficient;
  double log_beta;

  double log_density(double theta, double logit_pi) const
  {
    const double effect = mu + tau * theta;
    const double eta1   = logit_pi + 0.5 * effect;
    const double eta2   = logit_pi - 0.5 * effect;
    const double log_likelihood = log_coefficient +
      a * log_inv_logit(eta1) +
      (n1 - a) * log_one_minus_inv_logit(eta1) +
      c * log_inv_logit(eta2) +
      (n2 - c) * log_one_minus_inv_logit(eta2);
    const double log_prior =
      -0.5 * theta * theta - 0.5 * LOG_2PI +
      alpha * log_inv_logit(logit_pi) +
      beta * log_one_minus_inv_logit(logit_pi) - log_beta;

    return weight * log_likelihood + log_prior;
  }

  bool derivatives(double theta, double logit_pi, double *value,
                   double *g1, double *g2,
                   double *h11, double *h12, double *h22) const
  {
    *value = log_density(theta, logit_pi);
    if (!std::isfinite(*value)) {
      return false;
    }

    const double effect = mu + tau * theta;
    const double p1     = logistic(logit_pi + 0.5 * effect);
    const double p2     = logistic(logit_pi - 0.5 * effect);
    const double pi0    = logistic(logit_pi);
    const double r1     = a - n1 * p1;
    const double r2     = c - n2 * p2;
    const double v1     = n1 * p1 * (1.0 - p1);
    const double v2     = n2 * p2 * (1.0 - p2);

    *g1  = weight * 0.5 * tau * (r1 - r2) - theta;
    *g2  = weight * (r1 + r2) + alpha - (alpha + beta) * pi0;
    *h11 = 1.0 + weight * 0.25 * tau * tau * (v1 + v2);
    *h12 = weight * 0.5 * tau * (v1 - v2);
    *h22 = weight * (v1 + v2) +
      (alpha + beta) * pi0 * (1.0 - pi0);

    return std::isfinite(*g1) && std::isfinite(*g2) &&
      std::isfinite(*h11) && std::isfinite(*h12) && std::isfinite(*h22);
  }

  void initial(double *theta, double *nuisance) const
  {
    *theta = 0.0;
    const double events = static_cast<double>(a) + c;
    const double total  = static_cast<double>(n1) + n2;
    double pi0 = (alpha + weight * events) /
      (alpha + beta + weight * total);
    pi0 = std::max(1e-12, std::min(1.0 - 1e-12, pi0));
    *nuisance = std::log(pi0) - std::log1p(-pi0);
  }
};


struct PoissonProblem {
  int x1;
  int x2;
  double log_t1;
  double log_t2;
  double mu;
  double tau;
  double weight;
  double mean;
  double sd;
  double variance;
  double log_factorial;

  bool rates(double theta, double phi, double *log_lambda1,
             double *log_lambda2, double *lambda1, double *lambda2) const
  {
    const double effect = mu + tau * theta;
    *log_lambda1 = log_t1 + phi + 0.5 * effect;
    *log_lambda2 = log_t2 + phi - 0.5 * effect;
    if (!std::isfinite(*log_lambda1) || !std::isfinite(*log_lambda2) ||
        *log_lambda1 > LOG_DBL_MAX || *log_lambda2 > LOG_DBL_MAX) {
      return false;
    }

    *lambda1 = std::exp(*log_lambda1);
    *lambda2 = std::exp(*log_lambda2);
    return std::isfinite(*lambda1) && std::isfinite(*lambda2);
  }

  double log_density(double theta, double phi) const
  {
    double log_lambda1, log_lambda2, lambda1, lambda2;
    if (!rates(theta, phi, &log_lambda1, &log_lambda2,
               &lambda1, &lambda2)) {
      return -std::numeric_limits<double>::infinity();
    }

    const double log_likelihood =
      x1 * log_lambda1 - lambda1 +
      x2 * log_lambda2 - lambda2 - log_factorial;
    const double standardized = (phi - mean) / sd;
    const double log_prior =
      -0.5 * theta * theta - 0.5 * LOG_2PI -
      0.5 * standardized * standardized - std::log(sd) - 0.5 * LOG_2PI;

    return weight * log_likelihood + log_prior;
  }

  bool derivatives(double theta, double phi, double *value,
                   double *g1, double *g2,
                   double *h11, double *h12, double *h22) const
  {
    double log_lambda1, log_lambda2, lambda1, lambda2;
    if (!rates(theta, phi, &log_lambda1, &log_lambda2,
               &lambda1, &lambda2)) {
      return false;
    }
    *value = log_density(theta, phi);
    if (!std::isfinite(*value)) {
      return false;
    }

    const double r1 = x1 - lambda1;
    const double r2 = x2 - lambda2;
    *g1  = weight * 0.5 * tau * (r1 - r2) - theta;
    *g2  = weight * (r1 + r2) - (phi - mean) / variance;
    *h11 = 1.0 + weight * 0.25 * tau * tau * (lambda1 + lambda2);
    *h12 = weight * 0.5 * tau * (lambda1 - lambda2);
    *h22 = weight * (lambda1 + lambda2) + 1.0 / variance;

    return std::isfinite(*g1) && std::isfinite(*g2) &&
      std::isfinite(*h11) && std::isfinite(*h12) && std::isfinite(*h22);
  }

  void initial(double *theta, double *nuisance) const
  {
    *theta = 0.0;
    const double log_exposure = log_add_exp(
      log_t1 + 0.5 * mu,
      log_t2 - 0.5 * mu
    );
    *nuisance = std::log(static_cast<double>(x1) + x2 + 0.5) -
      log_exposure;
  }
};


template <typename Problem>
bool find_mode(const Problem &problem, double tolerance,
               double *theta, double *nuisance, Cholesky2 *chol,
               int *iterations)
{
  problem.initial(theta, nuisance);
  const int max_iterations = 60;

  for (int iteration = 1; iteration <= max_iterations; ++iteration) {
    double value, g1, g2, h11, h12, h22;
    if (!problem.derivatives(*theta, *nuisance, &value, &g1, &g2,
                             &h11, &h12, &h22) ||
        !cholesky2(h11, h12, h22, chol)) {
      return false;
    }

    double step1, step2;
    solve_hessian(*chol, g1, g2, &step1, &step2);
    const double improvement = g1 * step1 + g2 * step2;
    if (!std::isfinite(improvement) || improvement < 0.0) {
      return false;
    }
    if (0.5 * improvement <= tolerance) {
      *iterations = iteration;
      return true;
    }

    bool accepted     = false;
    double step_scale = 1.0;
    for (int backtrack = 0; backtrack < 60; ++backtrack) {
      const double candidate_theta    = *theta + step_scale * step1;
      const double candidate_nuisance = *nuisance + step_scale * step2;
      const double candidate_value = problem.log_density(
        candidate_theta,
        candidate_nuisance
      );
      if (std::isfinite(candidate_value) &&
          candidate_value >= value + 1e-4 * step_scale * improvement) {
        *theta    = candidate_theta;
        *nuisance = candidate_nuisance;
        accepted  = true;
        break;
      }
      step_scale *= 0.5;
    }
    if (!accepted) {
      return false;
    }
  }

  return false;
}


template <typename Problem>
double evaluate_rule(const Problem &problem, double theta, double nuisance,
                     const Cholesky2 &chol, const double *nodes,
                     const double *log_weights, int order)
{
  double term_max = -std::numeric_limits<double>::infinity();
  double term_sum = 0.0;

  for (int j = 0; j < order; ++j) {
    for (int i = 0; i < order; ++i) {
      double delta_theta, delta_nuisance;
      solve_cholesky(chol, nodes[i], nodes[j],
                     &delta_theta, &delta_nuisance);
      const double log_density = problem.log_density(
        theta + delta_theta,
        nuisance + delta_nuisance
      );
      if (!std::isfinite(log_density)) {
        continue;
      }

      const double term = log_density +
        0.5 * (nodes[i] * nodes[i] + nodes[j] * nodes[j]) +
        log_weights[i] + log_weights[j];
      if (!std::isfinite(term)) {
        continue;
      }

      if (term > term_max) {
        term_sum = std::isfinite(term_max) ?
          term_sum * std::exp(term_max - term) + 1.0 : 1.0;
        term_max = term;
      } else {
        term_sum += std::exp(term - term_max);
      }
    }
  }

  if (!std::isfinite(term_max) || !std::isfinite(term_sum) || term_sum <= 0.0) {
    return -std::numeric_limits<double>::infinity();
  }

  return LOG_2PI - std::log(chol.r11) - std::log(chol.r22) +
    term_max + std::log(term_sum);
}


template <typename Problem>
IntegralResult integrate_problem(const Problem &problem,
                                 const QuadratureRules &rules,
                                 double tolerance, int consecutive,
                                 double mode_tolerance)
{
  IntegralResult result;
  result.value                  = NA_REAL;
  result.change                 = std::numeric_limits<double>::infinity();
  result.order                  = 0;
  result.rule_index             = -1;
  result.mode_iterations        = 0;
  result.mode_converged         = false;
  result.quadrature_converged   = false;

  double theta, nuisance;
  Cholesky2 chol;
  result.mode_converged = find_mode(
    problem,
    mode_tolerance,
    &theta,
    &nuisance,
    &chol,
    &result.mode_iterations
  );
  if (!result.mode_converged) {
    return result;
  }

  double previous = NA_REAL;
  int below_count = 0;
  for (int r = 0; r < rules.n_rules; ++r) {
    const double current = evaluate_rule(
      problem,
      theta,
      nuisance,
      chol,
      rules.nodes[r],
      rules.log_weights[r],
      rules.orders[r]
    );
    if (!std::isfinite(current)) {
      return result;
    }

    if (r > 0) {
      result.change = std::fabs(current - previous);
      below_count = result.change <= tolerance ? below_count + 1 : 0;
    }
    previous          = current;
    result.value      = current;
    result.order      = rules.orders[r];
    result.rule_index = r;

    if (below_count >= consecutive) {
      result.quadrature_converged = true;
      return result;
    }
  }

  return result;
}


void update_diagnostics(const IntegralResult &result,
                        QuadratureDiagnostics *diagnostics)
{
  diagnostics->max_order = std::max(diagnostics->max_order, result.order);
  diagnostics->max_mode_iterations = std::max(
    diagnostics->max_mode_iterations,
    result.mode_iterations
  );
  diagnostics->max_change = std::max(
    diagnostics->max_change,
    result.change
  );
  diagnostics->order_counts[result.rule_index] += 1;
}


QuadratureDiagnostics initialize_diagnostics(int n_rules)
{
  QuadratureDiagnostics diagnostics;
  diagnostics.max_order           = 0;
  diagnostics.max_mode_iterations = 0;
  diagnostics.exact_count         = 0;
  diagnostics.max_change          = 0.0;
  diagnostics.order_counts = reinterpret_cast<int *>(R_alloc(
    static_cast<size_t>(n_rules),
    sizeof(int)
  ));
  std::fill(
    diagnostics.order_counts,
    diagnostics.order_counts + n_rules,
    0
  );
  return diagnostics;
}


SEXP make_result(SEXP value, const QuadratureDiagnostics &diagnostics,
                 const QuadratureRules &rules)
{
  const int n_fields = 6;
  SEXP out       = PROTECT(Rf_allocVector(VECSXP, n_fields));
  SEXP names     = PROTECT(Rf_allocVector(STRSXP, n_fields));
  SEXP max_order = PROTECT(Rf_ScalarInteger(diagnostics.max_order));
  SEXP max_change = PROTECT(Rf_ScalarReal(diagnostics.max_change));
  SEXP max_mode_iterations = PROTECT(
    Rf_ScalarInteger(diagnostics.max_mode_iterations)
  );
  SEXP exact_count = PROTECT(Rf_ScalarInteger(diagnostics.exact_count));
  SEXP order_counts = PROTECT(
    Rf_allocVector(INTSXP, rules.n_rules)
  );
  SEXP order_names = PROTECT(
    Rf_allocVector(STRSXP, rules.n_rules)
  );

  for (int r = 0; r < rules.n_rules; ++r) {
    INTEGER(order_counts)[r] = diagnostics.order_counts[r];
    char buffer[32];
    std::snprintf(buffer, sizeof(buffer), "%d", rules.orders[r]);
    SET_STRING_ELT(order_names, r, Rf_mkChar(buffer));
  }
  Rf_setAttrib(order_counts, R_NamesSymbol, order_names);

  SET_VECTOR_ELT(out, 0, value);
  SET_VECTOR_ELT(out, 1, max_order);
  SET_VECTOR_ELT(out, 2, max_change);
  SET_VECTOR_ELT(out, 3, max_mode_iterations);
  SET_VECTOR_ELT(out, 4, order_counts);
  SET_VECTOR_ELT(out, 5, exact_count);
  SET_STRING_ELT(names, 0, Rf_mkChar("value"));
  SET_STRING_ELT(names, 1, Rf_mkChar("max_order"));
  SET_STRING_ELT(names, 2, Rf_mkChar("max_change"));
  SET_STRING_ELT(names, 3, Rf_mkChar("max_mode_iterations"));
  SET_STRING_ELT(names, 4, Rf_mkChar("order_counts"));
  SET_STRING_ELT(names, 5, Rf_mkChar("exact_count"));
  Rf_setAttrib(out, R_NamesSymbol, names);

  UNPROTECT(8);
  return out;
}


SEXP run_binomial(SEXP ai, SEXP ci, SEXP n1i, SEXP n2i,
                  SEXP mu_samples, SEXP tau_within, SEXP weights,
                  SEXP alpha, SEXP beta, SEXP nodes, SEXP log_weights,
                  SEXP tolerance, SEXP consecutive, SEXP mode_tolerance,
                  bool row_sum)
{
  check_integer(ai, "ai");
  check_integer(ci, "ci");
  check_integer(n1i, "n1i");
  check_integer(n2i, "n2i");

  int S, K;
  validate_common(mu_samples, tau_within, weights, &S, &K);
  if (Rf_length(ai) != K || Rf_length(ci) != K ||
      Rf_length(n1i) != K || Rf_length(n2i) != K) {
    Rf_error("Binomial outcome vectors must have length K.");
  }

  const double alpha_value = scalar_real(alpha, "alpha", true);
  const double beta_value  = scalar_real(beta, "beta", true);
  const double tolerance_value = scalar_real(tolerance, "tolerance", true);
  const double mode_tolerance_value = scalar_real(
    mode_tolerance,
    "mode_tolerance",
    true
  );
  const int consecutive_value = scalar_integer(consecutive, "consecutive", 1);
  const QuadratureRules rules = validate_rules(
    nodes,
    log_weights,
    consecutive_value
  );

  const int *a_p         = INTEGER(ai);
  const int *c_p         = INTEGER(ci);
  const int *n1_p        = INTEGER(n1i);
  const int *n2_p        = INTEGER(n2i);
  const double *mu_p     = REAL(mu_samples);
  const double *tau_p    = REAL(tau_within);
  const double *weight_p = weights == R_NilValue ? NULL : REAL(weights);

  for (int k = 0; k < K; ++k) {
    if (n1_p[k] < 0 || n2_p[k] < 0 || a_p[k] < 0 || c_p[k] < 0 ||
        a_p[k] > n1_p[k] || c_p[k] > n2_p[k]) {
      Rf_error("Invalid binomial counts at observation %d.", k + 1);
    }
  }

  SEXP value = PROTECT(row_sum ?
    Rf_allocVector(REALSXP, S) : Rf_allocMatrix(REALSXP, S, K));
  double *value_p = REAL(value);
  if (row_sum) {
    std::fill(value_p, value_p + S, 0.0);
  }

  QuadratureDiagnostics diagnostics = initialize_diagnostics(rules.n_rules);
  for (int k = 0; k < K; ++k) {
    R_CheckUserInterrupt();

    const double log_coefficient =
      std::lgamma(n1_p[k] + 1.0) - std::lgamma(a_p[k] + 1.0) -
      std::lgamma(n1_p[k] - a_p[k] + 1.0) +
      std::lgamma(n2_p[k] + 1.0) - std::lgamma(c_p[k] + 1.0) -
      std::lgamma(n2_p[k] - c_p[k] + 1.0);

    for (int s = 0; s < S; ++s) {
      if ((s & 255) == 0) {
        R_CheckUserInterrupt();
      }
      const R_xlen_t cell = static_cast<R_xlen_t>(s) +
        static_cast<R_xlen_t>(S) * static_cast<R_xlen_t>(k);
      BinomialProblem problem;
      problem.a               = a_p[k];
      problem.c               = c_p[k];
      problem.n1              = n1_p[k];
      problem.n2              = n2_p[k];
      problem.mu              = mu_p[cell];
      problem.tau             = tau_p[cell];
      problem.weight          = weight_p == NULL ? 1.0 : weight_p[k];
      problem.alpha           = alpha_value;
      problem.beta            = beta_value;
      problem.log_coefficient = log_coefficient;
      problem.log_beta        = std::lgamma(alpha_value) +
        std::lgamma(beta_value) - std::lgamma(alpha_value + beta_value);

      IntegralResult result;
      const bool exact = problem.mu == 0.0 && problem.tau == 0.0;
      if (exact) {
        const double total      = static_cast<double>(problem.n1) + problem.n2;
        const double observed   = static_cast<double>(problem.a) + problem.c;
        const double events     = problem.weight * observed;
        const double nonevents  = problem.weight * (
          total - observed
        );
        const double alpha_post = problem.alpha + events;
        const double beta_post  = problem.beta + nonevents;

        result.value = problem.weight * problem.log_coefficient +
          std::lgamma(alpha_post) + std::lgamma(beta_post) -
          std::lgamma(alpha_post + beta_post) - problem.log_beta;
        if (!std::isfinite(result.value)) {
          Rf_error(
            "Exact powered beta-binomial value is nonfinite at sample %d, "
            "observation %d.",
            s + 1, k + 1
          );
        }
        result.change               = 0.0;
        result.order                = 0;
        result.rule_index           = -1;
        result.mode_iterations      = 0;
        result.mode_converged       = true;
        result.quadrature_converged = true;
      } else {
        result = integrate_problem(
          problem,
          rules,
          tolerance_value,
          consecutive_value,
          mode_tolerance_value
        );
      }
      if (!result.mode_converged) {
        Rf_error("Binomial AGHQ mode failed at sample %d, observation %d.",
                 s + 1, k + 1);
      }
      if (!result.quadrature_converged) {
        Rf_error(
          "Binomial AGHQ failed to converge at sample %d, observation %d "
          "(order %d, change %.8g).",
          s + 1, k + 1, result.order, result.change
        );
      }

      if (row_sum) {
        value_p[s] += result.value;
      } else {
        value_p[cell] = result.value;
      }
      if (exact) {
        diagnostics.exact_count += 1;
      } else {
        update_diagnostics(result, &diagnostics);
      }
    }
  }

  SEXP out = PROTECT(make_result(value, diagnostics, rules));
  UNPROTECT(2);
  return out;
}


SEXP run_poisson(SEXP x1i, SEXP x2i, SEXP t1i, SEXP t2i,
                 SEXP mu_samples, SEXP tau_within, SEXP weights,
                 SEXP mean, SEXP sd, SEXP nodes, SEXP log_weights,
                 SEXP tolerance, SEXP consecutive, SEXP mode_tolerance,
                 bool row_sum)
{
  check_integer(x1i, "x1i");
  check_integer(x2i, "x2i");
  check_real(t1i, "t1i");
  check_real(t2i, "t2i");

  int S, K;
  validate_common(mu_samples, tau_within, weights, &S, &K);
  if (Rf_length(x1i) != K || Rf_length(x2i) != K ||
      Rf_length(t1i) != K || Rf_length(t2i) != K) {
    Rf_error("Poisson outcome vectors must have length K.");
  }

  const double mean_value = scalar_real(mean, "mean", false);
  const double sd_value   = scalar_real(sd, "sd", true);
  const double tolerance_value = scalar_real(tolerance, "tolerance", true);
  const double mode_tolerance_value = scalar_real(
    mode_tolerance,
    "mode_tolerance",
    true
  );
  const int consecutive_value = scalar_integer(consecutive, "consecutive", 1);
  const QuadratureRules rules = validate_rules(
    nodes,
    log_weights,
    consecutive_value
  );

  const int *x1_p        = INTEGER(x1i);
  const int *x2_p        = INTEGER(x2i);
  const double *t1_p     = REAL(t1i);
  const double *t2_p     = REAL(t2i);
  const double *mu_p     = REAL(mu_samples);
  const double *tau_p    = REAL(tau_within);
  const double *weight_p = weights == R_NilValue ? NULL : REAL(weights);

  for (int k = 0; k < K; ++k) {
    if (x1_p[k] < 0 || x2_p[k] < 0 ||
        !std::isfinite(t1_p[k]) || !std::isfinite(t2_p[k]) ||
        t1_p[k] <= 0.0 || t2_p[k] <= 0.0) {
      Rf_error("Invalid Poisson data at observation %d.", k + 1);
    }
  }

  SEXP value = PROTECT(row_sum ?
    Rf_allocVector(REALSXP, S) : Rf_allocMatrix(REALSXP, S, K));
  double *value_p = REAL(value);
  if (row_sum) {
    std::fill(value_p, value_p + S, 0.0);
  }

  QuadratureDiagnostics diagnostics = initialize_diagnostics(rules.n_rules);
  for (int k = 0; k < K; ++k) {
    R_CheckUserInterrupt();

    const double log_factorial =
      std::lgamma(x1_p[k] + 1.0) + std::lgamma(x2_p[k] + 1.0);

    for (int s = 0; s < S; ++s) {
      if ((s & 255) == 0) {
        R_CheckUserInterrupt();
      }
      const R_xlen_t cell = static_cast<R_xlen_t>(s) +
        static_cast<R_xlen_t>(S) * static_cast<R_xlen_t>(k);
      PoissonProblem problem;
      problem.x1            = x1_p[k];
      problem.x2            = x2_p[k];
      problem.log_t1        = std::log(t1_p[k]);
      problem.log_t2        = std::log(t2_p[k]);
      problem.mu            = mu_p[cell];
      problem.tau           = tau_p[cell];
      problem.weight        = weight_p == NULL ? 1.0 : weight_p[k];
      problem.mean          = mean_value;
      problem.sd            = sd_value;
      problem.variance      = sd_value * sd_value;
      problem.log_factorial = log_factorial;

      const IntegralResult result = integrate_problem(
        problem,
        rules,
        tolerance_value,
        consecutive_value,
        mode_tolerance_value
      );
      if (!result.mode_converged) {
        Rf_error("Poisson AGHQ mode failed at sample %d, observation %d.",
                 s + 1, k + 1);
      }
      if (!result.quadrature_converged) {
        Rf_error(
          "Poisson AGHQ failed to converge at sample %d, observation %d "
          "(order %d, change %.8g).",
          s + 1, k + 1, result.order, result.change
        );
      }

      if (row_sum) {
        value_p[s] += result.value;
      } else {
        value_p[cell] = result.value;
      }
      update_diagnostics(result, &diagnostics);
    }
  }

  SEXP out = PROTECT(make_result(value, diagnostics, rules));
  UNPROTECT(2);
  return out;
}

} // namespace


extern "C" SEXP RoBMA_glmm_binom_aghq(
  SEXP ai, SEXP ci, SEXP n1i, SEXP n2i, SEXP mu_samples,
  SEXP tau_within, SEXP weights, SEXP alpha, SEXP beta,
  SEXP nodes, SEXP log_weights, SEXP tolerance,
  SEXP consecutive, SEXP mode_tolerance)
{
  return run_binomial(
    ai, ci, n1i, n2i, mu_samples, tau_within, weights,
    alpha, beta, nodes, log_weights, tolerance, consecutive,
    mode_tolerance, false
  );
}


extern "C" SEXP RoBMA_glmm_binom_aghq_row_sum(
  SEXP ai, SEXP ci, SEXP n1i, SEXP n2i, SEXP mu_samples,
  SEXP tau_within, SEXP weights, SEXP alpha, SEXP beta,
  SEXP nodes, SEXP log_weights, SEXP tolerance,
  SEXP consecutive, SEXP mode_tolerance)
{
  return run_binomial(
    ai, ci, n1i, n2i, mu_samples, tau_within, weights,
    alpha, beta, nodes, log_weights, tolerance, consecutive,
    mode_tolerance, true
  );
}


extern "C" SEXP RoBMA_glmm_pois_aghq(
  SEXP x1i, SEXP x2i, SEXP t1i, SEXP t2i, SEXP mu_samples,
  SEXP tau_within, SEXP weights, SEXP mean, SEXP sd,
  SEXP nodes, SEXP log_weights, SEXP tolerance,
  SEXP consecutive, SEXP mode_tolerance)
{
  return run_poisson(
    x1i, x2i, t1i, t2i, mu_samples, tau_within, weights,
    mean, sd, nodes, log_weights, tolerance, consecutive,
    mode_tolerance, false
  );
}


extern "C" SEXP RoBMA_glmm_pois_aghq_row_sum(
  SEXP x1i, SEXP x2i, SEXP t1i, SEXP t2i, SEXP mu_samples,
  SEXP tau_within, SEXP weights, SEXP mean, SEXP sd,
  SEXP nodes, SEXP log_weights, SEXP tolerance,
  SEXP consecutive, SEXP mode_tolerance)
{
  return run_poisson(
    x1i, x2i, t1i, t2i, mu_samples, tau_within, weights,
    mean, sd, nodes, log_weights, tolerance, consecutive,
    mode_tolerance, true
  );
}
