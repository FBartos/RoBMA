#include "selnorm-mv.h"

#include <R_ext/Lapack.h>
#include <JRmath.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "selnorm.h"

#ifndef FCONE
# define FCONE
#endif

namespace {

double log_add_exp(double x, double y)
{
  if (x == -std::numeric_limits<double>::infinity()) return y;
  if (y == -std::numeric_limits<double>::infinity()) return x;
  const double maximum = std::max(x, y);
  return maximum + std::log(std::exp(x - maximum) + std::exp(y - maximum));
}

double qmc_value(const double *values, int scrambles, int points,
                 int dimensions, int scramble, int point, int dimension)
{
  const std::size_t index = static_cast<std::size_t>(scramble) +
    static_cast<std::size_t>(scrambles) *
    (static_cast<std::size_t>(point) +
     static_cast<std::size_t>(points) * static_cast<std::size_t>(dimension));
  const std::size_t total = static_cast<std::size_t>(scrambles) *
    static_cast<std::size_t>(points) * static_cast<std::size_t>(dimensions);
  if (index >= total) return std::numeric_limits<double>::quiet_NaN();
  return values[index];
}

struct ClusterNormalizerContext {
  const double *mean;
  const double *residual_sd;
  const double *loading;
  int dimension;
  const double *selection_se;
  const double *omega;
  int kernel_mode;
  SelNormKernelData selection;
};

double cluster_log_integrand(const ClusterNormalizerContext &context,
                             double gamma)
{
  double value = -0.5 * gamma * gamma;
  for (int i = 0; i < context.dimension; ++i) {
    const double conditional_mean =
      context.mean[i] + context.loading[i] * gamma;
    const double local = cpp_selnorm_step_log_norm(
      conditional_mean, context.residual_sd[i], context.selection_se[i],
      context.omega, context.selection, 1, false
    );
    if (!std::isfinite(local)) return -std::numeric_limits<double>::infinity();
    value += local;
  }
  return value;
}

double cluster_rule_log_integral(const ClusterNormalizerContext &context,
                                 const double *nodes,
                                 const double *log_weights,
                                 int offset, int order, double mode,
                                 double curvature)
{
  const double negative_infinity = -std::numeric_limits<double>::infinity();
  const double inverse_scale = 1.0 / std::sqrt(curvature);
  const double log_jacobian = -0.5 * std::log(curvature);
  double out = negative_infinity;
  for (int j = 0; j < order; ++j) {
    const int index = offset + j;
    const double z = nodes[index];
    const double gamma = mode + z * inverse_scale;
    const double term = log_weights[index] +
      cluster_log_integrand(context, gamma) + 0.5 * z * z + log_jacobian;
    out = log_add_exp(out, term);
  }
  return out;
}

double cluster_relative_change(double coarse, double fine)
{
  const double difference = coarse - fine;
  if (!std::isfinite(difference) || std::fabs(difference) > 50.0) {
    return std::numeric_limits<double>::infinity();
  }
  return std::fabs(std::expm1(difference));
}

struct FactorNormalizerContext {
  const double *mean;
  const double *residual_sd;
  const double *loading;
  int dimension;
  int rank;
  const double *selection_se;
  const double *omega;
  int n_bins;
  int kernel_mode;
  SelNormKernelData selection;
};

double factor_log_integrand(const FactorNormalizerContext &context,
                            const double *latent, double *gradient)
{
  double value = 0.0;
  for (int factor = 0; factor < context.rank; ++factor) {
    value -= 0.5 * latent[factor] * latent[factor];
    if (gradient != nullptr) gradient[factor] = -latent[factor];
  }

  for (int i = 0; i < context.dimension; ++i) {
    double conditional_mean = context.mean[i];
    for (int factor = 0; factor < context.rank; ++factor) {
      conditional_mean +=
        context.loading[i + context.dimension * factor] * latent[factor];
    }
    const double local = cpp_selnorm_step_log_norm(
      conditional_mean, context.residual_sd[i], context.selection_se[i],
      context.omega, context.selection, 1, false
    );
    if (!std::isfinite(local)) {
      return -std::numeric_limits<double>::infinity();
    }
    value += local;

    if (gradient != nullptr) {
      double derivative_numerator = 0.0;
      const double signed_mean = context.selection.effect_sign *
        conditional_mean;
      for (int bin = 0; bin < context.n_bins; ++bin) {
        const double lower = context.selection.z_lower[bin] *
          context.selection_se[i];
        const double upper = context.selection.z_upper[bin] *
          context.selection_se[i];
        const double lower_score =
          (lower - signed_mean) / context.residual_sd[i];
        const double upper_score =
          (upper - signed_mean) / context.residual_sd[i];
        const double lower_density = std::isfinite(lower_score) ?
          dnorm(lower_score, 0.0, 1.0, false) : 0.0;
        const double upper_density = std::isfinite(upper_score) ?
          dnorm(upper_score, 0.0, 1.0, false) : 0.0;
        derivative_numerator += context.omega[bin] *
          (lower_density - upper_density);
      }
      const double local_probability = std::exp(local);
      if (!(local_probability > 0.0) ||
          !std::isfinite(derivative_numerator)) {
        return -std::numeric_limits<double>::infinity();
      }
      const double derivative = context.selection.effect_sign *
        derivative_numerator /
        (context.residual_sd[i] * local_probability);
      if (!std::isfinite(derivative)) {
        return -std::numeric_limits<double>::infinity();
      }
      for (int factor = 0; factor < context.rank; ++factor) {
        gradient[factor] +=
          context.loading[i + context.dimension * factor] * derivative;
      }
    }
  }
  return value;
}

bool factor_nested_rule_log_integral(
    const FactorNormalizerContext &context, const double *nodes,
    const double *log_weights, int offset, int order, double *result)
{
  std::vector<int> support_size(static_cast<std::size_t>(context.rank), 0);
  std::vector<int> permutation(static_cast<std::size_t>(context.rank));
  for (int factor = 0; factor < context.rank; ++factor) {
    permutation[factor] = factor;
    for (int i = 0; i < context.dimension; ++i) {
      if (context.loading[i + context.dimension * factor] != 0.0) {
        ++support_size[factor];
      }
    }
  }
  std::stable_sort(
    permutation.begin(), permutation.end(),
    [&support_size](int left, int right) {
      return support_size[left] > support_size[right];
    }
  );

  int effective_rank = context.rank;
  while (effective_rank > 0 &&
         support_size[permutation[effective_rank - 1]] == 0) {
    --effective_rank;
  }
  for (int position = 1; position < effective_rank; ++position) {
    const int previous = permutation[position - 1];
    const int current  = permutation[position];
    for (int i = 0; i < context.dimension; ++i) {
      if (context.loading[i + context.dimension * current] != 0.0 &&
          context.loading[i + context.dimension * previous] == 0.0) {
        return false;
      }
    }
  }

  std::vector<std::vector<int>> row_groups(
    static_cast<std::size_t>(effective_rank + 1)
  );
  for (int i = 0; i < context.dimension; ++i) {
    int active = 0;
    while (active < effective_rank &&
           context.loading[
             i + context.dimension * permutation[active]
           ] != 0.0) {
      ++active;
    }
    for (int position = active; position < effective_rank; ++position) {
      if (context.loading[
            i + context.dimension * permutation[position]
          ] != 0.0) {
        return false;
      }
    }
    row_groups[active].push_back(i);
  }
  if (effective_rank == context.rank &&
      row_groups[effective_rank].size() ==
        static_cast<std::size_t>(context.dimension)) {
    return false;
  }

  std::vector<std::size_t> sizes(
    static_cast<std::size_t>(effective_rank + 1), 1
  );
  for (int position = 1; position <= effective_rank; ++position) {
    sizes[position] = sizes[position - 1] * static_cast<std::size_t>(order);
  }
  std::vector<std::vector<long double>> factors(
    static_cast<std::size_t>(effective_rank + 1)
  );
  std::vector<double> latent(static_cast<std::size_t>(effective_rank));
  for (int group = 0; group <= effective_rank; ++group) {
    factors[group].assign(sizes[group], 1.0L);
    for (int row : row_groups[group]) {
      for (std::size_t point = 0; point < sizes[group]; ++point) {
        std::size_t remaining = point;
        double conditional_mean = context.mean[row];
        for (int position = 0; position < group; ++position) {
          const int node = static_cast<int>(remaining % order);
          remaining /= static_cast<std::size_t>(order);
          latent[position] = nodes[offset + node];
          conditional_mean += context.loading[
            row + context.dimension * permutation[position]
          ] * latent[position];
        }
        double omega_last = 0.0;
        double normalizer = 0.0;
        const bool ok = cpp_selnorm_step_cdf_telescope_plan(
          conditional_mean, context.residual_sd[row],
          context.selection_se[row], context.omega, context.selection,
          nullptr, nullptr, &omega_last, &normalizer, 1, false
        );
        if (!ok) return false;
        factors[group][point] *= static_cast<long double>(normalizer);
        if (!(factors[group][point] > 0.0L) ||
            !std::isfinite(factors[group][point])) {
          return false;
        }
      }
    }
  }

  std::vector<long double> weights(static_cast<std::size_t>(order));
  for (int node = 0; node < order; ++node) {
    weights[node] =
      std::exp(static_cast<long double>(log_weights[offset + node]));
  }
  std::vector<long double> current = std::move(factors[effective_rank]);
  for (int group = effective_rank; group > 0; --group) {
    const std::size_t lower_size = sizes[group - 1];
    std::vector<long double> integrated(lower_size, 0.0L);
    for (std::size_t lower = 0; lower < lower_size; ++lower) {
      for (int node = 0; node < order; ++node) {
        integrated[lower] += weights[node] * current[
          lower + static_cast<std::size_t>(node) * lower_size
        ];
      }
      integrated[lower] *= factors[group - 1][lower];
      if (!(integrated[lower] > 0.0L) ||
          !std::isfinite(integrated[lower])) {
        return false;
      }
    }
    current = std::move(integrated);
  }
  if (current.size() != 1 || !(current[0] > 0.0L) ||
      !std::isfinite(current[0])) {
    return false;
  }
  *result = static_cast<double>(std::log(current[0]));
  return std::isfinite(*result);
}

double factor_rule_log_integral(const FactorNormalizerContext &context,
                                const double *nodes,
                                const double *log_weights,
                                int offset, int order)
{
  double nested = 0.0;
  if (factor_nested_rule_log_integral(
        context, nodes, log_weights, offset, order, &nested
      )) {
    return nested;
  }
  std::size_t total = 1;
  for (int factor = 0; factor < context.rank; ++factor) {
    total *= static_cast<std::size_t>(order);
  }

  // Factor quadrature only reaches this kernel for an ordered step partition.
  // Accumulating its positive normalizers avoids a log/exp pair for every row
  // and quadrature point. Retain the log-scale path below for rare numerical
  // fallback cases.
  std::vector<long double> weights(static_cast<std::size_t>(order));
  for (int j = 0; j < order; ++j) {
    weights[static_cast<std::size_t>(j)] =
      std::exp(static_cast<long double>(log_weights[offset + j]));
  }
  long double direct_sum = 0.0L;
  bool direct_ok = true;
  std::vector<double> latent(static_cast<std::size_t>(context.rank));
  for (std::size_t point = 0; point < total && direct_ok; ++point) {
    std::size_t remaining = point;
    long double term = 1.0L;
    for (int factor = 0; factor < context.rank; ++factor) {
      const int index = static_cast<int>(remaining % order);
      remaining /= static_cast<std::size_t>(order);
      latent[factor] = nodes[offset + index];
      term *= weights[static_cast<std::size_t>(index)];
    }
    for (int i = 0; i < context.dimension; ++i) {
      double conditional_mean = context.mean[i];
      for (int factor = 0; factor < context.rank; ++factor) {
        conditional_mean +=
          context.loading[i + context.dimension * factor] * latent[factor];
      }
      double omega_last = 0.0;
      double normalizer = 0.0;
      direct_ok = cpp_selnorm_step_cdf_telescope_plan(
        conditional_mean, context.residual_sd[i], context.selection_se[i],
        context.omega, context.selection, nullptr, nullptr, &omega_last,
        &normalizer, 1, false
      );
      if (!direct_ok) break;
      term *= static_cast<long double>(normalizer);
    }
    direct_sum += term;
    direct_ok = direct_ok && std::isfinite(direct_sum) && direct_sum > 0.0L;
  }
  if (direct_ok) return static_cast<double>(std::log(direct_sum));

  double out = -std::numeric_limits<double>::infinity();
  for (std::size_t point = 0; point < total; ++point) {
    std::size_t remaining = point;
    double log_weight = 0.0;
    double normal_kernel = 0.0;
    for (int factor = 0; factor < context.rank; ++factor) {
      const int index = offset + static_cast<int>(remaining % order);
      remaining /= static_cast<std::size_t>(order);
      latent[factor] = nodes[index];
      normal_kernel += 0.5 * latent[factor] * latent[factor];
      log_weight += log_weights[index];
    }
    const double log_integrand = factor_log_integrand(
      context, latent.data(), nullptr
    );
    if (std::isfinite(log_integrand)) {
      out = log_add_exp(out, log_weight + log_integrand + normal_kernel);
    }
  }
  return out;
}

double vector_dot(const std::vector<double> &x,
                  const std::vector<double> &y)
{
  double out = 0.0;
  for (std::size_t i = 0; i < x.size(); ++i) out += x[i] * y[i];
  return out;
}

void factor_optimize_mode(const FactorNormalizerContext &context,
                          std::vector<double> *position)
{
  const int rank = context.rank;
  std::vector<double> gradient(static_cast<std::size_t>(rank));
  double value = factor_log_integrand(
    context, position->data(), gradient.data()
  );
  if (!std::isfinite(value)) return;

  std::vector<double> inverse_hessian(
    static_cast<std::size_t>(rank * rank), 0.0
  );
  for (int factor = 0; factor < rank; ++factor) {
    inverse_hessian[factor + rank * factor] = 1.0;
  }
  std::vector<double> best = *position;
  double best_value = value;
  const double gradient_tolerance =
    std::sqrt(std::numeric_limits<double>::epsilon());

  for (int iteration = 0; iteration < 100; ++iteration) {
    const double gradient_norm = std::sqrt(vector_dot(gradient, gradient));
    double position_norm = 0.0;
    for (int factor = 0; factor < rank; ++factor) {
      position_norm += (*position)[factor] * (*position)[factor];
    }
    position_norm = std::sqrt(position_norm);
    if (gradient_norm <= gradient_tolerance * (1.0 + position_norm)) break;

    std::vector<double> direction(static_cast<std::size_t>(rank), 0.0);
    for (int column = 0; column < rank; ++column) {
      for (int row = 0; row < rank; ++row) {
        direction[row] += inverse_hessian[row + rank * column] *
          gradient[column];
      }
    }
    double directional_derivative = vector_dot(gradient, direction);
    if (!(directional_derivative > 0.0) ||
        !std::isfinite(directional_derivative)) {
      direction = gradient;
      directional_derivative = vector_dot(gradient, direction);
    }

    std::vector<double> candidate(static_cast<std::size_t>(rank));
    std::vector<double> candidate_gradient(static_cast<std::size_t>(rank));
    double candidate_value = -std::numeric_limits<double>::infinity();
    double step = 1.0;
    bool accepted = false;
    for (int line_search = 0; line_search < 30; ++line_search) {
      for (int factor = 0; factor < rank; ++factor) {
        candidate[factor] = (*position)[factor] + step * direction[factor];
      }
      candidate_value = factor_log_integrand(
        context, candidate.data(), candidate_gradient.data()
      );
      if (std::isfinite(candidate_value) &&
          candidate_value >= value + 1e-4 * step * directional_derivative) {
        accepted = true;
        break;
      }
      step *= 0.5;
    }
    if (!accepted) break;

    std::vector<double> displacement(static_cast<std::size_t>(rank));
    std::vector<double> curvature(static_cast<std::size_t>(rank));
    for (int factor = 0; factor < rank; ++factor) {
      displacement[factor] = candidate[factor] - (*position)[factor];
      curvature[factor] = gradient[factor] - candidate_gradient[factor];
    }
    const double displacement_curvature =
      vector_dot(displacement, curvature);
    if (displacement_curvature > 0.0 &&
        std::isfinite(displacement_curvature)) {
      std::vector<double> h_curvature(static_cast<std::size_t>(rank), 0.0);
      for (int column = 0; column < rank; ++column) {
        for (int row = 0; row < rank; ++row) {
          h_curvature[row] += inverse_hessian[row + rank * column] *
            curvature[column];
        }
      }
      const double curvature_h_curvature =
        vector_dot(curvature, h_curvature);
      const double rho = 1.0 / displacement_curvature;
      const double scale =
        (1.0 + curvature_h_curvature * rho) * rho;
      for (int column = 0; column < rank; ++column) {
        for (int row = 0; row < rank; ++row) {
          inverse_hessian[row + rank * column] +=
            scale * displacement[row] * displacement[column] -
            rho * (displacement[row] * h_curvature[column] +
                   h_curvature[row] * displacement[column]);
        }
      }
    } else {
      std::fill(inverse_hessian.begin(), inverse_hessian.end(), 0.0);
      for (int factor = 0; factor < rank; ++factor) {
        inverse_hessian[factor + rank * factor] = 1.0;
      }
    }

    *position = candidate;
    gradient = candidate_gradient;
    value = candidate_value;
    if (value > best_value) {
      best = *position;
      best_value = value;
    }
  }
  *position = best;
}

bool factor_solve_cholesky(const std::vector<double> &cholesky, int rank,
                           std::vector<double> *value)
{
  for (int row = 0; row < rank; ++row) {
    for (int column = 0; column < row; ++column) {
      (*value)[row] -= cholesky[row + rank * column] * (*value)[column];
    }
    const double diagonal = cholesky[row + rank * row];
    if (!(diagonal > 0.0) || !std::isfinite(diagonal)) return false;
    (*value)[row] /= diagonal;
  }
  for (int row = rank - 1; row >= 0; --row) {
    for (int column = row + 1; column < rank; ++column) {
      (*value)[row] -= cholesky[column + rank * row] * (*value)[column];
    }
    (*value)[row] /= cholesky[row + rank * row];
  }
  return true;
}

double factor_log_mean(const double *log_values, int count)
{
  double out = -std::numeric_limits<double>::infinity();
  for (int i = 0; i < count; ++i) out = log_add_exp(out, log_values[i]);
  return out - std::log(static_cast<double>(count));
}

double factor_log_mean(const std::vector<double> &log_values, int count)
{
  return factor_log_mean(log_values.data(), count);
}

}

double cpp_selnorm_factor_step_lpdf(
    const double *x, const double *mean, const double *residual_sd,
    const double *loading, int dimension, int rank,
    const double *selection_se, const double *omega, int n_bins,
    const double *z_lower, const double *z_upper, const int *obs_bin,
    int effect_sign, bool telescope_probabilities, int kernel_mode,
    const double *quadrature_nodes, const double *quadrature_log_weights,
    const double *quadrature_orders, int quadrature_rule_count,
    const double *qmc, int initial_points, int max_points, int scrambles,
    double relative_tolerance, double *relative_mcse,
    double *relative_change)
{
  const double negative_infinity = -std::numeric_limits<double>::infinity();
  const double log_two_pi = std::log(6.283185307179586476925286766559);
  *relative_mcse = 0.0;
  *relative_change = 0.0;
  if (dimension < 1 || rank < 1 || quadrature_rule_count < 3 ||
      initial_points < 2 ||
      max_points < initial_points || scrambles < 2 ||
      !(relative_tolerance > 0.0)) {
    return negative_infinity;
  }

  double phack_z_zero[2] = {0, 0};
  double segment_bounds_zero[1] = {0};
  int segment_zero[1] = {0};
  SelNormKernelData selection;
  selection.n_bins = n_bins;
  selection.n_segments = 0;
  selection.effect_sign = effect_sign;
  selection.q = 0;
  selection.z_lower = z_lower;
  selection.z_upper = z_upper;
  selection.phack_z_source = phack_z_zero;
  selection.phack_z_dest = phack_z_zero;
  selection.segment_bounds = segment_bounds_zero;
  selection.segment_step_bin = segment_zero;
  selection.segment_phack_region = segment_zero;
  selection.segment_step_bin_real = 0;
  selection.segment_phack_region_real = 0;
  selection.trusted_step_partition = true;
  selection.telescope_probabilities = telescope_probabilities;

  std::vector<double> inverse_variance(static_cast<std::size_t>(dimension));
  std::vector<double> system(static_cast<std::size_t>(rank * rank), 0.0);
  std::vector<double> score(static_cast<std::size_t>(rank), 0.0);
  double log_det = 0.0;
  double diagonal_quadratic = 0.0;
  for (int i = 0; i < dimension; ++i) {
    if (!(residual_sd[i] > 0.0) || !std::isfinite(residual_sd[i]) ||
        !std::isfinite(x[i]) || !std::isfinite(mean[i])) {
      return negative_infinity;
    }
    const double inverse = 1.0 / (residual_sd[i] * residual_sd[i]);
    inverse_variance[i] = inverse;
    log_det += 2.0 * std::log(residual_sd[i]);
    const double residual = x[i] - mean[i];
    diagonal_quadratic += residual * residual * inverse;
    for (int factor = 0; factor < rank; ++factor) {
      const double coefficient = loading[i + dimension * factor];
      if (!std::isfinite(coefficient)) return negative_infinity;
      score[factor] += coefficient * residual * inverse;
      for (int other = 0; other <= factor; ++other) {
        system[factor + rank * other] += coefficient *
          loading[i + dimension * other] * inverse;
      }
    }
  }
  for (int factor = 0; factor < rank; ++factor) {
    system[factor + rank * factor] += 1.0;
    for (int other = 0; other < factor; ++other) {
      system[other + rank * factor] = system[factor + rank * other];
    }
  }
  std::vector<double> cholesky = system;
  int info = 0;
  F77_CALL(dpotrf)("L", &rank, cholesky.data(), &rank, &info FCONE);
  if (info != 0) return negative_infinity;
  for (int factor = 0; factor < rank; ++factor) {
    const double diagonal = cholesky[factor + rank * factor];
    if (!(diagonal > 0.0) || !std::isfinite(diagonal)) {
      return negative_infinity;
    }
    log_det += 2.0 * std::log(diagonal);
  }
  std::vector<double> solved_score = score;
  for (int row = 0; row < rank; ++row) {
    for (int column = 0; column < row; ++column) {
      solved_score[row] -= cholesky[row + rank * column] *
        solved_score[column];
    }
    solved_score[row] /= cholesky[row + rank * row];
  }
  double correction = vector_dot(solved_score, solved_score);
  double log_density = -0.5 * (
    static_cast<double>(dimension) * log_two_pi + log_det +
    diagonal_quadratic - correction
  );
  if (kernel_mode == SELKERNEL_NORMAL) return log_density;
  for (int i = 0; i < dimension; ++i) {
    const double log_weight = std::log(omega[obs_bin[i] - 1]);
    if (!std::isfinite(log_weight)) return negative_infinity;
    log_density += log_weight;
  }

  FactorNormalizerContext context;
  context.mean = mean;
  context.residual_sd = residual_sd;
  context.loading = loading;
  context.dimension = dimension;
  context.rank = rank;
  context.selection_se = selection_se;
  context.omega = omega;
  context.n_bins = n_bins;
  context.kernel_mode = kernel_mode;
  context.selection = selection;

  int quadrature_offset = 0;
  int quadrature_order = static_cast<int>(quadrature_orders[0]);
  double previous_quadrature = factor_rule_log_integral(
    context, quadrature_nodes, quadrature_log_weights,
    quadrature_offset, quadrature_order
  );
  double previous_quadrature_change =
    std::numeric_limits<double>::infinity();
  quadrature_offset += quadrature_order;
  for (int rule = 1; rule < quadrature_rule_count; ++rule) {
    quadrature_order = static_cast<int>(quadrature_orders[rule]);
    const double current_quadrature = factor_rule_log_integral(
      context, quadrature_nodes, quadrature_log_weights,
      quadrature_offset, quadrature_order
    );
    const double current_quadrature_change = cluster_relative_change(
      previous_quadrature, current_quadrature
    );
    if (std::isfinite(current_quadrature) &&
        current_quadrature_change <= relative_tolerance &&
        previous_quadrature_change <= relative_tolerance) {
      *relative_change = current_quadrature_change;
      return log_density - current_quadrature;
    }
    previous_quadrature = current_quadrature;
    previous_quadrature_change = current_quadrature_change;
    quadrature_offset += quadrature_order;
  }

  std::vector<std::vector<double> > starts;
  starts.push_back(std::vector<double>(static_cast<std::size_t>(rank), 0.0));
  std::vector<double> boundaries;
  for (int bin = 0; bin < n_bins; ++bin) {
    if (std::isfinite(z_lower[bin])) boundaries.push_back(z_lower[bin]);
    if (std::isfinite(z_upper[bin])) boundaries.push_back(z_upper[bin]);
  }
  std::sort(boundaries.begin(), boundaries.end());
  boundaries.erase(
    std::unique(boundaries.begin(), boundaries.end()), boundaries.end()
  );
  for (std::size_t boundary = 0; boundary < boundaries.size(); ++boundary) {
    std::vector<double> start(static_cast<std::size_t>(rank), 0.0);
    for (int i = 0; i < dimension; ++i) {
      const double target = effect_sign * boundaries[boundary] *
        selection_se[i] - mean[i];
      for (int factor = 0; factor < rank; ++factor) {
        start[factor] += loading[i + dimension * factor] *
          inverse_variance[i] * target;
      }
    }
    if (factor_solve_cholesky(cholesky, rank, &start)) {
      starts.push_back(start);
    }
  }

  std::vector<double> proposal = starts[0];
  double proposal_value = negative_infinity;
  for (std::size_t candidate = 0; candidate < starts.size(); ++candidate) {
    factor_optimize_mode(context, &starts[candidate]);
    const double value = factor_log_integrand(
      context, starts[candidate].data(), nullptr
    );
    if (std::isfinite(value) && value > proposal_value) {
      proposal = starts[candidate];
      proposal_value = value;
    }
  }
  if (!std::isfinite(proposal_value)) return negative_infinity;

  std::vector<double> scramble_log_mean(static_cast<std::size_t>(scrambles));
  std::vector<double> scramble_coarse_log_mean(
    static_cast<std::size_t>(scrambles)
  );
  const int proposal_count = 2;
  const int qmc_dimensions = proposal_count * rank;
  const double log_proposal_count = std::log(
    static_cast<double>(proposal_count)
  );
  std::vector<double> point_log_value(
    static_cast<std::size_t>(proposal_count * max_points * scrambles)
  );
  std::vector<double> latent(static_cast<std::size_t>(rank));
  std::vector<double> standard_normal(static_cast<std::size_t>(rank));
  std::vector<double> centered(static_cast<std::size_t>(rank));
  int evaluated_points = 0;
  int current_points = initial_points;
  int comparison_points = std::max(2, initial_points / 2);
  double log_normalizer = negative_infinity;
  while (true) {
    for (int scramble = 0; scramble < scrambles; ++scramble) {
      const std::size_t scramble_offset = static_cast<std::size_t>(
        proposal_count * max_points * scramble
      );
      for (int component = 0; component < proposal_count; ++component) {
        for (int point = evaluated_points; point < current_points; ++point) {
          for (int factor = 0; factor < rank; ++factor) {
            standard_normal[factor] = qnorm(
              qmc_value(
                qmc, scrambles, max_points, qmc_dimensions, scramble, point,
                factor + rank * component
              ),
              0.0, 1.0, true, false
            );
            if (!std::isfinite(standard_normal[factor])) {
              return negative_infinity;
            }
          }
          for (int factor = 0; factor < rank; ++factor) {
            latent[factor] = standard_normal[factor];
            if (component == 1) {
              latent[factor] += proposal[factor];
            }
          }
          double prior_log_density = 0.0;
          for (int factor = 0; factor < rank; ++factor) {
            prior_log_density -= 0.5 * latent[factor] * latent[factor];
            centered[factor] = latent[factor] - proposal[factor];
          }
          const double mode_log_density =
            -0.5 * vector_dot(centered, centered);
          const double target = factor_log_integrand(
            context, latent.data(), nullptr
          );
          if (!std::isfinite(target)) return negative_infinity;
          const double mixture_log_density = log_add_exp(
            prior_log_density, mode_log_density
          ) - log_proposal_count;
          point_log_value[
            scramble_offset + component + proposal_count * point
          ] = target - mixture_log_density;
        }
      }
      scramble_log_mean[scramble] = factor_log_mean(
        point_log_value.data() + scramble_offset,
        proposal_count * current_points
      );
      scramble_coarse_log_mean[scramble] = factor_log_mean(
        point_log_value.data() + scramble_offset,
        proposal_count * comparison_points
      );
    }
    log_normalizer = factor_log_mean(scramble_log_mean, scrambles);
    const double coarse_log_normalizer = factor_log_mean(
      scramble_coarse_log_mean, scrambles
    );
    if (!std::isfinite(log_normalizer) ||
        !std::isfinite(coarse_log_normalizer)) return negative_infinity;
    *relative_change = cluster_relative_change(
      coarse_log_normalizer, log_normalizer
    );

    double squared_relative = 0.0;
    for (int scramble = 0; scramble < scrambles; ++scramble) {
      const double ratio = std::exp(
        scramble_log_mean[scramble] - log_normalizer
      );
      const double difference = ratio - 1.0;
      squared_relative += difference * difference;
    }
    *relative_mcse = std::sqrt(
      squared_relative /
      (static_cast<double>(scrambles) * static_cast<double>(scrambles - 1))
    );
    if (std::max(*relative_mcse, *relative_change) <= relative_tolerance ||
        current_points == max_points) {
      break;
    }
    evaluated_points = current_points;
    comparison_points = current_points;
    current_points = std::min(max_points, 2 * current_points);
  }
  return log_density - log_normalizer;
}

double cpp_selnorm_cluster_step_lpdf(
    const double *x, const double *mean, const double *residual_sd,
    const double *loading, int dimension, const double *selection_se,
    const double *omega, int n_bins, const double *z_lower,
    const double *z_upper, const int *obs_bin, int effect_sign,
    bool telescope_probabilities, int kernel_mode,
    const double *quadrature_nodes, const double *quadrature_log_weights,
    const double *quadrature_orders, int quadrature_rule_count,
    double relative_tolerance, double *relative_change)
{
  const double negative_infinity = -std::numeric_limits<double>::infinity();
  const double log_two_pi = std::log(6.283185307179586476925286766559);
  *relative_change = 0.0;

  double phack_z_zero[2] = {0, 0};
  double segment_bounds_zero[1] = {0};
  int segment_zero[1] = {0};
  SelNormKernelData selection;
  selection.n_bins = n_bins;
  selection.n_segments = 0;
  selection.effect_sign = effect_sign;
  selection.q = 0;
  selection.z_lower = z_lower;
  selection.z_upper = z_upper;
  selection.phack_z_source = phack_z_zero;
  selection.phack_z_dest = phack_z_zero;
  selection.segment_bounds = segment_bounds_zero;
  selection.segment_step_bin = segment_zero;
  selection.segment_phack_region = segment_zero;
  selection.segment_step_bin_real = 0;
  selection.segment_phack_region_real = 0;
  selection.trusted_step_partition = true;
  selection.telescope_probabilities = telescope_probabilities;

  double log_det = 0.0;
  double denominator = 1.0;
  double score = 0.0;
  double diagonal_quadratic = 0.0;
  for (int i = 0; i < dimension; ++i) {
    if (!(residual_sd[i] > 0.0) || !std::isfinite(residual_sd[i]) ||
        !std::isfinite(loading[i]) || !std::isfinite(x[i]) ||
        !std::isfinite(mean[i])) return negative_infinity;
    const double variance = residual_sd[i] * residual_sd[i];
    const double residual = x[i] - mean[i];
    log_det += std::log(variance);
    denominator += loading[i] * loading[i] / variance;
    score += loading[i] * residual / variance;
    diagonal_quadratic += residual * residual / variance;
  }
  const double quadratic = diagonal_quadratic - score * score / denominator;
  double log_density = -0.5 * (
    static_cast<double>(dimension) * log_two_pi + log_det +
    std::log(denominator) + quadratic
  );
  if (kernel_mode == SELKERNEL_NORMAL) return log_density;

  for (int i = 0; i < dimension; ++i) {
    const double log_weight = std::log(omega[obs_bin[i] - 1]);
    if (!std::isfinite(log_weight)) return negative_infinity;
    log_density += log_weight;
  }

  ClusterNormalizerContext context;
  context.mean = mean;
  context.residual_sd = residual_sd;
  context.loading = loading;
  context.dimension = dimension;
  context.selection_se = selection_se;
  context.omega = omega;
  context.kernel_mode = kernel_mode;
  context.selection = selection;

  int offset = 0;
  int order = static_cast<int>(quadrature_orders[0]);
  double previous = cluster_rule_log_integral(
    context, quadrature_nodes, quadrature_log_weights,
    offset, order, 0.0, 1.0
  );
  offset += order;
  for (int rule = 1; rule < quadrature_rule_count; ++rule) {
    order = static_cast<int>(quadrature_orders[rule]);
    const double current = cluster_rule_log_integral(
      context, quadrature_nodes, quadrature_log_weights,
      offset, order, 0.0, 1.0
    );
    *relative_change = cluster_relative_change(previous, current);
    if (std::isfinite(current) && *relative_change <= relative_tolerance) {
      return log_density - current;
    }
    previous = current;
    offset += order;
  }

  return negative_infinity;
}

double cpp_selnorm_mnorm_step_lpdf(
    const double *x, const double *mean, const double *covariance_lower,
    int dimension, const double *selection_se, const double *omega,
    int n_bins, const double *z_lower, const double *z_upper,
    const int *obs_bin, int effect_sign, bool telescope_probabilities,
    int kernel_mode, const double *qmc, int points, int scrambles,
    double *relative_mcse)
{
  const int k = dimension;
  const int dimensions = 2 * k;
  const double two_pi = 6.283185307179586476925286766559;
  const double negative_infinity = -std::numeric_limits<double>::infinity();
  *relative_mcse = 0.0;

  double phack_z_zero[2] = {0, 0};
  double segment_bounds_zero[1] = {0};
  int segment_zero[1] = {0};
  SelNormKernelData selection;
  selection.n_bins = n_bins;
  selection.n_segments = 0;
  selection.effect_sign = effect_sign;
  selection.q = 0;
  selection.z_lower = z_lower;
  selection.z_upper = z_upper;
  selection.phack_z_source = phack_z_zero;
  selection.phack_z_dest = phack_z_zero;
  selection.segment_bounds = segment_bounds_zero;
  selection.segment_step_bin = segment_zero;
  selection.segment_phack_region = segment_zero;
  selection.segment_step_bin_real = 0;
  selection.segment_phack_region_real = 0;
  selection.trusted_step_partition = true;
  selection.telescope_probabilities = telescope_probabilities;

  if (k == 1) {
    const double variance = covariance_lower[0];
    if (!(variance > 0.0) || !std::isfinite(variance)) {
      return negative_infinity;
    }
    return cpp_selnorm_kernel_lpdf(
      x[0], mean[0], std::sqrt(variance), mean[0], std::sqrt(variance),
      selection_se[0], 1.0, omega, obs_bin[0], 0.0, 0, kernel_mode,
      selection, 1, false
    );
  }

  std::vector<double> covariance(static_cast<std::size_t>(k * k), 0.0);
  int position = 0;
  for (int column = 0; column < k; ++column) {
    for (int row = column; row < k; ++row) {
      const double value = covariance_lower[position++];
      covariance[static_cast<std::size_t>(row + k * column)] = value;
      covariance[static_cast<std::size_t>(column + k * row)] = value;
    }
  }
  std::vector<double> cholesky = covariance;
  int info = 0;
  F77_CALL(dpotrf)("L", &k, cholesky.data(), &k, &info FCONE);
  if (info != 0) return negative_infinity;

  double log_det = 0.0;
  std::vector<double> residual(static_cast<std::size_t>(k));
  for (int i = 0; i < k; ++i) {
    const double diagonal = cholesky[static_cast<std::size_t>(i + k * i)];
    if (!(diagonal > 0) || !std::isfinite(diagonal) ||
        !std::isfinite(x[i])) {
      return negative_infinity;
    }
    log_det += 2.0 * std::log(diagonal);
    residual[static_cast<std::size_t>(i)] = x[i] - mean[i];
  }
  for (int i = 0; i < k; ++i) {
    double value = residual[static_cast<std::size_t>(i)];
    for (int j = 0; j < i; ++j) {
      value -= cholesky[static_cast<std::size_t>(i + k * j)] *
        residual[static_cast<std::size_t>(j)];
    }
    residual[static_cast<std::size_t>(i)] =
      value / cholesky[static_cast<std::size_t>(i + k * i)];
  }
  double quadratic = 0.0;
  for (int i = 0; i < k; ++i) {
    const double value = residual[static_cast<std::size_t>(i)];
    if (!std::isfinite(value)) return negative_infinity;
    quadratic += value * value;
  }
  double log_density = -0.5 *
    (static_cast<double>(k) * std::log(two_pi) + log_det + quadratic);

  if (kernel_mode == SELKERNEL_NORMAL) return log_density;

  for (int i = 0; i < k; ++i) {
    const double observed_weight = omega[obs_bin[i] - 1];
    const double log_weight = std::log(observed_weight);
    if (!std::isfinite(log_weight)) return negative_infinity;
    log_density += log_weight;
  }

  bool diagonal_covariance = true;
  for (int column = 0; column < k && diagonal_covariance; ++column) {
    for (int row = column + 1; row < k; ++row) {
      if (covariance[static_cast<std::size_t>(row + k * column)] != 0) {
        diagonal_covariance = false;
        break;
      }
    }
  }
  if (diagonal_covariance) {
    double log_normalizer = 0.0;
    for (int i = 0; i < k; ++i) {
      log_normalizer += cpp_selnorm_kernel_log_norm(
        mean[i], std::sqrt(covariance[static_cast<std::size_t>(i + k * i)]),
        selection_se[i], omega, 0, 0, kernel_mode, selection, 1, false
      );
    }
    return log_density - log_normalizer;
  }

  std::vector<double> scramble_log_mean(static_cast<std::size_t>(scrambles));
  std::vector<double> latent(static_cast<std::size_t>(k));
  std::vector<double> mass(static_cast<std::size_t>(n_bins));
  std::vector<double> lower(static_cast<std::size_t>(n_bins));
  std::vector<double> upper(static_cast<std::size_t>(n_bins));
  for (int scramble = 0; scramble < scrambles; ++scramble) {
    double log_sum = negative_infinity;
    for (int point = 0; point < points; ++point) {
      double log_particle = 0.0;
      for (int i = 0; i < k; ++i) {
        double conditional_mean = mean[i];
        for (int j = 0; j < i; ++j) {
          conditional_mean +=
            cholesky[static_cast<std::size_t>(i + k * j)] *
            latent[static_cast<std::size_t>(j)];
        }
        const double conditional_sd =
          cholesky[static_cast<std::size_t>(i + k * i)];
        double log_local = 0.0;
        double sampled = 0.0;
        const double u_bin = qmc_value(
          qmc, scrambles, points, dimensions, scramble, point, 2 * i
        );
        const double u_interval = qmc_value(
          qmc, scrambles, points, dimensions, scramble, point, 2 * i + 1
        );
        if (selection.telescope_probabilities) {
          sampled = cpp_selnorm_step_log_norm_rng_workspace(
            conditional_mean, conditional_sd, selection_se[i], omega,
            u_bin, u_interval, selection, mass.data(), lower.data(),
            upper.data(), &log_local, 1, false
          );
        } else {
          log_local = cpp_selnorm_kernel_log_norm(
            conditional_mean, conditional_sd, selection_se[i], omega,
            0, 0, kernel_mode, selection, 1, false
          );
          sampled = cpp_selnorm_kernel_rng_workspace(
            conditional_mean, conditional_sd, selection_se[i], omega,
            u_bin, u_interval, 0, 0, kernel_mode, selection,
            mass.data(), lower.data(), upper.data(), 1, false
          );
        }
        if (!std::isfinite(log_local) || !std::isfinite(sampled)) {
          return negative_infinity;
        }
        log_particle += log_local;
        latent[static_cast<std::size_t>(i)] =
          (sampled - conditional_mean) / conditional_sd;
      }
      log_sum = log_add_exp(log_sum, log_particle);
    }
    scramble_log_mean[static_cast<std::size_t>(scramble)] =
      log_sum - std::log(static_cast<double>(points));
  }

  double log_normalizer = negative_infinity;
  for (int scramble = 0; scramble < scrambles; ++scramble) {
    log_normalizer = log_add_exp(
      log_normalizer,
      scramble_log_mean[static_cast<std::size_t>(scramble)]
    );
  }
  log_normalizer -= std::log(static_cast<double>(scrambles));
  if (!std::isfinite(log_normalizer)) return negative_infinity;

  double squared_relative = 0.0;
  for (int scramble = 0; scramble < scrambles; ++scramble) {
    const double ratio = std::exp(
      scramble_log_mean[static_cast<std::size_t>(scramble)] - log_normalizer
    );
    const double difference = ratio - 1.0;
    squared_relative += difference * difference;
  }
  *relative_mcse = std::sqrt(
    squared_relative /
    (static_cast<double>(scrambles) * static_cast<double>(scrambles - 1))
  );

  return log_density - log_normalizer;
}
