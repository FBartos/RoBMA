#include "wmnorm.h"
#include "mnorm.h"
#include "tools.h"

#include <util/nainf.h>
#include <JRmath.h>
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <limits>

namespace {

double clamp_weight(double w)
{
  if (w < 0) {
    return 0;
  }
  return w;
}

double maybe_clamp_weight(double w, bool clamp)
{
  return clamp ? clamp_weight(w) : w;
}

int wnorm_1s_weight_index(double x, double const *crit_x, const int J, bool rightmost_closed)
{
  const int n_cuts = J - 1;

  if (n_cuts <= 0) {
    return 0;
  }

  if (rightmost_closed) {
    if (x < crit_x[0]) {
      return 0;
    }
    if (x > crit_x[n_cuts - 1]) {
      return n_cuts;
    }
    if (x == crit_x[n_cuts - 1]) {
      return n_cuts - 1;
    }
    for (int j = 1; j < n_cuts; ++j) {
      if (x < crit_x[j]) {
        return j;
      }
    }
    return n_cuts - 1;
  }

  if (x < crit_x[0]) {
    return 0;
  }
  for (int j = 1; j < n_cuts; ++j) {
    if (x < crit_x[j]) {
      return j;
    }
  }
  return n_cuts;
}

double robma_logspace_add(double log_sum, double log_value)
{
  if (std::isnan(log_value)) {
    return log_sum;
  }
  if (log_sum == -INFINITY) {
    return log_value;
  }
  const double log_max = std::max(log_sum, log_value);
  return log_max + std::log(std::exp(log_sum - log_max) + std::exp(log_value - log_max));
}

double normal_cdf(double x, double mean, double sd)
{
  return pnorm(x, mean, sd, true, false);
}

double interval_log_term(double prob, double weight)
{
  if (prob > 0 && weight > 0) {
    return std::log(prob) + std::log(weight);
  }
  if (prob == 0 || weight == 0) {
    return -INFINITY;
  }
  return std::numeric_limits<double>::quiet_NaN();
}

}

static double cpp_wnorm_1s_log_norm(double const mean, double const sd, double const *crit_x, double const *omega, const int J, const bool clamp_omega)
{
  double cdf_sum = 0;
  double log_sum = -INFINITY;

  for (int j = 0; j < J; ++j) {
    double interval_prob;

    if (j == 0) {
      interval_prob = normal_cdf(crit_x[0], mean, sd);
    } else if (j < J - 1) {
      interval_prob = normal_cdf(crit_x[j], mean, sd) - cdf_sum;
    } else {
      interval_prob = 1 - cdf_sum;
    }

    if (interval_prob < 0) {
      interval_prob = 0;
    }

    cdf_sum += interval_prob;
    log_sum = robma_logspace_add(
      log_sum,
      interval_log_term(interval_prob, maybe_clamp_weight(omega[j], clamp_omega))
    );
  }

  return log_sum;
}

static double cpp_wnorm_1s_lpdf_core(double const x, double const mean, double const sd, double const *crit_x, double const *omega, const int J, const bool clamp_omega, const bool rightmost_closed)
{
  const int weight_index = wnorm_1s_weight_index(x, crit_x, J, rightmost_closed);
  const double weight = maybe_clamp_weight(omega[weight_index], clamp_omega);
  const double log_numerator = dnorm(x, mean, sd, true) + std::log(weight);
  const double log_std_const = cpp_wnorm_1s_log_norm(mean, sd, crit_x, omega, J, false);

  return log_numerator - log_std_const;
}

static double active_mix_crit(double const *crit_y_all, double const *crit_y_mapping, int j)
{
  return crit_y_all[static_cast<int>(crit_y_mapping[j]) - 1];
}

static double active_mix_weight(double const *omega_all, double const *crit_y_mapping, int j, const int omega_stride)
{
  if (j == 0) {
    return omega_all[0];
  }
  return omega_all[static_cast<int>(crit_y_mapping[j - 1]) * omega_stride];
}

static int wnorm_mix_weight_index(double x, double const *crit_y_all, double const *crit_y_mapping, const int crit_y_mapping_max)
{
  if (crit_y_mapping_max <= 0) {
    return 0;
  }

  if (x < active_mix_crit(crit_y_all, crit_y_mapping, 0)) {
    return 0;
  }

  for (int j = 1; j < crit_y_mapping_max; ++j) {
    if (x < active_mix_crit(crit_y_all, crit_y_mapping, j)) {
      return j;
    }
  }

  return crit_y_mapping_max;
}

double cpp_wnorm_mix_log_norm(double const mean, double const sd, double const *crit_y_all, double const *omega_all, double const *crit_y_mapping, const int crit_y_mapping_max, const int omega_stride)
{
  if (crit_y_mapping_max <= 0) {
    return 0;
  }

  double cdf_sum = 0;
  double log_sum = -INFINITY;

  for (int j = 0; j <= crit_y_mapping_max; ++j) {
    double interval_prob;

    if (j == 0) {
      interval_prob = normal_cdf(active_mix_crit(crit_y_all, crit_y_mapping, 0), mean, sd);
    } else if (j < crit_y_mapping_max) {
      interval_prob = normal_cdf(active_mix_crit(crit_y_all, crit_y_mapping, j), mean, sd) - cdf_sum;
    } else {
      interval_prob = 1 - cdf_sum;
    }

    if (interval_prob < 0) {
      interval_prob = 0;
    }

    cdf_sum += interval_prob;
    log_sum = robma_logspace_add(
      log_sum,
      interval_log_term(interval_prob, active_mix_weight(omega_all, crit_y_mapping, j, omega_stride))
    );
  }

  return log_sum;
}

double cpp_wnorm_mix_lpdf_precomputed(double const x, double const mu, double const sigma, double const *crit_y_all, double const *omega_all, double const *crit_y_mapping, const int crit_y_mapping_max, const double log_std_const, const int omega_stride)
{
  if (crit_y_mapping_max <= 0) {
    return dnorm(x, mu, sigma, true);
  }

  const int weight_index = wnorm_mix_weight_index(x, crit_y_all, crit_y_mapping, crit_y_mapping_max);
  const double log_numerator = dnorm(x, mu, sigma, true) +
    std::log(active_mix_weight(omega_all, crit_y_mapping, weight_index, omega_stride));

  return log_numerator - log_std_const;
}

double cpp_wnorm_mix_lpdf(double const x, double const mu, double const sigma, double const *crit_y_all, double const *omega_all, double const *crit_y_mapping, const int crit_y_mapping_max, const int omega_stride)
{
  const double log_std_const = cpp_wnorm_mix_log_norm(
    mu, sigma, crit_y_all, omega_all, crit_y_mapping,
    crit_y_mapping_max, omega_stride
  );

  return cpp_wnorm_mix_lpdf_precomputed(
    x, mu, sigma, crit_y_all, omega_all, crit_y_mapping,
    crit_y_mapping_max, log_std_const, omega_stride
  );
}

double cpp_wnorm_mix_cdf_precomputed(double const q, double const mu, double const sigma, double const *crit_y_all, double const *omega_all, double const *crit_y_mapping, const int crit_y_mapping_max, const double log_denom_total, const bool lower_tail, const bool log_p, const int omega_stride)
{
  if (crit_y_mapping_max <= 0) {
    return pnorm(q, mu, sigma, lower_tail, log_p);
  }

  double log_p_sum = -INFINITY;

  if (std::isinf(q)) {
    log_p_sum = q < 0 ? -INFINITY : 0;
  } else {
    double p_sum = 0;
    int j = 0;

    while (j < crit_y_mapping_max && q > active_mix_crit(crit_y_all, crit_y_mapping, j)) {
      double temp_p_add = normal_cdf(active_mix_crit(crit_y_all, crit_y_mapping, j), mu, sigma) - p_sum;
      if (temp_p_add < 0) {
        temp_p_add = 0;
      }
      p_sum += temp_p_add;
      log_p_sum = robma_logspace_add(
        log_p_sum,
        interval_log_term(temp_p_add, active_mix_weight(omega_all, crit_y_mapping, j, omega_stride))
      );
      ++j;
    }

    double temp_p_add = normal_cdf(q, mu, sigma) - p_sum;
    if (temp_p_add < 0) {
      temp_p_add = 0;
    }
    log_p_sum = robma_logspace_add(
      log_p_sum,
      interval_log_term(temp_p_add, active_mix_weight(omega_all, crit_y_mapping, j, omega_stride))
    );
  }

  double log_prob = log_p_sum - log_denom_total;
  if (log_prob > 0) {
    log_prob = 0;
  }

  if (!lower_tail) {
    log_prob = std::log1p(-std::exp(log_prob));
  }

  if (log_p) {
    return log_prob;
  }
  return std::exp(log_prob);
}

double cpp_wnorm_mix_cdf(double const q, double const mu, double const sigma, double const *crit_y_all, double const *omega_all, double const *crit_y_mapping, const int crit_y_mapping_max, const bool lower_tail, const bool log_p, const int omega_stride)
{
  const double log_denom_total = cpp_wnorm_mix_log_norm(
    mu, sigma, crit_y_all, omega_all, crit_y_mapping,
    crit_y_mapping_max, omega_stride
  );

  return cpp_wnorm_mix_cdf_precomputed(
    q, mu, sigma, crit_y_all, omega_all, crit_y_mapping,
    crit_y_mapping_max, log_denom_total, lower_tail, log_p, omega_stride
  );
}

void cpp_wnorm_mix_moments(double const mu, double const sigma, double const *crit_y_all, double const *omega_all, double const *crit_y_mapping, const int crit_y_mapping_max, double *moment_mean, double *moment_second, const int omega_stride)
{
  if (crit_y_mapping_max <= 0 || sigma < std::sqrt(DBL_EPSILON)) {
    *moment_mean   = mu;
    *moment_second = mu * mu + sigma * sigma;
    return;
  }

  double denom = 0;
  double weighted_mean = 0;
  double weighted_second = 0;

  for (int j = 0; j <= crit_y_mapping_max; ++j) {
    const double lower = j == 0 ? -INFINITY : active_mix_crit(crit_y_all, crit_y_mapping, j - 1);
    const double upper = j == crit_y_mapping_max ? INFINITY : active_mix_crit(crit_y_all, crit_y_mapping, j);

    const double alpha = (lower - mu) / sigma;
    const double beta  = (upper - mu) / sigma;

    double interval_prob = pnorm(beta, 0, 1, true, false) - pnorm(alpha, 0, 1, true, false);
    if (interval_prob < 0) {
      interval_prob = 0;
    }

    const double phi_alpha = dnorm(alpha, 0, 1, false);
    const double phi_beta  = dnorm(beta, 0, 1, false);

    double alpha_phi = alpha * phi_alpha;
    double beta_phi  = beta * phi_beta;
    if (!std::isfinite(alpha_phi)) {
      alpha_phi = 0;
    }
    if (!std::isfinite(beta_phi)) {
      beta_phi = 0;
    }

    const double tail_shift = phi_alpha - phi_beta;
    const double second_shift = alpha_phi - beta_phi;

    const double interval_mean = mu * interval_prob + sigma * tail_shift;
    const double interval_second = (mu * mu + sigma * sigma) * interval_prob +
      2 * mu * sigma * tail_shift + sigma * sigma * second_shift;
    const double weight = clamp_weight(active_mix_weight(omega_all, crit_y_mapping, j, omega_stride));

    denom           += weight * interval_prob;
    weighted_mean   += weight * interval_mean;
    weighted_second += weight * interval_second;
  }

  if (denom > 0) {
    *moment_mean   = weighted_mean / denom;
    *moment_second = weighted_second / denom;
  } else {
    *moment_mean   = std::numeric_limits<double>::quiet_NaN();
    *moment_second = std::numeric_limits<double>::quiet_NaN();
  }
}

double cpp_wnorm_1s_lpdf(double const *x, double const *mu, double const *sigma, double const *crit_x, double const *omega, const int J)
{
  return cpp_wnorm_1s_lpdf_core(*x, *mu, *sigma, crit_x, omega, J, false, false);
}

double cpp_wnorm_2s_lpdf(double const *x, double const *mu, double const *sigma, double const *crit_x, double const *omega, const int J)
{

  // obtain the weight (on log scale)
  double log_w = log_weight_twosided(&x[0], &crit_x[0], &omega[0], J);;

  // the log weighted normal likelihood
  double log_lik = dnorm(*x, *mu, *sigma, true) + log_w;

  // get the standardizing constant
  double log_std_const = log_std_constant_twosided(&x[0], &mu[0], &sigma[0], &crit_x[0], &omega[0], J);

  return log_lik - log_std_const;
}

double cpp_wmnorm_1s_lpdf(double const *x, double const *mu, double const *sigma, double const *crit_x, double const *omega, const int K, const int J)
{

  // obtain product of the weights (on log scale)
  double log_w = 0;
  for(int k = 0; k < K; k++){
    log_w += log_weight_onesided(&x[k], &crit_x[k * (J - 1)], &omega[0], J);
  }

  // the log weighted normal likelihood
  double log_lik = cpp_mnorm_lpdf(&x[0], &mu[0], &sigma[0], K) + log_w;

  // get the standardizing constant
  double log_std_const = log_std_m_constant_onesided(&x[0], &mu[0], &sigma[0], &crit_x[0], &omega[0], K, J);

  return log_lik - log_std_const;
}

double cpp_wmnorm_2s_lpdf(double const *x, double const *mu, double const *sigma, double const *crit_x, double const *omega, const int K, const int J)
{
  // obtain product of the weights (on log scale)
  double log_w = 0;
  for(int k = 0; k < K; k++){
    log_w += log_weight_twosided(&x[k], &crit_x[k * (J - 1)], &omega[0], J);
  }

  // the log weighted normal likelihood
  double log_lik = cpp_mnorm_lpdf(&x[0], &mu[0], &sigma[0], K) + log_w;

  // get the standardizing constant
  double log_std_const = log_std_m_constant_twosided(&x[0], &mu[0], &sigma[0], &crit_x[0], &omega[0], K, J);

  return log_lik - log_std_const;
}
