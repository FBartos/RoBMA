#include "selnorm-mv.h"

#include <R_ext/Lapack.h>

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
