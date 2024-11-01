#include "DMNMIXCD.h"

#include <rng/RNG.h>
#include <util/dim.h>
#include <util/nainf.h>
#include <cmath>
#include <vector>
#include <algorithm>

#include <JRmath.h>
#include "../source/mnorm.h"
#include "../source/tools.h"
#include "../matrix/matrix.h"

namespace jags {
namespace RoBMA {

DMNMIXCD::DMNMIXCD() : ArrayDist("dmnmix_cd", 3) {
  // Constructor for DMNMIXCD
}

std::vector<unsigned int> DMNMIXCD::dim(std::vector<std::vector<unsigned int>> const &dims) const {
  unsigned int mu_length = dims[0][0];      // Length of mu (M * K)
  unsigned int weights_length = dims[2][0]; // Length of weights (M)
  if (weights_length == 0) {
    return std::vector<unsigned int>(1, 0);
  }
  unsigned int K = mu_length / weights_length;
  return std::vector<unsigned int>(1, K);
}

bool DMNMIXCD::checkParameterDim(std::vector<std::vector<unsigned int>> const &dims) const {
  if (dims.size() != 3) {
    return false;
  }

  unsigned int mu_length = dims[0][0];       // Length of mu (M * K)
  unsigned int chol_length = dims[1][0];     // Length of chol (M * K * K)
  unsigned int weights_length = dims[2][0];  // Length of weights (M)

  if (weights_length == 0) {
    return false;
  }

  unsigned int M = weights_length;
  unsigned int K = mu_length / M;

  bool mu_OK = (mu_length == M * K);
  bool chol_OK = (chol_length == M * K * K);
  bool weights_OK = (weights_length == M);

  return mu_OK && chol_OK && weights_OK;
}

bool DMNMIXCD::checkParameterValue(std::vector<double const *> const &par,
                                   std::vector<std::vector<unsigned int>> const &dims) const {
  const double *chol = par[1];
  const double *weights = par[2];

  unsigned int mu_length = dims[0][0];
  unsigned int weights_length = dims[2][0];

  unsigned int M = weights_length;
  unsigned int K = mu_length / M;

  // Check if chol is upper triangular
  bool chol_OK = true;

  for (unsigned int m = 0; m < M; ++m) {
    const double *chol_m = chol + m * K * K;
    chol_OK = chol_OK && check_upper_triangular(chol_m, K);
  }

  // Check that weights are non-negative and sum to 1
  double sum_weights = 0.0;
  for (unsigned int m = 0; m < M; ++m) {
    if (weights[m] < 0) {
      return false;
    }
    sum_weights += weights[m];
  }

  if (std::fabs(sum_weights - 1.0) > 1e-6) {
    return false;
  }

  return chol_OK;
}

double DMNMIXCD::logDensity(double const *x, unsigned int length, PDFType type,
                            std::vector<double const *> const &par,
                            std::vector<std::vector<unsigned int>> const &dims,
                            double const *lower, double const *upper) const {
  const double *mu = par[0];
  const double *chol = par[1];
  const double *weights = par[2];

  unsigned int mu_length = dims[0][0];
  unsigned int weights_length = dims[2][0];

  unsigned int M = weights_length;
  unsigned int K = mu_length / M;

  std::vector<double> log_probs(M);

  for (unsigned int m = 0; m < M; ++m) {
    const double *mu_m = mu + m * K;
    const double *chol_m = chol + m * K * K;

    double log_weight = std::log(weights[m]);
    double log_density_m = cpp_mnorm_chol_lpdf(x, mu_m, chol_m, K);
    log_probs[m] = log_weight + log_density_m;
  }

  // Use log-sum-exp to compute the total log density
  double max_log_prob = *std::max_element(log_probs.begin(), log_probs.end());
  double sum_exp = 0.0;
  for (unsigned int m = 0; m < M; ++m) {
    sum_exp += std::exp(log_probs[m] - max_log_prob);
  }
  double log_density = max_log_prob + std::log(sum_exp);

  return log_density;
}

void DMNMIXCD::randomSample(double *x, unsigned int length,
                            std::vector<double const *> const &par,
                            std::vector<std::vector<unsigned int>> const &dims,
                            double const *lower, double const *upper,
                            RNG *rng) const {
  const double *mu = par[0];
  const double *chol = par[1];
  const double *weights = par[2];

  unsigned int mu_length = dims[0][0];
  unsigned int weights_length = dims[2][0];

  unsigned int M = weights_length;
  unsigned int K = mu_length / M;

  // Generate a component index m according to the mixture weights
  double u = rng->uniform();
  double cum_prob = 0.0;
  int m = M - 1; // Default to last component in case of rounding errors
  for (unsigned int i = 0; i < M; ++i) {
    cum_prob += weights[i];
    if (u <= cum_prob) {
      m = i;
      break;
    }
  }

  // Extract the corresponding mu and chol
  const double *mu_m = mu + m * K;
  const double *chol_m = chol + m * K * K;

  // Simulate from N(mu_m, chol_m)
  simulate_mnorm_chol(x, mu_m, chol_m, K, rng);
}

void DMNMIXCD::support(double *lower, double *upper, unsigned int length,
                       std::vector<double const *> const &par,
                       std::vector<std::vector<unsigned int>> const &dims) const {
  for (unsigned int i = 0; i < length; ++i) {
    lower[i] = JAGS_NEGINF;
    upper[i] = JAGS_POSINF;
  }
}

void DMNMIXCD::typicalValue(double *x, unsigned int length,
                            std::vector<double const *> const &par,
                            std::vector<std::vector<unsigned int>> const &dims,
                            double const *lower, double const *upper) const {
  // Not implemented
}

bool DMNMIXCD::isSupportFixed(std::vector<bool> const &fixmask) const {
  return true;
}

} // namespace RoBMA
} // namespace jags
