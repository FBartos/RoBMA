#include "DMNMIX.h"

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

DMNMIX::DMNMIX() : ArrayDist("dmnmix", 3) {
    // Constructor for DMNMIX
    // "DMNMIX" is the name of the distribution
    // 3 parameters: mu, sigma, weights
}

std::vector<unsigned int> DMNMIX::dim(std::vector<std::vector<unsigned int>> const &dims) const {
  unsigned int mu_length = dims[0][0];      // Length of mu (M * K)
  unsigned int weights_length = dims[2][0]; // Length of weights (M)
  if (weights_length == 0) {
    return std::vector<unsigned int>(1, 0);
  }
  unsigned int K = mu_length / weights_length;
  return std::vector<unsigned int>(1, K);
}


bool DMNMIX::checkParameterDim(std::vector<std::vector<unsigned int>> const &dims) const {
  if (dims.size() != 3) {
    return false;
  }

  unsigned int mu_length = dims[0][0];     // Length of mu (M * K)
  unsigned int sigma_length = dims[1][0];  // Length of sigma (M * K * K)
  unsigned int weights_length = dims[2][0]; // Length of weights (M)

  if (weights_length == 0) {
    return false;
  }

  unsigned int M = weights_length;
  unsigned int K = mu_length / M;

  // Check that mu and sigma lengths are consistent with M and K
  bool mu_OK = (mu_length == M * K);
  bool sigma_OK = (sigma_length == M * K * K);
  bool weights_OK = (weights_length == M);

  return mu_OK && sigma_OK && weights_OK;
}

bool DMNMIX::checkParameterValue(std::vector<double const *> const &par,
                                 std::vector<std::vector<unsigned int>> const &dims) const {
  const double *sigma = par[1];
  const double *weights = par[2];

  unsigned int mu_length = dims[0][0];
  unsigned int weights_length = dims[2][0];

  unsigned int M = weights_length;
  unsigned int K = mu_length / M;

  bool sigma_OK = true;

  for (unsigned int m = 0; m < M; ++m) {
    const double *sigma_m = sigma + m * K * K;

    // Check that sigma_m is symmetric positive definite
    sigma_OK = sigma_OK && check_symmetric_ispd(sigma_m, K);
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

  return sigma_OK;
}

double DMNMIX::logDensity(double const *x, unsigned int length, PDFType type,
                          std::vector<double const *> const &par,
                          std::vector<std::vector<unsigned int>> const &dims,
                          double const *lower, double const *upper) const {
  const double *mu = par[0];
  const double *sigma = par[1];
  const double *weights = par[2];

  unsigned int mu_length = dims[0][0];
  unsigned int weights_length = dims[2][0];

  unsigned int M = weights_length;
  unsigned int K = mu_length / M;

  std::vector<double> log_probs(M);

  for (unsigned int m = 0; m < M; ++m) {
    const double *mu_m = mu + m * K;
    const double *sigma_m = sigma + m * K * K;

    double log_weight = std::log(weights[m]);
    double log_density_m = cpp_mnorm_lpdf(x, mu_m, sigma_m, K);
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

void DMNMIX::randomSample(double *x, unsigned int length,
                          std::vector<double const *> const &par,
                          std::vector<std::vector<unsigned int>> const &dims,
                          double const *lower, double const *upper,
                          RNG *rng) const {
  const double *mu = par[0];
  const double *sigma = par[1];
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

  // Extract the corresponding mu and sigma
  const double *mu_m = mu + m * K;
  const double *sigma_m = sigma + m * K * K;

  // Simulate from N(mu_m, sigma_m)
  simulate_mnorm(x, mu_m, sigma_m, K, rng);
}


void DMNMIX::support(double *lower, double *upper, unsigned int length,
                     std::vector<double const *> const &par,
                     std::vector<std::vector<unsigned int>> const &dims) const {
    // The support is the entire real space for each dimension
    for (unsigned int i = 0; i < length; ++i) {
        lower[i] = JAGS_NEGINF;
        upper[i] = JAGS_POSINF;
    }
}

void DMNMIX::typicalValue(double *x, unsigned int length,
                          std::vector<double const *> const &par,
                          std::vector<std::vector<unsigned int>> const &dims,
                          double const *lower, double const *upper) const {
    // Not implemented
}

bool DMNMIX::isSupportFixed(std::vector<bool> const &fixmask) const {
    return true;
}

void DMNMIX::simulate_mnorm(double *x, const double *mu, const double *sigma, int K, RNG *rng) const {
  // Compute Cholesky decomposition of sigma
  // Allocate memory for U and compute Cholesky decomposition
  double *U = new double[K * K];
  cholesky_decomposition(U, sigma, K);
  // print_matrix(U, K, "U (Cholesky factor)");

  // Generate a vector of independent standard normal variates
  double *z = new double[K];
  for (int i = 0; i < K; ++i) {
    z[i] = rnorm(0.0, 1.0, rng);
  }

  // Compute x = mu + U^T * z
  for (int i = 0; i < K; ++i) {
    x[i] = mu[i];
    for (int j = 0; j <= i; ++j) {
      x[i] += U[j + i * K] * z[j];
    }
  }

  // Clean up
  delete[] U;
  delete[] z;
}

} // namespace RoBMA
} // namespace jags
