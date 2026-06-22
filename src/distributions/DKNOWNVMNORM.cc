#include "DKNOWNVMNORM.h"

#include <R_ext/Lapack.h>

#include <JRmath.h>
#include <rng/RNG.h>
#include <util/nainf.h>

#include <cmath>
#include <limits>
#include <vector>

#ifndef FCONE
# define FCONE
#endif

namespace jags {
namespace RoBMA {

DKNOWNVMNORM::DKNOWNVMNORM() : VectorDist("dknown_v_mnorm", 3) {}

bool DKNOWNVMNORM::checkParameterLength(std::vector<unsigned int> const &len) const
{
  if (len[0] == 0 || len[0] != len[1]) {
    return false;
  }

  const unsigned int k = len[0];
  return len[2] == k * (k + 1) / 2;
}

bool DKNOWNVMNORM::checkParameterValue(std::vector<double const *> const &par,
                                       std::vector<unsigned int> const &len) const
{
  const unsigned int k = len[0];

  for (unsigned int i = 0; i < k; ++i) {
    if (!std::isfinite(par[0][i]) || !std::isfinite(par[1][i]) ||
        par[1][i] < 0) {
      return false;
    }
  }

  for (unsigned int i = 0; i < len[2]; ++i) {
    if (!std::isfinite(par[2][i])) {
      return false;
    }
  }

  return true;
}

double DKNOWNVMNORM::logDensity(double const *x, unsigned int length,
                                PDFType type,
                                std::vector<double const *> const &par,
                                std::vector<unsigned int> const &len,
                                double const *lower, double const *upper) const
{
  const int k = static_cast<int>(length);
  const double two_pi = 6.283185307179586476925286766559;

  std::vector<double> covariance(static_cast<unsigned int>(k * k), 0.0);
  unsigned int pos = 0;
  for (int col = 0; col < k; ++col) {
    for (int row = col; row < k; ++row) {
      const double value = par[2][pos++];
      covariance[static_cast<unsigned int>(row + k * col)] = value;
    }
  }

  for (int i = 0; i < k; ++i) {
    covariance[static_cast<unsigned int>(i + k * i)] += par[1][i];
  }

  int info = 0;
  F77_CALL(dpotrf)("L", &k, covariance.data(), &k, &info FCONE);
  if (info != 0) {
    return JAGS_NEGINF;
  }

  double log_det = 0.0;
  std::vector<double> residual(static_cast<unsigned int>(k));
  for (int i = 0; i < k; ++i) {
    const double diag = covariance[static_cast<unsigned int>(i + k * i)];
    if (!(diag > 0) || !std::isfinite(diag) ||
        !std::isfinite(x[i]) || !std::isfinite(par[0][i])) {
      return JAGS_NEGINF;
    }
    log_det += 2.0 * std::log(diag);
    residual[static_cast<unsigned int>(i)] = x[i] - par[0][i];
  }

  for (int i = 0; i < k; ++i) {
    double value = residual[static_cast<unsigned int>(i)];
    for (int j = 0; j < i; ++j) {
      value -= covariance[static_cast<unsigned int>(i + k * j)] *
        residual[static_cast<unsigned int>(j)];
    }
    residual[static_cast<unsigned int>(i)] =
      value / covariance[static_cast<unsigned int>(i + k * i)];
  }

  double quad = 0.0;
  for (int i = 0; i < k; ++i) {
    const double z = residual[static_cast<unsigned int>(i)];
    if (!std::isfinite(z)) {
      return JAGS_NEGINF;
    }
    quad += z * z;
  }

  return -0.5 * (static_cast<double>(k) * std::log(two_pi) + log_det + quad);
}

void DKNOWNVMNORM::randomSample(double *x, unsigned int length,
                                std::vector<double const *> const &par,
                                std::vector<unsigned int> const &len,
                                double const *lower, double const *upper,
                                RNG *rng) const
{
  const int k = static_cast<int>(length);
  std::vector<double> covariance(static_cast<unsigned int>(k * k), 0.0);
  unsigned int pos = 0;
  for (int col = 0; col < k; ++col) {
    for (int row = col; row < k; ++row) {
      covariance[static_cast<unsigned int>(row + k * col)] = par[2][pos++];
    }
  }

  for (int i = 0; i < k; ++i) {
    covariance[static_cast<unsigned int>(i + k * i)] += par[1][i];
  }

  int info = 0;
  F77_CALL(dpotrf)("L", &k, covariance.data(), &k, &info FCONE);
  if (info != 0) {
    for (int i = 0; i < k; ++i) {
      x[i] = std::numeric_limits<double>::quiet_NaN();
    }
    return;
  }

  std::vector<double> z(static_cast<unsigned int>(k));
  for (int i = 0; i < k; ++i) {
    z[static_cast<unsigned int>(i)] = qnorm(rng->uniform(), 0.0, 1.0,
                                            true, false);
  }

  for (int i = 0; i < k; ++i) {
    double value = par[0][i];
    for (int j = 0; j <= i; ++j) {
      value += covariance[static_cast<unsigned int>(i + k * j)] *
        z[static_cast<unsigned int>(j)];
    }
    x[i] = value;
  }
}

void DKNOWNVMNORM::support(double *lower, double *upper, unsigned int length,
                           std::vector<double const *> const &par,
                           std::vector<unsigned int> const &len) const
{
  for (unsigned int i = 0; i < length; ++i) {
    lower[i] = JAGS_NEGINF;
    upper[i] = JAGS_POSINF;
  }
}

unsigned int DKNOWNVMNORM::length(std::vector<unsigned int> const &len) const
{
  return len[0];
}

void DKNOWNVMNORM::typicalValue(double *x, unsigned int length,
                                std::vector<double const *> const &par,
                                std::vector<unsigned int> const &len,
                                double const *lower, double const *upper) const
{
  for (unsigned int i = 0; i < length; ++i) {
    x[i] = par[0][i];
  }
}

bool DKNOWNVMNORM::isSupportFixed(std::vector<bool> const &fixmask) const
{
  return true;
}

}
}
