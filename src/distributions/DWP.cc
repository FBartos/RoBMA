#include "DWP.h"
#include <util/nainf.h>
#include <cmath>
#include <rng/RNG.h>
#include <JRmath.h>


namespace jags {
namespace RoBMA {

DWP::DWP() : VectorDist("dwpois", 2) {}


bool DWP::checkParameterLength(std::vector<unsigned int> const &len) const
{
  return true;
}

bool DWP::checkParameterValue(std::vector<double const *> const &par,
          std::vector<unsigned int> const &len) const
{
  bool lambda_OK = std::isfinite(*par[0]) && *par[0] >= 0.0;
  bool weight_OK = std::isfinite(*par[1]) && *par[1] > 0.0;

  return lambda_OK && weight_OK;
}

double DWP::logDensity(double const *x, unsigned int length, PDFType type,
          std::vector<double const *> const &par,
          std::vector<unsigned int> const &len,
          double const *lower, double const *upper) const
{
  double lambda = *par[0];
  double weight = *par[1];

  if(!std::isfinite(*x) || std::floor(*x) != *x || *x < 0.0){
    return JAGS_NEGINF;
  }

  return dpois(*x, lambda, true) * weight;
}

void DWP::randomSample(double *x, unsigned int length,
          std::vector<double const *> const &par,
          std::vector<unsigned int> const &len,
          double const *lower, double const *upper,
          RNG *rng) const
{
  double lambda = *par[0];
  x[0] = qpois(rng->uniform(), lambda, true, false);
}

void DWP::support(double *lower, double *upper, unsigned int length,
       std::vector<double const *> const &par,
       std::vector<unsigned int> const &len) const
{
  for (unsigned int i = 0; i < length; ++i) {
    lower[i] = 0.0;
    upper[i] = JAGS_POSINF;
  }
}

unsigned int DWP::length(std::vector<unsigned int> const &len) const
{
  return 1;
}


void DWP::typicalValue(double *x, unsigned int length,
          std::vector<double const *> const &par,
          std::vector<unsigned int> const &len,
          double const *lower, double const *upper) const
{
  x[0] = std::floor(*par[0]);
}


bool DWP::isSupportFixed(std::vector<bool> const &fixmask) const
{
  return true;
}

bool DWP::isDiscreteValued(std::vector<bool> const &mask) const
{
  return true;
}


}
}
