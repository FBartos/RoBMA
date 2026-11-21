#include "DWP.h"
#include <util/nainf.h>
#include <cmath>
#include <rng/RNG.h>
#include <JRmath.h>


namespace jags {
namespace RoBMA {

DWP::DWP() : VectorDist("dwpois", 2) {}


bool DWP::checkParameterLength(std::vector<unsigned long> const &len) const
{
  return true;
}

bool DWP::checkParameterValue(std::vector<double const *> const &par,
          std::vector<unsigned long> const &len) const
{
  bool lambda_OK = std::isfinite(*par[0]) && *par[0] >= 0.0;
  bool weight_OK = std::isfinite(*par[1]) && *par[1] > 0.0;

  return lambda_OK && weight_OK;
}

double DWP::logDensity(double const *x, PDFType type,
          std::vector<double const *> const &par,
          std::vector<unsigned long> const &len) const
{
  double lambda = *par[0];
  double weight = *par[1];

  if(!std::isfinite(*x) || std::floor(*x) != *x || *x < 0.0){
    return JAGS_NEGINF;
  }

  return dpois(*x, lambda, true) * weight;
}

void DWP::randomSample(double *x,
          std::vector<double const *> const &par,
          std::vector<unsigned long> const &len,
          RNG *rng) const
{
  double lambda = *par[0];
  x[0] = qpois(rng->uniform(), lambda, true, false);
}

void DWP::support(double *lower, double *upper,
       std::vector<double const *> const &par,
       std::vector<unsigned long> const &len) const
{
  *lower = 0.0;
  *upper = JAGS_POSINF;
}

unsigned long DWP::length(std::vector<unsigned long> const &len) const
{
  return 1;
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
