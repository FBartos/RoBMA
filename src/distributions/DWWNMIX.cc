#include "DWWNMIX.h"

#include <rng/RNG.h>
#include <util/dim.h>
#include <util/nainf.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <JRmath.h>
#include <numeric>

#include "../source/mnorm.h"
#include "../source/wmnorm.h"
#include "../source/tools.h"


namespace jags {
namespace RoBMA {

DWWNMIX::DWWNMIX() : VectorDist("dwwnorm_mix", 7) {}

bool DWWNMIX::checkParameterLength(std::vector<unsigned int> const &len) const
{
  // there is one less cut-point then weights
  return true;
}

bool DWWNMIX::checkParameterValue(std::vector<double const *> const &par,
			    std::vector<unsigned int> const &len) const
{
  return true;
}

// Log Density
double DWWNMIX::logDensity(double const *x, unsigned int length, PDFType type,
			  std::vector<double const *> const &par,
			  std::vector<unsigned int> const &len,
			  double const *lower, double const *upper) const
{
    // extract parameters
    const double *mu     = par[0];
    const double *sigma  = par[1];
    const double *crit_y_all = par[2];
    const double *omega_all  = par[3];

    const double *crit_y_mapping = par[4];
    const double crit_y_mapping_max = static_cast<int>(*par[5]);

    const double weight  = *par[6];

    double log_lik;

    if (crit_y_mapping_max == 0) {
        log_lik = dnorm(*x, *mu, *sigma, true) * weight;
    } else {
        // Map the input crit_y and omega to the actual crit_yalpha
        std::vector<double> crit_y(crit_y_mapping_max);
        std::vector<double> omega(crit_y_mapping_max + 1);
        const int J = crit_y_mapping_max + 1;

        // use first weight and each weight after the treshold change
        omega[0] = omega_all[0];
        for (int i = 0; i < crit_y_mapping_max; ++i) {
            crit_y[i] = crit_y_all[ static_cast<int>(crit_y_mapping[i]) - 1 ];
            omega[i + 1] = omega_all[ static_cast<int>(crit_y_mapping[i]) ];
        }

        log_lik = cpp_wnorm_1s_lpdf(x, mu, sigma, &crit_y[0], &omega[0], J) * weight;
    }

    return log_lik;
}

void DWWNMIX::randomSample(double *x, unsigned int length,
			  std::vector<double const *> const &par,
			  std::vector<unsigned int> const &len,
			  double const *lower, double const *upper,
			  RNG *rng) const
{
  // not implemented
}

void DWWNMIX::support(double *lower, double *upper, unsigned int length,
	     std::vector<double const *> const &par,
	     std::vector<unsigned int> const &len) const
{
  // no idea whether this is correct
  for (unsigned int i = 0; i < length; ++i) {
	  lower[i] = JAGS_NEGINF;
	  upper[i] = JAGS_POSINF;
  }
}

unsigned int DWWNMIX::length(std::vector<unsigned int> const &len) const
{
  // no idea how this works
  return 1;
}


void DWWNMIX::typicalValue(double *x, unsigned int length,
			  std::vector<double const *> const &par,
			  std::vector<unsigned int> const &len,
			  double const *lower, double const *upper) const
{
  // not implemented
}


bool DWWNMIX::isSupportFixed(std::vector<bool> const &fixmask) const
{
  return true;
}


}}
