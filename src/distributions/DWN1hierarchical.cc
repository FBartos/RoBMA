#include "DWN1hierarchical.h"
#include <util/nainf.h>
#include <cmath>
#include <rng/RNG.h>
#include <JRmath.h>

#include "../source/mnorm.h"
#include "../source/wmnorm.h"
#include "../source/tools.h"

// PARAMETER DEFINITIONS for dwnorm_hierarchical
// 1. theta   (Study-specific true effect)
// 2. prec_se (Precision of the study standard error: 1/se^2)
// 3. pop_mu  (Population Mean)
// 4. pop_tau (Population Heterogeneity SD)
// 5. crit_x  (Cutoffs vector)
// 6. omega   (Weights vector)

#define crit_x(par) (par[4])
#define omega(par)  (par[5])

// Dimensions
#define n_crit_x(len) (len[4])
#define n_omega(len)  (len[5])


namespace jags {
namespace RoBMA {

// Registering "dwnorm_hierarchical" with 6 parameters
DWN1hierarchical::DWN1hierarchical() : VectorDist("dwnorm_hierarchical", 6) {}


bool DWN1hierarchical::checkParameterLength(std::vector<unsigned int> const &len) const
{
  // there is one less cut-point then weights
  // Note: Indices shifted due to 4 scalar parameters preceding vectors
  return n_crit_x(len) == n_omega(len) - 1;
}

bool DWN1hierarchical::checkParameterValue(std::vector<double const *> const &par,
                           std::vector<unsigned int> const &len) const
{
  bool omega_OK  = true;

  // Check all weights are probabilities [0, 1]
  for(unsigned j = 0; j < (n_omega(len)); ++j){
    omega_OK = omega_OK && ( omega(par)[j] >= 0.0 ) && ( omega(par)[j] <= 1.0 );
  }

  // Check Precisions/SDs are positive
  // par[1] is prec_se (1/se^2)
  bool se_prec_OK = *par[1] > 0.0;
  
  // par[3] is pop_tau (SD). Strictly speaking SD >= 0.
  bool tau_OK = *par[3] >= 0.0;

  return omega_OK && se_prec_OK && tau_OK;
}

double DWN1hierarchical::logDensity(double const *x, unsigned int length, PDFType type,
                        std::vector<double const *> const &par,
                        std::vector<unsigned int> const &len,
                        double const *lower, double const *upper) const
{
  // 1. Extract Parameters
  const double *theta   = par[0];
  const double prec_se  = *par[1];
  const double *pop_mu  = par[2];
  const double pop_tau  = *par[3];
  
  const double *crit_x_ptr = par[4];
  const double *omega_ptr  = par[5];

  // Derived values
  const double sigma_se = std::sqrt(1.0 / prec_se);
  
  // 2. Information about the dimensions
  const int J = len[5]; // length of weights

  // 3. CALL THE CUSTOM C++ LPDF FUNCTION
  // We call the modified logic discussed previously:
  // Numerator uses (x | theta, sigma_se)
  // Denominator uses (pop_mu, sqrt(sigma_se^2 + pop_tau^2))
  
  // Note: cpp_wnorm_hierarchical_lpdf must be defined/included in your headers
  double log_lik = cpp_wnorm_hierarchical_lpdf(x, theta, &sigma_se, pop_mu, &pop_tau, crit_x_ptr, omega_ptr, J);

  return log_lik;
}

void DWN1hierarchical::randomSample(double *x, unsigned int length,
                        std::vector<double const *> const &par,
                        std::vector<unsigned int> const &len,
                        double const *lower, double const *upper,
                        RNG *rng) const
{
  // Not implemented for custom selection models usually
  // Requires rejection sampling logic
}

void DWN1hierarchical::support(double *lower, double *upper, unsigned int length,
                   std::vector<double const *> const &par,
                   std::vector<unsigned int> const &len) const
{
  for (unsigned int i = 0; i < length; ++i) {
    lower[i] = JAGS_NEGINF;
    upper[i] = JAGS_POSINF;
  }
}

unsigned int DWN1hierarchical::length(std::vector<unsigned int> const &len) const
{
  return 1;
}

void DWN1hierarchical::typicalValue(double *x, unsigned int length,
                        std::vector<double const *> const &par,
                        std::vector<unsigned int> const &len,
                        double const *lower, double const *upper) const
{
  // Not implemented
}

bool DWN1hierarchical::isSupportFixed(std::vector<bool> const &fixmask) const
{
  return true;
}

}
}