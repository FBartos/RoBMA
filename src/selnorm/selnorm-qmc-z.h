#ifndef ROBMA_SELNORM_QMC_Z_H
#define ROBMA_SELNORM_QMC_Z_H

#include <algorithm>

// Stopping multiplier of the randomized quasi-Monte Carlo (QMC) selection
// normalizers: a randomized estimate is accepted only when
//   max(selnorm_qmc_z * relative_mcse, relative_change) <= relative_tolerance
// (normalizers that report no change: selnorm_qmc_z * relative_mcse), so the
// Monte Carlo error is held to the tolerance at about the two-sided 95% level
// of the default 8 scrambles instead of at one standard error. The value is the
// 97.5% quantile of Student's t distribution with 7 degrees of freedom
// (2.3646) rounded to 2.36. It is fixed: a different number of scrambles does
// not change it. The native stops, the JAGS wrappers and the R judges of
// returned QMC diagnostics all use this one constant (R reads it through
// `RoBMA_selnorm_qmc_z`). The reported `relative_mcse` stays the raw standard
// error.
constexpr double selnorm_qmc_z = 2.36;

// The quantity compared with `relative_tolerance` by the factor, cluster and
// conditioned-sampling normalizers, which return both diagnostics.
inline double selnorm_qmc_diagnostic(double relative_mcse, double relative_change)
{
  return std::max(selnorm_qmc_z * relative_mcse, relative_change);
}

// The clause of the rejection messages that names the rule (one %g: selnorm_qmc_z).
// The R judges build the same wording in `.selection_qmc_criterion()`.
#define SELNORM_QMC_CRITERION_WITH_CHANGE \
  " (accepted when %g times the standard error and the change are within 'relative_tolerance')"
#define SELNORM_QMC_CRITERION_ERROR_ONLY \
  " (accepted when %g times the standard error is within 'relative_tolerance')"

#endif
