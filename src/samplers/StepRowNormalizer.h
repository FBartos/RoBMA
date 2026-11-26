#ifndef ROBMA_STEP_ROW_NORMALIZER_H
#define ROBMA_STEP_ROW_NORMALIZER_H

// The selected-normal normalizer of one step-selection row (the normalizer of
// a dselnorm_step_switch child in step mode) and its first two derivatives in
// the row's mean, for the location block's Newton-Laplace fits:
//   A(m) = sum_b omega_b P(lower_b < sign y / se < upper_b),  y ~ N(m, sd^2).
// Every interval probability, density and ratio is evaluated on the log
// scale with the tail-safe normal distribution functions, so that rows far in
// a tail keep finite derivatives. Bounds at +-1e300 or beyond are infinite
// (the package's JAGS encoding of open bins).

namespace jags {
namespace RoBMA {

struct StepRowLaw {
  double sd = 1;         // the row's marginal standard deviation
  double se = 1;         // its selection standard error
  int effect_sign = 1;   // the effect direction
  int bins = 0;
  double const *omega = nullptr;  // step weights
  double const *lower = nullptr;  // bin bounds on the z scale
  double const *upper = nullptr;
};

// log A(m), d log A / dm and d^2 log A / dm^2; false when A(m) is not
// positive or a value is not finite.
bool step_row_log_normalizer(StepRowLaw const &law, double mean, double *log_a, double *d1, double *d2);

// The row's own log density at observation x and mean m with its parameters
// read once (power, observed bin, kernel mode, telescoping): the one
// implementation of the row density, which DSELNORMSTEPSWITCH::logDensity
// calls and the location block evaluates directly. The caller has checked the
// parameters; the bounds are read only in step mode with an observed bin.
double step_row_log_density(StepRowLaw const &law, double power, int bin, int mode, bool telescope,
                            double x, double mean);

}
}
#endif
