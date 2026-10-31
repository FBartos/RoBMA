#ifndef ROBMA_PLOT_ROOT_H
#define ROBMA_PLOT_ROOT_H

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <limits>

namespace RoBMA {

// Smallest scale at which a root is resolved: DBL_EPSILON times the smallest
// positive component sd, so a mixture of narrow components is resolved at the
// width of its narrowest component rather than at the width of the bracket.
inline double plot_root_absolute_tolerance(double sigma_min)
{
  return std::max(
    2 * std::numeric_limits<double>::denorm_min(),
    DBL_EPSILON * sigma_min
  );
}

// Brent tolerance at the iterate `point`: 2 DBL_EPSILON |point| plus the
// absolute tolerance. The same local scale sizes the warm-start bracket step.
inline double plot_root_local_tolerance(double point, double absolute_tolerance)
{
  return 2 * DBL_EPSILON * std::fabs(point) + absolute_tolerance;
}

// Brent's method on the sign change of `function`, resolved to the local
// tolerance at the iterate (the root is not refined to adjacent doubles). On
// return `*root_value` is the function value at the returned point.
template <typename Function>
double plot_brent_root(Function function,
                       double lower,
                       double upper,
                       double lower_value,
                       double upper_value,
                       double absolute_tolerance,
                       double *root_value = NULL)
{
  double previous       = lower;
  double current        = upper;
  double bracket        = lower;
  double previous_value = lower_value;
  double current_value  = upper_value;
  double bracket_value  = lower_value;
  double previous_step  = current - previous;
  double current_step   = previous_step;

  for (int iteration = 0; iteration < 100; ++iteration) {
    if (previous_value != 0 && current_value != 0 &&
        std::signbit(previous_value) != std::signbit(current_value)) {
      bracket       = previous;
      bracket_value = previous_value;
      previous_step = current - previous;
      current_step  = previous_step;
    }

    if (std::fabs(bracket_value) < std::fabs(current_value)) {
      const double old_current       = current;
      const double old_current_value = current_value;
      previous       = current;
      previous_value = current_value;
      current        = bracket;
      current_value  = bracket_value;
      bracket        = old_current;
      bracket_value  = old_current_value;
    }

    const double delta = plot_root_local_tolerance(current, absolute_tolerance);
    const double bisection_step = 0.5 * (bracket - current);
    if (current_value == 0 || std::fabs(bisection_step) < delta) {
      if (root_value != NULL) {
        *root_value = current_value;
      }
      return current;
    }

    if (std::fabs(previous_step) > delta &&
        std::fabs(current_value) < std::fabs(previous_value)) {
      double interpolation_step;
      if (previous == bracket) {
        interpolation_step = -current_value * (current - previous) /
          (current_value - previous_value);
      } else {
        const double previous_divided_difference =
          (previous_value - current_value) / (previous - current);
        const double bracket_divided_difference =
          (bracket_value - current_value) / (bracket - current);
        interpolation_step = -current_value *
          (bracket_value * bracket_divided_difference -
           previous_value * previous_divided_difference) /
          (bracket_divided_difference * previous_divided_difference *
           (bracket_value - previous_value));
      }

      if (2 * std::fabs(interpolation_step) < std::min(
            std::fabs(previous_step),
            3 * std::fabs(bisection_step) - delta
          )) {
        previous_step = current_step;
        current_step  = interpolation_step;
      } else {
        previous_step = bisection_step;
        current_step  = bisection_step;
      }
    } else {
      previous_step = bisection_step;
      current_step  = bisection_step;
    }

    previous       = current;
    previous_value = current_value;
    current += std::fabs(current_step) > delta ?
      current_step : std::copysign(delta, bisection_step);
    current_value = function(current);
    if (!std::isfinite(current_value)) {
      return std::numeric_limits<double>::quiet_NaN();
    }
  }

  return std::numeric_limits<double>::quiet_NaN();
}

}

#endif
