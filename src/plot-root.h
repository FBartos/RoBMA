#ifndef ROBMA_PLOT_ROOT_H
#define ROBMA_PLOT_ROOT_H

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <limits>

namespace RoBMA {

inline double plot_root_tolerance(double lower, double upper, double step)
{
  const double tolerance_scale = std::max(
    std::max(std::fabs(lower), std::fabs(upper)),
    step
  );
  return std::max(
    2 * std::numeric_limits<double>::denorm_min(),
    DBL_EPSILON * tolerance_scale
  );
}

template <typename Function>
double plot_brent_root(Function function,
                       double lower,
                       double upper,
                       double lower_value,
                       double upper_value,
                       double tolerance)
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

    const double delta = 0.5 * (
      tolerance + 4 * DBL_EPSILON * std::fabs(current)
    );
    const double bisection_step = 0.5 * (bracket - current);
    if (current_value == 0 || std::fabs(bisection_step) < delta) {
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
