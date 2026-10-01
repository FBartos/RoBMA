#ifndef ROBMA_SELNORM_MIXTURE_COMPRESSION_H
#define ROBMA_SELNORM_MIXTURE_COMPRESSION_H

// Order-seven density-grid compression. The caller retains the
// original Gaussian evaluator, selection weights, SE Jacobian and mass checks.
// Bounds use conventional floating-point analysis, not an interval guarantee
// for library transcendental functions. No likelihood or RNG uses this plan.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <new>
#include <utility>
#include <type_traits>
#include <vector>

namespace selnorm_mixture_compression_detail {
constexpr long double phi_zero = 0.398942280401432677939946059934381868L;
constexpr long double phi_one = 0.241970724519143349797830192935560655L;
constexpr long double pi = 3.141592653589793238462643383279502884L;
constexpr long double eighth_constant = phi_zero / 384.0L;

inline long double round_up(long double value)
{
  return value == 0.0L ? value : std::nextafter(value,
    std::numeric_limits<long double>::infinity());
}

inline long double gamma(long double operations, long double epsilon)
{
  const long double product = operations * epsilon;
  return product < 1.0L ? product / (1.0L - product) :
    std::numeric_limits<long double>::infinity();
}

template<typename Value>
inline bool normal_or_zero(Value value)
{
  return std::isfinite(value) && (value == 0.0L || std::isnormal(value));
}

// sup |t^k phi(t)|. Prepared once, not in the grouping or grid loop.
inline std::array<long double, 8> const &monomial_bound()
{
  static const std::array<long double, 8> values = [] {
    std::array<long double, 8> out{};
    out[0] = phi_zero;
    for (int k = 1; k < 8; ++k) out[k] = phi_zero *
      std::exp(0.5L * k * (std::log(static_cast<long double>(k)) - 1.0L));
    return out;
  }();
  return values;
}

// sup |phi^(k)| / k!, using exact even bounds and Fourier odd bounds.
inline std::array<long double, 8> derivative_bound()
{
  return {phi_zero, phi_one, phi_zero / 2, 1 / (3 * pi),
    phi_zero / 8, 1 / (15 * pi), phi_zero / 48, 1 / (105 * pi)};
}

template<typename Accumulation>
inline std::array<long double, 8> polynomial(
    std::array<Accumulation, 9> const &input, bool absolute)
{
  // Coefficient construction retains the existing long-double arithmetic.
  std::array<long double, 9> m;
  for (int k = 0; k <= 8; ++k) m[k] = input[k];
  const long double sign = absolute ? 1.0L : -1.0L;
  return {m[0] + sign*m[2]/2 + m[4]/8 + sign*m[6]/48,
    m[1] + sign*m[3]/2 + m[5]/8 + sign*m[7]/48,
    m[2]/2 + sign*m[4]/4 + m[6]/16,
    m[3]/6 + sign*m[5]/12 + m[7]/48,
    m[4]/24 + sign*m[6]/48,
    m[5]/120 + sign*m[7]/240, m[6]/720, m[7]/5040};
}

template<typename Accumulation> struct Bin {
  std::array<Accumulation, 9> moment{}, absolute_moment{};
  double center = 0, minimum = 0, maximum = 0;
  std::size_t count = 0;
};
template<typename Accumulation> struct Assignment { std::size_t original, bin; Accumulation weight; };
}

struct SelNormGaussianMomentGroup {
  double mean = 0;
  std::array<double, 8> coefficient{};
  std::array<long double, 8> coefficient_long{};
  long double scaled_mass = 0, scaled_error = 0;
  long double minimum_offset = 0, maximum_offset = 0;
  long double tail_radius = std::numeric_limits<long double>::infinity();
  std::size_t source_begin = 0, source_count = 0;
  bool double_safe = false, omit_all = false, constant = false;
};

struct SelNormMixtureCompression {
  std::vector<SelNormGaussianMomentGroup> groups;
  // Each group's range indexes this array; entries index the ORIGINAL input.
  // Input arrays are not copied or owned. The core retains them for fallback.
  std::vector<std::size_t> source_index;
  long double log_scale = 0, absolute_error = 0;
  long double input_omission_mass = 0; // Incremental loss, not all compression L1 error.
  double sd = 0;
  std::size_t input_components = 0;
};

// Prepare fixed-width midpoint bins in linear time: determine the scale,
// compute each exponential once, accumulate bins/moments once, pack indices.
// False leaves *out untouched and selects the existing original-mixture path.
template<typename Accumulation>
inline bool selnorm_compress_gaussian_mixture_accumulate(
    double const *mean, double const *log_coefficient, std::size_t count,
    double sd, long double maximum_absolute_error, SelNormMixtureCompression *out)
{
  using namespace selnorm_mixture_compression_detail;
  if (out == nullptr || mean == nullptr || log_coefficient == nullptr || count < 8 ||
      !(sd > 0) || !std::isfinite(sd) || !(maximum_absolute_error > 0) ||
      !std::isfinite(maximum_absolute_error)) return false;
  try {
    const double negative_infinity = -std::numeric_limits<double>::infinity();
    const long double original_allowance = maximum_absolute_error;
    long double cutoff = -std::numeric_limits<long double>::infinity();
    long double omitted_coefficient_upper = 0, omitted_height_upper = 0;
    // This schedules at most one quarter of the unchanged PDF allowance.
    // Admission below uses the actual computed cutoff and its explicit bound.
    const long double proposed_cutoff = std::log(original_allowance) +
      std::log(static_cast<long double>(sd)) - std::log(4.0L) -
      std::log(static_cast<long double>(count)) - std::log(phi_zero);
    if (std::isfinite(proposed_cutoff)) {
      const long double coefficient = std::exp(proposed_cutoff);
      if (coefficient > 0 && std::isnormal(coefficient)) {
        const long double upper = round_up(coefficient);
        const long double product = upper * round_up(phi_zero);
        const long double height = round_up(product) / sd;
        if (std::isnormal(product) && height > 0 && std::isnormal(height)) {
          cutoff = proposed_cutoff;
          omitted_coefficient_upper = upper;
          omitted_height_upper = round_up(height);
        }
      }
    }
    std::size_t omitted = 0;
    double largest_log = negative_infinity;
    double minimum = std::numeric_limits<double>::infinity(), maximum = negative_infinity;
    std::size_t positive = 0;
    for (std::size_t i = 0; i < count; ++i) {
      if (!std::isfinite(mean[i]) || std::isnan(log_coefficient[i]) ||
          log_coefficient[i] == std::numeric_limits<double>::infinity()) return false;
      if (log_coefficient[i] == negative_infinity) continue;
      if (static_cast<long double>(log_coefficient[i]) <= cutoff) {
        ++omitted;
        continue;
      }
      ++positive;
      largest_log = std::max(largest_log, log_coefficient[i]);
      minimum = std::min(minimum, mean[i]);
      maximum = std::max(maximum, mean[i]);
    }
    const long double omission_error = omitted == 0 ? 0.0L :
      round_up(static_cast<long double>(omitted) * omitted_height_upper);
    const long double omission_mass = omitted == 0 ? 0.0L :
      round_up(static_cast<long double>(omitted) * omitted_coefficient_upper);
    if (!std::isfinite(omission_error) || !std::isfinite(omission_mass) ||
        !(omission_error < original_allowance)) return false;
    if (positive == 0 && omitted > 0) {
      // A numerically omitted mixture is not a structural zero. Its complete
      // PDF/mass loss remains visible to the unchanged outer accuracy gates.
      SelNormMixtureCompression empty;
      empty.sd = sd;
      empty.input_components = count;
      empty.absolute_error = omission_error;
      empty.input_omission_mass = omission_mass;
      *out = std::move(empty);
      return true;
    }
    if (positive < 8 || !std::isfinite(largest_log)) return false;
    maximum_absolute_error = original_allowance - omission_error;
    const long double log_scale = largest_log;
    std::vector<Assignment<Accumulation>> assignments;
    assignments.reserve(positive);
    long double total_mass = 0, weight_cost = 0;
    for (std::size_t i = 0; i < count; ++i) {
      if (log_coefficient[i] == negative_infinity ||
          static_cast<long double>(log_coefficient[i]) <= cutoff) continue;
      const long double shifted_log = static_cast<long double>(log_coefficient[i]) - log_scale;
      const Accumulation weight = std::exp(static_cast<Accumulation>(shifted_log));
      if (!(weight > 0) || !std::isnormal(weight)) return false;
      if constexpr (std::is_same<Accumulation, double>::value) {
        // Exponent-argument conversion changes relative weight accuracy with
        // |log(weight)|. Charge it explicitly, including ordinary exp rounding.
        weight_cost = std::max(weight_cost, 4.0L + 2.0L * std::abs(shifted_log));
      }
      total_mass += weight;
      assignments.push_back({i, 0, weight});
    }
    if (!(total_mass > 0) || !std::isfinite(total_mass)) return false;
    // Reserve half the allowance for arithmetic. Actual moments below, not
    // this scheduling formula, determine the returned error bound.
    const long double log_half_width = (std::log(maximum_absolute_error) + std::log(static_cast<long double>(sd)) -
      log_scale - std::log(total_mass) -
      std::log(2 * eighth_constant)) / 8;
    const long double width = 2 * static_cast<long double>(sd) * std::exp(log_half_width);
    if (!(width > 0) || !std::isfinite(width)) return false;
    const long double span = static_cast<long double>(maximum) - minimum;
    const long double last_bin = std::floor(span / width);
    if (!std::isfinite(last_bin) || last_bin < 0 ||
        last_bin >= static_cast<long double>(positive / 2)) return false;
    const std::size_t bin_count = static_cast<std::size_t>(last_bin) + 1;
    std::vector<Bin<Accumulation>> bins(bin_count);
    for (Assignment<Accumulation> &assignment : assignments) {
      const std::size_t i = assignment.original;
      const long double position = std::floor((static_cast<long double>(mean[i]) - minimum) / width);
      if (!(position >= 0) || position >= bin_count) return false;
      const std::size_t index = static_cast<std::size_t>(position);
      Bin<Accumulation> &bin = bins[index];
      if (bin.count == 0) {
        bin.center = static_cast<double>(static_cast<long double>(minimum) +
          (static_cast<long double>(index) + 0.5L) * width);
        if (!std::isfinite(bin.center)) return false;
        bin.minimum = bin.maximum = mean[i];
      }
      const Accumulation weight = assignment.weight;
      const long double original_offset = (static_cast<long double>(mean[i]) - bin.center) / sd;
      const Accumulation offset = static_cast<Accumulation>(original_offset);
      if (!normal_or_zero(offset) || (offset == 0 && original_offset != 0)) return false;
      Accumulation power = 1.0, term = 0.0;
      for (int k = 0; k <= 8; ++k) {
        term = weight * power;
        bin.moment[k] += term;
        bin.absolute_moment[k] += std::abs(term);
        if (k < 8) power *= offset;
      }
      // Magnitudes of successive powers/terms are monotone: for |offset|<=1
      // the last ones are smallest, otherwise largest. Starting weight and
      // power are normal, so these endpoints certify every intermediate.
      // A zero offset has exact zero higher powers and needs no exception.
      if (!normal_or_zero(power) || (power == 0 && offset != 0) ||
          !normal_or_zero(term) || (term == 0 && power != 0)) return false;
      bin.minimum = std::min(bin.minimum, mean[i]);
      bin.maximum = std::max(bin.maximum, mean[i]);
      ++bin.count;
      assignment.bin = index;
    }

    SelNormMixtureCompression plan;
    plan.log_scale = log_scale;
    plan.sd = sd;
    plan.input_components = count;
    plan.input_omission_mass = omission_mass;
    plan.groups.reserve(bin_count);
    std::vector<std::size_t> bin_group(bin_count, 0);
    const auto derivative = derivative_bound();
    const auto &monomial = monomial_bound();
    const long double epsilon_long = std::numeric_limits<long double>::epsilon();
    const long double epsilon_double = std::numeric_limits<double>::epsilon();
    const long double epsilon_accumulator = std::numeric_limits<Accumulation>::epsilon();
    long double total_error = 0.0L;
    std::size_t packed = 0;
    for (std::size_t index = 0; index < bin_count; ++index) {
      Bin<Accumulation> &bin = bins[index];
      if (bin.count == 0) continue;
      // Positive absolute sums cannot recover from overflow. An overflowed
      // signed sum remains nonfinite too; check once before any exact-mean
      // reset or coefficient construction, preserving the failure path.
      for (int k = 0; k <= 8; ++k) {
        if (!std::isfinite(bin.moment[k]) || !std::isfinite(bin.absolute_moment[k])) return false;
      }
      SelNormGaussianMomentGroup group;
      group.source_begin = packed;
      group.source_count = bin.count;
      packed += bin.count;
      // Identical original means need only a Gaussian mass sum. This is exact
      // structural equality, not a small-variance/rank threshold.
      if (bin.minimum == bin.maximum) {
        group.constant = true;
        bin.center = bin.minimum;
        for (int k = 1; k <= 8; ++k) bin.moment[k] = bin.absolute_moment[k] = 0;
      }
      group.mean = bin.center;
      group.scaled_mass = bin.moment[0];
      const long double low_offset = (static_cast<long double>(bin.minimum) - group.mean) / sd;
      const long double high_offset = (static_cast<long double>(bin.maximum) - group.mean) / sd;
      const long double offset_gamma = gamma(4, epsilon_long);
      group.minimum_offset = std::nextafter(low_offset - offset_gamma * std::abs(low_offset),
        -std::numeric_limits<long double>::infinity());
      group.maximum_offset = std::nextafter(high_offset + offset_gamma * std::abs(high_offset),
        std::numeric_limits<long double>::infinity());
      group.coefficient_long = polynomial(bin.moment, false);
      const auto absolute_coefficient = polynomial(bin.absolute_moment, true);
      group.double_safe = true;
      long double casting = 0.0L, polynomial_envelope = 0.0L, moment_error = 0.0L;
      for (int k = 0; k < 8; ++k) {
        if (!normal_or_zero(group.coefficient_long[k])) return false;
        group.coefficient[k] = static_cast<double>(group.coefficient_long[k]);
        group.double_safe = group.double_safe && std::isfinite(group.coefficient[k]) &&
          (group.coefficient[k] == 0 ? group.coefficient_long[k] == 0 : std::isnormal(group.coefficient[k]));
        if (std::isfinite(group.coefficient[k])) casting +=
          std::abs(group.coefficient_long[k] - static_cast<long double>(group.coefficient[k])) * monomial[k];
        polynomial_envelope += absolute_coefficient[k] * monomial[k];
        moment_error += bin.absolute_moment[k] * derivative[k];
      }
      long double sum_gamma = gamma(24.0L * (bin.count + 1) + weight_cost, epsilon_accumulator);
      if constexpr (std::is_same<Accumulation, double>::value) {
        // Express the forward bound using computed absolute moments. This
        // inflation also covers an underestimated positive absolute sum.
        if (!(sum_gamma < 1.0L)) return false;
        sum_gamma /= 1.0L - sum_gamma;
      }
      // Offset arithmetic, powers through eight and signed/absolute sums.
      long double error = eighth_constant * bin.absolute_moment[8] * (1 + sum_gamma) +
        sum_gamma * moment_error;
      // At most four terms form each coefficient. Evaluation includes seven
      // Horner stages, t casting, and the final Gaussian multiplication.
      error += gamma(16, epsilon_long) * polynomial_envelope;
      error += gamma(32, epsilon_long) * polynomial_envelope;
      if (group.double_safe) error += casting + gamma(32, epsilon_double) * polynomial_envelope;
      error = round_up(error / sd);
      if (!std::isfinite(error) || error < 0 || (bin.count > 1 && !std::isnormal(error))) return false;
      // A singleton always uses its original coefficient/primitive in core.
      group.scaled_error = bin.count == 1 ? 0.0L : error;
      if (group.scaled_error > 0) {
        const long double mass_upper = round_up(group.scaled_mass * (1 + sum_gamma));
        const long double radius_square = 2 * (std::log(mass_upper) + std::log(phi_zero) -
          std::log(static_cast<long double>(sd)) - std::log(group.scaled_error));
        if (!std::isfinite(radius_square)) return false;
        group.omit_all = radius_square <= 0;
        group.tail_radius = group.omit_all ? 0.0L :
          std::nextafter(std::sqrt(radius_square), std::numeric_limits<long double>::infinity());
      }
      total_error = round_up(total_error + group.scaled_error);
      bin_group[index] = plan.groups.size();
      plan.groups.push_back(group);
    }
    if (plan.groups.size() * 2 > positive || packed != positive ||
        !std::isfinite(total_error) || total_error < 0) return false;
    plan.absolute_error = total_error == 0 ? 0.0L :
      round_up(std::exp(log_scale + std::log(total_error)));
    if (!std::isfinite(plan.absolute_error) ||
        (total_error > 0 && !(plan.absolute_error > 0)) ||
        plan.absolute_error > maximum_absolute_error) return false;
    if (omitted > 0) plan.absolute_error = round_up(plan.absolute_error + omission_error);
    if (!std::isfinite(plan.absolute_error) || plan.absolute_error > original_allowance) return false;
    plan.source_index.resize(positive);
    std::vector<std::size_t> next(plan.groups.size());
    for (std::size_t group = 0; group < plan.groups.size(); ++group) next[group] = plan.groups[group].source_begin;
    for (Assignment<Accumulation> const &entry : assignments) plan.source_index[next[bin_group[entry.bin]]++] = entry.original;
    *out = std::move(plan);
    return true;
  } catch (std::bad_alloc const &) {
    return false;
  }
}

// The guarded double plan must satisfy the original absolute-error allowance.
// Any range, arithmetic or budget failure reruns the original LD computation;
// unsuccessful attempts never change *out. Platforms without extra LD precision
// retain the original implementation directly.
inline bool selnorm_compress_gaussian_mixture(
    double const *mean, double const *log_coefficient, std::size_t count,
    double sd, long double maximum_absolute_error, SelNormMixtureCompression *out)
{
  if constexpr (std::numeric_limits<long double>::digits > std::numeric_limits<double>::digits) {
    if (selnorm_compress_gaussian_mixture_accumulate<double>(mean, log_coefficient,
        count, sd, maximum_absolute_error, out)) return true;
  }
  return selnorm_compress_gaussian_mixture_accumulate<long double>(mean, log_coefficient,
    count, sd, maximum_absolute_error, out);
}

enum class SelNormGaussianMomentPoint { polynomial, omitted, exact };

// t uses the same emitted center and physical SD as the caller's Gaussian.
// On polynomial, multiply the returned SCALED polynomial by
// exp(plan.log_scale)*Normal(y|group.mean,plan.sd), then apply original SE/weight.
// On exact, core uses original components in group.source_begin/source_count.
// Both polynomial and omitted are covered by group.scaled_error*exp(log_scale).
inline SelNormGaussianMomentPoint selnorm_gaussian_moment_point(
    SelNormGaussianMomentGroup const &group, long double t,
    long double *value)
{
  if (value == nullptr || !std::isfinite(t) || group.source_count == 1) return SelNormGaussianMomentPoint::exact;
  long double distance = 0;
  if (t < group.minimum_offset) distance = group.minimum_offset - t;
  else if (t > group.maximum_offset) distance = t - group.maximum_offset;
  if (group.omit_all || distance > group.tail_radius) {
    *value = 0;
    return SelNormGaussianMomentPoint::omitted;
  }
  if (group.constant) {
    // The planner proved every original mean equal. Keep the tail decision
    // above and the same double/long-double coefficient choice as Horner.
    const double x = static_cast<double>(t);
    const bool ordinary = group.double_safe && std::isfinite(x) &&
      (x == 0 ? t == 0 : std::isnormal(x));
    *value = ordinary ? static_cast<long double>(group.coefficient[0]) : group.coefficient_long[0];
    return SelNormGaussianMomentPoint::polynomial;
  }
  if (group.double_safe) {
    const double x = static_cast<double>(t);
    bool safe = std::isfinite(x) && (x == 0 ? t == 0 : std::isnormal(x));
    double result = group.coefficient[7];
    for (int k = 6; k >= 0 && safe; --k) {
      const double previous = result;
      const double product = previous * x;
      safe = std::isfinite(product) && (product == 0 ? previous == 0 || x == 0 : std::isnormal(product));
      if (!safe) break;
      result = product + group.coefficient[k];
      safe = std::isfinite(result) && (result == 0 ?
        product == 0 && group.coefficient[k] == 0 : std::isnormal(result));
    }
    if (safe && result > 0) {
      *value = result;
      return SelNormGaussianMomentPoint::polynomial;
    }
  }
  long double result = group.coefficient_long[7];
  for (int k = 6; k >= 0; --k) {
    const long double product = result * t;
    if (!selnorm_mixture_compression_detail::normal_or_zero(product) ||
        (product == 0 && result != 0 && t != 0)) return SelNormGaussianMomentPoint::exact;
    result = product + group.coefficient_long[k];
    if (!selnorm_mixture_compression_detail::normal_or_zero(result) ||
        (result == 0 && (product != 0 || group.coefficient_long[k] != 0))) {
      return SelNormGaussianMomentPoint::exact;
    }
  }
  if (!(result > 0) || !std::isfinite(result)) return SelNormGaussianMomentPoint::exact;
  *value = result;
  return SelNormGaussianMomentPoint::polynomial;
}

#endif
