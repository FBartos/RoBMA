#include <chrono>
#include <cstdio>
#include <cstdlib>
#include "selnorm-mv.h"
#include "selnorm-mixture-compression.h"

#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include <JRmath.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstring>
#include <functional>
#include <iterator>
#include <list>
#include <memory_resource>
#include <mutex>
#include <new>
#include <optional>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <limits>
#include <stdexcept>
#include <vector>

#include "selnorm.h"

#if defined(_OPENMP)
#include <omp.h>
#endif

#ifndef FCONE
# define FCONE
#endif

namespace {

double log_add_exp(double x, double y)
{
  if (x == -std::numeric_limits<double>::infinity()) return y;
  if (y == -std::numeric_limits<double>::infinity()) return x;
  const double maximum = std::max(x, y);
  return maximum + std::log(std::exp(x - maximum) + std::exp(y - maximum));
}

double qmc_value(const double *values, int scrambles, int points,
                 int dimensions, int scramble, int point, int dimension)
{
  const std::size_t index = static_cast<std::size_t>(scramble) +
    static_cast<std::size_t>(scrambles) *
    (static_cast<std::size_t>(point) +
     static_cast<std::size_t>(points) * static_cast<std::size_t>(dimension));
  const std::size_t total = static_cast<std::size_t>(scrambles) *
    static_cast<std::size_t>(points) * static_cast<std::size_t>(dimensions);
  if (index >= total) return std::numeric_limits<double>::quiet_NaN();
  return values[index];
}

#include "selnorm-event.cc.inc"

struct ClusterNormalizerContext {
  int vector_rule = SELVECTOR_PRODUCT;
  const long double *quadrature_weights;
  const double *mean;
  const double *residual_sd;
  const double *loading;
  int dimension;
  const double *selection_se;
  const double *omega;
  int kernel_mode;
  SelNormKernelData selection;
};

double cluster_log_integrand(const ClusterNormalizerContext &context,
                             double gamma)
{
  double value = -0.5 * gamma * gamma;
  if (context.vector_rule != SELVECTOR_PRODUCT) {
    std::vector<double> means(context.dimension);
    for (int i = 0; i < context.dimension; ++i) {
      means[i] = context.mean[i] + context.loading[i] * gamma;
    }
    return value + independent_best_log_mass(means.data(), context.residual_sd,
      context.dimension, context.selection_se, context.omega,
      context.selection, context.vector_rule);
  }
  for (int i = 0; i < context.dimension; ++i) {
    const double conditional_mean =
      context.mean[i] + context.loading[i] * gamma;
    const double local = cpp_selnorm_step_log_norm(
      conditional_mean, context.residual_sd[i], context.selection_se[i],
      context.omega, context.selection, 1, false
    );
    if (!std::isfinite(local)) return -std::numeric_limits<double>::infinity();
    value += local;
  }
  return value;
}

// Optional projection state shared with the existing scalar/group rules.
struct ClusterRuleProjection {
  std::vector<double> means;
  std::vector<double> log_probability;
};


struct ClusterRuleWorkspace {
  std::vector<long double> products;
  std::vector<double> means;
  std::vector<long double> probability_cache;
  std::vector<double> projection_log_cache;
  std::vector<unsigned char> cache_ready;
  std::vector<long double> projection_probability;
};

double cluster_rule_log_integral(const ClusterNormalizerContext &context,
                                 const double *nodes,
                                 const double *log_weights,
                                 int offset, int order,
                                 ClusterRuleWorkspace *workspace = nullptr,
                                 const int *row_cache_slot = nullptr,
                                 int cache_slots = 0,
                                 ClusterRuleProjection *projection = nullptr)
{
  ClusterRuleWorkspace local;
  ClusterRuleWorkspace &scratch = workspace ? *workspace : local;
  std::vector<long double> &products = scratch.products;
  std::vector<double> &means = scratch.means;
  products.assign(static_cast<std::size_t>(order), 1.0L);
  means.resize(static_cast<std::size_t>(order));
  if (projection != nullptr) {
    projection->means.resize(static_cast<std::size_t>(context.dimension) * order);
    projection->log_probability.resize(static_cast<std::size_t>(context.dimension) * order);
    scratch.projection_probability.resize(order);
  }
  bool direct_ok = context.vector_rule == SELVECTOR_PRODUCT;
  const bool reuse = direct_ok && row_cache_slot != nullptr && cache_slots > 0;
  if (reuse) {
    scratch.probability_cache.resize(static_cast<std::size_t>(cache_slots) * order);
    if (projection != nullptr)
      scratch.projection_log_cache.resize(static_cast<std::size_t>(cache_slots) * order);
    scratch.cache_ready.assign(cache_slots, 0);
  }
  for (int i = 0; i < context.dimension && direct_ok; ++i) {
    const int slot = reuse ? row_cache_slot[i] : -1;
    const bool cached = slot >= 0 && scratch.cache_ready[slot];
    if (!cached || projection != nullptr) {
      for (int j = 0; j < order; ++j) {
        means[j] = context.mean[i] + context.loading[i] * nodes[offset + j];
      }
    }
    if (slot < 0 && projection == nullptr) {
      // No repetition: retain the original multiplication and fallback path.
      direct_ok = cpp_selnorm_step_normalizer_product(
        means.data(), means.size(), context.residual_sd[i],
        context.selection_se[i], context.omega, context.selection,
        products.data()
      );
    } else {
      long double *probability = slot < 0 ? scratch.projection_probability.data() :
        scratch.probability_cache.data() + static_cast<std::size_t>(slot) * order;
      if (!cached) {
        std::fill_n(probability, order, 1.0L);
        direct_ok = cpp_selnorm_step_normalizer_product(
          means.data(), means.size(), context.residual_sd[i],
          context.selection_se[i], context.omega, context.selection,
          probability
        );
        if (!direct_ok) break;
        if (slot >= 0) scratch.cache_ready[slot] = 1;
      }
      // Replay every original row in its original position. No power or
      // multiplicity regrouping changes the product's rounding order.
      for (int j = 0; j < order; ++j) {
        products[j] *= probability[j];
        if (projection != nullptr) {
          projection->means[static_cast<std::size_t>(i) * order + j] = means[j];
          // Equal row geometry already shares probability[]. Its logarithm
          // is identical too; retain the original product replay and row order.
          const double log_probability = cached ?
            scratch.projection_log_cache[static_cast<std::size_t>(slot) * order + j] :
            static_cast<double>(std::log(probability[j]));
          projection->log_probability[static_cast<std::size_t>(i) * order + j] = log_probability;
          if (slot >= 0 && !cached)
            scratch.projection_log_cache[static_cast<std::size_t>(slot) * order + j] = log_probability;
        }
        if (!(products[j] > 0.0L) || !std::isfinite(products[j])) {
          direct_ok = false;
          break;
        }
      }
    }
  }
  if (direct_ok) {
    long double integral = 0.0L;
    for (int j = 0; j < order; ++j) {
      const long double weight = context.quadrature_weights ?
        context.quadrature_weights[offset + j] :
        std::exp(static_cast<long double>(log_weights[offset + j]));
      integral += products[j] * weight;
    }
    if (integral > 0.0L && std::isfinite(integral)) {
      return static_cast<double>(std::log(integral));
    }
  }

  // A failed direct product discards the cache and uses the unchanged full
  // log-scale evaluator, including all rows and its existing tail safeguards.
  const double negative_infinity = -std::numeric_limits<double>::infinity();
  double out = negative_infinity;
  for (int j = 0; j < order; ++j) {
    const int index = offset + j;
    const double z = nodes[index];
    double term;
    if (projection == nullptr) {
      term = log_weights[index] + cluster_log_integrand(context, z) + 0.5 * z * z;
    } else {
      term = log_weights[index];
      for (int row = 0; row < context.dimension; ++row) {
        const double location = context.mean[row] + context.loading[row] * z;
        const double local = cpp_selnorm_step_log_norm(location, context.residual_sd[row],
          context.selection_se[row], context.omega, context.selection, 1, false);
        projection->means[static_cast<std::size_t>(row) * order + j] = location;
        projection->log_probability[static_cast<std::size_t>(row) * order + j] = local;
        term += local;
      }
    }
    out = log_add_exp(out, term);
  }
  return out;
}

double cluster_relative_change(double coarse, double fine)
{
  const double difference = coarse - fine;
  if (!std::isfinite(difference) || std::fabs(difference) > 50.0) {
    return std::numeric_limits<double>::infinity();
  }
  return std::fabs(std::expm1(difference));
}

double quadrature_log_gaussian_tail_bound(
    const double *nodes, int offset, int order, int rank,
    const double *omega, int n_bins, int dimension, int vector_rule)
{
  if (order < 1 || rank < 1) return std::numeric_limits<double>::infinity();
  const double *begin = nodes + offset;
  const double *end = begin + order;
  const double lower = *std::min_element(begin, end);
  const double upper = *std::max_element(begin, end);
  const double maximum_weight = *std::max_element(omega, omega + n_bins);
  if (!(maximum_weight > 0.0) || !std::isfinite(maximum_weight)) {
    return std::numeric_limits<double>::infinity();
  }

  // Successive Gaussian rules can agree while both miss a displaced selection
  // mode. Bound the entire omitted node-box region using the maximum possible
  // publication weight and a union bound over independent standard factors.
  // This is a coverage requirement in addition to the interior rule comparison;
  // it does not turn successive-rule differences into rigorous error bounds.
  const double log_tail = log_add_exp(
    pnorm(lower, 0.0, 1.0, true, true),
    pnorm(upper, 0.0, 1.0, false, true)
  ) + std::log(static_cast<double>(rank));
  const double log_weight_bound = std::log(maximum_weight) *
    static_cast<double>(vector_rule == SELVECTOR_PRODUCT ? dimension : 1);
  return log_tail + log_weight_bound;
}

bool quadrature_covers_gaussian_tails(
    const double *nodes, int offset, int order, int rank,
    const double *omega, int n_bins, int dimension, int vector_rule,
    double log_normalizer, double relative_tolerance)
{
  return std::isfinite(log_normalizer) &&
    quadrature_log_gaussian_tail_bound(nodes, offset, order, rank,
      omega, n_bins, dimension, vector_rule) - log_normalizer <=
      std::log(relative_tolerance);
}

struct SelNormCovarianceEnvelope {
  std::vector<int> type_index;
  int groups;
  double within_lower;
  double within_upper;
  double between_lower;
  double between_upper;
};

// One representable step encloses each correctly rounded basic operation.
// All interval operations using these have nonnegative operands. Handle exact
// zero products separately; underflowed positive products retain an upper
// bound of the least positive representable number. Overflow rejects the
// optional paths. These directed endpoints modify only auxiliary bounds.
// Shared by the covariance envelope and both covariance gap bounds so the
// outward-rounding convention has exactly one definition.
inline double mul_down(double x, double y)
{
  if (x == 0.0 || y == 0.0) return 0.0;
  const double value = x * y;
  return value == 0.0 ? 0.0 :
    std::nextafter(value, -std::numeric_limits<double>::infinity());
}

inline double mul_up(double x, double y)
{
  if (x == 0.0 || y == 0.0) return 0.0;
  return std::nextafter(x * y, std::numeric_limits<double>::infinity());
}

inline double add_down(double x, double y)
{
  if (x == 0.0) return y;
  if (y == 0.0) return x;
  return std::nextafter(x + y, -std::numeric_limits<double>::infinity());
}

inline double add_up(double x, double y)
{
  if (x == 0.0) return y;
  if (y == 0.0) return x;
  return std::nextafter(x + y, std::numeric_limits<double>::infinity());
}

// Equal means and diagonals with ordered off-diagonals order Gaussian
// expectations of nonnegative product weights monotone in one common direction.
// Along the covariance segment, the derivative is
// E[w_i' w_j' product_{l != i,j} w_l] >= 0; smooth monotone limits cover steps.
// The dense caller has already validated the target covariance by Cholesky.
// Positive residuals and common/type factors validate the auxiliary covariances.
bool covariance_envelope(const std::vector<double>& covariance, int dimension,
    const double* selection_se, SelNormCovarianceEnvelope* envelope)
{
  if (dimension < 2 || selection_se == nullptr || envelope == nullptr ||
      covariance.size() != static_cast<std::size_t>(dimension) * dimension) {
    return false;
  }
  const double infinity = std::numeric_limits<double>::infinity();

  std::vector<double> correlation(static_cast<std::size_t>(dimension) * dimension);
  std::vector<double> values;
  double common_lower = 0.0;
  double common_upper = infinity;
  for (int i = 0; i < dimension; ++i) {
    const double diagonal = covariance[i + dimension * i];
    if (!(selection_se[i] > 0.0) || !std::isfinite(selection_se[i]) ||
        !(diagonal > 0.0) || !std::isfinite(diagonal)) return false;
  }
  for (int j = 0; j < dimension; ++j) {
    for (int i = j + 1; i < dimension; ++i) {
      const double value = covariance[i + dimension * j];
      if (!(value >= 0.0) || !std::isfinite(value) ||
          value != covariance[j + dimension * i]) return false;
      const double scale = selection_se[i] * selection_se[j];
      if (!(scale > 0.0) || !std::isfinite(scale)) return false;
      const double ratio = value / scale;
      if (!std::isfinite(ratio)) return false;
      const double scale_lower = mul_down(selection_se[i], selection_se[j]);
      const double scale_upper = mul_up(selection_se[i], selection_se[j]);
      if (!(scale_lower > 0.0) || !std::isfinite(scale_upper)) return false;
      common_lower = std::max(common_lower, value == 0.0 ? 0.0 :
        std::nextafter(value / scale_upper, -infinity));
      common_upper = std::min(common_upper, value == 0.0 ? 0.0 :
        std::nextafter(value / scale_lower, infinity));
      correlation[i + dimension * j] = ratio;
      correlation[j + dimension * i] = ratio;
      values.push_back(ratio);
    }
  }

  // A largest-gap partition only proposes an efficient common/type factor
  // pattern. It never establishes an equality or discards a matrix residual.
  // Every off-diagonal is checked independently below. A failed or expensive
  // partition falls back to a one-factor covariance envelope.
  std::sort(values.begin(), values.end());
  double gap = 0.0;
  double split = values.front();
  for (std::size_t i = 1; i < values.size(); ++i) {
    const double difference = values[i] - values[i - 1];
    if (difference > gap) {
      gap = difference;
      split = values[i - 1];
    }
  }
  std::vector<int> types(static_cast<std::size_t>(dimension), -1);
  int groups = 0;
  // Overlapping arithmetic enclosures favour the cheaper compound proposal;
  // the final envelopes still retain and check the complete observed spread.
  if (gap > 0.0 && common_lower > common_upper) {
    for (int root = 0; root < dimension; ++root) {
      if (types[root] >= 0) continue;
      std::vector<int> pending(1, root);
      types[root] = groups;
      while (!pending.empty()) {
        const int row = pending.back();
        pending.pop_back();
        for (int column = 0; column < dimension; ++column) {
          if (column != row && types[column] < 0 &&
              correlation[row + dimension * column] > split) {
            types[column] = groups;
            pending.push_back(column);
          }
        }
      }
      ++groups;
    }
    for (int j = 0; j < dimension && groups <= 3; ++j) {
      for (int i = j + 1; i < dimension; ++i) {
        if ((types[i] == types[j]) !=
            (correlation[i + dimension * j] > split)) {
          groups = 4;
          break;
        }
      }
    }
  }
  if (groups < 2 || groups > 3) {
    groups = 1;
    std::fill(types.begin(), types.end(), 0);
  }

  // Trying the compound envelope second also handles an unsuitable proposed
  // type partition. Its potentially greater integration width is exposed to
  // the caller's unchanged numerical diagnostic, not silently accepted.
  for (int attempt = 0; attempt < 2; ++attempt) {
    double within_lower = infinity;
    double within_upper = 0.0;
    double between_lower = infinity;
    double between_upper = 0.0;
    for (int j = 0; j < dimension; ++j) {
      for (int i = j + 1; i < dimension; ++i) {
        const double value = covariance[i + dimension * j];
        const double scale_lower = mul_down(selection_se[i], selection_se[j]);
        const double scale_upper = mul_up(selection_se[i], selection_se[j]);
        if (!(scale_lower > 0.0) || !std::isfinite(scale_upper)) return false;
        const double lower = value == 0.0 ? 0.0 :
          std::nextafter(value / scale_upper, -infinity);
        const double upper = value == 0.0 ? 0.0 :
          std::nextafter(value / scale_lower, infinity);
        if (!(lower >= 0.0) || !std::isfinite(upper)) return false;
        if (types[i] == types[j]) {
          within_lower = std::min(within_lower, lower);
          within_upper = std::max(within_upper, upper);
        } else {
          between_lower = std::min(between_lower, lower);
          between_upper = std::max(between_upper, upper);
        }
      }
    }
    if (groups == 1) between_lower = between_upper = 0.0;
    bool valid = std::isfinite(within_lower) &&
      std::isfinite(between_lower) && between_lower <= within_lower &&
      between_upper <= within_upper;

    // The coefficient-defined envelope is already outward. Rounded square
    // roots and loading products can widen its numerical representation;
    // verify those products too. Move auxiliary coefficients outward by one
    // representable step per pass and accept only a verified result. The
    // bounded pass count is a work limit, never a numerical tolerance.
    for (int side = 0; side < 2 && valid; ++side) {
      double& within = side == 0 ? within_lower : within_upper;
      double& between = side == 0 ? between_lower : between_upper;
      bool certified = false;
      for (int pass = 0; pass < 32 && !certified; ++pass) {
        if (!(within >= between) || !(between >= 0.0) ||
            !std::isfinite(within)) break;
        const double common_root = std::sqrt(between);
        const double type_root = std::sqrt(within - between);
        std::vector<double> common(static_cast<std::size_t>(dimension));
        std::vector<double> type(static_cast<std::size_t>(dimension));
        certified = true;
        for (int i = 0; i < dimension; ++i) {
          common[i] = common_root * selection_se[i];
          type[i] = type_root * selection_se[i];
          const double nominal_variance = mul_up(
            mul_up(within, selection_se[i]), selection_se[i]);
          const double loading_variance = add_up(
            mul_up(common[i], common[i]), mul_up(type[i], type[i]));
          if (!std::isfinite(nominal_variance) ||
              !std::isfinite(loading_variance) ||
              !(covariance[i + dimension * i] > nominal_variance) ||
              !(covariance[i + dimension * i] > loading_variance)) {
            valid = false;
            certified = false;
            break;
          }
        }
        for (int j = 0; j < dimension && valid; ++j) {
          for (int i = j + 1; i < dimension; ++i) {
            const double coefficient = types[i] == types[j] ? within : between;
            const double nominal = side == 0 ?
              mul_up(mul_up(coefficient, selection_se[i]), selection_se[j]) :
              mul_down(mul_down(coefficient, selection_se[i]), selection_se[j]);
            double represented = side == 0 ?
              mul_up(common[i], common[j]) : mul_down(common[i], common[j]);
            if (types[i] == types[j]) {
              represented = side == 0 ?
                add_up(represented, mul_up(type[i], type[j])) :
                add_down(represented, mul_down(type[i], type[j]));
            }
            const double supplied = covariance[i + dimension * j];
            if (side == 0 ? (nominal > supplied || represented > supplied) :
                            (nominal < supplied || represented < supplied)) {
              certified = false;
            }
          }
        }
        if (!valid || certified) break;
        const double direction = side == 0 ? -infinity : infinity;
        if (within != 0.0 || side == 1) within = std::nextafter(within, direction);
        if (groups > 1 && (between != 0.0 || side == 1)) {
          between = std::nextafter(between, direction);
        }
      }
      valid = valid && certified;
    }
    if (valid) {
      envelope->type_index = types;
      envelope->groups = groups;
      envelope->within_lower = within_lower;
      envelope->within_upper = within_upper;
      envelope->between_lower = between_lower;
      envelope->between_upper = between_upper;
      return true;
    }
    if (groups == 1) break;
    groups = 1;
    std::fill(types.begin(), types.end(), 0);
  }
  return false;
}

struct EnvelopeGroupGeometry {
  std::vector<int> row_index;
  std::vector<double> mean;
  std::vector<double> common_loading;
  std::vector<double> loading;
  std::vector<double> residual_sd;
  std::vector<double> selection_se;
  std::vector<int> row_cache_slot;
  int cache_slots = 0;
  bool singleton_child = false;
};

struct EnvelopeRuleProjection {
  std::vector<std::vector<double>> means;
  std::vector<std::vector<double>> log_coefficient;
  double log_scale = 0.0; // context GH weight / full A / number of pooled rows
};

void append_cluster_projection(const EnvelopeGroupGeometry &group,
    const ClusterRuleProjection &rule, const double *log_weights, int offset,
    int order, double log_scale, EnvelopeRuleProjection *projection,
    std::vector<std::vector<int>> *common_indices = nullptr, int common = -1)
{
  const int size = group.mean.size();
  std::vector<double> prefix(size + 1), suffix(size + 1);
  for (int point = 0; point < order; ++point) {
    prefix[0] = 0.0;
    suffix[size] = 0.0;
    for (int row = 0; row < size; ++row) {
      prefix[row + 1] = prefix[row] + rule.log_probability[row * order + point];
    }
    for (int row = size - 1; row >= 0; --row) {
      suffix[row] = suffix[row + 1] + rule.log_probability[row * order + point];
    }
    for (int row = 0; row < size; ++row) {
      const int original = group.row_index[row];
      projection->means[original].push_back(rule.means[row * order + point]);
      projection->log_coefficient[original].push_back(
        log_scale + log_weights[offset + point] + prefix[row] + suffix[row + 1]);
      if (common_indices != nullptr) (*common_indices)[original].push_back(common);
    }
  }
}

struct EnvelopeRuleWorkspace {
  ClusterRuleWorkspace cluster;
  std::vector<double> conditional_mean;
  std::vector<double> log_products;
  std::vector<long double> singleton_probability;
};

// Repeat only the failure conditions of prepare_envelope_geometry for one
// envelope side without materializing its geometry. The residual arithmetic
// matches that function exactly, so feasibility here agrees with its success.
bool envelope_geometry_feasible(
    const std::vector<double> &covariance, const double* selection_se, int dimension,
    const SelNormCovarianceEnvelope &envelope, bool upper)
{
  const double within = upper ? envelope.within_upper : envelope.within_lower;
  const double between = upper ? envelope.between_upper : envelope.between_lower;
  const double common = std::sqrt(envelope.groups == 1 ? within : between);
  const double child = envelope.groups == 1 ? 0.0 : std::sqrt(within - between);
  const bool absorb_singleton = envelope.groups > 1;
  std::vector<int> group_sizes;
  if (absorb_singleton) {
    group_sizes.assign(envelope.groups, 0);
    for (int row = 0; row < dimension; ++row) ++group_sizes[envelope.type_index[row]];
  }
  for (int row = 0; row < dimension; ++row) {
    if (!(selection_se[row] > 0.0) || !std::isfinite(selection_se[row]) ||
        !(covariance[row + dimension * row] > 0.0) ||
        !std::isfinite(covariance[row + dimension * row])) return false;
    const bool singleton = absorb_singleton && group_sizes[envelope.type_index[row]] == 1;
    const double common_loading = common * selection_se[row];
    const double child_loading = singleton ? 0.0 : child * selection_se[row];
    const long double residual = static_cast<long double>(covariance[row + dimension * row]) -
      static_cast<long double>(common_loading) * common_loading -
      static_cast<long double>(child_loading) * child_loading;
    if (!(residual > 0.0) || !std::isfinite(residual)) return false;
    const double sd = std::sqrt(static_cast<double>(residual));
    if (!(sd > 0.0) || !std::isfinite(sd)) return false;
  }
  return true;
}

bool prepare_envelope_geometry(
    const std::vector<double> &covariance, const double *mean,
    const double *selection_se, int dimension, const SelNormCovarianceEnvelope &envelope,
    bool upper, std::vector<EnvelopeGroupGeometry> *geometry)
{
  const double within = upper ? envelope.within_upper : envelope.within_lower;
  const double between = upper ? envelope.between_upper : envelope.between_lower;
  const double common = std::sqrt(envelope.groups == 1 ? within : between);
  const double child = envelope.groups == 1 ? 0.0 : std::sqrt(within - between);
  const bool absorb_singleton = envelope.groups > 1;
  std::vector<int> group_sizes;
  if (absorb_singleton) {
    group_sizes.assign(envelope.groups, 0);
    for (int row = 0; row < dimension; ++row) ++group_sizes[envelope.type_index[row]];
  }
  geometry->resize(envelope.groups);
  for (int row = 0; row < dimension; ++row) {
    const bool singleton = absorb_singleton && group_sizes[envelope.type_index[row]] == 1;
    const double common_loading = common * selection_se[row];
    const double child_loading = singleton ? 0.0 : child * selection_se[row];
    // A unique-support child contributes only to this diagonal. Integrating it
    // analytically leaves the authoritative diagonal minus the common square.
    const long double residual = static_cast<long double>(covariance[row + dimension * row]) -
      static_cast<long double>(common_loading) * common_loading -
      static_cast<long double>(child_loading) * child_loading;
    if (!(residual > 0.0) || !std::isfinite(residual)) return false;
    const double sd = std::sqrt(static_cast<double>(residual));
    if (!(sd > 0.0) || !std::isfinite(sd)) return false;
    EnvelopeGroupGeometry &group = (*geometry)[envelope.type_index[row]];
    group.row_index.push_back(row);
    group.mean.push_back(mean[row]);
    group.common_loading.push_back(common_loading);
    group.loading.push_back(envelope.groups == 1 ? common_loading : child_loading);
    group.residual_sd.push_back(sd);
    group.selection_se.push_back(selection_se[row]);
    group.singleton_child = singleton;
  }
  // Exact complete geometry, including signed zero. Fixed cutoffs, weights,
  // direction, and kernel belong to the common call, never a sampled key.
  const auto equal = [](double x, double y) {
    return x == y && std::signbit(x) == std::signbit(y);
  };
  for (EnvelopeGroupGeometry &group : *geometry) {
    const int size = group.mean.size();
    group.row_cache_slot.assign(size, -1);
    for (int row = 1; row < size; ++row) {
      for (int previous = 0; previous < row; ++previous) {
        if (equal(group.mean[row], group.mean[previous]) &&
            equal(group.common_loading[row], group.common_loading[previous]) &&
            equal(group.loading[row], group.loading[previous]) &&
            equal(group.residual_sd[row], group.residual_sd[previous]) &&
            equal(group.selection_se[row], group.selection_se[previous])) {
          if (group.row_cache_slot[previous] < 0) {
            group.row_cache_slot[previous] = group.cache_slots++;
          }
          group.row_cache_slot[row] = group.row_cache_slot[previous];
          break;
        }
      }
    }
    if (group.cache_slots == 0) group.row_cache_slot.clear();
  }
  return true;
}

double envelope_rule_log_integral(
    const std::vector<EnvelopeGroupGeometry> &geometry, const double *omega,
    const SelNormKernelData &selection, const double *nodes, const double *log_weights,
    const long double *weights, int offset, int order, EnvelopeRuleWorkspace *workspace,
    EnvelopeRuleProjection *projection = nullptr)
{
  ClusterNormalizerContext context;
  context.vector_rule = SELVECTOR_PRODUCT;
  context.quadrature_weights = weights;
  context.omega = omega;
  context.kernel_mode = SELKERNEL_STEP;
  context.selection = selection;
  if (geometry.size() == 1) {
    const EnvelopeGroupGeometry &group = geometry[0];
    context.dimension = group.mean.size();
    context.mean = group.mean.data();
    context.residual_sd = group.residual_sd.data();
    context.loading = group.loading.data();
    context.selection_se = group.selection_se.data();
    ClusterRuleProjection row_rule;
    const double result = cluster_rule_log_integral(context, nodes, log_weights, offset, order,
      &workspace->cluster, group.row_cache_slot.empty() ? nullptr : group.row_cache_slot.data(),
      group.cache_slots, projection == nullptr ? nullptr : &row_rule);
    if (projection != nullptr && std::isfinite(result)) {
      append_cluster_projection(group, row_rule, log_weights, offset, order,
        projection->log_scale, projection);
    }
    return result;
  }

  std::vector<double> &log_products = workspace->log_products;
  std::vector<double> &conditional_mean = workspace->conditional_mean;
  log_products.assign(order, 0.0);
  std::vector<double> group_integrals;
  std::vector<std::size_t> starts;
  std::vector<std::vector<int>> common_indices;
  if (projection != nullptr) {
    group_integrals.resize(geometry.size() * order);
    starts.resize(projection->means.size());
    common_indices.resize(projection->means.size());
    for (std::size_t row = 0; row < starts.size(); ++row) starts[row] = projection->means[row].size();
  }
  for (std::size_t g = 0; g < geometry.size(); ++g) {
    const EnvelopeGroupGeometry &group = geometry[g];
    const int size = group.mean.size();
    if (group.singleton_child) {
      std::vector<double> &means = workspace->cluster.means;
      std::vector<long double> &probability = workspace->singleton_probability;
      means.resize(order);
      probability.assign(order, 1.0L);
      for (int point = 0; point < order; ++point) {
        means[point] = group.mean[0] + group.common_loading[0] * nodes[offset + point];
      }
      const bool direct = cpp_selnorm_step_normalizer_product(
        means.data(), means.size(), group.residual_sd[0], group.selection_se[0],
        omega, selection, probability.data());
      for (int point = 0; point < order; ++point) {
        const double local = direct ? static_cast<double>(std::log(probability[point])) :
          cpp_selnorm_step_log_norm(means[point], group.residual_sd[0],
            group.selection_se[0], omega, selection, 1, false);
        log_products[point] += local;
        if (projection != nullptr) {
          if (!std::isfinite(local)) return local;
          group_integrals[g * order + point] = local;
          const int row = group.row_index[0];
          projection->means[row].push_back(means[point]);
          projection->log_coefficient[row].push_back(projection->log_scale + log_weights[offset + point]);
          common_indices[row].push_back(point);
        }
      }
      continue;
    }
    conditional_mean.resize(size);
    context.dimension = size;
    context.mean = conditional_mean.data();
    context.residual_sd = group.residual_sd.data();
    context.loading = group.loading.data();
    context.selection_se = group.selection_se.data();
    ClusterRuleProjection row_rule;
    for (int point = 0; point < order; ++point) {
      for (int row = 0; row < size; ++row) {
        conditional_mean[row] = group.mean[row] +
          group.common_loading[row] * nodes[offset + point];
      }
      const double local = cluster_rule_log_integral(
        context, nodes, log_weights, offset, order, &workspace->cluster,
        group.row_cache_slot.empty() ? nullptr : group.row_cache_slot.data(), group.cache_slots,
        projection == nullptr ? nullptr : &row_rule);
      log_products[point] += local;
      if (projection != nullptr) {
        if (!std::isfinite(local)) return local;
        group_integrals[g * order + point] = local;
        append_cluster_projection(group, row_rule, log_weights, offset, order,
          projection->log_scale + log_weights[offset + point], projection, &common_indices, point);
      }
    }
  }
  if (projection != nullptr) {
    for (std::size_t g = 0; g < geometry.size(); ++g) {
      for (int row : geometry[g].row_index) {
        for (std::size_t component = starts[row]; component < projection->means[row].size(); ++component) {
          const int common = common_indices[row][component - starts[row]];
          for (std::size_t other = 0; other < geometry.size(); ++other) {
            if (other != g) projection->log_coefficient[row][component] += group_integrals[other * order + common];
          }
        }
      }
    }
  }
  double result = -std::numeric_limits<double>::infinity();
  for (int point = 0; point < order; ++point) {
    result = log_add_exp(result, log_weights[offset + point] + log_products[point]);
  }
  return result;
}

// Absolute Price-theorem bound |A(covariance) - A(lower envelope)|.
// The caller has certified equal diagonals, 0 <= lower_ij <= covariance_ij,
// positive definite covariance, and nonnegative product step weights.
// Off-diagonal differentiation has coefficient one (Price's theorem,
// Voigtlaender, 2021, Theorem 1, doi:10.1007/s10959-020-01017-w):
// each derivative measure has total variation
// sum |omega[b] - omega[b-1]|. Bound the remaining weights by their maximum
// and each bivariate density by its largest prefactor along the segment.
// Every actual lower loading uses the same double arithmetic as quadrature.
double covariance_price_gap_bound(const std::vector<double>& covariance,
    const double* selection_se, int dimension,
    const SelNormCovarianceEnvelope& envelope, const double* omega, int n_bins,
    bool monotone)
{
  const double infinity = std::numeric_limits<double>::infinity();
  if (dimension < 2 || n_bins < 1 || selection_se == nullptr || omega == nullptr ||
      covariance.size() != static_cast<std::size_t>(dimension) * dimension ||
      envelope.type_index.size() != static_cast<std::size_t>(dimension) ||
      envelope.groups < 1 || envelope.groups > 3) return infinity;
  double minimum = infinity, maximum = 0.0;
  for (int bin = 0; bin < n_bins; ++bin) {
    if (!(omega[bin] >= 0.0) || !std::isfinite(omega[bin])) return infinity;
    minimum = std::min(minimum, omega[bin]);
    maximum = std::max(maximum, omega[bin]);
  }
  if (!(maximum > 0.0)) return infinity;
  if (minimum == maximum) return 0.0;
  const double common_variance = envelope.groups == 1 ?
    envelope.within_lower : envelope.between_lower;
  const double child_variance = envelope.groups == 1 ? 0.0 :
    envelope.within_lower - envelope.between_lower;
  if (!(common_variance >= 0.0) || !(child_variance >= 0.0) ||
      !std::isfinite(common_variance) || !std::isfinite(child_variance)) return infinity;
  const double common = std::sqrt(common_variance);
  const double child = std::sqrt(child_variance);
  std::vector<double> common_loading(dimension), child_loading(dimension);
  for (int i = 0; i < dimension; ++i) {
    if (!(selection_se[i] > 0.0) || !std::isfinite(selection_se[i]) ||
        envelope.type_index[i] < 0 || envelope.type_index[i] >= envelope.groups ||
        !(covariance[i + dimension * i] > 0.0) ||
        !std::isfinite(covariance[i + dimension * i])) return infinity;
    common_loading[i] = common * selection_se[i];
    child_loading[i] = child * selection_se[i];
    if (!std::isfinite(common_loading[i]) || !std::isfinite(child_loading[i])) return infinity;
  }
  // This binary64 constant is strictly below mathematical 2*pi. Downward
  // products/square roots preserve a lower density-prefactor denominator.
  const double two_pi_lower = 6.283185307179586;
  double sum = 0.0;
  for (int j = 0; j < dimension; ++j) {
    for (int i = j + 1; i < dimension; ++i) {
      const double value = covariance[i + dimension * j];
      if (!(value >= 0.0) || !std::isfinite(value) ||
          value != covariance[j + dimension * i]) return infinity;
      double lower = mul_down(common_loading[i], common_loading[j]);
      if (envelope.type_index[i] == envelope.type_index[j]) {
        lower = add_down(lower, mul_down(child_loading[i], child_loading[j]));
      }
      if (!(lower <= value)) return infinity;
      if (lower == value) continue;
      const double gap = std::nextafter(value - lower, infinity);
      const double determinant = std::nextafter(
        mul_down(covariance[i + dimension * i], covariance[j + dimension * j]) -
        mul_up(value, value), -infinity);
      if (!(determinant > 0.0) || !std::isfinite(determinant)) return infinity;
      const double root = std::nextafter(std::sqrt(determinant), -infinity);
      const double denominator = mul_down(two_pi_lower, root);
      if (!(denominator > 0.0) || !std::isfinite(denominator)) return infinity;
      const double term = std::nextafter(gap / denominator, infinity);
      sum = add_up(sum, term);
      if (!std::isfinite(sum)) return infinity;
    }
  }
  if (sum == 0.0) return 0.0;
  // Both branches compute the same quantity: for a monotone weight sequence
  // the total variation sum |omega[b] - omega[b-1]| equals max - min exactly.
  // The closed form is kept only because it is one rounding step, whereas the
  // summed form accumulates outward rounding across every bin and is therefore
  // slightly looser. Dropping either branch changes only that slack, so do not
  // collapse them without accepting the wider monotone bound.
  double variation = std::nextafter(maximum - minimum, infinity);
  if (!monotone) {
    variation = 0.0;
    for (int bin = 1; bin < n_bins; ++bin) {
      if (omega[bin] == omega[bin - 1]) continue;
      variation = add_up(variation,
        std::nextafter(std::fabs(omega[bin] - omega[bin - 1]), infinity));
    }
  }
  double multiplier = mul_up(variation, variation);
  for (int i = 0; i < dimension - 2; ++i) multiplier = mul_up(multiplier, maximum);
  const double result = mul_up(sum, multiplier);
  return std::isfinite(result) ? result : infinity;
}

// Relative covariance perturbation for strictly positive product weights.
// The auxiliary covariance C0 = D + U U' uses the emitted SDs/loadings,
// including their diagonal rounding. Since C0 >= D, the maximum absolute row
// sum of D^(-1/2) (C-C0) D^(-1/2) bounds the whitened spectral error epsilon.
// Along the Gaussian covariance segment, Jensen and the chi-square MGF give
// E_selected[||Z||^2] <= 4*k*log(M/m) + 2*k*log(2). Integrating the Gaussian
// log-density derivative bounds |log A(C) - log A(C0)| without dividing an
// absolute error by a possibly tiny A. This is the bound used by zplot context
// integration. Zero weights retain the absolute Price bound instead.
double covariance_relative_gap_bound(const std::vector<double>& covariance,
    int dimension, const std::vector<EnvelopeGroupGeometry>& geometry,
    const double* omega, int n_bins)
{
  const double infinity = std::numeric_limits<double>::infinity();
  // Re-validate every input rather than relying on the caller's ordering:
  // the loops below index fixed-size vectors by geometry-supplied row indices.
  if (dimension < 2 || n_bins < 1 || omega == nullptr ||
      covariance.size() != static_cast<std::size_t>(dimension) * dimension ||
      geometry.empty()) return infinity;
  double minimum = infinity, maximum = 0.0;
  for (int bin = 0; bin < n_bins; ++bin) {
    if (!(omega[bin] > 0.0) || !std::isfinite(omega[bin])) return infinity;
    minimum = std::min(minimum, omega[bin]);
    maximum = std::max(maximum, omega[bin]);
  }
  if (!(minimum > 0.0) || !std::isfinite(minimum)) return infinity;
  if (minimum == maximum) return 0.0;
  // One group means the single common factor carries the whole loading, so
  // there is no separate child contribution. This matches the sibling bound's
  // envelope.groups == 1 test: prepare_envelope_geometry() resizes the
  // geometry to envelope.groups, so the two predicates are the same.
  const bool single_group = geometry.size() == 1;
  std::vector<double> sd(dimension), common(dimension), child(dimension);
  std::vector<int> groups(dimension);
  for (std::size_t group = 0; group < geometry.size(); ++group) {
    const EnvelopeGroupGeometry& part = geometry[group];
    if (part.residual_sd.size() != part.row_index.size() ||
        part.common_loading.size() != part.row_index.size() ||
        part.loading.size() != part.row_index.size()) return infinity;
    for (std::size_t local = 0; local < part.row_index.size(); ++local) {
      const int row = part.row_index[local];
      if (row < 0 || row >= dimension) return infinity;
      sd[row] = part.residual_sd[local];
      common[row] = part.common_loading[local];
      child[row] = single_group ? 0.0 : part.loading[local];
      groups[row] = static_cast<int>(group);
    }
  }
  std::vector<double> row_sum(dimension, 0.0);
  for (int column = 0; column < dimension; ++column) {
    for (int row = column; row < dimension; ++row) {
      double lower = row == column ? mul_down(sd[row], sd[row]) : 0.0;
      double upper = row == column ? mul_up(sd[row], sd[row]) : 0.0;
      lower = add_down(lower, mul_down(common[row], common[column]));
      upper = add_up(upper, mul_up(common[row], common[column]));
      if (groups[row] == groups[column]) {
        lower = add_down(lower, mul_down(child[row], child[column]));
        upper = add_up(upper, mul_up(child[row], child[column]));
      }
      const double value = covariance[row + dimension * column];
      // C0 = D + U U' can only represent nonnegative off-diagonals, so a
      // negative supplied entry is outside the envelope family. Reject it here
      // as the absolute bound does, instead of relying on epsilon reaching one.
      if (!std::isfinite(value) || (row != column && !(value >= 0.0))) return infinity;
      const double difference = std::nextafter(std::max(
        std::fabs(value - lower), std::fabs(value - upper)), infinity);
      const double scale = mul_down(sd[row], sd[column]);
      if (!(scale > 0.0) || !std::isfinite(scale) || !std::isfinite(difference)) return infinity;
      const double error = std::nextafter(difference / scale, infinity);
      row_sum[row] = add_up(row_sum[row], error);
      if (row != column) row_sum[column] = add_up(row_sum[column], error);
    }
  }
  const double epsilon = *std::max_element(row_sum.begin(), row_sum.end());
  if (!(epsilon >= 0.0) || !(epsilon < 1.0)) return infinity;
  const double complement = std::nextafter(1.0 - epsilon, -infinity);
  if (!(complement > 0.0)) return infinity;
  const double log_ratio = std::nextafter(
    std::nextafter(std::log(maximum), infinity) -
    std::nextafter(std::log(minimum), -infinity), infinity);
  const double log_two = std::nextafter(std::log(2.0), infinity);
  const double moment = add_up(mul_up(4.0 * dimension, log_ratio),
                              mul_up(2.0 * dimension, log_two));
  const double beta = std::nextafter(epsilon / (2.0 * complement), infinity);
  const double determinant = mul_up(0.5 * dimension,
    -std::nextafter(std::log1p(-epsilon), -infinity));
  const double log_bound = add_up(determinant, mul_up(beta, moment));
  if (!std::isfinite(log_bound)) return infinity;
  return std::nextafter(std::expm1(log_bound), infinity);
}

// Selection-envelope memoization. Hashes only select candidates: complete
// owned values, including the numerical controls, establish every cache hit.
struct EnvelopeMemoSpan {
  const double *data;
  std::size_t size;
};

struct EnvelopeMemoView {
  const EnvelopeMemoSpan *parts;
  std::size_t count;
  std::size_t size() const
  {
    std::size_t result = 0;
    for (std::size_t i = 0; i < count; ++i) result += parts[i].size;
    return result;
  }
  std::size_t hash() const
  {
    std::size_t value = 0;
    for (std::size_t i = 0; i < count; ++i) {
      const std::string_view bytes(reinterpret_cast<const char *>(parts[i].data),
        parts[i].size * sizeof(double));
      value ^= std::hash<std::string_view>{}(bytes) +
        static_cast<std::size_t>(0x9e3779b97f4a7c15ULL) +
        (value << 6) + (value >> 2) + parts[i].size;
    }
    return value;
  }
  bool equals(const std::pmr::vector<double> &stored) const
  {
    if (stored.size() != size()) return false;
    std::size_t offset = 0;
    for (std::size_t i = 0; i < count; ++i) {
      if (std::memcmp(stored.data() + offset, parts[i].data,
                      parts[i].size * sizeof(double)) != 0) return false;
      offset += parts[i].size;
    }
    return true;
  }
  void copy_to(std::pmr::vector<double> *stored) const
  {
    stored->reserve(size());
    for (std::size_t i = 0; i < count; ++i) {
      stored->insert(stored->end(), parts[i].data, parts[i].data + parts[i].size);
    }
  }
};

// Keys, controls, LRU nodes, and index blocks use this direct counting
// resource. A pool resource is deliberately not interposed here: its hidden
// pool/oversize state makes the byte cap indirect, and the actual crash was in
// a key copy on its exception-driven full-cache retry path. Direct requests
// are exactly countable and every eviction immediately returns usable bytes.
// Process RSS additionally includes allocator metadata outside these blocks.
class EnvelopeMemoResource final : public std::pmr::memory_resource {
  void *do_allocate(std::size_t bytes, std::size_t alignment) override
  {
    if (bytes > capacity - allocated) throw std::bad_alloc();
    void *result = std::pmr::new_delete_resource()->allocate(bytes, alignment);
    allocated += bytes;
    peak = std::max(peak, allocated);
    return result;
  }
  void do_deallocate(void *value, std::size_t bytes, std::size_t alignment) override
  {
    std::pmr::new_delete_resource()->deallocate(value, bytes, alignment);
    allocated -= bytes;
  }
  bool do_is_equal(const std::pmr::memory_resource &other) const noexcept override
  {
    return this == &other;
  }
public:
  std::size_t capacity, allocated = 0, peak = 0;
  explicit EnvelopeMemoResource(std::size_t bytes) : capacity(bytes) {}
};

struct EnvelopeMemoResult {
  bool success = false, used_envelope = false;
  double value = 0.0, change = 0.0, width = 0.0, tail = 0.0;
};

enum class EnvelopeMemoComponent : unsigned char { exact = 0, coarse = 1 };

// Portable fieldwise snapshot framing; double bytes are bound to this ABI.
// FNV-1a detects accidental corruption, not deliberate modification/authenticity.
constexpr std::uint64_t snapshot_magic = UINT64_C(0x315041434d424f52);
constexpr std::size_t snapshot_header = 48, snapshot_record = 48;
std::uint64_t snapshot_checksum(const unsigned char *data, std::size_t size)
{
  std::uint64_t value = UINT64_C(14695981039346656037);
  for (std::size_t i = 0; i < size; ++i) value = (value ^ data[i]) * UINT64_C(1099511628211);
  return value;
}
std::uint64_t snapshot_abi()
{
  const std::uint16_t one = 1;
  std::uint64_t value = sizeof(double) | (std::uint64_t(sizeof(long double)) << 8) |
    (std::uint64_t(std::numeric_limits<long double>::digits) << 16) |
    (std::uint64_t(sizeof(std::size_t)) << 24) | (std::uint64_t(sizeof(int)) << 32) |
    (std::uint64_t(*reinterpret_cast<const unsigned char *>(&one)) << 40);
#if defined(__MINGW32__) && defined(__GNUC__) && defined(__x86_64__)
  value |= std::uint64_t(selnorm_fma_detail::supported()) << 48;
#endif
  return value;
}
struct SnapshotWriter {
  unsigned char *data; std::size_t offset = 0;
  void integer(std::uint64_t value) {
    for (int i = 0; i < 8; ++i) data[offset++] = static_cast<unsigned char>(value >> (8 * i));
  }
  void bytes(const void *value, std::size_t size) {
    if (size) std::memcpy(data + offset, value, size);
    offset += size;
  }
  void real(double value) { bytes(&value, sizeof(value)); }
};
struct SnapshotReader {
  const unsigned char *data; std::size_t size, offset = 0;
  const unsigned char *bytes(std::size_t count) {
    if (count > size - offset) throw std::invalid_argument("Selection cache snapshot is truncated.");
    const unsigned char *value = data + offset; offset += count; return value;
  }
  std::uint64_t integer() {
    const unsigned char *value = bytes(8); std::uint64_t out = 0;
    for (int i = 0; i < 8; ++i) out |= std::uint64_t(value[i]) << (8 * i);
    return out;
  }
  std::size_t count() {
    const std::uint64_t value = integer();
    if (value > std::numeric_limits<std::size_t>::max()) throw std::invalid_argument("Selection cache snapshot size is unsupported.");
    return static_cast<std::size_t>(value);
  }
  const unsigned char *reals(std::size_t count) {
    if (count > (size - offset) / sizeof(double)) throw std::invalid_argument("Selection cache snapshot array is truncated.");
    return bytes(count * sizeof(double));
  }
  double real() { double value; std::memcpy(&value, bytes(sizeof(value)), sizeof(value)); return value; }
};
double snapshot_real(const unsigned char *data, std::size_t index)
{
  double value; std::memcpy(&value, data + index * sizeof(value), sizeof(value)); return value;
}
std::size_t snapshot_count(double value, std::size_t minimum)
{
  if (!std::isfinite(value) || value < minimum || value > std::numeric_limits<int>::max() || value != std::floor(value))
    throw std::invalid_argument("Selection cache snapshot has an invalid dimension.");
  return static_cast<std::size_t>(value);
}
void snapshot_controls(const unsigned char *data, std::size_t count)
{
  if (count == 0) return;
  if (count < 2) throw std::invalid_argument("Selection cache snapshot controls are incomplete.");
  const std::size_t rules = snapshot_count(snapshot_real(data, 0), 3);
  const std::size_t nodes = snapshot_count(snapshot_real(data, 1), 1);
  if (rules > count - 2 || nodes > (count - 2 - rules) / 2 || count != 2 + rules + 2 * nodes)
    throw std::invalid_argument("Selection cache snapshot controls have invalid lengths.");
  std::size_t total = 0, previous = 0;
  for (std::size_t i = 0; i < rules; ++i) {
    const std::size_t order = snapshot_count(snapshot_real(data, 2 + i), 1);
    if (order <= previous || order > nodes - total) throw std::invalid_argument("Selection cache snapshot quadrature orders are invalid.");
    total += order; previous = order;
  }
  if (total != nodes) throw std::invalid_argument("Selection cache snapshot quadrature length is invalid.");
  for (std::size_t i = 0; i < nodes; ++i) {
    const double node = snapshot_real(data, 2 + rules + i);
    const double weight = snapshot_real(data, 2 + rules + nodes + i);
    if (!std::isfinite(node) || std::isnan(weight) || weight == std::numeric_limits<double>::infinity())
      throw std::invalid_argument("Selection cache snapshot quadrature values are invalid.");
  }
}
struct SnapshotHeader {
  std::size_t entries = 0;
  std::array<std::size_t, 2> controls_size{{0, 0}};
  std::array<const unsigned char *, 2> controls{{nullptr, nullptr}};
};
SnapshotHeader snapshot_header_read(SnapshotReader &reader)
{
  SnapshotHeader out;
  if (reader.size == 0) return out;
  if (sizeof(double) != 8 || !std::numeric_limits<double>::is_iec559 ||
      reader.size < snapshot_header || reader.integer() != snapshot_magic ||
      reader.integer() != 1 || reader.integer() != snapshot_abi())
    throw std::invalid_argument("Selection cache snapshot format or numerical ABI is incompatible.");
  out.entries = reader.count();
  for (int component = 0; component < 2; ++component) out.controls_size[component] = reader.count();
  for (int component = 0; component < 2; ++component)
    out.controls[component] = reader.reals(out.controls_size[component]);
  if (out.entries > (reader.size - reader.offset) / snapshot_record)
    throw std::invalid_argument("Selection cache snapshot entry count is invalid.");
  return out;
}
struct SnapshotEntry {
  EnvelopeMemoComponent component;
  const unsigned char *key;
  std::size_t size;
  EnvelopeMemoResult result;
  EnvelopeMemoView view(EnvelopeMemoSpan (&parts)[7]) const {
    const std::size_t k = snapshot_count(snapshot_real(key, 0), 1);
    const std::size_t bins = snapshot_count(snapshot_real(key, 1), 1);
    if (size < 8 || k > (size - 8) / k) throw std::invalid_argument("Selection cache snapshot key shape is invalid.");
    const std::size_t counts[] = {8, k, k * k, k, bins, bins, bins};
    std::size_t offset = 0;
    for (int i = 0; i < 7; ++i) {
      if (counts[i] > size - offset) throw std::invalid_argument("Selection cache snapshot key shape is invalid.");
      // Hash/equality use byte representations only. Insertion memcpy-copies
      // into actual owned doubles, never dereferencing raw-buffer doubles.
      parts[i] = {reinterpret_cast<const double *>(key + offset * sizeof(double)), counts[i]};
      offset += counts[i];
    }
    if (offset != size) throw std::invalid_argument("Selection cache snapshot key length is invalid.");
    return {parts, 7};
  }
};
SnapshotEntry snapshot_entry_read(SnapshotReader &reader)
{
  SnapshotEntry out;
  const std::uint64_t flags = reader.integer();
  if (flags & ~UINT64_C(0x301)) throw std::invalid_argument("Selection cache snapshot flags are invalid.");
  out.component = static_cast<EnvelopeMemoComponent>(flags & 1);
  out.result.success = (flags & 0x100) != 0; out.result.used_envelope = (flags & 0x200) != 0;
  out.size = reader.count();
  out.result.value = reader.real(); out.result.change = reader.real();
  out.result.width = reader.real(); out.result.tail = reader.real();
  out.key = reader.reals(out.size);
  if (out.size < 8) throw std::invalid_argument("Selection cache snapshot key is incomplete.");
  return out;
}
void snapshot_entry_validate(const SnapshotEntry &entry, const SnapshotHeader &header)
{
  const auto valid_diagnostic = [&](double value) {
    return !std::isnan(value) && value >= 0.0 &&
      (!entry.result.success || std::isfinite(value));
  };
  if (!header.controls_size[static_cast<std::size_t>(entry.component)] ||
      (!entry.result.success && entry.component == EnvelopeMemoComponent::exact) ||
      !std::isfinite(entry.result.value) || (!entry.result.success && entry.result.value != 0.0) ||
      !valid_diagnostic(entry.result.change) || !valid_diagnostic(entry.result.width) ||
      !valid_diagnostic(entry.result.tail))
    throw std::invalid_argument("Selection cache snapshot result is invalid.");
  EnvelopeMemoSpan parts[7]; entry.view(parts);
  const double sign = snapshot_real(entry.key, 2), tolerance = snapshot_real(entry.key, 7);
  if (!(sign == -1 || sign == 1) || !(tolerance > 0) || !std::isfinite(tolerance))
    throw std::invalid_argument("Selection cache snapshot key controls are invalid.");
  for (std::size_t i = 3; i <= 6; ++i) snapshot_count(snapshot_real(entry.key, i), 0);
  for (std::size_t i = 8; i < entry.size; ++i) {
    const double value = snapshot_real(entry.key, i);
    if (std::isnan(value)) throw std::invalid_argument("Selection cache snapshot key contains missing values.");
  }
}

struct EnvelopeMemoEntry {
  std::pmr::vector<double> key;
  std::size_t hash;
  EnvelopeMemoComponent component;
  EnvelopeMemoResult result;
  EnvelopeMemoEntry(std::pmr::memory_resource *resource, const EnvelopeMemoView &view,
      std::size_t hash_value, EnvelopeMemoComponent owner, const EnvelopeMemoResult &value)
    : key(resource), hash(hash_value), component(owner), result(value)
  {
    view.copy_to(&key);
  }
  EnvelopeMemoEntry(std::pmr::memory_resource *resource, const unsigned char *data,
      std::size_t count, std::size_t hash_value, EnvelopeMemoComponent owner, const EnvelopeMemoResult &value)
    : key(resource), hash(hash_value), component(owner), result(value)
  {
    key.resize(count);
    if (count) std::memcpy(key.data(), data, count * sizeof(double));
  }
};

using EnvelopeMemoEntries = std::pmr::list<EnvelopeMemoEntry>;
using EnvelopeMemoIndex =
  std::pmr::unordered_multimap<std::size_t, EnvelopeMemoEntries::iterator>;

// Conservative direct-allocation allowance for one LRU node, one index node,
// and bucket storage (including the old-plus-new arrays during rehash).
constexpr std::size_t envelope_entry_overhead = 512;

struct EnvelopeMemoStorage {
  std::array<std::pmr::vector<double>, 2> controls;
  std::array<bool, 2> controls_ready{{false, false}};
  EnvelopeMemoEntries entries;
  EnvelopeMemoIndex index;
  explicit EnvelopeMemoStorage(std::pmr::memory_resource *resource)
    : controls{{std::pmr::vector<double>(resource), std::pmr::vector<double>(resource)}},
      entries(resource), index(resource)
  {
    // Rehash's old and new buckets both count against the one shared cap.
    index.max_load_factor(1.0f);
  }
};

bool envelope_memo_add(std::size_t left, std::size_t right, std::size_t *result)
{
  if (right > std::numeric_limits<std::size_t>::max() - left) return false;
  *result = left + right;
  return true;
}

bool envelope_memo_entry_bytes(std::size_t key_doubles, std::size_t controls_doubles,
    std::size_t *result)
{
  std::size_t bytes = envelope_entry_overhead;
  if (key_doubles > (std::numeric_limits<std::size_t>::max() - bytes) / sizeof(double) ||
      controls_doubles > (std::numeric_limits<std::size_t>::max() - bytes) / sizeof(double))
    return false;
  if (!envelope_memo_add(bytes, key_doubles * sizeof(double), &bytes) ||
      !envelope_memo_add(bytes, controls_doubles * sizeof(double), result)) return false;
  return true;
}

struct EnvelopeMemoTicket {
  std::uint64_t generation = 0;
  std::size_t hash = 0;
  EnvelopeMemoComponent component = EnvelopeMemoComponent::exact;
  bool cacheable = false;
};

// Exact and coarse entries share one budget and global recency order, without
// sharing numerical identity or admission policy. Evaluation and borrowed-key
// hashing stay outside the mutex; no entry or owned pointer escapes its lock.
class EnvelopeMemo {
  std::mutex mutex;
  std::atomic<std::size_t> capacity;
  EnvelopeMemoResource resource;
  std::optional<EnvelopeMemoStorage> storage;
  std::array<std::uint64_t, 2> generation{{0, 0}};
  std::array<SelNormCacheStats, 2> counters{};

  static std::size_t slot(EnvelopeMemoComponent component)
  {
    return static_cast<std::size_t>(component);
  }
  bool make_room(std::size_t required)
  {
    if (required > resource.capacity) return false;
    while (storage && resource.capacity - resource.allocated < required &&
           !storage->entries.empty()) {
      erase(std::prev(storage->entries.end()), true);
    }
    if (storage && storage->entries.empty() &&
        resource.capacity - resource.allocated < required) {
      // Unordered-map buckets do not shrink when nodes are erased. Return that
      // retained capacity before declaring an otherwise empty cache full.
      EnvelopeMemoIndex empty(&resource);
      storage->index.swap(empty);
    }
    return resource.capacity - resource.allocated >= required;
  }
  SelNormCacheInfo info_unlocked() const
  {
    SelNormCacheInfo out;
    out.capacity_bytes = resource.capacity;
    out.allocated_bytes = resource.allocated;
    out.peak_bytes = resource.peak;
    out.exact = counters[0];
    out.coarse = counters[1];
    return out;
  }
  void release_storage()
  {
    storage.reset();
    for (std::size_t component = 0; component < 2; ++component) {
      ++generation[component];
      counters[component].entries = 0;
    }
  }
  void erase(EnvelopeMemoEntries::iterator entry, bool eviction)
  {
    const std::size_t component = slot(entry->component);
    const auto range = storage->index.equal_range(entry->hash);
    for (auto found = range.first; found != range.second; ++found) {
      if (found->second == entry) {
        storage->index.erase(found);
        storage->entries.erase(entry);
        --counters[component].entries;
        if (eviction) ++counters[component].evictions;
        return;
      }
    }
    throw std::logic_error("Selection cache index invariant failed.");
  }
  void clear_component(EnvelopeMemoComponent component)
  {
    const std::size_t owner = slot(component);
    ++generation[owner];
    if (!storage) return;
    for (auto entry = storage->entries.begin(); entry != storage->entries.end();) {
      auto current = entry++;
      if (current->component == component) erase(current, false);
    }
    std::pmr::vector<double> empty(&resource);
    storage->controls[owner].swap(empty);
    storage->controls_ready[owner] = false;
  }
  void initialize_controls(EnvelopeMemoComponent component, const EnvelopeMemoView &controls)
  {
    if (!storage) storage.emplace(&resource);
    const std::size_t owner = slot(component);
    controls.copy_to(&storage->controls[owner]);
    storage->controls_ready[owner] = true;
  }
  bool ensure_controls(EnvelopeMemoComponent component, const EnvelopeMemoView &controls)
  {
    const std::size_t owner = slot(component);
    std::size_t required = 0;
    if (storage && storage->controls_ready[owner]) {
      if (controls.equals(storage->controls[owner])) return true;
      ++counters[owner].resets;
      clear_component(component);
    }
    if (!envelope_memo_entry_bytes(0, controls.size(), &required) ||
        !make_room(required)) return false;
    try {
      initialize_controls(component, controls);
      return true;
    } catch (const std::bad_alloc &) {
      ++counters[owner].allocation_failures;
      return false;
    }
  }
  void insert(const EnvelopeMemoView &key, const EnvelopeMemoView &controls,
      const EnvelopeMemoTicket &ticket, const EnvelopeMemoResult &result)
  {
    const std::size_t owner = slot(ticket.component);
    std::size_t required = 0;
    if (!envelope_memo_entry_bytes(key.size(), 0, &required) || !make_room(required)) {
      ++counters[owner].allocation_failures;
      return;
    }
    bool restarted = false;
    while (storage) {
      auto created = storage->entries.end();
      try {
        storage->entries.emplace_front(&resource, key, ticket.hash, ticket.component, result);
        created = storage->entries.begin();
        storage->index.emplace(ticket.hash, created);
        ++counters[owner].entries;
        return;
      } catch (const std::bad_alloc &) {
        ++counters[owner].allocation_failures;
        if (created != storage->entries.end()) storage->entries.erase(created);
        if (!storage->entries.empty()) {
          // Evict several LRU records before retrying. This is only a fallback
          // for a rehash's temporary old-plus-new bucket requirement or an
          // allocator failure; ordinary admission evicted proactively above.
          const std::size_t batch = std::max<std::size_t>(16, storage->entries.size() / 16);
          for (std::size_t evicted = 0; evicted < batch && !storage->entries.empty(); ++evicted)
            erase(std::prev(storage->entries.end()), true);
          continue;
        }
        if (restarted) return;
        restarted = true;
        ++counters[owner].resets;
        release_storage();
        try {
          initialize_controls(ticket.component, controls);
        } catch (const std::bad_alloc &) {
          ++counters[owner].allocation_failures;
          release_storage();
          return;
        }
      }
    }
  }
public:
  // R resolves the total automatic budget and assigns this process its share.
  EnvelopeMemo() : capacity(0), resource(0) {}
  std::size_t limit() const { return capacity.load(std::memory_order_relaxed); }

  SelNormCacheInfo control(bool set_capacity, std::size_t bytes, unsigned int clear_mask)
  {
    if (clear_mask > 3) throw std::invalid_argument("Invalid selection cache clear mask.");
    const std::lock_guard<std::mutex> lock(mutex);
    const std::size_t requested = set_capacity ? bytes : resource.capacity;
    if (requested != resource.capacity || clear_mask == 3) {
      release_storage();
      resource.capacity = requested;
      resource.peak = 0;
      counters = {};
      capacity.store(requested, std::memory_order_relaxed);
    } else {
      for (std::size_t component = 0; component < 2; ++component) {
        if (clear_mask & (1U << component)) {
          clear_component(static_cast<EnvelopeMemoComponent>(component));
          counters[component] = {};
        }
      }
    }
    return info_unlocked();
  }
  bool lookup(EnvelopeMemoComponent component, const EnvelopeMemoView &key,
      const EnvelopeMemoView &controls, EnvelopeMemoTicket *ticket, EnvelopeMemoResult *result)
  {
    const std::size_t owner = slot(component);
    ticket->component = component;
    ticket->hash = key.hash() ^ (static_cast<std::size_t>(0x9e3779b97f4a7c15ULL) * owner);
    const std::lock_guard<std::mutex> lock(mutex);
    if (resource.capacity == 0) return false;
    std::size_t required = 0;
    bool usable = envelope_memo_entry_bytes(key.size(), controls.size(), &required) &&
      required <= resource.capacity;
    if (usable) usable = ensure_controls(component, controls);
    ticket->cacheable = usable;
    ticket->generation = generation[owner];
    if (usable) {
      const auto range = storage->index.equal_range(ticket->hash);
      for (auto found = range.first; found != range.second; ++found) {
        if (found->second->component != component || !key.equals(found->second->key)) continue;
        *result = found->second->result;
        storage->entries.splice(storage->entries.begin(), storage->entries, found->second);
        ++counters[owner].hits;
        return true;
      }
    }
    ++counters[owner].misses;
    return false;
  }
  void remember(const EnvelopeMemoView &key, const EnvelopeMemoView &controls,
      const EnvelopeMemoTicket &ticket, const EnvelopeMemoResult &result)
  {
    const std::lock_guard<std::mutex> lock(mutex);
    const std::size_t owner = slot(ticket.component);
    if (!ticket.cacheable || !storage || !storage->controls_ready[owner] ||
        ticket.generation != generation[owner]) return;
    const auto range = storage->index.equal_range(ticket.hash);
    for (auto found = range.first; found != range.second; ++found) {
      if (found->second->component == ticket.component && key.equals(found->second->key)) {
        storage->entries.splice(storage->entries.begin(), storage->entries, found->second);
        return;
      }
    }
    insert(key, controls, ticket, result);
  }
  std::size_t snapshot_size_unlocked() const
  {
    if (!storage || storage->entries.empty()) return 0;
    std::size_t bytes = snapshot_header + 8;
    const auto add = [&](std::size_t count, std::size_t unit) {
      if (count > (std::numeric_limits<std::size_t>::max() - bytes) / unit)
        throw std::length_error("Selection cache snapshot is too large.");
      bytes += count * unit;
    };
    for (int component = 0; component < 2; ++component)
      add(storage->controls_ready[component] ? storage->controls[component].size() : 0, sizeof(double));
    for (const auto &entry : storage->entries) { add(1, snapshot_record); add(entry.key.size(), sizeof(double)); }
    if (bytes > resource.capacity) throw std::length_error("Selection cache snapshot exceeds the configured byte budget.");
    return bytes;
  }
  std::size_t snapshot_size()
  {
    const std::lock_guard<std::mutex> lock(mutex);
    return snapshot_size_unlocked();
  }
  bool snapshot_write(unsigned char *data, std::size_t size)
  {
    const std::lock_guard<std::mutex> lock(mutex);
    if (snapshot_size_unlocked() != size) return false;
    if (!size) return true;
    SnapshotWriter writer{data};
    writer.integer(snapshot_magic); writer.integer(1); writer.integer(snapshot_abi());
    writer.integer(storage->entries.size());
    for (int component = 0; component < 2; ++component)
      writer.integer(storage->controls_ready[component] ? storage->controls[component].size() : 0);
    for (int component = 0; component < 2; ++component)
      if (storage->controls_ready[component]) writer.bytes(storage->controls[component].data(), storage->controls[component].size() * sizeof(double));
    for (const auto &entry : storage->entries) {
      writer.integer(static_cast<std::uint64_t>(entry.component) |
        (std::uint64_t(entry.result.success) << 8) | (std::uint64_t(entry.result.used_envelope) << 9));
      writer.integer(entry.key.size());
      writer.real(entry.result.value); writer.real(entry.result.change);
      writer.real(entry.result.width); writer.real(entry.result.tail);
      writer.bytes(entry.key.data(), entry.key.size() * sizeof(double));
    }
    if (writer.offset != size - 8) throw std::logic_error("Selection cache snapshot size invariant failed.");
    writer.integer(snapshot_checksum(data, size - 8));
    return true;
  }
  SelNormCacheRestoreInfo restore(const SelNormCacheBlob *blobs, std::size_t count)
  {
    SelNormCacheRestoreInfo out; out.snapshots = count;
    // Validate every shard, including all later/cold records, before clearing.
    for (std::size_t shard = 0; shard < count; ++shard) {
      if (!blobs[shard].size) continue;
      if (blobs[shard].size < snapshot_header + 8) throw std::invalid_argument("Selection cache snapshot is truncated.");
      SnapshotReader checksum{blobs[shard].data + blobs[shard].size - 8, 8};
      if (checksum.integer() != snapshot_checksum(blobs[shard].data, blobs[shard].size - 8))
        throw std::invalid_argument("Selection cache snapshot checksum does not match.");
      SnapshotReader reader{blobs[shard].data, blobs[shard].size - 8};
      const auto header = snapshot_header_read(reader);
      for (int component = 0; component < 2; ++component) snapshot_controls(header.controls[component], header.controls_size[component]);
      for (std::size_t i = 0; i < header.entries; ++i) {
        const auto entry = snapshot_entry_read(reader); snapshot_entry_validate(entry, header);
        if (entry.component == EnvelopeMemoComponent::exact) ++out.available_exact; else ++out.available_coarse;
      }
      if (reader.offset != reader.size) throw std::invalid_argument("Selection cache snapshot has trailing records.");
    }
    const std::lock_guard<std::mutex> lock(mutex);
    release_storage(); resource.peak = 0; counters = {};
    bool full = resource.capacity == 0;
    for (std::size_t shard = 0; shard < count && !full; ++shard) {
      if (!blobs[shard].size) continue;
      SnapshotReader reader{blobs[shard].data, blobs[shard].size - 8};
      const auto header = snapshot_header_read(reader);
      for (std::size_t i = 0; i < header.entries && !full; ++i) {
        const auto entry = snapshot_entry_read(reader);
        const std::size_t owner = slot(entry.component);
        if (storage && storage->controls_ready[owner] &&
            (storage->controls[owner].size() != header.controls_size[owner] ||
             std::memcmp(storage->controls[owner].data(), header.controls[owner], header.controls_size[owner] * sizeof(double)) != 0)) continue;
        EnvelopeMemoSpan parts[7]; const auto key = entry.view(parts);
        const std::size_t hash = key.hash() ^ (static_cast<std::size_t>(0x9e3779b97f4a7c15ULL) * owner);
        if (storage) {
          bool duplicate = false;
          const auto range = storage->index.equal_range(hash);
          for (auto found = range.first; found != range.second; ++found)
            if (found->second->component == entry.component && key.equals(found->second->key)) { duplicate = true; break; }
          if (duplicate) continue;
        }
        std::size_t required = 0;
        if (!envelope_memo_entry_bytes(entry.size,
              storage && storage->controls_ready[owner] ? 0 : header.controls_size[owner],
              &required) ||
            !make_room(required)) {
          if (resource.capacity != 0) ++counters[owner].allocation_failures;
          full = true;
          break;
        }
        bool created = false;
        try {
          if (!storage) storage.emplace(&resource);
          if (!storage->controls_ready[owner]) {
            storage->controls[owner].resize(header.controls_size[owner]);
            std::memcpy(storage->controls[owner].data(), header.controls[owner], header.controls_size[owner] * sizeof(double));
            storage->controls_ready[owner] = true;
          }
          // Input order is MRU-first. Proactive admission evicts only from the
          // cold end, so a smaller budget retains the warm prefix and order.
          storage->entries.emplace_back(&resource, entry.key, entry.size, hash, entry.component, entry.result);
          created = true;
          storage->index.emplace(hash, std::prev(storage->entries.end()));
          ++counters[owner].entries;
          if (owner == 0) ++out.restored_exact; else ++out.restored_coarse;
        } catch (const std::bad_alloc &) {
          if (created) storage->entries.pop_back();
          ++counters[owner].allocation_failures; full = true;
        }
      }
    }
    if (storage && storage->entries.empty()) release_storage();
    out.skipped = out.available_exact + out.available_coarse - out.restored_exact - out.restored_coarse;
    return out;
  }
};

EnvelopeMemo &envelope_memo()
{
  static EnvelopeMemo memo;
  return memo;
}

double bounded_cdf_product_log_error(const SelNormKernelData &selection,
    const double *omega, int dimension, double tolerance)
{
  if (selection.n_bins != 2 || !selection.trusted_step_partition ||
      !selection.telescope_probabilities || !std::isfinite(selection.z_lower[0]) ||
      dimension < 1) return 0.0;
  const double minimum = std::min(omega[0], omega[1]);
  const double maximum = std::max(omega[0], omega[1]);
  if (!(minimum >= std::numeric_limits<double>::min()) ||
      !std::isfinite(maximum) || !(maximum > minimum)) return 0.0;
  const long double u = std::numeric_limits<double>::epsilon() / 2.0L;
  const long double gamma_blend = 5 * u / (1 - 5 * u);
  // Positive blending admits a relative operation-count bound independent
  // of cancellation. Five roundings also cover a subnormal difference*CDF
  // product, since the positive base is at least the normal minimum.
  const long double u_long = std::numeric_limits<long double>::epsilon() / 2.0L;
  const long double guard = 1 - 8 * u_long;
  const long double delta = ((static_cast<long double>(maximum) - minimum) /
    minimum * cpp_selnorm_bounded_cdf_absolute_error()) / guard;
  if (!(delta >= 0 && delta < 1)) return 0.0;
  const long double bound = -static_cast<long double>(dimension) *
    (std::log1p(-delta) + std::log1p(-gamma_blend)) / guard;
  const double upper = std::nextafter(static_cast<double>(bound),
    std::numeric_limits<double>::infinity());
  if (!std::isfinite(upper) || !(upper > 0)) return 0.0;
  // Even perfect quadrature cannot beat the two-envelope CDF-only score.
  // If it already exhausts the requested budget, use the ordinary primitive
  // immediately rather than refining an impossible approximate calculation.
  const double minimum_score = 3 * std::exp(upper) * std::expm1(2 * upper) + std::expm1(upper);
  if (!std::isfinite(minimum_score) || !(minimum_score < tolerance)) return 0.0;
  return upper;
}

bool dense_envelope_log_integral(
    const std::vector<double> &covariance, const double *mean,
    const double *selection_se, const double *omega, int dimension,
    const SelNormKernelData &selection, SelNormDenseIntegration *integration,
    double *log_normalizer)
{
  if (integration == nullptr || integration->rule_count < 3 ||
      !(integration->tolerance > 0.0)) return false;
  integration->cdf_log_error = 0.0;
  integration->cdf_relative_error = 0.0;
  const double cdf_log_error = integration->bounded_cdf ?
    bounded_cdf_product_log_error(selection, omega, dimension, integration->tolerance) : 0.0;
  SelNormKernelData evaluation_selection = selection;
  evaluation_selection.bounded_cdf = cdf_log_error > 0.0;
  const double cdf_scale = cdf_log_error > 0 ? std::exp(cdf_log_error) : 1.0;
  const double cdf_scale_twice = cdf_log_error > 0 ? std::exp(2 * cdf_log_error) : 1.0;
  const double cdf_relative_error = cdf_log_error > 0 ? std::expm1(cdf_log_error) : 0.0;
  const double cdf_change_error = cdf_log_error > 0 ? std::expm1(2 * cdf_log_error) : 0.0;
  const auto change_bound = [&](double raw) {
    // Exact-primitive sensitivity first, then convert its denominator to
    // the returned approximate integral. L=0 retains the original arithmetic.
    return cdf_log_error > 0 ? cdf_scale * (cdf_scale_twice * raw + cdf_change_error) : raw;
  };
  // Covariance ordering brackets normalizers only for monotone row weights.
  // Nonmonotone product weights instead require an absolute perturbation
  // bound; agreement of two auxiliary integrals cannot replace that bound.
  bool increasing = true;
  bool decreasing = true;
  for (int bin = 1; bin < selection.n_bins; ++bin) {
    increasing = increasing && omega[bin - 1] >= omega[bin];
    decreasing = decreasing && omega[bin - 1] <= omega[bin];
  }
  const bool monotone = increasing || decreasing;
  SelNormCovarianceEnvelope envelope;
  if (!covariance_envelope(covariance, dimension, selection_se, &envelope)) return false;
  // The absolute bound needs only the envelope, so evaluate it before any
  // geometry work: a nonmonotone case whose weights already rule out the
  // relative bound can then stop without preparing quadrature geometry.
  const double gap_bound = covariance_price_gap_bound(
    covariance, selection_se, dimension, envelope, omega, selection.n_bins, monotone
  );
  const double log_gap = gap_bound > 0.0 ? std::log(gap_bound) :
    -std::numeric_limits<double>::infinity();
  bool positive_weights = true;
  double maximum_weight = 0.0;
  for (int bin = 0; bin < selection.n_bins; ++bin) {
    positive_weights = positive_weights && omega[bin] > 0.0 && std::isfinite(omega[bin]);
    maximum_weight = std::max(maximum_weight, omega[bin]);
  }
  if (!monotone && !std::isfinite(gap_bound) && !positive_weights) return false;
  std::vector<EnvelopeGroupGeometry> lower_geometry, upper_geometry;
  if (!prepare_envelope_geometry(covariance, mean, selection_se, dimension,
                                 envelope, false, &lower_geometry)) return false;
  // The upper geometry is only consumed by a monotone direct comparison. The
  // covariance-bound route below accepts many states without ever evaluating
  // it, so its matrices stay lazy. The original eager rejection semantics are
  // preserved by the matching feasibility check below: it repeats the exact
  // residual arithmetic of prepare_envelope_geometry's upper pass.
  if (monotone && !envelope_geometry_feasible(covariance, selection_se, dimension,
                                              envelope, true)) return false;
  bool upper_prepared = false;
  const auto prepare_upper = [&]() {
    if (upper_prepared) return true;
    upper_prepared = prepare_envelope_geometry(covariance, mean, selection_se,
                                               dimension, envelope, true, &upper_geometry);
    return upper_prepared;
  };

  const double relative_gap_bound = monotone ?
    std::numeric_limits<double>::infinity() : covariance_relative_gap_bound(
      covariance, dimension, lower_geometry, omega, selection.n_bins);
  if (!monotone && !std::isfinite(gap_bound) && !std::isfinite(relative_gap_bound)) return false;
  // The product normalizer is at most max(omega)^dimension, so the absolute
  // bound can never yield a relative gap below this rule-independent value.
  // If neither that floor nor the relative bound can reach the tolerance, no
  // amount of refinement will accept, so stop before evaluating any rule.
  const bool use_relative_bound = !monotone && std::isfinite(relative_gap_bound);
  if (!monotone) {
    const double gap_floor = std::exp(log_gap -
      static_cast<double>(dimension) * std::log(maximum_weight));
    if (!(gap_floor <= integration->tolerance) &&
        !(use_relative_bound && relative_gap_bound <= integration->tolerance)) return false;
  }
  EnvelopeRuleWorkspace workspace;
  double previous_lower = std::numeric_limits<double>::quiet_NaN();
  double previous_upper = std::numeric_limits<double>::quiet_NaN();
  int previous_offset = 0;
  int previous_order = 0;
  std::vector<long double> weights;
  int offset = 0;
  for (int rule = 0; rule < integration->rule_count; ++rule) {
    const int order = static_cast<int>(integration->orders[rule]);
    // These weights are identical for both envelopes and every conditional
    // integral. Prepare only the rules actually visited by refinement.
    weights.resize(offset + order);
    for (int point = 0; point < order; ++point) {
      weights[offset + point] = std::exp(
        static_cast<long double>(integration->log_weights[offset + point])
      );
    }
    const double lower = envelope_rule_log_integral(lower_geometry, omega, evaluation_selection,
      integration->nodes, integration->log_weights, weights.data(), offset, order, &workspace);
    if (!std::isfinite(lower)) return false;
    const double raw_lower_change = cluster_relative_change(previous_lower, lower);
    const double lower_change = change_bound(raw_lower_change);
    const double log_tail = quadrature_log_gaussian_tail_bound(
      integration->nodes, offset, order,
      envelope.groups == 1 ? 1 : envelope.groups + 1,
      omega, selection.n_bins, dimension, SELVECTOR_PRODUCT);
    double relative_gap = cdf_log_error > 0 ?
      cdf_scale_twice * std::exp(log_gap - lower) : std::exp(log_gap - lower);
    const double lower_tail = cdf_log_error > 0 ?
      cdf_scale_twice * std::exp(log_tail - lower) : std::exp(log_tail - lower);
    // The relative covariance bound uses the exact auxiliary mass A(C0). The
    // returned denominator is the quadrature value A_hat, and
    // A(C0)/A_hat <= 1/(1 - quadrature_error), so divide by the complement
    // rather than multiplying by the first-order factor 1 + quadrature_error,
    // which understates the converted bound by O(quadrature_error^2). Round
    // every step outward, as elsewhere in this file. Nonmonotone weights never
    // use bounded CDFs: they need at least three bins, while the bounded-CDF
    // path requires exactly two.
    const double quadrature_error = add_up(mul_up(2.0, lower_change), lower_tail);
    if (use_relative_bound && std::isfinite(quadrature_error) && quadrature_error < 1.0) {
      const double complement = std::nextafter(1.0 - quadrature_error,
        -std::numeric_limits<double>::infinity());
      if (complement > 0.0) {
        relative_gap = std::min(relative_gap, std::nextafter(
          relative_gap_bound / complement, std::numeric_limits<double>::infinity()));
      }
    }
    const bool use_gap_bound = std::isfinite(relative_gap) && relative_gap >= 0.0 &&
      relative_gap <= integration->tolerance;
    if (use_gap_bound &&
        relative_gap + 2.0 * lower_change + lower_tail + cdf_relative_error <= integration->tolerance) {
      // The derivative bound measures uncertainty; adding a fraction of that
      // deliberately loose bound would introduce an artificial correction.
      *log_normalizer = lower;
      integration->used_envelope = true;
      integration->quadrature_change = lower_change;
      integration->covariance_width = relative_gap;
      integration->tail_bound = lower_tail;
      integration->cdf_log_error = cdf_log_error;
      integration->cdf_relative_error = cdf_relative_error;
      return true;
    }

    double upper = std::numeric_limits<double>::quiet_NaN();
    if (monotone && (!use_gap_bound || rule + 1 == integration->rule_count) &&
        prepare_upper()) {
      // A loose bound restores the two actual envelope quadratures. A deferred
      // previous upper value is evaluated at its original rule, not replaced
      // by a lower-covariance value. The final rule always tries this route.
      upper = envelope_rule_log_integral(upper_geometry, omega, evaluation_selection,
        integration->nodes, integration->log_weights, weights.data(), offset, order, &workspace);
      if (rule > 0 && std::isnan(previous_upper)) {
        previous_upper = envelope_rule_log_integral(upper_geometry, omega, evaluation_selection,
          integration->nodes, integration->log_weights, weights.data(),
          previous_offset, previous_order, &workspace);
      }
      if (!std::isfinite(upper)) return false;
      const double change = change_bound(std::max(raw_lower_change,
        cluster_relative_change(previous_upper, upper)));
      // Finite quadrature endpoints can cross. Count their absolute discrepancy
      // and estimated interior errors; do not manufacture a rigorous interval.
      const double width = cdf_log_error > 0 ?
        cdf_scale * std::expm1(std::fabs(upper - lower) + 2 * cdf_log_error) :
        std::expm1(std::fabs(upper - lower));
      const double tail = cdf_log_error > 0 ?
        cdf_scale_twice * std::exp(log_tail - std::min(lower, upper)) :
        std::exp(log_tail - std::min(lower, upper));
      if (width + 2.0 * change + tail + cdf_relative_error <= integration->tolerance) {
        *log_normalizer = log_add_exp(lower, upper) - std::log(2.0);
        integration->used_envelope = true;
        integration->quadrature_change = change;
        integration->covariance_width = width;
        integration->tail_bound = tail;
        integration->cdf_log_error = cdf_log_error;
        integration->cdf_relative_error = cdf_relative_error;
        return true;
      }
      if (2.0 * change + tail + cdf_relative_error <= integration->tolerance && width > integration->tolerance) {
        return false;
      }
    }
    previous_lower = lower;
    previous_upper = upper;
    previous_offset = offset;
    previous_order = order;
    offset += order;
  }
  return false;
}

// Component identity separates exact and surrogate values in the shared pool.
// Only coarse anchors may store the prescribed failed-envelope fallback.
bool memoized_envelope_log_integral(
    const std::vector<double> &covariance, const double *mean,
    const double *selection_se, const double *omega, int dimension,
    const SelNormKernelData &selection, SelNormDenseIntegration *integration,
    double *log_normalizer, EnvelopeMemoComponent component)
{
  EnvelopeMemo &memo = envelope_memo();
  const bool cache_failures = component == EnvelopeMemoComponent::coarse;
  const std::size_t capacity = memo.limit();
  if (capacity == 0 || integration == nullptr || integration->bounded_cdf || integration->rule_count < 3 ||
      !(integration->tolerance > 0.0)) {
    const bool success = dense_envelope_log_integral(covariance, mean, selection_se,
      omega, dimension, selection, integration, log_normalizer);
    if (!success && cache_failures) *log_normalizer = 0.0;
    return success;
  }
  std::size_t nodes = 0;
  for (int rule = 0; rule < integration->rule_count; ++rule) {
    const double order = integration->orders[rule];
    if (!(order >= 1) || !std::isfinite(order) || order != std::floor(order) ||
        order > capacity / sizeof(double) - nodes) {
      const bool success = dense_envelope_log_integral(covariance, mean, selection_se,
        omega, dimension, selection, integration, log_normalizer);
      if (!success && cache_failures) *log_normalizer = 0.0;
      return success;
    }
    nodes += static_cast<std::size_t>(order);
  }
  const double state_meta[] = {static_cast<double>(dimension), static_cast<double>(selection.n_bins),
    static_cast<double>(selection.effect_sign), static_cast<double>(selection.n_segments),
    static_cast<double>(selection.q), static_cast<double>(selection.trusted_step_partition),
    static_cast<double>(selection.telescope_probabilities), integration->tolerance};
  const double control_meta[] = {static_cast<double>(integration->rule_count),
    static_cast<double>(nodes)};
  const EnvelopeMemoSpan state_parts[] = {{state_meta, 8},
    {mean, static_cast<std::size_t>(dimension)}, {covariance.data(), covariance.size()},
    {selection_se, static_cast<std::size_t>(dimension)},
    {omega, static_cast<std::size_t>(selection.n_bins)},
    {selection.z_lower, static_cast<std::size_t>(selection.n_bins)},
    {selection.z_upper, static_cast<std::size_t>(selection.n_bins)}};
  const EnvelopeMemoSpan control_parts[] = {{control_meta, 2},
    {integration->orders, static_cast<std::size_t>(integration->rule_count)},
    {integration->nodes, nodes}, {integration->log_weights, nodes}};
  const EnvelopeMemoView key{state_parts, 7}, controls{control_parts, 4};
  EnvelopeMemoTicket ticket;
  EnvelopeMemoResult result;
  if (memo.lookup(component, key, controls, &ticket, &result)) {
    *log_normalizer = result.value;
    integration->used_envelope = result.used_envelope;
    integration->quadrature_change = result.change;
    integration->covariance_width = result.width;
    integration->tail_bound = result.tail;
    return result.success;
  }
  result.success = dense_envelope_log_integral(covariance, mean, selection_se,
    omega, dimension, selection, integration, log_normalizer);
  if (!result.success && cache_failures) *log_normalizer = 0.0;
  if (result.success || cache_failures) {
    result.value = *log_normalizer;
    result.used_envelope = integration->used_envelope;
    result.change = integration->quadrature_change;
    result.width = integration->covariance_width;
    result.tail = integration->tail_bound;
    memo.remember(key, controls, ticket, result);
  }
  return result.success;
}

bool exact_envelope_log_integral(
    const std::vector<double> &covariance, const double *mean,
    const double *selection_se, const double *omega, int dimension,
    const SelNormKernelData &selection, SelNormDenseIntegration *integration,
    double *log_normalizer)
{
  const bool retry_ordinary = integration != nullptr && integration->bounded_cdf &&
    bounded_cdf_product_log_error(selection, omega, dimension, integration->tolerance) > 0.0;
  if (memoized_envelope_log_integral(covariance, mean, selection_se, omega,
      dimension, selection, integration, log_normalizer, EnvelopeMemoComponent::exact)) return true;
  if (!retry_ordinary) return false;
  // Exhausting the bounded-CDF budget does not bypass the original envelope
  // algorithm. Retry it with its exact primitive, same rules and criterion,
  // before the caller's original QMC fallback. Neither attempt touches caches.
  SelNormDenseIntegration ordinary = *integration;
  ordinary.bounded_cdf = false;
  ordinary.used_envelope = false;
  ordinary.quadrature_change = ordinary.covariance_width = ordinary.tail_bound = 0.0;
  ordinary.cdf_log_error = ordinary.cdf_relative_error = 0.0;
  const bool success = dense_envelope_log_integral(covariance, mean, selection_se,
    omega, dimension, selection, &ordinary, log_normalizer);
  *integration = ordinary;
  integration->bounded_cdf = true;
  return success;
}

bool coarse_envelope_log_integral(
    const std::vector<double> &covariance, const double *mean,
    const double *selection_se, const double *omega, int dimension,
    const SelNormKernelData &selection, SelNormDenseIntegration *integration,
    double *log_normalizer)
{
  return memoized_envelope_log_integral(covariance, mean, selection_se, omega,
    dimension, selection, integration, log_normalizer, EnvelopeMemoComponent::coarse);
}

// Called only by the explicit corrected-sampler interface. Arithmetic-range or
// unsupported-anchor failures prescribe A=1. Allocation failures propagate;
// cache allocation cannot choose a different surrogate value.
double coarse_anchor_log_normalizer(
    const double *mean, const double *covariance_lower, int dimension,
    const double *selection_se, const double *omega, int n_bins,
    const double *z_lower, const double *z_upper, int effect_sign,
    bool telescope_probabilities, int kernel_mode, int vector_rule,
    const SelNormDenseIntegration &quadrature, const SelNormCoarseSettings &settings)
{
  if (!(settings.mean_step >= 0.0) || !std::isfinite(settings.mean_step) ||
      !(settings.diagonal_step >= 0.0) || !std::isfinite(settings.diagonal_step) ||
      !(settings.log_weight_step >= 0.0) || !std::isfinite(settings.log_weight_step)) {
    throw std::invalid_argument("Coarse selection-grid steps must be finite and nonnegative.");
  }
  if (dimension < 1 || n_bins < 1 || kernel_mode != SELKERNEL_STEP ||
      vector_rule != SELVECTOR_PRODUCT || settings.max_rules < 3 ||
      quadrature.rule_count < 3) return 0.0;
  double scale = 0.0;
  for (int row = 0; row < dimension; ++row) {
    if (!(selection_se[row] > 0.0) || !std::isfinite(selection_se[row]) ||
        !std::isfinite(mean[row])) return 0.0;
    scale = std::max(scale, selection_se[row]);
  }
  const long double mean_step = static_cast<long double>(settings.mean_step) * scale;
  const long double diagonal_step = static_cast<long double>(settings.diagonal_step) * scale * scale;
  std::vector<double> anchor_mean(mean, mean + dimension);
  std::vector<double> anchor_weight(omega, omega + n_bins);
  std::vector<double> anchor_covariance(static_cast<std::size_t>(dimension) * dimension);
  if (settings.mean_step > 0.0) {
    if (!(mean_step > 0.0L) || !std::isfinite(mean_step)) return 0.0;
    for (int row = 0; row < dimension; ++row) {
      const long double coordinate = static_cast<long double>(mean[row]) / mean_step;
      const double value = static_cast<double>(std::round(coordinate) * mean_step);
      if (!std::isfinite(value)) return 0.0;
      anchor_mean[row] = value;
    }
  }
  if (settings.diagonal_step > 0.0 &&
      (!(diagonal_step > 0.0L) || !std::isfinite(diagonal_step))) return 0.0;
  int position = 0;
  for (int column = 0; column < dimension; ++column) {
    for (int row = column; row < dimension; ++row) {
      const double original = covariance_lower[position++];
      if (!std::isfinite(original)) return 0.0;
      double value = original;
      if (row == column) {
        if (!(original > 0.0)) return 0.0;
        if (settings.diagonal_step > 0.0) {
          const long double origin = static_cast<long double>(selection_se[row]) * selection_se[row];
          const long double coordinate = (static_cast<long double>(original) - origin) / diagonal_step;
          value = static_cast<double>(origin + std::ceil(coordinate) * diagonal_step);
          // The returned matrix equals the actual covariance plus a
          // nonnegative diagonal. Never nudge a failed rounded candidate.
          if (!std::isfinite(value) || value < original) return 0.0;
        }
      }
      // Every off-diagonal assignment copies the supplied double unchanged.
      anchor_covariance[row + dimension * column] = value;
      anchor_covariance[column + dimension * row] = value;
    }
  }
  for (int bin = 0; bin < n_bins; ++bin) {
    if (!(omega[bin] >= 0.0) || !std::isfinite(omega[bin])) return 0.0;
    if (omega[bin] == 0.0 || settings.log_weight_step == 0.0) continue;
    const long double coordinate = std::log(static_cast<long double>(omega[bin])) /
      static_cast<long double>(settings.log_weight_step);
    const double value = static_cast<double>(std::exp(std::round(coordinate) *
      static_cast<long double>(settings.log_weight_step)));
    if (!(value > 0.0) || !std::isfinite(value)) return 0.0;
    anchor_weight[bin] = value;
  }
  SelNormDenseIntegration integration = quadrature;
  integration.rule_count = std::min<unsigned int>(settings.max_rules, quadrature.rule_count);
  integration.used_envelope = false;
  integration.quadrature_change = 0.0;
  integration.covariance_width = 0.0;
  integration.tail_bound = 0.0;
  SelNormKernelData selection = {};
  selection.n_bins = n_bins;
  selection.effect_sign = effect_sign;
  selection.z_lower = z_lower;
  selection.z_upper = z_upper;
  selection.trusted_step_partition = true;
  selection.telescope_probabilities = telescope_probabilities;
  double log_normalizer = 0.0;
  // The independent process-owned coarse memo stores both accepted anchors
  // and their prescribed deterministic unit-normalizer fallback. Its boolean
  // distinguishes those outcomes; either output is the surrogate definition.
  coarse_envelope_log_integral(anchor_covariance, anchor_mean.data(), selection_se,
    anchor_weight.data(), dimension, selection, &integration, &log_normalizer);
  return log_normalizer;
}

struct FactorNormalizerContext {
  int vector_rule = SELVECTOR_PRODUCT;
  const double *mean;
  const double *residual_sd;
  const double *loading;
  int dimension;
  int rank;
  const double *selection_se;
  const double *omega;
  int n_bins;
  int kernel_mode;
  SelNormKernelData selection;
};

double factor_log_integrand(const FactorNormalizerContext &context,
                            const double *latent, double *gradient)
{
  double value = 0.0;
  for (int factor = 0; factor < context.rank; ++factor) {
    value -= 0.5 * latent[factor] * latent[factor];
    if (gradient != nullptr) gradient[factor] = -latent[factor];
  }

  if (context.vector_rule != SELVECTOR_PRODUCT) {
    std::vector<double> means(context.dimension);
    std::vector<double> selection_gradient(context.rank);
    for (int i = 0; i < context.dimension; ++i) {
      means[i] = context.mean[i];
      for (int factor = 0; factor < context.rank; ++factor) {
        means[i] += context.loading[i + context.dimension * factor] * latent[factor];
      }
    }
    value += independent_best_log_mass(means.data(), context.residual_sd,
      context.dimension, context.selection_se, context.omega,
      context.selection, context.vector_rule, nullptr, nullptr,
      context.loading, context.rank,
      gradient == nullptr ? nullptr : selection_gradient.data());
    if (gradient != nullptr) {
      for (int factor = 0; factor < context.rank; ++factor) {
        gradient[factor] += selection_gradient[factor];
      }
    }
    return value;
  }

  for (int i = 0; i < context.dimension; ++i) {
    double conditional_mean = context.mean[i];
    for (int factor = 0; factor < context.rank; ++factor) {
      conditional_mean +=
        context.loading[i + context.dimension * factor] * latent[factor];
    }
    const double local = cpp_selnorm_step_log_norm(
      conditional_mean, context.residual_sd[i], context.selection_se[i],
      context.omega, context.selection, 1, false
    );
    if (!std::isfinite(local)) {
      return -std::numeric_limits<double>::infinity();
    }
    value += local;

    if (gradient != nullptr) {
      double derivative_numerator = 0.0;
      const double signed_mean = context.selection.effect_sign *
        conditional_mean;
      for (int bin = 0; bin < context.n_bins; ++bin) {
        const double lower = context.selection.z_lower[bin] *
          context.selection_se[i];
        const double upper = context.selection.z_upper[bin] *
          context.selection_se[i];
        const double lower_score =
          (lower - signed_mean) / context.residual_sd[i];
        const double upper_score =
          (upper - signed_mean) / context.residual_sd[i];
        const double lower_density = std::isfinite(lower_score) ?
          dnorm(lower_score, 0.0, 1.0, false) : 0.0;
        const double upper_density = std::isfinite(upper_score) ?
          dnorm(upper_score, 0.0, 1.0, false) : 0.0;
        derivative_numerator += context.omega[bin] *
          (lower_density - upper_density);
      }
      const double local_probability = std::exp(local);
      if (!(local_probability > 0.0) ||
          !std::isfinite(derivative_numerator)) {
        return -std::numeric_limits<double>::infinity();
      }
      const double derivative = context.selection.effect_sign *
        derivative_numerator /
        (context.residual_sd[i] * local_probability);
      if (!std::isfinite(derivative)) {
        return -std::numeric_limits<double>::infinity();
      }
      for (int factor = 0; factor < context.rank; ++factor) {
        gradient[factor] +=
          context.loading[i + context.dimension * factor] * derivative;
      }
    }
  }
  return value;
}

struct NormalProjectionGrid {
  const double *z;
  int size;
  std::vector<long double> gaps;
  std::vector<int> gap_index;
  std::vector<std::vector<int>> positions;
  std::vector<long double> delta, forward, backward, right, left;
  struct RepeatedDistance {
    int gap;
    long double value;
    std::vector<std::pair<int, int>> occurrences;
  };
  std::vector<RepeatedDistance> distances;

  NormalProjectionGrid(const double *grid, int length) : z(grid), size(length)
  {
    if (size < 2 || !std::isfinite(z[0])) return;
    for (int i = 1; i < size; ++i) {
      const long double gap = static_cast<long double>(z[i]) - z[i - 1];
      if (!(gap > 0.0L) || !std::isfinite(z[i])) {
        positions.clear();
        return;
      }
      const auto found = std::find(gaps.begin(), gaps.end(), gap);
      gap_index.push_back(static_cast<int>(found - gaps.begin()));
      if (found == gaps.end()) {
        gaps.push_back(gap);
        positions.emplace_back();
      }
      // Ratio initialization only pays for itself when spacings repeat.
      // Compare actual spacings exactly; never replace a rounded grid by seq.
      if (gaps.size() * 4 > static_cast<std::size_t>(size)) {
        positions.clear();
        return;
      }
      positions[gap_index.back()].push_back(i - 1);
    }
    delta.resize(gaps.size());
    right.resize(gaps.size());
    left.resize(gaps.size());
    forward.resize(size - 1);
    backward.resize(size - 1);
    for (std::size_t g = 0; g < positions.size(); ++g) {
      const std::size_t first_distance = distances.size();
      for (std::size_t j = 1; j < positions[g].size(); ++j) {
        const int previous = positions[g][j - 1];
        const int next = positions[g][j];
        const long double distance = static_cast<long double>(z[next]) - z[previous];
        auto found = std::find_if(distances.begin() + first_distance, distances.end(),
          [distance](RepeatedDistance const &entry) { return entry.value == distance; });
        if (found == distances.end()) {
          distances.push_back({static_cast<int>(g), distance, {}});
          found = distances.end() - 1;
        }
        found->occurrences.emplace_back(previous, next);
      }
    }
  }

  void prepare(double sei, double sd)
  {
    for (std::size_t g = 0; g < positions.size(); ++g) {
      delta[g] = gaps[g] * sei / sd;
    }
    for (RepeatedDistance const &entry : distances) {
      // Identical actual grid distances share the exact original arithmetic;
      // neither the grid nor its gaps are rounded to a regular sequence.
      const long double distance = entry.value * sei / sd;
      const long double multiplier = std::exp(-delta[entry.gap] * distance);
      for (auto const &occurrence : entry.occurrences) {
        forward[occurrence.first] = multiplier;
        backward[occurrence.second] = multiplier;
      }
    }
  }

  bool add(double mean, double sei, double sd, double log_scale,
           long double weight, const std::vector<double> &z_weight,
           std::vector<long double> &density)
  {
    if (positions.empty()) return false;
    if (weight == 0.0L) return true;
    if (!(weight > 0.0L) || !std::isfinite(weight)) return false;
    int anchor = static_cast<int>(std::lower_bound(z, z + size, mean / sei) - z);
    anchor = std::min(anchor, size - 1);
    const auto score = [this, sei, mean, sd](int i) {
      return selnorm_fma(z[i], sei, -mean) / sd;
    };
    if (anchor > 0 && std::abs(score(anchor - 1)) < std::abs(score(anchor))) --anchor;
    if (!std::isfinite(score(anchor))) return false;
    for (std::size_t g = 0; g < positions.size(); ++g) {
      const auto first = std::lower_bound(positions[g].begin(), positions[g].end(), anchor);
      if (first != positions[g].end()) {
        const double first_score = score(*first);
        if (!std::isfinite(first_score)) return false;
        right[g] = std::exp(-delta[g] * first_score - 0.5L * delta[g] * delta[g]);
        if (!std::isfinite(right[g]) || right[g] > 1.0L) return false;
      }
      if (first != positions[g].begin()) {
        const double first_score = score(*(first - 1) + 1);
        if (!std::isfinite(first_score)) return false;
        left[g] = std::exp(delta[g] * first_score - 0.5L * delta[g] * delta[g]);
        if (!std::isfinite(left[g]) || left[g] > 1.0L) return false;
      }
    }
    const long double initial = std::exp(std::log(weight) +
      cpp_selnorm_affine_normal_lpdf_log_scale(z[anchor], sei, mean, sd, log_scale));
    if (!std::isfinite(initial)) return false;
    density[anchor] += initial * z_weight[anchor];
    long double value = initial;
    // f(s + d) / f(s) = exp(-s*d - d*d/2). For another occurrence
    // of that exact gap, the ratio changes by exp(-d * distance).
    for (int i = anchor + 1; i < size; ++i) {
      const int g = gap_index[i - 1];
      value *= right[g];
      density[i] += value * z_weight[i];
      right[g] *= forward[i - 1];
    }
    value = initial;
    for (int i = anchor - 1; i >= 0; --i) {
      const int g = gap_index[i];
      value *= left[g];
      density[i] += value * z_weight[i];
      left[g] *= backward[i];
    }
    return true;
  }
};

// Project positive Gaussian mixtures without normalizing their coefficients.

struct SelNormGaussianMixtureRow {
  const double *mean;
  const double *log_coefficient;
  std::size_t count;
  double sd;
  double sei;
};

struct SelNormMixtureProjectionWorkspace {
  const double *z;
  int size;
  NormalProjectionGrid recurrence;
  std::vector<double> weight;
  std::vector<double> log_weight;
  double maximum_weight;
  std::vector<long double> row_density;

  static int checked_size(const double *grid, int length, const double *omega_on_grid)
  {
    if (grid == nullptr || omega_on_grid == nullptr || length < 1) {
      throw std::invalid_argument("Gaussian mixture projection requires a nonempty grid.");
    }
    return length;
  }

  SelNormMixtureProjectionWorkspace(const double *grid, int length,
                                   const double *omega_on_grid)
    : z(grid), size(checked_size(grid, length, omega_on_grid)), recurrence(grid, size), weight(size),
      log_weight(size), maximum_weight(0.0), row_density(size)
  {
    for (int point = 0; point < size; ++point) {
      if (!std::isfinite(z[point]) || !std::isfinite(omega_on_grid[point]) ||
          omega_on_grid[point] < 0.0) {
        throw std::invalid_argument("Gaussian mixture grid and step weights are invalid.");
      }
      weight[point] = omega_on_grid[point];
      log_weight[point] = std::log(weight[point]);
      maximum_weight = std::max(maximum_weight, weight[point]);
    }
  }

  // Worker copy: same grid and step weights, its own recurrence and buffers.
  SelNormMixtureProjectionWorkspace(const SelNormMixtureProjectionWorkspace &other)
    : SelNormMixtureProjectionWorkspace(other.z, other.size, other.weight.data()) {}
  SelNormMixtureProjectionWorkspace &operator=(
    const SelNormMixtureProjectionWorkspace &) = delete;
};

// Each view has physical-coordinate means and one fixed physical SD/original SE.
// Its coefficients already contain 1/K, quadrature/context/importance weights,
// inverse A, and any sample reduction required by the caller. Never divide by
// count, number of views, sum of weights, or number of proposals here.
// Returns a pooled z-density and an all-grid absolute compression-error bound.
// Allowance zero disables compression but still permits the exact Gaussian
// recurrence; a caller needing its historical raw-loop mode should branch before
// collecting mixtures. Probability/CDF projections do not call this helper.
// One mixture row's contribution, left in the workspace's row buffer, with its
// compression error and omitted mass reported separately. The caller reduces
// the rows in their original order, so a threaded evaluation accumulates
// exactly what a serial one does.
inline void selnorm_evaluate_gaussian_mixture_row(
    SelNormGaussianMixtureRow const &row, double compression_allowance,
    long double row_allowance, SelNormMixtureProjectionWorkspace &workspace,
    long double &error, long double &omitted_mass)
{
  error = 0.0L;
  omitted_mass = 0.0L;
  const double *means = row.mean;
  const double *logs = row.log_coefficient;
  std::size_t count = row.count;
  SelNormMixtureCompression compressed;
  const bool grouped = compression_allowance > 0.0 && workspace.size >= 16 &&
    selnorm_compress_gaussian_mixture(means, logs, count, row.sd,
      row_allowance / row.sei, &compressed);
  // The moment planner keeps the original input arrays and source indices.
  // Its error covers both the polynomial and explicit bounded tail omission;
  // unsafe points use the original Gaussian sum, never a clipped polynomial.
  const long double common_scale = grouped ? std::exp(compressed.log_scale) : 0.0L;
  if (grouped && common_scale > 0.0L && std::isfinite(common_scale)) {
    omitted_mass = compressed.input_omission_mass;
    error = compressed.absolute_error * row.sei * workspace.maximum_weight;
    std::fill(workspace.row_density.begin(), workspace.row_density.end(), 0.0L);
    workspace.recurrence.prepare(row.sei, row.sd);
    const double log_scale = std::log(row.sei) - std::log(row.sd);
    std::vector<double> unit_weight(workspace.size, 1.0);
    std::vector<long double> gaussian(workspace.size);
    for (SelNormGaussianMomentGroup const &group : compressed.groups) {
      std::fill(gaussian.begin(), gaussian.end(), 0.0L);
      // Singleton points always use the original component below; avoid
      // preparing a recurrence whose output cannot be consumed.
      bool recurrence_ok = group.source_count > 1 &&
        workspace.recurrence.add(group.mean, row.sei, row.sd,
          log_scale, common_scale, unit_weight, gaussian);
      if (recurrence_ok) {
        for (long double value : gaussian) {
          if (!std::isfinite(value) || value < 0.0L) recurrence_ok = false;
        }
      }
      for (int grid_point = 0; grid_point < workspace.size; ++grid_point) {
        if (!(workspace.weight[grid_point] > 0.0)) continue;
        const long double t = static_cast<long double>(selnorm_fma(
          workspace.z[grid_point], row.sei, -group.mean)) / row.sd;
        long double polynomial = 0.0L;
        SelNormGaussianMomentPoint action = selnorm_gaussian_moment_point(group, t, &polynomial);
        if (action == SelNormGaussianMomentPoint::omitted) continue;
        if (action == SelNormGaussianMomentPoint::polynomial) {
          const long double base_density = recurrence_ok ? gaussian[grid_point] :
            std::exp(compressed.log_scale + cpp_selnorm_affine_normal_lpdf_log_scale(
              workspace.z[grid_point], row.sei, group.mean, row.sd, log_scale));
          const long double value = base_density * polynomial * workspace.weight[grid_point];
          if (value > 0.0L && std::isfinite(value)) {
            workspace.row_density[grid_point] += value;
            continue;
          }
        }
        // Original group summation also handles extreme scaling and any
        // polynomial range/positivity failure. No extra coefficient factor.
        for (std::size_t index = group.source_begin;
             index < group.source_begin + group.source_count; ++index) {
          const std::size_t original = compressed.source_index[index];
          const long double log_value = static_cast<long double>(logs[original]) +
            workspace.log_weight[grid_point] + cpp_selnorm_affine_normal_lpdf_log_scale(
              workspace.z[grid_point], row.sei, means[original], row.sd, log_scale);
          workspace.row_density[grid_point] += std::exp(log_value);
        }
      }
    }
    return;
  }
  // Original outcome units retain the affine Gaussian FMA calculation. SE
  // appears once in the projected density and once in the error bound.
  std::fill(workspace.row_density.begin(), workspace.row_density.end(), 0.0L);
  workspace.recurrence.prepare(row.sei, row.sd);
  const double log_scale = std::log(row.sei) - std::log(row.sd);
  bool recurrence_ok = true;
  for (std::size_t point = 0; point < count; ++point) {
    const long double coefficient = std::exp(static_cast<long double>(logs[point]));
    if ((coefficient == 0.0L && std::isfinite(logs[point])) ||
        !workspace.recurrence.add(means[point], row.sei, row.sd, log_scale,
          coefficient, workspace.weight, workspace.row_density)) {
      recurrence_ok = false;
      break;
    }
  }
  if (recurrence_ok) {
    for (long double value : workspace.row_density) {
      if (!std::isfinite(value)) recurrence_ok = false;
    }
  }
  if (!recurrence_ok) {
    std::fill(workspace.row_density.begin(), workspace.row_density.end(), 0.0L);
    for (std::size_t point = 0; point < count; ++point) {
      for (int grid_point = 0; grid_point < workspace.size; ++grid_point) {
        if (!(workspace.weight[grid_point] > 0.0)) continue;
        const long double log_value = static_cast<long double>(logs[point]) +
          workspace.log_weight[grid_point] + cpp_selnorm_affine_normal_lpdf_log_scale(
            workspace.z[grid_point], row.sei, means[point], row.sd, log_scale);
        workspace.row_density[grid_point] += std::exp(log_value);
      }
    }
  }
}


inline void selnorm_evaluate_gaussian_mixtures(
    std::vector<SelNormGaussianMixtureRow> const &rows,
    double compression_allowance,
    SelNormMixtureProjectionWorkspace &workspace,
    std::vector<long double> &density,
    long double &absolute_error, long double &input_omitted_mass,
    int threads)
{
  if (!std::isfinite(compression_allowance) || compression_allowance < 0.0) {
    throw std::invalid_argument("Gaussian mixture compression allowance is invalid.");
  }
  for (SelNormGaussianMixtureRow const &row : rows) {
    if (!(row.sd > 0.0) || !std::isfinite(row.sd) ||
        !(row.sei > 0.0) || !std::isfinite(row.sei) ||
        (row.count > 0 && (row.mean == nullptr || row.log_coefficient == nullptr))) {
      throw std::invalid_argument("Gaussian mixture row parameters are invalid.");
    }
  }
  density.assign(workspace.size, 0.0L);
  absolute_error = 0.0L;
  input_omitted_mass = 0.0L;
  if (rows.empty() || !(workspace.maximum_weight > 0.0)) return;
  const long double row_allowance = static_cast<long double>(compression_allowance) /
    (static_cast<long double>(rows.size()) * workspace.maximum_weight);

  // Rows are independent: each one projects its own mixture onto the whole
  // grid. Only the reduction below is ordered, so the threaded result is the
  // serial result. A row's own buffer is kept and summed afterwards.
  const int count = static_cast<int>(rows.size());
  const int workers = threads > 1 && count > 1 ? std::min(threads, count) : 1;
  std::vector<long double> row_error(rows.size(), 0.0L);
  std::vector<long double> row_omitted(rows.size(), 0.0L);
  const auto accumulate = [&](int index) {
    if (row_omitted[index] > 0) {
      input_omitted_mass = selnorm_mixture_compression_detail::round_up(
        input_omitted_mass + row_omitted[index]);
    }
    if (row_error[index] > 0) {
      absolute_error = selnorm_mixture_compression_detail::round_up(
        absolute_error + row_error[index]);
    }
  };
  if (workers <= 1) {
    // Serial rows accumulate straight into the output, as before.
    for (int index = 0; index < count; ++index) {
      if (rows[index].count == 0) continue;
      selnorm_evaluate_gaussian_mixture_row(rows[index], compression_allowance,
        row_allowance, workspace, row_error[index], row_omitted[index]);
      accumulate(index);
      for (int grid_point = 0; grid_point < workspace.size; ++grid_point) {
        density[grid_point] += workspace.row_density[grid_point];
      }
    }
    return;
  }

  std::vector<std::vector<long double>> row_density(rows.size());
  std::atomic<bool> failed{false};
  const auto evaluate = [&](int index, SelNormMixtureProjectionWorkspace &scratch) {
    if (rows[index].count == 0) return;
    selnorm_evaluate_gaussian_mixture_row(rows[index], compression_allowance,
      row_allowance, scratch, row_error[index], row_omitted[index]);
    row_density[index] = scratch.row_density;
  };
#if defined(_OPENMP)
  {
    #pragma omp parallel num_threads(workers)
    {
      SelNormMixtureProjectionWorkspace scratch(workspace);
      #pragma omp for schedule(dynamic, 1)
      for (int index = 0; index < count; ++index) {
        if (failed.load(std::memory_order_relaxed)) continue;
        try {
          evaluate(index, scratch);
        } catch (...) {
          failed.store(true, std::memory_order_relaxed);
        }
      }
    }
    if (failed.load(std::memory_order_relaxed)) {
      throw std::runtime_error("Gaussian mixture row projection failed.");
    }
  }
#else
  for (int index = 0; index < count; ++index) evaluate(index, workspace);
#endif

  for (int index = 0; index < count; ++index) {
    if (row_density[index].empty()) continue;
    accumulate(index);
    for (int grid_point = 0; grid_point < workspace.size; ++grid_point) {
      density[grid_point] += row_density[index][grid_point];
    }
  }
}

// Deterministic quadrature over a forest of factor supports. The factors are
// sorted by decreasing support size; every pair of supports must then be
// nested or disjoint, which is what block-constant covariance structures on a
// nesting tree produce. A chain (each support inside the previous one) is the
// special case with one child per level and keeps its previous arithmetic.
//
// Each factor integrates over the tensor of its own axis and the axes of its
// ancestors only, so a root with disjoint children costs one rule per child
// rather than one full rank-dimensional tensor.
// Structural analysis of the factor supports, built once per likelihood state.
// When the supports form a forest the nested rule expands one axis per tree
// level, so its cost is `rows x order^(depth + 1)` and does not grow with the
// factor rank. The tensor fallback costs `order^rank`, which is what the rank
// cap used to encode.
struct FactorForestPlan {
  bool valid = false;
  int effective_rank = 0;
  int max_depth = 0;
  std::vector<int> permutation;
  std::vector<int> parent;
  std::vector<int> depth;
  std::vector<std::vector<int> > children;
  std::vector<int> roots;
  std::vector<std::vector<int> > axes;
  std::vector<std::vector<int> > node_rows;
  std::vector<int> free_rows;
};

bool factor_forest_plan(const FactorNormalizerContext &context,
                        const int *declared_support, FactorForestPlan *out)
{
  out->valid = false;
  if (context.vector_rule != SELVECTOR_PRODUCT) return false;
  const auto supported = [&context, declared_support](int row, int factor) {
    const int index = row + context.dimension * factor;
    return declared_support == nullptr ? context.loading[index] != 0.0 :
      declared_support[index] != 0;
  };
  std::vector<int> support_size(static_cast<std::size_t>(context.rank), 0);
  std::vector<int> &permutation = out->permutation;
  permutation.assign(static_cast<std::size_t>(context.rank), 0);
  for (int factor = 0; factor < context.rank; ++factor) {
    permutation[factor] = factor;
    for (int i = 0; i < context.dimension; ++i) {
      if (supported(i, factor)) {
        ++support_size[factor];
      }
    }
  }
  std::stable_sort(
    permutation.begin(), permutation.end(),
    [&support_size](int left, int right) {
      return support_size[left] > support_size[right];
    }
  );

  int effective_rank = context.rank;
  while (effective_rank > 0 &&
         support_size[permutation[effective_rank - 1]] == 0) {
    --effective_rank;
  }
  out->effective_rank = effective_rank;

  // Supports must form a forest: the closest containing support is a factor
  // parent, and a partial overlap is rejected.
  const std::size_t node_count = static_cast<std::size_t>(effective_rank);
  std::vector<int> &parent = out->parent;
  std::vector<int> &depth = out->depth;
  std::vector<std::vector<int> > &children = out->children;
  std::vector<int> &roots = out->roots;
  parent.assign(node_count, -1);
  depth.assign(node_count, 0);
  children.assign(node_count, std::vector<int>());
  roots.clear();
  for (int q = 0; q < effective_rank; ++q) {
    for (int p = 0; p < q; ++p) {
      bool inside = true;
      bool disjoint = true;
      for (int i = 0; i < context.dimension; ++i) {
        if (!supported(i, permutation[q])) continue;
        if (supported(i, permutation[p])) disjoint = false;
        else inside = false;
      }
      if (disjoint) continue;
      if (!inside) return false;
      parent[q] = p;
    }
    if (parent[q] < 0) {
      roots.push_back(q);
    } else {
      depth[q] = depth[parent[q]] + 1;
      children[parent[q]].push_back(q);
    }
  }

  // The axes of a factor are its ancestors, root first, then itself.
  std::vector<std::vector<int> > &axes = out->axes;
  axes.assign(node_count, std::vector<int>());
  for (int q = 0; q < effective_rank; ++q) {
    if (parent[q] >= 0) axes[q] = axes[parent[q]];
    axes[q].push_back(q);
  }

  // Every row belongs to the deepest factor supporting it, and that factor
  // axes must be exactly the factors supporting the row.
  std::vector<std::vector<int> > &node_rows = out->node_rows;
  std::vector<int> &free_rows = out->free_rows;
  node_rows.assign(node_count, std::vector<int>());
  free_rows.clear();
  for (int i = 0; i < context.dimension; ++i) {
    int owner = -1;
    int active = 0;
    for (int q = 0; q < effective_rank; ++q) {
      if (!supported(i, permutation[q])) continue;
      ++active;
      if (owner < 0 || depth[q] > depth[owner]) owner = q;
    }
    if (owner < 0) {
      free_rows.push_back(i);
      continue;
    }
    if (active != static_cast<int>(axes[owner].size())) return false;
    node_rows[owner].push_back(i);
  }

  out->max_depth = 0;
  for (int q = 0; q < effective_rank; ++q) {
    if (depth[q] > out->max_depth) out->max_depth = depth[q];
  }
  out->valid = true;
  return true;
}

bool factor_nested_rule_log_integral(
    const FactorNormalizerContext &context, const FactorForestPlan &plan,
    const double *nodes, const double *log_weights, int offset, int order,
    double *result, SelNormZProjection *projection = nullptr)
{
  if (!plan.valid) return false;
  const int effective_rank = plan.effective_rank;
  const std::vector<int> &permutation = plan.permutation;
  const std::vector<int> &parent = plan.parent;
  const std::vector<int> &depth = plan.depth;
  const std::vector<std::vector<int> > &children = plan.children;
  const std::vector<int> &roots = plan.roots;
  const std::vector<std::vector<int> > &axes = plan.axes;
  const std::vector<std::vector<int> > &node_rows = plan.node_rows;
  const std::vector<int> &free_rows = plan.free_rows;
  const std::size_t node_count = static_cast<std::size_t>(effective_rank);

  std::vector<std::size_t> power(node_count + 2, 1);
  for (std::size_t j = 1; j < power.size(); ++j) {
    power[j] = power[j - 1] * static_cast<std::size_t>(order);
  }

  std::vector<std::vector<long double>> own(node_count);
  std::vector<std::vector<long double>> integral(node_count);
  std::vector<std::vector<double>> row_means;
  std::vector<std::vector<long double>> row_normalizers;
  if (projection != nullptr) {
    row_means.resize(static_cast<std::size_t>(context.dimension));
    row_normalizers.resize(static_cast<std::size_t>(context.dimension));
  }
  const std::vector<int> no_axes;

  const auto accumulate_rows = [&](const std::vector<int> &rows,
                                   std::size_t size,
                                   const std::vector<int> &row_axes,
                                   std::vector<long double> &target) {
    std::vector<double> conditional_means(size);
    for (int row : rows) {
      conditional_means[0] = context.mean[row];
      for (std::size_t position = 0; position < row_axes.size(); ++position) {
        const double coefficient = context.loading[
          row + context.dimension * permutation[row_axes[position]]
        ];
        const std::size_t lower_size = power[position];
        // Expand one axis at a time. Visit node zero last so the prefix is
        // still available while filling the other slices of the same buffer.
        for (int point = order - 1; point >= 0; --point) {
          const double shift = coefficient * nodes[offset + point];
          for (std::size_t lower = 0; lower < lower_size; ++lower) {
            conditional_means[
              lower + static_cast<std::size_t>(point) * lower_size
            ] = conditional_means[lower] + shift;
          }
        }
      }
      if (projection != nullptr) {
        row_means[row] = conditional_means;
        row_normalizers[row].assign(size, 1.0L);
      }
      if (!cpp_selnorm_step_normalizer_product(
            conditional_means.data(), size, context.residual_sd[row],
            context.selection_se[row], context.omega, context.selection,
            projection == nullptr ? target.data() :
              row_normalizers[row].data())) return false;
      if (projection != nullptr) {
        for (std::size_t point = 0; point < size; ++point) {
          target[point] *= row_normalizers[row][point];
        }
      }
    }
    return true;
  };

  std::vector<long double> base(1, 1.0L);
  if (!accumulate_rows(free_rows, 1, no_axes, base)) return false;

  std::vector<long double> weights(static_cast<std::size_t>(order));
  for (int point = 0; point < order; ++point) {
    weights[point] =
      std::exp(static_cast<long double>(log_weights[offset + point]));
  }

  std::vector<int> by_depth(node_count);
  for (int q = 0; q < effective_rank; ++q) by_depth[q] = q;
  std::stable_sort(by_depth.begin(), by_depth.end(),
                   [&depth](int left, int right) {
                     return depth[left] > depth[right];
                   });

  for (int q : by_depth) {
    const std::size_t prefix = power[depth[q]];
    const std::size_t size = power[depth[q] + 1];
    own[q].assign(size, 1.0L);
    if (!accumulate_rows(node_rows[q], size, axes[q], own[q])) return false;
    for (int child : children[q]) {
      for (std::size_t point = 0; point < size; ++point) {
        own[q][point] *= integral[child][point];
      }
    }
    integral[q].assign(prefix, 0.0L);
    for (std::size_t lower = 0; lower < prefix; ++lower) {
      for (int point = 0; point < order; ++point) {
        integral[q][lower] += weights[point] * own[q][
          lower + static_cast<std::size_t>(point) * prefix
        ];
      }
      if (!(integral[q][lower] > 0.0L) ||
          !std::isfinite(integral[q][lower])) {
        return false;
      }
    }
  }

  long double total = base[0];
  for (int root : roots) total *= integral[root][0];
  if (!(total > 0.0L) || !std::isfinite(total)) return false;
  *result = static_cast<double>(std::log(total));

  if (projection != nullptr) {
    // The upward messages already integrate every deeper factor. Their
    // conditional quadrature weights project each row using only its own
    // axes, without rebuilding a full tensor for every density-grid point.
    std::vector<long double> density(projection->size, 0.0L);
    std::vector<double> z_weight(projection->size, context.omega[0]);
    if (!projection->probability) {
      for (int z = 0; z < projection->size; ++z) {
        const double signed_z = context.selection.effect_sign * projection->z[z];
        for (int bin = 0; bin < context.n_bins; ++bin) {
          if (signed_z >= context.selection.z_lower[bin] &&
              signed_z <= context.selection.z_upper[bin]) {
            z_weight[z] = context.omega[bin];
            break;
          }
        }
      }
    }
    NormalProjectionGrid recurrence(projection->z,
      projection->probability ? 0 : projection->size);
    std::vector<std::vector<long double>> marginal(node_count);
    const std::vector<long double> unit(1, 1.0L);

    const auto project_rows = [&](const std::vector<int> &rows,
                                  const std::vector<long double> &weight,
                                  std::size_t size) {
      for (int row : rows) {
        const double log_scale = std::log(context.selection_se[row]) -
          std::log(context.residual_sd[row]);
        recurrence.prepare(context.selection_se[row], context.residual_sd[row]);
        for (std::size_t point = 0; point < size; ++point) {
          const long double scaled = weight[point] / row_normalizers[row][point];
          if (recurrence.add(row_means[row][point], context.selection_se[row],
                context.residual_sd[row], log_scale, scaled, z_weight, density)) continue;
          for (int z = 0; z < projection->size; ++z) {
            if (projection->probability) {
              double inverse;
              density[z] += weight[point] * cpp_selnorm_kernel_threshold(
                projection->z[z], row_means[row][point], context.residual_sd[row],
                context.selection_se[row], context.omega, 0, 0,
                context.kernel_mode, context.selection, &inverse, 1, false
              );
            } else if (z_weight[z] > 0.0) {
              const double value = cpp_selnorm_affine_normal_pdf_log_scale(
                projection->z[z], context.selection_se[row], row_means[row][point],
                context.residual_sd[row], log_scale
              );
              const long double extended_value = value > 0.0 ?
                static_cast<long double>(value) :
                std::exp(static_cast<long double>(
                  cpp_selnorm_affine_normal_lpdf_log_scale(
                    projection->z[z], context.selection_se[row],
                    row_means[row][point], context.residual_sd[row], log_scale
                  )
                ));
              density[z] += scaled * z_weight[z] * extended_value;
            }
          }
        }
      }
    };

    project_rows(free_rows, unit, 1);
    std::vector<int> visit(node_count);
    for (int q = 0; q < effective_rank; ++q) visit[q] = q;
    std::stable_sort(visit.begin(), visit.end(),
                     [&depth](int left, int right) {
                       return depth[left] < depth[right];
                     });
    for (int q : visit) {
      const std::size_t prefix = power[depth[q]];
      const std::size_t size = power[depth[q] + 1];
      const std::vector<long double> &upper =
        parent[q] < 0 ? unit : marginal[parent[q]];
      marginal[q].assign(size, 0.0L);
      for (int point = 0; point < order; ++point) {
        for (std::size_t lower = 0; lower < prefix; ++lower) {
          const std::size_t index =
            lower + static_cast<std::size_t>(point) * prefix;
          marginal[q][index] = upper[lower] * weights[point] *
            own[q][index] / integral[q][lower];
        }
      }
      project_rows(node_rows[q], marginal[q], size);
    }
    for (int z = 0; z < projection->size; ++z) {
      projection->density[z] = static_cast<double>(density[z] / context.dimension);
      if (!std::isfinite(projection->density[z])) return false;
    }
    projection->relative_error = 0.0;
  }
  return std::isfinite(*result);
}

double factor_rule_log_integral(const FactorNormalizerContext &context,
                                const FactorForestPlan &plan,
                                const double *nodes,
                                const double *log_weights,
                                int offset, int order)
{
  double nested = 0.0;
  if (factor_nested_rule_log_integral(
        context, plan, nodes, log_weights, offset, order, &nested
      )) {
    return nested;
  }
  // A forest plan that fails numerically leaves only the tensor rule, which
  // the budget may refuse at this rank. Declining the rung keeps the ladder
  // and its randomized fallback in charge rather than spending the nodes.
  if (!selnorm_factor_rule_affordable(order, context.rank)) {
    return -std::numeric_limits<double>::infinity();
  }
  std::size_t total = 1;
  for (int factor = 0; factor < context.rank; ++factor) {
    total *= static_cast<std::size_t>(order);
  }

  // Factor quadrature only reaches this kernel for an ordered step partition.
  // Accumulating its positive normalizers avoids a log/exp pair for every row
  // and quadrature point. Retain the log-scale path below for rare numerical
  // fallback cases.
  std::vector<long double> weights(static_cast<std::size_t>(order));
  for (int j = 0; j < order; ++j) {
    weights[static_cast<std::size_t>(j)] =
      std::exp(static_cast<long double>(log_weights[offset + j]));
  }
  long double direct_sum = 0.0L;
  bool direct_ok = context.vector_rule == SELVECTOR_PRODUCT;
  std::vector<double> latent(static_cast<std::size_t>(context.rank));
  for (std::size_t point = 0; point < total && direct_ok; ++point) {
    std::size_t remaining = point;
    long double term = 1.0L;
    for (int factor = 0; factor < context.rank; ++factor) {
      const int index = static_cast<int>(remaining % order);
      remaining /= static_cast<std::size_t>(order);
      latent[factor] = nodes[offset + index];
      term *= weights[static_cast<std::size_t>(index)];
    }
    for (int i = 0; i < context.dimension; ++i) {
      double conditional_mean = context.mean[i];
      for (int factor = 0; factor < context.rank; ++factor) {
        conditional_mean +=
          context.loading[i + context.dimension * factor] * latent[factor];
      }
      double omega_last = 0.0;
      double normalizer = 0.0;
      direct_ok = cpp_selnorm_step_cdf_telescope_plan(
        conditional_mean, context.residual_sd[i], context.selection_se[i],
        context.omega, context.selection, nullptr, nullptr, &omega_last,
        &normalizer, 1, false
      );
      if (!direct_ok) break;
      term *= static_cast<long double>(normalizer);
    }
    direct_sum += term;
    direct_ok = direct_ok && std::isfinite(direct_sum) && direct_sum > 0.0L;
  }
  if (direct_ok) return static_cast<double>(std::log(direct_sum));

  double out = -std::numeric_limits<double>::infinity();
  for (std::size_t point = 0; point < total; ++point) {
    std::size_t remaining = point;
    double log_weight = 0.0;
    double normal_kernel = 0.0;
    for (int factor = 0; factor < context.rank; ++factor) {
      const int index = offset + static_cast<int>(remaining % order);
      remaining /= static_cast<std::size_t>(order);
      latent[factor] = nodes[index];
      normal_kernel += 0.5 * latent[factor] * latent[factor];
      log_weight += log_weights[index];
    }
    const double log_integrand = factor_log_integrand(
      context, latent.data(), nullptr
    );
    if (std::isfinite(log_integrand)) {
      out = log_add_exp(out, log_weight + log_integrand + normal_kernel);
    }
  }
  return out;
}

double vector_dot(const std::vector<double> &x,
                  const std::vector<double> &y)
{
  double out = 0.0;
  for (std::size_t i = 0; i < x.size(); ++i) out += x[i] * y[i];
  return out;
}

void factor_optimize_mode(const FactorNormalizerContext &context,
                          std::vector<double> *position,
                          std::vector<double> *proposal_covariance = nullptr)
{
  const int rank = context.rank;
  std::vector<double> gradient(static_cast<std::size_t>(rank));
  double value = factor_log_integrand(
    context, position->data(), gradient.data()
  );
  if (!std::isfinite(value)) return;

  std::vector<double> inverse_hessian(
    static_cast<std::size_t>(rank * rank), 0.0
  );
  for (int factor = 0; factor < rank; ++factor) {
    inverse_hessian[factor + rank * factor] = 1.0;
  }
  std::vector<double> best = *position;
  std::vector<double> best_covariance;
  if (proposal_covariance != nullptr) best_covariance = inverse_hessian;
  double best_value = value;
  const double gradient_tolerance =
    std::sqrt(std::numeric_limits<double>::epsilon());

  for (int iteration = 0; iteration < 100; ++iteration) {
    const double gradient_norm = std::sqrt(vector_dot(gradient, gradient));
    double position_norm = 0.0;
    for (int factor = 0; factor < rank; ++factor) {
      position_norm += (*position)[factor] * (*position)[factor];
    }
    position_norm = std::sqrt(position_norm);
    if (gradient_norm <= gradient_tolerance * (1.0 + position_norm)) break;

    std::vector<double> direction(static_cast<std::size_t>(rank), 0.0);
    for (int column = 0; column < rank; ++column) {
      for (int row = 0; row < rank; ++row) {
        direction[row] += inverse_hessian[row + rank * column] *
          gradient[column];
      }
    }
    double directional_derivative = vector_dot(gradient, direction);
    if (!(directional_derivative > 0.0) ||
        !std::isfinite(directional_derivative)) {
      direction = gradient;
      directional_derivative = vector_dot(gradient, direction);
    }

    std::vector<double> candidate(static_cast<std::size_t>(rank));
    std::vector<double> candidate_gradient(static_cast<std::size_t>(rank));
    double candidate_value = -std::numeric_limits<double>::infinity();
    double step = 1.0;
    bool accepted = false;
    for (int line_search = 0; line_search < 30; ++line_search) {
      for (int factor = 0; factor < rank; ++factor) {
        candidate[factor] = (*position)[factor] + step * direction[factor];
      }
      candidate_value = factor_log_integrand(
        context, candidate.data(), candidate_gradient.data()
      );
      if (std::isfinite(candidate_value) &&
          candidate_value >= value + 1e-4 * step * directional_derivative) {
        accepted = true;
        break;
      }
      step *= 0.5;
    }
    if (!accepted) break;

    std::vector<double> displacement(static_cast<std::size_t>(rank));
    std::vector<double> curvature(static_cast<std::size_t>(rank));
    for (int factor = 0; factor < rank; ++factor) {
      displacement[factor] = candidate[factor] - (*position)[factor];
      curvature[factor] = gradient[factor] - candidate_gradient[factor];
    }
    const double displacement_curvature =
      vector_dot(displacement, curvature);
    if (displacement_curvature > 0.0 &&
        std::isfinite(displacement_curvature)) {
      std::vector<double> h_curvature(static_cast<std::size_t>(rank), 0.0);
      for (int column = 0; column < rank; ++column) {
        for (int row = 0; row < rank; ++row) {
          h_curvature[row] += inverse_hessian[row + rank * column] *
            curvature[column];
        }
      }
      const double curvature_h_curvature =
        vector_dot(curvature, h_curvature);
      const double rho = 1.0 / displacement_curvature;
      const double scale =
        (1.0 + curvature_h_curvature * rho) * rho;
      for (int column = 0; column < rank; ++column) {
        for (int row = 0; row < rank; ++row) {
          inverse_hessian[row + rank * column] +=
            scale * displacement[row] * displacement[column] -
            rho * (displacement[row] * h_curvature[column] +
                   h_curvature[row] * displacement[column]);
        }
      }
    } else {
      std::fill(inverse_hessian.begin(), inverse_hessian.end(), 0.0);
      for (int factor = 0; factor < rank; ++factor) {
        inverse_hessian[factor + rank * factor] = 1.0;
      }
    }

    *position = candidate;
    gradient = candidate_gradient;
    value = candidate_value;
    if (value > best_value) {
      best = *position;
      best_value = value;
      if (proposal_covariance != nullptr) best_covariance = inverse_hessian;
    }
  }
  *position = best;
  if (proposal_covariance != nullptr) *proposal_covariance = best_covariance;
}

bool factor_solve_cholesky(const std::vector<double> &cholesky, int rank,
                           std::vector<double> *value)
{
  for (int row = 0; row < rank; ++row) {
    for (int column = 0; column < row; ++column) {
      (*value)[row] -= cholesky[row + rank * column] * (*value)[column];
    }
    const double diagonal = cholesky[row + rank * row];
    if (!(diagonal > 0.0) || !std::isfinite(diagonal)) return false;
    (*value)[row] /= diagonal;
  }
  for (int row = rank - 1; row >= 0; --row) {
    for (int column = row + 1; column < rank; ++column) {
      (*value)[row] -= cholesky[column + rank * row] * (*value)[column];
    }
    (*value)[row] /= cholesky[row + rank * row];
  }
  return true;
}

double factor_log_mean(const double *log_values, int count)
{
  double out = -std::numeric_limits<double>::infinity();
  for (int i = 0; i < count; ++i) out = log_add_exp(out, log_values[i]);
  return out - std::log(static_cast<double>(count));
}

double factor_log_mean(const std::vector<double> &log_values, int count)
{
  return factor_log_mean(log_values.data(), count);
}

}

struct SelNormNormalProjectionWorkspace::Impl {
  NormalProjectionGrid recurrence;
  std::vector<double> weights;
  std::vector<long double> output;

  Impl(const double *z, int size) : recurrence(z, size), weights(size, 1.0), output(size) {}
};

SelNormNormalProjectionWorkspace::SelNormNormalProjectionWorkspace(const double *z, int size)
  : impl(new Impl(z, size)) {}

SelNormNormalProjectionWorkspace::~SelNormNormalProjectionWorkspace() = default;

bool SelNormNormalProjectionWorkspace::density(
    const double *mean, int count, double sd, double sei,
    const double *log_coefficients, long double *density)
{
  impl->recurrence.prepare(sei, sd);
  std::fill(impl->output.begin(), impl->output.end(), 0.0L);
  const double log_scale = std::log(sei) - std::log(sd);
  for (int point = 0; point < count; ++point) {
    const long double weight = std::exp(static_cast<long double>(log_coefficients[point]));
    if (!impl->recurrence.add(mean[point], sei, sd, log_scale, weight,
                             impl->weights, impl->output)) {
      return false;
    }
  }
  std::copy(impl->output.begin(), impl->output.end(), density);
  return true;
}

int cpp_selnorm_factor_box_rng(
    const double *mean, const double *supplied, int n,
    const double *selection_se, const double *omega, int bins,
    const double *z_lower, const double *z_upper, int effect_sign,
    int max_attempts, double (*uniform)(), double *output)
{
  const bool profile_box = std::getenv("ROBMA_DEBUG_FACTOR_BOX") != nullptr;
  const auto profile_tick = []() { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); };
  const double profile_start = profile_box ? profile_tick() : 0.0;
  double profile_geometry = 0.0, profile_correction = 0.0, profile_integrands = 0.0, profile_intervals = 0.0, profile_totals = 0.0, profile_sampling = 0.0, profile_late_preparation = 0.0;
  std::uint64_t profile_queries = 0, profile_hits = 0, profile_interval_queries = 0;
  std::uint64_t preparation_work = 0;
  double log_lower = 0.0, log_upper = 0.0;
  const double infinity = std::numeric_limits<double>::infinity();
  const double negative_infinity = -infinity;
  const long double long_infinity = std::numeric_limits<long double>::infinity();
  if (n < 2 || bins < 1 || max_attempts < 1 || uniform == nullptr) return -1;
  bool increasing_signed = true, decreasing_signed = true;
  for (int bin = 1; bin < bins; ++bin) {
    increasing_signed = increasing_signed && omega[bin - 1] >= omega[bin];
    decreasing_signed = decreasing_signed && omega[bin - 1] <= omega[bin];
  }
  if (!increasing_signed && !decreasing_signed) return -1;
  const bool increasing = effect_sign > 0 ? increasing_signed : decreasing_signed;
  std::vector<double> covariance(supplied, supplied + static_cast<std::size_t>(n) * n);
  std::vector<double> target_cholesky = covariance;
  int info = 0;
  F77_CALL(dpotrf)("L", &n, target_cholesky.data(), &n, &info FCONE);
  if (info != 0) return -1;
  SelNormCovarianceEnvelope envelope;
  if (!covariance_envelope(covariance, n, selection_se, &envelope)) return -1;
  const int extracted_rank = envelope.groups == 1 ? 1 : envelope.groups + 1;
  std::vector<EnvelopeGroupGeometry> geometry;
  if (!prepare_envelope_geometry(covariance, mean, selection_se, n,
                                 envelope, false, &geometry)) return -1;
  std::vector<double> loading(static_cast<std::size_t>(n) * extracted_rank, 0.0), sd(n);
  std::vector<int> local(envelope.groups, 0);
  for (int row = 0; row < n; ++row) {
    const int group = envelope.type_index[row], position = local[group]++;
    sd[row] = geometry[group].residual_sd[position];
    loading[row] = geometry[group].common_loading[position];
    if (extracted_rank > 1) loading[row + n * (group + 1)] = geometry[group].loading[position];
  }
  // Exactly zero columns do not define a random coordinate. Removing them is
  // algebraic elimination, not numerical rank truncation or covariance repair.
  int rank = 0;
  for (int axis = 0; axis < extracted_rank; ++axis) {
    const double *column = loading.data() + n * axis;
    if (std::all_of(column, column + n, [](double value) { return value == 0.0; })) continue;
    if (rank != axis) std::copy_n(column, n, loading.data() + n * rank);
    ++rank;
  }
  if (rank < 1 || rank > 4) return -1;
  loading.resize(static_cast<std::size_t>(n) * rank);
  if (profile_box) profile_geometry = profile_tick() - profile_start;
  // Radius of each elementary long-double operation. This applies only to
  // the auxiliary covariance construction, never to modification of target C.
  const auto radius = [long_infinity](long double value) {
    if (!std::isfinite(value)) return long_infinity;
    return std::max(std::nextafter(value, long_infinity) - value,
                    value - std::nextafter(value, -long_infinity));
  };
  const auto up_add = [long_infinity](long double a, long double b) {
    return std::nextafter(a + b, long_infinity);
  };
  const auto represented_difference = [&](int row, int column) {
    long double represented = row == column ? static_cast<long double>(sd[row]) * sd[row] : 0.0L;
    long double error = row == column ? radius(represented) : 0.0L;
    for (int factor = 0; factor < rank; ++factor) {
      const long double product = static_cast<long double>(loading[row + n * factor]) *
        loading[column + n * factor];
      error = up_add(error, radius(product));
      represented += product;
      error = up_add(error, radius(represented));
    }
    const long double value = represented - covariance[row + n * column];
    return std::make_pair(value, up_add(error, radius(value)));
  };
  long double delta = 0.0L;
  for (int row = 0; row < n; ++row) {
    long double bound = 0.0L;
    for (int column = 0; column < n; ++column) {
      const auto residual = represented_difference(row, column);
      if (!std::isfinite(residual.first) || !std::isfinite(residual.second)) return -1;
      bound = up_add(bound, up_add(std::fabs(residual.first), residual.second));
    }
    delta = std::max(delta, bound);
  }
  if (!(delta > 0.0L) || !std::isfinite(delta)) return -1;
  for (int row = 0; row < n; ++row) {
    const long double old_square = static_cast<long double>(sd[row]) * sd[row];
    const long double required = up_add(std::nextafter(old_square, long_infinity),
                                       up_add(delta, delta));
    double candidate_sd = std::sqrt(static_cast<double>(required));
    bool adequate = false;
    for (int pass = 0; pass < 8; ++pass) {
      const long double square = static_cast<long double>(candidate_sd) * candidate_sd;
      if (std::isfinite(candidate_sd) && candidate_sd > 0.0 &&
          std::nextafter(square, -long_infinity) >= required) {
        adequate = true;
        break;
      }
      candidate_sd = std::nextafter(candidate_sd, infinity);
    }
    if (!adequate) return -1;
    sd[row] = candidate_sd;
  }
  // Verify domination for the actually emitted double SDs/loadings. A strict
  // diagonal-dominance margin also makes the Woodbury update factor positive.
  std::vector<double> update(static_cast<std::size_t>(n) * n);
  long double minimum_margin = long_infinity;
  for (int row = 0; row < n; ++row) {
    long double off_diagonal = 0.0L, diagonal_lower = 0.0L;
    for (int column = 0; column < n; ++column) {
      const auto value = represented_difference(row, column);
      if (!std::isfinite(value.first) || !std::isfinite(value.second)) return -1;
      update[row + n * column] = static_cast<double>(value.first);
      if (row == column) diagonal_lower = std::nextafter(value.first - value.second, -long_infinity);
      else off_diagonal = up_add(off_diagonal, up_add(std::fabs(value.first), value.second));
    }
    minimum_margin = std::min(minimum_margin, diagonal_lower - off_diagonal);
  }
  if (!(minimum_margin > 0.0L) || !std::isfinite(minimum_margin)) return -1;
  F77_CALL(dpotrf)("L", &n, update.data(), &n, &info FCONE);
  if (info != 0) return -1;
  for (int column = 0; column < n; ++column) {
    for (int row = 0; row < column; ++row) update[row + n * column] = 0.0;
  }
  const double one = 1.0, zero = 0.0;
  const int stride = 1;
  // J=L_C^-1 L_H, K=I+J'J. The correction is a sum of squares, not subtraction
  // of nearly equal inverse quadratic forms.
  F77_CALL(dtrsm)("L", "L", "N", "N", &n, &n, &one,
    target_cholesky.data(), &n, update.data(), &n FCONE FCONE FCONE FCONE);
  if (!std::all_of(update.begin(), update.end(), [](double value) { return std::isfinite(value); })) return -1;
  std::vector<double> correction(static_cast<std::size_t>(n) * n, 0.0);
  for (int row = 0; row < n; ++row) correction[row + n * row] = 1.0;
  F77_CALL(dsyrk)("L", "T", &n, &n, &one, update.data(), &n, &one,
    correction.data(), &n FCONE FCONE);
  F77_CALL(dpotrf)("L", &n, correction.data(), &n, &info FCONE);
  if (info != 0 || !std::all_of(correction.begin(), correction.end(),
      [](double value) { return std::isfinite(value); })) return -1;

  if (profile_box) profile_correction = profile_tick() - profile_start - profile_geometry;
  SelNormKernelData selection = {};
  selection.n_bins = bins;
  selection.effect_sign = effect_sign;
  selection.z_lower = z_lower;
  selection.z_upper = z_upper;
  selection.trusted_step_partition = true;
  selection.telescope_probabilities = true;
  struct Cell {
    std::array<double, 4> lower{}, upper{}, axis_mass{};
    double lower_corner = 0, upper_corner = 0, mass = 0, gap = 0;
  };
  const auto difference = [&](double upper, double lower) {
    if (upper < lower || std::isnan(upper) || std::isnan(lower)) {
      throw std::runtime_error("Factor-box bounds are not ordered");
    }
    if (upper == lower) return negative_infinity;
    if (lower == negative_infinity) return upper;
    return upper + std::log(-std::expm1(lower - upper));
  };
  std::vector<double> location(n);
  const auto integrand = [&](const std::array<double, 4> &point, bool preparation) {
    const double profile_integrand_start = profile_box ? profile_tick() : 0.0;
    for (int row = 0; row < n; ++row) {
      location[row] = mean[row];
      for (int axis = 0; axis < rank; ++axis) {
        if (loading[row + n * axis] != 0.0) location[row] += loading[row + n * axis] * point[axis];
      }
      if (std::isnan(location[row]) || (!preparation && !std::isfinite(location[row]))) {
        throw std::runtime_error("A factor-box mean is non-finite");
      }
    }
    long double product = 1.0L;
    bool direct = true;
    for (int row = 0; row < n && direct; ++row) {
      if (location[row] == infinity) product *= omega[effect_sign > 0 ? 0 : bins - 1];
      else if (location[row] == negative_infinity) product *= omega[effect_sign > 0 ? bins - 1 : 0];
      else {
        direct = cpp_selnorm_step_normalizer_product(location.data() + row, 1, sd[row],
          selection_se[row], omega, selection, &product);
        if (preparation) ++preparation_work;
      }
      direct = direct && product > 0.0L && std::isfinite(product);
    }
    if (direct) {
      const double value = static_cast<double>(std::log(product));
      if (profile_box && preparation) profile_integrands += profile_tick() - profile_integrand_start;
      return value;
    }
    // Unsafe direct products retain the complete existing scalar log path.
    // Preparation cost includes both attempted mass evaluations and fallback work.
    double total = 0.0;
    for (int row = 0; row < n; ++row) {
      double value;
      if (location[row] == infinity) value = std::log(omega[effect_sign > 0 ? 0 : bins - 1]);
      else if (location[row] == negative_infinity) value = std::log(omega[effect_sign > 0 ? bins - 1 : 0]);
      else {
        value = cpp_selnorm_step_log_norm(location[row], sd[row], selection_se[row], omega, selection, 1, false);
        if (preparation) ++preparation_work;
        if (!std::isfinite(value)) throw std::runtime_error("A factor-box normalizer is non-finite");
      }
      total += value;
    }
    if (std::isnan(total)) throw std::runtime_error("A factor-box corner is undefined");
    if (profile_box && preparation) profile_integrands += profile_tick() - profile_integrand_start;
    return total;
  };
  auto corner_hash = [rank](const std::array<double, 4> &point) {
    std::size_t hash = 0;
    for (int axis = 0; axis < rank; ++axis) {
      hash ^= std::hash<double>{}(point[axis]) + static_cast<std::size_t>(0x9e3779b9) +
        (hash << 6) + (hash >> 2);
    }
    return hash;
  };
  std::unordered_map<std::array<double, 4>, std::pair<double, std::uint64_t>,
                     decltype(corner_hash)> corner_cache(0, corner_hash);
  std::unordered_map<std::array<double, 4>, int, decltype(corner_hash)> profile_interval_keys(0, corner_hash);
  const auto corner = [&](const std::array<double, 4> &point) {
    if (profile_box) ++profile_queries;
    std::array<double, 4> key = point;
    std::fill(key.begin() + rank, key.end(), 0.0);
    const auto cached = corner_cache.find(key);
    if (cached != corner_cache.end()) {
      if (profile_box) ++profile_hits;
      // Preserve the original refinement decisions and RNG stream. This is
      // the same corner's logical preparation cost, even when its value is
      // already available from a neighbouring cell.
      preparation_work += cached->second.second;
      return cached->second.first;
    }
    const std::uint64_t before = preparation_work;
    const double value = integrand(point, true);
    corner_cache.emplace(key, std::make_pair(value, preparation_work - before));
    return value;
  };
  const auto low = [&](const Cell &cell) { return increasing ? cell.lower_corner : cell.upper_corner; };
  const auto high = [&](const Cell &cell) { return increasing ? cell.upper_corner : cell.lower_corner; };
  const auto set_mass_gap = [&](Cell &cell) {
    cell.mass = 0.0;
    for (int axis = 0; axis < rank; ++axis) cell.mass += cell.axis_mass[axis];
    cell.gap = cell.mass + difference(high(cell), low(cell));
  };
  std::vector<Cell> cells;
  cells.reserve(2000);
  Cell initial;
  initial.lower.fill(negative_infinity);
  initial.upper.fill(infinity);
  initial.lower_corner = corner(initial.lower);
  initial.upper_corner = corner(initial.upper);
  set_mass_gap(initial);
  cells.push_back(initial);
  const auto totals = [&]() {
    const double profile_total_start = profile_box ? profile_tick() : 0.0;
    log_lower = log_upper = negative_infinity;
    for (const Cell &cell : cells) {
      log_lower = log_add_exp(log_lower, cell.mass + low(cell));
      log_upper = log_add_exp(log_upper, cell.mass + high(cell));
    }
    if (profile_box) profile_totals += profile_tick() - profile_total_start;
  };
  const auto refine = [&]() {
    if (cells.size() >= 2000) return false;
    std::size_t selected = 0;
    for (std::size_t cell = 1; cell < cells.size(); ++cell) {
      if (cells[cell].gap > cells[selected].gap) selected = cell;
    }
    if (cells[selected].gap == negative_infinity) return false;
    const Cell parent = cells[selected];
    Cell best_left, best_right;
    double best_gap = infinity, widest = negative_infinity;
    bool found = false;
    for (int axis = 0; axis < rank; ++axis) {
      if (profile_box) {
        ++profile_interval_queries;
        ++profile_interval_keys[std::array<double, 4>{parent.lower[axis], parent.upper[axis], 0.0, 0.0}];
      }
      const double profile_interval_start = profile_box ? profile_tick() : 0.0;
      const double median = cpp_selnorm_normal_interval_quantile(
        parent.lower[axis], parent.upper[axis], 0.0, 1.0, 0.5);
      if (!std::isfinite(median) || median <= parent.lower[axis] || median >= parent.upper[axis]) continue;
      Cell left = parent, right = parent;
      left.upper[axis] = right.lower[axis] = median;
      left.axis_mass[axis] = cpp_selnorm_normal_interval_log_prob(parent.lower[axis], median, 0.0, 1.0);
      right.axis_mass[axis] = cpp_selnorm_normal_interval_log_prob(median, parent.upper[axis], 0.0, 1.0);
      if (profile_box) profile_intervals += profile_tick() - profile_interval_start;
      if (!std::isfinite(left.axis_mass[axis]) || !std::isfinite(right.axis_mass[axis])) continue;
      left.upper_corner = corner(left.upper);
      right.lower_corner = corner(right.lower);
      set_mass_gap(left);
      set_mass_gap(right);
      const double gap = log_add_exp(left.gap, right.gap);
      if (!found || gap < best_gap || (gap == best_gap && parent.axis_mass[axis] > widest)) {
        best_gap = gap; widest = parent.axis_mass[axis]; best_left = left; best_right = right; found = true;
      }
    }
    if (!found) return false;
    cells[selected] = best_left;
    cells.push_back(best_right);
    return true;
  };
  // Coarse groups of refinements avoid unbounded-corner zero-gain ties.
  // Stop preparation when its observed scalar-normalizer work plus the
  // conservative remaining proposal work ceases improving. This chooses an
  // exact finite proposal, not an integration accuracy or sample budget.
  double previous_cost = infinity;
  bool refinable = true;
  while (refinable && cells.size() < 2000) {
    for (int i = 0; i < 32 && cells.size() < 2000; ++i) {
      if (!refine()) { refinable = false; break; }
    }
    totals();
    const double cost = static_cast<double>(preparation_work) +
      n * std::exp(log_upper - log_lower);
    if (std::isfinite(previous_cost) && cost >= previous_cost) break;
    previous_cost = cost;
  }
  totals();
  const double profile_prepared = profile_box ? profile_tick() : 0.0;
  const auto profile_report = [&]() {
    if (!profile_box) return;
    static std::array<double, 7> times{};
    static std::array<std::uint64_t, 6> counts{};
    times[0] += profile_geometry; times[1] += profile_correction;
    times[2] += profile_integrands; times[3] += profile_intervals; times[4] += profile_totals;
    times[5] += profile_prepared - profile_start - profile_geometry - profile_correction + profile_late_preparation;
    times[6] += profile_sampling;
    ++counts[0]; counts[1] += profile_queries; counts[2] += profile_hits;
    counts[3] += profile_interval_queries; counts[4] += profile_interval_keys.size(); counts[5] += cells.size();
    if (counts[0] % 100 == 0) std::fprintf(stderr, "Box profile: calls %llu geometry %.6f correction %.6f integrands %.6f intervals %.6f totals %.6f preparation %.6f sampling %.6f corners %llu hits %llu interval_queries %llu unique_intervals %llu cells %llu\n", (unsigned long long) counts[0], times[0], times[1], times[2], times[3], times[4], times[5], times[6], (unsigned long long) counts[1], (unsigned long long) counts[2], (unsigned long long) counts[3], (unsigned long long) counts[4], (unsigned long long) counts[5]);
  };
  std::array<double, 4> factors{};
  std::vector<double> residual(n), transformed(n);
  std::vector<double> mass(bins), lower(bins), upper(bins);
  int completed_failures = 0;
  for (int attempt = 0; attempt < max_attempts; ++attempt) {
    // Adapt only before the next candidate, after previous rejections have
    // fully completed. An outer covariance rejection redraws both F and Y.
    if (completed_failures >= 64 && refinable && cells.size() < 2000) {
      const double profile_late_start = profile_box ? profile_tick() : 0.0;
      for (int i = 0; i < 32 && cells.size() < 2000; ++i) {
        if (!refine()) { refinable = false; break; }
      }
      totals();
      completed_failures = 0;
      if (profile_box) profile_late_preparation += profile_tick() - profile_late_start;
    }
    const double profile_sampling_start = profile_box ? profile_tick() : 0.0;
    const double threshold = std::log(uniform()) + log_upper;
    double cumulative = negative_infinity;
    std::size_t selected = cells.size();
    for (std::size_t cell = 0; cell < cells.size(); ++cell) {
      const double weight = cells[cell].mass + high(cells[cell]);
      if (!std::isfinite(weight)) continue;
      selected = cell;
      cumulative = log_add_exp(cumulative, weight);
      if (threshold <= cumulative) break;
    }
    if (selected == cells.size()) return -2;
    const Cell &cell = cells[selected];
    for (int axis = 0; axis < rank; ++axis) {
      factors[axis] = cpp_selnorm_normal_interval_quantile(
        cell.lower[axis], cell.upper[axis], 0.0, 1.0, uniform());
      if (!std::isfinite(factors[axis]) || factors[axis] < cell.lower[axis] || factors[axis] > cell.upper[axis]) return -2;
    }
    const double log_F = integrand(factors, false);
    if (log_F > high(cell) || std::isnan(log_F)) return -2;
    if (std::log(uniform()) > log_F - high(cell)) {
      if (profile_box) profile_sampling += profile_tick() - profile_sampling_start;
      ++completed_failures; continue;
    }
    for (int row = 0; row < n; ++row) {
      const double bin_uniform = uniform();
      const double interval_uniform = uniform();
      output[row] = cpp_selnorm_kernel_rng_workspace(location[row], sd[row], selection_se[row], omega,
        bin_uniform, interval_uniform, 0.0, 0, SELKERNEL_STEP, selection,
        mass.data(), lower.data(), upper.data(), 1, false);
      if (!std::isfinite(output[row])) return -2;
      residual[row] = output[row] - mean[row];
    }
    F77_CALL(dtrsv)("L", "N", "N", &n, target_cholesky.data(), &n, residual.data(), &stride FCONE FCONE FCONE);
    F77_CALL(dgemv)("T", &n, &n, &one, update.data(), &n, residual.data(), &stride,
      &zero, transformed.data(), &stride FCONE);
    F77_CALL(dtrsv)("L", "N", "N", &n, correction.data(), &n, transformed.data(), &stride FCONE FCONE FCONE);
    long double quadratic = 0.0L;
    for (double value : transformed) quadratic += static_cast<long double>(value) * value;
    if (!std::isfinite(quadratic)) return -2;
    if (std::log(uniform()) <= -0.5L * quadratic) {
      if (profile_box) profile_sampling += profile_tick() - profile_sampling_start;
      profile_report();
      return 1;
    }
    if (profile_box) profile_sampling += profile_tick() - profile_sampling_start;
    ++completed_failures;
  }
  return 0;
}

double cpp_selnorm_factor_step_lpdf(
    const double *x, const double *mean, const double *residual_sd,
    const double *loading, int dimension, int rank,
    const double *selection_se, const double *omega, int n_bins,
    const double *z_lower, const double *z_upper, const int *obs_bin,
    int effect_sign, bool telescope_probabilities, int kernel_mode,
    const double *quadrature_nodes, const double *quadrature_log_weights,
    const double *quadrature_orders, int quadrature_rule_count,
    const double *qmc, int initial_points, int max_points, int scrambles,
    double relative_tolerance, double *relative_mcse,
    double *relative_change, double *log_normalizer_out, int vector_rule)
{
  const double negative_infinity = -std::numeric_limits<double>::infinity();
  const double log_two_pi = std::log(6.283185307179586476925286766559);
  *relative_mcse = 0.0;
  *relative_change = 0.0;
  if (log_normalizer_out != nullptr) {
    *log_normalizer_out = std::numeric_limits<double>::quiet_NaN();
  }
  if (dimension < 1 || rank < 1 || rank > SELNORM_FACTOR_MAX_RANK ||
      (rank > 1 && quadrature_rule_count < 3) ||
      initial_points < 2 ||
      max_points < initial_points || scrambles < 2 ||
      !(relative_tolerance > 0.0)) {
    return negative_infinity;
  }

  double phack_z_zero[2] = {0, 0};
  double segment_bounds_zero[1] = {0};
  int segment_zero[1] = {0};
  SelNormKernelData selection;
  selection.n_bins = n_bins;
  selection.n_segments = 0;
  selection.effect_sign = effect_sign;
  selection.q = 0;
  selection.z_lower = z_lower;
  selection.z_upper = z_upper;
  selection.phack_z_source = phack_z_zero;
  selection.phack_z_dest = phack_z_zero;
  selection.segment_bounds = segment_bounds_zero;
  selection.segment_step_bin = segment_zero;
  selection.segment_phack_region = segment_zero;
  selection.segment_step_bin_real = 0;
  selection.segment_phack_region_real = 0;
  selection.trusted_step_partition = true;
  selection.telescope_probabilities = telescope_probabilities;

  std::vector<double> inverse_variance(static_cast<std::size_t>(dimension));
  std::vector<double> system(static_cast<std::size_t>(rank * rank), 0.0);
  std::vector<double> score(static_cast<std::size_t>(rank), 0.0);
  double log_det = 0.0;
  double diagonal_quadratic = 0.0;
  for (int i = 0; i < dimension; ++i) {
    if (!(residual_sd[i] > 0.0) || !std::isfinite(residual_sd[i]) ||
        !std::isfinite(x[i]) || !std::isfinite(mean[i])) {
      return negative_infinity;
    }
    const double inverse = 1.0 / (residual_sd[i] * residual_sd[i]);
    inverse_variance[i] = inverse;
    log_det += 2.0 * std::log(residual_sd[i]);
    const double residual = x[i] - mean[i];
    diagonal_quadratic += residual * residual * inverse;
    for (int factor = 0; factor < rank; ++factor) {
      const double coefficient = loading[i + dimension * factor];
      if (!std::isfinite(coefficient)) return negative_infinity;
      score[factor] += coefficient * residual * inverse;
      for (int other = 0; other <= factor; ++other) {
        system[factor + rank * other] += coefficient *
          loading[i + dimension * other] * inverse;
      }
    }
  }
  for (int factor = 0; factor < rank; ++factor) {
    system[factor + rank * factor] += 1.0;
    for (int other = 0; other < factor; ++other) {
      system[other + rank * factor] = system[factor + rank * other];
    }
  }
  std::vector<double> cholesky = system;
  int info = 0;
  F77_CALL(dpotrf)("L", &rank, cholesky.data(), &rank, &info FCONE);
  if (info != 0) return negative_infinity;
  for (int factor = 0; factor < rank; ++factor) {
    const double diagonal = cholesky[factor + rank * factor];
    if (!(diagonal > 0.0) || !std::isfinite(diagonal)) {
      return negative_infinity;
    }
    log_det += 2.0 * std::log(diagonal);
  }
  std::vector<double> solved_score = score;
  for (int row = 0; row < rank; ++row) {
    for (int column = 0; column < row; ++column) {
      solved_score[row] -= cholesky[row + rank * column] *
        solved_score[column];
    }
    solved_score[row] /= cholesky[row + rank * row];
  }
  double correction = vector_dot(solved_score, solved_score);
  double log_density = -0.5 * (
    static_cast<double>(dimension) * log_two_pi + log_det +
    diagonal_quadratic - correction
  );
  if (kernel_mode == SELKERNEL_NORMAL) {
    if (log_normalizer_out != nullptr) *log_normalizer_out = 0.0;
    return log_density;
  }
  if (vector_rule != SELVECTOR_PRODUCT) {
    const double weight = vector_observed_log_weight(x, dimension, selection_se,
      omega, z_lower, n_bins, effect_sign, vector_rule, obs_bin);
    if (!std::isfinite(weight) && log_normalizer_out == nullptr) return negative_infinity;
    log_density += weight;
  }
  for (int i = 0; i < dimension && vector_rule == SELVECTOR_PRODUCT; ++i) {
    const double log_weight = std::log(omega[obs_bin[i] - 1]);
    if (!std::isfinite(log_weight) &&
        (log_normalizer_out == nullptr || log_weight != negative_infinity)) {
      return negative_infinity;
    }
    log_density += log_weight;
  }

  FactorNormalizerContext context;
  context.vector_rule = vector_rule;
  context.mean = mean;
  context.residual_sd = residual_sd;
  context.loading = loading;
  context.dimension = dimension;
  context.rank = rank;
  context.selection_se = selection_se;
  context.omega = omega;
  context.n_bins = n_bins;
  context.kernel_mode = kernel_mode;
  context.selection = selection;

  FactorForestPlan forest;
  factor_forest_plan(context, nullptr, &forest);
  // The nested rule expands one axis per tree level, so a forest support pays
  // `order^(depth + 1)` and its ladder length does not depend on the rank.
  const int cost_exponent = forest.valid ? forest.max_depth + 1 : rank;

  // The cluster evaluator enters at rank one only after its own unchanged
  // quadrature sequence rejects. Share the QMC fallback below in that case.
  if (rank > 1) {
    int quadrature_offset = 0;
    int quadrature_order = static_cast<int>(quadrature_orders[0]);
    if (!selnorm_factor_rule_affordable(quadrature_order, cost_exponent)) {
      quadrature_order = 0;
    }
    double previous_quadrature = quadrature_order == 0 ? negative_infinity :
      factor_rule_log_integral(
        context, forest, quadrature_nodes, quadrature_log_weights,
        quadrature_offset, quadrature_order
      );
    double previous_quadrature_change =
      std::numeric_limits<double>::infinity();
    quadrature_offset += static_cast<int>(quadrature_orders[0]);
    for (int rule = 1; rule < quadrature_rule_count; ++rule) {
      quadrature_order = static_cast<int>(quadrature_orders[rule]);
      if (!selnorm_factor_rule_affordable(quadrature_order, cost_exponent)) break;
      const double current_quadrature = factor_rule_log_integral(
        context, forest, quadrature_nodes, quadrature_log_weights,
        quadrature_offset, quadrature_order
      );
      const double current_quadrature_change = cluster_relative_change(
        previous_quadrature, current_quadrature
      );
      if (std::isfinite(current_quadrature) &&
          current_quadrature_change <= relative_tolerance &&
          previous_quadrature_change <= relative_tolerance &&
          quadrature_covers_gaussian_tails(
            quadrature_nodes, quadrature_offset, quadrature_order, rank,
            omega, n_bins, dimension, vector_rule, current_quadrature,
            relative_tolerance
          )) {
        *relative_change = current_quadrature_change;
        if (log_normalizer_out != nullptr) {
          *log_normalizer_out = current_quadrature;
        }
        return log_density - current_quadrature;
      }
      previous_quadrature = current_quadrature;
      previous_quadrature_change = current_quadrature_change;
      quadrature_offset += quadrature_order;
    }
  }

  std::vector<std::vector<double> > starts;
  starts.push_back(std::vector<double>(static_cast<std::size_t>(rank), 0.0));
  std::vector<double> boundaries;
  for (int bin = 0; bin < n_bins; ++bin) {
    if (std::isfinite(z_lower[bin])) boundaries.push_back(z_lower[bin]);
    if (std::isfinite(z_upper[bin])) boundaries.push_back(z_upper[bin]);
  }
  std::sort(boundaries.begin(), boundaries.end());
  boundaries.erase(
    std::unique(boundaries.begin(), boundaries.end()), boundaries.end()
  );
  for (std::size_t boundary = 0; boundary < boundaries.size(); ++boundary) {
    std::vector<double> start(static_cast<std::size_t>(rank), 0.0);
    for (int i = 0; i < dimension; ++i) {
      const double target = effect_sign * boundaries[boundary] *
        selection_se[i] - mean[i];
      for (int factor = 0; factor < rank; ++factor) {
        start[factor] += loading[i + dimension * factor] *
          inverse_variance[i] * target;
      }
    }
    if (factor_solve_cholesky(cholesky, rank, &start)) {
      starts.push_back(start);
    }
  }

  std::vector<double> identity(static_cast<std::size_t>(rank * rank), 0.0);
  for (int factor = 0; factor < rank; ++factor) {
    identity[factor + rank * factor] = 1.0;
  }
  // Keep the prior and every finite optimized start: opposing selection modes
  // need not have the same location or scale. No mode deduplication is needed.
  std::vector<std::vector<double> > proposals(1, starts[0]);
  std::vector<std::vector<double> > proposal_factors(1, identity);
  std::vector<double> proposal_log_determinants(1, 0.0);
  for (std::size_t candidate = 0; candidate < starts.size(); ++candidate) {
    std::vector<double> covariance;
    factor_optimize_mode(context, &starts[candidate], &covariance);
    const double value = factor_log_integrand(
      context, starts[candidate].data(), nullptr
    );
    if (!std::isfinite(value)) continue;
    if (covariance.size() != identity.size() ||
        !std::all_of(covariance.begin(), covariance.end(),
                     [](double x) { return std::isfinite(x); })) {
      covariance = identity;
    }
    F77_CALL(dpotrf)("L", &rank, covariance.data(), &rank, &info FCONE);
    bool valid_factor = info == 0;
    for (int factor = 0; factor < rank && valid_factor; ++factor) {
      const double diagonal = covariance[factor + rank * factor];
      valid_factor = diagonal > 0.0 && std::isfinite(diagonal);
    }
    // An optimizer covariance only defines an importance proposal. A failed
    // approximation leaves the model covariance intact and uses unit scale.
    if (!valid_factor) covariance = identity;
    double proposal_log_determinant = 0.0;
    for (int factor = 0; factor < rank; ++factor) {
      proposal_log_determinant += std::log(
        covariance[factor + rank * factor]
      );
    }
    proposals.push_back(starts[candidate]);
    proposal_factors.push_back(covariance);
    proposal_log_determinants.push_back(proposal_log_determinant);
  }

  std::vector<double> scramble_log_mean(static_cast<std::size_t>(scrambles));
  std::vector<double> scramble_coarse_log_mean(
    static_cast<std::size_t>(scrambles)
  );
  const int proposal_count = static_cast<int>(proposals.size());
  const int qmc_dimensions = 2 * rank;
  const int max_per_proposal =
    (2 * max_points + proposal_count - 1) / proposal_count;
  const auto component_count = [proposal_count](int total, int component) {
    return total / proposal_count + (component < total % proposal_count ? 1 : 0);
  };
  std::vector<double> point_log_target(
    static_cast<std::size_t>(max_per_proposal) * proposal_count * scrambles
  );
  std::vector<double> point_latent(point_log_target.size() * rank);
  std::vector<int> fine_counts(static_cast<std::size_t>(proposal_count));
  std::vector<int> coarse_counts(static_cast<std::size_t>(proposal_count));
  std::vector<double> fine_log_counts(static_cast<std::size_t>(proposal_count));
  std::vector<double> coarse_log_counts(static_cast<std::size_t>(proposal_count));
  std::vector<double> latent(static_cast<std::size_t>(rank));
  std::vector<double> standard_normal(static_cast<std::size_t>(rank));
  std::vector<double> centered(static_cast<std::size_t>(rank));
  int evaluated_points = 0;
  int current_points = initial_points;
  int comparison_points = std::max(2, initial_points / 2);
  double log_normalizer = negative_infinity;
  while (true) {
    for (int component = 0; component < proposal_count; ++component) {
      fine_counts[component] = component_count(2 * current_points, component);
      coarse_counts[component] = component_count(
        2 * comparison_points, component
      );
      fine_log_counts[component] = std::log(
        static_cast<double>(fine_counts[component])
      );
      coarse_log_counts[component] = std::log(
        static_cast<double>(coarse_counts[component])
      );
    }
    for (int scramble = 0; scramble < scrambles; ++scramble) {
      scramble_log_mean[scramble] = negative_infinity;
      scramble_coarse_log_mean[scramble] = negative_infinity;
      for (int component = 0; component < proposal_count; ++component) {
        const int evaluated = component_count(2 * evaluated_points, component);
        const std::size_t component_offset = static_cast<std::size_t>(
          scramble * proposal_count + component
        ) * max_per_proposal;
        for (int point = evaluated; point < fine_counts[component]; ++point) {
          const int stream = proposal_count == 1 ? point % 2 : component % 2;
          const int design_point = proposal_count == 1 ? point / 2 : point;
          for (int factor = 0; factor < rank; ++factor) {
            standard_normal[factor] = qnorm(
              qmc_value(
                qmc, scrambles, max_points, qmc_dimensions, scramble,
                design_point, factor + rank * stream
              ),
              0.0, 1.0, true, false
            );
            if (!std::isfinite(standard_normal[factor])) {
              return negative_infinity;
            }
          }
          for (int factor = 0; factor < rank; ++factor) {
            latent[factor] = proposals[component][factor];
            for (int column = 0; column <= factor; ++column) {
              latent[factor] +=
                proposal_factors[component][factor + rank * column] *
                standard_normal[column];
            }
            point_latent[(component_offset + point) * rank + factor] =
              latent[factor];
          }
          const double target = factor_log_integrand(
            context, latent.data(), nullptr
          );
          if (!std::isfinite(target)) return negative_infinity;
          point_log_target[component_offset + point] = target;
        }
        for (int point = 0; point < fine_counts[component]; ++point) {
          const bool in_coarse = point < coarse_counts[component];
          double fine_mixture = negative_infinity;
          double coarse_mixture = negative_infinity;
          for (int proposal = 0; proposal < proposal_count; ++proposal) {
            if (fine_counts[proposal] == 0) continue;
            double kernel = -proposal_log_determinants[proposal];
            for (int factor = 0; factor < rank; ++factor) {
              centered[factor] =
                point_latent[(component_offset + point) * rank + factor] -
                proposals[proposal][factor];
              for (int column = 0; column < factor; ++column) {
                centered[factor] -=
                  proposal_factors[proposal][factor + rank * column] *
                  centered[column];
              }
              centered[factor] /=
                proposal_factors[proposal][factor + rank * factor];
              kernel -= 0.5 * centered[factor] * centered[factor];
            }
            fine_mixture = log_add_exp(
              fine_mixture, kernel + fine_log_counts[proposal]
            );
            if (in_coarse && coarse_counts[proposal] > 0) {
              coarse_mixture = log_add_exp(
                coarse_mixture, kernel + coarse_log_counts[proposal]
              );
            }
          }
          // For N points the mixture weights are n_j/N. Its N cancels the
          // pooled mean denominator. Coarse and fine counts define separate
          // valid mixtures, including when some proposals receive no points.
          const double target = point_log_target[component_offset + point];
          scramble_log_mean[scramble] = log_add_exp(
            scramble_log_mean[scramble], target - fine_mixture
          );
          if (in_coarse) {
            scramble_coarse_log_mean[scramble] = log_add_exp(
              scramble_coarse_log_mean[scramble], target - coarse_mixture
            );
          }
        }
      }
    }
    log_normalizer = factor_log_mean(scramble_log_mean, scrambles);
    const double coarse_log_normalizer = factor_log_mean(
      scramble_coarse_log_mean, scrambles
    );
    if (!std::isfinite(log_normalizer) ||
        !std::isfinite(coarse_log_normalizer)) return negative_infinity;
    *relative_change = cluster_relative_change(
      coarse_log_normalizer, log_normalizer
    );

    // Components reuse the two declared streams. Combine them before MCSE so
    // covariance between proposals stays within each independent scramble.
    double squared_relative = 0.0;
    for (int scramble = 0; scramble < scrambles; ++scramble) {
      const double ratio = std::exp(
        scramble_log_mean[scramble] - log_normalizer
      );
      const double difference = ratio - 1.0;
      squared_relative += difference * difference;
    }
    *relative_mcse = std::sqrt(
      squared_relative /
      (static_cast<double>(scrambles) * static_cast<double>(scrambles - 1))
    );
    if (std::max(*relative_mcse, *relative_change) <= relative_tolerance ||
        current_points == max_points) {
      break;
    }
    evaluated_points = current_points;
    comparison_points = current_points;
    current_points = std::min(max_points, 2 * current_points);
  }
  if (log_normalizer_out != nullptr) *log_normalizer_out = log_normalizer;
  return log_density - log_normalizer;
}

double cpp_selnorm_cluster_step_lpdf(
    const double *x, const double *mean, const double *residual_sd,
    const double *loading, int dimension, const double *selection_se,
    const double *omega, int n_bins, const double *z_lower,
    const double *z_upper, const int *obs_bin, int effect_sign,
    bool telescope_probabilities, int kernel_mode,
    const double *quadrature_nodes, const double *quadrature_log_weights,
    const double *quadrature_orders, int quadrature_rule_count,
    const double *qmc, int initial_points, int max_points, int scrambles,
    double relative_tolerance, double *relative_mcse, double *relative_change,
    const long double *quadrature_weights, double *log_normalizer_out,
    int vector_rule)
{
  const double negative_infinity = -std::numeric_limits<double>::infinity();
  const double log_two_pi = std::log(6.283185307179586476925286766559);
  *relative_mcse = 0.0;
  *relative_change = 0.0;
  if (log_normalizer_out != nullptr) {
    *log_normalizer_out = std::numeric_limits<double>::quiet_NaN();
  }

  double phack_z_zero[2] = {0, 0};
  double segment_bounds_zero[1] = {0};
  int segment_zero[1] = {0};
  SelNormKernelData selection;
  selection.n_bins = n_bins;
  selection.n_segments = 0;
  selection.effect_sign = effect_sign;
  selection.q = 0;
  selection.z_lower = z_lower;
  selection.z_upper = z_upper;
  selection.phack_z_source = phack_z_zero;
  selection.phack_z_dest = phack_z_zero;
  selection.segment_bounds = segment_bounds_zero;
  selection.segment_step_bin = segment_zero;
  selection.segment_phack_region = segment_zero;
  selection.segment_step_bin_real = 0;
  selection.segment_phack_region_real = 0;
  selection.trusted_step_partition = true;
  selection.telescope_probabilities = telescope_probabilities;

  double log_det = 0.0;
  double denominator = 1.0;
  double score = 0.0;
  double diagonal_quadratic = 0.0;
  for (int i = 0; i < dimension; ++i) {
    if (!(residual_sd[i] > 0.0) || !std::isfinite(residual_sd[i]) ||
        !std::isfinite(loading[i]) || !std::isfinite(x[i]) ||
        !std::isfinite(mean[i])) return negative_infinity;
    const double variance = residual_sd[i] * residual_sd[i];
    const double residual = x[i] - mean[i];
    log_det += std::log(variance);
    denominator += loading[i] * loading[i] / variance;
    score += loading[i] * residual / variance;
    diagonal_quadratic += residual * residual / variance;
  }
  const double quadratic = diagonal_quadratic - score * score / denominator;
  double log_density = -0.5 * (
    static_cast<double>(dimension) * log_two_pi + log_det +
    std::log(denominator) + quadratic
  );
  if (kernel_mode == SELKERNEL_NORMAL) {
    if (log_normalizer_out != nullptr) *log_normalizer_out = 0.0;
    return log_density;
  }

  if (vector_rule != SELVECTOR_PRODUCT) {
    const double weight = vector_observed_log_weight(x, dimension, selection_se,
      omega, z_lower, n_bins, effect_sign, vector_rule, obs_bin);
    if (!std::isfinite(weight) && log_normalizer_out == nullptr) return negative_infinity;
    log_density += weight;
  }
  for (int i = 0; i < dimension && vector_rule == SELVECTOR_PRODUCT; ++i) {
    const double log_weight = std::log(omega[obs_bin[i] - 1]);
    if (!std::isfinite(log_weight) &&
        (log_normalizer_out == nullptr || log_weight != negative_infinity)) {
      return negative_infinity;
    }
    log_density += log_weight;
  }

  ClusterNormalizerContext context;
  context.vector_rule = vector_rule;
  context.quadrature_weights = quadrature_weights;
  context.mean = mean;
  context.residual_sd = residual_sd;
  context.loading = loading;
  context.dimension = dimension;
  context.selection_se = selection_se;
  context.omega = omega;
  context.kernel_mode = kernel_mode;
  context.selection = selection;

  int offset = 0;
  int order = static_cast<int>(quadrature_orders[0]);
  double previous = cluster_rule_log_integral(
    context, quadrature_nodes, quadrature_log_weights,
    offset, order
  );
  offset += order;
  for (int rule = 1; rule < quadrature_rule_count; ++rule) {
    order = static_cast<int>(quadrature_orders[rule]);
    const double current = cluster_rule_log_integral(
      context, quadrature_nodes, quadrature_log_weights,
      offset, order
    );
    *relative_change = cluster_relative_change(previous, current);
    if (log_normalizer_out != nullptr) *log_normalizer_out = current;
    if (std::isfinite(current) && *relative_change <= relative_tolerance &&
        quadrature_covers_gaussian_tails(
          quadrature_nodes, offset, order, 1, omega, n_bins, dimension,
          vector_rule, current, relative_tolerance
        )) {
      return log_density - current;
    }
    previous = current;
    offset += order;
  }

  return cpp_selnorm_factor_step_lpdf(
    x, mean, residual_sd, loading, dimension, 1, selection_se, omega, n_bins,
    z_lower, z_upper, obs_bin, effect_sign, telescope_probabilities, kernel_mode,
    nullptr, nullptr, nullptr, 0, qmc, initial_points, max_points,
    scrambles, relative_tolerance, relative_mcse, relative_change,
    log_normalizer_out, vector_rule
  );
}

double cpp_selnorm_mnorm_step_lpdf(
    const double *x, const double *mean, const double *covariance_lower,
    int dimension, const double *selection_se, const double *omega,
    int n_bins, const double *z_lower, const double *z_upper,
    const int *obs_bin, int effect_sign, bool telescope_probabilities,
    int kernel_mode, const double *qmc, int points, int scrambles,
    double *relative_mcse, SelNormZProjection *projection,
    double *log_normalizer_out, int vector_rule,
    SelNormDenseIntegration *integration, const double *log_normalizer_override)
{
  const int k = dimension;
  const int dimensions = 2 * k;
  const double two_pi = 6.283185307179586476925286766559;
  const double negative_infinity = -std::numeric_limits<double>::infinity();
  *relative_mcse = 0.0;
  if (log_normalizer_override != nullptr &&
      (!std::isfinite(*log_normalizer_override) || projection != nullptr)) {
    throw std::invalid_argument("An explicit selection normalizer must be finite and cannot project outcomes.");
  }
  if (log_normalizer_out != nullptr) {
    *log_normalizer_out = std::numeric_limits<double>::quiet_NaN();
  }

  double phack_z_zero[2] = {0, 0};
  double segment_bounds_zero[1] = {0};
  int segment_zero[1] = {0};
  SelNormKernelData selection;
  selection.n_bins = n_bins;
  selection.n_segments = 0;
  selection.effect_sign = effect_sign;
  selection.q = 0;
  selection.z_lower = z_lower;
  selection.z_upper = z_upper;
  selection.phack_z_source = phack_z_zero;
  selection.phack_z_dest = phack_z_zero;
  selection.segment_bounds = segment_bounds_zero;
  selection.segment_step_bin = segment_zero;
  selection.segment_phack_region = segment_zero;
  selection.segment_step_bin_real = 0;
  selection.segment_phack_region_real = 0;
  selection.trusted_step_partition = true;
  selection.telescope_probabilities = telescope_probabilities;

  if (k == 1 && projection == nullptr) {
    const double variance = covariance_lower[0];
    if (!(variance > 0.0) || !std::isfinite(variance)) {
      return negative_infinity;
    }
    if (log_normalizer_override != nullptr) {
      double numerator = cpp_selnorm_kernel_lpdf(
        x[0], mean[0], std::sqrt(variance), mean[0], std::sqrt(variance),
        selection_se[0], 1.0, omega, obs_bin[0], 0.0, 0, SELKERNEL_NORMAL,
        selection, 1, false);
      if (kernel_mode != SELKERNEL_NORMAL) numerator += vector_observed_log_weight(
        x, 1, selection_se, omega, z_lower, n_bins, effect_sign, vector_rule, obs_bin);
      if (log_normalizer_out != nullptr) *log_normalizer_out = *log_normalizer_override;
      return numerator - *log_normalizer_override;
    }
    if (log_normalizer_out != nullptr) {
      *log_normalizer_out = cpp_selnorm_kernel_log_norm(
        mean[0], std::sqrt(variance), selection_se[0], omega,
        0, 0, kernel_mode, selection, 1, false
      );
    }
    return cpp_selnorm_kernel_lpdf(
      x[0], mean[0], std::sqrt(variance), mean[0], std::sqrt(variance),
      selection_se[0], 1.0, omega, obs_bin[0], 0.0, 0, kernel_mode,
      selection, 1, false
    );
  }

  std::vector<double> covariance(static_cast<std::size_t>(k * k), 0.0);
  int position = 0;
  for (int column = 0; column < k; ++column) {
    for (int row = column; row < k; ++row) {
      const double value = covariance_lower[position++];
      covariance[static_cast<std::size_t>(row + k * column)] = value;
      covariance[static_cast<std::size_t>(column + k * row)] = value;
    }
  }
  std::vector<double> cholesky = covariance;
  int info = 0;
  F77_CALL(dpotrf)("L", &k, cholesky.data(), &k, &info FCONE);
  if (info != 0) return negative_infinity;

  double log_det = 0.0;
  std::vector<double> residual(static_cast<std::size_t>(k));
  for (int i = 0; i < k; ++i) {
    const double diagonal = cholesky[static_cast<std::size_t>(i + k * i)];
    if (!(diagonal > 0) || !std::isfinite(diagonal) ||
        !std::isfinite(x[i])) {
      return negative_infinity;
    }
    log_det += 2.0 * std::log(diagonal);
    residual[static_cast<std::size_t>(i)] = x[i] - mean[i];
  }
  for (int i = 0; i < k; ++i) {
    double value = residual[static_cast<std::size_t>(i)];
    for (int j = 0; j < i; ++j) {
      value -= cholesky[static_cast<std::size_t>(i + k * j)] *
        residual[static_cast<std::size_t>(j)];
    }
    residual[static_cast<std::size_t>(i)] =
      value / cholesky[static_cast<std::size_t>(i + k * i)];
  }
  double quadratic = 0.0;
  for (int i = 0; i < k; ++i) {
    const double value = residual[static_cast<std::size_t>(i)];
    if (!std::isfinite(value)) return negative_infinity;
    quadratic += value * value;
  }
  double log_density = -0.5 *
    (static_cast<double>(k) * std::log(two_pi) + log_det + quadratic);

  if (kernel_mode == SELKERNEL_NORMAL && projection == nullptr) {
    if (log_normalizer_out != nullptr) *log_normalizer_out = 0.0;
    return log_density;
  }

  if (vector_rule != SELVECTOR_PRODUCT) {
    if (projection != nullptr) {
      *relative_mcse = std::numeric_limits<double>::infinity();
      projection->relative_error = *relative_mcse;
      return std::numeric_limits<double>::quiet_NaN();
    }
    const double weight = vector_observed_log_weight(x, k, selection_se,
      omega, z_lower, n_bins, effect_sign, vector_rule, obs_bin);
    if (!std::isfinite(weight) && log_normalizer_out == nullptr) return negative_infinity;
    if (log_normalizer_override != nullptr) {
      if (log_normalizer_out != nullptr) *log_normalizer_out = *log_normalizer_override;
      return log_density + weight - *log_normalizer_override;
    }
    const double normalizer = gaussian_event_log_mass(mean, cholesky, k,
      selection_se, omega, selection, kernel_mode, vector_rule,
      nullptr, nullptr, qmc, points, scrambles, relative_mcse);
    if (log_normalizer_out != nullptr) *log_normalizer_out = normalizer;
    return log_density + weight - normalizer;
  }

  for (int i = 0; i < k && projection == nullptr; ++i) {
    const double observed_weight = omega[obs_bin[i] - 1];
    const double log_weight = std::log(observed_weight);
    if (!std::isfinite(log_weight) &&
        (log_normalizer_out == nullptr || log_weight != negative_infinity)) {
      return negative_infinity;
    }
    log_density += log_weight;
  }

  if (log_normalizer_override != nullptr) {
    if (log_normalizer_out != nullptr) *log_normalizer_out = *log_normalizer_override;
    return log_density - *log_normalizer_override;
  }

  bool diagonal_covariance = true;
  for (int column = 0; column < k && diagonal_covariance; ++column) {
    for (int row = column + 1; row < k; ++row) {
      if (covariance[static_cast<std::size_t>(row + k * column)] != 0) {
        diagonal_covariance = false;
        break;
      }
    }
  }
  if (diagonal_covariance || kernel_mode == SELKERNEL_NORMAL) {
    double log_normalizer = 0.0;
    for (int i = 0; i < k; ++i) {
      log_normalizer += cpp_selnorm_kernel_log_norm(
        mean[i], std::sqrt(covariance[static_cast<std::size_t>(i + k * i)]),
        selection_se[i], omega, 0, 0, kernel_mode, selection, 1, false
      );
    }
    if (log_normalizer_out != nullptr) *log_normalizer_out = log_normalizer;
    if (projection != nullptr) {
      for (int z = 0; z < projection->size; ++z) {
        double value = 0.0;
        for (int i = 0; i < k; ++i) {
          const double sd = std::sqrt(covariance[i + k * i]);
          if (projection->probability) {
            double inverse;
            value += cpp_selnorm_kernel_threshold(projection->z[z], mean[i],
              sd, selection_se[i], omega, 0, 0, kernel_mode, selection,
              &inverse, 1, false);
          } else {
            int bin = 1;
            for (int b = 0; b < n_bins; ++b) {
              const double signed_z = effect_sign * projection->z[z];
              if (signed_z >= z_lower[b] && signed_z <= z_upper[b]) {
                bin = b + 1;
                break;
              }
            }
            value += selection_se[i] * std::exp(cpp_selnorm_kernel_lpdf(
              projection->z[z] * selection_se[i], mean[i], sd, mean[i], sd,
              selection_se[i], 1.0, omega, bin, 0.0, 0, kernel_mode,
              selection, 1, false
            ));
          }
        }
        projection->density[z] = value / k;
      }
      projection->relative_error = 0.0;
      return -log_normalizer;
    }
    return log_density - log_normalizer;
  }

  if (projection == nullptr && integration != nullptr && kernel_mode == SELKERNEL_STEP) {
    double log_normalizer = 0.0;
    if (exact_envelope_log_integral(covariance, mean, selection_se, omega, k,
                                    selection, integration, &log_normalizer)) {
      if (log_normalizer_out != nullptr) *log_normalizer_out = log_normalizer;
      return log_density - log_normalizer;
    }
  }

  std::vector<double> precision;
  std::vector<double> conditional_sd;
  std::vector<double> conditional_mean(static_cast<std::size_t>(k));
  std::vector<double> conditional_log_norm(static_cast<std::size_t>(k));
  std::vector<double> projection_log_sum;
  std::vector<double> projection_scale(static_cast<std::size_t>(scrambles), negative_infinity);
  const int factor_rank = projection == nullptr ? 0 : projection->rank;
  FactorNormalizerContext factor_context;
  std::vector<double> factor_mode(static_cast<std::size_t>(factor_rank), 0.0);
  std::vector<double> factor_latent(static_cast<std::size_t>(factor_rank));
  std::vector<double> factor_normal(static_cast<std::size_t>(factor_rank));
  std::vector<double> factor_centered(static_cast<std::size_t>(factor_rank));
  std::vector<double> proposal_factor(static_cast<std::size_t>(factor_rank * factor_rank), 0.0);
  double proposal_log_det = 0.0;
  const bool quadrature = projection != nullptr && projection->order > 0;
  if (factor_rank) {
    factor_context.mean = mean;
    factor_context.residual_sd = projection->residual_sd;
    factor_context.loading = projection->factor_loading;
    factor_context.dimension = k;
    factor_context.rank = factor_rank;
    factor_context.selection_se = selection_se;
    factor_context.omega = omega;
    factor_context.n_bins = n_bins;
    factor_context.kernel_mode = kernel_mode;
    factor_context.selection = selection;
  }
  if (quadrature && factor_rank > 0 && projection->loading_support != nullptr) {
    double log_normalizer = 0.0;
    FactorForestPlan forest;
    factor_forest_plan(factor_context, projection->loading_support, &forest);
    if (factor_nested_rule_log_integral(
          factor_context, forest, projection->nodes, projection->log_weights,
          0, projection->order, &log_normalizer, projection
        )) {
      *relative_mcse = 0.0;
      return -log_normalizer;
    }
  }
  if (factor_rank && !quadrature) {
    for (int i = 0; i < factor_rank; ++i) proposal_factor[i + factor_rank * i] = 1.0;
    factor_optimize_mode(factor_context, &factor_mode, &proposal_factor);
    F77_CALL(dpotrf)("L", &factor_rank, proposal_factor.data(), &factor_rank, &info FCONE);
    if (info != 0) {
      // Only the importance proposal changes here; the model covariance is
      // untouched. The prior component always retains full Gaussian support.
      std::fill(proposal_factor.begin(), proposal_factor.end(), 0.0);
      for (int i = 0; i < factor_rank; ++i) proposal_factor[i + factor_rank * i] = 1.0;
    }
    for (int i = 0; i < factor_rank; ++i) {
      proposal_log_det += std::log(proposal_factor[i + factor_rank * i]);
    }
  }
  if (projection != nullptr) {
    precision = cholesky;
    F77_CALL(dpotri)("L", &k, precision.data(), &k, &info FCONE);
    if (info != 0) return negative_infinity;
    conditional_sd.resize(k);
    for (int column = 0; column < k; ++column) {
      conditional_sd[column] = factor_rank ? projection->residual_sd[column] :
        1.0 / std::sqrt(precision[column + k * column]);
      for (int row = column + 1; row < k; ++row) {
        precision[column + k * row] = precision[row + k * column];
      }
    }
    projection_log_sum.assign(
      static_cast<std::size_t>(scrambles * projection->size), 0.0
    );
  }

  std::vector<double> scramble_log_mean(static_cast<std::size_t>(scrambles));
  std::vector<double> latent(static_cast<std::size_t>(k));
  std::vector<double> mass(static_cast<std::size_t>(n_bins));
  std::vector<double> lower(static_cast<std::size_t>(n_bins));
  std::vector<double> upper(static_cast<std::size_t>(n_bins));
  // The first conditional has no preceding latent coordinates. Prepare its
  // unchanged partition once, with the existing per-particle fallback intact.
  std::vector<double> first_mass, first_lower, first_upper;
  int first_groups = 0;
  double first_normalizer = 0.0;
  double first_log_normalizer = 0.0;
  bool first_prepared = false;
  if (projection == nullptr && kernel_mode == SELKERNEL_STEP &&
      vector_rule == SELVECTOR_PRODUCT && selection.telescope_probabilities) {
    first_mass.resize(n_bins);
    first_lower.resize(n_bins);
    first_upper.resize(n_bins);
    first_prepared = cpp_selnorm_step_rng_prepare_partition(
      mean[0], cholesky[0], selection_se[0], omega, selection,
      first_mass.data(), first_lower.data(), first_upper.data(),
      &first_groups, &first_normalizer, &first_log_normalizer
    );
  }
  for (int scramble = 0; scramble < scrambles; ++scramble) {
    // Deterministic quadrature has identical nodes in every scramble.
    if (quadrature && scramble > 0) {
      scramble_log_mean[scramble] = scramble_log_mean[0];
      for (int z = 0; z < projection->size; ++z) {
        projection_log_sum[scramble + scrambles * z] =
          projection_log_sum[scrambles * z];
      }
      continue;
    }
    double log_sum = negative_infinity;
    for (int point = 0; point < points; ++point) {
      for (int component = 0; component < (factor_rank && !quadrature ? 2 : 1); ++component) {
        double log_particle = 0.0;
        if (factor_rank) {
          double log_prior = 0.0;
          double log_mode = 0.0;
          int remaining = point;
          for (int factor = 0; factor < factor_rank; ++factor) {
            const int node = quadrature ? remaining % projection->order : 0;
            if (quadrature) remaining /= projection->order;
            const double normal = quadrature ? projection->nodes[node] :
              qnorm(qmc_value(qmc, scrambles, points, dimensions, scramble, point,
                factor + factor_rank * component), 0, 1, true, false);
            if (quadrature) log_particle += projection->log_weights[node];
            factor_normal[factor] = normal;
            factor_latent[factor] = normal;
            if (component == 1) {
              factor_latent[factor] = factor_mode[factor];
              for (int j = 0; j <= factor; ++j) {
                factor_latent[factor] += proposal_factor[factor + factor_rank * j] * factor_normal[j];
              }
            }
            log_prior -= .5 * factor_latent[factor] * factor_latent[factor];
            if (!quadrature) {
              factor_centered[factor] = factor_latent[factor] - factor_mode[factor];
              for (int j = 0; j < factor; ++j) {
                factor_centered[factor] -= proposal_factor[factor + factor_rank * j] * factor_centered[j];
              }
              factor_centered[factor] /= proposal_factor[factor + factor_rank * factor];
              log_mode -= .5 * factor_centered[factor] * factor_centered[factor];
            }
          }
          // Two equally allocated Gaussian proposals, including their mixture
          // importance correction. The outer reduction divides by 'points'.
          log_particle = quadrature ? log_particle + std::log(points) :
            log_prior - log_add_exp(log_prior, log_mode - proposal_log_det);
          for (int i = 0; i < k; ++i) {
            conditional_mean[i] = mean[i];
            for (int factor = 0; factor < factor_rank; ++factor) {
              conditional_mean[i] += projection->factor_loading[i + k * factor] *
                factor_latent[factor];
            }
            conditional_log_norm[i] = cpp_selnorm_step_log_norm(
              conditional_mean[i], conditional_sd[i], selection_se[i], omega,
              selection, 1, false);
            log_particle += conditional_log_norm[i];
          }
        }
        for (int i = 0; i < k && !factor_rank; ++i) {
          double conditional_mean = mean[i];
          for (int j = 0; j < i; ++j) {
            conditional_mean +=
              cholesky[static_cast<std::size_t>(i + k * j)] *
              latent[static_cast<std::size_t>(j)];
          }
          const double conditional_sd =
            cholesky[static_cast<std::size_t>(i + k * i)];
          // Only the final normalizer contributes to the path weight; no later
          // coordinate consumes a draw from this conditional.
          if (i == k - 1 && projection == nullptr) {
            const double log_local = cpp_selnorm_step_log_norm(
              conditional_mean, conditional_sd, selection_se[i], omega,
              selection, 1, false
            );
            if (!std::isfinite(log_local)) return negative_infinity;
            log_particle += log_local;
            break;
          }
          double log_local = 0.0;
          double sampled = 0.0;
          const double u_bin = qmc_value(
            qmc, scrambles, points, dimensions, scramble, point, 2 * i
          );
          const double u_interval = qmc_value(
            qmc, scrambles, points, dimensions, scramble, point, 2 * i + 1
          );
          if (selection.telescope_probabilities) {
            if (i == 0 && first_prepared) {
              log_local = first_log_normalizer;
              sampled = cpp_selnorm_step_rng_from_partition(
                conditional_mean, conditional_sd, effect_sign, u_bin, u_interval,
                first_mass.data(), first_lower.data(), first_upper.data(),
                first_groups, first_normalizer
              );
            } else {
              sampled = cpp_selnorm_step_log_norm_rng_workspace(
                conditional_mean, conditional_sd, selection_se[i], omega,
                u_bin, u_interval, selection, mass.data(), lower.data(),
                upper.data(), &log_local, 1, false
              );
            }
          } else {
            log_local = cpp_selnorm_kernel_log_norm(
              conditional_mean, conditional_sd, selection_se[i], omega,
              0, 0, kernel_mode, selection, 1, false
            );
            sampled = cpp_selnorm_kernel_rng_workspace(
              conditional_mean, conditional_sd, selection_se[i], omega,
              u_bin, u_interval, 0, 0, kernel_mode, selection,
              mass.data(), lower.data(), upper.data(), 1, false
            );
          }
          if (!std::isfinite(log_local) || !std::isfinite(sampled)) {
            return negative_infinity;
          }
          log_particle += log_local;
          latent[static_cast<std::size_t>(i)] =
            (sampled - conditional_mean) / conditional_sd;
          if (projection != nullptr) {
            // Preserve the full proposal vector for all-coordinate conditionals.
            // 'residual' is no longer needed by the Gaussian likelihood here.
            residual[i] = sampled - mean[i];
          }
        }
        log_sum = log_add_exp(log_sum, log_particle);
        if (projection != nullptr) {
          if (log_particle > projection_scale[scramble]) {
            const double scale = std::exp(projection_scale[scramble] - log_particle);
            for (int z = 0; z < projection->size; ++z) {
              projection_log_sum[scramble + scrambles * z] *= scale;
            }
            projection_scale[scramble] = log_particle;
          }
          const double particle_weight = std::exp(log_particle - projection_scale[scramble]);
          for (int i = 0; i < k && !factor_rank; ++i) {
            double score = 0.0;
            for (int j = 0; j < k; ++j) {
              if (j != i) score += precision[i + k * j] * residual[j];
            }
            conditional_mean[i] = mean[i] - score / precision[i + k * i];
            conditional_log_norm[i] = cpp_selnorm_step_log_norm(
              conditional_mean[i], conditional_sd[i], selection_se[i], omega,
              selection, 1, false
            );
          }
          for (int z = 0; z < projection->size; ++z) {
            double local_sum = 0.0;
            int bin = 1;
            for (int b = 0; b < n_bins; ++b) {
              const double signed_z = effect_sign * projection->z[z];
              if (signed_z >= z_lower[b] && signed_z <= z_upper[b]) {
                bin = b + 1;
                break;
              }
            }
            if (!projection->probability && !(omega[bin - 1] > 0.0)) continue;
            const double log_weight = std::log(omega[bin - 1]);
            for (int i = 0; i < k; ++i) {
              double local;
              if (projection->probability) {
                double inverse;
                local = std::log(cpp_selnorm_kernel_threshold(
                  projection->z[z], conditional_mean[i], conditional_sd[i],
                  selection_se[i], omega, 0, 0, kernel_mode, selection,
                  &inverse, 1, false));
              } else {
                local = std::log(selection_se[i]) + log_weight +
                  dnorm(projection->z[z] * selection_se[i], conditional_mean[i],
                        conditional_sd[i], true) - conditional_log_norm[i];
              }
              local_sum += std::exp(local);
            }
            const int index = scramble + scrambles * z;
            projection_log_sum[index] += particle_weight * local_sum / k;
          }
        }
      }
    }
    scramble_log_mean[static_cast<std::size_t>(scramble)] =
      log_sum - std::log(static_cast<double>(points));
    if (projection != nullptr) {
      for (int z = 0; z < projection->size; ++z) {
        const int index = scramble + scrambles * z;
        projection_log_sum[index] = std::log(projection_log_sum[index]) +
          projection_scale[scramble];
      }
    }
  }

  double log_normalizer = negative_infinity;
  for (int scramble = 0; scramble < scrambles; ++scramble) {
    log_normalizer = log_add_exp(
      log_normalizer,
      scramble_log_mean[static_cast<std::size_t>(scramble)]
    );
  }
  log_normalizer -= std::log(static_cast<double>(scrambles));
  if (log_normalizer_out != nullptr) *log_normalizer_out = log_normalizer;
  if (!std::isfinite(log_normalizer)) return negative_infinity;

  double squared_relative = 0.0;
  for (int scramble = 0; scramble < scrambles; ++scramble) {
    const double ratio = std::exp(
      scramble_log_mean[static_cast<std::size_t>(scramble)] - log_normalizer
    );
    const double difference = ratio - 1.0;
    squared_relative += difference * difference;
  }
  *relative_mcse = std::sqrt(
    squared_relative /
    (static_cast<double>(scrambles) * static_cast<double>(scrambles - 1))
  );

  if (projection != nullptr) {
    double peak = 0.0;
    double max_mcse = 0.0;
    for (int z = 0; z < projection->size; ++z) {
      double numerator = negative_infinity;
      for (int scramble = 0; scramble < scrambles; ++scramble) {
        numerator = log_add_exp(numerator, projection_log_sum[scramble + scrambles * z]);
      }
      const double value = std::exp(numerator - std::log(points * scrambles) - log_normalizer);
      projection->density[z] = value;
      peak = std::max(peak, value);
      double squared = 0.0;
      for (int scramble = 0; scramble < scrambles; ++scramble) {
        const double centered = std::exp(
          projection_log_sum[scramble + scrambles * z] - std::log(points) - log_normalizer
        ) - value * std::exp(scramble_log_mean[scramble] - log_normalizer);
        squared += centered * centered;
      }
      max_mcse = std::max(max_mcse, std::sqrt(squared / (scrambles * (scrambles - 1.0))));
    }
    projection->relative_error = peak > 0.0 ? max_mcse / peak : 0.0;
    return -log_normalizer;
  }

  return log_density - log_normalizer;
}

double cpp_selnorm_mnorm_step_surrogate_lpdf(
    const double *x, const double *mean, const double *covariance_lower,
    int dimension, const double *selection_se, const double *omega, int n_bins,
    const double *z_lower, const double *z_upper, const int *obs_bin,
    int effect_sign, bool telescope_probabilities, int kernel_mode, int vector_rule,
    const SelNormDenseIntegration &quadrature, const SelNormCoarseSettings &settings)
{
  const double unit_log_normalizer = 0.0;
  double ignored_mcse;
  // The existing body still evaluates the actual Gaussian density and actual
  // observed weights, validating actual covariance before anchor integration.
  const double numerator = cpp_selnorm_mnorm_step_lpdf(x, mean, covariance_lower, dimension,
    selection_se, omega, n_bins, z_lower, z_upper, obs_bin, effect_sign,
    telescope_probabilities, kernel_mode, nullptr, 1, 2, &ignored_mcse,
    nullptr, nullptr, vector_rule, nullptr, &unit_log_normalizer);
  if (!std::isfinite(numerator)) return numerator;
  const double log_normalizer = coarse_anchor_log_normalizer(mean, covariance_lower,
    dimension, selection_se, omega, n_bins, z_lower, z_upper, effect_sign,
    telescope_probabilities, kernel_mode, vector_rule, quadrature, settings);
  return numerator - log_normalizer;
}

double cpp_selnorm_gaussian_event_log_mass(
    const double *mean, const double *covariance_lower, int dimension,
    const double *selection_se, const double *omega, int n_bins,
    const double *z_lower, const double *z_upper, int effect_sign,
    int kernel_mode, int vector_rule, const double *lower,
    const double *upper, const double *qmc, int points, int scrambles,
    double *relative_mcse, const double *rank_one_loading,
    SelNormDenseIntegration *integration)
{
  *relative_mcse = std::numeric_limits<double>::infinity();
  const double invalid = std::numeric_limits<double>::quiet_NaN();
  if (dimension < 1 || points < 1 || scrambles < 2 ||
      vector_rule < SELVECTOR_PRODUCT || vector_rule > SELVECTOR_BEST_TWO_SIDED ||
      (kernel_mode != SELKERNEL_NORMAL && kernel_mode != SELKERNEL_STEP)) return invalid;
  SelNormKernelData selection = {};
  selection.n_bins = n_bins;
  selection.effect_sign = effect_sign;
  selection.z_lower = z_lower;
  selection.z_upper = z_upper;
  selection.trusted_step_partition = true;
  if (rank_one_loading != nullptr) {
    *relative_mcse = 0.0;
    return rank_one_event_log_mass(mean, rank_one_loading, dimension,
      selection_se, omega, selection, kernel_mode, vector_rule, lower, upper);
  }
  std::vector<double> cholesky(static_cast<std::size_t>(dimension * dimension));
  int position = 0;
  bool deterministic = true;
  bool diagonal = true;
  for (int column = 0; column < dimension; ++column) {
    for (int row = column; row < dimension; ++row) {
      const double value = covariance_lower[position++];
      if (!std::isfinite(value)) return invalid;
      deterministic = deterministic && value == 0.0;
      diagonal = diagonal && (row == column || value == 0.0);
      cholesky[row + dimension * column] = value;
      cholesky[column + dimension * row] = value;
    }
  }
  if (deterministic) {
    std::vector<double> loading(dimension, 0.0);
    *relative_mcse = 0.0;
    return rank_one_event_log_mass(mean, loading.data(), dimension,
      selection_se, omega, selection, kernel_mode, vector_rule, lower, upper);
  }
  bool full_space = true;
  for (int row = 0; row < dimension; ++row) {
    full_space = full_space &&
      (lower == nullptr || lower[row] == -std::numeric_limits<double>::infinity()) &&
      (upper == nullptr || upper[row] == std::numeric_limits<double>::infinity());
  }
  bool constant = true;
  for (int bin = 1; bin < n_bins; ++bin) constant = constant && omega[bin] == omega[0];
  const bool try_envelope = integration != nullptr && full_space && !diagonal &&
    !constant && kernel_mode == SELKERNEL_STEP && vector_rule == SELVECTOR_PRODUCT &&
    selnorm_is_descending_step_partition(z_lower, z_upper, n_bins);
  // Keep the supplied covariance before Cholesky overwrites it. The event
  // normalizer uses the same envelope and error gate as the dense likelihood.
  std::vector<double> covariance;
  if (try_envelope) covariance = cholesky;
  int info = 0;
  F77_CALL(dpotrf)("L", &dimension, cholesky.data(), &dimension, &info FCONE);
  if (info != 0) return invalid;
  if (try_envelope) {
    SelNormKernelData normalizer_selection = selection;
    normalizer_selection.telescope_probabilities = true;
    double log_mass;
    if (exact_envelope_log_integral(covariance, mean, selection_se, omega,
        dimension, normalizer_selection, integration, &log_mass)) {
      *relative_mcse = 0.0;
      return log_mass;
    }
  }
  return gaussian_event_log_mass(mean, cholesky, dimension, selection_se,
    omega, selection, kernel_mode, vector_rule, lower, upper, qmc,
    points, scrambles, relative_mcse);
}

#include "selnorm-sampling-conditioned.cc.inc"

std::size_t cpp_selnorm_cache_snapshot_size()
{
  return envelope_memo().snapshot_size();
}

bool cpp_selnorm_cache_snapshot_write(unsigned char *data, std::size_t size)
{
  return envelope_memo().snapshot_write(data, size);
}

SelNormCacheRestoreInfo cpp_selnorm_cache_restore(const SelNormCacheBlob *snapshots, std::size_t count)
{
  return envelope_memo().restore(snapshots, count);
}

SelNormCacheInfo cpp_selnorm_cache_control(
    bool set_capacity, std::size_t capacity_bytes, unsigned int clear_mask)
{
  return envelope_memo().control(set_capacity, capacity_bytes, clear_mask);
}


bool cpp_selnorm_covariance_envelope_components(
    const double *covariance_lower, int dimension, const double *selection_se,
    bool upper, double *residual_sd, double *loading, int *groups, int *type_index)
{
  if (dimension < 2 || covariance_lower == nullptr || selection_se == nullptr ||
      residual_sd == nullptr || loading == nullptr || groups == nullptr ||
      type_index == nullptr) return false;
  std::vector<double> covariance(static_cast<std::size_t>(dimension) * dimension);
  std::size_t position = 0;
  for (int column = 0; column < dimension; ++column) {
    for (int row = column; row < dimension; ++row) {
      const double value = covariance_lower[position++];
      if (!std::isfinite(value)) return false;
      covariance[row + dimension * column] = value;
      covariance[column + dimension * row] = value;
    }
  }
  // The existing envelope helper assumes the actual target passed Cholesky.
  // Validate an owned copy; never overwrite or repair the actual covariance.
  std::vector<double> cholesky = covariance;
  int info = 0;
  F77_CALL(dpotrf)("L", &dimension, cholesky.data(), &dimension, &info FCONE);
  if (info != 0) return false;
  SelNormCovarianceEnvelope envelope;
  if (!covariance_envelope(covariance, dimension, selection_se, &envelope)) {
    return false;
  }
  std::vector<double> mean(dimension, 0.0);
  std::vector<EnvelopeGroupGeometry> geometry;
  if (!prepare_envelope_geometry(covariance, mean.data(), selection_se,
                                 dimension, envelope, upper, &geometry)) return false;
  std::fill_n(loading, static_cast<std::size_t>(dimension) * 4, 0.0);
  std::vector<int> local_index(envelope.groups, 0);
  for (int row = 0; row < dimension; ++row) {
    const int group = envelope.type_index[row];
    const int local = local_index[group]++;
    const EnvelopeGroupGeometry &part = geometry[group];
    residual_sd[row] = part.residual_sd[local];
    const long double variance = static_cast<long double>(residual_sd[row]) * residual_sd[row];
    if (!(variance > 0.0L) || !std::isfinite(variance) ||
        !(static_cast<double>(variance) > 0.0)) return false;
    loading[row] = part.common_loading[local];
    if (envelope.groups > 1) {
      loading[row + dimension * (group + 1)] = part.loading[local];
    }
    type_index[row] = group;
  }
  *groups = envelope.groups;
  return true;
}



bool cpp_selnorm_context_star_projection(
  const double *means, int contexts, const double *covariance_lower, int dimension,
  const double *selection_se, const double *omega, int n_bins,
  const double *lower, const double *upper, int effect_sign, bool telescope,
  const double *log_normalizers, const double *context_log_weights,
  const double *nodes, const double *log_weights, int order,
  const double *z, int size, double allowance, double *density,
  double *mass, double *compression_error, std::size_t *components,
  double *compact_log_normalizers,
  bool bounded_cdf, double *cdf_factor_log_errors, double *input_omitted_mass_error,
  int threads)
{
  if (input_omitted_mass_error != nullptr) *input_omitted_mass_error = 0.0;
  std::vector<double> covariance(static_cast<std::size_t>(dimension) * dimension);
  std::size_t position = 0;
  for (int column = 0; column < dimension; ++column) {
    for (int row = column; row < dimension; ++row) {
      const double value = covariance_lower[position++];
      if (!std::isfinite(value)) return false;
      covariance[row + dimension * column] = value;
      covariance[column + dimension * row] = value;
    }
  }
  std::vector<double> cholesky = covariance;
  int info = 0;
  F77_CALL(dpotrf)("L", &dimension, cholesky.data(), &dimension, &info FCONE);
  if (info != 0) return false;
  SelNormCovarianceEnvelope envelope;
  if (!covariance_envelope(covariance, dimension, selection_se, &envelope)) return false;
  SelNormKernelData selection = {};
  selection.n_bins = n_bins;
  selection.effect_sign = effect_sign;
  selection.z_lower = lower;
  selection.z_upper = upper;
  selection.trusted_step_partition = true;
  selection.telescope_probabilities = telescope;
  const double cdf_factor_error = bounded_cdf ?
    bounded_cdf_product_log_error(selection, omega, 1,
      std::numeric_limits<double>::infinity()) : 0.0;
  selection.bounded_cdf = cdf_factor_error > 0.0;
  if (cdf_factor_log_errors != nullptr) {
    std::fill_n(cdf_factor_log_errors, contexts, cdf_factor_error);
  }
  std::size_t active_contexts = 0;
  for (int context = 0; context < contexts; ++context) {
    if (context_log_weights[context] != -std::numeric_limits<double>::infinity()) ++active_contexts;
  }
  EnvelopeRuleProjection projection;
  if (active_contexts > 0) {
    projection.means.resize(dimension);
    projection.log_coefficient.resize(dimension);
    std::vector<int> group_sizes(envelope.groups, 0);
    for (int row = 0; row < dimension; ++row) ++group_sizes[envelope.type_index[row]];
    for (int row = 0; row < dimension; ++row) {
      const bool child = envelope.groups > 1 && group_sizes[envelope.type_index[row]] > 1;
      const std::size_t count = active_contexts * order * (child ? order : 1);
      projection.means[row].reserve(count);
      projection.log_coefficient[row].reserve(count);
    }
  }
  std::vector<long double> weights(order);
  for (int node = 0; node < order; ++node) weights[node] = std::exp(static_cast<long double>(log_weights[node]));
  std::vector<double> residual_sd(active_contexts > 0 ? dimension : 0);
  double log_mass = -std::numeric_limits<double>::infinity();

  // Contexts are independent quadrature nodes of the retained axis. Each one
  // builds its own mixture fragment; the fragments are concatenated and the
  // masses reduced below in context order, so a threaded evaluation produces
  // the serial result element for element.
  const int context_workers = threads > 1 && contexts > 1 ?
    std::min(threads, contexts) : 1;
  std::vector<double> context_local(contexts,
    std::numeric_limits<double>::quiet_NaN());
  std::vector<EnvelopeRuleProjection> fragments(
    context_workers > 1 ? static_cast<std::size_t>(contexts) : 0);
  for (EnvelopeRuleProjection &fragment : fragments) {
    fragment.means.resize(dimension);
    fragment.log_coefficient.resize(dimension);
  }
  std::atomic<bool> context_failed{false};
  const auto evaluate_context = [&](int context, EnvelopeRuleWorkspace &workspace,
                                    EnvelopeRuleProjection *target) {
    std::vector<double> mean(dimension);
    for (int row = 0; row < dimension; ++row) mean[row] = means[context + contexts * row];
    std::vector<EnvelopeGroupGeometry> geometry;
    if (!prepare_envelope_geometry(covariance, mean.data(), selection_se, dimension,
        envelope, false, &geometry)) return false;
    if (active_contexts > 0) {
      for (const EnvelopeGroupGeometry &group : geometry) {
        for (std::size_t row = 0; row < group.row_index.size(); ++row) {
          residual_sd[group.row_index[row]] = group.residual_sd[row];
        }
      }
    }
    const bool collect = context_log_weights[context] != -std::numeric_limits<double>::infinity();
    const double log_scale = context_log_weights[context] -
      (log_normalizers == nullptr ? 0.0 : log_normalizers[context]);
    std::vector<std::size_t> starts(
      log_normalizers == nullptr && collect ? dimension : 0);
    if (log_normalizers == nullptr && collect) {
      for (int row = 0; row < dimension; ++row) starts[row] = target->log_coefficient[row].size();
    }
    target->log_scale = log_scale - std::log(static_cast<double>(dimension));
    const double local = envelope_rule_log_integral(geometry, omega, selection,
      nodes, log_weights, weights.data(), 0, order, &workspace, collect ? target : nullptr);
    if (!std::isfinite(local)) return false;
    context_local[context] = local;
    if (log_normalizers == nullptr && collect) {
      // Explicit private own-A0 mode. The caller must validate raw cross-rule
      // A0 changes/tails; this same-rule mass is algebraically normalized and
      // is not a diagnostic. The actual full-C law is covered separately.
      for (int row = 0; row < dimension; ++row) {
        for (std::size_t index = starts[row]; index < target->log_coefficient[row].size(); ++index) {
          target->log_coefficient[row][index] -= local;
        }
      }
    }
    return true;
  };

#if defined(_OPENMP)
  if (context_workers > 1) {
    #pragma omp parallel num_threads(context_workers)
    {
      EnvelopeRuleWorkspace workspace;
      #pragma omp for schedule(dynamic, 1)
      for (int context = 0; context < contexts; ++context) {
        if (context_failed.load(std::memory_order_relaxed)) continue;
        bool ok = false;
        try {
          ok = evaluate_context(context, workspace, &fragments[context]);
        } catch (...) {
          ok = false;
        }
        if (!ok) context_failed.store(true, std::memory_order_relaxed);
      }
    }
    if (context_failed.load(std::memory_order_relaxed)) return false;
    // Nothing was collected when no context carries weight, and the shared
    // projection's row vectors are then not even sized.
    for (int context = 0; active_contexts > 0 && context < contexts; ++context) {
      EnvelopeRuleProjection &fragment = fragments[context];
      for (int row = 0; row < dimension; ++row) {
        projection.means[row].insert(projection.means[row].end(),
          fragment.means[row].begin(), fragment.means[row].end());
        projection.log_coefficient[row].insert(projection.log_coefficient[row].end(),
          fragment.log_coefficient[row].begin(), fragment.log_coefficient[row].end());
      }
    }
  } else
#endif
  {
    EnvelopeRuleWorkspace workspace;
    for (int context = 0; context < contexts; ++context) {
      if (!evaluate_context(context, workspace, &projection)) return false;
    }
  }

  for (int context = 0; context < contexts; ++context) {
    const double local = context_local[context];
    if (compact_log_normalizers != nullptr) compact_log_normalizers[context] = local;
    const bool collect = context_log_weights[context] != -std::numeric_limits<double>::infinity();
    if (log_normalizers == nullptr && collect) {
      log_mass = log_add_exp(log_mass, context_log_weights[context]);
    } else if (log_normalizers != nullptr) {
      log_mass = log_add_exp(log_mass, context_log_weights[context] -
        log_normalizers[context] + local);
    }
  }
  if (active_contexts == 0) {
    // All raw anchors above were evaluated normally. Preserve the existing
    // zero-mass/unavailable result without building an unused mixture plan.
    std::fill_n(density, size, 0.0);
    *mass = 0.0;
    *compression_error = 0.0;
    *components = 0;
    return false;
  }
  std::vector<double> grid_weight(size);
  for (int point = 0; point < size; ++point) {
    int bin = 0;
    for (int candidate = 0; candidate < n_bins; ++candidate) {
      const double value = effect_sign * z[point];
      if (value >= lower[candidate] && value <= upper[candidate]) {
        bin = candidate;
        break;
      }
    }
    grid_weight[point] = omega[bin];
  }
  SelNormMixtureProjectionWorkspace mixture(z, size, grid_weight.data());
  std::vector<SelNormGaussianMixtureRow> rows;
  *components = 0;
  for (int row = 0; row < dimension; ++row) {
    const std::size_t count = projection.means[row].size();
    *components += count;
    rows.push_back({projection.means[row].data(), projection.log_coefficient[row].data(),
      count, residual_sd[row], selection_se[row]});
  }
  std::vector<long double> values;
  long double error = 0.0L, input_omitted_mass = 0.0L;
  selnorm_evaluate_gaussian_mixtures(rows, allowance, mixture, values, error,
    input_omitted_mass, threads);
  for (int point = 0; point < size; ++point) {
    density[point] = static_cast<double>(values[point]);
    if (!std::isfinite(density[point]) || density[point] < 0.0) return false;
  }
  // Mass uses the largest selection weight over ALL bins, not only grid
  // weights. The z-coordinate SE Jacobian cancels from the integrated loss.
  if (input_omitted_mass > 0) {
    const long double weighted = input_omitted_mass * *std::max_element(omega, omega + n_bins);
    if (!(weighted > 0) || !std::isfinite(weighted)) return false;
    input_omitted_mass = selnorm_mixture_compression_detail::round_up(weighted);
  }
  double represented_omission = static_cast<double>(input_omitted_mass);
  if (static_cast<long double>(represented_omission) < input_omitted_mass)
    represented_omission = std::nextafter(represented_omission, std::numeric_limits<double>::infinity());
  if (!std::isfinite(represented_omission)) return false;
  if (input_omitted_mass_error != nullptr) *input_omitted_mass_error = represented_omission;
  *mass = std::exp(log_mass);
  *compression_error = static_cast<double>(error);
  if (static_cast<long double>(*compression_error) < error) {
    *compression_error = std::nextafter(*compression_error, std::numeric_limits<double>::infinity());
  }
  return std::isfinite(*mass) && *mass > 0.0 && std::isfinite(*compression_error);
}
