// Private rank-one PRODUCT fallback. The analytic extrema enclose the complete
// Gaussian line, including both tails. Probability/libm evaluations carry the
// explicit working-precision allowance below; this is not a libm theorem.
// Included after FactorEventProposal, inside the implementation namespace.
namespace rank_one_bounds {
constexpr double infinity = std::numeric_limits<double>::infinity();
constexpr double zero = -infinity;
constexpr double allowance = 64.0 * std::numeric_limits<double>::epsilon();
struct Interval { double lower, upper; };
struct Refusal {};
// Structural infinities bypass these finite-arithmetic rounders. In particular,
// a finite -DBL_MAX must not become a zero-mass/one-sided structural -Inf.
inline double down(double x) {
  if (!std::isfinite(x)) throw Refusal{};
  const double result = std::nextafter(x, -infinity);
  if (!std::isfinite(result)) throw Refusal{};
  return result;
}
inline double up(double x) {
  if (!std::isfinite(x)) throw Refusal{};
  const double result = std::nextafter(x, infinity);
  if (!std::isfinite(result)) throw Refusal{};
  return result;
}
inline double finite(double x) {
  if (!std::isfinite(x)) throw Refusal{};
  return x;
}
inline Interval checked(Interval x) {
  if (std::isnan(x.lower) || std::isnan(x.upper) || x.upper == infinity ||
      x.lower == infinity || x.lower > x.upper) throw Refusal{};
  return x;
}
inline Interval rounded(double x) {
  finite(x);
  return checked({down(x), up(x)});
}
inline Interval log_value(double x) {
  if (x == 0.0) return {zero, zero};
  if (!(x > 0.0)) throw Refusal{};
  const double value = finite(std::log(x));
  const double delta = allowance * std::max(1.0, std::fabs(value));
  return checked({down(value - delta), up(value + delta)});
}
inline double log_add_bound(double a, double b, bool lower) {
  if (a == zero) return b;
  if (b == zero) return a;
  const double value = finite(log_add_exp(a, b));
  const double delta = allowance * std::max({1.0, std::fabs(a), std::fabs(b), std::fabs(value)});
  return finite(lower ? down(value - delta) : up(value + delta));
}
inline Interval log_add(Interval a, Interval b) {
  return checked({log_add_bound(a.lower, b.lower, true),
                  log_add_bound(a.upper, b.upper, false)});
}
inline Interval log_product(Interval a, Interval b) {
  return checked({a.lower == zero || b.lower == zero ? zero : down(finite(a.lower + b.lower)),
                  a.upper == zero || b.upper == zero ? zero : up(finite(a.upper + b.upper))});
}
inline Interval power(Interval x, int count) {
  if (count <= 0) throw Refusal{};
  return checked({x.lower == zero ? zero : down(finite(count * x.lower)),
                  x.upper == zero ? zero : up(finite(count * x.upper))});
}
inline Interval scaled_endpoint(double x, double se) {
  if (!(se > 0.0) || !std::isfinite(se) || std::isnan(x)) throw Refusal{};
  if (std::isinf(x) || x == 0.0 || se == 1.0) return {x, x};
  return rounded(finite(x * se));
}
inline Interval affine(double mean, double loading, double z) {
  if (loading == 0.0) return {mean, mean};
  if (std::isinf(z)) return {std::copysign(infinity, loading * z),
                           std::copysign(infinity, loading * z)};
  if (z == 0.0) return {mean, mean};
  const Interval product = rounded(finite(loading * z));
  return checked({down(finite(mean + product.lower)), up(finite(mean + product.upper))});
}
inline Interval mean_range(double mean, double loading, double l, double u) {
  const Interval a = affine(mean, loading, l), b = affine(mean, loading, u);
  // Infinite mean endpoints have exact limiting probabilities below.
  return {std::min(a.lower, b.lower), std::max(a.upper, b.upper)};
}
inline Interval midpoint(double l, double u) {
  const Interval a = rounded(finite(l / 2.0)), b = rounded(finite(u / 2.0));
  return checked({down(finite(a.lower + b.lower)), up(finite(a.upper + b.upper))});
}
inline double split_point(double lower, double upper) {
  double split;
  if (lower == -infinity) split = upper - std::max(1.0, std::fabs(upper));
  else if (upper == infinity) split = lower + std::max(1.0, std::fabs(lower));
  else split = lower / 2.0 + upper / 2.0;
  if (!std::isfinite(split) || !(lower < split && split < upper)) throw Refusal{};
  return split;
}
struct Result {
  bool accepted = false;
  double value = std::numeric_limits<double>::quiet_NaN();
  double radius = infinity;
  Interval total{zero, zero};
  std::size_t calls = 0, call_limit = 0, leaves = 0, leaf_limit = 0, bytes = 0;
};
inline std::size_t multiply(std::size_t a, std::size_t b) {
  if (b && a > std::numeric_limits<std::size_t>::max() / b) throw Refusal{};
  return a * b;
}
class Enclosure {
  struct Bin {
    double inner_lower, inner_upper, outer_lower, outer_upper;
    Interval excess, inner_center, center;
    double peak;
    bool finite_bin;
  };
  struct Group { int row, count; std::vector<Bin> bins; };
  struct Leaf { double lower, upper, gap; Interval value; std::size_t order; };
  struct Priority {
    double gap; std::size_t slot, order;
    bool operator<(const Priority &other) const {
      return gap < other.gap || (gap == other.gap && order > other.order);
    }
  };
  const FactorNormalizerContext &context;
  const FactorEventProposal &event;
  Result &result;
  std::vector<Group> groups;
  Interval minimum, maximum, constant{0.0, 0.0};
  void (*poll)();
  std::size_t until_poll = 1000000;

  // Every nonstructural scalar probability call, including latent mass,
  // peaks, inner/outer endpoint evaluations and constants, is charged here.
  Interval probability(double lower, double upper) {
    if (!(lower < upper) || std::isnan(lower) || std::isnan(upper)) throw Refusal{};
    if (lower == -infinity && upper == infinity) return {0.0, 0.0};
    if (result.calls >= result.call_limit) throw Refusal{};
    ++result.calls;
    if (poll && --until_poll == 0) { poll(); until_poll = 1000000; }
    const double value = finite(cpp_selnorm_normal_interval_log_prob(lower, upper, 0.0, 1.0));
    if (value > 0.0) throw Refusal{};
    const double delta = allowance * std::max(1.0, std::fabs(value));
    return checked({down(value - delta), up(value + delta)});
  }
  double at_mean(double lower, double upper, double mean, double sd, bool lower_bound) {
    if (mean == infinity) return upper == infinity ? 0.0 : zero;
    if (mean == -infinity) return lower == -infinity ? 0.0 : zero;
    const auto standardized = [&](double endpoint) {
      if (std::isinf(endpoint)) return Interval{endpoint, endpoint};
      const Interval residual = rounded(finite(endpoint - mean));
      return checked({down(finite(residual.lower / sd)), up(finite(residual.upper / sd))});
    };
    const Interval l = standardized(lower), u = standardized(upper);
    return lower_bound ? probability(l.upper, u.lower).lower : probability(l.lower, u.upper).upper;
  }
  Group prepare_group(int row, int count, std::vector<double> &knots) {
    Group group{row, count, {}};
    const double sd = context.residual_sd[row];
    const double se = context.selection_se[row];
    if (!(sd > 0.0) || !std::isfinite(sd) || !std::isfinite(context.mean[row]) ||
        !std::isfinite(context.loading[row])) throw Refusal{};
    for (int b = 0; b < context.n_bins; ++b) {
      if (!std::isfinite(event.log_excess[b])) continue;
      // event_interval defines the signed event. Carry the endpoint products
      // outward instead of treating its rounded finite endpoints as exact.
      const EventInterval signed_event = event_interval(context.selection.z_lower[b],
        context.selection.z_upper[b], 1.0, context.selection.effect_sign, -infinity, infinity);
      const Interval l = scaled_endpoint(signed_event.lower, se);
      const Interval u = scaled_endpoint(signed_event.upper, se);
      if (!(l.upper < u.lower) || !(l.lower < u.upper)) throw Refusal{};
      const double minimum_weight = *std::min_element(context.omega, context.omega + context.n_bins);
      const double excess = finite(context.omega[b] - minimum_weight);
      Interval excess_log = log_value(excess);
      // Sterbenz subtraction is exact, including the smallest positive
      // subnormal excess. Subtracting zero is exact as well. Outside that
      // proven range carry the subtraction outward rather than lose mass.
      const bool exact_excess = minimum_weight == 0.0 ||
        minimum_weight > std::numeric_limits<double>::max() / 2.0 ||
        context.omega[b] <= 2.0 * minimum_weight;
      if (!exact_excess) {
        const Interval excess_value = rounded(excess);
        if (!(excess_value.lower > 0.0)) throw Refusal{};
        excess_log = {log_value(excess_value.lower).lower, log_value(excess_value.upper).upper};
      }
      Bin bin{l.upper, u.lower, l.lower, u.upper, excess_log, {0.0, 0.0}, {0.0, 0.0}, zero,
              std::isfinite(l.lower) && std::isfinite(u.upper)};
      if (bin.finite_bin) {
        bin.inner_center = midpoint(bin.inner_lower, bin.inner_upper);
        bin.center = midpoint(bin.outer_lower, bin.outer_upper);
        // Width, not a rounded midpoint mean, determines the peak.
        const Interval width = rounded(finite(bin.outer_upper - bin.outer_lower));
        const double half = up(finite(up(finite(width.upper / 2.0)) / sd));
        if (!(half > 0.0)) throw Refusal{};
        bin.peak = probability(-half, half).upper;
      }
      for (double endpoint : {signed_event.lower, signed_event.upper}) {
        if (!std::isfinite(endpoint) || context.loading[row] == 0.0) continue;
        // These centers are refinement hints only. Overflow can omit a hint,
        // never an interval of the complete latent line.
        const double center = (endpoint * se - context.mean[row]) / context.loading[row];
        if (std::isfinite(center)) knots.push_back(center);
      }
      group.bins.push_back(bin);
    }
    return group;
  }
  Interval row_bounds(const Group &group, double lower, double upper) {
    const int row = group.row;
    const Interval means = mean_range(context.mean[row], context.loading[row], lower, upper);
    if (std::isnan(means.lower) || std::isnan(means.upper) || means.lower > means.upper) throw Refusal{};
    Interval total = minimum;
    double scale = 1.0;
    for (const Bin &bin : group.bins) {
      double lo, hi;
      if (!bin.finite_bin) {
        const bool increasing = bin.outer_upper == infinity;
        lo = at_mean(bin.inner_lower, bin.inner_upper,
          increasing ? means.lower : means.upper, context.residual_sd[row], true);
        hi = at_mean(bin.outer_lower, bin.outer_upper,
          increasing ? means.upper : means.lower, context.residual_sd[row], false);
      } else {
        // Inner and outer bins can have different rounded centers. Classify
        // their monotonic regions separately, and include any ambiguous peak.
        if (bin.inner_center.upper < means.lower) {
          lo = at_mean(bin.inner_lower, bin.inner_upper, means.upper, context.residual_sd[row], true);
        } else if (bin.inner_center.lower > means.upper) {
          lo = at_mean(bin.inner_lower, bin.inner_upper, means.lower, context.residual_sd[row], true);
        } else {
          lo = std::min(at_mean(bin.inner_lower, bin.inner_upper, means.lower, context.residual_sd[row], true),
                        at_mean(bin.inner_lower, bin.inner_upper, means.upper, context.residual_sd[row], true));
        }
        if (bin.center.upper < means.lower) {
          hi = at_mean(bin.outer_lower, bin.outer_upper, means.lower, context.residual_sd[row], false);
        } else if (bin.center.lower > means.upper) {
          hi = at_mean(bin.outer_lower, bin.outer_upper, means.upper, context.residual_sd[row], false);
        } else hi = bin.peak;
      }
      const Interval term = log_product({lo, hi}, bin.excess);
      for (double x : {term.lower, term.upper}) if (std::isfinite(x)) scale = std::max(scale, std::fabs(x));
      total = log_add(total, term);
    }
    for (double x : {total.lower, total.upper}) if (std::isfinite(x)) scale = std::max(scale, std::fabs(x));
    // Widen BEFORE intersecting the true global weight range: saturated sums
    // such as .9 + .1 must not create an inverted, falsely exact interval.
    if (total.lower != zero) total.lower = down(finite(total.lower - allowance * scale));
    if (total.upper != zero) total.upper = up(finite(total.upper + allowance * scale));
    total.lower = std::max(total.lower, minimum.lower);
    total.upper = std::min(total.upper, maximum.upper);
    return checked(total);
  }
  Interval segment(double lower, double upper) {
    Interval value = constant;
    for (const Group &group : groups) value = log_product(value, power(row_bounds(group, lower, upper), group.count));
    return log_product(value, probability(lower, upper));
  }
  static double gap(Interval x) {
    if (x.lower == x.upper) return zero;
    if (x.lower == zero) return x.upper;
    return finite(x.upper + std::log(-std::expm1(x.lower - x.upper)));
  }
  bool accept(Interval total, double tolerance) {
    result.total = total;
    if (!std::isfinite(total.lower) || !std::isfinite(total.upper)) return false;
    const double value = finite(log_add_exp(total.lower, total.upper) - std::log(2.0));
    if (value < total.lower || value > total.upper) return false;
    const double dl = up(finite(value - total.lower)), du = down(finite(value - total.upper));
    const double a = std::expm1(dl), b = -std::expm1(du);
    if (!std::isfinite(a) || !std::isfinite(b) || a < 0.0 || b < 0.0) return false;
    const double radius = up(std::max(up(a + allowance * a), up(b + allowance * b)));
    if (!std::isfinite(radius) || radius > tolerance) return false;
    result.value = value; result.radius = radius; result.accepted = true;
    return true;
  }
 public:
  Enclosure(const FactorNormalizerContext &value, const FactorEventProposal &proposal,
            Result &output, void (*poll_function)())
    : context(value), event(proposal), result(output), poll(poll_function) {}
  void run(int max_points, int scrambles, double tolerance) {
    const std::size_t points = multiply(multiply(2, static_cast<std::size_t>(max_points)), scrambles);
    std::size_t bins = 0;
    for (double x : event.log_excess) if (std::isfinite(x)) ++bins;
    result.call_limit = multiply(multiply(points, event.groups.size()), bins);
    result.leaf_limit = std::min(points, static_cast<std::size_t>(SELNORM_FACTOR_NODE_BUDGET));
    if (!result.call_limit || result.leaf_limit < 2) throw Refusal{};
    minimum = log_value(*std::min_element(context.omega, context.omega + context.n_bins));
    maximum = log_value(*std::max_element(context.omega, context.omega + context.n_bins));
    std::vector<double> knots{0.0};
    groups.reserve(event.groups.size());
    for (const auto &group : event.groups) groups.push_back(prepare_group(group.row, group.count, knots));
    // Factor-free rows remain constant; their enclosure arithmetic is charged.
    // Existing all-factor-free exact returns precede this helper.
    for (int row = 0; row < context.dimension; ++row) {
      if (context.loading[row] != 0.0) continue;
      const Group group = prepare_group(row, 1, knots);
      constant = log_product(constant, row_bounds(group, 0.0, 0.0));
    }
    std::sort(knots.begin(), knots.end());
    knots.erase(std::unique(knots.begin(), knots.end()), knots.end());
    if (knots.size() >= result.leaf_limit) throw Refusal{};
    knots.insert(knots.begin(), -infinity); knots.push_back(infinity);
    // leaf_limit caps ACTIVE slots; tree has exactly 2*next_power_of_two slots
    // and the heap at most leaf_limit entries. No per-leaf group/bin caches.
    std::size_t base = 1;
    while (base < result.leaf_limit) base = multiply(base, 2);
    const std::size_t tree_count = multiply(base, 2);
    result.bytes = multiply(result.leaf_limit, sizeof(Leaf) + sizeof(Priority)) + multiply(tree_count, sizeof(Interval));
    if (result.bytes < multiply(tree_count, sizeof(Interval)) || base > static_cast<std::size_t>(std::numeric_limits<int>::max())) throw Refusal{};
    std::vector<Interval> tree(tree_count, {zero, zero});
    std::vector<Leaf> leaves; leaves.reserve(result.leaf_limit);
    std::vector<Priority> heap; heap.reserve(result.leaf_limit);
    const auto update = [&](std::size_t slot, Interval value) {
      std::size_t at = base + slot; tree[at] = value;
      while (at > 1) { at /= 2; tree[at] = log_add(tree[2 * at], tree[2 * at + 1]); }
    };
    const auto push = [&](std::size_t slot) {
      const Leaf &leaf = leaves[slot];
      heap.push_back({leaf.gap, slot, leaf.order});
      std::push_heap(heap.begin(), heap.end());
    };
    std::size_t order = 0;
    for (std::size_t i = 1; i < knots.size(); ++i) {
      const Interval value = segment(knots[i - 1], knots[i]);
      leaves.push_back({knots[i - 1], knots[i], gap(value), value, order++});
      tree[base + leaves.size() - 1] = value; push(leaves.size() - 1);
    }
    for (std::size_t i = base - 1; i > 0; --i) tree[i] = log_add(tree[2 * i], tree[2 * i + 1]);
    result.leaves = leaves.size();
    while (!accept(tree[1], tolerance)) {
      if (tree[1].upper == zero || leaves.size() >= result.leaf_limit || heap.empty()) throw Refusal{};
      std::pop_heap(heap.begin(), heap.end());
      const std::size_t slot = heap.back().slot; heap.pop_back();
      const Leaf parent = leaves[slot];
      if (parent.gap == zero) throw Refusal{};
      const double split = split_point(parent.lower, parent.upper);
      const Interval left = segment(parent.lower, split), right = segment(split, parent.upper);
      leaves[slot] = {parent.lower, split, gap(left), left, order++};
      const std::size_t next = leaves.size();
      leaves.push_back({split, parent.upper, gap(right), right, order++});
      update(slot, left); update(next, right); push(slot); push(next);
      result.leaves = leaves.size();
    }
  }
};
inline Result integrate(const FactorNormalizerContext &context, const FactorEventProposal &event,
                        int max_points, int scrambles, double tolerance, void (*poll)() = nullptr) {
  Result result;
  try { Enclosure(context, event, result, poll).run(max_points, scrambles, tolerance); }
  catch (const Refusal &) {}
  catch (const std::bad_alloc &) {}
  catch (const std::length_error &) {}
  return result;
}
} // namespace rank_one_bounds
