// Native port of the zplot context projection controller; see zplot-context.h.
//
// Every expression below mirrors one in R/zplot-context.R (.zplot_context_projection,
// .zplot_context_geometry, .zplot_context_gaussian_bound,
// .selection_transition_rule_resolves_rows) or R/selection-likelihood.R
// (.selection_transition_rule_resolves). Binary arithmetic keeps R's
// left-to-right evaluation; reductions use the long-double accumulators of R's
// sum(), mean(), colSums(), colMeans() and rowSums(); linear algebra calls the
// routines R calls. Comments name the R statement where the mapping is not
// literal.
#include "zplot-context.h"

#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <map>
#include <memory>
#include <new>
#include <tuple>
#include <utility>

#include "selnorm.h"
#include "selnorm-mv.h"
#include "selnorm-compact-partition.h"
#include "selnorm-parallel.h"

#ifndef FCONE
# define FCONE
#endif

// Rmath, declared directly so that its unprefixed macros stay out of this
// file. Workers call only pnorm and dnorm, with a unit or positive finite
// scale, for which they cannot warn. pgamma and lgamma (the Gaussian
// resolution bounds) can reach non-convergence warnings, so they run on the
// main thread only: a draw that needs a bound it does not have suspends and
// the batch driver computes it between parallel rounds.
extern "C" {
double Rf_dnorm4(double, double, double, int);
double Rf_pnorm5(double, double, double, int, int);
double Rf_pgamma(double, double, double, int, int);
double Rf_lgammafn(double);
}

namespace {

constexpr double kInf = std::numeric_limits<double>::infinity();
constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
const double kNA = [] {
  // R's NA_real_: a quiet NaN whose low word is 1954.
  union { double value; std::uint32_t word[2]; } bits;
  bits.value = kNaN;
  bits.word[0] = 1954;
  return bits.value;
}();

bool r_is_na(double x)
{
  if (!std::isnan(x)) return false;
  union { double value; std::uint32_t word[2]; } bits;
  bits.value = x;
  return bits.word[0] == 1954;
}

// ---------------------------------------------------------------------------
// R reductions (R 4.6.0 src/main/summary.c and array.c)
// ---------------------------------------------------------------------------

// sum(x): rsum() with na.rm = FALSE.
double r_sum(const double *x, std::size_t n)
{
  long double s = 0.0L;
  for (std::size_t i = 0; i < n; ++i) s += x[i];
  if (s > DBL_MAX) return kInf;
  if (s < -DBL_MAX) return -kInf;
  return static_cast<double>(s);
}
double r_sum(const std::vector<double> &x) { return r_sum(x.data(), x.size()); }

// mean(x): real_mean().
double r_mean(const double *x, std::size_t n)
{
  long double s = 0.0L;
  for (std::size_t i = 0; i < n; ++i) s += x[i];
  const bool finite_s = std::isfinite(static_cast<double>(s));
  const long double count = static_cast<long double>(static_cast<double>(n));
  if (finite_s) {
    s /= count;
  } else {
    s = 0.0L;
    for (std::size_t i = 0; i < n; ++i) s += x[i] / count;
  }
  if (finite_s && std::isfinite(static_cast<double>(s))) {
    long double t = 0.0L;
    for (std::size_t i = 0; i < n; ++i) t += (x[i] - s);
    s += t / count;
  } else if (std::isfinite(static_cast<double>(s))) {
    long double t = 0.0L;
    for (std::size_t i = 0; i < n; ++i) t += (x[i] - s) / count;
    s += t;
  }
  return static_cast<double>(s);
}
double r_mean(const std::vector<double> &x) { return r_mean(x.data(), x.size()); }

// max(x) / min(x): rmax() / rmin() with na.rm = FALSE; x is nonempty.
double r_max(const double *x, std::size_t n)
{
  double s = 0.0;
  bool updated = false;
  for (std::size_t k = 0; k < n; ++k) {
    if (std::isnan(x[k])) {
      if (!r_is_na(s)) s = x[k];
      updated = true;
    } else if (!updated || x[k] > s) {
      s = x[k];
      updated = true;
    }
  }
  return s;
}
double r_max(const std::vector<double> &x) { return r_max(x.data(), x.size()); }
double r_min(const double *x, std::size_t n)
{
  double s = 0.0;
  bool updated = false;
  for (std::size_t k = 0; k < n; ++k) {
    if (std::isnan(x[k])) {
      if (!r_is_na(s)) s = x[k];
      updated = true;
    } else if (!updated || x[k] < s) {
      s = x[k];
      updated = true;
    }
  }
  return s;
}
double r_min(const std::vector<double> &x) { return r_min(x.data(), x.size()); }
double r_max2(double a, double b) { const double v[2] = {a, b}; return r_max(v, 2); }
double r_min2(double a, double b) { const double v[2] = {a, b}; return r_min(v, 2); }

// pmax(a, b) / pmin(a, b) for doubles (do_pmin with na.rm = FALSE).
double r_pmax(double a, double b)
{
  double r = a;
  if (b > r || std::isnan(b)) r = b;
  return r;
}
double r_pmin(double a, double b)
{
  double r = a;
  if (b < r || std::isnan(b)) r = b;
  return r;
}

// colSums / colMeans of an n x p column-major matrix.
std::vector<double> r_col_sums(const double *x, int n, int p, bool means = false)
{
  std::vector<double> out(p);
  for (int j = 0; j < p; ++j) {
    long double sum = 0.0L;
    const double *rx = x + static_cast<std::size_t>(n) * j;
    for (int i = 0; i < n; ++i) sum += *rx++;
    if (means) sum /= static_cast<long double>(n);
    out[j] = static_cast<double>(sum);
  }
  return out;
}

// rowSums of an n x p column-major matrix.
std::vector<double> r_row_sums(const double *x, int n, int p)
{
  std::vector<long double> sums(n, 0.0L);
  for (int j = 0; j < p; ++j) {
    const double *rx = x + static_cast<std::size_t>(n) * j;
    for (int i = 0; i < n; ++i) sums[i] += rx[i];
  }
  std::vector<double> out(n);
  for (int i = 0; i < n; ++i) out[i] = static_cast<double>(sums[i]);
  return out;
}

// which.min(x) / which.max(x): first extreme, NaN skipped; -1 if none.
int r_which_min(const std::vector<double> &x)
{
  int index = -1;
  double s = kInf;
  for (std::size_t i = 0; i < x.size(); ++i) {
    if (!std::isnan(x[i]) && (x[i] < s || index == -1)) {
      s = x[i];
      index = static_cast<int>(i);
    }
  }
  return index;
}
int r_which_max(const std::vector<double> &x)
{
  int index = -1;
  double s = -kInf;
  for (std::size_t i = 0; i < x.size(); ++i) {
    if (!std::isnan(x[i]) && (x[i] > s || index == -1)) {
      s = x[i];
      index = static_cast<int>(i);
    }
  }
  return index;
}

// ---------------------------------------------------------------------------
// R linear algebra (src/main/array.c, src/modules/lapack/Lapack.c)
// ---------------------------------------------------------------------------

bool may_have_nan_or_inf(const double *x, std::size_t n)
{
  if ((n & 1) != 0 && !std::isfinite(x[0])) return true;
  for (std::size_t i = n & 1; i < n; i += 2) {
    if (!std::isfinite(x[i] + x[i + 1])) return true;
  }
  return false;
}

// x %*% y for n x n matrices with n >= 2 under the default 'matprod'.
std::vector<double> r_matprod(const std::vector<double> &x, const std::vector<double> &y, int n)
{
  std::vector<double> z(static_cast<std::size_t>(n) * n);
  const std::size_t N = static_cast<std::size_t>(n);
  if (may_have_nan_or_inf(x.data(), N * n) || may_have_nan_or_inf(y.data(), N * n)) {
    for (int i = 0; i < n; ++i) {
      for (int k = 0; k < n; ++k) {
        double sum = 0.0;
        for (int j = 0; j < n; ++j) sum += x[i + j * N] * y[j + k * N];
        z[i + k * N] = sum;
      }
    }
    return z;
  }
  const char *trans = "N";
  const double one = 1.0, zero = 0.0;
  F77_CALL(dgemm)(trans, trans, &n, &n, &n, &one, x.data(), &n, y.data(), &n, &zero, z.data(), &n
    FCONE FCONE);
  return z;
}

// tcrossprod(x) for an n x n matrix under the default 'matprod'.
std::vector<double> r_tcrossprod(const std::vector<double> &x, int n)
{
  std::vector<double> z(static_cast<std::size_t>(n) * n);
  const std::size_t N = static_cast<std::size_t>(n);
  if (may_have_nan_or_inf(x.data(), N * n)) {
    for (int i = 0; i < n; ++i) {
      for (int k = 0; k < n; ++k) {
        double sum = 0.0;
        for (int j = 0; j < n; ++j) sum += x[i + j * N] * x[k + j * N];
        z[i + k * N] = sum;
      }
    }
    return z;
  }
  const char *trans = "N", *uplo = "U";
  const double one = 1.0, zero = 0.0;
  F77_CALL(dsyrk)(uplo, trans, &n, &n, &one, x.data(), &n, &zero, z.data(), &n FCONE FCONE);
  for (int i = 1; i < n; ++i)
    for (int j = 0; j < i; ++j) z[i + N * j] = z[j + N * i];
  return z;
}

// norm(x, "F").
double r_norm_frobenius(const std::vector<double> &x, int n)
{
  const char *type = "F";
  return F77_CALL(dlange)(type, &n, &n, x.data(), &n, nullptr FCONE);
}

std::vector<double> r_abs(const std::vector<double> &x)
{
  std::vector<double> out(x.size());
  for (std::size_t i = 0; i < x.size(); ++i) out[i] = std::fabs(x[i]);
  return out;
}

// ---------------------------------------------------------------------------
// Small R vector helpers used by the controller
// ---------------------------------------------------------------------------

// log_sum() of the controller.
double log_sum(const std::vector<double> &values)
{
  const double largest = r_max(values);
  if (std::isinf(largest) && largest < 0) return -kInf;
  std::vector<double> shifted(values.size());
  for (std::size_t i = 0; i < values.size(); ++i) shifted[i] = std::exp(values[i] - largest);
  return largest + std::log(r_sum(shifted));
}

// R's log(x) for doubles.
double r_log(double x) { return x > 0 ? std::log(x) : (x == 0 ? -kInf : kNaN); }

double gaussian_tail(const std::vector<double> &nodes)
{
  return log_sum({Rf_pnorm5(r_min(nodes), 0.0, 1.0, 1, 1),
                  Rf_pnorm5(r_max(nodes), 0.0, 1.0, 0, 1)});
}

std::vector<double> normal_interval(const std::vector<double> &lower, const std::vector<double> &upper)
{
  std::vector<double> out(lower.size());
  for (std::size_t i = 0; i < lower.size(); ++i) {
    double big, small;
    if (lower[i] >= 0) {
      big = Rf_pnorm5(lower[i], 0.0, 1.0, 0, 1);
      small = Rf_pnorm5(upper[i], 0.0, 1.0, 0, 1);
    } else {
      big = Rf_pnorm5(upper[i], 0.0, 1.0, 1, 1);
      small = Rf_pnorm5(lower[i], 0.0, 1.0, 1, 1);
    }
    out[i] = big + r_log(-std::expm1(small - big));
  }
  return out;
}

double log_abs_expm1(double value)
{
  if (value > std::log(2.0)) return value + std::log1p(-std::exp(-value));
  return r_log(std::fabs(std::expm1(value)));
}

double log_add(double first, double second)
{
  const double largest = r_pmax(first, second);
  double value = largest + std::log1p(std::exp(-std::fabs(first - second)));
  if (largest == -kInf) value = -kInf;
  return value;
}

[[noreturn]] void fail(const std::string &message) { throw ZplotNativeFailure(message); }
void stop_if_not(bool condition, const char *expression)
{
  if (!condition) fail(std::string(expression) + " is not TRUE");
}

// ---------------------------------------------------------------------------
// Transition-rule checks
// ---------------------------------------------------------------------------

// selection$z_lower[seq_len(length(omega) - 1L)] filtered to finite cuts at
// weight changes.
std::vector<double> transition_cuts(const ZplotProjectionSetup &setup, const double *omega, int bins)
{
  std::vector<double> cuts;
  for (int bin = 0; bin + 1 < bins; ++bin) {
    const double cut = bin < setup.rule_z_lower_length ? setup.rule_z_lower[bin] : kNA;
    if (std::isfinite(cut) && omega[bin] != omega[bin + 1]) cuts.push_back(cut);
  }
  return cuts;
}

// Shared preamble: sorted nodes, validity and the coupled gap. FALSE in
// 'valid' means the rule returns FALSE.
struct RuleNodes {
  std::vector<double> nodes;
  bool valid = false;
  double coupled_gap = kInf;
};

RuleNodes rule_nodes(const std::vector<double> &input)
{
  RuleNodes out;
  for (double node : input) if (!std::isnan(node)) out.nodes.push_back(node);
  std::stable_sort(out.nodes.begin(), out.nodes.end());
  const std::vector<double> &nodes = out.nodes;
  if (nodes.size() < 2) return out;
  for (double node : nodes) if (!std::isfinite(node)) return out;
  for (std::size_t i = 1; i < nodes.size(); ++i) if (nodes[i] - nodes[i - 1] <= 0) return out;
  bool any_coupled = false;
  double gap = -kInf;
  std::vector<double> gaps;
  for (std::size_t i = 0; i + 1 < nodes.size(); ++i) {
    if (nodes[i] <= 3 && nodes[i + 1] >= -3) gaps.push_back(nodes[i + 1] - nodes[i]);
  }
  if (!gaps.empty()) {
    any_coupled = true;
    gap = r_max(gaps);
  }
  out.coupled_gap = any_coupled ? gap : kInf;
  out.valid = true;
  return out;
}

// findInterval(x, nodes) for nodes[1] < x < nodes[n]: one-based i with
// nodes[i] <= x < nodes[i + 1].
int find_interval(const std::vector<double> &nodes, double x)
{
  return static_cast<int>(std::upper_bound(nodes.begin(), nodes.end(), x) - nodes.begin());
}

// loading: k x columns, column-major. axes: zero-based columns.
bool axis_coupled(const std::vector<double> &loading, int k, int columns,
                  const std::vector<int> &active, int axis)
{
  if (columns <= 1) return false;
  for (int row : active) {
    for (int other = 0; other < columns; ++other) {
      if (other != axis && loading[row + static_cast<std::size_t>(k) * other] != 0) return true;
    }
  }
  return false;
}

// .selection_transition_rule_resolves_rows(nodes, means, residual_sd,
// loading, sei, selection, axes); means: n_rows x k column-major.
bool transition_rule_resolves_rows(const ZplotProjectionSetup &setup, const double *omega,
    const std::vector<double> &input_nodes, const std::vector<double> &means, int n_rows,
    const std::vector<double> &residual_sd, const std::vector<double> &loading, int columns,
    const std::vector<int> &axes)
{
  const int k = setup.k;
  if (n_rows == 0) return true;
  const std::vector<double> cuts = transition_cuts(setup, omega, setup.bins);
  if (cuts.empty() || axes.empty()) return true;
  const RuleNodes rule = rule_nodes(input_nodes);
  if (!rule.valid) return false;
  const std::vector<double> &nodes = rule.nodes;
  const int n_nodes = static_cast<int>(nodes.size());
  for (int axis : axes) {
    std::vector<int> active;
    for (int row = 0; row < k; ++row) {
      if (loading[row + static_cast<std::size_t>(k) * axis] != 0) active.push_back(row);
    }
    if (active.empty()) continue;
    for (int row : active) {
      if (!std::isfinite(loading[row + static_cast<std::size_t>(k) * axis])) return false;
    }
    for (int row : active) {
      for (int r = 0; r < n_rows; ++r) {
        if (!std::isfinite(means[r + static_cast<std::size_t>(n_rows) * row])) return false;
      }
    }
    for (int row : active) {
      if (!std::isfinite(residual_sd[row]) || residual_sd[row] <= 0) return false;
    }
    if (axis_coupled(loading, k, columns, active, axis)) {
      for (int row : active) {
        if (rule.coupled_gap * std::fabs(loading[row + static_cast<std::size_t>(k) * axis]) >
            residual_sd[row]) return false;
      }
      continue;
    }
    for (double cut : cuts) {
      // center: n_rows x active column-major; anyNA() over all of it first.
      std::vector<double> center(static_cast<std::size_t>(n_rows) * active.size());
      for (std::size_t a = 0; a < active.size(); ++a) {
        const int row = active[a];
        const double shift = setup.rule_sign * cut * setup.sei[row];
        const double coefficient = loading[row + static_cast<std::size_t>(k) * axis];
        for (int r = 0; r < n_rows; ++r) {
          center[r + static_cast<std::size_t>(n_rows) * a] =
            (shift - means[r + static_cast<std::size_t>(n_rows) * row]) / coefficient;
        }
      }
      for (double value : center) if (std::isnan(value)) return false;
      for (std::size_t index = 0; index < center.size(); ++index) {
        const double value = center[index];
        if (!(value > nodes[0] && value < nodes[n_nodes - 1])) continue;
        const int lower = find_interval(nodes, value);
        double gap = nodes[lower] - nodes[lower - 1];
        if (nodes[lower - 1] == value && lower > 1) {
          gap = r_pmax(gap, nodes[lower - 1] - nodes[lower - 2]);
        }
        const int row = active[index / n_rows];
        if (gap * std::fabs(loading[row + static_cast<std::size_t>(k) * axis]) > residual_sd[row]) {
          return false;
        }
      }
    }
  }
  return true;
}

// .selection_transition_rule_resolves(nodes, mean, residual_sd, loading, sei,
// selection, axes).
bool transition_rule_resolves(const ZplotProjectionSetup &setup, const double *omega,
    const std::vector<double> &input_nodes, const double *mean,
    const std::vector<double> &residual_sd, const std::vector<double> &loading, int columns,
    const std::vector<int> &axes)
{
  const int k = setup.k;
  const std::vector<double> cuts = transition_cuts(setup, omega, setup.bins);
  if (cuts.empty() || axes.empty()) return true;
  const RuleNodes rule = rule_nodes(input_nodes);
  if (!rule.valid) return false;
  const std::vector<double> &nodes = rule.nodes;
  const int n_nodes = static_cast<int>(nodes.size());
  for (int axis : axes) {
    std::vector<int> active;
    for (int row = 0; row < k; ++row) {
      if (loading[row + static_cast<std::size_t>(k) * axis] != 0) active.push_back(row);
    }
    if (active.empty()) continue;
    for (int row : active) if (!std::isfinite(loading[row + static_cast<std::size_t>(k) * axis])) return false;
    for (int row : active) if (!std::isfinite(mean[row])) return false;
    for (int row : active) {
      if (!std::isfinite(residual_sd[row]) || residual_sd[row] <= 0) return false;
    }
    if (axis_coupled(loading, k, columns, active, axis)) {
      for (int row : active) {
        if (rule.coupled_gap * std::fabs(loading[row + static_cast<std::size_t>(k) * axis]) >
            residual_sd[row]) return false;
      }
      continue;
    }
    for (int row : active) {
      const double coefficient = loading[row + static_cast<std::size_t>(k) * axis];
      for (double cut : cuts) {
        const double center = (setup.rule_sign * cut * setup.sei[row] - mean[row]) / coefficient;
        if (std::isnan(center)) return false;
        if (center <= nodes[0] || center >= nodes[n_nodes - 1]) continue;
        const int lower = find_interval(nodes, center);
        double gap = nodes[lower] - nodes[lower - 1];
        if (nodes[lower - 1] == center && lower > 1) gap = r_max2(gap, nodes[lower - 1] - nodes[lower - 2]);
        if (gap * std::fabs(coefficient) > residual_sd[row]) return false;
      }
    }
  }
  return true;
}

// .zplot_context_gaussian_bound(beta, order). Main thread only (Rmath pgamma).
std::vector<double> gaussian_bound(const std::vector<double> &beta, int order)
{
  std::vector<double> result(beta.size(), 0.0);
  for (double value : beta) {
    if (std::isnan(value) || value < 0) fail("Gaussian resolution ratios must be nonnegative.");
  }
  bool any_active = false;
  for (std::size_t i = 0; i < beta.size(); ++i) {
    if (std::isinf(beta[i])) result[i] = kInf;
    if (std::isfinite(beta[i]) && beta[i] > 0) any_active = true;
  }
  if (!any_active) return result;
  const double o = static_cast<double>(order);
  const double scale = (std::log(2.0) + Rf_lgammafn(2 * o + 1) - Rf_lgammafn(o + 1)) / o -
    std::log(4.0);
  for (std::size_t i = 0; i < beta.size(); ++i) {
    if (!(std::isfinite(beta[i]) && beta[i] > 0)) continue;
    const double log_beta = std::log(beta[i]);
    const double log_x = scale - log_beta;
    const double x = std::exp(log_x);
    const double first = o * log_beta + Rf_pgamma(x, o + .5, 1.0, 1, 1);
    const double second = std::log(4.0) + Rf_pnorm5(-std::sqrt(2 * x), 0.0, 1.0, 1, 1);
    const double largest = r_pmax(first, second);
    result[i] = std::exp(largest + std::log(std::exp(first - largest) + std::exp(second - largest)));
  }
  return result;
}

// ---------------------------------------------------------------------------
// The controller
// ---------------------------------------------------------------------------

struct ContextTable {
  std::vector<double> means;   // rows x k
  int rows = 0;
  std::vector<double> nodes;   // the context rule
  std::vector<double> log_weights;
  std::vector<double> interval;  // empty for point contexts
  std::vector<int> prune_order;  // zero-based
  double mass = 0.0;
  double covariance_error = 0.0;
  double log_density_bound = 0.0;
  double outer_tail = 0.0;
};

struct NormalizerRecord {
  std::vector<double> log_A;
  std::vector<double> mass_error;
  std::vector<double> numerator_error;
};

struct Pruned {
  std::vector<int> keep;  // zero-based
  double mass = 0.0;
  double retained_mass = 0.0;
  double error = 0.0;
};

struct Evaluation {
  std::vector<double> density;
  double mass = 0.0;
  double compression_error = 0.0;
  double component_omission_mass = 0.0;
  double cdf_error = 0.0;
  double mass_cdf_error = 0.0;
  std::vector<double> denominator_error;
  double tail = 0.0;
  double mass_tail = 0.0;
  std::vector<double> log_D;
  const ContextTable *table = nullptr;
  std::vector<int> keep;
  double omitted_mass = 0.0;
  double retained_mass = 0.0;
  double pruning_error = 0.0;
};

struct InnerTable {
  std::vector<double> density;
  double error = 0.0;
  double mass_error = 0.0;
  double mass_omission_error = 0.0;
  double inner_change = 0.0;
  double normalizer = 0.0;
  double compressed = 0.0;
  double cdf = 0.0;
  double covariance = 0.0;
  double tail = 0.0;
  double pruning = 0.0;
  double outer_tail = 0.0;
  const Evaluation *fine = nullptr;
};

struct Geometry {
  std::vector<double> loading;  // k x columns (used columns only)
  int columns = 0;
  std::vector<double> residual_sd;
  std::vector<double> residual_variance;
  int groups = 0;
  std::vector<int> type_index;  // one-based
  double relative_error = 0.0;
};

// .zplot_context_geometry(); false for NULL.
bool context_geometry(int k, const std::vector<double> &covariance, const std::vector<double> &packed,
                      const double *sei, const std::vector<double> &cholesky, Geometry *out)
{
  ZplotEnvelopeComponents envelope;
  zplot_covariance_envelope_components(packed.data(), k, sei, false, &envelope);
  if (!envelope.available) return false;
  const int rank = envelope.rank;
  std::vector<int> used;
  for (int column = 0; column < rank; ++column) {
    int nonzero = 0;
    for (int row = 0; row < k; ++row) {
      if (envelope.loading[row + static_cast<std::size_t>(k) * column] != 0) ++nonzero;
    }
    if (nonzero > 0) used.push_back(column);
  }
  const int columns = static_cast<int>(used.size());
  if (columns == 0 || columns > k) return false;
  out->loading.assign(static_cast<std::size_t>(k) * columns, 0.0);
  for (int c = 0; c < columns; ++c) {
    std::copy_n(envelope.loading.data() + static_cast<std::size_t>(k) * used[c], k,
      out->loading.data() + static_cast<std::size_t>(k) * c);
  }
  out->columns = columns;
  const std::vector<double> &radius = envelope.covariance_residual_radius;
  for (double value : radius) if (!std::isfinite(value) || value < 0) return false;
  const std::size_t K = static_cast<std::size_t>(k);
  std::vector<double> sum(K * K);
  for (std::size_t i = 0; i < K * K; ++i) sum[i] = std::fabs(envelope.covariance_residual[i]) + radius[i];
  const double difference = r_norm_frobenius(sum, k);
  std::vector<double> identity(K * K, 0.0);
  for (std::size_t i = 0; i < K; ++i) identity[i + K * i] = 1.0;
  const auto gamma = [](int operations) {
    const double product = operations * DBL_EPSILON / 2;
    return product >= 1 ? kInf : product / (1 - product);
  };
  // factor <- t(cholesky); inverse <- forwardsolve(factor, identity)
  std::vector<double> factor(K * K);
  for (std::size_t i = 0; i < K; ++i)
    for (std::size_t j = 0; j < K; ++j) factor[i + K * j] = cholesky[j + K * i];
  for (std::size_t i = 0; i < K; ++i) {
    if (factor[i * (K + 1)] == 0.0) {
      fail("singular matrix in 'backsolve'. First zero in diagonal [" + std::to_string(i + 1) + "]");
    }
  }
  std::vector<double> inverse = identity;
  {
    const char *side = "L", *uplo = "L", *trans = "N", *diag = "N";
    const double one = 1.0;
    F77_CALL(dtrsm)(side, uplo, trans, diag, &k, &k, &one, factor.data(), &k, inverse.data(), &k
      FCONE FCONE FCONE FCONE);
  }
  const std::vector<double> product = r_matprod(factor, inverse, k);
  std::vector<double> residual(K * K);
  for (std::size_t i = 0; i < K * K; ++i) residual[i] = identity[i] - product[i];
  std::vector<double> absolute_product = r_matprod(r_abs(factor), r_abs(inverse), k);
  for (std::size_t i = 0; i < K * K; ++i) absolute_product[i] = absolute_product[i] + identity[i];
  const double inverse_residual = r_norm_frobenius(residual, k) +
    gamma(k + 2) * r_norm_frobenius(absolute_product, k);
  if (!std::isfinite(inverse_residual) || inverse_residual >= 1) return false;
  const double inverse_bound = r_norm_frobenius(inverse, k) / (1 - inverse_residual);
  std::vector<double> reconstruction = r_tcrossprod(factor, k);
  for (std::size_t i = 0; i < K * K; ++i) reconstruction[i] = reconstruction[i] - covariance[i];
  std::vector<double> absolute_reconstruction = r_tcrossprod(r_abs(factor), k);
  for (std::size_t i = 0; i < K * K; ++i) {
    absolute_reconstruction[i] = absolute_reconstruction[i] + std::fabs(covariance[i]);
  }
  const double reconstruction_bound = r_norm_frobenius(reconstruction, k) +
    gamma(k + 4) * r_norm_frobenius(absolute_reconstruction, k);
  const double relative_reconstruction = inverse_bound * inverse_bound * reconstruction_bound;
  if (!std::isfinite(relative_reconstruction) || relative_reconstruction >= 1) return false;
  const double epsilon = inverse_bound * inverse_bound * difference / (1 - relative_reconstruction);
  if (!std::isfinite(epsilon) || epsilon < 0 || epsilon >= 1) return false;
  out->residual_sd = envelope.residual_sd;
  out->residual_variance = envelope.residual_variance;
  out->groups = envelope.groups;
  out->type_index = envelope.type_index;
  out->relative_error = epsilon;
  return true;
}

const ZplotGaussHermiteRule &quadrature_rule(const ZplotProjectionSetup &setup, int order)
{
  for (const ZplotGaussHermiteRule &rule : *setup.rules) if (rule.order == order) return rule;
  fail("Internal error: Gauss-Hermite rule of order " + std::to_string(order) + " was not supplied.");
}


// Thrown at a resumable point when a Gaussian resolution bound is missing.
struct NeedBound {};

// One draw's controller. It runs in stages so that a parallel batch can stop
// a draw where it needs a Gaussian resolution bound (computed on the main
// thread) and resume it there: the bounds are requested only at the top of the
// final-resolution check and of a ladder iteration, before that step changes
// any state.
class ZplotContextController {
 public:
  ZplotContextController(const ZplotProjectionSetup &setup, const ZplotProjectionDraw &draw,
                         bool main_thread, ZplotProjectionOutcome *outcome)
    : s_(setup), d_(draw), main_thread_(main_thread), out_(outcome), k_(setup.k)
  {
    *out_ = ZplotProjectionOutcome();
  }

  // Runs until the draw finishes (TRUE) or needs bounds (FALSE). Never throws
  // except for an R unwind from a main-thread interrupt check.
  bool step();
  // Main thread: computes the requested bounds. FALSE when that fails, with
  // the failure recorded in the outcome.
  bool supply_bounds();

 private:
  enum Stage { stage_setup, stage_resolve, stage_ladder, stage_done };

  // Inputs
  const ZplotProjectionSetup &s_;
  const ZplotProjectionDraw &d_;
  bool main_thread_;
  ZplotProjectionOutcome *out_;
  int k_;
  Stage stage_ = stage_setup;
  const int context_count_ = static_cast<int>(s_.context_orders.size());
  const int inner_count_ = static_cast<int>(s_.inner_orders.size());
  std::map<int, std::vector<double>> factor_bounds_, context_bounds_;
  std::vector<int> factor_requests_, context_requests_;

  // Draw state
  std::vector<double> omega_;
  std::vector<double> covariance_;
  std::vector<double> cholesky_;
  std::vector<double> precision_diag_;
  std::vector<double> packed_;
  std::vector<double> point_means_;
  int point_rows_ = 1;
  std::vector<double> context_factor_;  // factor_rows x k
  bool point_context_ = false;
  bool cdf_enabled_ = true;
  Geometry geometry_;
  double epsilon_actual_ = 0, epsilon_ = 0, beta_ = 0, log_weight_ratio_ = 0;
  double covariance_log_A_bound_ = 0;
  std::vector<double> row_density_bound_;
  double compact_density_bound_ = 0;
  int integration_rank_ = 1;
  bool reuse_normalizer_ = false;
  double tolerance_ = 0;
  double grid_min_ = 0, grid_max_ = 0;
  std::vector<double> factor_beta_;  // k x columns
  std::vector<double> context_beta_;
  std::vector<double> height_ratio_;
  std::vector<double> context_factor_vector_;  // as.numeric(context_factor)
  int denominator_index_ = 3;  // one-based, as in R
  int context_index_ = 2;
  int inner_index_ = 2;
  int compression_refinements_ = 0;
  int pruning_refinements_ = 0;
  double allowance_ = 0;
  double pruning_allowance_ = 0;
  double reference_peak_ = 0;
  bool cache_raw_rules_ = false;
  double raw_rule_bytes_ = 0;
  const double raw_rule_budget_ = 128.0 * 1024 * 1024;

  // Memo tables
  std::map<std::pair<int, int>, std::vector<double>> resolutions_;
  std::map<int, ContextTable> tables_;
  std::map<std::tuple<int, int, bool>, NormalizerRecord> normalizers_;
  std::map<std::tuple<int, int, int, bool, double, double>, Evaluation> results_;
  std::map<std::tuple<int, int, bool, double>, std::vector<double>> raw_rules_;
  std::map<std::tuple<int, int, int>, std::vector<bool>> selection_resolutions_;
  std::map<std::pair<int, int>, bool> sweeps_;

  void setup();
  void resolve();
  void ladder();

  double determinant_bound(int dimension) const { return -.5 * dimension * std::log1p(-epsilon_); }
  double tilted_moment(int dimension) const
  {
    return 4 * dimension * log_weight_ratio_ + 2 * dimension * std::log(2.0);
  }
  double max_omega() const { return r_max(omega_); }
  double min_omega() const { return r_min(omega_); }

  const std::vector<double> &gaussian_resolution(int inner_order, int outer_order);
  const ContextTable *context_table(int order);
  double inverse_rule_bound(const ContextTable &table, const std::vector<double> &log_A,
                            const ZplotGaussHermiteRule &inner) const;
  Pruned prune_context(const ContextTable &table) const;
  void star(const ContextTable &table, const double *log_D, const std::vector<double> &weights,
            const ZplotGaussHermiteRule &rule, const double *z, int size, double allowance,
            const std::tuple<int, int, bool, double> *raw_key, ZplotStarResult *result);
  NormalizerRecord cdf_record(const ZplotStarResult &value, int contexts) const;
  const NormalizerRecord *raw_normalizer(int outer_order, int order);
  const Evaluation *evaluate(int outer_order, int inner_order, double allowance);
  bool inner_table(int outer_order, double allowance, InnerTable *out);
  double reference_at(int index, const std::vector<double> &reference_sd) const;
  const std::vector<bool> &selection_resolution(int inner_order, int denominator_order, int outer_order);
};


const std::vector<double> &ZplotContextController::gaussian_resolution(int inner_order, int outer_order)
{
  const auto key = std::make_pair(inner_order, outer_order);
  const auto found = resolutions_.find(key);
  if (found != resolutions_.end()) return found->second;
  // The bounds come from the main thread: in a parallel batch a missing one
  // suspends the draw here (nothing has changed in this step yet) together
  // with the orders the ladder usually asks for next.
  if (!main_thread_) {
    const bool need_factor = !factor_bounds_.count(inner_order);
    const bool need_context = !point_context_ && !context_bounds_.count(outer_order);
    if (need_factor || need_context) {
      factor_requests_ = {inner_order, 15};
      if (!point_context_) context_requests_ = {outer_order, 7, 15};
      throw NeedBound{};
    }
  }
  const int columns = geometry_.columns;
  if (!factor_bounds_.count(inner_order)) factor_bounds_[inner_order] = gaussian_bound(factor_beta_, inner_order);
  const std::vector<double> &bound = factor_bounds_[inner_order];
  const std::vector<double> row_bound = r_row_sums(bound.data(), k_, columns);
  std::vector<double> inner(k_), outer(k_, 0.0), total(k_);
  for (int row = 0; row < k_; ++row) inner[row] = height_ratio_[row] * row_bound[row];
  if (!point_context_) {
    if (!context_bounds_.count(outer_order)) context_bounds_[outer_order] = gaussian_bound(context_beta_, outer_order);
    const std::vector<double> &context = context_bounds_[outer_order];
    for (int row = 0; row < k_; ++row) outer[row] = height_ratio_[row] * context[row];
  }
  for (int row = 0; row < k_; ++row) total[row] = inner[row] + outer[row];
  return resolutions_[key] = {r_max(inner), r_max(outer), r_max(total)};
}

const ContextTable *ZplotContextController::context_table(int order)
{
  const auto found = tables_.find(order);
  if (found != tables_.end()) return &found->second;
  ContextTable table;
  const std::size_t K = static_cast<std::size_t>(k_);
  if (point_context_) {
    table.nodes.assign(point_rows_, 0.0);
    table.log_weights.assign(point_rows_, -std::log(static_cast<double>(point_rows_)));
    table.means = point_means_;
    table.rows = point_rows_;
  } else {
    const ZplotGaussHermiteRule &rule = quadrature_rule(s_, order);
    table.nodes = rule.nodes;
    table.log_weights = rule.log_weights;
    const int n = static_cast<int>(rule.nodes.size());
    table.rows = n;
    table.means.resize(static_cast<std::size_t>(n) * K);
    for (std::size_t j = 0; j < K; ++j)
      for (int i = 0; i < n; ++i)
        table.means[i + n * j] = rule.nodes[i] * context_factor_vector_[j] + d_.mean[j];
  }
  const int n = table.rows;
  const double log_context_mass = log_sum(table.log_weights);
  std::vector<double> lower_mean(K), upper_mean(K);
  for (std::size_t j = 0; j < K; ++j) {
    const double *column = table.means.data() + static_cast<std::size_t>(n) * j;
    if (point_context_) {
      lower_mean[j] = r_min(column, n);
      upper_mean[j] = r_max(column, n);
    } else {
      lower_mean[j] = r_pmin(column[0], column[n - 1]);
      upper_mean[j] = r_pmax(column[0], column[n - 1]);
    }
  }
  std::vector<double> standardized_distance(K);
  for (std::size_t j = 0; j < K; ++j) {
    standardized_distance[j] = r_pmax(std::fabs(grid_min_ * s_.sei[j] - upper_mean[j]),
      std::fabs(grid_max_ * s_.sei[j] - lower_mean[j])) /
      std::sqrt((1 - epsilon_actual_) * covariance_[j + K * j]);
  }
  std::vector<double> log_density_ratio(K);
  const double marginal_base = determinant_bound(1);
  std::vector<double> companion(K, 0.0);
  if (k_ > 1) {
    const double moment_bound = tilted_moment(k_ - 1);
    std::vector<double> mean_shift(K), inflated(K);
    for (std::size_t j = 0; j < K; ++j) {
      mean_shift[j] = epsilon_ / (1 - epsilon_) * standardized_distance[j];
      inflated[j] = mean_shift[j] / (1 - epsilon_);
    }
    for (std::size_t j = 0; j < K; ++j) if (!std::isfinite(inflated[j])) return nullptr;
    for (std::size_t j = 0; j < K; ++j) {
      if (beta_ + inflated[j] / (2 * std::sqrt(moment_bound)) > 1.0 / 4) return nullptr;
    }
    const double determinant = determinant_bound(k_ - 1);
    for (std::size_t j = 0; j < K; ++j) {
      companion[j] = determinant + mean_shift[j] * mean_shift[j] / (2 * (1 - epsilon_)) +
        beta_ * moment_bound + inflated[j] * std::sqrt(moment_bound);
    }
  }
  std::vector<double> error_terms(K);
  for (std::size_t j = 0; j < K; ++j) {
    const double marginal = marginal_base + beta_ * (standardized_distance[j] * standardized_distance[j]);
    log_density_ratio[j] = covariance_log_A_bound_ + marginal + companion[j];
    error_terms[j] = row_density_bound_[j] * std::expm1(log_density_ratio[j]);
  }
  const double covariance_error = r_mean(error_terms) * std::exp(log_context_mass);
  std::vector<double> scaled_precision(K);
  for (std::size_t j = 0; j < K; ++j) scaled_precision[j] = s_.sei[j] * std::sqrt(precision_diag_[j]);
  const double log_density_bound = r_log(max_omega()) - r_log(min_omega()) +
    r_log(r_mean(scaled_precision)) - .5 * std::log(2 * M_PI);
  if (!std::isfinite(covariance_error)) return nullptr;
  if (!point_context_) {
    std::vector<double> lower(n + 1), upper(n + 1);
    lower[0] = -kInf;
    for (int i = 0; i < n; ++i) lower[i + 1] = table.nodes[i];
    for (int i = 0; i < n; ++i) upper[i] = table.nodes[i];
    upper[n] = kInf;
    table.interval = normal_interval(lower, upper);
  }
  // order(log_weights, seq_along(log_weights)): ascending, ties by position.
  table.prune_order.resize(n);
  for (int i = 0; i < n; ++i) table.prune_order[i] = i;
  std::stable_sort(table.prune_order.begin(), table.prune_order.end(), [&table](int a, int b) {
    return table.log_weights[a] < table.log_weights[b];
  });
  table.mass = std::exp(log_context_mass);
  table.covariance_error = covariance_error;
  table.log_density_bound = log_density_bound;
  table.outer_tail = point_context_ ? 0 : std::exp(log_density_bound + gaussian_tail(table.nodes));
  return &(tables_[order] = std::move(table));
}

double ZplotContextController::inverse_rule_bound(const ContextTable &table, const std::vector<double> &log_A,
                                      const ZplotGaussHermiteRule &inner) const
{
  const double log_floor_inverse = -k_ * r_log(min_omega()) -
    integration_rank_ * log_sum(inner.log_weights);
  const std::size_t n = table.log_weights.size();
  std::vector<double> weighted(n);
  for (std::size_t i = 0; i < n; ++i) weighted[i] = table.log_weights[i] - log_A[i];
  double continuous = point_context_ ? log_sum(weighted) : log_floor_inverse;
  if (!point_context_) {
    bool increasing_weight = true, decreasing_weight = true;
    for (std::size_t b = 1; b < omega_.size(); ++b) {
      const double step = omega_[b] - omega_[b - 1];
      if (!(step <= 0)) increasing_weight = false;
      if (!(step >= 0)) decreasing_weight = false;
    }
    bool nonnegative = true, nonpositive = true;
    for (double value : context_factor_vector_) {
      const double signed_loading = s_.effect_sign * value;
      if (!(signed_loading >= 0)) nonnegative = false;
      if (!(signed_loading <= 0)) nonpositive = false;
    }
    const bool increasing_A = (increasing_weight && nonnegative) || (decreasing_weight && nonpositive);
    const bool decreasing_A = (decreasing_weight && nonnegative) || (increasing_weight && nonpositive);
    if (increasing_A || decreasing_A) {
      std::vector<double> rectangle(n + 1);
      if (increasing_A) {
        rectangle[0] = log_floor_inverse;
        for (std::size_t i = 0; i < n; ++i) rectangle[i + 1] = -log_A[i];
      } else {
        for (std::size_t i = 0; i < n; ++i) rectangle[i] = -log_A[i];
        rectangle[n] = log_floor_inverse;
      }
      std::vector<double> sum(n + 1);
      for (std::size_t i = 0; i <= n; ++i) sum[i] = table.interval[i] + rectangle[i];
      const double rectangle_bound = log_sum(sum);
      if (std::isfinite(rectangle_bound)) continuous = r_min2(continuous, rectangle_bound);
    }
  }
  return r_max2(continuous, log_sum(weighted));
}

Pruned ZplotContextController::prune_context(const ContextTable &table) const
{
  const int N = table.rows;
  Pruned out;
  double log_omitted = -kInf;
  std::vector<bool> removed(N, false);
  bool any_removed = false;
  if (!point_context_ && pruning_allowance_ > 0 && N > 1) {
    const double log_cap = r_log(pruning_allowance_) - table.log_density_bound;
    for (int c = 0; c + 1 < N; ++c) {
      const int index = table.prune_order[c];
      const double next_mass = log_sum({log_omitted, table.log_weights[index]});
      if (next_mass > log_cap) break;
      log_omitted = next_mass;
      removed[index] = true;
      any_removed = true;
    }
  }
  for (int i = 0; i < N; ++i) if (!any_removed || !removed[i]) out.keep.push_back(i);
  std::vector<double> kept(out.keep.size());
  for (std::size_t i = 0; i < out.keep.size(); ++i) kept[i] = table.log_weights[out.keep[i]];
  out.mass = std::exp(log_omitted);
  out.retained_mass = std::exp(log_sum(kept));
  out.error = (std::isinf(log_omitted) && log_omitted < 0) ? 0 :
    std::exp(table.log_density_bound + log_omitted);
  return out;
}

void ZplotContextController::star(const ContextTable &table, const double *log_D, const std::vector<double> &weights,
    const ZplotGaussHermiteRule &rule, const double *z, int size, double allowance,
    const std::tuple<int, int, bool, double> *raw_key, ZplotStarResult *result)
{
  ZplotStarCall call;
  call.means = table.means.data();
  call.contexts = table.rows;
  call.k = k_;
  call.covariance_lower = packed_.data();
  call.selection_se = s_.sei;
  call.omega = omega_.data();
  call.bins = s_.bins;
  call.lower = s_.z_lower;
  call.upper = s_.z_upper;
  call.effect_sign = s_.effect_sign;
  call.telescope = s_.telescope;
  call.log_normalizers = log_D;
  call.context_log_weights = weights.data();
  call.nodes = rule.nodes.data();
  call.log_weights = rule.log_weights.data();
  call.order = static_cast<int>(rule.nodes.size());
  call.z = z;
  call.size = size;
  call.allowance = allowance;
  call.bounded_cdf = cdf_enabled_;
  zplot_star_validate(call);
  // The R-owned raw-rule cache of the R controller: the same entries, keyed
  // and budgeted as there (object.size of the stored key/payload list).
  std::vector<double> saved;
  bool store = false;
  double stored_key_count = 0;
  if (raw_key != nullptr) {
    const auto found = raw_rules_.find(*raw_key);
    if (found != raw_rules_.end()) {
      call.cached_rules = found->second.data();
    } else {
      const double budget = std::min(raw_rule_budget_ - raw_rule_bytes_, 128.0 * 1024 * 1024);
      const double count = static_cast<double>(call.contexts) *
        (1.0 + static_cast<double>(k_) * call.order * (static_cast<double>(call.order) + 1.0));
      const double key_count = 8.0 + table.means.size() + packed_.size() + k_ + 3.0 * s_.bins +
        call.contexts + 2.0 * call.order;
      if ((count + key_count) * sizeof(double) + 1024 <= budget) {
        saved.assign(static_cast<std::size_t>(count), 0.0);
        call.saved_rules = saved.data();
        store = true;
        stored_key_count = key_count;
      }
    }
  }
  // R polls for interrupts between its .Call()s; serial batches do the same.
  if (main_thread_) robma_native::check_interrupt();
  zplot_star_evaluate(call, main_thread_ ? zplot_star_threads(call) : 1, result);
  if (store && result->available) {
    // object.size(list(key, payload)) for two long double vectors.
    const double bytes = 160.0 + 8.0 * (stored_key_count + saved.size());
    if (bytes <= raw_rule_budget_ - raw_rule_bytes_) {
      raw_rules_[*raw_key] = std::move(saved);
      raw_rule_bytes_ += bytes;
    }
  }
}

NormalizerRecord ZplotContextController::cdf_record(const ZplotStarResult &value, int contexts) const
{
  const std::vector<double> *bounds[] = {&value.cdf_factor_log_error, &value.cdf_log_error,
    &value.cdf_numerator_log_error};
  for (const std::vector<double> *bound : bounds) {
    bool valid = static_cast<int>(bound->size()) == contexts;
    for (double entry : *bound) if (!std::isfinite(entry) || !(entry >= 0)) valid = false;
    if (!valid) fail("Zplot CDF integration diagnostics are unavailable.");
  }
  if (!cdf_enabled_) {
    for (const std::vector<double> *bound : bounds)
      for (double entry : *bound)
        if (entry != 0) fail("Zplot CDF integration diagnostics are inconsistent.");
  }
  return {value.compact_log_normalizers, value.cdf_log_error, value.cdf_numerator_log_error};
}

// component_omission()
double component_omission(const ZplotStarResult &value)
{
  const double omitted = value.omitted_mass_error;
  if (!std::isfinite(omitted) || omitted < 0) fail("Zplot mixture omission diagnostics are unavailable.");
  return omitted;
}

const NormalizerRecord *ZplotContextController::raw_normalizer(int outer_order, int order)
{
  const auto key = std::make_tuple(outer_order, order, cdf_enabled_);
  const auto found = normalizers_.find(key);
  if (found != normalizers_.end()) return &found->second;
  const ContextTable *table = context_table(outer_order);
  if (table == nullptr) return nullptr;
  const ZplotGaussHermiteRule &quadrature = quadrature_rule(s_, order);
  const std::vector<double> weights(table->rows, -kInf);
  const double zero = 0.0;
  ZplotStarResult value;
  star(*table, nullptr, weights, quadrature, &zero, 1, 0.0, nullptr, &value);
  const double omitted_mass = component_omission(value);
  if (value.available || value.mass != 0 || value.components != 0 || value.compression_error != 0 ||
      omitted_mass != 0 || static_cast<int>(value.compact_log_normalizers.size()) != table->rows) {
    return nullptr;
  }
  for (double entry : value.compact_log_normalizers) if (!std::isfinite(entry)) return nullptr;
  return &(normalizers_[key] = cdf_record(value, table->rows));
}

const Evaluation *ZplotContextController::evaluate(int outer_order, int inner_order, double allowance)
{
  const int denominator_order = s_.inner_orders[denominator_index_ - 1];
  const auto key = std::make_tuple(outer_order, inner_order, denominator_order, cdf_enabled_,
    allowance_, pruning_allowance_);
  const auto found = results_.find(key);
  if (found != results_.end()) return &found->second;
  const ContextTable *table = context_table(outer_order);
  if (table == nullptr) return nullptr;
  const bool own_normalizer = reuse_normalizer_ && inner_order == denominator_order &&
    normalizers_.find(std::make_tuple(outer_order, denominator_order, cdf_enabled_)) == normalizers_.end();
  // R holds the denominator record by value: the numerator's raw A written
  // below under the same key must not change it.
  NormalizerRecord denominator;
  if (!own_normalizer) {
    const NormalizerRecord *found_denominator = raw_normalizer(outer_order, denominator_order);
    if (found_denominator == nullptr) return nullptr;
    denominator = *found_denominator;
  }
  std::vector<double> log_D;
  if (!own_normalizer) log_D = denominator.log_A;
  const ZplotGaussHermiteRule &inner = quadrature_rule(s_, inner_order);
  const ZplotGaussHermiteRule &denominator_rule = quadrature_rule(s_, denominator_order);
  const Pruned pruned = prune_context(*table);
  const std::vector<int> &keep = pruned.keep;
  std::vector<double> context_weights(table->rows, -kInf);
  for (int index : keep) context_weights[index] = table->log_weights[index];
  const auto raw_key = std::make_tuple(outer_order, inner_order, cdf_enabled_, pruning_allowance_);
  ZplotStarResult out;
  star(*table, own_normalizer ? nullptr : log_D.data(), context_weights, inner, s_.z, s_.size,
    allowance, cache_raw_rules_ ? &raw_key : nullptr, &out);
  if (!out.available) return nullptr;
  const double omitted_mass = component_omission(out);
  if (allowance == 0 && omitted_mass != 0) fail("Zplot mixture omission diagnostics are inconsistent.");
  const std::vector<double> &log_A = out.compact_log_normalizers;
  stop_if_not(static_cast<int>(log_A.size()) == table->rows, "length(log_A) == nrow(table$means)");
  bool finite = true;
  for (double entry : log_A) if (!std::isfinite(entry)) finite = false;
  stop_if_not(finite, "all(is.finite(log_A))");
  finite = true;
  for (double entry : out.density) if (!std::isfinite(entry)) finite = false;
  stop_if_not(finite, "all(is.finite(out$density))");
  bool nonnegative = true;
  for (double entry : out.density) if (!(entry >= 0)) nonnegative = false;
  stop_if_not(nonnegative, "all(out$density >= 0)");
  stop_if_not(std::isfinite(out.mass), "is.finite(out$mass)");
  stop_if_not(out.mass > 0, "out$mass > 0");
  stop_if_not(std::isfinite(out.compression_error), "is.finite(out$compression_error)");
  NormalizerRecord record = cdf_record(out, table->rows);
  normalizers_[std::make_tuple(outer_order, inner_order, cdf_enabled_)] = record;
  if (own_normalizer) {
    denominator = record;
    log_D = record.log_A;
  }
  const double log_tail_A = gaussian_tail(inner.nodes) + r_log(integration_rank_) +
    k_ * r_log(max_omega());
  std::vector<double> adjusted_D(log_D.size());
  for (std::size_t i = 0; i < log_D.size(); ++i) adjusted_D[i] = log_D[i] - denominator.mass_error[i];
  const double log_inverse = inverse_rule_bound(*table, adjusted_D, denominator_rule);
  const double numerator_mass_tail = std::exp(log_tail_A + log_inverse);
  std::vector<double> ratio(k_);
  for (int row = 0; row < k_; ++row) ratio[row] = s_.sei[row] / geometry_.residual_sd[row];
  const double numerator_tail = numerator_mass_tail * r_mean(ratio) / std::sqrt(2 * M_PI);
  const std::size_t kept = keep.size();
  std::vector<double> cdf_terms(kept), mass_terms(kept);
  for (std::size_t i = 0; i < kept; ++i) {
    const int index = keep[i];
    const double log_mass_ratio = log_A[index] - log_D[index];
    const double mass_log_error = record.mass_error[index] + denominator.mass_error[index];
    const double curve_log_error = record.numerator_error[index] + denominator.mass_error[index];
    const double log_weighted_ratio = table->log_weights[index] + log_mass_ratio;
    cdf_terms[i] = log_weighted_ratio + mass_log_error + log_abs_expm1(curve_log_error);
    mass_terms[i] = log_weighted_ratio + log_abs_expm1(mass_log_error);
  }
  Evaluation value;
  value.cdf_error = compact_density_bound_ * std::exp(log_sum(cdf_terms));
  value.mass_cdf_error = std::exp(log_sum(mass_terms));
  value.density = std::move(out.density);
  value.mass = out.mass;
  value.compression_error = out.compression_error;
  value.component_omission_mass = omitted_mass;
  value.denominator_error = denominator.mass_error;
  value.tail = numerator_tail;
  value.mass_tail = numerator_mass_tail;
  value.log_D = std::move(log_D);
  value.table = table;
  value.keep = keep;
  value.omitted_mass = pruned.mass;
  value.retained_mass = pruned.retained_mass;
  value.pruning_error = pruned.error;
  return &(results_[key] = std::move(value));
}

bool ZplotContextController::inner_table(int outer_order, double allowance, InnerTable *out)
{
  const Evaluation *fine = evaluate(outer_order, s_.inner_orders[inner_index_ - 1], allowance);
  const Evaluation *coarse = evaluate(outer_order, s_.inner_orders[inner_index_ - 2], allowance);
  if (fine == nullptr || coarse == nullptr) return false;
  stop_if_not(fine->keep == coarse->keep, "identical(fine$keep, coarse$keep)");
  bool same_D = fine->log_D.size() == coarse->log_D.size();
  for (std::size_t i = 0; same_D && i < fine->log_D.size(); ++i) {
    same_D = fine->log_D[i] == coarse->log_D[i] ||
      (std::isnan(fine->log_D[i]) && std::isnan(coarse->log_D[i]));
  }
  stop_if_not(same_D, "identical(fine$log_D, coarse$log_D)");
  bool same_error = fine->denominator_error.size() == coarse->denominator_error.size();
  for (std::size_t i = 0; same_error && i < fine->denominator_error.size(); ++i) {
    same_error = fine->denominator_error[i] == coarse->denominator_error[i] ||
      (std::isnan(fine->denominator_error[i]) && std::isnan(coarse->denominator_error[i]));
  }
  stop_if_not(same_error, "identical(fine$denominator_error, coarse$denominator_error)");
  const NormalizerRecord *found_previous = raw_normalizer(outer_order, s_.inner_orders[denominator_index_ - 2]);
  if (found_previous == nullptr) return false;
  const NormalizerRecord previous = *found_previous;
  const NormalizerRecord *previous_D = &previous;
  const std::vector<int> &keep = fine->keep;
  const std::size_t kept = keep.size();
  std::vector<double> log_weights(kept), log_difference(kept), adjacent_error(kept);
  bool any_adjacent = false;
  for (std::size_t i = 0; i < kept; ++i) {
    const int index = keep[i];
    log_weights[i] = fine->table->log_weights[index];
    log_difference[i] = log_abs_expm1(previous_D->log_A[index] - fine->log_D[index]);
    adjacent_error[i] = previous_D->mass_error[index] + fine->denominator_error[index];
    if (adjacent_error[i] > 0) any_adjacent = true;
  }
  if (any_adjacent) {
    for (std::size_t i = 0; i < kept; ++i) {
      log_difference[i] = log_add(adjacent_error[i] + log_difference[i], log_abs_expm1(adjacent_error[i]));
    }
  }
  std::vector<double> terms(kept);
  for (std::size_t i = 0; i < kept; ++i) terms[i] = log_weights[i] + log_difference[i];
  const double weighted_difference = std::exp(log_sum(terms));
  const ZplotGaussHermiteRule &denominator_rule = quadrature_rule(s_, s_.inner_orders[denominator_index_ - 1]);
  const double log_denominator_tail = gaussian_tail(denominator_rule.nodes) + r_log(integration_rank_) +
    k_ * r_log(max_omega());
  for (std::size_t i = 0; i < kept; ++i) {
    const int index = keep[i];
    terms[i] = log_weights[i] + log_denominator_tail - fine->log_D[index] + fine->denominator_error[index];
  }
  const double weighted_denominator_tail = std::exp(log_sum(terms));
  const double weighted_eta = 2 * weighted_difference + weighted_denominator_tail;
  std::vector<double> change(fine->density.size());
  for (std::size_t i = 0; i < change.size(); ++i) change[i] = std::fabs(fine->density[i] - coarse->density[i]);
  out->inner_change = 2 * r_max(change);
  out->normalizer = compact_density_bound_ * weighted_eta;
  out->compressed = 3 * fine->compression_error + 2 * coarse->compression_error;
  out->cdf = 3 * fine->cdf_error + 2 * coarse->cdf_error;
  out->error = out->inner_change + out->normalizer + out->compressed + out->cdf + fine->tail +
    fine->table->covariance_error + fine->pruning_error;
  const double mass_omission_error = 3 * fine->component_omission_mass + 2 * coarse->component_omission_mass;
  const double mass_error = std::fabs(fine->mass - fine->retained_mass) +
    2 * std::fabs(fine->mass - coarse->mass) + weighted_eta + fine->mass_tail + fine->omitted_mass +
    3 * fine->mass_cdf_error + 2 * coarse->mass_cdf_error + mass_omission_error;
  out->density = fine->density;
  out->mass_error = mass_error / fine->table->mass;
  out->mass_omission_error = mass_omission_error / fine->table->mass;
  out->covariance = fine->table->covariance_error;
  out->tail = fine->tail;
  out->pruning = fine->pruning_error;
  out->outer_tail = fine->table->outer_tail;
  out->fine = fine;
  return true;
}

double ZplotContextController::reference_at(int index, const std::vector<double> &reference_sd) const
{
  const std::size_t K = static_cast<std::size_t>(k_);
  const int m = point_rows_;
  std::vector<double> heights(static_cast<std::size_t>(m) * K);
  for (std::size_t j = 0; j < K; ++j) {
    const double x = s_.z[index] * s_.sei[j];
    for (int i = 0; i < m; ++i) {
      heights[i + m * j] = Rf_dnorm4(x, point_means_[i + m * j], reference_sd[j], 0);
    }
  }
  const std::vector<double> column_means = r_col_sums(heights.data(), m, k_, true);
  std::vector<double> scaled(K);
  for (std::size_t j = 0; j < K; ++j) scaled[j] = s_.sei[j] * column_means[j];
  return r_mean(scaled);
}

const std::vector<bool> &ZplotContextController::selection_resolution(int inner_order, int denominator_order,
                                                           int outer_order)
{
  const auto key = std::make_tuple(inner_order, denominator_order, outer_order);
  const auto found = selection_resolutions_.find(key);
  if (found != selection_resolutions_.end()) return found->second;
  const std::size_t K = static_cast<std::size_t>(k_);
  std::vector<double> outer_nodes;
  std::vector<double> means;
  int rows = 0;
  if (point_context_) {
    outer_nodes = {0.0};
    means = point_means_;
    rows = point_rows_;
  } else {
    outer_nodes = quadrature_rule(s_, outer_order).nodes_only;
    rows = static_cast<int>(outer_nodes.size());
    means.resize(static_cast<std::size_t>(rows) * K);
    for (std::size_t j = 0; j < K; ++j)
      for (int i = 0; i < rows; ++i)
        means[i + rows * j] = outer_nodes[i] * context_factor_vector_[j] + d_.mean[j];
  }
  const int columns = geometry_.columns + 1;
  std::vector<double> loading(K * columns);
  std::copy(geometry_.loading.begin(), geometry_.loading.end(), loading.begin());
  for (std::size_t j = 0; j < K; ++j) {
    loading[j + K * geometry_.columns] = point_context_ ? 0.0 : context_factor_vector_[j];
  }
  std::vector<int> inner_axes(geometry_.columns);
  for (int c = 0; c < geometry_.columns; ++c) inner_axes[c] = c;
  const auto check_inner = [&](int order) {
    const auto sweep = std::make_pair(order, outer_order);
    const auto hit = sweeps_.find(sweep);
    if (hit != sweeps_.end()) return hit->second;
    const bool value = transition_rule_resolves_rows(s_, omega_.data(), quadrature_rule(s_, order).nodes_only,
      means, rows, geometry_.residual_sd, loading, columns, inner_axes);
    sweeps_[sweep] = value;
    return value;
  };
  std::vector<bool> value(3);
  value[0] = check_inner(inner_order);
  value[1] = check_inner(denominator_order);
  value[2] = point_context_ || transition_rule_resolves(s_, omega_.data(), outer_nodes, d_.mean,
    geometry_.residual_sd, loading, columns, std::vector<int>{columns - 1});
  return selection_resolutions_[key] = value;
}

// The sorted refinement priority: sort(priority, decreasing = TRUE) is a
// stable decreasing sort that drops NaN.
std::vector<int> priority_order(const double (&priority)[7])
{
  std::vector<int> order;
  for (int i = 0; i < 7; ++i) if (!std::isnan(priority[i])) order.push_back(i);
  std::stable_sort(order.begin(), order.end(), [&priority](int a, int b) {
    return priority[a] > priority[b];
  });
  return order;
}

bool ZplotContextController::step()
{
  try {
    if (stage_ == stage_setup) {
      stage_ = stage_done;
      setup();
      if (stage_ == stage_done) return true;
    }
    if (stage_ == stage_resolve) {
      resolve();
      if (stage_ == stage_done) return true;
    }
    if (stage_ == stage_ladder) ladder();
    stage_ = stage_done;
    return true;
  } catch (const NeedBound &) {
    return false;
  } catch (const ZplotNativeFailure &failure) {
    out_->status = ZplotProjectionOutcome::failed;
    out_->message = failure.what();
  } catch (const std::bad_alloc &) {
    out_->status = ZplotProjectionOutcome::failed;
    out_->message = "Zplot context projection memory allocation failed.";
  } catch (const std::exception &exception) {
    out_->status = ZplotProjectionOutcome::failed;
    out_->message = exception.what();
  }
  stage_ = stage_done;
  return true;
}

bool ZplotContextController::supply_bounds()
{
  try {
    for (int order : factor_requests_) {
      if (!factor_bounds_.count(order)) factor_bounds_[order] = gaussian_bound(factor_beta_, order);
    }
    for (int order : context_requests_) {
      if (!context_bounds_.count(order)) context_bounds_[order] = gaussian_bound(context_beta_, order);
    }
  } catch (const ZplotNativeFailure &failure) {
    out_->status = ZplotProjectionOutcome::failed;
    out_->message = failure.what();
    stage_ = stage_done;
    return false;
  }
  factor_requests_.clear();
  context_requests_.clear();
  return true;
}

// Everything before the first Gaussian resolution bound. stage_done on any
// early return; stage_resolve when the ladder setup is complete.
void ZplotContextController::setup()
{
  const int k = k_;
  const std::size_t K = static_cast<std::size_t>(k);
  omega_.assign(d_.omega, d_.omega + s_.bins);
  for (double weight : omega_) if (weight <= 0) return;
  const bool cdf_enabled = true;
  cdf_enabled_ = cdf_enabled;
  // stopifnot(): only the draw-dependent conditions, in R's order. The
  // structural ones are checked once in R before the call.
  for (std::size_t j = 0; j < K; ++j) stop_if_not(std::isfinite(d_.mean[j]), "all(is.finite(mean))");
  const int factor_rows = d_.factor_rows;
  const std::size_t factor_size = static_cast<std::size_t>(factor_rows) * K;
  for (std::size_t i = 0; i < factor_size; ++i) {
    stop_if_not(std::isfinite(d_.context_factor[i]), "all(is.finite(context_factor))");
  }
  for (int i = 0; i < s_.size; ++i) stop_if_not(std::isfinite(s_.z[i]), "all(is.finite(z))");
  for (std::size_t j = 0; j < K; ++j) stop_if_not(std::isfinite(s_.sei[j]), "all(is.finite(sei))");
  for (std::size_t j = 0; j < K; ++j) stop_if_not(s_.sei[j] > 0, "all(sei > 0)");
  context_factor_.assign(d_.context_factor, d_.context_factor + factor_size);
  if (d_.point_means != nullptr) {
    bool valid = d_.point_rows >= 1;
    for (std::size_t i = 0; valid && i < static_cast<std::size_t>(d_.point_rows) * K; ++i) {
      if (!std::isfinite(d_.point_means[i])) valid = false;
    }
    for (std::size_t i = 0; valid && i < factor_size; ++i) if (context_factor_[i] != 0) valid = false;
    if (!valid) fail("Internal error: batched context projection inputs are invalid.");
    point_rows_ = d_.point_rows;
    point_means_.assign(d_.point_means, d_.point_means + static_cast<std::size_t>(point_rows_) * K);
  } else {
    point_rows_ = 1;
    point_means_.assign(d_.mean, d_.mean + K);
  }
  covariance_.assign(d_.covariance, d_.covariance + K * K);
  for (std::size_t i = 0; i < K; ++i) {
    for (std::size_t j = i + 1; j < K; ++j) {
      const double a = covariance_[i + K * j], b = covariance_[j + K * i];
      const bool same = a == b || (std::isnan(a) && std::isnan(b) && r_is_na(a) == r_is_na(b));
      if (!same) return;
    }
  }
  // chol(): dpotrf("U") on a copy with the lower triangle zeroed.
  cholesky_ = covariance_;
  for (std::size_t j = 0; j < K; ++j)
    for (std::size_t i = j + 1; i < K; ++i) cholesky_[i + K * j] = 0.0;
  {
    int info = 0;
    int n = k;
    F77_CALL(dpotrf)("U", &n, cholesky_.data(), &n, &info FCONE);
    if (info != 0) return;
  }
  // Uniform weights: the full Gaussian marginal (.zplot_normal_density_matrix
  // averaged over the point rows by colMeans).
  bool uniform = true;
  for (double weight : omega_) if (!(weight == omega_[0])) uniform = false;
  if (uniform) {
    const std::vector<double> factor_square_sums = [&] {
      std::vector<double> squares(factor_size);
      for (std::size_t i = 0; i < factor_size; ++i) squares[i] = context_factor_[i] * context_factor_[i];
      return r_col_sums(squares.data(), factor_rows, k);
    }();
    std::vector<double> sd(K), log_scale(K);
    for (std::size_t j = 0; j < K; ++j) {
      sd[j] = std::sqrt(covariance_[j + K * j] + factor_square_sums[j]);
      log_scale[j] = std::log(s_.sei[j]) - std::log(sd[j]);
    }
    const int m = point_rows_;
    std::vector<double> density(static_cast<std::size_t>(m) * s_.size);
    for (int row = 0; row < m; ++row) {
      for (int point = 0; point < s_.size; ++point) {
        double sum = 0;
        for (std::size_t j = 0; j < K; ++j) {
          sum += cpp_selnorm_affine_normal_pdf_log_scale(s_.z[point], s_.sei[j],
            point_means_[row + m * j], sd[j], log_scale[j]);
        }
        density[row + static_cast<std::size_t>(m) * point] = sum / k;
      }
    }
    out_->density = r_col_sums(density.data(), m, s_.size, true);
    out_->uniform = true;
    out_->status = ZplotProjectionOutcome::accepted;
    return;
  }
  // precision <- chol2inv(cholesky): dpotri("U"); only its diagonal is read.
  {
    std::vector<double> inverse(K * K, 0.0);
    for (std::size_t j = 0; j < K; ++j)
      for (std::size_t i = 0; i <= j; ++i) inverse[i + K * j] = cholesky_[i + K * j];
    int info = 0;
    int n = k;
    F77_CALL(dpotri)("U", &n, inverse.data(), &n, &info FCONE);
    if (info != 0) {
      if (info > 0) {
        fail("element (" + std::to_string(info) + ", " + std::to_string(info) +
          ") is zero, so the inverse cannot be computed");
      }
      fail("argument " + std::to_string(-info) + " of Lapack routine dpotri had invalid value");
    }
    precision_diag_.resize(K);
    for (std::size_t j = 0; j < K; ++j) precision_diag_[j] = inverse[j + K * j];
  }
  packed_.clear();
  for (std::size_t j = 0; j < K; ++j)
    for (std::size_t i = j; i < K; ++i) packed_.push_back(covariance_[i + K * j]);
  point_context_ = true;
  for (double value : context_factor_) if (!(value == 0)) point_context_ = false;
  if (s_.has_absolute_tolerance && !point_context_) {
    fail("An absolute context projection tolerance requires a point context.");
  }
  if (!point_context_ && factor_rows != 1) return;
  context_factor_vector_ = context_factor_;  // as.numeric(): one row when continuous
  if (!context_geometry(k, covariance_, packed_, s_.sei, cholesky_, &geometry_)) return;
  epsilon_actual_ = geometry_.relative_error;
  epsilon_ = epsilon_actual_ / (1 - epsilon_actual_);
  if (!std::isfinite(epsilon_) || epsilon_ < 0 || epsilon_ > 1.0 / 3) return;
  beta_ = epsilon_ / (2 * (1 - epsilon_));
  log_weight_ratio_ = r_log(max_omega()) - r_log(min_omega());
  covariance_log_A_bound_ = determinant_bound(k) + beta_ * tilted_moment(k);
  row_density_bound_.resize(K);
  for (std::size_t j = 0; j < K; ++j) {
    row_density_bound_[j] = max_omega() / min_omega() * s_.sei[j] / geometry_.residual_sd[j] /
      std::sqrt(2 * M_PI);
  }
  compact_density_bound_ = r_mean(row_density_bound_);
  if (geometry_.groups == 1) {
    integration_rank_ = 1;
  } else {
    std::vector<int> counts(geometry_.groups, 0);
    for (int type : geometry_.type_index) if (type >= 1 && type <= geometry_.groups) ++counts[type - 1];
    integration_rank_ = 1;
    for (int count : counts) if (count > 1) ++integration_rank_;
  }
  bool monotone = true;
  for (std::size_t b = 1; b < omega_.size(); ++b) if (!(omega_[b] - omega_[b - 1] <= 0)) monotone = false;
  reuse_normalizer_ = s_.telescope && monotone &&
    k * r_log(min_omega()) > r_log(DBL_MIN) / 2 &&
    k * r_log(max_omega()) < r_log(DBL_MAX) / 2;
  tolerance_ = s_.tolerance;
  grid_min_ = r_min(s_.z, s_.size);
  grid_max_ = r_max(s_.z, s_.size);
  const int columns = geometry_.columns;
  factor_beta_.resize(K * columns);
  for (int c = 0; c < columns; ++c) {
    for (std::size_t j = 0; j < K; ++j) {
      const double loading = geometry_.loading[j + K * c];
      factor_beta_[j + K * c] = loading * loading / (2 * geometry_.residual_variance[j]);
    }
  }
  std::vector<double> factor_squares(factor_size);
  for (std::size_t i = 0; i < factor_size; ++i) factor_squares[i] = context_factor_[i] * context_factor_[i];
  const std::vector<double> context_square_sums = r_col_sums(factor_squares.data(), factor_rows, k);
  std::vector<double> loading_squares(K * columns);
  for (std::size_t i = 0; i < K * columns; ++i) loading_squares[i] = geometry_.loading[i] * geometry_.loading[i];
  const std::vector<double> loading_square_sums = r_row_sums(loading_squares.data(), k, columns);
  context_beta_.resize(K);
  height_ratio_.resize(K);
  for (std::size_t j = 0; j < K; ++j) {
    context_beta_[j] = context_square_sums[j] / (2 * geometry_.residual_variance[j]);
    height_ratio_[j] = std::sqrt(1 + (loading_square_sums[j] + context_square_sums[j]) /
      geometry_.residual_variance[j]);
  }
  stage_ = stage_resolve;
}

// From the final-resolution check to the start of the ladder.
void ZplotContextController::resolve()
{
  const int k = k_;
  const std::size_t K = static_cast<std::size_t>(k);
  const int factor_rows = d_.factor_rows;
  const std::size_t factor_size = static_cast<std::size_t>(factor_rows) * K;
  const std::vector<double> &final_resolution = gaussian_resolution(s_.inner_orders[inner_count_ - 1],
    s_.context_orders[context_count_ - 1]);
  stage_ = stage_done;
  if (!std::isfinite(final_resolution[2]) || final_resolution[2] > tolerance_) return;
  denominator_index_ = std::min(3, inner_count_);
  context_index_ = inner_index_ = 2;
  compression_refinements_ = pruning_refinements_ = 0;
  cache_raw_rules_ = point_context_ && integration_rank_ > 1;
  // Reference peak.
  std::vector<double> reference_sd(K);
  {
    std::vector<double> squares(factor_size);
    for (std::size_t i = 0; i < factor_size; ++i) squares[i] = context_factor_[i] * context_factor_[i];
    const std::vector<double> sums = r_col_sums(squares.data(), factor_rows, k);
    for (std::size_t j = 0; j < K; ++j) reference_sd[j] = std::sqrt(covariance_[j + K * j] + sums[j]);
  }
  const std::vector<double> reference_mean = r_col_sums(point_means_.data(), point_rows_, k, true);
  std::vector<double> component_height(K);
  for (std::size_t j = 0; j < K; ++j) component_height[j] = s_.sei[j] / reference_sd[j];
  std::vector<int> reference_indices;
  bool all_finite = true, any_positive = false;
  for (double height : component_height) {
    if (!std::isfinite(height)) all_finite = false;
    if (height > 0) any_positive = true;
  }
  if (all_finite && any_positive) {
    const int highest = r_which_max(component_height);
    const double anchor = reference_mean[highest] / s_.sei[highest];
    std::vector<double> standardized(K);
    for (std::size_t j = 0; j < K; ++j) standardized[j] = reference_mean[j] / reference_sd[j];
    const double central = r_sum(standardized) / r_sum(component_height);
    if (std::isfinite(anchor) && std::isfinite(central)) {
      std::vector<double> distance(s_.size);
      for (int i = 0; i < s_.size; ++i) distance[i] = std::fabs(s_.z[i] - anchor);
      const int first = r_which_min(distance);
      for (int i = 0; i < s_.size; ++i) distance[i] = std::fabs(s_.z[i] - central);
      const int second = r_which_min(distance);
      if (first >= 0) reference_indices.push_back(first);
      if (second >= 0 && second != first) reference_indices.push_back(second);
    }
  }
  double reference_peak = kNA;
  if (!reference_indices.empty()) {
    std::vector<double> peaks;
    for (int index : reference_indices) peaks.push_back(reference_at(index, reference_sd));
    reference_peak = r_max(peaks);
  }
  if (!std::isfinite(reference_peak) || reference_peak <= 0) {
    std::vector<double> peaks(s_.size);
    for (int i = 0; i < s_.size; ++i) peaks[i] = reference_at(i, reference_sd);
    reference_peak = r_max(peaks);
  }
  reference_peak_ = reference_peak;
  allowance_ = tolerance_ * reference_peak / 128;
  pruning_allowance_ = point_context_ ? 0 : tolerance_ * reference_peak / 32;
  {
    const std::vector<bool> &resolved = selection_resolution(s_.inner_orders[inner_count_ - 1],
      s_.inner_orders[inner_count_ - 1], point_context_ ? 0 : s_.context_orders[context_count_ - 1]);
    if (!(resolved[0] && resolved[1] && resolved[2])) return;
  }
  stage_ = stage_ladder;
}

void ZplotContextController::ladder()
{
  const bool absolute = s_.has_absolute_tolerance;
  for (;;) {
    const int outer_order = point_context_ ? 0 : s_.context_orders[context_index_ - 1];
    const std::vector<double> resolution = gaussian_resolution(s_.inner_orders[inner_index_ - 1], outer_order);
    if (!std::isfinite(resolution[2]) || resolution[2] > tolerance_) {
      if (inner_index_ < inner_count_ && (point_context_ || resolution[0] >= resolution[1])) {
        ++inner_index_;
      } else if (!point_context_ && context_index_ < context_count_) {
        ++context_index_;
      } else if (inner_index_ < inner_count_) {
        ++inner_index_;
      } else return;
      continue;
    }
    const std::vector<bool> resolved = selection_resolution(s_.inner_orders[inner_index_ - 1],
      s_.inner_orders[denominator_index_ - 1], outer_order);
    if (!(resolved[0] && resolved[1] && resolved[2])) {
      if (!resolved[0] && inner_index_ < inner_count_) {
        ++inner_index_;
      } else if (!resolved[1] && denominator_index_ < inner_count_) {
        ++denominator_index_;
      } else if (!resolved[2] && !point_context_ && context_index_ < context_count_) {
        ++context_index_;
      } else return;
      continue;
    }
    InnerTable fine, coarse;
    const bool have_fine = inner_table(outer_order, allowance_, &fine);
    const bool have_coarse = point_context_ ? false :
      inner_table(s_.context_orders[context_index_ - 2], allowance_, &coarse);
    if (!have_fine || (!point_context_ && !have_coarse)) return;
    double outer_change = 0;
    if (!point_context_) {
      std::vector<double> change(fine.density.size());
      for (std::size_t i = 0; i < change.size(); ++i) change[i] = std::fabs(fine.density[i] - coarse.density[i]);
      outer_change = 2 * r_max(change);
    }
    const auto combine = [&](double InnerTable::*field) {
      return point_context_ ? fine.*field : 3 * (fine.*field) + 2 * (coarse.*field);
    };
    const double error = outer_change + combine(&InnerTable::error) + fine.outer_tail;
    const double peak_lower = r_max(fine.density) - fine.error;
    const double relative_error = peak_lower > 0 ? error / peak_lower : kInf;
    const InnerTable &limiting = (point_context_ || fine.mass_error >= coarse.mass_error) ? fine : coarse;
    const double mass_error = limiting.mass_error;
    const double mass_omission_error = limiting.mass_omission_error;
    const double normalizer_part = combine(&InnerTable::normalizer);
    const double compression_part = combine(&InnerTable::compressed);
    const double cdf_part = combine(&InnerTable::cdf);
    const double covariance_part = combine(&InnerTable::covariance);
    const double pruning_part = combine(&InnerTable::pruning);
    const double inner_part = combine(&InnerTable::inner_change) + combine(&InnerTable::tail);
    const bool density_accepted = !absolute ?
      (std::isfinite(relative_error) && relative_error <= tolerance_) :
      (std::isfinite(error) && error <= s_.absolute_tolerance);
    if (density_accepted && mass_error <= tolerance_) {
      out_->status = ZplotProjectionOutcome::accepted;
      out_->density = fine.density;
      out_->relative_error = relative_error;
      out_->mass_error = mass_error;
      out_->mass_omission_error = mass_omission_error;
      const double values[8] = {error, inner_part, normalizer_part, compression_part, cdf_part,
        covariance_part, pruning_part, outer_change + fine.outer_tail};
      std::copy(values, values + 8, out_->integration_error);
      out_->quadrature[0] = point_context_ ? 0 : s_.context_orders[context_index_ - 1];
      out_->quadrature[1] = s_.inner_orders[inner_index_ - 1];
      out_->quadrature[2] = s_.inner_orders[denominator_index_ - 1];
      out_->bounded_cdf = cdf_enabled_;
      return;
    }
    if (mass_error > tolerance_ && mass_omission_error >= mass_error / 2 && allowance_ > 0) {
      ++compression_refinements_;
      allowance_ = compression_refinements_ < 8 ? allowance_ / 2 : 0;
      continue;
    }
    const double routing_budget = !absolute ?
      tolerance_ * r_max2(r_max(fine.density), reference_peak_) : s_.absolute_tolerance;
    const double priority[7] = {inner_part, outer_change + fine.outer_tail, normalizer_part,
      compression_part, cdf_part, pruning_part, covariance_part};
    enum { inner = 0, outer, normalizer, compression, cdf, pruning, covariance };
    bool advanced = false;
    for (int action : priority_order(priority)) {
      if (action == normalizer && denominator_index_ < inner_count_) {
        ++denominator_index_;
      } else if (action == inner && inner_index_ < inner_count_) {
        ++inner_index_;
      } else if (action == outer && !point_context_ && context_index_ < context_count_) {
        ++context_index_;
      } else if (action == cdf && cdf_enabled_ && cdf_part > routing_budget / 4) {
        cdf_enabled_ = false;
      } else if (action == compression && allowance_ > 0 && compression_part > routing_budget / 4) {
        ++compression_refinements_;
        allowance_ = compression_refinements_ < 8 ? allowance_ / 2 : 0;
      } else if (action == pruning && pruning_allowance_ > 0 &&
          (pruning_part > routing_budget / 4 || fine.fine->omitted_mass > tolerance_ / 4)) {
        ++pruning_refinements_;
        pruning_allowance_ = pruning_refinements_ < 8 ? pruning_allowance_ / 2 : 0;
      } else if (action == covariance && covariance_part > routing_budget) {
        return;
      } else continue;
      advanced = true;
      break;
    }
    if (!advanced) return;
  }
}

}  // namespace

// ---------------------------------------------------------------------------
// Shared kernel and envelope cores
// ---------------------------------------------------------------------------

void zplot_star_validate(const ZplotStarCall &call)
{
  const auto finite_all = [](const double *values, std::size_t n) {
    for (std::size_t i = 0; i < n; ++i) if (!std::isfinite(values[i])) return false;
    return true;
  };
  if (call.contexts < 1 || call.k < 1 || call.bins < 1 || call.order < 1 || call.size < 1 ||
      !std::isfinite(call.allowance) || call.allowance < 0.0) {
    throw ZplotNativeFailure("Deterministic star projection dimensions or controls are invalid.");
  }
  if (!finite_all(call.means, static_cast<std::size_t>(call.contexts) * call.k) ||
      (call.log_normalizers != nullptr && !finite_all(call.log_normalizers, call.contexts)) ||
      !finite_all(call.nodes, call.order) || !finite_all(call.log_weights, call.order) ||
      !finite_all(call.z, call.size)) {
    throw ZplotNativeFailure("Deterministic star projection locations and quadrature values must be finite.");
  }
  for (int context = 0; context < call.contexts; ++context) {
    const double weight = call.context_log_weights[context];
    if (std::isnan(weight) || weight == kInf) {
      throw ZplotNativeFailure("Context log weights must be finite or negative infinity.");
    }
  }
  for (int row = 0; row < call.k; ++row) {
    if (!std::isfinite(call.selection_se[row]) || !(call.selection_se[row] > 0)) {
      throw ZplotNativeFailure("'selection_se' must be finite and positive.");
    }
  }
  for (int bin = 0; bin < call.bins; ++bin) {
    if (!std::isfinite(call.omega[bin]) || call.omega[bin] < 0) {
      throw ZplotNativeFailure("'omega' must contain finite non-negative values.");
    }
  }
  if (!selnorm_is_descending_step_partition(call.lower, call.upper, call.bins)) {
    throw ZplotNativeFailure("Selection bounds do not define a descending step partition.");
  }
}

int zplot_star_threads(const ZplotStarCall &call)
{
  int active_contexts = 0;
  for (int context = 0; context < call.contexts; ++context) {
    if (call.context_log_weights[context] != -kInf) ++active_contexts;
  }
  const int star_rows = std::max(call.contexts, call.k);
  const double star_work = static_cast<double>(call.k) * call.order * call.order *
    (0.4 * (call.contexts - active_contexts) + 10.0 * active_contexts);
  return robma_work_threads(star_rows, star_work / star_rows) > 1 ?
    robma_batch_threads(star_rows, 2) : 1;
}

void zplot_star_evaluate(const ZplotStarCall &call, int threads, ZplotStarResult *result)
{
  result->density.assign(call.size, 0.0);
  result->compact_log_normalizers.assign(call.contexts, kNA);
  result->cdf_factor_log_error.assign(call.contexts, 0.0);
  double mass = kNA, error = kNA, omitted_mass_error = 0.0;
  std::size_t components = 0;
  bool available = false;
  std::string message;
  try {
    CompactStepPartition partition;
    const bool compact_partition = call.bounded_cdf && call.telescope &&
      partition.prepare(call.omega, call.lower, call.upper, call.bins);
    available = cpp_selnorm_context_star_projection(call.means, call.contexts, call.covariance_lower,
      call.k, call.selection_se, compact_partition ? partition.omega.data() : call.omega,
      compact_partition ? static_cast<int>(partition.omega.size()) : call.bins,
      compact_partition ? partition.lower.data() : call.lower,
      compact_partition ? partition.upper.data() : call.upper, call.effect_sign, call.telescope,
      call.log_normalizers, call.context_log_weights, call.nodes, call.log_weights, call.order,
      call.z, call.size, call.allowance, result->density.data(), &mass, &error, &components,
      result->compact_log_normalizers.data(), call.bounded_cdf, result->cdf_factor_log_error.data(),
      &omitted_mass_error, threads, call.cached_rules, call.saved_rules);
  } catch (const std::exception &exception) {
    message = exception.what();
  }
  // The R entry reported the first 511 characters of the message, and an
  // exception without one as no failure.
  if (!message.empty()) {
    char buffer[512] = {};
    std::strncpy(buffer, message.c_str(), sizeof(buffer) - 1);
    throw ZplotNativeFailure(std::string("Deterministic star projection failed: ") + buffer + ".");
  }
  result->available = available;
  result->mass = mass;
  result->compression_error = error;
  result->omitted_mass_error = omitted_mass_error;
  result->components = static_cast<double>(components);
  const auto multiple_bound = [](double bound, int count) {
    return bound > 0 && count > 0 ? std::nextafter(bound * count, kInf) : 0.0;
  };
  result->cdf_log_error.resize(call.contexts);
  result->cdf_numerator_log_error.resize(call.contexts);
  for (int context = 0; context < call.contexts; ++context) {
    result->cdf_log_error[context] = multiple_bound(result->cdf_factor_log_error[context], call.k);
    result->cdf_numerator_log_error[context] = multiple_bound(result->cdf_factor_log_error[context], call.k - 1);
  }
}

void zplot_covariance_envelope_components(const double *covariance_lower, int k,
  const double *selection_se, bool upper, ZplotEnvelopeComponents *result)
{
  const std::size_t K = static_cast<std::size_t>(k);
  std::vector<double> sd(K), work(K * 4);
  std::vector<int> types(K);
  int groups = 0;
  bool available = false;
  try {
    available = cpp_selnorm_covariance_envelope_components(covariance_lower, k, selection_se, upper,
      sd.data(), work.data(), &groups, types.data());
  } catch (const std::bad_alloc &) {
    throw ZplotNativeFailure("Compact envelope allocation failed.");
  }
  result->available = available;
  if (!available) return;
  const int rank = groups == 1 ? 1 : groups + 1;
  result->groups = groups;
  result->rank = rank;
  result->residual_sd = sd;
  result->residual_variance.resize(K);
  result->loading.assign(work.begin(), work.begin() + K * rank);
  result->type_index.resize(K);
  for (std::size_t row = 0; row < K; ++row) {
    result->residual_variance[row] = static_cast<double>(static_cast<long double>(sd[row]) * sd[row]);
    result->type_index[row] = types[row] + 1;
  }
  const long double infinity = std::numeric_limits<long double>::infinity();
  const auto operation_radius = [infinity](long double value) {
    if (!std::isfinite(value)) return infinity;
    return std::max(std::nextafter(value, infinity) - value, value - std::nextafter(value, -infinity));
  };
  const auto add_radius = [infinity](long double left, long double right) {
    return std::nextafter(left + right, infinity);
  };
  result->covariance_residual.assign(K * K, 0.0);
  result->covariance_residual_radius.assign(K * K, 0.0);
  std::size_t position = 0;
  for (std::size_t column = 0; column < K; ++column) {
    for (std::size_t row = column; row < K; ++row) {
      long double represented = row == column ? static_cast<long double>(sd[row]) * sd[row] : 0.0L;
      long double error = row == column ? operation_radius(represented) : 0.0L;
      for (int factor = 0; factor < rank; ++factor) {
        const long double product = static_cast<long double>(result->loading[row + K * factor]) *
          result->loading[column + K * factor];
        error = add_radius(error, operation_radius(product));
        represented += product;
        error = add_radius(error, operation_radius(represented));
      }
      const long double exact_input_difference =
        static_cast<long double>(covariance_lower[position++]) - represented;
      error = add_radius(error, operation_radius(exact_input_difference));
      const double difference = static_cast<double>(exact_input_difference);
      const double cast_error = std::max(std::nextafter(difference, kInf) - difference,
        difference - std::nextafter(difference, -kInf));
      error = add_radius(error, static_cast<long double>(cast_error));
      const double returned_radius = std::nextafter(static_cast<double>(error), kInf);
      result->covariance_residual[row + K * column] = difference;
      result->covariance_residual[column + K * row] = difference;
      result->covariance_residual_radius[row + K * column] = returned_radius;
      result->covariance_residual_radius[column + K * row] = returned_radius;
    }
  }
}

void zplot_context_projection_batch(const ZplotProjectionSetup &setup,
  const std::vector<ZplotProjectionDraw> &draws, int threads, int chunk_rows,
  std::vector<ZplotProjectionOutcome> *outcomes)
{
  const int S = static_cast<int>(draws.size());
  outcomes->assign(S, ZplotProjectionOutcome());
  const bool serial = threads <= 1;
  std::vector<std::unique_ptr<ZplotContextController>> controllers(S);
  for (int s = 0; s < S; ++s) {
    controllers[s].reset(new ZplotContextController(setup, draws[s], serial, &(*outcomes)[s]));
  }
  if (serial) {
    for (int s = 0; s < S; ++s) {
      robma_native::check_interrupt();
      controllers[s]->step();
      controllers[s].reset();
    }
    return;
  }
  // Parallel rounds. A draw stops only where it needs Gaussian resolution
  // bounds, which the main thread computes between rounds in draw order; each
  // draw's computation is therefore the serial one at every thread count.
  std::vector<int> pending(S);
  for (int s = 0; s < S; ++s) pending[s] = s;
  while (!pending.empty()) {
    const int rows = static_cast<int>(pending.size());
    std::vector<char> finished(rows, 0);
    robma_parallel::for_rows(rows, threads, chunk_rows, [&](int index) {
      finished[index] = controllers[pending[index]]->step() ? 1 : 0;
    });
    std::vector<int> next;
    for (int index = 0; index < rows; ++index) {
      const int s = pending[index];
      if (!finished[index]) {
        robma_native::check_interrupt();
        if (controllers[s]->supply_bounds()) {
          next.push_back(s);
          continue;
        }
      }
      controllers[s].reset();
    }
    pending.swap(next);
  }
}
