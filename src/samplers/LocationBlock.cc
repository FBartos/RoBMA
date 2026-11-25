#include "LocationBlock.h"
#include "CoarseSurrogate.h"
#include "CorrectedSampler.h"
#include "StepRowNormalizer.h"
#include "LocationScale.h"
#include "../distributions/selnorm-jags-bounds.h"
#include "../selnorm/selnorm.h"

#include <BayesTools/factor_slice_kernel.h>

#include <distribution/ScalarDist.h>
#include <distribution/VectorDist.h>
#include <graph/AggNode.h>
#include <graph/DeterministicNode.h>
#include <graph/LogicalNode.h>
#include <graph/StochasticNode.h>
#include <module/ModuleError.h>
#include <rng/RNG.h>
#include <sampler/GraphView.h>
#include <sampler/MutableSampleMethod.h>
#include <sampler/MutableSampler.h>
#include <JRmath.h>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <mutex>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace jags {
namespace RoBMA {
namespace {

namespace slice = ::BayesTools::factor_slice;

using Vec = std::vector<double>;

constexpr char group_tag_name[] = "robma_location_group";
constexpr char latent_tag_name[] = "robma_location_latent";
constexpr char setup_tag_name[] = "robma_location_setup";
constexpr char sampler_name[] = "RoBMA::LocationBlock";
// The error of an untracked change of a likelihood mean (correct()); its
// bracketed identifier is the condition's stable class name.
constexpr char untracked_mean_message[] =
  "Location block sampler [RoBMA_location_untracked_mean]: the mean of this likelihood node changed in a sweep "
  "that did not move it, so the design the sampler verified does not hold at this state";
constexpr char tracked_mean_message[] =
  "Location block sampler [RoBMA_location_untracked_mean]: the mean of this likelihood node at the sweep's final "
  "state differs from the mean the sampler's verified design predicts, so the design does not hold at this state";
// The tolerance of that comparison: of the row's conditional SD, and of the
// size of its means (rounding).
constexpr double tracked_mean_tolerance = 1e-6;
constexpr double tracked_mean_rounding = 1e-12;
constexpr double log_two_pi = 1.8378770664093454836;
const double negative_infinity = -std::numeric_limits<double>::infinity();

// Partner and swap thresholds of the design (absolute correlations).
constexpr double partner_correlation = 0.3;
constexpr double swap_correlation = 0.5;
// Slice widths in proposal standard deviations.
constexpr double slice_width = 2.5;
// The damped Newton fits take a step unless it lowers the objective by more
// than this tolerance relative to its scale (1 + |f|): a step whose change is
// at the rounding level is always taken, so that the fits at the two ends of
// a move, whose inputs differ by rounding, take the same steps.
constexpr double fit_tolerance = 1e-10;
// The cross-family swap L2 (LocationBlockBiasMoves.h) rebuilds at its
// proposal what an attempt that starts there computes for the move back and
// is the identity when the two ends' stage-1 log ratios, or their log ratios
// of the proposal densities, are not inverses to this tolerance (log units):
// the bound of the frozen-state tests, about 1e4 times the residual of
// ordinary data (1e-13) and far below a missing Hastings term (at least 1e-3).
constexpr double swap_reciprocity_tolerance = 1e-9;
// L2 is the identity where its base point's working means are not
// representable for the swap: Q (LocationBlockBiasMoves.h,
// swap_working_criterion), the size of the working means' rounding in units
// of the rows' SDs weighted by their residuals, exceeds this tolerance. The
// rule above compares two ends that share their working means, while an
// attempt that starts at the proposal in a later sweep reaches the base point
// from other working means. Q is a representability criterion, not a bound
// on what re-read working means change (label_swap): it separates the
// regimes, at most 9.1e-11 on ordinary data (the real cases 2.3e-11, a
// 10,000-row fixture 9.1e-11) and at least 4e-6 at a row of 1000 / 1e-6;
// outcomes located far from zero relative to their SEs (from a median
// |y| / SE of about 1,000 on the tests' 10-row fixture) are skipped. A
// working mean's rounding is taken as swap_working_rounding eps
// (eps = 2^-52) times the size of its terms at the end's state.
constexpr double swap_working_tolerance = 1e-8;
constexpr double swap_working_rounding = 4.0;
// The refresh kernels of the active coefficients in selection labels (A3; see
// LocationBlock::refresh_variant): chosen by measurement.
constexpr int default_refresh = 1;

// ---- per-fit report and test hook

enum Component { flip_on, flip_off, refresh, study, swap, joint, ellipse_move, component_count };
const char *const component_names[component_count] = {
  "flip_on", "flip_off", "refresh", "study", "swap", "joint", "ellipse"};

// The bias-label components (LocationBlockBiasMoves.h): D, L0, L1, S, L2.
enum BiasComponent { bias_dormant, bias_label, bias_gaussian, bias_slope, bias_swap, bias_component_count };
const char *const bias_component_names[bias_component_count] = {
  "dormant", "label", "gaussian_family", "slope", "label_swap"};

// Where an update spends its work: child evaluations of the moves' acceptance
// ratios, of the Newton-Laplace fits, of the elliptical slice steps and of the
// latent updates (the exact correction has its own count), and, with the
// timing test switch, the wall time of each part.
enum Category { cat_prepare, cat_move, cat_fit, cat_ellipse, cat_latent, cat_correction,
                cat_dormant, cat_label, cat_gaussian, cat_slope, cat_swap, category_count };
const char *const category_names[category_count] = {
  "prepare", "move", "fit", "ellipse", "latent", "correction",
  "dormant", "label", "gaussian_family", "slope", "label_swap"};

struct Counters {
  std::atomic<std::uint64_t> samplers{0};
  std::atomic<std::uint64_t> sweeps{0};
  std::atomic<std::uint64_t> sweeps_selection{0};
  std::atomic<std::uint64_t> sweeps_evidence{0};
  std::atomic<std::uint64_t> corrections_accepted{0};
  std::atomic<std::uint64_t> corrections_selection_accepted{0};
  std::atomic<std::uint64_t> corrections_evidence_accepted{0};
  std::atomic<std::uint64_t> fallback_exact{0};
  std::atomic<std::uint64_t> surrogate_evaluations{0};
  std::atomic<std::uint64_t> exact_evaluations{0};
  std::atomic<std::uint64_t> correction_evaluations{0};
  std::atomic<std::uint64_t> fit_rule_fallback{0};
  std::atomic<std::uint64_t> ellipse_shrinks{0};
  // groups (per chain) whose refresh became the slice update, and of them those
  // of two or more coordinates (a composite component that reverses its
  // coordinates in reversed sweeps)
  std::atomic<std::uint64_t> slice_refresh_groups{0};
  std::atomic<std::uint64_t> slice_refresh_multi{0};
  // components skipped because a factorization of the working model failed
  // (a matrix not numerically positive definite)
  std::atomic<std::uint64_t> factorization_skips{0};
  // sweeps whose cross-family swap (L2) made no proposal because its base
  // point was not representable (LocationBlockBiasMoves.h, swap_unrepresentable:
  // M_0 or the working means' representability criterion Q), and after
  // adaptation the largest Q of the base points that pass that check (also
  // in attempts that end before a proposal; a Q above the tolerance at a
  // rebuilt end counts as swap_irreproducible)
  std::atomic<std::uint64_t> swap_unrepresentable{0};
  std::atomic<double> max_swap_working{0.0};
  // sweeps whose cross-family swap passed stage 1 and was then the identity
  // because its proposal end did not reproduce it (swap_irreproducible), and
  // after adaptation the largest reciprocity residuals of the attempts that
  // went on to stage 2: |stage 1 + reverse stage 1| and the sum of the two
  // ends' log ratios of the proposal densities
  std::atomic<std::uint64_t> swap_irreproducible{0};
  std::atomic<double> max_swap_stage1_residual{0.0};
  std::atomic<double> max_swap_full_residual{0.0};
  // after adaptation, the comparisons of an unmoved child's mean with its
  // entry mean that were not bitwise equal (correct()), and the largest such
  // discrepancy in the row's conditional SDs
  std::atomic<std::uint64_t> unmoved_mean_inexact{0};
  std::atomic<double> unmoved_mean_max{0.0};
  // after adaptation, the unmoved children whose graph mean changed and that
  // entered the exact correction with their own densities (correct()), and
  // the largest absolute change of such a child's log density
  std::atomic<std::uint64_t> unmoved_included{0};
  std::atomic<double> unmoved_density_max{0.0};
  std::atomic<std::uint64_t> proposed[component_count];
  std::atomic<std::uint64_t> stage0[component_count];
  std::atomic<std::uint64_t> stage1[component_count];
  std::atomic<std::uint64_t> accepted[component_count];
  std::atomic<std::uint64_t> bias_proposed[bias_component_count];
  std::atomic<std::uint64_t> bias_stage1[bias_component_count];
  std::atomic<std::uint64_t> bias_accepted[bias_component_count];
  std::atomic<std::uint64_t> evaluations[category_count];
  std::atomic<std::uint64_t> nanoseconds[category_count];
  std::atomic<double> max_correction{0.0};
  std::atomic<double> max_correction_gaussian{0.0};
  std::atomic<double> max_endpoint{0.0};
  std::atomic<double> max_endpoint_change{0.0};
  std::atomic<double> max_stage2_gaussian{0.0};

  Counters()
  {
    for (int i = 0; i < component_count; ++i) {
      proposed[i] = 0;
      stage0[i] = 0;
      stage1[i] = 0;
      accepted[i] = 0;
    }
    for (int i = 0; i < bias_component_count; ++i) {
      bias_proposed[i] = 0;
      bias_stage1[i] = 0;
      bias_accepted[i] = 0;
    }
    for (int i = 0; i < category_count; ++i) {
      evaluations[i] = 0;
      nanoseconds[i] = 0;
    }
  }
};

// Never destroyed: samplers of models that R frees late may still count.
Counters &counters()
{
  static Counters *store = new Counters();
  return *store;
}

void count(std::atomic<std::uint64_t> &counter, std::uint64_t n = 1)
{
  counter.fetch_add(n, std::memory_order_relaxed);
}

void record_max(std::atomic<double> &target, double value)
{
  if (!std::isfinite(value)) value = std::numeric_limits<double>::max();
  double current = target.load(std::memory_order_relaxed);
  while (value > current &&
         !target.compare_exchange_weak(current, value, std::memory_order_relaxed)) {}
}

std::atomic<int> test_mutant(0);
std::atomic<bool> test_probe(false);
// Test switches, frozen into the samplers built from now on (LocationBlockReport.h).
std::atomic<int> test_refresh(0);
std::atomic<int> test_fit_rule(0);
std::atomic<bool> test_timing(false);
// the schedule of L2 (0 the default, k: probability 1 / k) and its per-attempt trace
std::atomic<int> test_swap(0);
std::atomic<int> test_trace(0);
// the coin of the gate flips' and location swaps' first stage and of L0 (0 the
// default, a state-independent coin per attempt with probability
// two_state_coin or label_coin; 1 always the alternative: Barker, Gibbs; 2 never)
std::atomic<int> test_coin(0);
// The state-independent coins (spec: periodicity): Barker's first stage of the
// gate flips and location swaps with probability 1/4 per attempt, and the
// plain Gibbs draw of L0 with probability 1/8 (Liu's step otherwise). L2 needs
// no coin: it is attempted on a state-independent quarter of the sweeps.
constexpr double two_state_coin = 0.25;
constexpr double label_coin = 0.125;
// The probability of attempting the bias label's cross-family swap L2 in a sweep
// (wp3-spec.md, schedule): one unless L2 costs more than 5% of the sampling
// sweep on the Assink full-V moderator fit or the 18-moderator tutorial fit,
// then the largest of 1/2 and 1/4 within 5% on both. Measured at one (wp3
// profiles, L2's Laplace fits counted as L2): 15.7% and 10.6%; at 1/2: 8.5%
// and 5.6%; at 1/4: 4.4% and 2.9%.
constexpr double default_swap_probability = 0.25;

// The evidence record of the location gates (location_block_evidence): per
// gate slot and chain, over the updates after adaptation, the updates, the
// transitions of the gate between each update's entry state and the state it
// retains (on: off to on; a sweep the exact correction rejects changes
// nothing, a swap changes both of its gates), and of the inner stage (the
// sweep's components, before its exact correction) the flips proposed and
// accepted in each direction, the proposals the free screen removed, and for
// the proposals that reached stage 1 the count, sum and maximum of the log
// collapsed odds of the other state against the current one. Never
// destroyed, like the counters.
struct GateEvidence {
  double updates = 0, switched_on = 0, switched_off = 0;
  double on_proposed = 0, on_accepted = 0, off_proposed = 0, off_accepted = 0, screened_out = 0;
  double odds_count = 0, odds_sum = 0, odds_max = -std::numeric_limits<double>::infinity();

  void add(GateEvidence const &other)
  {
    updates += other.updates;
    switched_on += other.switched_on;
    switched_off += other.switched_off;
    on_proposed += other.on_proposed;
    on_accepted += other.on_accepted;
    off_proposed += other.off_proposed;
    off_accepted += other.off_accepted;
    screened_out += other.screened_out;
    odds_count += other.odds_count;
    odds_sum += other.odds_sum;
    odds_max = std::max(odds_max, other.odds_max);
  }
};

struct GateEvidenceStore {
  std::mutex mutex;
  std::map<std::pair<unsigned, unsigned>, GateEvidence> records;
};

GateEvidenceStore &gate_evidence_store()
{
  static GateEvidenceStore *store = new GateEvidenceStore();
  return *store;
}

struct ProbeStore {
  std::mutex mutex;
  std::vector<double> records;
};

// Never destroyed, like the sampler registry: samplers may run late.
ProbeStore &probe_store()
{
  static ProbeStore *store = new ProbeStore();
  return *store;
}

// The per-attempt record of L2 (test switch trace), never destroyed either.
ProbeStore &trace_store()
{
  static ProbeStore *store = new ProbeStore();
  return *store;
}

// The sweep sequences (test switch trace, bit 2), never destroyed either.
ProbeStore &sequence_store()
{
  static ProbeStore *store = new ProbeStore();
  return *store;
}

// The outcome of the last factory call that found a declaration (test hook):
// "claimed" or the first check that declined. Factories run on the main thread.
std::string &decision()
{
  static std::string *store = new std::string("none");
  return *store;
}

bool fail(char const *reason)
{
  decision() = reason;
  return false;
}

std::nullptr_t fail_null(char const *reason)
{
  decision() = reason;
  return nullptr;
}

// ---- small dense algebra (row-major)

// Lower Cholesky factor of a symmetric positive definite n x n matrix (its
// lower triangle is read), in place, upper triangle zeroed; false when a
// pivot is not positive and finite.
bool cholesky(Vec &a, std::size_t n)
{
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j <= i; ++j) {
      long double v = a[i * n + j];
      for (std::size_t k = 0; k < j; ++k) v -= static_cast<long double>(a[i * n + k]) * a[j * n + k];
      if (i == j) {
        if (!(v > 0) || !std::isfinite(static_cast<double>(v))) return false;
        a[i * n + i] = std::sqrt(static_cast<double>(v));
      } else {
        a[i * n + j] = static_cast<double>(v / a[j * n + j]);
      }
    }
    for (std::size_t j = i + 1; j < n; ++j) a[i * n + j] = 0.0;
  }
  return true;
}

// L x = b and L' x = b in place.
void solve_lower(Vec const &L, std::size_t n, double *b)
{
  for (std::size_t i = 0; i < n; ++i) {
    long double v = b[i];
    for (std::size_t k = 0; k < i; ++k) v -= static_cast<long double>(L[i * n + k]) * b[k];
    b[i] = static_cast<double>(v / L[i * n + i]);
  }
}

void solve_lower_transpose(Vec const &L, std::size_t n, double *b)
{
  for (std::size_t i = n; i-- > 0;) {
    long double v = b[i];
    for (std::size_t k = i + 1; k < n; ++k) v -= static_cast<long double>(L[k * n + i]) * b[k];
    b[i] = static_cast<double>(v / L[i * n + i]);
  }
}

Vec inverse_from_cholesky(Vec const &L, std::size_t n)
{
  Vec out(n * n, 0.0), column(n);
  for (std::size_t j = 0; j < n; ++j) {
    std::fill(column.begin(), column.end(), 0.0);
    column[j] = 1.0;
    solve_lower(L, n, column.data());
    solve_lower_transpose(L, n, column.data());
    for (std::size_t i = 0; i < n; ++i) out[i * n + j] = column[i];
  }
  // exact symmetry
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < i; ++j) {
      const double v = 0.5 * (out[i * n + j] + out[j * n + i]);
      out[i * n + j] = v;
      out[j * n + i] = v;
    }
  }
  return out;
}

double log_det_cholesky(Vec const &L, std::size_t n)
{
  double out = 0.0;
  for (std::size_t i = 0; i < n; ++i) out += 2.0 * std::log(L[i * n + i]);
  return out;
}

// log N(x; m, P^-1) from the lower Cholesky factor of the precision P.
double log_normal_precision(double const *x, Vec const &m, Vec const &Pchol, double log_det_P)
{
  const std::size_t n = m.size();
  // ||L'(x - m)||^2
  long double quadratic = 0;
  for (std::size_t j = 0; j < n; ++j) {
    long double v = 0;
    for (std::size_t i = j; i < n; ++i) v += static_cast<long double>(Pchol[i * n + j]) * (x[i] - m[i]);
    quadratic += v * v;
  }
  return 0.5 * log_det_P - 0.5 * static_cast<double>(n) * log_two_pi - 0.5 * static_cast<double>(quadratic);
}

// A draw from N(m, P^-1) given the lower Cholesky factor of P.
void draw_precision(Vec const &m, Vec const &Pchol, RNG *rng, double *out)
{
  const std::size_t n = m.size();
  Vec z(n);
  for (std::size_t i = 0; i < n; ++i) z[i] = rng->normal();
  solve_lower_transpose(Pchol, n, z.data());
  for (std::size_t i = 0; i < n; ++i) out[i] = m[i] + z[i];
}

double log1pexp(double x)
{
  return x > 0 ? x + std::log1p(std::exp(-x)) : std::log1p(std::exp(x));
}

// The largest difference of two vectors of graph values (the restoration
// probe, LocationBlock::restoration_probe): equal values and two NaN count as
// zero; a NaN or infinite value against another value, or vectors of
// different lengths, as an infinite difference (finiteness first: std::max
// would keep a NaN difference out).
double endpoint_difference(Vec const &a, Vec const &b)
{
  const double infinite = std::numeric_limits<double>::infinity();
  if (a.size() != b.size()) return infinite;
  double out = 0.0;
  for (std::size_t i = 0; i < a.size(); ++i) {
    if (a[i] == b[i] || (std::isnan(a[i]) && std::isnan(b[i]))) continue;
    out = std::isfinite(a[i]) && std::isfinite(b[i]) ? std::max(out, std::fabs(a[i] - b[i])) : infinite;
  }
  return out;
}

// The gate of the target cell of the cross-family swap L2: the log
// probabilities of the slab (on) and of the spike (off) from the cell's two
// log evidences, both from their difference d = on - off as -log(1 + e^-d)
// and -log(1 + e^d) (stable for either sign). No total of the two is formed:
// at evidences of size 1e9 its rounding drops the smaller cell, which the
// draw then never takes while a score from the total still gives it mass. A
// side with evidence -Inf (or none) has probability zero. The draw and the
// proposal density use these same two numbers.
void swap_gate_log_probabilities(double log_e_off, double log_e_on, double *log_on, double *log_off)
{
  const bool on_mass = log_e_on > negative_infinity, off_mass = log_e_off > negative_infinity;
  if (!on_mass || !off_mass) {
    *log_on = on_mass ? 0.0 : negative_infinity;
    *log_off = off_mass ? 0.0 : negative_infinity;
    return;
  }
  const double d = log_e_on - log_e_off;
  *log_on = -log1pexp(-d);
  *log_off = -log1pexp(d);
}

// Eigen-decomposition of a small symmetric n x n matrix by cyclic Jacobi
// rotations: values, and vectors as the columns of a row-major matrix.
void symmetric_eigen(Vec a, std::size_t n, Vec &values, Vec &vectors)
{
  vectors.assign(n * n, 0.0);
  for (std::size_t i = 0; i < n; ++i) vectors[i * n + i] = 1.0;
  for (unsigned sweep = 0; sweep < 100; ++sweep) {
    double off = 0.0, diagonal = 0.0;
    for (std::size_t p = 0; p < n; ++p) {
      diagonal += a[p * n + p] * a[p * n + p];
      for (std::size_t q = p + 1; q < n; ++q) off += a[p * n + q] * a[p * n + q];
    }
    if (!(off > 1e-30 * diagonal)) break;
    for (std::size_t p = 0; p < n; ++p) {
      for (std::size_t q = p + 1; q < n; ++q) {
        const double apq = a[p * n + q];
        if (apq == 0.0) continue;
        const double theta = (a[q * n + q] - a[p * n + p]) / (2.0 * apq);
        const double t = (theta >= 0 ? 1.0 : -1.0) / (std::fabs(theta) + std::sqrt(theta * theta + 1.0));
        const double c = 1.0 / std::sqrt(t * t + 1.0), s = t * c;
        for (std::size_t k = 0; k < n; ++k) {
          const double akp = a[k * n + p], akq = a[k * n + q];
          a[k * n + p] = c * akp - s * akq;
          a[k * n + q] = s * akp + c * akq;
        }
        for (std::size_t k = 0; k < n; ++k) {
          const double apk = a[p * n + k], aqk = a[q * n + k];
          a[p * n + k] = c * apk - s * aqk;
          a[q * n + k] = s * apk + c * aqk;
        }
        for (std::size_t k = 0; k < n; ++k) {
          const double vkp = vectors[k * n + p], vkq = vectors[k * n + q];
          vectors[k * n + p] = c * vkp - s * vkq;
          vectors[k * n + q] = s * vkp + c * vkq;
        }
      }
    }
  }
  values.resize(n);
  for (std::size_t i = 0; i < n; ++i) values[i] = a[i * n + i];
}

// The precision of a tilted law (spec: Selection tilt): H = P - ell'' with the
// eigenvalues of P^-1/2 H P^-1/2 floored at 0.25 (max(0.25 P, H) for a
// scalar). Pchol is the lower Cholesky factor of P.
bool floor_precision(Vec const &P, Vec const &Pchol, Vec const &hessian, std::size_t n, Vec &out)
{
  Vec H(n * n);
  for (std::size_t i = 0; i < n * n; ++i) H[i] = P[i] - hessian[i];
  if (n == 1) {
    out.assign(1, std::max(0.25 * P[0], H[0]));
    return std::isfinite(out[0]) && out[0] > 0;
  }
  // A = L^-1 H L^-T
  Vec A(H);
  for (std::size_t j = 0; j < n; ++j) {
    Vec column(n);
    for (std::size_t i = 0; i < n; ++i) column[i] = A[i * n + j];
    solve_lower(Pchol, n, column.data());
    for (std::size_t i = 0; i < n; ++i) A[i * n + j] = column[i];
  }
  for (std::size_t i = 0; i < n; ++i) {
    Vec row(A.begin() + i * n, A.begin() + (i + 1) * n);
    solve_lower(Pchol, n, row.data());
    for (std::size_t j = 0; j < n; ++j) A[i * n + j] = row[j];
  }
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < i; ++j) {
      const double v = 0.5 * (A[i * n + j] + A[j * n + i]);
      A[i * n + j] = v;
      A[j * n + i] = v;
    }
  }
  Vec values, vectors;
  symmetric_eigen(A, n, values, vectors);
  // out = L V diag(max(values, 0.25)) V' L'
  Vec M(n * n, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < n; ++j) {
      long double v = 0;
      for (std::size_t k = 0; k < n; ++k) {
        if (!std::isfinite(values[k])) return false;
        v += static_cast<long double>(vectors[i * n + k]) * std::max(values[k], 0.25) * vectors[j * n + k];
      }
      M[i * n + j] = static_cast<double>(v);
    }
  }
  out.assign(n * n, 0.0);
  Vec LM(n * n, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < n; ++j) {
      long double v = 0;
      for (std::size_t k = 0; k <= i; ++k) v += static_cast<long double>(Pchol[i * n + k]) * M[k * n + j];
      LM[i * n + j] = static_cast<double>(v);
    }
  }
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < n; ++j) {
      long double v = 0;
      for (std::size_t k = 0; k <= j; ++k) v += static_cast<long double>(LM[i * n + k]) * Pchol[j * n + k];
      out[i * n + j] = static_cast<double>(v);
    }
  }
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < i; ++j) {
      const double v = 0.5 * (out[i * n + j] + out[j * n + i]);
      out[i * n + j] = v;
      out[j * n + i] = v;
    }
  }
  return true;
}

// The row-marginal selected-normal normalizer of one row, A(mu) =
// sum_b omega_b P(lower_b < sign y / se < upper_b), y ~ N(mu, sd^2), and the
// first two derivatives of log A in mu (spec: the free screen), by the
// tail-safe evaluation of the step rows (StepRowNormalizer), so that a row far
// in a tail keeps its screen.
bool row_log_normalizer(double mu, double sd, double se, int sign, double const *omega,
                        double const *lower, double const *upper, int bins,
                        double *log_a, double *d1, double *d2)
{
  StepRowLaw law;
  law.sd = sd;
  law.se = se;
  law.effect_sign = sign;
  law.bins = bins;
  law.omega = omega;
  law.lower = lower;
  law.upper = upper;
  return step_row_log_normalizer(law, mu, log_a, d1, d2);
}

// ---- structure of the declared model

enum class ChildType { normal, step_switch, cluster, factor, mv };

bool binary(Node const *node, const char *syntax)
{
  auto const *logical = dynamic_cast<LogicalNode const *>(node);
  return logical != nullptr && node->parents().size() == 2 && logical->deparse({"a", "b"}) == syntax;
}

bool fixed_value(Node const *node, double value)
{
  if (!node->isFixed()) return false;
  for (unsigned long i = 0; i < node->length(); ++i) {
    if (node->value(0)[i] != value) return false;
  }
  return true;
}

struct GroupSpec {
  Node const *slab = nullptr;       // the slab as referenced: a node, or an aggregate of scalars
  std::vector<StochasticNode *> slab_nodes;  // its stochastic coordinates (members)
  StochasticNode *gate = nullptr;   // null: ungated (always on)
  Node const *physical = nullptr;   // the physical coefficient node
  double label = 0;                 // the slab's label of the gate
  unsigned d = 0;
  unsigned coord0 = 0;              // first coordinate in the coefficient vector
  unsigned gate_pos = 0, slab_pos = 0;  // positions in the view's value vector
  Vec m0, B, Bchol;                 // slab prior mean, precision and its factor
  double log_det_B = 0;
  double log_on = 0, log_off = 0;   // log prior probabilities of the slab and spike labels
  std::vector<Node const *> chart;  // comparisons and products of the chart
};

struct ChildSpec {
  StochasticNode *node = nullptr;
  ChildType type = ChildType::normal;
  ChildKind kind = ChildKind::other;  // the surrogate evaluator, if any
  unsigned n = 0, row0 = 0, woff = 0;
  int latent = -1;
  std::vector<unsigned long> lengths;
};

// The bias members (LocationBlockBias.h), when the block owns the bias label.
struct BiasSpec;

struct Shared {
  std::vector<GroupSpec> groups;
  std::vector<StochasticNode *> latents;
  std::vector<unsigned> latent_pos;
  std::vector<ChildSpec> children;
  std::vector<std::vector<unsigned>> latent_children;
  std::vector<std::vector<unsigned>> coord_children;
  std::vector<unsigned> coord_group;
  unsigned p = 0, rows = 0, wsize = 0, length = 0;
  Vec X;                                   // rows x p, row-major
  // rows x p: the size of the means at each entry's probe (column_check),
  // whose rounding 4 eps size bounds the entry's error per unit
  Vec X_size;
  std::unique_ptr<GraphView> latent_view;  // latents only, for the per-update probe
  SelNormCoarseSettings settings;
  SelNormCoarseSettings fit_settings;      // the coarsest rung of each schedule (Newton fits)
  bool corrected = false;                  // the surrogate is configured
  bool any_surrogate = false;              // some child has a surrogate evaluator
  bool coarse_fit = false;                 // some cluster or factor child has a coarse rung
  double min_acceptance = 0.0;
  int mutant = 0;
  bool probe = false;                      // test hook: dry-run records of the Gaussian moves
  int fit_rule = 0;                        // test hook: 0 coarse with safeguard, 1 coarse, 2 sampler's
  int refresh = 0;                         // test hook: the refresh variant (0 the default)
  bool timing = false;                     // test hook: wall time of the update's parts
  double swap_probability = 1.0;            // the schedule of L2
  int trace = 0;                           // test hook: 1 per-attempt record of L2, 2 sweep sequences (bits)
  int coin = 0;                            // test hook: the coin of the two-state moves and L0 (test_coin)
  std::shared_ptr<BiasSpec> bias;          // the bias indicator and slope sources, or null
};

// The stochastic coordinates of a slab reference: the node itself, or the
// scalar nodes of an aggregate (a factor emitted as independent scalars).
bool resolve_slab(Node const *slab, std::vector<StochasticNode *> &out)
{
  out.clear();
  if (auto *node = dynamic_cast<StochasticNode *>(const_cast<Node *>(slab))) {
    out.push_back(node);
    return true;
  }
  auto const *aggregate = dynamic_cast<AggNode const *>(slab);
  if (aggregate == nullptr) return false;
  for (std::size_t i = 0; i < aggregate->parents().size(); ++i) {
    auto *node = dynamic_cast<StochasticNode *>(const_cast<Node *>(aggregate->parents()[i]));
    if (node == nullptr || node->length() != 1 || aggregate->offsets()[i] != 0) return false;
    out.push_back(node);
  }
  return out.size() == slab->length();
}

bool same_slab(Node const *a, Node const *b)
{
  if (a == b) return true;
  auto const *x = dynamic_cast<AggNode const *>(a);
  auto const *y = dynamic_cast<AggNode const *>(b);
  return x != nullptr && y != nullptr && x->parents() == y->parents() && x->offsets() == y->offsets();
}

// The slab law with fixed parameters, unbounded: one scalar dnorm, one dmnorm,
// or independent scalar dnorm coordinates.
bool read_slab(GroupSpec &g)
{
  const unsigned long d = g.slab->length();
  g.d = static_cast<unsigned>(d);
  g.m0.assign(d, 0.0);
  g.B.assign(d * d, 0.0);
  if (g.slab_nodes.size() == 1 && g.slab_nodes[0]->distribution()->name() == "dmnorm") {
    StochasticNode const *slab = g.slab_nodes[0];
    auto const &par = slab->parents();
    if (slab->lowerBound() != nullptr || slab->upperBound() != nullptr || par.size() != 2 ||
        !par[0]->isFixed() || !par[1]->isFixed() || par[0]->length() != d || par[1]->length() != d * d) {
      return false;
    }
    g.m0.assign(par[0]->value(0), par[0]->value(0) + d);
    // JAGS stores the precision column-major; it must be symmetric.
    double const *P = par[1]->value(0);
    for (unsigned long i = 0; i < d; ++i) {
      for (unsigned long j = 0; j < d; ++j) {
        if (P[i + d * j] != P[j + d * i]) return false;
        g.B[i * d + j] = P[i + d * j];
      }
    }
  } else {
    if (g.slab_nodes.size() != d) return false;
    for (unsigned long i = 0; i < d; ++i) {
      StochasticNode const *slab = g.slab_nodes[i];
      auto const &par = slab->parents();
      if (slab->distribution()->name() != "dnorm" || slab->length() != 1 ||
          slab->lowerBound() != nullptr || slab->upperBound() != nullptr || par.size() != 2 ||
          !fixed_scalar(par[0]) || !fixed_scalar(par[1])) return false;
      g.m0[i] = par[0]->value(0)[0];
      g.B[i * d + i] = par[1]->value(0)[0];
    }
  }
  for (double v : g.m0) if (!std::isfinite(v)) return false;
  for (double v : g.B) if (!std::isfinite(v)) return false;
  g.Bchol = g.B;
  if (!cholesky(g.Bchol, d)) return false;
  g.log_det_B = log_det_cholesky(g.Bchol, d);
  return std::isfinite(g.log_det_B);
}

bool read_latent(StochasticNode const *z)
{
  auto const &par = z->parents();
  return z->distribution()->name() == "dnorm" && z->length() == 1 && z->lowerBound() == nullptr &&
    z->upperBound() == nullptr && par.size() == 2 && fixed_value(par[0], 0.0) && fixed_value(par[1], 1.0);
}

bool unobserved(StochasticNode const *node)
{
  for (unsigned long i = 0; i < node->length(); ++i) if (node->isObserved(i)) return false;
  return true;
}

// The check of a probed design column (spec: Eligibility, the factory's
// probes): the children's means mu0 at the baseline and mu1, mu2 with one
// coordinate perturbed by h and 2h. Each row is verified against its own
// response per unit (the column's scale is the largest): the perturbation is
// representable when the rounding of the row's baseline and perturbed means
// is below 1e-10 of its response times h, and the response is linear when the
// row's two responses agree to 1e-9 of it. Checked against the column's
// largest response, one row's entry could be lost to the rounding of its
// baseline while another row passed. A row whose means do not change (a zero
// entry) has no response of its own: its rounding is bounded by 1e-10 of the
// scale times h (a column lost entirely to the rounding of a large baseline
// would pass as zero), and an entry below that bound is the run-time
// comparison's (correct(): a mean the sweep did not move must not change, and
// a mean it moved must be the one the design predicts). Returns
// column_verified, column_unrepresentable or column_nonlinear (also for a mean
// that is not finite); *scale receives the scale.
enum ColumnCheck { column_verified, column_unrepresentable, column_nonlinear };

ColumnCheck column_check(Vec const &mu0, Vec const &mu1, Vec const &mu2, double h, double *scale)
{
  const double eps = std::numeric_limits<double>::epsilon();
  double s = 0.0;
  for (std::size_t r = 0; r < mu0.size(); ++r) {
    if (!std::isfinite(mu0[r]) || !std::isfinite(mu1[r]) || !std::isfinite(mu2[r])) return column_nonlinear;
    s = std::max(s, std::max(std::fabs(mu1[r] - mu0[r]) / h, std::fabs(mu2[r] - mu0[r]) / (2.0 * h)));
  }
  *scale = s;
  if (!std::isfinite(s)) return column_nonlinear;
  for (std::size_t r = 0; r < mu0.size(); ++r) {
    const double response = std::max(std::fabs(mu1[r] - mu0[r]) / h, std::fabs(mu2[r] - mu0[r]) / (2.0 * h));
    const double size = std::max(std::fabs(mu0[r]), std::max(std::fabs(mu1[r]), std::fabs(mu2[r])));
    if (!(4.0 * eps * size <= 1e-10 * (response > 0.0 ? response : s) * h)) return column_unrepresentable;
  }
  for (std::size_t r = 0; r < mu0.size(); ++r) {
    const double response = std::max(std::fabs(mu1[r] - mu0[r]) / h, std::fabs(mu2[r] - mu0[r]) / (2.0 * h));
    if (!(std::fabs((mu2[r] - mu0[r]) - 2.0 * (mu1[r] - mu0[r])) <= 2e-9 * response * h)) return column_nonlinear;
  }
  return column_verified;
}

#include "LocationBlockBias.h"

// The declared gate chart: physical = slab * (gate == label) + zero * (gate == other)
// in either order of the sum.
bool match_chart(GroupSpec &g)
{
  Node const *physical = g.physical;
  if (physical->length() != g.d || !binary(physical, "(a+b)")) return fail("chart: physical is not a sum");
  auto test = [&](Node const *node, double label) {
    return binary(node, "(a==b)") && node->parents()[0] == g.gate && fixed_value(node->parents()[1], label);
  };
  for (unsigned order = 0; order < 2; ++order) {
    Node const *on = physical->parents()[order];
    Node const *off = physical->parents()[1 - order];
    if (binary(on, "(a*b)") && binary(off, "(a*b)") && same_slab(on->parents()[0], g.slab) &&
        test(on->parents()[1], g.label) && fixed_value(off->parents()[0], 0.0) &&
        off->parents()[0]->length() == g.d && test(off->parents()[1], 3.0 - g.label)) {
      g.chart = {on->parents()[1], off->parents()[1], on, off};
      // an aggregate of scalar slab coordinates is part of the chart
      if (dynamic_cast<StochasticNode const *>(on->parents()[0]) == nullptr) g.chart.push_back(on->parents()[0]);
      return true;
    }
  }
  return fail("chart: no match");
}

ChildType child_type(StochasticNode const *child, bool *ok)
{
  const std::string name = child->distribution()->name();
  *ok = true;
  if (name == "dnorm") return ChildType::normal;
  if (name == "dselnorm_step_switch") return ChildType::step_switch;
  if (name == "dselnorm_cluster_step") return ChildType::cluster;
  if (name == "dselnorm_factor_step") return ChildType::factor;
  if (name == "dselnorm_mnorm_step") return ChildType::mv;
  *ok = false;
  return ChildType::normal;
}

// Structural admission (spec: Eligibility, native factory). The members enter
// their stochastic children only through the mean parameter; along every path
// that reaches a child, a node depending on the physical coefficients is a
// linear function of them with fixed coefficients (JAGS' closed-class check),
// a node depending on the latents is linear in them, the gates enter only
// through their declared charts and the gated slabs only through the charts.
// child_latent receives the latent of each child's mean (-1 none); latents_ok
// is false when a child's mean depends on two latents.
bool verify_structure(GraphView const &view, Graph const &graph, std::vector<GroupSpec> const &groups,
                      std::vector<StochasticNode *> const &latents, BiasSpec *bias,
                      std::vector<int> &child_latent, bool &latents_ok)
{
  std::set<Node const *> coef, zdep, gated, raw, members;
  std::map<Node const *, int> latent_id;
  std::set<Node const *> chart_nodes, physical_nodes;
  // an owned bias label (LocationBlockBias.h): the indicator and its
  // deterministic descendants outside the slope charts are label nodes; the
  // slope sources enter only through their charts, whose physical slopes are
  // coefficients of the means
  std::set<Node const *> labels, slope_raw;
  if (bias != nullptr) {
    members.insert(bias->indicator);
    labels.insert(bias->indicator);
    bias->label_nodes.clear();
    for (SlopeSpec const &spec : bias->slope) {
      if (spec.source == nullptr) continue;
      members.insert(spec.source);
      slope_raw.insert(spec.source);
      physical_nodes.insert(spec.physical);
    }
  }
  for (GroupSpec const &g : groups) {
    for (StochasticNode const *node : g.slab_nodes) members.insert(node);
    if (g.gate == nullptr) {
      coef.insert(g.slab);
      continue;
    }
    members.insert(g.gate);
    gated.insert(g.gate);
    for (StochasticNode const *node : g.slab_nodes) raw.insert(node);
    for (Node const *node : g.chart) chart_nodes.insert(node);
    physical_nodes.insert(g.physical);
  }
  for (unsigned s = 0; s < latents.size(); ++s) {
    members.insert(latents[s]);
    zdep.insert(latents[s]);
    latent_id[latents[s]] = static_cast<int>(s);
  }

  std::vector<DeterministicNode *> const &dnodes = view.deterministicChildren();
  std::map<Node const *, std::size_t> index;
  for (std::size_t i = 0; i < dnodes.size(); ++i) index[dnodes[i]] = i;
  // A node matters when a stochastic child depends on it.
  std::vector<char> reaches(dnodes.size(), 0);
  for (std::size_t i = dnodes.size(); i-- > 0;) {
    DeterministicNode *node = dnodes[i];
    if (!node->stochasticChildren()->empty()) reaches[i] = 1;
    for (DeterministicNode *child : *node->deterministicChildren()) {
      auto found = index.find(child);
      if (found != index.end() && reaches[found->second]) reaches[i] = 1;
    }
  }

  for (std::size_t i = 0; i < dnodes.size(); ++i) {
    DeterministicNode const *node = dnodes[i];
    if (chart_nodes.count(node)) {
      gated.insert(node);
      continue;
    }
    if (physical_nodes.count(node)) {
      coef.insert(node);
      continue;
    }
    bool from_coef = false, from_z = false, from_gate = false, from_raw = false;
    bool from_label = false, from_slope = false;
    int id = -1;
    for (Node const *parent : node->parents()) {
      if (coef.count(parent)) from_coef = true;
      if (gated.count(parent)) from_gate = true;
      if (raw.count(parent)) from_raw = true;
      if (labels.count(parent)) from_label = true;
      if (slope_raw.count(parent)) from_slope = true;
      if (zdep.count(parent)) {
        from_z = true;
        const int other = latent_id[parent];
        if (other == -2 || (id >= 0 && other != id)) id = -2;
        else if (id == -1) id = other;
      }
    }
    if (from_label) labels.insert(node);
    if (!reaches[i]) continue;
    if (from_gate || from_raw) return fail("structure: gate or raw slab outside its chart");
    if (from_slope) return fail("structure: slope source outside its chart");
    if (from_label && (from_coef || from_z)) return fail("structure: the bias label enters a mean outside the slope charts");
    if (from_label) bias->label_nodes.push_back(const_cast<DeterministicNode *>(node));
    if (from_coef) {
      if (!node->isClosed(coef, DNODE_LINEAR, true)) return fail("structure: not linear in the coefficients with fixed coefficients");
      coef.insert(node);
    }
    if (from_z) {
      if (!node->isClosed(zdep, DNODE_LINEAR, false)) return fail("structure: not linear in the latents");
      zdep.insert(node);
      latent_id[node] = id;
    }
  }

  // The gates have no other dependants, the gated slabs only the charts.
  for (GroupSpec const &g : groups) {
    if (g.gate == nullptr) continue;
    for (DeterministicNode *child : *g.gate->deterministicChildren()) {
      // the declaration tags of this sampler hang off the gates; a finite
      // readout tag or any other dependant declines
      if (chart_nodes.count(child) || terminal_tag(child, group_tag_name, 5, graph) ||
          terminal_tag(child, setup_tag_name, 3, graph)) continue;
      return fail("structure: gate has another dependant");
    }
    for (StochasticNode *node : g.slab_nodes) {
      for (DeterministicNode *child : *node->deterministicChildren()) {
        if (!chart_nodes.count(child) && index.count(child)) {
          if (reaches[index[child]]) return fail("structure: slab has another dependant");
        }
      }
    }
  }

  // The slope sources: their chart and the bias tag are their only
  // deterministic dependants, gate markers their only stochastic ones.
  if (bias != nullptr) {
    for (SlopeSpec const &spec : bias->slope) {
      if (spec.source == nullptr) continue;
      for (DeterministicNode *child : *spec.source->deterministicChildren()) {
        if (child == spec.physical || terminal_tag(child, bias_tag_name, 5, graph)) continue;
        return fail("structure: slope source has another dependant");
      }
      for (StochasticNode *child : *spec.source->stochasticChildren()) {
        if (!bias->skip.count(child)) return fail("structure: slope source has a stochastic dependant");
      }
    }
  }

  latents_ok = true;
  child_latent.clear();
  for (StochasticNode const *child : view.stochasticChildren()) {
    // gate markers and the nested weight node of an owned bias label are not
    // likelihood children (LocationBlockBias.h)
    if (bias != nullptr && bias->skip.count(child)) continue;
    if (members.count(child)) return fail("structure: a member is a child");
    bool ok;
    const ChildType type = child_type(child, &ok);
    if (!ok || !isObserved(child) || child->lowerBound() != nullptr || child->upperBound() != nullptr) {
      return fail("structure: child distribution, observation or bounds");
    }
    auto const &par = child->parents();
    for (std::size_t i = 1; i < par.size(); ++i) {
      Node const *parent = par[i];
      if (members.count(parent) || coef.count(parent) || zdep.count(parent) || gated.count(parent) ||
          raw.count(parent)) return fail("structure: member enters a non-mean parameter");
      if (labels.count(parent) && !label_slot_allowed(type, i)) {
        return fail("structure: the bias label enters a likelihood parameter the block reads");
      }
    }
    Node const *mean = par[0];
    if (gated.count(mean) || raw.count(mean) || labels.count(mean)) return fail("structure: mean depends on a gate or raw slab");
    if (!coef.count(mean) && !zdep.count(mean)) return fail("structure: mean does not depend on the members");
    int id = -1;
    if (zdep.count(mean)) id = latent_id[mean];
    if (id == -2) latents_ok = false;
    child_latent.push_back(id);
  }
  return true;
}

// ---- the sampler

struct Line {
  std::vector<unsigned> moving;    // coordinates of the moving group(s)
  std::vector<unsigned> partners;  // coordinates of the partners
  Vec m, Sigma, P, Pchol;          // marginal law of the moving coordinates given non-partners
  double log_det_P = 0;
  Vec KQ;                          // partners x moving
  std::vector<unsigned> lats;      // latents shifted by the line
  Vec KZ;                          // lats x moving
};

struct Shift {
  std::vector<unsigned> coords;
  Vec dcoef;
  std::vector<unsigned> lats;
  Vec dlat;
  std::vector<int> slopes;         // physical slope changes (owned bias label)
  Vec dslope;
  std::vector<unsigned> children;
  std::vector<unsigned> offset;    // into delta / wdelta per listed child
  Vec delta, wdelta, q1, q2, value;
  // the prior of shifted active slabs and latents: p1 t - p2 t^2 / 2
  double p1 = 0, p2 = 0;

  void clear()
  {
    coords.clear();
    dcoef.clear();
    lats.clear();
    dlat.clear();
    slopes.clear();
    dslope.clear();
  }
};

class LocationBlock : public MutableSampleMethod, public CorrectedSamplerWindow {
  GraphView const *_view;
  std::shared_ptr<Shared const> _s;
  const unsigned int _chain;
  std::vector<std::unique_ptr<SelectedChild>> _surrogates;  // per child, null without one
  bool _adapting = true;
  bool _exact = false;
  bool _registered = false;
  bool _injected = false;          // mutant 14 has thrown
  CorrectedTally _tally{corrected_kind_location};
  // the refresh kernel of each group in selection labels: the tilted
  // independence proposal, or (acceptance below 0.4 while adapting) the slice
  // update along the line, decided at the first adaptOff()
  std::vector<std::uint64_t> _refresh_proposed, _refresh_accepted;
  std::vector<char> _slice_refresh;
  // The Newton fits evaluate each cluster or factor child at the coarsest rung
  // of its schedule; while adapting, the stage-2 acceptance of the tilted flips
  // and swaps that passed stage 1 is tallied, and below 0.9 (50 proposals at
  // least) the fits use the sampler's rule from the first adaptOff() on.
  bool _coarse_fit = false;
  std::uint64_t _fit_proposed = 0, _fit_accepted = 0;
  // the evidence record of the gates (GateEvidence), flushed after each update
  std::vector<GateEvidence> _gate_evidence;
  bool _gate_dirty = false;

  // A component skipped because a factorization of the working model failed
  // (a matrix not numerically positive definite): counted in the report.
  void skipped_factorization() { count(counters().factorization_skips); }

  void evidence_proposed(unsigned j, bool on)
  {
    if (_adapting) return;
    GateEvidence &e = _gate_evidence[j];
    (on ? e.off_proposed : e.on_proposed) += 1;
    _gate_dirty = true;
  }

  void evidence_accepted(unsigned j, bool on)
  {
    if (_adapting) return;
    GateEvidence &e = _gate_evidence[j];
    (on ? e.off_accepted : e.on_accepted) += 1;
  }

  void evidence_screened(unsigned j)
  {
    if (!_adapting) _gate_evidence[j].screened_out += 1;
  }

  // log odds of the other state against the current one of a proposal that
  // reached stage 1
  void evidence_odds(unsigned j, double other)
  {
    if (_adapting) return;
    GateEvidence &e = _gate_evidence[j];
    e.odds_count += 1;
    e.odds_sum += other;
    e.odds_max = std::max(e.odds_max, other);
  }

  // After the exact correction: each gate's transition between the update's
  // entry and the state it retains (kept: the sweep's final state).
  void evidence_endpoint(bool kept)
  {
    if (_adapting) return;
    Shared const &S = s();
    for (unsigned j = 0; j < S.groups.size(); ++j) {
      GroupSpec const &g = S.groups[j];
      if (g.gate == nullptr) continue;
      GateEvidence &e = _gate_evidence[j];
      e.updates += 1;
      const bool entry = _entry[g.gate_pos] == g.label;
      const bool retained = kept ? _g[j] != 0 : entry;
      if (retained != entry) (retained ? e.switched_on : e.switched_off) += 1;
    }
    _gate_dirty = true;
  }

  void flush_evidence()
  {
    if (!_gate_dirty) return;
    GateEvidenceStore &store = gate_evidence_store();
    const std::lock_guard<std::mutex> lock(store.mutex);
    for (unsigned j = 0; j < _gate_evidence.size(); ++j) {
      if (s().groups[j].gate == nullptr) continue;
      store.records[std::make_pair(j + 1, _chain)].add(_gate_evidence[j]);
      _gate_evidence[j] = GateEvidence();
    }
    _gate_dirty = false;
  }
  // The part of the update that runs (its child evaluations are counted, and
  // with the timing switch its exclusive wall time).
  Category _category = cat_move;
  std::chrono::steady_clock::time_point _part_start;

  // Enters a part of the update and returns the previous one.
  Category enter(Category part)
  {
    // the bias label's swap L2 counts its Laplace fits as its own work
    if (_in_swap) part = cat_swap;
    if (s().timing) {
      const auto now = std::chrono::steady_clock::now();
      counters().nanoseconds[_category].fetch_add(static_cast<std::uint64_t>(
        std::chrono::duration_cast<std::chrono::nanoseconds>(now - _part_start).count()),
        std::memory_order_relaxed);
      _part_start = now;
    }
    const Category previous = _category;
    _category = part;
    return previous;
  }

  // ---- per-update state
  Vec _entry, _final;
  std::vector<char> _g;            // gate on
  Vec _beta, _b, _z;               // slabs, physical coefficients, latents
  Vec _y, _mean_entry, _mean, _r, _wr, _d, _W;
  // the update's entry: the physical coefficients and, per row, the size of
  // the terms of the means at the latent probe (prepare())
  Vec _b_entry, _d_size;
  std::vector<char> _selected, _value_ok, _dirty;
  Vec _value, _value_entry;        // inner-target terms of the selection children
  Vec _fit_value;                  // the coarse rung's terms at the current state
  std::vector<char> _fit_ok;
  std::vector<unsigned> _stamp;
  unsigned _stamp_id = 0;
  bool _selection = false;         // some child is a selection child in step mode
  bool _backward = false;          // this sweep runs the component list reversed
  // The surrogate replaces the exact selection densities: a property of the fit
  // (corrected and not fallen back to exact mode), so that every component of
  // a sweep targets one inner target whatever the label (sweep condition 2);
  // in a Gaussian label no surrogate is evaluated.
  bool _surrogate_mode = false;
  bool _not_exact = false;         // a surrogate evaluated this update was not exact
  long double _repr_delta = 0;     // accumulated log density change of the inner target
  Vec _correction_terms;           // scratch of correct(): the entry's exact terms
  std::vector<unsigned> _unmoved_changed;  // scratch of correct(): unmoved children whose graph mean changed
  Vec _unmoved_terms;              // scratch of correct(): their entry terms
  Vec _a, _lambda, _u, _P, _h;     // Gaussian working model
  std::vector<double const *> _params;
  Vec _scratch;

  // The parameters of each step-switch child, read once per update (the
  // analytic tilt and the direct row densities): valid when they pass the
  // distribution's parameter check at the entry mean.
  struct StepChild {
    bool valid = false;
    bool free_ok = false;          // the checks that do not depend on the bias label
    bool step_ok = false;          // the step-mode checks that do not (bin, bounds)
    double sd = 1, se = 1, power = 1;
    int bin = 0, sign = 1, mode = 0, telescope = 0, bins = 0;
    double const *omega = nullptr;
    std::size_t bounds = 0;        // offset of the decoded lower and upper bounds in _bounds
  };
  std::vector<StepChild> _step;
  Vec _bounds;
  std::vector<int> _tilt_position; // scratch of make_tilt, -1 outside it
  Vec _dc_scratch, _dz_scratch;    // scratch of prepare_shift

  Shared const &s() const { return *_s; }

  // ---- child densities

  void read_params(ChildSpec const &c)
  {
    auto const &par = c.node->parents();
    _params.resize(c.lengths.size());
    for (std::size_t i = 0; i < c.lengths.size(); ++i) _params[i] = par[i]->value(_chain);
  }

  // The child's own exact density at the given mean (PDF_LIKELIHOOD).
  double exact_with_mean(unsigned child, double const *mean)
  {
    ChildSpec const &c = s().children[child];
    count(counters().exact_evaluations);
    count(counters().evaluations[_category]);
    if (c.type == ChildType::step_switch) {
      // the row's density with its parameters read once per update (B1); of
      // the distribution's parameter check only the mean's finiteness depends
      // on the mean
      StepChild const &st = _step[child];
      if (!st.valid || !std::isfinite(mean[0])) return negative_infinity;
      return step_row_log_density(step_law(child), st.power, st.bin, st.mode, st.telescope != 0,
                                  c.node->value(_chain)[0], mean[0]);
    }
    read_params(c);
    _params[0] = mean;
    auto const *dist = static_cast<VectorDist const *>(c.node->distribution());
    if (c.type == ChildType::normal) {
      const double tau = *_params[1];
      const double x = c.node->value(_chain)[0];
      return 0.5 * (std::log(tau) - log_two_pi) - 0.5 * tau * (x - mean[0]) * (x - mean[0]);
    }
    if (!dist->checkParameterValue(_params, c.lengths)) return negative_infinity;
    return dist->logDensity(c.node->value(_chain), PDF_LIKELIHOOD, _params, c.lengths);
  }

  // The child's surrogate at a mean with the given settings (the sampler's,
  // or the coarsest rung for the Newton fits); not_exact flags an inexact
  // normalizer.
  double surrogate_with_mean(unsigned child, double const *mean, SelNormCoarseSettings const &settings,
                             bool *not_exact)
  {
    SelectedChild const &sc = *_surrogates[child];
    for (std::size_t i = 0; i < sc.parents.size(); ++i) sc.parameters[i] = sc.parents[i]->value(_chain);
    sc.parameters[0] = mean;
    bool valid;
    if (sc.controls_fixed) {
      if (!sc.controls_checked) {
        sc.controls_valid = sc.checkControl();
        sc.controls_checked = true;
      }
      valid = sc.controls_valid && sc.checkDynamic();
    } else {
      valid = static_cast<VectorDist const *>(sc.distribution)->checkParameterValue(sc.parameters, sc.lengths);
    }
    count(counters().surrogate_evaluations);
    count(counters().evaluations[_category]);
    if (!valid) return negative_infinity;
    return sc.surrogate(_chain, settings, not_exact);
  }

  double surrogate_with_mean(unsigned child, double const *mean)
  {
    return surrogate_with_mean(child, mean, s().settings, &_not_exact);
  }

  // Whether the Newton fits evaluate this child at its coarsest rung (A1): a
  // cluster or factor child with a surrogate, while the coarse fits are on.
  bool coarse_child(unsigned child) const
  {
    return _coarse_fit && _surrogates[child] &&
      (_surrogates[child]->kind == ChildKind::cluster || _surrogates[child]->kind == ChildKind::factor);
  }

  // The coarse rung's term of a child at a mean; it is never the inner target,
  // so it is no evidence about the surrogate.
  double coarse_value(unsigned child, double const *mean)
  {
    bool ignored = false;
    const double v = surrogate_with_mean(child, mean, s().fit_settings, &ignored);
    if (std::isnan(v) || v == std::numeric_limits<double>::infinity()) {
      throwNodeError(s().children[child].node,
        "Location block sampler: the log density of this likelihood node is not defined at a proposal");
    }
    return v;
  }

  void ensure_fit_value(unsigned child)
  {
    if (_fit_ok[child]) return;
    ChildSpec const &c = s().children[child];
    _fit_value[child] = coarse_value(child, &_mean[c.row0]);
    _fit_ok[child] = 1;
  }

  // The inner target's term of a selection child at a mean.
  double selected_value(unsigned child, double const *mean)
  {
    const double v = (_surrogate_mode && _surrogates[child]) ? surrogate_with_mean(child, mean) :
      exact_with_mean(child, mean);
    if (std::isnan(v) || v == std::numeric_limits<double>::infinity()) {
      throwNodeError(s().children[child].node,
        "Location block sampler: the log density of this likelihood node is not defined at a proposal");
    }
    return v;
  }

  bool gaussian_child(unsigned child) const { return !_selection || !_selected[child]; }

  void ensure_value(unsigned child)
  {
    if (_value_ok[child]) return;
    ChildSpec const &c = s().children[child];
    // A child first evaluated before any move changes it (its mean and the
    // label) gives its entry value.
    _value[child] = selected_value(child, &_mean[c.row0]);
    _value_ok[child] = 1;
    if (!_dirty[child] && _label == _label_entry) {
      _value_entry[child] = _value[child];
      _entry_valued[child] = 1;
    }
  }

  // ---- per-update preparation

  bool precision(ChildSpec const &c, double *W)
  {
    auto const &par = c.node->parents();
    const unsigned n = c.n;
    switch (c.type) {
    case ChildType::normal: {
      const double tau = par[1]->value(_chain)[0];
      if (!(tau > 0) || !std::isfinite(tau)) return false;
      W[0] = tau;
      return true;
    }
    case ChildType::step_switch: {
      const double sd = par[1]->value(_chain)[0];
      const double power = par[3]->value(_chain)[0];
      if (!(sd > 0) || !(power > 0) || !std::isfinite(power / (sd * sd))) return false;
      W[0] = power / (sd * sd);
      return true;
    }
    case ChildType::cluster: {
      double const *sd = par[1]->value(_chain);
      double const *loading = par[2]->value(_chain);
      Vec dinv(n);
      long double capacity = 1;
      for (unsigned i = 0; i < n; ++i) {
        if (!(sd[i] > 0) || !std::isfinite(sd[i]) || !std::isfinite(loading[i])) return false;
        dinv[i] = 1.0 / (sd[i] * sd[i]);
        capacity += static_cast<long double>(loading[i]) * loading[i] * dinv[i];
      }
      for (unsigned i = 0; i < n; ++i) {
        for (unsigned j = 0; j < n; ++j) {
          const double low_rank = static_cast<double>(dinv[i] * loading[i] * loading[j] * dinv[j] / capacity);
          W[i * n + j] = (i == j ? dinv[i] : 0.0) - low_rank;
        }
      }
      return true;
    }
    case ChildType::factor: {
      double const *sd = par[1]->value(_chain);
      double const *loading = par[2]->value(_chain);
      const unsigned r = c.n == 0 ? 0 : static_cast<unsigned>(c.lengths[2] / c.n);
      Vec dinv(n), C(r * r, 0.0), DL(n * r);
      for (unsigned i = 0; i < n; ++i) {
        if (!(sd[i] > 0) || !std::isfinite(sd[i])) return false;
        dinv[i] = 1.0 / (sd[i] * sd[i]);
        for (unsigned k = 0; k < r; ++k) DL[i * r + k] = dinv[i] * loading[i + n * k];
      }
      for (unsigned k = 0; k < r; ++k) {
        for (unsigned l = 0; l < r; ++l) {
          long double v = k == l ? 1.0 : 0.0;
          for (unsigned i = 0; i < n; ++i) v += static_cast<long double>(loading[i + n * k]) * DL[i * r + l];
          C[k * r + l] = static_cast<double>(v);
        }
      }
      if (!cholesky(C, r)) return false;
      const Vec Cinv = inverse_from_cholesky(C, r);
      for (unsigned i = 0; i < n; ++i) {
        for (unsigned j = 0; j < n; ++j) {
          long double v = 0;
          for (unsigned k = 0; k < r; ++k) {
            for (unsigned l = 0; l < r; ++l) v += static_cast<long double>(DL[i * r + k]) * Cinv[k * r + l] * DL[j * r + l];
          }
          W[i * n + j] = (i == j ? dinv[i] : 0.0) - static_cast<double>(v);
        }
      }
      return true;
    }
    case ChildType::mv: {
      double const *packed = par[1]->value(_chain);
      Vec S(n * n);
      std::size_t position = 0;
      for (unsigned col = 0; col < n; ++col) {
        for (unsigned row = col; row < n; ++row) {
          S[row * n + col] = packed[position];
          S[col * n + row] = packed[position];
          ++position;
        }
      }
      if (!cholesky(S, n)) return false;
      const Vec inv = inverse_from_cholesky(S, n);
      std::copy(inv.begin(), inv.end(), W);
      return true;
    }
    }
    return false;
  }

  bool is_selected(ChildSpec const &c) const { return kernel_selected(c, _chain); }

  // The parameters of a step-switch child for this update (StepChild).
  void read_step(unsigned c)
  {
    ChildSpec const &child = s().children[c];
    auto const &par = child.node->parents();
    StepChild &st = _step[c];
    st.sd = par[1]->value(_chain)[0];
    st.se = par[2]->value(_chain)[0];
    st.power = par[3]->value(_chain)[0];
    st.omega = par[4]->value(_chain);
    st.bins = static_cast<int>(child.lengths[4]);
    st.bin = static_cast<int>(par[7]->value(_chain)[0]);
    st.sign = static_cast<int>(par[8]->value(_chain)[0]);
    st.mode = static_cast<int>(par[9]->value(_chain)[0]);
    st.telescope = static_cast<int>(par[10]->value(_chain)[0]);
    st.bounds = _bounds.size();
    for (unsigned side = 5; side <= 6; ++side) {
      double const *values = par[side]->value(_chain);
      for (int b = 0; b < st.bins; ++b) _bounds.push_back(selnorm_jags_bound_to_ieee(values[b]));
    }
    read_params(child);
    st.valid = static_cast<VectorDist const *>(child.node->distribution())->checkParameterValue(_params, child.lengths);
    // the label-free parts of the distribution's parameter check, for the
    // validity at another bias label (LocationBlockBiasMoves.h, refresh_step)
    st.free_ok = st.sd > 0 && st.se > 0 && st.power > 0 && (st.sign == 1 || st.sign == -1) &&
      (st.telescope == 0 || st.telescope == 1) && par[8]->value(_chain)[0] == st.sign &&
      par[10]->value(_chain)[0] == st.telescope;
    st.step_ok = selnorm_jags_integer_in_range(par[7]->value(_chain)[0], 0, st.bins) &&
      selnorm_is_descending_step_partition(&_bounds[st.bounds], &_bounds[st.bounds + static_cast<std::size_t>(st.bins)], st.bins);
    for (int b = 0; st.step_ok && b < st.bins; ++b) {
      st.step_ok = _bounds[st.bounds + static_cast<std::size_t>(b)] < _bounds[st.bounds + static_cast<std::size_t>(st.bins + b)];
    }
  }

  StepRowLaw step_law(unsigned c) const
  {
    StepChild const &st = _step[c];
    StepRowLaw law;
    law.sd = st.sd;
    law.se = st.se;
    law.effect_sign = st.sign;
    law.bins = st.bins;
    law.omega = st.omega;
    law.lower = &_bounds[st.bounds];
    law.upper = &_bounds[st.bounds + static_cast<std::size_t>(st.bins)];
    return law;
  }

  // The working residuals r = y - mean and weighted residuals W r from the
  // working means (prepare(); the swap's probe, fresh_working_means()).
  void working_residuals()
  {
    Shared const &sh = s();
    for (unsigned i = 0; i < sh.rows; ++i) _r[i] = _y[i] - _mean[i];
    for (ChildSpec const &child : sh.children) {
      double const *W = &_W[child.woff];
      for (unsigned i = 0; i < child.n; ++i) {
        long double v = 0;
        for (unsigned j = 0; j < child.n; ++j) v += static_cast<long double>(W[i * child.n + j]) * _r[child.row0 + j];
        _wr[child.row0 + i] = static_cast<double>(v);
      }
    }
  }

  void prepare()
  {
    Shared const &sh = s();
    const unsigned G = static_cast<unsigned>(sh.groups.size());
    const unsigned S = static_cast<unsigned>(sh.latents.size());
    const unsigned p = sh.p;
    _view->getValue(_entry, _chain);
    for (unsigned j = 0; j < G; ++j) {
      GroupSpec const &g = sh.groups[j];
      _g[j] = g.gate == nullptr || _entry[g.gate_pos] == g.label;
      for (unsigned i = 0; i < g.d; ++i) {
        _beta[g.coord0 + i] = _entry[g.slab_pos + i];
        _b[g.coord0 + i] = _g[j] ? _beta[g.coord0 + i] : 0.0;
      }
    }
    for (unsigned k = 0; k < S; ++k) _z[k] = _entry[sh.latent_pos[k]];

    _selection = false;
    _bounds.clear();
    for (unsigned c = 0; c < sh.children.size(); ++c) {
      ChildSpec const &child = sh.children[c];
      double const *mean = child.node->parents()[0]->value(_chain);
      double const *y = child.node->value(_chain);
      for (unsigned i = 0; i < child.n; ++i) {
        _mean_entry[child.row0 + i] = mean[i];
        _y[child.row0 + i] = y[i];
      }
      if (!precision(child, &_W[child.woff])) {
        throwNodeError(child.node, "Location block sampler: the covariance of this likelihood node is not positive definite");
      }
      if (child.type == ChildType::step_switch) read_step(c);
      _selected[c] = is_selected(child);
      if (_selected[c]) _selection = true;
      _value_ok[c] = 0;
      _fit_ok[c] = 0;
      _dirty[c] = 0;
      _entry_valued[c] = 0;
    }
    _surrogate_mode = sh.corrected && !_exact;
    _selection_entry = _selection;
    _line.valid = false;
    _not_exact = false;
    _repr_delta = 0;

    // d: the change of each child's mean per unit of its latent (linear,
    // verified at construction), probed at the current non-members. The
    // probe's baseline is the sweep's own state, so its rounding is that of
    // the means themselves (unlike the factory's probe, whose baseline holds
    // the initial values); a mean change the block does not track is an
    // error in correct().
    std::fill(_d.begin(), _d.end(), 0.0);
    if (S > 0) {
      Vec shifted(_z);
      for (double &v : shifted) v += 1.0;
      sh.latent_view->setValue(shifted, _chain);
      for (ChildSpec const &child : sh.children) {
        if (child.latent < 0) continue;
        double const *mean = child.node->parents()[0]->value(_chain);
        for (unsigned i = 0; i < child.n; ++i) _d[child.row0 + i] = mean[i] - _mean_entry[child.row0 + i];
      }
      sh.latent_view->setValue(_z, _chain);
    }

    _mean = _mean_entry;
    working_residuals();

    // Gaussian working model with the latents integrated: P = B + X' S^-1 X,
    // h = B m0 + X' S^-1 r0, r0 = y - offset; per latent the within/between
    // form (spec: Numerics), never a difference of two large terms.
    std::fill(_P.begin(), _P.end(), 0.0);
    std::fill(_h.begin(), _h.end(), 0.0);
    std::fill(_u.begin(), _u.end(), 0.0);
    std::fill(_a.begin(), _a.end(), 0.0);
    // with the size of the means' terms at the latent probe (the entry and the
    // latents one above it; the bias slopes at the entry label): the rounding
    // of _d is at most 4 eps times it (correct(), the tracked comparison)
    _b_entry = _b;
    double slope_entry[2] = {0.0, 0.0};
    if (sh.bias != nullptr) {
      for (int b = 0; b < 2; ++b) {
        SlopeSpec const &spec = sh.bias->slope[b];
        if (spec.source != nullptr && spec.label == static_cast<int>(_entry[sh.bias->pos])) slope_entry[b] = _entry[spec.pos];
      }
    }
    Vec r0(sh.rows), WX, Wd, Wr0;
    for (ChildSpec const &child : sh.children) {
      for (unsigned i = 0; i < child.n; ++i) {
        const unsigned row = child.row0 + i;
        long double v = _r[row];
        double size = std::fabs(_mean_entry[row]) + std::fabs(_d[row]);
        for (unsigned k = 0; k < p; ++k) {
          v += static_cast<long double>(sh.X[row * p + k]) * _b[k];
          size += std::fabs(sh.X[row * p + k] * _b[k]);
        }
        if (child.latent >= 0) {
          v += static_cast<long double>(_d[row]) * _z[child.latent];
          size += std::fabs(_d[row] * _z[child.latent]);
        }
        for (int b = 0; b < 2 && sh.bias != nullptr; ++b) {
          if (slope_entry[b] != 0.0) size += std::fabs(sh.bias->slope[b].column[row] * slope_entry[b]);
        }
        r0[row] = static_cast<double>(v);
        _d_size[row] = size;
      }
    }
    std::vector<long double> M(p * p, 0.0L), v(p, 0.0L);
    auto child_products = [&](ChildSpec const &child, Vec &wx, Vec &wd, Vec &wr0) {
      const unsigned n = child.n;
      double const *W = &_W[child.woff];
      wx.assign(n * p, 0.0);
      wd.assign(n, 0.0);
      wr0.assign(n, 0.0);
      for (unsigned i = 0; i < n; ++i) {
        for (unsigned j = 0; j < n; ++j) {
          const double w = W[i * n + j];
          if (w == 0.0) continue;
          const unsigned row = child.row0 + j;
          for (unsigned k = 0; k < p; ++k) wx[i * p + k] += w * sh.X[row * p + k];
          wd[i] += w * _d[row];
          wr0[i] += w * r0[row];
        }
      }
    };
    // children without a latent
    for (ChildSpec const &child : sh.children) {
      if (child.latent >= 0) continue;
      child_products(child, WX, Wd, Wr0);
      for (unsigned i = 0; i < child.n; ++i) {
        const unsigned row = child.row0 + i;
        for (unsigned k = 0; k < p; ++k) {
          const double x = sh.X[row * p + k];
          if (x == 0.0) continue;
          v[k] += static_cast<long double>(x) * Wr0[i];
          for (unsigned l = 0; l < p; ++l) M[k * p + l] += static_cast<long double>(x) * WX[i * p + l];
        }
      }
    }
    for (unsigned lat = 0; lat < S; ++lat) {
      std::vector<unsigned> const &kids = sh.latent_children[lat];
      // a = d' W d, u = d' W X, rho = d' W r0 (the children's products are
      // recomputed in the second pass rather than stored per child)
      long double a = 0, rho = 0;
      std::vector<long double> u(p, 0.0L);
      for (unsigned c : kids) {
        ChildSpec const &child = sh.children[c];
        child_products(child, WX, Wd, Wr0);
        for (unsigned i = 0; i < child.n; ++i) {
          const unsigned row = child.row0 + i;
          a += static_cast<long double>(_d[row]) * Wd[i];
          rho += static_cast<long double>(Wd[i]) * r0[row];
          for (unsigned k = 0; k < p; ++k) u[k] += static_cast<long double>(Wd[i]) * sh.X[row * p + k];
        }
      }
      _a[lat] = static_cast<double>(a);
      _lambda[lat] = 1.0 + _a[lat];
      for (unsigned k = 0; k < p; ++k) _u[lat * p + k] = static_cast<double>(u[k]);
      const bool between = _a[lat] > 0;
      Vec xbar(p, 0.0);
      double rbar = 0;
      if (between) {
        for (unsigned k = 0; k < p; ++k) xbar[k] = static_cast<double>(u[k] / a);
        rbar = static_cast<double>(rho / a);
      }
      for (std::size_t q = 0; q < kids.size(); ++q) {
        ChildSpec const &child = sh.children[kids[q]];
        child_products(child, WX, Wd, Wr0);
        Vec const &wx = WX;
        Vec const &wd = Wd;
        Vec const &wr0 = Wr0;
        for (unsigned i = 0; i < child.n; ++i) {
          const unsigned row = child.row0 + i;
          const double di = _d[row];
          // centered rows and their products with W
          for (unsigned k = 0; k < p; ++k) {
            const double xc = sh.X[row * p + k] - (between ? di * xbar[k] : 0.0);
            if (xc == 0.0) continue;
            const double wrc = wr0[i] - (between ? wd[i] * rbar : 0.0);
            v[k] += static_cast<long double>(xc) * wrc;
            for (unsigned l = 0; l < p; ++l) {
              const double wxc = wx[i * p + l] - (between ? wd[i] * xbar[l] : 0.0);
              M[k * p + l] += static_cast<long double>(xc) * wxc;
            }
          }
        }
      }
      if (between) {
        const double scale = _a[lat] / _lambda[lat];
        for (unsigned k = 0; k < p; ++k) {
          v[k] += static_cast<long double>(scale) * xbar[k] * rbar;
          for (unsigned l = 0; l < p; ++l) M[k * p + l] += static_cast<long double>(scale) * xbar[k] * xbar[l];
        }
      }
    }
    for (unsigned k = 0; k < p; ++k) {
      for (unsigned l = 0; l < p; ++l) _P[k * p + l] = static_cast<double>(0.5L * (M[k * p + l] + M[l * p + k]));
      _h[k] = static_cast<double>(v[k]);
    }
    for (GroupSpec const &g : sh.groups) {
      for (unsigned i = 0; i < g.d; ++i) {
        long double bm = 0;
        for (unsigned j = 0; j < g.d; ++j) {
          _P[(g.coord0 + i) * p + g.coord0 + j] += g.B[i * g.d + j];
          bm += static_cast<long double>(g.B[i * g.d + j]) * g.m0[j];
        }
        _h[g.coord0 + i] += static_cast<double>(bm);
      }
    }
    bias_prepare();
  }

  // ---- lines

  bool active(unsigned group) const { return _g[group] != 0; }

  // The sheared line of the moving groups: their marginal Gaussian law given
  // the non-partner active coefficients, the regression of the partners on
  // them and the latent shear (spec: The sheared line of a group). Functions
  // of non-members and the other gates only.
  bool build_line(std::vector<unsigned> const &moving_groups, Line &L)
  {
    Shared const &sh = s();
    const unsigned p = sh.p;
    std::vector<char> moving(sh.groups.size(), 0);
    L.moving.clear();
    for (unsigned j : moving_groups) {
      moving[j] = 1;
      for (unsigned i = 0; i < sh.groups[j].d; ++i) L.moving.push_back(sh.groups[j].coord0 + i);
    }
    const std::size_t nm = L.moving.size();
    std::vector<unsigned> others;
    std::vector<unsigned> C(L.moving);
    for (unsigned j = 0; j < sh.groups.size(); ++j) {
      if (moving[j] || !active(j)) continue;
      others.push_back(j);
      for (unsigned i = 0; i < sh.groups[j].d; ++i) C.push_back(sh.groups[j].coord0 + i);
    }
    const std::size_t nc = C.size();
    Vec PC(nc * nc);
    for (std::size_t i = 0; i < nc; ++i) for (std::size_t j = 0; j < nc; ++j) PC[i * nc + j] = _P[C[i] * p + C[j]];
    if (!cholesky(PC, nc)) return false;
    const Vec SC = inverse_from_cholesky(PC, nc);
    // partners: an active group correlated with a moving coordinate
    L.partners.clear();
    std::vector<unsigned> nonpartners;
    std::size_t position = nm;
    for (unsigned j : others) {
      const unsigned d = sh.groups[j].d;
      double best = 0;
      for (unsigned i = 0; i < d; ++i) {
        for (std::size_t m = 0; m < nm; ++m) {
          const std::size_t a = position + i;
          const double corr = SC[a * nc + m] / std::sqrt(SC[a * nc + a] * SC[m * nc + m]);
          best = std::max(best, std::fabs(corr));
        }
      }
      for (unsigned i = 0; i < d; ++i) {
        (best > partner_correlation ? L.partners : nonpartners).push_back(sh.groups[j].coord0 + i);
      }
      position += d;
    }
    std::vector<unsigned> A(L.moving);
    A.insert(A.end(), L.partners.begin(), L.partners.end());
    const std::size_t na = A.size();
    Vec PA(na * na), hA(na);
    for (std::size_t i = 0; i < na; ++i) {
      long double h = _h[A[i]];
      for (unsigned k : nonpartners) h -= static_cast<long double>(_P[A[i] * p + k]) * _b[k];
      hA[i] = static_cast<double>(h);
      for (std::size_t j = 0; j < na; ++j) PA[i * na + j] = _P[A[i] * p + A[j]];
    }
    if (!cholesky(PA, na)) return false;
    const Vec SA = inverse_from_cholesky(PA, na);
    Vec muA(na, 0.0);
    for (std::size_t i = 0; i < na; ++i) {
      long double v = 0;
      for (std::size_t j = 0; j < na; ++j) v += static_cast<long double>(SA[i * na + j]) * hA[j];
      muA[i] = static_cast<double>(v);
    }
    L.m.assign(muA.begin(), muA.begin() + nm);
    L.Sigma.assign(nm * nm, 0.0);
    for (std::size_t i = 0; i < nm; ++i) for (std::size_t j = 0; j < nm; ++j) L.Sigma[i * nm + j] = SA[i * na + j];
    Vec Schol(L.Sigma);
    if (!cholesky(Schol, nm)) return false;
    L.P = inverse_from_cholesky(Schol, nm);
    L.Pchol = L.P;
    if (!cholesky(L.Pchol, nm)) return false;
    L.log_det_P = log_det_cholesky(L.Pchol, nm);
    const std::size_t nq = L.partners.size();
    L.KQ.assign(nq * nm, 0.0);
    for (std::size_t q = 0; q < nq; ++q) {
      for (std::size_t m = 0; m < nm; ++m) {
        long double v = 0;
        for (std::size_t l = 0; l < nm; ++l) v += static_cast<long double>(SA[(nm + q) * na + l]) * L.P[l * nm + m];
        L.KQ[q * nm + m] = static_cast<double>(v);
      }
    }
    // the latent shear: dZ_s = -(u_s[moving] + u_s[partners] KQ) db / lambda_s
    L.lats.clear();
    L.KZ.clear();
    for (unsigned lat = 0; lat < sh.latents.size(); ++lat) {
      bool any = false;
      Vec row(nm, 0.0);
      for (std::size_t m = 0; m < nm; ++m) {
        long double v = _u[lat * p + L.moving[m]];
        for (std::size_t q = 0; q < nq; ++q) v += static_cast<long double>(_u[lat * p + L.partners[q]]) * L.KQ[q * nm + m];
        row[m] = static_cast<double>(-v / _lambda[lat]);
        if (row[m] != 0.0) any = true;
      }
      if (any) {
        L.lats.push_back(lat);
        L.KZ.insert(L.KZ.end(), row.begin(), row.end());
      }
    }
    return true;
  }

  // The shift of a line for a change db of its moving coordinates; scale
  // multiplies the shear (the asymmetric-shear mutant).
  void line_shift(Line const &L, double const *db, Shift &sh, double scale = 1.0)
  {
    const std::size_t nm = L.moving.size();
    sh.clear();
    for (std::size_t m = 0; m < nm; ++m) {
      sh.coords.push_back(L.moving[m]);
      sh.dcoef.push_back(db[m]);
    }
    for (std::size_t q = 0; q < L.partners.size(); ++q) {
      long double v = 0;
      for (std::size_t m = 0; m < nm; ++m) v += static_cast<long double>(L.KQ[q * nm + m]) * db[m];
      sh.coords.push_back(L.partners[q]);
      sh.dcoef.push_back(scale * static_cast<double>(v));
    }
    for (std::size_t z = 0; z < L.lats.size(); ++z) {
      long double v = 0;
      for (std::size_t m = 0; m < nm; ++m) v += static_cast<long double>(L.KZ[z * nm + m]) * db[m];
      sh.lats.push_back(L.lats[z]);
      sh.dlat.push_back(scale * static_cast<double>(v));
    }
  }

  // ---- shifts: the children a change of the physical coefficients and the
  // latents moves, their mean changes, and the Gaussian quantities.

  void prepare_shift(Shift &sh, std::vector<unsigned> const &prior_coords)
  {
    Shared const &S = s();
    const unsigned p = S.p;
    ++_stamp_id;
    sh.children.clear();
    sh.offset.clear();
    Vec &dc = _dc_scratch;
    dc.assign(p, 0.0);
    for (std::size_t i = 0; i < sh.coords.size(); ++i) {
      if (sh.dcoef[i] == 0.0) continue;
      dc[sh.coords[i]] += sh.dcoef[i];
      for (unsigned c : S.coord_children[sh.coords[i]]) {
        if (_stamp[c] != _stamp_id) {
          _stamp[c] = _stamp_id;
          sh.children.push_back(c);
        }
      }
    }
    Vec &dz = _dz_scratch;
    dz.assign(S.latents.size(), 0.0);
    for (std::size_t i = 0; i < sh.lats.size(); ++i) {
      if (sh.dlat[i] == 0.0) continue;
      dz[sh.lats[i]] += sh.dlat[i];
      for (unsigned c : S.latent_children[sh.lats[i]]) {
        if (_stamp[c] != _stamp_id) {
          _stamp[c] = _stamp_id;
          sh.children.push_back(c);
        }
      }
    }
    // physical slope changes (an owned bias label): their columns
    for (std::size_t i = 0; i < sh.slopes.size(); ++i) {
      if (sh.dslope[i] == 0.0) continue;
      for (unsigned c : S.bias->slope[sh.slopes[i]].children) {
        if (_stamp[c] != _stamp_id) {
          _stamp[c] = _stamp_id;
          sh.children.push_back(c);
        }
      }
    }
    // the children in the order they are first reached (a function of the
    // shift's coordinates and latents)
    std::size_t total = 0;
    for (unsigned c : sh.children) {
      sh.offset.push_back(static_cast<unsigned>(total));
      total += S.children[c].n;
    }
    sh.delta.assign(total, 0.0);
    sh.wdelta.assign(total, 0.0);
    sh.q1.assign(sh.children.size(), 0.0);
    sh.q2.assign(sh.children.size(), 0.0);
    sh.value.assign(sh.children.size(), 0.0);
    for (std::size_t q = 0; q < sh.children.size(); ++q) {
      ChildSpec const &child = S.children[sh.children[q]];
      double *delta = &sh.delta[sh.offset[q]];
      for (unsigned i = 0; i < child.n; ++i) {
        const unsigned row = child.row0 + i;
        long double v = 0;
        for (unsigned k = 0; k < p; ++k) {
          if (dc[k] != 0.0) v += static_cast<long double>(S.X[row * p + k]) * dc[k];
        }
        if (child.latent >= 0) v += static_cast<long double>(_d[row]) * dz[child.latent];
        for (std::size_t b = 0; b < sh.slopes.size(); ++b) {
          v += static_cast<long double>(S.bias->slope[sh.slopes[b]].column[row]) * sh.dslope[b];
        }
        delta[i] = static_cast<double>(v);
      }
      double const *W = &_W[child.woff];
      double *wdelta = &sh.wdelta[sh.offset[q]];
      long double q1 = 0, q2 = 0;
      for (unsigned i = 0; i < child.n; ++i) {
        long double v = 0;
        for (unsigned j = 0; j < child.n; ++j) v += static_cast<long double>(W[i * child.n + j]) * delta[j];
        wdelta[i] = static_cast<double>(v);
        q1 += static_cast<long double>(delta[i]) * _wr[child.row0 + i];
        q2 += static_cast<long double>(delta[i]) * v;
      }
      sh.q1[q] = static_cast<double>(q1);
      sh.q2[q] = static_cast<double>(q2);
    }
    // priors of shifted active slabs (prior_coords) and latents
    long double p1 = 0, p2 = 0;
    std::vector<char> done(S.groups.size(), 0);
    for (unsigned k : prior_coords) {
      const unsigned j = S.coord_group[k];
      if (done[j]) continue;
      done[j] = 1;
      GroupSpec const &g = S.groups[j];
      for (unsigned i = 0; i < g.d; ++i) {
        const double vi = dc[g.coord0 + i];
        for (unsigned l = 0; l < g.d; ++l) {
          const double bil = g.B[i * g.d + l];
          p1 -= static_cast<long double>(vi) * bil * (_beta[g.coord0 + l] - g.m0[l]);
          p2 += static_cast<long double>(vi) * bil * dc[g.coord0 + l];
        }
      }
    }
    for (std::size_t lat = 0; lat < dz.size(); ++lat) {
      if (dz[lat] == 0.0) continue;
      p1 -= static_cast<long double>(dz[lat]) * _z[lat];
      p2 += static_cast<long double>(dz[lat]) * dz[lat];
    }
    sh.p1 = static_cast<double>(p1);
    sh.p2 = static_cast<double>(p2);
  }

  // The change of the inner target's likelihood terms at t times the shift
  // (each selection child evaluated, each Gaussian child exactly from its
  // quadratic form); the new values are kept in sh.value.
  long double children_delta(Shift &sh, double t)
  {
    Shared const &S = s();
    long double out = 0;
    for (std::size_t q = 0; q < sh.children.size(); ++q) {
      const unsigned c = sh.children[q];
      if (gaussian_child(c)) {
        out += static_cast<long double>(t) * sh.q1[q] - 0.5L * t * t * sh.q2[q];
        continue;
      }
      ChildSpec const &child = S.children[c];
      ensure_value(c);
      _scratch.resize(child.n);
      double const *delta = &sh.delta[sh.offset[q]];
      for (unsigned i = 0; i < child.n; ++i) _scratch[i] = _mean[child.row0 + i] + t * delta[i];
      const double v = selected_value(c, _scratch.data());
      sh.value[q] = v;
      if (v == negative_infinity) return negative_infinity;
      out += static_cast<long double>(v) - _value[c];
    }
    return out;
  }

  long double prior_delta(Shift const &sh, double t) const
  {
    return static_cast<long double>(t) * sh.p1 - 0.5L * t * t * sh.p2;
  }

  // Moves the state by t times the shift; the selection children take the
  // values of the last evaluation (children_delta at the same t).
  void apply_shift(Shift const &sh, double t)
  {
    Shared const &S = s();
    for (std::size_t i = 0; i < sh.coords.size(); ++i) {
      const unsigned k = sh.coords[i];
      _b[k] += t * sh.dcoef[i];
      if (active(S.coord_group[k])) _beta[k] = _b[k];
    }
    for (std::size_t i = 0; i < sh.lats.size(); ++i) _z[sh.lats[i]] += t * sh.dlat[i];
    for (std::size_t q = 0; q < sh.children.size(); ++q) {
      const unsigned c = sh.children[q];
      ChildSpec const &child = S.children[c];
      double const *delta = &sh.delta[sh.offset[q]];
      double const *wdelta = &sh.wdelta[sh.offset[q]];
      for (unsigned i = 0; i < child.n; ++i) {
        const unsigned row = child.row0 + i;
        _mean[row] += t * delta[i];
        _r[row] -= t * delta[i];
        _wr[row] -= t * wdelta[i];
      }
      if (!gaussian_child(c)) _value[c] = sh.value[q];
      _fit_ok[c] = 0;
      _dirty[c] = 1;
    }
  }

  // ---- slab priors

  double log_slab(GroupSpec const &g, double const *beta) const
  {
    return log_normal_precision(beta, g.m0, g.Bchol, g.log_det_B);
  }

  void draw_prior(GroupSpec const &g, RNG *rng, double *out)
  {
    draw_precision(g.m0, g.Bchol, rng, out);
  }

  bool all_finite(double const *x, std::size_t n) const
  {
    for (std::size_t i = 0; i < n; ++i) if (!std::isfinite(x[i])) return false;
    return true;
  }

  bool metropolis(RNG *rng, double log_ratio)
  {
    if (std::isnan(log_ratio)) return false;
    return log_ratio >= 0 || std::log(rng->uniform()) < log_ratio;
  }

  // A state-independent coin per attempt: true with the given probability
  // (the test switch coin forces it each way without a draw).
  bool alternative(RNG *rng, double probability)
  {
    if (s().coin == 1) return true;
    if (s().coin == 2) return false;
    return rng->uniform() < probability;
  }

  // The first stage of a two-state accept/reject (gate flips, location swaps):
  // Barker's acceptance r / (1 + r) after the coin, otherwise min(1, r). The
  // ratio of the two directions' acceptances is r under both, so the later
  // stages are unchanged; the mixture cannot alternate deterministically
  // between two states (a Metropolis step at ratio one would, and an evenly
  // thinned chain would keep one of them).
  bool first_stage(RNG *rng, double log_ratio)
  {
    if (std::isnan(log_ratio)) return false;
    if (!alternative(rng, two_state_coin)) return metropolis(rng, log_ratio);
    const double barker = log_ratio >= 0 ? 1.0 / (1.0 + std::exp(-log_ratio)) :
      std::exp(log_ratio) / (1.0 + std::exp(log_ratio));
    return rng->uniform() < barker;
  }

  // ---- component kernels

  // Dormant slab of an off group: an exact draw from its prior.
  void dormant_refresh(unsigned j, RNG *rng)
  {
    GroupSpec const &g = s().groups[j];
    Vec draw(g.d);
    draw_prior(g, rng, draw.data());
    if (!all_finite(draw.data(), g.d)) return;
    _repr_delta += log_slab(g, draw.data()) - log_slab(g, &_beta[g.coord0]);
    std::copy(draw.begin(), draw.end(), _beta.begin() + g.coord0);
  }

  // F_j in a Gaussian label: the collapsed flip along the sheared line with a
  // two-stage (factorized) acceptance.
  void gaussian_flip(unsigned j, RNG *rng, bool dry = false)
  {
    Shared const &S = s();
    GroupSpec const &g = S.groups[j];
    Line L;
    if (!build_line({j}, L)) { skipped_factorization(); return; }
    const bool on = active(j);
    const Component type = on ? flip_off : flip_on;
    // Savage-Dickey form of the collapsed odds of the line
    Vec zero(g.d, 0.0);
    const double log_odds = (g.log_on - g.log_off) + log_slab(g, zero.data()) -
      log_normal_precision(zero.data(), L.m, L.Pchol, L.log_det_P);
    const double stage1 = on ? -log_odds : log_odds;
    if (!dry) {
      count(counters().proposed[type]);
      count(counters().stage0[type]);
      // in a Gaussian label every proposal reaches stage 1 (no screen)
      evidence_proposed(j, on);
      evidence_odds(j, stage1);
      if (!first_stage(rng, stage1)) return;
      count(counters().stage1[type]);
    }
    Vec proposal(g.d), db(g.d);
    double const *old_beta = &_beta[g.coord0];
    long double proposal_terms;
    if (on) {
      draw_prior(g, rng, proposal.data());
      for (unsigned i = 0; i < g.d; ++i) db[i] = -old_beta[i];
      // reverse: the line law at the current slab; forward: the dormant prior
      proposal_terms = (S.mutant == 2 ? 0.0 : log_normal_precision(old_beta, L.m, L.Pchol, L.log_det_P)) -
        log_slab(g, proposal.data());
    } else {
      draw_precision(L.m, L.Pchol, rng, proposal.data());
      for (unsigned i = 0; i < g.d; ++i) db[i] = proposal[i];
      proposal_terms = (S.mutant == 2 ? 0.0 : log_slab(g, old_beta)) -
        log_normal_precision(proposal.data(), L.m, L.Pchol, L.log_det_P);
    }
    if (!all_finite(proposal.data(), g.d)) return;
    Shift sh;
    line_shift(L, db.data(), sh, (S.mutant == 1 && on) ? 1.5 : 1.0);
    prepare_shift(sh, L.partners);
    const long double children = children_delta(sh, 1.0);
    const long double rest = prior_delta(sh, 1.0) +
      (on ? g.log_off - g.log_on : g.log_on - g.log_off) +
      log_slab(g, proposal.data()) - log_slab(g, old_beta);
    const long double target = children + rest;
    const double full = static_cast<double>(target + proposal_terms);
    if (dry) {
      probe_record(type, stage1, gaussian_part(sh) + rest + proposal_terms, children - gaussian_part(sh), log_odds);
      return;
    }
    const double stage2 = full - stage1;
    if (!_selection) record_max(counters().max_stage2_gaussian, std::fabs(stage2));
    if (!std::isfinite(full) || !metropolis(rng, stage2)) return;
    count(counters().accepted[type]);
    evidence_accepted(j, on);
    apply_shift(sh, 1.0);
    _g[j] = on ? 0 : 1;
    std::copy(proposal.begin(), proposal.end(), _beta.begin() + g.coord0);
    for (unsigned i = 0; i < g.d; ++i) _b[g.coord0 + i] = on ? 0.0 : proposal[i];
    _repr_delta += target;
  }

  // The Gaussian part of the children terms of a shift at t = 1.
  long double gaussian_part(Shift const &sh) const
  {
    long double out = 0;
    for (std::size_t q = 0; q < sh.children.size(); ++q) out += sh.q1[q] - 0.5L * sh.q2[q];
    return out;
  }

  // Test hook (probe): one dry-run record of a collapsed Gaussian move: its
  // type, label, stage-1 log ratio, the full log ratio from the direct
  // Gaussian evaluation of priors, children and proposal densities (equal to
  // the stage-1 ratio when the line algebra is right), the selection part of
  // the inner target's ratio (log A~(old) - log A~(new)) and the collapsed log
  // odds.
  void probe_record(int type, double stage1, long double gaussian_full, long double selection, double log_odds)
  {
    ProbeStore &store = probe_store();
    const std::lock_guard<std::mutex> lock(store.mutex);
    store.records.insert(store.records.end(), {static_cast<double>(type), _selection ? 1.0 : 0.0, stage1,
      static_cast<double>(gaussian_full), static_cast<double>(selection), log_odds});
  }

  // R_j of an active group in a Gaussian label: an exact draw along the line.
  void gaussian_refresh(unsigned j, RNG *rng, bool dry = false)
  {
    GroupSpec const &g = s().groups[j];
    Line L;
    if (!build_line({j}, L)) { skipped_factorization(); return; }
    if (!dry) count(counters().proposed[refresh]);
    Vec proposal(g.d), db(g.d);
    draw_precision(L.m, L.Pchol, rng, proposal.data());
    if (!all_finite(proposal.data(), g.d)) return;
    double const *old_beta = &_beta[g.coord0];
    for (unsigned i = 0; i < g.d; ++i) db[i] = proposal[i] - old_beta[i];
    Shift sh;
    line_shift(L, db.data(), sh);
    // the moving group's slab is a shifted active slab too
    std::vector<unsigned> prior_coords(L.partners);
    prior_coords.push_back(g.coord0);
    prepare_shift(sh, prior_coords);
    const long double children = children_delta(sh, 1.0);
    const long double proposals = log_normal_precision(old_beta, L.m, L.Pchol, L.log_det_P) -
      log_normal_precision(proposal.data(), L.m, L.Pchol, L.log_det_P);
    const long double target = children + prior_delta(sh, 1.0);
    const double full = static_cast<double>(target + proposals);
    if (dry) {
      probe_record(refresh, 0.0, gaussian_part(sh) + prior_delta(sh, 1.0) + proposals,
                   children - gaussian_part(sh), 0.0);
      return;
    }
    if (!_selection) record_max(counters().max_stage2_gaussian, std::fabs(full));
    if (!std::isfinite(full) || !metropolis(rng, full)) return;
    count(counters().accepted[refresh]);
    apply_shift(sh, 1.0);
    _repr_delta += target;
  }

  // A univariate slice update of t -> pi~(state + t shift) (spec: Zs and the
  // stage-1 R_j in selection labels); stepping out and shrinkage as in
  // BayesTools' factor slice kernel.
  bool slice_shift(Shift &sh, double width, RNG *rng)
  {
    std::vector<double> u(1, 0.0), direction(1, 1.0);
    auto at = [&](std::vector<double> const &point) -> double {
      const double t = point[0];
      if (!std::isfinite(t)) return negative_infinity;
      const long double children = children_delta(sh, t);
      if (!(children > negative_infinity)) return negative_infinity;
      return static_cast<double>(children + prior_delta(sh, t));
    };
    const slice::SliceResult result = slice::slice_along(
      u, direction, width, 0.0, at, [rng]() { return rng->uniform(); },
      [rng]() { return rng->exponential(); });
    if (!result.accepted) return false;
    apply_shift(sh, u[0]);
    _repr_delta += result.log_density;
    return true;
  }

  // R_j of an active group in a selection label (stage 1): a slice update of
  // each coordinate along its sheared line, width 2.5 line SDs. A composite
  // component: its coordinates run in reverse order in the reversed sweep
  // (sweep condition 3), so that the whole sweep stays a palindrome.
  void slice_refresh(unsigned j, RNG *rng)
  {
    GroupSpec const &g = s().groups[j];
    Line L;
    if (!build_line({j}, L)) { skipped_factorization(); return; }
    const std::size_t nm = L.moving.size();
    std::vector<unsigned> prior_coords(L.partners);
    prior_coords.push_back(g.coord0);
    for (std::size_t step = 0; step < nm; ++step) {
      // mutant 13: the coordinates in the forward order in both directions
      const std::size_t i = (_backward && s().mutant != 13) ? nm - 1 - step : step;
      sequence_step(slot_coordinate, j, static_cast<unsigned>(i));
      count(counters().proposed[refresh]);
      Vec db(nm, 0.0);
      db[i] = 1.0;
      Shift sh;
      line_shift(L, db.data(), sh);
      prepare_shift(sh, prior_coords);
      if (slice_shift(sh, slice_width * std::sqrt(L.Sigma[i * nm + i]), rng)) count(counters().accepted[refresh]);
    }
  }

  // ---- stage 2: the selection tilt (spec: Selection tilt)

  // ell at a change db of the line's moving coordinates from the current
  // state: the inner target's terms of the children the line moves minus their
  // Gaussian part, i.e. -log A~ + const (zero for Gaussian children).
  // Only the Newton-Laplace fits evaluate ell: each cluster or factor child at
  // its coarsest rung while the coarse fits are on (A1), every other selection
  // child as the inner target does; the change from the current state, so
  // that the terms of children the shift leaves alone are zero.
  double line_ell(Line const &L, double const *db, Shift &sh)
  {
    const Category previous = enter(cat_fit);
    line_shift(L, db, sh);
    prepare_shift(sh, {});
    long double out = 0;
    for (std::size_t q = 0; q < sh.children.size() && out > negative_infinity; ++q) {
      const unsigned c = sh.children[q];
      if (gaussian_child(c)) continue;
      ChildSpec const &child = s().children[c];
      _scratch.resize(child.n);
      double const *delta = &sh.delta[sh.offset[q]];
      for (unsigned i = 0; i < child.n; ++i) _scratch[i] = _mean[child.row0 + i] + delta[i];
      double v, base;
      if (coarse_child(c)) {
        ensure_fit_value(c);
        base = _fit_value[c];
        v = coarse_value(c, _scratch.data());
      } else {
        ensure_value(c);
        base = _value[c];
        v = selected_value(c, _scratch.data());
      }
      out = v == negative_infinity ? negative_infinity : out + (static_cast<long double>(v) - base);
    }
    enter(previous);
    if (!(out > negative_infinity)) return negative_infinity;
    // the selection children's terms minus their Gaussian part (zero for
    // Gaussian children, whose terms are exactly their quadratic forms)
    long double gaussian = 0;
    for (std::size_t q = 0; q < sh.children.size(); ++q) {
      if (!gaussian_child(sh.children[q])) gaussian += sh.q1[q] - 0.5L * sh.q2[q];
    }
    return static_cast<double>(out - gaussian);
  }

  // The change of the line's moving coordinates that puts the free ones
  // [first, first + x.size()) at the physical values x and the others at the
  // physical values of `target`.
  void subline_db(Line const &L, std::size_t first, Vec const &x, Vec const &target, Vec &db) const
  {
    db.resize(L.moving.size());
    for (std::size_t m = 0; m < L.moving.size(); ++m) {
      const double value = (m >= first && m < first + x.size()) ? x[m - first] : target[m];
      db[m] = value - _b[L.moving[m]];
    }
  }

  struct Fit {
    bool ok = false;
    Vec xhat, Pt, Ptchol;
    double log_det_Pt = 0;
    double ell_hat = 0;
  };

  // The tilt of a line in closed form (A2): when every selection child the
  // line moves is a step-switch row, ell(x) = -sum_rows power log A_row(mean
  // at x) with A_row the row's selected-normal normalizer, and its gradient
  // and Hessian in the free coordinates V' g and V' diag(h) V from the rows'
  // closed-form derivatives in one pass (V the rows' mean changes per unit of
  // each moving coordinate). Otherwise the fits use central differences of
  // line_ell. Both are deterministic functions of the line.
  struct Tilt {
    bool analytic = false;
    std::size_t nm = 0;
    std::vector<unsigned> kids, rows;
    Vec V;                         // rows x moving coordinates
  };

  void make_tilt(Line const &L, Tilt &tilt, std::vector<Shift> const *given = nullptr)
  {
    Shared const &S = s();
    tilt = Tilt();
    tilt.nm = L.moving.size();
    std::vector<Shift> computed;
    if (given == nullptr) line_basis(L, computed);
    std::vector<Shift> const &basis = given != nullptr ? *given : computed;
    bool ok = true;
    for (std::size_t m = 0; m < tilt.nm && ok; ++m) {
      Shift const &b = basis[m];
      for (std::size_t q = 0; q < b.children.size(); ++q) {
        const unsigned c = b.children[q];
        if (gaussian_child(c)) continue;
        ChildSpec const &child = S.children[c];
        if (child.type != ChildType::step_switch || !_step[c].valid) {
          ok = false;
          break;
        }
        int &position = _tilt_position[c];
        if (position < 0) {
          position = static_cast<int>(tilt.kids.size());
          tilt.kids.push_back(c);
          tilt.rows.push_back(child.row0);
          tilt.V.resize(tilt.V.size() + tilt.nm, 0.0);
        }
        tilt.V[static_cast<std::size_t>(position) * tilt.nm + m] = b.delta[b.offset[q]];
      }
    }
    for (unsigned c : tilt.kids) _tilt_position[c] = -1;
    tilt.analytic = ok;
    if (!ok) {
      tilt.kids.clear();
      tilt.rows.clear();
      tilt.V.clear();
    }
  }

  // ell (absolute) at a point of the subline of the free coordinates [first,
  // first + d), with its gradient and Hessian there when asked.
  bool tilt_at(Tilt const &tilt, Line const &L, std::size_t first, Vec const &point, Vec const &target,
               double &ell, Vec *grad, Vec *hess)
  {
    const std::size_t d = point.size(), nm = tilt.nm;
    Vec db;
    subline_db(L, first, point, target, db);
    if (grad != nullptr) {
      grad->assign(d, 0.0);
      hess->assign(d * d, 0.0);
    }
    long double value = 0;
    for (std::size_t r = 0; r < tilt.kids.size(); ++r) {
      const unsigned c = tilt.kids[r];
      double const *v = &tilt.V[r * nm];
      long double mean = _mean[tilt.rows[r]];
      for (std::size_t m = 0; m < nm; ++m) mean += static_cast<long double>(db[m]) * v[m];
      double log_a, d1, d2;
      count(counters().evaluations[_category]);
      if (!step_row_log_normalizer(step_law(c), static_cast<double>(mean), &log_a, &d1, &d2)) return false;
      const double power = _step[c].power;
      value -= static_cast<long double>(power) * log_a;
      if (grad == nullptr) continue;
      for (std::size_t k = 0; k < d; ++k) {
        const double ak = v[first + k];
        if (ak == 0.0) continue;
        (*grad)[k] -= power * d1 * ak;
        for (std::size_t l = 0; l < d; ++l) (*hess)[k * d + l] -= power * d2 * ak * v[first + l];
      }
    }
    ell = static_cast<double>(value);
    return std::isfinite(ell);
  }

  // ell at a point of a subline by the tilt's own evaluator.
  double fit_ell(Tilt const &tilt, Line const &L, std::size_t first, Vec const &point, Vec const &target)
  {
    if (tilt.analytic) {
      double ell;
      return tilt_at(tilt, L, first, point, target, ell, nullptr, nullptr) ? ell : negative_infinity;
    }
    Vec db;
    Shift sh;
    subline_db(L, first, point, target, db);
    return line_ell(L, db.data(), sh);
  }

  // The damped Newton fit with the closed-form tilt (A2): the same iteration
  // as newton() with exact derivatives at every accepted point.
  Fit newton_analytic(Tilt const &tilt, Line const &L, std::size_t first, Vec const &target, Vec const &m,
                      Vec const &P, Vec const &Pchol, double log_det_P)
  {
    const std::size_t d = m.size();
    Fit fit;
    Vec x(m), grad, hess;
    double e;
    if (!tilt_at(tilt, L, first, x, target, e, &grad, &hess)) return fit;
    auto f_at = [&](Vec const &point, double ell) {
      return log_normal_precision(point.data(), m, Pchol, log_det_P) + ell;
    };
    for (int iteration = 0; iteration < 2; ++iteration) {
      Vec Pt;
      if (!floor_precision(P, Pchol, hess, d, Pt)) return fit;
      Vec step(d);
      for (std::size_t i = 0; i < d; ++i) {
        long double g = grad[i];
        for (std::size_t k = 0; k < d; ++k) g -= static_cast<long double>(P[i * d + k]) * (x[k] - m[k]);
        step[i] = static_cast<double>(g);
      }
      Vec chol(Pt);
      if (!cholesky(chol, d)) return fit;
      solve_lower(chol, d, step.data());
      solve_lower_transpose(chol, d, step.data());
      const double f0 = f_at(x, e);
      const double lowest = f0 - fit_tolerance * (1.0 + std::fabs(f0));
      bool moved = false;
      double scale = 1.0;
      for (int halving = 0; halving <= 5 && !moved; ++halving, scale *= 0.5) {
        Vec next(x), g_next, h_next;
        for (std::size_t i = 0; i < d; ++i) next[i] += scale * step[i];
        double en;
        if (tilt_at(tilt, L, first, next, target, en, &g_next, &h_next) && f_at(next, en) >= lowest) {
          x = next;
          e = en;
          grad = g_next;
          hess = h_next;
          moved = true;
        }
      }
      if (!moved) break;
    }
    if (!floor_precision(P, Pchol, hess, d, fit.Pt)) return fit;
    fit.Ptchol = fit.Pt;
    if (!cholesky(fit.Ptchol, d)) return fit;
    fit.log_det_Pt = log_det_cholesky(fit.Ptchol, d);
    fit.xhat = x;
    fit.ell_hat = e;
    fit.ok = std::isfinite(fit.log_det_Pt) && all_finite(x.data(), d);
    return fit;
  }

  // The damped Newton fit of f(x) = log N(x; m, P^-1) + ell(x) along the
  // subline of the free coordinates (spec: two iterations from m, curvature
  // floored at 0.25 P, at most five halvings, central differences with steps
  // of a quarter of the law's SD; the full Hessian up to four coordinates and
  // its diagonal beyond), or with the closed-form tilt when it applies. A
  // deterministic function of the line.
  Fit newton(Tilt const &tilt, Line const &L, std::size_t first, Vec const &target, Vec const &m, Vec const &P,
             Vec const &Pchol, double log_det_P)
  {
    if (tilt.analytic) return newton_analytic(tilt, L, first, target, m, P, Pchol, log_det_P);
    const std::size_t d = m.size();
    Fit fit;
    Vec db, x(m), h(d);
    Shift sh;
    const Vec Sigma = inverse_from_cholesky(Pchol, d);
    for (std::size_t i = 0; i < d; ++i) h[i] = 0.25 * std::sqrt(Sigma[i * d + i]);
    auto ell_at = [&](Vec const &point) {
      subline_db(L, first, point, target, db);
      return line_ell(L, db.data(), sh);
    };
    auto f_at = [&](Vec const &point, double ell) {
      return log_normal_precision(point.data(), m, Pchol, log_det_P) + ell;
    };
    auto derivatives = [&](Vec const &point, double e0, Vec &grad, Vec &hess) {
      grad.assign(d, 0.0);
      hess.assign(d * d, 0.0);
      Vec y(point), plus(d), minus(d);
      for (std::size_t i = 0; i < d; ++i) {
        y[i] = point[i] + h[i];
        plus[i] = ell_at(y);
        y[i] = point[i] - h[i];
        minus[i] = ell_at(y);
        y[i] = point[i];
        if (!std::isfinite(plus[i]) || !std::isfinite(minus[i])) return false;
        grad[i] = (plus[i] - minus[i]) / (2.0 * h[i]);
        hess[i * d + i] = (plus[i] - 2.0 * e0 + minus[i]) / (h[i] * h[i]);
      }
      if (d <= 4) {
        for (std::size_t i = 0; i < d; ++i) {
          for (std::size_t k = i + 1; k < d; ++k) {
            double corners[4];
            const double si[4] = {1, 1, -1, -1}, sk[4] = {1, -1, 1, -1};
            for (int c = 0; c < 4; ++c) {
              y[i] = point[i] + si[c] * h[i];
              y[k] = point[k] + sk[c] * h[k];
              corners[c] = ell_at(y);
              if (!std::isfinite(corners[c])) return false;
            }
            y[i] = point[i];
            y[k] = point[k];
            const double mixed = (corners[0] - corners[1] - corners[2] + corners[3]) / (4.0 * h[i] * h[k]);
            hess[i * d + k] = mixed;
            hess[k * d + i] = mixed;
          }
        }
      }
      return true;
    };
    double e = ell_at(x);
    if (!std::isfinite(e)) return fit;
    for (int iteration = 0; iteration < 2; ++iteration) {
      Vec grad, hess, Pt;
      if (!derivatives(x, e, grad, hess) || !floor_precision(P, Pchol, hess, d, Pt)) return fit;
      Vec step(d);
      for (std::size_t i = 0; i < d; ++i) {
        long double g = grad[i];
        for (std::size_t k = 0; k < d; ++k) g -= static_cast<long double>(P[i * d + k]) * (x[k] - m[k]);
        step[i] = static_cast<double>(g);
      }
      Vec chol(Pt);
      if (!cholesky(chol, d)) return fit;
      solve_lower(chol, d, step.data());
      solve_lower_transpose(chol, d, step.data());
      const double f0 = f_at(x, e);
      const double lowest = f0 - fit_tolerance * (1.0 + std::fabs(f0));
      bool moved = false;
      double scale = 1.0;
      for (int halving = 0; halving <= 5 && !moved; ++halving, scale *= 0.5) {
        Vec next(x);
        for (std::size_t i = 0; i < d; ++i) next[i] += scale * step[i];
        const double en = ell_at(next);
        if (std::isfinite(en) && f_at(next, en) >= lowest) {
          x = next;
          e = en;
          moved = true;
        }
      }
      if (!moved) break;
    }
    Vec grad, hess;
    if (!derivatives(x, e, grad, hess) || !floor_precision(P, Pchol, hess, d, fit.Pt)) return fit;
    fit.Ptchol = fit.Pt;
    if (!cholesky(fit.Ptchol, d)) return fit;
    fit.log_det_Pt = log_det_cholesky(fit.Ptchol, d);
    fit.xhat = x;
    fit.ell_hat = e;
    fit.ok = std::isfinite(fit.log_det_Pt) && all_finite(x.data(), d);
    return fit;
  }

  // The log of the tilt factor of a Laplace evidence along a subline:
  // log int N(x; m, P^-1) exp(ell(x) - ell(zero)) dx by the fit.
  double laplace_tilt(Fit const &fit, Vec const &m, Vec const &Pchol, double log_det_P, double ell_zero) const
  {
    const double d = static_cast<double>(m.size());
    return log_normal_precision(fit.xhat.data(), m, Pchol, log_det_P) + fit.ell_hat - ell_zero +
      0.5 * d * log_two_pi - 0.5 * fit.log_det_Pt;
  }

  // The row-marginal parameters of row i of a selection child: marginal SD,
  // selection SE, observed bin (0: exempt), sign, power, weights and bounds.
  struct RowModel {
    double sd, se, power;
    int bin, sign, bins;
    double const *omega;
    Vec lower, upper;
  };

  bool row_model(ChildSpec const &c, unsigned i, RowModel &out)
  {
    auto const &par = c.node->parents();
    auto bounds = [&](unsigned lo, unsigned hi, unsigned om) {
      const unsigned long nb = c.lengths[lo];
      const SelNormJagsBounds lower(par[lo]->value(_chain), static_cast<unsigned>(nb));
      const SelNormJagsBounds upper(par[hi]->value(_chain), static_cast<unsigned>(nb));
      out.lower.assign(lower.data(), lower.data() + nb);
      out.upper.assign(upper.data(), upper.data() + nb);
      out.omega = par[om]->value(_chain);
      out.bins = static_cast<int>(nb);
    };
    switch (c.type) {
    case ChildType::normal:
      return false;
    case ChildType::step_switch:
      out.sd = par[1]->value(_chain)[0];
      out.se = par[2]->value(_chain)[0];
      out.power = par[3]->value(_chain)[0];
      out.bin = static_cast<int>(par[7]->value(_chain)[0]);
      out.sign = static_cast<int>(par[8]->value(_chain)[0]);
      bounds(5, 6, 4);
      return true;
    case ChildType::cluster:
    case ChildType::factor: {
      const double rsd = par[1]->value(_chain)[i];
      double const *loading = par[2]->value(_chain);
      const unsigned r = c.type == ChildType::cluster ? 1 : static_cast<unsigned>(c.lengths[2] / c.n);
      long double variance = static_cast<long double>(rsd) * rsd;
      for (unsigned k = 0; k < r; ++k) variance += static_cast<long double>(loading[i + c.n * k]) * loading[i + c.n * k];
      out.sd = std::sqrt(static_cast<double>(variance));
      out.se = par[3]->value(_chain)[i];
      out.power = 1.0;
      out.bin = static_cast<int>(par[7]->value(_chain)[i]);
      out.sign = static_cast<int>(par[8]->value(_chain)[0]);
      bounds(5, 6, 4);
      return true;
    }
    case ChildType::mv: {
      double const *packed = par[1]->value(_chain);
      std::size_t position = 0;
      for (unsigned col = 0; col < i; ++col) position += c.n - col;
      out.sd = std::sqrt(packed[position]);
      out.se = par[2]->value(_chain)[i];
      out.power = 1.0;
      out.bin = static_cast<int>(par[6]->value(_chain)[i]);
      out.sign = static_cast<int>(par[7]->value(_chain)[0]);
      bounds(4, 5, 3);
      return true;
    }
    }
    return false;
  }

  // Unit shifts of each moving coordinate of a line (with the shear), prepared.
  void line_basis(Line const &L, std::vector<Shift> &basis)
  {
    const std::size_t nm = L.moving.size();
    basis.assign(nm, Shift());
    for (std::size_t m = 0; m < nm; ++m) {
      Vec unit(nm, 0.0);
      unit[m] = 1.0;
      line_shift(L, unit.data(), basis[m]);
      prepare_shift(basis[m], {});
    }
  }

  // ell0 of the row-marginal model, -sum_rows power log A_row, at a change db
  // of the line's moving coordinates, with its gradient and Hessian in the
  // free coordinates [first, first + d): closed form, no surrogate evaluation.
  bool row_marginal(std::vector<Shift> const &basis, Vec const &db, std::size_t first, std::size_t d,
                    double *value, Vec &grad, Vec &hess)
  {
    Shared const &S = s();
    const std::size_t nm = basis.size();
    // the selection children's rows x moving coordinates (the positions in
    // _tilt_position are reset before returning)
    std::vector<unsigned> kids;
    std::vector<std::size_t> offset;
    Vec deltas;
    for (std::size_t m = 0; m < nm; ++m) {
      Shift const &b = basis[m];
      for (std::size_t q = 0; q < b.children.size(); ++q) {
        const unsigned c = b.children[q];
        if (gaussian_child(c)) continue;
        int &position = _tilt_position[c];
        if (position < 0) {
          position = static_cast<int>(kids.size());
          kids.push_back(c);
          offset.push_back(deltas.size());
          deltas.resize(deltas.size() + S.children[c].n * nm, 0.0);
        }
        double *rows = &deltas[offset[static_cast<std::size_t>(position)]];
        for (unsigned i = 0; i < S.children[c].n; ++i) rows[i * nm + m] = b.delta[b.offset[q] + i];
      }
    }
    for (unsigned c : kids) _tilt_position[c] = -1;
    long double total = 0;
    grad.assign(d, 0.0);
    hess.assign(d * d, 0.0);
    RowModel row;
    for (std::size_t k = 0; k < kids.size(); ++k) {
      ChildSpec const &c = S.children[kids[k]];
      double const *rows = &deltas[offset[k]];
      for (unsigned i = 0; i < c.n; ++i) {
        long double mu = _mean[c.row0 + i];
        for (std::size_t m = 0; m < nm; ++m) mu += static_cast<long double>(db[m]) * rows[i * nm + m];
        double log_a, d1, d2, power;
        if (c.type == ChildType::step_switch && _step[kids[k]].valid) {
          // a step-switch row's own normalizer, its parameters read once per update
          StepChild const &st = _step[kids[k]];
          if (st.bin == 0) continue;
          power = st.power;
          if (!step_row_log_normalizer(step_law(kids[k]), static_cast<double>(mu), &log_a, &d1, &d2)) return false;
        } else {
          if (!row_model(c, i, row) || row.bin == 0) continue;
          power = row.power;
          if (!row_log_normalizer(static_cast<double>(mu), row.sd, row.se, row.sign, row.omega,
                                  row.lower.data(), row.upper.data(), row.bins, &log_a, &d1, &d2)) return false;
        }
        total -= static_cast<long double>(power) * log_a;
        for (std::size_t a = 0; a < d; ++a) {
          const double ak = rows[i * nm + first + a];
          if (ak == 0.0) continue;
          grad[a] -= power * d1 * ak;
          for (std::size_t l = 0; l < d; ++l) hess[a * d + l] -= power * d2 * ak * rows[i * nm + first + l];
        }
      }
    }
    *value = static_cast<double>(total);
    return std::isfinite(*value);
  }

  // The screen's tilt factor along a subline: the Laplace integral of
  // N(x; m, P^-1) exp(ell0(x) - ell0(zero)) with ell0 expanded to second
  // order at m (closed form).
  double screen_tilt(std::vector<Shift> const &basis, Line const &L, std::size_t first, Vec const &target,
                     Vec const &zero_target, Vec const &m, Vec const &P, Vec const &Pchol, double log_det_P)
  {
    const std::size_t d = m.size();
    Vec db, grad, hess, g0, h0;
    double at_m, at_zero;
    subline_db(L, first, m, target, db);
    if (!row_marginal(basis, db, first, d, &at_m, grad, hess)) return 0.0;
    Vec zero(d, 0.0);
    subline_db(L, first, zero, zero_target, db);
    if (!row_marginal(basis, db, first, d, &at_zero, g0, h0)) return 0.0;
    Vec Pf;
    if (!floor_precision(P, Pchol, hess, d, Pf)) return 0.0;
    Vec chol(Pf);
    if (!cholesky(chol, d)) return 0.0;
    Vec solved(grad);
    solve_lower(chol, d, solved.data());
    long double quadratic = 0;
    for (double v : solved) quadratic += static_cast<long double>(v) * v;
    const double out = at_m - at_zero + 0.5 * log_det_P - 0.5 * log_det_cholesky(chol, d) +
      0.5 * static_cast<double>(quadratic);
    return std::isfinite(out) ? out : 0.0;
  }

  static double clamp_screen(double log_odds)
  {
    const double limit = std::log(20.0);
    return std::max(-limit, std::min(limit, log_odds));
  }

  // The line, Gaussian odds, Laplace fit and screen of one group, shared by
  // F_j and R_j within one sweep (they are adjacent in either order and move
  // along the same line).
  struct GroupLine {
    unsigned group = 0;
    bool valid = false;
    Line L;
    double log_odds_gauss = 0;   // collapsed Gaussian log odds (prior odds included)
    bool screened = false;
    double log_screen = 0;       // clamped log C0
    bool fitted = false;
    Fit fit;
    bool have_basis = false;     // the unit shifts of the line, shared by the screen and the fit
    std::vector<Shift> basis;
    double log_odds = 0;         // Laplace log odds (prior odds included)
  };
  GroupLine _line;

  bool group_line(unsigned j)
  {
    if (_line.valid && _line.group == j) return true;
    _line = GroupLine();
    _line.group = j;
    if (!build_line({j}, _line.L)) { skipped_factorization(); return false; }
    GroupSpec const &g = s().groups[j];
    Vec zero(g.d, 0.0);
    _line.log_odds_gauss = (g.log_on - g.log_off) + log_slab(g, zero.data()) -
      log_normal_precision(zero.data(), _line.L.m, _line.L.Pchol, _line.L.log_det_P);
    _line.valid = true;
    return true;
  }

  void group_screen()
  {
    if (_line.screened) return;
    if (!_line.have_basis) line_basis(_line.L, _line.basis);
    _line.have_basis = true;
    Vec target(_line.L.m), zero(_line.L.m.size(), 0.0);
    _line.log_screen = clamp_screen(_line.log_odds_gauss + screen_tilt(
      _line.basis, _line.L, 0, target, zero, _line.L.m, _line.L.P, _line.L.Pchol, _line.L.log_det_P));
    _line.screened = true;
  }

  void group_fit()
  {
    if (_line.fitted) return;
    const Category previous = enter(cat_fit);
    Line const &L = _line.L;
    Vec target(L.m);
    Tilt tilt;
    if (!_line.have_basis) line_basis(L, _line.basis);
    _line.have_basis = true;
    make_tilt(L, tilt, &_line.basis);
    _line.fit = newton(tilt, L, 0, target, L.m, L.P, L.Pchol, L.log_det_P);
    // ell at the physical zero of the group (the off state)
    Vec zero(L.m.size(), 0.0);
    const double ell_zero = fit_ell(tilt, L, 0, zero, zero);
    if (!_line.fit.ok || !std::isfinite(ell_zero)) {
      // a deterministic fallback of the line: the Gaussian law, its ell at m
      Fit fallback;
      fallback.xhat = L.m;
      fallback.Pt = L.P;
      fallback.Ptchol = L.Pchol;
      fallback.log_det_Pt = L.log_det_P;
      fallback.ell_hat = fit_ell(tilt, L, 0, L.m, L.m);
      fallback.ok = std::isfinite(fallback.ell_hat) && std::isfinite(ell_zero);
      _line.fit = fallback;
    }
    _line.log_odds = _line.fit.ok ?
      _line.log_odds_gauss + laplace_tilt(_line.fit, L.m, L.Pchol, L.log_det_P, ell_zero) : negative_infinity;
    _line.fitted = true;
    enter(previous);
  }

  // F_j in a selection label (stage 2): the tilted flip with the free screen
  // (stage 0), the Laplace odds (stage 1) and the full ratio (stage 2).
  void tilted_flip(unsigned j, RNG *rng)
  {
    Shared const &S = s();
    GroupSpec const &g = S.groups[j];
    if (!group_line(j)) return;
    const bool on = active(j);
    const Component type = on ? flip_off : flip_on;
    count(counters().proposed[type]);
    evidence_proposed(j, on);
    group_screen();
    const double stage0 = on ? -_line.log_screen : _line.log_screen;
    if (!first_stage(rng, stage0)) {
      evidence_screened(j);
      return;
    }
    count(counters().stage0[type]);
    group_fit();
    if (!_line.fit.ok) return;
    const double log_odds = _line.log_odds;
    evidence_odds(j, on ? -log_odds : log_odds);
    const double stage1 = on ? -(log_odds - _line.log_screen) : (log_odds - _line.log_screen);
    if (!metropolis(rng, stage1)) return;
    count(counters().stage1[type]);
    if (_adapting) ++_fit_proposed;
    Fit const &fit = _line.fit;
    Vec proposal(g.d), db(g.d);
    double const *old_beta = &_beta[g.coord0];
    long double proposal_terms;
    if (on) {
      draw_prior(g, rng, proposal.data());
      for (unsigned i = 0; i < g.d; ++i) db[i] = -old_beta[i];
      proposal_terms = (S.mutant == 2 ? 0.0 : log_normal_precision(old_beta, fit.xhat, fit.Ptchol, fit.log_det_Pt)) -
        log_slab(g, proposal.data());
    } else {
      draw_precision(fit.xhat, fit.Ptchol, rng, proposal.data());
      for (unsigned i = 0; i < g.d; ++i) db[i] = proposal[i];
      proposal_terms = (S.mutant == 2 ? 0.0 : log_slab(g, old_beta)) -
        log_normal_precision(proposal.data(), fit.xhat, fit.Ptchol, fit.log_det_Pt);
    }
    if (!all_finite(proposal.data(), g.d)) return;
    Shift sh;
    line_shift(_line.L, db.data(), sh, (S.mutant == 1 && on) ? 1.5 : 1.0);
    prepare_shift(sh, _line.L.partners);
    const long double target = children_delta(sh, 1.0) + prior_delta(sh, 1.0) +
      (on ? g.log_off - g.log_on : g.log_on - g.log_off) +
      log_slab(g, proposal.data()) - log_slab(g, old_beta);
    const double full = static_cast<double>(target + proposal_terms);
    const double stage2 = full - (on ? -log_odds : log_odds);
    if (!std::isfinite(full) || !metropolis(rng, stage2)) return;
    count(counters().accepted[type]);
    evidence_accepted(j, on);
    if (_adapting) ++_fit_accepted;
    apply_shift(sh, 1.0);
    _g[j] = on ? 0 : 1;
    std::copy(proposal.begin(), proposal.end(), _beta.begin() + g.coord0);
    for (unsigned i = 0; i < g.d; ++i) _b[g.coord0 + i] = on ? 0.0 : proposal[i];
    _repr_delta += target;
  }

  // R_j of an active group in a selection label (stage 2): an independence
  // proposal from the tilted law along the line.
  void tilted_refresh(unsigned j, RNG *rng)
  {
    GroupSpec const &g = s().groups[j];
    if (!group_line(j)) return;
    count(counters().proposed[refresh]);
    group_fit();
    if (!_line.fit.ok) return;
    Fit const &fit = _line.fit;
    Vec proposal(g.d), db(g.d);
    draw_precision(fit.xhat, fit.Ptchol, rng, proposal.data());
    if (!all_finite(proposal.data(), g.d)) return;
    double const *old_beta = &_beta[g.coord0];
    for (unsigned i = 0; i < g.d; ++i) db[i] = proposal[i] - old_beta[i];
    Shift sh;
    line_shift(_line.L, db.data(), sh);
    std::vector<unsigned> prior_coords(_line.L.partners);
    prior_coords.push_back(g.coord0);
    prepare_shift(sh, prior_coords);
    const long double target = children_delta(sh, 1.0) + prior_delta(sh, 1.0);
    const double full = static_cast<double>(target +
      log_normal_precision(old_beta, fit.xhat, fit.Ptchol, fit.log_det_Pt) -
      log_normal_precision(proposal.data(), fit.xhat, fit.Ptchol, fit.log_det_Pt));
    const bool accepted = std::isfinite(full) && metropolis(rng, full);
    if (_adapting) {
      ++_refresh_proposed[j];
      if (accepted) ++_refresh_accepted[j];
    }
    if (!accepted) return;
    count(counters().accepted[refresh]);
    apply_shift(sh, 1.0);
    _repr_delta += target;
  }

  // S_jk in a selection label (stage 2): the swap with the group switched on
  // drawn from its tilted law along the pair line with the other off; the
  // pattern evidences by their own Laplace fits, screened like a flip.
  void tilted_swap(unsigned first, unsigned second, RNG *rng)
  {
    Shared const &S = s();
    if (active(first) == active(second)) return;
    const unsigned j = active(first) ? first : second;   // on -> off
    const unsigned k = active(first) ? second : first;   // off -> on
    GroupSpec const &gj = S.groups[j];
    GroupSpec const &gk = S.groups[k];
    Line L;
    if (!build_line({j, k}, L)) { skipped_factorization(); return; }
    count(counters().proposed[swap]);
    const std::size_t dj = gj.d, dk = gk.d, nm = dj + dk;
    // the conditional Gaussian laws of each group with the other at zero, and
    // the Gaussian log ratio of the pattern evidences
    Vec cj, Cj, ck, Ck, mj, Pj, mk, Pk;
    double lcj, lck, lj, lk;
    if (!swap_laws(L, dj, dk, mj, Pj, lj, mk, Pk, lk, cj, Cj, lcj, ck, Ck, lck)) { skipped_factorization(); return; }
    Vec zero_j(dj, 0.0), zero_k(dk, 0.0);
    const double prior_odds = gj.log_off + gk.log_on - gj.log_on - gk.log_off;
    const double gauss = prior_odds + log_slab(gk, zero_k.data()) + log_normal_precision(zero_j.data(), mj, Pj, lj) -
      log_slab(gj, zero_j.data()) - log_normal_precision(zero_k.data(), mk, Pk, lk);
    // targets on the pair line: A = (j free, k at zero), B = (j at zero, k free)
    Vec target_a(nm, 0.0), target_b(nm, 0.0), zeros(nm, 0.0);
    for (std::size_t i = 0; i < dj; ++i) target_a[i] = cj[i];
    for (std::size_t i = 0; i < dk; ++i) target_b[dj + i] = ck[i];
    Vec Cjp = precision_from_cholesky(Cj, dj), Ckp = precision_from_cholesky(Ck, dk);
    std::vector<Shift> basis;
    line_basis(L, basis);
    const double screen = clamp_screen(gauss +
      screen_tilt(basis, L, dj, target_b, zeros, ck, Ckp, Ck, lck) -
      screen_tilt(basis, L, 0, target_a, zeros, cj, Cjp, Cj, lcj));
    if (!first_stage(rng, screen)) return;
    count(counters().stage0[swap]);
    const Category previous = enter(cat_fit);
    Tilt tilt;
    make_tilt(L, tilt, &basis);
    const Fit fa = newton(tilt, L, 0, target_a, cj, Cjp, Cj, lcj);
    const Fit fb = newton(tilt, L, dj, target_b, ck, Ckp, Ck, lck);
    const double ell_zero = fit_ell(tilt, L, 0, zero_j, zeros);
    enter(previous);
    if (!fa.ok || !fb.ok || !std::isfinite(ell_zero)) return;
    const double log_odds = gauss + laplace_tilt(fb, ck, Ck, lck, ell_zero) - laplace_tilt(fa, cj, Cj, lcj, ell_zero);
    if (!metropolis(rng, log_odds - screen)) return;
    count(counters().stage1[swap]);
    if (_adapting) ++_fit_proposed;
    Vec beta_k(dk), dormant_j(dj);
    draw_precision(fb.xhat, fb.Ptchol, rng, beta_k.data());
    draw_prior(gj, rng, dormant_j.data());
    if (!all_finite(beta_k.data(), dk) || !all_finite(dormant_j.data(), dj)) return;
    double const *old_j = &_beta[gj.coord0];
    double const *old_k = &_beta[gk.coord0];
    Vec change(nm);
    for (std::size_t i = 0; i < dj; ++i) change[i] = -old_j[i];
    for (std::size_t i = 0; i < dk; ++i) change[dj + i] = beta_k[i];
    Shift move;
    line_shift(L, change.data(), move);
    prepare_shift(move, L.partners);
    const long double target = children_delta(move, 1.0) + prior_delta(move, 1.0) + prior_odds +
      log_slab(gj, dormant_j.data()) - log_slab(gj, old_j) + log_slab(gk, beta_k.data()) - log_slab(gk, old_k);
    const long double proposals = log_normal_precision(old_j, fa.xhat, fa.Ptchol, fa.log_det_Pt) +
      log_slab(gk, old_k) - log_normal_precision(beta_k.data(), fb.xhat, fb.Ptchol, fb.log_det_Pt) -
      log_slab(gj, dormant_j.data());
    const double full = static_cast<double>(target + proposals);
    if (!std::isfinite(full) || !metropolis(rng, full - log_odds)) return;
    count(counters().accepted[swap]);
    if (_adapting) ++_fit_accepted;
    apply_shift(move, 1.0);
    _g[j] = 0;
    _g[k] = 1;
    std::copy(dormant_j.begin(), dormant_j.end(), _beta.begin() + gj.coord0);
    std::copy(beta_k.begin(), beta_k.end(), _beta.begin() + gk.coord0);
    for (std::size_t i = 0; i < dj; ++i) _b[gj.coord0 + i] = 0.0;
    for (std::size_t i = 0; i < dk; ++i) _b[gk.coord0 + i] = beta_k[i];
    _repr_delta += target;
  }

  // The precision matrix from its lower Cholesky factor.
  static Vec precision_from_cholesky(Vec const &chol, std::size_t n)
  {
    Vec out(n * n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
      for (std::size_t j = 0; j < n; ++j) {
        long double v = 0;
        for (std::size_t k = 0; k <= std::min(i, j); ++k) v += static_cast<long double>(chol[i * n + k]) * chol[j * n + k];
        out[i * n + j] = static_cast<double>(v);
      }
    }
    return out;
  }

  // The marginal and the conditional (other group at zero) Gaussian laws of
  // the two groups of a pair line, as means and precision Cholesky factors.
  bool swap_laws(Line const &L, std::size_t dj, std::size_t dk, Vec &mj, Vec &Pj, double &lj, Vec &mk, Vec &Pk,
                 double &lk, Vec &cj, Vec &Cj, double &lcj, Vec &ck, Vec &Ck, double &lck)
  {
    const std::size_t nm = dj + dk;
    auto block = [&](std::size_t r0, std::size_t n0, std::size_t c0, std::size_t n1) {
      Vec out(n0 * n1);
      for (std::size_t i = 0; i < n0; ++i) for (std::size_t l = 0; l < n1; ++l) out[i * n1 + l] = L.Sigma[(r0 + i) * nm + c0 + l];
      return out;
    };
    auto marginal = [&](std::size_t r0, std::size_t n0, Vec &m, Vec &Pchol, double &log_det) {
      Vec Sg = block(r0, n0, r0, n0);
      m.assign(L.m.begin() + r0, L.m.begin() + r0 + n0);
      if (!cholesky(Sg, n0)) return false;
      Pchol = inverse_from_cholesky(Sg, n0);
      if (!cholesky(Pchol, n0)) return false;
      log_det = log_det_cholesky(Pchol, n0);
      return true;
    };
    auto conditional = [&](std::size_t r0, std::size_t n0, std::size_t c0, std::size_t n1,
                           Vec &m, Vec &Pchol, double &log_det) {
      Vec Soo = block(c0, n1, c0, n1);
      if (!cholesky(Soo, n1)) return false;
      const Vec Sinv = inverse_from_cholesky(Soo, n1);
      const Vec Sro = block(r0, n0, c0, n1);
      Vec Srr = block(r0, n0, r0, n0);
      m.assign(n0, 0.0);
      for (std::size_t i = 0; i < n0; ++i) {
        long double mi = L.m[r0 + i];
        for (std::size_t a = 0; a < n1; ++a) {
          long double coef = 0;
          for (std::size_t b = 0; b < n1; ++b) coef += static_cast<long double>(Sro[i * n1 + b]) * Sinv[b * n1 + a];
          mi -= coef * L.m[c0 + a];
          for (std::size_t l = 0; l < n0; ++l) Srr[i * n0 + l] -= static_cast<double>(coef * Sro[l * n1 + a]);
        }
        m[i] = static_cast<double>(mi);
      }
      if (!cholesky(Srr, n0)) return false;
      Pchol = inverse_from_cholesky(Srr, n0);
      if (!cholesky(Pchol, n0)) return false;
      log_det = log_det_cholesky(Pchol, n0);
      return true;
    };
    return marginal(0, dj, mj, Pj, lj) && marginal(dj, dk, mk, Pk, lk) &&
      conditional(0, dj, dj, dk, cj, Cj, lcj) && conditional(dj, dk, 0, dj, ck, Ck, lck);
  }

  // Z_s in a Gaussian label: an exact draw from its conditional.
  void gaussian_latent(unsigned lat, RNG *rng, bool dry = false)
  {
    Shared const &S = s();
    if (!dry) count(counters().proposed[study]);
    long double score = static_cast<long double>(_a[lat]) * _z[lat];
    for (unsigned c : S.latent_children[lat]) {
      ChildSpec const &child = S.children[c];
      for (unsigned i = 0; i < child.n; ++i) score += static_cast<long double>(_d[child.row0 + i]) * _wr[child.row0 + i];
    }
    const double lambda = _lambda[lat];
    const double mean = static_cast<double>(score / lambda);
    const double sd = 1.0 / std::sqrt(lambda);
    const double proposal = mean + sd * rng->normal();
    if (!std::isfinite(proposal)) return;
    Shift sh;
    sh.lats.push_back(lat);
    sh.dlat.push_back(proposal - _z[lat]);
    prepare_shift(sh, {});
    const long double children = children_delta(sh, 1.0);
    const long double target = children + prior_delta(sh, 1.0);
    const double zo = (_z[lat] - mean) / sd, zn = (proposal - mean) / sd;
    const double full = static_cast<double>(target) - 0.5 * zo * zo + 0.5 * zn * zn;
    if (dry) {
      probe_record(study, 0.0, gaussian_part(sh) + prior_delta(sh, 1.0) - 0.5 * zo * zo + 0.5 * zn * zn,
                   children - gaussian_part(sh), 0.0);
      return;
    }
    if (!_selection) record_max(counters().max_stage2_gaussian, std::fabs(full));
    if (!std::isfinite(full) || !metropolis(rng, full)) return;
    count(counters().accepted[study]);
    apply_shift(sh, 1.0);
    _repr_delta += target;
  }

  void slice_latent(unsigned lat, RNG *rng)
  {
    count(counters().proposed[study]);
    Shift sh;
    sh.lats.push_back(lat);
    sh.dlat.push_back(1.0);
    prepare_shift(sh, {});
    if (slice_shift(sh, slice_width / std::sqrt(_lambda[lat]), rng)) count(counters().accepted[study]);
  }

  // J: all active coefficients from their joint Gaussian law with the latent
  // shear (Gaussian labels only).
  void gaussian_joint(RNG *rng, bool dry = false)
  {
    Shared const &S = s();
    const unsigned p = S.p;
    std::vector<unsigned> A;
    for (unsigned j = 0; j < S.groups.size(); ++j) {
      if (!active(j)) continue;
      for (unsigned i = 0; i < S.groups[j].d; ++i) A.push_back(S.groups[j].coord0 + i);
    }
    const std::size_t na = A.size();
    if (na == 0) return;
    if (!dry) count(counters().proposed[joint]);
    Vec PA(na * na), hA(na);
    for (std::size_t i = 0; i < na; ++i) {
      hA[i] = _h[A[i]];
      for (std::size_t k = 0; k < na; ++k) PA[i * na + k] = _P[A[i] * p + A[k]];
    }
    Vec chol(PA);
    if (!cholesky(chol, na)) { skipped_factorization(); return; }
    Vec mean(hA);
    solve_lower(chol, na, mean.data());
    solve_lower_transpose(chol, na, mean.data());
    const double log_det = log_det_cholesky(chol, na);
    Vec proposal(na), old(na), db(na);
    draw_precision(mean, chol, rng, proposal.data());
    if (!all_finite(proposal.data(), na)) return;
    for (std::size_t i = 0; i < na; ++i) {
      old[i] = _b[A[i]];
      db[i] = proposal[i] - old[i];
    }
    Shift sh;
    for (std::size_t i = 0; i < na; ++i) {
      sh.coords.push_back(A[i]);
      sh.dcoef.push_back(db[i]);
    }
    for (unsigned lat = 0; lat < S.latents.size(); ++lat) {
      long double v = 0;
      for (std::size_t i = 0; i < na; ++i) v += static_cast<long double>(_u[lat * p + A[i]]) * db[i];
      if (v == 0) continue;
      sh.lats.push_back(lat);
      sh.dlat.push_back(static_cast<double>(-v / _lambda[lat]));
    }
    prepare_shift(sh, A);
    const long double children = children_delta(sh, 1.0);
    const long double proposals = log_normal_precision(old.data(), mean, chol, log_det) -
      log_normal_precision(proposal.data(), mean, chol, log_det);
    const long double target = children + prior_delta(sh, 1.0);
    const double full = static_cast<double>(target + proposals);
    if (dry) {
      probe_record(joint, 0.0, gaussian_part(sh) + prior_delta(sh, 1.0) + proposals, children - gaussian_part(sh), 0.0);
      return;
    }
    if (!_selection) record_max(counters().max_stage2_gaussian, std::fabs(full));
    if (!std::isfinite(full) || !metropolis(rng, full)) return;
    count(counters().accepted[joint]);
    apply_shift(sh, 1.0);
    _repr_delta += target;
  }

  // ---- E: the elliptical slice refresh of the active coefficients (A3)

  // The Gaussian working model of the active coefficients with the latents
  // integrated: the coordinates A, the law N(mu, P^-1) (Cholesky factor of P).
  bool active_law(std::vector<unsigned> &A, Vec &mu, Vec &chol, double &log_det)
  {
    Shared const &S = s();
    const unsigned p = S.p;
    A.clear();
    for (unsigned j = 0; j < S.groups.size(); ++j) {
      if (!active(j)) continue;
      for (unsigned i = 0; i < S.groups[j].d; ++i) A.push_back(S.groups[j].coord0 + i);
    }
    const std::size_t na = A.size();
    if (na == 0) return false;
    chol.assign(na * na, 0.0);
    mu.assign(na, 0.0);
    for (std::size_t i = 0; i < na; ++i) {
      mu[i] = _h[A[i]];
      for (std::size_t k = 0; k < na; ++k) chol[i * na + k] = _P[A[i] * p + A[k]];
    }
    if (!cholesky(chol, na)) return false;
    solve_lower(chol, na, mu.data());
    solve_lower_transpose(chol, na, mu.data());
    log_det = log_det_cholesky(chol, na);
    return true;
  }

  // The shift of the active coefficients by db with every latent sheared,
  // dZ_s = -u_s' db / lambda_s (the residuals Z_s - mu_s(beta) held fixed;
  // mutant 5 drops the shear).
  void active_shift(std::vector<unsigned> const &A, Vec const &db, Shift &sh)
  {
    Shared const &S = s();
    const unsigned p = S.p;
    sh.clear();
    for (std::size_t i = 0; i < A.size(); ++i) {
      sh.coords.push_back(A[i]);
      sh.dcoef.push_back(db[i]);
    }
    if (S.mutant != 5) {
      for (unsigned lat = 0; lat < S.latents.size(); ++lat) {
        long double v = 0;
        for (std::size_t i = 0; i < A.size(); ++i) v += static_cast<long double>(_u[lat * p + A[i]]) * db[i];
        if (v == 0) continue;
        sh.lats.push_back(lat);
        sh.dlat.push_back(static_cast<double>(-v / _lambda[lat]));
      }
    }
    prepare_shift(sh, A);
  }

  // One elliptical slice step (Murray, Adams and MacKay 2010) on the active
  // coefficients with the latents sheared. Along the shear the inner target is
  // N(beta_a; mu_a, P_a^-1) T, log T the sum over the children of the inner
  // target's term plus 1/2 r'Wr (zero for Gaussian children); the ellipse
  // through the current point and a draw nu ~ N(0, P_a^-1) around mu_a leaves
  // N(mu_a, P_a^-1) invariant and the slice on T makes the step reversible for
  // the inner target. Selection labels only (J is the exact draw otherwise).
  void ellipse(RNG *rng, bool dry = false)
  {
    Shared const &S = s();
    std::vector<unsigned> A;
    Vec mu, chol;
    double log_det = 0;
    if (!active_law(A, mu, chol, log_det)) { skipped_factorization(); return; }
    const std::size_t na = A.size();
    Vec old(na), nu(na), center(mu), db(na);
    for (std::size_t i = 0; i < na; ++i) old[i] = _b[A[i]];
    // mutant 4: the ellipse around the slab prior means instead of mu_a
    if (S.mutant == 4) {
      for (std::size_t i = 0; i < na; ++i) {
        GroupSpec const &g = S.groups[S.coord_group[A[i]]];
        center[i] = g.m0[A[i] - g.coord0];
      }
    }
    Vec zero(na, 0.0);
    draw_precision(zero, chol, rng, nu.data());
    if (!all_finite(nu.data(), na)) return;
    auto at = [&](double theta, Shift &sh) {
      const double c = std::cos(theta), sn = std::sin(theta);
      for (std::size_t i = 0; i < na; ++i) db[i] = center[i] + (old[i] - center[i]) * c + nu[i] * sn - old[i];
      active_shift(A, db, sh);
    };
    const double two_pi = 6.283185307179586477;
    if (dry) {
      Shift sh;
      at(two_pi * rng->uniform(), sh);
      ellipse_probe(A, mu, chol, log_det, sh);
      return;
    }
    count(counters().proposed[ellipse_move]);
    const double log_u = std::log(rng->uniform());
    double theta = two_pi * rng->uniform();
    double lower = theta - two_pi, upper = theta;
    for (;;) {
      Shift sh;
      at(theta, sh);
      const long double children = children_delta(sh, 1.0);
      const long double log_t = children - gaussian_part(sh);
      if (children > negative_infinity && std::isfinite(static_cast<double>(log_t)) && log_t > log_u) {
        count(counters().accepted[ellipse_move]);
        apply_shift(sh, 1.0);
        _repr_delta += children + prior_delta(sh, 1.0);
        return;
      }
      count(counters().ellipse_shrinks);
      if (theta < 0) lower = theta;
      else upper = theta;
      if (!(upper - lower >= 1e-10)) return;
      theta = lower + (upper - lower) * rng->uniform();
    }
  }

  // Test hook (probe) of E at the entry state: log T of a proposal as the block
  // computes it, and the changes of the inner target and of the Gaussian
  // working model evaluated separately (the slab and latent priors, every
  // child's term, N(beta_a; mu_a, P_a^-1) and each N(Z_s; mu_s(beta_a),
  // 1/lambda_s) from scratch); log T must equal their difference.
  void ellipse_probe(std::vector<unsigned> const &A, Vec const &mu, Vec const &chol, double log_det, Shift &sh)
  {
    Shared const &S = s();
    const long double children = children_delta(sh, 1.0);
    const long double log_t = children - gaussian_part(sh);
    const std::size_t na = A.size();
    Vec b_old(na), b_new(na);
    Vec dc(S.p, 0.0);
    for (std::size_t i = 0; i < sh.coords.size(); ++i) dc[sh.coords[i]] += sh.dcoef[i];
    for (std::size_t i = 0; i < na; ++i) {
      b_old[i] = _b[A[i]];
      b_new[i] = _b[A[i]] + dc[A[i]];
    }
    Vec dz(S.latents.size(), 0.0);
    for (std::size_t i = 0; i < sh.lats.size(); ++i) dz[sh.lats[i]] += sh.dlat[i];
    // the inner target: priors, then every changed child
    long double inner = 0;
    std::vector<char> done(S.groups.size(), 0);
    for (unsigned k : A) {
      const unsigned j = S.coord_group[k];
      if (done[j]) continue;
      done[j] = 1;
      GroupSpec const &g = S.groups[j];
      Vec next(g.d);
      for (unsigned i = 0; i < g.d; ++i) next[i] = _beta[g.coord0 + i] + dc[g.coord0 + i];
      inner += log_slab(g, next.data()) - log_slab(g, &_beta[g.coord0]);
    }
    for (std::size_t lat = 0; lat < dz.size(); ++lat) {
      inner += -0.5L * (static_cast<long double>(_z[lat] + dz[lat]) * (_z[lat] + dz[lat]) -
                        static_cast<long double>(_z[lat]) * _z[lat]);
    }
    // residuals of every child at the old and the new state
    Vec r_new(_r);
    for (std::size_t q = 0; q < sh.children.size(); ++q) {
      ChildSpec const &child = S.children[sh.children[q]];
      for (unsigned i = 0; i < child.n; ++i) r_new[child.row0 + i] -= sh.delta[sh.offset[q] + i];
    }
    auto quadratic = [&](ChildSpec const &child, Vec const &r) {
      double const *W = &_W[child.woff];
      long double v = 0;
      for (unsigned i = 0; i < child.n; ++i) {
        for (unsigned j = 0; j < child.n; ++j) v += static_cast<long double>(r[child.row0 + i]) * W[i * child.n + j] * r[child.row0 + j];
      }
      return v;
    };
    for (std::size_t q = 0; q < sh.children.size(); ++q) {
      const unsigned c = sh.children[q];
      ChildSpec const &child = S.children[c];
      if (gaussian_child(c)) {
        inner += -0.5L * (quadratic(child, r_new) - quadratic(child, _r));
      } else {
        ensure_value(c);
        inner += static_cast<long double>(sh.value[q]) - _value[c];
      }
    }
    // the Gaussian working model: N(beta_a; mu_a, P_a^-1) and, per latent, the
    // conditional N(Z_s; mu_s(beta_a), 1/lambda_s), mu_s = d_s' W r_s / lambda_s
    // with r_s the residuals without the latent
    long double gaussian = log_normal_precision(b_new.data(), mu, chol, log_det) -
      log_normal_precision(b_old.data(), mu, chol, log_det);
    for (unsigned lat = 0; lat < S.latents.size(); ++lat) {
      long double m_old = 0, m_new = 0;
      for (unsigned c : S.latent_children[lat]) {
        ChildSpec const &child = S.children[c];
        double const *W = &_W[child.woff];
        for (unsigned i = 0; i < child.n; ++i) {
          for (unsigned j = 0; j < child.n; ++j) {
            const unsigned rj = child.row0 + j;
            m_old += static_cast<long double>(_d[child.row0 + i]) * W[i * child.n + j] * (_r[rj] + _d[rj] * _z[lat]);
            m_new += static_cast<long double>(_d[child.row0 + i]) * W[i * child.n + j] *
              (r_new[rj] + _d[rj] * (_z[lat] + dz[lat]));
          }
        }
      }
      m_old /= _lambda[lat];
      m_new /= _lambda[lat];
      const long double e_old = _z[lat] - m_old, e_new = _z[lat] + dz[lat] - m_new;
      gaussian += -0.5L * _lambda[lat] * (e_new * e_new - e_old * e_old);
    }
    probe_record(ellipse_move, static_cast<double>(log_t), inner - gaussian, inner, static_cast<double>(gaussian));
  }

  // S_jk: exchange the gates of two correlated groups, the one switched on
  // drawn from its law along the line given the other off (Gaussian labels).
  void gaussian_swap(unsigned first, unsigned second, RNG *rng, bool dry = false)
  {
    Shared const &S = s();
    if (active(first) == active(second)) return;
    const unsigned j = active(first) ? first : second;   // on -> off
    const unsigned k = active(first) ? second : first;   // off -> on
    GroupSpec const &gj = S.groups[j];
    GroupSpec const &gk = S.groups[k];
    Line L;
    if (!build_line({j, k}, L)) { skipped_factorization(); return; }
    if (!dry) {
      count(counters().proposed[swap]);
      count(counters().stage0[swap]);
    }
    const std::size_t dj = gj.d, dk = gk.d, nm = dj + dk;
    // marginal laws of each group (from the joint of the pair) and the
    // conditionals of one given the other at zero
    Vec mj, Pj, mk, Pk, cj, Cj, ck, Ck;
    double lj, lk, lcj, lck;
    if (!swap_laws(L, dj, dk, mj, Pj, lj, mk, Pk, lk, cj, Cj, lcj, ck, Ck, lck)) { skipped_factorization(); return; }
    Vec zero_j(dj, 0.0), zero_k(dk, 0.0);
    const double stage1 = (gj.log_off + gk.log_on - gj.log_on - gk.log_off) +
      log_slab(gk, zero_k.data()) + log_normal_precision(zero_j.data(), mj, Pj, lj) -
      log_slab(gj, zero_j.data()) - log_normal_precision(zero_k.data(), mk, Pk, lk);
    if (!dry) {
      if (!first_stage(rng, stage1)) return;
      count(counters().stage1[swap]);
    }
    Vec beta_k(dk), dormant_j(dj);
    draw_precision(ck, Ck, rng, beta_k.data());
    draw_prior(gj, rng, dormant_j.data());
    if (!all_finite(beta_k.data(), dk) || !all_finite(dormant_j.data(), dj)) return;
    double const *old_j = &_beta[gj.coord0];
    double const *old_k = &_beta[gk.coord0];
    Vec db(nm);
    for (std::size_t i = 0; i < dj; ++i) db[i] = -old_j[i];
    for (std::size_t i = 0; i < dk; ++i) db[dj + i] = beta_k[i];
    Shift sh;
    line_shift(L, db.data(), sh);
    prepare_shift(sh, L.partners);
    const long double children = children_delta(sh, 1.0);
    const long double rest = prior_delta(sh, 1.0) + (gj.log_off + gk.log_on - gj.log_on - gk.log_off) +
      log_slab(gj, dormant_j.data()) - log_slab(gj, old_j) + log_slab(gk, beta_k.data()) - log_slab(gk, old_k);
    const long double target = children + rest;
    const long double proposals = log_normal_precision(old_j, cj, Cj, lcj) + log_slab(gk, old_k) -
      log_normal_precision(beta_k.data(), ck, Ck, lck) - log_slab(gj, dormant_j.data());
    const double full = static_cast<double>(target + proposals);
    if (dry) {
      probe_record(swap, stage1, gaussian_part(sh) + rest + proposals, children - gaussian_part(sh), stage1);
      return;
    }
    const double stage2 = full - stage1;
    if (!_selection) record_max(counters().max_stage2_gaussian, std::fabs(stage2));
    if (!std::isfinite(full) || !metropolis(rng, stage2)) return;
    count(counters().accepted[swap]);
    apply_shift(sh, 1.0);
    _g[j] = 0;
    _g[k] = 1;
    std::copy(dormant_j.begin(), dormant_j.end(), _beta.begin() + gj.coord0);
    std::copy(beta_k.begin(), beta_k.end(), _beta.begin() + gk.coord0);
    for (std::size_t i = 0; i < dj; ++i) _b[gj.coord0 + i] = 0.0;
    for (std::size_t i = 0; i < dk; ++i) _b[gk.coord0 + i] = beta_k[i];
    _repr_delta += target;
  }

  // Swap pairs: gated groups whose coordinates correlate above 0.5 in the
  // Gaussian posterior with every gate on (non-members only).
  std::vector<std::pair<unsigned, unsigned>> swap_pairs()
  {
    Shared const &S = s();
    const unsigned p = S.p;
    std::vector<std::pair<unsigned, unsigned>> pairs;
    Vec chol(_P);
    if (!cholesky(chol, p)) { skipped_factorization(); return pairs; }
    const Vec Sigma = inverse_from_cholesky(chol, p);
    for (unsigned j = 0; j < S.groups.size(); ++j) {
      if (S.groups[j].gate == nullptr) continue;
      for (unsigned k = j + 1; k < S.groups.size(); ++k) {
        if (S.groups[k].gate == nullptr) continue;
        double best = 0;
        for (unsigned a = 0; a < S.groups[j].d; ++a) {
          for (unsigned b = 0; b < S.groups[k].d; ++b) {
            const unsigned x = S.groups[j].coord0 + a, y = S.groups[k].coord0 + b;
            best = std::max(best, std::fabs(Sigma[x * p + y] / std::sqrt(Sigma[x * p + x] * Sigma[y * p + y])));
          }
        }
        if (best > swap_correlation) pairs.emplace_back(j, k);
      }
    }
    return pairs;
  }

#include "LocationBlockBiasMoves.h"

  // ---- the update

  enum Kind { kind_flip, kind_refresh, kind_swap, kind_joint, kind_ellipse, kind_latent, kind_bias };
  struct Step { Kind kind; unsigned a, b; };

  // Test hook (trace bit 2): the sequence of one sweep, the slots in the order
  // they run (kind, a, b) and within a composite slice refresh its
  // coordinates (slot_coordinate, group, coordinate); written to the sequence
  // store before the exact correction, after a header (-1, chain, reversed,
  // entry label is a selection label, entries).
  static constexpr unsigned slot_coordinate = 100;
  Vec _sequence;

  void sequence_step(unsigned kind, unsigned a, unsigned b)
  {
    if ((s().trace & 2) == 0) return;
    _sequence.insert(_sequence.end(), {static_cast<double>(kind), static_cast<double>(a), static_cast<double>(b)});
  }

  void sequence_flush()
  {
    if ((s().trace & 2) == 0) return;
    ProbeStore &store = sequence_store();
    const std::lock_guard<std::mutex> lock(store.mutex);
    store.records.insert(store.records.end(), {-1.0, static_cast<double>(_chain), _backward ? 1.0 : 0.0,
                                               _selection_entry ? 1.0 : 0.0, static_cast<double>(_sequence.size() / 3)});
    store.records.insert(store.records.end(), _sequence.begin(), _sequence.end());
    _sequence.clear();
  }

  // The refresh kernels of the active coefficients in selection labels (test
  // switch, 0 the default): 1 the per-group refreshes R_j, 2 one elliptical
  // slice step E, 3 two E steps, 4 E and the R_j, 5 the per-group refreshes
  // with every group's slice update forced from the first adaptOff() on.
  int refresh_variant() const { return s().refresh == 0 ? default_refresh : s().refresh; }

  void sweep(RNG *rng)
  {
    Shared const &S = s();
    if (S.timing) _part_start = std::chrono::steady_clock::now();
    _category = cat_prepare;
    prepare();
    enter(cat_move);
    count(counters().sweeps);
    std::vector<Step> steps;
    // the bias-label components first (D, L0, L1, S, L2; mirrored in the
    // reversed sweep), their slots fixed for the fit
    if (bias()) {
      for (unsigned k = 0; k < bias_component_count; ++k) steps.push_back({kind_bias, k, 0});
      bias_sweep_start(rng);
    }
    for (unsigned j = 0; j < S.groups.size(); ++j) {
      if (S.groups[j].gate != nullptr) steps.push_back({kind_flip, j, 0});
      steps.push_back({kind_refresh, j, 0});
    }
    const std::vector<std::pair<unsigned, unsigned>> pairs = swap_pairs();
    for (auto const &pair : pairs) steps.push_back({kind_swap, pair.first, pair.second});
    // The slot list is a function of the fit and of non-members only (sweep
    // condition 1): components whose eligibility depends on the label occupy
    // fixed slots and decide when they run whether they act. The joint
    // refresh acts in Gaussian labels only (the untilted joint law is rejected
    // under selection), E in selection labels only, after the flips and swaps
    // and before the latents.
    // (mutant 11: the joint slot from the entry label, not reversible once the
    // label can change inside the sweep)
    if (!S.groups.empty() && (S.mutant != 11 || !_selection)) steps.push_back({kind_joint, 0, 0});
    const int variant = refresh_variant();
    const int ellipses = variant == 3 ? 2 : (variant == 2 || variant == 4) ? 1 : 0;
    for (int e = 0; e < ellipses; ++e) steps.push_back({kind_ellipse, 0, 0});
    for (unsigned lat = 0; lat < S.latents.size(); ++lat) steps.push_back({kind_latent, lat, 0});
    if (S.probe) {
      // test hook: dry runs of every collapsed Gaussian move at the entry state
      for (unsigned j = 0; j < S.groups.size(); ++j) {
        if (S.groups[j].gate != nullptr) gaussian_flip(j, rng, true);
        if (active(j)) gaussian_refresh(j, rng, true);
      }
      for (auto const &pair : pairs) gaussian_swap(pair.first, pair.second, rng, true);
      if (!S.groups.empty()) gaussian_joint(rng, true);
      for (unsigned lat = 0; lat < S.latents.size(); ++lat) gaussian_latent(lat, rng, true);
      ellipse(rng, true);
    }
    // one coin per update: the whole list forward or reversed
    _backward = rng->uniform() < 0.5;
    _sequence.clear();
    for (std::size_t i = 0; i < steps.size(); ++i) {
      Step const &step = steps[_backward ? steps.size() - 1 - i : i];
      // mutants 27, 28 (test hooks): the dormant draws D and the swap L2 alone
      if ((S.mutant == 27 || S.mutant == 28) &&
          !(step.kind == kind_bias && (step.a == static_cast<unsigned>(bias_swap) ||
                                       step.a == static_cast<unsigned>(bias_dormant)))) continue;
      sequence_step(step.kind, step.a, step.b);
      // The line of a group serves its F_j and R_j, adjacent in either order;
      // any other component changes it.
      if ((step.kind != kind_flip && step.kind != kind_refresh) || _line.group != step.a) _line.valid = false;
      enter(step.kind == kind_latent ? cat_latent : step.kind == kind_ellipse ? cat_ellipse :
            step.kind == kind_bias ? bias_category(step.a) : cat_move);
      switch (step.kind) {
      case kind_bias:
        bias_component(step.a, rng);
        break;
      case kind_flip:
        if (_selection) tilted_flip(step.a, rng);
        else gaussian_flip(step.a, rng);
        break;
      case kind_refresh:
        if (!active(step.a)) dormant_refresh(step.a, rng);
        else if (!_selection) gaussian_refresh(step.a, rng);
        else if (variant == 2 || variant == 3) break;  // E refreshes the active coefficients
        else if (_slice_refresh[step.a]) slice_refresh(step.a, rng);
        else tilted_refresh(step.a, rng);
        break;
      case kind_ellipse:
        if (_selection) ellipse(rng);
        break;
      case kind_swap:
        if (_selection) tilted_swap(step.a, step.b, rng);
        else gaussian_swap(step.a, step.b, rng);
        break;
      case kind_joint:
        if (!_selection) gaussian_joint(rng);
        break;
      case kind_latent:
        if (_selection) slice_latent(step.a, rng);
        else gaussian_latent(step.a, rng);
        break;
      }
    }
    sequence_flush();
    if (selection_endpoint()) count(counters().sweeps_selection);
    enter(cat_correction);
    correct(rng);
    enter(cat_move);
    flush_evidence();
    flush_bias_evidence();
  }

  // Test hook (mutant 15 with the probe): the graph's endpoints, read
  // independently of the block's buffers: the members' values, the values of
  // the view's deterministic descendants and every parent value of every
  // stochastic child of the view (the likelihood children, the gate markers
  // and the nested node: means, weights, kernel modes, rules and the nested
  // node's parameters), and their densities (JAGS's own logDensity,
  // PDF_LIKELIHOOD) with the members' prior densities.
  struct Endpoints {
    Vec members, parents, densities;
  };
  Endpoints _restoration_entry;

  void entry_endpoints(Endpoints &out) const
  {
    out.members.assign(s().length, 0.0);
    _view->getValue(out.members, _chain);
    out.parents.clear();
    out.densities.clear();
    for (DeterministicNode const *node : _view->deterministicChildren()) {
      double const *value = node->value(_chain);
      out.parents.insert(out.parents.end(), value, value + node->length());
    }
    for (StochasticNode const *child : _view->stochasticChildren()) {
      for (Node const *parent : child->parents()) {
        double const *value = parent->value(_chain);
        out.parents.insert(out.parents.end(), value, value + parent->length());
      }
      out.densities.push_back(child->logDensity(_chain, PDF_LIKELIHOOD));
    }
    for (StochasticNode const *member : _view->nodes()) out.densities.push_back(member->logDensity(_chain, PDF_PRIOR));
  }

  // After the forced rejection and the restoring write: the endpoints read
  // again against those of the entry (probe record type 11: 1, the largest
  // differences of the members, the parent values and the densities, and
  // whether the counts agree). All must be zero.
  void restoration_probe()
  {
    Endpoints now;
    entry_endpoints(now);
    Endpoints const &entry = _restoration_entry;
    probe_raw(11, 1.0, endpoint_difference(now.members, entry.members),
              endpoint_difference(now.parents, entry.parents), endpoint_difference(now.densities, entry.densities),
              now.parents.size() == entry.parents.size() && now.densities.size() == entry.densities.size() ? 0.0 : 1.0);
  }

  // The graph's own log density of a child (PDF_LIKELIHOOD) at the graph's
  // current mean: a step-switch row through the distribution's evaluations
  // with its parameters read once per update (the same function and inputs as
  // its logDensity, B1), every other child through JAGS.
  double graph_density(unsigned c)
  {
    ChildSpec const &child = s().children[c];
    if (child.type == ChildType::step_switch && _step[c].valid) {
      StepChild const &st = _step[c];
      const double mean = child.node->parents()[0]->value(_chain)[0];
      if (!std::isfinite(mean)) return negative_infinity;
      return step_row_log_density(step_law(c), st.power, st.bin, st.mode, st.telescope != 0,
                                  child.node->value(_chain)[0], mean);
    }
    return child.node->logDensity(_chain, PDF_LIKELIHOOD);
  }

  // The exact correction: accept the sweep with [pi(new) / pi~(new)] /
  // [pi(old) / pi~(old)] over the changed members and the children whose
  // means changed, the exact terms from the graph's own densities.
  // The member values of the working state (the update's entry with the
  // gates, slabs, latents and the bias label and sources of the sweep).
  void working_members(Vec &out) const
  {
    Shared const &S = s();
    out = _entry;
    for (unsigned j = 0; j < S.groups.size(); ++j) {
      GroupSpec const &g = S.groups[j];
      if (g.gate != nullptr) out[g.gate_pos] = _g[j] ? g.label : 3.0 - g.label;
      for (unsigned i = 0; i < g.d; ++i) out[g.slab_pos + i] = _beta[g.coord0 + i];
    }
    for (unsigned lat = 0; lat < S.latents.size(); ++lat) out[S.latent_pos[lat]] = _z[lat];
    if (bias()) {
      out[S.bias->pos] = _label;
      for (int b = 0; b < 2; ++b) if (S.bias->slope[b].source != nullptr) out[S.bias->slope[b].pos] = _theta[b];
    }
  }

  void correct(RNG *rng)
  {
    Shared const &S = s();
    // the final member values
    working_members(_final);
    const bool label_changed = bias() && _label != _label_entry;
    if (_final == _entry) {
      count(counters().corrections_accepted);
      evidence_endpoint(true);
      bias_after_correction(true);
      return;
    }
    std::vector<StochasticNode *> const &members = _view->nodes();
    std::vector<unsigned> changed_members;
    {
      std::size_t position = 0;
      for (unsigned m = 0; m < members.size(); ++m) {
        const unsigned long len = members[m]->length();
        bool changed = false;
        for (unsigned long i = 0; i < len; ++i) changed = changed || _final[position + i] != _entry[position + i];
        if (changed) changed_members.push_back(m);
        position += len;
      }
    }
    std::vector<unsigned> kids;
    // every child counts as changed when the label changed (its weights,
    // kernel mode or mean may differ), and so does the nested node
    for (unsigned c = 0; c < S.children.size(); ++c) if (_dirty[c] || label_changed) kids.push_back(c);
    // the graph holds the entry state once the bias label's nodes are written
    // back (the block wrote them to evaluate other labels, never the means)
    if (bias()) set_graph_label(_label_entry);
    // mutant 15 with the probe (test hook): the entry's graph endpoints, read
    // independently of the block's buffers (restoration_probe)
    if (S.mutant == 15 && S.probe && label_changed) entry_endpoints(_restoration_entry);
    // The exact log ratio pairs each node's density at the final state with
    // its own entry value and sums the differences: totals of the two
    // endpoints would let one large term shared by both erase the others'
    // differences. An entry term that is not finite is an error naming its
    // node; a final term of -Inf (zero density) rejects the sweep.
    // The entry's terms (the graph still holds it) in the order changed
    // members, children, nested node, and the endpoint differences of the
    // selection children, sum (log A - log A~) up to sign: inner-target value
    // minus exact value.
    std::vector<double> &before = _correction_terms;
    before.clear();
    long double old_endpoint = 0, new_endpoint = 0;
    auto entry_term = [&](Node const *node, double v) {
      if (!std::isfinite(v)) {
        throwNodeError(node,
          "Location block sampler: the log density of this node is not finite at the update's entry state");
      }
      before.push_back(v);
    };
    for (unsigned m : changed_members) entry_term(members[m], members[m]->logDensity(_chain, PDF_PRIOR));
    // the precision guard over the entry's likelihood terms (CorrectedSampler.h)
    LikelihoodMagnitude magnitude;
    for (unsigned c : kids) {
      const double v = graph_density(c);
      count(counters().correction_evaluations);
      count(counters().evaluations[cat_correction]);
      entry_term(S.children[c].node, v);
      magnitude.add(S.children[c].node, v);
      if (_entry_valued[c]) old_endpoint += static_cast<long double>(_value_entry[c]) - v;
    }
    check_likelihood_magnitude(magnitude, _adapting, corrected_kind_location, "Location block sampler");
    if (label_changed) entry_term(bspec().nested, nested_term());
    _view->setValue(_final, _chain);
    // mutant 14: one error after this write, in the first sweep after adaptation
    // that changes the state (update() restores the entry)
    if (S.mutant == 14 && !_adapting && !_injected) {
      _injected = true;
      throw std::runtime_error("Location block sampler: an error injected after the final graph write (mutant 14)");
    }
    long double exact_change = 0;
    bool zero_density = false;
    std::size_t term = 0;
    auto final_term = [&](Node const *node, double v) {
      if (std::isnan(v) || v == std::numeric_limits<double>::infinity()) {
        throwNodeError(node,
          "Location block sampler: the log density of this node is not defined at the proposed state");
      }
      if (v == negative_infinity) zero_density = true;
      else exact_change += static_cast<long double>(v) - before[term];
      ++term;
    };
    if (label_changed) {
      // the graph's label is the final one: the rows' parameters follow it
      _graph_label = _label;
      for (unsigned c = 0; c < S.children.size(); ++c) if (_step_children[c]) refresh_step(c);
    }
    // A likelihood mean that the graph changed although no component of the
    // sweep moved it, or that differs from the mean the block tracked for the
    // final state where the sweep moved it: the verified design does not hold
    // at this state. An error naming the node (update() restores the entry),
    // never a silent rejection. A tracked mean is compared with a tolerance of
    // 1e-6 of the row's conditional SD (the precision matrix's diagonal) plus
    // 1e-12 of the means' size (the rounding of the block's accumulated
    // changes and of the graph's expression) and, for a moved row, the rounding
    // of the probed slopes times the coordinates' changes over the sweep, so
    // that a design entry the probe missed is not hidden by another coefficient
    // that moves the same child. The comparison runs while adapting too.
    // An unmoved mean is compared with its entry mean under the same
    // tolerance (a diagnostic of the design: a dependence the probe missed
    // changes the mean by far more than rounding). Its density is never left
    // out: the graph's recomputation of a mean that nothing moved need not be
    // bitwise identical on every platform (a Linux build with OpenBLAS changed
    // one by rounding), and a coefficient whose response rounded to zero in
    // the factory's probe (a mean offset + 2^-43 b at an offset of 32768)
    // still moves it, so every child outside `kids` with a row whose graph
    // mean after the final write is not bitwise its entry mean enters the
    // paired correction below with the graph's own densities at both
    // endpoints (mutant 21 leaves these children out). After adaptation the
    // comparisons of unmoved means that are not bitwise equal are counted
    // (unmoved_mean_inexact) with their largest discrepancy in the row's SDs
    // (unmoved_mean_max), and the children included this way
    // (unmoved_included) with their largest density change (unmoved_density_max).
    std::vector<unsigned> &unmoved = _unmoved_changed;
    unmoved.clear();
    for (unsigned c = 0; c < S.children.size(); ++c) {
      ChildSpec const &child = S.children[c];
      double const *mean = child.node->parents()[0]->value(_chain);
      const bool moved = _dirty[c] != 0;
      bool changed = false;
      for (unsigned i = 0; i < child.n; ++i) {
        const unsigned row = child.row0 + i;
        if (!moved && mean[i] != _mean_entry[row]) changed = true;
        const double sd = 1.0 / std::sqrt(_W[child.woff + i * child.n + i]);
        double expected = moved ? _mean[row] : _mean_entry[row];
        // mutants 17, 18 (test hook): the entry mean an unmoved mean is
        // compared with, one ulp or 1e-3 of the row's SD above it; 22, 23:
        // 0.1 or 10 times the row's tolerance above it (pinned below)
        const bool pinned = !moved && (S.mutant == 22 || S.mutant == 23);
        if (!moved) {
          if (S.mutant == 17) expected = std::nextafter(expected, std::numeric_limits<double>::infinity());
          if (S.mutant == 18) expected += 1e-3 * sd;
          if (!pinned && mean[i] == expected) continue;
        }
        // the size of the means and of the design's terms at the final state
        double size = std::max(std::fabs(_mean_entry[row]), std::max(std::fabs(mean[i]), std::fabs(expected)));
        for (unsigned k = 0; k < S.p; ++k) size += std::fabs(S.X[static_cast<std::size_t>(row) * S.p + k] * _b[k]);
        if (child.latent >= 0) size += std::fabs(_d[row] * _z[static_cast<unsigned>(child.latent)]);
        double tolerance = tracked_mean_tolerance * sd + tracked_mean_rounding * size;
        if (moved) {
          // A tracked mean is the entry mean plus the probed slopes times the
          // coordinates' changes, so it carries each slope's rounding times
          // its coordinate's change over the sweep: a slope is a difference of
          // two probed means at a step of one (the design's and the bias
          // slopes' at the factory's baseline, the latent's at the update's
          // entry), so its error is at most 4 eps times the size of those
          // means: for the design's and the bias slopes' columns the largest
          // |mu| of the probes (X_size, column_size, column_check's bound,
          // which does not cover a cancellation inside the graph's expression
          // at the factory's baseline), for the latent the sum of the entry's
          // absolute terms (_d_size), which can cancel: an intercept of 6,300
          // against a cluster effect of -6,300, then a latent change of 5e4,
          // missed the tolerance above by twice it while adapting.
          const double eps = std::numeric_limits<double>::epsilon();
          long double slopes = 0;
          for (unsigned k = 0; k < S.p; ++k) {
            const std::size_t at = static_cast<std::size_t>(row) * S.p + k;
            slopes += static_cast<long double>(S.X_size[at]) * std::fabs(_b[k] - _b_entry[k]);
          }
          if (child.latent >= 0) {
            const unsigned lat = static_cast<unsigned>(child.latent);
            slopes += static_cast<long double>(_d_size[row]) * std::fabs(_z[lat] - _entry[S.latent_pos[lat]]);
          }
          if (bias()) {
            for (int b = 0; b < 2; ++b) {
              SlopeSpec const &spec = S.bias->slope[b];
              if (spec.source == nullptr) continue;
              const double entry_slope = spec.label == _label_entry ? _entry[spec.pos] : 0.0;
              slopes += static_cast<long double>(spec.column_size[row]) * std::fabs(physical_at(_label, b) - entry_slope);
            }
          }
          tolerance += 4.0 * eps * static_cast<double>(slopes);
        }
        if (pinned) expected += (S.mutant == 22 ? 0.1 : 10.0) * tolerance;
        const double discrepancy = std::fabs(mean[i] - expected);
        if (!moved && !_adapting) {
          count(counters().unmoved_mean_inexact);
          record_max(counters().unmoved_mean_max, discrepancy / sd);
        }
        if (!(discrepancy <= tolerance)) {
          throwNodeError(child.node, moved ? tracked_mean_message : untracked_mean_message);
        }
      }
      if (changed && !moved && !label_changed && S.mutant != 21) unmoved.push_back(c);
    }
    for (unsigned m : changed_members) final_term(members[m], members[m]->logDensity(_chain, PDF_PRIOR));
    for (unsigned c : kids) {
      const double v = graph_density(c);
      count(counters().correction_evaluations);
      count(counters().evaluations[cat_correction]);
      final_term(S.children[c].node, v);
      if (!gaussian_child(c) && _value_ok[c] && v > negative_infinity) {
        new_endpoint += static_cast<long double>(_value[c]) - v;
      }
    }
    if (!unmoved.empty()) {
      // Their entry densities are no longer in the graph: the entry state is
      // written back once, those children evaluated (as entry terms of a
      // child in `kids`: finite or an error, the precision guard, the
      // endpoint record of a selection child; a step row's parameters are
      // the update's, the label being unchanged), and the final state
      // written again, where each is paired with its own entry term.
      std::vector<double> &entry_values = _unmoved_terms;
      entry_values.clear();
      _view->setValue(_entry, _chain);
      for (unsigned c : unmoved) {
        const double v = graph_density(c);
        count(counters().correction_evaluations);
        count(counters().evaluations[cat_correction]);
        if (!std::isfinite(v)) {
          throwNodeError(S.children[c].node,
            "Location block sampler: the log density of this node is not finite at the update's entry state");
        }
        entry_values.push_back(v);
        magnitude.add(S.children[c].node, v);
        if (_entry_valued[c]) old_endpoint += static_cast<long double>(_value_entry[c]) - v;
      }
      _view->setValue(_final, _chain);
      check_likelihood_magnitude(magnitude, _adapting, corrected_kind_location, "Location block sampler");
      for (std::size_t q = 0; q < unmoved.size(); ++q) {
        const unsigned c = unmoved[q];
        const double v = graph_density(c);
        count(counters().correction_evaluations);
        count(counters().evaluations[cat_correction]);
        if (std::isnan(v) || v == std::numeric_limits<double>::infinity()) {
          throwNodeError(S.children[c].node,
            "Location block sampler: the log density of this node is not defined at the proposed state");
        }
        if (v == negative_infinity) {
          zero_density = true;
          continue;
        }
        exact_change += static_cast<long double>(v) - entry_values[q];
        if (!gaussian_child(c) && _value_ok[c]) new_endpoint += static_cast<long double>(_value[c]) - v;
        if (!_adapting) {
          count(counters().unmoved_included);
          record_max(counters().unmoved_density_max, std::fabs(v - entry_values[q]));
        }
      }
    }
    if (label_changed) final_term(bspec().nested, nested_term());
    const double log_ratio = zero_density ? negative_infinity : static_cast<double>(exact_change - _repr_delta);
    bool accepted;
    if (zero_density) {
      accepted = false;
    } else if (S.mutant == 15 && label_changed) {
      // mutant 15: every sweep that changed the bias label is rejected
      accepted = false;
    } else if (S.mutant == 3) {
      accepted = true;
    } else {
      accepted = metropolis(rng, log_ratio);
    }
    if (std::isfinite(log_ratio)) {
      record_max(selection_endpoint() ? counters().max_correction : counters().max_correction_gaussian, std::fabs(log_ratio));
    }
    if (selection_endpoint()) {
      record_max(counters().max_endpoint, std::fabs(static_cast<double>(new_endpoint)));
      record_max(counters().max_endpoint_change, std::fabs(static_cast<double>(new_endpoint - old_endpoint)));
    }
    if (accepted) {
      count(counters().corrections_accepted);
      if (selection_endpoint()) count(counters().corrections_selection_accepted);
    } else {
      _view->setValue(_entry, _chain);
      if (S.mutant == 15 && S.probe && label_changed) restoration_probe();
    }
    evidence_endpoint(accepted);
    bias_after_correction(accepted);
    // evidence about the surrogate: a sweep that evaluated a surrogate that
    // was not exact (only selection labels evaluate surrogates)
    const bool evidence = _surrogate_mode && _not_exact && selection_endpoint();
    // the shared bookkeeping of the corrected samplers (process counters,
    // adaptive evidence, reporting window), for sweeps that visited a
    // selection label (a sweep entirely in Gaussian labels needs no surrogate)
    if (_surrogate_mode && selection_endpoint()) _tally.count(accepted, evidence);
    if (!evidence) return;
    count(counters().sweeps_evidence);
    if (accepted) count(counters().corrections_evidence_accepted);
  }

public:
  LocationBlock(GraphView const *view, std::shared_ptr<Shared const> shared, unsigned int chain)
    : _view(view), _s(std::move(shared)), _chain(chain)
  {
    Shared const &S = s();
    const unsigned G = static_cast<unsigned>(S.groups.size());
    _entry.resize(S.length);
    _g.assign(G, 0);
    _beta.assign(S.p, 0.0);
    _b.assign(S.p, 0.0);
    _z.assign(S.latents.size(), 0.0);
    _y.assign(S.rows, 0.0);
    _mean_entry.assign(S.rows, 0.0);
    _mean.assign(S.rows, 0.0);
    _r.assign(S.rows, 0.0);
    _wr.assign(S.rows, 0.0);
    _d.assign(S.rows, 0.0);
    _d_size.assign(S.rows, 0.0);
    _b_entry.assign(S.p, 0.0);
    _W.assign(S.wsize, 0.0);
    _selected.assign(S.children.size(), 0);
    _value_ok.assign(S.children.size(), 0);
    _dirty.assign(S.children.size(), 0);
    _value.assign(S.children.size(), 0.0);
    _value_entry.assign(S.children.size(), 0.0);
    _fit_value.assign(S.children.size(), 0.0);
    _fit_ok.assign(S.children.size(), 0);
    _step.assign(S.children.size(), StepChild());
    _entry_valued.assign(S.children.size(), 0);
    _step_children.assign(S.children.size(), 0);
    for (unsigned c = 0; c < S.children.size(); ++c) _step_children[c] = S.children[c].type == ChildType::step_switch;
    _tilt_position.assign(S.children.size(), -1);
    _gate_evidence.assign(G, GateEvidence());
    _coarse_fit = S.coarse_fit && S.fit_rule != 2;
    _stamp.assign(S.children.size(), 0);
    _a.assign(S.latents.size(), 0.0);
    _lambda.assign(S.latents.size(), 1.0);
    _u.assign(S.latents.size() * S.p, 0.0);
    _P.assign(S.p * S.p, 0.0);
    _h.assign(S.p, 0.0);
    _refresh_proposed.assign(G, 0);
    _refresh_accepted.assign(G, 0);
    _slice_refresh.assign(G, 0);
    for (ChildSpec const &child : S.children) {
      if (child.kind == ChildKind::other) _surrogates.emplace_back();
      else _surrogates.emplace_back(new SelectedChild(child.node, child.kind));
    }
    view->checkFinite(chain);
    count(counters().samplers);
    if (S.corrected && S.any_surrogate) {
      register_corrected_sampler(this, corrected_kind_location);
      _registered = true;
    }
  }

  ~LocationBlock() override
  {
    if (_registered) unregister_corrected_sampler(this);
  }

  std::uint64_t windowProposed() const override { return _tally.window_proposed(); }
  std::uint64_t windowAccepted() const override { return _tally.window_accepted(); }
  void clearWindow() override { _tally.clear_window(); }
  double rateExponent() const override { return 1.0; }

  void update(RNG *rng) override
  {
    // Any exception restores the update-entry values before it propagates.
    std::vector<double> entry(s().length);
    _view->getValue(entry, _chain);
    try {
      sweep(rng);
    } catch (...) {
      _view->setValue(entry, _chain);
      throw;
    }
  }

  bool isAdaptive() const override { return true; }

  // The fallback is decided once, at the first adaptOff() (spec: Fallback):
  // with at least 50 updates that are evidence about the surrogate and a
  // correction acceptance below the threshold, exact mode for the rest of
  // this object's life; with fewer, the sampler stays corrected.
  void adaptOff() override
  {
    if (!_adapting) return;
    _adapting = false;
    // the refresh kernel of each group, from its adaptive independence proposals
    // (refresh variant 5, a test switch: the slice update for every group)
    for (std::size_t j = 0; j < _slice_refresh.size(); ++j) {
      _slice_refresh[j] = s().refresh == 5 || (_refresh_proposed[j] > 0 &&
        static_cast<double>(_refresh_accepted[j]) < 0.4 * static_cast<double>(_refresh_proposed[j]));
      if (!_slice_refresh[j]) continue;
      count(counters().slice_refresh_groups);
      if (s().groups[j].d > 1) count(counters().slice_refresh_multi);
    }
    // the coarse rung of the Newton fits, from the stage-2 acceptance of the
    // tilted flips and swaps that passed stage 1 while adapting (A1)
    if (_coarse_fit && s().fit_rule == 0 && _fit_proposed >= 50 &&
        static_cast<double>(_fit_accepted) < 0.9 * static_cast<double>(_fit_proposed)) {
      _coarse_fit = false;
      count(counters().fit_rule_fallback);
    }
    if (!_registered) return;
    if (_tally.adapt_off(s().min_acceptance, 1.0)) {
      _exact = true;
      count(counters().fallback_exact);
    }
  }

  bool checkAdaptation() const override { return true; }
};

// ---- the factory

struct Declared {
  std::vector<GroupSpec> groups;
  std::vector<StochasticNode *> latents;
  bool has_bias = false;            // a valid bias declaration (LocationBlockBias.h)
  BiasDeclared bias;
  std::string bias_reason;          // why a bias tag was not taken
};

bool integer_value(Node const *node, double lower, double upper, double *out)
{
  if (!fixed_scalar(node)) return false;
  const double v = node->value(0)[0];
  if (!(v >= lower && v <= upper && v == std::floor(v))) return false;
  *out = v;
  return true;
}

bool collect(std::list<StochasticNode *> const &free, Graph const &graph, Declared &out)
{
  std::set<StochasticNode *> available(free.begin(), free.end());
  Node const *setup = nullptr;
  std::map<unsigned, GroupSpec> groups;
  std::map<unsigned, StochasticNode *> latents;
  Node const *bias_tag = nullptr;
  bool bias_tags_valid = true;
  for (StochasticNode *node : free) {
    for (Node const *tag : direct_terminal_tags(node, setup_tag_name, 3, graph)) {
      if (setup != nullptr) return fail("declaration: two setup tags");
      setup = tag;
    }
    // robma_location_group(anchor, slab, label, physical, slot): the anchor is
    // the gate, or for an ungated group (label 0) the scalar slab itself.
    for (Node const *tag : direct_terminal_tags(node, group_tag_name, 5, graph)) {
      auto const &par = tag->parents();
      double label = 0, slot = 0;
      if (!integer_value(par[2], 0, 2, &label) || !integer_value(par[4], 1, 1e6, &slot)) return fail("declaration: label or slot");
      GroupSpec g;
      g.slab = par[1];
      g.physical = par[3];
      if (!resolve_slab(g.slab, g.slab_nodes)) return fail("declaration: slab coordinates");
      if (label == 0) {
        if (g.slab != node || g.physical != node) return fail("declaration: ungated anchor");
      } else {
        StochasticNode *gate = node;
        if (!(label == 1 || label == 2)) return fail("declaration: label");
        const CategoricalWeights weights = categorical_weights(gate, 2);
        if (!weights) return fail("declaration: gate weights");
        g.gate = gate;
        g.label = label;
        const double log_total = std::log(weights.total);
        g.log_on = std::log(weights.values[static_cast<int>(label) - 1]) - log_total;
        g.log_off = std::log(weights.values[2 - static_cast<int>(label)]) - log_total;
      }
      if (groups.count(static_cast<unsigned>(slot))) return fail("declaration: duplicate group slot");
      groups[static_cast<unsigned>(slot)] = g;
    }
    for (Node const *tag : direct_terminal_tags(node, latent_tag_name, 2, graph)) {
      double slot = 0;
      if (!integer_value(tag->parents()[1], 1, 1e7, &slot) || latents.count(static_cast<unsigned>(slot))) return fail("declaration: latent slot");
      latents[static_cast<unsigned>(slot)] = node;
    }
    for (Node const *tag : direct_terminal_tags(node, bias_tag_name, 5, graph)) {
      if (bias_tag != nullptr) bias_tags_valid = false;
      bias_tag = tag;
    }
  }
  if (setup == nullptr) return false;  // not a location model: no decision recorded
  double n_groups = 0, n_latents = 0;
  if (!integer_value(setup->parents()[1], 1, 1e6, &n_groups) ||
      !integer_value(setup->parents()[2], 0, 1e7, &n_latents) ||
      groups.size() != n_groups || latents.size() != n_latents) return fail("declaration: counts");
  unsigned expected = 1;
  for (auto &entry : groups) {
    if (entry.first != expected++) return fail("declaration: group slots");
    out.groups.push_back(entry.second);
  }
  expected = 1;
  for (auto &entry : latents) {
    if (entry.first != expected++) return fail("declaration: latent slots");
    out.latents.push_back(entry.second);
  }
  if (out.groups.empty()) return fail("declaration: no group");
  Node const *anchor = out.groups[0].gate != nullptr ? out.groups[0].gate : out.groups[0].slab;
  if (setup->parents()[0] != anchor) return fail("declaration: setup anchor");
  // at least one gated group
  bool gated = false;
  for (GroupSpec &g : out.groups) {
    for (StochasticNode *node : g.slab_nodes) {
      if (!available.count(node) || !unobserved(node)) return fail("declaration: slab not free or observed");
    }
    if (!read_slab(g)) return fail("declaration: slab law");
    if (g.gate != nullptr) {
      gated = true;
      if (!available.count(g.gate) || !unobserved(g.gate) || !match_chart(g)) return fail("declaration: gate not free, observed or chart");
    }
  }
  if (!gated) return fail("declaration: no gated group");
  std::set<StochasticNode *> distinct;
  for (GroupSpec const &g : out.groups) {
    for (StochasticNode *node : g.slab_nodes) if (!distinct.insert(node).second) return fail("declaration: duplicate slab");
    if (g.gate != nullptr && !distinct.insert(g.gate).second) return fail("declaration: duplicate gate");
  }
  for (StochasticNode *z : out.latents) {
    if (!available.count(z) || !unobserved(z) || !read_latent(z) || !distinct.insert(z).second) {
      out.latents.clear();
      break;
    }
  }
  // the bias members: a refusal keeps the location members owned
  if (bias_tag != nullptr) {
    if (!bias_tags_valid) out.bias_reason = "two bias tags";
    else out.has_bias = read_bias_declaration(bias_tag, available, out.groups.size(), out.bias, out.bias_reason);
    if (out.has_bias && (distinct.count(out.bias.indicator) ||
                          (out.bias.sources[0] != nullptr && distinct.count(out.bias.sources[0])) ||
                          (out.bias.sources[1] != nullptr && distinct.count(out.bias.sources[1])))) {
      out.has_bias = false;
      out.bias_reason = "a bias member is a location member";
    }
  }
  return true;
}

// The member order of the view: gates, slabs, latents.
std::vector<StochasticNode *> member_list(std::vector<GroupSpec> &groups, std::vector<StochasticNode *> const &latents,
                                          std::vector<unsigned> &latent_pos, unsigned &length)
{
  std::vector<StochasticNode *> members;
  unsigned position = 0;
  for (GroupSpec &g : groups) {
    if (g.gate == nullptr) continue;
    members.push_back(g.gate);
    g.gate_pos = position++;
  }
  for (GroupSpec &g : groups) {
    g.slab_pos = position;
    for (StochasticNode *node : g.slab_nodes) members.push_back(node);
    position += g.d;
  }
  latent_pos.clear();
  for (StochasticNode *z : latents) {
    members.push_back(z);
    latent_pos.push_back(position++);
  }
  length = position;
  return members;
}

// The design X of the children's means in the physical coefficients, probed
// on the graph at construction (chain 0, restored): with every gate on, the
// slabs and latents at zero (the baseline, the non-members at their initial
// values), a unit change of each coordinate. A column counts as verified only
// when its perturbation is representable against the baseline means and its
// response is linear at one and two (column_check); otherwise the block is
// declined and decision() names the check. Guards: the gate off, the latents
// at one and two, and a mixed point agree with the affine prediction.
bool probe_design(GraphView const &view, Shared &S)
{
  const unsigned p = S.p;
  Vec saved(view.length());
  view.getValue(saved, 0);
  Vec v(saved.size());
  auto base = [&]() {
    v = saved;
    for (GroupSpec const &g : S.groups) {
      if (g.gate != nullptr) v[g.gate_pos] = g.label;
      for (unsigned i = 0; i < g.d; ++i) v[g.slab_pos + i] = 0.0;
    }
    for (unsigned pos : S.latent_pos) v[pos] = 0.0;
    // an owned bias label at a label without a slope (its first selection label)
    if (S.bias != nullptr) {
      for (int l = 1; l <= S.bias->levels; ++l) {
        if (S.bias->slope_of[static_cast<std::size_t>(l - 1)] < 0) {
          v[S.bias->pos] = l;
          break;
        }
      }
    }
  };
  auto means = [&](Vec &out) {
    view.setValue(v, 0);
    out.assign(S.rows, 0.0);
    for (ChildSpec const &c : S.children) {
      double const *m = c.node->parents()[0]->value(0);
      for (unsigned i = 0; i < c.n; ++i) out[c.row0 + i] = m[i];
    }
  };
  bool ok = true;
  const double eps = std::numeric_limits<double>::epsilon();
  try {
    Vec mu0, mu1, mu2;
    base();
    means(mu0);
    for (double x : mu0) if (!std::isfinite(x)) ok = fail("design probe: a baseline mean is not finite");
    S.X.assign(static_cast<std::size_t>(S.rows) * p, 0.0);
    S.X_size.assign(static_cast<std::size_t>(S.rows) * p, 0.0);
    // the sum of the columns' scales (column_check), the tolerance unit of the
    // gate-off and mixed points
    double total_scale = 0.0;
    auto agrees = [&](Vec const &observed, Vec const &predicted) {
      for (unsigned r = 0; r < S.rows; ++r) {
        const double rounding = 4.0 * eps * (std::fabs(observed[r]) + std::fabs(predicted[r]));
        if (!(std::fabs(observed[r] - predicted[r]) <= 1e-9 * total_scale + rounding)) return false;
      }
      return true;
    };
    for (GroupSpec const &g : S.groups) {
      for (unsigned i = 0; ok && i < g.d; ++i) {
        const unsigned k = g.coord0 + i;
        base();
        v[g.slab_pos + i] = 1.0;
        means(mu1);
        v[g.slab_pos + i] = 2.0;
        means(mu2);
        double scale = 0.0;
        const ColumnCheck check = column_check(mu0, mu1, mu2, 1.0, &scale);
        if (check == column_unrepresentable && scale == 0.0) {
          // no mean responds at all: a zero column against baseline means
          // whose rounding could hide a response
          ok = fail("design probe: a coordinate's perturbation changes no mean, against nonzero baseline means "
                    "(its column is zero or lost to their rounding)");
        } else if (check == column_unrepresentable) {
          ok = fail("design probe: a coordinate's perturbation is not representable against the baseline means");
        } else if (check == column_nonlinear) {
          ok = fail("design probe: the means are not linear in a coordinate");
        }
        total_scale += scale;
        for (unsigned r = 0; r < S.rows; ++r) {
          S.X[r * p + k] = mu1[r] - mu0[r];
          S.X_size[r * p + k] = std::max(std::fabs(mu0[r]), std::max(std::fabs(mu1[r]), std::fabs(mu2[r])));
        }
      }
      if (ok && g.gate != nullptr) {
        // the gate off hides the slab
        base();
        v[g.gate_pos] = 3.0 - g.label;
        for (unsigned i = 0; i < g.d; ++i) v[g.slab_pos + i] = 1.0;
        means(mu1);
        if (!agrees(mu1, mu0)) ok = fail("design probe: the gate off does not hide its slab");
      }
    }
    Vec dz(S.rows, 0.0);
    if (ok && !S.latents.empty()) {
      base();
      for (unsigned pos : S.latent_pos) v[pos] = 1.0;
      means(mu1);
      for (unsigned pos : S.latent_pos) v[pos] = 2.0;
      means(mu2);
      for (unsigned r = 0; r < S.rows; ++r) dz[r] = mu1[r] - mu0[r];
      for (ChildSpec const &c : S.children) {
        if (c.latent >= 0) continue;
        for (unsigned i = 0; i < c.n; ++i) {
          if (dz[c.row0 + i] != 0.0 || mu2[c.row0 + i] != mu0[c.row0 + i]) {
            ok = fail("design probe: a mean without a latent responds to the latents");
          }
        }
      }
      // The latents' response is a scale that depends on non-members (it is
      // measured at the sweep's own baseline in every update, prepare()):
      // linearity is required where the response is representable here, and
      // a response lost to the rounding of the baseline (or zero) checks
      // nothing that the structure does not.
      double scale = 0.0;
      const ColumnCheck check = column_check(mu0, mu1, mu2, 1.0, &scale);
      if (ok && check == column_nonlinear) ok = fail("design probe: the means are not linear in the latents");
      if (check == column_verified) total_scale += scale;
    }
    if (ok) {
      // mixed: every coordinate and latent at once
      base();
      for (GroupSpec const &g : S.groups) for (unsigned i = 0; i < g.d; ++i) v[g.slab_pos + i] = 1.0;
      for (unsigned pos : S.latent_pos) v[pos] = 1.0;
      means(mu1);
      Vec predicted(mu0);
      for (unsigned r = 0; r < S.rows; ++r) {
        long double x = mu0[r] + static_cast<long double>(dz[r]);
        for (unsigned k = 0; k < p; ++k) x += S.X[r * p + k];
        predicted[r] = static_cast<double>(x);
      }
      if (!agrees(mu1, predicted)) ok = fail("design probe: the means are not additive in the coordinates and latents");
    }
  } catch (...) {
    view.setValue(saved, 0);
    throw;
  }
  view.setValue(saved, 0);
  return ok;
}

std::shared_ptr<Shared> build_shared(Declared &declared, Graph const &graph, bool with_latents, bool with_bias,
                                     std::unique_ptr<GraphView> &view)
{
  auto S = std::make_shared<Shared>();
  S->groups = declared.groups;
  if (with_latents) S->latents = declared.latents;
  unsigned coord = 0;
  for (GroupSpec &g : S->groups) {
    g.coord0 = coord;
    coord += g.d;
    for (unsigned i = 0; i < g.d; ++i) S->coord_group.push_back(static_cast<unsigned>(&g - &S->groups[0]));
  }
  S->p = coord;
  std::vector<StochasticNode *> members = member_list(S->groups, S->latents, S->latent_pos, S->length);
  // the bias members after the location members (LocationBlockBias.h)
  if (with_bias) {
    auto B = std::make_shared<BiasSpec>();
    if (!init_bias(declared.bias, *S, *B)) return nullptr;
    B->pos = S->length++;
    members.push_back(B->indicator);
    for (SlopeSpec &spec : B->slope) {
      if (spec.source == nullptr) continue;
      spec.pos = S->length++;
      members.push_back(spec.source);
    }
    S->bias = B;
  }
  try {
    view.reset(new GraphView(members, graph, false));
  } catch (std::exception const &) {
    return fail_null("view construction");
  }
  if (with_bias && !classify_bias_children(*view, *S->bias)) return nullptr;
  std::vector<int> child_latent;
  bool latents_ok = true;
  if (!verify_structure(*view, graph, S->groups, S->latents, S->bias.get(), child_latent, latents_ok)) return nullptr;
  if (!latents_ok) return fail_null("latents: a child depends on two latents");
  // the likelihood children (gate markers and the nested node of an owned
  // bias label left out)
  std::vector<StochasticNode *> children;
  for (StochasticNode *child : view->stochasticChildren()) {
    if (S->bias == nullptr || !S->bias->skip.count(child)) children.push_back(child);
  }
  unsigned row = 0, woff = 0;
  S->latent_children.assign(S->latents.size(), {});
  for (std::size_t c = 0; c < children.size(); ++c) {
    ChildSpec child;
    bool ok;
    child.node = children[c];
    child.type = child_type(children[c], &ok);
    const unsigned int npar = children[c]->distribution()->npar();
    if (children[c]->parents().size() < npar) return fail_null("child parameters");
    for (unsigned int i = 0; i < npar; ++i) child.lengths.push_back(children[c]->parents()[i]->length());
    child.n = static_cast<unsigned>(children[c]->length());
    if (child.type == ChildType::normal && (child.n != 1 || npar != 2)) return fail_null("normal child shape");
    if (child.type == ChildType::step_switch && child.n != 1) return fail_null("switch child shape");
    if (child.lengths[0] != child.n) return fail_null("child mean length");
    child.row0 = row;
    child.woff = woff;
    child.latent = child_latent[c];
    row += child.n;
    woff += child.n * child.n;
    if (child.latent >= 0) S->latent_children[child.latent].push_back(static_cast<unsigned>(c));
    S->children.push_back(child);
  }
  S->rows = row;
  S->wsize = woff;
  if (!S->latents.empty()) {
    try {
      S->latent_view.reset(new GraphView(S->latents, graph, false));
    } catch (std::exception const &) {
      return fail_null("latent view");
    }
  }
  if (!probe_design(*view, *S)) return nullptr;  // decision() names the check
  if (with_bias && !probe_bias(*view, *S, *S->bias)) return nullptr;
  S->coord_children.assign(S->p, {});
  for (unsigned c = 0; c < S->children.size(); ++c) {
    ChildSpec const &child = S->children[c];
    for (unsigned k = 0; k < S->p; ++k) {
      for (unsigned i = 0; i < child.n; ++i) {
        if (S->X[(child.row0 + i) * S->p + k] != 0.0) {
          S->coord_children[k].push_back(c);
          break;
        }
      }
    }
  }
  return S;
}

class LocationBlockFactory : public SamplerFactory {
public:
  std::string name() const override { return sampler_name; }

  std::vector<Sampler *> makeSamplers(std::list<StochasticNode *> const &free, Graph const &graph) const override
  {
    Declared declared;
    if (!collect(free, graph, declared)) return {};
    std::unique_ptr<GraphView> view;
    // With the bias label first; a bias label that fails its checks leaves the
    // location members owned. If the latents fail their checks, own the gates
    // and slabs only.
    std::shared_ptr<Shared> shared;
    std::string bias_status = declared.bias_reason;
    for (const bool latents : {true, false}) {
      if (!latents && declared.latents.empty()) break;
      if (declared.has_bias) {
        shared = build_shared(declared, graph, latents, true, view);
        if (shared != nullptr) break;
        bias_status = decision();
      }
      shared = build_shared(declared, graph, latents, false, view);
      if (shared != nullptr) break;
    }
    if (shared == nullptr) return {};

    const CorrectedSamplerConfig config = corrected_sampler_configuration();
    shared->corrected = config.enabled;
    shared->settings = config.settings;
    shared->min_acceptance = config.min_acceptance;
    shared->fit_settings = config.settings;
    shared->fit_settings.rule = 1;
    shared->mutant = test_mutant.load(std::memory_order_relaxed);
    // mutant 12: the design column of the last coordinate dropped, so that the
    // graph changes means the block does not track (the untracked-mean error)
    if (shared->mutant == 12) {
      const unsigned k = shared->p - 1;
      for (unsigned r = 0; r < shared->rows; ++r) shared->X[r * shared->p + k] = 0.0;
      shared->coord_children[k].clear();
    }
    // mutant 16: the entry of the first coordinate dropped at the last row that
    // another coordinate also moves (a design entry the probe missed), the
    // children lists as the probe would have left them: the graph moves a
    // mean the block tracks to another value (the tracked-mean error)
    if (shared->mutant == 16) {
      Shared &S = *shared;
      for (unsigned r = S.rows; r-- > 0;) {
        bool other = false;
        for (unsigned k = 1; k < S.p; ++k) other = other || S.X[r * S.p + k] != 0.0;
        if (S.X[r * S.p] == 0.0 || !other) continue;
        S.X[r * S.p] = 0.0;
        S.coord_children[0].clear();
        for (unsigned c = 0; c < S.children.size(); ++c) {
          for (unsigned i = 0; i < S.children[c].n; ++i) {
            if (S.X[(S.children[c].row0 + i) * S.p] != 0.0) {
              S.coord_children[0].push_back(c);
              break;
            }
          }
        }
        break;
      }
    }
    shared->probe = test_probe.load(std::memory_order_relaxed);
    shared->fit_rule = test_fit_rule.load(std::memory_order_relaxed);
    shared->refresh = test_refresh.load(std::memory_order_relaxed);
    shared->timing = test_timing.load(std::memory_order_relaxed);
    {
      const int swap = test_swap.load(std::memory_order_relaxed);
      shared->swap_probability = swap > 0 ? 1.0 / swap : default_swap_probability;
    }
    shared->trace = test_trace.load(std::memory_order_relaxed);
    shared->coin = test_coin.load(std::memory_order_relaxed);
    for (ChildSpec &child : shared->children) {
      child.kind = child_kind(child.node, config.settings);
      if (child.kind != ChildKind::other) shared->any_surrogate = true;
      if (child.kind == ChildKind::cluster || child.kind == ChildKind::factor) shared->coarse_fit = true;
    }

    std::vector<std::unique_ptr<MutableSampleMethod>> methods;
    std::vector<MutableSampleMethod *> raw;
    std::shared_ptr<Shared const> constant = shared;
    for (unsigned int chain = 0; chain < view->nodes()[0]->nchain(); ++chain) {
      methods.push_back(std::make_unique<LocationBlock>(view.get(), constant, chain));
      raw.push_back(methods.back().get());
    }
    std::string name(sampler_name);
    std::unique_ptr<Sampler> sampler(new MutableSampler(view.release(), raw, name));
    for (auto &method : methods) method.release();
    decision() = shared->latents.empty() ? "claimed without latents" : "claimed";
    if (shared->bias != nullptr) decision() += " with the bias label";
    else if (!bias_status.empty()) decision() += "; bias label declined: " + bias_status;
    sampler = location_scale_owner(std::move(sampler), free, graph, shared->latents, this);
    std::vector<Sampler *> out{sampler.get()};
    sampler.release();
    return out;
  }
};

} // namespace

LocationTag::LocationTag(int kind)
  : JointSamplerTag(kind == 0 ? group_tag_name : kind == 1 ? latent_tag_name : kind == 2 ? setup_tag_name : bias_tag_name,
                    kind == 0 ? 5 : kind == 1 ? 2 : kind == 2 ? 3 : 5),
    _kind(kind) {}

bool LocationTag::checkParameterLength(std::vector<unsigned long> const &len) const
{
  switch (_kind) {
  case 0:
    return len.size() == 5 && len[0] == 1 && len[1] >= 1 && len[2] == 1 && len[3] == len[1] && len[4] == 1;
  case 1:
    return len.size() == 2 && len[0] == 1 && len[1] == 1;
  case 3:
    return len.size() == 5 && len[0] == 1 && len[1] >= 2 && len[2] == 1 && len[3] == 1 && len[4] == 1;
  default:
    return len.size() == 3 && len[0] == 1 && len[1] == 1 && len[2] == 1;
  }
}

SamplerFactory *location_block_factory()
{
  static LocationBlockFactory factory;
  return &factory;
}

std::vector<char const *> const &location_block_report_names()
{
  static const std::vector<char const *> names = [] {
    std::vector<char const *> out = {
      "samplers", "sweeps", "sweeps_selection", "sweeps_evidence", "corrections_accepted",
      "corrections_selection_accepted", "corrections_evidence_accepted", "fallback_exact",
      "surrogate_evaluations", "exact_evaluations", "correction_evaluations", "fit_rule_fallback", "ellipse_shrinks",
      "slice_refresh_groups", "slice_refresh_multi", "factorization_skips", "swap_unrepresentable",
      "swap_irreproducible", "unmoved_mean_inexact", "unmoved_included", "max_correction", "max_correction_gaussian",
      "max_endpoint", "max_endpoint_change", "max_stage2_gaussian", "max_swap_stage1_residual",
      "max_swap_full_residual", "max_swap_working", "unmoved_mean_max", "unmoved_density_max"};
    static std::vector<std::string> storage;
    for (int i = 0; i < component_count; ++i) {
      for (char const *stage : {"proposed", "stage0", "stage1", "accepted"}) {
        storage.push_back(std::string(component_names[i]) + "_" + stage);
      }
    }
    for (int i = 0; i < bias_component_count; ++i) {
      for (char const *stage : {"proposed", "stage1", "accepted"}) {
        storage.push_back(std::string("bias_") + bias_component_names[i] + "_" + stage);
      }
    }
    for (int i = 0; i < category_count; ++i) storage.push_back(std::string("evaluations_") + category_names[i]);
    for (int i = 0; i < category_count; ++i) storage.push_back(std::string("seconds_") + category_names[i]);
    for (std::string const &name : storage) out.push_back(name.c_str());
    return out;
  }();
  return names;
}

std::vector<double> location_block_report(bool clear)
{
  Counters &c = counters();
  auto take = [clear](std::atomic<std::uint64_t> &counter) {
    return static_cast<double>(clear ? counter.exchange(0, std::memory_order_relaxed) :
                               counter.load(std::memory_order_relaxed));
  };
  auto take_max = [clear](std::atomic<double> &value) {
    return clear ? value.exchange(0.0, std::memory_order_relaxed) : value.load(std::memory_order_relaxed);
  };
  std::vector<double> out = {
    take(c.samplers), take(c.sweeps), take(c.sweeps_selection), take(c.sweeps_evidence),
    take(c.corrections_accepted), take(c.corrections_selection_accepted),
    take(c.corrections_evidence_accepted), take(c.fallback_exact), take(c.surrogate_evaluations),
    take(c.exact_evaluations), take(c.correction_evaluations), take(c.fit_rule_fallback), take(c.ellipse_shrinks),
    take(c.slice_refresh_groups), take(c.slice_refresh_multi), take(c.factorization_skips),
    take(c.swap_unrepresentable), take(c.swap_irreproducible), take(c.unmoved_mean_inexact),
    take(c.unmoved_included),
    take_max(c.max_correction), take_max(c.max_correction_gaussian), take_max(c.max_endpoint),
    take_max(c.max_endpoint_change), take_max(c.max_stage2_gaussian), take_max(c.max_swap_stage1_residual),
    take_max(c.max_swap_full_residual), take_max(c.max_swap_working), take_max(c.unmoved_mean_max),
    take_max(c.unmoved_density_max)};
  for (int i = 0; i < component_count; ++i) {
    out.push_back(take(c.proposed[i]));
    out.push_back(take(c.stage0[i]));
    out.push_back(take(c.stage1[i]));
    out.push_back(take(c.accepted[i]));
  }
  for (int i = 0; i < bias_component_count; ++i) {
    out.push_back(take(c.bias_proposed[i]));
    out.push_back(take(c.bias_stage1[i]));
    out.push_back(take(c.bias_accepted[i]));
  }
  for (int i = 0; i < category_count; ++i) out.push_back(take(c.evaluations[i]));
  for (int i = 0; i < category_count; ++i) out.push_back(take(c.nanoseconds[i]) * 1e-9);
  return out;
}

int location_block_test_mutant(int mutant)
{
  return test_mutant.exchange(mutant, std::memory_order_relaxed);
}

std::string location_block_test_decision()
{
  return decision();
}

bool location_block_test_probe(int probe)
{
  if (probe < 0) return test_probe.load(std::memory_order_relaxed);
  return test_probe.exchange(probe != 0, std::memory_order_relaxed);
}

std::vector<double> location_block_test_take_probe()
{
  ProbeStore &store = probe_store();
  const std::lock_guard<std::mutex> lock(store.mutex);
  std::vector<double> out;
  out.swap(store.records);
  return out;
}

std::vector<char const *> const &location_block_evidence_names()
{
  static const std::vector<char const *> names = {
    "slot", "chain", "updates", "switched_on", "switched_off", "inner_flips_on_proposed", "inner_flips_on_accepted",
    "inner_flips_off_proposed", "inner_flips_off_accepted", "inner_screened_out", "inner_odds_evaluated",
    "inner_log_odds_sum", "inner_log_odds_max"};
  return names;
}

std::vector<double> location_block_evidence(bool clear)
{
  GateEvidenceStore &store = gate_evidence_store();
  const std::lock_guard<std::mutex> lock(store.mutex);
  std::vector<double> out;
  for (auto const &entry : store.records) {
    GateEvidence const &e = entry.second;
    out.insert(out.end(), {static_cast<double>(entry.first.first), static_cast<double>(entry.first.second),
                           e.updates, e.switched_on, e.switched_off,
                           e.on_proposed, e.on_accepted, e.off_proposed, e.off_accepted, e.screened_out,
                           e.odds_count, e.odds_sum, e.odds_max});
  }
  if (clear) store.records.clear();
  return out;
}


std::vector<char const *> const &location_bias_evidence_names()
{
  static const std::vector<char const *> names = {
    "chain", "label_proposed", "label_accepted", "family_proposed", "family_accepted", "slope_proposed",
    "slope_accepted", "swap_proposed", "swap_stage1", "swap_accepted", "swap_unrepresentable", "swap_irreproducible",
    "leave_selection_attempts", "leave_selection_log_ratio_sum", "leave_selection_log_ratio_max",
    "leave_selection_pass_probability_sum",
    "leave_gaussian_attempts", "leave_gaussian_log_ratio_sum", "leave_gaussian_log_ratio_max",
    "leave_gaussian_pass_probability_sum",
    "kept_leave_selection", "kept_leave_selection_log_ratio_sum", "kept_leave_selection_log_ratio_max",
    "kept_leave_gaussian", "kept_leave_gaussian_log_ratio_sum", "kept_leave_gaussian_log_ratio_max",
    "leave_selection_label_updates", "leave_selection_label_probability_sum", "leave_selection_label_crossings",
    "leave_selection_label_log_ratio_max",
    "leave_gaussian_label_updates", "leave_gaussian_label_probability_sum", "leave_gaussian_label_crossings",
    "leave_gaussian_label_log_ratio_max"};
  return names;
}

std::vector<double> location_bias_evidence(bool clear, std::vector<double> *labels)
{
  BiasEvidenceStore &store = bias_evidence_store();
  const std::lock_guard<std::mutex> lock(store.mutex);
  std::vector<double> out;
  labels->clear();
  for (auto const &entry : store.records) {
    BiasEvidence const &e = entry.second;
    const double chain = static_cast<double>(entry.first);
    out.push_back(chain);
    for (int i = 0; i < 4; ++i) {
      out.push_back(e.proposed[i]);
      if (i == 3) out.push_back(e.stage1[i]);
      out.push_back(e.accepted[i]);
    }
    out.push_back(e.swap_unrepresentable);
    out.push_back(e.swap_irreproducible);
    for (int d = 0; d < 2; ++d) out.insert(out.end(), {e.attempts[d], e.attempt_sum[d], e.attempt_max[d], e.pass_sum[d]});
    for (int d = 0; d < 2; ++d) out.insert(out.end(), {e.kept[d], e.kept_sum[d], e.kept_max[d]});
    for (int d = 0; d < 2; ++d) {
      out.insert(out.end(), {e.label_updates[d], e.label_cross_sum[d], e.label_crossings[d], e.label_log_ratio_max[d]});
    }
    for (std::size_t l = 0; l < e.label_sweeps.size(); ++l) {
      labels->insert(labels->end(), {chain, static_cast<double>(l + 1), e.label_sweeps[l]});
    }
  }
  if (clear) store.records.clear();
  return out;
}

std::vector<int> location_block_test_config(int refresh, int fit_rule, int timing, int swap, int trace, int coin)
{
  std::vector<int> previous = {test_refresh.load(std::memory_order_relaxed),
                               test_fit_rule.load(std::memory_order_relaxed),
                               test_timing.load(std::memory_order_relaxed) ? 1 : 0,
                               test_swap.load(std::memory_order_relaxed),
                               test_trace.load(std::memory_order_relaxed),
                               test_coin.load(std::memory_order_relaxed)};
  if (refresh >= 0) test_refresh.store(refresh, std::memory_order_relaxed);
  if (fit_rule >= 0) test_fit_rule.store(fit_rule, std::memory_order_relaxed);
  if (timing >= 0) test_timing.store(timing != 0, std::memory_order_relaxed);
  if (swap >= 0) test_swap.store(swap, std::memory_order_relaxed);
  if (trace >= 0) test_trace.store(trace, std::memory_order_relaxed);
  if (coin >= 0) test_coin.store(coin, std::memory_order_relaxed);
  return previous;
}

std::vector<double> location_block_test_take_trace()
{
  ProbeStore &store = trace_store();
  const std::lock_guard<std::mutex> lock(store.mutex);
  std::vector<double> out;
  out.swap(store.records);
  return out;
}

std::vector<double> location_block_test_take_sequence()
{
  ProbeStore &store = sequence_store();
  const std::lock_guard<std::mutex> lock(store.mutex);
  std::vector<double> out;
  out.swap(store.records);
  return out;
}

std::vector<double> location_block_test_swap_gate(double log_e_off, double log_e_on)
{
  double log_on, log_off;
  swap_gate_log_probabilities(log_e_off, log_e_on, &log_on, &log_off);
  return {log_on, log_off};
}

double location_block_test_endpoint_difference(std::vector<double> const &a, std::vector<double> const &b)
{
  return endpoint_difference(a, b);
}

}
}
