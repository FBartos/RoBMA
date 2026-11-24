#include "LocationBlock.h"
#include "CoarseSurrogate.h"
#include "CorrectedSampler.h"
#include "../distributions/selnorm-jags-bounds.h"

#include <BayesTools/factor_slice_kernel.h>

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

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
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
constexpr double log_two_pi = 1.8378770664093454836;
const double negative_infinity = -std::numeric_limits<double>::infinity();

// Partner and swap thresholds of the design (absolute correlations).
constexpr double partner_correlation = 0.3;
constexpr double swap_correlation = 0.5;
// Slice widths in proposal standard deviations.
constexpr double slice_width = 2.5;

// ---- per-fit report and test hook

enum Component { flip_on, flip_off, refresh, study, swap, joint, component_count };
const char *const component_names[component_count] = {
  "flip_on", "flip_off", "refresh", "study", "swap", "joint"};

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
  std::atomic<std::uint64_t> proposed[component_count];
  std::atomic<std::uint64_t> stage0[component_count];
  std::atomic<std::uint64_t> stage1[component_count];
  std::atomic<std::uint64_t> accepted[component_count];
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

double normal_cdf(double z)
{
  if (z == std::numeric_limits<double>::infinity()) return 1.0;
  if (z == -std::numeric_limits<double>::infinity()) return 0.0;
  return 0.5 * std::erfc(-z * 0.70710678118654752440);
}

double normal_density(double z)
{
  if (!std::isfinite(z)) return 0.0;
  return std::exp(-0.5 * z * z - 0.91893853320467274178);
}

// The row-marginal selected-normal normalizer of one row, A(mu) =
// sum_b omega_b P(lower_b < sign y / se < upper_b), y ~ N(mu, sd^2), and the
// first two derivatives of log A in mu (spec: the free screen).
bool row_log_normalizer(double mu, double sd, double se, int sign, double const *omega,
                        double const *lower, double const *upper, int bins,
                        double *log_a, double *d1, double *d2)
{
  long double a = 0, a1 = 0, a2 = 0;
  for (int b = 0; b < bins; ++b) {
    const double zu = (upper[b] * se - sign * mu) / sd;
    const double zl = (lower[b] * se - sign * mu) / sd;
    const double probability = normal_cdf(zu) - normal_cdf(zl);
    const double fu = normal_density(zu), fl = normal_density(zl);
    a += static_cast<long double>(omega[b]) * probability;
    a1 += static_cast<long double>(omega[b]) * (-sign * (fu - fl) / sd);
    a2 += static_cast<long double>(omega[b]) *
      (-((std::isfinite(zu) ? zu * fu : 0.0) - (std::isfinite(zl) ? zl * fl : 0.0)) / (sd * sd));
  }
  if (!(a > 0)) return false;
  *log_a = std::log(static_cast<double>(a));
  *d1 = static_cast<double>(a1 / a);
  *d2 = static_cast<double>(a2 / a) - (*d1) * (*d1);
  return std::isfinite(*log_a) && std::isfinite(*d1) && std::isfinite(*d2);
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
  std::unique_ptr<GraphView> latent_view;  // latents only, for the per-update probe
  SelNormCoarseSettings settings;
  bool corrected = false;                  // the surrogate is configured
  bool any_surrogate = false;              // some child has a surrogate evaluator
  double min_acceptance = 0.0;
  int mutant = 0;
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
                      std::vector<StochasticNode *> const &latents,
                      std::vector<int> &child_latent, bool &latents_ok)
{
  std::set<Node const *> coef, zdep, gated, raw, members;
  std::map<Node const *, int> latent_id;
  std::set<Node const *> chart_nodes, physical_nodes;
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
    int id = -1;
    for (Node const *parent : node->parents()) {
      if (coef.count(parent)) from_coef = true;
      if (gated.count(parent)) from_gate = true;
      if (raw.count(parent)) from_raw = true;
      if (zdep.count(parent)) {
        from_z = true;
        const int other = latent_id[parent];
        if (other == -2 || (id >= 0 && other != id)) id = -2;
        else if (id == -1) id = other;
      }
    }
    if (!reaches[i]) continue;
    if (from_gate || from_raw) return fail("structure: gate or raw slab outside its chart");
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

  latents_ok = true;
  child_latent.clear();
  for (StochasticNode const *child : view.stochasticChildren()) {
    if (members.count(child)) return fail("structure: a member is a child");
    bool ok;
    child_type(child, &ok);
    if (!ok || !isObserved(child) || child->lowerBound() != nullptr || child->upperBound() != nullptr) {
      return fail("structure: child distribution, observation or bounds");
    }
    auto const &par = child->parents();
    for (std::size_t i = 1; i < par.size(); ++i) {
      Node const *parent = par[i];
      if (members.count(parent) || coef.count(parent) || zdep.count(parent) || gated.count(parent) ||
          raw.count(parent)) return fail("structure: member enters a non-mean parameter");
    }
    Node const *mean = par[0];
    if (gated.count(mean) || raw.count(mean)) return fail("structure: mean depends on a gate or raw slab");
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
  std::uint64_t _adaptive_evidence = 0, _adaptive_accepted = 0;
  std::uint64_t _window_proposed = 0, _window_accepted = 0;
  // the refresh kernel of each group in selection labels: the tilted
  // independence proposal, or (acceptance below 0.4 while adapting) the slice
  // update along the line, decided at the first adaptOff()
  std::vector<std::uint64_t> _refresh_proposed, _refresh_accepted;
  std::vector<char> _slice_refresh;

  // ---- per-update state
  Vec _entry, _final;
  std::vector<char> _g;            // gate on
  Vec _beta, _b, _z;               // slabs, physical coefficients, latents
  Vec _y, _mean_entry, _mean, _r, _wr, _d, _W;
  std::vector<char> _selected, _value_ok, _dirty;
  Vec _value, _value_entry;        // inner-target terms of the selection children
  std::vector<unsigned> _stamp;
  unsigned _stamp_id = 0;
  bool _selection = false;         // some child is a selection child in step mode
  bool _surrogate_mode = false;    // the surrogate replaces the exact selection densities
  bool _not_exact = false;         // a surrogate evaluated this update was not exact
  long double _repr_delta = 0;     // accumulated log density change of the inner target
  Vec _a, _lambda, _u, _P, _h;     // Gaussian working model
  std::vector<double const *> _params;
  Vec _scratch;

  Shared const &s() const { return *_s; }

  // ---- child densities

  void read_params(ChildSpec const &c)
  {
    auto const &par = c.node->parents();
    _params.resize(c.lengths.size());
    for (std::size_t i = 0; i < c.lengths.size(); ++i) _params[i] = par[i]->value(_chain);
  }

  // The child's own exact density at the given mean (PDF_LIKELIHOOD).
  double exact_with_mean(ChildSpec const &c, double const *mean)
  {
    read_params(c);
    _params[0] = mean;
    auto const *dist = static_cast<VectorDist const *>(c.node->distribution());
    count(counters().exact_evaluations);
    if (c.type == ChildType::normal) {
      const double tau = *_params[1];
      const double x = c.node->value(_chain)[0];
      return 0.5 * (std::log(tau) - log_two_pi) - 0.5 * tau * (x - mean[0]) * (x - mean[0]);
    }
    if (!dist->checkParameterValue(_params, c.lengths)) return negative_infinity;
    return dist->logDensity(c.node->value(_chain), PDF_LIKELIHOOD, _params, c.lengths);
  }

  double surrogate_with_mean(unsigned child, double const *mean)
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
    if (!valid) return negative_infinity;
    return sc.surrogate(_chain, s().settings, &_not_exact);
  }

  // The inner target's term of a selection child at a mean.
  double selected_value(unsigned child, double const *mean)
  {
    const double v = (_surrogate_mode && _surrogates[child]) ? surrogate_with_mean(child, mean) :
      exact_with_mean(s().children[child], mean);
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
    // A child is first evaluated before any move changes it: its entry value.
    _value[child] = selected_value(child, &_mean[c.row0]);
    _value_entry[child] = _value[child];
    _value_ok[child] = 1;
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

  bool is_selected(ChildSpec const &c)
  {
    auto const &par = c.node->parents();
    switch (c.type) {
    case ChildType::normal:
      return false;
    case ChildType::step_switch:
      return par[9]->value(_chain)[0] != SELKERNEL_NORMAL && par[7]->value(_chain)[0] != 0;
    case ChildType::cluster:
    case ChildType::factor:
      return par[10]->value(_chain)[0] != SELKERNEL_NORMAL;
    case ChildType::mv:
      return par[9]->value(_chain)[0] != SELKERNEL_NORMAL;
    }
    return true;
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
      _selected[c] = is_selected(child);
      if (_selected[c]) _selection = true;
      _value_ok[c] = 0;
      _dirty[c] = 0;
    }
    _surrogate_mode = _selection && sh.corrected && !_exact;
    _line.valid = false;
    _not_exact = false;
    _repr_delta = 0;

    // d: the change of each child's mean per unit of its latent (linear,
    // verified at construction), probed at the current non-members.
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
    for (unsigned i = 0; i < sh.rows; ++i) _r[i] = _y[i] - _mean[i];
    for (ChildSpec const &child : sh.children) {
      double const *W = &_W[child.woff];
      for (unsigned i = 0; i < child.n; ++i) {
        long double v = 0;
        for (unsigned j = 0; j < child.n; ++j) v += static_cast<long double>(W[i * child.n + j]) * _r[child.row0 + j];
        _wr[child.row0 + i] = static_cast<double>(v);
      }
    }

    // Gaussian working model with the latents integrated: P = B + X' S^-1 X,
    // h = B m0 + X' S^-1 r0, r0 = y - offset; per latent the within/between
    // form (spec: Numerics), never a difference of two large terms.
    std::fill(_P.begin(), _P.end(), 0.0);
    std::fill(_h.begin(), _h.end(), 0.0);
    std::fill(_u.begin(), _u.end(), 0.0);
    std::fill(_a.begin(), _a.end(), 0.0);
    Vec r0(sh.rows), WX, Wd, Wr0;
    for (ChildSpec const &child : sh.children) {
      for (unsigned i = 0; i < child.n; ++i) {
        const unsigned row = child.row0 + i;
        long double v = _r[row];
        for (unsigned k = 0; k < p; ++k) v += static_cast<long double>(sh.X[row * p + k]) * _b[k];
        if (child.latent >= 0) v += static_cast<long double>(_d[row]) * _z[child.latent];
        r0[row] = static_cast<double>(v);
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
      // a = d' W d, u = d' W X, rho = d' W r0
      std::vector<std::pair<Vec, Vec>> wxs;
      std::vector<Vec> wds, wr0s;
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
        wxs.emplace_back(WX, Vec());
        wds.push_back(Wd);
        wr0s.push_back(Wr0);
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
        Vec const &wx = wxs[q].first;
        Vec const &wd = wds[q];
        Vec const &wr0 = wr0s[q];
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
    std::vector<double> dc(p, 0.0);
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
    std::vector<double> dz(S.latents.size(), 0.0);
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
    std::sort(sh.children.begin(), sh.children.end());
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
    if (s().mutant == 2) {
      // the reverse-density mutant: the dormant redraw is twice too wide
      for (unsigned i = 0; i < g.d; ++i) out[i] = g.m0[i] + 2.0 * (out[i] - g.m0[i]);
    }
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
  void gaussian_flip(unsigned j, RNG *rng)
  {
    Shared const &S = s();
    GroupSpec const &g = S.groups[j];
    Line L;
    if (!build_line({j}, L)) return;
    const bool on = active(j);
    const Component type = on ? flip_off : flip_on;
    count(counters().proposed[type]);
    count(counters().stage0[type]);
    // Savage-Dickey form of the collapsed odds of the line
    Vec zero(g.d, 0.0);
    const double log_odds = (g.log_on - g.log_off) + log_slab(g, zero.data()) -
      log_normal_precision(zero.data(), L.m, L.Pchol, L.log_det_P);
    const double stage1 = on ? -log_odds : log_odds;
    if (!metropolis(rng, stage1)) return;
    count(counters().stage1[type]);
    Vec proposal(g.d), db(g.d);
    double const *old_beta = &_beta[g.coord0];
    long double proposal_terms;
    if (on) {
      draw_prior(g, rng, proposal.data());
      for (unsigned i = 0; i < g.d; ++i) db[i] = -old_beta[i];
      // reverse: the line law at the current slab; forward: the dormant prior
      proposal_terms = log_normal_precision(old_beta, L.m, L.Pchol, L.log_det_P) -
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
    const long double target = children + prior_delta(sh, 1.0) +
      (on ? g.log_off - g.log_on : g.log_on - g.log_off) +
      log_slab(g, proposal.data()) - log_slab(g, old_beta);
    const double full = static_cast<double>(target + proposal_terms);
    const double stage2 = full - stage1;
    if (!_selection) record_max(counters().max_stage2_gaussian, std::fabs(stage2));
    if (!std::isfinite(full) || !metropolis(rng, stage2)) return;
    count(counters().accepted[type]);
    apply_shift(sh, 1.0);
    _g[j] = on ? 0 : 1;
    std::copy(proposal.begin(), proposal.end(), _beta.begin() + g.coord0);
    for (unsigned i = 0; i < g.d; ++i) _b[g.coord0 + i] = on ? 0.0 : proposal[i];
    _repr_delta += target;
  }

  // R_j of an active group in a Gaussian label: an exact draw along the line.
  void gaussian_refresh(unsigned j, RNG *rng)
  {
    GroupSpec const &g = s().groups[j];
    Line L;
    if (!build_line({j}, L)) return;
    count(counters().proposed[refresh]);
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
    const long double target = children_delta(sh, 1.0) + prior_delta(sh, 1.0);
    const double full = static_cast<double>(target +
      log_normal_precision(old_beta, L.m, L.Pchol, L.log_det_P) -
      log_normal_precision(proposal.data(), L.m, L.Pchol, L.log_det_P));
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
  // each coordinate along its sheared line, width 2.5 line SDs.
  void slice_refresh(unsigned j, RNG *rng)
  {
    GroupSpec const &g = s().groups[j];
    Line L;
    if (!build_line({j}, L)) return;
    const std::size_t nm = L.moving.size();
    std::vector<unsigned> prior_coords(L.partners);
    prior_coords.push_back(g.coord0);
    for (std::size_t i = 0; i < nm; ++i) {
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
  double line_ell(Line const &L, double const *db, Shift &sh)
  {
    line_shift(L, db, sh);
    prepare_shift(sh, {});
    const long double children = children_delta(sh, 1.0);
    if (!(children > negative_infinity)) return negative_infinity;
    long double gaussian = 0;
    for (std::size_t q = 0; q < sh.children.size(); ++q) gaussian += sh.q1[q] - 0.5L * sh.q2[q];
    return static_cast<double>(children - gaussian);
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

  // The damped Newton fit of f(x) = log N(x; m, P^-1) + ell(x) along the
  // subline of the free coordinates (spec: two iterations from m, curvature
  // floored at 0.25 P, at most five halvings, central differences with steps
  // of a quarter of the law's SD; the full Hessian up to four coordinates and
  // its diagonal beyond). A deterministic function of the line.
  Fit newton(Line const &L, std::size_t first, Vec const &target, Vec const &m, Vec const &P,
             Vec const &Pchol, double log_det_P)
  {
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
      bool moved = false;
      double scale = 1.0;
      for (int halving = 0; halving <= 5 && !moved; ++halving, scale *= 0.5) {
        Vec next(x);
        for (std::size_t i = 0; i < d; ++i) next[i] += scale * step[i];
        const double en = ell_at(next);
        if (std::isfinite(en) && f_at(next, en) >= f0) {
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
    std::map<unsigned, Vec> deltas;  // child -> rows x moving coordinates
    for (std::size_t m = 0; m < nm; ++m) {
      Shift const &b = basis[m];
      for (std::size_t q = 0; q < b.children.size(); ++q) {
        const unsigned c = b.children[q];
        if (gaussian_child(c)) continue;
        Vec &rows = deltas[c];
        if (rows.empty()) rows.assign(S.children[c].n * nm, 0.0);
        for (unsigned i = 0; i < S.children[c].n; ++i) rows[i * nm + m] = b.delta[b.offset[q] + i];
      }
    }
    long double total = 0;
    grad.assign(d, 0.0);
    hess.assign(d * d, 0.0);
    RowModel row;
    for (auto const &entry : deltas) {
      ChildSpec const &c = S.children[entry.first];
      Vec const &rows = entry.second;
      for (unsigned i = 0; i < c.n; ++i) {
        if (!row_model(c, i, row) || row.bin == 0) continue;
        long double mu = _mean[c.row0 + i];
        for (std::size_t m = 0; m < nm; ++m) mu += static_cast<long double>(db[m]) * rows[i * nm + m];
        double log_a, d1, d2;
        if (!row_log_normalizer(static_cast<double>(mu), row.sd, row.se, row.sign, row.omega,
                                row.lower.data(), row.upper.data(), row.bins, &log_a, &d1, &d2)) return false;
        total -= static_cast<long double>(row.power) * log_a;
        for (std::size_t k = 0; k < d; ++k) {
          const double ak = rows[i * nm + first + k];
          grad[k] -= row.power * d1 * ak;
          for (std::size_t l = 0; l < d; ++l) hess[k * d + l] -= row.power * d2 * ak * rows[i * nm + first + l];
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
    double log_odds = 0;         // Laplace log odds (prior odds included)
  };
  GroupLine _line;

  bool group_line(unsigned j)
  {
    if (_line.valid && _line.group == j) return true;
    _line = GroupLine();
    _line.group = j;
    if (!build_line({j}, _line.L)) return false;
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
    std::vector<Shift> basis;
    line_basis(_line.L, basis);
    Vec target(_line.L.m), zero(_line.L.m.size(), 0.0);
    _line.log_screen = clamp_screen(_line.log_odds_gauss + screen_tilt(
      basis, _line.L, 0, target, zero, _line.L.m, _line.L.P, _line.L.Pchol, _line.L.log_det_P));
    _line.screened = true;
  }

  void group_fit()
  {
    if (_line.fitted) return;
    Line const &L = _line.L;
    Vec target(L.m);
    _line.fit = newton(L, 0, target, L.m, L.P, L.Pchol, L.log_det_P);
    // ell at the physical zero of the group (the off state)
    Vec zero(L.m.size(), 0.0), db;
    Shift sh;
    subline_db(L, 0, zero, zero, db);
    const double ell_zero = line_ell(L, db.data(), sh);
    if (!_line.fit.ok || !std::isfinite(ell_zero)) {
      // a deterministic fallback of the line: the Gaussian law, its ell at m
      Fit fallback;
      fallback.xhat = L.m;
      fallback.Pt = L.P;
      fallback.Ptchol = L.Pchol;
      fallback.log_det_Pt = L.log_det_P;
      subline_db(L, 0, L.m, L.m, db);
      fallback.ell_hat = line_ell(L, db.data(), sh);
      fallback.ok = std::isfinite(fallback.ell_hat) && std::isfinite(ell_zero);
      _line.fit = fallback;
    }
    _line.log_odds = _line.fit.ok ?
      _line.log_odds_gauss + laplace_tilt(_line.fit, L.m, L.Pchol, L.log_det_P, ell_zero) : negative_infinity;
    _line.fitted = true;
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
    group_screen();
    const double stage0 = on ? -_line.log_screen : _line.log_screen;
    if (!metropolis(rng, stage0)) return;
    count(counters().stage0[type]);
    group_fit();
    if (!_line.fit.ok) return;
    const double log_odds = _line.log_odds;
    const double stage1 = on ? -(log_odds - _line.log_screen) : (log_odds - _line.log_screen);
    if (!metropolis(rng, stage1)) return;
    count(counters().stage1[type]);
    Fit const &fit = _line.fit;
    Vec proposal(g.d), db(g.d);
    double const *old_beta = &_beta[g.coord0];
    long double proposal_terms;
    if (on) {
      draw_prior(g, rng, proposal.data());
      for (unsigned i = 0; i < g.d; ++i) db[i] = -old_beta[i];
      proposal_terms = log_normal_precision(old_beta, fit.xhat, fit.Ptchol, fit.log_det_Pt) -
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
    if (!build_line({j, k}, L)) return;
    count(counters().proposed[swap]);
    const std::size_t dj = gj.d, dk = gk.d, nm = dj + dk;
    // the conditional Gaussian laws of each group with the other at zero, and
    // the Gaussian log ratio of the pattern evidences
    Vec cj, Cj, ck, Ck, mj, Pj, mk, Pk;
    double lcj, lck, lj, lk;
    if (!swap_laws(L, dj, dk, mj, Pj, lj, mk, Pk, lk, cj, Cj, lcj, ck, Ck, lck)) return;
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
    if (!metropolis(rng, screen)) return;
    count(counters().stage0[swap]);
    const Fit fa = newton(L, 0, target_a, cj, Cjp, Cj, lcj);
    const Fit fb = newton(L, dj, target_b, ck, Ckp, Ck, lck);
    Vec db;
    Shift sh;
    subline_db(L, 0, zero_j, zeros, db);
    const double ell_zero = line_ell(L, db.data(), sh);
    if (!fa.ok || !fb.ok || !std::isfinite(ell_zero)) return;
    const double log_odds = gauss + laplace_tilt(fb, ck, Ck, lck, ell_zero) - laplace_tilt(fa, cj, Cj, lcj, ell_zero);
    if (!metropolis(rng, log_odds - screen)) return;
    count(counters().stage1[swap]);
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
  void gaussian_latent(unsigned lat, RNG *rng)
  {
    Shared const &S = s();
    count(counters().proposed[study]);
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
    const long double target = children_delta(sh, 1.0) + prior_delta(sh, 1.0);
    const double zo = (_z[lat] - mean) / sd, zn = (proposal - mean) / sd;
    const double full = static_cast<double>(target) - 0.5 * zo * zo + 0.5 * zn * zn;
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
  void gaussian_joint(RNG *rng)
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
    count(counters().proposed[joint]);
    Vec PA(na * na), hA(na);
    for (std::size_t i = 0; i < na; ++i) {
      hA[i] = _h[A[i]];
      for (std::size_t k = 0; k < na; ++k) PA[i * na + k] = _P[A[i] * p + A[k]];
    }
    Vec chol(PA);
    if (!cholesky(chol, na)) return;
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
    const long double target = children_delta(sh, 1.0) + prior_delta(sh, 1.0);
    const double full = static_cast<double>(target + log_normal_precision(old.data(), mean, chol, log_det) -
      log_normal_precision(proposal.data(), mean, chol, log_det));
    if (!_selection) record_max(counters().max_stage2_gaussian, std::fabs(full));
    if (!std::isfinite(full) || !metropolis(rng, full)) return;
    count(counters().accepted[joint]);
    apply_shift(sh, 1.0);
    _repr_delta += target;
  }

  // S_jk: exchange the gates of two correlated groups, the one switched on
  // drawn from its law along the line given the other off (Gaussian labels).
  void gaussian_swap(unsigned first, unsigned second, RNG *rng)
  {
    Shared const &S = s();
    if (active(first) == active(second)) return;
    const unsigned j = active(first) ? first : second;   // on -> off
    const unsigned k = active(first) ? second : first;   // off -> on
    GroupSpec const &gj = S.groups[j];
    GroupSpec const &gk = S.groups[k];
    Line L;
    if (!build_line({j, k}, L)) return;
    count(counters().proposed[swap]);
    count(counters().stage0[swap]);
    const std::size_t dj = gj.d, dk = gk.d, nm = dj + dk;
    // marginal laws of each group (from the joint of the pair) and the
    // conditionals of one given the other at zero
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
    // conditional of block (r0, n0) given the other block (c0, n1) = 0
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
    Vec mj, Pj, mk, Pk, cj, Cj, ck, Ck;
    double lj, lk, lcj, lck;
    if (!marginal(0, dj, mj, Pj, lj) || !marginal(dj, dk, mk, Pk, lk) ||
        !conditional(0, dj, dj, dk, cj, Cj, lcj) || !conditional(dj, dk, 0, dj, ck, Ck, lck)) return;
    Vec zero_j(dj, 0.0), zero_k(dk, 0.0);
    const double stage1 = (gj.log_off + gk.log_on - gj.log_on - gk.log_off) +
      log_slab(gk, zero_k.data()) + log_normal_precision(zero_j.data(), mj, Pj, lj) -
      log_slab(gj, zero_j.data()) - log_normal_precision(zero_k.data(), mk, Pk, lk);
    if (!metropolis(rng, stage1)) return;
    count(counters().stage1[swap]);
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
    const long double target = children_delta(sh, 1.0) + prior_delta(sh, 1.0) +
      (gj.log_off + gk.log_on - gj.log_on - gk.log_off) +
      log_slab(gj, dormant_j.data()) - log_slab(gj, old_j) + log_slab(gk, beta_k.data()) - log_slab(gk, old_k);
    const long double proposals = log_normal_precision(old_j, cj, Cj, lcj) + log_slab(gk, old_k) -
      log_normal_precision(beta_k.data(), ck, Ck, lck) - log_slab(gj, dormant_j.data());
    const double full = static_cast<double>(target + proposals);
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
    if (!cholesky(chol, p)) return pairs;
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

  // ---- the update

  enum Kind { kind_flip, kind_refresh, kind_swap, kind_joint, kind_latent };
  struct Step { Kind kind; unsigned a, b; };

  void sweep(RNG *rng)
  {
    Shared const &S = s();
    prepare();
    count(counters().sweeps);
    if (_selection) count(counters().sweeps_selection);
    std::vector<Step> steps;
    for (unsigned j = 0; j < S.groups.size(); ++j) {
      if (S.groups[j].gate != nullptr) steps.push_back({kind_flip, j, 0});
      steps.push_back({kind_refresh, j, 0});
    }
    for (auto const &pair : swap_pairs()) steps.push_back({kind_swap, pair.first, pair.second});
    // the joint refresh is Gaussian only: the untilted joint law is rejected
    // under selection
    if (!_selection && !S.groups.empty()) steps.push_back({kind_joint, 0, 0});
    for (unsigned lat = 0; lat < S.latents.size(); ++lat) steps.push_back({kind_latent, lat, 0});
    // one coin per update: the whole list forward or reversed
    const bool backward = rng->uniform() < 0.5;
    for (std::size_t i = 0; i < steps.size(); ++i) {
      Step const &step = steps[backward ? steps.size() - 1 - i : i];
      // The line of a group serves its F_j and R_j, adjacent in either order;
      // any other component changes it.
      if ((step.kind != kind_flip && step.kind != kind_refresh) || _line.group != step.a) _line.valid = false;
      switch (step.kind) {
      case kind_flip:
        if (_selection) tilted_flip(step.a, rng);
        else gaussian_flip(step.a, rng);
        break;
      case kind_refresh:
        if (!active(step.a)) dormant_refresh(step.a, rng);
        else if (!_selection) gaussian_refresh(step.a, rng);
        else if (_slice_refresh[step.a]) slice_refresh(step.a, rng);
        else tilted_refresh(step.a, rng);
        break;
      case kind_swap:
        if (_selection) tilted_swap(step.a, step.b, rng);
        else gaussian_swap(step.a, step.b, rng);
        break;
      case kind_joint:
        gaussian_joint(rng);
        break;
      case kind_latent:
        if (_selection) slice_latent(step.a, rng);
        else gaussian_latent(step.a, rng);
        break;
      }
    }
    correct(rng);
  }

  // The exact correction: accept the sweep with [pi(new) / pi~(new)] /
  // [pi(old) / pi~(old)] over the changed members and the children whose
  // means changed, the exact terms from the graph's own densities.
  void correct(RNG *rng)
  {
    Shared const &S = s();
    // the final member values
    _final = _entry;
    for (unsigned j = 0; j < S.groups.size(); ++j) {
      GroupSpec const &g = S.groups[j];
      if (g.gate != nullptr) _final[g.gate_pos] = _g[j] ? g.label : 3.0 - g.label;
      for (unsigned i = 0; i < g.d; ++i) _final[g.slab_pos + i] = _beta[g.coord0 + i];
    }
    for (unsigned lat = 0; lat < S.latents.size(); ++lat) _final[S.latent_pos[lat]] = _z[lat];
    if (_final == _entry) {
      count(counters().corrections_accepted);
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
    for (unsigned c = 0; c < S.children.size(); ++c) if (_dirty[c]) kids.push_back(c);
    // the exact terms at the entry (the graph still holds it), and the
    // endpoint differences of the selection children, sum (log A - log A~)
    // up to sign: inner-target value minus exact value
    long double old_exact = 0, new_exact = 0, old_endpoint = 0, new_endpoint = 0;
    for (unsigned m : changed_members) old_exact += members[m]->logDensity(_chain, PDF_PRIOR);
    for (unsigned c : kids) {
      const double v = S.children[c].node->logDensity(_chain, PDF_LIKELIHOOD);
      count(counters().correction_evaluations);
      old_exact += v;
      if (!gaussian_child(c)) old_endpoint += static_cast<long double>(_value_entry[c]) - v;
    }
    _view->setValue(_final, _chain);
    bool hidden = false;
    for (unsigned c = 0; c < S.children.size() && !hidden; ++c) {
      if (_dirty[c]) continue;
      ChildSpec const &child = S.children[c];
      double const *mean = child.node->parents()[0]->value(_chain);
      for (unsigned i = 0; i < child.n; ++i) {
        if (mean[i] != _mean_entry[child.row0 + i]) hidden = true;
      }
    }
    for (unsigned m : changed_members) new_exact += members[m]->logDensity(_chain, PDF_PRIOR);
    for (unsigned c : kids) {
      const double v = S.children[c].node->logDensity(_chain, PDF_LIKELIHOOD);
      count(counters().correction_evaluations);
      new_exact += v;
      if (!gaussian_child(c)) new_endpoint += static_cast<long double>(_value[c]) - v;
    }
    if (std::isnan(static_cast<double>(new_exact))) {
      throwNodeError(members[0], "Location block sampler: the log density of the proposed state is not defined");
    }
    const double log_ratio = static_cast<double>((new_exact - old_exact) - _repr_delta);
    bool accepted;
    if (hidden || !(new_exact > negative_infinity)) {
      accepted = false;
    } else if (S.mutant == 3) {
      accepted = true;
    } else {
      accepted = metropolis(rng, log_ratio);
    }
    if (std::isfinite(log_ratio)) {
      record_max(_selection ? counters().max_correction : counters().max_correction_gaussian, std::fabs(log_ratio));
    }
    if (_selection) {
      record_max(counters().max_endpoint, std::fabs(static_cast<double>(new_endpoint)));
      record_max(counters().max_endpoint_change, std::fabs(static_cast<double>(new_endpoint - old_endpoint)));
    }
    if (accepted) {
      count(counters().corrections_accepted);
      if (_selection) count(counters().corrections_selection_accepted);
    } else {
      _view->setValue(_entry, _chain);
    }
    // evidence about the surrogate: a selection label with a surrogate that
    // was not exact
    const bool evidence = _surrogate_mode && _not_exact;
    if (_surrogate_mode) {
      corrected_sampler_counters.coarse_proposed.fetch_add(1, std::memory_order_relaxed);
      if (accepted) corrected_sampler_counters.coarse_accepted.fetch_add(1, std::memory_order_relaxed);
      else corrected_sampler_counters.coarse_correction_rejected.fetch_add(1, std::memory_order_relaxed);
    }
    if (!evidence) return;
    count(counters().sweeps_evidence);
    if (accepted) count(counters().corrections_evidence_accepted);
    if (_adapting) {
      ++_adaptive_evidence;
      if (accepted) ++_adaptive_accepted;
    } else if (_registered) {
      ++_window_proposed;
      corrected_sampler_counters.report_sampling_proposed.fetch_add(1, std::memory_order_relaxed);
      if (accepted) {
        ++_window_accepted;
        corrected_sampler_counters.report_sampling_accepted.fetch_add(1, std::memory_order_relaxed);
      }
    }
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
    _W.assign(S.wsize, 0.0);
    _selected.assign(S.children.size(), 0);
    _value_ok.assign(S.children.size(), 0);
    _dirty.assign(S.children.size(), 0);
    _value.assign(S.children.size(), 0.0);
    _value_entry.assign(S.children.size(), 0.0);
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
      register_corrected_sampler(this);
      _registered = true;
    }
  }

  ~LocationBlock() override
  {
    if (_registered) unregister_corrected_sampler(this);
  }

  std::uint64_t windowProposed() const override { return _window_proposed; }
  std::uint64_t windowAccepted() const override { return _window_accepted; }
  void clearWindow() override
  {
    _window_proposed = 0;
    _window_accepted = 0;
  }
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
    for (std::size_t j = 0; j < _slice_refresh.size(); ++j) {
      _slice_refresh[j] = _refresh_proposed[j] > 0 &&
        static_cast<double>(_refresh_accepted[j]) < 0.4 * static_cast<double>(_refresh_proposed[j]);
    }
    if (!_registered || !(s().min_acceptance > 0.0)) return;
    if (_adaptive_evidence >= corrected_min_adaptive_proposals &&
        static_cast<double>(_adaptive_accepted) / static_cast<double>(_adaptive_evidence) < s().min_acceptance) {
      _exact = true;
      corrected_sampler_counters.report_fallback_low_acceptance.fetch_add(1, std::memory_order_relaxed);
      count(counters().fallback_exact);
    }
  }

  bool checkAdaptation() const override { return true; }
};

// ---- the factory

struct Declared {
  std::vector<GroupSpec> groups;
  std::vector<StochasticNode *> latents;
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
// slabs and latents at zero, a unit change of each coordinate. Guards: a
// second point of each coordinate, the gate off, the latents at one and two,
// and a mixed point must agree with the affine prediction.
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
  try {
    Vec mu0, mu1, mu2;
    base();
    means(mu0);
    for (double x : mu0) if (!std::isfinite(x)) ok = false;
    S.X.assign(static_cast<std::size_t>(S.rows) * p, 0.0);
    Vec scale(S.rows, 1.0);
    for (unsigned i = 0; i < S.rows; ++i) scale[i] += std::fabs(mu0[i]);
    auto close = [&](Vec const &observed, Vec const &predicted) {
      for (unsigned i = 0; i < S.rows; ++i) {
        if (!(std::fabs(observed[i] - predicted[i]) <= 1e-9 * (scale[i] + std::fabs(predicted[i])))) return false;
      }
      return true;
    };
    for (GroupSpec const &g : S.groups) {
      for (unsigned i = 0; ok && i < g.d; ++i) {
        const unsigned k = g.coord0 + i;
        base();
        v[g.slab_pos + i] = 1.0;
        means(mu1);
        Vec predicted(mu0);
        for (unsigned r = 0; r < S.rows; ++r) {
          S.X[r * p + k] = mu1[r] - mu0[r];
          if (!std::isfinite(S.X[r * p + k])) ok = false;
          predicted[r] = mu0[r] + 2.0 * S.X[r * p + k];
        }
        v[g.slab_pos + i] = 2.0;
        means(mu2);
        if (!close(mu2, predicted)) ok = false;
      }
      if (ok && g.gate != nullptr) {
        // the gate off hides the slab
        base();
        v[g.gate_pos] = 3.0 - g.label;
        for (unsigned i = 0; i < g.d; ++i) v[g.slab_pos + i] = 1.0;
        means(mu1);
        if (!close(mu1, mu0)) ok = false;
      }
    }
    Vec dz(S.rows, 0.0);
    if (ok && !S.latents.empty()) {
      base();
      for (unsigned pos : S.latent_pos) v[pos] = 1.0;
      means(mu1);
      Vec predicted(mu0);
      for (unsigned r = 0; r < S.rows; ++r) {
        dz[r] = mu1[r] - mu0[r];
        predicted[r] = mu0[r] + 2.0 * dz[r];
      }
      for (ChildSpec const &c : S.children) {
        if (c.latent >= 0) continue;
        for (unsigned i = 0; i < c.n; ++i) if (dz[c.row0 + i] != 0.0) ok = false;
      }
      for (unsigned pos : S.latent_pos) v[pos] = 2.0;
      means(mu2);
      if (!close(mu2, predicted)) ok = false;
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
      if (!close(mu1, predicted)) ok = false;
    }
  } catch (...) {
    view.setValue(saved, 0);
    throw;
  }
  view.setValue(saved, 0);
  return ok;
}

std::shared_ptr<Shared> build_shared(Declared &declared, Graph const &graph, bool with_latents,
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
  const std::vector<StochasticNode *> members = member_list(S->groups, S->latents, S->latent_pos, S->length);
  try {
    view.reset(new GraphView(members, graph, false));
  } catch (std::exception const &) {
    return fail_null("view construction");
  }
  std::vector<int> child_latent;
  bool latents_ok = true;
  if (!verify_structure(*view, graph, S->groups, S->latents, child_latent, latents_ok)) return nullptr;
  if (!latents_ok) return fail_null("latents: a child depends on two latents");
  std::vector<StochasticNode *> const &children = view->stochasticChildren();
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
  if (!probe_design(*view, *S)) return fail_null("design probe");
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
    std::shared_ptr<Shared> shared = build_shared(declared, graph, true, view);
    // If the latents fail their checks, own the gates and slabs only.
    if (shared == nullptr && !declared.latents.empty()) shared = build_shared(declared, graph, false, view);
    if (shared == nullptr) return {};

    const CorrectedSamplerConfig config = corrected_sampler_configuration();
    shared->corrected = config.enabled;
    shared->settings = config.settings;
    shared->min_acceptance = config.min_acceptance;
    shared->mutant = test_mutant.load(std::memory_order_relaxed);
    for (ChildSpec &child : shared->children) {
      child.kind = child_kind(child.node, config.settings);
      if (child.kind != ChildKind::other) shared->any_surrogate = true;
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
    std::vector<Sampler *> out{sampler.get()};
    sampler.release();
    return out;
  }
};

} // namespace

LocationTag::LocationTag(int kind)
  : JointSamplerTag(kind == 0 ? group_tag_name : kind == 1 ? latent_tag_name : setup_tag_name,
                    kind == 0 ? 5 : kind == 1 ? 2 : 3),
    _kind(kind) {}

bool LocationTag::checkParameterLength(std::vector<unsigned long> const &len) const
{
  switch (_kind) {
  case 0:
    return len.size() == 5 && len[0] == 1 && len[1] >= 1 && len[2] == 1 && len[3] == len[1] && len[4] == 1;
  case 1:
    return len.size() == 2 && len[0] == 1 && len[1] == 1;
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
      "surrogate_evaluations", "exact_evaluations", "correction_evaluations",
      "max_correction", "max_correction_gaussian", "max_endpoint", "max_endpoint_change",
      "max_stage2_gaussian"};
    static std::vector<std::string> storage;
    for (int i = 0; i < component_count; ++i) {
      for (char const *stage : {"proposed", "stage0", "stage1", "accepted"}) {
        storage.push_back(std::string(component_names[i]) + "_" + stage);
      }
    }
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
    take(c.exact_evaluations), take(c.correction_evaluations), take_max(c.max_correction),
    take_max(c.max_correction_gaussian), take_max(c.max_endpoint), take_max(c.max_endpoint_change),
    take_max(c.max_stage2_gaussian)};
  for (int i = 0; i < component_count; ++i) {
    out.push_back(take(c.proposed[i]));
    out.push_back(take(c.stage0[i]));
    out.push_back(take(c.stage1[i]));
    out.push_back(take(c.accepted[i]));
  }
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

}
}
