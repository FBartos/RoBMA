#include "LocationScale.h"
#include "LocationScaleReport.h"
#include "GammaGroup.h"

#include <BayesTools/sampler_declarations.h>

#include <distribution/Distribution.h>
#include <distribution/ScalarDist.h>
#include <graph/AggNode.h>
#include <graph/DeterministicNode.h>
#include <graph/LogicalNode.h>
#include <graph/StochasticNode.h>
#include <model/Model.h>
#include <rng/RNG.h>
#include <sampler/GraphView.h>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <exception>
#include <limits>
#include <map>
#include <mutex>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace jags {
namespace RoBMA {
namespace {

using Vec = std::vector<double>;

const char *const tag_name = "robma_location_scale";
constexpr unsigned long tag_arity = 9;
constexpr double slice_width = 0.5;
constexpr int slice_budget = 40;
constexpr double probe_tolerance = 1e-9;
const double negative_infinity = -std::numeric_limits<double>::infinity();

// ---- counters, switches, decisions

struct Counters {
  std::atomic<std::uint64_t> owners{0}, declined{0}, moves{0}, identity{0}, skipped{0}, exhausted{0},
    evaluations{0};
  std::atomic<double> jump2{0.0}, seconds_embedded{0.0}, seconds_move{0.0};
};

Counters &counters()
{
  static Counters instance;
  return instance;
}

void count(std::atomic<std::uint64_t> &counter, std::uint64_t n = 1)
{
  counter.fetch_add(n, std::memory_order_relaxed);
}

void add(std::atomic<double> &target, double value)
{
  double old = target.load(std::memory_order_relaxed);
  while (!target.compare_exchange_weak(old, old + value, std::memory_order_relaxed)) {}
}

std::atomic<int> test_mutant(0);
std::atomic<int> test_decline(0);
std::atomic<int> test_timing(0);

struct DecisionStore {
  std::mutex mutex;
  std::vector<LocationScaleDecision> records;
};

DecisionStore &decision_store()
{
  static DecisionStore instance;
  return instance;
}

void record(LocationScaleDecision const &decision)
{
  DecisionStore &store = decision_store();
  const std::lock_guard<std::mutex> lock(store.mutex);
  store.records.push_back(decision);
}

// ---- the declared structure

enum class Form { none, gamma, share };

struct Spec {
  StochasticNode *a = nullptr;              // the source of a
  Node const *a_node = nullptr;             // the chart a
  Node const *sigma = nullptr;              // sigma_s
  Form form = Form::none;
  std::vector<StochasticNode *> weights;    // gamma: eta_1..eta_L; share: the share source
  std::vector<Node const *> w_nodes;        // gamma: the L ratios; share: the chart of w_s
  std::vector<StochasticNode *> indicators;
  std::vector<double> levels;
  StochasticNode *bias = nullptr;
  unsigned bias_levels = 0;
  unsigned L = 1, s = 0, K = 0;             // s is 0-based
  std::vector<StochasticNode *> latents;
  double a_lower = 0, a_upper = 0, share_lower = 0, share_upper = 1;
  std::unique_ptr<GraphView> view;          // the sources and latents: the move's writes
  std::vector<unsigned> offset;             // view offsets of the coordinates (latents, a, weights)
  std::vector<StochasticNode *> children;   // stochastic children of the view, markers excluded
  int mutant = 0;
  bool timing = false;
  std::string layout;

  unsigned coordinates() const { return K + 1 + static_cast<unsigned>(weights.size()); }

  std::vector<StochasticNode *> sources() const
  {
    std::vector<StochasticNode *> out(1, a);
    out.insert(out.end(), weights.begin(), weights.end());
    return out;
  }
};

// The move's inputs at one chain: the coordinates in the order latents, a,
// weight sources, with a, w_s, the gamma total T and sum(Z^2).
struct Point {
  Vec x;
  double a = 0, ws = 1, total = 0;
  long double z2 = 0;
};

// The prior log density (PDF_PRIOR) of a scalar source at x, as the node's own
// logDensity computes it from its fixed parameters and bounds.
double source_log_density(StochasticNode const *node, double x, unsigned int chain)
{
  auto const *dist = static_cast<ScalarDist const *>(node->distribution());
  const unsigned long npar = dist->npar();
  std::vector<double const *> parameters(npar);
  for (unsigned long i = 0; i < npar; ++i) parameters[i] = node->parents()[i]->value(chain);
  return dist->logDensity(x, PDF_PRIOR, parameters, node->lowerLimit(chain), node->upperLimit(chain));
}

void complete_point(Spec const &S, Point &p)
{
  long double z2 = 0;
  for (unsigned k = 0; k < S.K; ++k) z2 += static_cast<long double>(p.x[k]) * p.x[k];
  p.z2 = z2;
  p.a = p.x[S.K];
  p.ws = 1.0;
  p.total = 0.0;
  if (S.form == Form::gamma) {
    long double total = 0;
    for (unsigned j = 0; j < S.L; ++j) total += p.x[S.K + 1 + j];
    p.total = static_cast<double>(total);
    p.ws = p.x[S.K + 1 + S.s] / p.total;
  } else if (S.form == Form::share) {
    p.ws = p.x[S.K + 1];
  }
}

// T_c, c = exp(u), of the point (spec: The move); false when a transformed
// value is not strictly inside its node's support or not representable.
bool transform(Spec const &S, Point const &p, double u, int mutant, Vec &out)
{
  const double c = std::exp(u);
  const double c2 = std::exp(2.0 * u);
  const double d = 1.0 + p.ws * std::expm1(2.0 * u);
  if (!(c > 0.0) || !std::isfinite(c) || !(c2 > 0.0) || !std::isfinite(c2) || !(d > 0.0) || !std::isfinite(d)) {
    return false;
  }
  out.resize(p.x.size());
  for (unsigned k = 0; k < S.K; ++k) {
    out[k] = (mutant == 3 && k + 1 == S.K) ? p.x[k] : p.x[k] / c;
    if (!std::isfinite(out[k])) return false;
  }
  const double a1 = p.a * std::sqrt(d);
  if (!std::isfinite(a1) || !(a1 > S.a_lower && a1 < S.a_upper)) return false;
  out[S.K] = a1;
  if (S.form == Form::gamma) {
    for (unsigned j = 0; j < S.L; ++j) {
      const double eta = p.x[S.K + 1 + j];
      double value = j == S.s ? eta * c2 : eta;
      if (mutant != 4) value /= d;
      if (!std::isfinite(value) || !(value > 0.0)) return false;
      out[S.K + 1 + j] = value;
    }
  } else if (S.form == Form::share) {
    const double share = c2 * p.ws / d;
    const double rest = (1.0 - p.ws) / d;
    if (!std::isfinite(share) || !(share > S.share_lower && share < S.share_upper) || !(share < 1.0) ||
        !(rest > 0.0)) return false;
    out[S.K + 1] = share;
  }
  return true;
}

// log f(u) up to a constant (spec: The move); -Inf outside the support.
double log_f(Spec const &S, Point const &p, double u, unsigned int chain, int mutant, Vec &out)
{
  if (!transform(S, p, u, mutant, out)) return negative_infinity;
  long double lp = source_log_density(S.a, out[S.K], chain);
  for (std::size_t j = 0; j < S.weights.size(); ++j) lp += source_log_density(S.weights[j], out[S.K + 1 + j], chain);
  const long double zterm = p.z2 / static_cast<long double>(std::exp(2.0 * u));
  const double exponent = (mutant == 1 ? 1.0 : 2.0) - static_cast<double>(S.K);
  const double jacobian = mutant == 2 ? 0.0 : -(2.0 * S.L - 1.0) * std::log(out[S.K]);
  const double value = static_cast<double>(lp - 0.5L * zterm) + jacobian + exponent * u;
  return std::isfinite(value) ? value : negative_infinity;
}

// One slice update of t -> log f(t) from t = 0 (Neal 2003, Fig. 3): the slice
// level, an initial interval of width 0.5 placed at random around 0, a total
// budget of 40 expansions split at random between the sides (fixed_cap: 20
// per side, not a valid kernel; test hook only), then shrinkage. Returns the
// new point, or 0 with *exhausted when the interval is exhausted in floating
// point. Draws: one exponential, one uniform for the placement, one for the
// division and one per shrinkage step.
template <class LogF, class Uniform, class Exponential>
double slice_step(LogF &&log_density, double log_density0, Uniform &&uniform, Exponential &&exponential,
                  bool fixed_cap, bool *exhausted)
{
  *exhausted = false;
  const double logy = log_density0 - exponential();
  double left = -slice_width * uniform();
  double right = left + slice_width;
  int steps_left = static_cast<int>(std::floor(slice_budget * uniform()));
  int steps_right = slice_budget - 1 - steps_left;
  if (fixed_cap) steps_left = steps_right = slice_budget / 2;
  while (steps_left > 0 && log_density(left) > logy) {
    left -= slice_width;
    --steps_left;
  }
  while (steps_right > 0 && log_density(right) > logy) {
    right += slice_width;
    --steps_right;
  }
  for (;;) {
    const double middle = left + 0.5 * (right - left);
    if (!(middle > left && middle < right)) {
      *exhausted = true;
      return 0.0;
    }
    const double t = left + (right - left) * uniform();
    if (log_density(t) > logy) return t;
    if (t < 0.0) left = t;
    else right = t;
  }
}

// ---- the owner

class ScaleOwner;

// The live owners with their serial numbers (test hooks only).
struct OwnerRegistry {
  std::mutex mutex;
  std::vector<std::pair<int, ScaleOwner *>> owners;
  int next = 0;
};

OwnerRegistry &owner_registry()
{
  static OwnerRegistry instance;
  return instance;
}

// One sampler for the location members and the heterogeneity sources of all
// chains: the location block's own sampler, the sources' ordinary samplers
// (embedded unchanged) and the scale move.
class ScaleOwner : public Sampler {
  std::unique_ptr<Sampler> _location;
  std::vector<std::unique_ptr<Sampler>> _embedded;
  std::shared_ptr<Spec const> _spec;
  std::vector<Vec> _entry, _trial;
  std::vector<Point> _point;

public:
  ScaleOwner(GraphView *view, std::unique_ptr<Sampler> location, std::vector<std::unique_ptr<Sampler>> embedded,
             std::shared_ptr<Spec const> spec, unsigned int nchain)
    : Sampler(view), _location(std::move(location)), _embedded(std::move(embedded)), _spec(std::move(spec)),
      _entry(nchain), _trial(nchain), _point(nchain)
  {
    OwnerRegistry &registry = owner_registry();
    const std::lock_guard<std::mutex> lock(registry.mutex);
    registry.owners.emplace_back(registry.next++, this);
  }

  ~ScaleOwner() override
  {
    OwnerRegistry &registry = owner_registry();
    const std::lock_guard<std::mutex> lock(registry.mutex);
    registry.owners.erase(std::remove_if(registry.owners.begin(), registry.owners.end(),
                                         [this](std::pair<int, ScaleOwner *> const &entry) { return entry.second == this; }),
                          registry.owners.end());
  }

  Spec const &spec() const { return *_spec; }

  void update(unsigned int chain, RNG *rng) override
  {
    _location->update(chain, rng);
    if (_spec->timing) {
      const auto start = std::chrono::steady_clock::now();
      for (auto &sampler : _embedded) sampler->update(chain, rng);
      const auto middle = std::chrono::steady_clock::now();
      move(chain, rng);
      const auto end = std::chrono::steady_clock::now();
      add(counters().seconds_embedded, std::chrono::duration<double>(middle - start).count());
      add(counters().seconds_move, std::chrono::duration<double>(end - middle).count());
    } else {
      for (auto &sampler : _embedded) sampler->update(chain, rng);
      move(chain, rng);
    }
  }

  bool isAdaptive() const override
  {
    bool adaptive = _location->isAdaptive();
    for (auto const &sampler : _embedded) adaptive = sampler->isAdaptive() || adaptive;
    return adaptive;
  }

  void adaptOff() override
  {
    _location->adaptOff();
    for (auto &sampler : _embedded) sampler->adaptOff();
  }

  bool checkAdaptation() const override
  {
    bool done = _location->checkAdaptation();
    for (auto const &sampler : _embedded) done = sampler->checkAdaptation() && done;
    return done;
  }

  std::string name() const override { return _location->name(); }

  // whether every declared gate indicator is at its declared level
  bool active(unsigned int chain) const
  {
    Spec const &S = *_spec;
    for (std::size_t i = 0; i < S.indicators.size(); ++i) {
      if (S.indicators[i]->value(chain)[0] != S.levels[i]) return false;
    }
    return true;
  }

  void read(unsigned int chain, Point &p)
  {
    Spec const &S = *_spec;
    Vec &entry = _entry[chain];
    entry.resize(S.view->length());
    S.view->getValue(entry, chain);
    p.x.resize(S.coordinates());
    for (unsigned i = 0; i < S.coordinates(); ++i) p.x[i] = entry[S.offset[i]];
    complete_point(S, p);
  }

  // Writes the coordinates (with their deterministic descendants); restores
  // the entry values if the graph refuses them.
  void write(unsigned int chain, Vec const &x)
  {
    Spec const &S = *_spec;
    Vec entry(S.view->length());
    S.view->getValue(entry, chain);
    Vec values(entry);
    for (unsigned i = 0; i < S.coordinates(); ++i) values[S.offset[i]] = x[i];
    try {
      S.view->setValue(values, chain);
    } catch (...) {
      S.view->setValue(entry, chain);
      throw;
    }
  }

  double density(unsigned int chain, double u)
  {
    Point &p = _point[chain];
    read(chain, p);
    return log_f(*_spec, p, u, chain, _spec->mutant, _trial[chain]);
  }

  bool apply(unsigned int chain, double u)
  {
    Point &p = _point[chain];
    read(chain, p);
    if (!transform(*_spec, p, u, _spec->mutant, _trial[chain])) return false;
    write(chain, _trial[chain]);
    return true;
  }

private:
  void move(unsigned int chain, RNG *rng)
  {
    Spec const &S = *_spec;
    if (!active(chain)) {
      count(counters().identity);
      return;
    }
    Point &p = _point[chain];
    read(chain, p);
    Vec &trial = _trial[chain];
    std::uint64_t evaluations = 0;
    auto density_at = [&](double u) {
      ++evaluations;
      return log_f(S, p, u, chain, S.mutant, trial);
    };
    const double f0 = density_at(0.0);
    if (!std::isfinite(f0)) {
      count(counters().evaluations, evaluations);
      count(counters().skipped);
      return;
    }
    bool exhausted = false;
    const double u = slice_step(density_at, f0, [rng]() { return rng->uniform(); },
                                [rng]() { return rng->exponential(); }, false, &exhausted);
    count(counters().evaluations, evaluations);
    count(counters().moves);
    if (exhausted) {
      count(counters().exhausted);
      return;
    }
    if (u == 0.0) return;
    // an accepted point is inside the support, so the transform succeeds
    if (!transform(S, p, u, S.mutant, trial)) throw std::logic_error("location scale move outside the support");
    write(chain, trial);
    add(counters().jump2, u * u);
  }
};

ScaleOwner &live_owner(int serial)
{
  OwnerRegistry &registry = owner_registry();
  const std::lock_guard<std::mutex> lock(registry.mutex);
  for (auto const &entry : registry.owners) {
    if (entry.first == serial) return *entry.second;
  }
  throw std::runtime_error("No live location scale owner with this serial number.");
}

void check_chain(ScaleOwner &owner, unsigned int chain)
{
  if (chain >= owner.nodes()[0]->nchain()) throw std::runtime_error("The chain index is out of range.");
}

// ---- construction

// The deparse of a logical node with synthetic argument names x, y, ...
// equals one of the forms, with or without enclosing parentheses.
bool deparses_as(LogicalNode const *node, std::vector<std::string> const &forms)
{
  const std::size_t n = node->parents().size();
  std::vector<std::string> names;
  for (std::size_t i = 0; i < n; ++i) names.push_back(std::string(1, static_cast<char>('x' + i)));
  const std::string text = node->deparse(names);
  for (std::string const &form : forms) {
    if (text == form || text == "(" + form + ")") return true;
  }
  return false;
}

// The scalar nodes of a reference: the node itself, an aggregate of scalar
// nodes (offsets 0) or a concatenation c(...) of scalar nodes.
bool resolve_nodes(Node const *ref, std::vector<Node const *> &out)
{
  out.clear();
  if (auto const *aggregate = dynamic_cast<AggNode const *>(ref)) {
    for (std::size_t i = 0; i < aggregate->parents().size(); ++i) {
      if (aggregate->parents()[i]->length() != 1 || aggregate->offsets()[i] != 0) return false;
      out.push_back(aggregate->parents()[i]);
    }
  } else if (auto const *logical = dynamic_cast<LogicalNode const *>(ref)) {
    const std::size_t n = logical->parents().size();
    if (ref->length() == 1) {
      out.push_back(ref);
      return true;
    }
    std::string tight = "c(", spaced = "c(";
    for (std::size_t i = 0; i < n; ++i) {
      const std::string name(1, static_cast<char>('x' + i));
      tight += (i == 0 ? "" : ",") + name;
      spaced += (i == 0 ? "" : ", ") + name;
    }
    if (!deparses_as(logical, {tight + ")", spaced + ")"})) return false;
    for (Node const *parent : logical->parents()) {
      if (parent->length() != 1) return false;
      out.push_back(parent);
    }
  } else {
    if (ref->length() != 1) return false;
    out.push_back(ref);
  }
  return out.size() == ref->length();
}

bool resolve_stochastic(Node const *ref, std::vector<StochasticNode *> &out)
{
  std::vector<Node const *> nodes;
  out.clear();
  if (!resolve_nodes(ref, nodes)) return false;
  for (Node const *node : nodes) {
    auto *stochastic = dynamic_cast<StochasticNode *>(const_cast<Node *>(node));
    if (stochastic == nullptr) return false;
    out.push_back(stochastic);
  }
  return true;
}

bool unobserved(StochasticNode const *node)
{
  for (unsigned long i = 0; i < node->length(); ++i) if (node->isObserved(i)) return false;
  return true;
}

// A scalar continuous unobserved free node with fixed parameters and bounds
// and a fixed support (lower, upper) with lower < upper.
bool scalar_source(StochasticNode *node, std::set<StochasticNode *> const &available, double *lower, double *upper)
{
  if (node->length() != 1 || node->isDiscreteValued() || !available.count(node) || !unobserved(node) ||
      dynamic_cast<ScalarDist const *>(node->distribution()) == nullptr || !isSupportFixed(node)) return false;
  for (Node const *parent : node->parents()) if (!parent->isFixed()) return false;
  node->support(lower, upper, 1, 0);
  return *lower < *upper;
}

bool integer_vector(Node const *node, std::size_t length, std::vector<double> &out)
{
  if (!node->isFixed() || node->length() != length) return false;
  out.assign(node->value(0), node->value(0) + length);
  for (double v : out) if (!std::isfinite(v) || v != std::floor(v)) return false;
  return true;
}

bool close(double observed, double expected)
{
  return std::fabs(observed - expected) <= probe_tolerance * (1.0 + std::fabs(expected));
}

// The parent values of the children, concatenated.
void parent_values(std::vector<StochasticNode *> const &children, Vec &out)
{
  out.clear();
  for (StochasticNode const *child : children) {
    for (Node const *parent : child->parents()) {
      double const *value = parent->value(0);
      out.insert(out.end(), value, value + parent->length());
    }
  }
}

// The invariance check (spec: Ownership), on chain 0 with the state restored:
// with the indicators at their declared levels, the latents at nonzero values
// and every level of the bias indicator, the charts of a, w and sigma_s match
// the sources and, for c in {0.7, 1.3}, sigma_s' = c sigma_s and every parent
// value of every child is unchanged. Returns "passed" or "failed: <reason>".
std::string invariance_check(Spec const &S, Graph const &graph)
{
  std::vector<StochasticNode *> nodes = S.sources();
  nodes.insert(nodes.end(), S.latents.begin(), S.latents.end());
  nodes.insert(nodes.end(), S.indicators.begin(), S.indicators.end());
  if (S.bias != nullptr) nodes.push_back(S.bias);
  GraphView probe(nodes, graph, true);
  std::map<Node const *, unsigned> offset;
  {
    unsigned position = 0;
    for (StochasticNode *node : probe.nodes()) {
      offset[node] = position;
      position += static_cast<unsigned>(node->length());
    }
  }
  Vec saved(probe.length());
  probe.getValue(saved, 0);
  Vec base(saved);
  for (std::size_t i = 0; i < S.indicators.size(); ++i) base[offset[S.indicators[i]]] = S.levels[i];
  for (unsigned k = 0; k < S.K; ++k) {
    base[offset[S.latents[k]]] = (k % 2 == 0 ? 1.0 : -1.0) * (0.5 + 0.25 * (k % 4));
  }
  std::vector<double> bias_levels;
  if (S.bias == nullptr) bias_levels.push_back(0.0);
  else for (unsigned l = 1; l <= S.bias_levels; ++l) bias_levels.push_back(l);
  std::vector<StochasticNode *> coordinates = S.latents;
  coordinates.push_back(S.a);
  coordinates.insert(coordinates.end(), S.weights.begin(), S.weights.end());

  std::string outcome = "passed";
  try {
    for (double level : bias_levels) {
      Vec v(base);
      if (S.bias != nullptr) v[offset[S.bias]] = level;
      probe.setValue(v, 0);
      Point p;
      p.x.resize(coordinates.size());
      for (std::size_t i = 0; i < coordinates.size(); ++i) p.x[i] = v[offset[coordinates[i]]];
      complete_point(S, p);
      const double sigma0 = S.sigma->value(0)[0];
      if (!close(S.a_node->value(0)[0], p.a)) { outcome = "failed: the chart of a"; break; }
      if (S.form == Form::gamma) {
        for (unsigned j = 0; j < S.L; ++j) {
          if (!close(S.w_nodes[j]->value(0)[0], p.x[S.K + 1 + j] / p.total)) outcome = "failed: the weight charts";
        }
      } else if (S.form == Form::share) {
        if (!close(S.w_nodes[0]->value(0)[0], p.ws)) outcome = "failed: the share chart";
      }
      if (outcome != "passed") break;
      if (!(sigma0 > 0.0) || !close(sigma0, p.a * std::sqrt(p.ws))) { outcome = "failed: the chart of sigma_s"; break; }
      Vec before, after, moved;
      parent_values(S.children, before);
      for (double c : {0.7, 1.3}) {
        if (!transform(S, p, std::log(c), 0, moved)) { outcome = "failed: a probe point is outside the support"; break; }
        Vec w(v);
        for (std::size_t i = 0; i < coordinates.size(); ++i) w[offset[coordinates[i]]] = moved[i];
        probe.setValue(w, 0);
        if (!close(S.sigma->value(0)[0], c * sigma0)) { outcome = "failed: sigma_s does not scale by c"; break; }
        parent_values(S.children, after);
        if (after.size() != before.size()) { outcome = "failed: parent lengths"; break; }
        for (std::size_t i = 0; i < before.size(); ++i) {
          if (!close(after[i], before[i])) {
            outcome = "failed: a child's parent value changes";
            break;
          }
        }
        if (outcome != "passed") break;
        probe.setValue(v, 0);
      }
      if (outcome != "passed") break;
    }
  } catch (std::exception const &error) {
    outcome = std::string("failed: ") + error.what();
  }
  probe.setValue(saved, 0);
  return outcome;
}

std::string join(std::vector<std::string> const &parts)
{
  std::string out;
  for (std::size_t i = 0; i < parts.size(); ++i) out += (i ? "," : "") + parts[i];
  return out;
}

} // namespace

LocationScaleTag::LocationScaleTag() : JointSamplerTag(tag_name, tag_arity) {}

bool LocationScaleTag::checkParameterLength(std::vector<unsigned long> const &len) const
{
  return len.size() == tag_arity && len[0] == 1 && len[1] == 1 && len[2] == 1 && len[3] >= 1 && len[4] == 1 &&
    len[5] >= 1 && len[6] == len[5] && len[7] == 1 && len[8] == 3;
}

std::unique_ptr<Sampler> location_scale_owner(std::unique_ptr<Sampler> location,
                                              std::list<StochasticNode *> const &free, Graph const &graph,
                                              std::vector<StochasticNode *> const &latents,
                                              SamplerFactory const *location_factory)
{
  StochasticNode *anchor = nullptr;
  Node const *tag = nullptr;
  unsigned tags = 0;
  for (StochasticNode *node : free) {
    for (Node const *found : direct_terminal_tags(node, tag_name, tag_arity, graph)) {
      anchor = node;
      tag = found;
      ++tags;
    }
  }
  if (tags == 0) return location;

  auto spec = std::make_shared<Spec>();
  Spec &S = *spec;
  LocationScaleDecision decision;
  decision.check = "not run";
  decision.status = "declined";
  auto decline = [&](std::string const &reason) {
    decision.reason = reason;
    record(decision);
    count(counters().declined);
    return std::move(location);
  };
  if (tags != 1) return decline("declaration: more than one scale tag");
  auto const &par = tag->parents();
  std::vector<double> counts;
  if (!integer_vector(par[8], 3, counts)) return decline("declaration: counts");
  decision.L = static_cast<int>(counts[0]);
  decision.K = static_cast<int>(counts[2]);
  if (!(counts[0] >= 1 && counts[1] >= 1 && counts[1] <= counts[0] && counts[2] >= 1)) {
    return decline("declaration: counts");
  }
  S.L = static_cast<unsigned>(counts[0]);
  S.s = static_cast<unsigned>(counts[1]) - 1;
  S.K = static_cast<unsigned>(counts[2]);
  S.mutant = test_mutant.load(std::memory_order_relaxed);
  S.timing = test_timing.load(std::memory_order_relaxed) != 0;

  std::set<StochasticNode *> available(free.begin(), free.end());
  S.a = anchor;
  S.a_node = par[1];
  S.sigma = par[2];
  if (S.a_node->length() != 1 || S.sigma->length() != 1) return decline("declaration: charts");
  if (!scalar_source(S.a, available, &S.a_lower, &S.a_upper) || !(S.a_lower >= 0.0)) {
    return decline("the source of a: not a free continuous node with fixed parameters on (0, Inf)");
  }

  // the weights: gamma sources behind the declared ratios, a share, or none
  const bool no_weights = par[3]->isFixed();
  const bool no_share = par[4]->isFixed();
  decision.layout = no_weights ? "single" : !no_share ? "clustered" : "allocation";
  if (no_weights) {
    if (S.L != 1 || !no_share) return decline("declaration: weights");
    S.form = Form::none;
    S.layout = "single";
  } else if (!no_share) {
    std::vector<StochasticNode *> share;
    if (!resolve_stochastic(par[4], share) || share.size() != 1 || !resolve_nodes(par[3], S.w_nodes) ||
        S.w_nodes.size() != 1 || S.L != 2 || S.s != 0) return decline("declaration: share");
    if (!scalar_source(share[0], available, &S.share_lower, &S.share_upper) ||
        !(S.share_lower >= 0.0 && S.share_upper <= 1.0)) {
      return decline("the share: not a free continuous node with fixed parameters on (0, 1)");
    }
    S.form = Form::share;
    S.weights = share;
    S.layout = "clustered";
  } else {
    if (!resolve_nodes(par[3], S.w_nodes) || S.w_nodes.size() != S.L || S.L < 2) return decline("declaration: weights");
    Node const *denominator = nullptr;
    for (Node const *w : S.w_nodes) {
      auto const *ratio = dynamic_cast<LogicalNode const *>(w);
      if (ratio == nullptr || w->parents().size() != 2) return decline("the weights: not ratios");
      if (!deparses_as(ratio, {"x/y", "x / y"})) return decline("the weights: not ratios");
      auto *eta = dynamic_cast<StochasticNode *>(const_cast<Node *>(w->parents()[0]));
      if (eta == nullptr) return decline("the weights: not ratios of stochastic nodes");
      if (denominator != nullptr && w->parents()[1] != denominator) return decline("the weights: two denominators");
      denominator = w->parents()[1];
      S.weights.push_back(eta);
    }
    // the denominator is sum() of an aggregate of exactly the sources
    auto const *sum = dynamic_cast<LogicalNode const *>(denominator);
    if (sum == nullptr || denominator->parents().size() != 1 || !deparses_as(sum, {"sum(x)"})) {
      return decline("the weights: the denominator is not a sum");
    }
    std::vector<Node const *> parts;
    if (!resolve_nodes(denominator->parents()[0], parts) ||
        std::set<Node const *>(parts.begin(), parts.end()) != std::set<Node const *>(S.weights.begin(), S.weights.end()) ||
        parts.size() != S.L) return decline("the weights: the sum is not over the sources");
    std::vector<GammaGroup> groups = find_gamma_groups(S.weights);
    if (groups.size() != 1 || groups[0].members.size() != S.L) return decline("the weights: not a gamma group");
    for (StochasticNode *eta : S.weights) {
      double lower, upper;
      if (!scalar_source(eta, available, &lower, &upper)) return decline("the weights: a source is not free");
    }
    S.form = Form::gamma;
    S.layout = "allocation";
  }
  decision.layout = S.layout;

  // indicators, levels and the bias indicator
  if (!par[5]->isFixed()) {
    if (!resolve_stochastic(par[5], S.indicators)) return decline("declaration: indicators");
    if (!integer_vector(par[6], S.indicators.size(), S.levels)) return decline("declaration: levels");
    for (StochasticNode const *indicator : S.indicators) {
      if (indicator->length() != 1 || !indicator->isDiscreteValued()) return decline("declaration: indicators");
    }
  }
  if (!par[7]->isFixed()) {
    S.bias = dynamic_cast<StochasticNode *>(const_cast<Node *>(par[7]));
    const CategoricalWeights weights = categorical_weights(S.bias);
    if (S.bias == nullptr || !weights) return decline("declaration: bias indicator");
    S.bias_levels = static_cast<unsigned>(weights.size);
  }

  // the latents: every latent of the sampled term, owned by the block
  if (latents.size() != S.K) return decline("the latents: not all owned by the location block");
  S.latents = latents;
  std::vector<StochasticNode *> sources = S.sources();
  std::set<StochasticNode *> distinct(sources.begin(), sources.end());
  distinct.insert(latents.begin(), latents.end());
  if (distinct.size() != sources.size() + latents.size()) return decline("declaration: overlapping nodes");
  for (StochasticNode *indicator : S.indicators) {
    if (distinct.count(indicator)) return decline("declaration: overlapping nodes");
  }
  if (S.bias != nullptr && distinct.count(S.bias)) return decline("declaration: overlapping nodes");

  // ownership: no gate marker, and a block marker shared only among the sources
  std::set<unsigned int> blocks;
  for (StochasticNode *source : sources) {
    if (BayesTools::sampler_declarations::has_gate_marker(source)) return decline("a source has a gate marker");
    const unsigned int id = BayesTools::sampler_declarations::block_id(source);
    if (id != 0) blocks.insert(id);
  }
  if (!blocks.empty()) {
    for (StochasticNode *node : free) {
      if (distinct.count(node)) continue;
      if (blocks.count(BayesTools::sampler_declarations::block_id(node))) {
        return decline("a source shares its block with another node");
      }
    }
  }
  if (test_decline.load(std::memory_order_relaxed) != 0) return decline("test switch");

  try {
    std::vector<StochasticNode *> members(latents);
    members.insert(members.end(), sources.begin(), sources.end());
    S.view.reset(new GraphView(members, graph, true));
  } catch (std::exception const &) {
    return decline("the scale view");
  }
  {
    std::map<Node const *, unsigned> offset;
    unsigned position = 0;
    for (StochasticNode *node : S.view->nodes()) {
      offset[node] = position;
      position += static_cast<unsigned>(node->length());
    }
    for (StochasticNode *z : latents) S.offset.push_back(offset[z]);
    for (StochasticNode *source : sources) S.offset.push_back(offset[source]);
    for (StochasticNode *child : S.view->stochasticChildren()) {
      const std::string name = child->distribution()->name();
      if (distinct.count(child) || name == BayesTools::sampler_declarations::block_marker_name() ||
          name == BayesTools::sampler_declarations::gate_marker_name()) continue;
      S.children.push_back(child);
    }
  }

  decision.check = invariance_check(S, graph);
  if (decision.check != "passed") return decline("invariance check");

  // the sources' ordinary samplers: the factories after the location block's,
  // asked as the model asks them (each until it returns nothing)
  std::vector<std::unique_ptr<Sampler>> embedded;
  std::vector<std::string> names;
  try {
    std::list<StochasticNode *> remaining(sources.begin(), sources.end());
    std::list<SamplerFactory *> const &factories = Model::samplerFactories();
    auto it = std::find(factories.begin(), factories.end(), location_factory);
    if (it == factories.end()) return decline("the location factory is not registered");
    for (++it; it != factories.end() && !remaining.empty(); ++it) {
      if (!(*it)->isActive()) continue;
      for (;;) {
        std::vector<Sampler *> made = (*it)->makeSamplers(remaining, graph);
        if (made.empty()) break;
        bool inside = true;
        for (Sampler *sampler : made) {
          embedded.emplace_back(sampler);
          for (StochasticNode *node : sampler->nodes()) {
            if (std::find(remaining.begin(), remaining.end(), node) == remaining.end()) inside = false;
            remaining.remove(node);
          }
        }
        if (!inside) return decline("an embedded sampler claims another node");
        if (remaining.empty()) break;
      }
    }
    if (!remaining.empty()) return decline("no sampler for a source");
  } catch (std::exception const &error) {
    return decline(std::string("embedded samplers: ") + error.what());
  }
  for (auto const &sampler : embedded) names.push_back(sampler->name());

  const unsigned int nchain = location->nodes()[0]->nchain();
  std::vector<StochasticNode *> all(location->nodes());
  all.insert(all.end(), sources.begin(), sources.end());
  std::unique_ptr<GraphView> outer;
  try {
    outer.reset(new GraphView(all, graph, true));
  } catch (std::exception const &) {
    return decline("the owner view");
  }
  std::unique_ptr<Sampler> owner(new ScaleOwner(outer.release(), std::move(location), std::move(embedded),
                                                spec, nchain));
  decision.status = "owned";
  decision.embedded = join(names);
  record(decision);
  count(counters().owners);
  return owner;
}

std::vector<char const *> const &location_scale_report_names()
{
  static const std::vector<char const *> names = {
    "owners", "declined", "moves", "identity", "skipped", "exhausted", "evaluations", "jump2",
    "seconds_embedded", "seconds_move"};
  return names;
}

std::vector<double> location_scale_report(bool clear)
{
  Counters &c = counters();
  auto take = [clear](std::atomic<std::uint64_t> &counter) {
    return static_cast<double>(clear ? counter.exchange(0, std::memory_order_relaxed)
                                     : counter.load(std::memory_order_relaxed));
  };
  auto take_double = [clear](std::atomic<double> &value) {
    return clear ? value.exchange(0.0, std::memory_order_relaxed) : value.load(std::memory_order_relaxed);
  };
  return {take(c.owners), take(c.declined), take(c.moves), take(c.identity), take(c.skipped), take(c.exhausted),
          take(c.evaluations), take_double(c.jump2), take_double(c.seconds_embedded), take_double(c.seconds_move)};
}

std::vector<LocationScaleDecision> location_scale_decisions(bool clear)
{
  DecisionStore &store = decision_store();
  const std::lock_guard<std::mutex> lock(store.mutex);
  std::vector<LocationScaleDecision> out = store.records;
  if (clear) store.records.clear();
  return out;
}

std::vector<int> location_scale_test_config(int mutant, int decline, int timing)
{
  std::vector<int> previous = {test_mutant.load(), test_decline.load(), test_timing.load()};
  if (mutant >= 0) test_mutant.store(mutant);
  if (decline >= 0) test_decline.store(decline);
  if (timing >= 0) test_timing.store(timing);
  return previous;
}

int location_scale_test_last()
{
  OwnerRegistry &registry = owner_registry();
  const std::lock_guard<std::mutex> lock(registry.mutex);
  return registry.owners.empty() ? -1 : registry.owners.back().first;
}

LocationScaleState location_scale_test_state(int index, unsigned int chain)
{
  ScaleOwner &owner = live_owner(index);
  check_chain(owner, chain);
  Spec const &S = owner.spec();
  LocationScaleState state;
  state.layout = S.layout;
  state.L = static_cast<int>(S.L);
  state.s = static_cast<int>(S.s) + 1;
  state.K = static_cast<int>(S.K);
  state.active = owner.active(chain);
  Point p;
  owner.read(chain, p);
  state.values = p.x;
  if (S.form == Form::gamma) {
    for (unsigned j = 0; j < S.L; ++j) state.weights.push_back(p.x[S.K + 1 + j] / p.total);
  } else if (S.form == Form::share) {
    state.weights = {p.ws, 1.0 - p.ws};
  } else {
    state.weights = {1.0};
  }
  state.sigma = S.sigma->value(chain)[0];
  for (StochasticNode const *child : S.children) state.children.push_back(child->logDensity(chain, PDF_FULL));
  long double nodes = S.a->logDensity(chain, PDF_PRIOR);
  long double move = source_log_density(S.a, p.a, chain);
  for (std::size_t j = 0; j < S.weights.size(); ++j) {
    nodes += S.weights[j]->logDensity(chain, PDF_PRIOR);
    move += source_log_density(S.weights[j], p.x[S.K + 1 + j], chain);
  }
  state.prior_nodes = static_cast<double>(nodes);
  state.prior_move = static_cast<double>(move);
  return state;
}

void location_scale_test_set(int index, unsigned int chain, std::vector<double> const &values)
{
  ScaleOwner &owner = live_owner(index);
  check_chain(owner, chain);
  if (values.size() != owner.spec().coordinates()) throw std::runtime_error("The coordinates have the wrong length.");
  owner.write(chain, values);
}

bool location_scale_test_transform(int index, unsigned int chain, double u)
{
  ScaleOwner &owner = live_owner(index);
  check_chain(owner, chain);
  return owner.apply(chain, u);
}

std::vector<double> location_scale_test_density(int index, unsigned int chain, std::vector<double> const &u)
{
  ScaleOwner &owner = live_owner(index);
  check_chain(owner, chain);
  std::vector<double> out;
  for (double value : u) out.push_back(owner.density(chain, value));
  return out;
}

std::vector<double> location_scale_test_slice(double lower, double upper, std::vector<double> const &starts,
                                              bool fixed_cap, double (*uniform)(), double (*exponential)())
{
  std::vector<double> out;
  out.reserve(starts.size());
  for (double x0 : starts) {
    auto density = [&](double t) {
      const double x = x0 + t;
      return x > lower && x < upper ? 0.0 : negative_infinity;
    };
    bool exhausted = false;
    const double t = slice_step(density, 0.0, uniform, exponential, fixed_cap, &exhausted);
    out.push_back(x0 + t);
  }
  return out;
}

}
}
