#include "IndependentBiasJoint.h"
#include "IndependentBiasProposal.h"

#include "../distributions/selnorm-jags-bounds.h"

#include <graph/Graph.h>
#include <graph/StochasticNode.h>
#include <rng/RNG.h>
#include <sampler/GraphView.h>
#include <sampler/Sampler.h>
#include <sampler/SingletonGraphView.h>
#include <util/nainf.h>

#include <array>
#include <memory>
#include <numeric>
#include <set>
#include <stdexcept>
#include <utility>

namespace jags::RoBMA {

namespace {

constexpr char marker_name[] = "drobma_joint_bias";
constexpr char sampler_name[] = "RoBMA::IndependentBiasJoint";

bool fixed_scalar(Node const *node)
{
  return node->length() == 1 && node->isFixed();
}

bool observed_zero(StochasticNode const *node)
{
  return node->length() == 1 && isObserved(node) && node->value(0)[0] == 0;
}

std::vector<double> categorical_log_mass(StochasticNode const *node)
{
  const auto &parents = node->parents();
  if (node->length() != 1 || node->distribution()->name() != "dcat" ||
      parents.size() != 1 || !isSupportFixed(node) || !parents[0]->isFixed()) return {};
  const double *start = parents[0]->value(0);
  const double *end = start + parents[0]->length();
  const double sum = std::accumulate(start, end, 0.0);
  if (!std::isfinite(sum) || !(sum > 0)) return {};
  const double log_sum = std::log(sum);
  std::vector<double> out;
  out.reserve(parents[0]->length());
  for (; start != end; ++start) {
    if (!std::isfinite(*start) || !(*start > 0)) return {};
    out.push_back(std::log(*start) - log_sum);
  }
  return out;
}

bool normal_prior(StochasticNode const *node, double &mean, double &precision)
{
  const auto &parents = node->parents();
  if (node->length() != 1 || node->distribution()->name() != "dnorm" ||
      parents.size() != 2 || !isSupportFixed(node) ||
      !fixed_scalar(parents[0]) || !fixed_scalar(parents[1])) return false;
  double lower, upper;
  node->support(&lower, &upper, 1, 0);
  if (lower != JAGS_NEGINF || upper != JAGS_POSINF) return false;
  mean = parents[0]->value(0)[0];
  precision = parents[1]->value(0)[0];
  return std::isfinite(mean) && std::isfinite(precision) && precision > 0;
}

struct GeometryNodes {
  Node const *slopes = nullptr;
  Node const *omega = nullptr;
  std::vector<Node const *> total_sd;
  int slab_label = 0;
  int location_sign = 0;
};

// Stop at stochastic inputs. Their values stay fixed while this sampler
// changes other roots, even if their density depends on one of those roots.
std::set<Node const *> stochastic_inputs(Node const *node)
{
  std::vector<Node const *> pending = {node};
  std::set<Node const *> seen;
  std::set<Node const *> inputs;
  while (!pending.empty()) {
    const Node *current = pending.back();
    pending.pop_back();
    if (!seen.insert(current).second) continue;
    if (current->isStochastic()) {
      inputs.insert(current);
      continue;
    }
    const auto &parents = current->parents();
    pending.insert(pending.end(), parents.begin(), parents.end());
  }
  return inputs;
}

bool nested_weight_child(StochasticNode *node, const GraphView &view, const Graph &graph,
                         Node const *omega)
{
  const auto &parents = node->parents();
  if (node->distribution()->name() != "dbeta" || node->length() != 1 || isObserved(node) ||
      isBounded(node) || !isSupportFixed(node) || parents.size() != 2 ||
      stochastic_inputs(omega).count(node) == 0) return false;
  const auto *bias = view.nodes().front();
  bool bias_dependent = false;
  for (const Node *parent : parents) {
    if (parent->length() != 1) return false;
    for (const Node *input : stochastic_inputs(parent)) {
      if (input != bias) return false;
      bias_dependent = true;
    }
  }
  if (!bias_dependent) return false;

  // Downstream stochastic terms of the chart must already participate in the
  // joint view, where they are separately validated as rows or markers.
  SingletonGraphView chart(node, graph);
  const auto &known = view.stochasticChildren();
  for (const auto *child : chart.stochasticChildren()) {
    if (std::find(known.begin(), known.end(), child) == known.end()) return false;
  }
  return true;
}

bool add_likelihood_row(StochasticNode const *child, const GraphView &view,
                        joint_bias::Geometry &geometry, GeometryNodes &nodes)
{
  const auto &parents = child->parents();
  if (child->distribution()->name() != "dselnorm_step_switch" || !isObserved(child) ||
      child->length() != 1 || parents.size() != 11 || !view.isDependent(parents[0]) ||
      view.isDependent(parents[1]) || !fixed_scalar(parents[2]) ||
      !fixed_scalar(parents[3]) || !parents[5]->isFixed() || !parents[6]->isFixed() ||
      !fixed_scalar(parents[7]) || !fixed_scalar(parents[8]) ||
      parents[3]->value(0)[0] != 1) return false;

  const auto n_bins = parents[5]->length();
  if (n_bins < 2 || n_bins > static_cast<unsigned long>(std::numeric_limits<int>::max()) ||
      parents[6]->length() != n_bins ||
      !selnorm_jags_integer_in_range(parents[7]->value(0)[0], 0, static_cast<int>(n_bins))) return false;

  const double *z_lower = parents[5]->value(0);
  const double *z_upper = parents[6]->value(0);
  const double first = selnorm_jags_bound_to_ieee(z_upper[0]);
  const double last = selnorm_jags_bound_to_ieee(z_lower[n_bins - 1]);
  if (first != JAGS_POSINF || last != JAGS_NEGINF) return false;
  if (geometry.cuts.empty()) {
    geometry.cuts.resize(n_bins + 1);
    geometry.cuts[0] = first;
    for (unsigned long i = 0; i < n_bins; ++i) {
      const double lower = selnorm_jags_bound_to_ieee(z_lower[i]);
      if (!(lower < geometry.cuts[i]) ||
          selnorm_jags_bound_to_ieee(z_upper[i]) != geometry.cuts[i]) return false;
      geometry.cuts[i + 1] = lower;
    }
  } else {
    if (geometry.cuts.size() != n_bins + 1) return false;
    for (unsigned long i = 0; i < n_bins; ++i) {
      if (selnorm_jags_bound_to_ieee(z_upper[i]) != geometry.cuts[i] ||
          selnorm_jags_bound_to_ieee(z_lower[i]) != geometry.cuts[i + 1]) return false;
    }
  }

  const double outcome = child->value(0)[0];
  const double se = parents[2]->value(0)[0];
  const double sign = parents[8]->value(0)[0];
  const int bin = static_cast<int>(parents[7]->value(0)[0]);
  if (!std::isfinite(outcome) || !std::isfinite(se) || !(se > 0) ||
      (sign != 1 && sign != -1)) return false;
  geometry.rows.push_back({outcome, se, bin == 0 ? 0 : se, bin, static_cast<int>(sign)});
  nodes.total_sd.push_back(parents[1]);
  return true;
}

bool configure(StochasticNode const *marker, const GraphView &view, const Graph &graph,
               joint_bias::Geometry &geometry, GeometryNodes &nodes)
{
  const auto &par = marker->parents();
  const auto &owned = view.nodes();
  if (!fixed_scalar(par[3]) || !fixed_scalar(par[4]) || !par[5]->isFixed() ||
      view.isDependent(par[6]) || view.isDependent(par[7])) return false;
  const double slab = par[3]->value(0)[0];
  const double sign = par[4]->value(0)[0];
  if ((slab != 1 && slab != 2) || (sign != 1 && sign != -1)) return false;
  nodes.slab_label = static_cast<int>(slab);
  nodes.location_sign = static_cast<int>(sign);
  nodes.slopes = par[6];
  nodes.omega = par[7];

  const auto bias_mass = categorical_log_mass(owned[0]);
  const auto effect_mass = categorical_log_mass(owned[1]);
  const std::size_t n_bias = bias_mass.size();
  if (n_bias < 2 || n_bias > static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
      effect_mass.size() != 2 || par[5]->length() != n_bias ||
      nodes.slopes->length() != n_bias ||
      !normal_prior(owned[2], geometry.prior_mean, geometry.prior_precision)) return false;
  geometry.prior_mean *= nodes.location_sign;
  geometry.prior_sd = 1 / std::sqrt(geometry.prior_precision);

  bool has_selection = false;
  const double *family = par[5]->value(0);
  for (std::size_t i = 0; i < n_bias; ++i) {
    if (!selnorm_jags_integer_in_range(family[i], 0, 3)) return false;
    const auto code = static_cast<joint_bias::Family>(static_cast<int>(family[i]));
    geometry.families.push_back(code);
    has_selection = has_selection || code == joint_bias::Family::selection;
  }
  if (!has_selection) return false;
  for (int label : {3 - nodes.slab_label, nodes.slab_label}) {
    for (double log_mass : bias_mass) {
      const double log_cell_mass = log_mass + effect_mass[label - 1];
      if (!(std::exp(log_cell_mass) > 0)) return false;
      geometry.log_prior.push_back(log_cell_mass);
    }
  }

  bool nested_child = false;
  for (auto *child : view.stochasticChildren()) {
    if (child == marker) continue;
    // Nuisance branches may keep BayesTools' gated sampler. A marker naming
    // one of our own roots would declare an incompatible second owner.
    if (child->distribution()->name() == "dbt_gate" && observed_zero(child) &&
        child->parents().size() == 2 &&
        std::find(owned.begin(), owned.end(), child->parents().front()) == owned.end()) continue;
    // The shared weight is held fixed. GraphView evaluates its normalized
    // Beta density at each proposed bias label as a child likelihood term.
    if (!nested_child && nested_weight_child(child, view, graph, nodes.omega)) {
      nested_child = true;
      continue;
    }
    if (!add_likelihood_row(child, view, geometry, nodes)) return false;
  }
  return !geometry.rows.empty() && nodes.omega->length() % n_bias == 0 &&
    nodes.omega->length() / n_bias == geometry.cuts.size() - 1;
}

class IndependentBiasJoint : public Sampler {
  struct Chain {
    std::vector<double> total_sd;
    joint_bias::Proposal proposal;
    joint_bias::Scratch scratch;

    explicit Chain(const joint_bias::Geometry &geometry)
      : total_sd(geometry.rows.size()), proposal(geometry.log_prior.size()), scratch(geometry) {}
  };

  GraphView const *_view;
  joint_bias::Geometry _geometry;
  GeometryNodes _nodes;
  std::vector<Chain> _chains;

public:
  IndependentBiasJoint(std::unique_ptr<GraphView> view, joint_bias::Geometry geometry,
                       GeometryNodes nodes)
    : Sampler(view.get()), _view(view.release()), _geometry(std::move(geometry)),
      _nodes(std::move(nodes))
  {
    const unsigned int chains = _view->nodes()[0]->nchain();
    _chains.reserve(chains);
    for (unsigned int i = 0; i < chains; ++i) _chains.emplace_back(_geometry);
  }

  void update(unsigned int chain, RNG *rng) override
  {
    auto &work = _chains[chain];
    for (std::size_t i = 0; i < _nodes.total_sd.size(); ++i) {
      work.total_sd[i] = _nodes.total_sd[i]->value(chain)[0];
    }
    joint_bias::propose(_geometry, work.total_sd, _nodes.slopes->value(chain),
                        _nodes.omega->value(chain), work.proposal, work.scratch);
    const auto &q = work.proposal;

    const auto &roots = _view->nodes();
    const std::array<double, 3> previous = {
      roots[0]->value(chain)[0], roots[1]->value(chain)[0], roots[2]->value(chain)[0]
    };
    const auto n_bias = _geometry.families.size();
    if (!selnorm_jags_integer_in_range(previous[0], 1, static_cast<int>(n_bias)) ||
        (previous[1] != 1 && previous[1] != 2) || !std::isfinite(previous[2])) {
      throw std::runtime_error("Joint bias sampler encountered an invalid current state.");
    }
    const auto previous_cell = static_cast<std::size_t>(previous[0] - 1) +
      (previous[1] == _nodes.slab_label ? n_bias : 0);
    const bool use_prior = !q.local || rng->uniform() < joint_bias::prior_mixture;
    const std::size_t proposed_cell = q.cell(_geometry, use_prior, rng->uniform());
    const double mean = use_prior ? _geometry.prior_mean : q.mean[proposed_cell];
    const double sd = use_prior ? _geometry.prior_sd : q.sd[proposed_cell];
    const double value = mean + sd * rng->normal();
    if (!std::isfinite(value)) return;
    const std::array<double, 3> proposed = {
      static_cast<double>(proposed_cell % n_bias + 1),
      static_cast<double>(proposed_cell < n_bias ? 3 - _nodes.slab_label : _nodes.slab_label),
      _nodes.location_sign * value
    };

    const double previous_log_target = _view->logFullConditional(chain);
    const double previous_log_q = q.log_density(_geometry, previous_cell,
                                               _nodes.location_sign * previous[2]);
    if (!std::isfinite(previous_log_target) || !std::isfinite(previous_log_q)) {
      throw std::runtime_error("Joint bias sampler could not evaluate the current state.");
    }
    try {
      _view->setValue(proposed.data(), proposed.size(), chain);
      const double proposed_log_target = _view->logFullConditional(chain);
      if (proposed_log_target != JAGS_NEGINF) {
        const double proposed_log_q = q.log_density(_geometry, proposed_cell, value);
        const double ratio = (proposed_log_target - previous_log_target) +
          (previous_log_q - proposed_log_q);
        if (!std::isfinite(proposed_log_target) || !std::isfinite(proposed_log_q) || std::isnan(ratio)) {
          throw std::runtime_error("Joint bias sampler could not evaluate its acceptance probability.");
        }
        if (std::log(rng->uniform()) < ratio) return;
      }
    } catch (...) {
      _view->setValue(previous.data(), previous.size(), chain);
      throw;
    }
    _view->setValue(previous.data(), previous.size(), chain);
  }

  void adaptOff() override {}
  bool checkAdaptation() const override { return true; }
  bool isAdaptive() const override { return false; }
  std::string name() const override { return sampler_name; }
};

StochasticNode const *declared_marker(StochasticNode *node, const Graph &graph)
{
  StochasticNode const *marker = nullptr;
  for (StochasticNode const *child : *node->stochasticChildren()) {
    if (child->distribution()->name() != marker_name || !graph.contains(child) ||
        !observed_zero(child) || child->parents().size() != 8 || child->parents().front() != node) continue;
    if (marker) return nullptr;
    marker = child;
  }
  return marker;
}

} // namespace

DIndependentBiasJointMarker::DIndependentBiasJointMarker() : VectorDist(marker_name, 8) {}

double DIndependentBiasJointMarker::logDensity(
    double const *x, PDFType, std::vector<double const *> const &,
    std::vector<unsigned long> const &) const
{
  return *x == 0 ? 0 : JAGS_NEGINF;
}

void DIndependentBiasJointMarker::randomSample(double *x, std::vector<double const *> const &,
    std::vector<unsigned long> const &, RNG *) const
{
  *x = 0;
}

void DIndependentBiasJointMarker::support(double *lower, double *upper,
    std::vector<double const *> const &, std::vector<unsigned long> const &) const
{
  *lower = *upper = 0;
}

bool DIndependentBiasJointMarker::checkParameterLength(std::vector<unsigned long> const &len) const
{
  return len.size() == 8 && std::all_of(len.begin(), len.begin() + 5,
    [](unsigned long n) { return n == 1; });
}

bool DIndependentBiasJointMarker::checkParameterValue(std::vector<double const *> const &,
    std::vector<unsigned long> const &) const { return true; }
bool DIndependentBiasJointMarker::checkParameterDiscrete(std::vector<bool> const &) const { return true; }
bool DIndependentBiasJointMarker::isDiscreteValued(std::vector<bool> const &) const { return true; }
bool DIndependentBiasJointMarker::isSupportFixed(std::vector<bool> const &) const { return true; }
unsigned long DIndependentBiasJointMarker::length(std::vector<unsigned long> const &) const { return 1; }

std::vector<Sampler *> IndependentBiasJointFactory::makeSamplers(
    std::list<StochasticNode *> const &nodes, Graph const &graph) const
{
  for (auto *bias : nodes) {
    if (bias->distribution()->name() != "dcat") continue;
    auto *marker = declared_marker(bias, graph);
    if (!marker) continue;

    const auto &parents = marker->parents();
    std::vector<StochasticNode *> roots = {bias};
    for (unsigned int i = 1; i < 3; ++i) {
      auto entry = std::find(nodes.begin(), nodes.end(), parents[i]);
      if (entry == nodes.end() || std::find(roots.begin(), roots.end(), *entry) != roots.end()) break;
      roots.push_back(*entry);
    }
    if (roots.size() != 3 || roots[1]->distribution()->name() != "dcat" ||
        roots[2]->distribution()->name() != "dnorm") continue;

    // A multilevel view can classify a mismatched declaration without throwing;
    // configure() admits only roots with independent fixed-parameter priors.
    auto view = std::make_unique<GraphView>(roots, graph, true);
    joint_bias::Geometry geometry;
    GeometryNodes geometry_nodes;
    if (!configure(marker, *view, graph, geometry, geometry_nodes)) continue;
    auto sampler = std::make_unique<IndependentBiasJoint>(std::move(view), std::move(geometry),
                                                         std::move(geometry_nodes));
    std::vector<Sampler *> out = {sampler.get()};
    sampler.release();
    // JAGS calls this factory again with the claimed roots removed.
    return out;
  }
  return {};
}

std::string IndependentBiasJointFactory::name() const { return sampler_name; }

} // namespace jags::RoBMA
