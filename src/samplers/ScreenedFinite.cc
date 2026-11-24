#include "ScreenedFinite.h"
#include "CorrectedSampler.h"

#include <distribution/Distribution.h>
#include <graph/DeterministicNode.h>
#include <graph/LogicalNode.h>
#include <graph/StochasticNode.h>
#include <module/ModuleError.h>
#include <rng/RNG.h>
#include <sampler/MutableSampleMethod.h>
#include <sampler/MutableSampler.h>
#include <sampler/Sampler.h>
#include <sampler/SingletonGraphView.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <vector>

namespace jags {
namespace RoBMA {
namespace {

const char *const screened_finite_name = "RoBMA::ScreenedFinite";
const char *const screened_label = "Corrected finite sampler";
const double negative_infinity = -std::numeric_limits<double>::infinity();
// Larger supports keep the stock sampler.
const double maximum_levels = 10000;

}

ScreenedFiniteKernel::ScreenedFiniteKernel(StochasticNode *node, std::vector<StochasticNode *> const &children,
                                           unsigned int chain, SelNormCoarseSettings const &settings,
                                           double lower, unsigned int levels, int proposal)
  : _node(node), _children(children), _chain(chain), _settings(settings), _lower(lower), _levels(levels),
    _proposal(proposal)
{
  for (StochasticNode *child : _children) {
    const ChildKind kind = child_kind(child, settings);
    if (kind == ChildKind::other) {
      _child_selected.push_back(-1);
    } else {
      _child_selected.push_back(static_cast<int>(_selected.size()));
      _selected.emplace_back(child, kind);
    }
  }
  if (_selected.empty()) throw std::logic_error("A screened finite node has no selection child.");
}

double ScreenedFiniteKernel::prior_term() const
{
  return _node->logDensity(_chain, PDF_PRIOR);
}

// The children's terms of one level (the node is set): the surrogate of each
// coarse-capable child and the exact density of the others.
double ScreenedFiniteKernel::surrogate_score(bool *not_exact) const
{
  double out = prior_term();
  for (std::size_t i = 0; i < _children.size(); ++i) {
    const int selected = _child_selected[i];
    out += selected < 0 ? _children[i]->logDensity(_chain, PDF_LIKELIHOOD) :
      _selected[selected].logDensity(_chain, _settings, not_exact);
  }
  return out;
}

double ScreenedFiniteKernel::exact_score() const
{
  corrected_kind_counters(corrected_kind_finite).exact_scores.fetch_add(1, std::memory_order_relaxed);
  double out = prior_term();
  for (StochasticNode const *child : _children) out += child->logDensity(_chain, PDF_LIKELIHOOD);
  return out;
}

ScreenedFiniteKernel::Outcome ScreenedFiniteKernel::update(RNG *rng,
                                                           std::function<void(double)> const &set_value)
{
  Outcome outcome;
  const double current_value = _node->value(_chain)[0];
  const long current_index = static_cast<long>(current_value - _lower);
  if (current_index < 0 || current_index >= static_cast<long>(_levels)) {
    throwNodeError(_node, std::string(screened_label) + ": current value is outside the support");
  }
  // 1. the surrogate scores of every level with prior mass
  std::vector<double> score(_levels, negative_infinity);
  std::vector<char> inexact(_levels, 0);
  double maximum = negative_infinity;
  for (unsigned int l = 0; l < _levels; ++l) {
    set_value(_lower + l);
    if (prior_term() == negative_infinity) continue;
    bool not_exact = false;
    score[l] = surrogate_score(&not_exact);
    inexact[l] = not_exact;
    if (std::isnan(score[l]) || score[l] == std::numeric_limits<double>::infinity()) {
      throwNodeError(_node, std::string(screened_label) + ": a surrogate score is not finite or -Inf");
    }
    maximum = std::max(maximum, score[l]);
  }
  const unsigned int k = static_cast<unsigned int>(current_index);
  if (!std::isfinite(score[k])) {
    set_value(current_value);
    throwNodeError(_node, std::string(screened_label) + " requires a finite current surrogate score");
  }
  // 2. the proposal from the surrogate probabilities
  std::vector<double> mass(_levels, 0.0);
  for (unsigned int l = 0; l < _levels; ++l) mass[l] = std::exp(score[l] - maximum);
  unsigned int proposed = k;
  if (_proposal == 2) {
    // Liu's Metropolized Gibbs: a level other than k with probability
    // p_l / (1 - p_k), accepted on the surrogate with min(1, (1 - p_k) / (1 - p_l))
    double others = 0.0, total = 0.0;
    for (unsigned int l = 0; l < _levels; ++l) {
      total += mass[l];
      if (l != k) others += mass[l];
    }
    if (others > 0.0) {
      double draw = rng->uniform() * others;
      unsigned int l = 0;
      for (; l < _levels; ++l) {
        if (l == k || mass[l] == 0.0) continue;
        if (draw < mass[l]) break;
        draw -= mass[l];
      }
      if (l == _levels) {
        // rounding at the end of the cumulative sum: the last level with mass
        l = _levels;
        while (l-- > 0) if (l != k && mass[l] > 0.0) break;
      }
      double others_at_l = 0.0;
      for (unsigned int j = 0; j < _levels; ++j) if (j != l) others_at_l += mass[j];
      // (1 - p_k) / (1 - p_l) = others / others_at_l (the common total cancels)
      const double ratio = others / others_at_l;
      if (ratio >= 1.0 || rng->uniform() < ratio) proposed = l;
    }
  } else {
    std::vector<double> cumulative(_levels);
    double total = 0.0;
    for (unsigned int l = 0; l < _levels; ++l) {
      total += mass[l];
      cumulative[l] = total;
    }
    const double draw = rng->uniform() * total;
    proposed = static_cast<unsigned int>(
      std::upper_bound(cumulative.begin(), cumulative.end(), draw) - cumulative.begin());
    if (proposed >= _levels) proposed = _levels - 1;
    // a rounding draw of a level without mass goes to the last level with mass
    while (mass[proposed] == 0.0 && proposed > 0) --proposed;
  }
  // 3. a self-proposal ends the update without an exact evaluation
  if (proposed == k) {
    set_value(current_value);
    return outcome;
  }
  outcome.proposed = true;
  outcome.evidence = inexact[k] || inexact[proposed];
  const double proposed_value = _lower + proposed;
  if (_proposal == 1) {
    // test mutant: the surrogate draw without its correction
    set_value(proposed_value);
    outcome.accepted = true;
    return outcome;
  }
  // 4. the exact correction (the current level first, so that the last exact
  // evaluation is the proposed state)
  set_value(current_value);
  const double exact_current = exact_score();
  if (!std::isfinite(exact_current)) {
    throwNodeError(_node, std::string(screened_label) + " requires a finite current full density");
  }
  set_value(proposed_value);
  const double exact_proposed = exact_score();
  if (std::isnan(exact_proposed) || exact_proposed == std::numeric_limits<double>::infinity()) {
    set_value(current_value);
    throwNodeError(_node, std::string(screened_label) + " proposed a non-finite full density");
  }
  bool accepted = false;
  if (exact_proposed != negative_infinity) {
    const long double log_ratio =
      (static_cast<long double>(exact_proposed) - score[proposed]) -
      (static_cast<long double>(exact_current) - score[k]);
    if (std::isnan(log_ratio)) {
      set_value(current_value);
      throwNodeError(_node, std::string(screened_label) + " proposed an undefined acceptance ratio");
    }
    accepted = !(log_ratio < 0 && std::log(static_cast<long double>(rng->uniform())) > log_ratio);
  }
  if (!accepted) set_value(current_value);
  outcome.accepted = accepted;
  return outcome;
}

// The exact enumeration of JAGS 5's base::Finite: every level's full
// conditional, the same arithmetic order and one uniform draw.
void ScreenedFiniteKernel::enumerate(RNG *rng, std::function<void(double)> const &set_value)
{
  std::vector<double> likelihood(_levels);
  double maximum = negative_infinity;
  for (unsigned int i = 0; i < _levels; ++i) {
    set_value(_lower + i);
    likelihood[i] = exact_score();
    if (likelihood[i] > maximum) maximum = likelihood[i];
  }
  double total = 0;
  for (unsigned int i = 0; i < _levels; ++i) {
    total += std::exp(likelihood[i] - maximum);
    likelihood[i] = total;
  }
  if (!std::isfinite(total)) throwNodeError(_node, "Cannot normalize density");
  const double draw = rng->uniform() * total;
  const long index = std::upper_bound(likelihood.begin(), likelihood.end(), draw) - likelihood.begin();
  set_value(_lower + index);
}

void ScreenedFiniteKernel::scores(std::function<void(double)> const &set_value, double restore,
                                  std::vector<double> *surrogate, std::vector<double> *exact)
{
  surrogate->assign(_levels, negative_infinity);
  exact->assign(_levels, negative_infinity);
  for (unsigned int l = 0; l < _levels; ++l) {
    set_value(_lower + l);
    if (prior_term() == negative_infinity) continue;
    bool not_exact = false;
    (*surrogate)[l] = surrogate_score(&not_exact);
    (*exact)[l] = exact_score();
  }
  set_value(restore);
}

namespace {

class ScreenedFinite;

// The live screened samplers, for the test hook that reads their scores. Never
// destroyed, like the sampler registry: samplers may be destroyed late.
struct LiveScreened {
  std::mutex mutex;
  std::vector<ScreenedFinite *> live;
};

LiveScreened &live_screened()
{
  static LiveScreened *screened = new LiveScreened();
  return *screened;
}

class ScreenedFinite : public MutableSampleMethod, public CorrectedSamplerWindow {
  SingletonGraphView const *_graph;
  const unsigned int _chain;
  const double _min_acceptance;
  ScreenedFiniteKernel _kernel;
  std::function<void(double)> _set_value;
  bool _exact;
  CorrectedTally _tally;

  void leave_live()
  {
    LiveScreened &screened = live_screened();
    const std::lock_guard<std::mutex> lock(screened.mutex);
    std::vector<ScreenedFinite *> &live = screened.live;
    for (std::size_t i = 0; i < live.size(); ++i) {
      if (live[i] == this) {
        live.erase(live.begin() + static_cast<std::ptrdiff_t>(i));
        break;
      }
    }
  }

public:
  ScreenedFinite(SingletonGraphView const *graph, unsigned int chain, SelNormCoarseSettings const &settings,
                 double min_acceptance, double lower, unsigned int levels, int proposal)
    : _graph(graph), _chain(chain), _min_acceptance(min_acceptance),
      _kernel(graph->node(), graph->stochasticChildren(), chain, settings, lower, levels, proposal),
      _set_value([graph, chain](double value) { graph->setValue(&value, 1, chain); }),
      _exact(false), _tally(corrected_kind_finite)
  {
    graph->checkFinite(chain);
    {
      LiveScreened &screened = live_screened();
      const std::lock_guard<std::mutex> lock(screened.mutex);
      screened.live.push_back(this);
    }
    try {
      register_corrected_sampler(this, corrected_kind_finite);
    } catch (...) {
      leave_live();
      throw;
    }
  }

  ~ScreenedFinite() override
  {
    unregister_corrected_sampler(this);
    leave_live();
  }

  std::uint64_t windowProposed() const override { return _tally.window_proposed(); }
  std::uint64_t windowAccepted() const override { return _tally.window_accepted(); }
  void clearWindow() override { _tally.clear_window(); }
  double rateExponent() const override { return 1.0; }

  void update(RNG *rng) override
  {
    corrected_kind_counters(corrected_kind_finite).updates.fetch_add(1, std::memory_order_relaxed);
    // test hook: the indicator keeps its value
    if (corrected_test_freeze_finite()) return;
    const double entry = _graph->node()->value(_chain)[0];
    try {
      if (_exact) {
        _kernel.enumerate(rng, _set_value);
        return;
      }
      const ScreenedFiniteKernel::Outcome outcome = _kernel.update(rng, _set_value);
      if (outcome.proposed) _tally.count(outcome.accepted, outcome.evidence);
    }
    catch (...) {
      _graph->setValue(&entry, 1, _chain);
      throw;
    }
  }

  bool isAdaptive() const override
  {
    return true;
  }

  // The kernel is chosen once, at the first adaptOff(), by the evidence rule.
  void adaptOff() override
  {
    if (_tally.adapt_off(_min_acceptance, 1.0)) _exact = true;
  }

  bool checkAdaptation() const override
  {
    return true;
  }

  void append_scores(std::vector<double> &out)
  {
    std::vector<double> surrogate, exact;
    _kernel.scores(_set_value, _graph->node()->value(_chain)[0], &surrogate, &exact);
    out.push_back(_chain);
    out.push_back(static_cast<double>(surrogate.size()));
    out.insert(out.end(), surrogate.begin(), surrogate.end());
    out.insert(out.end(), exact.begin(), exact.end());
  }
};

// A finite readout tag (BayesTools' bt_finite_tag or bt_finite_witness_tag) on
// the node: the readout keeps the enumerating sampler.
bool has_finite_tag(StochasticNode *node)
{
  for (DeterministicNode const *child : *node->deterministicChildren()) {
    LogicalNode const *logical = dynamic_cast<LogicalNode const *>(child);
    if (logical == nullptr || logical->parents().empty() || logical->parents()[0] != node) continue;
    std::vector<std::string> names(logical->parents().size(), "w");
    names[0] = "x";
    const std::string text = logical->deparse(names);
    if (text.find("bt_finite_tag(") != std::string::npos ||
        text.find("bt_finite_witness_tag(") != std::string::npos) return true;
  }
  return false;
}

bool has_selection_child(SingletonGraphView const &view, SelNormCoarseSettings const &settings)
{
  for (StochasticNode *child : view.stochasticChildren()) {
    if (child_kind(child, settings) != ChildKind::other) return true;
  }
  return false;
}

}

std::vector<Sampler *> ScreenedFiniteFactory::makeSamplers(
    std::list<StochasticNode *> const &nodes, Graph const &graph) const
{
  std::vector<Sampler *> samplers;
  const CorrectedSamplerConfig config = corrected_sampler_configuration();
  if (!config.enabled) return samplers;
  const CorrectedTestFlags flags = corrected_test_flags();
  if (flags.disable_new || flags.disable_finite) return samplers;
  const std::string sampler_name(screened_finite_name);
  try {
    for (StochasticNode *node : nodes) {
      if (node->length() != 1 || !node->isDiscreteValued() || !node->fullRank()) continue;
      const std::string distribution = node->distribution()->name();
      if (distribution != "dcat" && distribution != "dbern") continue;
      if (node->lowerBound() != nullptr || node->upperBound() != nullptr || !isSupportFixed(node)) continue;
      bool fixed = true;
      for (Node const *parent : node->parents()) fixed = fixed && parent->isFixed();
      if (!fixed) continue;
      double lower, upper;
      node->support(&lower, &upper, 1, 0);
      if (!std::isfinite(lower) || !std::isfinite(upper) || upper < lower || lower != std::floor(lower) ||
          upper != std::floor(upper) || upper - lower + 1 > maximum_levels) continue;
      if (has_finite_tag(node)) continue;
      std::unique_ptr<SingletonGraphView> view(new SingletonGraphView(node, graph));
      if (!has_selection_child(*view, config.settings)) continue;
      const unsigned int levels = static_cast<unsigned int>(upper - lower + 1);
      std::vector<std::unique_ptr<MutableSampleMethod> > owned_methods;
      std::vector<MutableSampleMethod *> methods;
      for (unsigned int chain = 0; chain < node->nchain(); ++chain) {
        owned_methods.emplace_back(new ScreenedFinite(view.get(), chain, config.settings, config.min_acceptance,
                                                      lower, levels, flags.finite_proposal));
        methods.push_back(owned_methods.back().get());
      }
      // The Sampler base adopts the graph, and owned_methods keeps the methods
      // until construction succeeds.
      std::unique_ptr<Sampler> sampler(new MutableSampler(view.release(), methods, sampler_name));
      for (std::unique_ptr<MutableSampleMethod> &method : owned_methods) method.release();
      samplers.push_back(sampler.get());
      sampler.release();
    }
  }
  catch (...) {
    for (Sampler *sampler : samplers) delete sampler;
    throw;
  }
  return samplers;
}

std::string ScreenedFiniteFactory::name() const
{
  return screened_finite_name;
}

ScreenedFiniteFactory *screened_finite_factory()
{
  static ScreenedFiniteFactory factory;
  return &factory;
}

std::vector<double> screened_finite_test_scores()
{
  LiveScreened &screened = live_screened();
  const std::lock_guard<std::mutex> lock(screened.mutex);
  std::vector<double> out;
  for (ScreenedFinite *sampler : screened.live) sampler->append_scores(out);
  return out;
}

}
}
