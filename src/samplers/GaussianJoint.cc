#include "GaussianJoint.h"
#include "GaussianJointKernel.h"

#include <JRmath.h>
#include <graph/AggNode.h>
#include <graph/ScalarLogicalNode.h>
#include <graph/StochasticNode.h>
#include <rng/RNG.h>
#include <sampler/GraphView.h>
#include <sampler/ImmutableSampleMethod.h>
#include <sampler/ImmutableSampler.h>
#include <sampler/SamplerFactory.h>

#include <memory>
#include <set>
#include <stdexcept>
#include <utility>

namespace jags {
namespace RoBMA {
namespace {

namespace kernel = ::robma_gaussian_joint;
constexpr char tag_name[] = "robma_gaussian_joint";
constexpr char sampler_name[] = "RoBMA::GaussianJoint";
constexpr unsigned int root_count = 6;
constexpr unsigned int labels_parent = 6;
constexpr unsigned int moderator_parent = 7;
constexpr unsigned int outcome_parent = 8;
constexpr unsigned int se_parent = 9;
constexpr unsigned int sign_parent = 10;
constexpr unsigned int tag_arity = 11;

using Labels = std::array<double, 3>;

struct HalfNormalQuantile {
  double sd;
  double operator()(double upper_probability) const
  {
    return qnorm(upper_probability / 2, 0, sd, 0, 0);
  }
};

bool has_distribution(StochasticNode const &node, const char *name)
{
  return node.distribution()->name() == name && node.length() == 1;
}

bool read_prior(StochasticNode const &gate, StochasticNode const &slab,
                unsigned int role, double active_label, kernel::Model &model)
{
  if (!has_distribution(gate, "dcat") || !has_distribution(slab, "dnorm") ||
      gate.lowerBound() != nullptr || gate.upperBound() != nullptr ||
      gate.parents().size() != 1 || gate.parents()[0]->length() != 2 ||
      !gate.parents()[0]->isFixed() ||
      !(active_label == 1 || active_label == 2)) {
    return false;
  }
  const double *weights = gate.parents()[0]->value(0);
  if (!std::isfinite(weights[0]) || !std::isfinite(weights[1]) ||
      !(weights[0] > 0 && weights[1] > 0) || !std::isfinite(weights[0] + weights[1])) {
    return false;
  }
  const unsigned int on = static_cast<unsigned int>(active_label) - 1;
  const long double log_total = std::log(static_cast<long double>(weights[0]) + weights[1]);
  model.log_prior[role][0] = std::log(static_cast<long double>(weights[1 - on])) - log_total;
  model.log_prior[role][1] = std::log(static_cast<long double>(weights[on])) - log_total;

  const bool is_scale = role == 2;
  if (slab.parents().size() != (is_scale ? 3 : 2) ||
      slab.upperBound() != nullptr ||
      !fixed_scalar(slab.parents()[0]) || !fixed_scalar(slab.parents()[1]) ||
      (is_scale ? !fixed_scalar(slab.lowerBound()) || slab.lowerBound()->value(0)[0] != 0
                : slab.lowerBound() != nullptr)) {
    return false;
  }
  const double mean = slab.parents()[0]->value(0)[0];
  const double precision = slab.parents()[1]->value(0)[0];
  if (!std::isfinite(mean) || (is_scale && mean != 0) ||
      !std::isfinite(precision) || !(precision > 0)) {
    return false;
  }
  if (!is_scale) model.mean[role] = mean;
  model.sd[role] = 1 / std::sqrt(precision);
  return std::isfinite(model.sd[role]) && model.sd[role] > 0;
}

// The observed response vector of the marker is a JAGS aggregate of the
// likelihood nodes. Resolve those exact nodes in observation order.
std::vector<StochasticNode const *> read_observations(Node const *outcome)
{
  std::vector<StochasticNode const *> observed;
  if (outcome->length() == 1) {
    auto const *node = dynamic_cast<StochasticNode const *>(outcome);
    if (node != nullptr) observed.push_back(node);
    return observed;
  }
  auto const *aggregate = dynamic_cast<AggNode const *>(outcome);
  if (aggregate == nullptr) return observed;
  for (unsigned int i = 0; i < aggregate->length(); ++i) {
    auto const *node = dynamic_cast<StochasticNode const *>(aggregate->parents()[i]);
    if (node == nullptr || aggregate->offsets()[i] != 0) return {};
    observed.push_back(node);
  }
  return observed;
}

bool fixed_value(Node const *node, double value)
{
  return fixed_scalar(node) && node->value(0)[0] == value;
}

// JAGS exposes the operation of a scalar logical node through deparse().
// Supply placeholders instead of model names, then verify the actual
// parent pointers separately. This is a deliberately narrow match of the
// arithmetic emitted by the package, with no evaluation at trial states.
bool binary(Node const *node, const char *syntax)
{
  return dynamic_cast<ScalarLogicalNode const *>(node) != nullptr &&
    node->parents().size() == 2 && node->deparse({"a", "b"}) == syntax;
}

bool gated_component(Node const *node, Node const *gate, Node const *raw, double label)
{
  if (!binary(node, "(a*b)")) return false;
  const auto &parents = node->parents();
  const Node *test = parents[1];
  return (raw == nullptr ? fixed_value(parents[0], 0) : parents[0] == raw) &&
    binary(test, "(a==b)") && test->parents()[0] == gate &&
    fixed_value(test->parents()[1], label);
}

bool gated_parameter(Node const *node, Node const *gate, Node const *raw, double label)
{
  if (!binary(node, "(a+b)")) return false;
  const unsigned int on = static_cast<unsigned int>(label) - 1;
  return gated_component(node->parents()[on], gate, raw, label) &&
    gated_component(node->parents()[1 - on], gate, nullptr, 3 - label);
}

bool matching_location(Node const *mean, GraphView const &view,
                       Labels const &labels, double moderator, double effect_sign)
{
  if (effect_sign == -1) {
    if (dynamic_cast<ScalarLogicalNode const *>(mean) == nullptr ||
        mean->parents().size() != 1 || mean->deparse({"a"}) != "(-a)") {
      return false;
    }
    mean = mean->parents()[0];
  }
  if (!binary(mean, "(a+b)")) return false;
  const auto &roots = view.nodes();
  const auto &parents = mean->parents();
  const Node *slope = parents[1];
  return gated_parameter(parents[0], roots[0], roots[1], labels[0]) &&
    binary(slope, "(a*b)") &&
    gated_parameter(slope->parents()[0], roots[2], roots[3], labels[1]) &&
    fixed_value(slope->parents()[1], moderator);
}

bool matching_precision(Node const *precision, GraphView const &view,
                        Labels const &labels, double se)
{
  if (!binary(precision, "(a/b)") || !fixed_value(precision->parents()[0], 1)) return false;
  const Node *variance = precision->parents()[1];
  if (!binary(variance, "(a+b)")) return false;
  const Node *sampling = variance->parents()[0];
  const Node *heterogeneity = variance->parents()[1];
  const auto &roots = view.nodes();
  return binary(sampling, "(a^b)") && fixed_value(sampling->parents()[0], se) &&
    fixed_value(sampling->parents()[1], 2) && binary(heterogeneity, "(a^b)") &&
    gated_parameter(heterogeneity->parents()[0], roots[4], roots[5], labels[2]) &&
    fixed_value(heterogeneity->parents()[1], 2);
}

bool matching_observations(GraphView const &view,
                           std::vector<StochasticNode const *> const &observed,
                           kernel::Model const &model, Labels const &labels)
{
  std::set<StochasticNode const *> children(view.stochasticChildren().begin(),
                                            view.stochasticChildren().end());
  if (observed.empty() || children.size() != observed.size()) {
    return false;
  }
  for (std::size_t i = 0; i < observed.size(); ++i) {
    StochasticNode const *node = observed[i];
    const kernel::Observation &row = model.rows[i];
    if (children.erase(node) != 1 || !isObserved(node) ||
        !has_distribution(*node, "dnorm") ||
        node->lowerBound() != nullptr || node->upperBound() != nullptr ||
        node->parents().size() != 2 ||
        !matching_location(node->parents()[0], view, labels, row.moderator, model.effect_sign) ||
        !matching_precision(node->parents()[1], view, labels, row.se)) {
      return false;
    }
  }
  return true;
}

class GaussianJoint : public ImmutableSampleMethod {
  GraphView const *_view;
  const kernel::Model _model;
  const kernel::Proposal _proposal;
  const Labels _labels;

public:
  GaussianJoint(GraphView const *view, kernel::Model model, kernel::Proposal proposal,
                Labels labels)
    : _view(view), _model(std::move(model)), _proposal(std::move(proposal)), _labels(labels)
  {
  }

  void update(unsigned int chain, RNG *rng) const override
  {
    auto const &nodes = _view->nodes();
    const unsigned int coefficient_cell = (nodes[0]->value(chain)[0] == _labels[0] ? 1 : 0)
      + (nodes[2]->value(chain)[0] == _labels[1] ? 2 : 0);
    const bool scale_on = nodes[4]->value(chain)[0] == _labels[2];
    const double tau = scale_on ? nodes[5]->value(chain)[0] : 0;
    const unsigned int bin = scale_on ? _proposal.bucket(tau) : 0;
    const auto current = _proposal.state(_model, coefficient_cell, bin, tau);
    const auto candidate = _proposal.draw(_model, HalfNormalQuantile{_model.sd[2]}, *rng);
    if (!std::isfinite(current.log_ratio)) {
      throw std::runtime_error("RoBMA Gaussian joint sampler could not evaluate the retained state.");
    }
    const auto &retained = std::isfinite(candidate.log_ratio) &&
      std::log(rng->uniform()) < candidate.log_ratio - current.log_ratio ? candidate : current;

    const double intercept_normal = rng->normal();
    const double slope_normal = rng->normal();
    const auto coefficients = kernel::restore_coefficients(
      _model, retained.conditional, intercept_normal, slope_normal);
    const double raw_tau = retained.bin != 0 ? retained.tau : std::fabs(rng->normal()) * _model.sd[2];
    if (!std::isfinite(coefficients[0]) || !std::isfinite(coefficients[1]) ||
        !std::isfinite(raw_tau)) {
      throw std::runtime_error("RoBMA Gaussian joint sampler could not represent a raw slab draw.");
    }
    const std::array<double, root_count> restored {{
      kernel::active(retained.coefficients, 0) ? _labels[0] : 3 - _labels[0], coefficients[0],
      kernel::active(retained.coefficients, 1) ? _labels[1] : 3 - _labels[1], coefficients[1],
      retained.bin != 0 ? _labels[2] : 3 - _labels[2], raw_tau
    }};
    _view->setValue(restored.data(), restored.size(), chain);
  }
};

std::unique_ptr<Sampler> make_sampler(Node const *tag,
                                     std::list<StochasticNode *> const &available,
                                     Graph const &graph)
{
  auto const &parents = tag->parents();
  for (unsigned int i = root_count; i < tag_arity; ++i) {
    if (!parents[i]->isFixed()) return nullptr;
  }
  std::vector<StochasticNode *> nodes;
  for (unsigned int i = 0; i < root_count; ++i) {
    auto found = std::find(available.begin(), available.end(), parents[i]);
    if (found == available.end() || !(*found)->fullRank() ||
        std::find(nodes.begin(), nodes.end(), *found) != nodes.end()) {
      return nullptr;
    }
    nodes.push_back(*found);
  }
  Labels labels;
  kernel::Model model;
  for (unsigned int role = 0; role < labels.size(); ++role) {
    labels[role] = parents[labels_parent]->value(0)[role];
    if (!read_prior(*nodes[2 * role], *nodes[2 * role + 1], role, labels[role], model)) {
      return nullptr;
    }
  }
  model.effect_sign = parents[sign_parent]->value(0)[0];
  if (!(model.effect_sign == 1 || model.effect_sign == -1)) return nullptr;
  const auto observed = read_observations(parents[outcome_parent]);
  const double *moderator = parents[moderator_parent]->value(0);
  const double *se = parents[se_parent]->value(0);
  for (unsigned int i = 0; i < observed.size(); ++i) {
    const double y = observed[i]->value(0)[0];
    if (!std::isfinite(y) || !std::isfinite(moderator[i]) ||
        !std::isfinite(se[i]) || !(se[i] > 0)) {
      return nullptr;
    }
    model.rows.push_back({y, se[i], moderator[i]});
  }
  auto view = std::make_unique<GraphView>(nodes, graph);
  if (!matching_observations(*view, observed, model, labels)) return nullptr;
  kernel::Proposal proposal;
  const HalfNormalQuantile quantile{model.sd[2]};
  if (!proposal.prepare(model, quantile) || !proposal.coverage(model, quantile).certified) {
    return nullptr;
  }

  auto method = std::make_unique<GaussianJoint>(view.get(), std::move(model),
                                                 std::move(proposal), labels);
  std::string name(sampler_name);
  auto sampler = std::unique_ptr<Sampler>(new ImmutableSampler(view.release(), method.get(), name));
  method.release();
  return sampler;
}

class GaussianJointFactory : public SamplerFactory {
public:
  std::vector<Sampler *> makeSamplers(std::list<StochasticNode *> const &nodes,
                                      Graph const &graph) const override
  {
    for (StochasticNode *node : nodes) {
      const auto tags = direct_terminal_tags(node, tag_name, tag_arity, graph);
      if (tags.size() != 1) continue;
      auto sampler = make_sampler(tags.front(), nodes, graph);
      if (sampler != nullptr) {
        std::vector<Sampler *> claimed{sampler.get()};
        sampler.release();
        return claimed;
      }
    }
    return {};
  }

  std::string name() const override
  {
    return sampler_name;
  }
};

} // namespace

GaussianJointTag::GaussianJointTag() : JointSamplerTag(tag_name, tag_arity) {}

bool GaussianJointTag::checkParameterLength(const std::vector<unsigned long> &len) const
{
  return len.size() == tag_arity &&
    std::all_of(len.begin(), len.begin() + root_count, [](unsigned long n) { return n == 1; }) &&
    len[labels_parent] == 3 && len[moderator_parent] > 0 &&
    len[moderator_parent] == len[outcome_parent] && len[moderator_parent] == len[se_parent] &&
    len[sign_parent] == 1;
}

SamplerFactory *gaussian_joint_factory()
{
  static GaussianJointFactory factory;
  return &factory;
}

} // namespace RoBMA
} // namespace jags
