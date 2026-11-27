#include "../fp-contract-off.h"
#include "GammaGroup.h"

#include <distribution/Distribution.h>
#include <graph/AggNode.h>
#include <graph/DeterministicNode.h>
#include <graph/LogicalNode.h>
#include <graph/Node.h>

#include <algorithm>
#include <cmath>
#include <list>
#include <map>
#include <set>
#include <string>

namespace jags {
namespace RoBMA {
namespace {

// The children of a node; JAGS keeps the accessors non-const.
std::list<DeterministicNode *> const &deterministic_children(Node const *node)
{
  return *const_cast<Node *>(node)->deterministicChildren();
}

std::list<StochasticNode *> const &stochastic_children(Node const *node)
{
  return *const_cast<Node *>(node)->stochasticChildren();
}

bool is_sum_of(Node const *node, Node const *aggregate)
{
  LogicalNode const *logical = dynamic_cast<LogicalNode const *>(node);
  if (logical == nullptr || logical->length() != 1 || logical->parents().size() != 1 ||
      logical->parents()[0] != aggregate) return false;
  const std::string text = logical->deparse(std::vector<std::string>(1, "x"));
  return text == "sum(x)" || text == "(sum(x))";
}

bool is_ratio_of(Node const *node, Node const *sum)
{
  LogicalNode const *logical = dynamic_cast<LogicalNode const *>(node);
  if (logical == nullptr || logical->length() != 1 || logical->parents().size() != 2 ||
      logical->parents()[1] != sum) return false;
  std::vector<std::string> names;
  names.push_back("x");
  names.push_back("y");
  const std::string text = logical->deparse(names);
  return text == "x/y" || text == "(x/y)" || text == "x / y" || text == "(x / y)";
}

bool is_marker(StochasticNode const *node)
{
  const std::string name = node->distribution()->name();
  return name == "dbt_gate" || name == "dbt_block";
}

// A scalar dgamma member without truncation and with fixed, valid parameters.
bool gamma_member(StochasticNode const *node, double *shape, double *rate)
{
  if (node->length() != 1 || node->distribution()->name() != "dgamma" ||
      node->lowerBound() != nullptr || node->upperBound() != nullptr ||
      node->parents().size() != 2) return false;
  Node const *a = node->parents()[0];
  Node const *r = node->parents()[1];
  if (!a->isFixed() || !r->isFixed() || a->length() != 1 || r->length() != 1) return false;
  *shape = a->value(0)[0];
  *rate = r->value(0)[0];
  return std::isfinite(*shape) && *shape > 0 && std::isfinite(*rate) && *rate > 0;
}

}

std::vector<GammaGroup> find_gamma_groups(std::vector<StochasticNode *> const &members)
{
  std::map<Node const *, unsigned int> position;
  std::map<Node const *, std::pair<double, double> > parameters;
  for (unsigned int i = 0; i < members.size(); ++i) {
    double shape, rate;
    if (gamma_member(members[i], &shape, &rate)) {
      position[members[i]] = i;
      parameters[members[i]] = std::make_pair(shape, rate);
    }
  }
  std::vector<GammaGroup> groups;
  std::set<Node const *> grouped;
  for (StochasticNode *member : members) {
    if (position.count(member) == 0 || grouped.count(member) > 0) continue;
    // the one aggregate of the member
    AggNode const *aggregate = nullptr;
    unsigned int aggregates = 0;
    for (DeterministicNode const *child : deterministic_children(member)) {
      AggNode const *candidate = dynamic_cast<AggNode const *>(child);
      if (candidate != nullptr) {
        aggregate = candidate;
        ++aggregates;
      }
    }
    if (aggregates != 1) continue;
    // the aggregate is exactly a group of distinct candidate members
    std::vector<Node const *> const &parts = aggregate->parents();
    if (parts.size() < 2 || aggregate->length() != parts.size()) continue;
    std::set<Node const *> group_nodes(parts.begin(), parts.end());
    bool valid = group_nodes.size() == parts.size();
    for (std::size_t i = 0; valid && i < parts.size(); ++i) {
      valid = position.count(parts[i]) > 0 && grouped.count(parts[i]) == 0 &&
        aggregate->offsets()[i] == 0 && parameters[parts[i]].second == parameters[parts[0]].second;
    }
    // its only child is the sum, whose only children are ratios member / sum
    if (!valid || !stochastic_children(aggregate).empty() ||
        deterministic_children(aggregate).size() != 1) continue;
    Node const *sum = deterministic_children(aggregate).front();
    if (!is_sum_of(sum, aggregate) || !stochastic_children(sum).empty()) continue;
    std::set<Node const *> ratios;
    for (DeterministicNode const *child : deterministic_children(sum)) {
      if (!is_ratio_of(child, sum) || group_nodes.count(child->parents()[0]) == 0) {
        valid = false;
        break;
      }
      ratios.insert(child);
    }
    // the members' only children: markers, the aggregate and their own ratios
    for (std::size_t i = 0; valid && i < parts.size(); ++i) {
      for (StochasticNode const *child : stochastic_children(parts[i])) {
        if (!is_marker(child)) valid = false;
      }
      for (DeterministicNode const *child : deterministic_children(parts[i])) {
        if (child == aggregate) continue;
        if (ratios.count(child) == 0 || child->parents()[0] != parts[i]) valid = false;
      }
    }
    if (!valid) continue;
    GammaGroup group;
    for (Node const *part : parts) group.members.push_back(position[part]);
    std::sort(group.members.begin(), group.members.end());
    group.rate = parameters[parts[0]].second;
    group.shape_sum = 0.0;
    for (unsigned int index : group.members) {
      const double shape = parameters[members[index]].first;
      group.shapes.push_back(shape);
      group.shape_sum += shape;
    }
    for (Node const *part : parts) grouped.insert(part);
    groups.push_back(group);
  }
  return groups;
}

}
}
