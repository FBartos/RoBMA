#ifndef ROBMA_JOINT_SAMPLER_TAG_H
#define ROBMA_JOINT_SAMPLER_TAG_H

#include <distribution/Distribution.h>
#include <function/VectorFunction.h>
#include <graph/Graph.h>
#include <graph/StochasticNode.h>
#include <graph/VectorLogicalNode.h>

#include <cmath>
#include <string>
#include <vector>

namespace jags::RoBMA {

// A terminal declaration keeps its arguments visible to the sampler factory
// without becoming a density factor. Providers may check argument lengths.
class JointSamplerTag : public VectorFunction {
public:
  JointSamplerTag(std::string const &name, unsigned long arity)
    : VectorFunction(name, arity) {}

  unsigned long length(std::vector<unsigned long> const &,
                       std::vector<double const *> const &) const override
  {
    return 1;
  }

  void evaluate(double *value, std::vector<double const *> const &,
                std::vector<unsigned long> const &) const override
  {
    value[0] = 0.0;
  }
};

inline bool fixed_scalar(Node const *node)
{
  return node && node->length() == 1 && node->isFixed();
}

struct CategoricalWeights {
  double const *values = nullptr;
  unsigned long size = 0;
  double total = 0;

  explicit operator bool() const { return values != nullptr; }
};

// A zero required count admits any nonempty categorical support. The finite
// sum check uses double precision, independently of the caller's log masses.
inline CategoricalWeights categorical_weights(Node const *node, unsigned long required_count = 0)
{
  auto const *cat = dynamic_cast<StochasticNode const *>(node);
  if (!cat || cat->length() != 1 || cat->distribution()->name() != "dcat" ||
      cat->lowerBound() || cat->upperBound() || cat->parents().size() != 1) return {};
  Node const *weights = cat->parents().front();
  const auto size = weights->length();
  if (!weights->isFixed() || size == 0 || (required_count && size != required_count)) return {};
  double const *values = weights->value(0);
  double total = 0;
  for (unsigned long i = 0; i < size; ++i) {
    if (!(values[i] > 0 && std::isfinite(values[i]))) return {};
    total += values[i];
  }
  if (!std::isfinite(total)) return {};
  return {values, size, total};
}

inline bool terminal_tag(Node const *node, std::string const &name, unsigned long arity,
                         Graph const &graph)
{
  if (!dynamic_cast<VectorLogicalNode const *>(node) || graph.contains(node) ||
      node->length() != 1 || node->parents().size() != arity) return false;

  // JAGS exposes these read-only child lists through non-const accessors.
  auto *child = const_cast<Node *>(node);
  if (!child->stochasticChildren()->empty() || !child->deterministicChildren()->empty()) return false;

  // LogicalNode::deparse uses the compiled function's name. Use synthetic
  // argument names so recognition is independent of the caller's labels.
  std::string expected = "(" + name + "(";
  for (unsigned long i = 0; i < arity; ++i) expected += i == 0 ? "x" : ",x";
  expected += "))";
  return node->deparse(std::vector<std::string>(arity, "x")) == expected;
}

inline std::vector<Node const *> direct_terminal_tags(StochasticNode *first_parent,
                                                       std::string const &name,
                                                       unsigned long arity, Graph const &graph)
{
  std::vector<Node const *> tags;
  if (!first_parent) return tags;
  for (DeterministicNode const *child : *first_parent->deterministicChildren()) {
    auto const &parents = child->parents();
    if (!parents.empty() && parents.front() == first_parent && terminal_tag(child, name, arity, graph)) {
      tags.push_back(child);
    }
  }
  return tags;
}

}

#endif
