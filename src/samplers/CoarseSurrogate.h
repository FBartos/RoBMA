#ifndef ROBMA_COARSE_SURROGATE_H
#define ROBMA_COARSE_SURROGATE_H

// The surrogate density terms of RoBMA's corrected samplers: which stochastic
// children of a sampled node are selection blocks with a surrogate evaluator
// (child_kind), and the surrogate log density of one such child at a chain's
// current values (SelectedChild::logDensity, which also flags a surrogate that
// took an approximate route). They work on a child node and a chain, not on a
// GraphView, so a scalar and a block sampler use the same terms; the settings
// are the sampler's frozen SelNormCoarseSettings snapshot.

#include "../distributions/DSELNORMMVSTEP.h"
#include "../distributions/DSELNORMCLUSTERSTEP.h"
#include "../distributions/DSELNORMFACTORSTEP.h"
#include "../selnorm/selnorm.h"
#include "../selnorm/selnorm-mv.h"

#include <graph/StochasticNode.h>
#include <distribution/Distribution.h>

#include <algorithm>
#include <limits>
#include <stdexcept>
#include <vector>

namespace jags {
namespace RoBMA {

// How a stochastic child of the sampled node enters the surrogate density: a
// selection block with a surrogate evaluator (the dense multivariate block, or a
// cluster or factor block with a coarse rung), or an exact term.
enum class ChildKind { other, mv, cluster, factor };

// Parameter positions of the cluster and factor distributions (they share a
// layout up to the loading): kernel mode, vector rule, quadrature orders.
const unsigned int cluster_factor_mode = 10;
const unsigned int cluster_factor_orders = 13;
const unsigned int cluster_factor_rule = 19;

// A cluster or factor child is coarse-capable when its fixed controls allow a
// coarse rule: a vector rule not fixed at a best-weight value, a kernel mode not
// fixed at NORMAL, and, for a factor block, rule min(max_rules, rules) of its
// schedule affordable at the uncompacted rank (a compacted block or a forest
// only costs less). In model-averaged fits the vector rule is a node of the bias
// indicator that is product in every state of a product-only model, so a
// non-fixed rule counts as capable; in a state with a best-weight rule the
// child's surrogate is its exact density. The decision reads fixed parents only,
// so it is the same in every chain. Any other cluster or factor child stays an
// exact term.
inline ChildKind child_kind(StochasticNode const *child, SelNormCoarseSettings const &settings)
{
  Distribution const *distribution = child->distribution();
  if (dynamic_cast<DSELNORMMVSTEP const *>(distribution) != nullptr) return ChildKind::mv;
  const bool cluster = dynamic_cast<DSELNORMCLUSTERSTEP const *>(distribution) != nullptr;
  const bool factor = !cluster && dynamic_cast<DSELNORMFACTORSTEP const *>(distribution) != nullptr;
  if (!cluster && !factor) return ChildKind::other;
  std::vector<Node const *> const &parents = child->parents();
  if (parents.size() < distribution->npar()) return ChildKind::other;
  Node const *rule = parents[cluster_factor_rule];
  if (rule->isFixed() && rule->value(0)[0] != SELVECTOR_PRODUCT) return ChildKind::other;
  Node const *mode = parents[cluster_factor_mode];
  if (mode->isFixed() && mode->value(0)[0] == SELKERNEL_NORMAL) return ChildKind::other;
  Node const *orders = parents[cluster_factor_orders];
  if (!orders->isFixed() || orders->length() < 1) return ChildKind::other;
  if (factor) {
    const unsigned int rank = parents[2]->length() / parents[0]->length();
    const unsigned int rung = std::min(settings.max_rules, orders->length());
    if (rung < 1 || !selnorm_factor_rule_affordable(
          static_cast<int>(orders->value(0)[rung - 1]), static_cast<int>(rank))) {
      return ChildKind::other;
    }
  }
  return cluster ? ChildKind::cluster : ChildKind::factor;
}

struct SelectedChild {
  StochasticNode *node;
  ChildKind kind;
  // The distribution of the child, of the class the kind names (checked once).
  Distribution const *distribution;
  std::vector<Node const *> parents;
  std::vector<unsigned int> lengths;
  mutable std::vector<double const *> parameters;
  bool controls_fixed = true;
  mutable bool controls_checked = false;
  mutable bool controls_valid = false;

  // The parents that may change during sampling: mean, covariance (residual SD
  // and loadings), weights, kernel mode and, for cluster and factor children,
  // the vector rule (a mixture of weight rules composes it from the bias
  // indicator). Everything else is a control.
  static bool dynamic_parent(ChildKind kind, unsigned int i)
  {
    if (kind == ChildKind::mv) return i == 0 || i == 1 || i == 3 || i == 9;
    return i == 0 || i == 1 || i == 2 || i == 4 || i == cluster_factor_mode ||
      i == cluster_factor_rule;
  }

  SelectedChild(StochasticNode *child, ChildKind kind_of_child)
    : node(child), kind(kind_of_child), distribution(child->distribution())
  {
    const bool consistent =
      (kind == ChildKind::mv && dynamic_cast<DSELNORMMVSTEP const *>(distribution) != nullptr) ||
      (kind == ChildKind::cluster && dynamic_cast<DSELNORMCLUSTERSTEP const *>(distribution) != nullptr) ||
      (kind == ChildKind::factor && dynamic_cast<DSELNORMFACTORSTEP const *>(distribution) != nullptr);
    if (!consistent) {
      throw std::logic_error("Selection child has no surrogate evaluator.");
    }
    const unsigned int count = distribution->npar();
    if (child->parents().size() < count) {
      throw std::logic_error("Selection child is missing distribution parameters.");
    }
    // StochasticNode appends truncation bounds after its distribution parents.
    parents.assign(child->parents().begin(), child->parents().begin() + count);
    lengths.resize(count);
    parameters.resize(count);
    for (unsigned int i = 0; i < count; ++i) {
      lengths[i] = parents[i]->length();
      if (!dynamic_parent(kind, i) && !parents[i]->isFixed()) {
        controls_fixed = false;
      }
    }
  }

  bool checkControl() const
  {
    switch (kind) {
    case ChildKind::mv:
      return static_cast<DSELNORMMVSTEP const *>(distribution)->checkControlParameterValue(parameters, lengths);
    case ChildKind::cluster:
      return static_cast<DSELNORMCLUSTERSTEP const *>(distribution)->checkControlParameterValue(parameters, lengths);
    default:
      return static_cast<DSELNORMFACTORSTEP const *>(distribution)->checkControlParameterValue(parameters, lengths);
    }
  }

  bool checkDynamic() const
  {
    switch (kind) {
    case ChildKind::mv:
      return static_cast<DSELNORMMVSTEP const *>(distribution)->checkDynamicParameterValue(parameters, lengths);
    case ChildKind::cluster:
      return static_cast<DSELNORMCLUSTERSTEP const *>(distribution)->checkDynamicParameterValue(parameters, lengths);
    default:
      return static_cast<DSELNORMFACTORSTEP const *>(distribution)->checkDynamicParameterValue(parameters, lengths);
    }
  }

  double surrogate(unsigned int chain, SelNormCoarseSettings const &settings,
                   bool *not_exact) const
  {
    double const *x = node->value(chain);
    double const *lower = node->lowerLimit(chain);
    double const *upper = node->upperLimit(chain);
    switch (kind) {
    case ChildKind::mv:
      return static_cast<DSELNORMMVSTEP const *>(distribution)->surrogateLogDensity(
        x, node->length(), PDF_LIKELIHOOD, parameters, lengths, lower, upper, settings,
        not_exact);
    case ChildKind::cluster:
      return static_cast<DSELNORMCLUSTERSTEP const *>(distribution)->surrogateLogDensity(
        x, node->length(), PDF_LIKELIHOOD, parameters, lengths, lower, upper, settings,
        not_exact);
    default:
      return static_cast<DSELNORMFACTORSTEP const *>(distribution)->surrogateLogDensity(
        x, node->length(), PDF_LIKELIHOOD, parameters, lengths, lower, upper, settings,
        not_exact);
    }
  }

  // not_exact is set when the surrogate used a normalizer that is not exact.
  double logDensity(unsigned int chain, SelNormCoarseSettings const &settings,
                    bool *not_exact) const
  {
    for (unsigned int i = 0; i < parents.size(); ++i) {
      parameters[i] = parents[i]->value(chain);
    }
    bool valid;
    if (controls_fixed) {
      // The first surrogate call follows initialized-graph checkFinite().
      // This certificate belongs to one SelectedChild in one chain's sampler.
      if (!controls_checked) {
        controls_valid = checkControl();
        controls_checked = true;
      }
      valid = controls_valid && checkDynamic();
    } else {
      // Any nonfixed control retains the complete original node gate.
      valid = node->checkParentValues(chain);
    }
    if (!valid) return -std::numeric_limits<double>::infinity();
    return surrogate(chain, settings, not_exact);
  }
};

}
}
#endif
