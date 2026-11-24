#ifndef ROBMA_GAMMA_GROUP_H
#define ROBMA_GAMMA_GROUP_H

// The gamma-sum reparameterization of the corrected block samplers.
//
// Independent members eta_i ~ dgamma(a_i, r) with fixed shapes and one common
// fixed rate that reach every other node only through their normalized
// proportions p_i = eta_i / sum(eta) have a sum T = sum(eta) ~ Gamma(sum a, r)
// that is independent of the proportions, p ~ Dirichlet(a), and of the
// likelihood. A block holds such a group as T and the m - 1 additive log-ratio
// coordinates z_i = log(eta_i / eta_m) (i < m), and writes eta_i = T p_i with
// p_i = exp(z_i) / (1 + sum_k exp(z_k)), p_m = 1 / (1 + sum_k exp(z_k)).
// The absolute Jacobian of (T, z) -> eta is T^(m-1) prod_{i<=m} p_i, so the
// members' own densities at eta, plus (m - 1) log T + sum_i log p_i, are the
// joint density of (T, z): Gamma(T; sum a, r) times the Dirichlet density of
// p in z, prod p_i^a_i up to a constant.
//
// The grouping is decided on the compiled graph (find_gamma_groups): the
// members are scalar dgamma nodes of the block without truncation, with fixed
// shape and rate parameters and bitwise equal rates; the only stochastic
// children of each member are sampler markers (dbt_gate, dbt_block: log
// density 0); the only deterministic children of a member are one aggregate
// node of exactly the group (every member once) and division nodes
// member / S, where S = sum(aggregate) is the only child of the aggregate and
// has no children but such division nodes. Every path from the members to the
// rest of the model then passes through a ratio eta_i / sum(eta), which is
// invariant to the scale T. Anything else declines the group, and its members
// keep their own coordinates.

#include <graph/StochasticNode.h>

#include <vector>

namespace jags {
namespace RoBMA {

struct GammaGroup {
  // the members' positions in the block, ascending (the last is the reference
  // of the log-ratio coordinates)
  std::vector<unsigned int> members;
  std::vector<double> shapes;
  double rate;
  double shape_sum;
};

std::vector<GammaGroup> find_gamma_groups(std::vector<StochasticNode *> const &members);

}
}
#endif
