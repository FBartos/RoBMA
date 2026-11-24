#ifndef ROBMA_SCREENED_FINITE_H
#define ROBMA_SCREENED_FINITE_H

// The screened finite sampler (RoBMA::ScreenedFinite) and its kernel.
//
// The kernel (ScreenedFiniteKernel) updates one scalar finite-support node
// whose stochastic children include a coarse-capable selection child. For
// every level l of nonzero prior mass it scores
//   s_l = log prior(l) + sum over the stochastic children of the surrogate log
//         density with the node at l (children without a surrogate contribute
//         their exact density),
// proposes a level k' on these scores, p_l = exp(s_l) / sum exp(s), by one of
// two moves reversible for the surrogate law, chosen at each update by a coin
// that does not depend on the state:
//   - with probability 3/4 Liu's Metropolized Gibbs sampler (a level other
//     than the current level k with probability p_k' / (1 - p_k), accepted on
//     the surrogate with min(1, (1 - p_k) / (1 - p_k')); the identity when no
//     other level has mass);
//   - with probability 1/4 the plain surrogate full-conditional draw (k' with
//     probability p_k', the current level included: an independence proposal),
// and, when k' != k, accepts it with probability
//   min(1, exp((e_k' - s_k') - (e_k - s_k))),
// e the same score with the exact densities. Each move with its correction
// leaves the target invariant, and so does their mixture. The plain draw's
// holding probability makes the chain aperiodic: Liu's step alone between two
// equally scored levels swaps them at every update, so an evenly thinned chain
// would keep the level it started in. The update is exact for any
// surrogate that is one function of the evaluated state (s_l depends on the
// values of the graph with the node at l, never on the level the update
// started from, so that both ends of a move see the same scores) and is finite
// wherever the exact score is (the surrogate falls back to the exact density
// where its rule is unavailable); the row cache (StepRowCache.h) computes each
// row at the state it evaluates. A self-proposal needs no exact evaluation.
// The test hook selects either move alone. Every score uses the current values
// of all other nodes. The exact kernel enumerates the
// levels with the exact scores, with the arithmetic and the one uniform draw
// of JAGS's base::Finite.
// A sampler that owns an indicator jointly with other nodes can call the
// kernel with a setter that writes the indicator and its descendants.
//
// The sampler claims a scalar dcat (fixed probabilities) or dbern (fixed
// probability) node without truncation, with a fixed finite support, a
// coarse-capable selection child among its stochastic children and no finite
// readout tag (a readout that sampler.probabilities requests keeps the
// enumerating BayesTools sampler, whose recorded probabilities are the exact
// conditionals). RoBMA's joint samplers are asked first and keep the nodes they
// claim. Its kernel decision at the end of adaptation follows the evidence rule
// of CorrectedTally (changed proposals whose surrogate at either level was not
// exact).

#include "CoarseSurrogate.h"
#include "StepRowCache.h"

#include <sampler/SamplerFactory.h>

#include <functional>
#include <vector>

namespace jags {
struct RNG;
namespace RoBMA {

// The screened indicator's proposals (the test hook's finite_proposal): the
// plain surrogate full conditional, the same without the exact correction (a
// test mutant that is not exact), Liu's Metropolized Gibbs, and the default
// mixture of the plain draw (1/4) and Liu's step (3/4).
enum ScreenedProposal { proposal_plain = 0, proposal_uncorrected = 1, proposal_liu = 2, proposal_mixture = 3 };

class ScreenedFiniteKernel {
public:
  // proposal: a ScreenedProposal.
  // row_cache: singleton step-selection children enter the surrogate scores
  // through the row cache (StepRowCache.h), which evaluates each row at the
  // state of the level being scored.
  ScreenedFiniteKernel(StochasticNode *node, std::vector<StochasticNode *> const &children,
                       unsigned int chain, SelNormCoarseSettings const &settings, double lower,
                       unsigned int levels, int proposal, bool row_cache = false);

  // The outcome of one corrected update: whether another level was proposed
  // (with an exact correction), accepted, and whether the surrogate at either
  // level was not exact (evidence).
  struct Outcome {
    bool proposed = false;
    bool accepted = false;
    bool evidence = false;
    double leave = 0;   // the surrogate full-conditional probability of leaving the current level
  };

  // `set_value(v)` writes the node and its deterministic descendants. Both
  // updates leave the node at the drawn level; exceptions propagate (the caller
  // restores the node). enumerate() returns the exact full-conditional
  // probability of leaving the current level.
  Outcome update(RNG *rng, std::function<void(double)> const &set_value);
  double enumerate(RNG *rng, std::function<void(double)> const &set_value);

  // Test hook and frozen-state checks: the surrogate and exact scores of every
  // level at the current values of the other nodes (prior-free levels -Inf),
  // leaving the node at `restore`.
  void scores(std::function<void(double)> const &set_value, double restore,
              std::vector<double> *surrogate, std::vector<double> *exact);

private:
  StochasticNode *_node;
  std::vector<StochasticNode *> _children;
  const unsigned int _chain;
  const SelNormCoarseSettings _settings;
  const double _lower;
  const unsigned int _levels;
  const int _proposal;
  std::vector<SelectedChild> _selected;
  std::vector<int> _child_selected;
  StepRowCache _rows;

  double prior_term() const;
  double surrogate_score(bool *not_exact);
  double exact_score() const;
};

class ScreenedFiniteFactory : public SamplerFactory {
public:
  std::vector<Sampler*> makeSamplers(std::list<StochasticNode*> const &nodes,
                                     Graph const &graph) const override;
  std::string name() const override;
};

// Borrowed process-lifetime factory. Module registration must not delete it.
ScreenedFiniteFactory *screened_finite_factory();

// Test hook: the surrogate and exact scores of every level of the live
// screened samplers, at the current state of their chains: per sampler its
// chain, the number of levels L, then L surrogate and L exact scores.
std::vector<double> screened_finite_test_scores();

// The evidence record of the screened indicators of this process since the
// last clearing read, over the updates after adaptation: per indicator and
// chain the rows (serial, chain, levels and lower bound of the support, the
// node's value after its last recorded update, updates, changes of level, sum
// and maximum of the full-conditional probability of leaving the current
// level). The serial identifies a sampler (numbered as the factory creates
// them, across models); nodes have no names here, so R matches each serial to
// its node by the support and the chains' final values.
std::vector<char const *> const &screened_finite_evidence_names();
std::vector<double> screened_finite_evidence(bool clear);

}
}
#endif
