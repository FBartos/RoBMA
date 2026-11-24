#ifndef ROBMA_GATED_CORRECTED_BLOCK_H
#define ROBMA_GATED_CORRECTED_BLOCK_H

// The corrected gated block sampler (RoBMA::GatedCorrectedBlock): the selection
// weight parameters of one bias branch of a model-averaged fit, declared by
// BayesTools' gate markers (BayesTools_gate_<n> ~ dbt_gate(<node>, <active>)).
//
// Members: the still-unclaimed scalar continuous gated nodes with a fixed
// support that carry the identical <active> node (pointer identity), grouped
// in the order JAGS offers them; with sampler.nested.weights the shared
// coordinates (active under two labels) and the fine branch's own coordinate
// (active under one) form separate blocks. A group is claimed when its members
// do not depend on each other, it has a coarse-capable selection child, and
// every member reaches every coarse-capable child only through the child's
// weights (a weight parameter; a PET or PEESE slope reaches the mean and keeps
// BayesTools::GatedSlice, as does a branch whose rows have no coarse-capable
// child). A group of gamma members that enter the model only through their
// proportions (GammaGroup.h; the eta of a cumulative weight function) is held
// as its proportions, with its sum drawn from its gamma law.
//
// One update: while <active> is 0 every member is drawn from its own
// distribution, as BayesTools::GatedSlice does (the declaration guarantees
// that the likelihood is constant in the members), redrawn up to 100 times
// when a value is not strictly inside its support or a deterministic
// descendant is not finite and otherwise kept; while adapting, a change of the
// children's log density is the declaration error of GatedSlice. Otherwise one
// update of the corrected block kernel (CorrectedBlockKernel.h), whose axes
// adapt only on these updates; until they are learned the coordinate axes have
// widths 2.5 prior standard deviations of each coordinate (trigamma of the
// shapes for gamma and beta members and log-ratio coordinates, 1 otherwise).

#include <sampler/SamplerFactory.h>

namespace jags {
namespace RoBMA {

class GatedCorrectedBlockFactory : public SamplerFactory {
public:
  std::vector<Sampler*> makeSamplers(std::list<StochasticNode*> const &nodes,
                                     Graph const &graph) const override;
  std::string name() const override;
};

// Borrowed process-lifetime factory. Module registration must not delete it.
GatedCorrectedBlockFactory *gated_corrected_block_factory();

}
}
#endif
