#ifndef ROBMA_COARSE_CORRECTED_BLOCK_H
#define ROBMA_COARSE_CORRECTED_BLOCK_H

// The corrected block sampler (RoBMA::CoarseCorrectedBlock): a block of scalar
// parameters declared by BayesTools' dbt_block markers, with at least one
// coarse-capable selection child, gets one factor slice sweep on the surrogate
// density (BayesTools' kernel, BayesTools/factor_slice_kernel.h) in a random
// direction order, corrected by one Metropolis step against the full joint
// conditional.

#include <sampler/SamplerFactory.h>

#include <cstdint>
#include <vector>

namespace jags {
namespace RoBMA {

// Groups the still-unclaimed marked nodes by block id as BayesTools' factor
// slice factory does and claims a group when the shared corrected sampler
// configuration (CorrectedSampler.h) is enabled and the group has a
// coarse-capable selection child; every other group is left to BayesTools.
class CoarseCorrectedBlockFactory : public SamplerFactory {
public:
  std::vector<Sampler*> makeSamplers(std::list<StochasticNode*> const &nodes,
                                     Graph const &graph) const override;
  std::string name() const override;
};

// Borrowed process-lifetime factory. Module registration must not delete it.
CoarseCorrectedBlockFactory *coarse_corrected_block_factory();

// Test hooks, not package interface. fixed_order: samplers built from now on
// sweep the directions in their fixed order in corrected mode (a kernel that is
// not exact, for the power of the certification fixtures). trace: every
// density evaluation of a block sampler is recorded as
// (kind 0 surrogate / 1 full, chain, d, member values, log density without the
// Jacobian). The exhausted count is the number of slice intervals exhausted in
// floating point in corrected mode since the last take.
void corrected_block_test_set(bool fixed_order, bool trace);
bool corrected_block_test_fixed_order();
bool corrected_block_test_trace();
std::vector<double> corrected_block_test_take_trace();
std::uint64_t corrected_block_test_take_exhausted();
// The current widths of the live block samplers, in construction order:
// per sampler its chain, d and the d widths of its axes.
std::vector<double> corrected_block_test_widths();

}
}
#endif
