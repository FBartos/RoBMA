#ifndef ROBMA_COARSE_CORRECTED_SLICE_H
#define ROBMA_COARSE_CORRECTED_SLICE_H

#include <sampler/SamplerFactory.h>

namespace jags {
namespace RoBMA {

// The factory of the corrected scalar sampler. It reads the shared corrected
// sampler configuration (CorrectedSampler.h) at each model construction and
// builds nothing while that configuration is disabled.
class CoarseCorrectedSliceFactory : public SamplerFactory {
public:
  std::vector<Sampler*> makeSamplers(std::list<StochasticNode*> const &nodes,
                                     Graph const &graph) const override;
  std::string name() const override;
};

// Borrowed process-lifetime factory. Module registration must not delete it.
CoarseCorrectedSliceFactory *coarse_corrected_slice_factory();

}
}
#endif
