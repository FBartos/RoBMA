#ifndef ROBMA_COARSE_CORRECTED_SLICE_H
#define ROBMA_COARSE_CORRECTED_SLICE_H

#include <sampler/SamplerFactory.h>
#include <cstdint>
#include <mutex>
#include "../selnorm/selnorm-mv.h"

namespace jags {
namespace RoBMA {

struct CoarseCorrectedSliceConfig {
  bool enabled;
  SelNormCoarseSettings settings;
};

struct CoarseCorrectedSliceStats {
  std::uint64_t proposed;
  std::uint64_t accepted;
  std::uint64_t correction_rejected;
};

// The factory starts disabled. Configuration controls only future model
// construction; every existing method owns its immutable settings snapshot.
class CoarseCorrectedSliceFactory : public SamplerFactory {
  mutable std::mutex _configuration_mutex;
  CoarseCorrectedSliceConfig _configuration;
public:
  explicit CoarseCorrectedSliceFactory(
    SelNormCoarseSettings const &settings = SelNormCoarseSettings());
  void configure(bool enabled, SelNormCoarseSettings const &settings);
  CoarseCorrectedSliceConfig configuration() const;
  std::vector<Sampler*> makeSamplers(std::list<StochasticNode*> const &nodes,
                                     Graph const &graph) const;
  std::string name() const;
};

// Borrowed process-lifetime factory. Module registration must not delete it.
CoarseCorrectedSliceFactory *coarse_corrected_slice_factory();
// Counters cover all chains/samplers in this process; reset only between updates.
CoarseCorrectedSliceStats coarse_corrected_slice_stats(bool reset = false);

}
}
#endif
