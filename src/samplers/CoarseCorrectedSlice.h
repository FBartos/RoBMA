#ifndef ROBMA_COARSE_CORRECTED_SLICE_H
#define ROBMA_COARSE_CORRECTED_SLICE_H

#include <sampler/SamplerFactory.h>
#include <cstdint>
#include <mutex>
#include "../selnorm/selnorm-mv.h"

namespace jags {
namespace RoBMA {

// min_acceptance > 0 enables the fallback to the exact slice kernel that each
// sampler decides once, at the end of adaptation; 0 disables it.
struct CoarseCorrectedSliceConfig {
  bool enabled;
  SelNormCoarseSettings settings;
  double min_acceptance = 0.0;
};

struct CoarseCorrectedSliceStats {
  std::uint64_t proposed;
  std::uint64_t accepted;
  std::uint64_t correction_rejected;
};

// The per-fit counters, read and cleared by the runtime report of each fit:
// sampler objects constructed (node x chain, per model build), the fallbacks
// decided at adaptOff() by reason, the proposals and accepted corrections made
// after adaptOff() by samplers that stayed corrected, and the lowest acceptance
// rate of one such sampler over its proposals since the last clearing read (live
// samplers and those destroyed since; only samplers with at least 256 such
// proposals count, has_min_rate false when none does).
struct CoarseCorrectedSliceReport {
  std::uint64_t samplers;
  std::uint64_t fallback_few_proposals;
  std::uint64_t fallback_surrogate_exact;
  std::uint64_t fallback_low_acceptance;
  std::uint64_t sampling_proposed;
  std::uint64_t sampling_accepted;
  bool has_min_rate;
  double min_sampler_rate;
};

// The factory starts disabled. Configuration controls only future model
// construction; every existing method owns its immutable settings snapshot.
class CoarseCorrectedSliceFactory : public SamplerFactory {
  mutable std::mutex _configuration_mutex;
  CoarseCorrectedSliceConfig _configuration;
public:
  explicit CoarseCorrectedSliceFactory(
    SelNormCoarseSettings const &settings = SelNormCoarseSettings());
  void configure(bool enabled, SelNormCoarseSettings const &settings,
                 double min_acceptance = 0.0);
  CoarseCorrectedSliceConfig configuration() const;
  std::vector<Sampler*> makeSamplers(std::list<StochasticNode*> const &nodes,
                                     Graph const &graph) const;
  std::string name() const;
};

// Borrowed process-lifetime factory. Module registration must not delete it.
CoarseCorrectedSliceFactory *coarse_corrected_slice_factory();
// Counters cover all chains/samplers in this process; reset only between updates.
CoarseCorrectedSliceStats coarse_corrected_slice_stats(bool reset = false);
// Per-fit counters of this process: read, then cleared when clear is true.
CoarseCorrectedSliceReport coarse_corrected_slice_report(bool clear = false);

}
}
#endif
