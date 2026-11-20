#ifndef ROBMA_CORRECTED_SAMPLER_H
#define ROBMA_CORRECTED_SAMPLER_H

// What all of RoBMA's corrected samplers share: the frozen configuration their
// factories read, the process counters of selection_sampler_info() and of the
// per-fit report, and the registry of live samplers whose reporting windows the
// report reads (and clears).

#include <atomic>
#include <cstdint>
#include "../selnorm/selnorm-mv.h"

namespace jags {
namespace RoBMA {

// min_acceptance > 0 enables the fallback to the exact kernel that each sampler
// decides once, at the end of adaptation; 0 disables it.
struct CorrectedSamplerConfig {
  bool enabled = false;
  SelNormCoarseSettings settings;
  double min_acceptance = 0.0;
};

// The configuration starts disabled. It controls only future model
// construction; every existing sampler owns its immutable settings snapshot.
CorrectedSamplerConfig corrected_sampler_configuration();
void configure_corrected_samplers(CorrectedSamplerConfig const &config);

struct CorrectedSamplerStats {
  std::uint64_t proposed;
  std::uint64_t accepted;
  std::uint64_t correction_rejected;
};

// The per-fit counters, read and cleared by the runtime report of each fit:
// sampler objects constructed (sampled node or block x chain, per model build),
// the fallbacks decided at adaptOff() by reason, the proposals and accepted
// corrections made after adaptOff() by samplers that stayed corrected, and the
// lowest acceptance rate of one such sampler over its proposals since the last
// clearing read, judged as rate^rateExponent() (live samplers and those
// destroyed since; only samplers with at least 1,000 such proposals count,
// has_min_rate false when none does).
struct CorrectedSamplerReport {
  std::uint64_t samplers;
  std::uint64_t fallback_few_proposals;
  std::uint64_t fallback_surrogate_exact;
  std::uint64_t fallback_low_acceptance;
  std::uint64_t sampling_proposed;
  std::uint64_t sampling_accepted;
  // sum over the samplers of their proposals x window rate^rateExponent()
  // (the accepted proposals of a scalar sampler)
  double sampling_judged;
  bool has_min_rate;
  double min_sampler_rate;
};

// The process counters the samplers update (relaxed atomics): the three of
// selection_sampler_info() and the per-fit ones of CorrectedSamplerReport.
struct CorrectedSamplerCounters {
  std::atomic<std::uint64_t> coarse_proposed{0};
  std::atomic<std::uint64_t> coarse_accepted{0};
  std::atomic<std::uint64_t> coarse_correction_rejected{0};
  std::atomic<std::uint64_t> report_samplers{0};
  std::atomic<std::uint64_t> report_fallback_few_proposals{0};
  std::atomic<std::uint64_t> report_fallback_surrogate_exact{0};
  std::atomic<std::uint64_t> report_fallback_low_acceptance{0};
  std::atomic<std::uint64_t> report_sampling_proposed{0};
  std::atomic<std::uint64_t> report_sampling_accepted{0};
};
extern CorrectedSamplerCounters corrected_sampler_counters;

// Fallback rules decided at adaptOff(): fewer adaptive proposals than this
// cannot be judged.
const std::uint64_t corrected_min_adaptive_proposals = 50;

// A corrected sampler's reporting window: its proposals and accepted
// corrections after adaptOff() in corrected mode since the last clearing read
// of the report (the window of one run), and the exponent by which its
// acceptance rate is judged: 1 for a scalar, 1/d for a block of d members
// (rate^(1/d), the per-direction equivalent of a move that composes d
// direction moves).
class CorrectedSamplerWindow {
public:
  virtual ~CorrectedSamplerWindow() = default;
  virtual std::uint64_t windowProposed() const = 0;
  virtual std::uint64_t windowAccepted() const = 0;
  virtual void clearWindow() = 0;
  virtual double rateExponent() const = 0;
};

// Registration counts one sampler object of the report. A sampler unregisters
// in its own destructor (its window is still readable there), and a sampler
// destroyed before the report is read still counts for the worst rate of its run.
void register_corrected_sampler(CorrectedSamplerWindow *sampler);
void unregister_corrected_sampler(CorrectedSamplerWindow *sampler);

// Counters cover all chains/samplers in this process; reset only between updates.
CorrectedSamplerStats corrected_sampler_stats(bool reset = false);
// Per-fit counters of this process: read, then cleared when clear is true.
CorrectedSamplerReport corrected_sampler_report(bool clear = false);

}
}
#endif
