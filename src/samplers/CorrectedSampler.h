#ifndef ROBMA_CORRECTED_SAMPLER_H
#define ROBMA_CORRECTED_SAMPLER_H

// What all of RoBMA's corrected samplers share: the frozen configuration their
// factories read, the process counters of selection_sampler_info() and of the
// per-fit report, the registry of live samplers whose reporting windows the
// report reads (and clears), the evidence bookkeeping that decides each
// sampler's kernel at the end of adaptation, and the test hooks.

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
// the fallbacks decided at adaptOff(), the evidence proposals made after
// adaptOff() by samplers that stayed corrected and their accepted corrections,
// and the lowest acceptance rate of one such sampler over its evidence
// proposals since the last clearing read, judged as rate^rateExponent() (live
// samplers and those destroyed since; only samplers with at least 1,000 such
// proposals count, has_min_rate false when none does).
struct CorrectedSamplerReport {
  std::uint64_t samplers;
  std::uint64_t fallback_low_acceptance;
  std::uint64_t sampling_proposed;
  std::uint64_t sampling_accepted;
  // sum over the samplers of their evidence proposals x window
  // rate^rateExponent() (the accepted proposals of a scalar sampler)
  double sampling_judged;
  bool has_min_rate;
  double min_sampler_rate;
};

// The process counters the samplers update (relaxed atomics): the three of
// selection_sampler_info() (every completed proposal, adaptation included) and
// the per-fit ones of CorrectedSamplerReport.
struct CorrectedSamplerCounters {
  std::atomic<std::uint64_t> coarse_proposed{0};
  std::atomic<std::uint64_t> coarse_accepted{0};
  std::atomic<std::uint64_t> coarse_correction_rejected{0};
  std::atomic<std::uint64_t> report_samplers{0};
  std::atomic<std::uint64_t> report_fallback_low_acceptance{0};
  std::atomic<std::uint64_t> report_sampling_proposed{0};
  std::atomic<std::uint64_t> report_sampling_accepted{0};
};
extern CorrectedSamplerCounters corrected_sampler_counters;

// Evidence is a proposal made while at least one selection child evaluated a
// surrogate that is not its exact density (a coarse rung or an anchor
// normalizer under the selection kernel). A sampler falls back to its exact
// kernel at its first adaptOff() only when it made at least this many evidence
// proposals while adapting and their acceptance (rate^(1/d) for a block) is
// below the threshold.
const std::uint64_t corrected_min_evidence_proposals = 50;

// The kinds of corrected samplers, for the test hook's counters.
enum CorrectedKind { corrected_kind_slice = 0, corrected_kind_block = 1, corrected_kind_gated = 2,
                     corrected_kind_finite = 3, corrected_kind_location = 4, corrected_kind_count = 5 };

// Test hook counters per kind (process-wide, relaxed atomics, cleared by the
// hook): sampler objects, updates (every call, including gate-off draws and
// self-proposals of an indicator), proposals and accepted corrections,
// evidence proposals and their accepted corrections (adaptation included),
// fallbacks decided at adaptOff(), and exact evaluations of a whole full
// conditional (ScreenedFinite only).
struct CorrectedKindCounters {
  std::atomic<std::uint64_t> samplers{0};
  std::atomic<std::uint64_t> updates{0};
  std::atomic<std::uint64_t> proposed{0};
  std::atomic<std::uint64_t> accepted{0};
  std::atomic<std::uint64_t> evidence_proposed{0};
  std::atomic<std::uint64_t> evidence_accepted{0};
  std::atomic<std::uint64_t> fallbacks{0};
  std::atomic<std::uint64_t> exact_scores{0};
};
CorrectedKindCounters &corrected_kind_counters(int kind);

// Test flags (process-wide, read by the factories at model construction):
// the samplers of this round (gated weight blocks and screened indicators)
// left out, the gamma-sum Jacobian dropped (a kernel that is not exact), and
// the screened indicator's proposal: 2 (the default) Liu's Metropolized Gibbs on
// the surrogate with the exact correction, 0 the plain surrogate full
// conditional with the exact correction, 1 the same without the correction
// (not exact).
// Read at adaptOff() and update() instead: force_fallback switches every
// sampler that ends its adaptation to its exact kernel (counted as a
// low-acceptance fallback), and freeze_finite makes the screened indicators
// keep their values (no update).
struct CorrectedTestFlags {
  bool disable_new = false;
  // the same for one of them: gated weight blocks, screened indicators, gamma groups
  bool disable_gated = false;
  bool disable_finite = false;
  bool disable_sum = false;
  bool drop_sum_jacobian = false;
  int finite_proposal = 2;
  bool force_fallback = false;
  bool freeze_finite = false;
};
bool corrected_test_force_fallback();
bool corrected_test_freeze_finite();
CorrectedTestFlags corrected_test_flags();
void corrected_test_set_flags(CorrectedTestFlags const &flags);

// One completed proposal's bookkeeping and the kernel decision at the first
// adaptOff(), shared by the corrected samplers. `exponent` judges the rate
// (1/d for a block of d directions).
class CorrectedTally {
  int _kind;
  bool _adapting = true;
  std::uint64_t _adaptive_evidence_proposed = 0;
  std::uint64_t _adaptive_evidence_accepted = 0;
  std::uint64_t _window_proposed = 0;
  std::uint64_t _window_accepted = 0;
public:
  explicit CorrectedTally(int kind) : _kind(kind) {}
  bool adapting() const { return _adapting; }
  // One completed proposal: counted in the process counters and, when it is
  // evidence, in the adaptive tally or (after adaptation) in the reporting
  // window and the per-fit counters.
  void count(bool accepted, bool evidence);
  // The first call ends adaptation and returns true when the sampler is to
  // switch to its exact kernel (later calls return false and change nothing).
  bool adapt_off(double min_acceptance, double exponent);
  std::uint64_t window_proposed() const { return _window_proposed; }
  std::uint64_t window_accepted() const { return _window_accepted; }
  void clear_window() { _window_proposed = 0; _window_accepted = 0; }
  std::uint64_t adaptive_evidence_proposed() const { return _adaptive_evidence_proposed; }
};

// A corrected sampler's reporting window: its evidence proposals and accepted
// corrections after adaptOff() in corrected mode since the last clearing read
// of the report (the window of one run), and the exponent by which its
// acceptance rate is judged: 1 for a scalar, 1/d for a block of d directions
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
void register_corrected_sampler(CorrectedSamplerWindow *sampler, int kind);
void unregister_corrected_sampler(CorrectedSamplerWindow *sampler);

// Counters cover all chains/samplers in this process; reset only between updates.
CorrectedSamplerStats corrected_sampler_stats(bool reset = false);
// Per-fit counters of this process: read, then cleared when clear is true.
CorrectedSamplerReport corrected_sampler_report(bool clear = false);

}
}
#endif
