#include "../fp-contract-off.h"
#include "CorrectedSampler.h"

#include <module/ModuleError.h>

#include <cmath>
#include <cstdio>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <string>
#include <vector>

namespace jags {
namespace RoBMA {

CorrectedSamplerCounters corrected_sampler_counters;

namespace {

// The worst-sampler rate counts samplers with at least min_rate_proposals
// evidence proposals since the last clearing read of the report.
const std::uint64_t min_rate_proposals = 1000;

// The rate of one sampler's window, judged as rate^exponent, and whether it
// counts for the worst rate.
bool window_rate(std::uint64_t proposed, std::uint64_t accepted, double exponent, double *rate)
{
  if (proposed < min_rate_proposals) return false;
  *rate = static_cast<double>(accepted) / static_cast<double>(proposed);
  if (exponent != 1.0) *rate = std::pow(*rate, exponent);
  return true;
}

bool window_rate(CorrectedSamplerWindow const *sampler, double *rate)
{
  return window_rate(sampler->windowProposed(), sampler->windowAccepted(),
                     sampler->rateExponent(), rate);
}

// The judged acceptances of one sampler's window: its proposals times its rate
// to its exponent (rate^(1/d) for a block of d directions), so that the pooled
// judged rate is a proposal-weighted mean of per-direction rates; for a scalar
// sampler, its accepted proposals.
double judged_accepted(CorrectedSamplerWindow const *sampler)
{
  const std::uint64_t proposed = sampler->windowProposed();
  const std::uint64_t accepted = sampler->windowAccepted();
  const double exponent = sampler->rateExponent();
  if (proposed == 0) return 0.0;
  if (exponent == 1.0) return static_cast<double>(accepted);
  return static_cast<double>(proposed) *
    std::pow(static_cast<double>(accepted) / static_cast<double>(proposed), exponent);
}

struct ConfigurationStore {
  std::mutex mutex;
  CorrectedSamplerConfig configuration;
  CorrectedTestFlags test_flags;
};

ConfigurationStore &configuration_store()
{
  static ConfigurationStore store;
  return store;
}

// The live corrected samplers of this process, and the lowest window rate of
// those destroyed since the last clearing read of the report. Registration,
// destruction and report reads take the mutex. Each chain updates its own
// sampler's window, and the R caller reads the windows between JAGS runs.
struct SamplerRegistry {
  std::mutex mutex;
  std::vector<CorrectedSamplerWindow *> live;
  bool has_retired = false;
  double retired_min_rate = 0.0;
  // the judged acceptances of the samplers destroyed since the last clearing read
  double retired_judged = 0.0;
};

// Allocated once and never destroyed: samplers of models that R frees late
// (finalizers at unload or exit) unregister after static destructors may have
// run, so the registry must outlive every sampler.
SamplerRegistry &sampler_registry()
{
  static SamplerRegistry *registry = new SamplerRegistry();
  return *registry;
}

// test flags read during sampling (lock-free)
std::atomic<bool> test_force_fallback(false);
std::atomic<bool> test_freeze_finite(false);

CorrectedKindCounters *kind_counters()
{
  static CorrectedKindCounters *counters = new CorrectedKindCounters[corrected_kind_count];
  return counters;
}

}

namespace {

// The largest sum of absolute likelihood terms per guard site after
// adaptation (relaxed atomics; never destroyed, like the counters).
std::atomic<double> *likelihood_maxima()
{
  static std::atomic<double> *maxima = new std::atomic<double>[likelihood_site_count]();
  return maxima;
}

std::string format_number(double x)
{
  char buffer[32];
  std::snprintf(buffer, sizeof(buffer), "%.6g", x);
  return buffer;
}

}

void LikelihoodMagnitude::add(Node const *child, double term)
{
  const double size = std::fabs(term);
  if (!(size <= std::numeric_limits<double>::max())) return;
  sum += size;
  if (node == nullptr || size > largest) {
    largest = size;
    value = term;
    node = child;
  }
}

void check_likelihood_magnitude(LikelihoodMagnitude const &magnitude, bool adapting, int site,
                                std::string const &sampler)
{
  if (adapting) return;
  std::atomic<double> &maximum = likelihood_maxima()[site];
  double current = maximum.load(std::memory_order_relaxed);
  while (magnitude.sum > current &&
         !maximum.compare_exchange_weak(current, magnitude.sum, std::memory_order_relaxed)) {}
  if (!(magnitude.sum > likelihood_representable_limit) || magnitude.node == nullptr) return;
  throwNodeError(magnitude.node, sampler + " [RoBMA_likelihood_not_representable]: the log density of this "
                 "likelihood node is " + format_number(magnitude.value) + " at the chain's current state, and the "
                 "absolute log densities of the likelihood terms evaluated with it sum to " +
                 format_number(magnitude.sum) + ", above " + format_number(likelihood_representable_limit) +
                 ", where double precision cannot resolve the differences the sampler needs; check the data, "
                 "their units and the model");
}

std::vector<double> likelihood_magnitude_maxima(bool clear)
{
  std::vector<double> out(likelihood_site_count);
  for (int i = 0; i < likelihood_site_count; ++i) {
    out[static_cast<std::size_t>(i)] = clear ? likelihood_maxima()[i].exchange(0.0, std::memory_order_relaxed) :
      likelihood_maxima()[i].load(std::memory_order_relaxed);
  }
  return out;
}

CorrectedKindCounters &corrected_kind_counters(int kind)
{
  if (kind < 0 || kind >= corrected_kind_count) throw std::logic_error("Unknown corrected sampler kind.");
  return kind_counters()[kind];
}

CorrectedTestFlags corrected_test_flags()
{
  ConfigurationStore &store = configuration_store();
  const std::lock_guard<std::mutex> lock(store.mutex);
  return store.test_flags;
}

void corrected_test_set_flags(CorrectedTestFlags const &flags)
{
  if (flags.finite_proposal < 0 || flags.finite_proposal > 3) {
    throw std::invalid_argument("Unknown screened indicator proposal.");
  }
  ConfigurationStore &store = configuration_store();
  const std::lock_guard<std::mutex> lock(store.mutex);
  store.test_flags = flags;
  test_force_fallback.store(flags.force_fallback, std::memory_order_relaxed);
  test_freeze_finite.store(flags.freeze_finite, std::memory_order_relaxed);
}

bool corrected_test_force_fallback()
{
  return test_force_fallback.load(std::memory_order_relaxed);
}

bool corrected_test_freeze_finite()
{
  return test_freeze_finite.load(std::memory_order_relaxed);
}

void CorrectedTally::count(bool accepted, bool evidence)
{
  CorrectedSamplerCounters &counters = corrected_sampler_counters;
  counters.coarse_proposed.fetch_add(1, std::memory_order_relaxed);
  (accepted ? counters.coarse_accepted : counters.coarse_correction_rejected)
    .fetch_add(1, std::memory_order_relaxed);
  CorrectedKindCounters &kind = corrected_kind_counters(_kind);
  kind.proposed.fetch_add(1, std::memory_order_relaxed);
  if (accepted) kind.accepted.fetch_add(1, std::memory_order_relaxed);
  if (!evidence) return;
  kind.evidence_proposed.fetch_add(1, std::memory_order_relaxed);
  if (accepted) kind.evidence_accepted.fetch_add(1, std::memory_order_relaxed);
  if (_adapting) {
    ++_adaptive_evidence_proposed;
    if (accepted) ++_adaptive_evidence_accepted;
    return;
  }
  ++_window_proposed;
  counters.report_sampling_proposed.fetch_add(1, std::memory_order_relaxed);
  if (accepted) {
    ++_window_accepted;
    counters.report_sampling_accepted.fetch_add(1, std::memory_order_relaxed);
  }
}

bool CorrectedTally::adapt_off(double min_acceptance, double exponent)
{
  if (!_adapting) return false;
  _adapting = false;
  if (!corrected_test_force_fallback()) {
    if (!(min_acceptance > 0.0) || _adaptive_evidence_proposed < corrected_min_evidence_proposals) {
      return false;
    }
    double rate = static_cast<double>(_adaptive_evidence_accepted) /
      static_cast<double>(_adaptive_evidence_proposed);
    if (exponent != 1.0) rate = std::pow(rate, exponent);
    if (!(rate < min_acceptance)) return false;
  }
  corrected_sampler_counters.report_fallback_low_acceptance.fetch_add(1, std::memory_order_relaxed);
  corrected_kind_counters(_kind).fallbacks.fetch_add(1, std::memory_order_relaxed);
  return true;
}

CorrectedSamplerConfig corrected_sampler_configuration()
{
  ConfigurationStore &store = configuration_store();
  const std::lock_guard<std::mutex> lock(store.mutex);
  return store.configuration;
}

void configure_corrected_samplers(CorrectedSamplerConfig const &config)
{
  SelNormCoarseSettings const &settings = config.settings;
  const double min_acceptance = config.min_acceptance;
  if (!(settings.mean_step >= 0) || !std::isfinite(settings.mean_step) ||
      !(settings.diagonal_step >= 0) || !std::isfinite(settings.diagonal_step) ||
      !(settings.log_weight_step >= 0) || !std::isfinite(settings.log_weight_step) ||
      settings.max_rules < 3 || !(min_acceptance >= 0.0 && min_acceptance <= 1.0)) {
    throw std::invalid_argument("Invalid corrected selection surrogate settings.");
  }
  ConfigurationStore &store = configuration_store();
  const std::lock_guard<std::mutex> lock(store.mutex);
  store.configuration = config;
}

void register_corrected_sampler(CorrectedSamplerWindow *sampler, int kind)
{
  SamplerRegistry &registry = sampler_registry();
  const std::lock_guard<std::mutex> lock(registry.mutex);
  registry.live.push_back(sampler);
  corrected_sampler_counters.report_samplers.fetch_add(1, std::memory_order_relaxed);
  corrected_kind_counters(kind).samplers.fetch_add(1, std::memory_order_relaxed);
}

void unregister_corrected_sampler(CorrectedSamplerWindow *sampler)
{
  SamplerRegistry &registry = sampler_registry();
  const std::lock_guard<std::mutex> lock(registry.mutex);
  std::vector<CorrectedSamplerWindow *> &live = registry.live;
  for (std::size_t i = 0; i < live.size(); ++i) {
    if (live[i] == sampler) {
      live[i] = live.back();
      live.pop_back();
      break;
    }
  }
  double rate;
  if (window_rate(sampler, &rate) &&
      (!registry.has_retired || rate < registry.retired_min_rate)) {
    registry.retired_min_rate = rate;
    registry.has_retired = true;
  }
  registry.retired_judged += judged_accepted(sampler);
}

CorrectedSamplerStats corrected_sampler_stats(bool reset)
{
  CorrectedSamplerCounters &counters = corrected_sampler_counters;
  if (reset) {
    return {
      counters.coarse_proposed.exchange(0, std::memory_order_relaxed),
      counters.coarse_accepted.exchange(0, std::memory_order_relaxed),
      counters.coarse_correction_rejected.exchange(0, std::memory_order_relaxed)
    };
  }
  return {
    counters.coarse_proposed.load(std::memory_order_relaxed),
    counters.coarse_accepted.load(std::memory_order_relaxed),
    counters.coarse_correction_rejected.load(std::memory_order_relaxed)
  };
}

CorrectedSamplerReport corrected_sampler_report(bool clear)
{
  // Each counter is read, and cleared in the same operation when requested.
  const auto take = [clear](std::atomic<std::uint64_t> &counter) {
    return clear ? counter.exchange(0, std::memory_order_relaxed) :
      counter.load(std::memory_order_relaxed);
  };
  CorrectedSamplerCounters &counters = corrected_sampler_counters;
  CorrectedSamplerReport report;
  report.samplers = take(counters.report_samplers);
  report.fallback_low_acceptance = take(counters.report_fallback_low_acceptance);
  report.sampling_proposed = take(counters.report_sampling_proposed);
  report.sampling_accepted = take(counters.report_sampling_accepted);
  // The lowest rate of one sampler over its window: the destroyed samplers'
  // minimum and the live samplers with enough proposals (fallen-back samplers
  // have none).
  SamplerRegistry &registry = sampler_registry();
  const std::lock_guard<std::mutex> lock(registry.mutex);
  report.has_min_rate = registry.has_retired;
  report.min_sampler_rate = registry.has_retired ? registry.retired_min_rate : 0.0;
  report.sampling_judged = registry.retired_judged;
  for (CorrectedSamplerWindow *sampler : registry.live) {
    report.sampling_judged += judged_accepted(sampler);
    double rate;
    if (window_rate(sampler, &rate) &&
        (!report.has_min_rate || rate < report.min_sampler_rate)) {
      report.min_sampler_rate = rate;
      report.has_min_rate = true;
    }
    if (clear) sampler->clearWindow();
  }
  if (clear) {
    registry.has_retired = false;
    registry.retired_min_rate = 0.0;
    registry.retired_judged = 0.0;
  }
  return report;
}

}
}
