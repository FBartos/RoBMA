#include "CorrectedSampler.h"

#include <cmath>
#include <mutex>
#include <stdexcept>
#include <vector>

namespace jags {
namespace RoBMA {

CorrectedSamplerCounters corrected_sampler_counters;

namespace {

// The worst-sampler rate counts samplers with at least min_rate_proposals
// proposals since the last clearing read of the report.
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
// to its exponent (rate^(1/d) for a block of d members), so that the pooled
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
};

ConfigurationStore &configuration_store()
{
  static ConfigurationStore store;
  return store;
}

// The live corrected samplers of this process, and the lowest window rate of
// those destroyed since the last clearing read of the report. Registration,
// destruction and the report's read take the mutex; window counters change only
// inside updates, which JAGS runs on the thread that also reads the report
// between runs.
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

void register_corrected_sampler(CorrectedSamplerWindow *sampler)
{
  SamplerRegistry &registry = sampler_registry();
  const std::lock_guard<std::mutex> lock(registry.mutex);
  registry.live.push_back(sampler);
  corrected_sampler_counters.report_samplers.fetch_add(1, std::memory_order_relaxed);
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
  report.fallback_few_proposals = take(counters.report_fallback_few_proposals);
  report.fallback_surrogate_exact = take(counters.report_fallback_surrogate_exact);
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
