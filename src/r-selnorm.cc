#include <Rinternals.h>
#include <R_ext/Error.h>
#include <R_ext/Lapack.h>
#include <R_ext/Random.h>

#include <Rmath.h>

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "plot-root.h"
#include "selnorm/selnorm.h"
#include "selnorm/selnorm-mv.h"
#include "selnorm/selnorm-parallel.h"
#include "samplers/CoarseCorrectedSlice.h"

#ifndef FCONE
# define FCONE
#endif

extern "C" double Rf_dnorm4(double, double, double, int);

// Private implementation fragments share one anonymous namespace and are kept in
// this order so helper symbols stay internal to this translation unit.
namespace {
#include "r-selnorm-common.cc.inc"
#include "r-selnorm-funnel-common.cc.inc"
#include "r-selnorm-cluster.cc.inc"
}

#include "r-selnorm-loglik.cc.inc"
#include "r-selnorm-mv.cc.inc"
#include "r-selnorm-sampling-conditioned.cc.inc"
#include "r-selnorm-kernel.cc.inc"
#include "r-selnorm-funnel-zcurve.cc.inc"

extern "C" SEXP RoBMA_selnorm_set_native_threads(SEXP threads)
{
  if ((TYPEOF(threads) != INTSXP && TYPEOF(threads) != REALSXP) ||
      XLENGTH(threads) != 1 || Rf_inherits(threads, "factor")) {
    Rf_error("'threads' must be one positive integer or zero for the OpenMP default.");
  }
  const double value = Rf_asReal(threads);
  if (!std::isfinite(value) || value < 0.0 || value != std::floor(value) ||
      value > static_cast<double>(std::numeric_limits<int>::max())) {
    Rf_error("'threads' must be one nonnegative whole number of threads.");
  }
  robma_set_native_threads(static_cast<int>(value));
  return R_NilValue;
}

extern "C" SEXP RoBMA_selnorm_cache_snapshot()
{
  for (int attempt = 0; attempt < 2; ++attempt) {
    std::size_t size = 0; char error_message[512] = {};
    try { size = cpp_selnorm_cache_snapshot_size(); }
    catch (const std::exception &error) { std::strncpy(error_message, error.what(), sizeof(error_message) - 1); }
    if (error_message[0]) Rf_error("%s", error_message);
    if (size > static_cast<std::size_t>(R_XLEN_T_MAX)) Rf_error("Selection cache snapshot is too large for R.");
    SEXP output = PROTECT(Rf_allocVector(RAWSXP, static_cast<R_xlen_t>(size)));
    bool complete = false;
    try { complete = cpp_selnorm_cache_snapshot_write(RAW(output), size); }
    catch (const std::exception &error) { std::strncpy(error_message, error.what(), sizeof(error_message) - 1); }
    UNPROTECT(1);
    if (error_message[0]) Rf_error("%s", error_message);
    if (complete) return output;
  }
  Rf_error("Selection cache changed during snapshot capture.");
  return R_NilValue;
}

extern "C" SEXP RoBMA_selnorm_cache_restore(SEXP snapshots)
{
  if (TYPEOF(snapshots) != VECSXP) Rf_error("Selection cache snapshots must be a list of raw vectors.");
  const R_xlen_t count = XLENGTH(snapshots);
  if (static_cast<std::uint64_t>(count) > std::numeric_limits<std::size_t>::max() / sizeof(SelNormCacheBlob))
    Rf_error("Selection cache snapshot list is too large.");
  for (R_xlen_t i = 0; i < count; ++i)
    if (TYPEOF(VECTOR_ELT(snapshots, i)) != RAWSXP) Rf_error("Selection cache snapshots must be a list of raw vectors.");
  // R owns this small view array; no C++ object can be skipped by an R error.
  auto *blobs = reinterpret_cast<SelNormCacheBlob *>(R_alloc(static_cast<std::size_t>(count), sizeof(SelNormCacheBlob)));
  for (R_xlen_t i = 0; i < count; ++i) {
    SEXP value = VECTOR_ELT(snapshots, i);
    if (static_cast<std::uint64_t>(XLENGTH(value)) > std::numeric_limits<std::size_t>::max())
      Rf_error("Selection cache snapshot is too large.");
    blobs[i] = {RAW(value), static_cast<std::size_t>(XLENGTH(value))};
  }
  SelNormCacheRestoreInfo info; char error_message[512] = {};
  try { info = cpp_selnorm_cache_restore(blobs, static_cast<std::size_t>(count)); }
  catch (const std::exception &error) { std::strncpy(error_message, error.what(), sizeof(error_message) - 1); }
  if (error_message[0]) Rf_error("%s", error_message);
  const char *names[] = {"snapshots", "available_exact", "available_coarse", "restored_exact", "restored_coarse", "skipped"};
  const double values[] = {static_cast<double>(info.snapshots), static_cast<double>(info.available_exact),
    static_cast<double>(info.available_coarse), static_cast<double>(info.restored_exact),
    static_cast<double>(info.restored_coarse), static_cast<double>(info.skipped)};
  SEXP output = PROTECT(Rf_allocVector(VECSXP, 6));
  SEXP labels = PROTECT(Rf_allocVector(STRSXP, 6));
  for (int i = 0; i < 6; ++i) {
    SET_VECTOR_ELT(output, i, Rf_ScalarReal(values[i]));
    SET_STRING_ELT(labels, i, Rf_mkChar(names[i]));
  }
  Rf_setAttrib(output, R_NamesSymbol, labels);
  UNPROTECT(2);
  return output;
}

extern "C" SEXP RoBMA_selnorm_cache_control(SEXP capacity_bytes, SEXP clear_mask)
{
  if (TYPEOF(clear_mask) != INTSXP || XLENGTH(clear_mask) != 1 ||
      INTEGER(clear_mask)[0] < 0 || INTEGER(clear_mask)[0] > 3) {
    Rf_error("'clear_mask' must be one integer from 0 to 3.");
  }
  const bool set_capacity = capacity_bytes != R_NilValue;
  std::size_t capacity = 0;
  if (set_capacity) {
    if ((TYPEOF(capacity_bytes) != REALSXP && TYPEOF(capacity_bytes) != INTSXP) ||
        XLENGTH(capacity_bytes) != 1 || Rf_inherits(capacity_bytes, "factor")) {
      Rf_error("'capacity_bytes' must be NULL or one non-negative whole number of bytes.");
    }
    const double value = Rf_asReal(capacity_bytes);
    const long double size_limit = std::ldexp(1.0L, std::numeric_limits<std::size_t>::digits);
    if (!std::isfinite(value) || value < 0.0 || value != std::floor(value) ||
        static_cast<long double>(value) >= size_limit) {
      Rf_error("'capacity_bytes' must be a finite non-negative whole number of bytes supported on this platform.");
    }
    capacity = static_cast<std::size_t>(value);
  }
  SelNormCacheInfo info;
  char error_message[512] = {};
  try {
    info = cpp_selnorm_cache_control(set_capacity, capacity,
      static_cast<unsigned int>(INTEGER(clear_mask)[0]));
  } catch (const std::exception &error) {
    std::strncpy(error_message, error.what(), sizeof(error_message) - 1);
  }
  // Cache locks are released before any R error or allocation.
  if (error_message[0]) Rf_error("Selection cache operation failed: %s", error_message);
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 5));
  SEXP labels = PROTECT(Rf_allocVector(STRSXP, 5));
  const char *names[] = {"capacity_bytes", "allocated_bytes", "peak_bytes", "exact", "coarse"};
  for (int i = 0; i < 5; ++i) SET_STRING_ELT(labels, i, Rf_mkChar(names[i]));
  SET_VECTOR_ELT(out, 0, Rf_ScalarReal(static_cast<double>(info.capacity_bytes)));
  SET_VECTOR_ELT(out, 1, Rf_ScalarReal(static_cast<double>(info.allocated_bytes)));
  SET_VECTOR_ELT(out, 2, Rf_ScalarReal(static_cast<double>(info.peak_bytes)));
  const char *stat_names[] = {"entries", "hits", "misses", "evictions", "resets", "allocation_failures"};
  for (int component = 0; component < 2; ++component) {
    const SelNormCacheStats &stats = component == 0 ? info.exact : info.coarse;
    const double values[] = {static_cast<double>(stats.entries), static_cast<double>(stats.hits),
      static_cast<double>(stats.misses), static_cast<double>(stats.evictions),
      static_cast<double>(stats.resets), static_cast<double>(stats.allocation_failures)};
    SEXP part = PROTECT(Rf_allocVector(VECSXP, 6));
    SEXP part_names = PROTECT(Rf_allocVector(STRSXP, 6));
    for (int i = 0; i < 6; ++i) {
      SET_VECTOR_ELT(part, i, Rf_ScalarReal(values[i]));
      SET_STRING_ELT(part_names, i, Rf_mkChar(stat_names[i]));
    }
    Rf_setAttrib(part, R_NamesSymbol, part_names);
    SET_VECTOR_ELT(out, 3 + component, part);
    UNPROTECT(2);
  }
  Rf_setAttrib(out, R_NamesSymbol, labels);
  UNPROTECT(2);
  return out;
}

extern "C" SEXP RoBMA_selnorm_sampler_control(SEXP enabled, SEXP settings, SEXP clear)
{
  if (TYPEOF(clear) != LGLSXP || XLENGTH(clear) != 1 || LOGICAL(clear)[0] == NA_LOGICAL) {
    Rf_error("'clear' must be one non-missing logical value.");
  }
  if (enabled != R_NilValue && (TYPEOF(enabled) != LGLSXP || XLENGTH(enabled) != 1 ||
      LOGICAL(enabled)[0] == NA_LOGICAL)) {
    Rf_error("'enabled' must be NULL or one non-missing logical value.");
  }
  if (settings != R_NilValue) {
    if (TYPEOF(settings) != REALSXP || XLENGTH(settings) != 4) {
      Rf_error("'settings' must contain three grid steps and a rule count.");
    }
    const double *values = REAL(settings);
    for (int i = 0; i < 3; ++i) {
      if (!std::isfinite(values[i]) || values[i] < 0) {
        Rf_error("Coarse grid steps must be finite and nonnegative.");
      }
    }
    if (!std::isfinite(values[3]) || values[3] < 3 || values[3] != std::floor(values[3]) ||
        values[3] > std::numeric_limits<unsigned int>::max()) {
      Rf_error("The coarse rule count must be a supported integer of at least three.");
    }
  }
  jags::RoBMA::CoarseCorrectedSliceConfig config;
  jags::RoBMA::CoarseCorrectedSliceStats stats;
  char error_message[512] = {};
  try {
    auto *factory = jags::RoBMA::coarse_corrected_slice_factory();
    config = factory->configuration();
    if (enabled != R_NilValue) config.enabled = LOGICAL(enabled)[0] != 0;
    if (settings != R_NilValue) {
      config.settings.mean_step = REAL(settings)[0];
      config.settings.diagonal_step = REAL(settings)[1];
      config.settings.log_weight_step = REAL(settings)[2];
      config.settings.max_rules = static_cast<unsigned int>(REAL(settings)[3]);
    }
    if (enabled != R_NilValue || settings != R_NilValue) {
      factory->configure(config.enabled, config.settings);
    }
    if (LOGICAL(clear)[0]) jags::RoBMA::coarse_corrected_slice_stats(true);
    stats = jags::RoBMA::coarse_corrected_slice_stats();
  } catch (const std::exception &error) {
    std::strncpy(error_message, error.what(), sizeof(error_message) - 1);
  }
  if (error_message[0]) Rf_error("Selection sampler operation failed: %s", error_message);
  const char *names[] = {"enabled", "proposed", "accepted", "correction_rejected"};
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 4));
  SEXP labels = PROTECT(Rf_allocVector(STRSXP, 4));
  SET_VECTOR_ELT(out, 0, Rf_ScalarLogical(config.enabled));
  SET_VECTOR_ELT(out, 1, Rf_ScalarReal(static_cast<double>(stats.proposed)));
  SET_VECTOR_ELT(out, 2, Rf_ScalarReal(static_cast<double>(stats.accepted)));
  SET_VECTOR_ELT(out, 3, Rf_ScalarReal(static_cast<double>(stats.correction_rejected)));
  for (int i = 0; i < 4; ++i) SET_STRING_ELT(labels, i, Rf_mkChar(names[i]));
  Rf_setAttrib(out, R_NamesSymbol, labels);
  UNPROTECT(2);
  return out;
}
