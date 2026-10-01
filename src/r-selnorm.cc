#include <Rinternals.h>
#include <R_ext/Error.h>
#include <R_ext/Lapack.h>
#include <R_ext/Random.h>

#include <Rmath.h>

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <chrono>
#include <cstring>
#include <limits>
#include <memory>
#include <mutex>
#include <new>
#include <stdexcept>
#include <string>
#include <vector>

#include "plot-root.h"
#include "selnorm/selnorm.h"
#include "selnorm/selnorm-tail.h"
#include "selnorm/selnorm-mv.h"
#include "selnorm/selnorm-parallel.h"
#include "samplers/CoarseCorrectedSlice.h"

#ifndef FCONE
# define FCONE
#endif

extern "C" double Rf_dnorm4(double, double, double, int);

#include "r-native-api.h"

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
#include "r-native-probe.cc.inc"

// Standard-normal upper tail as the selection kernels evaluate it, exposed so
// the package tests can certify the batched approximation against pnorm over a
// dense argument grid, and so a length-one request and the same value inside a
// batch can be shown to agree.
extern "C" SEXP RoBMA_selnorm_normal_upper_tail(SEXP x, SEXP scalar)
{
  ROBMA_NATIVE_BEGIN(x, scalar)
  if ((TYPEOF(x) != REALSXP && TYPEOF(x) != INTSXP) || Rf_inherits(x, "factor")) {
    Rf_error("'x' must be a numeric vector.");
  }
  if (TYPEOF(scalar) != LGLSXP || Rf_length(scalar) != 1 ||
      LOGICAL(scalar)[0] == NA_LOGICAL) {
    Rf_error("'scalar' must be one non-missing logical value.");
  }
  SEXP values = PROTECT(Rf_coerceVector(x, REALSXP));
  const R_xlen_t count = XLENGTH(values);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, count));
  if (count > 0) {
    // The kernel's own affine form: score = 0 - x * (-1).
    selnorm_tail::upper_tail_affine(REAL(values),
      static_cast<std::size_t>(count), 0.0, -1.0, REAL(out),
      LOGICAL(scalar)[0] == TRUE);
  }
  UNPROTECT(2);
  return out;
  ROBMA_NATIVE_END
}

// The row schedule the batch kernels resolve for a shape, exposed so the
// package tests can certify which shapes stay serial and how large a parallel
// region is without timing anything.
extern "C" SEXP RoBMA_selnorm_row_schedule(SEXP rows, SEXP work_per_row)
{
  ROBMA_NATIVE_BEGIN(rows, work_per_row)
  if ((TYPEOF(rows) != INTSXP && TYPEOF(rows) != REALSXP) || XLENGTH(rows) != 1 ||
      (TYPEOF(work_per_row) != INTSXP && TYPEOF(work_per_row) != REALSXP) ||
      XLENGTH(work_per_row) != 1) {
    Rf_error("'rows' and 'work_per_row' must be single numbers.");
  }
  const double row_count = Rf_asReal(rows);
  const double work = Rf_asReal(work_per_row);
  if (!std::isfinite(row_count) || row_count < 0 ||
      row_count > static_cast<double>(std::numeric_limits<int>::max())) {
    Rf_error("'rows' must be one nonnegative whole number of rows.");
  }
  const RobmaRowSchedule schedule =
    robma_row_schedule(static_cast<int>(row_count), work);
  SEXP out = PROTECT(Rf_allocVector(INTSXP, 2));
  SEXP names = PROTECT(Rf_allocVector(STRSXP, 2));
  INTEGER(out)[0] = schedule.threads;
  INTEGER(out)[1] = schedule.chunk_rows;
  SET_STRING_ELT(names, 0, Rf_mkChar("threads"));
  SET_STRING_ELT(names, 1, Rf_mkChar("chunk_rows"));
  Rf_setAttrib(out, R_NamesSymbol, names);
  UNPROTECT(2);
  return out;
  ROBMA_NATIVE_END
}

// The rule-weight memo of the batch entry points, exposed so the package tests
// can certify that a remembered rule is served the same values a fresh
// exponentiation of its log weights produces, and that a changed rule is not
// served a remembered one. 'identifier' is the same for two calls that are
// served the same vector and differs whenever a new one was computed; the small
// registry behind it exists only for this certification entry point.
extern "C" SEXP RoBMA_selnorm_rule_weights_check(SEXP log_weights)
{
  ROBMA_NATIVE_BEGIN(log_weights)
  if (TYPEOF(log_weights) != REALSXP || XLENGTH(log_weights) < 1 ||
      XLENGTH(log_weights) > std::numeric_limits<int>::max()) {
    Rf_error("'log_weights' must be a nonempty numeric vector.");
  }
  const int count = static_cast<int>(XLENGTH(log_weights));
  const SelNormRuleWeights weights =
    selnorm_rule_weights(REAL(log_weights), count);
  if (!weights || static_cast<int>(weights->size()) != count) {
    Rf_error("The rule weight memo returned a vector of the wrong length.");
  }

  int mismatches = 0;
  for (int i = 0; i < count; ++i) {
    const long double fresh =
      std::exp(static_cast<long double>(REAL(log_weights)[i]));
    const long double served = (*weights)[static_cast<std::size_t>(i)];
    const bool both_missing = std::isnan(fresh) && std::isnan(served);
    if (!both_missing && !(served == fresh)) ++mismatches;
  }

  static std::vector<std::pair<SelNormRuleWeights, int> > registry;
  static int last_identifier = 0;
  int identifier = NA_INTEGER;
  for (std::size_t i = 0; i < registry.size(); ++i) {
    if (registry[i].first == weights) {
      identifier = registry[i].second;
      break;
    }
  }
  if (identifier == NA_INTEGER && registry.size() < 256) {
    registry.push_back(std::make_pair(weights, ++last_identifier));
    identifier = last_identifier;
  }

  SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
  SEXP names = PROTECT(Rf_allocVector(STRSXP, 3));
  SET_VECTOR_ELT(out, 0, Rf_ScalarInteger(mismatches));
  SET_VECTOR_ELT(out, 1, Rf_ScalarInteger(identifier));
  SET_VECTOR_ELT(out, 2, Rf_ScalarInteger(count));
  SET_STRING_ELT(names, 0, Rf_mkChar("mismatches"));
  SET_STRING_ELT(names, 1, Rf_mkChar("identifier"));
  SET_STRING_ELT(names, 2, Rf_mkChar("length"));
  Rf_setAttrib(out, R_NamesSymbol, names);
  UNPROTECT(2);
  return out;
  ROBMA_NATIVE_END
}

extern "C" SEXP RoBMA_selnorm_set_native_threads(SEXP threads)
{
  ROBMA_NATIVE_BEGIN(threads)
  if ((TYPEOF(threads) != INTSXP && TYPEOF(threads) != REALSXP) ||
      XLENGTH(threads) != 1 || Rf_inherits(threads, "factor")) {
    Rf_error("'threads' must be one positive integer or zero for the OpenMP default.");
  }
  const double value = Rf_asReal(threads);
  if (!std::isfinite(value) || value < 0.0 || value != std::floor(value) ||
      value > static_cast<double>(std::numeric_limits<int>::max())) {
    Rf_error("'threads' must be one nonnegative whole number of threads.");
  }
  // The replaced budget lets R-side scopes restore what they found.
  return Rf_ScalarInteger(robma_set_native_threads(static_cast<int>(value)));
  ROBMA_NATIVE_END
}

extern "C" SEXP RoBMA_selnorm_cache_snapshot()
{
  ROBMA_NATIVE_BEGIN()
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
  ROBMA_NATIVE_END
}

extern "C" SEXP RoBMA_selnorm_cache_restore(SEXP snapshots)
{
  ROBMA_NATIVE_BEGIN(snapshots)
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
  ROBMA_NATIVE_END
}

extern "C" SEXP RoBMA_selnorm_cache_control(SEXP capacity_bytes, SEXP clear_mask)
{
  ROBMA_NATIVE_BEGIN(capacity_bytes, clear_mask)
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
  ROBMA_NATIVE_END
}

extern "C" SEXP RoBMA_selnorm_sampler_control(SEXP enabled, SEXP settings, SEXP clear)
{
  ROBMA_NATIVE_BEGIN(enabled, settings, clear)
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
  ROBMA_NATIVE_END
}
