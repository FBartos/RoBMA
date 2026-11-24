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
#include "selnorm/selnorm-compact-partition.h"
#include "selnorm/zplot-context.h"
#include "selnorm/selnorm-parallel.h"
#include "selnorm/selnorm-qmc-z.h"
#include "samplers/CorrectedSampler.h"
#include "samplers/CoarseCorrectedBlock.h"
#include "samplers/ScreenedFinite.h"
#include "distributions/selnorm-node-memo.h"
#include "samplers/LocationBlockReport.h"
#include "samplers/StepRowNormalizer.h"
#include "distributions/DSELNORMCLUSTERSTEP.h"
#include "distributions/DSELNORMFACTORSTEP.h"
#include "distributions/DSELNORMSTEPSWITCH.h"
#include "distributions/selnorm-jags-bounds.h"

#ifndef FCONE
# define FCONE
#endif

extern "C" double Rf_dnorm4(double, double, double, int);

#include "r-native-api.h"

// Private implementation fragments share one anonymous namespace and are kept in
// this order so helper symbols stay internal to this translation unit.
namespace {
#include "r-native/r-selnorm-common.cc.inc"
#include "r-native/r-selnorm-funnel-common.cc.inc"
#include "r-native/r-selnorm-cluster.cc.inc"
}

#include "r-native/r-selnorm-loglik.cc.inc"
#include "r-native/r-selnorm-mv.cc.inc"
#include "r-native/r-selnorm-sampling-conditioned.cc.inc"
#include "r-native/r-selnorm-kernel.cc.inc"
#include "r-native/r-selnorm-funnel-zcurve.cc.inc"
#include "r-native/r-native-probe.cc.inc"

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

// The stopping multiplier of the randomized selection normalizers (selnorm-qmc-z.h).
// R judges returned QMC diagnostics with this value; it is the only copy R sees.
extern "C" SEXP RoBMA_selnorm_qmc_z()
{
  ROBMA_NATIVE_BEGIN()
  return Rf_ScalarReal(selnorm_qmc_z);
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
      row_count != std::floor(row_count) ||
      row_count > static_cast<double>(std::numeric_limits<int>::max())) {
    Rf_error("'rows' must be one nonnegative whole number of rows.");
  }
  if (!std::isfinite(work)) {
    Rf_error("'work_per_row' must be one finite number.");
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

// Private cache probe: all owners are local, so repeated test runs retain no
// test-only registry. Mutations deliberately reuse the same input allocations.
extern "C" SEXP RoBMA_selnorm_double_rule_weights_check(
    SEXP nodes, SEXP log_weights, SEXP shift, SEXP evict)
{
  ROBMA_NATIVE_BEGIN(nodes, log_weights, shift, evict)
  if (TYPEOF(nodes) != REALSXP || TYPEOF(log_weights) != REALSXP ||
      XLENGTH(nodes) < 1 || XLENGTH(nodes) != XLENGTH(log_weights) ||
      XLENGTH(nodes) > 100000) {
    Rf_error("'nodes' and 'log_weights' must be equal nonempty numeric vectors of at most 100000 values.");
  }
  if (TYPEOF(shift) != REALSXP || XLENGTH(shift) != 1) {
    Rf_error("'shift' must be one numeric value.");
  }
  if (TYPEOF(evict) != INTSXP || XLENGTH(evict) != 1 ||
      INTEGER(evict)[0] == NA_INTEGER || INTEGER(evict)[0] < 0 || INTEGER(evict)[0] > 64) {
    Rf_error("'evict' must be one integer between 0 and 64.");
  }
  const int count = static_cast<int>(XLENGTH(nodes));
  const double delta = REAL(shift)[0];
  std::vector<double> x(REAL(nodes), REAL(nodes) + count);
  std::vector<double> logs(REAL(log_weights), REAL(log_weights) + count);
  const auto get = [&](double d) {
    return cpp_selnorm_rule_weights_double(x.data(), logs.data(), count, d);
  };
  int mismatches = 0;
  const auto check = [&](const std::vector<double> &values, double d) {
    for (int i = 0; i < count; ++i) {
      const double fresh = std::exp(logs[i] - d * x[i] - 0.5 * d * d);
      const double served = values[i];
      if (!(std::isnan(fresh) && std::isnan(served)) &&
          (fresh != served || (fresh == 0.0 && std::signbit(fresh) != std::signbit(served)))) {
        ++mismatches;
      }
    }
  };
  const auto held = get(delta);
  const std::vector<double> saved = *held;
  check(*held, delta);
  // Distinct source allocations with the same bytes must hit a retained entry.
  const auto copied = cpp_selnorm_rule_weights_double(
    REAL(nodes), REAL(log_weights), count, delta);
  const bool reused = copied == held;
  bool distinct[3];
  x[0] = std::isfinite(x[0]) ? x[0] + 0.125 : 0.125;
  const auto changed_node = get(delta);
  distinct[0] = changed_node != held;
  check(*changed_node, delta);
  x[0] = REAL(nodes)[0];
  logs[0] = std::isfinite(logs[0]) ? logs[0] + 0.25 : -0.25;
  const auto changed_log = get(delta);
  distinct[1] = changed_log != held;
  check(*changed_log, delta);
  logs[0] = REAL(log_weights)[0];
  const double changed_delta = std::isfinite(delta) ? delta + 0.375 : 0.375;
  const auto changed_shift = get(changed_delta);
  distinct[2] = changed_shift != held;
  check(*changed_shift, changed_delta);
  bool evicted = false;
  if (INTEGER(evict)[0] > 0) {
    // Churn distinct keys; the original owner must remain readable.
    for (int i = 0; i < INTEGER(evict)[0]; ++i) {
      logs[0] = -1000.0 - i;
      const auto value = get(0.125);
      check(*value, 0.125);
    }
    logs[0] = REAL(log_weights)[0];
    evicted = get(delta) != held;
  }
  // Compare bytes: this also preserves the exact NaN payload in held storage.
  const bool retained = std::memcmp(saved.data(), held->data(),
    static_cast<std::size_t>(count) * sizeof(double)) == 0;
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 6));
  SEXP names = PROTECT(Rf_allocVector(STRSXP, 6));
  SEXP values = PROTECT(Rf_allocVector(REALSXP, count));
  std::copy(held->begin(), held->end(), REAL(values));
  SEXP mutations = PROTECT(Rf_allocVector(LGLSXP, 3));
  for (int i = 0; i < 3; ++i) LOGICAL(mutations)[i] = distinct[i];
  SET_VECTOR_ELT(out, 0, values);
  SET_VECTOR_ELT(out, 1, Rf_ScalarInteger(mismatches));
  SET_VECTOR_ELT(out, 2, Rf_ScalarLogical(reused));
  SET_VECTOR_ELT(out, 3, mutations);
  SET_VECTOR_ELT(out, 4, Rf_ScalarLogical(retained));
  SET_VECTOR_ELT(out, 5, Rf_ScalarLogical(evicted));
  const char *labels[] = {"values", "mismatches", "reused", "mutation_distinct", "retained", "evicted"};
  for (int i = 0; i < 6; ++i) SET_STRING_ELT(names, i, Rf_mkChar(labels[i]));
  Rf_setAttrib(out, R_NamesSymbol, names);
  UNPROTECT(4);
  return out;
  ROBMA_NATIVE_END
}

// Private direct-lane probe. Unsupported instructions are never called; guarded
// local buffers also verify that full blocks and padded remainders stay in bounds.
extern "C" SEXP RoBMA_selnorm_tail_lanes_check(
    SEXP locations, SEXP offset, SEXP inverse_sd, SEXP last, SEXP step)
{
  ROBMA_NATIVE_BEGIN(locations, offset, inverse_sd, last, step)
  if (TYPEOF(locations) != REALSXP || XLENGTH(locations) > 1000000) {
    Rf_error("'locations' must be a numeric vector of at most 1000000 values.");
  }
  for (SEXP argument : {offset, inverse_sd, last, step}) {
    if (TYPEOF(argument) != REALSXP || XLENGTH(argument) != 1) {
      Rf_error("'offset', 'inverse_sd', 'last' and 'step' must each be one numeric value.");
    }
  }
  const std::size_t count = static_cast<std::size_t>(XLENGTH(locations));
  const double intercept = REAL(offset)[0], scale = REAL(inverse_sd)[0];
  const double base = REAL(last)[0], change = REAL(step)[0];
  const double sentinel = 1234567.25;
  std::vector<double> input(count + 2, sentinel);
  std::copy(REAL(locations), REAL(locations) + count, input.begin() + 1);
  bool avx2 = false, avx512 = false;
#if SELNORM_TAIL_HAS_AVX2
  avx2 = selnorm_tail::avx2_available();
  avx512 = selnorm_tail::avx512_available();
#endif
  // lane 0 is automatic dispatch, 1 AVX2, 2 AVX512. Call only after gating.
  const auto evaluate = [&](int lane) -> SEXP {
    std::vector<double> tails(count + 2, sentinel), values(count + 2, sentinel);
    std::vector<double> scratch(count), magnitude(count);
    bool valid = false;
#if SELNORM_TAIL_HAS_AVX2
    if (lane == 1) {
      selnorm_tail::upper_tail_affine_avx2(input.data() + 1, count,
        intercept, scale, tails.data() + 1);
      valid = selnorm_tail::single_boundary_normalizers_avx2(input.data() + 1,
        count, intercept, scale, base, change, values.data() + 1);
    } else if (lane == 2) {
      selnorm_tail::upper_tail_affine_avx512(input.data() + 1, count,
        intercept, scale, tails.data() + 1);
      valid = selnorm_tail::single_boundary_normalizers_avx512(input.data() + 1,
        count, intercept, scale, base, change, values.data() + 1);
    } else
#endif
    {
      (void)lane;
      selnorm_tail::upper_tail_affine(input.data() + 1, count,
        intercept, scale, tails.data() + 1);
      valid = selnorm_tail::single_boundary_normalizers(input.data() + 1,
        count, intercept, scale, base, change, values.data() + 1,
        scratch.data(), magnitude.data());
    }
    const bool guards = tails.front() == sentinel && tails.back() == sentinel &&
      values.front() == sentinel && values.back() == sentinel &&
      input.front() == sentinel && input.back() == sentinel;
    SEXP result = PROTECT(Rf_allocVector(VECSXP, 4));
    SEXP names = PROTECT(Rf_allocVector(STRSXP, 4));
    SEXP tail_values = PROTECT(Rf_allocVector(REALSXP, count));
    SEXP normalizers = PROTECT(Rf_allocVector(REALSXP, count));
    std::copy(tails.begin() + 1, tails.end() - 1, REAL(tail_values));
    std::copy(values.begin() + 1, values.end() - 1, REAL(normalizers));
    SET_VECTOR_ELT(result, 0, tail_values);
    SET_VECTOR_ELT(result, 1, normalizers);
    SET_VECTOR_ELT(result, 2, Rf_ScalarLogical(valid));
    SET_VECTOR_ELT(result, 3, Rf_ScalarLogical(guards));
    const char *labels[] = {"tail", "normalizers", "valid", "guards"};
    for (int i = 0; i < 4; ++i) SET_STRING_ELT(names, i, Rf_mkChar(labels[i]));
    Rf_setAttrib(result, R_NamesSymbol, names);
    UNPROTECT(4);
    return result;
  };
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 6));
  SEXP names = PROTECT(Rf_allocVector(STRSXP, 6));
  SET_VECTOR_ELT(out, 0, Rf_ScalarLogical(avx2));
  SET_VECTOR_ELT(out, 1, Rf_ScalarLogical(avx512));
  if (avx2) SET_VECTOR_ELT(out, 2, evaluate(1));
  if (avx512) SET_VECTOR_ELT(out, 3, evaluate(2));
  SET_VECTOR_ELT(out, 4, evaluate(0));
  SEXP scalar = PROTECT(Rf_allocVector(REALSXP, count));
  selnorm_tail::upper_tail_affine(input.data() + 1, count, intercept, scale,
    REAL(scalar), true);
  SET_VECTOR_ELT(out, 5, scalar);
  const char *labels[] = {"avx2_available", "avx512_available", "avx2", "avx512", "automatic", "scalar"};
  for (int i = 0; i < 6; ++i) SET_STRING_ELT(names, i, Rf_mkChar(labels[i]));
  Rf_setAttrib(out, R_NamesSymbol, names);
  UNPROTECT(3);
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
  // Settings: three grid steps, a rule count and, optionally, the fallback
  // threshold. Settings frozen without a threshold (four values) mean 0: the
  // fallback is disabled.
  if (settings != R_NilValue) {
    if (TYPEOF(settings) != REALSXP || (XLENGTH(settings) != 4 && XLENGTH(settings) != 5)) {
      Rf_error("'settings' must contain three grid steps, a rule count and an optional acceptance threshold.");
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
    if (XLENGTH(settings) == 5 && !(values[4] >= 0 && values[4] <= 1)) {
      Rf_error("The coarse acceptance threshold must be a number between 0 and 1.");
    }
  }
  jags::RoBMA::CorrectedSamplerConfig config;
  jags::RoBMA::CorrectedSamplerStats stats;
  char error_message[512] = {};
  try {
    config = jags::RoBMA::corrected_sampler_configuration();
    if (enabled != R_NilValue) config.enabled = LOGICAL(enabled)[0] != 0;
    if (settings != R_NilValue) {
      config.settings.mean_step = REAL(settings)[0];
      config.settings.diagonal_step = REAL(settings)[1];
      config.settings.log_weight_step = REAL(settings)[2];
      config.settings.max_rules = static_cast<unsigned int>(REAL(settings)[3]);
      config.min_acceptance = XLENGTH(settings) == 5 ? REAL(settings)[4] : 0.0;
    }
    if (enabled != R_NilValue || settings != R_NilValue) {
      jags::RoBMA::configure_corrected_samplers(config);
    }
    if (LOGICAL(clear)[0]) jags::RoBMA::corrected_sampler_stats(true);
    stats = jags::RoBMA::corrected_sampler_stats();
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

// Whether the corrected sampler's factor rule of each order is affordable at the
// matching declared rank: the node budget the sampler applies when it decides
// whether a factor block has a coarse rule (selnorm_factor_rule_affordable).
extern "C" SEXP RoBMA_selnorm_factor_rule_affordable(SEXP order, SEXP rank)
{
  ROBMA_NATIVE_BEGIN(order, rank)
  if (TYPEOF(order) != REALSXP || TYPEOF(rank) != REALSXP || XLENGTH(order) != XLENGTH(rank)) {
    Rf_error("'order' and 'rank' must be numeric vectors of equal length.");
  }
  const R_xlen_t n = XLENGTH(order);
  const double *orders = REAL(order);
  const double *ranks = REAL(rank);
  for (R_xlen_t i = 0; i < n; ++i) {
    if (!std::isfinite(orders[i]) || orders[i] < 1 || orders[i] != std::floor(orders[i]) ||
        orders[i] > std::numeric_limits<int>::max() || !std::isfinite(ranks[i]) || ranks[i] < 1 ||
        ranks[i] != std::floor(ranks[i]) || ranks[i] > std::numeric_limits<int>::max()) {
      Rf_error("'order' and 'rank' must be positive whole numbers.");
    }
  }
  SEXP out = PROTECT(Rf_allocVector(LGLSXP, n));
  for (R_xlen_t i = 0; i < n; ++i) {
    LOGICAL(out)[i] = selnorm_factor_rule_affordable(static_cast<int>(orders[i]), static_cast<int>(ranks[i]));
  }
  UNPROTECT(1);
  return out;
  ROBMA_NATIVE_END
}

// An evidence record (rows of `columns.size()` values) as a numeric matrix with
// those column names.
static SEXP evidence_matrix(std::vector<double> const &values, std::vector<char const *> const &columns)
{
  const int width = static_cast<int>(columns.size());
  const int rows = static_cast<int>(values.size()) / width;
  SEXP out = PROTECT(Rf_allocMatrix(REALSXP, rows, width));
  for (int r = 0; r < rows; ++r) {
    for (int c = 0; c < width; ++c) REAL(out)[r + rows * c] = values[static_cast<std::size_t>(r * width + c)];
  }
  SEXP dimnames = PROTECT(Rf_allocVector(VECSXP, 2));
  SEXP names = PROTECT(Rf_allocVector(STRSXP, width));
  for (int c = 0; c < width; ++c) SET_STRING_ELT(names, c, Rf_mkChar(columns[static_cast<std::size_t>(c)]));
  SET_VECTOR_ELT(dimnames, 1, names);
  Rf_setAttrib(out, R_DimNamesSymbol, dimnames);
  UNPROTECT(3);
  return out;
}

// The per-fit counters of the corrected sampler in this process, read and then,
// when 'clear' is TRUE, cleared (the runtime report of each fit).
extern "C" SEXP RoBMA_selnorm_sampler_report(SEXP clear)
{
  ROBMA_NATIVE_BEGIN(clear)
  if (TYPEOF(clear) != LGLSXP || XLENGTH(clear) != 1 || LOGICAL(clear)[0] == NA_LOGICAL) {
    Rf_error("'clear' must be one non-missing logical value.");
  }
  const jags::RoBMA::CorrectedSamplerReport report =
    jags::RoBMA::corrected_sampler_report(LOGICAL(clear)[0] != 0);
  const char *names[] = {"samplers", "fallback_low_acceptance", "sampling_proposed",
    "sampling_accepted", "sampling_judged", "min_sampler_rate"};
  const double values[] = {static_cast<double>(report.samplers),
    static_cast<double>(report.fallback_low_acceptance),
    static_cast<double>(report.sampling_proposed),
    static_cast<double>(report.sampling_accepted),
    report.sampling_judged,
    report.has_min_rate ? report.min_sampler_rate : NA_REAL};
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 7));
  SEXP labels = PROTECT(Rf_allocVector(STRSXP, 7));
  for (int i = 0; i < 6; ++i) {
    SET_VECTOR_ELT(out, i, Rf_ScalarReal(values[i]));
    SET_STRING_ELT(labels, i, Rf_mkChar(names[i]));
  }
  // the screened indicators' evidence record (ScreenedFinite.h), one row per
  // indicator and chain
  SET_VECTOR_ELT(out, 6, evidence_matrix(jags::RoBMA::screened_finite_evidence(LOGICAL(clear)[0] != 0),
                                         jags::RoBMA::screened_finite_evidence_names()));
  SET_STRING_ELT(labels, 6, Rf_mkChar("screened"));
  Rf_setAttrib(out, R_NamesSymbol, labels);
  UNPROTECT(2);
  return out;
  ROBMA_NATIVE_END
}

// The evidence record of the location gates of this process (one row per gate
// slot and chain, LocationBlockReport.h), read and cleared when 'clear' is TRUE.
extern "C" SEXP RoBMA_location_block_evidence(SEXP clear)
{
  ROBMA_NATIVE_BEGIN(clear)
  if (TYPEOF(clear) != LGLSXP || XLENGTH(clear) != 1 || LOGICAL(clear)[0] == NA_LOGICAL) {
    Rf_error("'clear' must be one non-missing logical value.");
  }
  return evidence_matrix(jags::RoBMA::location_block_evidence(LOGICAL(clear)[0] != 0),
                         jags::RoBMA::location_block_evidence_names());
  ROBMA_NATIVE_END
}

// The per-fit counters of the location block samplers of this process (a
// named numeric vector), read and cleared when 'clear' is TRUE.
extern "C" SEXP RoBMA_location_block_report(SEXP clear)
{
  ROBMA_NATIVE_BEGIN(clear)
  if (TYPEOF(clear) != LGLSXP || XLENGTH(clear) != 1 || LOGICAL(clear)[0] == NA_LOGICAL) {
    Rf_error("'clear' must be one non-missing logical value.");
  }
  const std::vector<double> values = jags::RoBMA::location_block_report(LOGICAL(clear)[0] != 0);
  std::vector<char const *> const &names = jags::RoBMA::location_block_report_names();
  const R_xlen_t n = static_cast<R_xlen_t>(values.size());
  SEXP out = PROTECT(Rf_allocVector(REALSXP, n));
  SEXP labels = PROTECT(Rf_allocVector(STRSXP, n));
  for (R_xlen_t i = 0; i < n; ++i) {
    REAL(out)[i] = values[static_cast<std::size_t>(i)];
    SET_STRING_ELT(labels, i, Rf_mkChar(names[static_cast<std::size_t>(i)]));
  }
  Rf_setAttrib(out, R_NamesSymbol, labels);
  UNPROTECT(2);
  return out;
  ROBMA_NATIVE_END
}

// Test hook, not package interface: the location block sampler's mutant and
// probe (LocationBlockReport.h) for samplers built from now on; NULL keeps a
// setting. Returns the previous mutant and probe, the outcome of the last
// factory call that found a declaration, and the probe records since the last
// call (cleared).
extern "C" SEXP RoBMA_location_block_test(SEXP mutant, SEXP probe)
{
  ROBMA_NATIVE_BEGIN(mutant, probe)
  if (mutant != R_NilValue && (TYPEOF(mutant) != INTSXP || XLENGTH(mutant) != 1 ||
      INTEGER(mutant)[0] == NA_INTEGER || INTEGER(mutant)[0] < 0 || INTEGER(mutant)[0] > 5)) {
    Rf_error("The location block mutant must be NULL or one integer between 0 and 5.");
  }
  if (probe != R_NilValue && (TYPEOF(probe) != LGLSXP || XLENGTH(probe) != 1 ||
      LOGICAL(probe)[0] == NA_LOGICAL)) {
    Rf_error("The location block probe must be NULL or one non-missing logical value.");
  }
  int previous;
  if (mutant == R_NilValue) {
    previous = jags::RoBMA::location_block_test_mutant(0);
    jags::RoBMA::location_block_test_mutant(previous);
  } else {
    previous = jags::RoBMA::location_block_test_mutant(INTEGER(mutant)[0]);
  }
  const bool previous_probe = probe == R_NilValue ? jags::RoBMA::location_block_test_probe(-1) :
    jags::RoBMA::location_block_test_probe(LOGICAL(probe)[0] != 0 ? 1 : 0);
  const std::string decision = jags::RoBMA::location_block_test_decision();
  const std::vector<double> records = jags::RoBMA::location_block_test_take_probe();
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 4));
  SEXP labels = PROTECT(Rf_allocVector(STRSXP, 4));
  SET_VECTOR_ELT(out, 0, Rf_ScalarInteger(previous));
  SET_VECTOR_ELT(out, 1, Rf_ScalarLogical(previous_probe));
  SET_VECTOR_ELT(out, 2, Rf_mkString(decision.c_str()));
  SEXP values = PROTECT(Rf_allocVector(REALSXP, static_cast<R_xlen_t>(records.size())));
  for (std::size_t i = 0; i < records.size(); ++i) REAL(values)[i] = records[i];
  SET_VECTOR_ELT(out, 3, values);
  const char *names[] = {"mutant", "probe", "decision", "records"};
  for (int i = 0; i < 4; ++i) SET_STRING_ELT(labels, i, Rf_mkChar(names[i]));
  Rf_setAttrib(out, R_NamesSymbol, labels);
  UNPROTECT(3);
  return out;
  ROBMA_NATIVE_END
}

// Test hook, not package interface: the closed-form log density of one
// step-selection row (dselnorm_step_switch in step mode) and its first two
// derivatives in the mean, as the location block's analytic tilt computes
// them (StepRowNormalizer.h), at `mean`, and the distribution's own log density
// (DSELNORMSTEPSWITCH::logDensity) at every element of `means`, for central
// differences. Bounds use the JAGS encoding (+-1e300 open).
extern "C" SEXP RoBMA_step_row_derivative_probe(SEXP y, SEXP mean, SEXP row, SEXP omega, SEXP lower,
                                                SEXP upper, SEXP means)
{
  ROBMA_NATIVE_BEGIN(y, mean, row, omega, lower, upper, means)
  for (SEXP value : {y, mean, row, omega, lower, upper, means}) check_real(value, "argument");
  if (Rf_length(y) != 1 || Rf_length(mean) != 1 || Rf_length(row) != 6 || Rf_length(omega) < 1 ||
      Rf_length(lower) != Rf_length(omega) || Rf_length(upper) != Rf_length(omega)) {
    Rf_error("'row' must be c(sd, se, power, bin, sign, telescope) and the bins must match.");
  }
  const double *r = REAL(row);
  const int bins = Rf_length(omega);
  const double observation = REAL(y)[0];
  // the distribution's parameters: mean, sd, se, power, omega, lower, upper,
  // bin, sign, kernel mode (step), telescoping
  double mean_value = REAL(mean)[0], mode = 1.0;
  std::vector<double const *> par = {&mean_value, &r[0], &r[1], &r[2], REAL(omega), REAL(lower), REAL(upper),
                                     &r[3], &r[4], &mode, &r[5]};
  std::vector<unsigned long> len = {1, 1, 1, 1, static_cast<unsigned long>(bins), static_cast<unsigned long>(bins),
                                    static_cast<unsigned long>(bins), 1, 1, 1, 1};
  jags::RoBMA::DSELNORMSTEPSWITCH distribution;
  if (!distribution.checkParameterLength(len) || !distribution.checkParameterValue(par, len)) {
    Rf_error("The row's parameters are not valid for dselnorm_step_switch.");
  }
  const SelNormJagsBounds decoded_lower(REAL(lower), static_cast<unsigned int>(bins));
  const SelNormJagsBounds decoded_upper(REAL(upper), static_cast<unsigned int>(bins));
  jags::RoBMA::StepRowLaw law;
  law.sd = r[0];
  law.se = r[1];
  law.effect_sign = static_cast<int>(r[4]);
  law.bins = bins;
  law.omega = REAL(omega);
  law.lower = decoded_lower.data();
  law.upper = decoded_upper.data();
  double log_a = NA_REAL, d1 = NA_REAL, d2 = NA_REAL;
  const bool ok = jags::RoBMA::step_row_log_normalizer(law, mean_value, &log_a, &d1, &d2);
  // the row's log density power * [log N(y; m, sd) + log omega_bin - log A(m)]
  const double sd = r[0], power = r[2];
  const int bin = static_cast<int>(r[3]);
  const double residual = observation - mean_value;
  const double value = power * (-0.5 * std::log(2 * M_PI) - std::log(sd) - 0.5 * residual * residual / (sd * sd) +
                                std::log(REAL(omega)[bin - 1]) - log_a);
  const double first = power * (residual / (sd * sd) - d1);
  const double second = power * (-1.0 / (sd * sd) - d2);
  const R_xlen_t n = XLENGTH(means);
  SEXP exact = PROTECT(Rf_allocVector(REALSXP, n));
  SEXP direct = PROTECT(Rf_allocVector(REALSXP, n));
  for (R_xlen_t i = 0; i < n; ++i) {
    mean_value = REAL(means)[i];
    REAL(exact)[i] = distribution.logDensity(&observation, jags::PDF_LIKELIHOOD, par, len);
    // the location block's direct evaluation of the same density (B1)
    REAL(direct)[i] = jags::RoBMA::step_row_log_density(law, power, bin, 1, r[5] == 1, observation, mean_value);
  }
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 6));
  SEXP labels = PROTECT(Rf_allocVector(STRSXP, 6));
  SET_VECTOR_ELT(out, 0, Rf_ScalarLogical(ok));
  SET_VECTOR_ELT(out, 1, Rf_ScalarReal(value));
  SET_VECTOR_ELT(out, 2, Rf_ScalarReal(first));
  SET_VECTOR_ELT(out, 3, Rf_ScalarReal(second));
  SET_VECTOR_ELT(out, 4, exact);
  SET_VECTOR_ELT(out, 5, direct);
  const char *names[] = {"ok", "value", "first", "second", "exact", "direct"};
  for (int i = 0; i < 6; ++i) SET_STRING_ELT(labels, i, Rf_mkChar(names[i]));
  Rf_setAttrib(out, R_NamesSymbol, labels);
  UNPROTECT(4);
  return out;
  ROBMA_NATIVE_END
}

// Test hook, not package interface: the location block sampler's switches
// (LocationBlockReport.h, location_block_test_config) for samplers built from
// now on, as an integer vector c(refresh, fit_rule, timing) whose negative
// elements keep a setting, or NULL. Returns the previous settings.
extern "C" SEXP RoBMA_location_block_test_config(SEXP config)
{
  ROBMA_NATIVE_BEGIN(config)
  std::vector<int> wanted = {-1, -1, -1};
  if (config != R_NilValue) {
    if (TYPEOF(config) != INTSXP || XLENGTH(config) != 3) {
      Rf_error("The location block configuration must be NULL or three integers.");
    }
    const int upper[3] = {4, 2, 1};
    for (int i = 0; i < 3; ++i) {
      const int value = INTEGER(config)[i];
      if (value == NA_INTEGER || value > upper[i]) Rf_error("The location block configuration is out of range.");
      wanted[static_cast<std::size_t>(i)] = value;
    }
  }
  const std::vector<int> previous = jags::RoBMA::location_block_test_config(wanted[0], wanted[1], wanted[2]);
  SEXP out = PROTECT(Rf_allocVector(INTSXP, 3));
  SEXP labels = PROTECT(Rf_allocVector(STRSXP, 3));
  const char *names[] = {"refresh", "fit_rule", "timing"};
  for (int i = 0; i < 3; ++i) {
    INTEGER(out)[i] = previous[static_cast<std::size_t>(i)];
    SET_STRING_ELT(labels, i, Rf_mkChar(names[i]));
  }
  Rf_setAttrib(out, R_NamesSymbol, labels);
  UNPROTECT(2);
  return out;
  ROBMA_NATIVE_END
}

// Test hook, not package interface: the corrected block sampler's fixed-order
// mutant (for block samplers built from now on) and its evaluation trace. A
// NULL argument keeps the setting. Returns the settings, the slice intervals
// exhausted in corrected mode and the trace records since the last call
// (both cleared): per evaluation kind (0 surrogate, 1 full), chain, d, the d
// member values and the log density without the Jacobian; and the current
// widths of the live block samplers (per sampler chain, d, the d widths).
extern "C" SEXP RoBMA_corrected_block_test(SEXP fixed_order, SEXP trace)
{
  ROBMA_NATIVE_BEGIN(fixed_order, trace)
  for (SEXP flag : {fixed_order, trace}) {
    if (flag != R_NilValue && (TYPEOF(flag) != LGLSXP || XLENGTH(flag) != 1 ||
        LOGICAL(flag)[0] == NA_LOGICAL)) {
      Rf_error("The block sampler test flags must be NULL or one non-missing logical value.");
    }
  }
  const bool new_fixed_order = fixed_order == R_NilValue ?
    jags::RoBMA::corrected_block_test_fixed_order() : LOGICAL(fixed_order)[0] != 0;
  const bool new_trace = trace == R_NilValue ?
    jags::RoBMA::corrected_block_test_trace() : LOGICAL(trace)[0] != 0;
  jags::RoBMA::corrected_block_test_set(new_fixed_order, new_trace);
  const double exhausted = static_cast<double>(jags::RoBMA::corrected_block_test_take_exhausted());
  const std::vector<double> records = jags::RoBMA::corrected_block_test_take_trace();
  const std::vector<double> widths = jags::RoBMA::corrected_block_test_widths();
  const char *names[] = {"fixed_order", "trace", "exhausted", "records", "widths"};
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 5));
  SEXP labels = PROTECT(Rf_allocVector(STRSXP, 5));
  SET_VECTOR_ELT(out, 0, Rf_ScalarLogical(new_fixed_order));
  SET_VECTOR_ELT(out, 1, Rf_ScalarLogical(new_trace));
  SET_VECTOR_ELT(out, 2, Rf_ScalarReal(exhausted));
  SEXP values = PROTECT(Rf_allocVector(REALSXP, static_cast<R_xlen_t>(records.size())));
  for (std::size_t i = 0; i < records.size(); ++i) REAL(values)[i] = records[i];
  SET_VECTOR_ELT(out, 3, values);
  SEXP width_values = PROTECT(Rf_allocVector(REALSXP, static_cast<R_xlen_t>(widths.size())));
  for (std::size_t i = 0; i < widths.size(); ++i) REAL(width_values)[i] = widths[i];
  SET_VECTOR_ELT(out, 4, width_values);
  for (int i = 0; i < 5; ++i) SET_STRING_ELT(labels, i, Rf_mkChar(names[i]));
  Rf_setAttrib(out, R_NamesSymbol, labels);
  UNPROTECT(4);
  return out;
  ROBMA_NATIVE_END
}

// Test hook, not package interface: the flags of the corrected samplers'
// factories and of the node memo, and their process counters. `flags` is NULL
// (keep) or a list with any of disable_new (no gated weight blocks, screened
// indicators or gamma groups in models built from now on; disable_gated,
// disable_finite and disable_sum leave out one of them), drop_sum_jacobian
// (the gamma groups' mutant without sum log p), finite_proposal (0 corrected,
// 1 uncorrected mutant, 2 Liu's Metropolized Gibbs), force_fallback (every sampler
// switches to its exact kernel at the end of adaptation), freeze_finite (the
// screened indicators keep their values) and memo_slots (0, 1 or 2
// entries per selection node from now on). Returns the flags, the counters per
// sampler kind (rows slice, block, gated, finite), the memo counters
// (lookups and hits of the cluster and factor nodes) and, when `scores` is TRUE,
// the scores of the
// live screened samplers at their chains' current states (per sampler its
// chain, the number of levels L, L surrogate and L exact scores); `clear`
// resets the counters after reading them.
extern "C" SEXP RoBMA_corrected_sampler_test(SEXP flags, SEXP clear, SEXP scores)
{
  ROBMA_NATIVE_BEGIN(flags, clear, scores)
  if (TYPEOF(clear) != LGLSXP || XLENGTH(clear) != 1 || LOGICAL(clear)[0] == NA_LOGICAL) {
    Rf_error("'clear' must be one non-missing logical value.");
  }
  jags::RoBMA::CorrectedTestFlags current = jags::RoBMA::corrected_test_flags();
  int memo_slots = selnorm_node_memo_slots();
  if (flags != R_NilValue) {
    if (TYPEOF(flags) != VECSXP) Rf_error("'flags' must be NULL or a named list.");
    SEXP flag_names = Rf_getAttrib(flags, R_NamesSymbol);
    for (R_xlen_t i = 0; i < XLENGTH(flags); ++i) {
      const std::string name = flag_names == R_NilValue ? "" : CHAR(STRING_ELT(flag_names, i));
      SEXP value = VECTOR_ELT(flags, i);
      if (XLENGTH(value) != 1 || (TYPEOF(value) != LGLSXP && TYPEOF(value) != REALSXP)) {
        Rf_error("Each flag must be one logical or numeric value.");
      }
      const double number = TYPEOF(value) == LGLSXP ? LOGICAL(value)[0] : REAL(value)[0];
      if (!std::isfinite(number)) Rf_error("Each flag must be finite.");
      if (name == "disable_new") current.disable_new = number != 0;
      else if (name == "disable_gated") current.disable_gated = number != 0;
      else if (name == "disable_finite") current.disable_finite = number != 0;
      else if (name == "disable_sum") current.disable_sum = number != 0;
      else if (name == "drop_sum_jacobian") current.drop_sum_jacobian = number != 0;
      else if (name == "finite_proposal") current.finite_proposal = static_cast<int>(number);
      else if (name == "memo_slots") memo_slots = static_cast<int>(number);
      else if (name == "force_fallback") current.force_fallback = number != 0;
      else if (name == "freeze_finite") current.freeze_finite = number != 0;
      else Rf_error("Unknown corrected sampler test flag '%s'.", name.c_str());
    }
    if (memo_slots < 0 || memo_slots > 2) Rf_error("'memo_slots' must be 0, 1 or 2.");
    jags::RoBMA::corrected_test_set_flags(current);
    selnorm_node_memo_set_slots(memo_slots);
  }
  const bool reset = LOGICAL(clear)[0] != 0;
  const auto take = [reset](std::atomic<std::uint64_t> &counter) {
    return static_cast<double>(reset ? counter.exchange(0, std::memory_order_relaxed) :
                               counter.load(std::memory_order_relaxed));
  };
  const int kinds = jags::RoBMA::corrected_kind_count;
  SEXP counters = PROTECT(Rf_allocMatrix(REALSXP, kinds, 8));
  for (int k = 0; k < kinds; ++k) {
    jags::RoBMA::CorrectedKindCounters &kind = jags::RoBMA::corrected_kind_counters(k);
    double *column = REAL(counters) + k;
    column[0 * kinds] = take(kind.samplers);
    column[1 * kinds] = take(kind.updates);
    column[2 * kinds] = take(kind.proposed);
    column[3 * kinds] = take(kind.accepted);
    column[4 * kinds] = take(kind.evidence_proposed);
    column[5 * kinds] = take(kind.evidence_accepted);
    column[6 * kinds] = take(kind.fallbacks);
    column[7 * kinds] = take(kind.exact_scores);
  }
  SEXP dimnames = PROTECT(Rf_allocVector(VECSXP, 2));
  SEXP rows = PROTECT(Rf_allocVector(STRSXP, kinds));
  const char *row_names[] = {"slice", "block", "gated", "finite", "location"};
  for (int k = 0; k < kinds; ++k) SET_STRING_ELT(rows, k, Rf_mkChar(row_names[k]));
  SEXP columns = PROTECT(Rf_allocVector(STRSXP, 8));
  const char *column_names[] = {"samplers", "updates", "proposed", "accepted", "evidence_proposed",
                                "evidence_accepted", "fallbacks", "exact_scores"};
  for (int i = 0; i < 8; ++i) SET_STRING_ELT(columns, i, Rf_mkChar(column_names[i]));
  SET_VECTOR_ELT(dimnames, 0, rows);
  SET_VECTOR_ELT(dimnames, 1, columns);
  Rf_setAttrib(counters, R_DimNamesSymbol, dimnames);

  SelNormMemoCounters &memo = selnorm_node_memo_counters();
  SEXP memo_values = PROTECT(Rf_allocVector(REALSXP, 5));
  REAL(memo_values)[0] = take(memo.lookups[0]);
  REAL(memo_values)[1] = take(memo.hits[0]);
  REAL(memo_values)[2] = take(memo.lookups[1]);
  REAL(memo_values)[3] = take(memo.hits[1]);
  REAL(memo_values)[4] = selnorm_node_memo_slots();
  SEXP memo_names = PROTECT(Rf_allocVector(STRSXP, 5));
  const char *memo_labels[] = {"cluster_lookups", "cluster_hits", "factor_lookups", "factor_hits", "slots"};
  for (int i = 0; i < 5; ++i) SET_STRING_ELT(memo_names, i, Rf_mkChar(memo_labels[i]));
  Rf_setAttrib(memo_values, R_NamesSymbol, memo_names);

  // evaluating the live screened samplers sets their nodes and evaluates their
  // children (memo and counters included), so it runs only on request
  if (TYPEOF(scores) != LGLSXP || XLENGTH(scores) != 1 || LOGICAL(scores)[0] == NA_LOGICAL) {
    Rf_error("'scores' must be one non-missing logical value.");
  }
  const std::vector<double> screened = LOGICAL(scores)[0] ? jags::RoBMA::screened_finite_test_scores() :
    std::vector<double>();
  SEXP score_values = PROTECT(Rf_allocVector(REALSXP, static_cast<R_xlen_t>(screened.size())));
  for (std::size_t i = 0; i < screened.size(); ++i) REAL(score_values)[i] = screened[i];

  SEXP flag_values = PROTECT(Rf_allocVector(VECSXP, 9));
  SET_VECTOR_ELT(flag_values, 0, Rf_ScalarLogical(current.disable_new));
  SET_VECTOR_ELT(flag_values, 1, Rf_ScalarLogical(current.drop_sum_jacobian));
  SET_VECTOR_ELT(flag_values, 2, Rf_ScalarReal(current.finite_proposal));
  SET_VECTOR_ELT(flag_values, 3, Rf_ScalarReal(selnorm_node_memo_slots()));
  SET_VECTOR_ELT(flag_values, 4, Rf_ScalarLogical(current.force_fallback));
  SET_VECTOR_ELT(flag_values, 5, Rf_ScalarLogical(current.freeze_finite));
  SET_VECTOR_ELT(flag_values, 6, Rf_ScalarLogical(current.disable_gated));
  SET_VECTOR_ELT(flag_values, 7, Rf_ScalarLogical(current.disable_finite));
  SET_VECTOR_ELT(flag_values, 8, Rf_ScalarLogical(current.disable_sum));
  SEXP flag_labels = PROTECT(Rf_allocVector(STRSXP, 9));
  const char *flag_label_names[] = {"disable_new", "drop_sum_jacobian", "finite_proposal", "memo_slots",
                                    "force_fallback", "freeze_finite", "disable_gated", "disable_finite",
                                    "disable_sum"};
  for (int i = 0; i < 9; ++i) SET_STRING_ELT(flag_labels, i, Rf_mkChar(flag_label_names[i]));
  Rf_setAttrib(flag_values, R_NamesSymbol, flag_labels);

  SEXP out = PROTECT(Rf_allocVector(VECSXP, 4));
  SET_VECTOR_ELT(out, 0, flag_values);
  SET_VECTOR_ELT(out, 1, counters);
  SET_VECTOR_ELT(out, 2, memo_values);
  SET_VECTOR_ELT(out, 3, score_values);
  SEXP labels = PROTECT(Rf_allocVector(STRSXP, 4));
  const char *names[] = {"flags", "kinds", "memo", "screened_scores"};
  for (int i = 0; i < 4; ++i) SET_STRING_ELT(labels, i, Rf_mkChar(names[i]));
  Rf_setAttrib(out, R_NamesSymbol, labels);
  UNPROTECT(11);
  return out;
  ROBMA_NATIVE_END
}

// Test hook, not package interface: the surrogate and the exact density of one
// cluster or factor selection block given as its JAGS parameter list (in the
// parameter order of the distribution), evaluated through the distribution
// classes the sampler uses, with the validators and the surrogate's pieces
// (numerator, log normalizer of the rung, whether it fell back to the exact
// density). `settings` is the sampler's setting vector (three grid steps, rule
// count) and optionally a named rule (SelNormCoarseSettings::rule). An
// exception of an evaluation is returned as its message. With `unchecked` the
// evaluations also run on blocks whose validators refuse them only because the
// quadrature nodes or log weights have non-finite values, to reach the branches
// those values guard; every other control (lengths, orders that sum to the
// tables, bins, the QMC array) must pass, since the evaluators index with them.
extern "C" SEXP RoBMA_selnorm_surrogate_probe(SEXP kind, SEXP x, SEXP parameters, SEXP settings,
                                              SEXP unchecked)
{
  ROBMA_NATIVE_BEGIN(kind, x, parameters, settings, unchecked)
  if (TYPEOF(kind) != STRSXP || Rf_length(kind) != 1) {
    Rf_error("'kind' must be 'cluster' or 'factor'.");
  }
  const std::string name = CHAR(STRING_ELT(kind, 0));
  const bool factor = name == "factor";
  if (!factor && name != "cluster") Rf_error("'kind' must be 'cluster' or 'factor'.");
  check_real(x, "x");
  check_real(settings, "settings");
  check_logical(unchecked, "unchecked");
  if (Rf_length(settings) != 4 && Rf_length(settings) != 5) {
    Rf_error("'settings' must contain three grid steps, a rule count and optionally a named rule.");
  }
  if (TYPEOF(parameters) != VECSXP || Rf_length(parameters) != 20) {
    Rf_error("'parameters' must be a list of the 20 distribution parameters.");
  }
  SelNormCoarseSettings coarse;
  coarse.mean_step = REAL(settings)[0];
  coarse.diagonal_step = REAL(settings)[1];
  coarse.log_weight_step = REAL(settings)[2];
  coarse.max_rules = static_cast<unsigned int>(REAL(settings)[3]);
  if (Rf_length(settings) == 5) coarse.rule = static_cast<unsigned int>(REAL(settings)[4]);
  std::vector<double const *> par(20);
  std::vector<unsigned long> len(20);
  for (int i = 0; i < 20; ++i) {
    SEXP element = VECTOR_ELT(parameters, i);
    check_real(element, "parameter");
    par[static_cast<std::size_t>(i)] = REAL(element);
    len[static_cast<std::size_t>(i)] = static_cast<unsigned int>(Rf_length(element));
  }
  jags::RoBMA::DSELNORMCLUSTERSTEP cluster_distribution;
  jags::RoBMA::DSELNORMFACTORSTEP factor_distribution;
  const bool lengths_valid = factor ? factor_distribution.checkParameterLength(len) :
    cluster_distribution.checkParameterLength(len);
  if (!lengths_valid) Rf_error("The parameter lengths do not match the distribution.");
  const double *x_values = REAL(x);
  const unsigned int length = static_cast<unsigned int>(Rf_length(x));
  if (length != len[0]) Rf_error("'x' must have the length of the mean.");

  const bool dynamic_valid = factor ? factor_distribution.checkDynamicParameterValue(par, len) :
    cluster_distribution.checkDynamicParameterValue(par, len);
  const bool control_valid = factor ? factor_distribution.checkControlParameterValue(par, len) :
    cluster_distribution.checkControlParameterValue(par, len);
  const bool valid = factor ? factor_distribution.checkParameterValue(par, len) :
    cluster_distribution.checkParameterValue(par, len);
  const double nan = std::numeric_limits<double>::quiet_NaN();
  double surrogate = nan, exact = nan, numerator = nan, log_normalizer = nan;
  bool exact_fallback = false;
  std::string surrogate_error, exact_error;
  // Invalid controls are refused; invalid dynamic values are evaluated, and the
  // evaluators give them the exact kernel's -Inf.
  bool evaluate = control_valid;
  if (!control_valid && LOGICAL(unchecked)[0] == TRUE) {
    // the controls with the two tables replaced by finite ones
    const std::vector<double> finite_nodes(len[11], 0.0), finite_log_weights(len[12], 0.0);
    std::vector<double const *> finite_par(par);
    finite_par[11] = finite_nodes.data();
    finite_par[12] = finite_log_weights.data();
    const bool structure_valid = factor ? factor_distribution.checkControlParameterValue(finite_par, len) :
      cluster_distribution.checkControlParameterValue(finite_par, len);
    if (!structure_valid) {
      Rf_error("'unchecked' admits only non-finite values in the quadrature nodes and log weights.");
    }
    evaluate = true;
  }
  if (evaluate) {
    const int rows = static_cast<int>(length);
    const int rank = factor ? static_cast<int>(len[2] / len[0]) : 1;
    const SelNormJagsBounds lower_bounds(par[5], len[5]);
    const SelNormJagsBounds upper_bounds(par[6], len[6]);
    std::vector<int> obs_bin(static_cast<std::size_t>(rows));
    for (int i = 0; i < rows; ++i) obs_bin[static_cast<std::size_t>(i)] = static_cast<int>(par[7][i]);
    try {
      surrogate = factor ?
        factor_distribution.surrogateLogDensity(x_values, jags::PDF_LIKELIHOOD,
          par, len, coarse) :
        cluster_distribution.surrogateLogDensity(x_values, jags::PDF_LIKELIHOOD,
          par, len, coarse);
    } catch (const std::exception &failure) {
      surrogate_error = failure.what();
      if (surrogate_error.empty()) surrogate_error = "evaluation failed";
    }
    try {
      exact = factor ?
        factor_distribution.logDensity(x_values, jags::PDF_LIKELIHOOD,
          par, len) :
        cluster_distribution.logDensity(x_values, jags::PDF_LIKELIHOOD,
          par, len);
    } catch (const std::exception &failure) {
      exact_error = failure.what();
      if (exact_error.empty()) exact_error = "evaluation failed";
    }
    // The pieces of the surrogate: the exact kernel's numerator (its normalizer
    // overridden by zero) and the log normalizer of the rung from the native
    // evaluator the class calls.
    try {
      const int sign = static_cast<int>(*par[8]);
      const bool telescope = static_cast<int>(*par[9]) == 1;
      const int mode = static_cast<int>(*par[10]);
      const int rule = static_cast<int>(*par[19]);
      const int n_bins = static_cast<int>(len[4]);
      const int initial_points = static_cast<int>(*par[15]);
      const int maximum_points = static_cast<int>(*par[16]);
      const int scrambles = static_cast<int>(*par[17]);
      const double zero = 0.0;
      double mcse = 0.0, change = 0.0;
      if (factor) {
        numerator = cpp_selnorm_factor_step_lpdf(x_values, par[0], par[1], par[2], rows, rank,
          par[3], par[4], n_bins, lower_bounds.data(), upper_bounds.data(), obs_bin.data(),
          sign, telescope, mode, par[11], par[12], par[13], static_cast<int>(len[13]), par[14],
          initial_points, maximum_points, scrambles, *par[18], &mcse, &change, nullptr, rule,
          nullptr, nullptr, &zero);
        cpp_selnorm_factor_step_surrogate_lpdf(x_values, par[0], par[1], par[2], rows, rank,
          par[3], par[4], n_bins, lower_bounds.data(), upper_bounds.data(), obs_bin.data(),
          sign, telescope, mode, par[11], par[12], par[13], static_cast<int>(len[13]), rule,
          coarse, &exact_fallback, &log_normalizer);
      } else {
        numerator = cpp_selnorm_cluster_step_lpdf(x_values, par[0], par[1], par[2], rows,
          par[3], par[4], n_bins, lower_bounds.data(), upper_bounds.data(), obs_bin.data(),
          sign, telescope, mode, par[11], par[12], par[13], static_cast<int>(len[13]), par[14],
          initial_points, maximum_points, scrambles, *par[18], &mcse, &change, nullptr, nullptr,
          rule, nullptr, &zero);
        cpp_selnorm_cluster_step_surrogate_lpdf(x_values, par[0], par[1], par[2], rows,
          par[3], par[4], n_bins, lower_bounds.data(), upper_bounds.data(), obs_bin.data(),
          sign, telescope, mode, par[11], par[12], par[13], static_cast<int>(len[13]), rule,
          coarse, &exact_fallback, &log_normalizer);
      }
    } catch (const std::exception &) {
      numerator = nan;
    }
  }

  const int slots = 11;
  const char *labels[slots] = {"surrogate", "exact", "numerator", "log_normalizer", "exact_fallback",
    "valid", "dynamic_valid", "control_valid", "surrogate_error", "exact_error", "rule"};
  SEXP output = PROTECT(Rf_allocVector(VECSXP, slots));
  SEXP names = PROTECT(Rf_allocVector(STRSXP, slots));
  SET_VECTOR_ELT(output, 0, Rf_ScalarReal(surrogate));
  SET_VECTOR_ELT(output, 1, Rf_ScalarReal(exact));
  SET_VECTOR_ELT(output, 2, Rf_ScalarReal(numerator));
  SET_VECTOR_ELT(output, 3, Rf_ScalarReal(log_normalizer));
  SET_VECTOR_ELT(output, 4, Rf_ScalarLogical(exact_fallback));
  SET_VECTOR_ELT(output, 5, Rf_ScalarLogical(valid));
  SET_VECTOR_ELT(output, 6, Rf_ScalarLogical(dynamic_valid));
  SET_VECTOR_ELT(output, 7, Rf_ScalarLogical(control_valid));
  SET_VECTOR_ELT(output, 8, Rf_mkString(surrogate_error.c_str()));
  SET_VECTOR_ELT(output, 9, Rf_mkString(exact_error.c_str()));
  SET_VECTOR_ELT(output, 10, Rf_ScalarInteger(static_cast<int>(
    std::min(coarse.rule > 0 ? coarse.rule : coarse.max_rules, static_cast<unsigned int>(len[13])))));
  for (int i = 0; i < slots; ++i) SET_STRING_ELT(names, i, Rf_mkChar(labels[i]));
  Rf_setAttrib(output, R_NamesSymbol, names);
  UNPROTECT(2);
  return output;
  ROBMA_NATIVE_END
}
