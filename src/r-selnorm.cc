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
#include "selnorm/selnorm-qmc-z.h"
#include "samplers/CoarseCorrectedSlice.h"
#include "distributions/DSELNORMCLUSTERSTEP.h"
#include "distributions/DSELNORMFACTORSTEP.h"
#include "distributions/selnorm-jags-bounds.h"

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

// Test hook, not package interface: the surrogate and the exact density of one
// cluster or factor selection block given as its JAGS parameter list (in the
// parameter order of the distribution), evaluated through the distribution
// classes the sampler uses, with the validators and the surrogate's pieces
// (numerator, log normalizer of the rung, whether it fell back to the exact
// density). `settings` is the sampler's setting vector (three grid steps, rule
// count). An exception of an evaluation is returned as its message. With
// `unchecked` the evaluations also run on blocks whose controls the validators
// refuse (a quadrature table with a missing value, say), to reach the branches
// the validators otherwise guard.
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
  if (Rf_length(settings) != 4) Rf_error("'settings' must contain three grid steps and a rule count.");
  if (TYPEOF(parameters) != VECSXP || Rf_length(parameters) != 20) {
    Rf_error("'parameters' must be a list of the 20 distribution parameters.");
  }
  SelNormCoarseSettings coarse;
  coarse.mean_step = REAL(settings)[0];
  coarse.diagonal_step = REAL(settings)[1];
  coarse.log_weight_step = REAL(settings)[2];
  coarse.max_rules = static_cast<unsigned int>(REAL(settings)[3]);
  std::vector<double const *> par(20);
  std::vector<unsigned int> len(20);
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
  if (control_valid || LOGICAL(unchecked)[0] == TRUE) {
    const int rows = static_cast<int>(length);
    const int rank = factor ? static_cast<int>(len[2] / len[0]) : 1;
    const SelNormJagsBounds lower_bounds(par[5], len[5]);
    const SelNormJagsBounds upper_bounds(par[6], len[6]);
    std::vector<int> obs_bin(static_cast<std::size_t>(rows));
    for (int i = 0; i < rows; ++i) obs_bin[static_cast<std::size_t>(i)] = static_cast<int>(par[7][i]);
    try {
      surrogate = factor ?
        factor_distribution.surrogateLogDensity(x_values, length, jags::PDF_LIKELIHOOD,
          par, len, nullptr, nullptr, coarse) :
        cluster_distribution.surrogateLogDensity(x_values, length, jags::PDF_LIKELIHOOD,
          par, len, nullptr, nullptr, coarse);
    } catch (const std::exception &failure) {
      surrogate_error = failure.what();
      if (surrogate_error.empty()) surrogate_error = "evaluation failed";
    }
    try {
      exact = factor ?
        factor_distribution.logDensity(x_values, length, jags::PDF_LIKELIHOOD,
          par, len, nullptr, nullptr) :
        cluster_distribution.logDensity(x_values, length, jags::PDF_LIKELIHOOD,
          par, len, nullptr, nullptr);
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
    std::min(coarse.max_rules, static_cast<unsigned int>(len[13])))));
  for (int i = 0; i < slots; ++i) SET_STRING_ELT(names, i, Rf_mkChar(labels[i]));
  Rf_setAttrib(output, R_NamesSymbol, names);
  UNPROTECT(2);
  return output;
  ROBMA_NATIVE_END
}
