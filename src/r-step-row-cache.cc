// Test hook of the row cache of singleton step-selection rows
// (samplers/StepRowCache.h); not package interface.

#include <Rinternals.h>

#include <cmath>
#include <string>
#include <vector>

#include "r-native-boundary.h"
#include "samplers/StepRowCache.h"
#include "samplers/CorrectedBlockKernel.h"

// flags: NULL or a named list of `mutant` (0 none, 1 the power omitted, 2 the
// observed-bin term omitted, 3 a stale cache), `perturbation` and `disable`,
// read by the samplers built from now on. clear: reset the counters after reading them.
// probe: run the frozen-state probe of the live block samplers with a row
// cache. Returns the flags, the counters (refreshes, rows evaluated from the
// cache and by their node) and the probe records (chain, cache value, sum of
// the rows' graph log densities).
extern "C" SEXP RoBMA_step_row_cache_test(SEXP flags, SEXP clear, SEXP probe)
{
  ROBMA_NATIVE_BEGIN(flags, clear, probe)
  for (SEXP flag : {clear, probe}) {
    if (TYPEOF(flag) != LGLSXP || XLENGTH(flag) != 1 || LOGICAL(flag)[0] == NA_LOGICAL) {
      Rf_error("'clear' and 'probe' must be one non-missing logical value.");
    }
  }
  jags::RoBMA::StepRowCacheTestFlags current = jags::RoBMA::step_row_cache_test_flags();
  if (flags != R_NilValue) {
    if (TYPEOF(flags) != VECSXP) Rf_error("'flags' must be NULL or a named list.");
    SEXP names = Rf_getAttrib(flags, R_NamesSymbol);
    for (R_xlen_t i = 0; i < XLENGTH(flags); ++i) {
      const std::string name = names == R_NilValue ? "" : CHAR(STRING_ELT(names, i));
      SEXP value = VECTOR_ELT(flags, i);
      if (XLENGTH(value) != 1 || (TYPEOF(value) != LGLSXP && TYPEOF(value) != REALSXP &&
                                  TYPEOF(value) != INTSXP)) {
        Rf_error("Each flag must be one logical or numeric value.");
      }
      const double number = Rf_asReal(value);
      if (!std::isfinite(number)) Rf_error("Each flag must be finite.");
      if (name == "mutant") current.mutant = static_cast<int>(number);
      else if (name == "perturbation") current.perturbation = number;
      else if (name == "disable") current.disable = number != 0;
      else Rf_error("Unknown row cache test flag '%s'.", name.c_str());
    }
    if (current.mutant < 0 || current.mutant > 3) Rf_error("'mutant' must be 0, 1, 2 or 3.");
    jags::RoBMA::step_row_cache_test_set(current);
  }
  const std::vector<double> records = LOGICAL(probe)[0] ? jags::RoBMA::corrected_block_test_row_probe() :
    std::vector<double>();
  jags::RoBMA::StepRowCacheCounters &counters = jags::RoBMA::step_row_cache_counters();
  const bool reset = LOGICAL(clear)[0] != 0;
  const auto take = [reset](std::atomic<std::uint64_t> &counter) {
    return static_cast<double>(reset ? counter.exchange(0, std::memory_order_relaxed) :
                               counter.load(std::memory_order_relaxed));
  };
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 6));
  SEXP labels = PROTECT(Rf_allocVector(STRSXP, 6));
  const char *names[] = {"mutant", "perturbation", "counters", "probe", "tolerance", "disable"};
  SET_VECTOR_ELT(out, 0, Rf_ScalarReal(current.mutant));
  SET_VECTOR_ELT(out, 1, Rf_ScalarReal(current.perturbation));
  SET_VECTOR_ELT(out, 5, Rf_ScalarLogical(current.disable));
  SEXP counter_values = PROTECT(Rf_allocVector(REALSXP, 3));
  REAL(counter_values)[0] = take(counters.refreshes);
  REAL(counter_values)[1] = take(counters.cached);
  REAL(counter_values)[2] = take(counters.direct);
  SEXP counter_names = PROTECT(Rf_allocVector(STRSXP, 3));
  const char *counter_labels[] = {"refreshes", "cached", "direct"};
  for (int i = 0; i < 3; ++i) SET_STRING_ELT(counter_names, i, Rf_mkChar(counter_labels[i]));
  Rf_setAttrib(counter_values, R_NamesSymbol, counter_names);
  SET_VECTOR_ELT(out, 2, counter_values);
  SEXP probe_values = PROTECT(Rf_allocVector(REALSXP, static_cast<R_xlen_t>(records.size())));
  for (std::size_t i = 0; i < records.size(); ++i) REAL(probe_values)[i] = records[i];
  SET_VECTOR_ELT(out, 3, probe_values);
  SET_VECTOR_ELT(out, 4, Rf_ScalarReal(jags::RoBMA::step_row_cache_tolerance));
  for (int i = 0; i < 6; ++i) SET_STRING_ELT(labels, i, Rf_mkChar(names[i]));
  Rf_setAttrib(out, R_NamesSymbol, labels);
  UNPROTECT(5);
  return out;
  ROBMA_NATIVE_END
}
