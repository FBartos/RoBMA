// R entry points of the location block's heterogeneity scale ownership
// (samplers/LocationScale.h): the per-fit report and the test hooks.

#include <Rinternals.h>
#include <R_ext/Random.h>

#include <string>
#include <vector>

#include "samplers/LocationScaleReport.h"

#include "r-native-api.h"

namespace {

SEXP real_vector(std::vector<double> const &values)
{
  SEXP out = PROTECT(Rf_allocVector(REALSXP, static_cast<R_xlen_t>(values.size())));
  for (std::size_t i = 0; i < values.size(); ++i) REAL(out)[i] = values[i];
  UNPROTECT(1);
  return out;
}

SEXP named_list(std::vector<SEXP> const &values, std::vector<char const *> const &names)
{
  const R_xlen_t n = static_cast<R_xlen_t>(values.size());
  SEXP out = PROTECT(Rf_allocVector(VECSXP, n));
  SEXP labels = PROTECT(Rf_allocVector(STRSXP, n));
  for (R_xlen_t i = 0; i < n; ++i) {
    SET_VECTOR_ELT(out, i, values[static_cast<std::size_t>(i)]);
    SET_STRING_ELT(labels, i, Rf_mkChar(names[static_cast<std::size_t>(i)]));
  }
  Rf_setAttrib(out, R_NamesSymbol, labels);
  UNPROTECT(2);
  return out;
}

bool flag(SEXP value)
{
  if (TYPEOF(value) != LGLSXP || XLENGTH(value) != 1 || LOGICAL(value)[0] == NA_LOGICAL) {
    Rf_error("A flag must be one non-missing logical value.");
  }
  return LOGICAL(value)[0] != 0;
}

int integer_scalar(SEXP value, char const *what)
{
  if (TYPEOF(value) != INTSXP || XLENGTH(value) != 1 || INTEGER(value)[0] == NA_INTEGER) {
    Rf_error("'%s' must be one non-missing integer.", what);
  }
  return INTEGER(value)[0];
}

std::vector<double> real_values(SEXP value, char const *what)
{
  if (TYPEOF(value) != REALSXP) Rf_error("'%s' must be a numeric vector.", what);
  return std::vector<double>(REAL(value), REAL(value) + XLENGTH(value));
}

}

// The per-fit counters (named numeric) and the factory decisions (columns
// layout, L, K, status, reason, check, embedded) of this process, read and
// cleared when 'clear' is TRUE.
extern "C" SEXP RoBMA_location_scale_report(SEXP clear)
{
  ROBMA_NATIVE_BEGIN(clear)
  const bool clearing = flag(clear);
  const std::vector<double> values = jags::RoBMA::location_scale_report(clearing);
  std::vector<char const *> const &names = jags::RoBMA::location_scale_report_names();
  SEXP counters = PROTECT(real_vector(values));
  SEXP labels = PROTECT(Rf_allocVector(STRSXP, static_cast<R_xlen_t>(names.size())));
  for (std::size_t i = 0; i < names.size(); ++i) SET_STRING_ELT(labels, static_cast<R_xlen_t>(i), Rf_mkChar(names[i]));
  Rf_setAttrib(counters, R_NamesSymbol, labels);
  const std::vector<jags::RoBMA::LocationScaleDecision> decisions = jags::RoBMA::location_scale_decisions(clearing);
  const R_xlen_t n = static_cast<R_xlen_t>(decisions.size());
  SEXP layout = PROTECT(Rf_allocVector(STRSXP, n));
  SEXP status = PROTECT(Rf_allocVector(STRSXP, n));
  SEXP reason = PROTECT(Rf_allocVector(STRSXP, n));
  SEXP check = PROTECT(Rf_allocVector(STRSXP, n));
  SEXP embedded = PROTECT(Rf_allocVector(STRSXP, n));
  SEXP L = PROTECT(Rf_allocVector(INTSXP, n));
  SEXP K = PROTECT(Rf_allocVector(INTSXP, n));
  for (R_xlen_t i = 0; i < n; ++i) {
    jags::RoBMA::LocationScaleDecision const &d = decisions[static_cast<std::size_t>(i)];
    SET_STRING_ELT(layout, i, Rf_mkChar(d.layout.c_str()));
    SET_STRING_ELT(status, i, Rf_mkChar(d.status.c_str()));
    SET_STRING_ELT(reason, i, Rf_mkChar(d.reason.c_str()));
    SET_STRING_ELT(check, i, Rf_mkChar(d.check.c_str()));
    SET_STRING_ELT(embedded, i, Rf_mkChar(d.embedded.c_str()));
    INTEGER(L)[i] = d.L;
    INTEGER(K)[i] = d.K;
  }
  SEXP table = PROTECT(named_list({layout, L, K, status, reason, check, embedded},
                                  {"layout", "L", "K", "status", "reason", "check", "embedded"}));
  SEXP out = named_list({counters, table}, {"counters", "decisions"});
  UNPROTECT(10);
  return out;
  ROBMA_NATIVE_END
}

// Test hook, not package interface: the mutant, decline, timing and check
// switches (LocationScaleReport.h) for scale owners built from now on (three
// integers keep the check); NULL or negative values keep a setting. Returns the
// previous four.
extern "C" SEXP RoBMA_location_scale_test_config(SEXP config)
{
  ROBMA_NATIVE_BEGIN(config)
  std::vector<int> wanted = {-1, -1, -1, -1};
  if (config != R_NilValue) {
    if (TYPEOF(config) != INTSXP || (XLENGTH(config) != 3 && XLENGTH(config) != 4)) {
      Rf_error("The location scale configuration must be NULL or three or four integers.");
    }
    const int upper[4] = {4, 1, 1, 1000000};
    for (R_xlen_t i = 0; i < XLENGTH(config); ++i) {
      const int value = INTEGER(config)[i];
      if (value == NA_INTEGER || value > upper[i]) Rf_error("The location scale configuration is out of range.");
      wanted[static_cast<std::size_t>(i)] = value;
    }
  }
  const std::vector<int> previous =
    jags::RoBMA::location_scale_test_config(wanted[0], wanted[1], wanted[2], wanted[3]);
  SEXP out = PROTECT(Rf_allocVector(INTSXP, 4));
  SEXP labels = PROTECT(Rf_allocVector(STRSXP, 4));
  const char *names[] = {"mutant", "decline", "timing", "check"};
  for (int i = 0; i < 4; ++i) {
    INTEGER(out)[i] = previous[static_cast<std::size_t>(i)];
    SET_STRING_ELT(labels, i, Rf_mkChar(names[i]));
  }
  Rf_setAttrib(out, R_NamesSymbol, labels);
  UNPROTECT(2);
  return out;
  ROBMA_NATIVE_END
}

// Test hook, not package interface (LocationScaleReport.h). 'what':
//   "last"       the serial number of the most recently built live owner;
//   "state"      the state of chain 'chain' (0-based) of owner 'index';
//   "set"        writes the coordinates 'values';
//   "transform"  applies T_c with u = values[1]; FALSE outside the support;
//   "density"    log f(u) at every u in 'values';
//   "slice"      one slice update from each start values[-(1:2)] of a uniform
//                target on (values[1], values[2]); 'index' 1 the fixed-cap
//                mutant (R's random number generator).
extern "C" SEXP RoBMA_location_scale_test(SEXP what, SEXP index, SEXP chain, SEXP values)
{
  ROBMA_NATIVE_BEGIN(what, index, chain, values)
  if (TYPEOF(what) != STRSXP || XLENGTH(what) != 1) Rf_error("'what' must be one string.");
  const std::string action = CHAR(STRING_ELT(what, 0));
  if (action == "last") return Rf_ScalarInteger(jags::RoBMA::location_scale_test_last());
  const int owner = integer_scalar(index, "index");
  const int chain_index = integer_scalar(chain, "chain");
  if (chain_index < 0) Rf_error("'chain' must be nonnegative.");
  const unsigned int ch = static_cast<unsigned int>(chain_index);
  if (action == "state") {
    const jags::RoBMA::LocationScaleState state = jags::RoBMA::location_scale_test_state(owner, ch);
    SEXP layout = PROTECT(Rf_mkString(state.layout.c_str()));
    SEXP counts = PROTECT(Rf_allocVector(INTSXP, 3));
    INTEGER(counts)[0] = state.L;
    INTEGER(counts)[1] = state.s;
    INTEGER(counts)[2] = state.K;
    SEXP active = PROTECT(Rf_ScalarLogical(state.active ? 1 : 0));
    SEXP x = PROTECT(real_vector(state.values));
    SEXP w = PROTECT(real_vector(state.weights));
    SEXP sigma = PROTECT(Rf_ScalarReal(state.sigma));
    SEXP children = PROTECT(real_vector(state.children));
    SEXP prior = PROTECT(real_vector({state.prior_nodes, state.prior_move}));
    SEXP out = named_list({layout, counts, active, x, w, sigma, children, prior},
                          {"layout", "counts", "active", "values", "weights", "sigma", "children", "prior"});
    UNPROTECT(8);
    return out;
  }
  const std::vector<double> numbers = real_values(values, "values");
  if (action == "set") {
    jags::RoBMA::location_scale_test_set(owner, ch, numbers);
    return R_NilValue;
  }
  if (action == "transform") {
    if (numbers.size() != 1) Rf_error("'values' must be one number.");
    return Rf_ScalarLogical(jags::RoBMA::location_scale_test_transform(owner, ch, numbers[0]) ? 1 : 0);
  }
  if (action == "density") return real_vector(jags::RoBMA::location_scale_test_density(owner, ch, numbers));
  if (action == "slice") {
    if (numbers.size() < 2) Rf_error("'values' must hold the bounds and the starts.");
    const std::vector<double> starts(numbers.begin() + 2, numbers.end());
    GetRNGstate();
    const std::vector<double> out =
      jags::RoBMA::location_scale_test_slice(numbers[0], numbers[1], starts, owner == 1, &unif_rand, &exp_rand);
    PutRNGstate();
    return real_vector(out);
  }
  Rf_error("Unknown location scale test action.");
  return R_NilValue;
  ROBMA_NATIVE_END
}
