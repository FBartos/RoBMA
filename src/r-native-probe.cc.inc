// Registered development probes exercise the same boundary as real kernels.
// Their counters belong only to probe-owned workspaces, not model evaluation.
namespace {
std::atomic<int> native_probe_alive{0};
std::atomic<int> native_probe_destroyed{0};

struct NativeProbeWorkspace {
  std::vector<double> values;
  NativeProbeWorkspace() : values(2048, 1.0) { ++native_probe_alive; }
  ~NativeProbeWorkspace() { --native_probe_alive; ++native_probe_destroyed; }
};
}

extern "C" SEXP RoBMA_native_unwind_probe(SEXP action, SEXP callback, SEXP threads)
{
  // Observe prior state before this entry installs its own context. These
  // nonallocating comparisons never dereference a possibly stale pointer.
  const bool prior_context = robma_native::active_context != nullptr;
  const bool prior_rng_owner = robma_native::rng_owner != nullptr;
  const bool prior_context_is_owner = prior_context &&
    robma_native::active_context == robma_native::rng_owner;
  ROBMA_NATIVE_BEGIN(action, callback, threads)
  if (TYPEOF(action) != STRSXP || XLENGTH(action) != 1 ||
      TYPEOF(threads) != INTSXP || XLENGTH(threads) != 1 ||
      INTEGER(threads)[0] < 1 || INTEGER(threads)[0] > 2) {
    Rf_error("Invalid native unwind probe arguments.");
  }
  const char *mode = CHAR(STRING_ELT(action, 0));
  if (std::strcmp(mode, "state") == 0) {
    SEXP state = PROTECT(Rf_allocVector(LGLSXP, 3));
    LOGICAL(state)[0] = prior_context;
    LOGICAL(state)[1] = prior_rng_owner;
    LOGICAL(state)[2] = prior_context_is_owner;
    UNPROTECT(1);
    return state;
  }
  if (std::strcmp(mode, "counts") == 0) {
    SEXP counts = PROTECT(Rf_allocVector(INTSXP, 2));
    INTEGER(counts)[0] = native_probe_alive.load();
    INTEGER(counts)[1] = native_probe_destroyed.load();
    UNPROTECT(1);
    return counts;
  }
  if (std::strcmp(mode, "openmp") == 0) {
#if defined(_OPENMP)
    return Rf_ScalarLogical(TRUE);
#else
    return Rf_ScalarLogical(FALSE);
#endif
  }
  NativeProbeWorkspace outer;
  {
    NativeProbeWorkspace inner;
    const bool nil_length = std::strcmp(mode, "nil_length") == 0;
    const bool nil_xlength = std::strcmp(mode, "nil_xlength") == 0;
    int nil_size = -1;
    if (nil_length || nil_xlength) {
      GetRNGstate();
      unif_rand();
      nil_size = nil_length ? Rf_length(R_NilValue) :
        static_cast<int>(XLENGTH(R_NilValue));
    }
    if (std::strcmp(mode, "cpp") == 0) throw std::runtime_error("native probe C++ failure");
    if (std::strcmp(mode, "bad_alloc") == 0) throw std::bad_alloc();
    if (std::strcmp(mode, "allocation") == 0) return Rf_allocVector(REALSXP, -1);
    if (std::strcmp(mode, "malformed_data") == 0) return Rf_ScalarReal(REAL(callback)[0]);
    if (std::strcmp(mode, "protect_stack") == 0) {
      for (int i = 0; i < 1000000; ++i) PROTECT(R_NilValue);
      UNPROTECT(1000000);
      Rf_error("Native probe did not reach the protection-stack limit.");
    }
    if (std::strcmp(mode, "warning") == 0) {
      robma_native::call(Rf_warning, "native probe warning");
    }
    if (std::strcmp(mode, "r_error") == 0) Rf_error("native probe validation failure");
    if (std::strcmp(mode, "altrep") == 0) {
      if (TYPEOF(callback) != REALSXP || XLENGTH(callback) != 16) Rf_error("Expected 16 numeric values.");
      SEXP out = PROTECT(Rf_allocVector(REALSXP, 16));
      robma_parallel::for_rows(16, INTEGER(threads)[0], 4, [&](int row) {
        REAL(out)[row] = REAL(callback)[row] + 1.0;
      });
      UNPROTECT(1);
      return out;
    }
    if (std::strcmp(mode, "rng_success") == 0) {
      GetRNGstate();
      unif_rand();
    }
    if (std::strcmp(mode, "callback") == 0 || std::strcmp(mode, "rng_callback") == 0 ||
        std::strcmp(mode, "time_limit") == 0 ||
        ((nil_length || nil_xlength) && callback != R_NilValue)) {
      if (std::strcmp(mode, "rng_callback") == 0) {
        GetRNGstate();
        unif_rand();
      }
      robma_native::Context *const parent = robma_native::active_context;
      robma_native::Context *const owner = robma_native::rng_owner;
      SEXP invocation = PROTECT(robma_native::call(Rf_lang1, callback));
      robma_native::call(Rf_eval, invocation, R_GlobalEnv);
      UNPROTECT(1);
      if (robma_native::active_context != parent || robma_native::rng_owner != owner) {
        Rf_error("Native probe callback did not restore its parent context and RNG owner.");
      }
      if (std::strcmp(mode, "time_limit") == 0) {
        // The R callback installs a short time limit. This bounded fallback
        // fails the test instead of hanging if interrupts stop being checked.
        const auto start = std::chrono::steady_clock::now();
        while (std::chrono::steady_clock::now() - start < std::chrono::seconds(2)) {
          R_CheckUserInterrupt();
        }
        Rf_error("Native interrupt probe did not observe its time limit.");
      }
    }
    if (nil_length || nil_xlength) return Rf_ScalarInteger(nil_size);
    if (std::strcmp(mode, "worker") == 0 || std::strcmp(mode, "worker_api") == 0) {
      const bool api = std::strcmp(mode, "worker_api") == 0;
      robma_parallel::for_rows(16, INTEGER(threads)[0], 4, [&](int row) {
        NativeProbeWorkspace worker;
        if (row == 2) {
          if (api) robma_native::call(&(Rf_allocVector), REALSXP, 1);
          else throw std::runtime_error("native probe worker failure");
        }
      });
    }
  }
  return Rf_ScalarInteger(42);
  ROBMA_NATIVE_END
}
