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
    if (std::strcmp(mode, "rank_one_poll") == 0) {
      const int rows = 64;
      std::vector<double> values(rows, 0.0), means(rows), sd(rows), se(rows), loading(rows);
      std::vector<int> bins(rows, 1);
      for (int row = 0; row < rows; ++row) {
        const int j = row + 1;
        const double x = (j - (rows + 1) / 2.0) / rows;
        means[row] = 2 * x + .15 * std::sin(j);
        sd[row] = .03 + .02 * (j % 4);
        se[row] = .5 + .25 * (j % 5);
        loading[row] = (j % 3 == 1 ? -1.0 : 1.0) * (.5 + .1 * (j % 7));
      }
      const double omega[] = {.4, 1.0, .4};
      const double lower[] = {1.0, -1.0, -std::numeric_limits<double>::infinity()};
      const double upper[] = {std::numeric_limits<double>::infinity(), 1.0, -1.0};
      double mcse, change;
      const auto failure = +[]() { throw std::runtime_error("native rank-one poll"); };
      // No QMC data are read on rank-one PRODUCT fallback. The unchanged
      // scalar-work poll threshold must throw through the real entry boundary.
      cpp_selnorm_factor_step_lpdf(values.data(), means.data(), sd.data(), loading.data(),
        rows, 1, se.data(), omega, 3, lower, upper, bins.data(), 1, true, SELKERNEL_STEP,
        nullptr, nullptr, nullptr, 0, nullptr, 4, 8192, 8, 1e-12,
        &mcse, &change, nullptr, SELVECTOR_PRODUCT, nullptr, nullptr, nullptr, failure);
      Rf_error("The rank-one work poll was not reached.");
    }
    if (std::strncmp(mode, "projection_", 11) == 0) {
      std::vector<double> grid(64), means(1024, 0.0), logs(1024, -std::log(1024.0));
      std::vector<long double> density(64);
      for (int i = 0; i < 64; ++i) grid[i] = (i - 32) / 8.0;
      const bool disabled = std::strcmp(mode, "projection_disabled") == 0;
      const bool checked = std::strcmp(mode, "projection_checked") == 0;
      const bool constructor = std::strcmp(mode, "projection_constructor") == 0;
      const bool zero = std::strcmp(mode, "projection_zero") == 0;
      if (zero) std::fill(logs.begin(), logs.end(), -std::numeric_limits<double>::infinity());
      const auto quiet = +[]() {};
      const auto failure = +[]() { throw std::runtime_error("native projection poll"); };
      // No hook argument at all on the default-disabled path. The other
      // modes exercise actual construction, recurrence and zero components.
      std::unique_ptr<SelNormNormalProjectionWorkspace> projection;
      if (disabled) projection.reset(new SelNormNormalProjectionWorkspace(grid.data(), 64));
      else projection.reset(new SelNormNormalProjectionWorkspace(grid.data(), 64,
        checked ? quiet : failure, constructor || checked ? 1 : 512));
      if (!projection->density(means.data(), 1024, 1.0, 1.0, logs.data(), density.data())) {
        Rf_error("Native projection probe could not evaluate its regular grid.");
      }
      SEXP out = PROTECT(Rf_allocVector(REALSXP, 64));
      for (int i = 0; i < 64; ++i) REAL(out)[i] = static_cast<double>(density[i]);
      UNPROTECT(1);
      return out;
    }
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
