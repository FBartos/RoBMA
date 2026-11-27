# Platform properties behind exact numerical checks. Where long double is wider
# than double, R's accumulators (sum, rowSums, the internal matrix product, numeric
# literals) and the native long double products round in extended precision, so
# many results are exact to the last bit. Apple silicon has no such type: there
# they round in double and these results differ in the last bits.
.test_long_double_extended <- function() .Machine$sizeof.longdouble > 8

# Whether the package's native code fuses a * b + c into one rounding. Compilers
# do where FMA is part of the instruction set they target: always on aarch64, and
# on x86-64 from -march=x86-64-v3 (the GCC default on the Ubuntu 26.04 runners).
# The native probe reports what this build does.
.test_fused_multiply_add <- function() {
  .Call("RoBMA_native_unwind_probe", "fused_multiply_add", NULL, 1L, PACKAGE = "RoBMA")
}

# Bit for bit where long double is extended, otherwise within an absolute bound
# measured at the call site.
.test_expect_exact_or_within <- function(object, expected, bound, info = NULL) {

  if (.test_long_double_extended()) {
    testthat::expect_identical(object, expected, info = info)
  } else {
    testthat::expect_lte(max(abs(object - expected)), bound)
  }
}

# covr (>= 4.0.0) instruments package code with '.Call(covr:::covr_count, id)'.
# A test that puts its own '.Call' into a package function's environment gets
# those counting calls too: pass them, which name a native symbol rather than a
# string, to the real '.Call' and the package's routines to 'handler'.
.native_call_mock <- function(handler) {

  force(handler)
  function(.NAME, ...) {
    if (!is.character(.NAME)) {
      return(base::.Call(.NAME, ...))
    }
    handler(.NAME, ...)
  }
}

# Whether the tests run under covr, which rewrites package function bodies.
.test_in_covr <- function() identical(Sys.getenv("R_COVR"), "true")

# Whether seeded draws are those of the build that made the locked references:
# not a build that fuses multiply-adds, not covr's unoptimized build. The
# location block sampler and the selection samplers added with it
# (RoBMA::LocationBlock, ::GatedCorrectedBlock, ::ScreenedFinite) turn last-bit
# differences into different chains, so the same seed gives different, equally
# valid draws on those builds (RoBMA, BayesTools and JAGS fuse on the Linux and
# macOS runners): there the short fit dat.lehmann2018_RoBMA gave BF01 = 4.495
# against the reference 3.435. With `fit`, also TRUE for a fit none of them
# samples.
.test_reference_build_draws <- function(fit = NULL) {

  if (!.test_in_covr() && !.test_fused_multiply_add()) return(TRUE)
  if (is.null(fit)) return(FALSE)
  jags     <- fit[["fit"]]
  tables   <- c(list(jags[["samplers"]]), lapply(jags[["samplers_by_worker"]], `[[`, "samplers"))
  samplers <- unlist(lapply(tables, `[[`, "Sampler"))
  return(!any(c("RoBMA::LocationBlock", "RoBMA::GatedCorrectedBlock", "RoBMA::ScreenedFinite") %in% samplers))
}

# Skips, after the test's other expectations, what was not compared because
# .test_reference_build_draws() is FALSE (`what` names the fits or fixtures).
.skip_reference_build_draws <- function(what) {

  if (length(what) == 0L) return(invisible(FALSE))
  if (exists(".announce_existing_visual_snapshots", mode = "function")) {
    .announce_existing_visual_snapshots()
  }
  testthat::skip(paste0(
    "Seeded draws of the location block and its selection samplers differ from the reference build's ",
    "on this build (fused multiply-adds or covr); not checked here: ", paste(unique(what), collapse = ", "), "."
  ))
}
