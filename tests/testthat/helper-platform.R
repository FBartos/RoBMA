# Platform properties behind exact numerical checks. Where long double is wider
# than double, R's accumulators (sum, rowSums, the internal matrix product, numeric
# literals) and the native long double products round in extended precision, so
# many results are exact to the last bit. Apple silicon has no such type: there
# they round in double and these results differ in the last bits.
.test_long_double_extended <- function() .Machine$sizeof.longdouble > 8

# Compilers fuse a * b + c in native code where FMA is part of the baseline
# instruction set (aarch64); x86-64 builds keep the two roundings.
.test_fused_multiply_add <- function() R.version$arch %in% c("aarch64", "arm64")

# Bit for bit where long double is extended, otherwise within an absolute bound
# measured at the call site.
.test_expect_exact_or_within <- function(object, expected, bound, info = NULL) {

  if (.test_long_double_extended()) {
    testthat::expect_identical(object, expected, info = info)
  } else {
    testthat::expect_lte(max(abs(object - expected)), bound)
  }
}
