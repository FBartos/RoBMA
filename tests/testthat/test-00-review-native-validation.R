context("Native input validation regressions")

test_that("native GLMM outcomes are validated before likelihood arithmetic", {

  binomial <- list(3L, 2L, 20L, 22L, matrix(0), matrix(.1), NULL,
                   0, 0, matrix(0), matrix(0))
  poisson <- list(3L, 2L, 20, 22, matrix(0), matrix(.1), NULL,
                  0, 0, matrix(0), matrix(0))
  call_native <- function(symbol, args) {

    do.call(.Call, c(list(symbol), args, list(PACKAGE = "RoBMA")))
  }
  for (row_sum in c(FALSE, TRUE)) {
    suffix <- if (row_sum) "_row_sum" else ""
    bin_name <- paste0("RoBMA_glmm_binom_marginal_loglik", suffix)
    pois_name <- paste0("RoBMA_glmm_pois_marginal_loglik", suffix)
    expect_true(all(is.finite(call_native(bin_name, binomial))))
    expect_true(all(is.finite(call_native(pois_name, poisson))))
    for (bad in c(NA_integer_, -1L, 21L)) {
      args <- binomial
      args[[1L]] <- bad
      expect_error(call_native(bin_name, args), "Binomial counts")
    }
    for (bad in c(NA_integer_, -1L)) {
      args <- poisson
      args[[1L]] <- bad
      expect_error(call_native(pois_name, args), "Poisson counts")
    }
    for (bad in c(NA_real_, NaN, Inf, 0, -1)) {
      args <- poisson
      args[[3L]] <- bad
      expect_error(call_native(pois_name, args), "Poisson exposures")
    }
  }
  conditional_binomial <- c(binomial[1:5], list(matrix(0), NULL))
  conditional_poisson <- c(poisson[1:5], list(matrix(0), NULL))
  conditional_binomial[[1L]] <- NA_integer_
  conditional_poisson[[3L]] <- 0
  expect_error(call_native("RoBMA_glmm_binom_conditional_loglik_sum",
                          conditional_binomial), "Binomial counts")
  expect_error(call_native("RoBMA_glmm_pois_conditional_loglik_sum",
                          conditional_poisson), "Poisson exposures")
})
