test_that("conditional Poisson sums preserve zero and infinite rate limits", {

  skip_if_not(.has_native_glmm_row_sum("pois", conditional = TRUE))

  log_phi <- matrix(c(-Inf, 0, 750, Inf), ncol = 1L)
  mu <- matrix(0, nrow(log_phi), 1L)
  for (counts in list(c(0L, 0L), c(2L, 0L), c(0L, 3L))) {
    expected <- stats::dpois(counts[[1L]], exp(log_phi[, 1L]), log = TRUE) +
      stats::dpois(counts[[2L]], 2 * exp(log_phi[, 1L]), log = TRUE)
    for (weight in c(.5, 1, 2)) {
      actual <- .outcome_pdf_sum.pois_conditional(
        x1i = counts[[1L]], x2i = counts[[2L]], t1i = 1, t2i = 2,
        mu_samples = mu, log_phi = log_phi, weights = weight
      )
      expect_equal(actual, weight * expected, tolerance = 1e-12)
    }
  }
})
