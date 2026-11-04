context("Native selection event regressions")

test_that("rank-one two-sided best events split both threshold crossings", {

  evaluate <- function(bounds, weights, mean = 0, loading = 1,
                       lower = -Inf, upper = Inf) {

    .Call("RoBMA_selnorm_gaussian_event_mass_batch",
      matrix(mean, 1L), NULL, 1, matrix(weights, 1L),
      tail(bounds, -1L), head(bounds, -1L), 1L, 1L, 2L,
      matrix(lower, 1L), matrix(upper, 1L), numeric(), 8L, 2L,
      matrix(loading, 1L), NULL, .005, PACKAGE = "RoBMA")
  }
  for (mean in c(-.4, 0, .7)) {
    for (loading in c(-1.2, .8)) {
      # The event is |Y| >= 1. The reference integrates its two tails
      # directly, independent of the native latent-line partition.
      outside <- pnorm(-1, mean, abs(loading)) +
        pnorm(1, mean, abs(loading), lower.tail = FALSE)
      expected <- .2 + .8 * outside
      for (bounds in list(c(Inf, 1, -Inf), c(Inf, 1, 0, -1, -Inf))) {
        weights <- if (length(bounds) == 3L) c(1, .2) else c(1, .2, .2, 1)
        actual <- evaluate(bounds, weights, mean, loading)
        expect_equal(exp(actual$log_mass), expected, tolerance = 1e-13)
        expect_identical(actual$relative_mcse, 0)
      }
      negative <- evaluate(c(Inf, 1, -Inf), c(1, .2), mean, loading,
                           upper = 0)
      expected_negative <- .2 * pnorm(0, mean, abs(loading)) +
        .8 * pnorm(-1, mean, abs(loading))
      expect_equal(exp(negative$log_mass), expected_negative, tolerance = 1e-13)
    }
  }
})

