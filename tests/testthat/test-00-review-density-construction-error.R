context("Density construction failures on public density-only routes")

test_that("density-only routes signal a public class with the KDE remedy", {

  object <- review_factor_replay(c("A", "B", "C"), "meandif")
  means  <- marginal_means(object, bf = FALSE)
  stopper <- function(...) {
    .iwmde_stop_construction_failure(
      "q_grid_cmde", "mu_g[A]", c(3L, 7L), "normalization", "forced failure"
    )
  }
  capture <- function(expr) {
    tryCatch(
      withCallingHandlers(expr, warning = function(w) invokeRestart("muffleWarning")),
      error = identity
    )
  }
  grDevices::pdf(NULL)
  withr::defer(grDevices::dev.off())
  testthat::local_mocked_bindings(.iwmde_plan = stopper, .package = "RoBMA")

  errors <- list(
    marginal_means = capture(marginal_means(object, bf = FALSE, density_method = "qCMDE")),
    plot           = capture(plot(means, "g", density_method = "qCMDE"))
  )
  for (route in names(errors)) {
    error <- errors[[route]]
    expect_s3_class(error, "RoBMA_density_construction_error")
    # a member of the density condition family
    expect_s3_class(error, "RoBMA_density_error")
    # the internal class stays for the callers that match it
    expect_s3_class(error, "iwmde_construction_error")
    expect_s3_class(error, "error")
    if (inherits(error, "error")) {
      expect_match(conditionMessage(error), "construction failed for target 'mu_g[A]'",
                   fixed = TRUE, info = route)
      expect_match(conditionMessage(error), "density_method = 'KDE'", fixed = TRUE,
                   info = route)
    }
  }

  # the hypothesis ordinate route keeps its own class
  error <- capture(hypothesis(object, "g[A] = 0", density_method = "qCMDE"))
  expect_s3_class(error, "RoBMA_density_ordinate_error")
  expect_s3_class(error, "iwmde_construction_error")
  expect_false(inherits(error, "RoBMA_density_construction_error"))
})
