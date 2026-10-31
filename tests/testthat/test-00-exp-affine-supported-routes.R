source(testthat::test_path("common-functions.R"))

test_that("exp-affine coefficients refuse qCMDE/IWMDE with the original-scale classes", {

  skip_on_cran()
  fit <- load_fit("lane_exp_affine")
  plan <- engine_plans(fit, "intercept = 0.5", component = "scale")[[1L]]
  expect_identical(plan[["kind"]], "exp_affine")

  classes <- c("RoBMA_density_method_original_scale",
               "RoBMA_density_method_unavailable")
  for (method in c("qCMDE", "IWMDE")) {
    error <- tryCatch(
      plot(fit, "intercept", component = "scale", density_method = method,
           density_control = list(n_points = 20, samples = 50)),
      error = identity
    )
    expect_identical(class(error), c(classes, "error", "condition"), info = method)
    expect_identical(
      conditionMessage(error),
      paste0(
        "qCMDE/IWMDE does not support the fitted nonlinear joint transform ",
        "for 'log_tau_intercept'. Use density_method = 'KDE' or ",
        "standardized_coefficients = TRUE."
      ),
      info = method
    )
    error <- tryCatch(
      hypothesis(fit, "intercept = 0.5", component = "scale",
                 density_method = method),
      error = identity
    )
    expect_identical(
      class(error),
      c("BayesTools_hypothesis_method", "BayesTools_hypothesis_unavailable", classes,
        "error", "condition"),
      info = method
    )
  }
  # The normal approximation is no qCMDE/IWMDE request: its refusal keeps the
  # hypothesis classes only.
  error <- tryCatch(
    hypothesis(fit, "intercept = 0.5", component = "scale",
               density_method = "normal"),
    error = identity
  )
  expect_identical(
    class(error),
    c("BayesTools_hypothesis_method", "BayesTools_hypothesis_unavailable", "error",
      "condition")
  )
  kde <- suppressWarnings(hypothesis(
    fit, "intercept = 0.5", component = "scale", density_method = "KDE"
  ))
  expect_true(is.finite(attr(kde, "raw_BF")))
})
