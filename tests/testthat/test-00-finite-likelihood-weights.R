test_that("shared outcome inputs reject infinite weights and preserve fractional weights", {

  constructors <- list(
    normal = function(weights) brma.norm(
      yi = c(.1, .2, .3), vi = c(.04, .09, .16), weights = weights,
      measure = "GEN", only_data = TRUE),
    binomial = function(weights) brma.glmm(
      ai = c(1, 2, 3), ci = c(2, 3, 4), n1i = rep(10, 3), n2i = rep(12, 3),
      weights = weights, measure = "OR", only_data = TRUE),
    poisson = function(weights) brma.glmm(
      x1i = c(1, 2, 3), x2i = c(2, 3, 4), t1i = rep(10, 3), t2i = rep(12, 3),
      weights = weights, measure = "IRR", only_data = TRUE)
  )
  for (construct in constructors) {
    error <- tryCatch(construct(c(1, Inf, 1)), error = identity)
    expect_identical(conditionMessage(error),
                     "The 'weights' argument must contain only finite values.")
    expect_null(conditionCall(error))
    # Existing positivity and missing-value diagnostics retain their precedence.
    expect_error(construct(c(1, -Inf, 1)),
                 "The 'weights' must be higher than 0.", fixed = TRUE)
    for (missing in c(NA_real_, NaN)) {
      expect_error(construct(c(1, missing, 1)),
                   "The 'weights' argument must not contain missing values.", fixed = TRUE)
    }
    for (weights in list(c(1, 2, 3), c(.5, 1.25, 2))) {
      object <- construct(weights)
      expect_identical(object[["data"]][["outcome"]][["weights"]], weights)
    }
  }
})

test_that("missing likelihood weights keep their existing pre-filter rejection", {

  expect_error(brma.norm(
    yi = c(.1, NA_real_, .3), vi = c(.04, .09, .16),
    weights = c(1, NA_real_, 1), measure = "GEN", only_data = TRUE
  ), "The 'weights' argument must not contain missing values.", fixed = TRUE)
})
