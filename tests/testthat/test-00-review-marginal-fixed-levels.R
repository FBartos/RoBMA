context("Marginal-means hypotheses compare variable levels with a fixed level")

# A three-level treatment factor with a fixed-zero intercept: the reference
# level 'a' is a fixed marginal mean (exactly zero), 'b' and 'c' are variable.
.review_fixed_reference_replay <- function() {

  review_factor_replay(
    c("a", "b", "c"), "treatment", variable_intercept = FALSE,
    independent = TRUE, distribution = "normal"
  )
}

# The posterior probability of the event of a region statement.
.review_event_probability <- function(object, statement) {

  result <- hypothesis(object, statement, density_method = "KDE", columns = "all")
  odds <- as.numeric(result[["posterior"]])

  odds / (1 + odds)
}

test_that("a comparison with a fixed marginal mean is evaluated like on the fitted object", {

  object <- .review_fixed_reference_replay()
  mm     <- marginal_means(object)
  levels <- mm[["inference"]][["averaged"]][["mu_g"]]
  expect_equal(max(abs(as.numeric(levels[["a"]]))), 0)

  # the events of the comparisons with the fixed level are those of the
  # variable level against zero, on the draws of the marginal means
  expected <- mean(as.numeric(levels[["b"]]) > as.numeric(levels[["a"]]))
  expect_gt(expected, 0)
  expect_lt(expected, 1)
  for (statement in c("g[b] > g[a]", "g[b] - g[a] > 0", "g[b] > 0")) {
    observed <- tryCatch(.review_event_probability(mm, statement), error = identity)
    expect_type(observed, "double")
    if (is.numeric(observed)) {
      expect_equal(observed, expected, tolerance = 1e-12, info = statement)
    }
  }
  # the comparison of two variable levels is unaffected
  expect_equal(
    .review_event_probability(mm, "g[b] > g[c]"),
    mean(as.numeric(levels[["b"]]) > as.numeric(levels[["c"]])),
    tolerance = 1e-12
  )
  # the fitted object evaluates the same comparison on its own draws
  expect_equal(
    .review_event_probability(object, "g[b] > g[a]"), expected, tolerance = 1e-12
  )
})

test_that("tests of the fixed marginal mean itself stay refused", {

  mm <- marginal_means(.review_fixed_reference_replay())

  for (statement in c("g[a] > 0", "g[a] = 0", "g > 0", "g = 0")) {
    expect_error(
      hypothesis(mm, statement, density_method = "KDE"),
      class = "BayesTools_hypothesis_fixed", info = statement
    )
  }
})
