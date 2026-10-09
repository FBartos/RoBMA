context("Hypothesis statements share a group exactly when BayesTools accepts one linear target")

# The Bayes factors of the rows of one call are those of the statements run
# one at a time.
.review_expect_combined_equals_single <- function(object, statements, info) {

  call <- function(statement) {
    hypothesis(object, statement, density_method = "KDE", columns = "all", seed = 1)
  }
  combined <- tryCatch(call(statements), error = identity)
  expect_s3_class(combined, "data.frame")
  if (!is.data.frame(combined)) {
    return(invisible(NULL))
  }
  expect_equal(nrow(combined), length(statements), info = info)
  for (i in seq_along(statements)) {
    expect_equal(
      as.numeric(combined[["BF"]][[i]]), as.numeric(call(statements[[i]])[["BF"]]),
      tolerance = 1e-8, info = paste(info, statements[[i]])
    )
  }
  invisible(combined)
}

test_that("different level expressions with equal fitted weights are separate groups", {

  # g[A] - g[B] and 2 * g[A] + g[C] have the same fitted weights in the
  # sum-to-zero parameterization but are different linear combinations of levels
  three <- review_factor_replay(
    c("A", "B", "C"), "meandif", variable_intercept = FALSE, independent = TRUE
  )
  seed_pair <- c("g[A] - g[B] = 0.1", "2 * g[A] + g[C] = 0.1")
  # for two levels the weights are even identical
  two <- review_factor_replay(
    c("A", "B"), "meandif", variable_intercept = FALSE, independent = TRUE
  )
  exact_pair <- c("g[A] - g[B] = 0.1", "2 * g[A] = 0.1")

  groups <- function(object, statements) {
    vapply(.hypothesis_plans(object, statements), `[[`, character(1), "group")
  }
  expect_length(unique(groups(three, seed_pair)), 2L)
  expect_length(unique(groups(two, exact_pair)), 2L)
  # one level expression, two thresholds: one group
  expect_length(unique(groups(three, c("g[A] - g[B] = 0.1", "g[A] - g[B] = 0.2"))), 1L)

  .review_expect_combined_equals_single(three, seed_pair, "three levels")
  .review_expect_combined_equals_single(two, exact_pair, "two levels")
  same <- .review_expect_combined_equals_single(
    three, c("g[A] - g[B] = 0.1", "g[A] - g[B] = 0.2"), "same expression"
  )
  expect_equal(nrow(same), 2L)

  # the marginal means of the same model
  .review_expect_combined_equals_single(marginal_means(three), seed_pair,
                                        "three levels, marginal means")
  .review_expect_combined_equals_single(marginal_means(two), exact_pair,
                                        "two levels, marginal means")
})

# The group key of a level combination reads the level coefficients of a
# statement side from this unexported BayesTools helper
# (.hypothesis_plan_combination_group). A rename or a changed return value
# there must fail here, with the reason, instead of in every grouped call.
test_that("the BayesTools linear-target side helper the group key uses is unchanged", {

  side_function <- tryCatch(
    utils::getFromNamespace(".hypothesis_linear_target_side", "BayesTools"),
    error = identity
  )
  expect_true(
    is.function(side_function),
    info = paste(
      "BayesTools:::.hypothesis_linear_target_side is gone; RoBMA's",
      ".hypothesis_plan_combination_group() depends on it."
    )
  )
  if (!is.function(side_function)) {
    return(invisible(NULL))
  }

  coefficients <- function(statement) {
    hypothesis <- BayesTools::hypothesis_parse(statement)
    side <- side_function(hypothesis[["statements"]][[1L]][["left"]], "g")
    expect_true(
      is.list(side) && "coefficients" %in% names(side),
      info = paste("no 'coefficients' returned for", statement)
    )
    side[["coefficients"]]
  }
  difference <- coefficients("g[A] - g[B] = 0.1")
  expect_true(is.numeric(difference) && !is.null(names(difference)))
  expect_identical(names(difference), c("A", "B"))
  expect_equal(unname(difference[["A"]]), -unname(difference[["B"]]))
  # the key does not depend on the threshold, and differs between expressions
  expect_identical(coefficients("g[A] - g[B] = 0.2"), difference)
  weighted <- coefficients("2 * g[A] + g[C] = 0.1")
  expect_identical(names(weighted), c("A", "C"))
  expect_equal(unname(weighted[["A"]]), 2 * unname(weighted[["C"]]))
})
