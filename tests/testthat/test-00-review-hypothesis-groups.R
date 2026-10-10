context("Hypothesis statements share a group exactly when BayesTools accepts one linear target")

source(testthat::test_path("common-functions.R"))

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

# Statements on the marginal means of different parameters are tested in one
# call, each on the parameter it names: the rows are those of the statements
# run alone, in the order of the statements. The fit has two factor
# moderators, g1 (levels 5, 10, 20) and g2 (levels 1 to 4).
test_that("marginal-means statements on several parameters equal the statements run alone", {

  skip_on_cran()
  means <- marginal_means(load_fit("hypothesis_plan_treatment"), bf = FALSE)
  g1_10 <- "g1[10] > 0 vs g1[10] = 0"
  g1_20 <- "g1[20] > 0 vs g1[20] = 0"
  g2_2  <- "g2[2] > 0 vs g2[2] = 0"

  forward <- .review_expect_combined_equals_single(means, c(g1_10, g2_2), "g1, g2")
  expect_identical(rownames(forward), c("mu_g1", "mu_g2"))
  reversed <- .review_expect_combined_equals_single(means, c(g2_2, g1_10), "g2, g1")
  expect_identical(rownames(reversed), c("mu_g2", "mu_g1"))
  interleaved <- .review_expect_combined_equals_single(
    means, c(g1_10, g2_2, g1_20), "g1, g2, g1"
  )
  expect_identical(rownames(interleaved), c("mu_g1 (1)", "mu_g2", "mu_g1 (3)"))
  # A linear combination of one parameter's levels next to a level of another
  # parameter, and region statements.
  .review_expect_combined_equals_single(
    means, c("g1[10] - g1[20] > 0", g2_2), "combination of g1, level of g2"
  )
  .review_expect_combined_equals_single(
    means, c("g2[2] > g2[3]", "g1[10] > 0", g2_2), "regions and points"
  )

  # qCMDE, the default density method of a fitted brma, gives the rows of the
  # statements run alone as well, also for points on one parameter next to a
  # region on another.
  requests <- list(
    c(g1_10, g2_2),
    c(g2_2, g1_10, g1_20),
    c("g1[10] = 0.1", "g2[2] > 0", "g2[3] = 0.1")
  )
  for (statements in requests) {
    combined <- hypothesis(
      means, statements, density_method = "qCMDE", columns = "all", seed = 1
    )
    alone <- lapply(statements, function(statement) {
      hypothesis(
        means, statement, density_method = "qCMDE", columns = "all", seed = 1
      )
    })
    expect_identical(
      combined[["method"]],
      unlist(lapply(alone, `[[`, "method"), use.names = FALSE)
    )
    for (column in c("prior", "posterior")) {
      expect_equal(
        combined[[column]],
        unlist(lapply(alone, `[[`, column), use.names = FALSE),
        tolerance = 1e-12
      )
    }
    expect_equal(
      attr(combined, "raw_log_BF"),
      unlist(lapply(alone, attr, "raw_log_BF"), use.names = FALSE),
      tolerance = 1e-12
    )
    expect_length(attr(combined, "numerical_diagnostics"), length(statements))
  }
})


test_that("a marginal-means parameter is refused for statements of another parameter", {

  skip_on_cran()
  means <- marginal_means(load_fit("hypothesis_plan_treatment"), bf = FALSE)
  mismatch <- c(
    "RoBMA_hypothesis_parameter_mismatch", "RoBMA_hypothesis_statement",
    "error", "condition"
  )
  requests <- list(
    list(c("g1[10] > 0 vs g1[10] = 0", "g2[2] > 0 vs g2[2] = 0"), "g1"),
    list("g2[2] > 0 vs g2[2] = 0", "g1"),
    list("g2[2] > 0", "g1"),
    list("g2 > 0", "mu_g1"),
    list("g1[10] > 0", "g2"),
    list("g1[10] > g2[2]", "g1")
  )
  for (request in requests) {
    error <- tryCatch(
      hypothesis(means, request[[1L]], parameter = request[[2L]], seed = 1),
      error = identity
    )
    expect_identical(
      class(error), mismatch,
      info = paste(request[[1L]], request[[2L]])
    )
    expect_match(
      conditionMessage(error),
      "' refers to the marginal means of 'mu_g[12]', but 'parameter' selects 'mu_g[12]'. ",
      info = paste(request[[1L]], request[[2L]])
    )
    expect_match(
      conditionMessage(error),
      "omit 'parameter' to test statements on different parameters in one call",
      fixed = TRUE
    )
    expect_false(grepl("prior density", conditionMessage(error), fixed = TRUE))
  }
  expect_match(
    tryCatch(
      hypothesis(means, "g2[2] > 0 vs g2[2] = 0", parameter = "g1"),
      error = conditionMessage
    ),
    "Statement '`g2[2]` > 0 vs `g2[2]` = 0' refers to the marginal means of 'mu_g2', but 'parameter' selects 'mu_g1'.",
    fixed = TRUE
  )

  # The statements of the selected parameter, and unknown names, as before.
  statement <- "g1[10] > 0 vs g1[10] = 0"
  expect_identical(
    hypothesis(means, statement, parameter = "g1", density_method = "KDE", seed = 1),
    hypothesis(means, statement, density_method = "KDE", seed = 1)
  )
  expect_error(
    hypothesis(means, "zz[10] > 0", parameter = "g1"),
    class = "BayesTools_parameter_not_found"
  )
  expect_error(
    hypothesis(means, c("g1[10] > 0", "zz[10] > 0")),
    class = "BayesTools_parameter_not_found"
  )
})


test_that("one marginal-means statement names the marginal means of one parameter", {

  skip_on_cran()
  means <- marginal_means(load_fit("hypothesis_plan_treatment"), bf = FALSE)
  for (statement in c("g1[10] > g2[2]", "g1[10] > g2[2] vs g1[10] = g2[2]")) {
    error <- tryCatch(hypothesis(means, statement), error = identity)
    expect_identical(
      class(error),
      c("RoBMA_hypothesis_ambiguous", "RoBMA_hypothesis_statement",
        "error", "condition"),
      info = statement
    )
    expect_match(
      conditionMessage(error),
      "one statement must refer to the marginal means of one parameter",
      fixed = TRUE
    )
    expect_false(grepl("Specify 'parameter'", conditionMessage(error), fixed = TRUE))
  }
})
