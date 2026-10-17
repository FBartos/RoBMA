.expect_random_prior_rejected <- function(expr) {

  reporter <- testthat::ListReporter$new()
  testthat::with_reporter(reporter, {
    reporter$start_file("random-prior-negative-control")
    testthat::test_that("the supplied curve is rejected", force(expr))
    reporter$end_file()
  })
  result <- as.data.frame(reporter$get_results())
  expect_gt(sum(result[["failed"]]), 0L)
  expect_false(any(result[["error"]] | result[["skipped"]]))
}

test_that("random prior comparisons retain and check every actual line ordinate", {

  x <- seq_len(200L) / 201
  reference <- rep(1 / 3, length(x))
  line <- function(y) {
    ggplot2::ggplot(data.frame(x = x, y = y), ggplot2::aes(x = x, y = y)) +
      ggplot2::geom_line()
  }
  curve <- .random_prior_curve(line(reference))
  expect_equal(curve, data.frame(x = x, y = reference))
  .expect_random_prior_curve(curve, reference, 1e-6)
  for (invalid in c(0, -1, NA_real_)) {
    values <- reference
    values[[50L]] <- invalid
    curve <- .random_prior_curve(line(values))
    expect_equal(nrow(curve), length(x))
    .expect_random_prior_rejected(.expect_random_prior_curve(curve, reference, 1e-6))
  }
  .expect_random_prior_rejected(.expect_random_prior_curve(.random_prior_curve(line(rep(0, length(x)))), reference, 1e-6))
  .expect_random_prior_rejected(.expect_random_prior_curve(curve, reference[-1L], 1e-6))
  .expect_random_prior_rejected(.expect_random_prior_curve(curve[FALSE, ], numeric(), 1e-6))
})

test_that("random prior comparisons check declared boundaries and outside zeros", {

  x <- c(-.1, 0, seq_len(200L) / 201, 1, 1.1)
  expected <- stats::dbeta(x, 2, 3)
  .expect_random_prior_curve(data.frame(x = x, y = expected), expected, 1e-6,
    support = c(0, 1))
  x <- c(0, 0, seq_len(200L) / 201, 1, 1)
  y <- c(0, rep(1 / 3, length(x) - 2L), 0)
  .expect_random_prior_curve(data.frame(x = x, y = y), rep(1 / 3, length(x)),
    1e-6, support = c(0, 1))
  y[[1L]] <- -1
  .expect_random_prior_rejected(.expect_random_prior_curve(data.frame(x = x, y = y),
    rep(1 / 3, length(x)), 1e-6, support = c(0, 1)))
  x <- c(0, seq_len(200L) / 201)
  expected <- stats::dbeta(x, .5, 2)
  .expect_random_prior_curve(data.frame(x = x, y = expected), expected, 1e-6,
    support = c(0, 1))
})

test_that("plan assertions distinguish numerical refusals from construction errors", {

  expressions <- as.list(parse(testthat::test_path("helper-contracts.R")))
  helper <- Filter(function(node) is.call(node) && identical(node[[1L]], as.name("<-")) &&
    identical(node[[2L]], as.name(".expect_plans_consistent")), expressions)[[1L]]
  assertion <- refusal <- NULL
  visit <- function(node) {
    if (missing(node)) return(invisible(NULL))
    if (is.call(node)) {
      if (identical(node[[1L]], as.name("expect_false"))) assertion <<- node
      if (identical(node[[1L]], as.name("<-")) && identical(node[[2L]], as.name("permitted_refusal"))) refusal <<- node
      for (child in as.list(node)) visit(child)
    }
  }
  visit(helper)
  expect_false(is.null(assertion) || is.null(refusal))
  condition <- function(classes) structure(list(message = "explicit synthetic plan condition"),
    class = c(classes, "error", "condition"))
  controls <- list(numeric = 1, plain = simpleError("explicit synthetic plain error"),
    ordinate = condition("RoBMA_density_ordinate_error"),
    numerical = condition(c("RoBMA_density_ordinate_error", "iwmde_ordinate_numerical_error")),
    construction = condition(c("RoBMA_density_ordinate_error", "iwmde_construction_error")),
    chain = condition(c("RoBMA_density_ordinate_error", "RoBMA_density_chain_coverage", "iwmde_construction_error")))
  for (precomputed in c(TRUE, FALSE)) for (name in names(controls)) {
    input <- list2env(list(out = controls[[name]], precomputed = precomputed,
      case = paste("synthetic", name), expect_false = testthat::expect_false))
    eval(refusal, input)
    reporter <- testthat::ListReporter$new()
    testthat::with_reporter(reporter, {
      reporter$start_file("plan-condition-control")
      testthat::test_that("the actual plan assertion evaluates the synthetic condition", eval(assertion, input))
      reporter$end_file()
    })
    result <- as.data.frame(reporter$get_results())
    rejected <- name != "numeric" && (!precomputed || name %in% c("plain", "construction"))
    expect_equal(sum(result[["failed"]]), as.integer(rejected), info = paste(precomputed, name))
    expect_false(any(result[["error"]] | result[["skipped"]]))
  }
})
