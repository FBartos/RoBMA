test_that("interpret includes multivariate heterogeneity inclusion evidence", {

  inclusion <- function(parameters, BF, bounds = rep(NA_character_, length(BF))) {

    out <- data.frame(
      prior_prob = .5, post_prob = BF / (1 + BF), inclusion_BF = BF,
      row.names = parameters
    )
    attr(BF, "bound_operator") <- bounds
    out[["inclusion_BF"]] <- BayesTools::format_BF(BF, inclusion = TRUE)
    class(out) <- c("BayesTools_table", "data.frame")
    attr(out, "type") <- c("prior_prob", "post_prob", "inclusion_BF")
    attr(out, "parameters") <- parameters
    out
  }
  summary_object <- list(
    name = "Bayesian Multivariate Model-Averaged Meta-Analysis",
    inclusion_components = inclusion(c("Effect", "Publication Bias"), c(2, .5)),
    inclusion_random = inclusion(c("Study: tau", "tau_total"), c(9, 100), c(NA, ">"))
  )
  original <- summary_object
  object <- structure(
    list(fit = list(1), data = structure(list(), measure = "GEN")),
    class = c("RoBMA", "brma")
  )
  testthat::local_mocked_bindings(
    summary.brma = function(...) summary_object,
    .package = "RoBMA"
  )

  output <- interpret(object, scope = "components")
  records <- attr(output, "records")
  evidence <- records[records[["kind"]] == "evidence", , drop = FALSE]
  expect_identical(evidence[["row"]], c(
    "Effect", "Heterogeneity: Study: tau", "Heterogeneity: tau_total", "Publication Bias"
  ))
  expect_equal(evidence[["BF_canonical_value"]], c(2, 9, 100, .5))
  expect_identical(evidence[["BF_canonical_bound_operator"]], c(NA, NA, ">", NA))
  expect_true(any(grepl("Heterogeneity: Study: tau inclusion", output, fixed = TRUE)))
  expect_true(any(grepl("BF > 100.000", output, fixed = TRUE)))
  expect_identical(summary_object, original)

  summary_object[["inclusion_components"]] <- list()
  random_only <- attr(interpret(object, scope = "components"), "records")
  expect_equal(sum(random_only[["kind"]] == "evidence"), 2L)

  bias_only <- attr(interpret(object, scope = "bias"), "records")
  expect_equal(sum(bias_only[["kind"]] == "evidence"), 0L)
})
