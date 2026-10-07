test_that("variance point-region identities retain canonical logs and row diagnostics", {

  object <- gated_random_object()
  tau2 <- "`(mu) study: tau2(intercept)`"
  run <- function(statement) hypothesis(object, statement, density_method = "KDE",
    columns = "all", seed = 1, logBF = TRUE)
  point <- run(paste0(tau2, " = 0.09 vs ", tau2, " != 0.09"))
  region <- run(paste0(tau2, " > 0.09 vs ", tau2, " >= 0"))
  mixed <- run(paste0(tau2, " = 0.09 vs ", tau2, " > 0.09"))
  reversed <- run(paste0(tau2, " > 0.09 vs ", tau2, " = 0.09"))
  expected <- attr(point, "raw_log_BF") - attr(region, "raw_log_BF")
  expect_length(expected, 1L)
  expect_equal(attr(mixed, "raw_log_BF"), expected, tolerance = 1e-12)
  expect_equal(attr(reversed, "raw_log_BF"), -expected, tolerance = 1e-12)
  expect_equal(as.numeric(mixed$BF), expected, tolerance = 1e-12)
  expect_length(attr(mixed, "numerical_diagnostics"), 1L)
  expect_length(attr(mixed, "prior_numerical_diagnostics"), 1L)
})

test_that("parameter groups bind producer logs in statement order and export visible values", {

  statements <- c("a = 356 vs a != 356", "b > 0", "a = 0")
  make_table <- function(parameter, hypothesis) BayesTools::hypothesis_BF(
    c(1, 2, 3), BayesTools::prior("normal", list(0, 1)), hypothesis,
    parameter = parameter, density_method = "normal", columns = "all")
  first <- make_table("a", statements[c(1L, 3L)])
  second <- make_table("b", statements[2L])
  out <- .hypothesis_brma_bind_parameter_results(list(first, second),
    list(c(1L, 3L), 2L), BayesTools::hypothesis_parse(statements))
  expected <- c(attr(first, "raw_log_BF")[1L], attr(second, "raw_log_BF"),
    attr(first, "raw_log_BF")[2L])
  expect_equal(expected[1L], 710, tolerance = 1e-12)
  expect_identical(attr(out, "raw_log_BF"), expected)
  expect_length(attr(out, "numerical_diagnostics"), 3L)
  expect_length(attr(out, "prior_numerical_diagnostics"), 3L)
  inverse <- update(out, BF01 = TRUE)
  expect_equal(as.numeric(inverse$BF)[1L], exp(-710), tolerance = 1e-14)
  plain <- .output_table_as_long_data_frame(out, "location")
  bound <- .output_bind_long_data_frames(list(plain, plain[2L, , drop = FALSE]))
  expect_identical(as.numeric(bound$BF), c(as.numeric(out$BF), as.numeric(out$BF)[2L]))
  expect_null(attr(plain$BF, "canonical_log_BF"))
  expect_null(attr(bound$BF, "canonical_log_BF"))
  expect_null(attr(bound$BF, "producer_values"))
})

test_that("declared numerical measure causes preserve detail through the RoBMA reader", {

  reason <- "Numerical scale arithmetic is unavailable for this contribution."
  sample <- with_draw_metadata(c(1, 2), measure_unavailable =
    data.frame(column = "mu", measure = "prior_density", reason = reason,
      cause = "numerical_scale_unavailable"))
  condition <- tryCatch(.brma_formula_measure_check(sample, "prior_density"), error = identity)
  expect_s3_class(condition, "BayesTools_formula_measure_unavailable")
  expect_identical(condition$reason, "numerical_scale_unavailable")
  expect_identical(condition$detail, reason)
  unknown <- with_draw_metadata(c(1, 2), measure_unavailable =
    data.frame(column = "mu", measure = "prior_density", reason = "Missing source metadata."))
  condition <- tryCatch(.brma_formula_measure_check(unknown, "prior_density"), error = identity)
  expect_s3_class(condition, "BayesTools_formula_prior_density_unavailable")
  expect_false(inherits(condition, "BayesTools_formula_measure_unavailable"))
})
