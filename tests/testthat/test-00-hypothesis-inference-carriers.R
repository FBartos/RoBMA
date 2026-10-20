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
  expect_identical(attr(out$BF_error, "name"), "error%(BF)")
  expect_true(any(grepl("error%(BF)", capture.output(print(out)), fixed = TRUE)))
  expect_identical(as.numeric(out$BF_error),
    c(as.numeric(first$BF_error)[1L], as.numeric(second$BF_error), as.numeric(first$BF_error)[2L]))
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
test_that("variance point-region diagnostics align actual repeated point producers", {

  object <- gated_random_object()
  tau2 <- "`(mu) study: tau2(intercept)`"
  points <- c(.09, .16, .09)
  run <- function(statement) hypothesis(object, statement, density_method = "KDE", columns = "all", seed = 1, logBF = TRUE)
  direct <- run(sprintf("%s = %.17g vs %s != %.17g", tau2, points, tau2, points))
  expected <- attr(direct, "numerical_diagnostics")
  for(reverse in c(FALSE, TRUE)){
    statements <- if(reverse) sprintf("%s > %.17g vs %s = %.17g", tau2, points, tau2, points) else
      sprintf("%s = %.17g vs %s > %.17g", tau2, points, tau2, points)
    mixed <- run(statements)
    expect_identical(attr(mixed, "numerical_diagnostics"), expected)
    expect_identical(attr(mixed$BF, "numerical_diagnostics"), expected)
    expect_length(attr(mixed, "prior_numerical_diagnostics"), length(points))
  }
  for(i in seq_along(points)){
    diagnostic <- expected[[i]]
    expect_equal(diagnostic$log_prior_height, log(direct$prior[i]), tolerance = 1e-12)
    expect_equal(diagnostic$log_posterior_height, log(direct$posterior[i]), tolerance = 1e-12)
    expect_identical(diagnostic$requested_value, points[i])
    expect_equal(diagnostic$evaluation_value, sqrt(points[i]), tolerance = 1e-14)
  }
})


test_that("equivalent SD and variance zero-point zero-region BFs are canonical undefined", {

  object <- gated_random_object()
  sd <- "`(mu) study: tau(intercept)`"
  variance <- "`(mu) study: tau2(intercept)`"
  for(reverse in c(FALSE, TRUE)){
    statement <- function(name, value){
      if(reverse) sprintf("%s > %.17g vs %s = %.17g", name, value, name, value) else
        sprintf("%s = %.17g vs %s > %.17g", name, value, name, value)
    }
    sd_result <- hypothesis(object, statement(sd, 10), density_method = "KDE", columns = "all", seed = 1)
    var_result <- hypothesis(object, statement(variance, 100), density_method = "KDE", columns = "all", seed = 1)
    expect_identical(attr(var_result, "raw_log_BF"), attr(sd_result, "raw_log_BF"))
    expect_identical(attr(var_result, "raw_log_BF"), NA_real_)
    expect_false(any(is.nan(attr(var_result, "raw_log_BF"))))
    expect_true(any(grepl("both numerically zero; transitive BF is undefined",
      attr(var_result, "warnings"), fixed = TRUE)))
    expect_identical(as.numeric(update(var_result, logBF = TRUE)$BF), NA_real_)
    expect_identical(as.numeric(update(var_result, BF01 = TRUE)$BF), NA_real_)
    expect_identical(attr(var_result[1L, , drop = FALSE], "raw_log_BF"), NA_real_)
  }
})

test_that("single and gated variance composition retains indexed diagnostics and valid boundaries", {

  objects <- list(gated = gated_random_object(),
    single = single_sd_random_object(BayesTools::prior("gamma", list(2, 2)), varying_draws = TRUE))
  undefined_message <- "The posterior point ordinate and region mass are both numerically zero; transitive BF is undefined."
  boundary_message <- "Posterior region mass is zero; reported BF is boundary-valued."
  cuts <- c(10, .3, 3, 10)
  for(name in names(objects)){
    object <- objects[[name]]
    original <- serialize(object, NULL, version = 2)
    sd <- if(name == "gated") "`(mu) study: tau(intercept)`" else "`(mu) tau(intercept)`"
    variance <- if(name == "gated") "`(mu) study: tau2(intercept)`" else "`(mu) tau2(intercept)`"
    run <- function(statements, ...) hypothesis(object, statements,
      density_method = "KDE", columns = "all", seed = 1, ...)
    real_point <- run(sprintf("%s = %.17g vs %s != %.17g", variance, cuts^2, variance, cuts^2))
    point_diagnostics <- attr(real_point, "numerical_diagnostics")
    for(reverse in c(FALSE, TRUE)){
      statement <- function(parameter, values){
        if(reverse) sprintf("%s > %.17g vs %s = %.17g", parameter, values, parameter, values) else
          sprintf("%s = %.17g vs %s > %.17g", parameter, values, parameter, values)
      }
      sd_result <- run(statement(sd, cuts))
      var_result <- run(statement(variance, cuts^2))
      logs <- attr(var_result, "raw_log_BF")
      expect_identical(logs[c(1L, 4L)], rep(NA_real_, 2L))
      expect_equal(logs, attr(sd_result, "raw_log_BF"), tolerance = 1e-12)
      expect_true(is.infinite(logs[[3L]]))
      expect_identical(sign(logs[[3L]]), if(reverse) -1 else 1)
      expect_identical(attr(var_result, "numerical_diagnostics"), point_diagnostics)
      expect_identical(attr(var_result$BF, "numerical_diagnostics"), point_diagnostics)
      prior_diagnostics <- attr(var_result, "prior_numerical_diagnostics")
      expect_length(prior_diagnostics, length(cuts))
      expect_true(all(vapply(prior_diagnostics, Negate(is.null), logical(1))))
      warnings <- attr(var_result, "warnings")
      for(row in rownames(var_result)[c(1L, 4L)]){
        selected <- warnings[names(warnings) == row]
        expect_true(any(grepl(undefined_message, selected, fixed = TRUE)))
        expect_false(any(grepl(boundary_message, selected, fixed = TRUE)))
      }
      expect_true(any(grepl(boundary_message,
        warnings[names(warnings) == rownames(var_result)[[3L]]], fixed = TRUE)))
      for(logBF in c(FALSE, TRUE)) for(BF01 in c(FALSE, TRUE)){
        viewed <- update(var_result, logBF = logBF, BF01 = BF01)
        expect_identical(attr(viewed, "raw_log_BF"), logs)
        expect_identical(attr(viewed$BF, "canonical_log_BF"), logs)
        expect_identical(as.numeric(viewed$BF)[c(1L, 4L)], rep(NA_real_, 2L))
        expect_identical(attr(update(viewed, logBF = FALSE, BF01 = FALSE), "raw_log_BF"), logs)
      }
      order <- c(4L, 2L, 1L, 3L, 2L)
      subset <- var_result[order, c("BF", "BF_error"), drop = FALSE]
      expect_identical(attr(subset, "raw_log_BF"), logs[order])
      expect_identical(attr(subset$BF, "numerical_diagnostics"), point_diagnostics[order])
      without_bf <- var_result[order, "BF_error", drop = FALSE]
      expect_null(attr(without_bf, "raw_BF"))
      expect_null(attr(without_bf, "raw_log_BF"))
      expect_identical(attr(without_bf, "numerical_diagnostics"), point_diagnostics[order])
      expect_identical(attr(without_bf, "prior_numerical_diagnostics"), prior_diagnostics[order])
      for(i in c(1L, 4L)){
        expect_identical(point_diagnostics[[i]]$log_posterior_height, -Inf)
        expect_true(is.finite(point_diagnostics[[i]]$log_prior_height))
        expect_identical(point_diagnostics[[i]]$evaluation_value, cuts[[i]])
        expect_identical(point_diagnostics[[i]]$requested_value, cuts[[i]]^2)
      }
    }
    expect_identical(serialize(object, NULL, version = 2), original)
  }
})
