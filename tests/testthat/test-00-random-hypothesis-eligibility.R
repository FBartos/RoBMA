test_that("random hypothesis eligibility uses structural rather than sampled constancy", {

  for (fixed in c(FALSE, TRUE)) {
    prior <- if (fixed) {
      BayesTools::prior("spike", list(location = .2))
    } else {
      BayesTools::prior("gamma", list(shape = 2, rate = 2))
    }
    object <- single_sd_random_object(prior)
    quantities <- hypothesis_quantities(object)
    random <- quantities[quantities[["component"]] == "random", , drop = FALSE]

    expect_gt(nrow(random), 0L)
    expect_identical(unique(random[["direction_test"]]), !fixed)
    expect_identical(unique(random[["point_test"]]), !fixed)
    expect_identical(any(grepl("fixed by the fitted model", random[["reason"]])), fixed)
  }
})


test_that("sampled random quantities take region prior odds from the exact prior", {

  # The posterior draws of the study SD are constant (0.2), yet the quantity
  # is sampled: its region hypotheses use the prior odds of its canonical
  # prior density, Gamma(2, 2), not prior draws.
  object <- single_sd_random_object(
    BayesTools::prior("gamma", list(shape = 2, rate = 2))
  )
  out <- suppressWarnings(hypothesis(
    object, "`(mu) tau(intercept)` > .5 vs `(mu) tau(intercept)` <= .5",
    density_method = "KDE", seed = 1, n_samples = 20L, columns = "all"
  ))
  expect_equal(
    as.numeric(out[["prior"]]),
    stats::pgamma(.5, 2, 2, lower.tail = FALSE) / stats::pgamma(.5, 2, 2),
    tolerance = 1e-10
  )
})


test_that("constant continuous fixtures retain the exact zero-spread KDE refusal", {

  single <- single_sd_random_object(BayesTools::prior("gamma", list(2, 2)))
  nested <- nested_allocation_random_object()
  metadata <- .brma_parameter_catalog_metadata(nested)
  quantities <- metadata[["catalog"]][["quantities"]]
  entry <- metadata[["entries"]]
  intercept <- entry[entry[["parameter"]] == "mu_intercept" &
    entry[["component"]] == "mods", , drop = FALSE]
  expect_equal(nrow(intercept), 1L)
  selector <- BayesTools::parameter_labels(quantities[
    match(intercept[["quantity_id"]], quantities[["quantity_id"]]), , drop = FALSE
  ], style = "selector")
  cases <- list(
    single = list(object = single, statement = "`(mu) tau(intercept)` = .3",
      finite_draws = 20L, value = .2),
    nested = list(object = nested, statement = paste0("`", selector, "` = 0"),
      finite_draws = 4L, value = .5)
  )
  for (name in names(cases)) {
    case <- cases[[name]]
    original <- serialize(case[["object"]], NULL, version = 2)
    error <- tryCatch(hypothesis(case[["object"]], case[["statement"]],
      density_method = "KDE", seed = 1), error = identity)
    expect_identical(class(error), c("BayesTools_hypothesis_numerical_unavailable",
      "BayesTools_hypothesis_ordinate", "BayesTools_numerical_unavailable",
      "BayesTools_numerical_condition", "error", "condition"), info = name)
    expect_identical(error[["reason"]], "zero_spread_continuous_draws", info = name)
    expect_identical(error[["operation"]], "posterior kernel density ordinate", info = name)
    expect_identical(error[["family"]], "hypothesis", info = name)
    expect_identical(error[["requested_scale"]], "log", info = name)
    expect_identical(error[["indices"]], 1L, info = name)
    expect_identical(error[["inputs"]], list(finite_draws = case[["finite_draws"]],
      value = case[["value"]]), info = name)
    expect_null(error[["diagnostics"]], info = name)
    expect_null(conditionCall(error), info = name)
    expect_identical(conditionMessage(error), paste0(
      "The Savage-Dickey Bayes factor is unavailable: zero spread continuous draws in ",
      "posterior kernel density ordinate. Use more posterior draws from the fitted model, ",
      "or check whether the parameter is fixed (identical draws have no kernel density ordinate)."
    ), info = name)
    expect_identical(serialize(case[["object"]], NULL, version = 2), original, info = name)
  }
})
