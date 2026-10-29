test_that("cross-level hypotheses keep all rows and build one shared context", {

  calls    <- 0L
  attached <- NULL
  testthat::local_mocked_bindings(
    hypothesis_BF.default = function(...) data.frame(BF = c(2, 3)),
    .hypothesis_model_append_ordinate_warnings = function(table, ...) table,
    .package = "BayesTools"
  )
  testthat::local_mocked_bindings(
    .iwmde_check_point_ordinate_supported = function(...) NULL,
    .iwmde_context = function(...) { calls <<- calls + 1L; list() },
    .hypothesis_density_attach_scalar = function(posterior, request, ...) {
      attached <<- list(value = request[["values"]], spec = request[["spec"]])
      # The engine requires a valid ordinate at every requested value.
      BayesTools::posterior_metadata(posterior, "posterior_ordinate") <-
        BayesTools::posterior_ordinate_attribute(
          value          = request[["values"]],
          ordinate       = rep(1, length(request[["values"]])),
          method         = "qCMDE",
          density_method = "qCMDE"
        )
      posterior
    },
    .package = "RoBMA"
  )
  spec <- list(type = "linear", weights = c(a = 1, b = -1), prior_density = NULL)
  plan <- list(
    parameter     = "mu_group",
    label         = "group",
    point         = TRUE,
    conditional   = FALSE,
    linear_target = list(
      posterior = 1:4, hypothesis = c("contrast = 0", "contrast = 1"),
      parameter = "contrast", weights = c(a = 1, b = -1)
    ),
    targets = list(
      list(value = 0, spec = spec),
      list(value = 1, spec = spec)
    )
  )
  object <- structure(list(fit = list(TRUE)), class = "brma")
  result <- BayesTools:::.hypothesis_model_plan_execute_combination(
    plans = list(plan), hypothesis = NULL,
    provider = list(object = object, capability_object = object),
    backend = BayesTools:::.hypothesis_model_backend("qCMDE"),
    logBF = FALSE, BF01 = FALSE, seed = NULL,
    density_control = list(n_points = 50L), columns = NULL
  )
  expect_identical(result[["BF"]], c(2, 3))
  expect_identical(rownames(result), c("mu_group (1)", "mu_group (2)"))
  expect_identical(calls, 1L)
  # Both values share one estimate of the linear target with the plan's
  # fitted-coordinate weights.
  expect_identical(attached[["value"]], c(0, 1))
  expect_identical(attached[["spec"]], spec)
})

test_that("BayesTools refusals are matched by their condition class", {

  status <- function(condition, reason) {
    data.frame(
      value = 0, eligible = FALSE, condition = condition, reason = reason,
      continuous_behavior = "regular", stringsAsFactors = FALSE
    )
  }
  # The providers of a model-averaged object and of a single model; only the
  # former words a prior point mass as an inclusion Bayes factor.
  robma  <- engine_provider(gated_random_object())
  single_model <- engine_provider(single_sd_random_object(BayesTools::prior("gamma", list(2, 2))))

  # A prior point mass at the value names the inclusion Bayes factor for
  # model-averaged objects whatever the wording of the BayesTools message;
  # the refusal keeps the BayesTools class.
  density <- BayesTools::prior("normal", list(0, 1))
  target <- list(
    prior_density     = density,
    status            = status("BayesTools_point_mass_at_null", "Reworded refusal."),
    point_mass_reason = BayesTools:::.hypothesis_model_plan_point_mass_reason(robma)
  )
  refusal <- BayesTools:::.hypothesis_model_plan_target_refusal(target)
  expect_identical(
    refusal[["class"]],
    c("BayesTools_point_mass_at_null", "BayesTools_hypothesis_ordinate")
  )
  expect_match(
    refusal[["reason"]],
    paste0(
      "This parameter has a null component, so its evidence against the ",
      "null is the inclusion Bayes factor reported by 'summary()' and ",
      "'summary_models()'."
    ),
    fixed = TRUE
  )
  expect_true(refusal[["value_specific"]])
  expect_error(BayesTools:::.hypothesis_model_stop(refusal), class = "BayesTools_point_mass_at_null")

  # Other classes, and point masses of single models, keep the BayesTools
  # message.
  infinite <- BayesTools:::.hypothesis_model_plan_target_refusal(list(
    prior_density     = density,
    status            = status("BayesTools_infinite_ordinate", "Infinite ordinate."),
    point_mass_reason = BayesTools:::.hypothesis_model_plan_point_mass_reason(robma)
  ))
  expect_identical(infinite[["reason"]], "Infinite ordinate.")
  expect_identical(infinite[["class"]][[1L]], "BayesTools_infinite_ordinate")
  single <- BayesTools:::.hypothesis_model_plan_target_refusal(list(
    prior_density     = density,
    status            = status("BayesTools_point_mass_at_null", "Point mass."),
    point_mass_reason = BayesTools:::.hypothesis_model_plan_point_mass_reason(single_model)
  ))
  expect_identical(single[["reason"]], "Point mass.")

  # Missing monitored columns of a prior are recognized by class only.
  missing_columns <- structure(
    class = c(
      "BayesTools_missing_monitored_columns", "BayesTools_marglik_input",
      "error", "condition"
    ),
    list(message = "Reworded missing columns.", call = NULL)
  )
  expect_true(.iwmde_marglik_parameters_missing(missing_columns))
  expect_false(.iwmde_marglik_parameters_missing(simpleError(
    "'samples' does not contain all monitored parameters."
  )))
})

test_that("qCMDE and IWMDE are RoBMA's precomputed density methods", {

  expect_true(.density_method_uses_precomputed("qCMDE"))
  expect_true(.density_method_uses_precomputed("iwmde"))
  expect_false(.density_method_uses_precomputed("KDE"))
  expect_false(.density_method_uses_precomputed("normal", allow_normal = TRUE))
})
