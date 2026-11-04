test_that("formula prior arguments keep named NULL slots", {

  object <- brma(yi = c(.1, .2, .3), sei = rep(.1, 3),
                  mods = ~ x, scale = ~ x, data = data.frame(x = 1:3),
                  measure = "SMD", only_priors = TRUE, silent = TRUE)
  priors <- object[["priors"]]
  priors[["mods"]] <- NULL
  priors[["scale"]] <- NULL
  args <- .create_jags_formula_args(object[["data"]], priors)
  expect_identical(names(args[["formula_prior_list"]]),
                   names(args[["formula_list"]]))
  expect_true(all(vapply(args[["formula_prior_list"]], is.null, logical(1))))
})


test_that("formula compilation emits actual sampled parameter monitors", {

  object <- brma(yi = c(.1, .2, .3), sei = rep(.2, 3), mods = ~ x,
    data = data.frame(x = 1:3), measure = "SMD", only_priors = TRUE)
  args <- .create_jags_formula_args(object[["data"]], object[["priors"]])
  compiled <- BayesTools::JAGS_formula(formula = args[["formula_list"]][["mu"]],
    parameter = "mu", data = args[["formula_data_list"]][["mu"]],
    prior_list = args[["formula_prior_list"]][["mu"]],
    formula_scale = args[["formula_scale_list"]][["mu"]])
  monitors <- BayesTools::JAGS_to_monitor(compiled[["prior_list"]])
  expect_true(all(c("mu_intercept", "mu_x") %in% monitors))
  expect_match(compiled[["formula_syntax"]], "mu_x", fixed = TRUE)
})
