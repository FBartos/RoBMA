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

