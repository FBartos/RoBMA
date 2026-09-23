test_that("missing random scale-block metadata has a targeted error", {

  object <- brma.mv(yi = c(.1, .2, .3), vi = rep(.01, 3),
                     random = ~ 1 | id, scale = ~ x,
                     data = data.frame(id = 1:3, x = 1:3),
                     measure = "SMD", only_priors = TRUE, silent = TRUE)
  design <- .fitted_formula_design(object, "mu")
  design[["random_effects"]][[1L]][["block_name"]] <- "absent"
  expect_error(.validate_prior_random_scale_sources(object[["data"]], design),
               "^Random-effect scale metadata contain an unknown block\\.$")
})
