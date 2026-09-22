test_that("marginal means retain rejected ordinate diagnostics", {

  ordinate <- BayesTools::posterior_ordinate_attribute(
    value = 0, ordinate = 1, method = "q_grid_cmde", density_method = "qCMDE",
    diagnostics = list(estimator = "q_grid_cmde", ordinate_relative_change = .5)
  )
  object <- list(inference = list(conditional = list(mu = list(a = 1:3))))
  result <- .marginal_means_attach_iwmde_ordinate_type(
    object, "conditional", list(a = list(parameter = "mu", level = "a")),
    list(a = list(diagnostics = list(ordinate = list(status = "ok")),
                  rejected_posterior_ordinate = ordinate))
  )
  attached <- attr(result[["inference"]][["conditional"]][["mu"]][["a"]],
                   "posterior_ordinate", exact = TRUE)
  expect_identical(attached, ordinate)
  expect_match(.marginal_means_iwmde_bf_warning(attached),
               "was rejected by diagnostics", fixed = TRUE)
  expect_false(.iwmde_posterior_ordinate_supports_bf(attached))
})

test_that("single-level mixture priors color both plotted components", {

  style <- .set_dots_prior_marginal_means(NULL, 1L, 2L)
  expect_identical(style[["col"]], c("black", "black"))
  expect_identical(style[["lty"]], c(2, 2))
})
