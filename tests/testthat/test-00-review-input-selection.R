test_that("sensitivity rejects partitions that split normalization events", {

  expect_error(.selection_sensitivity_run(
    means = matrix(0, 1L, 2L), covariance = .block_covariance_zero(1L, 2L),
    context_covariance = .block_covariance_zero(1L, 2L),
    selection_sei = c(1, 1), selection_context = list(use_normal = TRUE),
    normalization_units = list(1:2), row_blocks = list(1L, 2L),
    execution_plan = list(), latent_samples = 5L
  ), "^Selection sensitivity blocks must contain complete normalization events\\.$")
})


test_that("inactive RoBMA selection branches still validate supplied controls", {

  args <- list(yi = c(.1, .2, .3), measure = "SMD", model_type = "PP",
               only_priors = TRUE, silent = TRUE)
  for (constructor in list(RoBMA, RoBMA.mv)) {
    model_args <- args
    if (identical(constructor, RoBMA.mv)) {
      model_args[["V"]] <- diag(rep(.04, 3))
    } else {
      model_args[["sei"]] <- rep(.2, 3)
    }
    expect_error(do.call(constructor, c(model_args, list(selection_control = list()))),
      "'selection_control' must be created by set_selection_likelihood_control().", fixed = TRUE)
    expect_silent(object <- do.call(constructor, c(model_args,
      list(selection_control = set_selection_likelihood_control()))))
    expect_false(.is_priors_weightfunction(object[["priors"]]))
  }
})
