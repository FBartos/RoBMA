skip_on_cran()

.robma_mv_input_data <- function() {

  data.frame(
    yi    = c(0.05, 0.16, 0.22, 0.34),
    study = factor(rep(c("s1", "s2"), each = 2L)),
    obs   = factor(seq_len(4L)),
    x     = rep(c(0, 1), 2L)
  )
}


.robma_mv_input_V <- function() {

  kronecker(
    diag(2L),
    matrix(c(0.010, 0.004, 0.004, 0.014), nrow = 2L)
  )
}


test_that("RoBMA.mv stages the complete exact product space", {

  data <- .robma_mv_input_data()
  object <- RoBMA.mv(
    yi = yi, V = .robma_mv_input_V(), mods = ~ x,
    random = list(study = ~ 1 | study, observation = ~ 1 | obs),
    data = data, measure = "GEN", prior_unit_information_sd = 1,
    only_priors = TRUE
  )

  expect_identical(
    class(object),
    c(
      "only_priors.brma", "RoBMA.mv", "RoBMA", "brma.mv", "brma.norm",
      "brma"
    )
  )
  expect_identical(object[["selection_likelihood"]][["type"]], "exact")
  expect_identical(
    object[["selection_likelihood"]][["target"]],
    "finite_vector_product_selection"
  )
  expect_identical(.data_selection_likelihood(object[["data"]]), "exact")
  expect_true(.is_priors_weightfunction(object[["priors"]]))
  expect_true(.is_priors_PET(object[["priors"]]))
  expect_true(.is_priors_PEESE(object[["priors"]]))
  expect_length(object[["priors"]][["outcome"]][["bias"]], 9L)

  compile_modes <- vapply(
    object[["formula_design"]][["mu"]][["random_effects"]],
    `[[`,
    character(1),
    "compile_mode"
  )
  expect_true(all(compile_modes == "marginalized"))

  allocation <- object[["priors"]][["random"]][["allocation"]][[1L]]
  expect_equal(
    unname(vapply(allocation[["inclusion"]], mean, numeric(1))),
    c(0.5, 0.5)
  )
  expect_length(unique(vapply(
    object[["formula_design"]][["mu"]][["random_allocations"]][[1L]][["inclusion"]],
    `[[`,
    character(1),
    "indicator_name"
  )), 2L)
})


test_that("exact selection data retain PET and PEESE predictors", {

  data <- .robma_mv_input_data()
  object <- RoBMA.mv(
    yi = yi, V = .robma_mv_input_V(), random = ~ 1 | study,
    data = data, measure = "GEN", prior_unit_information_sd = 1,
    only_priors = TRUE
  )
  fit_data <- .create_fit_data(
    data   = object[["data"]],
    priors = object[["priors"]]
  )

  expect_equal(fit_data[["sei"]], sqrt(diag(.robma_mv_input_V())))
  expect_true(any(grepl("sel_exact_block_", names(fit_data), fixed = TRUE)))
})


test_that("RoBMA.mv bypasses selection preparation for PP ensembles", {

  data <- .robma_mv_input_data()
  object <- RoBMA.mv(
    yi = yi, V = .robma_mv_input_V(),
    data = data, measure = "GEN", model_type = "PP",
    selection_likelihood = "exact",
    known_v_parameterization = "block_mvn",
    prior_unit_information_sd = 1,
    only_priors = TRUE
  )

  expect_false(.is_priors_weightfunction(object[["priors"]]))
  expect_true(.is_priors_PET(object[["priors"]]))
  expect_true(.is_priors_PEESE(object[["priors"]]))
  expect_null(object[["selection_likelihood"]])
  expect_null(attr(
    object[["data"]],
    "selection_likelihood",
    exact = TRUE
  ))
  expect_identical(
    .data_known_v_effective_backend(object[["data"]]),
    "block_mvn"
  )
})


test_that("approximate RoBMA.mv selection uses only the latent known-V target", {

  data <- .robma_mv_input_data()
  automatic <- RoBMA.mv(
    yi = yi, V = .robma_mv_input_V(), random = ~ 1 | study,
    data = data, measure = "GEN",
    selection_likelihood = "approximate",
    prior_unit_information_sd = 1,
    only_priors = TRUE
  )
  known_V <- .data_known_v_data(automatic[["data"]])

  expect_identical(.known_v_effective_backend(known_V), "latent")
  expect_identical(.known_v_requested_parameterization(known_V), "auto")
  expect_identical(
    automatic[["selection_likelihood"]][["target"]],
    "row_selected_normal_conditional_on_random_effects"
  )

  for (backend in c("whitened", "block_mvn")) {
    expect_error(
      RoBMA.mv(
        yi = yi, V = .robma_mv_input_V(), random = ~ 1 | study,
        data = data, measure = "GEN",
        selection_likelihood = "approximate",
        known_v_parameterization = backend,
        prior_unit_information_sd = 1,
        only_priors = TRUE
      ),
      paste0(
        "'selection_likelihood = \"approximate\"' requires ",
        "'known_v_parameterization = \"latent\"'."
      ),
      fixed = TRUE
    )
  }
})


test_that("approximate PP ensembles retain non-latent known-V backends", {

  data <- .robma_mv_input_data()
  object <- RoBMA.mv(
    yi = yi, V = .robma_mv_input_V(),
    data = data, measure = "GEN", model_type = "PP",
    selection_likelihood = "approximate",
    known_v_parameterization = "whitened",
    prior_unit_information_sd = 1,
    only_priors = TRUE
  )

  expect_false(.is_priors_weightfunction(object[["priors"]]))
  expect_identical(
    .data_known_v_effective_backend(object[["data"]]),
    "whitened"
  )
})
