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


.robma_mv_input_bias_priors <- function(other_random_effects = "condition",
                                        known_sampling_variance = "condition") {

  .default_prior.bias_alt(
    model_type                = "PSMA",
    measure                   = "GEN",
    data                      = NULL,
    prior_unit_information_sd = 1,
    weightfunction_model = BayesTools::selection_model(
      other_random_effects    = other_random_effects,
      known_sampling_variance = known_sampling_variance,
      group                   = "study"
    )
  )
}


test_that("omitted and NULL mv random terms give fixed effects without a heterogeneity mixture", {

  args <- list(yi = c(.1, .2), vi = c(.04, .09), data = data.frame(x = c(0, 1)),
    measure = "GEN", prior_unit_information_sd = 1, only_priors = TRUE, silent = TRUE)
  no_random <- NULL
  for (constructor in list(brma.mv, BMA.mv, bselmodel.mv, RoBMA.mv, bPET.mv, bPEESE.mv)) {
    object <- do.call(constructor, args)
    literal <- do.call(constructor, c(args, list(random = NULL)))
    evaluated <- do.call(constructor, c(args, list(random = quote(no_random))))
    expect_identical(object[["data"]], literal[["data"]])
    expect_identical(object[["data"]], evaluated[["data"]])
    expect_identical(object[["priors"]], literal[["priors"]])
    expect_identical(object[["priors"]], evaluated[["priors"]])
    expect_identical(.fixed_tau_prior_value(object[["priors"]]), 0)
    expect_false(BayesTools::is.prior.mixture(object[["priors"]][["outcome"]][["tau"]]))
    expect_null(object[["priors"]][["random"]])
    expect_false(.is_data_random(object[["data"]]))
    expect_false("theta" %in% names(.create_fit_priors(object[["data"]], object[["priors"]])))
    if (.is_data_joint_selection(object[["data"]])) {
      model <- .data_selection_model(object[["data"]])
      expect_length(model[["sources"]][["random"]], 0L)
      expect_false(model[["applicability"]][["estimate_random_effects"]])
      expect_false(model[["applicability"]][["other_random_effects"]])
    }
    for (prior in list(NULL, FALSE, BayesTools::prior("spike", parameters = list(0)),
      BayesTools::prior("normal", parameters = list(0, 1), truncation = list(0, Inf)))) {
      expect_error(do.call(constructor, c(args, list(prior_heterogeneity = prior))),
        "The 'prior_heterogeneity' argument can be used only when 'random' is specified.",
        fixed = TRUE)
    }
    expect_error(do.call(constructor, c(args, list(scale = ~ x))),
      "The 'scale' argument requires 'random' in multivariate models.", fixed = TRUE)
    if (identical(constructor, BMA.mv) || identical(constructor, RoBMA.mv)) {
      expect_error(do.call(constructor, c(args, list(prior_heterogeneity_null = NULL))),
        "The 'prior_heterogeneity_null' argument can be used only when 'random' is specified.",
        fixed = TRUE)
    }
  }
  ordinary <- do.call(brma, args)
  expect_false(BayesTools::is.prior.point(ordinary[["priors"]][["outcome"]][["tau"]]))
  fixed <- do.call(brma, c(args, list(prior_heterogeneity = NULL)))
  expect_identical(.fixed_tau_prior_value(fixed[["priors"]]), 0)
})


test_that("constructor defaults agree and explicit child targets remain authoritative", {

  fields <- c("estimate_random_effects", "other_random_effects", "known_sampling_variance")
  for (constructor in list(bselmodel, RoBMA, bselmodel.mv, RoBMA.mv)) {
    expect_error(constructor(selection = list()),
      "'selection' must be a specification from 'selection_model()'.", fixed = TRUE)
    multivariate <- identical(constructor, bselmodel.mv) || identical(constructor, RoBMA.mv)
    args <- list(yi = c(.1, .2, .3, .4), measure = "GEN",
      prior_unit_information_sd = 1, only_priors = TRUE, silent = TRUE)
    if (multivariate) {
      args$V <- diag(rep(.01, 4L))
      args$data <- .robma_mv_input_data()
      args$random <- ~ 1 | study
    } else {
      args$sei <- rep(.1, 4L)
      args$cluster <- c("a", "a", "b", "b")
    }
    default <- BayesTools::selection_model()
    expect_identical(default[fields], list(
      estimate_random_effects = "integrate", other_random_effects = "condition",
      known_sampling_variance = "integrate"))
    if (multivariate) default <- BayesTools::selection_model(group = "study")
    explicit <- do.call(constructor, c(args, list(selection = default)))
    generated <- .data_selection_model(explicit[["data"]])
    expect_identical(generated[fields], default[fields])
    expect_true(all(vapply(generated[["branches"]][generated[["active_branches"]]],
      identical, logical(1), y = default)))
    if (!multivariate) {
      implicit <- do.call(constructor, args)
      expect_identical(implicit[["data"]], explicit[["data"]])
      expect_identical(implicit[["priors"]], explicit[["priors"]])
      expect_identical(
        .create_model_syntax(implicit[["data"]], implicit[["priors"]]),
        .create_model_syntax(explicit[["data"]], explicit[["priors"]]))
    }

    requested <- BayesTools::selection_model(estimate_random_effects = "condition",
      other_random_effects = "integrate", known_sampling_variance = "integrate",
      group = "study")
    if (!multivariate) requested <- BayesTools::selection_model(
      estimate_random_effects = "condition", other_random_effects = "integrate",
      known_sampling_variance = "integrate")
    args[["prior_bias"]] <- BayesTools::prior_weightfunction(
      "one-sided", steps = .025, weights = BayesTools::wf_fixed(c(1, .5)), model = requested)
    supplied <- do.call(constructor, c(args, list(selection = default)))
    model <- .data_selection_model(supplied[["data"]])
    expect_identical(model[fields], requested[fields])
    expect_true(all(vapply(model[["branches"]][model[["active_branches"]]],
      identical, logical(1), y = requested)))
  }
})


test_that("all conditioning cells and weight rules retain their declared source roles", {

  data <- .robma_mv_input_data()
  V <- known_v_factor(c(.006, .010, .006, .010), cbind(
    c(sqrt(.004), sqrt(.004), 0, 0),
    c(0, 0, sqrt(.004), sqrt(.004))
  ))
  cells <- expand.grid(estimate_random_effects = c("condition", "integrate"),
    other_random_effects = c("condition", "integrate"),
    known_sampling_variance = c("condition", "integrate"), stringsAsFactors = FALSE)
  for (i in seq_len(nrow(cells))) {
    cell <- as.list(cells[i, , drop = FALSE])
    for (weight_rule in c("product", "best")) {
      requested <- do.call(BayesTools::selection_model,
        c(cell, list(weight_rule = weight_rule, group = "study")))
      prior <- BayesTools::prior_weightfunction(
        "one-sided", steps = .025, weights = BayesTools::wf_fixed(c(1, .5)),
        model = requested
      )
      object <- bselmodel.mv(
        yi = yi, V = V, data = data,
        random = list(study = ~ 1 | study, observation = ~ 1 | obs),
        prior_bias = prior, measure = "GEN", prior_unit_information_sd = 1,
        only_priors = TRUE, silent = TRUE
      )
      model <- .data_selection_model(object[["data"]])
      expect_identical(model[["branches"]][[1L]], requested)
      expect_identical(model[["groups"]][["group_index"]], c(1L, 1L, 2L, 2L))
      expect_identical(.selection_retains_any_random(object[["data"]]),
        cell$estimate_random_effects == "condition" || cell$other_random_effects == "condition")
      expect_identical(.selection_retains_sampling(object[["data"]]),
        cell$known_sampling_variance == "condition")
      expect_identical(vapply(model[["sources"]][["random"]], `[[`,
        character(1), "role"), c("other", "estimate"))
      expect_identical(vapply(model[["sources"]][["random"]], `[[`,
        logical(1), "retained"),
        c(cell$other_random_effects == "condition", cell$estimate_random_effects == "condition"))
      stored_prior <- .selection_bias_priors(object[["priors"]])[[1L]]
      expect_identical(stored_prior, prior)
    }
  }
})


test_that("independent random slopes remain integrated estimate-level sources", {

  data <- .robma_mv_input_data()
  data[["x"]] <- c(-1, .5, 2, .8)
  object <- bselmodel.mv(
    yi = yi, vi = rep(.02, 4L), random = ~ diag(0 + x | obs), data = data,
    prior_bias = BayesTools::prior_weightfunction(
      "one-sided", steps = .025, weights = BayesTools::wf_fixed(c(1, .5)),
      model = BayesTools::selection_model(group = "study")
    ),
    measure = "GEN", prior_unit_information_sd = 1,
    only_priors = TRUE, silent = TRUE
  )
  model <- .data_selection_model(object[["data"]])
  expect_identical(model[["sources"]][["random"]],
                   list(list(name = "obs", role = "estimate", retained = FALSE)))
  expect_true(model[["applicability"]][["estimate_random_effects"]])
  expect_false(model[["applicability"]][["other_random_effects"]])
  expect_identical(.data_selection_execution_plan(object[["data"]])[["row_blocks"]],
                   as.list(seq_len(4L)))
  expect_identical(model[["groups"]][["group_index"]], c(1L, 1L, 2L, 2L))
})


test_that("active selection branches require a common cell and publication partition", {

  data <- .robma_mv_input_data()
  data[["publication"]] <- c("p1", "p1", "p2", "p2")
  data[["crossed"]] <- c("p1", "p2", "p1", "p2")
  prior <- function(model) BayesTools::prior_weightfunction(
    "one-sided", steps = .025, weights = BayesTools::wf_fixed(c(1, .5)),
    model = model
  )
  stage <- function(models) RoBMA.mv(
    yi = yi, V = .robma_mv_input_V(), random = ~ 1 | study,
    data = data, prior_bias = lapply(models, prior), measure = "GEN",
    prior_unit_information_sd = 1, only_priors = TRUE, silent = TRUE
  )
  common <- BayesTools::selection_model(group = "study")
  other_cell <- BayesTools::selection_model(
    other_random_effects = "integrate", group = "study"
  )
  expect_error(stage(list(common, other_cell)), paste0(
    "Active weightfunction branches must use the same 'estimate_random_effects', ",
    "'other_random_effects', and 'known_sampling_variance' settings."
  ), fixed = TRUE)
  expect_error(stage(list(common, BayesTools::selection_model(group = "crossed"))),
               "Active weightfunction branches must use the same publication partition.",
               fixed = TRUE)
  equivalent <- BayesTools::selection_model(group = "publication")
  object <- stage(list(common, equivalent))
  model <- .data_selection_model(object[["data"]])
  expect_identical(model[["branches"]][model[["active_branches"]]],
                   list(common, equivalent))
  expect_identical(model[["branch_groups"]][[1L]][["group_index"]],
                   model[["branch_groups"]][[2L]][["group_index"]])
})


test_that("RoBMA.mv stages the complete marginal product space", {

  data <- .robma_mv_input_data()
  object <- RoBMA.mv(
    yi = yi, V = .robma_mv_input_V(), mods = ~ x,
    random = list(study = ~ 1 | study, observation = ~ 1 | obs),
    data = data, measure = "GEN", prior_unit_information_sd = 1,
    prior_bias = .robma_mv_input_bias_priors("integrate", "integrate"),
    only_priors = TRUE
  )

  expect_identical(
    class(object),
    c(
      "only_priors.brma", "RoBMA.mv", "RoBMA", "brma.mv", "brma.norm",
      "brma"
    )
  )
  expect_identical(
    .data_selection_model(object[["data"]])[c("estimate_random_effects", "other_random_effects", "known_sampling_variance")],
    list(estimate_random_effects = "integrate", other_random_effects = "integrate", known_sampling_variance = "integrate")
  )
  expect_identical(
    .data_selection_execution_plan(object[["data"]])[["statistical_target"]],
    "conditional_gaussian_vector_selection"
  )
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


test_that("joint selection data retain PET and PEESE predictors", {

  data <- .robma_mv_input_data()
  object <- RoBMA.mv(
    yi = yi, V = .robma_mv_input_V(), random = ~ 1 | study,
    data = data, measure = "GEN", prior_unit_information_sd = 1,
    prior_bias = .robma_mv_input_bias_priors("integrate", "integrate"),
    only_priors = TRUE
  )
  fit_data <- .create_fit_data(
    data   = object[["data"]],
    priors = object[["priors"]]
  )

  expect_equal(fit_data[["sei"]], sqrt(diag(.robma_mv_input_V())))
  expect_true(any(grepl("sel_joint_block_", names(fit_data), fixed = TRUE)))
})


test_that("RoBMA.mv bypasses selection preparation for PP ensembles", {

  data <- .robma_mv_input_data()
  object <- RoBMA.mv(
    yi = yi, V = .robma_mv_input_V(),
    data = data, measure = "GEN", model_type = "PP",
    selection = BayesTools::selection_model(
      other_random_effects = "integrate", known_sampling_variance = "integrate"),
    known_v_parameterization = "block_mvn",
    prior_unit_information_sd = 1,
    only_priors = TRUE
  )

  expect_false(.is_priors_weightfunction(object[["priors"]]))
  expect_true(.is_priors_PET(object[["priors"]]))
  expect_true(.is_priors_PEESE(object[["priors"]]))
  expect_null(.data_selection_model(object[["data"]]))
  expect_null(.data_selection_execution_plan(object[["data"]]))
  expect_identical(
    .data_known_v_effective_backend(object[["data"]]),
    "block_mvn"
  )
})


test_that("conditional RoBMA.mv sampling structure is independent of the Gaussian backend", {

  data <- .robma_mv_input_data()
  automatic <- RoBMA.mv(
    yi = yi, V = .robma_mv_input_V(), random = ~ 1 | study,
    data = data, measure = "GEN",
    prior_bias = .robma_mv_input_bias_priors(),
    prior_unit_information_sd = 1,
    only_priors = TRUE
  )
  known_V <- .data_known_v_data(automatic[["data"]])

  expect_identical(.known_v_requested_parameterization(known_V), "auto")
  expect_identical(
    .data_selection_model(automatic[["data"]])[c("estimate_random_effects", "other_random_effects", "known_sampling_variance")],
    list(estimate_random_effects = "integrate", other_random_effects = "condition", known_sampling_variance = "condition")
  )
  sampling <- .selection_sampling_structure(automatic[["data"]])

  for (backend in c("whitened", "block_mvn")) {
    object <- RoBMA.mv(
        yi = yi, V = .robma_mv_input_V(), random = ~ 1 | study,
        data = data, measure = "GEN",
        prior_bias = .robma_mv_input_bias_priors(),
        known_v_parameterization = backend,
        prior_unit_information_sd = 1,
        only_priors = TRUE
      )
    expect_identical(
      .data_selection_model(object[["data"]]),
      .data_selection_model(automatic[["data"]])
    )
    expect_equal(
      .selection_sampling_structure(object[["data"]]),
      sampling,
      tolerance = 0
    )
  }
})


test_that("PP ensembles retain non-latent known-V backends", {

  data <- .robma_mv_input_data()
  object <- RoBMA.mv(
    yi = yi, V = .robma_mv_input_V(),
    data = data, measure = "GEN", model_type = "PP",
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
