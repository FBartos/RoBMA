skip_on_cran()

.bma_mv_input_data <- function() {

  data.frame(
    yi    = c(0.10, 0.25, 0.15, 0.30),
    study = factor(c("s1", "s1", "s2", "s2")),
    esid  = factor(seq_len(4L)),
    obs   = factor(seq_len(4L)),
    group = factor(c("a", "b", "a", "b")),
    x     = c(0, 1, 0, 1),
    z     = c(-1, -1, 1, 1)
  )
}

.bma_mv_input_V <- function() {

  diag(c(0.04, 0.06, 0.05, 0.08))
}

.bma_mv_gate_priors <- function(prior_random) {

  allocations <- prior_random[["allocation"]]
  unlist(lapply(allocations, function(allocation) {
    allocation[["inclusion"]]
  }), recursive = FALSE)
}


test_that("BMA.mv stages multivariate product-space objects", {

  data <- .bma_mv_input_data()
  V    <- .bma_mv_input_V()
  object <- BMA.mv(
    yi = yi, V = V, random = ~ 1 | obs, data = data,
    measure = "GEN", prior_unit_information_sd = 1,
    only_priors = TRUE
  )

  expect_identical(
    class(object),
    c("only_priors.brma", "BMA.mv", "RoBMA", "brma.mv", "brma.norm", "brma")
  )
  expect_null(object[["priors"]][["outcome"]][["bias"]])
  expect_s3_class(object[["priors"]][["random"]], "prior_random")

  allocation <- object[["priors"]][["random"]][["allocation"]][[1L]]
  expect_length(allocation[["terms"]], 1L)
  expect_null(allocation[["weights"]])
  expect_equal(mean(allocation[["inclusion"]][[1L]]), 0.5)
})


test_that("BMA.mv uses independent default gates without renormalizing allocations", {

  data <- .bma_mv_input_data()
  object <- BMA.mv(
    yi = yi, V = .bma_mv_input_V(),
    random = list(study = ~ 1 | study, observation = ~ 1 | obs),
    data = data, measure = "GEN", prior_unit_information_sd = 1,
    only_priors = TRUE
  )
  allocation <- object[["priors"]][["random"]][["allocation"]][[1L]]
  gates      <- .bma_mv_gate_priors(object[["priors"]][["random"]])

  expect_equal(unname(allocation[["terms"]]), c("study", "observation"))
  expect_s3_class(allocation[["weights"]], "prior.simplex")
  expect_equal(allocation[["weights"]][["parameters"]][["alpha"]], c(1, 1))
  expect_equal(unname(vapply(gates, mean, numeric(1))), c(0.5, 0.5))

  formula_design <- .compile_prior_random_formula_design(
    data           = object[["data"]],
    prior_location = object[["priors"]][["location"]],
    prior_random   = object[["priors"]][["random"]]
  )
  allocation_compiled <- formula_design[["random_allocations"]][[1L]]
  indicators <- vapply(
    allocation_compiled[["inclusion"]],
    `[[`,
    character(1),
    "indicator_name"
  )
  expect_length(unique(indicators), 2L)
  expect_setequal(names(indicators), c("study", "observation"))
})


test_that("BMA.mv decomposes heterogeneity mixture weights into gate and slab", {

  data <- .bma_mv_input_data()
  alternatives <- list(
    BayesTools::prior(
      "normal", parameters = list(mean = 0, sd = 0.25),
      truncation = list(0, Inf), prior_weights = 3
    ),
    BayesTools::prior(
      "normal", parameters = list(mean = 0, sd = 0.75),
      truncation = list(0, Inf), prior_weights = 2
    )
  )
  null <- BayesTools::prior(
    "spike", parameters = list(location = 0), prior_weights = 5
  )
  object <- BMA.mv(
    yi = yi, V = .bma_mv_input_V(), random = ~ 1 | obs,
    prior_heterogeneity = alternatives,
    prior_heterogeneity_null = null,
    data = data, measure = "GEN", prior_unit_information_sd = 1,
    only_priors = TRUE
  )
  allocation <- object[["priors"]][["random"]][["allocation"]][[1L]]
  slab       <- allocation[["sd"]]
  gate       <- allocation[["inclusion"]][[1L]]

  expect_s3_class(slab, "prior.mixture")
  expect_equal(attr(slab, "prior_weights"), c(3, 2))
  expect_equal(mean(gate), 0.5)
})


test_that("BMA.mv supports fixed inclusion and exclusion gates", {

  data <- .bma_mv_input_data()
  included <- BMA.mv(
    yi = yi, V = .bma_mv_input_V(), random = ~ 1 | obs,
    prior_heterogeneity_null = FALSE,
    data = data, measure = "GEN", prior_unit_information_sd = 1,
    only_priors = TRUE
  )
  excluded <- BMA.mv(
    yi = yi, V = .bma_mv_input_V(), random = ~ 1 | obs,
    prior_heterogeneity = FALSE,
    data = data, measure = "GEN", prior_unit_information_sd = 1,
    only_priors = TRUE
  )

  expect_equal(mean(.bma_mv_gate_priors(included[["priors"]][["random"]])[[1L]]), 1)
  expect_equal(mean(.bma_mv_gate_priors(excluded[["priors"]][["random"]])[[1L]]), 0)

  included_null <- BMA.mv(
    yi = yi, V = .bma_mv_input_V(), random = ~ 1 | obs,
    prior_heterogeneity_null = NULL,
    data = data, measure = "GEN", prior_unit_information_sd = 1,
    only_priors = TRUE
  )
  excluded_null <- BMA.mv(
    yi = yi, V = .bma_mv_input_V(), random = ~ 1 | obs,
    prior_heterogeneity = NULL,
    data = data, measure = "GEN", prior_unit_information_sd = 1,
    only_priors = TRUE
  )

  expect_equal(mean(.bma_mv_gate_priors(included_null[["priors"]][["random"]])[[1L]]), 1)
  expect_equal(mean(.bma_mv_gate_priors(excluded_null[["priors"]][["random"]])[[1L]]), 0)
})


test_that("BMA.mv gates nested structures before their within-component split", {

  data <- .bma_mv_input_data()
  object <- BMA.mv(
    yi = yi, V = .bma_mv_input_V(), random = ~ 1 | study / esid,
    data = data, measure = "GEN", prior_unit_information_sd = 1,
    only_priors = TRUE
  )
  allocations <- object[["priors"]][["random"]][["allocation"]]

  expect_length(allocations, 2L)
  expect_null(allocations[[1L]][["weights"]])
  expect_length(allocations[[1L]][["inclusion"]], 1L)
  expect_s3_class(allocations[[2L]][["weights"]], "prior.simplex")
  expect_length(allocations[[2L]][["terms"]], 2L)
  expect_equal(
    allocations[[2L]][["parent"]][["allocation"]],
    allocations[[1L]][["name"]]
  )

  inclusion_map <- .random_component_inclusion_map(object)
  gate <- unique(unlist(inclusion_map, use.names = FALSE))
  expect_length(gate, 1L)
  expect_identical(inclusion_map[["study"]], gate)
  expect_identical(inclusion_map[["esid_study"]], gate)

  conditioning <- .random_component_conditioning_parameters(
    object     = object,
    components = c("study", "esid_study", "total")
  )
  expect_identical(conditioning[["study"]], gate)
  expect_identical(conditioning[["esid_study"]], gate)
  expect_identical(conditioning[["total"]], gate)
  expect_error(
    .random_component_conditioning_parameters(
      object     = object,
      components = "missing"
    ),
    "No inclusion gate is available for random-effect component 'missing'.",
    fixed = TRUE
  )
})


test_that("BMA.mv gates a random-coefficient block before its SD split", {

  data <- .bma_mv_input_data()
  object <- BMA.mv(
    yi = yi, V = .bma_mv_input_V(),
    random = ~ us(0 + group | study),
    prior_heterogeneity = BayesTools::prior_random(
      study = BayesTools::random_block(
        contrasts = c(group = "independent")
      )
    ),
    data = data, measure = "GEN", prior_unit_information_sd = 1,
    only_priors = TRUE
  )
  allocations <- object[["priors"]][["random"]][["allocation"]]

  expect_length(allocations, 2L)
  expect_length(allocations[[1L]][["inclusion"]], 1L)
  expect_null(allocations[[1L]][["weights"]])
  expect_identical(allocations[[2L]][["target"]], "sd_component")
  expect_identical(allocations[[2L]][["scale"]], "mean_variance")
  expect_equal(
    allocations[[2L]][["weights"]][["parameters"]][["alpha"]],
    c(1, 1)
  )
})


test_that("BMA.mv gates component-specific scale regressions independently", {

  data <- .bma_mv_input_data()
  object <- BMA.mv(
    yi = yi, V = .bma_mv_input_V(),
    random = list(study = ~ 1 | study, observation = ~ 1 | obs),
    scale = list(study = ~ x, observation = ~ 1),
    data = data, measure = "GEN", prior_unit_information_sd = 1,
    only_priors = TRUE
  )
  allocations <- object[["priors"]][["random"]][["allocation"]]

  expect_identical(
    .data_scale_formula_sources(object[["data"]]),
    c(study = "tau_study", observation = "tau_observation")
  )
  expect_length(allocations, 2L)
  expect_true(all(vapply(
    allocations,
    function(allocation) length(allocation[["inclusion"]]) == 1L,
    logical(1)
  )))
  expect_true(all(vapply(
    allocations,
    function(allocation) is.null(allocation[["weights"]]),
    logical(1)
  )))
  expect_named(
    .random_component_inclusion_map(object),
    c("study", "observation")
  )
})


test_that("BMA.mv preserves structural prior_random overrides", {

  data <- .bma_mv_input_data()
  override <- BayesTools::prior_random(
    study = BayesTools::random_block(parameterization = "noncentered")
  )
  object <- BMA.mv(
    yi = yi, V = .bma_mv_input_V(), random = ~ us(1 + x | study),
    prior_heterogeneity = override,
    data = data, measure = "GEN", prior_unit_information_sd = 1,
    only_priors = TRUE
  )

  expect_identical(
    object[["priors"]][["random"]][["blocks"]][["study"]][["parameterization"]],
    "noncentered"
  )
  expect_length(.bma_mv_gate_priors(object[["priors"]][["random"]]), 1L)

  expect_error(
    BMA.mv(
      yi = yi, V = .bma_mv_input_V(), random = ~ 1 | obs,
      prior_heterogeneity = BayesTools::prior_random(
        sd = BayesTools::prior(
          "normal", list(mean = 0, sd = 1), truncation = list(0, Inf)
        )
      ),
      data = data, measure = "GEN", prior_unit_information_sd = 1,
      only_priors = TRUE
    ),
    "cannot add component gates to a user-supplied random-effect scale architecture",
    fixed = TRUE
  )
})


test_that("BMA.mv singular known V requires structurally active regularization", {

  data <- data.frame(
    yi  = c(0.10, 0.25),
    obs = seq_len(2L)
  )
  singular_V <- matrix(0.04, nrow = 2L, ncol = 2L)

  suppressWarnings(expect_error(
    BMA.mv(
      yi = yi, V = singular_V, random = ~ 1 | obs,
      known_v_parameterization = "block_mvn",
      data = data, measure = "GEN", prior_unit_information_sd = 1,
      only_priors = TRUE
    ),
    "singular",
    ignore.case = TRUE
  ))
  suppressWarnings(expect_no_error(BMA.mv(
    yi = yi, V = singular_V, random = ~ 1 | obs,
    known_v_parameterization = "block_mvn",
    prior_heterogeneity_null = FALSE,
    data = data, measure = "GEN", prior_unit_information_sd = 1,
    only_priors = TRUE
  )))
})


test_that("BMA.mv rejects publication-bias arguments", {

  data <- .bma_mv_input_data()
  expect_error(
    BMA.mv(
      yi = yi, V = .bma_mv_input_V(), random = ~ 1 | obs,
      prior_bias = BayesTools::prior_none(),
      data = data, measure = "GEN", prior_unit_information_sd = 1,
      only_priors = TRUE
    ),
    "Selection/publication-bias priors are not supported in BMA.mv().",
    fixed = TRUE
  )
})
