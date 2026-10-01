test_that("post-fit selection controls are explicit and validated", {

  control <- set_selection_likelihood_control(max_points_per_scramble = 32768)
  for (method in c("qCMDE", "IWMDE")) {
    defaults <- .density_control_normalize(method)
    explicit_null <- .density_control_normalize(
      method, list(integration_control = NULL)
    )
    requested <- .density_control_normalize(
      method, list(integration_control = control)
    )
    expect_null(defaults[["integration_control"]])
    expect_identical(
      .iwmde_density_control_provenance(defaults),
      .iwmde_density_control_provenance(explicit_null)
    )
    expect_identical(requested[["integration_control"]], control)
    expect_identical(
      .iwmde_density_control_provenance(requested)[["integration_control"]],
      unclass(control)
    )
    expect_error(
      .density_control_normalize(method, list(integration_control = list())),
      "'density_control$integration_control' must be created by set_selection_likelihood_control().",
      fixed = TRUE
    )
  }
  expect_error(
    .density_control_normalize("KDE", list(integration_control = control)),
    "'density_control' is only used when 'density_method' is 'qCMDE' or 'IWMDE'.",
    fixed = TRUE
  )
  invalid <- control
  invalid[["max_points_per_scramble"]] <- 0L
  expect_error(
    .density_control_normalize("qCMDE", list(integration_control = invalid)),
    "The 'max_points_per_scramble' must be equal or higher than 512.",
    fixed = TRUE
  )
  expect_error(
    .iwmde_context_with_integration_control(list(data = list()), control),
    "'density_control$integration_control' is unavailable for this model: post-fit integration controls require a bound Gaussian selection model.",
    fixed = TRUE
  )
})


test_that("post-fit budgets preserve fitted metadata and isolate evaluation caches", {

  control <- set_selection_likelihood_control(
    points_per_scramble = 16L, max_points_per_scramble = 64L, scrambles = 2L
  )
  requested <- control
  requested[["max_points_per_scramble"]] <- 256L
  sampling <- structure(list(
    representation = "diagonal_factor",
    diagonal       = rep(.2, 8L),
    loading        = cbind(c(.1, .2, .1, rep(0, 5L)))
  ), class = c("RoBMA_selection_sampling_plan", "list"))
  row_blocks <- list(1:3, 4:5, 6:7, 8L)
  factor_blocks <- .selection_joint_sampling_factor_blocks(sampling, row_blocks)
  plan <- .selection_joint_execution_plan(
    row_blocks             = row_blocks,
    block_methods          = c("factor", "rank_one", "dense", "singleton"),
    factor_ranks           = c(3L, 1L, NA_integer_, 0L),
    selection_control      = control,
    sampling               = sampling,
    sampling_factor_blocks = factor_blocks,
    random_covariance      = NULL
  )
  data <- structure(
    list(measure = "SMD"), known_V = TRUE,
    selection_model = structure(list(
      schema_version = 3L,
      estimate_random_effects = "integrate", other_random_effects = "integrate", known_sampling_variance = "integrate"
    ), class = c("RoBMA_selection_model", "list")),
    selection_execution_plan = plan
  )
  object <- list(data = data, fit = list(), likelihood = list(family = "normal"))
  context <- .iwmde_context_ensure_caches(list(
    object            = object,
    data              = data,
    posterior_samples = matrix(seq_len(24L) / 10, ncol = 1L,
      dimnames = list(NULL, "mu")),
    priors            = list()
  ))
  assign("old_likelihood", -123, context[["likelihood_cache"]])
  assign("old_row", -456, context[["row_cache"]])
  source_bytes <- serialize(object, NULL)
  set.seed(72)
  seed <- .Random.seed
  kind <- RNGkind()

  expect_identical(.iwmde_context_with_integration_control(context, NULL), context)
  expect_identical(.iwmde_context_with_integration_control(context, control), context)
  updated <- .iwmde_context_with_integration_control(context, requested)
  refreshed <- .data_selection_execution_plan(updated[["data"]])
  expect_identical(.Random.seed, seed)
  expect_identical(RNGkind(), kind)
  expect_identical(serialize(object, NULL), source_bytes)
  expect_identical(updated[["object"]], object)
  expect_identical(context[["data"]], data)
  expect_identical(updated[["posterior_samples"]], context[["posterior_samples"]])
  expect_identical(updated[["selection_spec"]], context[["selection_spec"]])
  expect_identical(updated[["formula_inputs"]], context[["formula_inputs"]])
  expect_false(identical(updated[["source_fingerprint"]], context[["source_fingerprint"]]))
  for (name in names(context)[vapply(context, is.environment, logical(1))]) {
    expect_false(identical(updated[[name]], context[[name]]))
    expect_length(ls(updated[[name]], all.names = TRUE), 0L)
  }
  expect_identical(get("old_likelihood", context[["likelihood_cache"]]), -123)
  expect_identical(get("old_row", context[["row_cache"]]), -456)
  unchanged <- setdiff(names(plan), c(
    "max_points_per_scramble", "factor_max_points_per_proposal", "designs"
  ))
  expect_identical(refreshed[unchanged], plan[unchanged])
  for (name in names(plan[["designs"]])) {
    original <- plan[["designs"]][[name]]
    expect_identical(
      refreshed[["designs"]][[name]][, seq_len(dim(original)[2L]), , drop = FALSE],
      original
    )
  }
  expect_identical(
    .iwmde_context_with_integration_control(updated, requested), updated
  )
  cluster_context <- context
  attr(cluster_context[["data"]], "known_V") <- FALSE
  cluster_updated <- .iwmde_context_with_integration_control(cluster_context, requested)
  expect_identical(
    .data_selection_execution_plan(cluster_updated[["data"]]),
    .data_selection_execution_plan(updated[["data"]])
  )
})


test_that("plot and hypothesis helpers pass the explicit post-fit control", {

  control <- set_selection_likelihood_control(max_points_per_scramble = 32768)
  captured <- list()
  testthat::local_mocked_bindings(
    .iwmde_context = function(object, integration_control = NULL) {
      expect_identical(integration_control, control)
      list()
    },
    .iwmde_estimate = function(context, parameter, density_method, density_control,
                              outputs, ...) {
      captured[[outputs]] <<- density_control[["integration_control"]]
      stop("Captured estimator input.", call. = FALSE)
    },
    .package = "RoBMA"
  )
  values <- seq(-1, 1, length.out = 30L)
  expect_error(
    .plot_brma_attach_iwmde(
      object = list(), samples = list(mu = values), parameter = "mu",
      sample_parameter = "mu", conditional = FALSE,
      n_points = 20L, sample_budget = 20L,
      normalization_points = 20L, normalization_prob = .99,
      density_method = "qCMDE", display_grid = "adaptive",
      integration_control = control
    ),
    "Captured estimator input.", fixed = TRUE
  )
  expect_error(
    .hypothesis_brma_attach_iwmde_scalar(
      posterior = values, raw_posterior = values, context = list(),
      estimate_cache = NULL, parameter = "mu", parameter_label = "mu",
      value = 0, conditional = FALSE, n_points = 20L, samples = 20L,
      target_relative_mcse = .05, normalization_points = 20L,
      normalization_prob = .99, density_method = "qCMDE",
      integration_control = control
    ),
    "Captured estimator input.", fixed = TRUE
  )
  expect_identical(captured, list(density = control, ordinate = control))
})
