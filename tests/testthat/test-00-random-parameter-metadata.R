make_random_metadata_prior <- function(
    quantity, label, parameter = "mu", block = NA_character_,
    grouping = NA_character_, allocation = NA_character_,
    component = NA_character_) {

  prior <- BayesTools::prior(
    distribution = "normal",
    parameters   = list(mean = 0, sd = 1)
  )
  attr(prior, "random_summary")         <- quantity
  attr(prior, "random_summary_label")   <- label
  attr(prior, "parameter")              <- parameter
  attr(prior, "random_block")           <- block
  attr(prior, "random_grouping_factor") <- grouping
  attr(prior, "random_allocation")      <- allocation
  attr(prior, "random_component")       <- component
  prior
}


make_random_metadata_catalog <- function(summaries, mappings) {

  quantities <- lapply(seq_along(summaries), function(i) {
    mapping <- mappings[[i]]
    BayesTools:::.bt_parameter_catalog_quantity(
      canonical_name = names(summaries)[i],
      namespace = "mu",
      role = mapping[["role"]],
      formula_parameter = "mu",
      owner_type = mapping[["owner_type"]],
      owner_name = mapping[["owner_name"]],
      quantity = mapping[["quantity"]],
      arguments = mapping[["arguments"]],
      source_type = mapping[["source_type"]],
      extraction_key = list(
        type = "random_summary",
        dependencies = character(),
        source_type = mapping[["source_type"]],
        source_parameter = mapping[["source_parameter"]],
        source_prior = mapping[["source_prior"]],
        source_transform = mapping[["source_transform"]],
        source_scale = mapping[["source_scale"]]
      )
    )
  })
  list(quantities = do.call(rbind, quantities))
}


test_that("the random semantic interface exposes every public quantity", {

  expect_identical(
    .brma_random_parameter_supported_quantities(),
    c(
      "sd", "var", "sd_total", "var_total", "sd_common", "var_common",
      "cor", "var_prop", "var_ratio", "sd_ratio"
    )
  )
})


test_that("random-parameter sources are consumed from the BayesTools catalog", {

  summaries <- list(
    "(mu) total: sd_total" = make_random_metadata_prior(
      "sd_total", "total: sd_total", allocation = "total"
    ),
    "(mu) total: var_prop(a)" = make_random_metadata_prior(
      "var_prop", "total: var_prop(a)", allocation = "total", component = "a"
    ),
    "(mu) study: sd(a)" = make_random_metadata_prior(
      "sd", "study: sd(a)", block = "study", grouping = "study", component = "a"
    ),
    "(mu) study: cor" = make_random_metadata_prior(
      "cor", "study: cor", block = "study", grouping = "study"
    ),
    "(mu) study: cor(a,b)" = make_random_metadata_prior(
      "cor", "study: cor(a,b)", block = "study", grouping = "study"
    )
  )
  mapping <- function(role, owner_type, owner_name, quantity,
                      arguments = character(), source_type = "identity",
                      source_parameter = "", source_prior = "",
                      source_transform = "identity", source_scale = 1,
                      allocation_derived = FALSE) {
    list(
      role = role, owner_type = owner_type, owner_name = owner_name,
      quantity = quantity, arguments = arguments, source_type = source_type,
      source_parameter = source_parameter, source_prior = source_prior,
      source_transform = source_transform, source_scale = source_scale,
      allocation_derived = allocation_derived
    )
  }
  mappings <- list(
    mapping("random_sd_total", "variance_allocation", "total", "sd_total",
            source_parameter = "total_source", source_prior = "total_source"),
    mapping("random_var_prop", "variance_allocation", "total", "var_prop", "a",
            source_parameter = "weight", source_prior = "weight"),
    mapping("random_sd", "random_block", "study", "sd", "a",
            source_parameter = "sd_a", source_prior = "sd_a"),
    mapping("random_correlation", "random_block", "study", "cor",
            source_parameter = "study_rho_z", source_prior = "study_rho_z",
            source_type = "one_to_one_transform",
            source_transform = "fisher_z", source_scale = NA_real_),
    mapping("random_correlation", "random_block", "study", "cor", c("a", "b"),
            source_type = "composite", source_transform = "lkj",
            source_scale = NA_real_)
  )
  catalog <- make_random_metadata_catalog(summaries, mappings)
  specs <- .brma_random_parameter_specs(catalog[["quantities"]])

  expect_equal(
    specs[["source_parameter"]],
    c("total_source", "weight", "sd_a", "study_rho_z", "")
  )
  expect_equal(
    specs[["source_type"]],
    c("identity", "identity", "identity", "one_to_one_transform", "composite")
  )
  expect_equal(specs[["owner_name"]], c("total", "total", rep("study", 3L)))
  expect_identical(specs[["source_transform"]][4L], "fisher_z")
  expect_false(any(specs[["allocation_derived"]]))
})


test_that("selected random extraction materializes only its catalog row", {

  summaries <- list(
    "(mu) study: sd(a)" = make_random_metadata_prior(
      "sd", "study: sd(a)", block = "study", component = "a"
    ),
    "(mu) study: sd(b)" = make_random_metadata_prior(
      "sd", "study: sd(b)", block = "study", component = "b"
    )
  )
  mappings <- lapply(c("sd_a", "sd_b"), function(source) {
    list(
      role             = "random_sd",
      owner_type       = "random_block",
      owner_name       = "study",
      quantity         = "sd",
      arguments        = sub("sd_", "", source),
      source_type      = "identity",
      source_parameter = source,
      source_prior     = source,
      source_transform = "identity",
      source_scale     = 1
    )
  })
  quantities <- make_random_metadata_catalog(
    summaries,
    mappings
  )[["quantities"]]
  selection <- structure(
    list(quantities = quantities[2L, , drop = FALSE]),
    class = c("BayesTools_parameter_selection", "list")
  )
  fit <- coda::mcmc.list(coda::mcmc(matrix(1:6, ncol = 2L)))
  extracted_ids <- character()
  testthat::local_mocked_bindings(
    parameter_catalog = function(...) {
      stop("the full catalog must not be traversed")
    },
    parameter_draws = function(object, selection, ...) {
      extracted_ids <<- c(
        extracted_ids,
        selection[["quantities"]][["canonical_name"]]
      )
      matrix(c(0.2, 0.3, 0.4), ncol = 1L)
    },
    parameter_transform = function(...) list(type = "identity"),
    .package = "BayesTools"
  )

  out <- .brma_random_parameter_extract_fit(
    fit,
    selections = list(selection)
  )

  expect_identical(extracted_ids, "(mu) study: sd(b)")
  expect_identical(colnames(out[["samples"]]), "(mu) study: sd(b)")
  expect_identical(out[["specs"]][["source_parameter"]], "sd_b")
})


test_that("random-parameter metadata preserves semantic support", {

  simplex <- BayesTools::prior(
    "dirichlet",
    parameters = list(alpha = c(1, 1))
  )
  expect_equal(
    .brma_random_parameter_support(list(quantity = "var_prop"),
                                   source_prior = simplex),
    c(0, 1)
  )

  expect_equal(
    .brma_random_parameter_support(
      list(quantity = "sd_ratio"),
      allocation = list(scale = "mean_variance", n_targets = 4L)
    ),
    c(0, 2)
  )
  expect_equal(
    .brma_random_parameter_support(
      list(quantity = "sd_ratio"),
      allocation = list(scale = "total_variance", n_targets = 4L)
    ),
    c(0, 1)
  )
  expect_equal(
    .brma_random_parameter_support(list(quantity = "sd_ratio")),
    c(0, Inf)
  )
})


test_that("random-parameter support consumes BayesTools transform descriptors", {

  beta <- BayesTools::prior("beta", list(alpha = 1, beta = 1))
  expect_equal(
    .brma_random_parameter_support(
      list(
        quantity = "cor",
        display_transform = list(type = "affine", offset = -1, scale = 2)
      ),
      source_prior = beta
    ),
    c(-1, 1)
  )
  uniform_z <- BayesTools::prior("uniform", list(a = -1, b = 1))
  expect_equal(
    .brma_random_parameter_support(
      list(
        quantity = "cor",
        display_transform = list(type = "tanh")
      ),
      source_prior = uniform_z
    ),
    tanh(c(-1, 1))
  )
})


test_that("independent factor coefficients retain separate plot levels", {

  group <- matrix(rnorm(20), ncol = 2L)
  class(group) <- c("mixed_posteriors.factor", class(group))
  attr(group, "independent") <- TRUE
  attr(group, "level_names") <- c("sensitivity", "specificity")
  samples <- list(group = group)

  n_levels <- .get_samples_n_levels(samples, "group")
  dots     <- .set_dots_plot(n_levels = n_levels)

  expect_identical(n_levels, 2L)
  expect_gte(length(unique(dots[["col"]])), 2L)
})


test_that("random qCMDE targets retain stored semantic coordinates", {

  z <- c(-0.5, 0, 0.5)
  fit <- structure(
    list(mcmc = matrix(z, ncol = 1L, dimnames = list(NULL, "study_rho_z"))),
    prior_list = list(
      study_rho_z = BayesTools::prior(
        "normal",
        parameters = list(mean = 0, sd = 1)
      )
    )
  )
  selected <- list(
    spec = list(
      quantity     = "cor",
      source_type      = "one_to_one_transform",
      source_parameter = "study_rho_z",
      source_transform = "fisher_z",
      display_transform = list(type = "tanh"),
      label            = "study: cor"
    ),
    samples = matrix(NA_real_, nrow = length(z), ncol = 1L)
  )
  testthat::local_mocked_bindings(
    .brma_random_parameter_select = function(...) selected,
    .package = "RoBMA"
  )

  target <- .brma_random_parameter_density_target(list(fit = fit), "study: cor")
  expect_identical(target[["parameter"]], "study_rho_z")
  expect_identical(target[["parameter_spec"]][["type"]], "primitive")
  expect_identical(target[["display_transform"]], list(type = "tanh"))
})


test_that("bivariate LKJ correlations expose their scalar qCMDE coordinate", {

  probability <- c(0.25, 0.5, 0.75)
  fit <- structure(
    list(mcmc = matrix(
      probability,
      ncol = 1L,
      dimnames = list(NULL, "study_lkj_probability")
    )),
    prior_list = list(
      study_lkj_probability = BayesTools::prior(
        "beta",
        parameters = list(alpha = 1, beta = 1)
      )
    )
  )
  selected <- list(
    spec = list(
      quantity     = "cor",
      source_type      = "one_to_one_transform",
      source_parameter = "study_lkj_probability",
      source_transform = "lkj2",
      display_transform = list(type = "affine", offset = -1, scale = 2),
      label            = "study: cor(group[sensitivity],group[specificity])"
    ),
    samples = matrix(NA_real_, nrow = length(probability), ncol = 1L)
  )
  testthat::local_mocked_bindings(
    .brma_random_parameter_select = function(...) selected,
    .package = "RoBMA"
  )

  target <- .brma_random_parameter_density_target(
    list(fit = fit),
    "study: cor(group[sensitivity],group[specificity])"
  )
  expect_identical(target[["parameter"]], "study_lkj_probability")
  expect_identical(target[["parameter_spec"]][["type"]], "primitive")
  expect_identical(
    target[["display_transform"]],
    list(type = "affine", offset = -1, scale = 2)
  )
})


test_that("variance aggregates expose squared-SD qCMDE coordinates", {

  common_sd <- c(0.25, 0.5, 1)
  fit <- structure(
    list(mcmc = matrix(
      common_sd,
      ncol = 1L,
      dimnames = list(NULL, "heterogeneity_common_sd")
    )),
    prior_list = list(
      heterogeneity_common_sd = BayesTools::prior(
        "normal",
        parameters = list(mean = 0, sd = 1),
        truncation = list(lower = 0, upper = Inf)
      )
    )
  )
  selected <- list(
    spec = list(
      quantity         = "var_common",
      source_type      = "one_to_one_transform",
      source_parameter = "heterogeneity_common_sd",
      source_transform = "square",
      display_transform = list(type = "square"),
      label            = "var_common"
    ),
    samples = matrix(NA_real_, nrow = length(common_sd), ncol = 1L)
  )
  testthat::local_mocked_bindings(
    .brma_random_parameter_select = function(...) selected,
    .package = "RoBMA"
  )

  target <- .brma_random_parameter_density_target(
    list(fit = fit),
    "var_common"
  )
  expect_identical(target[["parameter"]], "heterogeneity_common_sd")
  expect_identical(target[["parameter_spec"]][["type"]], "primitive")
  expect_identical(target[["display_transform"]], list(type = "square"))
})


test_that("random qCMDE targets support general simplex allocations", {

  eta <- rbind(
    c(1, 2, 3),
    c(3, 1, 2),
    c(2, 3, 1)
  )
  weights <- eta / rowSums(eta)
  colnames(weights) <- paste0("allocation[", 1:3, "]")
  colnames(eta) <- .iwmde_simplex_auxiliary_columns("allocation", 3L)
  posterior <- cbind(weights, eta)
  source_prior <- BayesTools::prior(
    "dirichlet",
    parameters = list(alpha = c(1, 2, 3))
  )
  summary_prior <- make_random_metadata_prior(
    "var_ratio",
    "allocation: var_ratio(study)",
    allocation = "allocation",
    component  = "study"
  )
  attr(summary_prior, "random_allocation_metadata") <- list(
    scale       = "mean_variance",
    n_targets   = 3L,
    weight_name = "allocation"
  )
  attr(summary_prior, "random_allocation_index") <- 2L
  fit <- structure(
    list(mcmc = posterior),
    prior_list = list(allocation = source_prior)
  )
  selected <- list(
    spec = list(
      quantity          = "var_ratio",
      evaluator         = "allocation",
      source_type       = "one_to_one_transform",
      source_parameter = "allocation",
      source_transform = "var_ratio",
      display_transform = list(type = "affine", offset = 0, scale = 3),
      label             = "allocation: var_ratio(study)",
      allocation_index  = 2L
    ),
    samples      = matrix(NA_real_, nrow = nrow(weights), ncol = 1L),
    prior        = summary_prior,
    source_prior = source_prior,
    allocation_definition = list(
      scale       = "mean_variance",
      n_targets   = 3L,
      weight_name = "allocation"
    )
  )
  testthat::local_mocked_bindings(
    .brma_random_parameter_select = function(...) selected,
    .package = "RoBMA"
  )

  target <- .brma_random_parameter_density_target(
    list(fit = fit),
    "allocation: var_ratio(study)"
  )
  expect_identical(target[["parameter"]], "allocation[2]")
  expect_identical(target[["parameter_spec"]][["type"]], "simplex_pair")
  expect_identical(target[["parameter_spec"]][["n_targets"]], 3L)
  expect_identical(
    target[["display_transform"]],
    list(type = "affine", offset = 0, scale = 3)
  )
})


test_that("allocated component SD targets use declared catalog provenance", {

  posterior <- cbind(
    tau = c(0.5, 0.8),
    `weight[1]` = c(0.25, 0.4),
    `weight[2]` = c(0.75, 0.6),
    `prior_par_eta_weight[1]` = c(1, 2),
    `prior_par_eta_weight[2]` = c(3, 3)
  )
  allocation <- list(
    target               = "sd_component",
    source               = list(name = "tau", shape = "scalar"),
    parent_factors       = list(),
    weight_name          = "weight",
    scale                = "mean_variance",
    n_targets            = 2L,
    leaf_index_by_column = 1:2
  )
  term <- list(
    block_name         = "study",
    group_label        = "study",
    structure          = "diag",
    sd_component_terms = c("a", "b"),
    sd_binding = list(
      true_allocation = TRUE,
      allocations     = list(allocation)
    )
  )
  fit <- structure(
    list(mcmc = posterior),
    formula_design = list(list(
      parameter      = "mu",
      random_effects = list(term)
    ))
  )
  selected <- list(
    spec = list(
      quantity           = "sd",
      evaluator          = "sd",
      source_type        = "composite",
      allocation_derived = TRUE,
      formula_parameter  = "mu",
      block              = "study",
      grouping           = "",
      random_component   = "a",
      label              = "study: sd(a)",
      display_transform  = NULL
    ),
    samples = matrix(NA_real_, nrow = nrow(posterior), ncol = 1L)
  )
  testthat::local_mocked_bindings(
    .brma_random_parameter_select = function(...) selected,
    .package = "RoBMA"
  )

  target <- .brma_random_parameter_density_target(
    list(fit = fit),
    "study: sd(a)"
  )
  expect_identical(target[["parameter"]], "tau")
  expect_identical(
    target[["parameter_spec"]][["type"]],
    "random_component_sd"
  )
  expect_identical(
    target[["parameter_spec"]][["factor_columns"]],
    "weight[1]"
  )

  selected[["spec"]][["allocation_derived"]] <- FALSE
  unsupported <- .brma_random_parameter_density_target(
    list(fit = fit),
    "study: sd(a)"
  )
  expect_match(unsupported[["reason"]], "no supported scalar")
})


test_that("compiled allocation metadata takes precedence over definitions", {

  compiled <- list(
    label       = "heterogeneity",
    weight_name = "heterogeneity_weight",
    scale       = "mean_variance",
    n_targets   = 2L
  )
  definition <- list(
    label       = "heterogeneity",
    weight_name = "heterogeneity_weight",
    scale       = "mean_variance"
  )
  formula_design <- list(list(
    parameter          = "mu",
    random_allocations = list(definition),
    random_effects     = list(list(
      sd_binding = list(allocations = list(compiled))
    ))
  ))

  actual <- .brma_random_parameter_design_allocation(
    formula_design,
    list(allocation = "heterogeneity", formula_parameter = "mu")
  )
  expect_identical(actual, compiled)
})
