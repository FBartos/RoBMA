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


# Label parts of a random-effect quantity '(mu) <owner>: <quantity>(<args>)'.
random_metadata_parts <- function(selector, owner, quantity,
                                  arguments = character(),
                                  display_arguments = arguments) {

  catalog_label_parts(
    selector, if (nzchar(owner)) owner else quantity, "mu",
    random = list(
      owner             = owner,
      quantity          = quantity,
      arguments         = arguments,
      display_arguments = display_arguments
    )
  )
}


make_random_metadata_catalog <- function(summaries, mappings) {

  quantities <- lapply(seq_along(summaries), function(i) {
    mapping <- mappings[[i]]
    BayesTools:::.bt_parameter_catalog_quantity(
      canonical_name = names(summaries)[i],
      namespace = "mu",
      role = mapping[["role"]],
      formula_parameter = "mu",
      label_parts = random_metadata_parts(
        names(summaries)[i],
        mapping[["owner_name"]],
        mapping[["quantity"]],
        mapping[["arguments"]]
      ),
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


mock_random_marginal_update_plan <- function() {

  plan <- structure(
    list(
      family = "unsupported",
      reason = "metadata_target_fixture"
    ),
    class = c(
      "BayesTools_random_effects_marginal_update_plan",
      "list"
    )
  )
  testthat::local_mocked_bindings(
    random_effects_marginal_update_plan = function(...) plan,
    .package = "BayesTools",
    .env = parent.frame()
  )
}


test_that("the random semantic interface exposes every public quantity", {

  expect_identical(
    .brma_random_parameter_supported_quantities(),
    c(
      "sd", "var", "sd_total", "var_total", "sd_common", "var_common",
      "cor", "var_prop", "var_mult", "sd_mult"
    )
  )
})


test_that("RoBMA renders BayesTools random quantities by its I/O names", {

  expect_identical(
    .brma_random_parameter_io_quantity_map(),
    c(
      sd         = "tau",
      var        = "tau2",
      sd_total   = "tau_total",
      var_total  = "tau2_total",
      sd_common  = "tau_common",
      var_common = "tau2_common",
      cor        = "rho",
      var_prop   = "tau2_prop",
      sd_mult    = "tau_mult",
      var_mult   = "tau2_mult"
    )
  )
  quantities <- data.frame(
    canonical_name = c(
      "(mu) study: sd(intercept)", "(mu) var_total", "(mu) cor(a,b)",
      "(mu) allocation: var_prop(study)", "(mu) sd_common"
    ),
    stringsAsFactors = FALSE
  )
  quantities[["label_parts"]] <- I(list(
    random_metadata_parts(
      "(mu) study: sd(intercept)", "study", "sd", "intercept", character()
    ),
    random_metadata_parts("(mu) var_total", "", "var_total"),
    random_metadata_parts("(mu) cor(a,b)", "", "cor", c("a", "b")),
    random_metadata_parts(
      "(mu) allocation: var_prop(study)", "allocation", "var_prop", "study"
    ),
    random_metadata_parts("(mu) sd_common", "", "sd_common")
  ))

  expect_identical(
    .brma_random_parameter_io_labels(quantities, "selector"),
    c(
      "(mu) study: tau(intercept)", "(mu) tau2_total", "(mu) rho(a,b)",
      "(mu) allocation: tau2_prop(study)", "(mu) tau_common"
    )
  )
  expect_identical(
    .brma_random_parameter_io_labels(quantities, "label"),
    c(
      "study: tau", "tau2_total", "rho(a,b)",
      "allocation: tau2_prop(study)", "tau_common"
    )
  )
  expect_identical(
    .brma_random_parameter_io_aliases(as.list(quantities[1L, , drop = FALSE])),
    data.frame(
      alias      = c(
        "study: tau(intercept)", "(mu) study: tau", "study: tau", "tau"
      ),
      simplified = c(FALSE, TRUE, TRUE, TRUE),
      stringsAsFactors = FALSE
    )
  )
})


test_that("public draws replace backend random coordinates with RoBMA names", {

  coordinate_values <- cbind(
    mu                 = 1:3,
    backend_random_sd  = c(.2, .3, .4),
    random_inclusion   = c(1, 0, 1),
    backend_random_z   = c(-1, 0, 1)
  )
  semantic_values <- cbind(
    backend_tau  = coordinate_values[, "backend_random_sd"],
    backend_tau2 = coordinate_values[, "backend_random_sd"]^2
  )
  coordinate_chain <- coda::mcmc(coordinate_values, start = 5, thin = 2)
  semantic_chain   <- coda::mcmc(semantic_values, start = 5, thin = 2)
  object <- structure(list(fit = list()), class = "brma")

  testthat::local_mocked_bindings(
    parameter_coordinates = function(...) data.frame(
      coordinate_name = colnames(coordinate_values),
      role = c(
        "fixed_coefficient", "random_sd", "random_inclusion", "random_latent"
      )
    ),
    .package = "BayesTools"
  )
  # No random-inclusion gates in this synthetic catalog, so the gate columns
  # are a no-op and the coordinate replacement is what is under test.
  testthat::local_mocked_bindings(
    parameter_catalog = function(...) list(
      quantities = data.frame(
        role             = character(),
        stringsAsFactors = FALSE
      )
    ),
    .package = "BayesTools"
  )
  testthat::local_mocked_bindings(
    .brma_random_parameter_bundle = function(...) list(
      samples = coda::mcmc.list(semantic_chain),
      specs   = data.frame(label = c("tau", "tau2"))
    ),
    .package = "RoBMA"
  )

  out <- .brma_replace_random_coordinates(
    object,
    coda::mcmc.list(coordinate_chain)
  )
  values <- as.matrix(out[[1L]])

  expect_identical(
    colnames(values),
    c("mu", "random_inclusion", "backend_random_z", "tau", "tau2")
  )
  expect_equal(values[, "tau"], coordinate_values[, "backend_random_sd"])
  expect_equal(values[, "tau2"], coordinate_values[, "backend_random_sd"]^2)
  expect_identical(coda::mcpar(out[[1L]]), coda::mcpar(coordinate_chain))
})


test_that("random summary rows use the RoBMA I/O map", {

  estimates <- matrix(
    c(.2, .1, .4, .3),
    nrow = 2L,
    dimnames = list(c("sd", "cor(a,b)"), c("Mean", "SD"))
  )
  attr(estimates, "parameters") <- c("(mu) sd", "(mu) cor(a,b)")
  quantities <- data.frame(
    canonical_name = c("(mu) sd", "(mu) cor(a,b)"),
    role           = c("random_sd", "random_correlation"),
    quantity       = c("sd", "cor")
  )
  quantities[["label_parts"]] <- I(list(
    random_metadata_parts("(mu) sd", "", "sd"),
    random_metadata_parts("(mu) cor(a,b)", "", "cor", c("a", "b"))
  ))
  testthat::local_mocked_bindings(
    parameter_catalog = function(...) list(quantities = quantities),
    .package = "BayesTools"
  )

  out <- .summary_random_repair_parameter_names(
    estimates,
    list(fit = list())
  )

  expect_identical(rownames(out), c("tau", "rho(a,b)"))
  expect_identical(
    attr(out, "parameters"),
    c("(mu) tau", "(mu) rho(a,b)")
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


test_that("random extraction avoids redundant catalog and dependency work", {

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
  extraction_keys <- quantities[["extraction_key"]]
  for (i in seq_along(extraction_keys)) {
    extraction_keys[[i]][["dependencies"]] <- c("sd_a", "sd_b")[i]
  }
  quantities[["extraction_key"]] <- I(extraction_keys)
  selections <- lapply(seq_len(nrow(quantities)), function(i) {
    structure(
      list(quantities = quantities[i, , drop = FALSE]),
      class = c("BayesTools_parameter_selection", "list")
    )
  })
  fit <- coda::mcmc.list(coda::mcmc(matrix(1:6, ncol = 2L)))
  extracted_ids             <- character()
  supplied_model_samples    <- logical()
  materialized_dependencies <- list()
  testthat::local_mocked_bindings(
    parameter_catalog = function(...) {
      stop("the full catalog must not be traversed")
    },
    JAGS_materialize_draws = function(object, parameters, ...) {
      materialized_dependencies[[length(materialized_dependencies) + 1L]] <<-
        parameters
      values <- matrix(
        seq_len(3L * length(parameters)),
        nrow = 3L,
        dimnames = list(NULL, parameters)
      )
      coda::mcmc.list(coda::mcmc(values))
    },
    parameter_draws = function(object, selection, model_samples = NULL, ...) {
      extracted_ids <<- c(
        extracted_ids,
        selection[["quantities"]][["canonical_name"]]
      )
      supplied_model_samples <<- c(
        supplied_model_samples,
        !is.null(model_samples)
      )
      matrix(c(0.2, 0.3, 0.4), ncol = 1L)
    },
    parameter_transform = function(...) list(type = "identity"),
    .package = "BayesTools"
  )

  out <- .brma_random_parameter_extract_fit(
    fit,
    selections = selections[2L]
  )

  expect_identical(extracted_ids, "(mu) study: sd(b)")
  expect_length(materialized_dependencies, 0L)
  expect_identical(supplied_model_samples, FALSE)
  expect_identical(colnames(out[["samples"]]), "(mu) study: tau(b)")
  expect_identical(out[["specs"]][["source_parameter"]], "sd_b")

  .brma_random_parameter_extract_fit(fit, selections = selections)

  expect_identical(
    materialized_dependencies,
    list(c("sd_a", "sd_b"))
  )
  expect_identical(
    extracted_ids,
    c("(mu) study: sd(b)", "(mu) study: sd(a)", "(mu) study: sd(b)")
  )
  expect_identical(supplied_model_samples, c(FALSE, TRUE, TRUE))
})


test_that("random-parameter support is the catalog's declared support", {

  selected <- function(support) {
    quantities <- data.frame(quantity_id = "q", stringsAsFactors = FALSE)
    quantities[["support"]] <- I(list(support))
    list(entry = list(
      parameter = "study: tau",
      selection = list(quantities = quantities)
    ))
  }

  exact <- BayesTools::posterior_support_attribute(c(0, 2))
  expect_identical(.brma_random_parameter_catalog_support(selected(exact)), exact)
  expect_equal(.brma_random_parameter_support(selected(exact)), c(0, 2))

  # Plotting limits and underivable supports claim no exact bounds.
  limits <- BayesTools::posterior_support_attribute(c(0, 2), exact = FALSE)
  expect_equal(.brma_random_parameter_support(selected(limits)), c(-Inf, Inf))
  expect_null(.brma_random_parameter_catalog_support(selected(NULL)))
  expect_equal(.brma_random_parameter_support(selected(NULL)), c(-Inf, Inf))

  expect_error(
    .brma_random_parameter_support(list(entry = list(parameter = "study: tau"))),
    "Random-effect quantity 'study: tau' has no catalog support metadata.",
    fixed = TRUE
  )
})


test_that("random-effect diagnostics bound densities by the catalog's plotting limits", {

  # A component SD allocated from a scale prior that is bounded above: the
  # catalog declares its hull [0, Inf) as plotting limits (exact = FALSE),
  # which bound the diagnostic density but no hypothesis.
  quantities <- data.frame(quantity_id = "q", stringsAsFactors = FALSE)
  quantities[["support"]] <- I(list(
    BayesTools::posterior_support_attribute(c(0, Inf), exact = FALSE)
  ))
  selected <- list(
    entry   = list(
      parameter = "(mu) study: tau(intercept)",
      selection = list(quantities = quantities)
    ),
    spec    = list(label = "study: tau"),
    samples = matrix(c(0.1, 0.2), ncol = 1L)
  )
  testthat::local_mocked_bindings(
    .brma_random_parameter_select           = function(...) selected,
    .brma_random_parameter_fit_with_samples = function(fit, samples) list(),
    .package = "RoBMA"
  )

  diagnostic <- .brma_random_parameter_diagnostic_fit(
    list(fit = list()),
    "(mu) study: tau(intercept)"
  )
  prior <- attr(diagnostic[["fit"]], "prior_list")[[diagnostic[["parameter"]]]]
  expect_equal(unlist(prior[["truncation"]]), c(lower = 0, upper = Inf))

  expect_equal(.brma_random_parameter_support(selected), c(-Inf, Inf))
  expect_equal(.brma_random_parameter_support(selected, limits = TRUE), c(0, Inf))
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

  mock_random_marginal_update_plan()
  z <- c(-0.5, 0, 0.5)
  fit <- structure(
    coda::mcmc(matrix(z, ncol = 1L, dimnames = list(NULL, "study_rho_z"))),
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

  mock_random_marginal_update_plan()
  probability <- c(0.25, 0.5, 0.75)
  fit <- structure(
    coda::mcmc(matrix(
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

  mock_random_marginal_update_plan()
  common_sd <- c(0.25, 0.5, 1)
  fit <- structure(
    coda::mcmc(matrix(
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

  mock_random_marginal_update_plan()
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
    "var_mult",
    "allocation: var_mult(study)",
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
    coda::mcmc(posterior),
    prior_list = list(allocation = source_prior)
  )
  selected <- list(
    spec = list(
      quantity          = "var_mult",
      evaluator         = "allocation",
      source_type       = "one_to_one_transform",
      source_parameter = "allocation",
      source_transform = "var_mult",
      display_transform = list(type = "affine", offset = 0, scale = 3),
      label             = "allocation: var_mult(study)",
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
    "allocation: var_mult(study)"
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

  mock_random_marginal_update_plan()
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
    leaf_index_by_column = 1:2,
    leaf_names           = c("mu__xREx__study_a", "mu__xREx__study_b")
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
    coda::mcmc(posterior),
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
  expect_identical(target[["parameter_spec"]][["node"]], "mu__xREx__study_a")

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
