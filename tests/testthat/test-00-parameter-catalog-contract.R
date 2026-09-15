context("BayesTools parameter catalog contract")


.test_parameter_catalog <- function() {

  quantity <- function(canonical_name, namespace, role,
                       formula_parameter = "", term = "", component = "",
                       display_label = canonical_name, status = "sampled",
                       fixed_value = NA_real_, extraction_key = NULL,
                       owner_type = "",
                       owner_name = "", quantity_name = "",
                       arguments = character(), source_type = "none") {

    if (is.null(extraction_key)) {
      extraction_key <- list(
        type         = "coordinate",
        dependencies = canonical_name
      )
    }
    BayesTools:::.bt_parameter_catalog_quantity(
      canonical_name    = canonical_name,
      namespace         = namespace,
      role              = role,
      formula_parameter = formula_parameter,
      term              = term,
      component         = component,
      display_label     = display_label,
      display_scale     = "original",
      status            = status,
      fixed_value       = fixed_value,
      owner_type        = owner_type,
      owner_name        = owner_name,
      quantity          = quantity_name,
      arguments         = arguments,
      source_type       = source_type,
      extraction_key    = extraction_key
    )
  }
  quantities <- do.call(rbind, list(
    quantity("mu_x", "mu", "fixed_coefficient", "mu", "x", "mu_x"),
    quantity(
      "mu_f[1]", "mu", "fixed_coefficient", "mu", "f", "A",
      "(mu) f[A]"
    ),
    quantity(
      "mu_f[2]", "mu", "fixed_coefficient", "mu", "f", "B",
      "(mu) f[B]"
    ),
    quantity(
      "mu_f[dif: C]", "mu", "fixed_coefficient", "mu", "f", "C",
      "(mu) f[C]", "structural", 0,
      list(type = "factor_level", dependencies = character(), weights = numeric())
    ),
    quantity(
      "log_tau_x", "log_tau", "fixed_coefficient", "log_tau", "x",
      "log_tau_x"
    ),
    quantity("PET", "model", "parameter"),
    quantity(
      "random_sd_study_intercept", "mu", "random_sd", "mu", "intercept",
      "intercept", "(mu) study: sd",
      owner_type    = "random_block",
      owner_name    = "study",
      quantity_name = "sd",
      arguments     = "intercept",
      source_type   = "identity"
    ),
    quantity(
      "random_cor_study_group", "mu", "random_correlation",
      formula_parameter = "mu",
      term              = "study",
      display_label     = paste0(
        "(mu) cor(group[sensitivity],group[specificity])"
      ),
      owner_type        = "random_block",
      owner_name        = "study",
      quantity_name     = "cor",
      arguments         = c(
        "group[sensitivity]", "group[specificity]"
      ),
      source_type       = "identity"
    )
  ))
  out <- BayesTools:::.bt_parameter_catalog_new(
    quantities = quantities,
    aliases    = BayesTools:::.bt_parameter_catalog_aliases(quantities)
  )
  return(out)
}


.test_formula_name_map <- function(parameter) {

  out <- data.frame(
    encoded_name      = c("mu_x", "mu_f", "log_tau_x"),
    jags_name         = c("mu_x", "mu_f", "log_tau_x"),
    kind              = "fixed",
    formula_parameter = c("mu", "mu", "log_tau"),
    term              = c("x", "f", "x"),
    role              = "coefficient",
    stringsAsFactors  = FALSE
  )
  out <- out[out[["formula_parameter"]] == parameter, , drop = FALSE]
  class(out) <- c("BayesTools_formula_name_map", "data.frame")
  attr(out, "schema_version") <- 1L
  return(out)
}


test_that("fitted parameter discovery is metadata-only and component-aware", {

  scale <- data.frame(x = c(0, 1))
  attr(scale, "parameter") <- "log_tau"
  attr(scale, "source")    <- "tau"
  attr(scale, "aliases")   <- "tau"
  data <- list(
    outcome  = data.frame(yi = c(0, 1), sei = c(1, 1)),
    location = data.frame(x = c(0, 1)),
    scale    = scale
  )
  attr(data, "mods")   <- FALSE
  attr(data, "random") <- TRUE
  attr(data, "scale")  <- TRUE
  object <- list(
    data   = data,
    fit    = structure(list(sentinel = TRUE), class = "BayesTools_fit"),
    priors = list(outcome = list(
      bias = BayesTools::prior_PET("normal", list(mean = 0, sd = 1))
    ))
  )
  checked <- NULL
  testthat::local_mocked_bindings(
    JAGS_validate_fit_contract = function(fit, requires) {
      checked <<- unique(c(checked, requires))
      invisible(TRUE)
    },
    parameter_catalog = function(object, ...) .test_parameter_catalog(),
    # Catalog metadata requires the fitted parameter map unconditionally; the
    # sentinel fit carries none, so stand in a map that exposes the runtime
    # cache the metadata builder reuses.
    parameter_map = function(object, ...) {
      structure(list(), runtime_cache = new.env(parent = emptyenv()))
    },
    JAGS_formula_name_map = function(fit, parameter) {
      .test_formula_name_map(parameter)
    },
    .package = "BayesTools"
  )
  testthat::local_mocked_bindings(
    .brma_random_parameter_bundle = function(...) {
      stop("posterior discovery is forbidden")
    },
    .package = "RoBMA"
  )

  catalog <- .brma_parameter_catalog(object)
  mods    <- .brma_parameter_select_entry(object, "x", component = "mods")
  scale   <- .brma_parameter_select_entry(object, "x", component = "scale")
  factor  <- .brma_parameter_select_entry(object, "f", component = "mods")
  pet     <- .brma_parameter_select_entry(object, "PET", component = "bias")
  random  <- .brma_parameter_select_entry(
    object,
    "study: tau",
    component = "random"
  )
  random_short <- .brma_parameter_select_entry(
    object,
    "tau",
    component = "random"
  )

  expect_setequal(
    checked,
    c(
      "name_encoding", "formula_name_map", "formula_design",
      "parameter_map"
    )
  )
  expect_true(any(catalog[["component"]] == "random"))
  expect_identical(mods[["parameter"]], "mu_x")
  expect_identical(scale[["parameter"]], "log_tau_x")
  expect_identical(factor[["parameter"]], "mu_f")
  expect_identical(
    factor[["selection"]][["quantities"]][["provider"]],
    "RoBMA"
  )
  expect_identical(
    factor[["selection"]][["quantities"]][["role"]],
    "formula_coefficient_group"
  )
  native_factor <- factor[["selection"]][["quantities"]][["extraction_key"]][[1L]]
  expect_identical(native_factor[["dependencies"]], c("mu_f[1]", "mu_f[2]"))
  factor_hypothesis <- .hypothesis_brma_select_parameter(
    object     = object,
    hypothesis = "f[B] > f[C]",
    component  = "auto"
  )
  factor_hypothesis_explicit <- .hypothesis_brma_select_parameter(
    object     = object,
    hypothesis = "f[B] > f[C]",
    component  = "mods"
  )
  expect_identical(factor_hypothesis[["parameter"]], "mu_f")
  expect_identical(
    factor_hypothesis_explicit[["parameter"]],
    factor_hypothesis[["parameter"]]
  )
  expect_setequal(
    unique(factor_hypothesis[["resolution"]][["occurrences"]][["level"]]),
    c("B", "C")
  )
  expect_true(all(
    factor_hypothesis[["resolution"]][["occurrences"]][["quantity_id"]] %in%
      factor[["member_quantity_ids"]][[1L]]
  ))
  resolved_catalog <- .brma_parameter_catalog_metadata(object)[["catalog"]]
  level_b <- BayesTools::parameter_catalog_resolve(
    resolved_catalog,
    alias     = "f",
    component = "B"
  )
  level_c <- BayesTools::parameter_catalog_resolve(
    resolved_catalog,
    alias     = "f",
    component = "C"
  )
  expect_identical(level_b[["quantities"]][["provider"]], "BayesTools")
  expect_identical(level_b[["quantities"]][["status"]], "sampled")
  expect_identical(level_c[["quantities"]][["provider"]], "BayesTools")
  expect_identical(level_c[["quantities"]][["status"]], "structural")
  expect_identical(
    pet[["selection"]][["quantities"]][["provider"]],
    "BayesTools"
  )
  expect_identical(
    random[["selection"]][["quantities"]][["provider"]],
    "BayesTools"
  )
  expect_identical(random[["quantity"]], "sd")
  expect_identical(random[["status"]], "sampled")
  expect_identical(random_short[["canonical_name"]], random[["canonical_name"]])
  expect_error(
    .brma_parameter_select_entry(
      object,
      "study: sd",
      component = "random"
    ),
    class = "BayesTools_parameter_not_found"
  )
  random_hypothesis <- BayesTools::hypothesis_parse(
    "rho(group[sensitivity],group[specificity]) = 0",
    catalog        = resolved_catalog,
    simplify_names = TRUE
  )
  selected_random_hypothesis <- .hypothesis_brma_select_parameter(
    object     = object,
    hypothesis = random_hypothesis,
    component  = "random"
  )
  expect_identical(
    selected_random_hypothesis[["parameter"]],
    "random_cor_study_group"
  )
  expect_error(
    BayesTools::hypothesis_parse(
      "cor(group[sensitivity],group[specificity]) = 0",
      catalog        = resolved_catalog,
      simplify_names = TRUE
    ),
    "Unsupported hypothesis expression operator or function 'cor'",
    fixed = TRUE
  )
  expect_s3_class(mods[["selection"]], "BayesTools_parameter_selection")
  expect_error(
    .brma_parameter_select_entry(object, "x"),
    class = "BayesTools_parameter_ambiguous"
  )
  expect_error(
    .brma_parameter_select_entry(object, "missing", component = "mods"),
    class = "BayesTools_parameter_not_found"
  )

  materialized <- NULL
  testthat::local_mocked_bindings(
    JAGS_materialize_draws = function(fit, parameters) {
      materialized <<- list(fit = fit, parameters = parameters)
      "draws"
    },
    .package = "BayesTools"
  )
  class(object) <- "brma"
  expect_identical(
    BayesTools::parameter_draws(object, factor[["selection"]]),
    "draws"
  )
  expect_identical(materialized[["fit"]], object[["fit"]])
  expect_identical(materialized[["parameters"]], c("mu_f[1]", "mu_f[2]"))
})


test_that("indexed factor plots retain semantic cells without cached fits", {

  make_object <- function(contrast, levels) {

    dat <- data.frame(yi = c(-.1, .2, .05, .3, -.2, .1),
      study = factor(rep(letters[1:3], each = 2L)),
      group = factor(rep(levels, length.out = 6L), levels = levels))
    bias <- BayesTools::prior_weightfunction("one-sided", .025,
      BayesTools::wf_fixed(c(1, .5)), model = selection_model(
        group = "study", other_random_effects = "condition",
        known_sampling_variance = "integrate"))
    object <- bselmodel.mv(yi = yi, V = diag(.04, 6L), mods = ~ group,
      random = ~ 1 | study, data = dat, measure = "GEN",
      prior_unit_information_sd = 1,
      prior_effect = BayesTools::prior("normal", list(0, .7)),
      prior_mods = list(group = BayesTools::prior_factor(
        if (contrast == "orthonormal") "mnormal" else "normal",
        list(0, .6), contrast = contrast)),
      prior_heterogeneity = BayesTools::prior("point", list(location = .2)),
      prior_bias = bias, only_priors = TRUE, silent = TRUE)
    design <- object$formula_design$mu
    map <- design$name_map
    intercept <- map$jags_name[map$kind == "fixed" & map$term == "intercept"]
    factor_name <- map$jags_name[map$kind == "fixed" & map$term == "group"]
    factor_columns <- BayesTools:::.JAGS_prior_factor_names(
      factor_name, design$prior_list[[factor_name]])
    samples <- matrix(seq(-.1, .2, length.out = 24L), ncol = 1L,
      dimnames = list(NULL, intercept))
    for (column in factor_columns) {
      samples <- cbind(samples, seq(-.4, .5, length.out = 24L))
      colnames(samples)[ncol(samples)] <- column
    }
    term <- design$random_effects[[1L]]
    samples <- cbind(samples, rep(.2, 24L))
    colnames(samples)[ncol(samples)] <- term$sd_parameter_names
    latent <- as.vector(BayesTools:::.bt_random_effect_latent_names(
      term, n_groups = term$n_groups, n_columns = 1L))
    for (column in latent) {
      samples <- cbind(samples, rep(0, 24L))
      colnames(samples)[ncol(samples)] <- column
    }
    fit <- coda::mcmc.list(coda::mcmc(samples))
    class(fit) <- c("BayesTools_fit", class(fit))
    attr(fit, "formula_design") <- object$formula_design
    attr(fit, "prior_list") <- c(design$prior_list,
      .create_fit_priors(object$data, object$priors))
    fit <- BayesTools:::.bt_attach_parameter_map(fit)
    fit <- BayesTools:::.bt_attach_draw_geometry(fit)
    fit <- BayesTools:::.bt_attach_fit_contract(fit)
    object$fit <- fit
    class(object) <- setdiff(class(object), "only_priors.brma")
    list(object = object, samples = samples, factor_name = factor_name,
      factor_columns = factor_columns)
  }

  rendered <- NULL
  attached <- NULL
  testthat::local_mocked_bindings(
    .plot_brma_attach_iwmde = function(object, samples, parameter,
        sample_parameter, parameter_spec, ...) {
      attached <<- list(parameter = parameter, spec = parameter_spec)
      attr(samples[[sample_parameter]], "posterior_density") <- list(test_sentinel = TRUE)
      samples
    }, .package = "RoBMA")
  testthat::local_mocked_bindings(
    plot_marginal = function(samples, parameter, ...) {
      rendered <<- list(samples = samples, parameter = parameter, args = list(...), renderer = "cell")
      structure(list(), class = "mock_plot")
    },
    plot_posterior = function(samples, parameter, ...) {
      rendered <<- list(samples = samples, parameter = parameter, args = list(...), renderer = "term")
      structure(list(), class = "mock_plot")
    }, .package = "BayesTools")

  for (contrast in c("treatment", "orthonormal")) {
    fixture <- make_object(contrast, if (contrast == "treatment") letters[1:3] else letters[1:2])
    object <- fixture$object
    for (level in c("a", "b")) {
      selector <- paste0("group[", level, "]")
      entry <- .brma_parameter_select_entry(object, selector, allow_factor_cells = TRUE)
      key <- entry$selection$quantities$extraction_key[[1L]]
      # Treatment a is structural zero and b has weight +1. Two-level
      # orthonormal cells have weights +/-1/sqrt(2); metadata owns which
      # label receives each orientation, independently of posterior values.
      manual_weight <- if (contrast == "treatment") {
        if (level == "a") 0 else 1
      } else sign(key$weights[[1L]]) / sqrt(2)
      expected <- manual_weight * fixture$samples[, fixture$factor_columns[[1L]]]
      attached <- NULL
      plot(object, selector, component = "mods", prior = TRUE,
        standardized_coefficients = TRUE, density_method = "qCMDE", plot_type = "ggplot")
      expect_identical(rendered$renderer, "cell")
      expect_identical(names(rendered$samples), entry$parameter)
      expect_identical(rendered$parameter, entry$parameter)
      sample <- rendered$samples[[entry$parameter]]
      expect_identical(attr(sample, "level_name", exact = TRUE), level)
      expect_equal(as.numeric(sample), as.numeric(expected), tolerance = 1e-14)
      prior <- attr(sample, "prior_density", exact = TRUE)
      expect_s3_class(prior, "prior_linear_density")
      if (manual_weight == 0) {
        expect_null(attached)
        expect_identical(BayesTools::prior_density_ordinate(prior, 0)$point_mass, 1)
      } else {
        expect_identical(names(attached$spec$weights), fixture$factor_columns[[1L]])
        expect_equal(unname(attached$spec$weights), manual_weight, tolerance = 1e-14)
        expect_identical(rendered$args$density_method, "precomputed")
        expect_equal(BayesTools::prior_density_ordinate(prior, .3)$log_density,
          stats::dnorm(.3, sd = .6 * abs(manual_weight), log = TRUE), tolerance = 1e-12)
      }
      error <- tryCatch(plot_prior(object, selector, component = "mods"), error = function(e) e)
      expect_identical(conditionMessage(error), paste0(
        "Individual factor-cell selection is unavailable for this method. ",
        "Use 'parameter = \"group\"' to select the whole factor term."))
      expect_null(conditionCall(error))
    }
    plot(object, "group", prior = TRUE, plot_type = "ggplot")
    expect_identical(rendered$renderer, "term")
    expect_identical(rendered$parameter, fixture$factor_name)
    expect_s3_class(rendered$samples[[fixture$factor_name]], "mixed_posteriors.factor")
  }
})
