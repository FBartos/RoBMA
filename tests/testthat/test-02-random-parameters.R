context("Semantic random-effect parameters")

source(testthat::test_path("common-functions.R"))

.random_parameter_fit_names <- c(
  "brma.mv_block_mvn_random",
  "brma.mv_block_mvn_random_scale",
  "brma.mv_block_mvn_known_R",
  "brma.mv_v14_konstantopoulos2011_cs",
  "brma.mv_v14_ishak2007_har",
  "brma.mv_v14_begg1989_study_treatment"
)

.active_random_parameter_fit_names <- function() {

  names <- intersect(
    .random_parameter_fit_names,
    active_fit_catalog()[["name"]]
  )
  if (length(names) == 0L) {
    testthat::skip("No semantic random-parameter fixtures are active.")
  }

  return(names)
}

.random_parameter_hypothesis <- function(label, operator_left, value,
                                         operator_right) {

  parameter <- paste0("`", label, "`")
  paste(
    parameter, operator_left, format(value, digits = 17),
    "vs", parameter, operator_right, format(value, digits = 17)
  )
}

.random_parameter_weights <- function(S, K) {

  positions <- seq_len(S)
  weights <- vapply(seq_len(K), function(i) {
    center <- 1L + ((17L * i) %% S)
    exp(-abs(positions - center) / max(S / 5, 1))
  }, numeric(S))
  sweep(weights, 2L, colSums(weights), "/")
}

test_that("random catalog matches summaries across structures", {

  fit_names <- .active_random_parameter_fit_names()
  skip_if_missing_fits(fit_names)

  for (name in fit_names) {
    fit      <- load_fit(name, validate = FALSE)
    bundle   <- .brma_random_parameter_bundle(fit)
    summary  <- summary(fit)
    quantity <- hypothesis_quantities(fit)
    random_quantity <- quantity[quantity[["component"]] == "random", , drop = FALSE]

    expect_true(nrow(bundle[["specs"]]) > 0L, info = name)
    expect_true(all(
      bundle[["specs"]][["label"]] %in%
        rownames(summary[["estimates_random"]])
    ), info = name)
    expect_setequal(
      bundle[["specs"]][["parameter"]],
      unique(random_quantity[["parameter"]])
    )
    fixed <- vapply(seq_len(ncol(bundle[["samples"]])), function(i) {
      diff(range(bundle[["samples"]][, i])) <= sqrt(.Machine$double.eps)
    }, logical(1))
    names(fixed) <- bundle[["specs"]][["parameter"]]
    expected_direction <- unname(!fixed[random_quantity[["parameter"]]])
    expect_equal(random_quantity[["direction_test"]], expected_direction,
                 info = name)
    expect_equal(
      random_quantity[["direction_test_methods"]][expected_direction],
      rep("KDE, normal", sum(expected_direction)),
      info = name
    )
  }
})

test_that("standard random fixture crosses every semantic consumer boundary", {

  skip_if_missing_fits("brma.mv_block_mvn_random")
  fit <- load_fit("brma.mv_block_mvn_random", validate = FALSE)
  bundle <- .brma_random_parameter_bundle(fit)
  parameter <- "sd(intercept)"
  index <- match(parameter, bundle[["specs"]][["label"]])

  expect_false(is.na(index))
  expect_identical(
    bundle[["specs"]][["display_transform"]][[index]],
    list(type = "identity")
  )
  expect_true(parameter %in% rownames(summary(fit)[["estimates_random"]]))
  expect_true(parameter %in% hypothesis_quantities(fit)[["alias"]])
  expect_identical(
    .brma_random_parameter_density_target(fit, parameter)[["parameter_spec"]][["type"]],
    "primitive"
  )
  expect_s3_class(
    plot(
      fit,
      parameter = parameter,
      component = "random",
      prior = TRUE,
      plot_type = "ggplot"
    ),
    "ggplot"
  )
  value <- stats::median(bundle[["samples"]][, index])
  expect_s3_class(
    hypothesis(
      fit,
      .random_parameter_hypothesis(parameter, ">", value, "<"),
      component = "random",
      density_method = "KDE",
      n_samples = 500L,
      seed = 29
    ),
    "BayesTools_hypothesis_BF"
  )
})

test_that("ordinary plots and hypotheses ignore unrelated random priors", {

  skip_if_missing_fits("brma.mv_block_mvn_random_scale")
  fit <- load_fit("brma.mv_block_mvn_random_scale", validate = FALSE)

  expect_s3_class(
    plot(fit, parameter = "mu", component = "mods", prior = TRUE,
         plot_type = "ggplot"),
    "ggplot"
  )
  expect_s3_class(
    plot(fit, parameter = "intercept", component = "scale", prior = TRUE,
         plot_type = "ggplot"),
    "ggplot"
  )
  expect_s3_class(
    hypothesis(fit, "mu > 0 vs mu < 0", component = "mods",
               n_samples = 1000, seed = 11),
    "BayesTools_hypothesis_BF"
  )
  expect_s3_class(
    hypothesis(fit, "intercept > 0.2 vs intercept < 0.2",
               component = "scale", n_samples = 1000, seed = 12),
    "BayesTools_hypothesis_BF"
  )
})

test_that("random plots and MCMC diagnostics use semantic draws", {

  fit_names <- .active_random_parameter_fit_names()
  skip_if_missing_fits(fit_names)

  for (name in fit_names) {
    fit   <- load_fit(name, validate = FALSE)
    specs <- .brma_random_parameter_bundle(fit)[["specs"]]
    label <- specs[["label"]][1L]

    expect_s3_class(
      plot(fit, parameter = label, component = "random", prior = TRUE,
           plot_type = "ggplot"),
      "ggplot"
    )
    expect_s3_class(
      plot_diagnostic_trace(
        fit, parameter = label, component = "random", plot_type = "ggplot"
      ),
      "ggplot"
    )
  }

  if ("brma.mv_v14_konstantopoulos2011_cs" %in% fit_names) {
    fit <- load_fit("brma.mv_v14_konstantopoulos2011_cs", validate = FALSE)
    for (type in c("density", "autocorrelation")) {
      expect_s3_class(
        plot_diagnostic(
          fit, parameter = "district: cor", component = "random",
          type = type, plot_type = "ggplot"
        ),
        "ggplot"
      )
    }
  }
})

test_that("random directional hypotheses use induced joint priors", {

  fit_names <- .active_random_parameter_fit_names()
  skip_if_missing_fits(fit_names)

  for (i in seq_along(fit_names)) {
    name     <- fit_names[[i]]
    fit      <- load_fit(name, validate = FALSE)
    selected <- .brma_random_parameter_bundle(fit)
    label    <- selected[["specs"]][["label"]][1L]
    value    <- stats::median(selected[["samples"]][, 1L])
    statement <- .random_parameter_hypothesis(label, ">", value, "<")

    out <- hypothesis(
      fit,
      statement,
      component      = "random",
      density_method = "KDE",
      n_samples      = 1000,
      seed           = 100 + i
    )
    expect_s3_class(out, "BayesTools_hypothesis_BF")
  }
})

test_that("random point hypotheses follow quantity-specific policy", {

  fit_names <- .active_random_parameter_fit_names()
  skip_if_missing_fits(fit_names)

  if ("brma.mv_v14_konstantopoulos2011_cs" %in% fit_names) {
    fit_cs <- load_fit("brma.mv_v14_konstantopoulos2011_cs", validate = FALSE)
    expect_s3_class(
      hypothesis(
        fit_cs,
        .random_parameter_hypothesis("district: cor", "=", 0.5, "!="),
        component = "random", n_samples = 2000, seed = 21
      ),
      "BayesTools_hypothesis_BF"
    )
  }

  if ("brma.mv_v14_ishak2007_har" %in% fit_names) {
    fit_har <- load_fit("brma.mv_v14_ishak2007_har", validate = FALSE)
    expect_error(
      hypothesis(
        fit_har,
        .random_parameter_hypothesis(
          "study: sd(time[1])", "=", 0.2, "!="
        ),
        component = "random", n_samples = 1000, seed = 22
      ),
      "derived component SD"
    )
  }

  if ("brma.mv_v14_begg1989_study_treatment" %in% fit_names) {
    fit_mixed <- load_fit(
      "brma.mv_v14_begg1989_study_treatment",
      validate = FALSE
    )
    mixed_quantities <- hypothesis_quantities(fit_mixed)
    mixed_rho <- mixed_quantities[
      mixed_quantities[["alias"]] == "treatment: cor" &
        mixed_quantities[["component"]] == "random",
      ,
      drop = FALSE
    ]
    expect_false(any(mixed_rho[["point_test"]]))
    expect_false(any(mixed_rho[["direction_test"]]))
    expect_true(all(grepl("fixed by the fitted model", mixed_rho[["reason"]],
                          fixed = TRUE)))
    expect_error(
      hypothesis(
        fit_mixed,
        .random_parameter_hypothesis("treatment: cor", "=", 0, "!="),
        component = "random", n_samples = 1000, seed = 23
      ),
      "contains a point mass"
    )
  }

  if ("brma.mv_block_mvn_random_scale" %in% fit_names) {
    fit_alloc <- load_fit("brma.mv_block_mvn_random_scale", validate = FALSE)
    expect_error(
      hypothesis(
        fit_alloc,
        .random_parameter_hypothesis(
          "var_prop(study)", "=", 0, ">"
        ),
        component = "random", n_samples = 1000, seed = 24
      ),
      "support boundary"
    )
    expect_s3_class(
      hypothesis(
        fit_alloc,
        .random_parameter_hypothesis(
          "var_prop(study)", ">", 0.5, "<"
        ),
        component = "random", density_method = "qCMDE",
        n_samples = 1000, seed = 25
      ),
      "BayesTools_hypothesis_BF"
    )
  }
})

test_that("random influence matches weighted scalar moment oracles", {

  skip_if_missing_fits("brma.mv_block_mvn_random_scale")
  fit     <- load_fit("brma.mv_block_mvn_random_scale", validate = FALSE)
  bundle  <- .brma_random_parameter_bundle(fit)
  samples <- bundle[["samples"]]
  weights <- .random_parameter_weights(nrow(samples), nobs(fit))

  observed <- dfbetas(
    fit,
    component = "random",
    .weights  = weights
  )
  loo_mean <- crossprod(weights, samples)
  expected <- matrix(NA_real_, nrow = nobs(fit), ncol = ncol(samples))
  for (j in seq_len(ncol(samples))) {
    centered <- outer(samples[, j], loo_mean[, j], "-")
    loo_sd   <- sqrt(colSums(weights * centered^2))
    expected[, j] <- (mean(samples[, j]) - loo_mean[, j]) / loo_sd
  }
  colnames(expected) <- colnames(samples)
  expect_equal(unname(as.matrix(observed)), unname(expected), tolerance = 1e-12)
  expect_equal(colnames(observed), bundle[["specs"]][["parameter"]])

  parameter <- bundle[["specs"]][["label"]][1L]
  x         <- samples[, 1L]
  full_var  <- mean((x - mean(x))^2)
  loo_var <- vapply(seq_len(ncol(weights)), function(i) {
    loo_mean_i <- sum(weights[, i] * x)
    sum(weights[, i] * (x - loo_mean_i)^2)
  }, numeric(1))
  observed_covratio <- covratio(
    fit,
    component = "random",
    parameter = parameter,
    .weights  = weights
  )
  expect_equal(unname(observed_covratio), loo_var / full_var,
               tolerance = 1e-12)
  expect_error(
    covratio(fit, component = "random", .weights = weights),
    "requires one explicit 'parameter'"
  )
  expect_error(
    dfbetas(
      fit, type = "scale", component = "random", .weights = weights
    ),
    "select different parameter namespaces"
  )
})

test_that("fixed and unavailable random influence targets are explicit", {

  skip_if_missing_fits("brma.mv_v14_begg1989_study_treatment")
  fit     <- load_fit("brma.mv_v14_begg1989_study_treatment", validate = FALSE)
  bundle  <- .brma_random_parameter_bundle(fit)
  weights <- .random_parameter_weights(nrow(bundle[["samples"]]), nobs(fit))

  dfb <- dfbetas(
    fit,
    component = "random",
    parameter = "treatment: cor",
    .weights  = weights
  )
  covr <- covratio(
    fit,
    component = "random",
    parameter = "treatment: cor",
    .weights  = weights
  )
  expect_true(all(is.nan(as.matrix(dfb))))
  expect_match(attr(dfb, "note"), "LOO posterior variance is zero")
  expect_true(all(is.nan(covr)))
  expect_match(attr(covr, "note"), "no parameters with non-zero posterior variance")
  expect_error(
    dfbetas(
      fit, component = "random", parameter = "latent_z", .weights = weights
    ),
    "not available|Could not select|No public parameter quantity"
  )
})

test_that("no-intercept random formulas expose only fitted location terms", {

  skip_if_missing_fits("brma.mv_v14_ishak2007_har")
  fit     <- load_fit("brma.mv_v14_ishak2007_har", validate = FALSE)
  catalog <- .brma_parameter_catalog(fit)
  mods    <- catalog[catalog[["component"]] == "mods", , drop = FALSE]

  expect_false(any(mods[["alias"]] %in% c("mu", "mu_intercept", "intercept")))
  expect_identical(.brma_parameter_default_formula(fit, "mods"), "mu_time_factor")
  expect_s3_class(
    plot(fit, component = "mods", plot_type = "ggplot"),
    "ggplot"
  )
  expect_error(
    hypothesis(fit, "mu > 0 vs mu < 0", component = "mods"),
    "No public parameter quantity matches 'mu'"
  )
})


test_that("IWMDE disables focal prior delta for sampled random SD rows", {

  skip_if_missing_fits("brma.mv_block_mvn_random_scale")

  context <- .iwmde_context(load_fit(
    "brma.mv_block_mvn_random_scale",
    validate = FALSE
  ))
  parameter <- "log_tau_intercept"
  if (!parameter %in% colnames(context[["posterior_samples"]])) {
    skip("brma.mv random-scale fixture does not contain log_tau_intercept.")
  }
  rows <- which(is.finite(context[["posterior_samples"]][, parameter]))
  rows <- head(rows, 3L)

  expect_gt(length(rows), 0L)
  for (row in rows) {
    state <- .iwmde_row_state(context, row, parameter)
    expect_false(state[["use_focal_prior_delta"]])
  }
})


test_that("nonlinear transformed scale intercepts fail qCMDE and IWMDE closed", {

  skip_if_missing_fits("brma.mv_block_mvn_random_scale")
  fit <- load_fit("brma.mv_block_mvn_random_scale", validate = FALSE)
  entry <- .brma_parameter_select_entry(
    fit,
    parameter = "intercept",
    component = "scale"
  )
  transform <- BayesTools::JAGS_formula_coefficient_transform(
    fit[["fit"]],
    parameter = entry[["formula_parameter"]]
  )
  target_i <- match(entry[["parameter"]], transform[["target_names"]])
  dependencies <- transform[["dependencies"]][
    transform[["dependencies"]][["target"]] == entry[["parameter"]],
    ,
    drop = FALSE
  ]
  nonlinear_joint <- nrow(dependencies) > 1L &&
    (transform[["output_transforms"]][[target_i]] != "identity" ||
       any(transform[["source_transforms"]][dependencies[["source"]]] !=
             "identity"))
  expect_true(nonlinear_joint)

  kde <- hypothesis(
    fit,
    "intercept = 0.2 vs intercept != 0.2",
    component      = "scale",
    density_method = "KDE",
    n_samples = 1000L,
    seed      = 31
  )

  expect_s3_class(kde, "BayesTools_hypothesis_BF")
  for (method in c("qCMDE", "IWMDE")) {
    expect_error(
      hypothesis(
        fit,
        "intercept = 0.2 vs intercept != 0.2",
        component       = "scale",
        density_method  = method,
        density_control = list(
          n_points             = 30L,
          samples              = 40L,
          normalization_points = 40L
        ),
        n_samples = 1000L,
        seed      = 32
      ),
      "supported only with density_method = 'KDE'"
    )
  }
})

test_that("random prior overlays and diagnostic labels are semantic", {

  fit_names <- intersect(c(
    "brma.mv_v14_begg1989_study_treatment",
    "brma.mv_v14_konstantopoulos2011_cs"
  ), active_fit_catalog()[["name"]])
  if (length(fit_names) == 0L) {
    testthat::skip("No random-prior overlay fixtures are active.")
  }
  skip_if_missing_fits(fit_names)

  if ("brma.mv_v14_begg1989_study_treatment" %in% fit_names) {
    fixed <- load_fit(
      "brma.mv_v14_begg1989_study_treatment",
      validate = FALSE
    )
    expect_s3_class(
      plot(
        fixed,
        parameter = "treatment: cor", component = "random", prior = TRUE,
        plot_type = "ggplot"
      ),
      "ggplot"
    )
  }

  if ("brma.mv_v14_konstantopoulos2011_cs" %in% fit_names) {
    cs <- load_fit("brma.mv_v14_konstantopoulos2011_cs", validate = FALSE)
    diagnostic <- plot_diagnostic_density(
      cs,
      parameter = "district: cor", component = "random",
      plot_type = "ggplot"
    )
    expect_identical(diagnostic[["labels"]][["title"]], "district: cor")
  }
})

test_that("Dirichlet allocation priors use their exact beta marginals", {

  selected <- list(
    spec         = list(quantity = "var_prop", allocation_index = 2L),
    source_prior = BayesTools::prior(
      "dirichlet",
      parameters = list(alpha = c(2, 3, 5))
    ),
    prior        = NULL
  )
  prior <- .brma_random_parameter_exact_prior(selected)

  expect_true(BayesTools::is.prior(prior))
  expect_equal(
    BayesTools::lpdf(prior, c(0.2, 0.4, 0.7)),
    stats::dbeta(c(0.2, 0.4, 0.7), 3, 7, log = TRUE)
  )
})

test_that("simplex density replacements preserve auxiliary-gamma coordinates", {

  source            <- "mu_allocation"
  columns           <- paste0(source, "[", 1:2, "]")
  auxiliary_columns <- .iwmde_simplex_auxiliary_columns(source, 2L)
  samples <- matrix(
    c(
      0.25, 0.75, 1, 3,
      0.50, 0.50, 2, 2
    ),
    nrow = 2L,
    byrow = TRUE,
    dimnames = list(NULL, c(columns, auxiliary_columns))
  )
  context <- list(posterior_samples = samples)
  spec <- .iwmde_parameter_spec(
    context,
    columns[[1L]],
    list(type = "simplex_pair", parameter = source, index = 1L)
  )
  inconsistent <- samples
  inconsistent[, columns[[1L]]] <- 0.4
  inconsistent[, columns[[2L]]] <- 0.6
  expect_identical(
    .iwmde_parameter_spec(
      list(posterior_samples = inconsistent),
      columns[[1L]],
      list(type = "simplex_pair", parameter = source, index = 1L)
    )[["status"]],
    "ok"
  )
  outside_support <- samples
  outside_support[, columns[[1L]]] <- 1.1
  expect_identical(
    .iwmde_parameter_spec(
      list(posterior_samples = outside_support),
      columns[[1L]],
      list(type = "simplex_pair", parameter = source, index = 1L)
    )[["status"]],
    "unsupported"
  )
  replacement <- .iwmde_replacement_spec(context, columns[[1L]], spec)
  replaced <- .iwmde_replace_row_for_value(
    context,
    list(row = samples[1L, ]),
    columns[[1L]],
    0.4,
    replacement
  )

  expect_equal(unname(replaced[["row"]][columns]), c(0.4, 0.6))
  expect_equal(unname(replaced[["row"]][auxiliary_columns]), c(1.6, 2.4))

  prior <- BayesTools::prior(
    "dirichlet",
    parameters = list(alpha = c(2, 3))
  )
  expected <- stats::dgamma(1.6, 2, 1, log = TRUE) +
    stats::dgamma(2.4, 3, 1, log = TRUE)
  expect_equal(
    .iwmde_log_prior_row(replaced[["row"]], stats::setNames(list(prior), source)),
    expected
  )
  prior_list <- stats::setNames(list(prior), source)
  expect_equal(
    .iwmde_replacement_log_prior(
      parameter       = columns[[1L]],
      values          = 0.4,
      valid_samples   = samples,
      valid_positions = 1:2,
      candidates      = list(
        state_index = 1:2,
        grid_index  = c(1L, 1L)
      ),
      row_states = list(
        list(use_focal_prior_delta = FALSE, prior_list = prior_list),
        list(use_focal_prior_delta = FALSE, prior_list = prior_list)
      ),
      replacement = list(type = "simplex_pair")
    ),
    c(
      stats::dgamma(1, 2, 1, log = TRUE) +
        stats::dgamma(3, 3, 1, log = TRUE),
      stats::dgamma(2, 2, 1, log = TRUE) +
        stats::dgamma(2, 3, 1, log = TRUE)
    )
  )
})

test_that("allocation replacements synchronize derived random-effect SDs", {

  dat <- data.frame(
    yi    = c(.10, .20, .30, .40),
    study = c("a", "a", "b", "b"),
    esid  = c("a1", "a2", "b1", "b2")
  )
  object <- brma.mv(
    yi                        = yi,
    V                         = diag(rep(.04, 4L)),
    random                    = ~ 1 | study / esid,
    data                      = dat,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  design     <- .fitted_formula_design(object, "mu")
  allocation <- design[["random_allocations"]][["heterogeneity"]]
  source     <- allocation[["source"]][["name"]]
  weights    <- paste0(allocation[["weight_name"]], "[", 1:2, "]")
  components <- vapply(
    design[["random_effects"]],
    function(term) term[["sd_parameter_names"]][[1L]],
    character(1)
  )
  samples <- matrix(
    c(.8, .25, .75, 9, 9),
    nrow = 1L,
    dimnames = list(NULL, c(source, weights, components))
  )
  context <- list(object = object, data = object[["data"]])

  synced <- .iwmde_sync_random_allocation_sd_matrix(
    context    = context,
    samples    = samples,
    parameters = source
  )
  expect_true(synced[["valid"]])
  expect_equal(
    unname(synced[["samples"]][1L, components]),
    .8 * sqrt(c(.25, .75))
  )

  samples[1L, weights] <- c(.4, .6)
  synced <- .iwmde_sync_random_allocation_sd_matrix(
    context    = context,
    samples    = samples,
    parameters = weights
  )
  expect_equal(
    unname(synced[["samples"]][1L, components]),
    .8 * sqrt(c(.4, .6))
  )
})

test_that("allocated component SD targets retain their replacement structure", {

  source  <- "mu_total_sd"
  weight  <- "mu_allocation"
  columns <- c(
    source,
    paste0(weight, "[", 1:2, "]"),
    .iwmde_simplex_auxiliary_columns(weight, 2L)
  )
  samples <- matrix(
    c(0.8, 0.25, 0.75, 1, 3),
    nrow = 1L,
    dimnames = list(NULL, columns)
  )
  context <- list(
    posterior_samples = samples,
    indicator_names   = character(),
    flat_prior_list   = stats::setNames(
      list(BayesTools::prior(
        "dirichlet",
        parameters = list(alpha = c(1, 1))
      )),
      weight
    ),
    selection_spec = NULL
  )
  input_spec <- list(
    type                 = "random_component_sd",
    source_parameter     = source,
    factors              = list(list(
      weight_name = weight,
      index       = 1L,
      scale       = "mean_variance",
      n_targets   = 2L
    )),
    target_columns       = source,
    conditioning_exclude = paste0(weight, "[2]")
  )
  spec <- .iwmde_parameter_spec(context, source, input_spec)
  replacement <- .iwmde_replacement_spec(context, source, spec)
  replaced <- .iwmde_replace_row_for_value(
    context,
    list(row = samples[1L, ]),
    source,
    0.7,
    replacement
  )

  expect_equal(
    .iwmde_parameter_values(context, source, spec),
    0.8 * sqrt(2 * 0.25)
  )
  expect_equal(replaced[["row"]][[source]], 0.7 / sqrt(2 * 0.25))
  expect_equal(
    .iwmde_chen_conditioning_columns(context, source, spec),
    paste0(weight, "[1]")
  )
  expect_named(
    .iwmde_plan_parameter_spec(spec),
    c(
      "type", "source_parameter", "factors", "target_columns",
      "factor_columns", "auxiliary_columns", "conditioning_exclude", "status"
    )
  )
})

test_that("direct multivariate random quantities expose density targets", {

  skip_if_missing_fits("brma.mv_v14_assink2016_nested")
  fit <- load_fit("brma.mv_v14_assink2016_nested", validate = FALSE)

  total <- .brma_random_parameter_density_target(
    fit,
    "sd_total"
  )
  nested_sd <- .brma_random_parameter_density_target(
    fit,
    "esid_study: sd(intercept)"
  )
  study_sd <- .brma_random_parameter_density_target(
    fit,
    "study: sd(intercept)"
  )
  nested <- .brma_random_parameter_density_target(
    fit,
    "var_prop(esid_study)"
  )
  study <- .brma_random_parameter_density_target(
    fit,
    "var_prop(study)"
  )

  expect_identical(total[["parameter_spec"]][["type"]], "primitive")
  expect_identical(
    nested_sd[["parameter_spec"]][["type"]],
    "random_component_sd"
  )
  expect_identical(nested_sd[["parameter_spec"]][["factors"]][[1L]][["index"]], 1L)
  expect_identical(study_sd[["parameter_spec"]][["factors"]][[1L]][["index"]], 2L)
  expect_identical(nested[["parameter_spec"]][["type"]], "simplex_pair")
  expect_identical(nested[["parameter_spec"]][["index"]], 1L)
  expect_identical(study[["parameter_spec"]][["index"]], 2L)

  total_plot <- plot(
    fit,
    "sd_total",
    component = "random",
    plot_type = "ggplot"
  )
  expect_identical(
    total_plot$scales$get_scales("x")$name,
    "sd_total"
  )

  samples <- .brma_random_parameter_mixed_posterior(
    fit,
    "var_prop(esid_study)",
    prior = TRUE
  )
  parameter <- names(samples)[[1L]]
  prior     <- attr(samples, "prior_list", exact = TRUE)[[parameter]]
  expect_equal(BayesTools::lpdf(prior, c(0.2, 0.5, 0.8)), rep(0, 3L))
  expect_null(attr(samples[[parameter]], "prior_density", exact = TRUE))

  expect_s3_class(
    plot(
      fit,
      "study: sd(intercept)",
      component       = "random",
      density_method  = "qCMDE",
      density_control = list(
        n_points             = 20L,
        normalization_points = 200L
      ),
      plot_type = "ggplot"
    ),
    "ggplot"
  )

  context <- .iwmde_context(fit)
  columns <- paste0(nested[["parameter_spec"]][["parameter"]], "[", 1:2, "]")
  state <- .iwmde_row_state(
    context,
    1L,
    nested[["parameter"]],
    nested[["parameter_spec"]]
  )
  replacement <- .iwmde_replacement_spec(
    context,
    nested[["parameter"]],
    nested[["parameter_spec"]]
  )
  replaced <- .iwmde_replace_row_for_value(
    context,
    state,
    nested[["parameter"]],
    0.37,
    replacement
  )
  expect_equal(unname(unlist(replaced[["row"]][columns])), c(0.37, 0.63))
  eta_sum <- sum(state[["row"]][nested[["parameter_spec"]][["auxiliary_columns"]]])
  expect_equal(
    unname(unlist(replaced[["row"]][nested[["parameter_spec"]][["auxiliary_columns"]]])),
    eta_sum * c(0.37, 0.63)
  )
})

test_that("qCMDE plots a two-component multivariate allocation proportion", {

  skip_if_missing_fits("brma.mv_v14_assink2016_nested")
  fit <- load_fit("brma.mv_v14_assink2016_nested", validate = FALSE)

  expect_s3_class(
    plot(
      fit,
      "var_prop(esid_study)",
      component       = "random",
      prior           = TRUE,
      density_method  = "qCMDE",
      density_control = list(
        n_points             = 20L,
        samples              = 300L,
        normalization_points = 60L
      ),
      plot_type = "ggplot"
    ),
    "ggplot"
  )
})

test_that("random DFBETAS zero-variance handling is cellwise", {

  skip_if_missing_fits("brma.mv_block_mvn_random_scale")
  fit     <- load_fit("brma.mv_block_mvn_random_scale", validate = FALSE)
  samples <- .brma_random_parameter_bundle(fit)[["samples"]]
  weights <- .random_parameter_weights(nrow(samples), nobs(fit))
  weights[, 1L] <- 0
  weights[1L, 1L] <- 1

  observed <- suppressWarnings(dfbetas(
    fit,
    component = "random",
    .weights  = weights
  ))
  expect_true(all(is.nan(as.matrix(observed)[1L, ])))
  expect_true(all(is.finite(as.matrix(observed)[-1L, ])))
  expect_match(attr(observed, "note"), "LOO posterior variance is zero")
})

test_that("bounded induced-prior densities are normalized and corrected", {

  samples <- seq(0.0005, 0.9995, length.out = 2000L)
  density <- .brma_random_parameter_prior_density(
    samples = samples,
    support = c(0, 1),
    n_points = 2048L
  )[["density"]]
  integral <- sum(diff(density[["x"]]) *
    (density[["y"]][-1L] + density[["y"]][-length(density[["y"]])]) / 2)
  center <- which.min(abs(density[["x"]] - 0.5))

  expect_equal(integral, 1, tolerance = 1e-12)
  expect_gt(density[["y"]][1L] / density[["y"]][center], 0.8)
  expect_gt(density[["y"]][length(density[["y"]])] /
    density[["y"]][center], 0.8)
})
