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

# The label of the fit's one random-effect quantity of catalog type
# 'quantity' (e.g. "cor" -> "rho" or "treatment: rho").
.random_parameter_catalog_label <- function(fit, quantity) {

  specs <- .brma_random_parameter_bundle(fit)[["specs"]]
  label <- specs[["label"]][specs[["quantity"]] == quantity]
  if (length(label) != 1L) {
    stop(
      "The fit has ", length(label), " random-effect catalog quantities of ",
      "type '", quantity, "'; expected one.",
      call. = FALSE
    )
  }

  return(label)
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
  skip_if_not(length(fit_names) > 0L, "No random-parameter fixtures are active.")
  skip_if_missing_fits(fit_names)

  for (name in fit_names) {
    fit      <- load_fit(name, validate = FALSE)
    bundle   <- .brma_random_parameter_bundle(fit)
    summary  <- summary(fit)
    quantity <- hypothesis_quantities(fit)
    random_quantity <- quantity[quantity[["component"]] == "random", , drop = FALSE]

    expect_true(nrow(bundle[["specs"]]) > 0L, info = name)
    summary_labels <- rownames(summary[["estimates_random"]])
    summary_parameters <- attr(summary[["estimates_random"]], "parameters")
    matched <- summary_parameters %in% bundle[["specs"]][["parameter"]]
    expect_true(any(matched), info = name)
    expect_equal(
      summary_labels[matched],
      bundle[["specs"]][["label"]][match(
        summary_parameters[matched],
        bundle[["specs"]][["parameter"]]
      )],
      info = name
    )
    expect_setequal(
      bundle[["specs"]][["parameter"]],
      unique(random_quantity[["parameter"]])
    )
    # Region tests are available for every quantity that the catalog does
    # not declare structurally fixed; they use no density method.
    fixed <- bundle[["specs"]][["status"]] == "structural"
    names(fixed) <- bundle[["specs"]][["parameter"]]
    expected_direction <- unname(!fixed[random_quantity[["parameter"]]])
    expect_equal(random_quantity[["direction_test"]], expected_direction,
                 info = name)
    expect_false("direction_test_methods" %in% names(random_quantity))
  }
})

test_that("standard random fixture crosses every semantic consumer boundary", {

  skip_if_missing_fits("brma.mv_block_mvn_random")
  fit <- load_fit("brma.mv_block_mvn_random", validate = FALSE)
  bundle <- .brma_random_parameter_bundle(fit)
  parameter <- "tau"
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
          fit, parameter = "rho", component = "random",
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

  fit_names <- intersect(
    .active_random_parameter_fit_names(),
    c(
      "brma.mv_v14_konstantopoulos2011_cs",
      "brma.mv_v14_ishak2007_har",
      "brma.mv_v14_begg1989_study_treatment",
      "brma.mv_block_mvn_random_scale"
    )
  )
  skip_if_not(length(fit_names) > 0L, "No random-parameter fixtures are active.")
  skip_if_missing_fits(fit_names)

  if ("brma.mv_v14_konstantopoulos2011_cs" %in% fit_names) {
    fit_cs <- load_fit("brma.mv_v14_konstantopoulos2011_cs", validate = FALSE)
    expect_s3_class(
      hypothesis(
        fit_cs,
        .random_parameter_hypothesis("rho", "=", 0.5, "!="),
        component = "random", n_samples = 2000, seed = 21
      ),
      "BayesTools_hypothesis_BF"
    )
  }

  if ("brma.mv_v14_ishak2007_har" %in% fit_names) {
    fit_har <- load_fit("brma.mv_v14_ishak2007_har", validate = FALSE)
    # The derived component SD tau(time[1]) = allocation SD * sqrt(4 * weight)
    # has an exact prior ordinate, so its point hypothesis is eligible. The
    # null 0.2 lies far below every posterior draw (all above 1.4): qCMDE
    # returns an extreme Bayes factor whose relative precision target is
    # missed while its conclusion is certain, and the attached diagnostics
    # report the collapsed importance weights.
    har <- hypothesis(
      fit_har,
      .random_parameter_hypothesis(
        "tau(time[1])", "=", 0.2, "!="
      ),
      component = "random", n_samples = 1000, seed = 22
    )
    expect_s3_class(har, "BayesTools_hypothesis_BF")
    expect_true(is.finite(har[["BF"]]) && har[["BF"]] < 1e-10)
    har_diagnostics <- density_diagnostics(har)
    expect_identical(nrow(har_diagnostics), 1L)
    expect_false(har_diagnostics[["precision_target_met"]])
    expect_true(har_diagnostics[["bf_grade_met"]])
    expect_lt(har_diagnostics[["ess"]], har_diagnostics[["warning_min_ess"]])
    expect_gt(
      har_diagnostics[["max_weight_share"]],
      har_diagnostics[["warning_max_weight_share"]]
    )
    expect_true(nzchar(har_diagnostics[["warnings"]]))
  }

  if ("brma.mv_v14_begg1989_study_treatment" %in% fit_names) {
    fit_mixed <- load_fit(
      "brma.mv_v14_begg1989_study_treatment",
      validate = FALSE
    )
    mixed_label      <- .random_parameter_catalog_label(fit_mixed, "cor")
    mixed_quantities <- hypothesis_quantities(fit_mixed)
    mixed_rho <- mixed_quantities[
      mixed_quantities[["alias"]] == mixed_label &
        mixed_quantities[["component"]] == "random",
      ,
      drop = FALSE
    ]
    expect_identical(nrow(mixed_rho), 1L)
    expect_false(any(mixed_rho[["point_test"]]))
    expect_false(any(mixed_rho[["direction_test"]]))
    expect_true(all(grepl("fixed by the fitted model", mixed_rho[["reason"]],
                          fixed = TRUE)))
    expect_error(
      hypothesis(
        fit_mixed,
        .random_parameter_hypothesis(mixed_label, "=", 0, "!="),
        component = "random", n_samples = 1000, seed = 23
      ),
      class = "RoBMA_hypothesis_fixed"
    )
  }

  if ("brma.mv_block_mvn_random_scale" %in% fit_names) {
    fit_alloc <- load_fit("brma.mv_block_mvn_random_scale", validate = FALSE)
    # The variance proportion at its lower bound: the one-sided Dirichlet
    # share ordinate is exact, finite, and positive, so the boundary point
    # null has a Savage-Dickey Bayes factor with one-sided densities.
    boundary <- hypothesis(
      fit_alloc,
      .random_parameter_hypothesis(
        "tau2_prop(study)", "=", 0, ">"
      ),
      component = "random", n_samples = 1000, seed = 24
    )
    expect_true(is.finite(attr(boundary, "raw_BF")))
    # The qCMDE ordinate at the bound is continuous with the ordinates inside
    # the support: its jump from the ordinate at 0.005 is no larger than the
    # change over the next 0.005 plus four combined Monte Carlo standard
    # errors (BF_error is the relative error of the ordinate).
    ordinates <- lapply(c(0, 0.005, 0.01), function(value) {
      out <- hypothesis(
        fit_alloc,
        paste0("`tau2_prop(study)` = ", value),
        component = "random", n_samples = 1000, seed = 24, columns = "all"
      )
      c(ordinate = as.numeric(out[["posterior"]]),
        mcse     = as.numeric(out[["posterior"]]) * as.numeric(out[["BF_error"]]) / 100)
    })
    ordinate <- vapply(ordinates, `[[`, numeric(1), "ordinate")
    mcse     <- vapply(ordinates, `[[`, numeric(1), "mcse")
    expect_true(all(is.finite(ordinate) & ordinate > 0))
    expect_lte(
      abs(ordinate[[1L]] - ordinate[[2L]]),
      abs(ordinate[[3L]] - ordinate[[2L]]) + 4 * sqrt(mcse[[1L]]^2 + mcse[[2L]]^2)
    )
    expect_s3_class(
      hypothesis(
        fit_alloc,
        .random_parameter_hypothesis(
          "tau2_prop(study)", ">", 0.5, "<"
        ),
        component = "random", density_method = "qCMDE",
        n_samples = 1000, seed = 25
      ),
      "BayesTools_hypothesis_BF"
    )
  }
})

test_that("random-effect variances are tested through their standard deviations", {

  skip_if_missing_fits("BMA.mv_random_components")
  fit <- load_fit("BMA.mv_random_components", validate = FALSE)
  # 'tau2 = v^2' is 'tau = v' with the prior and posterior densities divided
  # by the derivative 2 * v of the square map, so every density method gives
  # the same Bayes factor for the gated SD and its variance.
  value <- 0.06
  for (method in c("KDE", "qCMDE", "IWMDE")) {
    control <- if (!identical(method, "KDE")) list(n_points = 20L, samples = 50L)
    run <- function(statement) {
      suppressWarnings(hypothesis(
        fit, statement, density_method = method, density_control = control,
        columns = "all", seed = 1
      ))
    }
    sd  <- run(paste0("`(mu) study: tau(intercept)` = ", value))
    var <- run(paste0("`(mu) study: tau2(intercept)` = ", value^2))
    expect_true(is.finite(attr(sd, "raw_BF")), info = method)
    expect_equal(attr(var, "raw_BF"), attr(sd, "raw_BF"), tolerance = 1e-10,
                 info = method)
    expect_equal(as.numeric(var[["prior"]]),
                 as.numeric(sd[["prior"]]) / (2 * value),
                 tolerance = 1e-10, info = method)
    expect_equal(as.numeric(var[["posterior"]]),
                 as.numeric(sd[["posterior"]]) / (2 * value),
                 tolerance = 1e-10, info = method)
    # A variance point against a region: the point part (the point statement
    # against its complement, through the SD) over the region part (the
    # region against the encompassing model, on the variance draws), with
    # their errors combined on the log scale.
    tau2   <- "`(mu) study: tau2(intercept)`"
    point  <- run(paste0(tau2, " = ", value^2, " vs ", tau2, " != ", value^2))
    region <- run(paste0(tau2, " > ", value^2, " vs ", tau2, " >= 0"))
    mixed  <- run(paste0(tau2, " = ", value^2, " vs ", tau2, " > ", value^2))
    expect_equal(attr(mixed, "raw_BF"), attr(point, "raw_BF") / attr(region, "raw_BF"),
                 tolerance = 1e-10, info = method)
    expect_equal(
      as.numeric(mixed[["BF_error"]]),
      sqrt(as.numeric(point[["BF_error"]])^2 + as.numeric(region[["BF_error"]])^2),
      tolerance = 1e-10, info = method
    )
  }
  # Draws conditional on the inclusion of the gated SD (KDE, the method of
  # conditional random-effect hypotheses): the conditional point part differs
  # from the averaged one by the inclusion Bayes factor, and the identity holds
  # on the conditional draws.
  run_conditional <- function(statement) {
    suppressWarnings(hypothesis(
      fit, statement, density_method = "KDE", conditional = TRUE,
      columns = "all", seed = 1
    ))
  }
  point  <- run_conditional(paste0(tau2, " = ", value^2, " vs ", tau2, " != ", value^2))
  region <- run_conditional(paste0(tau2, " > ", value^2, " vs ", tau2, " >= 0"))
  mixed  <- run_conditional(paste0(tau2, " = ", value^2, " vs ", tau2, " > ", value^2))
  expect_equal(attr(mixed, "raw_BF"), attr(point, "raw_BF") / attr(region, "raw_BF"),
               tolerance = 1e-10)
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
  label   <- .random_parameter_catalog_label(fit, "cor")

  dfb <- dfbetas(
    fit,
    component = "random",
    parameter = label,
    .weights  = weights
  )
  covr <- covratio(
    fit,
    component = "random",
    parameter = label,
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


test_that("nonlinear transformed scale intercepts test points with KDE and fail qCMDE and IWMDE closed", {

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

  # The exp(affine) target admits KDE only; qCMDE/IWMDE are refused by the
  # plan's method refusal.
  plan <- .hypothesis_plans(
    fit, "intercept = 0.2 vs intercept != 0.2", component = "scale"
  )[[1L]]
  expect_identical(plan[["kind"]], "exp_affine")
  for (method in c("qCMDE", "IWMDE", "normal")) {
    expect_match(
      plan[["method_refusals"]][[method]][["reason"]],
      "supported only with density_method = 'KDE'",
      fixed = TRUE,
      info  = method
    )
  }
  # Its certified prior density is exact: the original-scale intercept is the
  # fitted intercept (a log source) times the exp of the Gaussian slope part,
  # a scale product with an exact ordinate. Point hypotheses therefore run
  # with KDE, and qCMDE/IWMDE stop with the plan's method refusal.
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
      "supported only with density_method = 'KDE'",
      fixed = TRUE,
      class = "RoBMA_hypothesis_method",
      info  = method
    )
  }
  kde <- hypothesis(
    fit,
    "intercept = 0.2 vs intercept != 0.2",
    component      = "scale",
    density_method = "KDE",
    n_samples      = 1000L,
    seed           = 32,
    columns        = "all"
  )
  # Independent Savage-Dickey: the prior ordinate by numerical integration of
  # the scale product over the slope prior, and the posterior ordinate as the
  # Gaussian kernel sum of the fitted draws of the intercept reflected at 0
  # (helper-contracts.R); the transform's weight on the slope is the fitted
  # standardization shift -mean(x) / sd(x).
  reference <- exp_affine_scale_intercept_reference(fit, slope = "x", value = 0.2)
  expect_equal(unname(plan[["weights"]][["log_tau_x"]]), -reference[["shift"]],
               tolerance = 1e-12)
  expect_equal(as.numeric(plan[["draws"]][["posterior"]]), reference[["draws"]],
               tolerance = 1e-12)
  expect_equal(as.numeric(kde[["prior"]]), reference[["prior"]], tolerance = 1e-10)
  expect_equal(as.numeric(kde[["posterior"]]), reference[["posterior"]], tolerance = 1e-12)
  expect_equal(attr(kde, "raw_BF"), reference[["posterior"]] / reference[["prior"]],
               tolerance = 1e-10)
  expect_s3_class(
    hypothesis(
      fit,
      "intercept > 0.2",
      component      = "scale",
      density_method = "KDE",
      n_samples      = 1000L,
      seed           = 31
    ),
    "BayesTools_hypothesis_BF"
  )
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
        parameter = .random_parameter_catalog_label(fixed, "cor"),
        component = "random", prior = TRUE, plot_type = "ggplot"
      ),
      "ggplot"
    )
  }

  if ("brma.mv_v14_konstantopoulos2011_cs" %in% fit_names) {
    cs <- load_fit("brma.mv_v14_konstantopoulos2011_cs", validate = FALSE)
    diagnostic <- plot_diagnostic_density(
      cs,
      parameter = "rho", component = "random",
      plot_type = "ggplot"
    )
    expect_identical(diagnostic[["labels"]][["title"]], "rho")
  }
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
      replacement = list(type = "simplex_pair", index = 1L)
    ),
    stats::dbeta(c(.25, .5), 2, 3, log = TRUE)
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
  # The SDs are recomputed by the BayesTools nodes of the fit: the prior-only
  # object gets a fit with its BayesTools formula and these draws.
  formula <- .object_bayestools_formula(
    object    = object,
    parameter = "mu",
    source    = .fitted_formula_source(parameter = "mu", data = object[["data"]])
  )
  fit <- coda::mcmc.list(coda::mcmc(samples))
  attr(fit, "prior_list")     <- formula[["prior_list"]]
  attr(fit, "formula_design") <- list(mu = formula[["formula_design"]])
  object[["fit"]] <- as_bayestools_fit(fit)
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

  object <- sd_component_allocation_object(list(
    source = 0.8,
    weight = cbind(0.25, 0.75),
    eta    = cbind(1, 3)
  ))
  samples    <- as.matrix(object[["fit"]][[1L]])
  term       <- attr(object[["fit"]], "formula_design")[["mu"]][["random_effects"]][[1L]]
  allocation <- term[["sd_binding"]][["allocations"]][[1L]]
  source     <- allocation[["source"]][["name"]]
  weight     <- allocation[["weight_name"]]
  context <- list(
    object            = object,
    posterior_samples = samples,
    indicator_names   = character(),
    flat_prior_list   = attr(object[["fit"]], "prior_list"),
    selection_spec    = NULL,
    evaluator_cache   = new.env(parent = emptyenv())
  )
  input_spec <- list(
    type                 = "random_component_sd",
    source_parameter     = source,
    node                 = "mu__xREx__g_intercept",
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
  # The target's source and the excluded second weight are not conditioned on.
  expect_equal(
    .iwmde_chen_conditioning_columns(context, source, spec),
    c("mu_intercept", paste0(weight, "[1]"))
  )
  expect_named(
    .iwmde_plan_parameter_spec(spec),
    c(
      "type", "source_parameter", "node", "factors", "target_columns",
      "factor_columns", "auxiliary_columns", "conditioning_exclude", "status"
    )
  )
})

test_that("direct multivariate random quantities expose density targets", {

  skip_if_missing_fits("brma.mv_v14_assink2016_nested")
  fit <- load_fit("brma.mv_v14_assink2016_nested", validate = FALSE)

  total <- .brma_random_parameter_density_target(
    fit,
    "tau_total"
  )
  nested_sd <- .brma_random_parameter_density_target(
    fit,
    "esid_study: tau(intercept)"
  )
  study_sd <- .brma_random_parameter_density_target(
    fit,
    "study: tau(intercept)"
  )
  nested <- .brma_random_parameter_density_target(
    fit,
    "tau2_prop(esid_study)"
  )
  study <- .brma_random_parameter_density_target(
    fit,
    "tau2_prop(study)"
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
    "tau_total",
    component = "random",
    plot_type = "ggplot"
  )
  expect_identical(
    total_plot$scales$get_scales("x")$name,
    "tau_total"
  )

  # The two-component proportion's canonical prior is the Beta(1, 1) share
  # marginal of its Dirichlet(1, 1) weights.
  samples <- .brma_random_parameter_mixed_posterior(
    fit,
    "tau2_prop(esid_study)"
  )
  parameter <- names(samples)[[1L]]
  prior     <- BayesTools::posterior_metadata(samples[[parameter]], "prior_density")
  expect_equal(
    vapply(c(0.2, 0.5, 0.8), function(value) {
      BayesTools::prior_density_ordinate(prior, value)[["log_density"]]
    }, numeric(1)),
    rep(0, 3L),
    tolerance = 1e-10
  )

  component_samples <- .brma_random_parameter_mixed_posterior(
    fit,
    "study: tau(intercept)"
  )
  component_parameter <- names(component_samples)[[1L]]
  expect_s3_class(
    BayesTools::posterior_metadata(
      component_samples[[component_parameter]],
      "prior_density"
    ),
    "prior_linear_density"
  )

  expect_s3_class(
    plot(
      fit,
      "study: tau(intercept)",
      component       = "random",
      density_method  = "qCMDE",
      density_control = list(
        n_points             = 20L,
        samples              = 1500L,
        normalization_points = 200L
      ),
      plot_type = "ggplot"
    ),
    "ggplot"
  )

  expect_error(
    hypothesis(
      fit,
      "`study: tau(intercept)` = 0",
      density_method = "qCMDE"
    ),
    paste0(
      "Point-null Bayes factors are unavailable for allocation-derived ",
      "random-effect quantity 'study: tau' at 0 because zero is a ",
      "nonregular product boundary of the common scale and allocation ",
      "weight. Test 'tau2_prop(study) = 0' to compare omission of this ",
      "component."
    ),
    fixed = TRUE
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
      "tau2_prop(esid_study)",
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


