context("Semantic random-effect parameters")

source(testthat::test_path("common-functions.R"))

.random_parameter_fit_names <- c(
  "brma.mv_block_mvn_random_scale",
  "brma.mv_block_mvn_known_R",
  "brma.mv_v14_konstantopoulos2011_cs",
  "brma.mv_v14_ishak2007_har",
  "brma.mv_v14_begg1989_study_treatment"
)

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

  skip_if_missing_fits(.random_parameter_fit_names)

  for (name in .random_parameter_fit_names) {
    fit      <- load_fit(name, validate = FALSE)
    bundle   <- .brma_random_parameter_bundle(fit)
    summary  <- summary(fit, random_effects = "standard")
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

  skip_if_missing_fits(.random_parameter_fit_names)

  for (name in .random_parameter_fit_names) {
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

  fit <- load_fit("brma.mv_v14_konstantopoulos2011_cs", validate = FALSE)
  for (type in c("density", "autocorrelation")) {
    expect_s3_class(
      plot_diagnostic(
        fit, parameter = "rho(district)", component = "random",
        type = type, plot_type = "ggplot"
      ),
      "ggplot"
    )
  }
})

test_that("random directional hypotheses use induced joint priors", {

  skip_if_missing_fits(.random_parameter_fit_names)

  for (i in seq_along(.random_parameter_fit_names)) {
    name     <- .random_parameter_fit_names[[i]]
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

  skip_if_missing_fits(.random_parameter_fit_names)

  fit_cs <- load_fit("brma.mv_v14_konstantopoulos2011_cs", validate = FALSE)
  expect_s3_class(
    hypothesis(
      fit_cs,
      .random_parameter_hypothesis("rho(district)", "=", 0.5, "!="),
      component = "random", n_samples = 2000, seed = 21
    ),
    "BayesTools_hypothesis_BF"
  )

  fit_har <- load_fit("brma.mv_v14_ishak2007_har", validate = FALSE)
  expect_error(
    hypothesis(
      fit_har,
      .random_parameter_hypothesis(
        "sd(time[1] | study)", "=", 0.2, "!="
      ),
      component = "random", n_samples = 1000, seed = 22
    ),
    "derived component SD"
  )

  fit_mixed <- load_fit(
    "brma.mv_v14_begg1989_study_treatment",
    validate = FALSE
  )
  mixed_quantities <- hypothesis_quantities(fit_mixed)
  mixed_rho <- mixed_quantities[
    mixed_quantities[["alias"]] == "rho(study)" &
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
      .random_parameter_hypothesis("rho(study)", "=", 0, "!="),
      component = "random", n_samples = 1000, seed = 23
    ),
    "derived pairwise correlation"
  )

  fit_alloc <- load_fit("brma.mv_block_mvn_random_scale", validate = FALSE)
  expect_error(
    hypothesis(
      fit_alloc,
      .random_parameter_hypothesis(
        "var_frac(total: study)", "=", 0, ">"
      ),
      component = "random", n_samples = 1000, seed = 24
    ),
    "support boundary"
  )
  expect_error(
    hypothesis(
      fit_alloc,
      .random_parameter_hypothesis(
        "var_frac(total: study)", ">", 0.5, "<"
      ),
      component = "random", density_method = "qCMDE",
      n_samples = 1000, seed = 25
    ),
    "not available for semantic random-effect quantities"
  )
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
    parameter = "rho(study)",
    .weights  = weights
  )
  covr <- covratio(
    fit,
    component = "random",
    parameter = "rho(study)",
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
    "not available|Could not select"
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
    "Could not infer a model parameter"
  )
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
  if (!nonlinear_joint) {
    skip("Cached fit does not contain a nonlinear joint scale transform.")
  }

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
          max_samples          = 40L,
          normalization_points = 40L
        ),
        n_samples = 1000L,
        seed      = 32
      ),
      "nonlinear joint transform"
    )
  }
})

test_that("random prior overlays and diagnostic labels are semantic", {

  skip_if_missing_fits(c(
    "brma.mv_v14_begg1989_study_treatment",
    "brma.mv_v14_konstantopoulos2011_cs"
  ))
  fixed <- load_fit(
    "brma.mv_v14_begg1989_study_treatment",
    validate = FALSE
  )
  expect_error(
    plot(
      fixed,
      parameter = "rho(study)", component = "random", prior = TRUE,
      plot_type = "ggplot"
    ),
    "prior-density overlay is not available for fixed random-effect quantity"
  )

  cs <- load_fit("brma.mv_v14_konstantopoulos2011_cs", validate = FALSE)
  diagnostic <- plot_diagnostic_density(
    cs,
    parameter = "rho(district)", component = "random",
    plot_type = "ggplot"
  )
  expect_identical(diagnostic[["labels"]][["title"]], "rho(district)")
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
