context("brma.mv heterogeneity")

source(testthat::test_path("common-functions.R"))
source(testthat::test_path("helper-contracts.R"))

fit_names <- list_fits()
fits      <- lazy_fits(fit_names, validate = FALSE)

.brma_mv_heterogeneity_named_object <- function() {

  dat <- data.frame(
    yi     = c(0.10, 0.20, 0.30, 0.40),
    study  = c("s1", "s1", "s2", "s2"),
    effect = paste0("e", 1:4)
  )

  brma.mv(
    yi                        = yi,
    V                         = diag(rep(0.04, 4)),
    data                      = dat,
    random                    = list(
      "Study effects"  = ~ 1 | study,
      "Effect effects" = ~ 1 | effect
    ),
    known_v_parameterization  = "block_mvn",
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
}


.brma_mv_heterogeneity_slope_object <- function() {

  dat <- data.frame(
    yi    = c(0.10, 0.20, 0.30, 0.40),
    study = c("s1", "s1", "s2", "s2"),
    x     = c(0, 1, 0, 1)
  )

  brma.mv(
    yi                        = yi,
    V                         = diag(rep(0.04, 4)),
    data                      = dat,
    random                    = ~ diag(1 + x | study),
    known_v_parameterization  = "block_mvn",
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
}


.brma_mv_heterogeneity_nested_object <- function() {

  dat <- data.frame(
    yi    = c(0.10, 0.20, 0.30, 0.40),
    study = c("s1", "s1", "s2", "s2"),
    esid  = paste0("e", 1:4)
  )

  brma.mv(
    yi                        = yi,
    V                         = diag(rep(0.04, 4)),
    data                      = dat,
    random                    = ~ 1 | study / esid,
    known_v_parameterization  = "block_mvn",
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
}


.expect_mv_heterogeneity_summary <- function(x, info = NULL) {

  expect_true(inherits(x, "summary_heterogeneity.brma"), info = info)
  expect_equal(sort(rownames(x[["estimates"]])), c("tau", "tau2"), info = info)
  estimates <- as.matrix(x[["estimates"]][c("tau", "tau2"), , drop = FALSE])
  expect_true(all(is.finite(estimates)), info = info)
  expect_true(all(estimates >= 0), info = info)
  expect_false("I2" %in% rownames(x[["estimates"]]), info = info)
  expect_false("H2" %in% rownames(x[["estimates"]]), info = info)
}


.brma_mv_nested_allocation_posterior <- function() {

  posterior_samples <- matrix(
    c(
      0.50, 0.64, 0.36, 0.40, 0.30,
      0.80, 0.25, 0.75, 0.40, 0.80 * sqrt(0.75)
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(NULL, c(
      "mu__xRE_ALLOCx_random_total__total_sd",
      "mu__xRE_ALLOCx_random_total__weight[1]",
      "mu__xRE_ALLOCx_random_total__weight[2]",
      "mu__xREx__esid_study_intercept",
      "mu__xREx__study_intercept"
    ))
  )

  posterior_samples
}


.brma_mv_two_component_allocation_posterior <- function(first_sd,
                                                        second_sd,
                                                        first_column,
                                                        second_column) {

  total_sd <- sqrt(first_sd^2 + second_sd^2)
  out <- cbind(
    total_sd,
    first_sd^2 / total_sd^2,
    second_sd^2 / total_sd^2,
    first_sd,
    second_sd
  )
  colnames(out) <- c(
    "mu__xRE_ALLOCx_random_total__total_sd",
    "mu__xRE_ALLOCx_random_total__weight[1]",
    "mu__xRE_ALLOCx_random_total__weight[2]",
    first_column,
    second_column
  )

  out
}


test_that("brma.mv heterogeneity resolves aliases and component errors", {

  object <- .brma_mv_heterogeneity_named_object()
  posterior_samples <- .brma_mv_two_component_allocation_posterior(
    first_sd      = c(0.20, 0.25),
    second_sd     = c(0.30, 0.35),
    first_column  = "mu__xREx__Study_effects_intercept",
    second_column = "mu__xREx__Effect_effects_intercept"
  )
  components <- .brma_mv_heterogeneity_components(
    object            = object,
    posterior_samples = posterior_samples
  )
  allocation_samples <- .brma_mv_allocation_sample_lists(
    object            = object,
    posterior_samples = posterior_samples
  )

  all_components <- pooled_heterogeneity(
    object,
    .posterior_samples = posterior_samples
  )
  study <- pooled_heterogeneity(
    object,
    component          = "Study effects",
    .posterior_samples = posterior_samples
  )
  total <- pooled_heterogeneity(
    object,
    component          = "total",
    .posterior_samples = posterior_samples
  )

  expect_named(all_components, c(names(allocation_samples), names(components)))
  expect_equal(
    unname(as.matrix(study)),
    matrix(c(0.20, 0.25), ncol = 1),
    tolerance = 1e-12
  )
  expect_equal(
    unname(as.matrix(total)),
    matrix(c(sqrt(0.20^2 + 0.30^2), sqrt(0.25^2 + 0.35^2)), ncol = 1),
    tolerance = 1e-12
  )
  expect_error(
    pooled_heterogeneity(
      object,
      component          = "missing",
      .posterior_samples = posterior_samples
    ),
    "Unknown heterogeneity component"
  )
})


test_that("brma.mv pooled heterogeneity uses RMS row aggregation", {

  samples <- matrix(c(1, 3, 5, 7), nrow = 2, byrow = TRUE)
  expected <- matrix(
    c(
      sqrt(mean(c(1^2, 3^2))),
      sqrt(mean(c(5^2, 7^2)))
    ),
    ncol = 1,
    dimnames = list(NULL, "tau")
  )

  expect_equal(
    .pooled_brma_mv_heterogeneity_samples(samples),
    expected,
    tolerance = 1e-12
  )
})


test_that("brma.mv pooled allocation heterogeneity selects total SD by name", {

  samples <- list(
    "var_frac(allocation: first)" = c(.25, .50),
    "sd_total(allocation)"        = c(2, 4),
    "var_frac(allocation: second)" = c(.75, .50)
  )
  expected <- matrix(c(2, 4), ncol = 1, dimnames = list(NULL, "tau"))

  expect_equal(
    .pooled_brma_mv_heterogeneity_samples(samples),
    expected
  )
})


test_that("brma.mv summary heterogeneity reports shared allocation nodes", {

  object            <- .brma_mv_heterogeneity_nested_object()
  posterior_samples <- .brma_mv_nested_allocation_posterior()

  summaries <- summary_heterogeneity(
    object,
    .posterior_samples = posterior_samples
  )
  allocation <- summary_heterogeneity(
    object,
    component          = "study/esid",
    .posterior_samples = posterior_samples
  )
  allocation_by_node <- summary_heterogeneity(
    object,
    component          = "random_total",
    .posterior_samples = posterior_samples
  )
  pooled_allocation <- pooled_heterogeneity(
    object,
    component          = "study/esid",
    .posterior_samples = posterior_samples
  )
  pooled_by_node <- pooled_heterogeneity(
    object,
    component          = "random_total",
    .posterior_samples = posterior_samples
  )
  study <- summary_heterogeneity(
    object,
    component          = "study",
    .posterior_samples = posterior_samples
  )

  expect_named(summaries, c("study/esid", "esid_study", "study"))
  expect_equal(
    rownames(summaries[["study/esid"]][["estimates"]]),
    c(
      "tau",
      "tau2",
      "var_frac(random_total: esid_study)",
      "var_frac(random_total: study)"
    )
  )
  expect_equal(
    summaries[["study/esid"]][["estimates"]]["tau", "Mean"],
    mean(c(0.50, 0.80)),
    tolerance = 1e-12
  )
  expect_equal(
    summaries[["study/esid"]][["estimates"]]["tau2", "Mean"],
    mean(c(0.50^2, 0.80^2)),
    tolerance = 1e-12
  )
  expect_equal(
    summaries[["study/esid"]][["estimates"]][
      "var_frac(random_total: esid_study)", "Mean"
    ],
    mean(c(0.64, 0.25)),
    tolerance = 1e-12
  )
  expect_s3_class(allocation, "summary_heterogeneity.brma")
  expect_equal(allocation[["estimates"]], allocation_by_node[["estimates"]])
  expect_equal(
    unname(as.matrix(pooled_allocation)),
    matrix(c(0.50, 0.80), ncol = 1),
    tolerance = 1e-12
  )
  expect_equal(
    unname(as.matrix(pooled_allocation)),
    unname(as.matrix(pooled_by_node)),
    tolerance = 1e-12
  )
  expect_error(
    pooled_heterogeneity(
      object,
      component          = "missing",
      .posterior_samples = posterior_samples
    ),
    "study/esid"
  )
  .expect_mv_heterogeneity_summary(study, "study leaf summary")
})


test_that("brma.mv heterogeneity reports ambiguous component aliases", {

  dat <- data.frame(
    yi       = c(0.10, 0.20, 0.30, 0.40),
    study    = c("s1", "s1", "s2", "s2"),
    estimate = paste0("e", 1:4)
  )
  object <- brma.mv(
    yi                        = yi,
    V                         = diag(rep(0.04, 4)),
    data                      = dat,
    random                    = ~ 1 | study / estimate,
    known_v_parameterization  = "block_mvn",
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  posterior_samples <- .brma_mv_two_component_allocation_posterior(
    first_sd      = c(0.20, 0.25),
    second_sd     = c(0.30, 0.35),
    first_column  = "mu__xREx__estimate_study_intercept",
    second_column = "mu__xREx__study_intercept"
  )

  expect_error(
    pooled_heterogeneity(
      object,
      component          = "Component 1",
      .posterior_samples = posterior_samples
    ),
    "ambiguous"
  )
})


test_that("brma.mv heterogeneity falls back to row-marginal random SDs", {

  object <- .brma_mv_heterogeneity_slope_object()
  posterior_samples <- matrix(
    c(
      0.20, 0.30,
      0.25, 0.35
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(NULL, c(
      "mu__xREx__study_intercept",
      "mu__xREx__study_x"
    ))
  )

  formula_design <- .fitted_formula_design(object, "mu", required = TRUE)
  term           <- formula_design[["random_effects"]][[1L]]

  allocation_samples <- matrix(
    c(0.50, 0.70),
    nrow = 2,
    dimnames = list(NULL, "mu__xRE_ALLOCx_allocation__total_sd")
  )
  expect_error(
    .random_effect_term_sd_samples(
      term              = term,
      posterior_samples = allocation_samples,
      K                 = nobs(object)
    ),
    "SD-component allocation"
  )

  components <- .brma_mv_heterogeneity_components(
    object            = object,
    posterior_samples = posterior_samples
  )

  formula_fit <- .posterior_formula_fit(
    fit               = object[["fit"]],
    posterior_samples = posterior_samples,
    formula_design    = TRUE
  )
  attr(formula_fit, "formula_design") <- list(mu = formula_design)
  random_vcov <- BayesTools::random_effects_marginal_vcov(
    fit               = formula_fit,
    parameter         = "mu",
    data              = object[["data"]][["location"]],
    posterior_samples = posterior_samples,
    prior_list        = object[["priors"]][["location"]],
    blocks            = "study"
  )
  expected <- matrix(NA_real_, nrow = nrow(posterior_samples),
                     ncol = nobs(object))
  for (s in seq_len(nrow(expected))) {
    expected[s, ] <- sqrt(pmax(diag(random_vcov[["samples"]][s, , ]), 0))
  }

  expect_named(components, "study")
  expect_equal(components[["study"]], expected, tolerance = 1e-12)
})


test_that("brma.mv summary heterogeneity reports SD-component allocation tables", {

  dat <- data.frame(
    yi    = c(0.10, 0.20, 0.30, 0.40),
    study = c("s1", "s1", "s2", "s2"),
    time  = c(1, 2, 1, 2)
  )
  object <- brma.mv(
    yi                        = yi,
    V                         = diag(rep(0.04, 4)),
    data                      = dat,
    random                    = ~ har(time | study),
    known_v_parameterization  = "block_mvn",
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  posterior_samples <- matrix(
    c(
      2.0, 0.25, 0.75,
      4.0, 0.50, 0.50
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(NULL, c(
      "mu__xRE_ALLOCx_allocation__total_sd",
      "mu__xRE_ALLOCx_allocation__weight[1]",
      "mu__xRE_ALLOCx_allocation__weight[2]"
    ))
  )

  allocation_summaries <- .summary_heterogeneity_brma_mv_allocations(
    object            = object,
    posterior_samples = posterior_samples,
    probs             = c(.025, .975)
  )

  expect_named(allocation_summaries, "study")
  expect_true(
    any(grepl("sd(time[1] | study)",
              rownames(allocation_summaries[["study"]][["estimates"]]),
              fixed = TRUE))
  )
  expect_true(
    any(grepl("var_ratio(allocation: time[2])",
              rownames(allocation_summaries[["study"]][["estimates"]]),
              fixed = TRUE))
  )
  expect_equal(
    allocation_summaries[["study"]][["estimates"]]["sd_total(allocation)", "Mean"],
    3
  )
  expect_equal(
    allocation_summaries[["study"]][["estimates"]]["var_ratio(allocation: time[1])", "Mean"],
    0.75
  )
})


test_that("brma.mv SD-component allocation summaries map row SD sources through observations", {

  dat <- data.frame(
    yi    = c(0.10, 0.20, 0.30, 0.40),
    study = c("s1", "s1", "s2", "s2"),
    time  = c(0, 1, 0, 1)
  )
  object <- brma.mv(
    yi                        = yi,
    V                         = diag(rep(0.04, 4)),
    data                      = dat,
    random                    = ~ har(time | study),
    scale                     = ~ time,
    known_v_parameterization  = "block_mvn",
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  design     <- .fitted_formula_design(object, "mu", required = TRUE)
  term       <- design[["random_effects"]][[1L]]
  allocation <- term[["sd_binding"]][["allocations"]][[1L]]
  posterior_samples <- matrix(
    c(.25, .75),
    nrow     = 1L,
    dimnames = list(NULL, paste0(allocation[["weight_name"]], "[", 1:2, "]"))
  )
  source_samples <- list(tau = matrix(c(1, 10, 2, 20), nrow = 1L))

  samples <- .brma_mv_sd_component_allocation_summary_samples(
    allocation        = allocation,
    term              = term,
    posterior_samples = posterior_samples,
    source_samples    = source_samples,
    K                 = 4L
  )

  expect_equal(
    unname(samples[["sd(time[0] | study)"]]),
    sqrt(mean(c(1, 2)^2) * (2 * .25)),
    tolerance = 1e-12
  )
  expect_equal(
    unname(samples[["sd(time[1] | study)"]]),
    sqrt(mean(c(10, 20)^2) * (2 * .75)),
    tolerance = 1e-12
  )
})


test_that("brma.mv heterogeneity supports single known-V component", {

  name <- "brma.mv_block_mvn"
  skip_if_missing_fits(name)

  fit_brma <- fits[[name]]

  pooled <- pooled_heterogeneity(fit_brma)
  direct <- predict(fit_brma, type = "terms.scale", quiet = TRUE)
  direct <- matrix(sqrt(rowMeans(as.matrix(direct)^2)), ncol = 1L)
  total  <- pooled_heterogeneity(fit_brma, component = "total")
  tau    <- pooled_heterogeneity(fit_brma, component = "tau")
  summary <- summary_heterogeneity(fit_brma)

  expect_brma_samples_matrix(pooled, 1, paste(name, "pooled heterogeneity"))
  expect_equal(unname(as.matrix(pooled)), unname(as.matrix(direct)),
               tolerance = 1e-12)
  expect_equal(unname(as.matrix(total)), unname(as.matrix(pooled)),
               tolerance = 1e-12)
  expect_equal(unname(as.matrix(tau)), unname(as.matrix(pooled)),
               tolerance = 1e-12)
  .expect_mv_heterogeneity_summary(summary, paste(name, "summary"))
})


test_that("brma.mv heterogeneity decomposes random-formula SD components", {

  name <- "brma.mv_block_mvn_random_scale"
  skip_if_missing_fits(name)

  fit_brma          <- fits[[name]]
  posterior_samples <- .get_posterior_samples(fit_brma[["fit"]])
  components        <- .brma_mv_heterogeneity_components(
    object            = fit_brma,
    posterior_samples = posterior_samples
  )
  allocation_samples <- .brma_mv_allocation_sample_lists(
    object            = fit_brma,
    posterior_samples = posterior_samples
  )
  pooled_component_names <- c(
    names(allocation_samples),
    names(components)[
      !names(components) %in%
        .brma_mv_allocation_replaced_components(allocation_samples)
    ]
  )

  pooled <- pooled_heterogeneity(fit_brma)
  total  <- pooled_heterogeneity(fit_brma, component = "total")
  study  <- pooled_heterogeneity(fit_brma, component = "study")

  expect_type(pooled, "list")
  expect_named(pooled, pooled_component_names)
  for (component in intersect(names(components), names(pooled))) {
    expected <- matrix(
      .brma_mv_rms_sd_samples(components[[component]]),
      ncol = 1L
    )
    colnames(expected) <- "tau"
    expect_brma_samples_matrix(
      pooled[[component]],
      1,
      paste(name, component, "pooled heterogeneity")
    )
    expect_equal(unname(as.matrix(pooled[[component]])), unname(expected),
                 tolerance = 1e-12)
  }
  for (component in names(allocation_samples)) {
    expected <- .pooled_brma_mv_heterogeneity_samples(
      allocation_samples[[component]]
    )
    expect_equal(unname(as.matrix(pooled[[component]])), unname(expected),
                 tolerance = 1e-12)
  }

  total_expected <- matrix(
    .brma_mv_rms_sd_samples(.total_brma_mv_heterogeneity_samples(components)),
    ncol = 1L
  )
  colnames(total_expected) <- "tau"
  expect_equal(unname(as.matrix(total)), unname(total_expected),
               tolerance = 1e-12)
  expect_equal(unname(as.matrix(study)),
               unname(as.matrix(pooled[["study"]])),
               tolerance = 1e-12)
})


test_that("brma.mv summary heterogeneity returns absolute component summaries", {

  name <- "brma.mv_block_mvn_random_scale"
  skip_if_missing_fits(name)

  fit_brma   <- fits[[name]]
  components <- .brma_mv_heterogeneity_components(fit_brma)
  allocation_summaries <- .summary_heterogeneity_brma_mv_allocations(
    object            = fit_brma,
    posterior_samples = .get_posterior_samples(fit_brma[["fit"]]),
    probs             = c(.025, .975)
  )
  summaries  <- summary_heterogeneity(fit_brma)
  total      <- summary_heterogeneity(fit_brma, component = "total")

  expect_type(summaries, "list")
  expect_named(summaries, c(names(allocation_summaries), names(components)))
  for (component in names(components)) {
    .expect_mv_heterogeneity_summary(
      summaries[[component]],
      paste(name, component, "summary")
    )
  }
  for (component in names(allocation_summaries)) {
    expect_true(
      any(grepl("var_frac(", rownames(summaries[[component]][["estimates"]]),
                fixed = TRUE)),
      info = paste(name, component, "allocation summary")
    )
  }
  .expect_mv_heterogeneity_summary(total, paste(name, "total summary"))
})


test_that("brma.mv summary heterogeneity print omits list headers", {

  name <- "brma.mv_block_mvn_random_scale"
  skip_if_missing_fits(name)

  output <- capture.output(print(summary_heterogeneity(fits[[name]])))

  expect_false(any(grepl("^\\$", output)))
  expect_true(any(grepl("Heterogeneity Estimates (study/effect):", output,
                        fixed = TRUE)))
  expect_true(any(grepl("Heterogeneity Estimates (study):", output,
                        fixed = TRUE)))
  expect_true(any(grepl("Heterogeneity Estimates (effect_study):", output,
                        fixed = TRUE)))
})


test_that("brma.mv random-only one-component heterogeneity simplifies", {

  name <- "brma.mv_block_mvn_random"
  skip_if_missing_fits(name)

  fit_brma <- fits[[name]]

  pooled    <- pooled_heterogeneity(fit_brma)
  estimate  <- pooled_heterogeneity(fit_brma, component = "estimate")
  total     <- pooled_heterogeneity(fit_brma, component = "total")
  summary   <- summary_heterogeneity(fit_brma)

  expect_brma_samples_matrix(pooled, 1, paste(name, "pooled heterogeneity"))
  expect_equal(unname(as.matrix(estimate)), unname(as.matrix(pooled)),
               tolerance = 1e-12)
  expect_equal(unname(as.matrix(total)), unname(as.matrix(pooled)),
               tolerance = 1e-12)
  .expect_mv_heterogeneity_summary(summary, paste(name, "summary"))
})
