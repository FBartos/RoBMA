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


.expect_mv_heterogeneity_summary <- function(x,
                                             info = NULL,
                                             expected_rows = c("sd", "var")) {

  expect_true(inherits(x, "summary_heterogeneity.brma"), info = info)
  expect_equal(sort(rownames(x[["estimates"]])), sort(expected_rows), info = info)
  estimates <- as.matrix(x[["estimates"]][expected_rows, , drop = FALSE])
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
      "mu__xRE_ALLOCx_heterogeneity__allocation_sd",
      "mu__xRE_ALLOCx_heterogeneity__weight[1]",
      "mu__xRE_ALLOCx_heterogeneity__weight[2]",
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
    "mu__xRE_ALLOCx_heterogeneity__allocation_sd",
    "mu__xRE_ALLOCx_heterogeneity__weight[1]",
    "mu__xRE_ALLOCx_heterogeneity__weight[2]",
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

  expect_named(all_components, names(components))
  for (component in all_components) {
    component_df <- as.data.frame(component)
    expect_identical(
      names(component_df),
      c("component", "parameter", "Mean", "Median", "CI_0.025", "CI_0.975")
    )
    expect_equal(
      unname(as.matrix(component_df)),
      unname(as.matrix(as.data.frame(summary(component)))),
      tolerance = 0
    )
  }
  component_table  <- as.data.frame(all_components)
  component_tables <- as.data.frame(all_components, format = "list")
  expect_s3_class(all_components, "brma_samples_list")
  expect_s3_class(component_table, "data.frame")
  expect_named(
    component_table[1:2],
    c("component", "parameter")
  )
  expect_setequal(component_table[["component"]], names(all_components))
  expect_true(all(component_table[["parameter"]] == "sd"))
  expect_type(component_tables, "list")
  expect_named(component_tables, names(all_components))
  expect_true(all(vapply(component_tables, is.data.frame, logical(1))))
  expect_identical(
    data.frame(all_components),
    data.frame(component_table)
  )
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


test_that("brma.mv pooled heterogeneity validates forwarded arguments", {

  object <- .brma_mv_heterogeneity_named_object()
  posterior_samples <- .brma_mv_two_component_allocation_posterior(
    first_sd      = c(0.20, 0.25),
    second_sd     = c(0.30, 0.35),
    first_column  = "mu__xREx__Study_effects_intercept",
    second_column = "mu__xREx__Effect_effects_intercept"
  )

  expect_no_error(pooled_heterogeneity(
    object,
    component          = "Study effects",
    .posterior_samples = posterior_samples
  ))
  expect_error(
    pooled_heterogeneity(
      object,
      component          = "Study effects",
      posterior_sample   = posterior_samples
    ),
    "Unused argument.*posterior_sample"
  )
})


test_that("brma.mv pooled heterogeneity evaluates the average random design", {

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
  row_components <- .brma_mv_heterogeneity_components(
    object            = object,
    posterior_samples = posterior_samples
  )
  row_rms <- matrix(
    .brma_mv_rms_sd_samples(row_components[["study"]]),
    ncol = 1L
  )
  expected <- matrix(
    posterior_samples[, "mu__xREx__study_intercept"],
    ncol = 1L
  )

  pooled <- pooled_heterogeneity(
    object,
    component          = "study",
    .posterior_samples = posterior_samples
  )

  expect_equal(unname(as.matrix(pooled)), unname(expected), tolerance = 1e-12)
  expect_false(isTRUE(all.equal(unname(expected), unname(row_rms))))
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
  study <- summary_heterogeneity(
    object,
    component          = "study",
    .posterior_samples = posterior_samples
  )

  expect_named(summaries, c("study/esid", "esid_study", "study"))
  expect_equal(
    rownames(summaries[["study/esid"]][["estimates"]]),
    c(
      "sd_total",
      "var_total",
      "var_prop(esid_study)",
      "var_prop(study)"
    )
  )
  expect_equal(
    summaries[["study/esid"]][["estimates"]]["sd_total", "Mean"],
    mean(c(0.50, 0.80)),
    tolerance = 1e-12
  )
  expect_equal(
    summaries[["study/esid"]][["estimates"]]["var_total", "Mean"],
    mean(c(0.50^2, 0.80^2)),
    tolerance = 1e-12
  )
  expect_equal(
    summaries[["study/esid"]][["estimates"]][
      "var_prop(esid_study)", "Mean"
    ],
    mean(c(0.64, 0.25)),
    tolerance = 1e-12
  )
  expect_s3_class(allocation, "summary_heterogeneity.brma")
  expect_error(
    summary_heterogeneity(
      object,
      component          = "heterogeneity",
      .posterior_samples = posterior_samples
    ),
    "Available components"
  )
  expect_error(
    pooled_heterogeneity(
      object,
      component          = "missing",
      .posterior_samples = posterior_samples
    ),
    "esid_study"
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
      component          = "component 1",
      .posterior_samples = posterior_samples
    ),
    "ambiguous"
  )
})


test_that("brma.mv heterogeneity falls back to row-marginal random SDs", {

  old_options <- options(RoBMA.known_v_covariance_max_bytes = 1)
  on.exit(options(old_options), add = TRUE)

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
    dimnames = list(NULL, "mu__xRE_ALLOCx_heterogeneity__allocation_sd")
  )
  expect_error(
    .random_effect_term_sd_samples(
      term              = term,
      posterior_samples = allocation_samples,
      K                 = nobs(object)
    ),
    "missing posterior column"
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
    blocks            = "study",
    diagonal_only     = TRUE
  )
  expected <- sqrt(random_vcov[["samples"]])

  expect_named(components, "study")
  expect_equal(components[["study"]], expected, tolerance = 1e-12)
})


test_that("brma.mv validates marginal random variance schemas", {

  valid <- list(
    samples = matrix(
      c(0.04, 0.09, 0.16, 0.25),
      nrow = 2,
      dimnames = list(NULL, c("1", "2"))
    ),
    metadata = list(
      representation = "diagonal",
      quantity       = "variance",
      diagonal_only  = TRUE,
      dense          = FALSE,
      n_draws        = 2L,
      n_rows         = 2L,
      row_order      = 1:2,
      row_names      = c("1", "2")
    )
  )

  expect_equal(
    .brma_mv_validate_random_marginal_variance_samples(valid, 2L, 2L),
    valid[["samples"]]
  )

  invalid <- list(
    representation = function(x) {
      x[["metadata"]][["representation"]] <- "dense"
      x
    },
    row_order = function(x) {
      x[["metadata"]][["row_order"]] <- 2:1
      x
    },
    row_names = function(x) {
      colnames(x[["samples"]]) <- c("2", "1")
      x
    },
    negative = function(x) {
      x[["samples"]][1L, 1L] <- -1e-12
      x
    }
  )
  for (mutate in invalid) {
    expect_error(
      .brma_mv_validate_random_marginal_variance_samples(
        mutate(valid), 2L, 2L
      ),
      "invalid diagonal variance schema",
      fixed = TRUE
    )
  }
})


test_that("brma.mv random covariance adapter prefers compiled priors", {

  object            <- .brma_mv_heterogeneity_slope_object()
  posterior_samples <- matrix(
    c(0.20, 0.30),
    nrow = 1L,
    dimnames = list(NULL, c(
      "mu__xREx__study_intercept",
      "mu__xREx__study_x"
    ))
  )
  formula_design  <- .fitted_formula_design(object, "mu", required = TRUE)
  compiled_priors <- formula_design[["prior_list"]]
  expect_false(is.null(compiled_priors))

  attr(object[["fit"]], "prior_list") <- NULL
  object[["priors"]][["location"]] <- list(raw_fallback = TRUE)
  captured <- NULL
  sentinel <- list(samples = matrix(1, nrow = 1L, ncol = nobs(object)))

  testthat::local_mocked_bindings(
    .fitted_formula_design = function(object, parameter, required = FALSE) {
      formula_design
    },
    .package = "RoBMA"
  )
  testthat::local_mocked_bindings(
    random_effects_marginal_vcov = function(
        fit, parameter, data, posterior_samples, prior_list, blocks,
        new_levels, diagonal_only) {

      captured <<- list(
        prior_list    = prior_list,
        blocks        = blocks,
        new_levels    = new_levels,
        diagonal_only = diagonal_only
      )
      sentinel
    },
    .package = "BayesTools"
  )

  out <- .brma_mv_random_effects_marginal_vcov(
    object            = object,
    posterior_samples = posterior_samples,
    blocks            = "study",
    diagonal_only     = TRUE
  )

  expect_identical(out, sentinel)
  expect_identical(captured[["prior_list"]], compiled_priors)
  expect_identical(captured[["blocks"]], "study")
  expect_null(captured[["new_levels"]])
  expect_true(captured[["diagonal_only"]])
})


test_that("HAR row-marginal heterogeneity scales to all posterior draws", {

  skip_on_cran()

  n_draws <- 40000L
  n_rows  <- 82L
  dat <- data.frame(
    yi    = seq_len(n_rows) / n_rows,
    study = factor(rep("s1", n_rows)),
    time  = seq_len(n_rows)
  )
  object <- brma.mv(
    yi                        = yi,
    V                         = diag(rep(0.04, n_rows)),
    data                      = dat,
    random                    = ~ har(time | study),
    known_v_parameterization  = "block_mvn",
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  posterior_samples <- cbind(
    mu__xRE_ALLOCx_heterogeneity__allocation_sd = rep(2, n_draws),
    matrix(
      1 / n_rows,
      nrow = n_draws,
      ncol = n_rows,
      dimnames = list(
        NULL,
        paste0(
          "mu__xRE_ALLOCx_heterogeneity__weight[",
          seq_len(n_rows),
          "]"
        )
      )
    )
  )
  old_options <- options(RoBMA.known_v_covariance_max_bytes = 1)
  on.exit(options(old_options), add = TRUE)

  samples <- .brma_mv_random_block_row_sd_samples(
    object            = object,
    posterior_samples = posterior_samples,
    block             = "study",
    K                 = n_rows,
    original_error    = simpleError("forced fallback")
  )

  expect_equal(dim(samples), c(n_draws, n_rows))
  expect_true(all(is.finite(samples)))
  expect_equal(
    unname(samples[c(1L, n_draws), c(1L, n_rows)]),
    matrix(2, nrow = 2L, ncol = 2L),
    tolerance = 1e-12
  )
  expect_lt(as.numeric(utils::object.size(samples)), 40 * 1024^2)
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
      "mu__xRE_ALLOCx_heterogeneity__allocation_sd",
      "mu__xRE_ALLOCx_heterogeneity__weight[1]",
      "mu__xRE_ALLOCx_heterogeneity__weight[2]"
    ))
  )

  allocation_summaries <- .summary_heterogeneity_brma_mv_allocations(
    object            = object,
    posterior_samples = posterior_samples,
    probs             = c(.025, .975)
  )
  pooled <- pooled_heterogeneity(
    object,
    .posterior_samples = posterior_samples
  )

  expect_named(allocation_summaries, "study")
  expect_true(
    any(grepl("sd(time[1])",
              rownames(allocation_summaries[["study"]][["estimates"]]),
              fixed = TRUE))
  )
  expect_true(
    any(grepl("var(time[1])",
              rownames(allocation_summaries[["study"]][["estimates"]]),
              fixed = TRUE))
  )
  expect_true(
    any(grepl("var_ratio(time[2])",
              rownames(allocation_summaries[["study"]][["estimates"]]),
              fixed = TRUE))
  )
  expect_equal(
    allocation_summaries[["study"]][["estimates"]]["sd_common", "Mean"],
    3
  )
  expect_equal(
    allocation_summaries[["study"]][["estimates"]]["var_common", "Mean"],
    10
  )
  expect_identical(colnames(pooled), "sd_common")
  expect_equal(as.numeric(pooled[, 1L]), c(2, 4))
  expect_equal(
    allocation_summaries[["study"]][["estimates"]]["var(time[1])", "Mean"],
    9
  )
  expect_equal(
    allocation_summaries[["study"]][["estimates"]]["var_ratio(time[1])", "Mean"],
    0.75
  )
  expect_true(
    "sd_ratio(time[1])" %in%
      rownames(allocation_summaries[["study"]][["estimates"]])
  )
  expect_equal(
    rownames(allocation_summaries[["study"]][["estimates"]]),
    c(
      "sd_common", "var_common", "sd(time[1])", "sd(time[2])",
      "var(time[1])", "var(time[2])", "sd_ratio(time[1])",
      "sd_ratio(time[2])", "var_ratio(time[1])", "var_ratio(time[2])"
    )
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
    unname(samples[["study: sd(time[0])"]]),
    sqrt(mean(c(1, 2)^2) * (2 * .25)),
    tolerance = 1e-12
  )
  expect_equal(
    unname(samples[["study: sd(time[1])"]]),
    sqrt(mean(c(10, 20)^2) * (2 * .75)),
    tolerance = 1e-12
  )
})


test_that("brma.mv heterogeneity supports single known-V component", {

  name <- "brma.mv_block_mvn"
  skip_if_missing_fits(name)

  fit_brma <- fits[[name]]

  pooled <- pooled_heterogeneity(fit_brma)
  direct <- .pooled_brma_mv_heterogeneity_components(fit_brma)[["estimate"]]
  total  <- pooled_heterogeneity(fit_brma, component = "total")
  estimate <- pooled_heterogeneity(fit_brma, component = "estimate")
  summary <- summary_heterogeneity(fit_brma)

  expect_brma_samples_matrix(pooled, 1, paste(name, "pooled heterogeneity"))
  expect_equal(unname(as.matrix(pooled)), unname(as.matrix(direct)),
               tolerance = 1e-12)
  expect_equal(unname(as.matrix(total)), unname(as.matrix(pooled)),
               tolerance = 1e-12)
  expect_equal(unname(as.matrix(estimate)), unname(as.matrix(pooled)),
               tolerance = 1e-12)
  .expect_mv_heterogeneity_summary(summary, paste(name, "summary"))
})


test_that("brma.mv heterogeneity decomposes random-formula SD components", {

  name <- "brma.mv_block_mvn_random_scale"
  skip_if_missing_fits(name)

  fit_brma          <- fits[[name]]
  posterior_samples <- .get_posterior_samples(fit_brma[["fit"]])
  components        <- .pooled_brma_mv_heterogeneity_components(
    object            = fit_brma,
    posterior_samples = posterior_samples
  )

  pooled <- pooled_heterogeneity(fit_brma)
  total  <- pooled_heterogeneity(fit_brma, component = "total")
  study  <- pooled_heterogeneity(fit_brma, component = "study")

  expect_type(pooled, "list")
  expect_named(pooled, names(components))
  for (component in names(components)) {
    expected <- components[[component]]
    expect_brma_samples_matrix(
      pooled[[component]],
      1,
      paste(name, component, "pooled heterogeneity")
    )
    expect_equal(unname(as.matrix(pooled[[component]])), unname(expected),
                 tolerance = 1e-12)
  }

  total_expected <- matrix(
    sqrt(Reduce(`+`, lapply(components, function(samples) samples[, 1L]^2))),
    ncol = 1L
  )
  colnames(total_expected) <- "sd_total"
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
    expect_true(all(
      c("var_prop(effect_study)", "var_prop(study)") %in%
        rownames(summaries[[component]][["estimates"]])
    ),
      info = paste(name, component, "allocation summary")
    )
  }
  .expect_mv_heterogeneity_summary(
    total,
    paste(name, "total summary"),
    expected_rows = c("sd_total", "var_total")
  )
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
