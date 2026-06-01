# ============================================================================ #
# test-02-iwmde-fast-paths.R
# ============================================================================ #

context("IWMDE fast-path equivalence")

source(testthat::test_path("common-functions.R"))
source(testthat::test_path("helper-iwmde.R"))


test_that("IWMDE density aggregation matches row-wise reference", {

  log_terms <- matrix(c(
    0, -1, -Inf, NA,
    -Inf, -Inf, NA, NA,
    2, 1, 0, -1
  ), nrow = 3, byrow = TRUE)
  active_mass <- .75
  denominator <- ncol(log_terms)

  fast <- .iwmde_density_aggregate(
    log_terms   = log_terms,
    active_mass = active_mass,
    denominator = denominator
  )

  y                <- numeric(nrow(log_terms))
  finite_terms     <- integer(nrow(log_terms))
  max_log_ratio    <- rep(Inf, nrow(log_terms))
  ess              <- numeric(nrow(log_terms))
  max_weight_share <- rep(1, nrow(log_terms))
  contributions    <- matrix(0, nrow = nrow(log_terms), ncol = ncol(log_terms))

  for (g in seq_len(nrow(log_terms))) {
    finite <- is.finite(log_terms[g, ])
    finite_terms[g] <- sum(finite)
    if (any(finite)) {
      max_term            <- max(log_terms[g, finite])
      scaled_terms        <- exp(log_terms[g, finite] - max_term)
      sum_scaled_terms    <- sum(scaled_terms)
      y[g]                <- active_mass * exp(max_term) *
        sum_scaled_terms / denominator
      contributions[g, finite] <- active_mass * exp(log_terms[g, finite])
      max_log_ratio[g]    <- max_term - stats::median(log_terms[g, finite])
      ess[g]              <- sum_scaled_terms^2 / sum(scaled_terms^2)
      max_weight_share[g] <- max(scaled_terms) / sum_scaled_terms
    }
  }

  expect_equal(fast[["y"]], y)
  expect_equal(fast[["finite_terms"]], finite_terms)
  expect_equal(fast[["max_log_ratio"]], max_log_ratio)
  expect_equal(fast[["ess"]], ess)
  expect_equal(fast[["max_weight_share"]], max_weight_share)
  expect_equal(fast[["contributions"]], contributions)
})

test_that("IWMDE omega matrix collapse matches row-wise collapse", {

  omega       <- matrix(c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6), nrow = 3, byrow = TRUE)
  global_cuts <- c(0, .025, .05, .50, 1)
  active_cuts <- c(0, .05, 1)

  fast <- .iwmde_collapse_omega_matrix(
    omega       = omega,
    global_cuts = global_cuts,
    active_cuts = active_cuts
  )
  ref <- t(apply(omega, 1, .iwmde_collapse_omega,
    global_cuts = global_cuts,
    active_cuts = active_cuts
  ))

  expect_equal(fast, ref)
  expect_equal(fast[, 1], c(1, 3, 5))
})


test_that("IWMDE omega collapse rejects unequal merged bins", {

  global_cuts <- c(0, .025, .05, .50, 1)
  active_cuts <- c(0, .05, .50, 1)
  omega       <- c(1, .5, .25, .25)

  collapsed <- .iwmde_collapse_omega(
    omega       = omega,
    global_cuts = global_cuts,
    active_cuts = active_cuts
  )

  expect_true(is.na(collapsed[1]))
  expect_equal(collapsed[2:3], c(.25, .25))
})


test_that("IWMDE omega helpers use selection-specific JAGS names", {

  context <- list(
    selection_spec = list(
      n_bins     = 4L,
      p_cuts     = c(0, .025, .05, .50, 1),
      jags_omega = "custom.omega+beta"
    )
  )
  active_setup <- list(
    is_weightfunction = TRUE,
    selection_spec    = list(
      n_bins     = 2L,
      p_cuts     = c(0, .05, 1),
      jags_omega = "active.omega+beta"
    ),
    priors = list(outcome = list(bias = NULL))
  )
  samples <- matrix(
    c(1, 1, .5, .5, 2, 10,
      2, 2, .4, .4, 2, 11),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(NULL, c(
      paste0("custom.omega+beta[", 1:4, "]"),
      "bias_indicator",
      "mu"
    ))
  )

  out <- .iwmde_likelihood_posterior_samples(
    context      = context,
    samples      = samples,
    active_setup = active_setup
  )
  row_omega <- .iwmde_active_omega(
    context      = context,
    row          = samples[1, ],
    active_setup = list(
      selection_spec = context[["selection_spec"]],
      priors         = active_setup[["priors"]]
    )
  )

  expect_equal(samples[, "bias_indicator"], c(2, 2))
  expect_equal(out[, "bias_indicator"], c(1, 1))
  expect_false(any(grepl("^custom\\.omega\\+beta\\[", colnames(out))))
  expect_true(all(paste0("active.omega+beta[", 1:2, "]") %in% colnames(out)))
  expect_equal(out[, "active.omega+beta[1]"], samples[, "custom.omega+beta[1]"])
  expect_equal(out[, "active.omega+beta[2]"], samples[, "custom.omega+beta[3]"])
  expect_equal(out[, "mu"], samples[, "mu"])
  expect_equal(row_omega, as.numeric(samples[1, paste0("custom.omega+beta[", 1:4, "]")]))
})


test_that("IWMDE parameter filters use selection-specific JAGS names", {

  context <- list(
    selection_spec = list(jags_omega = "custom.omega+beta"),
    posterior_samples = matrix(
      0,
      nrow = 2,
      ncol = 4,
      dimnames = list(NULL, c(
        "mu",
        "custom.omega+beta[1]",
        "custom.omega+beta[2]",
        "eta[1]"
      ))
    ),
    indicator_names    = character(),
    flat_prior_list    = list(),
    focal_prior_cache  = new.env(parent = emptyenv()),
    support_cache      = new.env(parent = emptyenv())
  )

  expect_equal(.iwmde_default_parameters(context), "mu")
  expect_true(.iwmde_parameter_is_weightfunction_coordinate(
    parameter = "custom.omega+beta[1]",
    context   = context
  ))
  expect_true(.iwmde_parameter_is_weightfunction_coordinate(
    parameter = "eta[1]",
    context   = context
  ))
  expect_false(.iwmde_parameter_is_weightfunction_coordinate(
    parameter = "log_omega[1]",
    context   = context
  ))
  expect_true(.iwmde_chen_is_global_conditioning_column(
    context = context,
    column  = "custom.omega+beta[1]"
  ))
  expect_false(.iwmde_chen_is_global_conditioning_column(
    context = context,
    column  = "eta[1]"
  ))
  expect_false(.iwmde_chen_is_global_conditioning_column(
    context = context,
    column  = "log_omega[1]"
  ))
  expect_equal(
    .iwmde_parameter_spec(context, "custom.omega+beta[1]")[["status"]],
    "unsupported"
  )
})


test_that("IWMDE omega localization renames same-cut active branches", {

  context <- list(
    selection_spec = list(
      n_bins     = 2L,
      p_cuts     = c(0, .05, 1),
      jags_omega = "custom.omega+beta"
    )
  )
  active_setup <- list(
    is_weightfunction = TRUE,
    selection_spec    = list(
      n_bins     = 2L,
      p_cuts     = c(0, .05, 1),
      jags_omega = "active.omega+beta"
    ),
    priors = list(outcome = list(bias = NULL))
  )
  samples <- matrix(
    c(1, .5, 2, 10,
      2, .4, 2, 11),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(NULL, c(
      paste0("custom.omega+beta[", 1:2, "]"),
      "bias_indicator",
      "mu"
    ))
  )

  out <- .iwmde_likelihood_posterior_samples(
    context      = context,
    samples      = samples,
    active_setup = active_setup
  )
  out_again <- .iwmde_likelihood_posterior_samples(
    context      = context,
    samples      = out,
    active_setup = active_setup
  )
  with_active <- cbind(
    samples,
    "active.omega+beta[1]" = c(.9, .8),
    "active.omega+beta[2]" = c(.7, .6)
  )
  active_first <- .iwmde_likelihood_posterior_samples(
    context      = context,
    samples      = with_active,
    active_setup = active_setup
  )

  expect_false(any(grepl("^custom\\.omega\\+beta\\[", colnames(out))))
  expect_true(all(paste0("active.omega+beta[", 1:2, "]") %in% colnames(out)))
  expect_equal(out[, "active.omega+beta[1]"], samples[, "custom.omega+beta[1]"])
  expect_equal(out[, "active.omega+beta[2]"], samples[, "custom.omega+beta[2]"])
  expect_equal(out[, "bias_indicator"], c(1, 1))
  expect_equal(out_again, out)
  expect_equal(active_first[, "active.omega+beta[1]"], c(.9, .8))
  expect_equal(active_first[, "active.omega+beta[2]"], c(.7, .6))
  expect_false(any(grepl("^custom\\.omega\\+beta\\[", colnames(active_first))))
})


test_that("IWMDE active omega does not truncate unexpected omega lengths", {

  context <- list(
    selection_spec = list(
      n_bins     = 4L,
      p_cuts     = c(0, .025, .05, .50, 1),
      jags_omega = "omega"
    )
  )
  active_setup <- list(
    selection_spec = list(
      n_bins     = 2L,
      p_cuts     = c(0, .05, 1),
      jags_omega = "omega"
    ),
    priors = list(outcome = list(bias = NULL))
  )
  row <- c("omega[1]" = 1, "omega[2]" = .5, "omega[3]" = .25)

  omega <- .iwmde_active_omega(
    context      = context,
    row          = row,
    active_setup = active_setup
  )

  expect_true(all(is.na(omega)))
  expect_length(omega, 2L)
})

test_that("IWMDE active omega fails closed when weights are unavailable", {

  context <- list(
    selection_spec = list(
      n_bins     = 2L,
      p_cuts     = c(0, .05, 1),
      jags_omega = "global.omega"
    )
  )
  active_setup <- list(
    selection_spec = list(
      n_bins     = 2L,
      p_cuts     = c(0, .05, 1),
      jags_omega = "active.omega"
    ),
    priors = list(outcome = list(bias = NULL))
  )

  omega <- .iwmde_active_omega(
    context      = context,
    row          = c(mu = 0),
    active_setup = active_setup
  )

  expect_equal(omega, rep(NA_real_, 2L))

  active_setup[["selection_spec"]][["jags_data"]] <- list(
    active.omega = c(1, .5)
  )
  omega <- .iwmde_active_omega(
    context      = context,
    row          = c(mu = 0),
    active_setup = active_setup
  )

  expect_equal(omega, c(1, .5))
})

test_that("IWMDE selection context uses localized active-branch indicators", {

  prior <- BayesTools::prior_weightfunction(
    side    = "one-sided",
    steps   = c(.025, .05),
    weights = BayesTools::wf_fixed(c(1, .7, .35))
  )
  data <- list(
    outcome = data.frame(
      yi  = c(.1, .2),
      sei = c(.1, .2)
    )
  )
  attr(data, "effect_direction") <- "positive"

  context <- list(
    object = list(fit = list()),
    data   = data
  )
  active_setup <- list(
    priors            = list(outcome = list(bias = prior)),
    is_weightfunction = TRUE
  )
  samples <- matrix(
    c(2, 2),
    nrow     = 2,
    dimnames = list(NULL, "bias_indicator")
  )

  selection_context <- .iwmde_selection_context_active_branch(
    context           = context,
    active_setup      = active_setup,
    posterior_samples = samples
  )

  expect_equal(selection_context[["bias_indicator"]], c(1L, 1L))
  expect_false(any(selection_context[["use_normal"]]))
  expect_true(all(selection_context[["kernel_mode"]] != SELKERNEL_NORMAL))
  expect_equal(samples[, "bias_indicator"], c(2, 2))
})

test_that("IWMDE active-branch samples reject uncollapsed prior mixtures", {

  prior <- BayesTools::prior_mixture(
    prior_list = list(
      BayesTools::prior("spike", parameters = list(location = 0)),
      BayesTools::prior("normal", parameters = list(mean = 0, sd = 1))
    ),
    is_null = c(TRUE, FALSE)
  )
  samples <- matrix(
    2,
    nrow     = 1,
    dimnames = list(NULL, "mu_indicator")
  )

  expect_error(
    .iwmde_likelihood_posterior_samples(
      context      = list(),
      samples      = samples,
      active_setup = list(
        priors            = list(outcome = list(mu = prior)),
        is_weightfunction = FALSE
      )
    ),
    "requires priors localized to a single mixture component"
  )
})

test_that("IWMDE batched q evaluation warns on invalid likelihood length", {

  context <- list(
    posterior_samples = matrix(
      0,
      nrow     = 1,
      dimnames = list(NULL, "mu")
    ),
    row_cache = new.env(parent = emptyenv())
  )
  row_states <- list(
    list(
      row_index       = 1L,
      active_key      = "same",
      active_setup    = list(),
      likelihood_mode = "marginal",
      prior_list      = list()
    )
  )

  expect_warning(
    out <- .iwmde_log_q_grid_from_samples(
      context         = context,
      parameter       = "mu",
      values          = c(0, 1),
      row_states      = row_states,
      replacement     = list(type = "scalar"),
      likelihood_mode = "marginal",
      log_lik_fun     = function(samples, active_setup) 0
    ),
    "invalid length"
  )
  expect_null(out)
})


test_that("IWMDE replacement keeps inverse-gamma auxiliary parameters synced", {

  prior <- BayesTools::prior(
    "invgamma",
    parameters = list(2, .5)
  )
  context <- list(flat_prior_list = list(tau = prior))
  samples <- matrix(
    c(.5, 99,
      2, 99),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(NULL, c("tau", "inv_tau"))
  )

  synced <- .iwmde_sync_invgamma_auxiliary_matrix(
    context    = context,
    samples    = samples,
    parameters = "tau"
  )
  bad <- .iwmde_sync_invgamma_auxiliary_matrix(
    context    = context,
    samples    = samples,
    parameters = "mu"
  )
  row <- .iwmde_sync_invgamma_auxiliary_row(
    context    = context,
    row        = c(tau = 4, inv_tau = 99),
    parameters = "tau"
  )
  invalid <- .iwmde_sync_invgamma_auxiliary_row(
    context    = context,
    row        = c(tau = 0, inv_tau = 99),
    parameters = "tau"
  )

  expect_false(.iwmde_can_use_focal_prior_delta(prior))
  expect_true(all(synced[["valid"]]))
  expect_equal(synced[["samples"]][, "inv_tau"], 1 / samples[, "tau"])
  expect_equal(bad[["samples"]], samples)
  expect_true(isTRUE(row[["valid"]]))
  expect_equal(row[["row"]][["inv_tau"]], .25)
  expect_false(invalid[["valid"]])
  expect_true(is.na(invalid[["row"]][["inv_tau"]]))
})


test_that("IWMDE batched q evaluation matches scalar fallback", {

  cases <- list(
    list("bcg_meta-regression", "mu_ablat", NULL),
    list("bangertdrowns2004_location-scale", "log_tau_intercept", NULL),
    list("dat.lehmann2018-3PSM", "mu", NULL),
    list("dat.lehmann2018_RoBMA_3lvl_mods_scale", "mu_Preregistered", NULL),
    list("bcg_BMA.glmm_3lvl_location_scale", "mu_year", NULL),
    list("konstantopoulos2011_3lvl", "gamma[1]", NULL),
    list("bcg_glmm", "theta[1]", NULL),
    list("bcg_glmm", "pi[1]", NULL),
    list("nielweise2008_glmm", "phi[1]", NULL),
    list(
      "bcg_meta-regression",
      "mu_intercept + mu_ablat",
      list(type = "linear", weights = c(mu_intercept = 1, mu_ablat = 1))
    )
  )

  for (case in cases) {
    .expect_iwmde_batch_equals_scalar(
      fit_name       = case[[1L]],
      parameter      = case[[2L]],
      parameter_spec = case[[3L]]
    )
  }
})

test_that("IWMDE predictor fast path matches scalar fallback", {

  cases <- list(
    list("bcg_meta-analysis", "mu", NULL),
    list("bcg_meta-analysis", "tau", NULL),
    list("konstantopoulos2011_3lvl", "rho", NULL),
    list("dat.lehmann2018-PET", "PET", NULL),
    list("dat.lehmann2018-PEESE", "PEESE", NULL),
    list("bcg_meta-regression", "mu_ablat", NULL),
    list("bangertdrowns2004_location-scale", "log_tau_intercept", NULL),
    list(
      "bcg_meta-regression",
      "mu_intercept + mu_ablat",
      list(type = "linear", weights = c(mu_intercept = 1, mu_ablat = 1))
    )
  )

  for (case in cases) {
    .expect_iwmde_predictor_fast_equals_scalar(
      fit_name       = case[[1L]],
      parameter      = case[[2L]],
      parameter_spec = case[[3L]]
    )
  }
})

test_that("IWMDE normal quadratic fast path matches scalar fallback", {

  cases <- list(
    list("bcg_meta-analysis", "mu", NULL),
    list("dat.lehmann2018-3PSM", "mu", NULL),
    list("dat.lehmann2018-PET", "PET", NULL),
    list("dat.lehmann2018-PEESE", "PEESE", NULL),
    list("bcg_meta-regression", "mu_ablat", NULL),
    list("konstantopoulos2011_3lvl", "mu", NULL),
    list(
      "bcg_meta-regression",
      "mu_intercept + mu_ablat",
      list(type = "linear", weights = c(mu_intercept = 1, mu_ablat = 1))
    )
  )

  for (case in cases) {
    .expect_iwmde_normal_quadratic_equals_scalar(
      fit_name       = case[[1L]],
      parameter      = case[[2L]],
      parameter_spec = case[[3L]]
    )
  }
})
