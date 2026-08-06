# ============================================================================ #
# test-02-iwmde-api.R
# ============================================================================ #

context("IWMDE public API")

source(testthat::test_path("common-functions.R"))
source(testthat::test_path("helper-iwmde.R"))


.load_raw_fit_or_skip <- function(name) {

  if (!name %in% list_fits(validate = FALSE, active_only = TRUE)) {
    skip(paste0("Raw cached fit unavailable for the active profile: ", name))
  }
  fit <- try(load_fit(name, validate = FALSE), silent = TRUE)
  if (inherits(fit, "try-error")) {
    skip(paste0("Raw cached fit unavailable: ", name))
  }

  fit
}


.mock_random_non_known_v_iwmde_object <- function() {

  data <- structure(list(), random = TRUE)
  structure(
    list(
      fit    = list(dummy = TRUE),
      data   = data,
      priors = list()
    ),
    class = c("brma.mv", "brma")
  )
}


test_that("density methods are case-insensitive", {

  expect_equal(.density_method_normalize("kde"), "KDE")
  expect_equal(.density_method_normalize("QcMdE"), "qCMDE")
  expect_equal(.density_method_normalize("iwmde"), "IWMDE")
  expect_equal(.density_method_normalize("NORMAL", allow_normal = TRUE), "normal")
  expect_equal(.density_method_normalize_precomputed("qcmde"), "qCMDE")
  expect_equal(.density_method_normalize_precomputed("Iwmde"), "IWMDE")
})


test_that("IWMDE plan marks non-known-V random-formula contexts unsupported", {

  context <- list(
    data              = .mock_random_non_known_v_iwmde_object()[["data"]],
    posterior_samples = matrix(
      seq_len(25),
      ncol     = 1L,
      dimnames = list(NULL, "mu")
    ),
    flat_prior_list = list()
  )
  plan <- .iwmde_plan(
    context         = context,
    parameter       = "mu",
    density_method  = "qCMDE",
    density_control = list(n_points = 20, max_samples = 20),
    outputs         = "density",
    parameter_spec  = list(type = "primitive")
  )

  expect_equal(plan[["status"]], "unsupported")
  expect_match(plan[["reason"]], "random-formula")

  known_v_context <- context
  attr(known_v_context[["data"]], "known_V") <- TRUE
  expect_null(.iwmde_context_unavailable_reason(known_v_context))
})

test_that("IWMDE histogram uses global active-sample mass", {

  values <- c(-2, -1, -.5, 0, .5, 1, 2)
  hist_data <- .iwmde_histogram(values, xlim = c(-1, 1), mass = .8)
  hist_mass <- sum(hist_data[["density"]] * diff(hist_data[["breaks"]]))

  expect_equal(hist_mass, .8 * 5 / 7, tolerance = 1e-12)

  empty <- .iwmde_histogram(values, xlim = c(3, 4), mass = .8)
  expect_equal(sum(empty[["density"]] * diff(empty[["breaks"]])), 0)
  expect_true(all(empty[["density"]] == 0))

  all_missing <- .iwmde_histogram(c(NA_real_, NaN), xlim = c(0, 1), mass = .8)
  expect_equal(sum(all_missing[["density"]] * diff(all_missing[["breaks"]])), 0)
  expect_true(all(all_missing[["density"]] == 0))
})


test_that("plot qCMDE/IWMDE matrix samples use exact parameter columns", {

  sample <- matrix(
    c(1, 2, 3, 4),
    nrow     = 2,
    dimnames = list(NULL, c("theta_extra", "theta"))
  )
  samples <- list(theta = sample)

  expect_equal(
    .plot_brma_plotted_samples(samples, "theta", "theta"),
    c(3, 4)
  )

  sample_missing <- sample[, 1L, drop = FALSE]
  colnames(sample_missing) <- "theta_extra"
  samples <- list(theta = sample_missing)

  expect_null(.plot_brma_plotted_samples(samples, "theta", "theta"))
})


test_that("fixed nonzero priors fill missing evaluator columns", {

  priors <- list(
    outcome = list(
      mu  = BayesTools::prior("spike", list(0.20)),
      tau = BayesTools::prior("spike", list(0.05))
    )
  )
  samples <- matrix(
    seq_len(4) / 10,
    ncol     = 1L,
    dimnames = list(NULL, "theta")
  )

  expect_equal(.fixed_mu_prior_value(priors), 0.20)
  expect_equal(.fixed_tau_prior_value(priors), 0.05)

  factor_prior <- BayesTools::prior_factor(
    "spike",
    parameters = list(0.10),
    contrast   = "treatment"
  )
  expect_equal(.point_prior_value(factor_prior), 0.10)
  expect_equal(
    .fixed_prior_list_values(list(mod_factor = factor_prior)),
    c(mod_factor = 0.10)
  )

  mu_result <- .evaluate.brma.mu(
    fit               = NULL,
    outcome_data      = data.frame(sei = c(0.10, 0.20, 0.30)),
    mods_data         = NULL,
    mods_formula      = NULL,
    mods_priors       = NULL,
    is_mods           = FALSE,
    is_PET            = FALSE,
    is_PEESE          = FALSE,
    effect_direction  = "positive",
    bias_adjusted     = FALSE,
    K                 = 3L,
    posterior_samples = samples,
    priors            = priors
  )
  expect_equal(mu_result, matrix(0.20, nrow = 4L, ncol = 3L))

  expect_error(
    .evaluate.brma.tau(
      fit               = NULL,
      scale_data        = NULL,
      scale_formula     = NULL,
      scale_priors      = NULL,
      is_scale          = FALSE,
      is_multilevel     = FALSE,
      K                 = 3L,
      posterior_samples = samples
    ),
    "Missing posterior tau columns"
  )

  tau_result <- .evaluate.brma.tau(
    fit               = NULL,
    scale_data        = NULL,
    scale_formula     = NULL,
    scale_priors      = NULL,
    is_scale          = FALSE,
    is_multilevel     = FALSE,
    K                 = 3L,
    posterior_samples = samples,
    fixed_tau         = .fixed_tau_prior_value(priors)
  )
  expect_equal(tau_result[["tau_total"]], matrix(0.05, nrow = 4L, ncol = 3L))

  active_setup <- list(
    priors            = priors,
    fit_priors        = list(
      mu  = BayesTools::prior("spike", list(0.20)),
      tau = BayesTools::prior("spike", list(0.05))
    ),
    is_PET            = FALSE,
    is_PEESE          = FALSE,
    is_weightfunction = FALSE
  )
  localized <- .iwmde_likelihood_posterior_samples(
    context      = list(indicator_names = character(), selection_spec = NULL),
    samples      = samples,
    active_setup = active_setup
  )
  expect_equal(localized[, "mu"], rep(0.20, nrow(samples)))
  expect_equal(localized[, "tau"], rep(0.05, nrow(samples)))

  data <- list(
    outcome = data.frame(yi = seq_len(3)),
    measure = "GEN"
  )
  attr(data, "outcome_type") <- "norm"
  row_parameters <- .iwmde_row_parameters(
    context      = list(data = data),
    row          = samples[1L, ],
    active_setup = active_setup
  )
  expect_equal(row_parameters[["mu"]], 0.20)
  expect_equal(row_parameters[["tau"]], 0.05)

  raw_fixed_samples <- matrix(
    c(9, 8),
    nrow     = 1L,
    dimnames = list(NULL, c("mu", "tau"))
  )
  raw_row_parameters <- .iwmde_row_parameters(
    context      = list(data = data),
    row          = raw_fixed_samples[1L, ],
    active_setup = active_setup
  )
  expect_equal(raw_row_parameters[["mu"]], 0.20)
  expect_equal(raw_row_parameters[["tau"]], 0.05)
  expect_equal(
    .iwmde_log_prior_row(raw_fixed_samples[1L, ], active_setup[["fit_priors"]]),
    .iwmde_log_prior_row(c(mu = 0.20, tau = 0.05), active_setup[["fit_priors"]])
  )
})


test_that("fixed PET and PEESE priors fill missing bias columns", {

  samples <- matrix(
    seq_len(4) / 10,
    ncol     = 1L,
    dimnames = list(NULL, "mu")
  )
  outcome_data <- data.frame(sei = c(0.10, 0.20, 0.30))

  pet_priors <- list(
    outcome = list(
      bias = BayesTools::prior_PET("spike", list(location = 0.30))
    )
  )
  pet_mu <- .evaluate.brma.mu(
    fit               = NULL,
    outcome_data      = outcome_data,
    mods_data         = NULL,
    mods_formula      = NULL,
    mods_priors       = NULL,
    is_mods           = FALSE,
    is_PET            = TRUE,
    is_PEESE          = FALSE,
    effect_direction  = "positive",
    bias_adjusted     = FALSE,
    K                 = 3L,
    posterior_samples = samples,
    priors            = pet_priors
  )
  expect_equal(pet_mu, matrix(samples[, "mu"], nrow = 4L, ncol = 3L) +
                 outer(rep(0.30, 4L), outcome_data[["sei"]]))

  data <- list(
    outcome = data.frame(yi = seq_len(3), sei = outcome_data[["sei"]]),
    measure = "GEN"
  )
  attr(data, "outcome_type") <- "norm"
  row_parameters <- .iwmde_row_parameters(
    context      = list(data = data),
    row          = samples[1L, ],
    active_setup = list(
      priors            = pet_priors,
      fit_priors        = list(PET = pet_priors[["outcome"]][["bias"]]),
      is_PET            = TRUE,
      is_PEESE          = FALSE,
      is_weightfunction = FALSE
    )
  )
  expect_equal(row_parameters[["PET"]], 0.30)

  raw_pet_row <- c(mu = 0.10, PET = 9)
  row_parameters <- .iwmde_row_parameters(
    context      = list(data = data),
    row          = raw_pet_row,
    active_setup = list(
      priors            = pet_priors,
      fit_priors        = list(PET = pet_priors[["outcome"]][["bias"]]),
      is_PET            = TRUE,
      is_PEESE          = FALSE,
      is_weightfunction = FALSE
    )
  )
  expect_equal(row_parameters[["PET"]], 0.30)

  peese_priors <- list(
    outcome = list(
      bias = BayesTools::prior_PEESE("spike", list(location = 0.40))
    )
  )
  peese_offset <- .evaluate.brma.bias_offset(
    fit               = NULL,
    outcome_data      = outcome_data,
    is_PET            = FALSE,
    is_PEESE          = TRUE,
    effect_direction  = "positive",
    K                 = 3L,
    posterior_samples = samples,
    priors            = peese_priors
  )
  expect_equal(peese_offset, outer(rep(0.40, 4L), outcome_data[["sei"]]^2))

  raw_peese_row <- c(mu = 0.10, PEESE = 9)
  row_parameters <- .iwmde_row_parameters(
    context      = list(data = data),
    row          = raw_peese_row,
    active_setup = list(
      priors            = peese_priors,
      fit_priors        = list(PEESE = peese_priors[["outcome"]][["bias"]]),
      is_PET            = FALSE,
      is_PEESE          = TRUE,
      is_weightfunction = FALSE
    )
  )
  expect_equal(row_parameters[["PEESE"]], 0.40)
})


test_that("fixed priors are constants, not IWMDE geometry dimensions", {

  samples <- matrix(
    seq_len(4) / 10,
    ncol     = 1L,
    dimnames = list(NULL, "mu")
  )
  context <- .iwmde_context_ensure_caches(list(
    posterior_samples = samples,
    flat_prior_list   = list(
      tau      = BayesTools::prior("spike", list(0.05)),
      mu_fixed = BayesTools::prior("spike", list(0.20))
    ),
    object = list(
      fit        = list(),
      likelihood = list(family = "normal"),
      data       = list(measure = "GEN")
    ),
    data  = list(measure = "GEN"),
    priors = list()
  ))

  fixed_spec <- .iwmde_parameter_spec(context, "mu_fixed")
  expect_equal(fixed_spec[["status"]], "ok")
  expect_equal(
    .iwmde_parameter_values(context, "mu_fixed", fixed_spec),
    rep(0.20, nrow(samples))
  )
  fixed_component <- .iwmde_parameter_components(context, "mu_fixed", fixed_spec)
  expect_false(any(fixed_component[["active"]]))
  expect_equal(fixed_component[["point_masses"]][["x"]], 0.20)

  linear_spec <- .iwmde_parameter_spec(
    context,
    "mu_plus_tau",
    list(type = "linear", weights = c(mu = 1, tau = 2))
  )
  expect_equal(linear_spec[["status"]], "ok")
  expect_equal(
    .iwmde_parameter_values(context, "mu_plus_tau", linear_spec),
    as.numeric(samples[, "mu"]) + 0.10
  )
  expect_equal(
    .iwmde_linear_value_row(context, samples[1L, ], linear_spec[["weights"]]),
    as.numeric(samples[1L, "mu"]) + 0.10
  )
  expect_equal(
    .iwmde_linear_active_columns(context, samples[1L, ], linear_spec[["weights"]]),
    "mu"
  )

  varying_samples <- matrix(
    c(seq_len(4) / 10, c(0.10, 0.20, 0.30, 0.40)),
    ncol     = 2L,
    dimnames = list(NULL, c("mu", "tau"))
  )
  conditioning_context <- context
  conditioning_context[["posterior_samples"]] <- varying_samples
  conditioning <- .iwmde_chen_conditioning_matrix(
    context        = conditioning_context,
    parameter      = "mu",
    parameter_spec = list(type = "primitive", parameter = "mu"),
    active_rows    = seq_len(nrow(varying_samples)),
    weight_rows    = seq_len(nrow(varying_samples))
  )
  expect_false("tau" %in% conditioning[["columns"]])
})


test_that("factor qCMDE/IWMDE densities attach only when all plotted columns resolve", {

  expect_true(.plot_brma_factor_density_complete(
    density_columns  = c("b", "a"),
    expected_columns = c("a", "b")
  ))
  expect_false(.plot_brma_factor_density_complete(
    density_columns  = "a",
    expected_columns = c("a", "b")
  ))
  expect_false(.plot_brma_factor_density_complete(
    density_columns  = character(),
    expected_columns = "a"
  ))
})


test_that("global IWMDE state excludes local latent parameters", {

  samples <- matrix(
    c(.10, .20, .40, .30, -.20),
    nrow     = 1L,
    dimnames = list(NULL, c("mu", "tau", "rho", "gamma[1]", "gamma[2]"))
  )
  data <- list(
    outcome = data.frame(
      yi      = c(.1, .2, .3),
      sei     = c(.1, .1, .1),
      cluster = c(1L, 1L, 2L)
    ),
    measure = "GEN"
  )
  attr(data, "outcome_type") <- "norm"
  attr(data, "cluster")      <- TRUE
  active_setup <- list(
    priors            = list(outcome = list()),
    fit_priors        = list(
      mu    = BayesTools::prior("normal", parameters = list(mean = 0, sd = 1)),
      tau   = BayesTools::prior("normal", parameters = list(mean = 0, sd = 1)),
      rho   = BayesTools::prior("beta", parameters = list(alpha = 1, beta = 1)),
      gamma = BayesTools::prior("normal", parameters = list(mean = 0, sd = 1))
    ),
    is_PET            = FALSE,
    is_PEESE          = FALSE,
    is_weightfunction = FALSE
  )
  context <- list(
    data            = data,
    flat_prior_list = active_setup[["fit_priors"]]
  )

  global_parameters <- .iwmde_row_parameters(
    context      = context,
    row          = samples[1L, ],
    active_setup = active_setup,
    state_scope  = "global"
  )
  local_parameters <- .iwmde_row_parameters(
    context      = context,
    row          = samples[1L, ],
    active_setup = active_setup,
    state_scope  = "local"
  )
  global_priors <- .iwmde_active_flat_prior_list(
    context     = context,
    row         = samples[1L, ],
    parameter   = "mu",
    state_scope = "global"
  )
  local_priors <- .iwmde_active_flat_prior_list(
    context     = context,
    row         = samples[1L, ],
    parameter   = "gamma[1]",
    state_scope = "local"
  )

  expect_null(global_parameters[["gamma"]])
  expect_equal(local_parameters[["gamma"]], c(.30, -.20))
  expect_false("gamma" %in% names(global_priors))
  expect_true("gamma" %in% names(local_priors))
})


test_that("known-V marginalized random SDs stay in global IWMDE state", {

  fit     <- .load_raw_fit_or_skip("brma.mv_block_mvn_random")
  context <- .iwmde_context(fit)
  term    <- .data_marginalized_random_effects(fit[["data"]])[[1L]]
  row     <- context[["posterior_samples"]][1L, ]
  setup   <- .iwmde_active_setup(context, row)

  parameters <- .iwmde_row_parameters(
    context      = context,
    row          = row,
    active_setup = setup,
    state_scope  = "global"
  )

  expect_true(term[["sd_parameter_names"]] %in% names(parameters))
  expect_silent(.iwmde_log_likelihood_parameters_marginal_internal(
    context      = context,
    parameters   = parameters,
    active_setup = setup
  ))

  parameter <- .check_and_select_plot_parameter(
    parameter       = "mu",
    parameter_mods  = NULL,
    parameter_scale = NULL,
    component       = NULL,
    object          = fit
  )
  sample_parameter <- .as_mixed_posteriors_parameters(fit, parameter)
  samples <- BayesTools::as_mixed_posteriors(
    model            = fit[["fit"]],
    parameters       = sample_parameter,
    conditional      = NULL,
    transform_scaled = TRUE
  )
  density_sample_parameter <- .plot_brma_density_sample_parameter(
    samples          = samples,
    parameter        = parameter,
    sample_parameter = sample_parameter
  )
  samples <- .plot_brma_attach_iwmde(
    object               = fit,
    samples              = samples,
    parameter            = parameter,
    sample_parameter     = density_sample_parameter,
    conditional          = NULL,
    n_points             = 20,
    max_samples          = 20,
    normalization_points = 20,
    normalization_prob   = .999,
    density_method       = "qCMDE",
    display_grid         = "adaptive"
  )
  diagnostic <- attr(samples, "iwmde_diagnostics")[[1L]]
  diagnostic_reason <- if (is.null(diagnostic[["reason"]])) {
    ""
  } else {
    diagnostic[["reason"]]
  }
  density_failure <- .iwmde_diagnostics_density_failure_reason(
    diagnostic[["diagnostics"]]
  )
  density_failure <- if (is.null(density_failure)) {
    ""
  } else {
    density_failure
  }

  expect_equal(diagnostic[["status"]], "ok")
  expect_false(grepl("Missing posterior column", diagnostic_reason))
  expect_false(grepl("Missing posterior column", density_failure))
})


test_that("known-V marginalized allocation weights stay in global IWMDE state", {

  fit     <- .load_raw_fit_or_skip("brma.mv_block_mvn_random_scale")
  context <- .iwmde_context(fit)
  row     <- context[["posterior_samples"]][1L, ]
  setup   <- .iwmde_active_setup(context, row)

  parameters <- .iwmde_row_parameters(
    context      = context,
    row          = row,
    active_setup = setup,
    state_scope  = "global"
  )

  expect_true(any(grepl("__xRE_ALLOCx_", names(parameters), fixed = TRUE)))
  expect_silent(.iwmde_log_likelihood_parameters_marginal_internal(
    context      = context,
    parameters   = parameters,
    active_setup = setup
  ))
})


test_that("density diagnostics gate unstable qCMDE/IWMDE attributes", {

  diagnostics <- list(
    estimator                         = "q_grid_cmde",
    plot_scale_relative_mcse          = .01,
    bulk_max_relative_mcse            = .01,
    bulk_min_ess                      = 80,
    bulk_max_weight_share             = .05,
    max_relative_mcse                 = .01,
    min_finite_terms                  = 80L,
    min_ess                           = 80,
    max_weight_share                  = .05,
    active_mass                       = 1,
    final_normalization_integral      = 1,
    normalization_relative_error      = 0,
    normalization_mass_ratio          = 1,
    max_ordinate_relative_change      = .06,
    max_normalizer_relative_change    = .06,
    max_quadrature_relative_change    = NA_real_
  )
  diagnostic <- list(
    status       = "ok",
    iwmde        = list(
      x         = c(0, 1),
      y         = c(1, 1),
      estimator = "q_grid_cmde"
    ),
    diagnostics  = diagnostics,
    point_masses = NULL
  )

  expect_match(
    .iwmde_diagnostics_density_failure_reason(diagnostics),
    "qCMDE.*ordinate"
  )
  expect_null(.iwmde_posterior_density_attribute(
    diagnostic     = diagnostic,
    density_method = "qCMDE"
  ))

  diagnostics[["max_ordinate_relative_change"]] <- .01
  diagnostics[["max_normalizer_relative_change"]] <- .01
  diagnostics[["n_estimator_rows"]] <- 80
  diagnostics[["n_active_state_keys"]] <- 1
  expect_null(.iwmde_diagnostics_density_failure_reason(diagnostics))
  expect_true(any(grepl(
    "estimator rows",
    .iwmde_diagnostics_density_warning(diagnostics)
  )))
  diagnostics[["n_active_state_keys"]] <- 2
  expect_match(
    .iwmde_diagnostics_density_failure_reason(diagnostics),
    "model-averaged density"
  )
  diagnostics[["n_estimator_rows"]] <- 600
  diagnostics[["n_active_state_keys"]] <- 2
  diagnostics[["bulk_min_ess"]] <- 120

  diagnostics[["max_quadrature_relative_change"]] <- .03
  expect_null(.iwmde_diagnostics_density_failure_reason(diagnostics))
  expect_true(any(grepl(
    "quadrature",
    .iwmde_diagnostics_density_warning(diagnostics)
  )))

  diagnostics[["max_quadrature_relative_change"]] <- .06
  expect_match(
    .iwmde_diagnostics_density_failure_reason(diagnostics),
    "quadrature"
  )
})


test_that("density gates use bulk and 5/95 tail checkpoints", {

  density <- list(
    x                = c(-4, -1.645, 0, 1.645, 4),
    y                = c(.001, .10, 1, .10, .001),
    mcse             = c(.03, .008, .01, .008, .03),
    relative_mcse    = c(30, .08, .01, .08, 30),
    sampling_mcse    = c(.03, .008, .01, .008, .03),
    sampling_relative_mcse = c(30, .08, .01, .08, 30),
    ess              = c(2, 120, 400, 120, 2),
    max_weight_share = c(.70, .04, .01, .04, .70)
  )
  curve <- .iwmde_density_curve_diagnostics(
    density            = density,
    tail_probabilities = c(.05, .95),
    tail_values        = c(-1.645, 1.645)
  )

  expect_equal(.iwmde_density_tail_probabilities(), c(.05, .95))
  expect_equal(curve[["plot_scale_relative_mcse"]], .03)
  expect_equal(curve[["plot_scale_sampling_relative_mcse"]], .03)
  expect_equal(curve[["bulk_max_relative_mcse"]], .08)
  expect_equal(curve[["bulk_max_sampling_relative_mcse"]], .08)
  expect_equal(curve[["bulk_min_ess"]], 120)
  expect_equal(curve[["bulk_max_weight_share"]], .04)
  expect_equal(curve[["tail_evaluation_x"]], c(-1.645, 1.645))
  expect_equal(curve[["tail_relative_mcse"]], c(.08, .08))
  expect_equal(curve[["tail_ess"]], c(120, 120))
  expect_equal(curve[["tail_max_weight_share"]], c(.04, .04))

  diagnostics <- c(list(
    estimator                       = "q_grid_cmde",
    max_relative_mcse               = 30,
    max_sampling_relative_mcse      = 30,
    min_finite_terms                = 500L,
    min_ess                         = 2,
    max_weight_share                = .70,
    n_estimator_rows                = 500L,
    n_active_state_keys             = 1L,
    active_mass                     = 1,
    final_normalization_integral    = 1,
    normalization_relative_error    = 0,
    normalization_mass_ratio        = 1,
    max_ordinate_relative_change    = 0,
    max_normalizer_relative_change  = 0,
    max_quadrature_relative_change  = NA_real_
  ), curve)

  expect_null(.iwmde_diagnostics_density_failure_reason(diagnostics))
  expect_length(.iwmde_diagnostics_density_warning(diagnostics), 0L)

  diagnostics[["plot_scale_sampling_relative_mcse"]] <- .12
  expect_match(
    .iwmde_diagnostics_density_warning(diagnostics),
    "plot-scale finite-population sampling MCSE"
  )
  diagnostics[["plot_scale_sampling_relative_mcse"]] <- .30
  expect_match(
    .iwmde_diagnostics_density_failure_reason(diagnostics),
    "plot-scale finite-population sampling MCSE"
  )
  diagnostics[["plot_scale_sampling_relative_mcse"]] <- .03

  diagnostics[["bulk_max_sampling_relative_mcse"]] <- .12
  expect_match(
    .iwmde_diagnostics_density_warning(diagnostics),
    "bulk finite-population sampling relative MCSE"
  )
  diagnostics[["bulk_max_sampling_relative_mcse"]] <- .30
  expect_match(
    .iwmde_diagnostics_density_failure_reason(diagnostics),
    "bulk finite-population sampling relative MCSE"
  )
  diagnostics[["bulk_max_sampling_relative_mcse"]] <- .08

  diagnostics[["plot_scale_relative_mcse"]] <- .12
  expect_match(
    .iwmde_diagnostics_density_warning(diagnostics),
    "plot-scale MCSE"
  )
  diagnostics[["plot_scale_relative_mcse"]] <- .30
  expect_match(
    .iwmde_diagnostics_density_failure_reason(diagnostics),
    "plot-scale MCSE"
  )
  diagnostics[["plot_scale_relative_mcse"]] <- .03

  diagnostics[["bulk_max_relative_mcse"]] <- .12
  expect_match(
    .iwmde_diagnostics_density_warning(diagnostics),
    "bulk relative MCSE"
  )
  diagnostics[["bulk_max_relative_mcse"]] <- .30
  expect_match(
    .iwmde_diagnostics_density_failure_reason(diagnostics),
    "bulk relative MCSE"
  )
  diagnostics[["bulk_max_relative_mcse"]] <- .08

  diagnostics[["bulk_min_ess"]] <- 80
  expect_match(
    .iwmde_diagnostics_density_warning(diagnostics),
    "bulk effective sample size"
  )
  diagnostics[["bulk_min_ess"]] <- 40
  expect_match(
    .iwmde_diagnostics_density_failure_reason(diagnostics),
    "bulk effective sample size"
  )
  diagnostics[["bulk_min_ess"]] <- 120

  diagnostics[["bulk_max_weight_share"]] <- .12
  expect_match(
    .iwmde_diagnostics_density_warning(diagnostics),
    "bulk density importance weight"
  )
  diagnostics[["bulk_max_weight_share"]] <- .30
  expect_match(
    .iwmde_diagnostics_density_failure_reason(diagnostics),
    "bulk density importance weight"
  )
})


test_that("primitive AST symbols match qCMDE/IWMDE warning records", {

  record <- data.frame(
    label     = "mu_intercept",
    parameter = "mu_intercept",
    value     = 0,
    stringsAsFactors = FALSE
  )
  side <- list(
    type  = "point",
    expr  = as.name("mu_intercept"),
    value = 0
  )

  expect_true(.hypothesis_brma_warning_side_matches(side, record))
  side[["expr"]] <- as.name("tau")
  expect_false(.hypothesis_brma_warning_side_matches(side, record))
})


test_that("IWMDE BF diagnostics retain the exact boundary ordinate", {

  evaluation_value  <- 0
  expected_ordinate <- stats::dnorm(evaluation_value)
  density <- list(
    x                = 0,
    evaluation_x     = evaluation_value,
    y                = expected_ordinate,
    mcse             = .1,
    relative_mcse    = .05,
    finite_terms     = 40L,
    ess              = 25,
    max_weight_share = .2,
    max_log_ratio    = 1,
    pilot_y          = 2,
    ordinate_relative_change = 0,
    ordinate_log_change      = 0
  )

  bf_diagnostics <- .iwmde_density_bf_diagnostics(density, 0)
  diagnostic <- list(
    status      = "ok",
    diagnostics = c(
      bf_diagnostics,
      list(
        estimator              = "q_grid_cmde",
        weight_method          = "conditional_grid",
        active_mass            = 1,
        final_normalization_integral = 1,
        normalization_relative_error = 0,
        normalization_mass_ratio = 1,
        max_ordinate_relative_change = 0,
        max_normalizer_relative_change = 0,
        normalization_range    = c(evaluation_value / 2, evaluation_value * 2)
      )
    )
  )
  ordinate <- .iwmde_posterior_ordinate_attribute(diagnostic, "qCMDE")

  expect_true(bf_diagnostics[["bf_included"]])
  expect_equal(bf_diagnostics[["bf_value"]], 0)
  expect_equal(bf_diagnostics[["bf_evaluation_value"]], evaluation_value)
  expect_equal(ordinate[["value"]], 0)
  expect_equal(ordinate[["evaluation_value"]], evaluation_value)
  expect_equal(ordinate[["ordinate"]], expected_ordinate)
  expect_equal(
    ordinate[["iwmde_provenance"]][["result"]][["evaluation_value"]],
    evaluation_value
  )
  expect_true(.iwmde_posterior_ordinate_supports_bf(ordinate))
})


test_that("IWMDE delegates exact prior ordinates to BayesTools", {

  beta_zero <- BayesTools::prior(
    "beta",
    parameters = list(alpha = 2, beta = 2)
  )
  beta_infinite <- BayesTools::prior(
    "beta",
    parameters = list(alpha = .5, beta = 2)
  )
  beta_regular <- BayesTools::prior(
    "beta",
    parameters = list(alpha = 1, beta = 1)
  )
  normal_regular <- BayesTools::prior(
    "normal",
    parameters = list(mean = 0, sd = 1)
  )

  behavior <- function(prior, value) {

    .iwmde_prior_ordinate_classifications(prior, value)[[1L]][["behavior"]]
  }

  expect_identical(behavior(beta_zero, 0), "zero")
  expect_identical(behavior(beta_infinite, 0), "infinite")
  expect_identical(behavior(beta_regular, 0), "regular")
  expect_identical(behavior(normal_regular, 40), "regular")
  expect_identical(behavior(normal_regular, 1e200), "regular")
  expect_identical(behavior(
    BayesTools::prior("lognormal", parameters = list(meanlog = 0, sdlog = 1)),
    0
  ), "zero")
  expect_identical(behavior(
    BayesTools::prior("gamma", parameters = list(shape = .5, rate = 1)),
    0
  ), "infinite")
  expect_identical(behavior(
    BayesTools::prior(
      "moment",
      parameters = list(location = 0, tau = 1, order = 1)
    ),
    0
  ), "zero")
  expect_identical(behavior(
    BayesTools::prior("point", parameters = list(location = 0)),
    0
  ), "point_mass")
  expect_identical(
    .iwmde_prior_ordinate_classifications(NULL, 0)[[1L]][["behavior"]],
    "unknown"
  )
  unknown <- .iwmde_prior_ordinate_classifications(NULL, 0)[[1L]]
  expect_identical(unknown[["method"]], "unsupported_provenance")
  expect_equal(unknown[["point_mass"]], 0)

  induced <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(beta = normal_regular),
    weights    = c(beta = 2),
    n_grid     = 128
  )
  induced_ordinate <- .iwmde_prior_ordinate_classifications(induced, 0)[[1L]]
  expect_identical(induced_ordinate[["behavior"]], "regular")
  expect_identical(induced_ordinate[["method"]], "scalar_affine")

  adjacent <- .iwmde_prior_ordinate_classifications(
    normal_regular,
    c(1, 1 + .Machine$double.eps)
  )
  expect_length(adjacent, 2L)
})


test_that("IWMDE prepares compact prior ordinates idempotently", {

  prior <- BayesTools::prior(
    "beta",
    parameters = list(alpha = 2, beta = 2)
  )
  context <- list(flat_prior_list = list(rho = prior))

  primitive <- .iwmde_prepare_prior_ordinates(
    context        = context,
    parameter      = "rho",
    parameter_spec = list(type = "primitive"),
    values         = 0
  )
  expect_identical(
    primitive[["prior_ordinates"]][[1L]][["behavior"]],
    "zero"
  )

  supplied <- .iwmde_prepare_prior_ordinates(
    context        = list(flat_prior_list = list()),
    parameter      = "contrast",
    parameter_spec = list(type = "linear", prior_density = prior),
    values         = 0
  )
  expect_null(supplied[["prior_density"]])
  expect_identical(supplied[["prior_ordinates"]][[1L]][["behavior"]], "zero")

  conditional <- .iwmde_prepare_prior_ordinates(
    context        = context,
    parameter      = "rho",
    parameter_spec = list(type = "primitive", conditional = "rho"),
    values         = 0
  )
  expect_identical(
    conditional[["prior_ordinates"]][[1L]][["behavior"]],
    "unknown"
  )

  n_classifications <- 0L
  testthat::local_mocked_bindings(
    .iwmde_prior_ordinate_classifications = function(...) {
      n_classifications <<- n_classifications + 1L
      ordinate <- list(.iwmde_unknown_prior_ordinate(0))
      names(ordinate) <- .iwmde_key_number(0)
      ordinate
    },
    .package = "RoBMA"
  )
  prepared <- .iwmde_prepare_prior_ordinates(
    context        = context,
    parameter      = "rho",
    parameter_spec = list(type = "primitive"),
    values         = 0
  )
  .iwmde_prepare_prior_ordinates(
    context        = context,
    parameter      = "rho",
    parameter_spec = prepared,
    values         = 0
  )
  expect_equal(n_classifications, 1L)

  stale <- primitive
  stale[["prior_ordinates"]][[1L]][["value"]] <- 1
  expect_error(
    .iwmde_prepare_prior_ordinates(
      context        = context,
      parameter      = "rho",
      parameter_spec = stale,
      values         = 0
    ),
    "do not match the current request"
  )
})


test_that("IWMDE warns about nonregular priors before estimator failure", {

  ordinates <- .iwmde_prior_ordinate_classifications(
    BayesTools::prior("beta", parameters = list(alpha = 2, beta = 2)),
    0
  )
  warnings <- .iwmde_ordinate_prior_warnings("rho", ordinates)
  expect_length(warnings, 1L)
  expect_match(warnings, "exact requested value is retained")
  expect_match(
    .iwmde_diagnostics_bf_warning(list(
      estimator         = "iwmde",
      ordinate_warnings = warnings
    )),
    "prior density tends to zero"
  )

  expect_match(
    .iwmde_ordinate_prior_warnings(
      "rho",
      .iwmde_prior_ordinate_classifications(
        BayesTools::prior("beta", parameters = list(alpha = .5, beta = 2)),
        0
      )
    ),
    "singular (tends to infinity)",
    fixed = TRUE
  )
  expect_match(
    .iwmde_ordinate_prior_warnings(
      "rho",
      .iwmde_prior_ordinate_classifications(
        BayesTools::prior("point", parameters = list(location = 0)),
        0
      )
    ),
    "positive point mass"
  )
  expect_match(
    .iwmde_ordinate_prior_warnings(
      "rho",
      .iwmde_prior_ordinate_classifications(NULL, 0)
    ),
    "deterministic provenance"
  )
  expect_false(grepl(
    "is nonregular",
    .iwmde_ordinate_prior_warnings(
      "rho",
      .iwmde_prior_ordinate_classifications(NULL, 0)
    ),
    fixed = TRUE
  ))

  events <- character()
  testthat::local_mocked_bindings(
    .iwmde_context_ensure_caches = function(context) context,
    .iwmde_check_context_density_method_supported = function(...) {
      events <<- c(events, "error")
      stop("estimator failed", call. = FALSE)
    },
    .package = "RoBMA"
  )
  expect_error(
    withCallingHandlers(
      .iwmde_estimate(
        context         = list(flat_prior_list = list()),
        parameter       = "rho",
        density_method  = "qCMDE",
        density_control = list(),
        outputs         = "ordinate",
        values          = 0,
        parameter_spec  = list(
          type          = "primitive",
          prior_density = BayesTools::prior(
            "beta",
            parameters = list(alpha = 2, beta = 2)
          )
        )
      ),
      warning = function(warning) {
        events <<- c(events, "warning")
        invokeRestart("muffleWarning")
      }
    ),
    "estimator failed"
  )
  expect_identical(events, c("warning", "error"))
})


test_that("IWMDE preflights only ordinate outputs", {

  captured_spec <- NULL
  testthat::local_mocked_bindings(
    .iwmde_context_ensure_caches = function(context) context,
    .iwmde_check_context_density_method_supported = function(...) NULL,
    .iwmde_plan = function(..., parameter_spec) {
      captured_spec <<- parameter_spec
      list(prior_ordinates = list(), ordinate_warnings = character())
    },
    .iwmde_estimate_from_plan = function(context, plan, cache) {
      list(plan = plan)
    },
    .package = "RoBMA"
  )

  expect_no_warning(.iwmde_estimate(
    context         = list(flat_prior_list = list()),
    parameter       = "rho",
    density_method  = "qCMDE",
    density_control = list(),
    outputs         = "density",
    values          = 0,
    parameter_spec  = list(
      type          = "primitive",
      prior_density = BayesTools::prior(
        "beta",
        parameters = list(alpha = 2, beta = 2)
      )
    )
  ))
  expect_length(captured_spec[["prior_ordinates"]], 0L)

  expect_warning(
    .iwmde_estimate(
      context         = list(flat_prior_list = list()),
      parameter       = "rho",
      density_method  = "qCMDE",
      density_control = list(),
      outputs         = c("density", "ordinate"),
      values          = c(0, 1),
      parameter_spec  = list(
        type          = "primitive",
        prior_density = BayesTools::prior(
          "beta",
          parameters = list(alpha = 2, beta = 2)
        )
      )
    ),
    "prior density tends to zero"
  )
  expect_length(captured_spec[["prior_ordinates"]], 1L)
  expect_identical(captured_spec[["prior_ordinates"]][[1L]][["value"]], 0)
})


test_that("IWMDE keeps adjacent representable coordinates distinct", {

  value          <- 1
  adjacent_value <- value + .Machine$double.eps
  ordinate       <- list(value = value)
  density        <- list(
    x           = c(value, adjacent_value),
    y           = c(2, 3),
    finite_terms = c(10, 20),
    ess           = c(5, 10)
  )

  expect_true(.iwmde_ordinate_value_matches(ordinate, value))
  expect_false(.iwmde_ordinate_value_matches(ordinate, adjacent_value))

  diagnostics <- .iwmde_density_bf_diagnostics(density, adjacent_value)
  expect_true(diagnostics[["bf_included"]])
  expect_equal(diagnostics[["bf_grid_index"]], 2L)
  expect_equal(diagnostics[["bf_ordinate"]], 3)

  grid <- .iwmde_include_display_values(
    grid    = c(0, value, 2),
    values  = adjacent_value,
    xlim    = c(0, 2),
    support = c(0, 2)
  )
  expect_true(any(grid == adjacent_value))
  expect_false(any(grid == value))

  clustered_values <- c(.Machine$double.eps, 2 * .Machine$double.eps)
  clustered_grid <- .iwmde_include_display_values(
    grid    = c(0, 1, 2),
    values  = clustered_values,
    xlim    = c(0, 2),
    support = c(0, 2)
  )
  expect_length(clustered_grid, 3L)
  expect_true(all(clustered_values %in% clustered_grid))
})


test_that("qCMDE/IWMDE posterior attributes carry RoBMA provenance", {

  density_control <- list(
    n_points             = 20,
    max_samples          = 30,
    normalization_points = 40,
    normalization_prob   = .95,
    display_grid         = "adaptive"
  )
  metadata <- list(parameter = "mu", conditional_rule = "AND")
  diagnostic <- list(
    status       = "ok",
    target_key   = "primitive|mu",
    active_mass  = 1,
    point_masses = data.frame(x = numeric(), mass = numeric()),
    iwmde        = list(
      x                       = c(-1, 0, 1),
      y                       = c(.1, .4, .1),
      estimator               = "q_grid_cmde",
      final_normalization_integral = 1,
      normalization_mass_ratio = 1,
      ordinate_relative_change = c(0, 0, 0)
    ),
    diagnostics = list(
      bf_included         = TRUE,
      bf_value            = 0,
      bf_evaluation_value = 0,
      bf_ordinate         = .4,
      bf_mcse             = .01,
      bf_relative_mcse    = .01,
      bf_error_percent    = 1,
      bf_finite_terms     = 500,
      bf_ess              = 120,
      bf_max_weight_share = .05,
      bf_max_log_ratio    = 0,
      n_estimator_rows    = 500,
      plot_scale_relative_mcse    = .01,
      bulk_probability_range      = c(.05, .95),
      bulk_x_range                = c(-1, 1),
      bulk_max_relative_mcse      = .01,
      bulk_min_ess                = 120,
      bulk_max_weight_share       = .05,
      tail_probabilities          = c(.05, .95),
      tail_target_x               = c(-1, 1),
      tail_evaluation_x           = c(-1, 1),
      tail_relative_mcse          = c(.01, .01),
      tail_ess                    = c(120, 120),
      tail_max_weight_share       = c(.05, .05),
      active_mass         = 1,
      final_normalization_integral = 1,
      normalization_relative_error = 0,
      bf_ordinate_relative_change = 0,
      max_ordinate_relative_change = 0,
      max_normalizer_relative_change = 0,
      normalization_range = c(-1, 1),
      estimator           = "q_grid_cmde",
      weight_method       = "mock"
    )
  )

  density_attr <- .iwmde_posterior_density_attribute(
    diagnostic      = diagnostic,
    density_method  = "qCMDE",
    metadata        = metadata,
    density_control = density_control
  )
  ordinate_attr <- .iwmde_posterior_ordinate_attribute(
    diagnostic      = diagnostic,
    density_method  = "qCMDE",
    metadata        = metadata,
    density_control = density_control
  )
  expected_provenance <- .iwmde_provenance_request(
    density_method   = "qCMDE",
    method           = "q_grid_cmde",
    metadata         = metadata,
    density_control  = density_control,
    value            = 0,
    attribute        = "ordinate",
    target_key       = "primitive|mu"
  )
  expected_density_provenance <- .iwmde_provenance_request(
    density_method   = "qCMDE",
    method           = "q_grid_cmde",
    metadata         = metadata,
    density_control  = density_control,
    attribute        = "density",
    target_key       = "primitive|mu"
  )

  expect_true("iwmde_provenance" %in% names(density_attr))
  expect_true("iwmde_provenance" %in% names(ordinate_attr))
  expect_true(.iwmde_posterior_density_matches_request(
    posterior_density = density_attr,
    provenance        = expected_density_provenance
  ))
  no_provenance_density <- density_attr
  no_provenance_density[["iwmde_provenance"]] <- NULL
  expect_false(.iwmde_posterior_density_matches_request(
    posterior_density = no_provenance_density,
    provenance        = expected_density_provenance
  ))
  expect_false(.iwmde_posterior_density_matches_request(
    posterior_density = density_attr,
    provenance        = .iwmde_provenance_request(
      density_method   = "qCMDE",
      method           = "q_grid_cmde",
      metadata         = metadata,
      density_control  = modifyList(density_control, list(n_points = 25L)),
      attribute        = "density",
      target_key       = "primitive|mu"
    )
  ))
  expect_false(.iwmde_posterior_density_matches_request(
    posterior_density = density_attr,
    provenance        = .iwmde_provenance_request(
      density_method   = "qCMDE",
      method           = "q_grid_cmde",
      metadata         = metadata,
      density_control  = density_control,
      attribute        = "density",
      target_key       = "primitive|tau"
    )
  ))

  provenance <- ordinate_attr[["iwmde_provenance"]]
  expect_equal(provenance[["schema_version"]], "4")
  expect_equal(provenance[["algorithm_version"]], "11")
  expect_equal(provenance[["provenance_level"]], "diagnostic_adapter")
  expect_equal(provenance[["density_method"]], "qCMDE")
  expect_equal(provenance[["internal_method"]], "q_grid_cmde")
  expect_equal(provenance[["target"]][["parameter"]], "mu")
  expect_equal(provenance[["target"]][["conditional_rule"]], "AND")
  expect_equal(provenance[["requested_value"]], .iwmde_key_number(0))
  expect_equal(
    provenance[["result"]][["evaluation_value"]],
    0
  )
  expect_equal(
    provenance[["result"]][["evaluation_value_key"]],
    .iwmde_key_number(0)
  )
  expect_true(isTRUE(provenance[["result"]][["bf_grade"]]))
  expect_true(.iwmde_posterior_ordinate_matches_request(
    posterior_ordinate = ordinate_attr,
    value              = 0,
    provenance         = expected_provenance
  ))
  expect_false(.iwmde_posterior_ordinate_matches_request(
    posterior_ordinate = ordinate_attr,
    value              = 0,
    provenance         = .iwmde_provenance_request(
      density_method   = "IWMDE",
      method           = "iwmde",
      metadata         = metadata,
      density_control  = density_control,
      value            = 0,
      attribute        = "ordinate",
      target_key       = "primitive|mu"
    )
  ))
})


test_that("iwmde_estimate returns plan-backed attributes and caches by provenance", {

  density_control <- list(
    n_points             = 20,
    max_samples          = 20,
    normalization_points = 20,
    normalization_prob   = .95,
    display_grid         = "adaptive"
  )
  prior_density <- BayesTools::prior(
    "normal",
    parameters = list(mean = 0, sd = 1)
  )
  context <- list(
    posterior_samples = matrix(
      stats::rnorm(60),
      ncol     = 1L,
      dimnames = list(NULL, "mu")
    ),
    object = list(
      fit        = list(),
      likelihood = list(family = "normal"),
      data       = list(measure = "yi")
    )
  )
  mock_diagnostic <- function(parameter, target_key) {

    list(
      parameter    = parameter,
      status       = "ok",
      target_key   = target_key,
      active_mass  = 1,
      point_masses = data.frame(x = numeric(), mass = numeric()),
      iwmde        = list(
        x                         = c(-1, 0, 1),
        y                         = c(.1, .4, .1),
        estimator                 = "q_grid_cmde",
        final_normalization_integral = 1,
        normalization_mass_ratio  = 1,
        ordinate_relative_change  = c(0, 0, 0)
      ),
      diagnostics = list(
        bf_included         = TRUE,
        bf_value            = 0,
        bf_evaluation_value = 0,
        bf_ordinate         = .4,
        bf_mcse             = .01,
        bf_relative_mcse    = .01,
        bf_error_percent    = 1,
        bf_finite_terms     = 500,
        bf_ess              = 120,
        bf_max_weight_share = .05,
        bf_max_log_ratio    = 0,
        n_estimator_rows    = 500,
        plot_scale_relative_mcse    = .01,
        bulk_probability_range      = c(.05, .95),
        bulk_x_range                = c(-1, 1),
        bulk_max_relative_mcse      = .01,
        bulk_min_ess                = 120,
        bulk_max_weight_share       = .05,
        tail_probabilities          = c(.05, .95),
        tail_target_x               = c(-1, 1),
        tail_evaluation_x           = c(-1, 1),
        tail_relative_mcse          = c(.01, .01),
        tail_ess                    = c(120, 120),
        tail_max_weight_share       = c(.05, .05),
        active_mass         = 1,
        final_normalization_integral = 1,
        normalization_relative_error = 0,
        bf_ordinate_relative_change = 0,
        max_ordinate_relative_change = 0,
        max_normalizer_relative_change = 0,
        normalization_range = c(-1, 1),
        estimator           = "q_grid_cmde",
        weight_method       = "mock"
      )
    )
  }
  density_calls  <- 0L
  ordinate_calls <- 0L

  testthat::local_mocked_bindings(
    .iwmde_row_states = function(context, rows, parameter, parameter_spec,
                                 estimator = NULL) {

      lapply(rows, function(row) {
        .iwmde_new_row_state(list(
          baseline_log_q = 0,
          active_key     = "all",
          row_index      = row
        ))
      })
    },
    .iwmde_execute_plan_diagnostic = function(context, plan, output,
                                              execution_cache = NULL,
                                              diagnostic_cache = NULL) {

      if (identical(output, "density")) {
        density_calls <<- density_calls + 1L
      } else {
        ordinate_calls <<- ordinate_calls + 1L
      }
      mock_diagnostic(
        parameter  = plan[["target"]][["parameter"]],
        target_key = plan[["target"]][["target_key"]]
      )
    },
    .package = "RoBMA"
  )

  cache <- .iwmde_estimate_cache()
  estimate <- .iwmde_estimate(
    context         = context,
    parameter       = "mu",
    density_method  = "qCMDE",
    density_control = density_control,
    outputs         = c("density", "ordinate"),
    values          = 0,
    parameter_spec  = list(type = "primitive", prior_density = prior_density),
    metadata        = list(parameter = "mu"),
    cache           = cache
  )
  estimate_again <- .iwmde_estimate(
    context         = context,
    parameter       = "mu",
    density_method  = "qCMDE",
    density_control = density_control,
    outputs         = c("density", "ordinate"),
    values          = 0,
    parameter_spec  = list(type = "primitive", prior_density = prior_density),
    metadata        = list(parameter = "mu"),
    cache           = cache
  )

  expect_s3_class(estimate, "iwmde_estimate")
  expect_s3_class(estimate[["plan"]], "iwmde_plan")
  expect_null(estimate[["plan"]][["parameter_spec"]][["prior_density"]])
  expect_identical(
    estimate[["plan"]][["prior_ordinates"]][[1L]][["behavior"]],
    "regular"
  )
  expect_equal(estimate_again[["plan"]][["plan_key"]], estimate[["plan"]][["plan_key"]])
  expect_equal(density_calls, 1L)
  expect_equal(ordinate_calls, 1L)
  expect_true("iwmde_provenance" %in% names(estimate[["posterior_density"]]))
  expect_true("iwmde_provenance" %in% names(estimate[["posterior_ordinate"]]))

  provenance <- estimate[["posterior_ordinate"]][["iwmde_provenance"]]
  expect_equal(provenance[["provenance_level"]], "iwmde_plan")
  expect_equal(provenance[["plan_key"]], estimate[["plan"]][["plan_key"]])
  expect_identical(provenance[["prior_ordinates"]],
                   estimate[["plan"]][["prior_ordinates"]])
  expect_true(is.list(provenance[["source_fingerprint"]]))
  expect_true(.iwmde_posterior_ordinate_matches_request(
    posterior_ordinate = estimate[["posterior_ordinate"]],
    value              = 0,
    provenance         = .iwmde_request_provenance(
      context         = context,
      parameter       = "mu",
      density_method  = "qCMDE",
      density_control = density_control,
      attribute       = "ordinate",
      value           = 0,
      parameter_spec  = list(type = "primitive", prior_density = prior_density),
      metadata        = list(parameter = "mu")
    )
  ))
  density_request <- .iwmde_request_provenance(
    context         = context,
    parameter       = "mu",
    density_method  = "qCMDE",
    density_control = density_control,
    attribute       = "density",
    parameter_spec  = list(type = "primitive", prior_density = prior_density),
    metadata        = list(parameter = "mu")
  )
  expect_true(.iwmde_posterior_density_matches_request(
    posterior_density = estimate[["posterior_density"]],
    provenance        = density_request
  ))
  expect_null(
    estimate[["posterior_density"]][["iwmde_provenance"]][["prior_ordinates"]]
  )
  expect_identical(
    density_request[["request_key"]],
    .iwmde_request_provenance(
      context         = context,
      parameter       = "mu",
      density_method  = "qCMDE",
      density_control = density_control,
      attribute       = "density",
      parameter_spec  = list(
        type          = "primitive",
        prior_density = BayesTools::prior(
          "normal",
          parameters = list(mean = 0, sd = 2)
        )
      ),
      metadata        = list(parameter = "mu")
    )[["request_key"]]
  )

  context_changed <- context
  context_changed[["posterior_samples"]][1, 1] <- 99
  expect_false(.iwmde_posterior_ordinate_matches_request(
    posterior_ordinate = estimate[["posterior_ordinate"]],
    value              = 0,
    provenance         = .iwmde_request_provenance(
      context         = context_changed,
      parameter       = "mu",
      density_method  = "qCMDE",
      density_control = density_control,
      attribute       = "ordinate",
      value           = 0,
      parameter_spec  = list(type = "primitive", prior_density = prior_density),
      metadata        = list(parameter = "mu")
    )
  ))
})


test_that("qCMDE ordinate and IWMDE mass thresholds warn before failing", {

  ordinate <- function(estimator, normalization_integral,
                       ordinate_relative_change = NA_real_) {

    BayesTools::posterior_ordinate_attribute(
      value          = 0,
      ordinate       = .4,
      method         = estimator,
      density_method = if (identical(estimator, "q_grid_cmde")) {
        "qCMDE"
      } else {
        "IWMDE"
      },
      diagnostics    = list(
        evaluation_value       = 0,
        relative_mcse          = .01,
        finite_terms           = 200,
        ess                    = 150,
        max_weight_share       = .05,
        active_mass            = 1,
        normalization_relative_error = abs(normalization_integral - 1),
        ordinate_relative_change = ordinate_relative_change,
        max_normalizer_relative_change = ordinate_relative_change,
        normalization_range    = c(-1, 1),
        estimator              = estimator,
        weight_method          = "test"
      )
    )
  }

  qcmde_ok <- ordinate("q_grid_cmde", 1, .02)
  expect_true(.iwmde_posterior_ordinate_supports_bf(qcmde_ok))
  expect_length(.iwmde_posterior_ordinate_warnings(qcmde_ok), 0L)

  qcmde_warn <- ordinate("q_grid_cmde", 1, .03)
  expect_true(.iwmde_posterior_ordinate_supports_bf(qcmde_warn))
  expect_true(any(grepl(
    "qCMDE.*ordinate.*3%.*2.5%.*5%",
    .iwmde_posterior_ordinate_warnings(qcmde_warn)
  )))

  qcmde_fail <- ordinate("q_grid_cmde", 1, .06)
  expect_false(.iwmde_posterior_ordinate_supports_bf(qcmde_fail))
  expect_match(
    .iwmde_posterior_ordinate_failure_reasons(qcmde_fail),
    "qCMDE.*ordinate.*6%.*5%"
  )

  iwmde_warn <- ordinate("iwmde", .94)
  expect_true(.iwmde_posterior_ordinate_supports_bf(iwmde_warn))
  expect_true(any(grepl(
    "IWMDE.*6%.*5%.*10%",
    .iwmde_posterior_ordinate_warnings(iwmde_warn)
  )))

  iwmde_fail <- ordinate("iwmde", .89)
  expect_false(.iwmde_posterior_ordinate_supports_bf(iwmde_fail))
  expect_match(
    .iwmde_posterior_ordinate_failure_reasons(iwmde_fail),
    "IWMDE.*11%.*10%"
  )

  bf <- .iwmde_bf_append_warning(1, iwmde_warn)
  expect_true(any(grepl("IWMDE.*6%", attr(bf, "warnings"))))

  posterior <- structure(
    stats::rnorm(20),
    posterior_ordinate = qcmde_warn
  )
  table <- data.frame(BF = "1", row.names = "mu")
  table <- .hypothesis_brma_append_iwmde_warnings(
    table     = table,
    posterior = posterior
  )
  expect_true(any(grepl("qCMDE.*ordinate.*3%", attr(table, "warnings"))))

  multi_ordinate <- list(
    status    = "ok",
    ordinates = list(
      list(
        value       = 0,
        ordinate    = 1,
        diagnostics = list(warning = "warning for zero"),
        parameter   = "theta"
      ),
      list(
        value       = 1,
        ordinate    = 1,
        diagnostics = list(warning = "warning for one"),
        parameter   = "theta"
      )
    )
  )
  posterior <- structure(
    stats::rnorm(20),
    posterior_ordinate = multi_ordinate
  )
  table <- data.frame(
    Alternative = c("theta != 0", "theta != 1"),
    Null        = c("theta = 0", "theta = 1"),
    BF          = c("1", "1"),
    row.names   = c("theta", "theta1"),
    check.names = FALSE
  )
  attr(table, "parsed") <- list(
    list(
      left  = list(type = "point", expr = "theta", value = 0, label = "theta = 0"),
      right = list(type = "not_point", expr = "theta", value = 0, label = "theta != 0")
    ),
    list(
      left  = list(type = "point", expr = "theta", value = 1, label = "theta = 1"),
      right = list(type = "not_point", expr = "theta", value = 1, label = "theta != 1")
    )
  )
  table <- .hypothesis_brma_append_iwmde_warnings(
    table     = table,
    posterior = posterior
  )
  expect_equal(
    attr(table, "warnings"),
    c(theta = "warning for zero", theta1 = "warning for one")
  )
  compact_table <- data.frame(
    BF          = c("1", "1"),
    row.names   = c("theta", "theta1"),
    check.names = FALSE
  )
  attr(compact_table, "parsed") <- attr(table, "parsed", exact = TRUE)
  compact_table <- .hypothesis_brma_append_iwmde_warnings(
    table     = compact_table,
    posterior = posterior
  )
  expect_equal(
    attr(compact_table, "warnings"),
    c(theta = "warning for zero", theta1 = "warning for one")
  )

  factor_posterior <- list(
    a = structure(
      stats::rnorm(20),
      posterior_ordinate = list(
        value       = 0,
        ordinate    = 1,
        diagnostics = list(warning = "warning for a"),
        parameter   = "mu",
        level       = "a"
      )
    ),
    b = structure(
      stats::rnorm(20),
      posterior_ordinate = list(
        value       = 0,
        ordinate    = 1,
        diagnostics = list(warning = "warning for b"),
        parameter   = "mu",
        level       = "b"
      )
    )
  )
  table <- data.frame(
    Alternative = c("mu != 0", "mu != 0"),
    Null        = c("mu = 0", "mu = 0"),
    BF          = c("1", "1"),
    row.names   = c("mu[a]", "mu[b]"),
    check.names = FALSE
  )
  attr(table, "parsed") <- list(
    list(
      left  = list(type = "point", expr = "mu", value = 0, label = "mu = 0"),
      right = list(type = "not_point", expr = "mu", value = 0, label = "mu != 0")
    )
  )
  table <- .hypothesis_brma_append_iwmde_warnings(
    table     = table,
    posterior = factor_posterior
  )
  expect_equal(
    attr(table, "warnings"),
    c("mu[a]" = "warning for a", "mu[b]" = "warning for b")
  )

  raw_diagnostic <- list(
    status      = "ok",
    diagnostics = list(
      bf_ordinate            = .4,
      bf_relative_mcse       = .1,
      bf_finite_terms        = 60,
      bf_ess                 = 30,
      bf_max_weight_share    = .2,
      active_mass            = 1,
      final_normalization_integral = 1,
      normalization_relative_error = 0,
      bf_ordinate_relative_change = .06,
      max_normalizer_relative_change = .06,
      normalization_range    = c(-1, 1),
      estimator              = "q_grid_cmde"
    )
  )
  expect_match(
    .hypothesis_brma_diagnostic_reason(raw_diagnostic),
    "qCMDE.*ordinate.*6%.*5%"
  )
})


test_that("qCMDE and IWMDE BF warnings cover Monte Carlo reliability", {

  diagnostics <- list(
    evaluation_value       = 0,
    relative_mcse          = .173,
    finite_terms           = 60,
    ess                    = 67.8,
    max_weight_share       = .30,
    active_mass            = 1,
    normalization_relative_error = 0,
    ordinate_relative_change = 0,
    max_normalizer_relative_change = 0,
    normalization_range    = c(-1, 1),
    estimator              = "q_grid_cmde",
    weight_method          = "test"
  )
  ordinate <- BayesTools::posterior_ordinate_attribute(
    value          = 0,
    ordinate       = .4,
    method         = "q_grid_cmde",
    density_method = "qCMDE",
    diagnostics    = diagnostics
  )
  warnings <- .iwmde_posterior_ordinate_warnings(ordinate)

  expect_true(.iwmde_posterior_ordinate_supports_bf(ordinate))
  expect_true("qCMDE relative MCSE is 17.3%" %in% warnings)
  expect_true(any(grepl("uses only\\s+60.*finite importance terms.*100", warnings)))
  expect_true("qCMDE effective sample size is 67.8" %in% warnings)
  expect_true(any(grepl("largest importance weight.*30%.*20%", warnings)))

  diagnostics[["estimator"]] <- "iwmde"
  warnings <- .iwmde_diagnostics_bf_warning(diagnostics)
  expect_true("IWMDE relative MCSE is 17.3%" %in% warnings)
  expect_true("IWMDE effective sample size is 67.8" %in% warnings)
})


test_that("old IWMDE variant labels are rejected", {

  .skip_if_missing_raw_fits(c("bcg_meta-analysis", "bcg_meta-regression2"))

  expect_error(
    plot(
      load_fit("bcg_meta-analysis", validate = FALSE),
      parameter      = "mu",
      density_method = "IWMDE-CMDE"
    ),
    "must be one of"
  )
  expect_error(
    marginal_means(
      load_fit("bcg_meta-regression2", validate = FALSE),
      n_samples      = 1000,
      density_method = "IWMDE-Chen"
    ),
    "must be one of"
  )
})

test_that("unused dots warn in cleaned density interfaces", {

  .skip_if_missing_raw_fits(c("bcg_meta-analysis", "bcg_meta-regression2"))

  expect_warning(
    plot(
      load_fit("bcg_meta-analysis", validate = FALSE),
      parameter       = "mu",
      density_method  = "qCMDE",
      density_control = list(n_points = 20, max_samples = 20),
      iwmde_n_points  = 20
    ),
    "Unused argument.*iwmde_n_points"
  )
  expect_warning(
    marginal_means(
      load_fit("bcg_meta-regression2", validate = FALSE),
      n_samples          = 1000,
      density_method     = "qCMDE",
      density_control    = list(n_points = 20, max_samples = 20),
      iwmde_max_samples  = 20
    ),
    "Unused argument.*iwmde_max_samples"
  )
})

test_that("density_control validates public density settings", {

  valid <- .density_control_normalize(
    density_method  = "qCMDE",
    density_control = list(
      n_points             = 20,
      max_samples          = 30,
      normalization_points = 40,
      normalization_prob   = .95,
      display_grid         = "uniform"
    )
  )

  expect_equal(valid[["n_points"]], 20)
  expect_equal(valid[["max_samples"]], 30)
  expect_equal(valid[["normalization_points"]], 40)
  expect_equal(valid[["normalization_prob"]], .95)
  expect_equal(valid[["display_grid"]], "uniform")
  expect_error(
    .density_control_normalize("qCMDE", list(unknown = 1)),
    "unrecognized"
  )
  expect_error(
    .density_control_normalize("qCMDE", list(20)),
    "fully named"
  )
  expect_error(
    .density_control_normalize(
      "qCMDE",
      structure(list(20, 30), names = c("n_points", "n_points"))
    ),
    "duplicate"
  )
  expect_error(
    .density_control_normalize("KDE", list(n_points = 20)),
    "only used"
  )
  valid_iwmde <- .density_control_normalize(
    density_method  = "IWMDE",
    density_control = list(
      normalization_points = 40,
      normalization_prob   = .95
    )
  )
  expect_equal(valid_iwmde[["n_points"]], 100)
  expect_equal(valid_iwmde[["max_samples"]], 1000)
  expect_equal(valid_iwmde[["normalization_points"]], 40)
  expect_equal(valid_iwmde[["normalization_prob"]], .95)
  mm_iwmde <- .hypothesis_marginal_means_density_control(
    object = list(density_settings = list(
      n_points             = 20,
      max_samples          = 30,
      normalization_points = 40,
      normalization_prob   = .95,
      display_grid         = "uniform"
    )),
    density_method  = "IWMDE",
    density_control = NULL
  )
  expect_equal(mm_iwmde[["normalization_points"]], 40)
  expect_equal(mm_iwmde[["normalization_prob"]], .95)
  expect_error(
    .density_control_normalize("qCMDE", list(normalization_prob = 0)),
    "higher than 0"
  )
  expect_error(
    .density_control_normalize("qCMDE", list(normalization_prob = NA_real_)),
    "cannot contain"
  )
})

test_that("IWMDE row sampling is reproducible, nested, and RNG-local", {

  rows <- seq(2L, 80L, by = 2L)
  set.seed(123)
  old_seed <- .Random.seed
  old_kind <- RNGkind()

  out10 <- .nested_srs_rows(rows = rows, max_samples = 10)
  out20 <- .nested_srs_rows(rows = rows, max_samples = 20)

  expect_length(out10, 10L)
  expect_length(out20, 20L)
  expect_true(all(out10 %in% out20))
  expect_true(is.unsorted(out10) == FALSE)
  expect_identical(.Random.seed, old_seed)
  expect_identical(RNGkind(), old_kind)
  expect_equal(
    .nested_srs_rows(rows = rows, max_samples = 10),
    out10
  )
})


test_that("nested row sampling is independent of product-state labels", {

  rows     <- seq_len(10000L)
  state    <- rep("common", length(rows))
  selected <- .nested_srs_rows(rows = rows, max_samples = 1000L)
  omitted  <- setdiff(rows, selected)[[1L]]
  state[omitted] <- "rare"

  expect_equal(length(selected), 1000L)
  expect_true(all(selected %in% rows))
  expect_false(omitted %in% selected)
  expect_false("rare" %in% state[selected])
  expect_identical(
    .nested_srs_rows(rows = rows, max_samples = 1000L),
    selected
  )
  expect_false("group" %in% names(formals(.nested_srs_rows)))
})


test_that("IWMDE fails when an eligible population omits a fitted chain", {

  context <- list(
    posterior_samples = matrix(
      seq_len(50L),
      ncol     = 1L,
      dimnames = list(NULL, "mu")
    ),
    chain_id = rep(1:2, each = 25L),
    object = list(
      fit        = list(),
      likelihood = list(family = "normal"),
      data       = list(measure = "GEN")
    ),
    data            = list(),
    priors          = list(),
    flat_prior_list = list()
  )
  testthat::local_mocked_bindings(
    .iwmde_parameter_components = function(...) {
      list(
        active         = c(rep(TRUE, 25L), rep(FALSE, 25L)),
        point_masses   = data.frame(x = numeric(), mass = numeric()),
        point_location = rep(NA_real_, 50L)
      )
    },
    .package = "RoBMA"
  )

  failure <- tryCatch(
    .iwmde_plan(
      context         = context,
      parameter       = "mu",
      density_method  = "qCMDE",
      density_control = list(n_points = 20, max_samples = 20),
      outputs         = "ordinate",
      values          = 0,
      parameter_spec  = list(type = "primitive")
    ),
    iwmde_construction_error = identity
  )

  expect_s3_class(failure, "iwmde_construction_error")
  expect_equal(failure[["stage"]], "chain-coverage validation")
  expect_equal(
    failure[["chain_coverage"]][["expected_chain_ids"]],
    1:2
  )
  expect_equal(
    unname(failure[["chain_coverage"]][["eligible_chain_counts"]]),
    c(25L, 0L)
  )
  expect_match(conditionMessage(failure), "fitted chain\\(s\\): 2")
})


test_that("IWMDE plan fails the target on a non-finite baseline row", {

  context <- list(
    posterior_samples = matrix(
      seq_len(25),
      ncol     = 1L,
      dimnames = list(NULL, "mu")
    ),
    object = list(
      fit        = list(),
      likelihood = list(family = "normal"),
      data       = list(measure = "GEN")
    ),
    data            = list(),
    priors          = list(),
    flat_prior_list = list()
  )
  drop_positions <- 1L
  condition_rows <- rep(TRUE, nrow(context[["posterior_samples"]]))

  testthat::local_mocked_bindings(
    .iwmde_row_states = function(context, rows, parameter = NULL,
                                 parameter_spec = NULL, estimator = NULL) {

      lapply(seq_along(rows), function(i) {
        baseline_log_q <- if (i %in% drop_positions) -Inf else 0
        if (is.finite(baseline_log_q)) {
          .iwmde_new_row_state(list(baseline_log_q = baseline_log_q))
        } else {
          list(baseline_log_q = baseline_log_q)
        }
      })
    },
    .iwmde_parameter_condition_rows = function(context, parameter_spec) {
      condition_rows
    },
    .package = "RoBMA"
  )

  make_plan <- function() {
    .iwmde_plan(
      context         = context,
      parameter       = "mu",
      density_method  = "qCMDE",
      density_control = list(
        n_points             = 20,
        samples              = 21,
        normalization_points = 20
      ),
      outputs         = "ordinate",
      values          = 0,
      parameter_spec  = list(type = "primitive"),
      metadata        = NULL
    )
  }

  expect_error(
    make_plan(),
    "qCMDE construction failed.*target 'mu'.*posterior row.*baseline joint-density"
  )

  drop_positions <- integer()
  context[["posterior_samples"]][3L, "mu"] <- NA_real_
  expect_error(
    make_plan(),
    "target 'mu'.*posterior row 3.*target-draw validation"
  )
  condition_rows[[3L]] <- FALSE
  expect_equal(make_plan()[["status"]], "ok")
  context[["posterior_samples"]][3L, "mu"] <- 3
  condition_rows[[3L]] <- TRUE

  drop_positions <- integer()
  plan <- make_plan()
  expect_equal(plan[["status"]], "ok")
  expect_equal(plan[["rows"]][["n_denominator_rows"]], 21L)
  expect_equal(plan[["rows"]][["n_candidate_rows"]], 21L)
  expect_equal(plan[["rows"]][["n_estimator_rows"]], 21L)
  expect_equal(length(plan[["rows"]][["row_states"]]), 21L)
  expect_equal(
    plan[["rows"]][["estimator_rows"]],
    plan[["rows"]][["candidate_rows"]]
  )

  execution <- .iwmde_plan_row_execution(
    context = context,
    plan    = plan
  )
  expect_equal(execution[["active_rows"]], plan[["rows"]][["estimator_rows"]])
  expect_equal(execution[["n_candidate_rows"]], 21L)

  density_plan <- .iwmde_plan(
    context         = context,
    parameter       = "mu",
    density_method  = "qCMDE",
    density_control = list(
      n_points             = 20,
      max_samples          = 21,
      normalization_points = 20
    ),
    outputs         = "density",
    parameter_spec  = list(type = "primitive")
  )
  expected_tail <- stats::quantile(
    seq_len(25),
    probs = c(.05, .95),
    names = FALSE,
    type  = 8
  )
  expect_equal(density_plan[["grids"]][["tail_probabilities"]], c(.05, .95))
  expect_equal(density_plan[["grids"]][["tail_values"]], expected_tail)
  expect_true(all(
    expected_tail %in% density_plan[["grids"]][["display_grid"]]
  ))

  drop_positions <- 1:2
  expect_error(make_plan(), "posterior rows .*baseline joint-density")
})


test_that("IWMDE row plans keep population mass separate from SRS rows", {

  component <- list(
    active       = c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE),
    point_masses = data.frame(x = numeric(), mass = numeric())
  )
  rows <- .iwmde_plan_rows(
    posterior_values = seq_len(6L),
    finite_rows      = rep(TRUE, 6L),
    component        = component,
    continuous_rows  = c(1L, 3L, 5L),
    candidate_rows   = c(1L, 5L),
    candidate_values = c(1, 5),
    continuous_values = c(1, 3, 5)
  )

  expect_equal(rows[["active_mass"]], .5)
  expect_equal(rows[["continuous_rows"]], c(1L, 3L, 5L))
  expect_equal(rows[["n_denominator_rows"]], 2L)
  expect_equal(rows[["candidate_rows"]], c(1L, 5L))
})


test_that("IWMDE condition rows keep one-row shape and ignore missing values", {

  context <- list(
    posterior_samples = data.frame(mu = 0),
    flat_prior_list   = list()
  )
  rows <- .iwmde_parameter_condition_rows(
    context = context,
    parameter_spec = list(
      conditional      = c("missing_a", NA, "missing_b"),
      conditional_rule = "OR"
    )
  )

  expect_length(rows, 1L)
  expect_false(rows)
})

test_that("IWMDE marginal-mean specs preserve child condition metadata", {

  condition_key <- paste(c("OR", 2L, "mu_alloc", "mu_intercept"), collapse = "\r")
  condition_event <- structure(
    list(
      conditional      = c("mu_intercept", "mu_alloc"),
      conditional_rule = "OR",
      condition_key    = condition_key
    ),
    class = "BayesTools_condition_event"
  )
  prior_density_context <- structure(
    list(condition_key = condition_key),
    class = "prior_density_context"
  )
  prior_density <- BayesTools:::.prior_linear_density_point(0)
  samples <- structure(
    stats::rnorm(20),
    linear_weights              = c(mu_intercept = 1, mu_alloc = 1),
    conditional                 = "stale_parent_condition",
    conditional_rule            = "AND",
    effective_conditional       = c("mu_intercept", "mu_alloc"),
    effective_conditional_rule  = "OR",
    condition_key               = condition_key,
    condition_event             = condition_event,
    resolved_condition_event    = condition_event,
    prior_density_context       = prior_density_context,
    prior_density               = prior_density
  )
  marginal_means_object <- list(
    inference        = list(
      conditional = list(mu_alloc = list(alternate = samples))
    ),
    parameters       = "mu_alloc",
    term_map         = data.frame(
      term             = "alloc",
      parameter        = "mu_alloc",
      label            = "alloc",
      stringsAsFactors = FALSE
    ),
    conditional_rule = "AND"
  )

  specs <- .iwmde_marginal_means_specs(
    marginal_means_object = marginal_means_object,
    parameter             = "alloc",
    type                  = "conditional",
    levels                = "alternate"
  )
  spec <- specs[["alloc: alternate"]]

  expect_equal(spec[["conditional"]], c("mu_intercept", "mu_alloc"))
  expect_equal(spec[["conditional_rule"]], "OR")
  expect_equal(spec[["condition_key"]], condition_key)
  expect_identical(spec[["condition_event"]], condition_event)
  expect_identical(spec[["resolved_condition_event"]], condition_event)
  expect_identical(spec[["prior_density_context"]], prior_density_context)
  expect_identical(spec[["prior_density"]], prior_density)
  expect_equal(.iwmde_parameter_condition_key(spec), condition_key)
  expect_equal(
    .iwmde_posterior_metadata(samples, "mu_alloc", "alternate")[["condition_key"]],
    condition_key
  )

  stale_samples <- samples
  attr(stale_samples, "effective_conditional") <- NULL
  attr(stale_samples, "effective_conditional_rule") <- NULL
  attr(stale_samples, "condition_key") <- NULL
  marginal_means_object[["inference"]][["conditional"]][["mu_alloc"]][["alternate"]] <- stale_samples

  expect_error(
    .iwmde_marginal_means_specs(
      marginal_means_object = marginal_means_object,
      parameter             = "alloc",
      type                  = "conditional",
      levels                = "alternate"
    ),
    "metadata is incomplete"
  )
})

test_that("restricted Chen fallback weights are proper densities", {

  gamma_weight <- .iwmde_chen_gamma_log_weight_single(
    active_values = seq(.001, 3, length.out = 500),
    weight_values = c(.25, .30, .35, .40, .45, .50),
    support       = c(0, Inf)
  )
  gamma_mass <- .iwmde_trapz(
    seq(.001, 3, length.out = 500),
    exp(gamma_weight[["log_weight"]])
  )

  beta_grid <- seq(.001, .999, length.out = 500)
  beta_weight <- .iwmde_chen_beta_log_weight(
    active_values = beta_grid,
    weight_values = c(.25, .30, .35, .40, .45, .50),
    support       = c(0, 1)
  )
  beta_mass <- .iwmde_trapz(beta_grid, exp(beta_weight[["log_weight"]]))

  expect_equal(gamma_weight[["method"]], "chen_gamma")
  expect_equal(beta_weight[["method"]], "chen_beta")
  expect_equal(gamma_mass, 1, tolerance = .01)
  expect_equal(beta_mass, 1, tolerance = .01)
})

test_that("IWMDE Chen conditional-normal weights error on failed conditioning", {

  context <- list(
    posterior_samples = data.frame(
      mu  = c(-0.2, -0.1, 0.1, 0.2),
      tau = c(0.3, 0.4, 0.5, 0.6)
    ),
    indicator_names = character(),
    selection_spec  = NULL
  )

  expect_error(
    .iwmde_chen_conditional_normal_log_weight(
      context        = context,
      parameter      = "mu",
      parameter_spec = list(type = "primitive"),
      active_rows    = 1:4,
      active_values  = context[["posterior_samples"]][["mu"]],
      weight_rows    = 1:4,
      weight_values  = c(NA_real_, NA_real_, 0.1, NA_real_)
    ),
    "conditioning fit rows contain non-finite values"
  )

  expect_error(
    .iwmde_chen_logit_conditional_normal_log_weight(
      context        = context,
      parameter      = "mu",
      parameter_spec = list(type = "primitive"),
      active_rows    = 1:4,
      active_values  = c(0.2, 0.3, 0.4, 0.5),
      weight_rows    = 1:4,
      weight_values  = c(NA_real_, NA_real_, 0.4, NA_real_),
      support        = c(0, 1)
    ),
    "conditioning fit rows contain non-finite values"
  )
})
