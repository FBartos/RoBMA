context("brma.mv known-V backend preparation")

test_that("brma.mv honors forced known-V backend requests", {

  V <- matrix(c(0.04, 0.015, 0.015, 0.09), nrow = 2)

  for (parameterization in c("latent", "whitened", "block_mvn")) {
    object <- brma.mv(
      yi                        = c(0.10, 0.20),
      V                         = V,
      known_v_parameterization  = parameterization,
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_data                 = TRUE
    )
    known_V <- attr(object[["data"]], "known_V_data")

    expect_equal(known_V[["parameterization_requested"]], parameterization)
    expect_equal(known_V[["parameterization"]], parameterization)
    expect_equal(.data_known_v_parameterization(object[["data"]]), parameterization)
  }
})


test_that("LOO target comparison rejects missing likelihood target labels", {

  metadata <- list(
    unit             = "estimate",
    retained_context = "remaining_data",
    data_hash        = "same-data",
    target           = "estimate_log_score"
  )
  loo_a <- structure(list(), class = "loo")
  loo_b <- structure(list(), class = "loo")
  attr(loo_a, "RoBMA_target") <- metadata
  attr(loo_b, "RoBMA_target") <- metadata[names(metadata) != "target"]

  expect_error(
    .check_loo_compare_targets(list(loo_a, loo_b)),
    "without likelihood target labels"
  )
})


test_that("brma.mv adds latent known-V factors to fit syntax", {

  object <- brma.mv(
    yi                        = c(0.10, 0.20),
    V                         = matrix(c(0.04, 0.01, 0.01, 0.09), nrow = 2),
    known_v_parameterization  = "latent",
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )

  fit_priors <- .create_fit_priors(object[["data"]], object[["priors"]])
  fit_data   <- .create_fit_data(object[["data"]], object[["priors"]])
  syntax     <- .create_model_syntax(object[["data"]], object[["priors"]])

  expect_true("sampling_z" %in% names(fit_priors))
  expect_true("sampling_var" %in% names(fit_data))
  expect_true("sampling_B_1" %in% names(fit_data))
  expect_false("sampling_B" %in% names(fit_data))
  expect_match(syntax, "sampling_dependency", fixed = TRUE)
  expect_match(syntax, "sampling_var\\[i\\]")
  expect_silent(
    .check_log_lik_target_available(object, "estimate", "test()")
  )
  expect_error(
    .check_log_lik_target_available(object, "cluster", "test()"),
    "unit = 'cluster'"
  )
  expect_silent(
    .check_log_lik_target_available(object, "estimate", "add_waic()")
  )
})


test_that("diagonal brma.mv estimate likelihood matches brma.norm factorization", {

  yi <- c(0.10, 0.20, -0.05)
  vi <- c(0.04, 0.09, 0.16)
  mu <- matrix(
    c(
      0.00, 0.10, -0.10,
      0.05, 0.15, -0.02
    ),
    nrow  = 2,
    byrow = TRUE
  )
  tau <- matrix(
    c(
      0.10, 0.20, 0.05,
      0.05, 0.15, 0.08
    ),
    nrow  = 2,
    byrow = TRUE
  )

  mv_object <- brma.mv(
    yi                        = yi,
    V                         = diag(vi),
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  uni_object <- brma.norm(
    yi                        = yi,
    vi                        = vi,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )

  make_setup <- function(object) {
    list(
      fit               = NULL,
      priors            = object[["priors"]],
      data              = object[["data"]],
      yi                = object[["data"]][["outcome"]][["yi"]],
      sei               = object[["data"]][["outcome"]][["sei"]],
      selection_sei     = object[["data"]][["outcome"]][["sei"]],
      K                 = length(yi),
      S                 = nrow(mu),
      mu                = mu,
      tau_within        = tau,
      tau_between       = matrix(0, nrow = nrow(mu), ncol = length(yi)),
      cluster           = NULL,
      weights           = NULL,
      data_hash         = .get_outcome_hash(object),
      is_weightfunction = FALSE,
      outcome_type      = "norm",
      effect_direction  = "positive",
      posterior_samples = matrix(numeric(0), nrow = nrow(mu), ncol = 0L)
    )
  }

  mv_setup  <- make_setup(mv_object)
  uni_setup <- make_setup(uni_object)

  expect_false(.known_v_estimate_target_uses_backend(mv_object[["data"]]))
  expect_equal(
    .log_lik_estimate_from_setup(mv_setup),
    .log_lik_estimate_from_setup(uni_setup),
    tolerance = 1e-12
  )

  target <- .estimate_log_lik_target_metadata(
    setup     = mv_setup,
    data_hash = .get_outcome_hash(mv_object)
  )
  expect_equal(target[["target"]], "estimate_log_score")
  expect_false(target[["known_v_estimate_backend"]])
  expect_false(target[["dependency_conditioning"]])
})


test_that("diagonal brma.mv random estimate target remains factorized", {

  object <- brma.mv(
    yi                         = c(0.10, 0.20, -0.05),
    V                          = diag(c(0.04, 0.09, 0.16)),
    random                     = ~ 1 | study,
    data                       = data.frame(study = c("a", "b", "c")),
    measure                    = "GEN",
    prior_unit_information_sd  = 1,
    marginalize_estimate_level = FALSE,
    only_priors                = TRUE
  )
  setup <- list(
    fit               = NULL,
    priors            = object[["priors"]],
    data              = object[["data"]],
    yi                = object[["data"]][["outcome"]][["yi"]],
    sei               = object[["data"]][["outcome"]][["sei"]],
    selection_sei     = object[["data"]][["outcome"]][["sei"]],
    K                 = 3L,
    S                 = 1L,
    mu                = matrix(c(0, 0, 0), nrow = 1),
    tau_within        = matrix(0, nrow = 1, ncol = 3),
    tau_between       = matrix(0, nrow = 1, ncol = 3),
    cluster           = NULL,
    weights           = NULL,
    data_hash         = .get_outcome_hash(object),
    is_weightfunction = FALSE,
    outcome_type      = "norm",
    effect_direction  = "positive",
    posterior_samples = matrix(numeric(0), nrow = 1, ncol = 0L)
  )

  target <- .estimate_log_lik_target_metadata(
    setup             = setup,
    data_hash         = .get_outcome_hash(object),
    dependency_blocks = .random_effect_dependency_blocks(
      sampling_covariance = .known_v_dependency_covariance(
        object[["data"]],
        sampling_latent_marginalized = TRUE
      ),
      formula_design = .fitted_formula_design(
        object,
        "mu",
        required = TRUE
      ),
      blocks = .data_sampled_random_effect_blocks(
        object[["data"]]
      )
    )
  )
  expect_true(target[["known_v_estimate_backend"]])
  expect_false(target[["dependency_conditioning"]])
  expect_equal(target[["target"]], "estimate_log_score")
  expect_equal(target[["dependency_component_sizes"]], rep(1L, 3L))
  expect_equal(target[["random_effect_representation"]], "sampled")
  expect_equal(target[["latent_effect_handling"]], "integrated")
})


test_that("brma.mv omits latent factors for diagonal known V", {

  object <- brma.mv(
    yi                        = c(0.10, 0.20),
    V                         = diag(c(0.04, 0.09)),
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )

  expect_equal(
    attr(object[["data"]], "known_V_data")[["parameterization"]],
    "whitened"
  )
  expect_equal(
    attr(object[["data"]], "known_V_data")[["parameterization_requested"]],
    "auto"
  )

  fit_priors <- .create_fit_priors(object[["data"]], object[["priors"]])
  fit_data   <- .create_fit_data(object[["data"]], object[["priors"]])
  syntax     <- .create_model_syntax(object[["data"]], object[["priors"]])

  expect_false("sampling_z" %in% names(fit_priors))
  expect_false("sampling_rank" %in% names(fit_data))
  expect_false("sampling_B" %in% names(fit_data))
  expect_false(grepl("sampling_dependency", syntax, fixed = TRUE))
  expect_silent(
    .check_log_lik_target_available(object, "estimate", "test()")
  )
})


test_that("brma.mv auto selects exact known-V backend when feasible", {

  dat <- data.frame(
    yi    = c(0.10, 0.20, 0.30),
    x     = c(0, 1, 2),
    study = c("s1", "s2", "s3")
  )
  V <- matrix(
    c(
      0.04, 0.02, 0.01,
      0.02, 0.09, 0.03,
      0.01, 0.03, 0.16
    ),
    nrow  = 3,
    byrow = TRUE
  )

  no_scale <- brma.mv(
    yi                        = yi,
    V                         = V,
    data                      = dat,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_data                 = TRUE
  )
  no_scale_known_V <- attr(no_scale[["data"]], "known_V_data")
  expect_equal(no_scale_known_V[["parameterization_requested"]], "auto")
  expect_equal(no_scale_known_V[["parameterization"]], "whitened")

  random_scalar <- brma.mv(
    yi                        = yi,
    V                         = V,
    random                    = ~ 1 | study,
    data                      = dat,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  random_scalar_known_V <- attr(random_scalar[["data"]], "known_V_data")
  expect_equal(random_scalar_known_V[["parameterization_requested"]], "auto")
  expect_equal(random_scalar_known_V[["parameterization"]], "whitened")
  expect_true(.data_has_marginalized_random_effects(random_scalar[["data"]]))
  expect_false(.marginalized_random_effects_row_varying(
    .data_marginalized_random_effects(random_scalar[["data"]])
  ))

  scale_small <- brma.mv(
    yi                        = yi,
    V                         = V,
    scale                     = ~ x,
    data                      = dat,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_data                 = TRUE
  )
  scale_small_known_V <- attr(scale_small[["data"]], "known_V_data")
  expect_equal(scale_small_known_V[["parameterization_requested"]], "auto")
  expect_equal(scale_small_known_V[["parameterization"]], "block_mvn")

  old_max_block <- getOption("RoBMA.known_v_block_mvn_max_block_size", NULL)
  on.exit({
    options(RoBMA.known_v_block_mvn_max_block_size = old_max_block)
  }, add = TRUE)

  options(RoBMA.known_v_block_mvn_max_block_size = 3L)
  scale_boundary <- brma.mv(
    yi                        = yi,
    V                         = V,
    scale                     = ~ x,
    data                      = dat,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_data                 = TRUE
  )
  expect_equal(
    attr(scale_boundary[["data"]], "known_V_data")[["parameterization"]],
    "block_mvn"
  )

  options(RoBMA.known_v_block_mvn_max_block_size = Inf)
  scale_inf <- brma.mv(
    yi                        = yi,
    V                         = V,
    scale                     = ~ x,
    data                      = dat,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_data                 = TRUE
  )
  expect_equal(
    attr(scale_inf[["data"]], "known_V_data")[["parameterization"]],
    "block_mvn"
  )

  random_scale_small <- brma.mv(
    yi                        = yi,
    V                         = V,
    scale                     = ~ x,
    random                    = ~ 1 | study,
    data                      = dat,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  random_scale_small_known_V <- attr(random_scale_small[["data"]], "known_V_data")
  expect_equal(random_scale_small_known_V[["parameterization_requested"]], "auto")
  expect_equal(random_scale_small_known_V[["parameterization"]], "block_mvn")
  expect_true(.data_has_marginalized_random_effects(random_scale_small[["data"]]))
  expect_true(.marginalized_random_effects_row_varying(
    .data_marginalized_random_effects(random_scale_small[["data"]])
  ))

  options(RoBMA.known_v_block_mvn_max_block_size = 2L)

  scale_large <- brma.mv(
    yi                        = yi,
    V                         = V,
    scale                     = ~ x,
    data                      = dat,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_data                 = TRUE
  )
  scale_large_known_V <- attr(scale_large[["data"]], "known_V_data")
  expect_equal(scale_large_known_V[["parameterization_requested"]], "auto")
  expect_equal(scale_large_known_V[["parameterization"]], "latent")

  random_scale_large <- brma.mv(
    yi                        = yi,
    V                         = V,
    scale                     = ~ x,
    random                    = ~ 1 | study,
    data                      = dat,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  random_scale_large_known_V <- attr(random_scale_large[["data"]], "known_V_data")
  expect_equal(random_scale_large_known_V[["parameterization_requested"]], "auto")
  expect_equal(random_scale_large_known_V[["parameterization"]], "latent")

  explicit_block <- brma.mv(
    yi                        = yi,
    V                         = V,
    scale                     = ~ x,
    data                      = dat,
    known_v_parameterization  = "block_mvn",
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_data                 = TRUE
  )
  explicit_known_V <- attr(explicit_block[["data"]], "known_V_data")
  expect_equal(explicit_known_V[["parameterization_requested"]], "block_mvn")
  expect_equal(explicit_known_V[["parameterization"]], "block_mvn")

  options(RoBMA.known_v_block_mvn_max_block_size = NA_real_)
  expect_error(
    brma.mv(
      yi                        = yi,
      V                         = V,
      scale                     = ~ x,
      data                      = dat,
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_data                 = TRUE
    ),
    "single positive"
  )
})


test_that("brma.mv prepares whitened known-V backend", {

  V <- matrix(c(0.04, 0.03, 0.03, 0.09), nrow = 2)

  expect_warning(
    object <- brma.mv(
      yi                        = c(0.10, 0.20),
      V                         = V,
      known_v_parameterization  = "whitened",
      known_v_residual_fraction = 0.50,
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_data                 = TRUE
    ),
    "disregarded"
  )

  known_V <- attr(object[["data"]], "known_V_data")
  expect_equal(known_V[["parameterization"]], "whitened")
  expect_true(known_V[["correlated"]])
  expect_equal(known_V[["rank"]], 0L)
  expect_equal(known_V[["residual_fraction_requested"]], 0.50)

  whitening_block <- known_V[["whitening_blocks"]][[1L]]
  rotation        <- whitening_block[["rotation"]]
  rotated_V       <- rotation %*%
    V[
      whitening_block[["index"]],
      whitening_block[["index"]],
      drop = FALSE
    ] %*%
    t(rotation)
  expect_equal(
    rotated_V,
    diag(whitening_block[["variance"]], nrow = 2),
    tolerance = 1e-10
  )
  expect_null(known_V[["whitening_matrix"]])
  expect_null(known_V[["sampling_factor"]])
})


test_that("brma.mv auto known-V update preserves residual fraction", {

  old_options <- options(RoBMA.known_v_block_mvn_max_block_size = 1L)
  on.exit(options(old_options))

  dat <- data.frame(
    yi       = c(0.10, 0.20, 0.30),
    estimate = factor(c("e1", "e2", "e3"), levels = c("e1", "e2", "e3"))
  )
  V <- matrix(
    c(
      0.040, 0.018, 0.012,
      0.018, 0.090, 0.024,
      0.012, 0.024, 0.160
    ),
    nrow  = 3,
    byrow = TRUE
  )
  R <- diag(c(1, 4, 9))
  dimnames(R) <- list(levels(dat[["estimate"]]), levels(dat[["estimate"]]))

  expect_warning(
    object <- brma.mv(
      yi                          = yi,
      V                           = V,
      random                      = ~ 1 | estimate,
      R                           = R,
      Rscale                      = "none",
      data                        = dat,
      known_v_parameterization    = "auto",
      known_v_residual_fraction   = 0.35,
      measure                     = "GEN",
      prior_unit_information_sd   = 1,
      only_priors                 = TRUE
    ),
    NA
  )

  known_V <- attr(object[["data"]], "known_V_data")
  expect_equal(known_V[["parameterization_requested"]], "auto")
  expect_equal(known_V[["parameterization"]], "latent")
  expect_equal(known_V[["residual_fraction_requested"]], 0.35)

  data <- brma.mv(
    yi                        = yi,
    V                         = diag(c(0.04, 0.09, 0.16)),
    data                      = data.frame(yi = c(0.10, 0.20, 0.30)),
    known_v_parameterization  = "auto",
    known_v_residual_fraction = 0.40,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_data                 = TRUE
  )[["data"]]
  term <- list(
    block_name           = "study",
    sd_parameter_names   = "tau",
    row_multiplier       = c(1, 4, 9),
    row_multiplier_name  = "known_r_multiplier"
  )

  updated_known_V <- attr(.known_v_auto_update_for_marginalized_random(
    data  = data,
    terms = list(term)
  ), "known_V_data")

  expect_equal(updated_known_V[["parameterization_requested"]], "auto")
  expect_equal(updated_known_V[["parameterization"]], "whitened")
  expect_equal(updated_known_V[["effective_backend"]], "diagonal")
  expect_equal(updated_known_V[["residual_fraction_requested"]], 0.40)
})


test_that("brma.mv marginalized known R contributes known-V row variance", {

  dat <- data.frame(
    yi       = c(0.10, 0.20, 0.30),
    estimate = factor(c("e1", "e2", "e3"), levels = c("e1", "e2", "e3"))
  )
  R <- diag(c(4, 9, 16))
  dimnames(R) <- list(levels(dat[["estimate"]]), levels(dat[["estimate"]]))

  object <- brma.mv(
    yi                        = yi,
    V                         = diag(0.01, nrow(dat)),
    random                    = ~ 1 | estimate,
    R                         = R,
    Rscale                    = "none",
    data                      = dat,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  known_V  <- attr(object[["data"]], "known_V_data")
  term     <- .data_marginalized_random_effects(object[["data"]])[[1L]]
  sd_name  <- term[["sd_parameter_names"]]
  row_name <- term[["row_multiplier_name"]]
  fit_data <- .create_fit_data(object[["data"]], object[["priors"]])
  syntax   <- .create_model_syntax(object[["data"]], object[["priors"]])

  expect_equal(known_V[["parameterization_requested"]], "auto")
  expect_equal(known_V[["parameterization"]], "whitened")
  expect_equal(known_V[["effective_backend"]], "diagonal")
  expect_equal(unname(term[["row_multiplier"]]), c(4, 9, 16))
  expect_equal(unname(fit_data[[row_name]]), c(4, 9, 16))
  expect_match(
    .data_marginalized_random_variance_expression(object[["data"]]),
    paste0("pow\\(", sd_name, ",2\\) \\* ", row_name, "\\[i\\]")
  )
  expect_match(syntax, paste0("sampling_var\\[i\\] \\+ pow\\(", sd_name))
  expect_match(syntax, paste0(row_name, "\\[i\\]"))

  posterior <- matrix(
    0.20,
    nrow     = 1L,
    dimnames = list(NULL, sd_name)
  )
  expected_extra_variance <- matrix(
    0.20^2 * c(4, 9, 16),
    nrow = 1L
  )
  expect_equal(
    .evaluate_marginalized_random_variance(
      data              = object[["data"]],
      posterior_samples = posterior
    ),
    expected_extra_variance
  )

  setup <- list(
    fit               = NULL,
    priors            = object[["priors"]],
    data              = object[["data"]],
    yi                = dat[["yi"]],
    sei               = object[["data"]][["outcome"]][["sei"]],
    selection_sei     = object[["data"]][["outcome"]][["sei"]],
    K                 = 3L,
    S                 = 1L,
    mu                = matrix(0, nrow = 1L, ncol = 3L),
    tau_within        = matrix(0, nrow = 1L, ncol = 3L),
    tau_between       = matrix(0, nrow = 1L, ncol = 3L),
    cluster           = NULL,
    weights           = NULL,
    data_hash         = .get_outcome_hash(object),
    is_weightfunction = FALSE,
    outcome_type      = "norm",
    effect_direction  = "positive",
    posterior_samples = posterior
  )
  expected_log_lik <- matrix(
    stats::dnorm(
      dat[["yi"]],
      mean = 0,
      sd   = sqrt(0.01 + as.numeric(expected_extra_variance)),
      log  = TRUE
    ),
    nrow = 1L
  )
  parameters <- list(mu = 0)
  parameters[[sd_name]] <- 0.20

  expect_equal(
    .known_v_extra_variance_from_setup(setup),
    expected_extra_variance
  )
  expect_equal(
    .log_lik_estimate_from_setup(setup),
    expected_log_lik,
    tolerance = 1e-12
  )
  expect_equal(
    .marglik_known_v_extra_variance(
      parameters         = parameters,
      model_data         = object[["data"]],
      bridge_context     = NULL,
      tau_within_samples = matrix(0, nrow = 1L, ncol = 3L),
      is_random          = TRUE,
      K                  = 3L
    ),
    expected_extra_variance
  )
  expect_equal(
    .log_posterior(
      parameters                 = parameters,
      data                       = fit_data,
      is_scale                   = FALSE,
      is_random                  = TRUE,
      is_multilevel              = FALSE,
      is_weights                 = FALSE,
      is_known_v                 = TRUE,
      is_PET                     = FALSE,
      is_PEESE                   = FALSE,
      is_weightfunction          = FALSE,
      effect_direction           = "positive",
      outcome_type               = "norm",
      model_data                 = object[["data"]]
    ),
    sum(expected_log_lik),
    tolerance = 1e-12
  )

  target <- .estimate_log_lik_target_metadata(
    setup             = setup,
    data_hash         = .get_outcome_hash(object),
    dependency_blocks = .known_v_dependency_blocks(
      object[["data"]],
      nrow(dat)
    )
  )
  expect_true(target[["known_r"]])
  expect_equal(target[["known_r_blocks"]], "estimate")
  expect_true(grepl("BayesTools-prepared row variance",
                    target[["known_r_semantics"]], fixed = TRUE))
  expect_equal(target[["random_effect_representation"]], "marginalized")
  expect_equal(target[["latent_effect_handling"]], "integrated")
})


test_that("brma.mv known R row multipliers enter known-V backends consistently", {

  dat <- data.frame(
    yi       = c(0.10, 0.20, 0.30),
    estimate = factor(c("e1", "e2", "e3"), levels = c("e1", "e2", "e3")),
    x        = c(0, 1, 2)
  )
  R <- diag(c(4, 9, 16))
  dimnames(R) <- list(levels(dat[["estimate"]]), levels(dat[["estimate"]]))
  V <- matrix(
    c(
      0.04, 0.01, 0.005,
      0.01, 0.09, 0.002,
      0.005, 0.002, 0.16
    ),
    nrow  = 3,
    byrow = TRUE
  )

  latent <- brma.mv(
    yi                        = yi,
    V                         = V,
    random                    = ~ 1 | estimate,
    R                         = R,
    Rscale                    = "none",
    data                      = dat,
    known_v_parameterization  = "latent",
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  latent_term   <- .data_marginalized_random_effects(latent[["data"]])[[1L]]
  latent_syntax <- .create_model_syntax(latent[["data"]], latent[["priors"]])
  expect_match(latent_syntax, "sampling_var\\[i\\]")
  expect_match(latent_syntax, latent_term[["row_multiplier_name"]],
               fixed = TRUE)

  whitened <- brma.mv(
    yi                        = yi,
    V                         = V,
    random                    = ~ 1 | estimate,
    R                         = R,
    Rscale                    = "cor",
    data                      = dat,
    known_v_parameterization  = "whitened",
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  whitened_term   <- .data_marginalized_random_effects(whitened[["data"]])[[1L]]
  whitened_syntax <- .create_model_syntax(whitened[["data"]], whitened[["priors"]])
  expect_false(.marginalized_random_effects_row_varying(list(whitened_term)))
  expect_equal(unname(whitened_term[["row_multiplier"]]), c(1, 1, 1))
  expect_match(whitened_syntax, "whitening_var_1\\[j\\]")
  expect_match(whitened_syntax, whitened_term[["row_multiplier_name"]],
               fixed = TRUE)

  expect_error(
    brma.mv(
      yi                        = yi,
      V                         = V,
      random                    = ~ 1 | estimate,
      R                         = R,
      Rscale                    = "none",
      data                      = dat,
      known_v_parameterization  = "whitened",
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_priors               = TRUE
    ),
    "row-varying marginalized estimate-level random effects"
  )

  expect_error(
    brma.mv(
      yi                        = yi,
      V                         = diag(0.01, nrow(dat)),
      scale                     = ~ x,
      random                    = ~ 1 | estimate,
      R                         = R,
      Rscale                    = "none",
      data                      = dat,
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_priors               = TRUE
    ),
    "row-indexed external SD sources"
  )
})


test_that("brma.mv adds whitened known-V likelihood to fit syntax", {

  object <- brma.mv(
    yi                        = c(0.10, 0.20),
    V                         = matrix(c(0.04, 0.03, 0.03, 0.09), nrow = 2),
    known_v_parameterization  = "whitened",
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )

  fit_priors <- .create_fit_priors(object[["data"]], object[["priors"]])
  fit_data   <- .create_fit_data(object[["data"]], object[["priors"]])
  syntax     <- .create_model_syntax(object[["data"]], object[["priors"]])

  expect_false("sampling_z" %in% names(fit_priors))
  expect_true(all(c("whitening_y_1", "whitening_var_1", "whitening_matrix_1") %in%
                  names(fit_data)))
  expect_false(any(c("whitening_y", "whitening_var", "whitening_matrix") %in%
                   names(fit_data)))
  expect_false("sampling_var" %in% names(fit_data))
  expect_false("sampling_B" %in% names(fit_data))
  expect_match(syntax, "mu_observed", fixed = TRUE)
  expect_match(syntax, "whitening_mu", fixed = TRUE)
  expect_match(syntax, "whitening_y", fixed = TRUE)
  expect_false(grepl("sampling_dependency", syntax, fixed = TRUE))
})


test_that("brma.mv prepares block-MVN known-V backend", {

  V <- matrix(
    c(
      0.04, 0.03, 0.00,
      0.03, 0.09, 0.00,
      0.00, 0.00, 0.16
    ),
    nrow = 3,
    byrow = TRUE
  )

  object <- brma.mv(
    yi                       = c(0.10, 0.20, 0.30),
    V                        = V,
    known_v_parameterization = "block_mvn",
    measure                  = "GEN",
    prior_unit_information_sd = 1,
    only_data                = TRUE
  )

  known_V <- attr(object[["data"]], "known_V_data")
  expect_equal(known_V[["parameterization"]], "block_mvn")
  expect_true(known_V[["correlated"]])
  expect_equal(known_V[["rank"]], 0L)
  expect_equal(known_V[["residual_variance"]], diag(V))
  expect_null(known_V[["sampling_factor"]])
  expect_null(known_V[["residual_fraction_requested"]])

  expect_length(known_V[["block_mvn_blocks"]], 1L)
  expect_equal(known_V[["block_mvn_blocks"]][[1]][["index"]], 1:2)
  expect_equal(
    known_V[["block_mvn_blocks"]][[1]][["v_lower"]],
    V[1:2, 1:2][lower.tri(V[1:2, 1:2], diag = TRUE)]
  )
  expect_equal(.known_v_independent_indices(known_V), 3L)
})


test_that("brma.mv adds block-MVN known-V likelihood to fit syntax", {

  object <- brma.mv(
    yi                       = c(0.10, 0.20, 0.30),
    V                        = matrix(
      c(
        0.04, 0.03, 0.00,
        0.03, 0.09, 0.00,
        0.00, 0.00, 0.16
      ),
      nrow = 3,
      byrow = TRUE
    ),
    known_v_parameterization = "block_mvn",
    measure                  = "GEN",
    prior_unit_information_sd = 1,
    only_priors              = TRUE
  )

  fit_priors <- .create_fit_priors(object[["data"]], object[["priors"]])
  fit_data   <- .create_fit_data(object[["data"]], object[["priors"]])
  syntax     <- .create_model_syntax(object[["data"]], object[["priors"]])

  expect_false("sampling_z" %in% names(fit_priors))
  expect_false(any(c("sampling_var", "sampling_B") %in% names(fit_data)))
  expect_false(any(c("whitening_y", "whitening_var", "whitening_matrix") %in% names(fit_data)))
  expect_true(all(c(
    "known_v_independent_y", "known_v_independent_var",
    "known_v_y_1", "known_v_lower_1"
  ) %in% names(fit_data)))
  expect_equal(fit_data[["known_v_y_1"]], c(0.10, 0.20))
  expect_equal(length(fit_data[["known_v_lower_1"]]), 3L)

  expect_match(syntax, "mu_observed", fixed = TRUE)
  expect_match(syntax, "tau2_observed", fixed = TRUE)
  expect_match(syntax, "dknown_v_mnorm", fixed = TRUE)
  expect_match(syntax, "known_v_y_1\\[1:2\\] ~ dknown_v_mnorm")
  expect_match(syntax, "known_v_independent_y\\[j\\] ~ dnorm")
  expect_false(grepl("sampling_dependency", syntax, fixed = TRUE))
})


test_that("brma.mv omits scalar known-V data for all-correlated block-MVN fits", {

  object <- brma.mv(
    yi                       = c(0.10, 0.20),
    V                        = matrix(c(0.04, 0.03, 0.03, 0.09), nrow = 2),
    known_v_parameterization = "block_mvn",
    measure                  = "GEN",
    prior_unit_information_sd = 1,
    only_priors              = TRUE
  )

  fit_data <- .create_fit_data(object[["data"]], object[["priors"]])
  syntax   <- .create_model_syntax(object[["data"]], object[["priors"]])

  expect_false("yi" %in% names(fit_data))
  expect_false("known_v_var" %in% names(fit_data))
  expect_true(all(c("known_v_y_1", "known_v_lower_1") %in% names(fit_data)))
  expect_false(grepl("yi\\[", syntax))
  expect_match(syntax, "known_v_y_1\\[1:2\\] ~ dknown_v_mnorm")
})


test_that("brma.mv block-MVN packing handles mixed-sign covariance blocks", {

  R <- matrix(
    c(
      1.00, -0.65,  0.45,
     -0.65,  1.00, -0.40,
      0.45, -0.40,  1.00
    ),
    nrow = 3,
    byrow = TRUE
  )
  V_block <- diag(c(0.18, 0.24, 0.30)) %*% R %*% diag(c(0.18, 0.24, 0.30))
  V_block <- (V_block + t(V_block)) / 2
  V <- rbind(
    cbind(V_block, matrix(0, 3, 1)),
    cbind(matrix(0, 1, 3), matrix(0.025, 1, 1))
  )

  object <- brma.mv(
    yi                        = c(0.10, 0.20, 0.30, 0.40),
    V                         = V,
    known_v_parameterization  = "block_mvn",
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_data                 = TRUE
  )

  blocks <- attr(object[["data"]], "known_V_data")[["block_mvn_blocks"]]

  expect_length(blocks, 1L)
  expect_equal(blocks[[1]][["index"]], 1:3)
  expect_equal(blocks[[1]][["v_lower"]], V_block[lower.tri(V_block, diag = TRUE)])
  expect_equal(
    .known_v_independent_indices(.data_known_v_data(object[["data"]])),
    4L
  )
})


test_that("brma.mv block-MVN accepts near-singular positive known-V blocks", {

  rho <- 0.99995
  R   <- matrix(rho, nrow = 3, ncol = 3)
  diag(R) <- 1
  V <- diag(c(0.10, 0.20, 0.40)) %*% R %*% diag(c(0.10, 0.20, 0.40))

  object <- brma.mv(
    yi                        = c(0.10, 0.20, 0.30),
    V                         = V,
    known_v_parameterization  = "block_mvn",
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_data                 = TRUE
  )

  known_V <- attr(object[["data"]], "known_V_data")

  expect_equal(known_V[["block_mvn_blocks"]][[1]][["index"]], 1:3)
  expect_true(all(is.finite(known_V[["block_mvn_blocks"]][[1]][["v_lower"]])))
  expect_gt(known_V[["diagnostics"]][["min_eigenvalue"]], 0)

  rho_tight <- 1 - 1e-10
  V_tight <- 0.04 * matrix(c(
    1, rho_tight,
    rho_tight, 1
  ), nrow = 2, byrow = TRUE)

  tight_object <- brma.mv(
    yi                        = c(0.10, 0.20),
    V                         = V_tight,
    known_v_parameterization  = "block_mvn",
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_data                 = TRUE
  )
  tight_known_V <- attr(tight_object[["data"]], "known_V_data")

  expect_equal(tight_known_V[["block_mvn_blocks"]][[1]][["index"]], 1:2)
  expect_gt(tight_known_V[["diagnostics"]][["min_eigenvalue"]], 0)
})


test_that("brma.mv allows block-MVN known-V scale models", {

  dat <- data.frame(
    yi = c(0.10, 0.20, 0.30),
    x  = c(0, 1, 2)
  )

  expect_silent(
    object <- brma.mv(
      yi                       = yi,
      V                        = matrix(
        c(
          0.04, 0.03, 0.00,
          0.03, 0.09, 0.00,
          0.00, 0.00, 0.16
        ),
        nrow = 3,
        byrow = TRUE
      ),
      scale                    = ~ x,
      data                     = dat,
      known_v_parameterization = "block_mvn",
      measure                  = "GEN",
      prior_unit_information_sd = 1,
      only_priors              = TRUE
    )
  )

  syntax <- .create_model_syntax(object[["data"]], object[["priors"]])
  expect_match(syntax, "tau\\[i\\] = exp\\(log_tau\\[i\\]\\)")
  expect_match(syntax, "tau2_observed\\[i\\] = pow\\(tau\\[i\\],2\\)")
  expect_match(syntax, "dknown_v_mnorm", fixed = TRUE)

  random_object <- brma.mv(
    yi                       = yi,
    V                        = matrix(
      c(
        0.04, 0.03, 0.00,
        0.03, 0.09, 0.00,
        0.00, 0.00, 0.16
      ),
      nrow = 3,
      byrow = TRUE
    ),
    scale                    = ~ x,
    random                   = ~ 1 | study,
    data                     = data.frame(dat, study = c("a", "a", "b")),
    known_v_parameterization = "block_mvn",
    measure                  = "GEN",
    prior_unit_information_sd = 1,
    only_priors              = TRUE
  )
  random_syntax <- .create_model_syntax(random_object[["data"]], random_object[["priors"]])
  random_design <- .fitted_formula_design(random_object, "mu", required = TRUE)
  expect_match(random_syntax, "tau\\[i\\] = exp\\(log_tau\\[i\\]\\)")
  expect_match(random_syntax, "tau2_observed\\[i\\] = 0")
  expect_match(random_syntax, "dknown_v_mnorm", fixed = TRUE)
  expect_equal(
    random_design[["random_effects"]][[1]][["sd_binding"]][["source"]][["shape"]],
    "row"
  )
})


test_that("brma.mv rejects unsupported whitened known-V combinations", {

  dat <- data.frame(
    yi = c(0.10, 0.20),
    x  = c(0, 1)
  )

  expect_error(
    brma.mv(
      yi                        = yi,
      V                         = matrix(c(0.04, 0.03, 0.03, 0.09), nrow = 2),
      scale                     = ~ x,
      data                      = dat,
      known_v_parameterization  = "whitened",
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_data                 = TRUE
    ),
    "scale regression"
  )
  expect_error(
    .known_v_prepare(
      V                                   = matrix(c(0.04, 0.03, 0.03, 0.09), nrow = 2),
      keep_rows                           = c(TRUE, TRUE),
      known_v_parameterization            = "whitened",
      known_v_residual_fraction           = NULL,
      known_v_residual_fraction_specified = FALSE,
      known_v_is_scale                    = TRUE
    ),
    "scale regression"
  )
})


test_that("known-V preparation can suppress duplicate singular warnings", {

  V <- matrix(0.04, nrow = 2L, ncol = 2L)

  expect_warning(
    .known_v_prepare(
      V                                   = V,
      keep_rows                           = c(TRUE, TRUE),
      known_v_parameterization            = "whitened",
      known_v_residual_fraction           = NULL,
      known_v_residual_fraction_specified = FALSE,
      warn_singular                       = TRUE
    ),
    "rank-deficient correlation structure"
  )
  expect_silent(
    .known_v_prepare(
      V                                   = V,
      keep_rows                           = c(TRUE, TRUE),
      known_v_parameterization            = "whitened",
      known_v_residual_fraction           = NULL,
      known_v_residual_fraction_specified = FALSE,
      warn_singular                       = FALSE
    )
  )
})
