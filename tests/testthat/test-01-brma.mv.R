context("Model fitting for brma.mv")

source(testthat::test_path("common-functions.R"))
skip_on_cran()
skip_refit_if_cached("brma.mv")

.brma_mv_smoke_data <- function() {

  dat <- data.frame(
    yi     = c(0.08, 0.13, 0.18, 0.20, 0.01, 0.05),
    study  = rep(c("s1", "s2", "s3"), each = 2),
    effect = rep(c("a", "b"), 3),
    x      = c(0, 1, 0, 1, 0, 1),
    z      = c(-1, -1, 0, 0, 1, 1)
  )
  V <- matrix(0, nrow = 6, ncol = 6)
  block <- matrix(c(0.04, 0.018, 0.018, 0.05), nrow = 2)
  for (g in seq_len(3)) {
    index <- ((g - 1L) * 2L + 1L):(g * 2L)
    V[index, index] <- block
  }

  list(dat = dat, V = V)
}

.brma_mv_fit_args <- function() {

  list(
    measure                   = "GEN",
    chains                    = 2,
    sample                    = 120,
    burnin                    = 60,
    adapt                     = 500,
    seed                      = 1,
    silent                    = TRUE,
    prior_unit_information_sd = 1,
    convergence_checks        = set_convergence_checks(
      max_Rhat = NULL,
      min_ESS  = NULL
    )
  )
}

.brma_mv_metafor_fit_args <- function(seed) {

  list(
    chains             = 2,
    sample             = 1200,
    burnin             = 1200,
    adapt              = 2000,
    seed               = seed,
    silent             = TRUE,
    convergence_checks = set_convergence_checks(
      max_Rhat = NULL,
      min_ESS  = NULL
    )
  )
}

.brma_mv_add_cached_diagnostics <- function(fit, seed) {

  set.seed(seed)
  fit <- add_marglik(fit)
  fit <- suppressWarnings(add_loo(fit))

  return(fit)
}

.brma_mv_model_syntax <- function(fit) {

  syntax <- fit[["fit"]][["model_syntax"]]
  if (is.null(syntax)) {
    syntax <- attr(fit[["fit"]], "model_syntax", exact = TRUE)
  }
  if (is.null(syntax)) {
    syntax <- .create_model_syntax(fit[["data"]], fit[["priors"]])
  }

  syntax
}

test_that("brma.mv fits known-V backend smoke models", {

  smoke <- .brma_mv_smoke_data()
  dat      <- smoke[["dat"]]
  V        <- smoke[["V"]]
  V_blocks <- lapply(seq_len(3L), function(g) {
    index <- ((g - 1L) * 2L + 1L):(g * 2L)
    V[index, index, drop = FALSE]
  })
  args     <- .brma_mv_fit_args()

  fit_latent <- brma.mv(
    yi                        = yi,
    V                         = V,
    data                      = dat,
    known_v_parameterization  = "latent",
    measure                   = args[["measure"]],
    chains                    = args[["chains"]],
    sample                    = args[["sample"]],
    burnin                    = args[["burnin"]],
    adapt                     = 1000,
    seed                      = args[["seed"]],
    silent                    = args[["silent"]],
    prior_unit_information_sd = args[["prior_unit_information_sd"]],
    convergence_checks        = args[["convergence_checks"]]
  )
  latent_syntax <- .brma_mv_model_syntax(fit_latent)
  expect_s3_class(fit_latent, "brma.mv")
  expect_true("mu" %in% rownames(fit_latent[["summary"]]))
  expect_match(latent_syntax, "sampling_z", fixed = TRUE)
  fit_latent <- suppressWarnings(add_loo(fit_latent))
  expect_s3_class(fit_latent[["loo"]][["estimate"]], "loo")
  save_fit("brma.mv_latent", fit_latent, info = list(data = dat, V = V))

  fit_whitened <- brma.mv(
    yi                        = yi,
    V                         = V_blocks,
    data                      = dat,
    known_v_parameterization  = "whitened",
    measure                   = args[["measure"]],
    chains                    = args[["chains"]],
    sample                    = args[["sample"]],
    burnin                    = args[["burnin"]],
    adapt                     = 1000,
    seed                      = args[["seed"]],
    silent                    = args[["silent"]],
    prior_unit_information_sd = args[["prior_unit_information_sd"]],
    convergence_checks        = args[["convergence_checks"]]
  )
  whitened_syntax <- .brma_mv_model_syntax(fit_whitened)
  expect_s3_class(fit_whitened, "brma.mv")
  expect_true("mu" %in% rownames(fit_whitened[["summary"]]))
  expect_match(whitened_syntax, "whitening_y", fixed = TRUE)
  whitened_known_V <- .data_known_v_data(fit_whitened[["data"]])
  expect_identical(whitened_known_V[["storage"]], "blocks")
  expect_null(whitened_known_V[["V"]])
  expect_null(whitened_known_V[["whitening_matrix"]])
  fit_whitened <- suppressWarnings(add_loo(fit_whitened))
  expect_s3_class(fit_whitened[["loo"]][["estimate"]], "loo")
  save_fit("brma.mv_whitened", fit_whitened, info = list(data = dat, V = V))

  fit_block_mvn <- brma.mv(
    yi                        = yi,
    V                         = V,
    data                      = dat,
    known_v_parameterization  = "block_mvn",
    measure                   = args[["measure"]],
    chains                    = args[["chains"]],
    sample                    = args[["sample"]],
    burnin                    = args[["burnin"]],
    adapt                     = args[["adapt"]],
    seed                      = args[["seed"]],
    silent                    = args[["silent"]],
    prior_unit_information_sd = args[["prior_unit_information_sd"]],
    convergence_checks        = args[["convergence_checks"]]
  )
  block_syntax <- .brma_mv_model_syntax(fit_block_mvn)
  expect_s3_class(fit_block_mvn, "brma.mv")
  expect_true("mu" %in% rownames(fit_block_mvn[["summary"]]))
  expect_match(block_syntax, "dknown_v_mnorm", fixed = TRUE)
  fit_block_mvn <- suppressWarnings(add_loo(fit_block_mvn))
  expect_s3_class(fit_block_mvn[["loo"]][["estimate"]], "loo")
  save_fit("brma.mv_block_mvn", fit_block_mvn, info = list(data = dat, V = V))

  fit_block_random <- brma.mv(
    yi                        = yi,
    V                         = V,
    data                      = transform(dat, estimate = seq_len(nrow(dat))),
    random                    = ~ 1 | estimate,
    known_v_parameterization  = "block_mvn",
    measure                   = args[["measure"]],
    chains                    = args[["chains"]],
    sample                    = args[["sample"]],
    burnin                    = args[["burnin"]],
    adapt                     = args[["adapt"]],
    seed                      = args[["seed"]],
    silent                    = args[["silent"]],
    prior_unit_information_sd = args[["prior_unit_information_sd"]],
    convergence_checks        = args[["convergence_checks"]]
  )
  block_random_syntax <- .brma_mv_model_syntax(fit_block_random)
  expect_s3_class(fit_block_random, "brma.mv")
  expect_true(.is_random(fit_block_random))
  expect_match(block_random_syntax, "dknown_v_mnorm", fixed = TRUE)
  expect_match(block_random_syntax, "mu__xREx__estimate_intercept", fixed = TRUE)
  block_random_summary <- summary(
    fit_block_random,
    include_mcmc_diagnostics = FALSE
  )
  block_random_output <- capture.output(print(block_random_summary))
  expect_true(length(block_random_summary[["estimates_random"]]) > 0L)
  expect_true(any(grepl("Bayesian Multivariate", block_random_output,
                        fixed = TRUE)))
  expect_true(any(grepl("Known-V backend", block_random_output,
                        fixed = TRUE)))
  expect_false(any(grepl("Random Components", block_random_output,
                         fixed = TRUE)))
  expect_true(any(block_random_output == "Random"))
  expect_true(any(grepl("sd(", block_random_output, fixed = TRUE)))
  fit_block_random <- suppressWarnings(add_loo(fit_block_random))
  fit_block_random <- suppressWarnings(add_waic(fit_block_random))
  expect_s3_class(fit_block_random[["loo"]][["estimate"]], "loo")
  expect_s3_class(fit_block_random[["waic"]][["estimate"]], "waic")
  fit_block_random_cache <- fit_block_random
  fit_block_random_cache[["waic"]] <- NULL
  save_fit(
    "brma.mv_block_mvn_random",
    fit_block_random_cache,
    info = list(data = dat, V = V, random = "estimate")
  )

  fit_block_random_sampled <- brma.mv(
    yi                         = yi,
    V                          = V,
    data                       = transform(dat, estimate = seq_len(nrow(dat))),
    random                     = ~ 1 | estimate,
    marginalize_estimate_level = FALSE,
    known_v_parameterization   = "block_mvn",
    measure                    = args[["measure"]],
    chains                     = args[["chains"]],
    sample                     = args[["sample"]],
    burnin                     = args[["burnin"]],
    adapt                      = args[["adapt"]],
    seed                       = args[["seed"]],
    silent                     = args[["silent"]],
    prior_unit_information_sd  = args[["prior_unit_information_sd"]],
    convergence_checks         = args[["convergence_checks"]]
  )
  fit_block_random_sampled <- suppressWarnings(add_loo(
    fit_block_random_sampled
  ))
  fit_block_random_sampled <- suppressWarnings(add_waic(
    fit_block_random_sampled
  ))

  marginalized_target <- attr(
    log_lik(fit_block_random),
    "RoBMA_target",
    exact = TRUE
  )
  sampled_target <- attr(
    log_lik(fit_block_random_sampled),
    "RoBMA_target",
    exact = TRUE
  )
  comparison_fields <- c(
    "unit", "conditioning_depth", "target", "data_hash"
  )
  expect_identical(
    marginalized_target[comparison_fields],
    sampled_target[comparison_fields]
  )
  expect_identical(
    marginalized_target[["estimate_level_random"]],
    "marginalized"
  )
  expect_identical(sampled_target[["random_effects"]], "conditioned")

  loo_comparison <- loo_compare(
    fit_block_random,
    fit_block_random_sampled
  )
  waic_comparison <- loo_compare(
    waic(fit_block_random),
    waic(fit_block_random_sampled)
  )
  expect_true(all(is.finite(loo_comparison)))
  expect_true(all(is.finite(waic_comparison)))
  expect_lt(
    abs(loo_comparison[2L, "elpd_diff"]),
    0.25
  )
  expect_lt(
    abs(waic_comparison[2L, "elpd_diff"]),
    0.25
  )

  fit_block_random_sampled_cache <- fit_block_random_sampled
  fit_block_random_sampled_cache[["waic"]] <- NULL
  save_fit(
    "brma.mv_block_mvn_random_sampled",
    fit_block_random_sampled_cache,
    info = list(data = dat, V = V, random = "estimate")
  )

  study_R <- matrix(
    c(
      4.0, 1.0, 0.5,
      1.0, 9.0, 2.0,
      0.5, 2.0, 16.0
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("s1", "s2", "s3"), c("s1", "s2", "s3"))
  )
  fit_known_R_metafor <- metafor::rma.mv(
    yi      = yi,
    V       = V,
    mods    = ~ x + z,
    random  = ~ 1 | study,
    R       = list(study = study_R),
    Rscale  = "none",
    data    = dat,
    method  = "ML"
  )
  fit_known_R <- brma.mv(
    yi                        = yi,
    V                         = V,
    mods                      = ~ x + z,
    data                      = dat,
    random                    = ~ 1 | study,
    R                         = list(study = study_R),
    Rscale                    = list(study = "none"),
    known_v_parameterization  = "block_mvn",
    measure                   = args[["measure"]],
    chains                    = args[["chains"]],
    sample                    = args[["sample"]],
    burnin                    = args[["burnin"]],
    adapt                     = args[["adapt"]],
    seed                      = 11,
    silent                    = args[["silent"]],
    prior_unit_information_sd = args[["prior_unit_information_sd"]],
    convergence_checks        = args[["convergence_checks"]]
  )
  known_R_syntax <- .brma_mv_model_syntax(fit_known_R)
  known_R_term <- .fitted_formula_design(fit_known_R, "mu")[["random_effects"]][[1L]]
  known_R_summary <- summary(
    fit_known_R,
    include_mcmc_diagnostics = FALSE
  )
  known_R_output <- capture.output(print(known_R_summary))
  expect_s3_class(fit_known_R, "brma.mv")
  expect_true(.is_mods(fit_known_R))
  expect_true(.is_random(fit_known_R))
  expect_equal(known_R_term[["compile_mode"]], "sampled")
  expect_equal(known_R_term[["group_covariance"]][["scale"]], "none")
  expect_match(known_R_syntax, "_xRE_GROUP_Zx", fixed = TRUE)
  expect_true(any(grepl("sd_multiplier(", known_R_output, fixed = TRUE)))
  fit_known_R <- suppressWarnings(add_loo(fit_known_R))
  expect_s3_class(fit_known_R[["loo"]][["estimate"]], "loo")
  save_fit(
    "brma.mv_block_mvn_known_R",
    fit_known_R,
    info = list(
      metafor = fit_known_R_metafor,
      data    = dat,
      V       = V,
      R       = study_R,
      Rscale  = "none",
      mods    = c("x", "z"),
      random  = "study"
    )
  )

  estimate_dat <- transform(dat, estimate = seq_len(nrow(dat)))
  fit_latent_estimate_scale <- brma.mv(
    yi                        = yi,
    V                         = V,
    data                      = estimate_dat,
    random                    = ~ 1 | estimate,
    scale                     = ~ x,
    known_v_parameterization  = "latent",
    measure                   = args[["measure"]],
    chains                    = args[["chains"]],
    sample                    = args[["sample"]],
    burnin                    = args[["burnin"]],
    adapt                     = 1000,
    seed                      = args[["seed"]],
    silent                    = args[["silent"]],
    prior_unit_information_sd = args[["prior_unit_information_sd"]],
    convergence_checks        = args[["convergence_checks"]]
  )
  latent_estimate_scale_syntax <- .brma_mv_model_syntax(fit_latent_estimate_scale)
  expect_s3_class(fit_latent_estimate_scale, "brma.mv")
  expect_true(.is_random(fit_latent_estimate_scale))
  expect_true(.is_scale(fit_latent_estimate_scale))
  expect_match(latent_estimate_scale_syntax, "sampling_dependency", fixed = TRUE)
  expect_match(latent_estimate_scale_syntax, "sampling_var\\[i\\] \\+ pow\\(tau\\[i\\],2\\)")
  fit_latent_estimate_scale <- suppressWarnings(add_loo(fit_latent_estimate_scale))
  expect_s3_class(fit_latent_estimate_scale[["loo"]][["estimate"]], "loo")
  save_fit(
    "brma.mv_latent_estimate_scale",
    fit_latent_estimate_scale,
    info = list(data = estimate_dat, V = V, random = "estimate", scale = "x")
  )

  fit_block_estimate_scale <- brma.mv(
    yi                        = yi,
    V                         = V,
    data                      = estimate_dat,
    random                    = ~ 1 | estimate,
    scale                     = ~ x,
    known_v_parameterization  = "block_mvn",
    measure                   = args[["measure"]],
    chains                    = args[["chains"]],
    sample                    = args[["sample"]],
    burnin                    = args[["burnin"]],
    adapt                     = args[["adapt"]],
    seed                      = args[["seed"]],
    silent                    = args[["silent"]],
    prior_unit_information_sd = args[["prior_unit_information_sd"]],
    convergence_checks        = args[["convergence_checks"]]
  )
  block_estimate_scale_syntax <- .brma_mv_model_syntax(fit_block_estimate_scale)
  expect_s3_class(fit_block_estimate_scale, "brma.mv")
  expect_true(.is_random(fit_block_estimate_scale))
  expect_true(.is_scale(fit_block_estimate_scale))
  expect_match(block_estimate_scale_syntax, "dknown_v_mnorm", fixed = TRUE)
  expect_match(block_estimate_scale_syntax, "tau2_observed\\[i\\] = pow\\(tau\\[i\\],2\\)")
  fit_block_estimate_scale <- suppressWarnings(add_loo(fit_block_estimate_scale))
  expect_s3_class(fit_block_estimate_scale[["loo"]][["estimate"]], "loo")
  save_fit(
    "brma.mv_block_mvn_estimate_scale",
    fit_block_estimate_scale,
    info = list(data = estimate_dat, V = V, random = "estimate", scale = "x")
  )

  fit_block_random_scale <- brma.mv(
    yi                        = yi,
    V                         = V,
    data                      = dat,
    random                    = ~ 1 | study / effect,
    scale                     = ~ x,
    known_v_parameterization  = "block_mvn",
    measure                   = args[["measure"]],
    chains                    = args[["chains"]],
    sample                    = args[["sample"]],
    burnin                    = args[["burnin"]],
    adapt                     = args[["adapt"]],
    seed                      = args[["seed"]],
    silent                    = args[["silent"]],
    prior_unit_information_sd = args[["prior_unit_information_sd"]],
    convergence_checks        = args[["convergence_checks"]]
  )
  random_scale_syntax <- .brma_mv_model_syntax(fit_block_random_scale)
  expect_s3_class(fit_block_random_scale, "brma.mv")
  expect_true(.summary_random_components_enabled(fit_block_random_scale))
  expect_true(any(grepl("log_tau", rownames(fit_block_random_scale[["summary"]]),
                        fixed = TRUE)))
  expect_match(random_scale_syntax, "dknown_v_mnorm", fixed = TRUE)
  expect_match(random_scale_syntax, "tau\\[i\\] = exp\\(log_tau\\[i\\]\\)")
  random_scale_summary <- summary(
    fit_block_random_scale,
    include_mcmc_diagnostics = FALSE
  )
  random_scale_output <- capture.output(print(random_scale_summary))
  expect_true(length(random_scale_summary[["estimates_random"]]) > 0L)
  expect_false(any(grepl("Random Components", random_scale_output,
                         fixed = TRUE)))
  expect_true(any(grepl("Scale", random_scale_output, fixed = TRUE)))
  expect_true(any(random_scale_output == "Random"))
  expect_true(any(grepl("var_frac(", random_scale_output, fixed = TRUE)))
  fit_block_random_scale <- suppressWarnings(add_loo(fit_block_random_scale))
  expect_s3_class(fit_block_random_scale[["loo"]][["estimate"]], "loo")
  save_fit(
    "brma.mv_block_mvn_random_scale",
    fit_block_random_scale,
    info = list(data = dat, V = V, scale = "x")
  )

  fit_3lvl_scale_total <- brma.mv(
    yi                        = yi,
    V                         = V,
    data                      = dat,
    random                    = list(total = ~ 1 | study / effect),
    scale                     = ~ x,
    known_v_parameterization  = "block_mvn",
    measure                   = args[["measure"]],
    chains                    = args[["chains"]],
    sample                    = args[["sample"]],
    burnin                    = args[["burnin"]],
    adapt                     = args[["adapt"]],
    seed                      = args[["seed"]],
    silent                    = args[["silent"]],
    prior_unit_information_sd = args[["prior_unit_information_sd"]],
    convergence_checks        = args[["convergence_checks"]]
  )
  scale_total_syntax <- .brma_mv_model_syntax(fit_3lvl_scale_total)
  expect_s3_class(fit_3lvl_scale_total, "brma.mv")
  expect_true(.is_random(fit_3lvl_scale_total))
  expect_true(.is_scale(fit_3lvl_scale_total))
  expect_true(.summary_random_components_enabled(fit_3lvl_scale_total))
  expect_equal(
    .data_scale_formula_parameters(fit_3lvl_scale_total[["data"]]),
    c(tau = "log_tau")
  )
  expect_equal(
    .data_scale_formula_sources(fit_3lvl_scale_total[["data"]]),
    c(tau = "tau")
  )
  expect_true(any(grepl("(log_tau) x", rownames(fit_3lvl_scale_total[["summary"]]),
                        fixed = TRUE)))
  expect_match(scale_total_syntax, "tau\\[i\\] = exp\\(log_tau\\[i\\]\\)")
  expect_match(scale_total_syntax, "mu__xRE_ALLOCx_total__weight", fixed = TRUE)
  fit_3lvl_scale_total <- suppressWarnings(add_loo(fit_3lvl_scale_total))
  expect_s3_class(fit_3lvl_scale_total[["loo"]][["estimate"]], "loo")
  save_fit(
    "brma.mv_block_mvn_3lvl_scale_total",
    fit_3lvl_scale_total,
    info = list(data = dat, V = V, random = "study/effect", scale = "total")
  )

  fit_3lvl_scale_top <- brma.mv(
    yi                        = yi,
    V                         = V,
    data                      = dat,
    random                    = list(study = ~ 1 | study, effect = ~ 1 | study:effect),
    scale                     = list(study = ~ x, effect = ~ 1),
    known_v_parameterization  = "block_mvn",
    measure                   = args[["measure"]],
    chains                    = args[["chains"]],
    sample                    = args[["sample"]],
    burnin                    = args[["burnin"]],
    adapt                     = args[["adapt"]],
    seed                      = args[["seed"]],
    silent                    = args[["silent"]],
    prior_unit_information_sd = args[["prior_unit_information_sd"]],
    convergence_checks        = args[["convergence_checks"]]
  )
  scale_top_syntax <- .brma_mv_model_syntax(fit_3lvl_scale_top)
  expect_s3_class(fit_3lvl_scale_top, "brma.mv")
  expect_true(.is_random(fit_3lvl_scale_top))
  expect_true(.is_scale(fit_3lvl_scale_top))
  expect_equal(
    .data_scale_formula_parameters(fit_3lvl_scale_top[["data"]]),
    c(study = "log_tau_study", effect = "log_tau_effect")
  )
  expect_equal(
    .data_scale_formula_sources(fit_3lvl_scale_top[["data"]]),
    c(study = "tau_study", effect = "tau_effect")
  )
  expect_true(any(grepl("(log_tau_study) x", rownames(fit_3lvl_scale_top[["summary"]]),
                        fixed = TRUE)))
  expect_true(any(grepl("(log_tau_effect) intercept", rownames(fit_3lvl_scale_top[["summary"]]),
                        fixed = TRUE)))
  expect_false(any(grepl("(log_tau_effect) x", rownames(fit_3lvl_scale_top[["summary"]]),
                         fixed = TRUE)))
  expect_match(scale_top_syntax, "tau_study\\[i\\] = exp\\(log_tau_study\\[i\\]\\)")
  expect_match(scale_top_syntax, "tau_effect\\[i\\] = exp\\(log_tau_effect\\[i\\]\\)")
  fit_3lvl_scale_top <- suppressWarnings(add_loo(fit_3lvl_scale_top))
  expect_s3_class(fit_3lvl_scale_top[["loo"]][["estimate"]], "loo")
  save_fit(
    "brma.mv_block_mvn_3lvl_scale_top",
    fit_3lvl_scale_top,
    info = list(data = dat, V = V, random = c("study", "study:effect"), scale = "study")
  )

  fit_3lvl_scale_bottom <- brma.mv(
    yi                        = yi,
    V                         = V,
    data                      = dat,
    random                    = list(study = ~ 1 | study, effect = ~ 1 | study:effect),
    scale                     = list(study = ~ 1, effect = ~ x),
    known_v_parameterization  = "block_mvn",
    measure                   = args[["measure"]],
    chains                    = args[["chains"]],
    sample                    = args[["sample"]],
    burnin                    = args[["burnin"]],
    adapt                     = args[["adapt"]],
    seed                      = args[["seed"]],
    silent                    = args[["silent"]],
    prior_unit_information_sd = args[["prior_unit_information_sd"]],
    convergence_checks        = args[["convergence_checks"]]
  )
  scale_bottom_syntax <- .brma_mv_model_syntax(fit_3lvl_scale_bottom)
  expect_s3_class(fit_3lvl_scale_bottom, "brma.mv")
  expect_true(.is_random(fit_3lvl_scale_bottom))
  expect_true(.is_scale(fit_3lvl_scale_bottom))
  expect_equal(
    .data_scale_formula_parameters(fit_3lvl_scale_bottom[["data"]]),
    c(study = "log_tau_study", effect = "log_tau_effect")
  )
  expect_equal(
    .data_scale_formula_sources(fit_3lvl_scale_bottom[["data"]]),
    c(study = "tau_study", effect = "tau_effect")
  )
  expect_true(any(grepl("(log_tau_study) intercept", rownames(fit_3lvl_scale_bottom[["summary"]]),
                        fixed = TRUE)))
  expect_false(any(grepl("(log_tau_study) x", rownames(fit_3lvl_scale_bottom[["summary"]]),
                         fixed = TRUE)))
  expect_true(any(grepl("(log_tau_effect) x", rownames(fit_3lvl_scale_bottom[["summary"]]),
                        fixed = TRUE)))
  expect_match(scale_bottom_syntax, "tau_study\\[i\\] = exp\\(log_tau_study\\[i\\]\\)")
  expect_match(scale_bottom_syntax, "tau_effect\\[i\\] = exp\\(log_tau_effect\\[i\\]\\)")
  fit_3lvl_scale_bottom <- suppressWarnings(add_loo(fit_3lvl_scale_bottom))
  expect_s3_class(fit_3lvl_scale_bottom[["loo"]][["estimate"]], "loo")
  save_fit(
    "brma.mv_block_mvn_3lvl_scale_bottom",
    fit_3lvl_scale_bottom,
    info = list(data = dat, V = V, random = c("study", "study:effect"), scale = "effect")
  )

  fit_block_mods <- brma.mv(
    yi                        = yi,
    V                         = V,
    mods                      = ~ x + z,
    data                      = dat,
    known_v_parameterization  = "block_mvn",
    measure                   = args[["measure"]],
    chains                    = args[["chains"]],
    sample                    = args[["sample"]],
    burnin                    = args[["burnin"]],
    adapt                     = args[["adapt"]],
    seed                      = args[["seed"]],
    silent                    = args[["silent"]],
    prior_unit_information_sd = args[["prior_unit_information_sd"]],
    convergence_checks        = args[["convergence_checks"]]
  )
  expect_s3_class(fit_block_mods, "brma.mv")
  expect_true(.is_mods(fit_block_mods))
  fit_block_mods <- suppressWarnings(add_loo(fit_block_mods))
  expect_s3_class(fit_block_mods[["loo"]][["estimate"]], "loo")
  save_fit(
    "brma.mv_block_mvn_mods",
    fit_block_mods,
    info = list(data = dat, V = V, mods = c("x", "z"))
  )

  fit_block_random_mods_scale <- brma.mv(
    yi                        = yi,
    V                         = V,
    mods                      = ~ x + z,
    data                      = dat,
    random                    = ~ 1 | study / effect,
    scale                     = ~ x,
    known_v_parameterization  = "block_mvn",
    measure                   = args[["measure"]],
    chains                    = args[["chains"]],
    sample                    = args[["sample"]],
    burnin                    = args[["burnin"]],
    adapt                     = 1000,
    seed                      = args[["seed"]],
    silent                    = args[["silent"]],
    prior_unit_information_sd = args[["prior_unit_information_sd"]],
    convergence_checks        = args[["convergence_checks"]]
  )
  expect_s3_class(fit_block_random_mods_scale, "brma.mv")
  expect_true(.is_mods(fit_block_random_mods_scale))
  expect_true(.is_random(fit_block_random_mods_scale))
  expect_true(.is_scale(fit_block_random_mods_scale))
  save_fit(
    "brma.mv_block_mvn_random_mods_scale",
    fit_block_random_mods_scale,
    info = list(data = dat, V = V, mods = c("x", "z"), random = "study/effect", scale = "x")
  )
})


test_that("brma.mv fits structurally regularized singular V", {

  skip_if_fit_not_active("brma.mv_singular_regularized_whitened")

  dat <- data.frame(
    yi    = c(0.08, 0.13, 0.18),
    study = rep("s1", 3)
  )
  sei <- c(0.20, 0.30, 0.40)
  V   <- tcrossprod(sei)
  A         <- matrix(c(1, 0, 1, 0, 1, 1), nrow = 3, ncol = 2)
  V_general <- tcrossprod(A)
  args <- .brma_mv_fit_args()

  prior_zero <- BayesTools::prior(
    distribution = "spike",
    parameters   = list(location = 0)
  )
  prior_positive <- BayesTools::prior(
    distribution = "spike",
    parameters   = list(location = 0.10)
  )
  expect_error(
    suppressWarnings(brma.mv(
      yi                        = yi,
      V                         = V,
      data                      = dat,
      prior_heterogeneity       = prior_zero,
      known_v_parameterization  = "whitened",
      measure                   = args[["measure"]],
      prior_unit_information_sd = args[["prior_unit_information_sd"]],
      chains                    = args[["chains"]],
      sample                    = args[["sample"]],
      burnin                    = args[["burnin"]],
      adapt                     = args[["adapt"]],
      seed                      = args[["seed"]],
      silent                    = args[["silent"]],
      convergence_checks        = args[["convergence_checks"]]
    )),
    "not structurally regularized"
  )
  expect_error(
    suppressWarnings(brma.mv(
      yi                        = yi,
      V                         = V,
      random                    = ~ 1 | study,
      data                      = dat,
      known_v_parameterization  = "block_mvn",
      measure                   = args[["measure"]],
      prior_unit_information_sd = args[["prior_unit_information_sd"]],
      chains                    = args[["chains"]],
      sample                    = args[["sample"]],
      burnin                    = args[["burnin"]],
      adapt                     = args[["adapt"]],
      seed                      = args[["seed"]],
      silent                    = args[["silent"]],
      convergence_checks        = args[["convergence_checks"]]
    )),
    "Sampled random effects change the conditional mean"
  )
  expect_error(
    suppressWarnings(brma.mv(
      yi                        = yi,
      V                         = V_general,
      data                      = dat,
      prior_heterogeneity       = prior_zero,
      known_v_parameterization  = "block_mvn",
      measure                   = args[["measure"]],
      prior_unit_information_sd = args[["prior_unit_information_sd"]],
      chains                    = args[["chains"]],
      sample                    = args[["sample"]],
      burnin                    = args[["burnin"]],
      adapt                     = args[["adapt"]],
      seed                      = args[["seed"]],
      silent                    = args[["silent"]],
      convergence_checks        = args[["convergence_checks"]]
    )),
    "not structurally regularized"
  )

  general_fit <- suppressWarnings(brma.mv(
    yi                        = yi,
    V                         = V_general,
    data                      = dat,
    prior_heterogeneity       = prior_positive,
    known_v_parameterization  = "block_mvn",
    measure                   = args[["measure"]],
    prior_unit_information_sd = args[["prior_unit_information_sd"]],
    chains                    = args[["chains"]],
    sample                    = args[["sample"]],
    burnin                    = args[["burnin"]],
    adapt                     = args[["adapt"]],
    seed                      = args[["seed"]],
    silent                    = args[["silent"]],
    convergence_checks        = args[["convergence_checks"]]
  ))
  expect_s3_class(general_fit, "brma.mv")
  expect_true(isTRUE(general_fit[["fit"]][["has_posterior"]]))

  V_tolerance <- matrix(c(1, 1 + 1e-9, 1 + 1e-9, 1), nrow = 2)
  expect_error(
    suppressWarnings(brma.mv(
      yi                       = yi,
      V                        = V_tolerance,
      data                     = dat[1:2, , drop = FALSE],
      known_v_parameterization = "block_mvn",
      measure                  = args[["measure"]],
      chains                   = args[["chains"]],
      sample                   = args[["sample"]],
      burnin                   = args[["burnin"]],
      adapt                    = args[["adapt"]],
      seed                     = args[["seed"]],
      silent                   = args[["silent"]],
      convergence_checks       = args[["convergence_checks"]]
    )),
    "must be positive semidefinite"
  )

  estimate_dat <- transform(dat, estimate = paste0("e", seq_len(nrow(dat))))
  marginalized_fit <- suppressWarnings(brma.mv(
    yi                        = yi,
    V                         = V,
    random                    = ~ 1 | estimate,
    data                      = estimate_dat,
    prior_heterogeneity       = BayesTools::prior_random(
      estimate = BayesTools::random_block(sd = prior_positive)
    ),
    known_v_parameterization  = "block_mvn",
    measure                   = args[["measure"]],
    prior_unit_information_sd = args[["prior_unit_information_sd"]],
    chains                    = args[["chains"]],
    sample                    = args[["sample"]],
    burnin                    = args[["burnin"]],
    adapt                     = args[["adapt"]],
    seed                      = args[["seed"]],
    silent                    = args[["silent"]],
    convergence_checks        = args[["convergence_checks"]]
  ))
  expect_s3_class(marginalized_fit, "brma.mv")
  expect_true(isTRUE(marginalized_fit[["fit"]][["has_posterior"]]))

  fixed_scale_fit <- suppressWarnings(brma.mv(
    yi                        = yi,
    V                         = V,
    scale                     = ~ 1,
    data                      = dat,
    prior_scale               = list(
      intercept = BayesTools::prior(
        distribution = "spike",
        parameters   = list(location = 0.10)
      )
    ),
    known_v_parameterization  = "block_mvn",
    measure                   = args[["measure"]],
    prior_unit_information_sd = args[["prior_unit_information_sd"]],
    chains                    = args[["chains"]],
    sample                    = args[["sample"]],
    burnin                    = args[["burnin"]],
    adapt                     = args[["adapt"]],
    seed                      = args[["seed"]],
    silent                    = args[["silent"]],
    convergence_checks        = args[["convergence_checks"]]
  ))
  expect_s3_class(fixed_scale_fit, "brma.mv")
  expect_true(isTRUE(fixed_scale_fit[["fit"]][["has_posterior"]]))

  fits <- list()
  for (parameterization in c("whitened", "block_mvn")) {
    fit <- suppressWarnings(brma.mv(
      yi                        = yi,
      V                         = V,
      data                      = dat,
      prior_heterogeneity       = prior_positive,
      known_v_parameterization  = parameterization,
      measure                   = args[["measure"]],
      prior_unit_information_sd = args[["prior_unit_information_sd"]],
      chains                    = args[["chains"]],
      sample                    = args[["sample"]],
      burnin                    = args[["burnin"]],
      adapt                     = args[["adapt"]],
      seed                      = args[["seed"]],
      silent                    = args[["silent"]],
      convergence_checks        = args[["convergence_checks"]]
    ))
    expect_s3_class(fit, "brma.mv")
    expect_true(isTRUE(fit[["fit"]][["has_posterior"]]))
    fits[[parameterization]] <- suppressWarnings(add_loo(fit))
  }

  save_fit(
    "brma.mv_singular_regularized_whitened",
    fits[["whitened"]],
    info = list(data = dat, V = V)
  )
  save_fit(
    "brma.mv_singular_regularized_block_mvn",
    fits[["block_mvn"]],
    info = list(data = dat, V = V)
  )
})

test_that("brma.mv fits fixed-effect known-V model with random = NULL", {

  skip_if_not_installed("metafor")

  smoke <- .brma_mv_smoke_data()
  dat   <- smoke[["dat"]]
  V     <- smoke[["V"]]
  args  <- .brma_mv_fit_args()

  fit_metafor <- metafor::rma.mv(
    yi,
    V      = V,
    random = NULL,
    data   = dat,
    method = "ML"
  )
  fit_brma <- brma.mv(
    yi                        = yi,
    V                         = V,
    data                      = dat,
    random                    = NULL,
    prior_heterogeneity       = BayesTools::prior(
      distribution = "spike",
      parameters   = list(location = 0)
    ),
    known_v_parameterization  = "block_mvn",
    measure                   = args[["measure"]],
    chains                    = args[["chains"]],
    sample                    = args[["sample"]],
    burnin                    = args[["burnin"]],
    adapt                     = args[["adapt"]],
    seed                      = args[["seed"]],
    silent                    = args[["silent"]],
    prior_unit_information_sd = args[["prior_unit_information_sd"]],
    convergence_checks        = args[["convergence_checks"]]
  )
  fixed_syntax <- .brma_mv_model_syntax(fit_brma)
  expect_s3_class(fit_brma, "brma.mv")
  expect_false(.is_random(fit_brma))
  expect_false(any(grepl("__xRE", rownames(fit_brma[["summary"]]), fixed = TRUE)))
  expect_match(fixed_syntax, "dknown_v_mnorm", fixed = TRUE)
  fit_brma <- suppressWarnings(add_loo(fit_brma))
  expect_s3_class(fit_brma[["loo"]][["estimate"]], "loo")
  save_fit(
    "brma.mv_block_mvn_fixed_random_null",
    fit_brma,
    info = list(
      metafor = fit_metafor,
      data    = dat,
      V       = V
    )
  )
})

test_that("brma.mv fits v14 metafor parity models", {

  skip_if_not_certification("These high-draw parity fits are release evidence.")
  skip_if_not_installed("metadat")
  skip_if_not_installed("metafor")

  data("dat.konstantopoulos2011", package = "metadat")
  dat_konstantopoulos <- dat.konstantopoulos2011
  fit1_metafor <- metafor::rma.mv(
    yi,
    vi,
    random = ~ school | district,
    struct = "CS",
    data   = dat_konstantopoulos
  )
  args <- .brma_mv_metafor_fit_args(seed = 1)
  fit1_brma <- brma.mv(
    yi                 = yi,
    V                  = vi,
    measure            = "SMD",
    random             = ~ cs(school | district),
    data               = dat_konstantopoulos,
    chains             = args[["chains"]],
    sample             = args[["sample"]],
    burnin             = args[["burnin"]],
    adapt              = args[["adapt"]],
    seed               = args[["seed"]],
    silent             = args[["silent"]],
    convergence_checks = args[["convergence_checks"]]
  )
  expect_s3_class(fit1_brma, "brma.mv")
  fit1_brma <- .brma_mv_add_cached_diagnostics(fit1_brma, seed = 101)
  save_fit(
    "brma.mv_v14_konstantopoulos2011_cs",
    fit1_brma,
    info = list(
      metafor = fit1_metafor,
      data    = dat_konstantopoulos
    )
  )

  data("dat.assink2016", package = "metadat")
  dat_assink <- dat.assink2016
  V_assink   <- metafor::vcalc(
    vi,
    cluster = study,
    type    = deltype,
    obs     = esid,
    rho     = c(0.7, 0.5),
    data    = dat_assink
  )
  fit2_metafor <- metafor::rma.mv(
    yi,
    V_assink,
    random = ~ 1 | study / esid,
    data   = dat_assink
  )
  args <- .brma_mv_metafor_fit_args(seed = 2)
  fit2_brma <- brma.mv(
    yi                 = yi,
    V                  = V_assink,
    measure            = "SMD",
    random             = ~ 1 | study / esid,
    data               = dat_assink,
    chains             = args[["chains"]],
    sample             = args[["sample"]],
    burnin             = args[["burnin"]],
    adapt              = args[["adapt"]],
    seed               = args[["seed"]],
    silent             = args[["silent"]],
    convergence_checks = args[["convergence_checks"]]
  )
  expect_s3_class(fit2_brma, "brma.mv")
  fit2_brma <- .brma_mv_add_cached_diagnostics(fit2_brma, seed = 102)
  save_fit(
    "brma.mv_v14_assink2016_nested",
    fit2_brma,
    info = list(
      metafor = fit2_metafor,
      data    = dat_assink,
      V       = V_assink
    )
  )

  data("dat.ishak2007", package = "metadat")
  dat_ishak <- dat.ishak2007
  dat_ishak <- reshape(
    dat_ishak,
    direction = "long",
    idvar     = "study",
    v.names   = c("yi", "vi"),
    varying   = list(c(2, 4, 6, 8), c(3, 5, 7, 9))
  )
  dat_ishak <- dat_ishak[order(dat_ishak$study, dat_ishak$time), ]
  dat_ishak <- dat_ishak[!is.na(dat_ishak$yi), ]
  rownames(dat_ishak) <- NULL
  dat_ishak$time_factor <- factor(dat_ishak$time)
  V_ishak <- metafor::vcalc(
    vi,
    cluster = study,
    time1   = time,
    phi     = 0.97,
    data    = dat_ishak
  )
  fit3_metafor <- metafor::rma.mv(
    yi,
    V_ishak,
    mods   = ~ 0 + time_factor,
    random = ~ time | study,
    struct = "HAR",
    data   = dat_ishak
  )
  ishak_unit_information_sd <- estimate_unit_information_sd(
    sei = sqrt(14.3),
    ni  = 15
  )
  ishak_time_prior <- BayesTools::prior_factor(
    distribution = "normal",
    parameters   = list(mean = 0, sd = ishak_unit_information_sd),
    contrast     = "independent"
  )
  args <- .brma_mv_metafor_fit_args(seed = 3)
  # HAR mixes slowly in this sparse repeated-measures example. Keep this
  # reference fit large enough that parity is not determined by sampler noise.
  args[["chains"]] <- 4
  args[["sample"]] <- 10000
  args[["burnin"]] <- 4000
  fit3_brma <- brma.mv(
    yi                        = yi,
    V                         = V_ishak,
    measure                   = "GEN",
    mods                      = ~ 0 + time_factor,
    random                    = ~ har(time | study),
    prior_effect              = NULL,
    prior_mods                = ishak_time_prior,
    prior_unit_information_sd = ishak_unit_information_sd,
    data                      = dat_ishak,
    chains                    = args[["chains"]],
    sample                    = args[["sample"]],
    burnin                    = args[["burnin"]],
    adapt                     = args[["adapt"]],
    seed                      = args[["seed"]],
    silent                    = args[["silent"]],
    convergence_checks        = args[["convergence_checks"]]
  )
  expect_s3_class(fit3_brma, "brma.mv")
  fit3_brma <- .brma_mv_add_cached_diagnostics(fit3_brma, seed = 103)
  save_fit(
    "brma.mv_v14_ishak2007_har",
    fit3_brma,
    info = list(
      metafor = fit3_metafor,
      data    = dat_ishak,
      V       = V_ishak
    )
  )

  data("dat.begg1989", package = "metadat")
  dat_begg <- dat.begg1989
  dat_begg$trt <- relevel(factor(dat_begg$trt), ref = "CMO")
  dat_begg$ni  <- as.numeric(round(with(dat_begg, yi * (1 - yi) / sei^2)))
  fit4_metafor <- metafor::rma.mv(
    yi,
    vi,
    mods   = ~ trt,
    random = list(
      ~ 1 | study,
      ~ trt | study
    ),
    struct = "CS",
    rho    = 0,
    data   = dat_begg,
    method = "ML"
  )
  begg_unit_information_sd <- estimate_unit_information_sd(
    sei = dat_begg$sei,
    ni  = dat_begg$ni
  )
  sd_prior <- BayesTools::prior(
    distribution = "normal",
    parameters   = list(mean = 0, sd = begg_unit_information_sd),
    truncation   = list(lower = 0, upper = Inf)
  )
  rho_zero <- BayesTools::prior(
    distribution = "spike",
    parameters   = list(location = 0)
  )
  cs_zero <- BayesTools::random_covariance(
    rho       = rho_zero,
    rho_scale = "rho"
  )
  args <- .brma_mv_metafor_fit_args(seed = 4)
  fit4_brma <- brma.mv(
    yi                        = yi,
    V                         = vi,
    ni                        = ni,
    measure                   = "GEN",
    mods                      = ~ trt,
    random                    = list(
      study     = ~ 1 | study,
      treatment = ~ cs(trt | study)
    ),
    prior_heterogeneity       = BayesTools::prior_random(
      study     = BayesTools::random_block(sd = sd_prior),
      treatment = BayesTools::random_block(
        sd         = sd_prior,
        covariance = cs_zero
      )
    ),
    prior_unit_information_sd = begg_unit_information_sd,
    data                      = dat_begg,
    chains                    = args[["chains"]],
    sample                    = args[["sample"]],
    burnin                    = args[["burnin"]],
    adapt                     = args[["adapt"]],
    seed                      = args[["seed"]],
    silent                    = args[["silent"]],
    convergence_checks        = args[["convergence_checks"]]
  )
  expect_s3_class(fit4_brma, "brma.mv")
  fit4_brma <- .brma_mv_add_cached_diagnostics(fit4_brma, seed = 104)
  save_fit(
    "brma.mv_v14_begg1989_study_treatment",
    fit4_brma,
    info = list(
      metafor = fit4_metafor,
      data    = dat_begg
    )
  )
})
