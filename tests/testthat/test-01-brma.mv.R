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
  dat   <- smoke[["dat"]]
  V     <- smoke[["V"]]
  args  <- .brma_mv_fit_args()

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
    V                         = V,
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
  expect_false(any(grepl("Known-V backend", block_random_output,
                         fixed = TRUE)))
  expect_false(any(grepl("Random Components", block_random_output,
                         fixed = TRUE)))
  expect_true(any(block_random_output == "Random"))
  expect_true(any(grepl("sd(", block_random_output, fixed = TRUE)))
  fit_block_random <- suppressWarnings(add_loo(fit_block_random))
  expect_s3_class(fit_block_random[["loo"]][["estimate"]], "loo")
  save_fit(
    "brma.mv_block_mvn_random",
    fit_block_random,
    info = list(data = dat, V = V, random = "estimate")
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
