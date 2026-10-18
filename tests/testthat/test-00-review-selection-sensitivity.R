.review_sampling_sensitivity_replay <- function(case, direction = "positive") {

  model <- if (case == "A") BayesTools::selection_model(weight_rule = "best", group = "paper",
    known_sampling_variance = "condition") else BayesTools::selection_model(weight_rule = "product",
      other_random_effects = "integrate", known_sampling_variance = "condition")
  args <- list(yi = if (case == "A") c(1, -1) else c(.1, .2),
    V = if (case == "A") matrix(c(1, .5, .5, 1), 2L) else diag(2L),
    data = if (case == "A") data.frame(paper = c("a", "b")) else data.frame(study = c("s", "s")),
    random = if (case == "A") NULL else ~1 | study,
    prior_effect = BayesTools::prior("point", list(location = 0)),
    prior_bias = BayesTools::prior_weightfunction("one-sided", .5, BayesTools::wf_fixed(c(1, .2)), model = model),
    effect_direction = direction, measure = "GEN", prior_unit_information_sd = 1, only_priors = TRUE, silent = TRUE)
  if (case == "B") args$prior_heterogeneity <- BayesTools::prior("point", list(location = 1))
  setup <- review_compile_object(do.call(bselmodel.mv, args))
  priors <- unlist(lapply(setup$compiled, `[[`, "prior_list"), recursive = FALSE)
  priors <- c(priors, .create_fit_priors(setup$object$data, setup$object$priors))
  columns <- setdiff(review_prior_columns(priors), "omega")
  if (case == "B") {
    term <- setup$object$formula_design$mu$random_effects[[1L]]
    columns <- c(columns, paste0(term$parameter_stem, "_xRE_Zx[1,1]"))
  }
  samples <- matrix(0, 24L, length(columns), dimnames = list(NULL, columns))
  for (name in names(priors)) if (BayesTools::is.prior.point(priors[[name]])) samples[, name] <- priors[[name]]$parameters$location
  for (block in .selection_sampling_structure(setup$object$data)$latent_blocks) {
    z <- solve(block$B, setup$object$data$outcome$yi[block$index])
    samples[, paste0("sampling_z[", seq.int(block$z_start, block$z_end), "]")] <- matrix(z, 24L, length(z), byrow = TRUE)
  }
  initial <- review_attach_replay(setup, samples)
  generated <- BayesTools::JAGS_evaluate_deterministic(initial$fit, draws = samples)
  for (name in colnames(generated)) {
    if (name %in% colnames(samples)) samples[, name] <- generated[, name] else samples <- cbind(samples, generated[, name, drop = FALSE])
  }
  object <- review_attach_replay(setup, samples)
  list(object = object, samples = samples)
}

test_that("F3-1 public sensitivity preserves separate publication events and whole-context RNG", {

  replay <- .review_sampling_sensitivity_replay("A")
  object <- replay$object
  withr::local_seed(2666)
  before <- .Random.seed; kind <- RNGkind()
  actual <- suppressWarnings(selection_sensitivity_diagnostics(object, max_posterior_samples = 10,
    latent_samples = 64, seed = 2601))
  expect_identical(.Random.seed, before); expect_identical(RNGkind(), kind)
  set.seed(2601)
  metrics <- t(vapply(seq_len(10L), function(i) {
    context <- .selection_sensitivity_covariance_draws(matrix(c(1, .5, .5, 1), 2L), 64L)
    w <- apply(ifelse(context >= 0, 1, .2), 1L, prod)
    c(TV = .5 * mean(abs(w / mean(w) - 1)), ESS = mean(w)^2 / mean(w^2))
  }, numeric(2L)))
  expect_lte(abs(actual$total_variation_median - median(metrics[, "TV"])), 1e-12)
  expect_lte(abs(actual$ess_fraction_median - median(metrics[, "ESS"])), 1e-12)
  expect_identical(attr(actual, "target")$publication_groups$row_blocks, list(1L, 2L))
  expect_identical(actual$integration_relative_mcse_max, 0)
})

test_that("F3-2 declared singular random factors remain available to sensitivity", {

  replay <- .review_sampling_sensitivity_replay("B")
  object <- replay$object
  posterior <- replay$samples[1L, , drop = FALSE]
  scales <- .zplot_tau_samples(object, posterior)
  setup <- list(fit = object$fit, data = object$data, priors = object$priors,
    posterior_samples = posterior, S = 1L, K = 2L,
    tau_within = scales$tau_within, tau_between = scales$tau_between)
  normalizers <- .selection_conditioned_sampling_normalizer(setup, rbind(c(0, 0), c(1, -1)))
  oracle <- c(.52, .04 * pnorm(-1) + .2 * (pnorm(1) - pnorm(-1)) + 1 - pnorm(1))
  expect_lte(max(abs(exp(normalizers$log_mass) - oracle)), 1e-12)
  withr::local_seed(2666); before <- .Random.seed
  actual <- suppressWarnings(selection_sensitivity_diagnostics(object, max_posterior_samples = 10,
    latent_samples = 8, seed = 2601))
  expect_identical(.Random.seed, before)
  expect_true(all(is.finite(as.matrix(actual[, vapply(actual, is.numeric, logical(1)), drop = FALSE]))))
  expect_identical(actual$integration_relative_mcse_max, 0)
  expect_identical(actual$integration_relative_quadrature_error_max, 0)
})

test_that("conditioned sensitivity uses the validated override in its actual state setup", {

  object <- .review_sampling_sensitivity_replay("B")$object
  original_plan <- .data_selection_execution_plan(object$data)
  control <- set_selection_likelihood_control(points_per_scramble = 64L,
    max_points_per_scramble = 128L, scrambles = 3L, seed = 43L, relative_tolerance = .01)
  observed <- list(); producer <- .selection_conditioned_sampling_factors
  testthat::local_mocked_bindings(.selection_conditioned_sampling_factors = function(setup) {
    observed[[length(observed) + 1L]] <<- .data_selection_execution_plan(setup$data)
    producer(setup)
  }, .package = "RoBMA")
  suppressWarnings(selection_sensitivity_diagnostics(object, max_posterior_samples = 10,
    latent_samples = 8, integration_control = control))
  expect_gt(length(observed), 0L)
  for (plan in observed) expect_identical(plan[names(control)], unclass(control))
  expect_identical(.data_selection_execution_plan(object$data), original_plan)
})
