# Genuine compiler/fit-contract supplied-row replays; these are transport and
# declared-target regressions, not MCMC or posterior calibration fixtures.
.review_h_empty_condition <- function() {

  dat <- data.frame(yi = c(.1, .2, .15, .25), sei = .2, x = c(-1, 1, -1, 1))
  setup <- review_compile_object(BMA(yi = yi, sei = sei, mods = ~ x, data = dat,
    measure = "GEN", prior_mods = list(intercept = BayesTools::prior("normal", list(.4, .2))),
    prior_heterogeneity = BayesTools::prior("point", list(0)), prior_heterogeneity_null = NULL,
    prior_unit_information_sd = 1, standardize_continuous_predictors = FALSE,
    only_priors = TRUE, silent = TRUE))
  slope <- setup$compiled$mu$prior_list$mu_x
  components <- vapply(slope, function(prior) attr(prior, "component"), character(1L))
  alternative <- which(components == "alternative")
  null <- which(components == "null")
  samples <- cbind(mu_intercept = .2 + .025 * stats::qnorm((seq_len(120) - .5) / 120),
    mu_x = c(rep(0, 60), seq(.01, .15, length.out = 60)),
    mu_x_indicator = rep(c(null, alternative), each = 60), tau = 0)
  review_attach_replay(setup, samples)
}

test_that("E2-1 declared empty marginal cells retain their actual event", {

  object <- .review_h_empty_condition()
  means <- marginal_means(object, bf = FALSE, density_method = "KDE", n_samples = 128)
  cells <- means$inference$conditional$mu_x
  for (parameter in c("intercept", "x")) {
    specs <- .iwmde_marginal_means_specs(means, parameter, "conditional", NULL)
    expect_gt(length(specs), 0L)
    for (spec in specs) {
      cell <- spec$marginal_samples
      declared <- BayesTools::posterior_metadata(cell, "condition")
      expect_identical(spec$condition_key, declared$condition_key)
      expect_identical(spec$condition_event, declared$condition_event)
      expect_identical(spec$prior_density_context, BayesTools::posterior_metadata(cell, "prior_context"))
      if (isTRUE(declared$averaged)) expect_length(spec$conditional, 0L)
    }
  }
  empty <- means$inference$conditional$mu_intercept$intercept
  for (defect in c("missing", "conflict")) {
    changed <- empty
    metadata <- BayesTools::posterior_metadata(changed, "condition")
    if (defect == "missing") metadata$condition_event <- NULL else metadata$condition_key <- "stale"
    BayesTools::posterior_metadata(changed, "condition") <- metadata
    expect_error(.iwmde_sample_condition_metadata(changed, require_child_condition = TRUE), "metadata")
  }
  # A public point reaches estimation; a subsequent diagnostic refusal is kept.
  result <- tryCatch(hypothesis(means, "intercept = 0", density_method = "qCMDE", seed = 19), error = identity)
  expect_false(inherits(result, "error") && grepl("metadata is incomplete", conditionMessage(result), fixed = TRUE))
})
