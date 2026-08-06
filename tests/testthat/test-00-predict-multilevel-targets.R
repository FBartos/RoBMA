context("Multilevel prediction targets")

test_that("same-data cluster and response targets retain latent cluster draws", {

  n_draws <- 4L
  n_rows  <- 2L

  fixed_location <- matrix(10, nrow = n_draws, ncol = n_rows)
  latent_cluster <- matrix(
    rep(c(-2, -1, 1, 2), n_rows),
    nrow = n_draws,
    ncol = n_rows
  )
  blup_cluster  <- matrix(0.25, nrow = n_draws, ncol = n_rows)
  blup_estimate <- matrix(0.50, nrow = n_draws, ncol = n_rows)

  context <- list(
    object              = list(fit = structure(list(), class = "mock_fit")),
    type                = "estimate",
    same_data           = TRUE,
    new_data            = list(mods = NULL),
    known_V_new         = NULL,
    priors              = list(mods = list()),
    is_mods             = FALSE,
    is_multilevel       = TRUE,
    is_random           = FALSE,
    is_PET              = FALSE,
    is_PEESE            = FALSE,
    is_weightfunction   = FALSE,
    is_weights          = FALSE,
    is_known_v          = FALSE,
    outcome_type        = "norm",
    effect_direction    = "positive",
    posterior_samples   = matrix(0, nrow = n_draws, ncol = 1L),
    outcome_data        = list(
      yi  = c(0.1, -0.2),
      sei = c(0.2, 0.3)
    ),
    fit_data            = list(cluster = c(1L, 1L)),
    K_original          = n_rows,
    K                   = n_rows,
    random_mv           = FALSE,
    bias_adjusted       = FALSE
  )
  scale_state <- list(
    within  = matrix(0.4, nrow = n_draws, ncol = n_rows),
    between = matrix(0.6, nrow = n_draws, ncol = n_rows)
  )

  testthat::local_mocked_bindings(
    .evaluate.brma.mu = function(...) fixed_location,
    .evaluate.brma.cluster_effects = function(...) latent_cluster,
    .evaluate.brma.multilevel_blup.norm = function(...) {
      list(cluster = blup_cluster, estimate = blup_estimate)
    },
    .package = "RoBMA"
  )

  estimate_state <- .predict_brma_location_state(context, scale_state)
  expect_equal(
    estimate_state[["mu"]],
    fixed_location + blup_cluster,
    info = "estimate predictions use the conditional cluster BLUP"
  )
  expect_equal(
    estimate_state[["multilevel_blup"]][["estimate"]],
    blup_estimate
  )

  for (target in c("cluster", "response")) {
    context[["type"]] <- target
    latent_state <- .predict_brma_location_state(context, scale_state)

    expect_null(
      latent_state[["multilevel_blup"]],
      info = paste(target, "does not use the BLUP shortcut")
    )
    expect_equal(
      latent_state[["mu"]],
      fixed_location + latent_cluster,
      info = paste(target, "retains fitted latent cluster draws")
    )
    expect_true(
      stats::var(latent_state[["mu"]][, 1L]) > 0,
      info = paste(target, "retains cluster uncertainty across draws")
    )
  }
})
