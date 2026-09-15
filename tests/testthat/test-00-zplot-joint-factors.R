test_that("the shared context fallback preserves block weighting", {

  z <- c(-2, 0, 2)
  sei <- c(.5, .7, .8)
  mean <- matrix(c(.1, -.2, .3, .1, -.1, .4), 2L)
  sigma <- matrix(c(.5, .1, 0, .1, .7, 0, 0, 0, .9), 3L)
  prior <- BayesTools::prior_weightfunction("one-sided", .025, BayesTools::wf_fixed(c(1, .2)))
  context <- .selection_spec(list(outcome = list(bias = prior)), rep(0, 3L), sei,
    effect_direction = "positive", signed_data = FALSE)
  context$omega <- matrix(rep(c(1, .2), 2L), 2L, byrow = TRUE)
  context$alpha <- numeric(2L)
  context$phack_kind <- integer(2L)
  context$kernel_mode <- rep(SELKERNEL_NORMAL, 2L)
  context$vector_rule <- integer(2L)
  context$use_normal <- rep(TRUE, 2L)
  result <- .zplot_full_event_context_mixture(z, mean,
    array(rep(sigma, each = 2L), c(2L, 3L, 3L)), array(0, c(2L, 3L, 3L)), sei,
    context, FALSE, set_selection_likelihood_control(), list(row_blocks = list(1:2, 3L)))
  expected <- vapply(z, function(value) vapply(seq_len(2L), function(row) {
    mean(sei * stats::dnorm(value * sei, mean[row, ], sqrt(diag(sigma))))
  }, numeric(1L)), numeric(2L))
  expect_equal(result, expected, tolerance = 1e-12)
})


test_that("product publication inputs preserve the streamed zplot factor path", {

  dat <- data.frame(yi = c(.1, .2, -.1, .3), study = c("a", "a", "b", "b"),
                    esid = 1:4, paper = "all")
  V <- kronecker(diag(2L), matrix(c(.04, .01, .01, .06), 2L))
  samples <- matrix(c(.1, -.2), 2L, dimnames = list(NULL, "mu_intercept"))
  prepared <- lapply(list(selection_model(), selection_model(group = paper)), function(model) {

    prior <- BayesTools::prior_weightfunction("one-sided", .025,
      BayesTools::wf_fixed(c(1, .2)), model = model)
    object <- bselmodel.mv(yi = yi, V = V, random = ~ 1 | study / esid, data = dat,
      prior_bias = prior, measure = "GEN", prior_unit_information_sd = 1,
      prior_heterogeneity = BayesTools::prior_random(
        sd = BayesTools::prior("point", list(location = .4))),
      only_priors = TRUE, silent = TRUE)
    object$fit <- structure(list(), formula_design = object$formula_design,
      prior_list = c(object$formula_design$mu$prior_list,
                     .create_fit_priors(object$data, object$priors)))
    selection <- .selection_spec(object$priors, dat$yi, sqrt(diag(V)), "positive")
    selection$omega <- matrix(c(1, .2), 2L, 2L, byrow = TRUE)
    selection$alpha <- numeric(2L)
    selection$phack_kind <- integer(2L)
    selection$kernel_mode <- rep(SELKERNEL_STEP, 2L)
    selection$vector_rule <- integer(2L)
    selection$use_normal <- rep(FALSE, 2L)
    .zplot_joint_factor_preparation(object, samples, selection, FALSE,
      set_selection_likelihood_control())
  })
  for (result in prepared) {
    expect_false(is.null(result))
    expect_identical(result$execution_plan$row_blocks, list(1:2, 3:4))
    expect_equal(result$integrated$diagonal, matrix(.4^2, 2L, 4L))
    expect_identical(result$retained$ranks, c(1L, 1L))
  }
  expect_identical(prepared[[1L]]$integrated, prepared[[2L]]$integrated)
  expect_identical(prepared[[1L]]$retained, prepared[[2L]]$retained)
})


test_that("streamed contexts keep row weights and current-cell fallback aligned", {

  dat <- data.frame(yi = c(.1, .2, .3), study = c("a", "a", "b"), esid = 1:3)
  V <- matrix(c(.25, .1, 0, .1, .49, 0, 0, 0, .64), 3L)
  object <- bselmodel.mv(yi = yi, V = V, random = ~ 1 | study / esid, data = dat,
    selection = selection_model(group = study), measure = "GEN", prior_unit_information_sd = 1,
    only_priors = TRUE, silent = TRUE)
  plan <- .data_selection_execution_plan(object$data)
  expect_identical(unname(plan$row_blocks), list(1:2, 3L))
  samples <- matrix(1:3, 3L, 1L, dimnames = list(NULL, "draw_id"))
  active <- c(3L, 1L)
  tau_e <- c(.1, .2, .3)
  tau_b <- c(0, .4, .5)
  sei <- sqrt(diag(V))
  predictive <- list(mu = matrix(c(.1, .2, .3), 3L, 3L), sei = sei,
    tau_within = matrix(sqrt(tau_e^2 + tau_b^2), 3L, 3L))
  predictive$mu_extrapolated <- predictive$mu
  prior <- BayesTools::prior_weightfunction("one-sided", .025, BayesTools::wf_fixed(c(1, .2)))
  selection <- .selection_spec(list(outcome = list(bias = prior)), dat$yi, sei,
    effect_direction = "positive", signed_data = FALSE)
  selection$omega <- matrix(c(1, .2, 1, .4, 1, .6), 3L, byrow = TRUE)
  selection$alpha <- numeric(3L)
  selection$phack_kind <- integer(3L)
  selection$kernel_mode <- c(SELKERNEL_STEP, SELKERNEL_NORMAL, SELKERNEL_STEP)
  selection$vector_rule <- integer(3L)
  selection$use_normal <- c(FALSE, TRUE, FALSE)
  integrated <- list(diagonal = matrix(tau_e[active]^2, 2L, 3L), ranks = c(0L, 0L),
    loadings = list(array(numeric(), c(2L, 2L, 0L)), array(numeric(), c(2L, 1L, 0L))),
    loading_supports = list(matrix(logical(), 2L, 0L), matrix(logical(), 1L, 0L)))
  retained <- list(diagonal = cbind(0, 0, tau_b[active]^2), ranks = c(1L, 0L),
    loadings = list(array(matrix(tau_b[active], 2L, 2L), c(2L, 2L, 1L)), array(numeric(), c(2L, 1L, 0L))))
  prepared <- list(active = active, setup = list(data = object$data, S = 2L, K = 3L),
    execution_plan = plan, integrated = integrated, retained = retained,
    selection = BayesTools::selection_context_subset_rows(selection, active))
  calls <- list()
  fallback <- list()
  scalar <- list()
  z <- c(-1, 0, 2)
  testthat::local_mocked_bindings(
    .zplot_joint_factor_preparation = function(...) prepared,
    .predict_joint_selection_gaussian_parts = function(...) stop("Dense covariance cubes must not be requested."),
    .zplot_context_projection = function(z, mean, covariance, context_factor, sei, selection, control) {
      stopifnot(exists("static", envir = selection$native_cache, inherits = FALSE))
      calls[[length(calls) + 1L]] <<- list(mean = mean, C = covariance, L = context_factor,
        omega = selection$omega, static = BayesTools::selection_native_static_args(selection))
      if (selection$omega[1L, 2L] == .6) return(NULL)
      list(density = matrix(.4, 1L, length(z)), relative_error = 0, mass_error = 0)
    },
    .zplot_full_event_context_draw = function(z, mean, sigma, latent, sei, context, probability,
                                              control, execution_plan, ...) {
      fallback[[length(fallback) + 1L]] <<- list(mean = mean, C = sigma, G = latent,
        omega = context$omega, plan = execution_plan)
      matrix(.8, 1L, length(z))
    },
    .zplot_latent_mixture = function(z, mean, sd, latent_sd, sei, selection, ...) {
      scalar <<- list(mean = mean, sd = sd, retained = latent_sd, omega = selection$omega)
      list(fitted = matrix(c(.15, .25), 2L, length(z)))
    },
    .package = "RoBMA"
  )
  result <- .zplot_joint_marginal(object, samples, predictive, selection, z, FALSE,
    set_selection_likelihood_control())
  expect_length(calls, 2L)
  expect_equal(vapply(calls, function(call) call$omega[1L, 2L], numeric(1L)), c(.6, .2))
  expect_identical(calls[[1L]]$static, calls[[2L]]$static)
  expect_equal(calls[[1L]]$C, V[1:2, 1:2] + diag(tau_e[3L]^2, 2L))
  expect_equal(calls[[2L]]$C, V[1:2, 1:2] + diag(tau_e[1L]^2, 2L))
  expect_length(fallback, 1L)
  expect_equal(fallback[[1L]]$G, matrix(tau_b[3L]^2, 2L, 2L))
  expect_identical(fallback[[1L]]$plan, plan)
  expect_equal(as.numeric(scalar$retained), tau_b[active])
  expect_equal(as.numeric(scalar$sd), sqrt(V[3L, 3L] + tau_e[active]^2))
  expect_equal(result$fitted[3L, ], rep(.8 * 2 / 3 + .15 / 3, length(z)))
  expect_equal(result$fitted[1L, ], rep(.4 * 2 / 3 + .25 / 3, length(z)))
  expect_equal(result$fitted[2L, ], result$extrapolated[2L, ])
})
