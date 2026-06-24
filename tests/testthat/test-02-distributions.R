context("Distribution helpers")
skip_on_cran()

test_that(".rowLogSumExps handles all -Inf rows", {

  x <- matrix(c(-Inf, -Inf, -1000, -1001), nrow = 2, byrow = TRUE)

  expect_equal(.rowLogSumExps(x)[1], -Inf)
  expect_equal(.rowLogSumExps(x)[2], -1000 + log1p(exp(-1)))
})

test_that(".rowLogSumExps propagates missing quadrature failures", {

  x <- matrix(
    c(0, NaN,
      NA_real_, -1,
      NaN, NaN,
      -Inf, -Inf,
      -1000, -1001),
    nrow  = 5,
    byrow = TRUE
  )
  out <- .rowLogSumExps(x)

  expect_true(is.na(out[1]))
  expect_true(is.na(out[2]))
  expect_true(is.na(out[3]))
  expect_equal(out[4], -Inf)
  expect_equal(out[5], -1000 + log1p(exp(-1)))
})


test_that("analytic normal cluster likelihood matches mvtnorm oracle", {

  S <- 3L
  K <- 6L
  setup <- list(
    S           = S,
    mu          = matrix(c(
      -.10, .00, .05, .20, .30, .35,
       .02, .08, .12, .18, .22, .28,
      -.15, .04, .09, .13, .19, .25
    ), nrow = S, byrow = TRUE),
    tau_within  = matrix(c(
      .08, .10, .12, .14, .16, .18,
      .09, .11, .13, .15, .17, .19,
      .07, .12, .10, .18, .14, .20
    ), nrow = S, byrow = TRUE),
    tau_between = matrix(c(
      .05, .07, .09, .11, .13, .15,
      .06, .08, .10, .12, .14, .16,
      .04, .09, .06, .13, .08, .17
    ), nrow = S, byrow = TRUE),
    cluster     = list(single = 1L, pair = c(2L, 4L), triple = c(3L, 5L, 6L))
  )
  yi <- c(.04, .15, -.05, .30, .08, .18)
  vi <- c(.10, .12, .08, .20, .15, .11)^2

  analytic <- .log_lik_cluster_norm_analytic(
    setup = setup,
    yi    = yi,
    vi    = vi
  )
  oracle <- matrix(NA_real_, nrow = S, ncol = length(setup[["cluster"]]))
  for (s in seq_len(S)) {
    for (g in seq_along(setup[["cluster"]])) {
      idx        <- setup[["cluster"]][[g]]
      diagonal   <- vi[idx] + setup[["tau_within"]][s, idx]^2
      covariance <- diag(diagonal, nrow = length(idx), ncol = length(idx)) +
        tcrossprod(setup[["tau_between"]][s, idx])
      oracle[s, g] <- mvtnorm::dmvnorm(
        x     = yi[idx],
        mean  = setup[["mu"]][s, idx],
        sigma = covariance,
        log   = TRUE
      )
    }
  }

  expect_equal(analytic, oracle, tolerance = 1e-12)
})


test_that("analytic normal cluster likelihood agrees with unweighted quadrature", {

  S <- 2L
  setup <- list(
    S           = S,
    mu          = matrix(c(
      -.04, .05, .12, .18,
       .02, .08, .15, .22
    ), nrow = S, byrow = TRUE),
    tau_within  = matrix(c(
      .06, .09, .12, .15,
      .08, .10, .14, .16
    ), nrow = S, byrow = TRUE),
    tau_between = matrix(c(
      .04, .07, .10, .12,
      .05, .08, .09, .13
    ), nrow = S, byrow = TRUE),
    cluster     = list(a = c(1L, 3L), b = c(2L, 4L)),
    weights     = NULL
  )
  yi  <- c(.02, .15, -.05, .30)
  sei <- c(.10, .12, .08, .20)

  analytic <- .log_lik_cluster_norm_analytic(
    setup = setup,
    yi    = yi,
    vi    = sei^2
  )
  quadrature <- .log_lik_cluster_norm_quadrature_r(
    setup             = setup,
    yi                = yi,
    sei               = sei,
    is_weightfunction = FALSE,
    n_gamma           = 61L
  )

  expect_equal(quadrature, analytic, tolerance = 1e-8)
})


test_that("GLMM prior grid keys reuse cached prior hashes", {

  clear_cache <- function() {
    rm(
      list  = ls(envir = .glmm_prior_grid_hash_cache),
      envir = .glmm_prior_grid_hash_cache
    )
  }
  clear_cache()
  on.exit(clear_cache(), add = TRUE)

  prior <- BayesTools::prior("normal", list(0, 1))
  other <- BayesTools::prior("normal", list(0, 2))
  calls <- 0L

  testthat::local_mocked_bindings(
    .glmm_prior_grid_hash_compute = function(prior) {

      calls <<- calls + 1L
      c("11111111", "22222222")
    },
    .package = "RoBMA"
  )

  expect_equal(.glmm_prior_grid_key(prior, 9), "9|11111111|22222222|")
  expect_equal(.glmm_prior_grid_key(prior, 11, "bin|4"),
               "11|11111111|22222222|bin|4")
  expect_equal(calls, 1L)

  .glmm_prior_grid_key(other, 9)
  expect_equal(calls, 2L)
})


test_that("native coercion helpers preserve matrix and vector shapes", {

  x <- matrix(1:4, nrow = 2)

  expect_true(is.double(.native_numeric_matrix(x)))
  expect_equal(dim(.native_numeric_matrix(1:4)), c(4L, 1L))
  expect_true(is.double(.native_numeric_vector(1:4)))
  expect_true(is.integer(.native_integer_vector(c(1, 2))))
})

test_that("native normal cluster likelihood matches R quadrature", {

  skip_if_not(.has_native_norm_cluster_quadrature(selection = FALSE))

  set.seed(20260514)
  S <- 4
  K <- 5
  setup <- list(
    S           = S,
    mu          = matrix(rnorm(S * K, 0, .25), nrow = S, ncol = K),
    tau_within  = matrix(runif(S * K, .03, .25), nrow = S, ncol = K),
    tau_between = matrix(runif(S * K, .02, .20), nrow = S, ncol = K),
    cluster     = list(a = c(1L, 3L), b = c(2L, 4L, 5L)),
    weights     = c(1, .5, 1.25, 2, .75)
  )
  yi  <- c(.02, .15, -.05, .30, .08)
  sei <- c(.10, .12, .08, .20, .15)

  native <- .log_lik_cluster_norm_quadrature_native(
    setup             = setup,
    yi                = yi,
    sei               = sei,
    is_weightfunction = FALSE,
    n_gamma           = 5
  )
  ref <- .log_lik_cluster_norm_quadrature_r(
    setup             = setup,
    yi                = yi,
    sei               = sei,
    is_weightfunction = FALSE,
    n_gamma           = 5
  )

  expect_equal(native, ref, tolerance = 1e-12)
  expect_equal(
    .log_lik_cluster_norm_quadrature_sum(
      setup             = setup,
      yi                = yi,
      sei               = sei,
      is_weightfunction = FALSE,
      n_gamma           = 5
    ),
    rowSums(ref),
    tolerance = 1e-12
  )
})

test_that("native selected-normal cluster likelihood matches high-order R quadrature", {

  skip_if_not(.has_native_norm_cluster_quadrature(selection = TRUE))

  prior <- BayesTools::prior_weightfunction(
    side    = "one-sided",
    steps   = c(.025, .05),
    weights = BayesTools::wf_fixed(c(1, .5, .25))
  )
  yi  <- c(.02, .15, -.05, .30, .08)
  sei <- c(.10, .12, .08, .20, .15)
  spec <- .selection_spec(
    priors           = list(outcome = list(bias = prior)),
    yi               = yi,
    sei              = sei,
    effect_direction = "positive"
  )

  S <- 4
  K <- length(yi)
  setup <- list(
    S           = S,
    mu          = matrix(c(
      -.05, .00, .05, .10, .15,
       .02, .08, .12, .18, .22,
      -.10, .04, .09, .13, .19,
       .05, .10, .15, .20, .25
    ), nrow = S, byrow = TRUE),
    tau_within  = matrix(runif(S * K, .04, .22), nrow = S, ncol = K),
    tau_between = matrix(runif(S * K, .02, .18), nrow = S, ncol = K),
    cluster     = list(a = c(1L, 3L), b = c(2L, 4L, 5L)),
    weights     = c(1, .5, 1.25, 2, .75)
  )
  selection_context <- spec
  selection_context[["omega"]] <- matrix(c(
    1, .5, .25,
    1, .7, .40,
    1, .4, .20,
    1, .8, .35
  ), nrow = S, byrow = TRUE)
  selection_context[["alpha"]]       <- rep(0, S)
  selection_context[["phack_kind"]]  <- rep(0L, S)
  selection_context[["kernel_mode"]] <- c(
    SELKERNEL_STEP,
    SELKERNEL_NORMAL,
    SELKERNEL_STEP,
    SELKERNEL_STEP
  )

  native <- .log_lik_cluster_norm_quadrature_native(
    setup             = setup,
    yi                = yi,
    sei               = sei,
    is_weightfunction = TRUE,
    selection_context = selection_context,
    n_gamma           = 5
  )
  ref <- .log_lik_cluster_norm_quadrature_r(
    setup             = setup,
    yi                = yi,
    sei               = sei,
    is_weightfunction = TRUE,
    selection_context = selection_context,
    n_gamma           = 61
  )

  expect_equal(native, ref, tolerance = 1e-3)
  expect_equal(
    .log_lik_cluster_norm_quadrature_sum(
      setup             = setup,
      yi                = yi,
      sei               = sei,
      is_weightfunction = TRUE,
      selection_context = selection_context,
      n_gamma           = 5
    ),
    rowSums(ref),
    tolerance = 1e-3
  )
})

test_that("native estimate row-sum kernels match matrix likelihoods", {

  skip_if_not(.has_native_norm_loglik_row_sum(selection = FALSE))
  skip_if_not(.has_native_norm_loglik_row_sum(selection = TRUE))

  yi  <- c(.02, .15, -.05, .30)
  sei <- c(.10, .12, .08, .20)
  S   <- 3
  K   <- length(yi)

  mu <- matrix(c(
    -.05, .00, .05, .10,
     .02, .08, .12, .18,
    -.10, .04, .09, .13
  ), nrow = S, byrow = TRUE)
  tau <- matrix(c(
    .04, .08, .12, .16,
    .05, .09, .13, .17,
    .06, .10, .14, .18
  ), nrow = S, byrow = TRUE)
  weights <- c(1, .5, 1.25, 2)

  expect_equal(
    .outcome_pdf_sum.norm(
      yi         = yi,
      mu_samples = mu,
      tau_within = tau,
      sei        = sei,
      weights    = weights
    ),
    rowSums(.apply_log_lik_weights(
      .outcome_pdf.norm(
        yi         = yi,
        mu_samples = mu,
        tau_within = tau,
        sei        = sei
      ),
      weights
    )),
    tolerance = 1e-12
  )

  prior <- BayesTools::prior_weightfunction(
    side    = "one-sided",
    steps   = c(.025, .05),
    weights = BayesTools::wf_fixed(c(1, .5, .25))
  )
  selection_context <- .selection_spec(
    priors           = list(outcome = list(bias = prior)),
    yi               = yi,
    sei              = sei,
    effect_direction = "positive"
  )
  selection_context[["omega"]] <- matrix(c(
    1, .5, .25,
    1, .7, .40,
    1, .4, .20
  ), nrow = S, byrow = TRUE)
  selection_context[["alpha"]]       <- rep(0, S)
  selection_context[["phack_kind"]]  <- rep(0L, S)
  selection_context[["kernel_mode"]] <- c(
    SELKERNEL_STEP,
    SELKERNEL_NORMAL,
    SELKERNEL_STEP
  )
  selection_context <- .selection_reset_native_cache(selection_context)
  total_sd          <- sqrt(tau^2 + matrix(sei, nrow = S, ncol = K, byrow = TRUE)^2)

  expect_equal(
    .outcome_pdf_sum.selnorm(
      yi                = yi,
      mu_samples        = mu,
      tau_within        = tau,
      sei               = sei,
      selection_context = selection_context,
      weights           = weights
    ),
    rowSums(.selnorm_kernel_loglik_matrix(
      yi             = yi,
      mu_num         = mu,
      sigma_num      = total_sd,
      mu_norm        = mu,
      sigma_norm     = total_sd,
      sei            = sei,
      omega          = selection_context[["omega"]],
      selection_spec = selection_context,
      alpha          = selection_context[["alpha"]],
      phack_kind     = selection_context[["phack_kind"]],
      kernel_mode    = selection_context[["kernel_mode"]],
      weights        = weights
    )),
    tolerance = 1e-12
  )
})

test_that("selected-normal helpers separate likelihood and selection SEs", {

  yi             <- c(2.0, 1.7)
  likelihood_sei <- c(.20, .25)
  selection_sei  <- c(1.00, .80)
  S              <- 2
  K              <- length(yi)
  mu             <- matrix(c(0, .1, .2, .3), nrow = S, byrow = TRUE)
  tau            <- matrix(c(.05, .06, .07, .08), nrow = S, byrow = TRUE)
  prior          <- BayesTools::prior_weightfunction(
    side    = "one-sided",
    steps   = c(.025, .05),
    weights = BayesTools::wf_fixed(c(1, .2, .05))
  )
  selection_context <- .selection_spec(
    priors           = list(outcome = list(bias = prior)),
    yi               = yi,
    sei              = selection_sei,
    effect_direction = "positive"
  )
  selection_context[["omega"]] <- matrix(c(
    1, .2, .05,
    1, .4, .10
  ), nrow = S, byrow = TRUE)
  selection_context[["alpha"]]       <- rep(0, S)
  selection_context[["phack_kind"]]  <- rep(0L, S)
  selection_context[["kernel_mode"]] <- rep(SELKERNEL_STEP, S)
  selection_context <- .selection_reset_native_cache(selection_context)

  total_sd <- sqrt(
    tau^2 + matrix(likelihood_sei, nrow = S, ncol = K, byrow = TRUE)^2
  )
  separated <- .outcome_pdf.selnorm(
    yi                = yi,
    mu_samples        = mu,
    tau_within        = tau,
    sei               = likelihood_sei,
    selection_sei     = selection_sei,
    selection_context = selection_context
  )
  direct <- .selnorm_kernel_loglik_matrix(
    yi             = yi,
    mu_num         = mu,
    sigma_num      = total_sd,
    mu_norm        = mu,
    sigma_norm     = total_sd,
    sei            = selection_sei,
    omega          = selection_context[["omega"]],
    selection_spec = selection_context,
    alpha          = selection_context[["alpha"]],
    phack_kind     = selection_context[["phack_kind"]],
    kernel_mode    = selection_context[["kernel_mode"]],
    weights        = NULL
  )
  old_scale <- .outcome_pdf.selnorm(
    yi                = yi,
    mu_samples        = mu,
    tau_within        = tau,
    sei               = likelihood_sei,
    selection_context = selection_context
  )

  expect_equal(separated, direct, tolerance = 1e-12)
  expect_false(isTRUE(all.equal(separated, old_scale)))
  expect_equal(
    .outcome_pdf_sum.selnorm(
      yi                = yi,
      mu_samples        = mu,
      tau_within        = tau,
      sei               = likelihood_sei,
      selection_sei     = selection_sei,
      selection_context = selection_context
    ),
    rowSums(separated),
    tolerance = 1e-12
  )
})

test_that("estimate selected-normal likelihood forwards selection SEs", {

  captured <- list()
  setup <- list(
    fit                = list(),
    data               = list(outcome = data.frame(yi = c(.1, .2), sei = c(1, 1.5))),
    priors             = list(),
    yi                 = c(.1, .2),
    sei                = c(.2, .3),
    selection_sei      = c(1, 1.5),
    K                  = 2L,
    S                  = 1L,
    mu                 = matrix(0, nrow = 1, ncol = 2),
    tau_within         = matrix(0, nrow = 1, ncol = 2),
    weights            = NULL,
    is_weightfunction  = TRUE,
    outcome_type       = "norm",
    effect_direction   = "positive",
    posterior_samples  = matrix(numeric(0), nrow = 1, ncol = 0)
  )

  testthat::local_mocked_bindings(
    .selection_context_from_parts = function(...) {

      list(marker = TRUE)
    },
    .outcome_pdf.selnorm = function(yi, mu_samples, tau_within, sei,
                                    selection_context, weights = NULL,
                                    selection_sei = sei) {

      captured[["matrix"]] <<- list(sei = sei, selection_sei = selection_sei)
      matrix(0, nrow = nrow(mu_samples), ncol = ncol(mu_samples))
    },
    .outcome_pdf_sum.selnorm = function(yi, mu_samples, tau_within, sei,
                                        selection_context, weights = NULL,
                                        selection_sei = sei) {

      captured[["sum"]] <<- list(sei = sei, selection_sei = selection_sei)
      rep(0, nrow(mu_samples))
    },
    .package = "RoBMA"
  )

  expect_equal(.log_lik_estimate_from_setup(setup), matrix(0, nrow = 1, ncol = 2))
  expect_equal(.log_lik_estimate_sum_from_setup(setup), 0)
  expect_equal(captured[["matrix"]][["sei"]], setup[["sei"]])
  expect_equal(captured[["matrix"]][["selection_sei"]], setup[["selection_sei"]])
  expect_equal(captured[["sum"]][["sei"]], setup[["sei"]])
  expect_equal(captured[["sum"]][["selection_sei"]], setup[["selection_sei"]])
})

test_that("marglik selected-normal known-V path separates SE roles", {

  captured <- NULL
  data <- list(
    K             = 2L,
    yi            = c(.1, .2),
    sei           = c(1, 1.5),
    sampling_var  = c(.04, .09),
    sampling_rank = 0L
  )
  parameters <- list(mu = 0, tau = 0)

  testthat::local_mocked_bindings(
    .marglik_selection_context = function(parameters, data) {

      list(marker = TRUE)
    },
    .outcome_pdf.selnorm = function(yi, mu_samples, tau_within, sei,
                                    selection_context, weights = NULL,
                                    selection_sei = sei) {

      captured <<- list(sei = sei, selection_sei = selection_sei)
      matrix(0, nrow = nrow(mu_samples), ncol = ncol(mu_samples))
    },
    .package = "RoBMA"
  )

  expect_equal(.log_posterior(
    parameters        = parameters,
    data              = data,
    is_mods           = FALSE,
    is_scale          = FALSE,
    is_multilevel     = FALSE,
    is_weights        = FALSE,
    is_known_v        = TRUE,
    is_PET            = FALSE,
    is_PEESE          = FALSE,
    is_weightfunction = TRUE,
    effect_direction  = "positive",
    outcome_type      = "norm"
  ), 0)
  expect_equal(captured[["sei"]], sqrt(data[["sampling_var"]]))
  expect_equal(captured[["selection_sei"]], data[["sei"]])
})

test_that("native GLMM row-sum kernels match matrix likelihoods", {

  skip_if_not(.has_native_glmm_row_sum("bin"))
  skip_if_not(.has_native_glmm_row_sum("pois"))
  skip_if_not(.has_native_glmm_row_sum("bin", conditional = TRUE))
  skip_if_not(.has_native_glmm_row_sum("pois", conditional = TRUE))

  ai  <- c(3L, 0L, 8L, 12L)
  ci  <- c(2L, 4L, 0L, 10L)
  n1i <- c(20L, 18L, 15L, 30L)
  n2i <- c(22L, 17L, 16L, 31L)
  x1i <- c(0L, 3L, 10L, 12L)
  x2i <- c(1L, 0L, 8L, 9L)
  t1i <- c(12, 30, 45, 50)
  t2i <- c(10, 28, 43, 48)

  mu <- matrix(
    c(-0.4, 0.1, 0.7, 0.2,
      -0.2, 0.3, 0.5, 0.0,
       0.1, 0.5, 0.3, -0.2),
    nrow = 3,
    ncol = 4
  )
  tau <- matrix(
    c(0.0, 0.1, 0.4, 0.2,
      0.2, 0.3, 0.1, 0.5,
      0.4, 0.2, 0.0, 0.3),
    nrow = 3,
    ncol = 4
  )
  weights   <- c(1, .5, 1.25, 2)
  prior_pi  <- BayesTools::prior("beta", list(1.5, 2.5))
  prior_phi <- BayesTools::prior("normal", list(-1, 1.5))

  expect_equal(
    .outcome_pdf_sum.binom(
      ai, ci, n1i, n2i, mu, tau, prior_pi,
      weights = weights, n_theta = 5, n_pi = 7
    ),
    rowSums(.outcome_pdf.binom(
      ai, ci, n1i, n2i, mu, tau, prior_pi,
      weights = weights, n_theta = 5, n_pi = 7
    )),
    tolerance = 1e-10
  )
  expect_equal(
    .outcome_pdf_sum.pois(
      x1i, x2i, t1i, t2i, mu, tau, prior_phi,
      weights = weights, n_theta = 5, n_phi = 7
    ),
    rowSums(.outcome_pdf.pois(
      x1i, x2i, t1i, t2i, mu, tau, prior_phi,
      weights = weights, n_theta = 5, n_phi = 7
    )),
    tolerance = 1e-10
  )

  logit_baserate <- matrix(
    qlogis(c(.20, .30, .40, .25,
             .25, .35, .45, .30,
             .30, .40, .50, .35)),
    nrow = 3,
    byrow = TRUE
  )
  log_phi <- matrix(
    log(c(.10, .15, .20, .25,
          .12, .18, .22, .28,
          .14, .20, .24, .30)),
    nrow = 3,
    byrow = TRUE
  )

  expect_equal(
    .outcome_pdf_sum.binom_conditional(
      ai, ci, n1i, n2i, mu, logit_baserate, weights = weights
    ),
    rowSums(.apply_log_lik_weights(
      .outcome_pdf.binom_conditional(
        ai, ci, n1i, n2i, mu, logit_baserate
      ),
      weights
    )),
    tolerance = 1e-12
  )
  expect_equal(
    .outcome_pdf_sum.pois_conditional(
      x1i, x2i, t1i, t2i, mu, log_phi, weights = weights
    ),
    rowSums(.apply_log_lik_weights(
      .outcome_pdf.pois_conditional(
        x1i, x2i, t1i, t2i, mu, log_phi
      ),
      weights
    )),
    tolerance = 1e-12
  )
})

test_that("native GLMM cluster row-sum kernels match matrix kernels", {

  skip_if_not(.has_native_glmm_row_sum("bin", cluster = TRUE))
  skip_if_not(.has_native_glmm_row_sum("pois", cluster = TRUE))

  set.seed(20260514)
  S <- 4
  K <- 5
  setup <- list(
    mu          = matrix(rnorm(S * K, 0, .25), nrow = S, ncol = K),
    tau_within  = matrix(runif(S * K, .05, .25), nrow = S, ncol = K),
    tau_between = matrix(runif(S * K, .02, .18), nrow = S, ncol = K),
    cluster     = list(a = c(1L, 3L), b = c(2L, 4L, 5L)),
    weights     = c(1, .5, 1.25, 2, .75)
  )
  bin_data <- list(outcome = data.frame(
    ai  = c(3L, 0L, 8L, 12L, 2L),
    ci  = c(2L, 4L, 0L, 10L, 1L),
    n1i = c(20L, 18L, 15L, 30L, 22L),
    n2i = c(22L, 17L, 16L, 31L, 21L)
  ))
  pois_data <- list(outcome = data.frame(
    x1i = c(0L, 3L, 10L, 12L, 1L),
    x2i = c(1L, 0L, 8L, 9L, 2L),
    t1i = c(12, 30, 45, 50, 25),
    t2i = c(10, 28, 43, 48, 24)
  ))
  bin_priors <- list(outcome = list(
    pi = BayesTools::prior("beta", list(1.5, 2.5))
  ))
  pois_priors <- list(outcome = list(
    phi = BayesTools::prior("normal", list(-1, 1.5))
  ))

  expect_equal(
    .log_lik_cluster_glmm_native_sum(
      setup, bin_data, bin_priors, "bin",
      n_theta = 5, n_gamma = 5, n_pi = 7
    ),
    rowSums(.log_lik_cluster_glmm_native(
      setup, bin_data, bin_priors, "bin",
      n_theta = 5, n_gamma = 5, n_pi = 7
    )),
    tolerance = 1e-10
  )
  expect_equal(
    .log_lik_cluster_glmm_native_sum(
      setup, pois_data, pois_priors, "pois",
      n_theta = 5, n_gamma = 5, n_phi = 7
    ),
    rowSums(.log_lik_cluster_glmm_native(
      setup, pois_data, pois_priors, "pois",
      n_theta = 5, n_gamma = 5, n_phi = 7
    )),
    tolerance = 1e-10
  )
})

test_that("native log-sum-exp kernels propagate bad quadrature nodes", {

  skip_if_not(.has_native_glmm("bin"))
  skip_if_not(.has_native_norm_cluster_quadrature(selection = FALSE))

  binom_out <- .Call(
    "RoBMA_glmm_binom_marginal_loglik",
    .native_integer_vector(3L),
    .native_integer_vector(2L),
    .native_integer_vector(20L),
    .native_integer_vector(22L),
    matrix(0, nrow = 1, ncol = 1),
    matrix(.1, nrow = 1, ncol = 1),
    NULL,
    c(-1, 0, 1),
    c(log(.3), NaN, log(.7)),
    matrix(0, nrow = 1, ncol = 1),
    matrix(0, nrow = 1, ncol = 1),
    PACKAGE = "RoBMA"
  )
  expect_true(is.na(binom_out[1, 1]))

  norm_out <- .Call(
    "RoBMA_norm_cluster_loglik",
    .native_numeric_vector(.1),
    .native_numeric_vector(.2),
    matrix(0, nrow = 1, ncol = 1),
    matrix(.1, nrow = 1, ncol = 1),
    matrix(.1, nrow = 1, ncol = 1),
    .native_integer_vector(1L),
    .native_integer_vector(1L),
    NULL,
    c(-1, 0, 1),
    c(log(.3), NaN, log(.7)),
    PACKAGE = "RoBMA"
  )
  expect_true(is.na(norm_out[1, 1]))
})

test_that("selection native static cache is reused and reset for subsets", {

  prior <- BayesTools::prior_weightfunction(
    side    = "one-sided",
    steps   = c(.025, .05),
    weights = BayesTools::wf_fixed(c(1, .5, .25))
  )
  selection_context <- .selection_spec(
    priors           = list(outcome = list(bias = prior)),
    yi               = c(.02, .15, -.05),
    sei              = c(.10, .12, .08),
    effect_direction = "positive"
  )
  selection_context[["omega"]]       <- matrix(1, nrow = 3, ncol = 3)
  selection_context[["alpha"]]       <- rep(0, 3)
  selection_context[["phack_kind"]]  <- rep(0L, 3)
  selection_context[["kernel_mode"]] <- rep(SELKERNEL_STEP, 3)
  selection_context <- .selection_reset_native_cache(selection_context)

  first  <- BayesTools::selection_native_static_args(selection_context)
  second <- BayesTools::selection_native_static_args(selection_context)
  expect_identical(first, second)

  row_subset <- BayesTools::selection_context_subset_rows(selection_context, 1:2)
  obs_subset <- BayesTools::selection_context_subset_observations(selection_context, 1:2)

  expect_true(is.environment(row_subset[["native_cache"]]))
  expect_true(is.environment(obs_subset[["native_cache"]]))
  expect_false(identical(
    selection_context[["native_cache"]],
    row_subset[["native_cache"]]
  ))
  expect_false(identical(
    selection_context[["native_cache"]],
    obs_subset[["native_cache"]]
  ))
  expect_equal(
    BayesTools::selection_native_static_args(obs_subset)[["z_lower"]],
    .native_numeric_vector(obs_subset[["z_lower"]])
  )
})

test_that("native selected-normal log-normalizer delta matches matrix reference", {

  skip_if_not(.has_native_selnorm_log_norm_delta())

  prior <- BayesTools::prior_weightfunction(
    side    = "one-sided",
    steps   = c(.025, .05),
    weights = BayesTools::wf_fixed(c(1, .5, .25))
  )
  yi  <- c(.02, .15, -.05)
  sei <- c(.10, .12, .08)
  spec <- .selection_spec(
    priors           = list(outcome = list(bias = prior)),
    yi               = yi,
    sei              = sei,
    effect_direction = "positive"
  )

  mean <- matrix(c(
    -.05, .00, .05,
     .02, .08, .12
  ), nrow = 2, byrow = TRUE)
  sd <- matrix(c(
    .10, .12, .14,
    .11, .13, .16
  ), nrow = 2, byrow = TRUE)
  basis <- matrix(c(
    1, .5, .25,
    1, .2, .10
  ), nrow = 2, byrow = TRUE)
  omega <- matrix(c(
    1, .5, .25,
    1, .7, .40
  ), nrow = 2, byrow = TRUE)
  current     <- c(.01, .03)
  values      <- c(-.05, .00, .10)
  weights     <- c(1, .5, 2)
  kernel_mode <- c(SELKERNEL_STEP, SELKERNEL_NORMAL)

  current_log_norm <- .selnorm_kernel_log_norm_matrix(
    mean           = mean,
    sd             = sd,
    sei            = sei,
    omega          = omega,
    selection_spec = spec,
    kernel_mode    = kernel_mode
  )
  native <- .selnorm_kernel_log_norm_delta_grid(
    mean             = mean,
    sd               = sd,
    basis            = basis,
    current_log_norm = current_log_norm,
    current          = current,
    values           = values,
    sei              = sei,
    weights          = weights,
    omega            = omega,
    selection_spec   = spec,
    kernel_mode      = kernel_mode
  )

  G          <- length(values)
  S          <- nrow(mean)
  expect_true(is.matrix(native))
  expect_equal(dim(native), c(G, S))

  row_index  <- rep(seq_len(S), each = G)
  grid_index <- rep(seq_len(G), times = S)
  delta      <- values[grid_index] - current[row_index]
  candidate_context <- spec
  candidate_context[["omega"]]       <- omega[row_index, , drop = FALSE]
  candidate_context[["alpha"]]       <- rep(0, length(row_index))
  candidate_context[["phack_kind"]]  <- rep(0L, length(row_index))
  candidate_context[["kernel_mode"]] <- kernel_mode[row_index]
  candidate_log_norm <- .selnorm_kernel_log_norm_matrix(
    mean           = mean[row_index, , drop = FALSE] +
      basis[row_index, , drop = FALSE] * delta,
    sd             = sd[row_index, , drop = FALSE],
    sei            = sei,
    omega          = candidate_context[["omega"]],
    selection_spec = candidate_context,
    kernel_mode    = candidate_context[["kernel_mode"]]
  )
  ref <- rowSums(sweep(
    candidate_log_norm - current_log_norm[row_index, , drop = FALSE],
    2L,
    weights,
    "*"
  ))

  expect_equal(as.numeric(native), ref, tolerance = 1e-12)
})

test_that("selected-normal step kernel handles boundary and extreme p-bins", {

  skip_if_not(.has_native_selnorm_kernel())

  prior <- BayesTools::prior_weightfunction(
    side    = "one-sided",
    steps   = c(.001, .025, .05, .50),
    weights = BayesTools::wf_fixed(c(1, .8, .4, .2, .05))
  )
  sei <- c(.05, .10, .20, .40)
  yi  <- stats::qnorm(1 - c(.001, .025, .05, .50)) * sei
  spec <- .selection_spec(
    priors           = list(outcome = list(bias = prior)),
    yi               = yi,
    sei              = sei,
    effect_direction = "positive"
  )

  mu <- matrix(c(
    -.05, .00, .05, .10,
     .10, .05, .00, -.05
  ), nrow = 2, byrow = TRUE)
  sigma <- matrix(c(
    .06, .12, .25, .45,
    .08, .15, .30, .50
  ), nrow = 2, byrow = TRUE)
  omega <- matrix(c(
    1, .8, .4, .2, .05,
    1, .1, .9, .3, .01
  ), nrow = 2, byrow = TRUE)

  out <- .selnorm_kernel_loglik_matrix(
    yi             = yi,
    mu_num         = mu,
    sigma_num      = sigma,
    sei            = sei,
    omega          = omega,
    selection_spec = spec
  )

  expect_true(all(is.finite(out)))
  expect_equal(dim(out), c(2L, 4L))
  expect_equal(out, .test_step_reference(yi, mu, sigma, sei, omega, spec),
               tolerance = 1e-12)
})

test_that("selected-normal log-likelihood matrix honors observation weights", {

  skip_if_not(.has_native_selnorm_kernel())

  prior <- BayesTools::prior_weightfunction(
    side    = "one-sided",
    steps   = c(.025, .05),
    weights = BayesTools::wf_fixed(c(1, .5, .25))
  )
  yi <- c(.05, .15, .30)
  sei <- c(.10, .12, .14)
  spec <- .selection_spec(
    priors           = list(outcome = list(bias = prior)),
    yi               = yi,
    sei              = sei,
    effect_direction = "positive"
  )
  mu      <- matrix(c(.00, .10, .20, .03, .12, .24), nrow = 2, byrow = TRUE)
  sigma   <- matrix(.20, nrow = 2, ncol = 3)
  omega   <- matrix(c(1, .5, .25, 1, .7, .4), nrow = 2, byrow = TRUE)
  weights <- c(1, 2, .25)

  out_unweighted <- .selnorm_kernel_loglik_matrix(
    yi, mu, sigma, sei = sei, omega = omega, selection_spec = spec
  )
  out_weighted <- .selnorm_kernel_loglik_matrix(
    yi, mu, sigma, sei = sei, omega = omega, selection_spec = spec,
    weights = weights
  )

  expect_equal(out_weighted, sweep(out_unweighted, 2, weights, "*"),
               tolerance = 1e-12)
})

test_that("native selected-normal CDF and moments honor normal subdispatch", {

  skip_if_not(.has_native_selnorm_kernel())

  prior <- BayesTools::prior_weightfunction(
    side    = "one-sided",
    steps   = c(.001, .025, .05, .50),
    weights = BayesTools::wf_fixed(c(1, .8, .4, .2, .05))
  )
  yi  <- c(.05, .15, .30, .60)
  sei <- c(.05, .10, .20, .40)
  spec <- .selection_spec(
    priors           = list(outcome = list(bias = prior)),
    yi               = yi,
    sei              = sei,
    effect_direction = "positive"
  )

  mu <- matrix(c(
    -.05, .00, .05, .10,
     .10, .05, .00, -.05,
     .20, .10, .05, .00
  ), nrow = 3, byrow = TRUE)
  sigma <- matrix(c(
    .06, .12, .25, .45,
    .08, .15, .30, .50,
    .10, .20, .35, .60
  ), nrow = 3, byrow = TRUE)
  omega <- matrix(c(
    2, 3, 4, 5, 6,
    1, .1, .9, .3, .01,
    1, .8, .4, .2, .05
  ), nrow = 3, byrow = TRUE)
  kernel_mode <- c(SELKERNEL_NORMAL, SELKERNEL_STEP, SELKERNEL_STEP)
  context <- spec
  context[["omega"]]       <- omega
  context[["alpha"]]       <- rep(0, 3)
  context[["phack_kind"]]  <- rep(0L, 3)
  context[["kernel_mode"]] <- kernel_mode
  context[["use_normal"]]  <- kernel_mode == SELKERNEL_NORMAL

  q <- c(-.10, .05, .20, .80)
  cdf_lower <- .selection_step_cdf_matrix(
    q                 = q,
    mean              = mu,
    sd                = sigma,
    sei               = sei,
    selection_context = context,
    lower.tail        = TRUE
  )
  cdf_upper <- .selection_step_cdf_matrix(
    q                 = q,
    mean              = mu,
    sd                = sigma,
    sei               = sei,
    selection_context = context,
    lower.tail        = FALSE
  )
  moments <- .selection_step_moments_matrix(
    mean              = mu,
    sd                = sigma,
    sei               = sei,
    selection_context = context
  )

  expect_equal(
    cdf_lower,
    .test_step_cdf_reference(q, mu, sigma, sei, omega, spec, kernel_mode),
    tolerance = 1e-12
  )
  expect_equal(cdf_lower + cdf_upper, matrix(1, nrow = 3, ncol = 4),
               tolerance = 1e-10)
  expect_equal(
    moments,
    .test_step_moments_reference(mu, sigma, sei, omega, spec, kernel_mode),
    tolerance = 1e-12
  )
  expect_equal(moments[["mean"]][1, ], mu[1, ], tolerance = 1e-12)
  expect_equal(moments[["second"]][1, ], sigma[1, ]^2 + mu[1, ]^2,
               tolerance = 1e-12)
  expect_equal(
    .selection_step_log_norm_matrix(mu, sigma, sei, context)[1, ],
    rep(0, ncol(mu)),
    tolerance = 1e-12
  )
})

test_that("native selected-normal CDF mirrors negative direction and preserves extreme upper tails", {

  skip_if_not(.has_native_selnorm_kernel())

  prior <- BayesTools::prior_weightfunction(
    side    = "one-sided",
    steps   = c(.025, .05),
    weights = BayesTools::wf_fixed(c(1, .5, .25))
  )
  yi    <- c(.05, .20)
  sei   <- c(.10, .20)
  q     <- c(.10, .50)
  mu    <- matrix(c(.02, .25, .05, .30), nrow = 2, byrow = TRUE)
  sigma <- matrix(c(.15, .30, .18, .35), nrow = 2, byrow = TRUE)
  omega <- matrix(c(1, .5, .25, 1, .7, .4), nrow = 2, byrow = TRUE)

  spec_pos <- .selection_spec(
    priors           = list(outcome = list(bias = prior)),
    yi               = yi,
    sei              = sei,
    effect_direction = "positive"
  )
  context_pos <- spec_pos
  context_pos[["omega"]]       <- omega
  context_pos[["alpha"]]       <- rep(0, 2)
  context_pos[["phack_kind"]]  <- rep(0L, 2)
  context_pos[["kernel_mode"]] <- rep(SELKERNEL_STEP, 2)

  spec_neg <- .selection_spec(
    priors           = list(outcome = list(bias = prior)),
    yi               = -yi,
    sei              = sei,
    effect_direction = "negative"
  )
  context_neg <- spec_neg
  context_neg[["omega"]]       <- omega
  context_neg[["alpha"]]       <- rep(0, 2)
  context_neg[["phack_kind"]]  <- rep(0L, 2)
  context_neg[["kernel_mode"]] <- rep(SELKERNEL_STEP, 2)

  pos_upper <- .selection_step_cdf_matrix(
    q, mu, sigma, sei, context_pos, lower.tail = FALSE
  )
  neg_lower <- .selection_step_cdf_matrix(
    -q, -mu, sigma, sei, context_neg, lower.tail = TRUE
  )
  expect_equal(neg_lower, pos_upper, tolerance = 1e-12)

  context_extreme <- spec_pos
  context_extreme[["omega"]]       <- matrix(1, nrow = 1, ncol = spec_pos[["n_bins"]])
  context_extreme[["alpha"]]       <- 0
  context_extreme[["phack_kind"]]  <- 0L
  context_extreme[["kernel_mode"]] <- SELKERNEL_STEP

  extreme <- .selection_step_cdf_matrix(
    q                 = 10,
    mean              = matrix(0, nrow = 1, ncol = 1),
    sd                = matrix(1, nrow = 1, ncol = 1),
    sei               = 1,
    selection_context = context_extreme,
    lower.tail        = FALSE
  )
  expect_gt(extreme[1, 1], 0)
  expect_equal(extreme[1, 1], stats::pnorm(10, lower.tail = FALSE),
               tolerance = 1e-20)
})
