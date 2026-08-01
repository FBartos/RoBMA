context("Summary heterogeneity helpers")

test_that("rho allocation retains endpoints and rejects invalid values", {

  tau <- matrix(.3, nrow = 3L, ncol = 2L)
  out <- RoBMA:::.heterogeneity_components(
    tau_total     = tau,
    rho           = c(0, .5, 1),
    is_multilevel = TRUE
  )

  expect_equal(out[["rho"]], c(0, .5, 1))
  expect_equal(out[["tau_within"]][1, ], c(.3, .3))
  expect_equal(out[["tau_between"]][1, ], c(0, 0))
  expect_equal(out[["tau_within"]][3, ], c(0, 0))
  expect_equal(out[["tau_between"]][3, ], c(.3, .3))
  expect_equal(
    out[["tau_within"]]^2 + out[["tau_between"]]^2,
    tau^2
  )

  expect_error(
    RoBMA:::.heterogeneity_components(
      tau[1, , drop = FALSE],
      -1e-12,
      TRUE
    ),
    "within \\[0, 1\\]"
  )
  expect_error(
    RoBMA:::.heterogeneity_components(
      tau[1, , drop = FALSE],
      1 + 1e-12,
      TRUE
    ),
    "within \\[0, 1\\]"
  )

  posterior <- matrix(
    c(.2, .4),
    ncol     = 1L,
    dimnames = list(NULL, "tau")
  )
  fixed <- RoBMA:::.evaluate.brma.tau(
    fit               = NULL,
    scale_data        = NULL,
    scale_formula     = NULL,
    scale_priors      = NULL,
    is_scale          = FALSE,
    is_multilevel     = TRUE,
    K                 = 1L,
    posterior_samples = posterior,
    fixed_rho         = 1
  )
  expect_equal(fixed[["tau_within"]], matrix(0, nrow = 2L))
  expect_equal(as.numeric(fixed[["tau_between"]]), c(.2, .4))
})


test_that("fixed rho summaries use evaluated allocation samples", {

  observed_rho <- NULL
  testthat::local_mocked_bindings(
    .outcome_data_vi = function(object) c(1, 1),
    .get_model_matrix = function(object) matrix(1, nrow = 2L, ncol = 1L),
    .get_posterior_samples = function(fit) matrix(
      c(.2, .4),
      ncol     = 1L,
      dimnames = list(NULL, "tau")
    ),
    .evaluate.brma.tau = function(...) list(
      tau_within  = matrix(0, nrow = 2L, ncol = 2L),
      tau_between = matrix(c(.2, .2, .4, .4), nrow = 2L),
      rho         = c(1, 1)
    ),
    .summary_heterogeneity_samples = function(
        tau_within_samples, tau_between_samples, rho_samples, ...) {
      observed_rho <<- rho_samples
      list(rho = rho_samples)
    },
    .package = "RoBMA"
  )

  data <- list(scale = NULL)
  attr(data, "cluster") <- TRUE
  attr(data, "scale")   <- FALSE
  object <- structure(list(
    fit    = NULL,
    data   = data,
    priors = list(
      outcome = list(
        rho = BayesTools::prior(
          "point",
          parameters = list(location = 1)
        )
      ),
      scale = NULL
    )
  ), class = "brma")

  summary_heterogeneity.brma(object)

  expect_identical(observed_rho, c(1, 1))
})

test_that("scale heterogeneity summaries aggregate variances before tau", {

  tau_within <- matrix(
    c(1, 3,
      2, 4),
    nrow = 2,
    byrow = TRUE
  )

  samples <- RoBMA:::.summary_heterogeneity_samples(
    tau_within_samples  = tau_within,
    tau_between_samples = matrix(0, nrow = 2, ncol = 2),
    v_tilde             = 1,
    is_multilevel       = FALSE
  )

  expected_tau2 <- rowMeans(tau_within^2)
  old_tau2      <- rowMeans(tau_within)^2

  expect_equal(samples[["tau2"]], expected_tau2)
  expect_equal(samples[["tau"]], sqrt(expected_tau2))
  expect_false(isTRUE(all.equal(samples[["tau2"]], old_tau2)))
  expect_equal(samples[["I2"]], rowMeans(100 * tau_within^2 / (tau_within^2 + 1)))
  expect_equal(samples[["H2"]], rowMeans(tau_within^2 + 1))
})

test_that("multilevel scale heterogeneity partitions variance and I2", {

  tau_within <- matrix(
    c(1, 2,
      3, 4),
    nrow = 2,
    byrow = TRUE
  )
  tau_between <- matrix(
    c(2, 1,
      1, 3),
    nrow = 2,
    byrow = TRUE
  )

  samples <- RoBMA:::.summary_heterogeneity_samples(
    tau_within_samples  = tau_within,
    tau_between_samples = tau_between,
    rho_samples         = c(0.8, 0.2),
    v_tilde             = 2,
    is_multilevel       = TRUE
  )

  sigma2_within  <- tau_within^2
  sigma2_between <- tau_between^2
  sigma2_total   <- sigma2_within + sigma2_between
  denominator    <- sigma2_total + 2

  expect_equal(samples[["tau2 [within]"]], rowMeans(sigma2_within))
  expect_equal(samples[["tau2 [between]"]], rowMeans(sigma2_between))
  expect_equal(samples[["tau2"]], rowMeans(sigma2_total))
  expect_equal(samples[["rho"]], c(0.8, 0.2))
  expect_equal(samples[["I2"]], rowMeans(100 * sigma2_total / denominator))
  expect_equal(samples[["I2 [within]"]], rowMeans(100 * sigma2_within / denominator))
  expect_equal(samples[["I2 [between]"]], rowMeans(100 * sigma2_between / denominator))
  expect_equal(samples[["I2"]], samples[["I2 [within]"]] + samples[["I2 [between]"]])
})


test_that("heterogeneity summaries reject invalid rho without projection", {

  tau <- matrix(0, nrow = 2L, ncol = 1L)
  unidentified <- RoBMA:::.summary_heterogeneity_samples(
    tau_within_samples  = tau,
    tau_between_samples = tau,
    v_tilde             = 1,
    is_multilevel       = TRUE
  )

  expect_true(all(is.na(unidentified[["rho"]])))
  expect_error(
    RoBMA:::.summary_heterogeneity_samples(
      tau_within_samples  = tau,
      tau_between_samples = tau,
      rho_samples         = c(-1e-12, 1),
      v_tilde             = 1,
      is_multilevel       = TRUE
    ),
    "within \\[0, 1\\]"
  )
})
