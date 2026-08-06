.low_rank_known_v_test_matrix <- function() {

  factor <- matrix(
    c(
      -0.62645381074233242, 0.18364332422208224,
      -0.83562861241004716, 1.5952808021377916,
      0.32950777181536051, -0.82046838411801526,
      0.48742905242848528, 0.73832470512921733
    ),
    nrow = 4L,
    ncol = 2L
  )
  tcrossprod(factor)
}


test_that("low-rank known V uses a non-negative spectral representation", {

  V             <- .low_rank_known_v_test_matrix()
  factorization <- .covariance_factorization(V)

  expect_identical(factorization[["covariance"]], V)
  expect_identical(factorization[["status"]], "positive_semidefinite")
  expect_equal(sum(factorization[["spectral_values"]] == 0), 2L)
  expect_true(all(factorization[["spectral_values"]] >= 0))

  sampling_factor <- .covariance_sampling_factor(factorization)
  expect_equal(crossprod(sampling_factor), V, tolerance = 1e-14)

  mu_samples         <- matrix(c(0.10, 0.20, -0.05, 0.15), nrow = 1L)
  covariance_samples <- array(V, dim = c(1L, 4L, 4L))
  set.seed(412)
  response <- .outcome_rng.norm_known_v_covariance(
    mu_samples         = mu_samples,
    covariance_samples = covariance_samples
  )
  expect_identical(dim(response), c(1L, 4L))
  expect_true(all(is.finite(response)))

  indefinite <- matrix(-0.75, nrow = 3L, ncol = 3L)
  diag(indefinite) <- 1
  invalid <- .covariance_factorization(indefinite)
  expect_identical(invalid[["status"]], "indefinite")
  expect_true(any(invalid[["spectral_values"]] < 0))
  expect_null(.covariance_sampling_factor(invalid))
  expect_error(
    .known_v_as_matrix(indefinite),
    "positive semidefinite"
  )
})


test_that("auto-whitened low-rank known V retains the intended likelihood", {

  yi <- c(0.10, 0.20, -0.05, 0.15)
  V  <- .low_rank_known_v_test_matrix()
  tau <- 0.10
  prior_tau <- BayesTools::prior(
    distribution = "spike",
    parameters   = list(location = tau)
  )

  expect_warning(
    object <- brma.mv(
      yi                        = yi,
      V                         = V,
      prior_heterogeneity       = prior_tau,
      known_v_parameterization  = "auto",
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_priors               = TRUE
    ),
    "positive semidefinite"
  )

  known_V <- .data_known_v_data(object[["data"]])
  expect_identical(.known_v_materialize(known_V), V)
  expect_identical(.known_v_effective_backend(known_V), "whitened")
  expect_true(all(vapply(
    .known_v_backend_blocks(known_V, "whitened"),
    function(block) all(block[["variance"]] >= 0),
    logical(1)
  )))

  setup <- list(
    fit               = NULL,
    priors            = object[["priors"]],
    data              = object[["data"]],
    yi                = yi,
    K                 = length(yi),
    S                 = 1L,
    mu                = matrix(0, nrow = 1L, ncol = length(yi)),
    tau_within        = matrix(tau, nrow = 1L, ncol = length(yi)),
    posterior_samples = matrix(numeric(0), nrow = 1L, ncol = 0L),
    weights           = NULL,
    is_weightfunction = FALSE,
    outcome_type      = "norm",
    effect_direction  = "positive"
  )
  expected <- mvtnorm::dmvnorm(
    x     = yi,
    mean  = rep(0, length(yi)),
    sigma = V + diag(tau^2, nrow(V)),
    log   = TRUE
  )
  expect_equal(
    .log_lik_known_v_joint_sum_from_setup(setup),
    expected,
    tolerance = 1e-12
  )
})
