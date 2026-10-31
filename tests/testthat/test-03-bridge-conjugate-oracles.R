# Independent conjugate posterior rows; these replays are not MCMC samples.
.cleanup_bridge_conjugate_fixture <- function(family) {

  yi <- c(.2, -.1, .3, .4)
  sei <- c(.2, .3, .25, .2)
  prior_mu <- BayesTools::prior("normal", list(mean = 0, sd = .5))
  point <- function(location) BayesTools::prior("point", list(location = location))
  Z <- cbind(c(1, 1, 0, 0), c(0, 0, 1, 1))
  if (family == "normal") {
    M <- diag(sei^2 + .15^2)
    object <- brma(yi = yi, sei = sei, measure = "GEN", prior_effect = prior_mu,
      prior_heterogeneity = point(.15), prior_unit_information_sd = 1,
      seed = 198, only_priors = TRUE, silent = TRUE)
  } else if (family == "cluster") {
    M <- diag(sei^2 + .08) + .08 * tcrossprod(Z)
    object <- brma(yi = yi, sei = sei, cluster = c("a", "a", "b", "b"), measure = "GEN",
      prior_effect = prior_mu, prior_heterogeneity = point(.4),
      prior_heterogeneity_allocation = point(.5), prior_unit_information_sd = 1,
      seed = 198, only_priors = TRUE, silent = TRUE)
  } else {
    M <- diag(sei^2)
    M[1, 2] <- M[2, 1] <- .01
    M[3, 4] <- M[4, 3] <- .01
    object <- brma.mv(yi = yi, V = M, random = NULL, measure = "GEN",
      prior_effect = prior_mu, prior_unit_information_sd = 1,
      known_v_parameterization = "block_mvn", seed = 198,
      only_priors = TRUE, silent = TRUE)
  }
  setup <- review_compile_object(object)
  priors <- unlist(lapply(setup$compiled, `[[`, "prior_list"), recursive = FALSE)
  priors <- c(priors, .create_fit_priors(setup$object$data, setup$object$priors))
  parameter <- if (is.null(setup$object$formula_design$mu)) "mu" else
    review_coefficient_name(setup, "intercept")
  gamma_columns <- if (family == "cluster") review_prior_columns(priors["gamma"]) else character()
  if (family == "cluster") stopifnot(length(gamma_columns) == 2L)
  if (family == "cluster") {
    A <- cbind(1, sqrt(.08) * Z)
    D <- diag(sei^2 + .08)
    prior_variance <- c(.25, 1, 1)
  } else {
    A <- matrix(1, length(yi), 1L)
    D <- M
    prior_variance <- .25
  }
  D_chol <- chol(D)
  inverse_D_A <- backsolve(D_chol, forwardsolve(t(D_chol), A))
  inverse_D_y <- backsolve(D_chol, forwardsolve(t(D_chol), yi))
  precision <- diag(1 / prior_variance, length(prior_variance)) + crossprod(A, inverse_D_A)
  inverse_chol <- backsolve(chol(precision), diag(ncol(A)))
  covariance <- tcrossprod(inverse_chol)
  mean <- as.vector(covariance %*% crossprod(A, inverse_D_y))
  withr::local_seed(198)
  chains <- lapply(seq_len(2L), function(chain) {
    draws <- sweep(matrix(stats::rnorm(2000L * ncol(A)), 2000L, ncol(A)) %*%
      chol(covariance), 2L, mean, "+")
    colnames(draws) <- c(parameter, gamma_columns)
    for (name in names(priors)) if (BayesTools::is.prior.point(priors[[name]]) &&
      !name %in% colnames(draws)) {
      draws <- cbind(draws, stats::setNames(data.frame(rep(priors[[name]]$parameters$location,
        nrow(draws))), name))
    }
    as.matrix(draws)
  })
  monitors <- BayesTools::JAGS_to_monitor(priors)
  additional <- unique(unlist(lapply(setup$compiled, `[[`, "add_parameters")))
  for (monitor in c(monitors, additional)) stopifnot(monitor %in% colnames(chains[[1L]]) ||
    any(startsWith(colnames(chains[[1L]]), paste0(monitor, "["))))
  object <- review_attach_replay(setup, chains)
  contract <- BayesTools::JAGS_validate_fit_contract(object$fit, requires = c("name_encoding",
    "formula_name_map", "formula_design", "parameter_map", "draw_geometry"))
  oracle_chol <- chol(M + .25 * matrix(1, length(yi), length(yi)))
  standardized <- forwardsolve(t(oracle_chol), yi)
  oracle <- -.5 * length(yi) * log(2 * pi) - sum(log(diag(oracle_chol))) -
    .5 * sum(standardized^2)
  list(object = object, oracle = unname(oracle), monitors = monitors,
    columns = colnames(chains[[1L]]), gamma_columns = gamma_columns, contract = contract)
}

test_that("normal clustered and known-V bridge logml agree with conjugate oracles", {

  fixtures <- lapply(c("normal", "cluster", "known_v"), .cleanup_bridge_conjugate_fixture)
  names(fixtures) <- c("normal", "cluster", "known_v")
  for (family in names(fixtures)) {
    fixture <- fixtures[[family]]
    object <- add_marglik(fixture$object, method = "normal", repetitions = 1L,
      parallel = FALSE, cores = 1L)
    repetitions <- object$marglik$repetitions
    expect_identical(nrow(repetitions), 1L, info = family)
    expect_true(all(repetitions$finite), info = family)
    expect_true(all(repetitions$success), info = family)
    expect_true(all(repetitions$within_maxiter), info = family)
    expect_length(repetitions$mcse, 1L)
    expect_true(all(is.finite(repetitions$mcse) & repetitions$mcse >= 0), info = family)
    roundoff <- 64 * .Machine$double.eps * (1 + abs(fixture$oracle))
    expect_true(abs(logml(object) - fixture$oracle) <= 4 * max(repetitions$mcse) + roundoff,
      info = family)
  }
})
