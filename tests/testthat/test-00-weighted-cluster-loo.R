test_that("cached weighted criteria track their integration source", {

  expect_true("R/log-lik-weighted-cluster.R" %in% .fit_cache_required_source_files())
})

test_that("weighted Gaussian deletion matches independent integrals and the fitted joint likelihood", {

  r <- c(2, -.3, .7)
  d <- c(.3, .9, .5)
  for(a in list(c(.5, 1.3, .4), c(0, 0, 0))){
    for(w in list(c(1, 1, 1), c(1.000001, 1, 1), c(.4, 2, 1.3), c(1e-5, 10, .25))){
      value <- .weighted_normal_cluster_deletion(matrix(r, 1L), matrix(sqrt(d), 1L),
        matrix(a, 1L), w, list(1:3))
      full <- .log_lik_cluster_norm_stable(r, sqrt(d), a, w)
      for(i in seq_along(r)){
        keep <- seq_along(r) != i
        v <- 1 / (1 + sum(w[keep] * a[keep]^2 / d[keep]))
        m <- v * sum(w[keep] * a[keep] * r[keep] / d[keep])
        oracle <- log(stats::integrate(function(g){
          exp(w[i] * stats::dnorm(r[i], a[i] * g, sqrt(d[i]), log = TRUE) +
                stats::dnorm(g, m, sqrt(v), log = TRUE))
        }, -Inf, Inf, rel.tol = 1e-11)$value)
        deleted <- .log_lik_cluster_norm_stable(r, sqrt(d), a, replace(w, i, 0))
        expect_equal(value$log_score[1L, i], oracle, tolerance = 1e-10)
        expect_equal(value$log_score[1L, i], full - deleted, tolerance = 1e-10)
        expect_equal(value$residual[1L, i], r[i] - a[i] * m, tolerance = 1e-12)
        expect_equal(exp(2 * value$log_sd[1L, i]), d[i] + a[i]^2 * v, tolerance = 1e-12)
      }
    }
  }
})

test_that("weighted deletion preserves singleton priors and whole-row replication", {

  r <- c(.2, -.4, 1)
  sd <- c(.1, .3, .2)
  a <- c(.4, .8, 0)
  w <- c(2, 1, 3)
  value <- .weighted_normal_cluster_deletion(matrix(r, 1L), matrix(sd, 1L),
    matrix(a, 1L), w, as.list(1:3))
  expect_equal(exp(2 * value$log_sd), matrix(sd^2 + a^2, 1L), tolerance = 1e-13)
  expect_equal(value$shift, matrix(0, 1L, 3L))
  for(i in 1:3){
    expect_equal(value$log_score[1L, i],
      .log_lik_cluster_norm_stable(rep(r[i], w[i]), rep(sd[i], w[i]), rep(a[i], w[i])),
      tolerance = 1e-12)
  }
})

test_that("weighted deletion refuses representability limits and preserves ordinary cancellation", {

  for(power in c(100, 150)){
    a <- rep(10^power, 3L)
    sd <- rep(10^(-power), 3L)
    expect_error(.weighted_normal_cluster_deletion(matrix(a, 1L), matrix(sd, 1L),
      matrix(a, 1L), c(.5, 1, 2), list(1:3)),
      class = "RoBMA_weighted_gaussian_unavailable")
  }
  for(w in c(2^120, 1e200)){
    expect_error(.weighted_normal_cluster_deletion(matrix(c(1, 2, 3), 1L),
      matrix(1, 1L, 3L), matrix(1, 1L, 3L), c(w, 1, 2), list(1:3)),
      class = "RoBMA_weighted_gaussian_unavailable")
  }
  expect_error(.weighted_normal_cluster_deletion(matrix(c(1, 1e300, -1e300, 1), 1L),
    matrix(1, 1L, 4L), matrix(1, 1L, 4L), rep(1, 4L), list(1:4)),
    class = "RoBMA_weighted_gaussian_unavailable")
  expect_error(.weighted_normal_cluster_deletion(matrix(c(0, 2^52 + 1, -2^52), 1L),
    matrix(c(.01, 1, 1), 1L), matrix(1, 1L, 3L), c(2, 1, 1), list(1:3)),
    class = "RoBMA_weighted_gaussian_unavailable")
  equivalent <- .weighted_normal_cluster_deletion(matrix(c(0, 1, -1), 1L),
    matrix(c(1, 1, 2), 1L), matrix(1, 1L, 3L), c(2, 1, 4), list(1:3))
  expect_equal(equivalent$shift[1L, 1L], 0, tolerance = 0)
  unequal <- .weighted_normal_cluster_deletion(matrix(c(0, 1, -4), 1L),
    matrix(c(1, 1, 2), 1L), matrix(1, 1L, 3L), c(2, 1, 1), list(1:3))
  expect_equal(unequal$shift[1L, 1L], 0, tolerance = 0)
})

test_that("large locations are refused before deletion loses observation-scale information", {

  setup <- list(S = 1L, K = 2L, yi = c(1, 0), mu = matrix(-2^100, 1L, 2L),
    sei = c(2^-40, 1), tau_within = matrix(0, 1L, 2L), tau_between = matrix(1, 1L, 2L),
    weights = c(1, 2^120), cluster = list(1:2), data = list(), effect_direction = "positive")
  expect_error(.weighted_normal_cluster_summary_from_setup(setup),
    class = "RoBMA_weighted_gaussian_unavailable")
  setup$yi <- rep(2^52 + 1, 2L)
  setup$mu[,] <- 2^52
  setup$sei <- c(1, 1)
  setup$weights <- c(2, 1)
  expect_error(.weighted_normal_cluster_summary_from_setup(setup),
    class = "RoBMA_weighted_gaussian_unavailable")
  setup$yi <- rep(101, 2L)
  setup$mu[,] <- 100
  actual <- .weighted_normal_cluster_summary_from_setup(setup)
  expected <- .weighted_normal_cluster_deletion(matrix(1, 1L, 2L), matrix(1, 1L, 2L),
    matrix(1, 1L, 2L), setup$weights, list(1:2))
  expect_equal(actual$log_density, expected$log_score, tolerance = 1e-12)
  expect_equal(actual$cdf, stats::pnorm(expected$residual, sd = exp(expected$log_sd)), tolerance = 1e-12)
  object <- brma(yi = setup$yi, sei = setup$sei, cluster = c(1, 1),
    weights = setup$weights, only_data = TRUE, measure = "GEN")
  object$priors <- list(outcome = list(bias = NULL))
  setup$data <- object$data
  setup$priors <- object$priors
  setup$outcome_type <- "norm"
  setup$is_weightfunction <- FALSE
  moments <- .loo_predictive_moments_estimate(object, setup, NULL, matrix(1, 1L, 2L),
    predictive_summary = .weighted_normal_cluster_summary_from_setup(setup))
  expect_equal(moments$resid, as.numeric(expected$residual), tolerance = 1e-12)
  expect_equal(moments$se, as.numeric(exp(expected$log_sd)), tolerance = 1e-12)
})

test_that("weighted deletion handles flagged draws and dominating retained powers", {

  residual <- rbind(c(.2, .3, -.1), c(0, 1, -1 + 1e-10), c(1, 1.5, 1.5))
  weights <- c(2, 1, 1)
  value <- .weighted_normal_cluster_deletion(residual, matrix(1, 3L, 3L),
    matrix(1, 3L, 3L), weights, list(1:3))
  for(i in 1:3){
    precision <- 1 + sum(weights[-i])
    shift <- rowSums(sweep(residual[, -i, drop = FALSE], 2L, weights[-i], `*`)) / precision
    expect_equal(value$shift[, i], shift, tolerance = 1e-14)
    expect_equal(value$residual[, i], residual[, i] - shift, tolerance = 1e-14)
  }
  value <- .weighted_normal_cluster_deletion(matrix(c(1, 2, 3), 1L),
    matrix(1, 1L, 3L), matrix(1, 1L, 3L), c(1e12, 1, 2), list(1:3))
  expect_equal(value$shift[1L, 1L], 2, tolerance = 1e-12)
  expect_equal(value$log_sd[1L, 1L], .5 * log(1.25), tolerance = 1e-12)
})

test_that("ordinary weighted deletion does not invoke the fallback", {

  testthat::local_mocked_bindings(
    .weighted_cluster_deletion_fallback = function(...) stop("Unexpected fallback."),
    .package = "RoBMA")
  value <- .weighted_normal_cluster_deletion(matrix(c(.2, .3, .5), 2L, 3L, byrow = TRUE),
    matrix(1, 2L, 3L), matrix(.4, 2L, 3L), c(.5, 1, 2), list(1:3))
  expect_true(all(is.finite(value$log_score)))
})

.weighted_cluster_loo_replay <- function(direction = "positive", weights = c(.5, 1.5, 2, 1)) {

  object <- brma(yi = c(.2, -.3, .5, .1), sei = c(.2, .3, .2, .4),
    cluster = c(1, 1, 2, 3), weights = weights, measure = "GEN",
    prior_unit_information_sd = 1,
    only_priors = TRUE, silent = TRUE)
  setup <- review_compile_object(object)
  attr(setup$object$data, "effect_direction") <- direction
  samples <- cbind(mu = seq(-.1, .1, length.out = 100L), tau = .4, rho = .5,
    `gamma[1]` = 5, `gamma[2]` = -7, `gamma[3]` = 9)
  object <- review_attach_replay(setup, samples)
  class(object) <- c("brma.norm", "brma")
  object
}

test_that("public weighted LOO and diagnostics integrate rather than retain sampled cluster effects", {

  for(direction in c("positive", "negative")){
    object <- .weighted_cluster_loo_replay(direction)
    changed <- object
    draws <- as.matrix(object$fit$mcmc)
    draws[, grep("gamma", colnames(draws))] <- -100
    changed$fit <- BayesTools::JAGS_with_draws(object$fit, coda::mcmc.list(coda::mcmc(draws)))
    expect_false(isTRUE(all.equal(
      .estimate_likelihood_setup.brma(object, condition_local_effects = TRUE)$mu,
      .estimate_likelihood_setup.brma(changed, condition_local_effects = TRUE)$mu)))
    first <- log_lik(object)
    expect_equal(log_lik(changed), first, tolerance = 0)
    target <- attr(first, "RoBMA_target")
    expect_identical(target$latent_effect_handling, "integrated")
    expect_identical(target$weighted_cluster_integration, 1L)
    expect_equal(target$dependency_component_sizes, c(2L, 1L, 1L))
    object <- add_loo(object, r_eff = 1)
    changed <- add_loo(changed, r_eff = 1)
    expect_equal(rstudent(changed), rstudent(object), tolerance = 0)
    expect_equal(influence(changed), influence(object), tolerance = 0)
    stale <- object
    attr(stale$loo$estimate, "RoBMA_target")$weighted_cluster_integration <- NULL
    expect_error(loo(stale), "Stored LOO predates integrated weighted Gaussian cluster deletion. Recompute with add_loo().", fixed = TRUE)
    expect_error(rstudent(stale), "Recompute with add_loo()", fixed = TRUE)
    object <- add_waic(object)
    attr(object$waic$estimate, "RoBMA_target")$weighted_cluster_integration <- NULL
    expect_error(waic(object), "Stored WAIC predates integrated weighted Gaussian cluster deletion. Recompute with add_waic().", fixed = TRUE)
  }
})

test_that("weighted observation CDF uses original variance and unit weights meet the unweighted limit", {

  object <- .weighted_cluster_loo_replay()
  setup <- .estimate_likelihood_setup.brma(object, condition_local_effects = FALSE)
  summary <- .weighted_normal_cluster_summary_from_setup(setup)
  changed <- setup
  changed$weights[3L] <- 100
  changed <- .weighted_normal_cluster_summary_from_setup(changed)
  expect_identical(summary$cdf[, 3L], changed$cdf[, 3L])
  expect_equal(summary$variance[, 3L], rep(.2^2 + .4^2, setup$S), tolerance = 1e-12)
  expect_false(isTRUE(all.equal(summary$log_density[, 3L], changed$log_density[, 3L])))
  setup$weights[] <- 1
  integrated <- .weighted_normal_cluster_summary_from_setup(setup)
  unweighted <- .weighted_cluster_loo_replay(weights = rep(1, 4L))
  unweighted <- .estimate_likelihood_setup.brma(unweighted, condition_local_effects = FALSE)
  expect_equal(unname(integrated$log_density), unname(.log_lik_estimate_from_setup(unweighted)), tolerance = 1e-12)
  reference <- .normal_covariance_estimate_target_summary_from_setup(unweighted)
  for(component in names(reference)) expect_equal(integrated[[component]], reference[[component]], tolerance = 1e-12)
})

test_that("weighted location-scale regression scores use the evaluated row-specific design", {

  data <- data.frame(yi = c(.2, -.3, .5, .1), sei = c(.2, .3, .2, .4),
                     w = c(.5, 1.5, 2, 1), study = c(1, 1, 2, 3), x = c(-1, 0, 1, 2), z = c(2, -1, 0, 1))
  object <- brma(yi = yi, sei = sei, weights = w, cluster = study, mods = ~ x, scale = ~ z,
    data = data, measure = "GEN", prior_unit_information_sd = 1, only_priors = TRUE, silent = TRUE)
  compiled <- review_compile_object(object)
  names <- unique(c(unlist(lapply(compiled$compiled, function(item) review_prior_columns(item$prior_list))),
                     "rho", paste0("gamma[", 1:3, "]")))
  samples <- matrix(.1, 4L, length(names), dimnames = list(NULL, names))
  samples[, review_coefficient_name(compiled, "intercept", "log_tau")] <- c(.2, .3, .4, .5)
  samples[, "rho"] <- .4
  object <- review_attach_replay(compiled, samples)
  class(object) <- c("brma.norm", "brma")
  setup <- .estimate_likelihood_setup.brma(object, condition_local_effects = FALSE)
  actual <- log_lik(object)
  for(s in seq_len(setup$S)) for(i in seq_len(setup$K)){
    group <- setup$cluster[[which(vapply(setup$cluster, function(rows) i %in% rows, logical(1)))]]
    r <- setup$yi[group] - setup$mu[s, group]
    sd <- sqrt(setup$sei[group]^2 + setup$tau_within[s, group]^2)
    a <- setup$tau_between[s, group]
    w <- setup$weights[group]
    oracle <- .log_lik_cluster_norm_stable(r, sd, a, w) -
      .log_lik_cluster_norm_stable(r, sd, a, replace(w, match(i, group), 0))
    expect_equal(unname(actual[s, i]), unname(oracle), tolerance = 1e-12)
  }
})

test_that("weighted predictive residual moments refuse extreme location draws", {

  object <- brma(yi = c(2^53, 1), sei = c(1, 1), weights = c(2, 1), cluster = c(1, 1),
    measure = "GEN", only_data = TRUE)
  object$priors <- list(outcome = list(bias = NULL))
  setup <- list(S = 2L, K = 2L, yi = c(2^53, 1), mu = rbind(c(2^53 - 1, 0), c(-2^53, 0)),
    sei = c(1, 1), tau_within = matrix(0, 2L, 2L), tau_between = rbind(c(1, 1), c(0, 0)),
    weights = c(2, 1), cluster = list(1:2), data = object$data, priors = object$priors,
    effect_direction = "positive", outcome_type = "norm", is_weightfunction = FALSE)
  expect_error(.weighted_normal_cluster_summary_from_setup(setup),
    class = "RoBMA_weighted_gaussian_unavailable")
})

test_that("weighted scores and tails refuse unrepresentable covariance variances", {

  for(power in c(-600, 600)){
    setup <- list(S = 1L, K = 2L, yi = c(0, 0), mu = matrix(0, 1L, 2L),
      sei = rep(2^power, 2L), tau_within = matrix(0, 1L, 2L),
      tau_between = matrix(2^600, 1L, 2L), weights = c(1, 2), cluster = list(1:2),
      data = list(), effect_direction = "positive")
    expect_error(.weighted_normal_cluster_summary_from_setup(setup, c("log_density", "cdf")),
      class = "RoBMA_weighted_gaussian_unavailable")
    expect_error(.weighted_normal_cluster_summary_from_setup(setup, "variance"),
      class = "RoBMA_weighted_gaussian_unavailable")
  }
  expect_error(.weighted_normal_cluster_deletion(matrix(c(0, 1e300, -1e300), 1L),
    matrix(c(1, 1.1, 1.2), 1L), matrix(1, 1L, 3L), c(2, 1.21, 1.44), list(1:3)),
    class = "RoBMA_weighted_gaussian_unavailable")
})

test_that("arbitrarily small positive powers retain their determinant contribution", {

  for(w in c(1e-20, 1e-100, 1e-200)){
    actual <- .weighted_normal_cluster_deletion(matrix(0, 1L), matrix(2, 1L),
      matrix(3, 1L), w, list(1L))$log_score[1L, 1L]
    reference <- -.5 * w * log(2 * pi * 4) - .5 * log1p(w * 9 / 4)
    expect_equal(actual / w, reference / w, tolerance = 1e-12)
  }
})
