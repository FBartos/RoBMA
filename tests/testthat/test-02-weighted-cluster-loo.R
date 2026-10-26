source(testthat::test_path("common-functions.R"))

test_that("cached weighted cluster LOO agrees with exact Gaussian deletion posteriors", {

  skip_if_missing_fits("weighted_gaussian_cluster_loo")
  object <- load_fit("weighted_gaussian_cluster_loo")
  y <- .outcome_data_yi(object)
  se <- .outcome_data_sei(object)
  weights <- .outcome_data_weights(object)
  cluster <- object$data$outcome$cluster
  Z <- stats::model.matrix(~ 0 + factor(cluster))
  X <- cbind(1, sqrt(.08) * Z)
  d <- se^2 + .08
  exact_score <- exact_cdf <- numeric(length(y))
  for(i in seq_along(y)){
    keep <- seq_along(y) != i
    precision <- diag(ncol(X)) + crossprod(X[keep, , drop = FALSE],
      (weights[keep] / d[keep]) * X[keep, , drop = FALSE])
    variance <- solve(precision)
    mean <- variance %*% crossprod(X[keep, , drop = FALSE], weights[keep] * y[keep] / d[keep])
    location <- as.numeric(X[i, , drop = FALSE] %*% mean)
    latent_variance <- as.numeric(X[i, , drop = FALSE] %*% variance %*% X[i, ])
    constant <- .5 * (1 - weights[i]) * log(2 * pi * d[i]) - .5 * log(weights[i])
    exact_score[i] <- constant + stats::dnorm(y[i], location,
      sqrt(d[i] / weights[i] + latent_variance), log = TRUE)
    exact_cdf[i] <- stats::pnorm(y[i], location, sqrt(d[i] + latent_variance))
  }
  score <- log_lik(object)
  setup <- .estimate_likelihood_setup.brma(object, condition_local_effects = FALSE)
  cdf <- .weighted_normal_cluster_summary_from_setup(setup, "cdf")$cdf
  # Independent raw importance estimates in contiguous within-chain batches
  # estimate Monte Carlo error, rather than demanding equality to exact refits.
  batches <- split(seq_len(nrow(score)), rep(seq_len(12L), each = nrow(score) / 12L))
  batch_score <- batch_cdf <- matrix(NA_real_, length(batches), ncol(score))
  for(b in seq_along(batches)) for(i in seq_len(ncol(score))){
    rows <- batches[[b]]
    log_weight <- -score[rows, i]
    offset <- max(log_weight)
    w <- exp(log_weight - offset)
    batch_score[b, i] <- -offset - log(mean(w))
    batch_cdf[b, i] <- sum(w * cdf[rows, i]) / sum(w)
  }
  score_mcse <- apply(batch_score, 2L, stats::sd) / sqrt(length(batches))
  cdf_mcse <- apply(batch_cdf, 2L, stats::sd) / sqrt(length(batches))
  actual_score <- loo(object)$pointwise[, "elpd_loo"]
  actual_cdf <- stats::pnorm(rstudent(object)$z)
  expect_true(all(abs(actual_score - exact_score) < 6 * score_mcse + 1e-5))
  expect_true(all(abs(actual_cdf - exact_cdf) < 6 * cdf_mcse + 1e-5))
  expect_identical(attr(loo(object), "RoBMA_target")$weighted_cluster_integration, 1L)
  expect_identical(attr(waic(object), "RoBMA_target")$weighted_cluster_integration, 1L)
})
