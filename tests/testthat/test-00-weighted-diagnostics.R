context("Weighted diagnostics")

skip_on_cran()


test_that("Outcome likelihood weights default to one and preserve supplied weights", {

  object_unweighted <- brma.norm(
    yi        = c(.10, .20, .30),
    sei       = c(.10, .20, .30),
    only_data = TRUE
  )

  expect_equal(.outcome_data_weights(object_unweighted), c(1, 1, 1))

  object_weighted <- brma.norm(
    yi        = c(.10, .20, .30),
    sei       = c(.10, .20, .30),
    weights   = c(.5, 1, 2),
    only_data = TRUE
  )

  expect_equal(.outcome_data_weights(object_weighted), c(.5, 1, 2))
})


test_that("VIF uses likelihood weights in the meta-regression precision", {

  dat <- data.frame(
    yi  = c(.10, .25, .05, .30, .18, .12),
    sei = c(.12, .18, .14, .20, .16, .22),
    w   = c(3.0, .4, 1.5, .7, 2.2, .9),
    x1  = c(-2, -1, 0, 1, 2, 3),
    x2  = c(3, 1, 2, 5, 4, 7)
  )

  object <- brma.norm(
    yi        = yi,
    sei       = sei,
    weights   = w,
    mods      = ~ x1 + x2,
    data      = dat,
    only_data = TRUE
  )
  object[["priors"]] <- list(outcome = list(bias = NULL))
  object[["fit"]] <- coda::mcmc.list(coda::mcmc(
    matrix(
      rep(.20, 4),
      ncol     = 1,
      dimnames = list(NULL, "tau")
    )
  ))
  class(object) <- c(
    "brma.norm",
    "brma"
  )

  brma_vif <- .compute_vif(object)

  X       <- .get_model_matrix(object)
  assign  <- attr(X, "assign")
  keep    <- assign != 0
  vi      <- dat[["sei"]]^2
  tau2    <- .20^2
  W       <- diag(dat[["w"]] / (vi + tau2), nrow = nrow(dat))
  vcov    <- solve(crossprod(X, W) %*% X)
  vcov    <- vcov[keep, keep, drop = FALSE]
  assign  <- assign[keep]
  R       <- stats::cov2cor(vcov)
  det_R   <- det(R)
  terms_i <- sort(unique(assign))

  expected <- vapply(terms_i, function(term_i) {
    cols <- which(assign == term_i)
    det(R[cols, cols, drop = FALSE]) *
      det(R[-cols, -cols, drop = FALSE]) / det_R
  }, numeric(1))

  W_unweighted    <- diag(1 / (vi + tau2), nrow = nrow(dat))
  vcov_unweighted <- solve(crossprod(X, W_unweighted) %*% X)
  vcov_unweighted <- vcov_unweighted[keep, keep, drop = FALSE]
  R_unweighted    <- stats::cov2cor(vcov_unweighted)
  det_unweighted  <- det(R_unweighted)
  expected_unweighted <- vapply(terms_i, function(term_i) {
    cols <- which(assign == term_i)
    det(R_unweighted[cols, cols, drop = FALSE]) *
      det(R_unweighted[-cols, -cols, drop = FALSE]) / det_unweighted
  }, numeric(1))

  expect_equal(brma_vif[["GVIF"]], expected, tolerance = 1e-12)
  expect_false(isTRUE(all.equal(expected, expected_unweighted, tolerance = 1e-8)))
})


test_that("VIF covariance helper uses observation-specific tau samples", {

  X <- cbind(
    intercept = 1,
    x1        = c(-2, -1, 0, 1, 2),
    x2        = c(2, 1, 3, 5, 4)
  )
  vi <- c(.02, .04, .03, .05, .06)
  weights <- c(1.0, 2.0, .5, 1.5, .8)
  tau_within <- rbind(
    c(.05, .10, .15, .20, .25),
    c(.20, .18, .16, .14, .12)
  )

  brma_vcov <- .vif_vcov_from_tau_samples(
    X                  = X,
    vi                 = vi,
    weights            = weights,
    tau_within_samples = tau_within
  )

  expected <- Reduce(`+`, lapply(seq_len(nrow(tau_within)), function(s) {
    W <- diag(weights / (vi + tau_within[s, ]^2), nrow = length(vi))
    solve(crossprod(X, W) %*% X)
  })) / nrow(tau_within)

  scalar_tau <- sqrt(mean(tau_within^2))
  scalar_W <- diag(weights / (vi + scalar_tau^2), nrow = length(vi))
  scalar_vcov <- solve(crossprod(X, scalar_W) %*% X)

  expect_equal(unname(brma_vcov), unname(expected), tolerance = 1e-12)
  expect_false(isTRUE(all.equal(brma_vcov, scalar_vcov, tolerance = 1e-8)))
})


test_that("VIF covariance helper uses multilevel block precision", {

  X <- cbind(
    intercept = 1,
    x1        = c(-1, 0, 1, -2, 2),
    x2        = c(2, 4, 3, 1, 5)
  )
  vi <- c(.03, .05, .04, .06, .02)
  weights <- c(1.0, .7, 1.4, .9, 1.2)
  cluster <- c(1, 1, 2, 2, 2)
  tau_within <- rbind(
    c(.10, .12, .08, .11, .09),
    c(.15, .14, .13, .12, .11)
  )
  tau_between <- rbind(
    c(.20, .20, .15, .15, .15),
    c(.10, .10, .18, .18, .18)
  )

  brma_vcov <- .vif_vcov_from_tau_samples(
    X                   = X,
    vi                  = vi,
    weights             = weights,
    tau_within_samples  = tau_within,
    tau_between_samples = tau_between,
    cluster             = cluster
  )

  block_indices <- split(seq_along(cluster), cluster)
  expected <- Reduce(`+`, lapply(seq_len(nrow(tau_within)), function(s) {
    M <- diag((vi + tau_within[s, ]^2) / weights, nrow = length(vi))
    for (idx in block_indices) {
      M[idx, idx] <- M[idx, idx] + tcrossprod(tau_between[s, idx])
    }
    W <- solve(M)
    solve(crossprod(X, W) %*% X)
  })) / nrow(tau_within)

  diagonal_only <- .vif_vcov_from_tau_samples(
    X                   = X,
    vi                  = vi,
    weights             = weights,
    tau_within_samples  = tau_within,
    tau_between_samples = tau_between
  )

  expect_equal(unname(brma_vcov), unname(expected), tolerance = 1e-12)
  expect_false(isTRUE(all.equal(brma_vcov, diagonal_only, tolerance = 1e-8)))
})
