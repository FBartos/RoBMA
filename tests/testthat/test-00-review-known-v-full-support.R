test_that("declared singleton generators certify full-row support before null projection", {

  V <- tcrossprod(cbind(1, 3500 + 0:2)) / 3500^2
  factorization <- .covariance_factorization(V)
  permutations <- list(1:3, c(3L, 1L, 2L), c(2L, 3L, 1L))
  for (order in permutations) {
    f <- .covariance_factorization(V[order, order])
    for (scale in c(1, 1e-30, 1e30)) {
      generators <- diag(c(1, -2, 3) * scale)[order, c(3L, 1L, 2L)]
      expect_true(.known_v_nullspace_is_covered(f, generators))
      expect_true(.known_v_nullspace_is_covered(f,
        cbind(generators, generators[, 1L], c(1, 2, 3))))
    }
  }
  expect_false(.known_v_nullspace_is_covered(factorization,
    cbind(diag(3)[, 1:2], diag(3)[, 1:2])))
  expect_false(.known_v_nullspace_is_covered(factorization,
    cbind(1, 3500 + 0:2)))
  expect_error(.known_v_nullspace_is_covered(factorization,
    cbind(diag(3), c(Inf, 0, 0))))
})

test_that("positive independent integrated effects regularize the required ill-conditioned V", {

  V <- tcrossprod(cbind(1, 3500 + 0:2)) / 3500^2
  dat <- data.frame(estimate = factor(1:3))
  args <- list(yi = rep(0, 3), V = V, random = ~ 1 | estimate, data = dat,
    prior_effect = prior("point", list(0)),
    prior_heterogeneity = prior_random(sd = prior("point", list(.2))),
    measure = "GEN", prior_unit_information_sd = 1, only_priors = TRUE, silent = TRUE)
  ordinary <- tryCatch(suppressWarnings(do.call(brma.mv, args)), error = identity)
  expect_s3_class(ordinary, "brma")
  args$prior_bias <- prior_weightfunction("one-sided", .5, wf_fixed(c(1, .5)),
    model = selection_model(estimate_random_effects = "integrate"))
  selected <- tryCatch(suppressWarnings(do.call(bselmodel.mv, args)), error = identity)
  expect_s3_class(selected, "brma")
  args$prior_bias <- prior_weightfunction("one-sided", .5, wf_fixed(c(1, .5)),
    model = selection_model(estimate_random_effects = "condition",
      known_sampling_variance = "integrate"))
  expect_error(suppressWarnings(do.call(bselmodel.mv, args)),
    "do not guarantee a nondegenerate covariance", fixed = TRUE)
  args$prior_bias <- prior_weightfunction("one-sided", .5, wf_fixed(c(1, .5)),
    model = selection_model(estimate_random_effects = "condition",
      known_sampling_variance = "condition"))
  expect_s3_class(suppressWarnings(do.call(bselmodel.mv, args)), "brma")
  expect_gt(min(eigen(V + diag(.04, 3), symmetric = TRUE, only.values = TRUE)$values), .039)
})
