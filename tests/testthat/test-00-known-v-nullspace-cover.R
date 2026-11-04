context("(00) structural null-space coverage of singular known-V blocks")

# The integrated variation of a singular sampling covariance V is given by
# generators (the columns of a factor X of its covariance, X X'); it regularizes
# V when range(X) contains null(V). The null basis of V carries the error of
# its eigensolve, so a direction that only rounding noise places in range(X)
# must not count as covered, however the magnitudes of V and X happen to fall.

.covered <- function(V, generators) {

  .known_v_nullspace_is_covered(.covariance_factorization(V), generators)
}


test_that("the correlation condition of the null basis is recorded with the factorization", {

  expect_identical(.covariance_factorization(matrix(1 / 7, 2, 2))[["correlation_condition"]], 1)
  expect_identical(.covariance_factorization(diag(c(1, 1e-20, 0)))[["correlation_condition"]], 1)
  expect_true(is.na(.covariance_factorization(matrix(0, 2, 2))[["correlation_condition"]]))
  # An integrated sampling covariance [1 x][1 x]' has the correlation spectrum
  # of two nearly collinear columns: the condition grows with the offset of x.
  condition <- vapply(list(c(0, 1, 2), c(10, 11, 12), c(2000, 2001, 2002)), function(x) {
    .covariance_factorization(tcrossprod(cbind(1, x)))[["correlation_condition"]]
  }, numeric(1))
  expect_true(all(diff(condition) > 0))
  expect_gt(condition[[3L]], 1e6)
})


test_that("integrated variation on the range of V does not cover its null space", {

  # An integrated paper intercept (support 11') beside V = v 11': the paper
  # intercept spans range(V), and the direction (1, -1) is uncovered. The
  # projected support N'(11')N is rounding noise, 1e-32 for v = 1/97.
  for (v in c(1, 1 / 97, .04, 1 / seq_len(200))) {
    expect_false(.covered(matrix(v, 2, 2), matrix(1, 2, 1)), info = format(v))
  }

  # One paper, three estimates, V = v 11' (null space of dimension two) and an
  # integrated random intercept and slope on x: range(X) = span(1, x) contains
  # a single null direction at most.
  for (case in list(c(1, .1, .2, .3), c(1 / 97, -.861, .377, -.398), c(.04, .1, .2, .3))) {
    expect_false(.covered(case[[1L]] * matrix(1, 3, 3), cbind(1, case[-1L])))
  }

  # V = S = UU' with the sampling covariance as the only integrated variation.
  U <- cbind(1:4, c(2, 1, 0, 1))
  expect_false(.covered(tcrossprod(U), U))
})


test_that("integrated variation on the range of a weakly conditioned V does not cover its null space", {

  # V = [1 x][1 x]' (null space (1, -2, 1)) beside an integrated intercept and
  # slope on x: span(1, x) = range(V). The null basis of a V whose correlation
  # spectrum is this ill-conditioned is off by eps lambda_1 / lambda_r, far
  # beyond eps relative to the singular values of its factor.
  for (x in list(c(0, 1, 2), c(10, 11, 12), c(100, 101, 102), c(1000, 1001, 1002),
                 c(2000, 2001, 2002))) {
    V <- tcrossprod(cbind(1, x))
    expect_false(.covered(V, cbind(1, x)), info = paste(x, collapse = ","))
    expect_false(.covered(V, cbind(1, x - mean(x))), info = paste(x, collapse = ","))
  }
})


test_that("integrated variation that spans the null space covers it, however small the variance", {

  small <- diag(c(1, 1e-20, 0))
  expect_true(.covered(small, diag(3)[, 3, drop = FALSE]))
  expect_false(.covered(small, diag(3)[, 2, drop = FALSE]))
  expect_true(.known_v_nullspace_is_regularized(small, c(FALSE, FALSE, TRUE)))
  expect_false(.known_v_nullspace_is_regularized(small, c(FALSE, TRUE, FALSE)))
  # A source whose variation is tiny but exactly representable still counts, in
  # any scale of the generators.
  expect_true(.covered(small, 1e-30 * diag(3)[, 3, drop = FALSE]))
  expect_true(.covered(small, 1e30 * diag(3)[, 3, drop = FALSE]))

  # Row-level variation beside the intercept and slope covers the null space
  # of the weakly conditioned V.
  for (x in list(c(10, 11, 12), c(2000, 2001, 2002))) {
    expect_true(.covered(tcrossprod(cbind(1, x)), cbind(1, x - mean(x), diag(3))))
  }
  # Row-level generators (an estimate-level random intercept) cover the null
  # space of any block; none, or only zero columns, cover nothing.
  expect_true(.covered(matrix(.04, 2, 2), diag(2)))
  expect_false(.covered(matrix(.04, 2, 2), matrix(0, 2, 1)))
  expect_false(.covered(matrix(.04, 2, 2), matrix(0, 2, 0)))
})


test_that("an independent regularized row never covers the null direction of a dependent block", {

  # V = v 11' on two rows and an independent row, in every order of the rows:
  # only the independent row is regularized, and the null direction (1, -1) of
  # the dependent pair is not covered. The unpermuted order was always refused;
  # in the others the null basis put rounding noise (1e-16) on the regularized
  # row.
  permutations <- list(1:3, c(1, 3, 2), c(2, 1, 3), c(2, 3, 1), c(3, 1, 2), c(3, 2, 1))
  for (v in c(.7, 1e-3, 123.456)) {
    for (p in permutations) {
      V <- matrix(0, 3, 3)
      V[1:2, 1:2] <- v
      V[3, 3] <- 1
      expect_false(.known_v_nullspace_is_regularized(V[p, p], p == 3L),
                   info = paste(v, paste(p, collapse = "")))
    }
  }
})


# ---------------------------------------------------------------------------- #
# Public constructor
# ---------------------------------------------------------------------------- #

.selection_integrated_model <- function(V, random, data) {

  suppressWarnings(bselmodel.mv(
    yi = yi, V = V, random = random, measure = "GEN", data = data,
    prior_unit_information_sd = 1,
    prior_bias = prior_weightfunction("one-sided", .025, wf_fixed(c(1, .2)),
      model = selection_model(other_random_effects = "integrate",
                              known_sampling_variance = "integrate", group = paper)),
    only_priors = TRUE, silent = TRUE
  ))
}


test_that("the constructor refuses integrated variation on the range of a singular V", {

  message <- "do not guarantee a nondegenerate covariance"
  # Two perfectly correlated sampling errors and an integrated paper intercept.
  for (v in c(1, .1, .04, 1 / 97)) {
    expect_error(.selection_integrated_model(matrix(v, 2, 2), ~ 1 | paper,
      data.frame(yi = c(.2, .1), paper = "paper")), message, info = format(v))
  }
  # Three estimates, an integrated intercept and a slope on x.
  for (case in list(c(1, .1, .2, .3), c(1 / 97, -.861, .377, -.398), c(.04, .1, .2, .3))) {
    expect_error(.selection_integrated_model(case[[1L]] * matrix(1, 3, 3), ~ x | paper,
      data.frame(yi = c(.2, .1, .3), x = case[-1L], paper = "paper")), message)
  }
  # V = [1 x][1 x]', whose null direction (1, -2, 1) only an ill-conditioned
  # solve could place in span(1, x).
  for (x in list(c(10, 11, 12), c(100, 101, 102), c(2000, 2001, 2002))) {
    expect_error(.selection_integrated_model(tcrossprod(cbind(1, x)), ~ x | paper,
      data.frame(yi = c(.2, .1, .3), x = x, paper = "paper")), message,
      info = paste(x, collapse = ","))
  }
})


test_that("the constructor accepts integrated variation that covers the null space", {

  # Perfectly correlated sampling errors and an estimate-level random
  # intercept, whose integrated support covers every row.
  dat <- data.frame(yi = c(.2, .1), estimate = factor(1:2), paper = "paper")
  expect_s3_class(
    .selection_integrated_model(matrix(.04, 2, 2), ~ 1 | estimate, dat),
    "brma"
  )
  # The weakly conditioned V with an intercept, a slope and an estimate-level
  # intercept.
  x <- c(2000, 2001, 2002)
  dat <- data.frame(yi = c(.2, .1, .3), x = x, estimate = factor(1:3), paper = "paper")
  expect_s3_class(
    .selection_integrated_model(tcrossprod(cbind(1, x)),
      list(~ x | paper, ~ 1 | estimate), dat),
    "brma"
  )
})
