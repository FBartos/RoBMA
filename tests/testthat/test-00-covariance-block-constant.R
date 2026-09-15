# Exact recovery of diagonal-plus-block-constant structure from a plain
# covariance block. Acceptance is an algebraic identity certified by
# reconstruction, never a numerical-rank decision, so every case below states
# the exact structure it expects rather than a tolerance on eigenvalues.

.block_constant_covariance <- function(variance, membership, rho) {

  correlation <- matrix(0, length(variance), length(variance))
  for (level in rev(seq_along(rho))) {
    same <- outer(membership[[level]], membership[[level]], `==`)
    correlation[same] <- rho[[level]]
  }
  diag(correlation) <- 1
  covariance <- correlation * tcrossprod(sqrt(variance))
  .known_v_exact_symmetrize(covariance)
}


.reconstruct_block_constant <- function(factor) {

  out <- tcrossprod(factor[["loading"]])
  diag(out) <- diag(out) + factor[["diagonal"]]
  out
}


.assink_covariance <- function(rho = c(.7, .5)) {

  data("dat.assink2016", package = "metadat", envir = environment())
  V <- suppressWarnings(metafor::vcalc(vi, cluster = study, type = deltype,
    obs = esid, rho = rho, data = dat.assink2016))
  V <- matrix(as.numeric(V), nrow = nrow(V))
  .known_v_exact_symmetrize(V)
}


test_that("compound symmetry recovers exactly as diagonal plus rank one", {

  variance <- c(.04, .09, .05, .08)
  for (rho in c(.1, .5, .999)) {
    covariance <- .block_constant_covariance(
      variance, list(rep(1L, 4L)), rho)
    factor <- .covariance_block_constant_factor(covariance)

    expect_false(is.null(factor))
    expect_identical(factor[["rank"]], 1L)
    expect_identical(factor[["support"]], "chain")
    expect_equal(factor[["levels"]], rho, tolerance = 1e-14)
    expect_true(all(factor[["diagonal"]] > 0))
    expect_lte(
      max(abs(.reconstruct_block_constant(factor) - covariance)) /
        max(abs(covariance)),
      8 * nrow(covariance) * .Machine$double.eps
    )
    # The residual diagonal is the complement of the shared level.
    expect_equal(factor[["diagonal"]], (1 - rho) * variance, tolerance = 1e-14)
  }
})


test_that("every positive two-row correlation is exactly rank one", {

  for (rho in c(.05, .4, .95)) {
    covariance <- .block_constant_covariance(c(.04, .09), list(c(1L, 1L)), rho)
    factor <- .covariance_block_constant_factor(covariance)
    expect_identical(factor[["rank"]], 1L)
    expect_lte(max(abs(.reconstruct_block_constant(factor) - covariance)),
               16 * .Machine$double.eps * max(abs(covariance)))
  }
})


test_that("singleton groups are absorbed and one-type blocks stay rank one", {

  # Two rows of one type, one row of its own type: the singleton contributes
  # only to the diagonal, so the rank is two, not three.
  covariance <- .block_constant_covariance(
    c(.04, .05, .06), list(c(1L, 1L, 2L), rep(1L, 3L)), c(.7, .5))
  factor <- .covariance_block_constant_factor(covariance)
  expect_identical(factor[["rank"]], 2L)
  expect_identical(factor[["support"]], "chain")
  expect_identical(sort(as.integer(colSums(factor[["loading"]] != 0))), c(2L, 3L))

  # One type only: the two levels collapse into a single column.
  one_type <- .block_constant_covariance(
    c(.04, .05, .06), list(rep(1L, 3L), rep(1L, 3L)), c(.7, .7))
  collapsed <- .covariance_block_constant_factor(one_type)
  expect_identical(collapsed[["rank"]], 1L)
  expect_equal(collapsed[["levels"]], .7, tolerance = 1e-14)
})


test_that("nested vcalc types recover with the documented per-block ranks", {

  V <- .assink_covariance()
  blocks <- Filter(function(index) length(index) > 1L, .known_v_block_indices(V))
  observed <- vapply(blocks, function(index) {
    factor <- .covariance_block_constant_factor(V[index, index, drop = FALSE])
    if (is.null(factor)) return(NA_integer_)
    factor[["rank"]]
  }, integer(1))
  sizes <- lengths(blocks)

  expect_false(anyNA(observed))
  expect_identical(observed[sizes == 22L], 3L)
  expect_identical(observed[sizes == 16L], 4L)
  expect_identical(observed[sizes == 5L & observed > 1L], 2L)
  expect_identical(sort(observed), sort(c(3L, 4L, 2L, rep(1L, 11L))))

  for (index in blocks) {
    covariance <- V[index, index, drop = FALSE]
    factor <- .covariance_block_constant_factor(covariance)
    expect_lte(
      max(abs(.reconstruct_block_constant(factor) - covariance)) /
        max(abs(covariance)),
      8 * length(index) * .Machine$double.eps
    )
    expect_true(all(factor[["diagonal"]] > 0))
    shape <- if (factor[["rank"]] >= 3L) "tree" else "chain"
    expect_identical(factor[["support"]], shape)
  }
})


test_that("recovery does not depend on row order", {

  V <- .assink_covariance()
  blocks <- Filter(function(index) length(index) > 1L, .known_v_block_indices(V))
  set.seed(4016)
  for (index in blocks) {
    covariance <- V[index, index, drop = FALSE]
    reference  <- .covariance_block_constant_factor(covariance)
    for (replicate in seq_len(30L)) {
      permutation <- sample(length(index))
      permuted <- .covariance_block_constant_factor(
        covariance[permutation, permutation, drop = FALSE])
      expect_identical(permuted[["rank"]], reference[["rank"]])
      expect_identical(permuted[["support"]], reference[["support"]])
      expect_equal(permuted[["levels"]], reference[["levels"]], tolerance = 1e-14)
      expect_equal(permuted[["diagonal"]], reference[["diagonal"]][permutation],
                   tolerance = 1e-14)
    }
  }
})


test_that("recovery absorbs vcalc rounding but not a perturbed covariance", {

  V <- .assink_covariance()
  blocks <- Filter(function(index) length(index) > 1L, .known_v_block_indices(V))
  reference <- vapply(blocks, function(index) {
    .covariance_block_constant_factor(V[index, index, drop = FALSE])[["rank"]]
  }, integer(1))

  set.seed(11)
  perturb <- function(V, scale) {
    lower <- lower.tri(V)
    V[lower] <- V[lower] * (1 + scale)
    V[upper.tri(V)] <- t(V)[upper.tri(V)]
    V
  }
  ulps <- perturb(V, sample(c(-2, -1, 1, 2), sum(lower.tri(V)), TRUE) *
                    .Machine$double.eps)
  observed <- vapply(blocks, function(index) {
    factor <- .covariance_block_constant_factor(ulps[index, index, drop = FALSE])
    if (is.null(factor)) NA_integer_ else factor[["rank"]]
  }, integer(1))
  expect_identical(observed, reference)

  coarse <- perturb(V, 1e-6 * stats::rnorm(sum(lower.tri(V))))
  declined <- vapply(blocks, function(index) {
    is.null(.covariance_block_constant_factor(coarse[index, index, drop = FALSE]))
  }, logical(1))
  # Only two-row blocks survive a perturbation of that size: any positive
  # two-row correlation is still exactly diagonal plus rank one.
  expect_identical(declined, lengths(blocks) > 2L)
})


test_that("structures that are not block constant on nested partitions decline", {

  # A negative increment: the across-type correlation exceeds the within-type
  # correlation, so no non-negative level decomposition exists.
  expect_null(.covariance_block_constant_factor(.block_constant_covariance(
    c(.04, .05, .06, .07), list(c(1L, 1L, 2L, 2L), rep(1L, 4L)), c(.5, .7))))

  # Three levels that do not form nested partitions.
  general <- diag(c(.04, .05, .06))
  general[1L, 2L] <- general[2L, 1L] <- .01
  general[1L, 3L] <- general[3L, 1L] <- .01
  general[2L, 3L] <- general[3L, 2L] <- .02
  expect_null(.covariance_block_constant_factor(general))

  # AR(1) in time: three distinct levels on a path, not on a partition.
  variance <- c(.03, .07, .05, .04)
  autoregressive <- outer(seq_along(variance), seq_along(variance),
                          function(i, j) .9^abs(i - j)) * tcrossprod(sqrt(variance))
  expect_null(.covariance_block_constant_factor(
    .known_v_exact_symmetrize(autoregressive)))

  # A shared-control structure: the covariance follows the control arm, which
  # is not constant within any partition of the rows.
  control <- matrix(c(
    .05, .02, .02, 0,
    .02, .06, .02, 0,
    .02, .02, .07, .03,
    0, 0, .03, .08), 4L, byrow = TRUE)
  expect_null(.covariance_block_constant_factor(control))

  # Negative correlation.
  negative <- diag(c(.04, .09))
  negative[1L, 2L] <- negative[2L, 1L] <- -.006
  expect_null(.covariance_block_constant_factor(negative))
})


test_that("an exactly singular block stays with the exact rank-one detector", {

  variance <- c(.04, .09, .05)
  singular <- tcrossprod(sqrt(variance))
  expect_null(.covariance_block_constant_factor(singular))
  expect_false(is.null(.covariance_exact_rank_one_factor(singular)))

  # The same holds for the matrix vcalc(rho = 1) produces.
  V <- .assink_covariance(rho = c(1, 1))
  for (index in Filter(function(i) length(i) > 1L, .known_v_block_indices(V))) {
    block <- V[index, index, drop = FALSE]
    expect_null(.covariance_block_constant_factor(block))
  }
})


test_that("a recovered factor is attached per block and leaves V authoritative", {

  V <- .assink_covariance()
  known_V <- .known_v_canonicalize(V)

  expect_identical(.known_v_storage(known_V), "blocks")
  expect_true(.known_v_has_certified_factor(known_V))
  expect_identical(.known_v_certified_factor_status(known_V),
                   "recovered_block_constant")
  expect_length(.known_v_certified_factor_dense_rows(known_V), 0L)
  # V stays the stored covariance, bit for bit.
  expect_identical(.known_v_covariance_matrix(known_V), V)

  loading <- .known_v_certified_factor_loading(known_V)
  expect_identical(ncol(loading), 20L)
  reconstruction <- tcrossprod(loading)
  diag(reconstruction) <- diag(reconstruction) +
    .known_v_certified_factor_diagonal(known_V)
  expect_lte(max(abs(reconstruction - V)),
             8 * 22 * .Machine$double.eps * max(abs(V)))

  # The compact representation does not carry a dense K x p loading.
  expect_null(known_V[["certified_factor"]][["loading"]])
  expect_true(all(vapply(.known_v_certified_factor_blocks(known_V),
    function(block) nrow(block[["loading"]]) == length(block[["index"]]),
    logical(1))))
})


test_that("a mixed matrix recovers block by block", {

  variance <- c(.04, .05, .06, .03, .07, .05)
  compound <- .block_constant_covariance(variance[1:3], list(rep(1L, 3L)), .6)
  autoregressive <- .known_v_exact_symmetrize(
    outer(1:3, 1:3, function(i, j) .8^abs(i - j)) * tcrossprod(sqrt(variance[4:6])))
  V <- matrix(0, 6L, 6L)
  V[1:3, 1:3] <- compound
  V[4:6, 4:6] <- autoregressive

  known_V <- .known_v_canonicalize(V)
  expect_identical(.known_v_certified_factor_status(known_V),
                   "recovered_block_constant")
  expect_length(.known_v_certified_factor_blocks(known_V), 1L)
  expect_identical(.known_v_certified_factor_dense_rows(known_V), 4:6)
  expect_identical(.known_v_covariance_matrix(known_V), V)

  object <- bselmodel.mv(
    yi = c(.1, -.2, .3, .05, .25, -.15), V = V,
    data = data.frame(id = 1:6), measure = "GEN",
    prior_unit_information_sd = 1, only_priors = TRUE, silent = TRUE)
  plan <- .data_selection_execution_plan(object[["data"]])
  expect_identical(plan[["row_blocks"]], list(1:3, 4:6))
  expect_identical(plan[["block_methods"]], c("rank_one", "dense"))
  # The dense block keeps the supplied entries, not a reconstruction.
  expect_identical(
    .selection_joint_sampling_block(plan[["sampling"]], 4:6), autoregressive)
  expect_identical(
    .selection_joint_sampling_block(plan[["sampling"]], 1:3), compound)
})
