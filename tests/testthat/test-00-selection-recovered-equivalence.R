# A recovered factor is a representation of the supplied covariance, so the
# selected law it produces must match the dense route it replaces within the
# integration diagnostics both routes report. The block evaluators are compared
# directly on the same states, and the fitted plumbing is checked separately.

.recovered_equivalence_blocks <- function(rho = c(.7, .5)) {

  # Two dependency blocks. The first has three effect-size types of two rows
  # each, so it recovers as a rank-four tree; the second has one type and
  # recovers as rank one.
  study <- rep(c("a", "b"), c(6L, 3L))
  type  <- c("x", "x", "y", "y", "z", "z", "w", "w", "w")
  variance <- seq(0.03, 0.11, length.out = length(study))
  correlation <- matrix(0, length(study), length(study))
  correlation[outer(study, study, `==`)] <- rho[[2L]]
  correlation[outer(study, study, `==`) & outer(type, type, `==`)] <- rho[[1L]]
  diag(correlation) <- 1

  list(
    yi = c(-.30, .12, .41, -.08, .22, .05, .18, -.21, .33),
    V  = .known_v_exact_symmetrize(correlation * tcrossprod(sqrt(variance)))
  )
}


.recovered_equivalence_context <- function(sei, draws, weights = c(1, .7, .35)) {

  prior <- BayesTools::prior_weightfunction("one-sided", steps = c(.025, .05),
    weights = BayesTools::wf_fixed(weights))
  spec <- .selection_spec(priors = list(outcome = list(bias = prior)),
    yi = rep(0, length(sei)), sei = sei, effect_direction = "positive",
    signed_data = FALSE)
  spec[["omega"]]        <- matrix(weights, draws, length(weights), byrow = TRUE)
  spec[["kernel_mode"]]  <- rep(SELKERNEL_STEP, draws)
  spec[["vector_rule"]]  <- rep(0L, draws)
  spec[["native_cache"]] <- new.env(parent = emptyenv())
  spec
}


.recovered_equivalence_plan <- function(k, method, rank) {

  .selection_joint_execution_plan(
    row_blocks = list(seq_len(k)), block_methods = method,
    factor_ranks = rank, selection_control = set_selection_likelihood_control(),
    sampling = NULL, sampling_factor_blocks = NULL, random_covariance = NULL)
}


test_that("recovered and dense block evaluators agree on the selected law", {

  fixture <- .recovered_equivalence_blocks()
  known_V <- .known_v_canonicalize(fixture[["V"]])
  expect_identical(.known_v_certified_factor_status(known_V),
                   "recovered_block_constant")
  blocks <- .known_v_certified_factor_blocks(known_V)
  expect_identical(vapply(blocks, function(b) ncol(b[["loading"]]), integer(1)),
                   c(4L, 1L))

  set.seed(4021)
  draws     <- 100L
  tolerance <- set_selection_likelihood_control()[["relative_tolerance"]]

  for (block in blocks) {
    index <- block[["index"]]
    k     <- length(index)
    rank  <- ncol(block[["loading"]])
    yi    <- fixture[["yi"]][index]
    sei   <- sqrt(diag(fixture[["V"]])[index])
    context <- .recovered_equivalence_context(sei, draws)

    # The same states for both routes: a posterior-like mean sweep and an
    # extra integrated variance that both representations must carry.
    means <- matrix(stats::rnorm(draws * k, 0.15, 0.10), draws, k)
    extra <- stats::runif(draws, 0.01, 0.06)

    residual_sd <- sqrt(matrix(block[["diagonal"]], draws, k, byrow = TRUE) +
                          matrix(extra, draws, k))
    loading <- matrix(as.numeric(block[["loading"]]), draws, k * rank,
                      byrow = TRUE)
    observed <- if (rank == 1L) {
      .selection_joint_cluster_loglik_block(
        yi = yi, means = means, residual_sd = residual_sd, loading = loading,
        sei = sei, selection_context = context,
        execution_plan = .recovered_equivalence_plan(k, "rank_one", 1L),
        return_normalizer = TRUE)
    } else {
      .selection_joint_factor_loglik_block(
        yi = yi, means = means, residual_sd = residual_sd, loading = loading,
        sei = sei, selection_context = context,
        execution_plan = .recovered_equivalence_plan(k, "factor", rank),
        block_index = 1L, return_normalizer = TRUE)
    }

    # The dense route receives the supplied entries plus the same extra
    # variance, packed in the order the dense kernel expects.
    supplied <- fixture[["V"]][index, index, drop = FALSE]
    pairs <- .selection_joint_lower_pairs(
      .recovered_equivalence_plan(k, "dense", NA_integer_), seq_len(k))
    packed <- matrix(supplied[cbind(pairs[["row_1"]], pairs[["row_2"]])],
                     draws, length(pairs[["row_1"]]), byrow = TRUE)
    diagonal <- pairs[["row_1"]] == pairs[["row_2"]]
    packed[, diagonal] <- packed[, diagonal] + extra
    reference <- .selection_joint_dense_loglik_block(
      yi = yi, means = means, covariance_lower = packed, sei = sei,
      selection_context = context,
      execution_plan = .recovered_equivalence_plan(k, "dense", NA_integer_),
      block_size = k, return_normalizer = TRUE)

    expect_true(all(is.finite(observed[["log_density"]])))
    expect_true(all(is.finite(reference[["log_density"]])))
    # Both routes accept a state only when their own diagnostics are inside
    # the fitted relative tolerance, so agreement is judged on that scale.
    budget <- tolerance + max(observed[["relative_mcse"]],
                              observed[["relative_change"]],
                              reference[["relative_mcse"]])
    expect_lt(max(abs(observed[["log_density"]] - reference[["log_density"]])),
              budget)
  }
})


test_that("recovered and dense zplot block densities agree", {

  fixture <- .recovered_equivalence_blocks()
  known_V <- .known_v_canonicalize(fixture[["V"]])
  blocks  <- .known_v_certified_factor_blocks(known_V)
  z       <- seq(-3, 3, by = .25)
  control <- set_selection_likelihood_control()

  set.seed(9110)
  draws <- 12L
  for (block in blocks) {
    index <- block[["index"]]
    k     <- length(index)
    rank  <- ncol(block[["loading"]])
    sei   <- sqrt(diag(fixture[["V"]])[index])
    context <- .recovered_equivalence_context(sei, draws)
    means <- matrix(stats::rnorm(draws * k, 0.15, 0.10), draws, k)
    extra <- stats::runif(draws, 0.01, 0.06)

    supplied <- fixture[["V"]][index, index, drop = FALSE]
    pairs <- .selection_joint_lower_pairs(
      .recovered_equivalence_plan(k, "dense", NA_integer_), seq_len(k))
    packed <- matrix(supplied[cbind(pairs[["row_1"]], pairs[["row_2"]])],
                     draws, length(pairs[["row_1"]]), byrow = TRUE)
    diagonal <- pairs[["row_1"]] == pairs[["row_2"]]
    packed[, diagonal] <- packed[, diagonal] + extra

    factors <- list(
      residual_sd = sqrt(matrix(block[["diagonal"]], draws, k, byrow = TRUE) +
                           matrix(extra, draws, k)),
      loading = matrix(as.numeric(block[["loading"]]), draws, k * rank,
                       byrow = TRUE),
      loading_support = block[["loading"]] != 0
    )

    observed <- .zplot_joint_block(z, means, packed, sei, context, FALSE,
      control, new.env(parent = emptyenv()), factors)
    reference <- .zplot_joint_block(z, means, packed, sei, context, FALSE,
      control, new.env(parent = emptyenv()), NULL)

    expect_true(all(is.finite(observed[["density"]])))
    expect_true(all(observed[["density"]] >= 0))
    peak <- max(reference[["density"]])
    expect_gt(peak, 0)
    budget <- control[["relative_tolerance"]] +
      max(observed[["relative_mcse"]], reference[["relative_mcse"]])
    expect_lt(max(abs(observed[["density"]] - reference[["density"]])) / peak,
              budget)
  }
})


test_that("recovery leaves Gaussian known-V models byte-identical", {

  fixture <- .recovered_equivalence_blocks()
  dat <- data.frame(
    yi    = fixture[["yi"]],
    study = factor(rep(c("a", "b"), c(6L, 3L))),
    esid  = 1:9
  )

  gaussian_data <- function(recovered) {
    build <- function() {
      brma.mv(yi = yi, V = fixture[["V"]], random = ~ 1 | study / esid,
        data = dat, measure = "GEN", prior_unit_information_sd = 1,
        only_priors = TRUE, silent = TRUE)
    }
    if (!recovered) {
      testthat::local_mocked_bindings(
        .covariance_block_constant_factor = function(...) NULL,
        .package = "RoBMA")
    }
    object <- build()
    list(fit_data = .create_fit_data(object$data, object$priors),
         syntax   = .create_model_syntax(object$data, object$priors),
         known_V  = .data_known_v_data(object$data))
  }

  recovered <- gaussian_data(TRUE)
  dense     <- gaussian_data(FALSE)

  expect_identical(recovered$fit_data, dense$fit_data)
  expect_identical(recovered$syntax, dense$syntax)
  # V stays the stored covariance and the supplied blocks are unchanged.
  expect_identical(.known_v_blocks(recovered$known_V), .known_v_blocks(dense$known_V))
  expect_identical(.known_v_covariance_matrix(recovered$known_V), fixture[["V"]])
  expect_true(.known_v_has_certified_factor(recovered$known_V))
  expect_false(.known_v_has_certified_factor(dense$known_V))
})


test_that("a recovered selection plan routes its blocks and keeps V exact", {

  fixture <- .recovered_equivalence_blocks()
  dat <- data.frame(yi = fixture[["yi"]], id = seq_along(fixture[["yi"]]))
  prior <- BayesTools::prior_weightfunction("one-sided", steps = c(.025, .05),
    weights = BayesTools::wf_fixed(c(1, .7, .35)))

  build <- function(recovered) {
    if (!recovered) {
      testthat::local_mocked_bindings(
        .covariance_block_constant_factor = function(...) NULL,
        .package = "RoBMA")
    }
    bselmodel.mv(yi = yi, V = fixture[["V"]], data = dat, measure = "GEN",
      prior_unit_information_sd = 1, prior_bias = prior,
      effect_direction = "positive", only_priors = TRUE, silent = TRUE)
  }

  recovered <- .data_selection_execution_plan(build(TRUE)[["data"]])
  dense     <- .data_selection_execution_plan(build(FALSE)[["data"]])

  expect_identical(recovered[["row_blocks"]], list(1:6, 7:9))
  expect_identical(recovered[["block_methods"]], c("factor", "rank_one"))
  expect_identical(recovered[["factor_ranks"]], c(4L, 1L))
  expect_identical(dense[["row_blocks"]], recovered[["row_blocks"]])
  expect_identical(dense[["block_methods"]], c("dense", "dense"))

  # Both plans materialize the supplied entries, not a reconstruction.
  for (plan in list(recovered, dense)) {
    for (rows in plan[["row_blocks"]]) {
      expect_identical(.selection_joint_sampling_block(plan[["sampling"]], rows),
                       fixture[["V"]][rows, rows, drop = FALSE])
    }
  }
})
