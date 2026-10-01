# The nested factor quadrature integrates a forest of supports, not only a
# chain. Every case is checked against an independent oracle: the selected
# block normalizer is the expected product of step weights under the block
# Gaussian law, so it is the sum of weighted bin-rectangle probabilities.

.tree_rule_oracle <- function(yi, mean, covariance, sei, omega, spec) {

  bins  <- as.matrix(expand.grid(rep(list(seq_len(spec[["n_bins"]])),
                                     length(yi))))
  total <- 0
  for (point in seq_len(nrow(bins))) {
    index <- bins[point, ]
    total <- total + prod(omega[index]) * suppressWarnings(mvtnorm::pmvnorm(
      lower     = spec[["z_lower"]][index] * sei,
      upper     = spec[["z_upper"]][index] * sei,
      mean      = spec[["sign"]] * mean,
      sigma     = covariance,
      algorithm = mvtnorm::Miwa(steps = 512L)
    )[[1L]])
  }

  mvtnorm::dmvnorm(spec[["sign"]] * yi, spec[["sign"]] * mean, covariance,
                   log = TRUE) +
    sum(log(omega[spec[["obs_bin"]]])) - log(total)
}


.tree_rule_evaluate <- function(loading, weights = c(1, .7, .35),
                                steps = c(.025, .05)) {

  k    <- nrow(loading)
  rank <- ncol(loading)
  yi   <- seq(-0.35, 0.45, length.out = k)
  sei  <- seq(0.18, 0.34, length.out = k)
  residual_sd <- sqrt(seq(0.02, 0.05, length.out = k))
  mean_vector <- rep(0.12, k)

  prior <- BayesTools::prior_weightfunction("one-sided", steps = steps,
    weights = BayesTools::wf_fixed(weights))
  spec <- .selection_spec(priors = list(outcome = list(bias = prior)),
    yi = yi, sei = sei, effect_direction = "positive", signed_data = FALSE)
  context <- spec
  context[["omega"]]        <- matrix(weights, nrow = 1L)
  context[["kernel_mode"]]  <- SELKERNEL_STEP
  context[["vector_rule"]]  <- 0L
  context[["native_cache"]] <- new.env(parent = emptyenv())

  plan <- .selection_joint_execution_plan(
    row_blocks = list(seq_len(k)), block_methods = "factor",
    factor_ranks = as.integer(rank),
    selection_control = set_selection_likelihood_control(),
    sampling = NULL, sampling_factor_blocks = NULL, random_covariance = NULL)

  observed <- .selection_joint_factor_loglik_block(
    yi = yi, means = matrix(mean_vector, 1L, k),
    residual_sd = matrix(residual_sd, 1L, k),
    loading = matrix(as.numeric(loading), 1L, k * rank),
    sei = sei, selection_context = context, execution_plan = plan,
    block_index = 1L, return_normalizer = TRUE)

  covariance <- diag(residual_sd^2, k) + tcrossprod(loading)
  list(
    observed  = observed[["log_density"]][[1L]],
    reference = .tree_rule_oracle(yi, mean_vector, covariance, sei, weights, spec),
    mcse      = observed[["relative_mcse"]][[1L]],
    change    = observed[["relative_change"]][[1L]]
  )
}


test_that("chain supports keep their deterministic quadrature", {

  result <- .tree_rule_evaluate(cbind(
    c(.12, .12, .12, .12), c(.09, .09, .09, 0), c(.11, .11, 0, 0)))

  expect_identical(result[["mcse"]], 0)
  expect_lt(result[["change"]], set_selection_likelihood_control()[["relative_tolerance"]])
  expect_equal(result[["observed"]], result[["reference"]], tolerance = 1e-5)
})


test_that("a root with disjoint children integrates as a tree", {

  # One column over the whole block and two disjoint child columns: the shape
  # vcalc() produces for a study with two effect-size types.
  result <- .tree_rule_evaluate(cbind(
    c(.12, .12, .12, .12), c(.09, .09, 0, 0), c(0, 0, .11, .11)))

  expect_identical(result[["mcse"]], 0)
  expect_lt(result[["change"]], set_selection_likelihood_control()[["relative_tolerance"]])
  expect_equal(result[["observed"]], result[["reference"]], tolerance = 1e-5)
})


test_that("a rank-four tree and a two-root forest integrate deterministically", {

  tree <- .tree_rule_evaluate(cbind(
    rep(.10, 6), c(.08, .08, 0, 0, 0, 0), c(0, 0, .07, .07, 0, 0),
    c(0, 0, 0, 0, .09, .09)), weights = c(1, .5), steps = .025)
  expect_identical(tree[["mcse"]], 0)
  expect_equal(tree[["observed"]], tree[["reference"]], tolerance = 1e-5)

  forest <- .tree_rule_evaluate(cbind(c(.13, .13, 0, 0), c(0, 0, .10, .10)))
  expect_identical(forest[["mcse"]], 0)
  expect_equal(forest[["observed"]], forest[["reference"]], tolerance = 1e-5)
})


test_that("partially overlapping supports are not a forest and fall back", {

  # The nested rule must decline this shape; the tensor or importance route
  # still has to return the same integral.
  result <- .tree_rule_evaluate(cbind(c(.12, .12, .12, 0), c(0, .10, .10, .10)))
  expect_equal(result[["observed"]], result[["reference"]], tolerance = 1e-4)
})


test_that("recovered Assink tree blocks reach the deterministic factor route", {

  data("dat.assink2016", package = "metadat", envir = environment())
  V <- suppressWarnings(metafor::vcalc(vi, cluster = study, type = deltype,
    obs = esid, rho = c(.7, .5), data = dat.assink2016))
  V <- .known_v_exact_symmetrize(matrix(as.numeric(V), nrow = nrow(V)))
  known_V <- .known_v_canonicalize(V)

  prior <- BayesTools::prior_weightfunction("one-sided", steps = c(.025, .05),
    weights = BayesTools::wf_fixed(c(1, .7, .35)))
  blocks <- Filter(function(block) ncol(block[["loading"]]) >= 2L,
                   .known_v_certified_factor_blocks(known_V))
  expect_length(blocks, 3L)

  for (block in blocks) {
    index <- block[["index"]]
    k     <- length(index)
    rank  <- ncol(block[["loading"]])
    yi    <- dat.assink2016$yi[index]
    sei   <- sqrt(diag(V)[index])
    spec <- .selection_spec(priors = list(outcome = list(bias = prior)),
      yi = yi, sei = sei, effect_direction = "positive", signed_data = FALSE)
    context <- spec
    context[["omega"]]        <- matrix(c(1, .7, .35), nrow = 1L)
    context[["kernel_mode"]]  <- SELKERNEL_STEP
    context[["vector_rule"]]  <- 0L
    context[["native_cache"]] <- new.env(parent = emptyenv())
    plan <- .selection_joint_execution_plan(
      row_blocks = list(seq_len(k)), block_methods = "factor",
      factor_ranks = as.integer(rank),
      selection_control = set_selection_likelihood_control(),
      sampling = NULL, sampling_factor_blocks = NULL, random_covariance = NULL)

    result <- .selection_joint_factor_loglik_block(
      yi = yi, means = matrix(rep(.18, k), 1L, k),
      residual_sd = matrix(sqrt(block[["diagonal"]]), 1L, k),
      loading = matrix(as.numeric(block[["loading"]]), 1L, k * rank),
      sei = sei, selection_context = context, execution_plan = plan,
      block_index = 1L, return_normalizer = TRUE)

    # Deterministic quadrature, not the mode-proposal importance sampler.
    expect_identical(result[["relative_mcse"]][[1L]], 0)
    expect_lt(result[["relative_change"]][[1L]],
              set_selection_likelihood_control()[["relative_tolerance"]])
    expect_true(is.finite(result[["log_density"]][[1L]]))
  }
})
