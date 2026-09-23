.pr82_global_conditioned_setup <- function(S = 5L, K = 4L) {

  object <- bselmodel(
    yi = seq(-.3, .5, length.out = K), sei = seq(.7, 1, length.out = K),
    cluster = rep("shared", K), measure = "GEN", effect_direction = "positive",
    prior_unit_information_sd = 1,
    prior_bias = BayesTools::prior_weightfunction("one-sided", .5,
      BayesTools::wf_fixed(c(1, .4)), model = BayesTools::selection_model(
        other_random_effects = "condition", known_sampling_variance = "condition")),
    only_priors = TRUE
  )
  samples <- cbind(mu = seq(-.1, .2, length.out = S),
                   tau = seq(.2, .5, length.out = S), rho = .3,
                   "gamma[1]" = seq(-.2, .2, length.out = S))
  for (name in c("theta", "sampling_z")) {
    values <- matrix(seq(-.4, .4, length.out = S * K), S, K)
    colnames(values) <- paste0(name, "[", seq_len(K), "]")
    samples <- cbind(samples, values)
  }
  .log_lik_posterior_setup(object$fit, samples, object$data, object$priors, "estimate", NULL)
}

.pr82_dense_sampling_setup <- function(S = 5L, K = 6L) {

  covariance <- .3^abs(outer(seq_len(K), seq_len(K), `-`))
  object <- bselmodel.mv(
    yi = seq(-.3, .5, length.out = K), V = covariance, random = ~ 1 | esid,
    data = data.frame(esid = factor(seq_len(K))), measure = "GEN",
    effect_direction = "positive", prior_unit_information_sd = 1,
    prior_heterogeneity = BayesTools::prior_random(sd = BayesTools::prior("point", list(.3))),
    prior_bias = BayesTools::prior_weightfunction("one-sided", .5,
      BayesTools::wf_fixed(c(1, .4)), model = BayesTools::selection_model(
        known_sampling_variance = "condition")), only_priors = TRUE, silent = TRUE
  )
  fit <- structure(list(), formula_design = object[["formula_design"]],
    prior_list = c(object[["formula_design"]][["mu"]][["prior_list"]],
                   .create_fit_priors(object[["data"]], object[["priors"]])))
  samples <- matrix(seq(-.1, .2, length.out = S), S, 1L,
                     dimnames = list(NULL, "mu_intercept"))
  for (term in .formula_design_random_effects_by_mode(object[["formula_design"]][["mu"]], "sampled")) {
    values <- matrix(seq(-.2, .2, length.out = S * term[["n_groups"]]), S)
    colnames(values) <- paste0(term[["parameter_stem"]], "_xRE_COEFx[",
                               seq_len(term[["n_groups"]]), ",1]")
    samples <- cbind(samples, values)
  }
  rank <- .selection_sampling_structure(object[["data"]])[["rank"]]
  auxiliary <- matrix(seq(-.2, .2, length.out = S * rank), S, rank)
  colnames(auxiliary) <- paste0("sampling_z[", seq_len(rank), "]")
  samples <- cbind(samples, auxiliary)
  .log_lik_posterior_setup(fit, samples, object[["data"]], object[["priors"]], "estimate", NULL)
}

test_that("dense sampling covariance keeps global dependence across draw chunks", {

  setup <- .pr82_dense_sampling_setup()
  plan <- .data_selection_execution_plan(setup[["data"]])
  expect_identical(plan[["row_blocks"]], list(seq_len(setup[["K"]])))
  expect_true(all(plan[["sampling_covariance"]] != 0))
  state <- .selection_conditioned_sampling_state(setup)
  normalizer <- .selection_conditioned_sampling_normalizer(setup, state[["baseline_mu"]])
  posterior <- .selection_conditioned_sampling_posterior(setup)
  withr::local_options(RoBMA.known_v_covariance_max_bytes =
                        4 * .known_v_covariance_peak_bytes(2L, setup[["K"]]))
  expect_identical(.selection_joint_block_loglik_from_setup(setup), state[["block_log_lik"]])
  expect_identical(.selection_conditioned_sampling_normalizer(setup, state[["baseline_mu"]]), normalizer)
  expect_identical(.selection_conditioned_sampling_posterior(setup), posterior)
})

test_that("global conditioned covariances stream without changing scores or source draws", {

  setup <- .pr82_global_conditioned_setup()
  state <- .selection_conditioned_sampling_state(setup)
  expect_equal(max(lengths(state[["total_covariance"]][["blocks"]])), setup[["K"]])
  expected_sources <- .selection_random_source_posterior(setup, state)
  expected_means <- .selection_random_source_conditional_means(setup, state, expected_sources)
  expected_estimate <- .selection_conditioned_sampling_estimate_targets(setup, c("log_density", "mean"))
  units <- list(seq_len(setup[["K"]]))
  ordinary_setup <- setup
  ordinary_setup[["priors"]][["outcome"]][["bias"]] <- BayesTools::prior_weightfunction(
    "one-sided", .5, BayesTools::wf_fixed(c(1, 1)))
  expected_cluster <- .selection_conditioned_sampling_deletion_loglik(ordinary_setup, units)
  expected_normalizer <- .selection_conditioned_sampling_normalizer(setup, state[["baseline_mu"]])
  withr::local_seed(82)
  seed <- .Random.seed
  for (chunk_size in c(1L, 2L, 3L)) {
    withr::local_options(RoBMA.known_v_covariance_max_bytes =
                          4 * .known_v_covariance_peak_bytes(chunk_size, setup[["K"]]))
    expect_lte(max(lengths(.selection_conditioned_sampling_chunks(setup[["S"]], setup[["K"]]))), chunk_size)
    expect_identical(.selection_joint_block_loglik_from_setup(setup), state[["block_log_lik"]])
    sources <- .selection_conditioned_sampling_posterior(setup)
    expect_equal(sources[["e"]], state[["e"]], tolerance = 0)
    expect_equal(sources[["sources"]], expected_sources, tolerance = 0)
    expect_equal(sources[["source_means"]], expected_means, tolerance = 0)
    expect_identical(.selection_conditioned_sampling_estimate_targets(setup, c("log_density", "mean")),
                     expected_estimate)
    expect_identical(.selection_conditioned_sampling_deletion_loglik(ordinary_setup, units), expected_cluster)
    expect_identical(.selection_conditioned_sampling_normalizer(setup, state[["baseline_mu"]]),
                     expected_normalizer)
  }
  expect_identical(.Random.seed, seed)
})

test_that("conditioned setup slices broadcast scalars and retain cached-factor row identity", {

  setup <- .pr82_global_conditioned_setup()
  setup[["tau_between"]] <- setup[["tau_between"]][1L, , drop = FALSE]
  setup[["selection_conditioned_factors"]] <- .selection_conditioned_sampling_factors(setup)
  setup[["selection_random_factor_samples"]] <- setup[["selection_conditioned_factors"]]
  rows <- c(5L, 1L, 3L)
  sliced <- .selection_conditioned_sampling_subset_setup(setup, rows)
  expect_identical(sliced[["posterior_samples"]], setup[["posterior_samples"]][rows, , drop = FALSE])
  expect_identical(sliced[["tau_between"]], setup[["tau_between"]][rep(1L, 3L), , drop = FALSE])
  expect_identical(sliced[["selection_conditioned_factors"]][["diagonal"]],
                   setup[["selection_conditioned_factors"]][["diagonal"]][rows, , drop = FALSE])
  expect_identical(sliced[["selection_random_factor_samples"]][["loadings"]],
                   lapply(setup[["selection_random_factor_samples"]][["loadings"]], function(x) x[rows, , , drop = FALSE]))
})

test_that("streaming preserves product-space branch indicators at chunk boundaries", {

  setup <- .pr82_global_conditioned_setup()
  alternative <- BayesTools::prior_weightfunction("one-sided", .5,
    BayesTools::wf_fixed(c(1, .8)), model = BayesTools::selection_model(
      other_random_effects = "condition", known_sampling_variance = "condition"))
  setup[["priors"]][["outcome"]][["bias"]] <- BayesTools::prior_mixture(
    list(setup[["priors"]][["outcome"]][["bias"]], alternative))
  indicators <- c(1, 2, 2, 1, 2)
  setup[["posterior_samples"]] <- cbind(setup[["posterior_samples"]], bias_indicator = indicators)
  context <- .selection_conditioned_sampling_context(setup)
  expect_equal(as.numeric(context[["omega"]][, 2L]), ifelse(indicators == 1, .4, .8), tolerance = 0)
  reference <- .selection_conditioned_sampling_state(setup)[["block_log_lik"]]
  withr::local_options(RoBMA.known_v_covariance_max_bytes =
                        4 * .known_v_covariance_peak_bytes(2L, setup[["K"]]))
  expect_identical(.selection_joint_block_loglik_from_setup(setup), reference)
})

test_that("streamed normalizer diagnostics retain every row and reject a later chunk", {

  setup <- .pr82_global_conditioned_setup()
  withr::local_options(RoBMA.known_v_covariance_max_bytes =
                        4 * .known_v_covariance_peak_bytes(2L, setup[["K"]]))
  plan <- .data_selection_execution_plan(setup[["data"]])
  testthat::local_mocked_bindings(
    .selection_conditioned_sampling_normalizer_chunk = function(setup, means) {
      diagnostics <- list(relative_mcse = means[, 1L] / 10000,
                           relative_change = means[, 1L] / 20000)
      .selection_conditioned_sampling_diagnostics(diagnostics, plan)
      list(log_mass = means[, 1L],
           diagnostics = rep(list(diagnostics), length(plan[["row_blocks"]])))
    }, .package = "RoBMA"
  )
  means <- matrix(seq_len(setup[["S"]]), setup[["S"]], setup[["K"]])
  result <- .selection_conditioned_sampling_normalizer(setup, means)
  for (diagnostics in result[["diagnostics"]]) {
    expect_identical(diagnostics[["relative_mcse"]], seq_len(setup[["S"]]) / 10000)
    expect_identical(diagnostics[["relative_change"]], seq_len(setup[["S"]]) / 20000)
  }
  means[setup[["S"]], 1L] <- 10000
  expect_error(.selection_conditioned_sampling_normalizer(setup, means),
               "relative Monte Carlo standard error was 1", fixed = TRUE)
})

test_that("global-block public score computations allocate only bounded covariance chunks", {

  skip_if_not(capabilities("profmem"))
  setup <- .pr82_dense_sampling_setup(S = 80L, K = 20L)
  withr::local_options(RoBMA.known_v_covariance_max_bytes =
                        4 * .known_v_covariance_peak_bytes(3L, setup[["K"]]))
  # Warm metadata/bytecode paths before measuring the computation's allocations.
  .selection_joint_block_loglik_from_setup(.selection_conditioned_sampling_subset_setup(setup, 1:3))
  allocation_file <- tempfile()
  on.exit(unlink(allocation_file), add = TRUE)
  utils::Rprofmem(allocation_file, threshold = 1000)
  on.exit(utils::Rprofmem(NULL), add = TRUE)
  result <- .selection_joint_block_loglik_from_setup(setup)
  utils::Rprofmem(NULL)
  allocations <- readLines(allocation_file, warn = FALSE)
  allocations <- allocations[grepl("^[0-9]+ :", allocations)]
  sizes <- as.numeric(sub(" .*", "", allocations))
  expect_true(length(sizes) > 0L)
  expect_lt(max(sizes), .known_v_covariance_bytes(setup[["S"]], setup[["K"]]))
  expect_equal(dim(result), c(setup[["S"]], length(.data_selection_execution_plan(setup[["data"]])[["row_blocks"]])))
  expect_true(all(is.finite(result)))
  expect_error(.selection_conditioned_sampling_state(setup),
               "Conditioned sampling state would require approximately", fixed = TRUE)
})
