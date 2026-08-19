# ============================================================================ #
# test-02-iwmde-fast-paths.R
# ============================================================================ #

context("IWMDE fast-path equivalence")

source(testthat::test_path("common-functions.R"))
source(testthat::test_path("helper-iwmde.R"))


test_that("IWMDE predictor cache keys remain compact", {

  key <- .iwmde_predictor_cache_key(
    prefix     = "test",
    active_key = paste(rep("branch", 100), collapse = "|"),
    rows       = seq_len(5000),
    unit       = "estimate"
  )
  cache <- new.env(parent = emptyenv())

  expect_lt(nchar(key), 64)
  assign(key, TRUE, envir = cache)
  expect_true(get(key, envir = cache, inherits = FALSE))
})


test_that("IWMDE numeric cache keys preserve represented coordinates", {

  values <- c(
    0,
    .Machine$double.xmin,
    1e-20,
    2e-20,
    1,
    1 + .Machine$double.eps
  )
  keys <- .iwmde_key_number(values)

  expect_length(unique(keys), length(values))
  expect_identical(.iwmde_key_number(-0), .iwmde_key_number(0))
  expect_identical(
    .iwmde_key_number(c(0, 1, -1)),
    c("0000000000000000", "3ff0000000000000", "bff0000000000000")
  )
})


test_that("IWMDE extracts only statically continuous stored columns directly", {

  samples <- cbind(
    theta           = c(10, 20),
    theta_indicator = c(1, 2)
  )
  context <- .iwmde_context_ensure_caches(list(
    posterior_samples = samples,
    flat_prior_list   = list(
      theta = BayesTools::prior("normal", list(mean = 0, sd = 1))
    ),
    indicator_names = "theta_indicator"
  ))

  testthat::local_mocked_bindings(
    .iwmde_parameter_value_row = function(...) {
      stop("row dispatch should not be used", call. = FALSE)
    },
    .package = "RoBMA"
  )
  expect_equal(
    .iwmde_parameter_column_values(context, samples, "theta"),
    samples[, "theta"]
  )

  context[["flat_prior_list"]][["theta"]] <- BayesTools::prior_mixture(list(
    BayesTools::prior("point", list(location = 0)),
    BayesTools::prior("normal", list(mean = 0, sd = 1))
  ))
  expect_error(
    .iwmde_parameter_column_values(context, samples, "theta"),
    "row dispatch should not be used"
  )
})


test_that("IWMDE evaluates focal support once per active key", {

  samples <- cbind(
    theta           = seq_len(6L),
    theta_indicator = c(1, 1, 2, 2, 1, 2)
  )
  context <- list(
    posterior_samples = samples,
    indicator_names   = "theta_indicator"
  )
  calls <- character()
  testthat::local_mocked_bindings(
    .iwmde_focal_prior_cached = function(context, parameter, row,
                                          active_key) {

      calls <<- c(calls, active_key)
      BayesTools::prior("normal", list(mean = 0, sd = 1))
    },
    .package = "RoBMA"
  )

  supports <- .iwmde_parameter_row_supports(
    context   = context,
    parameter = "theta",
    rows      = c(6L, 1L, 4L, 2L)
  )

  expect_length(calls, 2L)
  expect_setequal(calls, c("theta_indicator=1", "theta_indicator=2"))
  expect_equal(supports, matrix(c(-Inf, Inf), nrow = 4L, ncol = 2L,
                                byrow = TRUE))
})


test_that("IWMDE linear supports vectorize the exact row-wise intersection", {

  samples <- cbind(
    theta = c(-.8, -.1, .4, 1.2),
    phi   = c(.3, 1.1, 1.8, .6)
  )
  context <- .iwmde_context_ensure_caches(list(
    posterior_samples = samples,
    flat_prior_list = list(
      theta = BayesTools::prior(
        "normal",
        list(mean = 0, sd = 1),
        truncation = list(-1, 2)
      ),
      phi = BayesTools::prior(
        "normal",
        list(mean = 0, sd = 1),
        truncation = list(0, 2.5)
      )
    ),
    indicator_names = character()
  ))
  weights <- c(theta = 2, phi = -1)

  expected <- matrix(NA_real_, nrow = nrow(samples), ncol = 2L)
  for (row_i in seq_len(nrow(samples))) {
    row <- samples[row_i, ]
    current <- sum(row[names(weights)] * weights)
    coefficients <- weights / sum(weights^2)
    lower_bound <- -Inf
    upper_bound <- Inf
    for (parameter in names(weights)) {
      support <- .iwmde_prior_support(
        context[["flat_prior_list"]][[parameter]]
      )
      coefficient <- coefficients[[parameter]]
      if (coefficient > 0) {
        lower <- current + (support[1L] - row[[parameter]]) / coefficient
        upper <- current + (support[2L] - row[[parameter]]) / coefficient
      } else {
        lower <- current + (support[2L] - row[[parameter]]) / coefficient
        upper <- current + (support[1L] - row[[parameter]]) / coefficient
      }
      lower_bound <- max(lower_bound, lower)
      upper_bound <- min(upper_bound, upper)
    }
    expected[row_i, ] <- c(lower_bound, upper_bound)
  }

  expect_equal(
    .iwmde_linear_row_supports(context, seq_len(nrow(samples)), weights),
    expected,
    tolerance = 0
  )
})


test_that("IWMDE focal-prior deltas preserve row-wise arithmetic", {

  prior  <- BayesTools::prior("normal", list(mean = .2, sd = 1.3))
  values <- c(-.5, 0, .75)
  states <- lapply(seq_len(4L), function(i) {
    list(
      use_focal_prior_delta    = TRUE,
      focal_prior              = prior,
      baseline_log_prior       = i / 7,
      baseline_focal_log_prior = -i / 11
    )
  })
  expected <- unlist(lapply(states, function(state) {
    state[["baseline_log_prior"]] +
      BayesTools::lpdf(prior, values) -
      state[["baseline_focal_log_prior"]]
  }), use.names = FALSE)

  observed <- .iwmde_predictor_log_prior(
    context     = list(),
    parameter   = "theta",
    values      = values,
    row_states  = states,
    replacement = list(type = "scalar")
  )

  expect_identical(observed, expected)
})


test_that("known-V structured location changes match Cholesky reference", {

  set.seed(27)
  S        <- 5L
  yi       <- c(.2, -.1, .4, .3, -.2)
  mu       <- matrix(stats::rnorm(S * length(yi)), nrow = S)
  mu_basis <- matrix(stats::rnorm(S * length(yi)), nrow = S)
  context  <- list(predictor_cache = new.env(parent = emptyenv()))

  diagonal_blocks <- lapply(seq_along(yi), function(i) {
    list(index = i, covariance = matrix(.05 + i / 100, 1L, 1L))
  })
  diagonal_extra <- matrix(
    stats::runif(S * length(yi), .01, .10),
    nrow = S
  )
  diagonal <- .iwmde_normal_location_likelihood_change_known_v_structured(
    context        = context,
    yi             = yi,
    mu             = mu,
    mu_basis       = mu_basis,
    block_data     = diagonal_blocks,
    extra_variance = diagonal_extra
  )
  diagonal_reference <-
    .iwmde_normal_location_likelihood_change_known_v_cholesky(
      yi             = yi,
      mu             = mu,
      mu_basis       = mu_basis,
      block_data     = diagonal_blocks,
      extra_variance = diagonal_extra
    )
  expect_equal(diagonal, diagonal_reference, tolerance = 1e-12)

  full_blocks <- list(
    list(index = 1:3, covariance = matrix(
      c(.20, .04, .02, .04, .18, .03, .02, .03, .16),
      nrow = 3L
    )),
    list(index = 4:5, covariance = matrix(
      c(.12, .03, .03, .14),
      nrow = 2L
    ))
  )
  shifts <- cbind(
    matrix(seq(.01, .05, length.out = S), nrow = S, ncol = 3L),
    matrix(seq(.02, .06, length.out = S), nrow = S, ncol = 2L)
  )
  full <- .iwmde_normal_location_likelihood_change_known_v_structured(
    context        = context,
    yi             = yi,
    mu             = mu,
    mu_basis       = mu_basis,
    block_data     = full_blocks,
    extra_variance = shifts
  )
  full_reference <- .iwmde_normal_location_likelihood_change_known_v_cholesky(
    yi             = yi,
    mu             = mu,
    mu_basis       = mu_basis,
    block_data     = full_blocks,
    extra_variance = shifts
  )
  expect_equal(full, full_reference, tolerance = 1e-12)

  noncommon <- shifts
  noncommon[, 2L] <- noncommon[, 2L] + .01
  expect_null(.iwmde_normal_location_likelihood_change_known_v_structured(
    context        = context,
    yi             = yi,
    mu             = mu,
    mu_basis       = mu_basis,
    block_data     = full_blocks,
    extra_variance = noncommon
  ))
})


test_that("known-V common-shift likelihood matches the joint reference", {

  set.seed(28)
  S  <- 5L
  yi <- c(.2, -.1, .4, .3, -.2)
  V  <- matrix(0, nrow = length(yi), ncol = length(yi))
  V[1:3, 1:3] <- matrix(
    c(.20, .04, .02, .04, .18, .03, .02, .03, .16),
    nrow = 3L
  )
  V[4:5, 4:5] <- matrix(c(.12, .03, .03, .14), nrow = 2L)
  data <- data.frame(yi = yi)
  attr(data, "known_V")      <- TRUE
  attr(data, "known_V_data") <- .known_v_prepare(
    V                         = V,
    keep_rows                 = rep(TRUE, length(yi)),
    known_v_parameterization  = "block_mvn",
    known_v_residual_fraction = NULL
  )

  shifts <- cbind(
    matrix(seq(.01, .05, length.out = S), nrow = S, ncol = 3L),
    matrix(seq(.02, .06, length.out = S), nrow = S, ncol = 2L)
  )
  setup <- list(
    outcome_type       = "norm",
    is_weightfunction  = FALSE,
    weights            = NULL,
    data               = data,
    posterior_samples  = matrix(numeric(), nrow = S, ncol = 0L),
    K                  = length(yi),
    S                  = S,
    yi                 = yi,
    mu                 = matrix(stats::rnorm(S * length(yi)), nrow = S),
    tau_within         = sqrt(shifts),
    effect_direction   = "positive"
  )
  context <- list(
    data            = data,
    predictor_cache = new.env(parent = emptyenv())
  )

  expect_equal(
    .iwmde_log_lik_known_v_joint_sum_common_shift(context, setup),
    .log_lik_known_v_joint_sum_from_setup(setup),
    tolerance = 1e-12
  )

  invariant_shifts <- matrix(
    c(.03, .03, .03, .04, .04),
    nrow = S,
    ncol = length(yi),
    byrow = TRUE
  )
  setup[["tau_within"]] <- sqrt(invariant_shifts)
  expect_equal(
    .iwmde_log_lik_known_v_joint_sum_common_shift(context, setup),
    .log_lik_known_v_joint_sum_from_setup(setup),
    tolerance = 1e-12
  )

  setup[["tau_within"]][, 2L] <- setup[["tau_within"]][, 2L] + .01
  expect_null(.iwmde_log_lik_known_v_joint_sum_common_shift(context, setup))
})


test_that("IWMDE adaptive-grid scaling preserves valid numeric states", {

  tiny <- c(
    .Machine$double.xmin * .Machine$double.eps,
    2 * .Machine$double.xmin * .Machine$double.eps
  )

  expect_equal(.iwmde_rescale_positive(numeric()), numeric())
  expect_equal(.iwmde_rescale_positive(c(0, 0)), c(0, 0))
  expect_equal(.iwmde_rescale_positive(tiny), c(.5, 1))
  expect_error(.iwmde_rescale_positive(c(1, NA_real_)), "must be finite")
  expect_error(.iwmde_rescale_positive(c(1, Inf)), "must be finite")
  expect_error(.iwmde_rescale_positive(c(1, -1e-300)), "must be nonnegative")
})


test_that("IWMDE identifies sampled random SD focal parameters", {

  dat <- data.frame(
    yi    = c(.10, .20, .30),
    x     = c(0, 1, 2),
    study = c("s1", "s2", "s3")
  )
  object <- brma.mv(
    yi                         = yi,
    V                          = diag(rep(.04, 3L)),
    random                     = ~ 1 | study,
    mods                       = ~ x,
    scale                      = ~ x,
    data                       = dat,
    measure                    = "GEN",
    prior_unit_information_sd  = 1,
    marginalize_estimate_level = FALSE,
    only_priors                = TRUE
  )
  prior <- BayesTools::prior("normal", list(mean = 0, sd = 1))
  context <- list(
    object          = object,
    data            = object[["data"]],
    flat_prior_list = list(
      tau               = prior,
      log_tau_intercept = prior,
      log_tau_x         = prior,
      mu                = prior
    )
  )

  expect_true(.iwmde_parameter_controls_sampled_random_sd(
    context,
    "log_tau_intercept"
  ))
  expect_true(.iwmde_parameter_controls_sampled_random_sd(context, "log_tau_x"))
  expect_true(.iwmde_parameter_controls_sampled_random_sd(context, "tau[1]"))
  expect_false(.iwmde_parameter_controls_sampled_random_sd(context, "mu"))
  expect_equal(
    .iwmde_predictor_column_basis(
      context = context,
      column  = "mu_intercept",
      state   = list(active_setup = list())
    )[["mu_basis"]],
    rep(1, 3L)
  )
  context_no_mods <- context
  attr(context_no_mods[["data"]], "mods") <- FALSE
  expect_equal(
    .iwmde_predictor_column_basis(
      context = context_no_mods,
      column  = "mu_intercept",
      state   = list(active_setup = list())
    )[["mu_basis"]],
    rep(1, 3L)
  )
  expect_null(.iwmde_predictor_column_basis(
    context = context,
    column  = "tau",
    state   = list(active_setup = list())
  ))
})


test_that("IWMDE density aggregation matches row-wise reference", {

  log_terms <- matrix(c(
    0, -1, -Inf, -Inf,
    -Inf, -Inf, -Inf, -Inf,
    2, 1, 0, -1
  ), nrow = 3, byrow = TRUE)
  active_mass <- .75
  denominator <- ncol(log_terms)

  fast <- .iwmde_density_aggregate(
    log_terms   = log_terms,
    active_mass = active_mass,
    denominator = denominator
  )

  y                <- numeric(nrow(log_terms))
  finite_terms     <- integer(nrow(log_terms))
  max_log_ratio    <- rep(Inf, nrow(log_terms))
  ess              <- numeric(nrow(log_terms))
  max_weight_share <- rep(1, nrow(log_terms))
  contributions    <- matrix(0, nrow = nrow(log_terms), ncol = ncol(log_terms))

  for (g in seq_len(nrow(log_terms))) {
    finite <- is.finite(log_terms[g, ])
    finite_terms[g] <- sum(finite)
    if (any(finite)) {
      max_term            <- max(log_terms[g, finite])
      scaled_terms        <- exp(log_terms[g, finite] - max_term)
      sum_scaled_terms    <- sum(scaled_terms)
      y[g]                <- active_mass * exp(max_term) *
        sum_scaled_terms / denominator
      contributions[g, finite] <- active_mass * exp(log_terms[g, finite])
      max_log_ratio[g]    <- max_term - stats::median(log_terms[g, finite])
      ess[g]              <- sum_scaled_terms^2 / sum(scaled_terms^2)
      max_weight_share[g] <- max(scaled_terms) / sum_scaled_terms
    }
  }

  expect_equal(fast[["y"]], y)
  expect_equal(fast[["finite_terms"]], finite_terms)
  expect_equal(fast[["max_log_ratio"]], max_log_ratio)
  expect_equal(fast[["ess"]], ess)
  expect_equal(fast[["max_weight_share"]], max_weight_share)
  expect_equal(as.numeric(fast[["contributions"]]), as.numeric(contributions))
  expect_equal(fast[["sampling_mcse"]], rep(0, nrow(log_terms)))

  padded <- .iwmde_density_aggregate(
    log_terms   = matrix(log(2), nrow = 1L),
    active_mass = 1,
    denominator = 2L
  )
  expect_equal(padded[["y"]], 1)
  expect_equal(rowMeans(padded[["contributions"]]), padded[["y"]])

  expect_error(
    .iwmde_density_aggregate(
      log_terms   = matrix(c(0, Inf), nrow = 1L),
      active_mass = 1,
      denominator = 2L
    ),
    "positive-infinite or undefined"
  )
  expect_error(
    .iwmde_density_aggregate(
      log_terms   = matrix(c(0, NaN), nrow = 1L),
      active_mass = 1,
      denominator = 2L
    ),
    "positive-infinite or undefined"
  )
})


test_that("IWMDE row sampling uncertainty is separate and compact", {

  population_rows   <- seq_len(100000L)
  contribution_rows <- seq(100L, 10000L, length.out = 100L)
  contribution_rows <- unique(as.integer(round(contribution_rows)))
  values            <- seq_len(length(contribution_rows))
  chain_id          <- rep(1:2, length.out = length(contribution_rows))

  density <- .iwmde_density_aggregate(
    log_terms         = matrix(log(values), nrow = 1L),
    active_mass       = 1,
    denominator       = length(contribution_rows),
    contribution_rows = contribution_rows,
    sampling_population_rows = population_rows,
    chain_id         = chain_id
  )
  contributions <- density[["contributions"]]
  expected_sampling_mcse <- sqrt(
    (1 - length(contribution_rows) / length(population_rows)) *
      stats::var(values) / length(contribution_rows)
  )

  expect_equal(as.numeric(contributions), as.numeric(values))
  expect_equal(rowMeans(contributions), density[["y"]])
  expect_equal(density[["sampling_mcse"]], expected_sampling_mcse)
  expect_equal(
    density[["sampling_relative_mcse"]],
    expected_sampling_mcse / mean(values)
  )
  expect_equal(density[["sampling_fraction"]], .001)
  expect_lt(as.numeric(object.size(contributions)), 10000)

  census <- .iwmde_density_aggregate(
    log_terms         = matrix(0, nrow = 1L, ncol = 20L),
    active_mass       = .4,
    denominator       = 20L,
    contribution_rows = seq_len(20L),
    sampling_population_rows = seq_len(20L),
    chain_id         = rep(1L, 20L)
  )
  expect_equal(census[["y"]], .4)
  expect_equal(census[["sampling_mcse"]], 0)
  expect_equal(census[["sampling_fraction"]], 1)

  one_row_census <- .iwmde_density_aggregate(
    log_terms         = matrix(log(2), nrow = 1L),
    active_mass       = .4,
    denominator       = 1L,
    contribution_rows = 7L,
    sampling_population_rows = 7L,
    chain_id          = 1L
  )
  expect_equal(one_row_census[["y"]], .8)
  expect_equal(one_row_census[["sampling_mcse"]], 0)
  expect_equal(one_row_census[["sampling_relative_mcse"]], 0)

  for (m in c(20L, 80L)) {
    constant <- .iwmde_density_aggregate(
      log_terms         = matrix(log(3), nrow = 1L, ncol = m),
      active_mass       = .4,
      denominator       = m,
      contribution_rows = seq_len(m),
      sampling_population_rows = seq_len(100L),
      chain_id          = rep(1L, m)
    )
    expect_equal(constant[["y"]], 1.2)
    expect_equal(constant[["sampling_mcse"]], 0)
  }

  finite_population <- c(1, 3, 6, 10)
  srs_estimates <- apply(
    utils::combn(seq_along(finite_population), 2L),
    2L,
    function(sample_rows) {
      sampled <- .iwmde_density_aggregate(
        log_terms         = matrix(
          log(finite_population[sample_rows]),
          nrow = 1L
        ),
        active_mass       = 1,
        denominator       = length(sample_rows),
        contribution_rows = sample_rows,
        sampling_population_rows = seq_along(finite_population),
        chain_id          = rep(1L, length(sample_rows))
      )

      sampled[["y"]]
    }
  )
  expect_equal(mean(srs_estimates), mean(finite_population))

  expect_error(
    .iwmde_density_aggregate(
      log_terms         = matrix(0, nrow = 1L, ncol = 20L),
      active_mass       = 1,
      denominator       = 100L,
      contribution_rows = seq_len(20L),
      sampling_population_rows = seq_len(100L),
      chain_id          = rep(1L, 20L)
    ),
    "selected-row denominator"
  )
  expect_error(
    .iwmde_density_aggregate(
      log_terms         = matrix(0, nrow = 1L, ncol = 20L),
      active_mass       = 1,
      denominator       = 10L,
      contribution_rows = seq_len(20L),
      sampling_population_rows = seq_len(100L),
      chain_id          = rep(1L, 20L)
    ),
    "selected-row denominator"
  )

  alternating <- matrix(rep(c(0, 1), 50L), nrow = 1L)
  alternating_mcse <- .iwmde_batch_mcse(alternating)
  expect_equal(alternating_mcse[["mcse"]], 0)
  expect_equal(alternating_mcse[["ess"]], 50)
  expect_equal(
    alternating_mcse[["uncertainty_scope"]],
    "selected_continuous_rows_only"
  )
  expect_equal(alternating_mcse[["uncertainty_status"]], "available")
  expect_null(alternating_mcse[["uncertainty_reason"]])

  missed_chain <- .iwmde_density_aggregate(
    log_terms         = matrix(log(seq_len(20L)), nrow = 1L),
    active_mass       = 1,
    denominator       = 20L,
    contribution_rows = seq_len(20L),
    sampling_population_rows = seq_len(100L),
    chain_id          = rep(1L, 20L),
    expected_chain_ids = c(1L, 2L)
  )
  missed_chain_mcse <- .iwmde_batch_mcse(missed_chain[["contributions"]])
  expect_true(all(is.na(missed_chain_mcse[["mcse"]])))
  expect_true(all(is.na(missed_chain_mcse[["relative_mcse"]])))
  expect_true(all(is.na(missed_chain_mcse[["ess"]])))
  expect_equal(
    missed_chain_mcse[["uncertainty_scope"]],
    "unavailable_missing_selected_chain"
  )
  expect_equal(missed_chain_mcse[["uncertainty_status"]], "unavailable")
  expect_match(missed_chain_mcse[["uncertainty_reason"]], "chain\\(s\\): 2")
  expect_true(is.finite(missed_chain[["y"]]))
  expect_true(is.finite(missed_chain[["sampling_mcse"]]))

  sparse_chain <- matrix(seq_len(101L), nrow = 1L)
  attr(sparse_chain, "chain_id") <- c(1L, rep(2L, 100L))
  attr(sparse_chain, "expected_chain_ids") <- c(1L, 2L)
  attr(sparse_chain, "target") <- rowMeans(sparse_chain)
  sparse_chain_mcse <- .iwmde_batch_mcse(sparse_chain)

  expect_true(all(is.na(sparse_chain_mcse[["mcse"]])))
  expect_true(all(is.na(sparse_chain_mcse[["relative_mcse"]])))
  expect_true(all(is.na(sparse_chain_mcse[["ess"]])))
  expect_equal(
    sparse_chain_mcse[["uncertainty_scope"]],
    "unavailable_insufficient_chain_batches"
  )
  expect_equal(sparse_chain_mcse[["uncertainty_status"]], "unavailable")
  expect_match(
    sparse_chain_mcse[["uncertainty_reason"]],
    "fewer than two complete batches.*chain\\(s\\): 1"
  )
  expect_equal(sparse_chain_mcse[["batch_size"]], 2L)
  expect_equal(sparse_chain_mcse[["n_batches"]], 50L)
})


test_that("mixed estimates batch the full conditioned chain sequence", {

  set.seed(20260806)
  conditioned_rows <- seq_len(4000L)
  conditioned_chain_id <- rep(1:4, each = 1000L)
  active_rows <- unlist(lapply(0:3, function(chain) {
    chain * 1000L + sort(sample.int(1000L, 50L))
  }))

  density <- .iwmde_density_aggregate(
    log_terms                = matrix(0, nrow = 1L, ncol = 200L),
    active_mass              = .05,
    denominator              = 200L,
    contribution_rows        = active_rows,
    sampling_population_rows = active_rows,
    chain_id                 = conditioned_chain_id[active_rows],
    expected_chain_ids       = 1:4,
    conditioned_rows         = conditioned_rows,
    conditioned_chain_id     = conditioned_chain_id
  )
  conditional_mcse <- .iwmde_batch_mcse(density[["contributions"]])
  mixture_contributions <- density[["mcmc_contributions"]]
  mixture_mcse <- .iwmde_batch_mcse(mixture_contributions)

  expect_equal(density[["y"]], .05)
  expect_equal(rowMeans(mixture_contributions), density[["y"]])
  expect_equal(
    as.numeric(mixture_contributions[, active_rows]),
    rep(1, 200L)
  )
  expect_equal(
    as.numeric(mixture_contributions[, -active_rows]),
    rep(0, 3800L)
  )
  expect_equal(attr(mixture_contributions, "chain_id"), conditioned_chain_id)
  expect_equal(conditional_mcse[["mcse"]], 0)
  expect_equal(
    mixture_mcse[["relative_mcse"]],
    sqrt(.95 / 200),
    tolerance = .005
  )
  expect_equal(
    density[["active_mass_error"]][["relative_mcse"]],
    mixture_mcse[["relative_mcse"]]
  )
  expect_equal(mixture_mcse[["uncertainty_scope"]], "full_conditioned_rows")
  expect_equal(mixture_mcse[["uncertainty_status"]], "available")

  partial_rows <- active_rows[c(
    seq_len(25L),
    50L + seq_len(25L),
    100L + seq_len(25L),
    150L + seq_len(25L)
  )]
  partial_log_terms <- log(rep(c(.5, 1, 1.5, 2), 25L))
  partial <- .iwmde_density_aggregate(
    log_terms                = matrix(partial_log_terms, nrow = 1L),
    active_mass              = .05,
    denominator              = 100L,
    contribution_rows        = partial_rows,
    sampling_population_rows = active_rows,
    chain_id                 = conditioned_chain_id[partial_rows],
    expected_chain_ids       = 1:4,
    conditioned_rows         = conditioned_rows,
    conditioned_chain_id     = conditioned_chain_id
  )
  partial_mcse <- .iwmde_mixture_mcse(
    contributions      = partial[["contributions"]],
    mcmc_contributions = partial[["mcmc_contributions"]],
    active_mass_error  = partial[["active_mass_error"]],
    active_mass        = .05
  )

  expect_true(all(is.finite(partial_mcse[["mcse"]])))
  expect_true(all(is.finite(partial_mcse[["relative_mcse"]])))
  expect_true(all(is.finite(partial_mcse[["ess"]])))
  expect_gt(partial_mcse[["active_branch_mcse"]], 0)
  expect_gt(partial_mcse[["active_mass_component_mcse"]], 0)
  expect_equal(
    partial_mcse[["mcse"]],
    partial_mcse[["active_branch_mcse"]] +
      partial_mcse[["active_mass_component_mcse"]]
  )
  expect_equal(
    partial_mcse[["uncertainty_scope"]],
    "selected_active_rows_with_mass_bound"
  )
  expect_equal(partial_mcse[["uncertainty_status"]], "available")
  expect_null(partial_mcse[["uncertainty_reason"]])
  expect_equal(
    partial_mcse[["mixture_mcse_type"]],
    "worst_correlation_delta_upper_bound"
  )
})


test_that("IWMDE vectorizes a common focal-prior delta exactly", {

  prior  <- BayesTools::prior("normal", list(mean = 0, sd = 1))
  values <- c(-0.4, 0.1, 0.7)
  baseline_values <- c(-0.2, 0.5)
  baseline_log_prior <- c(-3.1, -2.4)
  baseline_focal <- BayesTools::lpdf(prior, baseline_values)
  row_states <- lapply(seq_along(baseline_values), function(i) {
    list(
      prior_list               = list(theta = prior),
      focal_prior              = prior,
      baseline_log_prior       = baseline_log_prior[[i]],
      baseline_focal_log_prior = baseline_focal[[i]],
      use_focal_prior_delta    = TRUE
    )
  })
  candidates <- list(
    state_index = rep(seq_along(row_states), each = length(values)),
    grid_index  = rep(seq_along(values), times = length(row_states))
  )
  valid_positions <- seq_along(candidates[["state_index"]])
  actual <- .iwmde_replacement_log_prior(
    parameter       = "theta",
    values          = values,
    valid_samples   = matrix(0, nrow = length(valid_positions), ncol = 0L),
    valid_positions = valid_positions,
    candidates      = candidates,
    row_states      = row_states,
    replacement     = list(type = "scalar")
  )
  expected <- vapply(valid_positions, function(position) {
    state <- candidates[["state_index"]][[position]]
    grid  <- candidates[["grid_index"]][[position]]
    baseline_log_prior[[state]] +
      BayesTools::lpdf(prior, values[[grid]]) - baseline_focal[[state]]
  }, numeric(1))

  expect_equal(actual, expected, tolerance = 0)
})


test_that("IWMDE scalar latent random effects match formula reconstruction", {

  dat <- data.frame(
    yi    = c(.10, .20, .30, .40),
    study = c("a", "a", "b", "c")
  )
  object <- brma.mv(
    yi                         = yi,
    V                          = diag(rep(.04, 4L)),
    random                     = ~ 1 | study,
    data                       = dat,
    measure                    = "GEN",
    prior_unit_information_sd  = 1,
    marginalize_estimate_level = FALSE,
    only_priors                = TRUE
  )
  term <- .fitted_formula_design(
    object,
    "mu",
    required = TRUE
  )[["random_effects"]][[1L]]
  z_names <- paste0(
    term[["parameter_stem"]],
    "_xRE_Zx[",
    seq_len(term[["n_groups"]]),
    ",1]"
  )
  samples <- matrix(
    c(
      0.0, .5, 1, 2, 3,
      0.1, 1.0, 4, 5, 6
    ),
    nrow = 2L,
    byrow = TRUE,
    dimnames = list(
      NULL,
      c("mu_intercept", term[["sd_parameter_names"]], z_names)
    )
  )
  priors <- object[["priors"]]
  attr(priors[["location"]][["intercept"]], "parameter") <- "mu"
  context <- list(
    object          = object,
    data            = object[["data"]],
    predictor_cache = new.env(parent = emptyenv())
  )

  formula_priors <- .repair_formula_prior_list(
    prior_list = priors[["location"]],
    parameter  = "mu"
  )
  direct_mu <- .evaluate.brma.mu(
    fit               = object[["fit"]],
    outcome_data      = object[["data"]][["outcome"]],
    mods_data         = object[["data"]][["mods"]],
    mods_formula      = NULL,
    mods_priors       = priors[["location"]],
    priors            = priors,
    is_mods           = FALSE,
    is_PET            = FALSE,
    is_PEESE          = FALSE,
    effect_direction  = "positive",
    bias_adjusted     = TRUE,
    K                 = nrow(dat),
    posterior_samples = samples
  )
  formula_mu <- t(BayesTools::JAGS_evaluate_formula(
    fit            = .posterior_formula_fit(
      fit               = object[["fit"]],
      posterior_samples = samples,
      formula_design    = FALSE
    ),
    formula        = stats::as.formula("~ 1"),
    parameter      = "mu",
    data           = data.frame(row.names = seq_len(nrow(dat))),
    prior_list     = formula_priors,
    formula_target = "fixed"
  ))
  expect_identical(direct_mu, formula_mu)

  fast <- .iwmde_conditioned_random_effects_from_latent(
    context           = context,
    posterior_samples = samples,
    unit              = "estimate"
  )
  reference <- .evaluate.brma.random_effects(
    fit               = object[["fit"]],
    data              = object[["data"]],
    priors            = priors,
    posterior_samples = samples,
    same_data         = TRUE,
    required          = TRUE,
    object            = object
  )
  expect_equal(fast, reference, tolerance = 1e-12)

  fast_setup <- .log_lik_posterior_setup(
    fit                        = object[["fit"]],
    posterior_samples          = samples,
    data                       = object[["data"]],
    priors                     = priors,
    unit                       = "estimate",
    data_hash                  = NULL,
    conditioned_random_effects = fast
  )
  reference_setup <- .log_lik_posterior_setup(
    fit               = object[["fit"]],
    posterior_samples = samples,
    data              = object[["data"]],
    priors            = priors,
    unit              = "estimate",
    data_hash         = NULL
  )
  expect_equal(fast_setup[["mu"]], reference_setup[["mu"]], tolerance = 1e-12)
  expect_equal(
    .log_lik_known_v_joint_sum_from_setup(fast_setup),
    .log_lik_known_v_joint_sum_from_setup(reference_setup),
    tolerance = 1e-12
  )

  coefficient_name <- paste0(term[["parameter_stem"]], "_xRE_COEFx[1,1]")
  coefficient <- matrix(
    0,
    nrow     = nrow(samples),
    ncol     = 1L,
    dimnames = list(NULL, coefficient_name)
  )
  samples_with_coefficient <- cbind(samples, coefficient)
  expect_null(.iwmde_conditioned_random_effects_from_latent(
    context           = context,
    posterior_samples = samples_with_coefficient,
    unit              = "estimate"
  ))
})


test_that("partial mixed qCMDE curves receive conservative diagnostics", {

  grid <- seq(0, 1, length.out = 21L)
  conditioned_rows <- seq_len(80L)
  conditioned_chain_id <- rep(1:2, each = 40L)
  active_rows <- c(seq(1L, 39L, by = 2L), seq(41L, 79L, by = 2L))
  estimator_rows <- active_rows[c(seq_len(10L), 20L + seq_len(10L))]

  testthat::local_mocked_bindings(
    .iwmde_log_q_grid = function(context, parameter, values, row_states,
                                 replacement) {
      matrix(0, nrow = length(values), ncol = length(row_states))
    },
    .package = "RoBMA"
  )

  density <- .iwmde_density_grid(
    context            = list(),
    parameter          = "mu",
    display_grid       = grid,
    normalization_grid = list(
      x            = grid,
      z            = grid,
      log_jacobian = rep(0, length(grid))
    ),
    transform          = .iwmde_parameter_transform(c(0, 1)),
    row_states         = rep(list(list()), 20L),
    active_mass        = .5,
    replacement        = list(type = "scalar"),
    estimator_rows     = estimator_rows,
    population_rows    = active_rows,
    chain_id           = conditioned_chain_id[estimator_rows],
    expected_chain_ids = 1:2,
    conditioned_rows   = conditioned_rows,
    conditioned_chain_id = conditioned_chain_id
  )
  density <- .iwmde_new_density_result(
    fields    = density,
    estimator = "q_grid_cmde"
  )

  continuous_density <- .iwmde_density_grid(
    context            = list(),
    parameter          = "mu",
    display_grid       = grid,
    normalization_grid = list(
      x            = grid,
      z            = grid,
      log_jacobian = rep(0, length(grid))
    ),
    transform          = .iwmde_parameter_transform(c(0, 1)),
    row_states         = rep(list(list()), 20L),
    active_mass        = 1,
    replacement        = list(type = "scalar"),
    estimator_rows     = estimator_rows,
    population_rows    = active_rows,
    chain_id           = conditioned_chain_id[estimator_rows],
    expected_chain_ids = 1:2
  )
  # At an atom location the continuous ordinate remains alpha * f_c; the
  # complementary atom probability is not added to density height.
  expect_equal(density[["y"]], .5 * continuous_density[["y"]])
  expect_true(all(is.finite(density[["mcse"]])))
  expect_true(all(is.finite(density[["ess"]])))
  expect_equal(
    density[["mcse"]],
    density[["active_branch_mcse"]] +
      density[["active_mass_component_mcse"]]
  )
  expect_equal(
    density[["mcmc_uncertainty_scope"]],
    "selected_active_rows_with_mass_bound"
  )
  expect_equal(density[["mcmc_uncertainty_status"]], "available")
  expect_equal(
    density[["mixture_mcse_type"]],
    "worst_correlation_delta_upper_bound"
  )
  expect_equal(density[["sampling_fraction"]], .5)
})


test_that("mixed density estimators receive the full conditioned chain", {

  captured <- NULL
  conditioned_rows <- seq_len(4000L)
  conditioned_chain_id <- rep(1:4, each = 1000L)

  testthat::local_mocked_bindings(
    .iwmde_density_grid = function(...) {
      captured <<- list(...)
      list(y = 1)
    },
    .iwmde_new_density_result = function(fields, estimator) fields,
    .package = "RoBMA"
  )

  plan <- list(
    method      = "q_grid_cmde",
    rows        = list(
      point_mass_total = .95,
      active_mass      = .05
    ),
    grids       = list(
      display_grid       = 0,
      normalization_grid = list(x = c(0, 1))
    ),
    target      = list(parameter = "mu"),
    support     = list(transform = NULL),
    replacement = NULL
  )
  execution <- list(
    row_states               = list(list()),
    active_rows              = 1L,
    sampling_population_rows = 1L,
    population_rows          = conditioned_rows,
    chain_id                 = 1L,
    expected_chain_ids       = 1:4,
    conditioned_chain_id     = conditioned_chain_id,
    n_candidate_rows         = 1L
  )

  .iwmde_execute_plan_estimator(
    context   = list(),
    plan      = plan,
    output    = "density",
    execution = execution
  )

  expect_equal(captured[["conditioned_rows"]], conditioned_rows)
  expect_equal(
    captured[["conditioned_chain_id"]],
    conditioned_chain_id
  )
})


test_that("IWMDE rejects invalid Chen log weights", {

  expect_error(
    .iwmde_density_iwmde(
      context          = list(),
      parameter        = "mu",
      display_grid     = 0,
      row_states       = list(list(), list()),
      active_rows      = 1:2,
      active_values    = c(0, 1),
      proposal_weight  = list(log_weight = c(0, Inf), method = "test"),
      active_mass       = 1,
      replacement       = NULL,
      n_candidate_rows = 2L
    ),
    "target 'mu'.*posterior row 2.*proposal-density construction"
  )
})


test_that("qCMDE density normalization is row-wise and mass preserving", {

  grid <- seq(0, 1, length.out = 101)
  normalization_grid <- list(
    x            = grid,
    z            = grid,
    log_jacobian = rep(0, length(grid))
  )
  row_states <- list(list(), list())

  testthat::local_mocked_bindings(
    .iwmde_log_q_grid = function(context, parameter, values, row_states,
                                 replacement) {

      cbind(
        rep(log(2), length(values)),
        rep(log(8), length(values))
      )
    },
    .package = "RoBMA"
  )

  density <- .iwmde_density_grid(
    context            = list(),
    parameter          = "mu",
    display_grid       = grid,
    normalization_grid = normalization_grid,
    transform          = .iwmde_parameter_transform(c(-Inf, Inf)),
    row_states         = row_states,
    active_mass        = .6,
    replacement        = list(type = "scalar")
  )

  expect_equal(density[["estimator"]], "q_grid_cmde")
  expect_equal(density[["normalization_mass_ratio"]], 1)
  final_width <- diff(density[["normalization_range"]])
  pilot_width <- diff(density[["normalization_initial_range"]])
  expect_equal(density[["pilot_normalization_integral"]],
               .6 * pilot_width / final_width,
               tolerance = 1e-10)
  expect_equal(density[["final_normalization_integral"]], .6,
               tolerance = 1e-10)
  expect_equal(density[["normalization_relative_error"]], 0,
               tolerance = 1e-10)
  expect_equal(.iwmde_trapz(density[["x"]], density[["y"]]), .6 / final_width,
               tolerance = 1e-10)
  expect_equal(density[["y"]], rep(.6 / final_width, length(grid)),
               tolerance = 1e-10)
  expect_equal(density[["pilot_y"]], rep(.6 / pilot_width, length(grid)),
               tolerance = 1e-10)
  expect_equal(density[["n_evaluated_rows"]], 2L)
})


test_that("qCMDE density matches an analytic row-normalized normal mixture", {

  display_grid <- seq(-3, 3, length.out = 101)
  norm_grid    <- seq(-6, 6, length.out = 601)
  means        <- c(-.75, 1.25)
  sds          <- c(.5, 1.1)
  active_mass  <- .8
  normalization_grid <- list(
    x            = norm_grid,
    z            = norm_grid,
    log_jacobian = rep(0, length(norm_grid))
  )

  testthat::local_mocked_bindings(
    .iwmde_log_q_grid = function(context, parameter, values, row_states,
                                 replacement) {

      cbind(
        stats::dnorm(values, means[[1L]], sds[[1L]], log = TRUE),
        stats::dnorm(values, means[[2L]], sds[[2L]], log = TRUE)
      )
    },
    .package = "RoBMA"
  )

  density <- .iwmde_density_grid(
    context            = list(),
    parameter          = "mu",
    display_grid       = display_grid,
    normalization_grid = normalization_grid,
    transform          = .iwmde_parameter_transform(c(-Inf, Inf)),
    row_states         = list(list(), list()),
    active_mass        = active_mass,
    replacement        = list(type = "scalar")
  )
  normalizer <- stats::pnorm(
    density[["normalization_range"]][[2L]],
    means,
    sds
  ) - stats::pnorm(
    density[["normalization_range"]][[1L]],
    means,
    sds
  )
  expected <- active_mass / 2 * (
    stats::dnorm(display_grid, means[[1L]], sds[[1L]]) / normalizer[[1L]] +
      stats::dnorm(display_grid, means[[2L]], sds[[2L]]) / normalizer[[2L]]
  )

  expect_equal(density[["y"]], expected, tolerance = 2e-4)
  expect_equal(.iwmde_trapz(density[["x"]], density[["y"]]),
               .iwmde_trapz(display_grid, expected),
               tolerance = 1e-6)
})


test_that("qCMDE and IWMDE match conjugate normal-normal posterior oracle", {

  yi <- c(-.25, .10, .35, .60)
  sei <- c(.30, .45, .40, .50)
  prior_mean <- .15
  prior_sd   <- 1.20

  posterior_var <- 1 / (1 / prior_sd^2 + sum(1 / sei^2))
  posterior_sd  <- sqrt(posterior_var)
  posterior_mean <- posterior_var * (
    prior_mean / prior_sd^2 + sum(yi / sei^2)
  )
  log_q <- function(mu) {
    vapply(mu, function(value) {
      sum(stats::dnorm(yi, mean = value, sd = sei, log = TRUE)) +
        stats::dnorm(value, mean = prior_mean, sd = prior_sd, log = TRUE)
    }, numeric(1))
  }

  display_grid <- seq(
    posterior_mean - 3 * posterior_sd,
    posterior_mean + 3 * posterior_sd,
    length.out = 91
  )
  norm_grid <- seq(
    posterior_mean - 8 * posterior_sd,
    posterior_mean + 8 * posterior_sd,
    length.out = 501
  )
  normalization_grid <- list(
    x            = norm_grid,
    z            = norm_grid,
    log_jacobian = rep(0, length(norm_grid))
  )
  active_values <- stats::qnorm(
    stats::ppoints(80),
    mean = posterior_mean,
    sd   = posterior_sd
  )
  row_states <- lapply(active_values, function(value) {
    list(baseline_log_q = log_q(value))
  })

  testthat::local_mocked_bindings(
    .iwmde_log_q_grid = function(context, parameter, values, row_states,
                                 replacement) {

      matrix(
        log_q(values),
        nrow = length(values),
        ncol = length(row_states)
      )
    },
    .package = "RoBMA"
  )

  qcmde <- .iwmde_density_grid(
    context            = list(),
    parameter          = "mu",
    display_grid       = display_grid,
    normalization_grid = normalization_grid,
    transform          = .iwmde_parameter_transform(c(-Inf, Inf)),
    row_states         = row_states,
    active_mass        = 1,
    replacement        = list(type = "scalar")
  )
  qcmde_mass <- diff(stats::pnorm(
    qcmde[["normalization_range"]],
    mean = posterior_mean,
    sd   = posterior_sd
  ))
  expect_equal(
    qcmde[["y"]],
    stats::dnorm(display_grid, posterior_mean, posterior_sd) / qcmde_mass,
    tolerance = 1e-3
  )

  iwmde <- .iwmde_density_iwmde(
    context            = list(),
    parameter          = "mu",
    display_grid       = display_grid,
    row_states         = row_states,
    active_rows        = seq_along(active_values),
    active_values      = active_values,
    proposal_weight    = list(
      log_weight = stats::dnorm(
        active_values,
        mean = posterior_mean,
        sd   = posterior_sd,
        log  = TRUE
      ),
      method = "oracle_posterior"
    ),
    active_mass        = 1,
    replacement        = list(type = "scalar"),
    normalization_grid = normalization_grid
  )
  expect_equal(
    iwmde[["y"]],
    stats::dnorm(display_grid, posterior_mean, posterior_sd),
    tolerance = 1e-10
  )
  expect_equal(iwmde[["weight_method"]], "oracle_posterior")
})


test_that("qCMDE and IWMDE match correlated linear-contrast normal oracle", {

  mean_beta <- c(mu_intercept = .20, mu_ablat = -.35)
  sigma <- matrix(
    c(.50^2, .50 * .30 * .65,
      .50 * .30 * .65, .30^2),
    nrow = 2,
    dimnames = list(names(mean_beta), names(mean_beta))
  )
  weights <- c(mu_intercept = 1, mu_ablat = -.5)
  target_mean <- sum(weights * mean_beta[names(weights)])
  target_sd   <- sqrt(as.numeric(t(weights) %*% sigma[names(weights),
                                                        names(weights)] %*%
                                   weights))
  log_q <- function(value) {
    stats::dnorm(value, mean = target_mean, sd = target_sd, log = TRUE)
  }

  display_grid <- seq(target_mean - 3 * target_sd,
                      target_mean + 3 * target_sd,
                      length.out = 81)
  norm_grid <- seq(target_mean - 8 * target_sd,
                   target_mean + 8 * target_sd,
                   length.out = 501)
  normalization_grid <- list(
    x            = norm_grid,
    z            = norm_grid,
    log_jacobian = rep(0, length(norm_grid))
  )
  active_values <- stats::qnorm(
    stats::ppoints(80),
    mean = target_mean,
    sd   = target_sd
  )
  row_states <- lapply(active_values, function(value) {
    list(baseline_log_q = log_q(value))
  })
  spec <- list(type = "linear", weights = weights)

  testthat::local_mocked_bindings(
    .iwmde_log_q_grid = function(context, parameter, values, row_states,
                                 replacement) {

      matrix(
        log_q(values),
        nrow = length(values),
        ncol = length(row_states)
      )
    },
    .package = "RoBMA"
  )

  qcmde <- .iwmde_density_grid(
    context            = list(),
    parameter          = "mu_intercept - .5 * mu_ablat",
    display_grid       = display_grid,
    normalization_grid = normalization_grid,
    transform          = .iwmde_parameter_transform(c(-Inf, Inf)),
    row_states         = row_states,
    active_mass        = 1,
    replacement        = spec
  )
  qcmde_mass <- diff(stats::pnorm(
    qcmde[["normalization_range"]],
    mean = target_mean,
    sd   = target_sd
  ))
  expect_equal(
    qcmde[["y"]],
    stats::dnorm(display_grid, target_mean, target_sd) / qcmde_mass,
    tolerance = 1e-3
  )

  iwmde <- .iwmde_density_iwmde(
    context            = list(),
    parameter          = "mu_intercept - .5 * mu_ablat",
    display_grid       = display_grid,
    row_states         = row_states,
    active_rows        = seq_along(active_values),
    active_values      = active_values,
    proposal_weight    = list(
      log_weight = stats::dnorm(
        active_values,
        mean = target_mean,
        sd   = target_sd,
        log  = TRUE
      ),
      method = "oracle_linear_contrast"
    ),
    active_mass        = 1,
    replacement        = spec,
    normalization_grid = normalization_grid
  )
  expect_equal(
    iwmde[["y"]],
    stats::dnorm(display_grid, target_mean, target_sd),
    tolerance = 1e-10
  )
  expect_equal(iwmde[["weight_method"]], "oracle_linear_contrast")
})


test_that("qCMDE point ordinates are independent of unrelated display values", {

  norm_grid <- seq(-6, 6, length.out = 301)
  normalization_grid <- list(
    x            = norm_grid,
    z            = norm_grid,
    log_jacobian = rep(0, length(norm_grid))
  )

  testthat::local_mocked_bindings(
    .iwmde_log_q_grid = function(context, parameter, values, row_states,
                                 replacement) {

      matrix(
        stats::dnorm(values, log = TRUE),
        nrow = length(values),
        ncol = length(row_states)
      )
    },
    .package = "RoBMA"
  )

  base <- .iwmde_density_grid(
    context            = list(),
    parameter          = "mu",
    display_grid       = 0,
    normalization_grid = normalization_grid,
    transform          = .iwmde_parameter_transform(c(-Inf, Inf)),
    row_states         = rep(list(list()), 30),
    active_mass        = 1,
    replacement        = list(type = "scalar")
  )
  with_far_value <- .iwmde_density_grid(
    context            = list(),
    parameter          = "mu",
    display_grid       = c(0, 100),
    normalization_grid = normalization_grid,
    transform          = .iwmde_parameter_transform(c(-Inf, Inf)),
    row_states         = rep(list(list()), 30),
    active_mass        = 1,
    replacement        = list(type = "scalar")
  )

  expect_equal(with_far_value[["y"]][[1L]], base[["y"]][[1L]],
               tolerance = 1e-12)
  expect_equal(with_far_value[["normalization_range"]],
               base[["normalization_range"]])
})


test_that("qCMDE point ordinates are stable under doubled integration budget", {

  normalizer_grid <- function(n) {
    z <- seq(-6, 6, length.out = n)
    list(
      x            = z,
      z            = z,
      log_jacobian = rep(0, length(z))
    )
  }

  testthat::local_mocked_bindings(
    .iwmde_log_q_grid = function(context, parameter, values, row_states,
                                 replacement) {

      matrix(
        stats::dnorm(values, mean = .2, sd = 1.1, log = TRUE),
        nrow = length(values),
        ncol = length(row_states)
      )
    },
    .package = "RoBMA"
  )

  base <- .iwmde_density_grid(
    context            = list(),
    parameter          = "mu",
    display_grid       = 0,
    normalization_grid = normalizer_grid(151),
    transform          = .iwmde_parameter_transform(c(-Inf, Inf)),
    row_states         = rep(list(list()), 40),
    active_mass        = 1,
    replacement        = list(type = "scalar")
  )
  doubled <- .iwmde_density_grid(
    context            = list(),
    parameter          = "mu",
    display_grid       = 0,
    normalization_grid = normalizer_grid(301),
    transform          = .iwmde_parameter_transform(c(-Inf, Inf)),
    row_states         = rep(list(list()), 40),
    active_mass        = 1,
    replacement        = list(type = "scalar")
  )

  expect_equal(doubled[["y"]][[1L]], base[["y"]][[1L]], tolerance = 1e-4)
  expect_lt(doubled[["ordinate_relative_change"]][[1L]], .01)
})


test_that("qCMDE fails the target when refinement cannot normalize a row", {

  grid <- seq(0, 1, length.out = 101)
  normalization_grid <- list(
    x            = grid,
    z            = grid,
    log_jacobian = rep(0, length(grid))
  )
  row_states <- list(list(), list())

  testthat::local_mocked_bindings(
    .iwmde_log_q_grid = function(context, parameter, values, row_states,
                                 replacement) {

      cbind(
        rep(log(2), length(values)),
        rep(-Inf, length(values))
      )
    },
    .package = "RoBMA"
  )

  failure <- tryCatch(
    .iwmde_density_grid(
      context            = list(),
      parameter          = "mu",
      display_grid       = grid,
      normalization_grid = normalization_grid,
      transform          = .iwmde_parameter_transform(c(-Inf, Inf)),
      row_states         = row_states,
      active_mass        = .6,
      replacement        = list(type = "scalar")
    ),
    iwmde_construction_error = function(e) e
  )
  expect_s3_class(failure, "iwmde_construction_error")
  expect_equal(failure[["target"]], "mu")
  expect_equal(failure[["posterior_rows"]], 2L)
  expect_equal(failure[["stage"]], "conditional-density normalization")
  expect_match(failure[["detail"]], "no finite positive normalizer")
})


test_that("qCMDE retains certified zero-density ordinates", {

  normalization_values <- seq(-3, 3, length.out = 101)
  normalization_grid <- list(
    x            = normalization_values,
    z            = normalization_values,
    log_jacobian = rep(0, length(normalization_values))
  )

  testthat::local_mocked_bindings(
    .iwmde_log_q_grid = function(context, parameter, values, row_states,
                                 replacement) {

      out <- matrix(
        stats::dnorm(values, log = TRUE),
        nrow = length(values),
        ncol = length(row_states)
      )
      out[values == 5, ] <- -Inf
      out
    },
    .package = "RoBMA"
  )

  density <- .iwmde_density_grid(
    context            = list(),
    parameter          = "mu",
    display_grid       = 5,
    normalization_grid = normalization_grid,
    transform          = .iwmde_parameter_transform(c(-Inf, Inf)),
    row_states         = rep(list(list()), 20),
    active_mass        = 1,
    replacement        = list(type = "scalar")
  )

  expect_identical(density[["y"]], 0)
  expect_equal(density[["n_evaluated_rows"]], 20L)
})


test_that("qCMDE rejects malformed joint-density grids atomically", {

  normalization_values <- seq(-1, 1, length.out = 21)
  normalization_grid <- list(
    x            = normalization_values,
    z            = normalization_values,
    log_jacobian = rep(0, length(normalization_values))
  )

  testthat::local_mocked_bindings(
    .iwmde_log_q_grid = function(context, parameter, values, row_states,
                                 replacement) {
      matrix(0, nrow = length(values), ncol = 1L)
    },
    .package = "RoBMA"
  )

  expect_error(
    .iwmde_density_grid(
      context            = list(),
      parameter          = "mu",
      display_grid       = 0,
      normalization_grid = normalization_grid,
      transform          = .iwmde_parameter_transform(c(-Inf, Inf)),
      row_states         = list(list(), list()),
      active_mass        = 1,
      replacement        = list(type = "scalar")
    ),
    "target 'mu'.*joint-density grid evaluation.*invalid matrix"
  )
})


test_that("qCMDE refinement waits for an all-finite certified pair", {

  refinement <- .iwmde_qcmde_select_refinement(
    log_q_display = matrix(0, nrow = 1L, ncol = 2L),
    log_normalizer_sequence = list(
      c(0, 0),
      c(0, -Inf),
      c(0, 0),
      c(0, 0)
    ),
    active_mass = 1,
    denominator = 2L
  )

  expect_equal(refinement[["final_index"]], 3L)
  expect_equal(refinement[["validation_index"]], 4L)
})


test_that("qCMDE refinement grids are nested", {

  z <- seq(-2, 2, length.out = 100L)
  transform <- .iwmde_parameter_transform(c(-Inf, Inf))
  sequence <- .iwmde_qcmde_grid_sequence(
    normalization_grid = .iwmde_qcmde_grid_from_z(z, transform),
    transform          = transform
  )

  expect_length(sequence, 4L)
  for (index in seq_len(length(sequence) - 1L)) {
    expect_true(all(sequence[[index]][["z"]] %in%
                      sequence[[index + 1L]][["z"]]))
  }

  plan <- .iwmde_qcmde_normalizer_plan(
    normalization_grid = .iwmde_qcmde_grid_from_z(z, transform),
    transform          = transform
  )
  expect_identical(plan[["all_grid"]][["z"]],
                   plan[["final_grid"]][["z"]])
  expect_lt(length(plan[["all_grid"]][["z"]]), 400L)
})


test_that("qCMDE evaluates only refinement nodes needed for certification", {

  transform <- .iwmde_parameter_transform(c(-Inf, Inf))
  normalization_grid <- .iwmde_qcmde_grid_from_z(
    z         = seq(-2, 2, length.out = 50L),
    transform = transform
  )
  plan <- .iwmde_qcmde_normalizer_plan(
    normalization_grid = normalization_grid,
    transform          = transform
  )
  display_grid <- 7
  row_states   <- list(list(), list())
  evaluated_values <- list()

  testthat::local_mocked_bindings(
    .iwmde_log_q_grid = function(context, parameter, values, row_states,
                                 replacement) {

      evaluated_values[[length(evaluated_values) + 1L]] <<- values
      matrix(
        stats::dnorm(values, log = TRUE),
        nrow = length(values),
        ncol = length(row_states)
      )
    },
    .package = "RoBMA"
  )

  incremental <- .iwmde_qcmde_evaluate_grid_sequence(
    context         = list(),
    parameter       = "mu",
    display_grid    = display_grid,
    normalizer_plan = plan,
    row_states      = row_states,
    replacement     = list(type = "scalar"),
    estimator_rows  = seq_along(row_states),
    active_mass     = 1,
    denominator     = length(row_states)
  )

  eager_log_q_display <- matrix(
    stats::dnorm(display_grid, log = TRUE),
    nrow = length(display_grid),
    ncol = length(row_states)
  )
  eager_log_q_all <- matrix(
    stats::dnorm(plan[["all_grid"]][["x"]], log = TRUE),
    nrow = length(plan[["all_grid"]][["x"]]),
    ncol = length(row_states)
  )
  eager_log_q_sequence <- lapply(plan[["grid_sequence"]], function(grid) {
    eager_log_q_all[grid[["all_index"]], , drop = FALSE]
  })
  eager_log_normalizer_sequence <- lapply(
    seq_along(eager_log_q_sequence),
    function(index) {
      grid <- plan[["grid_sequence"]][[index]]
      .iwmde_log_trapz_columns(
        x     = grid[["z"]],
        log_y = eager_log_q_sequence[[index]] + grid[["log_jacobian"]]
      )
    }
  )
  eager_refinement <- .iwmde_qcmde_select_refinement(
    log_q_display           = eager_log_q_display,
    log_normalizer_sequence = eager_log_normalizer_sequence,
    active_mass             = 1,
    denominator             = length(row_states)
  )
  incremental_refinement <- .iwmde_qcmde_select_refinement(
    log_q_display           = incremental[["log_q_display"]],
    log_normalizer_sequence = incremental[["log_normalizer_sequence"]],
    active_mass             = 1,
    denominator             = length(row_states)
  )

  grid_sizes <- vapply(
    plan[["grid_sequence"]],
    function(grid) length(grid[["x"]]),
    integer(1)
  )
  expect_identical(grid_sizes, c(50L, 113L, 133L, 175L))
  expect_identical(
    vapply(evaluated_values, length, integer(1)),
    c(length(display_grid) + grid_sizes[[1L]], diff(grid_sizes[1:3]))
  )
  expect_length(incremental[["log_q_sequence"]], 3L)
  expect_equal(
    incremental[["log_q_sequence"]],
    eager_log_q_sequence[seq_len(3L)],
    tolerance = 0
  )
  expect_equal(
    incremental[["log_normalizer_sequence"]],
    eager_log_normalizer_sequence[seq_len(3L)],
    tolerance = 0
  )
  expect_identical(incremental_refinement, eager_refinement)

  unused_final_x <- setdiff(
    plan[["grid_sequence"]][[4L]][["x"]],
    plan[["grid_sequence"]][[3L]][["x"]]
  )
  expect_length(unused_final_x, grid_sizes[[4L]] - grid_sizes[[3L]])
  expect_false(any(unused_final_x %in% unlist(evaluated_values)))

  incremental_final_index <- incremental_refinement[["final_index"]]
  eager_final_index       <- eager_refinement[["final_index"]]
  incremental_normalizer  <- incremental[["log_normalizer_sequence"]]
  incremental_normalizer  <- incremental_normalizer[[incremental_final_index]]
  incremental_y <- .iwmde_qcmde_density_from_normalizer(
    log_q_display  = incremental[["log_q_display"]],
    log_normalizer = incremental_normalizer,
    active_mass    = 1,
    denominator    = length(row_states)
  )
  eager_y <- .iwmde_qcmde_density_from_normalizer(
    log_q_display  = eager_log_q_display,
    log_normalizer = eager_log_normalizer_sequence[[eager_final_index]],
    active_mass    = 1,
    denominator    = length(row_states)
  )
  expect_equal(incremental_y, eager_y, tolerance = 0)
})


test_that("qCMDE exhausts uncertified grids and preserves infinite diagnostics", {

  transform <- .iwmde_parameter_transform(c(-Inf, Inf))
  normalization_grid <- .iwmde_qcmde_grid_from_z(
    z         = seq(-2, 2, length.out = 20L),
    transform = transform
  )
  plan <- .iwmde_qcmde_normalizer_plan(
    normalization_grid = normalization_grid,
    transform          = transform
  )
  evaluated_values <- list()

  testthat::local_mocked_bindings(
    .iwmde_log_q_grid = function(context, parameter, values, row_states,
                                 replacement) {

      evaluated_values[[length(evaluated_values) + 1L]] <<- values
      out <- matrix(
        stats::dnorm(values, log = TRUE),
        nrow = length(values),
        ncol = length(row_states)
      )
      if (length(evaluated_values) == 2L) {
        attr(out, "max_quadrature_relative_change") <- Inf
      }
      out
    },
    .iwmde_qcmde_refinement_pair_converged = function(...) FALSE,
    .package = "RoBMA"
  )

  evaluation <- .iwmde_qcmde_evaluate_grid_sequence(
    context         = list(),
    parameter       = "mu",
    display_grid    = 7,
    normalizer_plan = plan,
    row_states      = list(list()),
    replacement     = list(type = "scalar"),
    estimator_rows  = 1L,
    active_mass     = 1,
    denominator     = 1L
  )

  normalization_values <- unlist(evaluated_values, use.names = FALSE)[-1L]
  expect_length(evaluation[["log_q_sequence"]], length(plan[["grid_sequence"]]))
  expect_equal(
    sort(normalization_values),
    plan[["all_grid"]][["x"]],
    tolerance = 0
  )
  expect_length(unique(normalization_values), length(normalization_values))
  expect_identical(evaluation[["quadrature_change"]], Inf)
})


test_that("qCMDE refinement preserves infinite relative changes", {

  density_call <- 0L
  testthat::local_mocked_bindings(
    .iwmde_qcmde_density_from_normalizer = function(...) {
      density_call <<- density_call + 1L
      switch(
        as.character(density_call),
        "1" = c(0, 100),
        "2" = c(1, 101),
        c(1, 101)
      )
    },
    .package = "RoBMA"
  )

  refinement <- .iwmde_qcmde_select_refinement(
    log_q_display          = matrix(0, nrow = 2L, ncol = 1L),
    log_normalizer_sequence = rep(list(0), 4L),
    active_mass             = 1,
    denominator             = 1L
  )

  change <- .iwmde_qcmde_ordinate_change(
    pilot_y = c(0, 100),
    final_y = c(1, 101)
  )
  expect_identical(.iwmde_max_or_na(change[["relative"]]), Inf)
  expect_equal(refinement[["final_index"]], 3L)
  expect_equal(refinement[["validation_index"]], 4L)
})


test_that("qCMDE fails when the validation normalizer is non-finite", {

  normalization_values <- seq(-1, 1, length.out = 21)
  normalization_grid <- list(
    x            = normalization_values,
    z            = normalization_values,
    log_jacobian = rep(0, length(normalization_values))
  )
  normalizer_call <- 0L

  testthat::local_mocked_bindings(
    .iwmde_log_q_grid = function(context, parameter, values, row_states,
                                 replacement) {
      matrix(0, nrow = length(values), ncol = length(row_states))
    },
    .iwmde_log_trapz_columns = function(x, log_y) {
      normalizer_call <<- normalizer_call + 1L
      if (normalizer_call == 3L) c(0, -Inf) else c(0, 0)
    },
    .iwmde_qcmde_select_refinement = function(...) {
      list(
        pilot_index        = 1L,
        final_index        = 2L,
        validation_index   = 3L,
        n_refinement_steps = 1L
      )
    },
    .package = "RoBMA"
  )

  expect_error(
    .iwmde_density_grid(
      context            = list(),
      parameter          = "mu",
      display_grid       = 0,
      normalization_grid = normalization_grid,
      transform          = .iwmde_parameter_transform(c(-Inf, Inf)),
      row_states         = list(list(), list()),
      active_mass        = 1,
      replacement        = list(type = "scalar")
    ),
    "target 'mu'.*posterior row 2.*normalization validation"
  )
})


test_that("qCMDE uses refined normalizers and diagnoses pilot-grid impact", {

  display_grid <- 0
  norm_grid    <- seq(-.5, .5, length.out = 21)
  row_states   <- rep(list(list()), 20)
  normalization_grid <- list(
    x            = norm_grid,
    z            = norm_grid,
    log_jacobian = rep(0, length(norm_grid))
  )

  testthat::local_mocked_bindings(
    .iwmde_log_q_grid = function(context, parameter, values, row_states,
                                 replacement) {

      matrix(
        stats::dnorm(values, log = TRUE),
        nrow = length(values),
        ncol = length(row_states)
      )
    },
    .package = "RoBMA"
  )

  density <- .iwmde_density_grid(
    context            = list(),
    parameter          = "mu",
    display_grid       = display_grid,
    normalization_grid = normalization_grid,
    transform          = .iwmde_parameter_transform(c(-Inf, Inf)),
    row_states         = row_states,
    active_mass        = 1,
    replacement        = list(type = "scalar")
  )

  expect_lt(density[["pilot_normalization_integral"]], 1)
  expect_equal(density[["final_normalization_integral"]], 1, tolerance = 1e-10)
  expect_gt(density[["ordinate_relative_change"]][[1]], .05)
  expect_gt(density[["max_normalizer_relative_change"]], .05)
  expect_match(
    .iwmde_diagnostics_bf_failure_reason(list(
      estimator              = "q_grid_cmde",
      ordinate               = density[["y"]][[1]],
      relative_mcse          = .1,
      finite_terms           = 20,
      ess                    = 20,
      max_weight_share       = .2,
      evaluation_value       = 0,
      normalization_range    = range(norm_grid),
      normalization_relative_error =
        density[["normalization_relative_error"]],
      ordinate_relative_change =
        density[["ordinate_relative_change"]][[1]],
      max_normalizer_relative_change =
        density[["max_normalizer_relative_change"]]
    )),
    "qCMDE.*ordinate"
  )
})


test_that("qCMDE BF gate rejects validation movement despite perfect mass", {

  diagnostics <- list(
    estimator              = "q_grid_cmde",
    ordinate               = 1,
    relative_mcse          = .01,
    finite_terms           = 80L,
    ess                    = 40,
    max_weight_share       = .10,
    evaluation_value       = 0,
    final_normalization_integral = 1,
    normalization_relative_error = 0,
    active_mass            = 1,
    normalization_mass_ratio = 1,
    ordinate_relative_change = .06,
    max_normalizer_relative_change = 0
  )

  expect_match(
    .iwmde_diagnostics_bf_failure_reason(diagnostics),
    "qCMDE.*ordinate"
  )
})


test_that("IWMDE density reports raw support mass without scaling the curve", {

  grid <- seq(0, 1, length.out = 101)
  normalization_grid <- list(
    x            = grid,
    z            = grid,
    log_jacobian = rep(0, length(grid))
  )
  row_states <- list(
    list(baseline_log_q = 0),
    list(baseline_log_q = 0)
  )

  testthat::local_mocked_bindings(
    .iwmde_log_q_grid = function(context, parameter, values, row_states,
                                 replacement) {

      matrix(0, nrow = length(values), ncol = length(row_states))
    },
    .package = "RoBMA"
  )

  density <- .iwmde_density_iwmde(
    context            = list(),
    parameter          = "mu",
    display_grid       = grid,
    row_states         = row_states,
    active_rows        = 1:2,
    active_values      = c(.25, .75),
    proposal_weight    = list(
      log_weight = rep(log(.5), 2L),
      method     = "mock"
    ),
    active_mass        = 1,
    replacement        = list(type = "scalar"),
    normalization_grid = normalization_grid
  )

  expect_equal(density[["estimator"]], "iwmde")
  expect_equal(density[["weight_method"]], "mock")
  expect_equal(
    density[["support_grid_normalization_integral"]],
    .5,
    tolerance = 1e-10
  )
  expect_equal(density[["normalization_mass_ratio"]], 2, tolerance = 1e-10)
  expect_equal(.iwmde_trapz(density[["x"]], density[["y"]]), .5,
               tolerance = 1e-10)
  expect_equal(density[["y"]], rep(.5, length(grid)), tolerance = 1e-10)
})


test_that("IWMDE fails the target on a non-finite proposal density", {

  grid <- seq(0, 1, length.out = 101)
  normalization_grid <- list(
    x            = grid,
    z            = grid,
    log_jacobian = rep(0, length(grid))
  )
  row_states <- list(
    list(baseline_log_q = 0),
    list(baseline_log_q = 0)
  )

  testthat::local_mocked_bindings(
    .iwmde_log_q_grid = function(context, parameter, values, row_states,
                                 replacement) {

      matrix(0, nrow = length(values), ncol = length(row_states))
    },
    .package = "RoBMA"
  )

  failure <- tryCatch(
    .iwmde_density_iwmde(
      context            = list(),
      parameter          = "mu",
      display_grid       = grid,
      row_states         = row_states,
      active_rows        = 1:2,
      active_values      = c(.25, .75),
      proposal_weight    = list(
        log_weight = c(log(.5), -Inf),
        method     = "mock"
      ),
      active_mass        = 1,
      replacement        = list(type = "scalar"),
      normalization_grid = normalization_grid
    ),
    iwmde_construction_error = function(e) e
  )
  expect_s3_class(failure, "iwmde_construction_error")
  expect_equal(failure[["target"]], "mu")
  expect_equal(failure[["posterior_rows"]], 2L)
  expect_equal(failure[["stage"]], "proposal-density construction")
  expect_match(failure[["detail"]], "zero or non-finite")
})


test_that("restricted Chen fallback weights are moment matched", {

  p <- stats::ppoints(500)

  distance_fit <- stats::qgamma(p, shape = 18, scale = .04)
  active       <- seq(.02, 2, length.out = 100)
  gamma_weight <- .iwmde_chen_gamma_log_weight_single(
    active_values = active,
    weight_values = distance_fit,
    support       = c(0, Inf)
  )

  shape <- mean(distance_fit)^2 / stats::var(distance_fit)
  rate  <- mean(distance_fit) / stats::var(distance_fit)

  expect_equal(gamma_weight[["method"]], "chen_gamma")
  expect_equal(
    gamma_weight[["log_weight"]],
    stats::dgamma(active, shape = shape, rate = rate, log = TRUE),
    tolerance = 1e-12
  )

  prob_fit <- stats::qbeta(p, shape1 = 6, shape2 = 3)
  active   <- seq(.01, .99, length.out = 100)
  beta_weight <- .iwmde_chen_beta_log_weight(
    active_values = active,
    weight_values = prob_fit,
    support       = c(0, 1)
  )

  location <- mean(prob_fit)
  variance <- stats::var(prob_fit)
  common   <- location * (1 - location) / variance - 1

  expect_equal(beta_weight[["method"]], "chen_beta")
  expect_equal(
    beta_weight[["log_weight"]],
    stats::dbeta(
      active,
      shape1 = location * common,
      shape2 = (1 - location) * common,
      log    = TRUE
    ),
    tolerance = 1e-12
  )
})


test_that("Chen weights dispatch by row-specific support", {

  supports <- matrix(
    c(0, 1,
      0, 1,
      1, 2,
      1, 2),
    ncol  = 2,
    byrow = TRUE
  )
  rownames(supports) <- as.character(seq_len(4L))
  calls <- list()

  testthat::local_mocked_bindings(
    .iwmde_parameter_row_supports = function(context, parameter, rows,
                                             parameter_spec) {

      context[["supports"]][as.character(rows), , drop = FALSE]
    },
    .iwmde_chen_logit_conditional_normal_log_weight = function(
        context, parameter, parameter_spec, active_rows, active_values,
        weight_rows, weight_values, support) {

      calls[[length(calls) + 1L]] <<- list(
        active_rows = active_rows,
        weight_rows = weight_rows,
        support     = support
      )

      list(
        log_weight = rep(sum(support), length(active_values)),
        method     = paste0("bounded_", paste(support, collapse = "_"))
      )
    },
    .package = "RoBMA"
  )

  weight <- .iwmde_chen_log_weight(
    context        = list(supports = supports),
    parameter      = "theta",
    parameter_spec = list(type = "primitive"),
    active_rows    = seq_len(4L),
    active_values  = c(.25, .75, 1.25, 1.75),
    weight_rows    = seq_len(4L),
    weight_values  = c(.30, .70, 1.30, 1.70),
    support        = c(0, 2)
  )

  expect_length(calls, 2L)
  expect_equal(calls[[1L]][["active_rows"]], 1:2)
  expect_equal(calls[[1L]][["weight_rows"]], 1:2)
  expect_equal(calls[[1L]][["support"]], c(0, 1))
  expect_equal(calls[[2L]][["active_rows"]], 3:4)
  expect_equal(calls[[2L]][["weight_rows"]], 3:4)
  expect_equal(calls[[2L]][["support"]], c(1, 2))
  expect_equal(weight[["log_weight"]], c(1, 1, 3, 3))
  expect_equal(weight[["method"]], "chen_mixed(bounded_0_1,bounded_1_2)")
  expect_length(weight[["partitions"]], 2L)
  expect_equal(
    vapply(weight[["partitions"]], `[[`, character(1), "method"),
    c("bounded_0_1", "bounded_1_2")
  )
  expect_equal(
    vapply(weight[["partitions"]], `[[`, integer(1), "n_eval_rows"),
    c(2L, 2L)
  )
})


test_that("Chen proposal pooling ignores nuisance product states", {

  supports <- matrix(
    rep(c(0, 1), 4L),
    ncol  = 2,
    byrow = TRUE
  )
  rownames(supports) <- as.character(seq_len(4L))
  samples <- data.frame(
    bias_indicator = c(1, 1, 2, 2),
    `omega[1]`     = c(.4, .4, .4, .4),
    check.names    = FALSE
  )
  calls   <- list()

  testthat::local_mocked_bindings(
    .iwmde_parameter_row_supports = function(context, parameter, rows,
                                             parameter_spec) {

      context[["supports"]][as.character(rows), , drop = FALSE]
    },
    .iwmde_chen_logit_conditional_normal_log_weight = function(
        context, parameter, parameter_spec, active_rows, active_values,
        weight_rows, weight_values, support) {

      calls[[length(calls) + 1L]] <<- list(
        active_rows = active_rows,
        weight_rows = weight_rows,
        support     = support
      )

      list(
        log_weight = rep(active_rows[[1L]], length(active_values)),
        method     = paste0("branch_", active_rows[[1L]])
      )
    },
    .package = "RoBMA"
  )

  weight <- .iwmde_chen_log_weight(
    context        = list(
      supports          = supports,
      posterior_samples = samples,
      indicator_names   = "bias_indicator",
      selection_spec    = list(jags_omega = "omega")
    ),
    parameter      = "theta",
    parameter_spec = list(type = "primitive"),
    active_rows    = seq_len(4L),
    active_values  = c(.25, .75, .25, .75),
    weight_rows    = seq_len(4L),
    weight_values  = c(.30, .70, .30, .70),
    support        = c(0, 1)
  )

  expect_length(calls, 1L)
  expect_equal(calls[[1L]][["active_rows"]], 1:4)
  expect_equal(calls[[1L]][["weight_rows"]], 1:4)
  expect_equal(weight[["log_weight"]], rep(1, 4L))
  expect_equal(weight[["method"]], "branch_1")
  expect_length(weight[["partitions"]], 1L)
})


test_that("Chen weight dispatcher falls back when conditioning fails", {

  supports <- matrix(
    rep(c(-Inf, Inf), 3L),
    ncol  = 2L,
    byrow = TRUE
  )
  rownames(supports) <- as.character(seq_len(3L))

  testthat::local_mocked_bindings(
    .iwmde_parameter_row_supports = function(context, parameter, rows,
                                             parameter_spec) {

      context[["supports"]][as.character(rows), , drop = FALSE]
    },
    .iwmde_chen_conditional_normal_log_weight = function(...) {

      .iwmde_chen_conditional_stop("conditioning failed")
    },
    .iwmde_chen_marginal_normal_log_weight = function(active_values,
                                                      weight_values) {

      list(log_weight = rep(log(.25), length(active_values)),
           method = "chen_marginal_normal")
    },
    .package = "RoBMA"
  )

  weight <- .iwmde_chen_log_weight(
    context        = list(supports = supports),
    parameter      = "mu",
    parameter_spec = list(type = "primitive"),
    active_rows    = seq_len(3L),
    active_values  = c(-1, 0, 1),
    weight_rows    = seq_len(3L),
    weight_values  = c(-1, 0, 1),
    support        = c(-Inf, Inf)
  )

  expect_equal(weight[["method"]], "chen_marginal_normal")
  expect_equal(weight[["log_weight"]], rep(log(.25), 3L))

  supports[,] <- rep(c(0, 1), each = 3L)
  testthat::local_mocked_bindings(
    .iwmde_parameter_row_supports = function(context, parameter, rows,
                                             parameter_spec) {

      context[["supports"]][as.character(rows), , drop = FALSE]
    },
    .iwmde_chen_logit_conditional_normal_log_weight = function(...) {

      .iwmde_chen_conditional_stop("conditioning failed")
    },
    .iwmde_chen_beta_log_weight = function(active_values, weight_values,
                                           support) {

      list(log_weight = rep(log(.5), length(active_values)),
           method = "chen_beta")
    },
    .package = "RoBMA"
  )

  weight <- .iwmde_chen_log_weight(
    context        = list(supports = supports),
    parameter      = "rho",
    parameter_spec = list(type = "primitive"),
    active_rows    = seq_len(3L),
    active_values  = c(.25, .5, .75),
    weight_rows    = seq_len(3L),
    weight_values  = c(.25, .5, .75),
    support        = c(0, 1)
  )

  expect_equal(weight[["method"]], "chen_beta")
  expect_equal(weight[["log_weight"]], rep(log(.5), 3L))
})


test_that("Chen conditional-normal weights match bivariate normal oracle", {

  n <- 500L
  x <- seq(-2, 2, length.out = n)
  z <- .3 + 1.2 * x + stats::qnorm(stats::ppoints(n)) * .4
  active_values <- z + .15 * sin(seq_len(n) / 17)
  samples <- cbind(mu = z, PET = x)
  context <- list(
    posterior_samples = samples,
    indicator_names   = character(),
    selection_spec    = NULL,
    flat_prior_list   = list(mu = TRUE, PET = TRUE)
  )

  weight <- .iwmde_chen_conditional_normal_log_weight(
    context        = context,
    parameter      = "mu",
    parameter_spec = list(type = "primitive"),
    active_rows    = seq_len(n),
    active_values  = active_values,
    weight_rows    = seq_len(n),
    weight_values  = z
  )

  x_center <- mean(x)
  x_scale  <- stats::sd(x)
  x_fit    <- (x - x_center) / x_scale
  cov_mat  <- stats::cov(cbind(z, x_fit))
  beta     <- cov_mat[1L, 2L] / cov_mat[2L, 2L]
  means    <- mean(z) + ((x - x_center) / x_scale) * beta
  variance <- sum((z - means)^2) / (length(z) - 1L)
  expected <- stats::dnorm(
    active_values,
    mean = means,
    sd   = sqrt(variance),
    log  = TRUE
  )

  expect_equal(weight[["method"]], "chen_conditional_normal")
  expect_equal(weight[["log_weight"]], expected, tolerance = 1e-12)
})


test_that("Chen Gaussian weights preserve represented small scales", {

  x      <- 2^-400 * seq(-2, 2, length.out = 200L)
  z      <- 2^-300 * (seq(-1, 1, length.out = 200L) +
    stats::qnorm(stats::ppoints(200L)) / 10)
  active <- z + 2^-305 * sin(seq_along(z))

  gaussian <- .iwmde_chen_conditional_gaussian(
    z_fit  = z,
    x_fit  = matrix(x),
    x_eval = matrix(x)
  )
  expected <- stats::dnorm(
    active / gaussian[["focal_scale"]],
    mean = gaussian[["means"]],
    sd   = gaussian[["sd"]],
    log  = TRUE
  ) - log(gaussian[["focal_scale"]])

  expect_true(all(is.finite(expected)))
  expect_true(gaussian[["sd"]] > 0)
})


test_that("Chen conditioning transforms follow declared prior support", {

  context <- list(
    flat_prior_list = list(
      arbitrary_beta = BayesTools::prior(
        "beta",
        parameters = list(alpha = 2, beta = 3)
      ),
      arbitrary_gamma = BayesTools::prior(
        "gamma",
        parameters = list(shape = 2, rate = 3)
      ),
      arbitrary_cor = BayesTools::prior(
        "uniform",
        parameters = list(a = -1, b = 1)
      ),
      upper = BayesTools::prior(
        "normal",
        parameters = list(mean = 0, sd = 1),
        truncation = list(lower = -Inf, upper = 3)
      ),
      weight = BayesTools::prior(
        "dirichlet",
        parameters = list(alpha = c(1, 1))
      ),
      unrestricted = BayesTools::prior(
        "normal",
        parameters = list(mean = 0, sd = 1)
      )
    ),
    selection_spec = NULL
  )
  unit <- .iwmde_chen_transform_conditioning_column(
    context,
    c(0, 1, -.Machine$double.eps, 1 + .Machine$double.eps),
    c(0, 1),
    "arbitrary_beta"
  )
  nonnegative <- .iwmde_chen_transform_conditioning_column(
    context,
    c(-.Machine$double.xmin, 0, .Machine$double.xmin),
    0,
    "arbitrary_gamma"
  )
  correlation <- .iwmde_chen_transform_conditioning_column(
    context,
    c(-0.5, 0, 0.5),
    0,
    "arbitrary_cor"
  )
  upper <- .iwmde_chen_transform_conditioning_column(
    context,
    c(2, 3, 4),
    3,
    "upper"
  )
  simplex <- .iwmde_chen_transform_conditioning_column(
    context,
    c(0.25, 0.75),
    0.5,
    "weight[1]"
  )
  unrestricted <- .iwmde_chen_transform_conditioning_column(
    context,
    c(-2, 2),
    0,
    "unrestricted"
  )

  expect_identical(unit[["fit"]], c(-Inf, Inf, NA_real_, NA_real_))
  expect_identical(
    nonnegative[["fit"]],
    c(NA_real_, 0, log1p(.Machine$double.xmin))
  )
  expect_equal(correlation[["fit"]], stats::qlogis(c(0.25, 0.5, 0.75)))
  expect_identical(upper[["fit"]], c(-log(2), 0, NA_real_))
  expect_equal(simplex[["fit"]], stats::qlogis(c(0.25, 0.75)))
  expect_identical(unrestricted[["fit"]], c(-2, 2))
})


test_that("logit conditional Chen weights are proper on original scale", {

  n       <- 500L
  x       <- seq(-2, 2, length.out = n)
  z       <- -.2 + .8 * x + stats::qnorm(stats::ppoints(n)) * .15
  rho     <- stats::plogis(z)
  samples <- cbind(rho = rho, mu = x)
  context <- list(
    posterior_samples = samples,
    indicator_names   = character(),
    selection_spec    = NULL,
    flat_prior_list   = list(rho = TRUE, mu = TRUE)
  )

  grid <- seq(.001, .999, length.out = 2000)
  weight <- .iwmde_chen_logit_conditional_normal_log_weight(
    context        = context,
    parameter      = "rho",
    parameter_spec = list(type = "primitive"),
    active_rows    = rep(250L, length(grid)),
    active_values  = grid,
    weight_rows    = seq_len(n),
    weight_values  = rho,
    support        = c(0, 1)
  )

  expect_equal(weight[["method"]], "chen_logit_conditional_normal")
  expect_equal(
    .iwmde_trapz(grid, exp(weight[["log_weight"]])),
    1,
    tolerance = 5e-3
  )
})

test_that("IWMDE omega matrix collapse matches row-wise collapse", {

  omega       <- matrix(c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6), nrow = 3, byrow = TRUE)
  global_cuts <- c(0, .025, .05, .50, 1)
  active_cuts <- c(0, .05, 1)

  fast <- .iwmde_collapse_omega_matrix(
    omega       = omega,
    global_cuts = global_cuts,
    active_cuts = active_cuts
  )
  ref <- t(apply(omega, 1, .iwmde_collapse_omega,
    global_cuts = global_cuts,
    active_cuts = active_cuts
  ))

  expect_equal(fast, ref)
  expect_equal(fast[, 1], c(1, 3, 5))
})


test_that("IWMDE omega collapse rejects unequal merged bins", {

  global_cuts <- c(0, .025, .05, .50, 1)
  active_cuts <- c(0, .05, .50, 1)
  omega       <- c(1, .5, .25, .25)

  collapsed <- .iwmde_collapse_omega(
    omega       = omega,
    global_cuts = global_cuts,
    active_cuts = active_cuts
  )

  expect_true(is.na(collapsed[1]))
  expect_equal(collapsed[2:3], c(.25, .25))
})


test_that("IWMDE structural identities are exact", {

  expect_error(
    .iwmde_indicator_index(1 + .Machine$double.eps, "model_indicator", 2L),
    "integer-valued"
  )
  expect_false(.iwmde_same_p_cuts(
    c(0, .5, 1),
    c(0, .5 + .Machine$double.eps, 1)
  ))

  collapsed <- .iwmde_collapse_omega(
    omega       = c(1, 1 + .Machine$double.eps),
    global_cuts = c(0, .25, .5),
    active_cuts = c(0, .5)
  )
  expect_true(is.na(collapsed[[1L]]))

  tiny <- .Machine$double.eps / 2
  expect_identical(
    .iwmde_linear_weights(c(mu = tiny, tau = 0)),
    c(mu = tiny)
  )

  points <- .iwmde_point_mass_table(
    c(1, 1 + .Machine$double.eps, 1),
    denominator = 3L
  )
  expect_identical(points[["x"]], c(1, 1 + .Machine$double.eps))
  expect_identical(points[["mass"]], c(2 / 3, 1 / 3))
})


test_that("IWMDE omega helpers use selection-specific JAGS names", {

  context <- list(
    selection_spec = list(
      n_bins     = 4L,
      p_cuts     = c(0, .025, .05, .50, 1),
      jags_omega = "custom.omega+beta"
    )
  )
  active_setup <- list(
    is_weightfunction = TRUE,
    selection_spec    = list(
      n_bins     = 2L,
      p_cuts     = c(0, .05, 1),
      jags_omega = "active.omega+beta"
    ),
    priors = list(outcome = list(bias = NULL))
  )
  samples <- matrix(
    c(1, 1, .5, .5, 2, 10,
      2, 2, .4, .4, 2, 11),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(NULL, c(
      paste0("custom.omega+beta[", 1:4, "]"),
      "bias_indicator",
      "mu"
    ))
  )

  out <- .iwmde_likelihood_posterior_samples(
    context      = context,
    samples      = samples,
    active_setup = active_setup
  )
  row_omega <- .iwmde_active_omega(
    context      = context,
    row          = samples[1, ],
    active_setup = list(
      selection_spec = context[["selection_spec"]],
      priors         = active_setup[["priors"]]
    )
  )

  expect_equal(samples[, "bias_indicator"], c(2, 2))
  expect_equal(out[, "bias_indicator"], c(1, 1))
  expect_false(any(grepl("^custom\\.omega\\+beta\\[", colnames(out))))
  expect_true(all(paste0("active.omega+beta[", 1:2, "]") %in% colnames(out)))
  expect_equal(out[, "active.omega+beta[1]"], samples[, "custom.omega+beta[1]"])
  expect_equal(out[, "active.omega+beta[2]"], samples[, "custom.omega+beta[3]"])
  expect_equal(out[, "mu"], samples[, "mu"])
  expect_equal(row_omega, as.numeric(samples[1, paste0("custom.omega+beta[", 1:4, "]")]))
})


test_that("IWMDE parameter filters use selection-specific JAGS names", {

  context <- list(
    selection_spec = list(jags_omega = "custom.omega+beta"),
    posterior_samples = matrix(
      0,
      nrow = 2,
      ncol = 4,
      dimnames = list(NULL, c(
        "mu",
        "custom.omega+beta[1]",
        "custom.omega+beta[2]",
        "eta[1]"
      ))
    ),
    indicator_names    = character(),
    flat_prior_list    = list(),
    focal_prior_cache  = new.env(parent = emptyenv()),
    support_cache      = new.env(parent = emptyenv())
  )

  expect_true(.iwmde_parameter_is_weightfunction_coordinate(
    parameter = "custom.omega+beta[1]",
    context   = context
  ))
  expect_true(.iwmde_parameter_is_weightfunction_coordinate(
    parameter = "eta[1]",
    context   = context
  ))
  expect_false(.iwmde_parameter_is_weightfunction_coordinate(
    parameter = "log_omega[1]",
    context   = context
  ))
  expect_true(.iwmde_chen_is_global_conditioning_column(
    context = context,
    column  = "custom.omega+beta[1]"
  ))
  expect_false(.iwmde_chen_is_global_conditioning_column(
    context = context,
    column  = "eta[1]"
  ))
  expect_false(.iwmde_chen_is_global_conditioning_column(
    context = context,
    column  = "log_omega[1]"
  ))
  expect_equal(
    .iwmde_parameter_spec(context, "custom.omega+beta[1]")[["status"]],
    "unsupported"
  )
})


test_that("Chen conditioning uses independent fitted prior coordinates", {

  columns <- c(
    "mu_intercept",
    "allocation_sd",
    "weight[1]",
    "weight[2]",
    "prior_par_eta_weight[1]",
    "prior_par_eta_weight[2]",
    "derived_random_sd",
    "random_latent[1]"
  )
  context <- list(
    posterior_samples = matrix(
      seq_len(24L),
      nrow     = 3L,
      dimnames = list(NULL, columns)
    ),
    flat_prior_list = list(
      mu_intercept = BayesTools::prior(
        "normal",
        parameters = list(mean = 0, sd = 1)
      ),
      allocation_sd = BayesTools::prior(
        "normal",
        parameters = list(mean = 0, sd = 1)
      ),
      weight = BayesTools::prior(
        "dirichlet",
        parameters = list(alpha = c(1, 1))
      )
    ),
    indicator_names = character(),
    selection_spec  = NULL
  )

  expect_identical(
    .iwmde_chen_global_conditioning_columns(context),
    c("mu_intercept", "allocation_sd", "weight[1]")
  )
})


test_that("IWMDE omega localization renames same-cut active branches", {

  context <- list(
    selection_spec = list(
      n_bins     = 2L,
      p_cuts     = c(0, .05, 1),
      jags_omega = "custom.omega+beta"
    )
  )
  active_setup <- list(
    is_weightfunction = TRUE,
    selection_spec    = list(
      n_bins     = 2L,
      p_cuts     = c(0, .05, 1),
      jags_omega = "active.omega+beta"
    ),
    priors = list(outcome = list(bias = NULL))
  )
  samples <- matrix(
    c(1, .5, 2, 10,
      2, .4, 2, 11),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(NULL, c(
      paste0("custom.omega+beta[", 1:2, "]"),
      "bias_indicator",
      "mu"
    ))
  )

  out <- .iwmde_likelihood_posterior_samples(
    context      = context,
    samples      = samples,
    active_setup = active_setup
  )
  out_again <- .iwmde_likelihood_posterior_samples(
    context      = context,
    samples      = out,
    active_setup = active_setup
  )
  with_active <- cbind(
    samples,
    "active.omega+beta[1]" = c(.9, .8),
    "active.omega+beta[2]" = c(.7, .6)
  )
  active_first <- .iwmde_likelihood_posterior_samples(
    context      = context,
    samples      = with_active,
    active_setup = active_setup
  )

  expect_false(any(grepl("^custom\\.omega\\+beta\\[", colnames(out))))
  expect_true(all(paste0("active.omega+beta[", 1:2, "]") %in% colnames(out)))
  expect_equal(out[, "active.omega+beta[1]"], samples[, "custom.omega+beta[1]"])
  expect_equal(out[, "active.omega+beta[2]"], samples[, "custom.omega+beta[2]"])
  expect_equal(out[, "bias_indicator"], c(1, 1))
  expect_equal(out_again, out)
  expect_equal(active_first[, "active.omega+beta[1]"], c(.9, .8))
  expect_equal(active_first[, "active.omega+beta[2]"], c(.7, .6))
  expect_false(any(grepl("^custom\\.omega\\+beta\\[", colnames(active_first))))
})


test_that("IWMDE active omega does not truncate unexpected omega lengths", {

  context <- list(
    selection_spec = list(
      n_bins     = 4L,
      p_cuts     = c(0, .025, .05, .50, 1),
      jags_omega = "omega"
    )
  )
  active_setup <- list(
    selection_spec = list(
      n_bins     = 2L,
      p_cuts     = c(0, .05, 1),
      jags_omega = "omega"
    ),
    priors = list(outcome = list(bias = NULL))
  )
  row <- c("omega[1]" = 1, "omega[2]" = .5, "omega[3]" = .25)

  omega <- .iwmde_active_omega(
    context      = context,
    row          = row,
    active_setup = active_setup
  )

  expect_true(all(is.na(omega)))
  expect_length(omega, 2L)
})

test_that("IWMDE active omega fails closed when weights are unavailable", {

  context <- list(
    selection_spec = list(
      n_bins     = 2L,
      p_cuts     = c(0, .05, 1),
      jags_omega = "global.omega"
    )
  )
  active_setup <- list(
    selection_spec = list(
      n_bins     = 2L,
      p_cuts     = c(0, .05, 1),
      jags_omega = "active.omega"
    ),
    priors = list(outcome = list(bias = NULL))
  )

  omega <- .iwmde_active_omega(
    context      = context,
    row          = c(mu = 0),
    active_setup = active_setup
  )

  expect_equal(omega, rep(NA_real_, 2L))

  active_setup[["selection_spec"]][["jags_data"]] <- list(
    active.omega = c(1, .5)
  )
  omega <- .iwmde_active_omega(
    context      = context,
    row          = c(mu = 0),
    active_setup = active_setup
  )

  expect_equal(omega, c(1, .5))
})

test_that("IWMDE selection context uses localized active-branch indicators", {

  prior <- BayesTools::prior_weightfunction(
    side    = "one-sided",
    steps   = c(.025, .05),
    weights = BayesTools::wf_fixed(c(1, .7, .35))
  )
  data <- list(
    outcome = data.frame(
      yi  = c(.1, .2),
      sei = c(.1, .2)
    )
  )
  attr(data, "effect_direction") <- "positive"

  context <- list(
    object = list(fit = list()),
    data   = data
  )
  active_setup <- list(
    priors            = list(outcome = list(bias = prior)),
    is_weightfunction = TRUE
  )
  samples <- matrix(
    c(2, 2),
    nrow     = 2,
    dimnames = list(NULL, "bias_indicator")
  )

  selection_context <- .iwmde_selection_context_active_branch(
    context           = context,
    active_setup      = active_setup,
    posterior_samples = samples
  )

  expect_equal(selection_context[["bias_indicator"]], c(1L, 1L))
  expect_false(any(selection_context[["use_normal"]]))
  expect_true(all(selection_context[["kernel_mode"]] != SELKERNEL_NORMAL))
  expect_equal(samples[, "bias_indicator"], c(2, 2))
})


test_that("negative-direction PET and PEESE location fast paths match flipped likelihoods", {

  yi  <- c(.20, -.10, .35)
  sei <- c(.25, .40, .55)
  data <- list(outcome = data.frame(yi = yi, sei = sei))
  attr(data, "outcome_type")      <- "norm"
  attr(data, "effect_direction")  <- "negative"

  mu <- matrix(
    c(.10, -.05, .20,
      -.15, .05, .25),
    nrow  = 2L,
    byrow = TRUE
  )
  tau <- matrix(
    c(.30, .35, .40,
      .32, .38, .44),
    nrow  = 2L,
    byrow = TRUE
  )
  current <- c(.15, -.20)
  values  <- c(-.30, .05, .40)
  active_setup <- list(is_weightfunction = FALSE)
  context <- list(data = data)
  setup <- list(
    yi               = yi,
    mu               = mu,
    tau_within       = tau,
    sei              = sei,
    weights          = NULL,
    effect_direction = "negative"
  )
  explicit_log_lik <- function(row, value, mu_basis) {

    candidate_mu <- mu[row, ] + (value - current[[row]]) * mu_basis
    sum(stats::dnorm(
      -yi,
      mean = -candidate_mu,
      sd   = sqrt(tau[row, ]^2 + sei^2),
      log  = TRUE
    ))
  }

  cases <- list(
    list(parameter = "PET",   mu_basis = -sei),
    list(parameter = "PEESE", mu_basis = -sei^2)
  )
  for (case in cases) {
    mu_basis <- case[["mu_basis"]]
    row_states <- lapply(seq_len(nrow(mu)), function(row) {
      list(
        active_setup     = active_setup,
        baseline_log_lik = explicit_log_lik(row, current[[row]], mu_basis)
      )
    })
    basis <- list(
      formula_mu      = FALSE,
      formula_logtau  = FALSE,
      scale_update    = "none",
      log_tau_basis   = NULL,
      mu_basis        = matrix(mu_basis, nrow = nrow(mu), ncol = length(sei),
                               byrow = TRUE),
      current         = current
    )

    testthat::local_mocked_bindings(
      .iwmde_predictor_log_prior = function(context, parameter, values,
                                            row_states, replacement) {

        rep(0, length(values) * length(row_states))
      },
      .package = "RoBMA"
    )

    fast <- .iwmde_log_q_grid_normal_location_group(
      context     = context,
      parameter   = case[["parameter"]],
      values      = values,
      row_states  = row_states,
      replacement = list(type = "scalar"),
      setup       = setup,
      basis       = basis
    )
    reference <- matrix(
      vapply(seq_along(row_states), function(row) {
        vapply(values, explicit_log_lik, numeric(1), row = row,
               mu_basis = mu_basis)
      }, numeric(length(values))),
      nrow = length(values)
    )

    expect_equal(fast, reference, tolerance = 1e-12)
  }
})


test_that("unmaterialized formula predictors bypass the quadratic location path", {

  data <- list(outcome = data.frame(yi = 0, sei = 1))
  attr(data, "outcome_type") <- "norm"
  row_states <- list(list(active_setup = list(is_weightfunction = FALSE)))

  for (formula_component in c("formula_mu", "formula_logtau")) {
    basis <- list(
      formula_mu     = FALSE,
      formula_logtau = FALSE,
      scale_update   = "none",
      log_tau_basis  = NULL,
      mu_basis       = NULL
    )
    basis[[formula_component]] <- TRUE

    expect_null(.iwmde_log_q_grid_normal_location_group(
      context     = list(data = data),
      parameter   = "mu_x",
      values      = 0,
      row_states  = row_states,
      replacement = list(type = "scalar"),
      setup       = list(),
      basis       = basis
    ))
  }
})


test_that("fixed formula coordinates materialize an exact location basis", {

  mods <- data.frame(
    group = factor(
      c("sensitivity", "specificity", "sensitivity", "specificity"),
      levels = c("sensitivity", "specificity")
    )
  )
  formula_result <- BayesTools::JAGS_formula(
    formula = ~ 0 + group,
    parameter = "mu",
    data = mods,
    prior_list = list(
      group = BayesTools::prior_factor(
        "normal",
        list(mean = 0, sd = 1),
        contrast = "independent"
      )
    )
  )
  columns <- c("mu_group[1]", "mu_group[2]")
  posterior_samples <- matrix(
    c(-0.2, 0.5, 1.1, 0.1, 0.8, 1.2),
    nrow = 2L,
    byrow = TRUE,
    dimnames = list(NULL, c(columns, "irrelevant"))
  )
  fit <- coda::mcmc(posterior_samples)
  class(fit) <- c("BayesTools_fit", class(fit))
  attr(fit, "prior_list") <- c(
    formula_result$prior_list,
    list(irrelevant = BayesTools::prior("normal", list(mean = 0, sd = 1)))
  )
  attr(fit, "formula_design") <- list(mu = formula_result$formula_design)
  fit <- BayesTools:::.bt_attach_parameter_map(fit)
  fit <- BayesTools:::.bt_attach_draw_geometry(fit)
  fit <- BayesTools:::.bt_attach_fit_contract(fit)

  data <- list(
    outcome = data.frame(yi = rep(0, nrow(mods)), sei = rep(1, nrow(mods))),
    mods    = mods
  )
  attr(data, "outcome_type") <- "norm"
  context <- list(
    data        = data,
    formula_fit = fit
  )
  setup <- list(
    posterior_samples = posterior_samples,
    mu = matrix(0, nrow = nrow(posterior_samples), ncol = nrow(mods))
  )
  row_states <- lapply(seq_len(nrow(posterior_samples)), function(i) {
    list(row_index = i, active_setup = list())
  })

  basis <- .iwmde_predictor_update_basis(
    context     = context,
    parameter   = "mu_group[1]",
    row_states  = row_states,
    replacement = list(type = "scalar"),
    setup       = setup
  )

  expected <- matrix(
    formula_result$formula_design$model_matrix[, "group1"],
    nrow = nrow(posterior_samples),
    ncol = nrow(mods),
    byrow = TRUE
  )
  expect_false(isTRUE(basis[["formula_mu"]]))
  expect_true(isTRUE(basis[["formula_mu_affine"]]))
  expect_identical(basis[["current"]], posterior_samples[, "mu_group[1]"])
  expect_identical(unname(basis[["mu_basis"]]), unname(expected))
  expect_false(.iwmde_q_grid_fixed_mu_is_invariant(
    context,
    list(changed_parameters = "mu_group[1]")
  ))
  expect_true(.iwmde_q_grid_fixed_mu_is_invariant(
    context,
    list(changed_parameters = "irrelevant")
  ))
})


test_that("fixed-mu reuse follows persisted indirect formula dependencies", {

  x_prior <- BayesTools::prior("normal", list(mean = 0, sd = 1))
  attr(x_prior, "multiply_by") <- "tau"
  formula_result <- BayesTools::JAGS_formula(
    formula = ~ 0 + x,
    parameter = "mu",
    data = data.frame(x = c(-1, 1)),
    prior_list = list(x = x_prior)
  )
  posterior_samples <- matrix(
    c(0.2, 0.5, 1.1, -0.1, 0.8, 1.2),
    nrow = 2L,
    byrow = TRUE,
    dimnames = list(NULL, c("mu_x", "tau", "irrelevant"))
  )
  fit <- coda::mcmc(posterior_samples)
  class(fit) <- c("BayesTools_fit", class(fit))
  attr(fit, "prior_list") <- c(
    formula_result$prior_list,
    list(
      tau = BayesTools::prior("gamma", list(shape = 2, rate = 1)),
      irrelevant = BayesTools::prior("normal", list(mean = 0, sd = 1))
    )
  )
  attr(fit, "formula_design") <- list(mu = formula_result$formula_design)
  fit <- BayesTools:::.bt_attach_parameter_map(fit)
  fit <- BayesTools:::.bt_attach_draw_geometry(fit)
  fit <- BayesTools:::.bt_attach_fit_contract(fit)
  context <- list(formula_fit = fit)

  expect_false(.iwmde_q_grid_fixed_mu_is_invariant(
    context,
    list(changed_parameters = "tau")
  ))
  expect_true(.iwmde_q_grid_fixed_mu_is_invariant(
    context,
    list(changed_parameters = "irrelevant")
  ))
})


test_that("known-V random predictor fallbacks remain marginal", {

  testthat::local_mocked_bindings(
    .iwmde_log_q_grid_known_v_random_group_iid = function(...) NULL,
    .iwmde_predictor_setup = function(...) list(),
    .iwmde_predictor_update_basis = function(...) {
      list(formula_mu = TRUE)
    },
    .iwmde_log_q_grid_normal_location_group = function(...) NULL,
    .iwmde_predictor_candidates = function(...) {
      stop("conditional predictor fallback was evaluated", call. = FALSE)
    },
    .iwmde_uses_known_v_random_marginal_likelihood = function(...) TRUE,
    .package = "RoBMA"
  )

  expect_null(.iwmde_log_q_grid_predictor_group(
    context     = list(data = list()),
    parameter   = "derived_fixed_target",
    values      = 0,
    row_states  = list(list(active_setup = list())),
    replacement = list(type = "fallback")
  ))
})


test_that("known-V formula random locations use the exact marginal quadratic", {

  S <- 2L
  K <- 4L
  model_matrix <- rbind(
    c(1, 0), c(0, 1), c(1, 0), c(0, 1)
  )
  group_map <- c(1L, 1L, 2L, 2L)
  coefficient_factors <- list(
    matrix(c(.7, .2, 0, .5), nrow = 2L, byrow = TRUE),
    matrix(c(.4, -.1, 0, .8), nrow = 2L, byrow = TRUE)
  )
  factor_plans <- list(list(
    type                  = "group",
    model_matrix          = model_matrix,
    group_map             = group_map,
    coefficient_structure = "dense"
  ))
  factor_states <- lapply(coefficient_factors, function(factor) {
    list(list(coefficient_factor = factor))
  })
  sampling_covariance <- matrix(0, nrow = K, ncol = K)
  sampling_covariance[1:2, 1:2] <- matrix(c(.4, .1, .1, .5), 2L)
  sampling_covariance[3:4, 3:4] <- matrix(c(.6, -.05, -.05, .3), 2L)
  dependency_blocks <- list(1:2, 3:4)
  fixed_mu <- rbind(
    c(-.2, .5, .1, .4),
    c(.3, -.1, .2, .6)
  )
  mu_basis <- rbind(
    c(1, 0, 1, 0),
    c(.8, .1, .9, -.1)
  )
  yi <- c(.1, .7, -.2, .9)
  setup <- list(
    S                 = S,
    K                 = K,
    effect_direction  = "positive",
    posterior_samples = matrix(0, nrow = S, ncol = 0L)
  )
  random_setup <- list(
    blocks              = "study",
    dependency_blocks   = dependency_blocks,
    sampling_covariance = sampling_covariance
  )
  random_factors <- list(
    row_blocks    = dependency_blocks,
    factor_plans  = factor_plans,
    factor_states = factor_states
  )

  testthat::local_mocked_bindings(
    .iwmde_known_v_random_marginal_setup = function(context) random_setup,
    .brma_mv_random_effects_marginal_factor_states = function(...) {
      random_factors
    },
    .package = "RoBMA"
  )

  actual <- .iwmde_normal_location_likelihood_change_known_v_random(
    context  = list(object = list(), data = list()),
    setup    = setup,
    yi       = yi,
    mu       = fixed_mu,
    mu_basis = mu_basis
  )
  expected_linear    <- numeric(S)
  expected_quadratic <- numeric(S)
  for (s in seq_len(S)) {
    coefficient_covariance <- tcrossprod(coefficient_factors[[s]])
    for (index in dependency_blocks) {
      design <- model_matrix[index, , drop = FALSE]
      covariance <- sampling_covariance[index, index, drop = FALSE] +
        design %*% coefficient_covariance %*% t(design)
      residual <- yi[index] - fixed_mu[s, index]
      expected_linear[[s]] <- expected_linear[[s]] + as.numeric(
        crossprod(mu_basis[s, index], solve(covariance, residual))
      )
      expected_quadratic[[s]] <- expected_quadratic[[s]] + as.numeric(
        crossprod(mu_basis[s, index], solve(covariance, mu_basis[s, index]))
      )
    }
  }

  expect_equal(actual[["linear"]], expected_linear, tolerance = 1e-12)
  expect_equal(actual[["quadratic"]], expected_quadratic, tolerance = 1e-12)
})


test_that("known-V location quadratics cover every random factor contract", {

  K <- 4L
  S <- 2L
  y <- c(.1, .7, -.2, .9)
  means <- rbind(
    c(-.2, .5, .1, .4),
    c(.3, -.1, .2, .6)
  )
  bases <- rbind(
    c(1, 0, 1, 0),
    c(.8, .1, .9, -.1)
  )
  sampling_covariance <- matrix(
    c(
      .50, .08, .02, 0,
      .08, .60, 0, .03,
      .02, 0, .45, -.04,
      0, .03, -.04, .55
    ),
    nrow = K,
    byrow = TRUE
  )
  extra_variances <- rbind(
    c(.02, .01, .03, .04),
    c(.01, .03, .02, .01)
  )
  model_matrix <- rbind(
    c(1, 0), c(0, 1), c(1, .5), c(.3, 1)
  )
  group_map <- as.integer(c(1, 1, 2, 2))
  crossed_group_map <- as.integer(c(1, 2, 1, 2))

  coefficient_factor <- function(scale, rho) {

    correlation <- matrix(c(1, rho, rho, 1), nrow = 2L)
    sweep(t(chol(correlation)), 1L, scale, "*")
  }
  markov_state <- function(scale, rho) {

    list(
      coefficient_factor              = coefficient_factor(scale, rho),
      coefficient_scale               = scale,
      markov_transition               = rho,
      markov_innovation_variance       = 1 - rho^2
    )
  }
  plan <- function(type, coefficient_structure = "dense",
                   map = group_map, group_covariance = NULL) {

    out <- list(
      type                  = type,
      model_matrix          = model_matrix,
      group_map             = map,
      coefficient_structure = coefficient_structure
    )
    if (!is.null(group_covariance)) {
      out[["group_covariance"]] <- group_covariance
    }
    out
  }

  dense_roots <- list(
    matrix(c(.15, .02, 0, .1, .03, 0, .12, -.01), nrow = K),
    matrix(c(.08, 0, .02, .14, 0, .04, -.01, .1), nrow = K)
  )
  known_group_covariance <- matrix(c(1, .3, .3, 1.4), nrow = 2L)
  cases <- list(
    dense = list(
      plans = list(list(type = "dense")),
      states = lapply(dense_roots, function(root) {
        list(list(covariance = tcrossprod(root)))
      })
    ),
    id_diag = list(
      plans = list(plan("group", "diagonal")),
      states = list(
        list(list(coefficient_factor = diag(c(.4, .7)))),
        list(list(coefficient_factor = diag(c(.6, .3))))
      )
    ),
    us_cs_hcs = list(
      plans = list(plan("group", "dense")),
      states = list(
        list(list(coefficient_factor = coefficient_factor(c(.5, .7), .25))),
        list(list(coefficient_factor = coefficient_factor(c(.4, .8), -.2)))
      )
    ),
    ar_har_car = list(
      plans = list(plan("group", "markov")),
      states = list(
        list(markov_state(c(.5, .7), .25)),
        list(markov_state(c(.4, .8), -.2))
      )
    ),
    row_scaled = list(
      plans = list(plan("row_group", "dense")),
      states = list(
        list(list(
          coefficient_factor = coefficient_factor(c(.5, .7), .25),
          row_scale = c(.8, 1.1, .9, 1.2)
        )),
        list(list(
          coefficient_factor = coefficient_factor(c(.4, .8), -.2),
          row_scale = c(1.2, .7, 1, .9)
        ))
      )
    ),
    known_group = list(
      plans = list(plan(
        "known_group",
        "dense",
        map = crossed_group_map,
        group_covariance = known_group_covariance
      )),
      states = list(
        list(list(coefficient_factor = coefficient_factor(c(.5, .7), .25))),
        list(list(coefficient_factor = coefficient_factor(c(.4, .8), -.2)))
      )
    )
  )

  covariance_contribution <- function(factor_plan, factor_state) {

    if (identical(factor_plan[["type"]], "dense")) {
      return(factor_state[["covariance"]])
    }
    root <- factor_plan[["model_matrix"]] %*%
      factor_state[["coefficient_factor"]]
    if (identical(factor_plan[["type"]], "row_group")) {
      root <- root * factor_state[["row_scale"]]
    }
    covariance <- tcrossprod(root)
    map <- factor_plan[["group_map"]]
    if (identical(factor_plan[["type"]], "known_group")) {
      covariance <- covariance * factor_plan[["group_covariance"]][
        map,
        map,
        drop = FALSE
      ]
    } else {
      covariance <- covariance * outer(map, map, "==")
    }
    covariance
  }

  for (case_name in names(cases)) {
    case <- cases[[case_name]]
    actual <- .marglik_covariance_plan_location_quadratic_batch(
      cache                    = NULL,
      y                        = y,
      means                    = means,
      bases                    = bases,
      sampling_covariance      = sampling_covariance,
      random_covariance_plans  = case[["plans"]],
      random_covariance_states = case[["states"]],
      block_indices            = list(seq_len(K)),
      extra_variances          = extra_variances
    )
    expected_linear    <- numeric(S)
    expected_quadratic <- numeric(S)
    for (draw in seq_len(S)) {
      covariance <- sampling_covariance + diag(extra_variances[draw, ])
      for (factor_i in seq_along(case[["plans"]])) {
        covariance <- covariance + covariance_contribution(
          case[["plans"]][[factor_i]],
          case[["states"]][[draw]][[factor_i]]
        )
      }
      residual <- y - means[draw, ]
      expected_linear[[draw]] <- as.numeric(crossprod(
        bases[draw, ],
        solve(covariance, residual)
      ))
      expected_quadratic[[draw]] <- as.numeric(crossprod(
        bases[draw, ],
        solve(covariance, bases[draw, ])
      ))
    }

    expect_equal(
      actual[["linear"]],
      expected_linear,
      tolerance = 2e-12,
      info = case_name
    )
    expect_equal(
      actual[["quadratic"]],
      expected_quadratic,
      tolerance = 2e-12,
      info = case_name
    )
  }
})


test_that("metadata-declared non-affine formula locations retain evaluation", {

  formula <- ~ 1
  attr(formula, "log(intercept)") <- TRUE
  formula_result <- BayesTools::JAGS_formula(
    formula   = formula,
    parameter = "mu",
    data      = data.frame(row = 1:2),
    prior_list = list(
      intercept = BayesTools::prior(
        "gamma",
        list(shape = 2, rate = 1)
      )
    )
  )
  posterior_samples <- matrix(
    c(0.5, 1.5),
    ncol     = 1L,
    dimnames = list(NULL, "mu_intercept")
  )
  fit <- coda::mcmc(posterior_samples)
  class(fit) <- c("BayesTools_fit", class(fit))
  attr(fit, "prior_list") <- formula_result$prior_list
  attr(fit, "formula_design") <- list(mu = formula_result$formula_design)
  fit <- BayesTools:::.bt_attach_parameter_map(fit)
  fit <- BayesTools:::.bt_attach_draw_geometry(fit)
  fit <- BayesTools:::.bt_attach_fit_contract(fit)

  data <- list(
    outcome = data.frame(yi = c(0, 0), sei = c(1, 1))
  )
  attr(data, "outcome_type") <- "norm"
  context <- list(data = data, formula_fit = fit)
  setup <- list(
    posterior_samples = posterior_samples,
    mu                = matrix(0, nrow = 2L, ncol = 2L)
  )
  row_states <- lapply(seq_len(2L), function(i) {
    list(row_index = i, active_setup = list())
  })

  basis <- .iwmde_predictor_update_basis(
    context     = context,
    parameter   = "mu_intercept",
    row_states  = row_states,
    replacement = list(type = "fallback"),
    setup       = setup
  )

  expect_true(isTRUE(basis[["formula_mu"]]))
  expect_identical(basis[["formula_mu_status"]], "non_affine")
  expect_null(basis[["mu_basis"]])
})


test_that("normal cluster rho grid preserves boundaries and prior rows", {

  S      <- 2L
  values <- c(-.1, 0, .4, 1, 1.1)
  data <- list(outcome = data.frame(
    yi  = c(.02, .15, -.05, .30),
    sei = c(.10, .12, .08, .20)
  ))
  attr(data, "outcome_type") <- "norm"
  attr(data, "cluster")      <- TRUE
  setup <- list(
    mu = matrix(c(
      -.04, .05, .12, .18,
       .02, .08, .15, .22
    ), nrow = S, byrow = TRUE),
    tau_total = matrix(c(
      .10, .14, .18, .22,
      .12, .16, .20, .24
    ), nrow = S, byrow = TRUE),
    yi                = data[["outcome"]][["yi"]],
    sei               = data[["outcome"]][["sei"]],
    cluster           = list(a = c(1L, 3L), b = c(2L, 4L)),
    is_weightfunction = FALSE,
    weights           = NULL
  )
  row_states <- rep(list(list(active_setup = list())), S)
  basis <- list(
    scale_update   = "rho",
    formula_mu     = FALSE,
    formula_logtau = FALSE,
    mu_basis       = NULL,
    log_tau_basis  = NULL
  )
  prior <- matrix(seq_len(length(values) * S) / 100,
                  nrow = length(values), ncol = S)

  testthat::local_mocked_bindings(
    .iwmde_predictor_log_prior = function(...) as.numeric(prior),
    .package = "RoBMA"
  )
  observed <- .iwmde_log_q_grid_normal_cluster_rho_group(
    context     = list(data = data),
    parameter   = "rho",
    values      = values,
    row_states  = row_states,
    replacement = list(type = "scalar"),
    setup       = setup,
    basis       = basis
  )
  valid <- values >= 0 & values <= 1
  expected <- matrix(-Inf, nrow = length(values), ncol = S)
  expected[valid, ] <- .log_lik_cluster_norm_analytic_rho_grid_sum(
    setup = setup,
    yi    = setup[["yi"]],
    vi    = setup[["sei"]]^2,
    rho   = values[valid]
  )
  expected <- expected + prior

  expect_equal(observed, expected, tolerance = 1e-12)
})


test_that("negative-direction selected-normal normalizer change matches matrix reference", {

  prior <- BayesTools::prior_weightfunction(
    side    = "one-sided",
    steps   = c(.025, .05),
    weights = BayesTools::wf_fixed(c(1, .7, .35))
  )
  yi  <- c(.10, -.20, .30)
  sei <- c(.20, .30, .40)
  data <- list(outcome = data.frame(yi = yi, sei = sei))
  attr(data, "outcome_type")     <- "norm"
  attr(data, "effect_direction") <- "negative"

  S      <- 2L
  G      <- 3L
  values <- c(-.30, .05, .40)
  current <- c(.10, -.15)

  setup <- list(
    mu                = matrix(c(.10, .15, .20, -.05, .00, .05),
                               nrow = S, byrow = TRUE),
    tau_within        = matrix(c(.25, .30, .35, .28, .32, .36),
                               nrow = S, byrow = TRUE),
    sei               = sei,
    weights           = c(1, .5, 2),
    posterior_samples = matrix(
      0,
      nrow     = S,
      ncol     = 1L,
      dimnames = list(NULL, "mu")
    )
  )
  active_setup <- list(
    priors            = list(outcome = list(bias = prior)),
    is_weightfunction = TRUE
  )
  context <- list(
    object          = list(fit = list()),
    data            = data,
    predictor_cache = new.env(parent = emptyenv())
  )
  selection_context <- .iwmde_selection_context_active_branch(
    context           = context,
    active_setup      = active_setup,
    posterior_samples = setup[["posterior_samples"]]
  )
  basis <- list(
    formula_logtau = FALSE,
    scale_update   = "none",
    log_tau_basis  = NULL,
    mu_basis       = matrix(-sei, nrow = S, ncol = length(sei), byrow = TRUE),
    current        = current
  )
  row_states <- lapply(seq_len(S), function(row) {
    list(row_index = row, active_key = "negative-selection")
  })

  fast <- .iwmde_selected_normal_location_normalizer_change(
    context      = context,
    setup        = setup,
    basis        = basis,
    values       = values,
    row_states   = row_states,
    active_setup = active_setup
  )

  row_index  <- rep(seq_len(S), each = G)
  grid_index <- rep(seq_len(G), times = S)
  delta      <- values[grid_index] - current[row_index]
  sd <- sqrt(setup[["tau_within"]]^2 +
    matrix(sei^2, nrow = S, ncol = length(sei), byrow = TRUE))
  current_log_norm <- .selection_step_log_norm_matrix(
    mean              = setup[["mu"]],
    sd                = sd,
    sei               = sei,
    selection_context = selection_context
  )
  candidate_context <- BayesTools::selection_context_subset_rows(
    context = selection_context,
    rows    = row_index
  )
  candidate_log_norm <- .selection_step_log_norm_matrix(
    mean              = setup[["mu"]][row_index, , drop = FALSE] +
      basis[["mu_basis"]][row_index, , drop = FALSE] * delta,
    sd                = sd[row_index, , drop = FALSE],
    sei               = sei,
    selection_context = candidate_context
  )
  reference <- rowSums(.apply_log_lik_weights(
    candidate_log_norm - current_log_norm[row_index, , drop = FALSE],
    setup[["weights"]]
  ))

  expect_equal(fast, reference, tolerance = 1e-12)
})


test_that("multilevel weightfunction location fast path dispatches to cluster selected-normal grid", {

  yi  <- c(.10, .20)
  sei <- c(.20, .25)
  data <- list(outcome = data.frame(yi = yi, sei = sei, cluster = c(1L, 1L)))
  attr(data, "outcome_type")     <- "norm"
  attr(data, "effect_direction") <- "positive"
  attr(data, "cluster")          <- TRUE

  values <- c(-.10, .30)
  setup <- list(
    yi                = yi,
    sei               = sei,
    mu                = matrix(c(.05, .10, .15, .20), nrow = 2L, byrow = TRUE),
    tau_within        = matrix(.20, nrow = 2L, ncol = 2L),
    tau_between       = matrix(.10, nrow = 2L, ncol = 2L),
    cluster           = list(`1` = 1:2),
    posterior_samples = matrix(0, nrow = 2L, ncol = 1L)
  )
  basis <- list(
    formula_mu      = FALSE,
    formula_logtau  = FALSE,
    scale_update    = "none",
    log_tau_basis   = NULL,
    mu_basis        = matrix(1, nrow = 2L, ncol = 2L),
    current         = c(0, .10)
  )
  row_states <- lapply(seq_len(2L), function(row) {
    list(
      row_index    = row,
      active_setup = list(is_weightfunction = TRUE)
    )
  })
  log_lik <- matrix(c(1, 2, 3, 4), nrow = length(values))
  log_prior <- c(.1, .2, .3, .4)
  captured <- NULL

  testthat::local_mocked_bindings(
    .iwmde_selection_context_active_branch = function(context, active_setup,
                                                      posterior_samples) {

      list(marker = TRUE)
    },
    .log_lik_cluster_selnorm_location_grid = function(setup, yi, sei, basis,
                                                      current, values,
                                                      selection_context) {

      captured <<- list(
        yi                = yi,
        sei               = sei,
        basis             = basis,
        current           = current,
        values            = values,
        selection_context = selection_context
      )
      log_lik
    },
    .iwmde_predictor_log_prior = function(context, parameter, values,
                                          row_states, replacement) {

      log_prior
    },
    .package = "RoBMA"
  )

  fast <- .iwmde_log_q_grid_normal_location_group(
    context     = list(data = data),
    parameter   = "mu",
    values      = values,
    row_states  = row_states,
    replacement = list(type = "scalar"),
    setup       = setup,
    basis       = basis
  )

  expect_equal(fast, log_lik + matrix(log_prior, nrow = length(values)))
  expect_equal(captured[["yi"]], yi)
  expect_equal(captured[["sei"]], sei)
  expect_equal(captured[["basis"]], basis[["mu_basis"]])
  expect_equal(captured[["current"]], basis[["current"]])
  expect_equal(captured[["values"]], values)
  expect_true(isTRUE(captured[["selection_context"]][["marker"]]))
})


test_that("multilevel weightfunction formula path matches scalar fallback", {

  skip_if_not_certification(
    "The exact multilevel scalar-vs-grid equivalence matrix is certification coverage."
  )

  fit_name  <- "dat.lehmann2018_RoBMA_3lvl_mods_scale"
  parameter <- "mu_Preregistered"
  .skip_if_missing_raw_fits(fit_name)

  context <- .iwmde_context(load_fit(fit_name, validate = FALSE))
  spec    <- .iwmde_parameter_spec(context, parameter, NULL)
  values  <- .iwmde_parameter_values(context, parameter, spec)
  finite  <- is.finite(values)
  component <- .iwmde_parameter_components(context, parameter, spec)
  rows      <- which(component[["active"]] & finite)

  expect_gt(length(rows), 0L)

  row_states <- .iwmde_row_states(context, rows, parameter, spec)
  keep <- vapply(row_states, function(state) {
    isTRUE(state[["active_setup"]][["is_weightfunction"]]) &&
      identical(state[["likelihood_mode"]], "marginal") &&
      is.finite(state[["baseline_log_q"]])
  }, logical(1))
  row_states <- row_states[keep]

  expect_gt(length(row_states), 0L)

  keys <- vapply(row_states, function(state) {
    .iwmde_state_active_key(context, state)
  }, character(1))
  active_key <- names(sort(table(keys), decreasing = TRUE))[[1L]]
  row_states <- head(row_states[keys == active_key], 3L)

  expect_gt(length(row_states), 0L)

  active_values <- values[component[["active"]] & finite]
  grid_values   <- as.numeric(stats::quantile(
    active_values,
    probs = seq(.25, .75, length.out = 3L),
    names = FALSE,
    type  = 8
  ))
  grid_values <- unique(grid_values[is.finite(grid_values)])

  expect_gt(length(grid_values), 0L)

  replacement  <- .iwmde_replacement_spec(context, parameter, spec)
  active_setup <- row_states[[1L]][["active_setup"]]
  setup <- .iwmde_predictor_setup(
    context      = context,
    row_states   = row_states,
    active_setup = active_setup,
    unit         = "cluster"
  )
  basis <- .iwmde_predictor_update_basis(
    context     = context,
    parameter   = parameter,
    row_states  = row_states,
    replacement = replacement,
    setup       = setup
  )

  expect_false(is.null(basis))

  expect_false(isTRUE(basis[["formula_mu"]]))
  expect_true(is.matrix(basis[["mu_basis"]]))
  location_fast <- .iwmde_log_q_grid_normal_location_group(
    context     = context,
    parameter   = parameter,
    values      = grid_values,
    row_states  = row_states,
    replacement = replacement,
    setup       = setup,
    basis       = basis
  )
  expect_true(is.matrix(location_fast))
  fast <- .iwmde_log_q_grid_predictor_group(
    context     = context,
    parameter   = parameter,
    values      = grid_values,
    row_states  = row_states,
    replacement = replacement
  )
  scalar <- .iwmde_log_q_grid_scalar(
    context     = context,
    parameter   = parameter,
    values      = grid_values,
    row_states  = row_states,
    replacement = replacement
  )

  expect_true(is.matrix(fast))
  expect_equal(dim(fast), dim(scalar))
  expect_equal(is.finite(fast), is.finite(scalar))
  finite_terms <- is.finite(fast) & is.finite(scalar)
  expect_true(any(finite_terms))
  expect_equal(fast[finite_terms], scalar[finite_terms], tolerance = 1e-8)
})


test_that("negative-direction multilevel selected-normal location grid matches quadrature reference", {

  skip_if_not(.has_native_selnorm_cluster_location_grid())

  prior <- BayesTools::prior_weightfunction(
    side    = "one-sided",
    steps   = c(.025, .05),
    weights = BayesTools::wf_fixed(c(1, .65, .30))
  )
  yi  <- c(.12, -.18, .31, -.05)
  sei <- c(.20, .28, .24, .18)
  selection_context <- .selection_spec(
    priors           = list(outcome = list(bias = prior)),
    yi               = yi,
    sei              = sei,
    effect_direction = "negative"
  )

  expect_equal(selection_context[["sign"]], -1L)

  S <- 2L
  K <- length(yi)
  setup <- list(
    S           = S,
    mu          = matrix(c(
       .05, .08, .12, .16,
      -.04, .02, .09, .13
    ), nrow = S, byrow = TRUE),
    tau_within  = matrix(c(
      .09, .11, .13, .15,
      .08, .10, .12, .14
    ), nrow = S, byrow = TRUE),
    tau_between = matrix(c(
      .05, .07, .09, .11,
      .04, .06, .08, .10
    ), nrow = S, byrow = TRUE),
    cluster     = list(a = c(1L, 3L), b = c(2L, 4L)),
    weights     = NULL
  )
  selection_context[["omega"]] <- matrix(c(
    1, .65, .30,
    1, .55, .25
  ), nrow = S, byrow = TRUE)
  selection_context[["alpha"]]       <- rep(0, S)
  selection_context[["phack_kind"]]  <- rep(0L, S)
  selection_context[["kernel_mode"]] <- rep(SELKERNEL_STEP, S)

  basis <- matrix(c(
     .04, -.03, .05, -.02,
    -.02,  .06, .01,  .04
  ), nrow = S, byrow = TRUE)
  current <- c(-.10, .15)
  values  <- c(-.20, .05, .30)
  n_gamma <- 11L

  fast <- .log_lik_cluster_selnorm_location_grid(
    setup             = setup,
    yi                = yi,
    sei               = sei,
    basis             = basis,
    current           = current,
    values            = values,
    selection_context = selection_context
  )
  reference <- matrix(NA_real_, nrow = length(values), ncol = S)
  for (s in seq_len(S)) {
    context_s <- BayesTools::selection_context_subset_rows(
      context = selection_context,
      rows    = s
    )
    for (g in seq_along(values)) {
      candidate_setup <- setup
      delta <- values[[g]] - current[[s]]
      candidate_setup[["S"]] <- 1L
      candidate_setup[["mu"]] <- setup[["mu"]][s, , drop = FALSE] +
        basis[s, , drop = FALSE] * delta
      candidate_setup[["tau_within"]]  <- setup[["tau_within"]][s, , drop = FALSE]
      candidate_setup[["tau_between"]] <- setup[["tau_between"]][s, , drop = FALSE]
      reference[g, s] <- sum(.log_lik_cluster_norm_quadrature_r(
        setup             = candidate_setup,
        yi                = yi,
        sei               = sei,
        is_weightfunction = TRUE,
        selection_context = context_s,
        n_gamma           = n_gamma
      ))
    }
  }

  expect_true(is.matrix(fast))
  expect_equal(dim(fast), dim(reference))
  quadrature_change <- attr(fast, "quadrature_relative_change", exact = TRUE)
  quadrature_order <- attr(fast, "quadrature_order", exact = TRUE)
  expect_equal(dim(quadrature_change), dim(fast))
  expect_true(all(is.finite(quadrature_change)))
  expect_equal(dim(quadrature_order), dim(fast))
  expect_true(all(quadrature_order %in% c(15L, 31L)))
  attr(fast, "quadrature_relative_change") <- NULL
  attr(fast, "quadrature_order") <- NULL
  expect_equal(fast, reference, tolerance = 1e-7)
})


test_that("IWMDE selected-normal native normalizer skips fallback candidates", {

  G <- 3L
  S <- 2L
  K <- 2L
  setup <- list(
    mu                = matrix(c(.1, .2, .3, .4), nrow = S),
    tau_within        = c(.05, .10),
    sei               = c(.20, .30),
    weights           = NULL,
    posterior_samples = matrix(0, nrow = S, ncol = 1)
  )
  basis <- list(
    log_tau_basis  = NULL,
    scale_update   = "none",
    formula_logtau = FALSE,
    mu_basis       = NULL,
    current        = c(.10, .30)
  )
  row_states <- list(
    list(row_index = 1L),
    list(row_index = 2L)
  )
  selection_context <- list(
    omega       = matrix(1, nrow = S, ncol = 1),
    alpha       = rep(NA_real_, S),
    phack_kind  = rep(NA_integer_, S),
    kernel_mode = rep(SELKERNEL_NORMAL, S)
  )

  testthat::local_mocked_bindings(
    .iwmde_selection_context_active_branch = function(context, active_setup,
                                                      posterior_samples) {

      selection_context
    },
    .iwmde_selected_normal_current_log_norm = function(context, setup, sd,
                                                       selection_context,
                                                       row_states) {

      matrix(0, nrow = S, ncol = K)
    },
    .has_native_selnorm_log_norm_delta = function() TRUE,
    .selnorm_kernel_log_norm_delta_grid = function(
        mean, sd, basis, current_log_norm, current, values, sei, weights,
        omega, selection_spec, alpha, phack_kind, kernel_mode) {

      expect_equal(mean, setup[["mu"]])
      expect_null(basis)
      matrix(.25, nrow = length(values), ncol = length(current))
    },
    .package = "RoBMA"
  )

  out <- .iwmde_selected_normal_location_normalizer_change(
    context      = list(),
    setup        = setup,
    basis        = basis,
    values       = seq(-1, 1, length.out = G),
    row_states   = row_states,
    active_setup = list()
  )

  expect_equal(out, rep(.25, G * S))
})


test_that("IWMDE active-branch samples reject uncollapsed prior mixtures", {

  prior <- BayesTools::prior_mixture(
    prior_list = list(
      BayesTools::prior("spike", parameters = list(location = 0)),
      BayesTools::prior("normal", parameters = list(mean = 0, sd = 1))
    ),
    is_null = c(TRUE, FALSE)
  )
  samples <- matrix(
    2,
    nrow     = 1,
    dimnames = list(NULL, "mu_indicator")
  )

  expect_error(
    .iwmde_likelihood_posterior_samples(
      context      = list(),
      samples      = samples,
      active_setup = list(
        priors            = list(outcome = list(mu = prior)),
        is_weightfunction = FALSE
      )
    ),
    "requires priors localized to a single mixture component"
  )
})

test_that("IWMDE batched q evaluation warns on invalid likelihood length", {

  context <- list(
    posterior_samples = matrix(
      0,
      nrow     = 1,
      dimnames = list(NULL, "mu")
    ),
    row_cache = new.env(parent = emptyenv())
  )
  row_states <- list(
    list(
      row_index       = 1L,
      active_key      = "same",
      active_setup    = list(),
      likelihood_mode = "marginal",
      prior_list      = list()
    )
  )

  expect_warning(
    out <- .iwmde_log_q_grid_from_samples(
      context         = context,
      parameter       = "mu",
      values          = c(0, 1),
      row_states      = row_states,
      replacement     = list(type = "scalar"),
      likelihood_mode = "marginal",
      log_lik_fun     = function(samples, active_setup, batch) 0
    ),
    "invalid length"
  )
  expect_null(out)
})


test_that("IWMDE replacement keeps inverse-gamma auxiliary parameters synced", {

  prior <- BayesTools::prior(
    "invgamma",
    parameters = list(2, .5)
  )
  context <- list(flat_prior_list = list(tau = prior))
  samples <- matrix(
    c(.5, 99,
      2, 99),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(NULL, c("tau", "inv_tau"))
  )

  synced <- .iwmde_sync_invgamma_auxiliary_matrix(
    context    = context,
    samples    = samples,
    parameters = "tau"
  )
  bad <- .iwmde_sync_invgamma_auxiliary_matrix(
    context    = context,
    samples    = samples,
    parameters = "mu"
  )
  row <- .iwmde_sync_invgamma_auxiliary_row(
    context    = context,
    row        = c(tau = 4, inv_tau = 99),
    parameters = "tau"
  )
  invalid <- .iwmde_sync_invgamma_auxiliary_row(
    context    = context,
    row        = c(tau = 0, inv_tau = 99),
    parameters = "tau"
  )

  expect_false(.iwmde_can_use_focal_prior_delta(prior))
  expect_true(all(synced[["valid"]]))
  expect_equal(synced[["samples"]][, "inv_tau"], 1 / samples[, "tau"])
  expect_equal(bad[["samples"]], samples)
  expect_true(isTRUE(row[["valid"]]))
  expect_equal(row[["row"]][["inv_tau"]], .25)
  expect_false(invalid[["valid"]])
  expect_true(is.na(invalid[["row"]][["inv_tau"]]))
})


test_that("IWMDE inverse-gamma auxiliary matrix sync reuses active states", {

  inv_prior <- BayesTools::prior(
    "invgamma",
    parameters = list(2, .5)
  )
  normal_prior <- BayesTools::prior(
    "normal",
    parameters = list(0, 1)
  )
  context <- list(
    flat_prior_list = list(
      tau = BayesTools::prior_mixture(list(inv_prior, normal_prior))
    ),
    indicator_names = "tau_indicator"
  )
  samples <- matrix(
    c(.5, 99, 1,
      2, 99, 1,
      0, 99, 1,
      0, 99, 2,
      4, 99, 2,
      5, 99, 1),
    nrow = 6,
    byrow = TRUE,
    dimnames = list(NULL, c("tau", "inv_tau", "tau_indicator"))
  )

  calls <- 0L
  testthat::local_mocked_bindings(
    .iwmde_row_uses_invgamma_prior = function(context, parameter, row) {
      calls <<- calls + 1L
      row[["tau_indicator"]] == 1
    },
    .package = "RoBMA"
  )
  synced <- .iwmde_sync_invgamma_auxiliary_matrix(
    context    = context,
    samples    = samples,
    parameters = "tau"
  )

  expect_equal(calls, 2L)
  expect_equal(
    synced[["samples"]][c(1, 2, 6), "inv_tau"],
    1 / samples[c(1, 2, 6), "tau"]
  )
  expect_true(is.na(synced[["samples"]][3, "inv_tau"]))
  expect_equal(synced[["samples"]][c(4, 5), "inv_tau"], samples[c(4, 5), "inv_tau"])
  expect_equal(synced[["valid"]], c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE))
})


.iwmde_active_cases <- function(cases) {

  active_names <- active_fit_catalog()[["name"]]
  cases <- Filter(function(case) case[[1L]] %in% active_names, cases)
  if (length(cases) == 0L) {
    testthat::skip("No fixtures for this certification matrix are active.")
  }

  return(cases)
}


test_that("IWMDE batched q evaluation matches scalar fallback", {

  skip_if_not_certification(
    "The exact batched-vs-scalar likelihood matrix is certification coverage."
  )

  cases <- list(
    list("bcg_meta-regression", "mu_ablat", NULL),
    list("bangertdrowns2004_location-scale", "log_tau_intercept", NULL),
    list("dat.lehmann2018-3PSM", "mu", NULL),
    list("dat.lehmann2018-3PSM_neg", "mu", NULL),
    list("dat.lehmann2018-PET_neg", "PET", NULL),
    list("dat.lehmann2018-PEESE_neg", "PEESE", NULL),
    list("dat.lehmann2018_RoBMA_3lvl_mods_scale", "mu_Preregistered", NULL),
    list("bcg_BMA.glmm_3lvl_location_scale", "mu_year", NULL),
    list("konstantopoulos2011_3lvl", "gamma[1]", NULL),
    list("bcg_glmm", "theta[1]", NULL),
    list("bcg_glmm", "pi[1]", NULL),
    list("nielweise2008_glmm", "phi[1]", NULL),
    list(
      "bcg_meta-regression",
      "mu_intercept + mu_ablat",
      list(type = "linear", weights = c(mu_intercept = 1, mu_ablat = 1))
    )
  )
  cases <- .iwmde_active_cases(cases)

  for (case in cases) {
    .expect_iwmde_batch_equals_scalar(
      fit_name       = case[[1L]],
      parameter      = case[[2L]],
      parameter_spec = case[[3L]]
    )
  }
})

test_that("IWMDE predictor fast path matches scalar fallback", {

  skip_if_not_certification(
    "The exhaustive predictor fast-path equivalence matrix is certification coverage."
  )

  cases <- list(
    list("bcg_meta-analysis", "mu", NULL),
    list("brma.mv_block_mvn_fixed_random_null", "mu", NULL),
    list("brma.mv_block_mvn", "mu", NULL),
    list("bcg_meta-analysis", "tau", NULL),
    list("konstantopoulos2011_3lvl", "rho", NULL),
    list("dat.lehmann2018-PET", "PET", NULL),
    list("dat.lehmann2018-PEESE", "PEESE", NULL),
    list("dat.lehmann2018-PET_neg", "PET", NULL),
    list("dat.lehmann2018-PEESE_neg", "PEESE", NULL),
    list("dat.lehmann2018-3PSM_neg", "mu", NULL),
    list("bcg_meta-regression", "mu_ablat", NULL),
    list("bangertdrowns2004_location-scale", "log_tau_intercept", NULL),
    list(
      "bcg_meta-regression",
      "mu_intercept + mu_ablat",
      list(type = "linear", weights = c(mu_intercept = 1, mu_ablat = 1))
    )
  )
  cases <- .iwmde_active_cases(cases)

  for (case in cases) {
    .expect_iwmde_predictor_fast_equals_scalar(
      fit_name       = case[[1L]],
      parameter      = case[[2L]],
      parameter_spec = case[[3L]]
    )
  }
})

test_that("IWMDE normal quadratic fast path matches scalar fallback", {

  skip_if_not_certification(
    "The exhaustive quadratic fast-path equivalence matrix is certification coverage."
  )

  cases <- list(
    list("bcg_meta-analysis", "mu", NULL),
    list("brma.mv_block_mvn_fixed_random_null", "mu", NULL),
    list("brma.mv_block_mvn", "mu", NULL),
    list("dat.lehmann2018-3PSM", "mu", NULL),
    list("dat.lehmann2018-PET", "PET", NULL),
    list("dat.lehmann2018-PEESE", "PEESE", NULL),
    list("dat.lehmann2018-PET_neg", "PET", NULL),
    list("dat.lehmann2018-PEESE_neg", "PEESE", NULL),
    list("dat.lehmann2018-3PSM_neg", "mu", NULL),
    list("konstantopoulos2011_3lvl", "mu", NULL)
  )
  cases <- .iwmde_active_cases(cases)

  for (case in cases) {
    .expect_iwmde_normal_quadratic_equals_scalar(
      fit_name       = case[[1L]],
      parameter      = case[[2L]],
      parameter_spec = case[[3L]]
    )
  }
})


test_that("known-V normal q-grid caches invariant spectral blocks", {

  fit_name <- "brma.mv_block_mvn_fixed_random_null"
  .skip_if_missing_raw_fits(fit_name)

  context   <- .iwmde_context(load_fit(fit_name, validate = FALSE))
  parameter <- "mu"
  spec      <- .iwmde_parameter_spec(context, parameter, NULL)
  values    <- .iwmde_parameter_values(context, parameter, spec)
  component <- .iwmde_parameter_components(context, parameter, spec)
  rows      <- which(component[["active"]] & is.finite(values))
  rows      <- head(rows, 3L)
  row_states <- .iwmde_row_states(context, rows, parameter, spec)
  keep       <- vapply(row_states, function(state) {
    is.finite(state[["baseline_log_q"]])
  }, logical(1))
  row_states <- row_states[keep]

  expect_gt(length(row_states), 0L)

  active_setup <- row_states[[1L]][["active_setup"]]
  setup <- .iwmde_predictor_setup(
    context      = context,
    row_states   = row_states,
    active_setup = active_setup,
    unit         = "estimate"
  )
  replacement <- .iwmde_replacement_spec(context, parameter, spec)
  basis <- .iwmde_predictor_update_basis(
    context     = context,
    parameter   = parameter,
    row_states  = row_states,
    replacement = replacement,
    setup       = setup
  )
  grid <- seq(min(values[rows]), max(values[rows]), length.out = 20L)

  predictor_cache <- context[["predictor_cache"]]
  spectral_keys <- grep(
    "^known_v_location_spectral_blocks\\|",
    ls(envir = predictor_cache),
    value = TRUE
  )
  if (length(spectral_keys) > 0L) {
    rm(list = spectral_keys, envir = predictor_cache)
  }

  out <- .iwmde_log_q_grid_normal_location_group(
    context     = context,
    parameter   = parameter,
    values      = grid,
    row_states  = row_states,
    replacement = replacement,
    setup       = setup,
    basis       = basis
  )
  cached_keys <- grep(
    "^known_v_location_spectral_blocks\\|",
    ls(envir = predictor_cache),
    value = TRUE
  )
  reused <- .iwmde_log_q_grid_normal_location_group(
    context     = context,
    parameter   = parameter,
    values      = rev(grid),
    row_states  = row_states,
    replacement = replacement,
    setup       = setup,
    basis       = basis
  )
  reused_keys <- grep(
    "^known_v_location_spectral_blocks\\|",
    ls(envir = predictor_cache),
    value = TRUE
  )

  expect_true(is.matrix(out))
  expect_equal(dim(out), c(length(grid), length(row_states)))
  expect_equal(reused, out[nrow(out):1L, , drop = FALSE], tolerance = 1e-12)
  expect_length(cached_keys, 1L)
  expect_identical(reused_keys, cached_keys)
})


test_that("known-V grouped baseline states match scalar construction", {

  fit_name <- "brma.mv_block_mvn_random"
  .skip_if_missing_raw_fits(fit_name)

  fit             <- load_fit(fit_name, validate = FALSE)
  parameter       <- "mu"
  grouped_context <- .iwmde_context(fit)
  scalar_context  <- .iwmde_context(fit)
  spec            <- .iwmde_parameter_spec(grouped_context, parameter, NULL)
  values          <- .iwmde_parameter_values(grouped_context, parameter, spec)
  component       <- .iwmde_parameter_components(
    grouped_context,
    parameter,
    spec
  )
  rows <- head(which(component[["active"]] & is.finite(values)), 5L)

  grouped <- .iwmde_row_states_grouped_marginal(
    context        = grouped_context,
    rows           = rows,
    parameter      = parameter,
    parameter_spec = spec,
    estimator      = "q_grid_cmde"
  )
  scalar <- .iwmde_row_states(
    context        = scalar_context,
    rows           = rows,
    parameter      = parameter,
    parameter_spec = spec,
    estimator      = "q_grid_cmde"
  )

  fields <- c(
    "row_index", "active_key", "baseline_log_lik", "baseline_log_prior",
    "baseline_focal_log_prior", "use_focal_prior_delta", "baseline_log_q",
    "likelihood_mode", "state_scope"
  )
  expect_equal(
    lapply(grouped, `[`, fields),
    lapply(scalar, `[`, fields),
    tolerance = 1e-10
  )
})


test_that("ordinary marginal baseline states match scalar construction", {

  fit_names  <- c("bcg_meta-analysis", "bcg_meta-regression")
  parameters <- c("mu", "mu_ablat")
  .skip_if_missing_raw_fits(fit_names)

  for (i in seq_along(fit_names)) {
    fit             <- load_fit(fit_names[[i]], validate = FALSE)
    parameter       <- parameters[[i]]
    grouped_context <- .iwmde_context(fit)
    scalar_context  <- .iwmde_context(fit)
    spec            <- .iwmde_parameter_spec(grouped_context, parameter, NULL)
    values          <- .iwmde_parameter_values(
      grouped_context,
      parameter,
      spec
    )
    component <- .iwmde_parameter_components(
      grouped_context,
      parameter,
      spec
    )
    rows <- head(which(component[["active"]] & is.finite(values)), 5L)

    grouped <- .iwmde_row_states_grouped_marginal(
      context        = grouped_context,
      rows           = rows,
      parameter      = parameter,
      parameter_spec = spec,
      estimator      = "q_grid_cmde"
    )
    scalar <- .iwmde_row_states(
      context        = scalar_context,
      rows           = rows,
      parameter      = parameter,
      parameter_spec = spec,
      estimator      = "q_grid_cmde"
    )

    fields <- c(
      "row_index", "active_key", "baseline_log_lik", "baseline_log_prior",
      "baseline_focal_log_prior", "use_focal_prior_delta", "baseline_log_q",
      "likelihood_mode", "state_scope"
    )
    expect_equal(
      lapply(grouped, `[`, fields),
      lapply(scalar, `[`, fields),
      tolerance = 1e-10
    )
  }
})
