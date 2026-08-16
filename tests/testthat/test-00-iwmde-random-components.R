test_that("simplex endpoint replacement preserves row-wise prior densities", {

  parameter <- "rho"
  eta <- .iwmde_simplex_auxiliary_columns(parameter, 2L)
  samples <- matrix(
    c(0, 0, 4, 1, 2, 0),
    nrow = 2L,
    byrow = TRUE,
    dimnames = list(NULL, c("theta", eta))
  )
  prior_list <- list(
    theta = BayesTools::prior("normal", list(mean = 0, sd = 1)),
    rho   = BayesTools::prior("dirichlet", list(alpha = c(1, 1)))
  )
  replacement <- list(
    type              = "simplex_pair",
    parameter         = parameter,
    auxiliary_columns = eta
  )

  actual <- .iwmde_replacement_log_prior_rows(
    samples     = samples,
    prior_list  = prior_list,
    replacement = replacement
  )
  expected <- stats::dnorm(samples[, "theta"], log = TRUE) +
    stats::dgamma(samples[, eta[[1L]]], shape = 1, log = TRUE) +
    stats::dgamma(samples[, eta[[2L]]], shape = 1, log = TRUE)

  expect_equal(actual, expected, tolerance = 0)
  expect_true(all(is.finite(actual)))
})

test_that("diagonal allocation grid equals the full marginal covariance", {

  yi <- c(0.1, -0.2, 0.3, 0.05)
  vi <- c(0.04, 0.05, 0.06, 0.07)
  cluster_map <- c(1L, 1L, 2L, 2L)
  posterior_samples <- matrix(
    c(0.2, 0.45),
    ncol = 1L,
    dimnames = list(NULL, "tau")
  )
  mu_samples <- rbind(
    c(0.02, 0.02, -0.01, -0.01),
    c(-0.03, 0.01, 0.04, 0.02)
  )
  values <- c(0, 0.4, 1)
  row_states <- list(list(row_index = 1L), list(row_index = 2L))
  context <- list(
    posterior_samples = posterior_samples,
    data = list(outcome = data.frame(yi = yi))
  )
  testthat::local_mocked_bindings(
    .iwmde_known_v_random_allocation_cluster_plan = function(...) {
      list(
        source_parameter           = "tau",
        cluster_map                = cluster_map,
        target_is_cluster_fraction = TRUE,
        sampling_variance          = vi
      )
    },
    .iwmde_predictor_evaluate_fixed_mu = function(...) mu_samples,
    .iwmde_predictor_log_prior = function(...) {
      numeric(length(values) * length(row_states))
    },
    .data_effect_direction = function(...) "positive",
    .package = "RoBMA"
  )

  actual <- .iwmde_log_q_grid_known_v_random_allocation_group(
    context      = context,
    parameter    = "rho[2]",
    values       = values,
    row_states   = row_states,
    replacement  = list(),
    active_setup = list()
  )
  expected <- matrix(NA_real_, nrow = length(values), ncol = nrow(mu_samples))
  same_cluster <- outer(cluster_map, cluster_map, "==")
  for (value_i in seq_along(values)) {
    for (draw in seq_len(nrow(mu_samples))) {
      tau <- posterior_samples[draw, "tau"]
      covariance <- diag(vi + tau^2 * (1 - values[[value_i]])) +
        tau^2 * values[[value_i]] * same_cluster
      expected[value_i, draw] <- .marglik_mvn_log_density(
        y          = yi,
        mean       = mu_samples[draw, ],
        covariance = covariance
      )
    }
  }

  expect_equal(actual, expected, tolerance = 1e-12)
})

test_that("semantic random qCMDE hypotheses use the plotting density target", {

  selected <- list(
    entry = list(term = "study variance fraction"),
    spec = list(
      summary_type     = "var_frac",
      source_parameter = "rho",
      label            = "var_frac(random_total: study)"
    ),
    samples      = matrix(seq(0.01, 0.99, length.out = 100L), ncol = 1L),
    prior        = BayesTools::prior("beta", list(alpha = 1, beta = 1)),
    source_prior = BayesTools::prior(
      "dirichlet",
      list(alpha = c(1, 1))
    )
  )
  attached_values <- numeric()
  target_spec <- list(
    type                 = "simplex_pair",
    parameter            = "rho",
    index                = 2L,
    n_targets            = 2L,
    target_columns       = paste0("rho[", 1:2, "]"),
    auxiliary_columns    = paste0("rho_eta[", 1:2, "]"),
    conditioning_exclude = paste0("rho[", 1:2, "]")
  )
  used_density_method <- NULL
  attachment_calls    <- 0L
  testthat::local_mocked_bindings(
    .brma_random_parameter_select = function(...) selected,
    .brma_random_parameter_density_target = function(...) {
      list(parameter = "rho[2]", parameter_spec = target_spec)
    },
    .brma_random_parameter_mixed_posterior = function(...) list(theta = 1:3),
    .iwmde_context = function(...) list(),
    .iwmde_estimate_cache = function(...) new.env(parent = emptyenv()),
    .hypothesis_brma_attach_iwmde_scalar = function(
        posterior, value, parameter, parameter_spec, ...) {
      attachment_calls <<- attachment_calls + 1L
      attached_values <<- c(attached_values, value)
      expect_identical(parameter, "rho[2]")
      expect_identical(parameter_spec, target_spec)
      posterior
    },
    .hypothesis_brma_append_iwmde_warnings = function(table, ...) table,
    .package = "RoBMA"
  )
  testthat::local_mocked_bindings(
    marginal_posterior = function(...) seq(0.01, 0.99, length.out = 100L),
    hypothesis_BF = function(..., density_method) {
      used_density_method <<- density_method
      structure(data.frame(BF = 1), class = c(
        "BayesTools_hypothesis_BF",
        "data.frame"
      ))
    },
    .package = "BayesTools"
  )

  out <- .hypothesis_brma_random(
    object                    = list(),
    parameter                 = "theta",
    hypothesis                = BayesTools::hypothesis_parse(c(
      "theta != 0 vs theta = 0",
      "theta != 1 vs theta = 1"
    )),
    standardized_coefficients = FALSE,
    conditional               = FALSE,
    logBF                     = FALSE,
    BF01                      = FALSE,
    seed                      = 1,
    density_method            = "qCMDE",
    density_control           = list(
      n_points             = 20L,
      samples              = 50L,
      target_relative_mcse = 0.05,
      normalization_points = 50L,
      normalization_prob   = 0.999
    ),
    n_samples = 100L,
    columns   = "default"
  )

  expect_s3_class(out, "BayesTools_hypothesis_BF")
  expect_identical(attachment_calls, 1L)
  expect_identical(attached_values, c(0, 1))
  expect_identical(used_density_method, "precomputed")
})

test_that("batched ordinates retain value-specific diagnostics", {

  values <- c(0, 1)
  prior_ordinates <- .iwmde_prior_ordinate_classifications(
    prior_density = BayesTools::prior("beta", list(alpha = 1, beta = 1)),
    values        = values
  )
  density <- list(
    x                                = values,
    evaluation_x                     = values,
    y                                = c(0.25, 0.5),
    pilot_y                          = c(0.24, 0.48),
    validation_y                     = c(0.26, 0.52),
    ordinate_relative_change         = c(0.01, 0.02),
    ordinate_log_change              = c(0.01, 0.02),
    pilot_ordinate_relative_change   = c(0.01, 0.02),
    pilot_ordinate_log_change        = c(0.01, 0.02),
    mcse                             = c(0.01, 0.02),
    relative_mcse                    = c(0.04, 0.08),
    active_branch_mcse               = c(0.01, 0.02),
    active_branch_relative_mcse      = c(0.04, 0.08),
    active_mass_component_mcse       = c(0, 0),
    sampling_mcse                    = c(0.03, 0.04),
    sampling_relative_mcse           = c(0.12, 0.16),
    finite_terms                     = c(500L, 450L),
    ess                              = c(300, 40),
    max_weight_share                 = c(0.02, 0.25),
    max_log_ratio                    = c(1, 2),
    estimator                        = "q_grid_cmde"
  )
  diagnostics <- .mock_iwmde_good_diagnostics(
    estimator  = "q_grid_cmde",
    rows       = 500L,
    value      = 0,
    include_bf = TRUE
  )
  diagnostics[["prior_ordinates"]]   <- prior_ordinates
  diagnostics[["ordinate_warnings"]] <- character()
  diagnostic <- list(
    status       = "ok",
    parameter    = "rho[2]",
    target_key   = "simplex_pair|rho|2",
    point_masses = .mock_iwmde_empty_point_masses(),
    iwmde        = density,
    diagnostics  = diagnostics,
    plan = list(
      plan_key           = "shared-plan",
      source_fingerprint = list(source = "test"),
      prior_ordinates    = prior_ordinates,
      target = list(
        target_key = "simplex_pair|rho|2",
        metadata   = list(parameter = "rho[2]")
      )
    )
  )

  attributes <- .iwmde_posterior_ordinate_attributes(
    diagnostic      = diagnostic,
    density_method  = "qCMDE",
    density_control = list(samples = 500L)
  )
  entries <- .iwmde_posterior_ordinate_entries(attributes[["accepted"]])

  expect_equal(vapply(entries, `[[`, numeric(1), "value"), values)
  expect_equal(vapply(entries, `[[`, numeric(1), "ordinate"), density[["y"]])
  expect_equal(vapply(entries, function(entry) {
    entry[["diagnostics"]][["relative_mcse"]]
  }, numeric(1)), density[["relative_mcse"]])
  expect_equal(vapply(entries, function(entry) {
    entry[["diagnostics"]][["sampling_relative_mcse"]]
  }, numeric(1)), density[["sampling_relative_mcse"]])
  expect_equal(vapply(entries, function(entry) {
    entry[["diagnostics"]][["ess"]]
  }, numeric(1)), density[["ess"]])
  expect_true(all(vapply(entries, function(entry) {
    length(entry[["iwmde_provenance"]][["prior_ordinates"]]) == 1L
  }, logical(1))))
  expect_identical(attributes[["rejected"]], NULL)
})

test_that("separate random targets share one hypothesis result", {

  hypothesis <- BayesTools::hypothesis_parse(c(
    "`var_frac(random_total: study)` != 0 vs `var_frac(random_total: study)` = 0",
    "`sd_total(random_total)` = 0",
    "`var_frac(random_total: study)` != 1 vs `var_frac(random_total: study)` = 1"
  ))
  selections <- list(
    list(component = "random", parameter = "rho[2]"),
    list(component = "random", parameter = "tau"),
    list(component = "random", parameter = "rho[2]")
  )
  calls <- list()
  testthat::local_mocked_bindings(
    hypothesis.brma = function(object, hypothesis, ...) {

      n <- length(hypothesis)
      is_allocation <- grepl("var_frac", hypothesis[[1L]], fixed = TRUE)
      out <- BayesTools::hypothesis_BF(
        posterior  = if (is_allocation) {
          seq(-1, 2, length.out = 101L)
        } else {
          seq(-2, 1, length.out = 101L)
        },
        prior      = seq(-2, 2, length.out = 101L),
        hypothesis = paste("theta >", seq_len(n) / 10),
        parameter  = "theta"
      )
      attr(out, "warnings") <- stats::setNames(
        paste("warning", seq_len(n)),
        rownames(out)
      )
      calls[[length(calls) + 1L]] <<- list(
        hypothesis = hypothesis,
        result     = out
      )
      out
    },
    .package = "RoBMA"
  )

  out <- .hypothesis_brma_multiple_parameters(
    object                    = list(),
    hypothesis                = hypothesis,
    selections                = selections,
    standardized_coefficients = FALSE,
    conditional               = FALSE,
    conditional_omitted       = TRUE,
    logBF                     = FALSE,
    BF01                      = FALSE,
    seed                      = 1,
    density_method            = "qCMDE",
    density_control           = list(),
    n_samples                 = 100L,
    columns                   = "default"
  )

  expect_length(calls, 2L)
  expect_identical(lengths(lapply(calls, `[[`, "hypothesis")), c(2L, 1L))
  allocation_BF <- attr(calls[[1L]][["result"]], "raw_BF", exact = TRUE)
  total_BF      <- attr(calls[[2L]][["result"]], "raw_BF", exact = TRUE)
  expect_identical(
    attr(out, "raw_BF"),
    c(allocation_BF[[1L]], total_BF, allocation_BF[[2L]])
  )
  expect_identical(
    names(attr(out, "warnings")),
    c("theta", "theta2", "theta1")
  )
  expect_identical(attr(out, "hypothesis_ast"), hypothesis)
  expect_s3_class(out, "BayesTools_hypothesis_BF")
})
