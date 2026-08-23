# ============================================================================ #
# test-02-iwmde-scalar-baseline.R
# ============================================================================ #

context("IWMDE scalar baseline")

source(testthat::test_path("common-functions.R"))
source(testthat::test_path("helper-iwmde.R"))


test_that("IWMDE mixture indicators are validated before active-prior selection", {

  prior <- BayesTools::prior_mixture(list(
    BayesTools::prior("normal", parameters = list(mean = 0, sd = 1)),
    BayesTools::prior_none()
  ))

  expect_identical(
    .iwmde_select_prior_component(
      prior     = prior,
      row       = c(mu_indicator = 1),
      indicator = "mu_indicator"
    ),
    prior[[1L]]
  )
  expect_error(
    .iwmde_select_prior_component(
      prior     = prior,
      row       = c(mu_indicator = 1.25),
      indicator = "mu_indicator"
    ),
    "integer-valued"
  )
  expect_error(
    .iwmde_select_prior_component(
      prior     = prior,
      row       = c(mu_indicator = 3),
      indicator = "mu_indicator"
    ),
    "outside the available prior components"
  )
})

test_that("IWMDE vectorized prior states preserve product-space components", {

  samples <- cbind(
    mu            = seq_len(6) / 10,
    tau           = seq_len(6) / 20,
    mu_indicator  = c(1, 2, 1, 2, 2, 1),
    tau_indicator = c(1, 1, 2, 2, 1, 2)
  )
  context <- .iwmde_context_ensure_caches(list(
    posterior_samples = samples,
    indicator_names   = c("mu_indicator", "tau_indicator"),
    flat_prior_list   = list(
      mu = BayesTools::prior_mixture(list(
        BayesTools::prior("normal", parameters = list(mean = 0, sd = 1)),
        BayesTools::prior_none()
      )),
      tau = BayesTools::prior_mixture(list(
        BayesTools::prior("spike", parameters = list(location = .5)),
        BayesTools::prior("normal", parameters = list(mean = 0, sd = 1))
      ))
    )
  ))

  mu_component <- .iwmde_parameter_components(context, "mu")
  expect_identical(mu_component[["active"]], samples[, "mu_indicator"] == 1)
  expect_equal(
    mu_component[["point_location"]],
    ifelse(samples[, "mu_indicator"] == 2, 0, NA_real_)
  )

  linear_component <- .iwmde_parameter_components(
    context,
    "mu_plus_tau",
    list(type = "linear", weights = c(mu = 1, tau = 2))
  )
  expect_identical(
    linear_component[["active"]],
    samples[, "mu_indicator"] == 1 | samples[, "tau_indicator"] == 2
  )
  expect_equal(
    linear_component[["point_location"]],
    ifelse(
      samples[, "mu_indicator"] == 2 & samples[, "tau_indicator"] == 1,
      1,
      NA_real_
    )
  )
  expected_keys <- c(
    "mu_indicator=1|tau_indicator=1",
    "mu_indicator=2|tau_indicator=1",
    "mu_indicator=1|tau_indicator=2",
    "mu_indicator=2|tau_indicator=2",
    "mu_indicator=2|tau_indicator=1",
    "mu_indicator=1|tau_indicator=2"
  )
  expect_identical(.iwmde_active_keys(context), expected_keys)
  expect_identical(.iwmde_active_keys(context), expected_keys)

  invalid_context <- context
  invalid_context[["posterior_samples"]][1, "mu_indicator"] <- 1.5
  invalid_context[["row_cache"]] <- new.env(parent = emptyenv())
  expect_error(.iwmde_active_keys(invalid_context), "integer-valued")
})

test_that("IWMDE likelihood does not mask internal errors", {

  data <- list(outcome = data.frame(yi = 0, sei = 1))
  attr(data, "outcome_type")     <- "norm"
  attr(data, "effect_direction") <- "positive"

  context <- list(
    object = list(fit = list()),
    data   = data
  )
  active_setup <- list(
    fit_data          = list(),
    priors            = list(),
    is_PET            = FALSE,
    is_PEESE          = FALSE,
    is_weightfunction = FALSE
  )

  testthat::local_mocked_bindings(
    .log_posterior = function(...) {
      stop("sentinel conditional likelihood error", call. = FALSE)
    }
  )

  expect_error(
    .iwmde_log_likelihood_parameters(
      context      = context,
      parameters   = list(mu = 0, tau = 1),
      active_setup = active_setup
    ),
    "sentinel conditional likelihood error"
  )
})

test_that("IWMDE marginal likelihood does not mask internal errors", {

  data <- list(outcome = data.frame(yi = 0, sei = 1))
  attr(data, "outcome_type")     <- "norm"
  attr(data, "effect_direction") <- "positive"

  context <- list(
    object = list(fit = list()),
    data   = data
  )
  active_setup <- list(
    fit_data          = list(K = 1),
    priors            = list(),
    is_PET            = FALSE,
    is_PEESE          = FALSE,
    is_weightfunction = FALSE
  )

  testthat::local_mocked_bindings(
    .iwmde_log_likelihood_parameters_marginal_internal = function(...) {
      stop("sentinel marginal likelihood error", call. = FALSE)
    }
  )

  expect_error(
    .iwmde_log_likelihood_parameters(
      context         = context,
      parameters      = list(mu = 0, tau = 1),
      active_setup    = active_setup,
      likelihood_mode = "marginal"
    ),
    "sentinel marginal likelihood error"
  )
})

test_that("IWMDE scalar evaluators preserve zeros and expose invalid results", {

  data <- list(outcome = data.frame(yi = 0, sei = 1))
  attr(data, "outcome_type")     <- "norm"
  attr(data, "effect_direction") <- "positive"
  context <- list(
    object = list(fit = list()),
    data   = data
  )
  active_setup <- list(
    fit_data          = list(),
    priors            = list(),
    is_PET            = FALSE,
    is_PEESE          = FALSE,
    is_weightfunction = FALSE
  )
  conditional_value <- 0
  marginal_value    <- 0
  full_prior_value  <- 0
  focal_prior_value <- 0

  testthat::local_mocked_bindings(
    .log_posterior = function(...) conditional_value,
    .iwmde_log_likelihood_parameters_marginal_internal = function(...) {
      marginal_value
    },
    .resolve_fixed_prior_row = function(row, prior_list) row
  )
  testthat::local_mocked_bindings(
    JAGS_marglik_priors = function(...) full_prior_value,
    is.prior.none       = function(...) FALSE,
    lpdf                = function(...) focal_prior_value,
    .package            = "BayesTools"
  )

  evaluate <- function(value) {
    conditional_value <<- value
    marginal_value    <<- value
    full_prior_value  <<- value
    focal_prior_value <<- value

    list(
      conditional = .iwmde_log_likelihood_parameters(
        context      = context,
        parameters   = list(mu = 0, tau = 1),
        active_setup = active_setup
      ),
      marginal = .iwmde_log_likelihood_parameters(
        context         = context,
        parameters      = list(mu = 0, tau = 1),
        active_setup    = active_setup,
        likelihood_mode = "marginal"
      ),
      full_prior = .iwmde_log_prior_row(
        row        = c(mu = 0),
        prior_list = list(mu = structure(list(), class = "mock_prior"))
      ),
      focal_prior = .iwmde_focal_log_prior(
        prior = structure(list(), class = "mock_prior"),
        value = 0
      )
    )
  }

  exact_zero <- evaluate(-Inf)
  expect_true(all(vapply(exact_zero, identical, logical(1), -Inf)))

  missing <- evaluate(NA_real_)
  expect_true(all(vapply(missing, function(value) {
    is.na(value) && !is.nan(value)
  }, logical(1))))

  undefined <- evaluate(NaN)
  expect_true(all(vapply(undefined, is.nan, logical(1))))

  positive_infinite <- evaluate(Inf)
  expect_true(all(vapply(positive_infinite, identical, logical(1), Inf)))

  for (invalid in list("invalid", c(0, 1), matrix(0, nrow = 1L))) {
    result <- evaluate(invalid)
    expect_true(all(vapply(result, identical, logical(1), NA_real_)))
  }
})

test_that("scalar joint-density failures identify the exact posterior row", {

  source_kind <- "likelihood"
  bad_value   <- 0

  testthat::local_mocked_bindings(
    .iwmde_log_q_grid_predictor_batch = function(...) NULL,
    .iwmde_log_q_grid_glmm_conditional_batch = function(...) NULL,
    .iwmde_log_q_grid_batch = function(...) NULL,
    .iwmde_replace_parameters = function(context, state, parameter, value,
                                         row, replacement) {

      parameters              <- state[["parameters"]]
      parameters[[parameter]] <- value
      list(
        valid                 = TRUE,
        row                   = row,
        parameters            = parameters,
        use_focal_prior_delta = state[["use_focal_prior_delta"]]
      )
    },
    .iwmde_log_likelihood_parameters = function(context, parameters,
                                                active_setup,
                                                likelihood_mode,
                                                row) {

      if (identical(source_kind, "likelihood") &&
          parameters[["marker"]] == 2 && parameters[["mu"]] == 5) {
        return(bad_value)
      }
      return(0)
    },
    .iwmde_replacement_log_prior_rows = function(samples, prior_list,
                                                  replacement) {

      if (identical(source_kind, "full prior") &&
          samples[1L, "marker"] == 2 && samples[1L, "mu"] == 5) {
        return(bad_value)
      }
      return(0)
    },
    .iwmde_focal_log_prior = function(prior, value, parameter = NULL) {

      if (identical(source_kind, "focal prior") &&
          prior[["marker"]] == 2 && value == 5) {
        return(bad_value)
      }
      return(0)
    }
  )

  normalization_values <- seq(-3, 3, length.out = 101)
  normalization_grid <- list(
    x            = normalization_values,
    z            = normalization_values,
    log_jacobian = rep(0, length(normalization_values))
  )
  make_states <- function(use_focal_prior_delta) {
    lapply(seq_len(2L), function(i) {
      list(
        row_index                = c(7L, 19L)[[i]],
        row                      = c(mu = 0, marker = i),
        parameters               = list(mu = 0, marker = i),
        active_setup             = list(),
        likelihood_mode          = "conditional",
        state_scope              = "local",
        prior_list               = list(),
        focal_prior              = list(marker = i),
        baseline_log_prior       = 0,
        baseline_focal_log_prior = 0,
        use_focal_prior_delta    = use_focal_prior_delta
      )
    })
  }
  evaluate <- function() {
    row_states <- make_states(identical(source_kind, "focal prior"))
    .iwmde_density_grid(
      context            = list(),
      parameter          = "mu",
      display_grid       = 5,
      normalization_grid = normalization_grid,
      transform          = .iwmde_parameter_transform(c(-Inf, Inf)),
      row_states         = row_states,
      active_mass        = 1,
      replacement        = list(type = "scalar"),
      estimator_rows     = c(7L, 19L),
      population_rows    = c(7L, 19L),
      chain_id           = c(1L, 1L)
    )
  }

  for (source in c("likelihood", "full prior", "focal prior")) {
    source_kind <- source
    for (invalid in list(NA_real_, NaN, Inf, "invalid", c(0, 1))) {
      bad_value <- invalid
      failure <- tryCatch(
        evaluate(),
        iwmde_construction_error = function(e) e
      )
      expect_s3_class(failure, "iwmde_construction_error")
      expect_equal(failure[["target"]], "mu")
      expect_equal(failure[["posterior_rows"]], 19L)
      expect_equal(failure[["stage"]], "joint-density grid evaluation")
    }

    bad_value <- -Inf
    density   <- evaluate()
    expect_true(is.finite(density[["y"]][[1L]]))
    expect_equal(density[["n_evaluated_rows"]], 2L)
  }
})

test_that("IWMDE active-branch wrappers call localized likelihood helpers", {

  data <- list(outcome = data.frame(yi = 0, sei = 1))
  attr(data, "effect_direction") <- "positive"

  context <- list(
    object = list(fit = list()),
    data   = data
  )
  active_setup <- list(priors = list())
  samples      <- matrix(0, nrow = 1, ncol = 0)

  calls <- list()

  testthat::local_mocked_bindings(
    .log_lik_posterior_setup = function(fit, posterior_samples,
                                         data, priors, unit,
                                          data_hash = NULL) {
      calls[["setup"]] <<- TRUE
      return(list(posterior_samples = posterior_samples))
    },
    .log_lik_from_posterior_samples_sum = function(fit, posterior_samples,
                                                    data, priors, unit,
                                                    data_hash = NULL) {
      calls[["sum"]] <<- TRUE
      return(rep(0, nrow(posterior_samples)))
    },
    .log_lik_from_evaluated_predictors_sum = function(fit, data, priors,
                                                      mu_samples,
                                                      tau_within_samples,
                                                      tau_between_samples = NULL,
                                                      posterior_samples = NULL,
                                                      unit = "estimate",
                                                      data_hash = NULL) {
      calls[["evaluated_sum"]] <<- TRUE
      return(rep(0, nrow(mu_samples)))
    },
    .selection_context_from_parts = function(fit, data, priors,
                                             posterior_samples,
                                             effect_direction,
                                             newdata = NULL) {
      calls[["selection"]] <<- TRUE
      return(list())
    }
  )

  expect_equal(
    .iwmde_log_lik_posterior_setup_active_branch(
      context           = context,
      posterior_samples = samples,
      active_setup      = active_setup,
      unit              = "estimate"
    ),
    list(posterior_samples = samples)
  )
  expect_equal(
    .iwmde_selection_context_active_branch(
      context           = context,
      active_setup      = active_setup,
      posterior_samples = samples
    ),
    list()
  )
  expect_equal(
    .iwmde_log_lik_from_posterior_samples_sum_active_branch(
      context           = context,
      posterior_samples = samples,
      active_setup      = active_setup,
      unit              = "estimate"
    ),
    0
  )
  expect_equal(
    .iwmde_log_lik_from_evaluated_predictors_sum_active_branch(
      context            = context,
      active_setup       = active_setup,
      mu_samples         = matrix(0, nrow = 1L, ncol = 1L),
      tau_within_samples = matrix(1, nrow = 1L, ncol = 1L),
      posterior_samples  = samples,
      unit               = "estimate"
    ),
    0
  )
  expect_equal(
    names(calls),
    c("setup", "selection", "sum", "evaluated_sum")
  )
})

test_that("IWMDE support transforms preserve density mass", {

  transforms <- list(
    .iwmde_parameter_transform(c(-Inf, Inf)),
    .iwmde_parameter_transform(c(0, Inf)),
    .iwmde_parameter_transform(c(-Inf, 1)),
    .iwmde_parameter_transform(c(0, 1))
  )

  for (transform in transforms) {
    z <- if (identical(transform[["type"]], "reverse_log")) {
      seq(5, -5, length.out = 5000)
    } else {
      seq(-5, 5, length.out = 5000)
    }
    x <- .iwmde_from_internal(z, transform)
    y <- stats::dnorm(z) / exp(.iwmde_log_jacobian(z, transform))

    expect_true(all(diff(x) > 0))
    expect_equal(.iwmde_to_internal(x, transform), z, tolerance = 1e-10)
    expect_equal(.iwmde_trapz(x, y), diff(stats::pnorm(c(-5, 5))), tolerance = 1e-6)
  }

  grid <- .iwmde_display_grid(
    xlim        = c(.001, 5),
    n_points    = 100,
    transform   = .iwmde_parameter_transform(c(0, Inf)),
    grid_method = "uniform"
  )
  expect_equal(
    diff(grid),
    rep(diff(range(grid)) / 99, 99),
    tolerance = 1e-12
  )
})


test_that("IWMDE support interiors respect narrow parameter scales", {

  bounded <- .iwmde_open_finite_support(
    xlim    = c(0, 1e-10),
    support = c(0, 1e-10)
  )
  asymmetric <- .iwmde_open_finite_support(
    xlim    = c(.9, 1),
    support = c(-1e100, 1)
  )
  subnormal <- .iwmde_open_finite_support(
    xlim    = c(0, 1e-320),
    support = c(0, Inf)
  )

  expect_true(0 < bounded[1L] && bounded[1L] < bounded[2L])
  expect_true(bounded[2L] < 1e-10)
  expect_true(.8 < asymmetric[1L] && asymmetric[1L] < asymmetric[2L])
  expect_true(asymmetric[2L] < 1)
  expect_true(0 < subnormal[1L] && subnormal[1L] < subnormal[2L])
})


test_that("IWMDE integration grid ignores display and evaluation ordinates", {

  values       <- seq(-1, 1, length.out = 101)
  display_grid <- seq(-1, 1, length.out = 21)
  transform    <- .iwmde_parameter_transform(c(-Inf, Inf))

  base <- .iwmde_normalization_grid(
    values               = values,
    display_grid         = display_grid,
    support              = c(-Inf, Inf),
    transform            = transform,
    normalization_points = 41,
    normalization_prob   = .90
  )
  with_far <- .iwmde_normalization_grid(
    values               = values,
    display_grid         = c(display_grid, 100),
    support              = c(-Inf, Inf),
    transform            = transform,
    normalization_points = 41,
    normalization_prob   = .90
  )

  expect_equal(with_far[["x"]], base[["x"]], tolerance = 1e-12)
  expect_equal(with_far[["z"]], base[["z"]], tolerance = 1e-12)
  expect_lt(max(with_far[["x"]]), 5)
})


test_that("IWMDE normalization grids reject saturated support transforms", {

  bounded_transform <- .iwmde_parameter_transform(c(0, 1))
  lower_saturation  <- c(2^-1074, 2^-1073, 2^-1072)
  upper_saturation  <- 1 - c(.Machine$double.eps, 2^-45, 2^-40)

  for (values in list(lower_saturation, upper_saturation)) {
    expect_null(.iwmde_normalization_grid(
      values               = values,
      display_grid         = values,
      support              = c(0, 1),
      transform            = bounded_transform,
      normalization_points = 41L,
      normalization_prob   = .90
    ))
  }

  cases <- list(
    list(values = seq(.1, .9, length.out = 20), support = c(0, 1)),
    list(values = exp(seq(-2, 2, length.out = 20)), support = c(0, Inf)),
    list(values = -exp(seq(-2, 2, length.out = 20)), support = c(-Inf, 0))
  )
  for (case in cases) {
    grid <- .iwmde_normalization_grid(
      values               = case[["values"]],
      display_grid         = case[["values"]],
      support              = case[["support"]],
      transform            = .iwmde_parameter_transform(case[["support"]]),
      normalization_points = 41L,
      normalization_prob   = .90
    )
    expect_false(is.null(grid))
    expect_true(all(diff(grid[["z"]]) > 0))
    expect_true(all(is.finite(grid[["log_jacobian"]])))
  }
})


test_that("IWMDE adaptive display grid concentrates support points", {

  set.seed(1)
  values <- c(
    stats::rnorm(800, mean = 0, sd = .30),
    stats::rnorm(200, mean = 2, sd = .12)
  )
  xlim <- c(-2, 3)

  adaptive <- .iwmde_display_grid(
    xlim        = xlim,
    n_points    = 60,
    transform   = .iwmde_parameter_transform(c(-Inf, Inf)),
    values      = values,
    grid_method = "adaptive"
  )
  uniform <- .iwmde_display_grid(
    xlim        = xlim,
    n_points    = 60,
    transform   = .iwmde_parameter_transform(c(-Inf, Inf)),
    values      = values,
    grid_method = "uniform"
  )

  expect_length(adaptive, 60L)
  expect_true(all(diff(adaptive) > 0))
  expect_equal(range(adaptive), xlim)
  expect_gt(stats::sd(diff(adaptive)), 0)
  expect_gt(sum(abs(adaptive) < .5), sum(abs(uniform) < .5))
  expect_gt(sum(abs(adaptive - 2) < .25), sum(abs(uniform - 2) < .25))
})

test_that("Chen conditioning rejects structural boundary mixtures", {

  .skip_if_missing_raw_fits("dat.lehmann2018_RoBMA")

  context <- .iwmde_context(load_fit("dat.lehmann2018_RoBMA", validate = FALSE))
  spec    <- .iwmde_parameter_spec(context, "mu")
  columns <- .iwmde_chen_conditioning_columns(
    context        = context,
    parameter      = "mu",
    parameter_spec = spec
  )

  expect_true(all(c(
    "tau", "omega[1]", "omega[2]", "omega[3]", "omega[4]",
    "omega[5]", "omega[6]", "PET", "PEESE"
  ) %in% columns))
  expect_gt(length(columns), 8L)

  component <- .iwmde_parameter_components(context, "mu", spec)
  values    <- .iwmde_parameter_values(context, "mu", spec)
  rows      <- which(component[["active"]] & is.finite(values))
  expect_error(
    .iwmde_chen_conditioning_matrix(
      context        = context,
      parameter      = "mu",
      parameter_spec = spec,
      active_rows    = head(rows, 50L),
      weight_rows    = rows
    ),
    "conditioning fit columns contain non-finite transformed values"
  )
})

test_that("IWMDE replacement maps preserve the baseline posterior ordinate", {

  fit_names <- c(
    "bcg_meta-analysis",
    "dat.lehmann2018-3PSM",
    "dat.lehmann2018_RoBMA",
    "bcg_glmm"
  )
  .skip_if_missing_raw_fits(fit_names)

  cases <- list(
    "bcg_meta-analysis"     = c("mu", "tau"),
    "dat.lehmann2018-3PSM"  = c("mu", "tau"),
    "dat.lehmann2018_RoBMA" = c("PET", "PEESE"),
    "bcg_glmm"              = c("pi[1]", "theta[1]")
  )

  for (name in names(cases)) {
    context <- .iwmde_context(load_fit(name, validate = FALSE))
    samples <- context[["posterior_samples"]]

    for (parameter in cases[[name]]) {
      component <- .iwmde_parameter_components(context, parameter)
      rows      <- which(component[["active"]] & is.finite(samples[, parameter]))
      rows      <- head(rows, 5L)

      row_states  <- .iwmde_row_states(context, rows, parameter)
      replacement <- .iwmde_replacement_spec(context, parameter)
      finite      <- vapply(row_states, function(state) {
        is.finite(state[["baseline_log_q"]])
      }, logical(1))
      rows        <- rows[finite]
      row_states  <- row_states[finite]
      expect_gt(length(rows), 0L)

      for (i in seq_along(rows)) {
        state <- row_states[[i]]
        log_q <- .iwmde_log_q_replacement(
          context     = context,
          state       = state,
          parameter   = parameter,
          value       = samples[rows[i], parameter],
          replacement = replacement
        )

        expect_equal(log_q, state[["baseline_log_q"]], tolerance = 1e-8)
      }
    }
  }
})

test_that("IWMDE scalar prior delta matches full prior ordinate", {

  .skip_if_missing_raw_fits("bcg_meta-analysis")

  context   <- .iwmde_context(load_fit("bcg_meta-analysis", validate = FALSE))
  samples   <- context[["posterior_samples"]]
  parameter <- "mu"
  rows      <- head(which(is.finite(samples[, parameter])), 3L)

  expect_gt(length(rows), 0L)

  values <- as.numeric(stats::quantile(
    samples[rows, parameter],
    probs = c(.25, .5, .75),
    names = FALSE,
    type  = 8
  ))
  values <- unique(values[is.finite(values)])

  expect_gt(length(values), 0L)

  for (row in rows) {
    state <- .iwmde_row_state(context, row, parameter)

    expect_true(state[["use_focal_prior_delta"]])

    focal_log_prior <- .iwmde_focal_log_prior_values(
      prior     = state[["focal_prior"]],
      values    = values,
      parameter = parameter
    )

    for (i in seq_along(values)) {
      candidate <- state[["row"]]
      candidate[[parameter]] <- values[[i]]
      delta_log_prior <- state[["baseline_log_prior"]] +
        focal_log_prior[[i]] - state[["baseline_focal_log_prior"]]
      full_log_prior <- .iwmde_log_prior_row(candidate, state[["prior_list"]])

      expect_equal(delta_log_prior, full_log_prior, tolerance = 1e-10)
    }
  }
})

test_that("IWMDE disables focal prior delta for vector priors", {

  .skip_if_missing_raw_fits("dat.lehmann2018_RoBMA_3lvl_mods_scale")

  context   <- .iwmde_context(load_fit("dat.lehmann2018_RoBMA_3lvl_mods_scale", validate = FALSE))
  samples   <- context[["posterior_samples"]]
  parameter <- "mu_Preregistered"
  rows      <- which(is.finite(samples[, parameter]) & samples[, parameter] != 0)
  rows      <- head(rows, 3L)

  expect_gt(length(rows), 0L)

  for (row in rows) {
    state <- .iwmde_row_state(context, row, parameter)

    expect_false(state[["use_focal_prior_delta"]])
  }
})

test_that("IWMDE conditions GLMM global parameters on sampled local states", {

  fit_names <- c("bcg_glmm", "nielweise2008_glmm")
  .skip_if_missing_raw_fits(fit_names)

  for (name in fit_names) {
    context <- .iwmde_context(load_fit(name, validate = FALSE))
    samples <- context[["posterior_samples"]]

    for (parameter in c("mu", "tau")) {
      rows        <- head(which(is.finite(samples[, parameter])), 2L)
      row_states  <- .iwmde_row_states(context, rows, parameter)
      replacement <- .iwmde_replacement_spec(context, parameter)

      expect_gt(length(rows), 0L)

      for (i in seq_along(rows)) {
        state <- row_states[[i]]
        expect_equal(state[["likelihood_mode"]], "conditional")
        expect_equal(state[["state_scope"]], "local")
        expect_true("theta" %in% names(state[["prior_list"]]))
        expect_true(any(c("pi", "phi") %in% names(state[["prior_list"]])))
        expect_true(is.finite(state[["baseline_log_lik"]]))

        log_q <- .iwmde_log_q_replacement(
          context     = context,
          state       = state,
          parameter   = parameter,
          value       = samples[rows[i], parameter],
          replacement = replacement
        )
        expect_equal(log_q, state[["baseline_log_q"]], tolerance = 1e-8)
      }
    }
  }
})

test_that("IWMDE uses marginal cluster likelihood for multilevel selection rows", {

  .skip_if_missing_raw_fits("dat.lehmann2018_RoBMA_3lvl_mods_scale")

  context   <- .iwmde_context(load_fit("dat.lehmann2018_RoBMA_3lvl_mods_scale", validate = FALSE))
  samples   <- context[["posterior_samples"]]
  parameter <- "mu_Preregistered"
  rows      <- which(is.finite(samples[, parameter]) & samples[, parameter] != 0)
  states    <- list()

  for (row in rows) {
    state <- .iwmde_row_state(context, row, parameter)
    if (isTRUE(state[["active_setup"]][["is_weightfunction"]])) {
      states[[length(states) + 1L]] <- state
    }
    if (length(states) >= 3L) {
      break
    }
  }

  expect_gt(length(states), 0L)

  replacement <- .iwmde_replacement_spec(context, parameter)
  for (state in states) {
    expect_equal(state[["likelihood_mode"]], "marginal")
    expect_true(is.finite(state[["baseline_log_lik"]]))

    log_q <- .iwmde_log_q_replacement(
      context     = context,
      state       = state,
      parameter   = parameter,
      value       = state[["row"]][[parameter]],
      replacement = replacement
    )
    expect_equal(log_q, state[["baseline_log_q"]], tolerance = 1e-8)
  }
})

test_that("IWMDE linear replacement maps preserve derived ordinates", {

  .skip_if_missing_raw_fits("bcg_meta-regression")

  context <- .iwmde_context(load_fit("bcg_meta-regression", validate = FALSE))
  samples <- context[["posterior_samples"]]
  spec    <- list(
    type    = "linear",
    weights = c(mu_intercept = 1, mu_ablat = 1)
  )
  parameter <- "mu_intercept + mu_ablat"
  values    <- .iwmde_parameter_values(
    context        = context,
    parameter      = parameter,
    parameter_spec = .iwmde_parameter_spec(context, parameter, spec)
  )
  component <- .iwmde_parameter_components(context, parameter, spec)
  rows      <- which(component[["active"]] & is.finite(values))
  rows      <- head(rows, 5L)

  row_states  <- .iwmde_row_states(context, rows, parameter, spec)
  replacement <- .iwmde_replacement_spec(context, parameter, spec)

  for (i in seq_along(rows)) {
    state <- row_states[[i]]
    log_q <- .iwmde_log_q_replacement(
      context     = context,
      state       = state,
      parameter   = parameter,
      value       = values[rows[i]],
      replacement = replacement
    )

    expect_equal(log_q, state[["baseline_log_q"]], tolerance = 1e-8)

    shifted <- .iwmde_replace_parameters(
      context     = context,
      state       = state,
      parameter   = parameter,
      value       = values[rows[i]] + .1,
      row         = samples[rows[i], ],
      replacement = replacement
    )
    expect_true(shifted[["valid"]])
    expect_equal(
      .iwmde_linear_value_row(context, shifted[["row"]], spec[["weights"]]),
      values[rows[i]] + .1,
      tolerance = 1e-10
    )
  }
})
