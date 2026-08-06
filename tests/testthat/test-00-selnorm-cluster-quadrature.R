.selnorm_cluster_quadrature_case <- function(direction = "positive",
                                             omega = NULL) {

  base_omega <- c(
    1,
    0.00033241074238402009,
    9.9224406360828319e-06,
    7.5355328624205985e-06,
    3.8179196456412814e-06,
    2.7288332895299089e-07
  )
  if (is.null(omega)) {
    omega <- base_omega
  }

  sign <- if (identical(direction, "positive")) 1 else -1
  yi   <- sign * 0.080077558738760934
  sei  <- 0.056401467112344769
  mu   <- sign * -0.054812987377585123
  steps <- c(.005, .025, .05, .10, .50)
  prior <- BayesTools::prior_weightfunction(
    side    = "one-sided",
    steps   = steps,
    weights = BayesTools::wf_fixed(omega)
  )
  context <- .selection_spec(
    priors           = list(outcome = list(bias = prior)),
    yi               = yi,
    sei              = sei,
    effect_direction = direction
  )
  context[["omega"]]      <- matrix(omega, nrow = 1L)
  context[["alpha"]]      <- 0
  context[["phack_kind"]] <- 0L
  context[["kernel_mode"]] <- SELKERNEL_STEP

  return(list(
    yi       = yi,
    sei      = sei,
    mu       = mu,
    steps    = steps,
    omega    = omega,
    sign     = sign,
    context  = context,
    setup    = list(
      S           = 1L,
      mu          = matrix(mu, nrow = 1L),
      tau_within  = matrix(0.0021728030140613199, nrow = 1L),
      tau_between = matrix(0.14061504344827008, nrow = 1L),
      cluster     = list(cluster_1 = 1L),
      weights     = 2.839324103295803
    )
  ))
}


.selnorm_cluster_interval_probability <- function(lower, upper, mean, sd) {

  lower_cdf <- stats::pnorm(lower, mean = mean, sd = sd)
  upper_cdf <- stats::pnorm(upper, mean = mean, sd = sd)
  lower_sf  <- stats::pnorm(lower, mean = mean, sd = sd, lower.tail = FALSE)
  upper_sf  <- stats::pnorm(upper, mean = mean, sd = sd, lower.tail = FALSE)

  return(pmax(upper_cdf - lower_cdf, lower_sf - upper_sf))
}


.selnorm_cluster_log_integrand <- function(gamma, case) {

  z_bound <- stats::qnorm(1 - case[["steps"]])
  lower   <- c(z_bound, -Inf)
  upper   <- c(Inf, z_bound)
  sd      <- sqrt(
    case[["setup"]][["tau_within"]][1, 1]^2 + case[["sei"]]^2
  )
  directed_y <- case[["sign"]] * case[["yi"]]

  vapply(gamma, function(gamma_value) {
    directed_mean <- case[["sign"]] * (
      case[["mu"]] + gamma_value *
        case[["setup"]][["tau_between"]][1, 1]
    )
    bin_probability <- .selnorm_cluster_interval_probability(
      lower = lower,
      upper = upper,
      mean  = directed_mean / case[["sei"]],
      sd    = sd / case[["sei"]]
    )
    selection_probability <- sum(case[["omega"]] * bin_probability)
    selected_log_density <- stats::dnorm(
      directed_y,
      mean = directed_mean,
      sd   = sd,
      log  = TRUE
    ) + log(case[["omega"]][4L]) - log(selection_probability)

    return(
      case[["setup"]][["weights"]] * selected_log_density +
        stats::dnorm(gamma_value, log = TRUE)
    )
  }, numeric(1L))
}


.selnorm_cluster_piecewise_reference <- function(case) {

  optimum <- stats::optimize(
    f       = function(gamma) .selnorm_cluster_log_integrand(gamma, case),
    interval = c(-12, 12),
    maximum  = TRUE,
    tol      = 1e-12
  )
  breakpoints <- seq(-12, 12, by = .25)
  pieces <- vapply(seq_len(length(breakpoints) - 1L), function(i) {
    stats::integrate(
      f = function(gamma) {
        exp(.selnorm_cluster_log_integrand(gamma, case) - optimum[["objective"]])
      },
      lower        = breakpoints[i],
      upper        = breakpoints[i + 1L],
      subdivisions = 200L,
      rel.tol      = 1e-10,
      abs.tol      = 0
    )[["value"]]
  }, numeric(1L))

  return(optimum[["objective"]] + log(sum(pieces)))
}


test_that("selected-normal cluster quadrature finds a displaced rare-selection mode", {

  skip_if_not(.has_native_norm_cluster_quadrature(selection = TRUE))
  skip_if_not(.has_native_selnorm_cluster_location_grid())

  case      <- .selnorm_cluster_quadrature_case()
  reference <- .selnorm_cluster_piecewise_reference(case)
  matrix_value <- .log_lik_cluster_norm_quadrature_native(
    setup             = case[["setup"]],
    yi                = case[["yi"]],
    sei               = case[["sei"]],
    is_weightfunction = TRUE,
    selection_context = case[["context"]]
  )
  row_sum_value <- .log_lik_cluster_norm_quadrature_sum(
    setup             = case[["setup"]],
    yi                = case[["yi"]],
    sei               = case[["sei"]],
    is_weightfunction = TRUE,
    selection_context = case[["context"]]
  )
  grid_value <- .log_lik_cluster_selnorm_location_grid(
    setup             = case[["setup"]],
    yi                = case[["yi"]],
    sei               = case[["sei"]],
    basis             = matrix(0, nrow = 1L),
    current           = 0,
    values            = 0,
    selection_context = case[["context"]]
  )

  expect_equal(matrix_value[1, 1], reference, tolerance = 1e-4)
  expect_equal(row_sum_value[1], reference, tolerance = 1e-4)
  expect_equal(grid_value[1, 1], reference, tolerance = 1e-4)
  expect_equal(attr(grid_value, "quadrature_order")[1, 1], 31L)
  expect_lt(attr(grid_value, "quadrature_relative_change")[1, 1], 5e-4)
})


test_that("selected-normal cluster quadrature mirrors the selection direction", {

  skip_if_not(.has_native_norm_cluster_quadrature(selection = TRUE))

  positive <- .selnorm_cluster_quadrature_case("positive")
  negative <- .selnorm_cluster_quadrature_case("negative")
  positive_value <- .log_lik_cluster_norm_quadrature_native(
    setup             = positive[["setup"]],
    yi                = positive[["yi"]],
    sei               = positive[["sei"]],
    is_weightfunction = TRUE,
    selection_context = positive[["context"]]
  )
  negative_value <- .log_lik_cluster_norm_quadrature_native(
    setup             = negative[["setup"]],
    yi                = negative[["yi"]],
    sei               = negative[["sei"]],
    is_weightfunction = TRUE,
    selection_context = negative[["context"]]
  )

  expect_equal(negative_value, positive_value, tolerance = 1e-12)
})


test_that("benign selected-normal clusters stop at the 15-node rule", {

  skip_if_not(.has_native_selnorm_cluster_location_grid())

  omega <- c(1, .65, .30)
  yi    <- .12
  sei   <- .20
  prior <- BayesTools::prior_weightfunction(
    side    = "one-sided",
    steps   = c(.025, .05),
    weights = BayesTools::wf_fixed(omega)
  )
  context <- .selection_spec(
    priors           = list(outcome = list(bias = prior)),
    yi               = yi,
    sei              = sei,
    effect_direction = "positive"
  )
  context[["omega"]]       <- matrix(omega, nrow = 1L)
  context[["alpha"]]       <- 0
  context[["phack_kind"]]  <- 0L
  context[["kernel_mode"]] <- SELKERNEL_STEP
  setup <- list(
    S           = 1L,
    mu          = matrix(.05, nrow = 1L),
    tau_within  = matrix(.09, nrow = 1L),
    tau_between = matrix(.11, nrow = 1L),
    cluster     = list(cluster_1 = 1L),
    weights     = NULL
  )
  out <- .log_lik_cluster_selnorm_location_grid(
    setup             = setup,
    yi                = yi,
    sei               = sei,
    basis             = matrix(0, nrow = 1L),
    current           = 0,
    values            = 0,
    selection_context = context
  )

  expect_true(is.finite(out[1, 1]))
  expect_equal(attr(out, "quadrature_order")[1, 1], 15L)
  expect_lt(attr(out, "quadrature_relative_change")[1, 1], 5e-4)
})


test_that("selected-normal cluster quadrature fails explicitly after order 31", {

  skip_if_not(.has_native_norm_cluster_quadrature(selection = TRUE))

  base_case <- .selnorm_cluster_quadrature_case()
  case <- .selnorm_cluster_quadrature_case(
    omega = base_case[["omega"]]^1.5
  )

  expect_error(
    .log_lik_cluster_norm_quadrature_native(
      setup             = case[["setup"]],
      yi                = case[["yi"]],
      sei               = case[["sei"]],
      is_weightfunction = TRUE,
      selection_context = case[["context"]]
    ),
    "Selected-normal cluster quadrature failed to converge.*after 31 nodes",
    fixed = FALSE
  )
})
