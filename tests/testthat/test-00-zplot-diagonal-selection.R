.diagonal_product_reference <- function(z, mean, integrated, retained, sei,
                                         omega = .3, probability = FALSE) {

  cutoff <- stats::qnorm(.025, lower.tail = FALSE) * sei
  total_sd <- sqrt(integrated + retained)
  if (integrated == 0) {
    if (probability) return(stats::pnorm(-z * sei, mean, total_sd) +
      stats::pnorm(z * sei, mean, total_sd, lower.tail = FALSE))
    return(sei * stats::dnorm(z * sei, mean, total_sd))
  }
  sd <- sqrt(integrated)
  A <- function(location) omega + (1 - omega) *
    stats::pnorm(cutoff, location, sd, lower.tail = FALSE)
  if (!probability) {
    y <- z * sei
    shifted <- retained / (integrated + retained) * (y - mean)
    conditional_sd <- sqrt(integrated * retained / (integrated + retained))
    inverse <- if (retained == 0) 1 / A(mean) else stats::integrate(function(u) {
      stats::dnorm(u) / A(mean + shifted + conditional_sd * u)
    }, -Inf, Inf, rel.tol = 1e-10)$value
    return(sei * stats::dnorm(y, mean, total_sd) *
      ifelse(y >= cutoff, 1, omega) * inverse)
  }
  conditional_tail <- function(location) {
    (omega * (stats::pnorm(-z * sei, location, sd) +
      stats::pnorm(z * sei, location, sd, lower.tail = FALSE)) +
      (1 - omega) * stats::pnorm(max(z * sei, cutoff), location, sd,
        lower.tail = FALSE)) / A(location)
  }
  if (retained == 0) return(conditional_tail(mean))
  stats::integrate(function(u) stats::dnorm(u) *
    conditional_tail(mean + sqrt(retained) * u), -Inf, Inf, rel.tol = 1e-10)$value
}

test_that("diagonal product zplots reuse scalar integration for every source axis", {

  S <- 3L
  K <- 4L
  samples <- matrix(seq_len(S), S, 1L, dimnames = list(NULL, "draw_id"))
  sei <- c(.1, .2, .15, .25)
  mean <- matrix(c(-.1, .3, .1, .2, -.1, .25, .05, .2, -.15, .3, .1, .05), S)
  estimate_sd <- c(.06, .05, 0)
  study_sd <- c(.12, 0, .08)
  control <- set_selection_likelihood_control(relative_tolerance = 1e-7)
  z <- c(-2, -.3, .8, stats::qnorm(.025, lower.tail = FALSE), 3)
  testthat::local_mocked_bindings(
    .brma_mv_heterogeneity_components = function(object, posterior_samples, ...) {
      rows <- as.integer(posterior_samples[, "draw_id"])
      sources <- .data_selection_model(object$data)$sources$random
      stats::setNames(lapply(sources, function(source) {
        value <- if (source$role == "estimate") estimate_sd else study_sd
        matrix(value[rows], length(rows), K)
      }), vapply(sources, `[[`, character(1L), "name"))
    },
    .predict_joint_selection_gaussian_parts = function(...) stop("Dense covariance cubes must not be requested."),
    .zplot_full_event_context_mixture = function(...) stop("Nested context projection must not be requested."),
    .package = "RoBMA"
  )
  make_case <- function(estimate, sampling, study, correlated = FALSE, omega = .3,
                       other = "condition") {
    dat <- data.frame(yi = 3 * sei, study = factor(study), esid = seq_len(K), paper = "p")
    V <- if (correlated) .4 * outer(sei, sei) else matrix(0, K, K)
    diag(V) <- sei^2
    prior <- BayesTools::prior_weightfunction("one-sided", .025,
      BayesTools::wf_fixed(c(1, omega)), model = selection_model(
        estimate_random_effects = estimate, other_random_effects = other,
        known_sampling_variance = sampling, group = paper))
    object <- bselmodel.mv(yi = yi, V = V, random = ~ 1 | study/esid, data = dat,
      prior_bias = prior, measure = "GEN", prior_unit_information_sd = 1,
      effect_direction = "positive", only_priors = TRUE, silent = TRUE)
    selection <- .selection_spec(list(outcome = list(bias = prior)), dat$yi, sei,
      "positive", signed_data = FALSE)
    selection$omega <- matrix(rep(c(1, omega), S), S, 2L, byrow = TRUE)
    selection$alpha <- numeric(S)
    selection$phack_kind <- integer(S)
    selection$kernel_mode <- rep(SELKERNEL_STEP, S)
    selection$vector_rule <- integer(S)
    selection$use_normal <- rep(FALSE, S)
    list(object = object, selection = selection)
  }
  predictive <- list(mu = mean, mu_extrapolated = mean + .05,
    tau_within = sqrt(matrix(estimate_sd^2 + study_sd^2, S, K)), sei = sei)
  for (estimate in c("integrate", "condition")) {
    for (sampling in c("integrate", "condition")) {
      integrated <- matrix(if (estimate == "integrate") estimate_sd^2 else 0, S, K)
      if (sampling == "integrate") integrated <- sweep(integrated, 2L, sei^2, "+")
      retained <- matrix(study_sd^2 + if (estimate == "condition") estimate_sd^2 else 0, S, K)
      if (sampling == "condition") retained <- sweep(retained, 2L, sei^2, "+")
      expected <- vapply(z, function(value) vapply(seq_len(S), function(row) {
        mean(vapply(seq_len(K), function(column) .diagonal_product_reference(
          value, mean[row, column], integrated[row, column], retained[row, column], sei[column]), numeric(1L)))
      }, numeric(1L)), numeric(S))
      first <- NULL
      for (study in list(rep("a", K), rep(c("a", "b"), each = 2L))) {
        case <- make_case(estimate, sampling, study)
        actual <- .zplot_joint_marginal(case$object, samples, predictive,
          case$selection, z, FALSE, control)
        expect_equal(actual$fitted, expected, tolerance = 2e-6)
        expect_equal(actual$extrapolated, .zplot_normal_density_matrix(z,
          predictive$mu_extrapolated, sqrt(integrated + retained), sei), tolerance = 1e-13)
        expect_identical(actual$weights, rep(1, S))
        expect_null(actual$EDR)
        if (is.null(first)) first <- actual else
          expect_equal(actual, first, tolerance = 1e-13)
      }
      expected_tail <- vapply(seq_len(S), function(row) {
        mean(vapply(seq_len(K), function(column) .diagonal_product_reference(
          1.5, mean[row, column], integrated[row, column], retained[row, column],
          sei[column], probability = TRUE), numeric(1L)))
      }, numeric(1L))
      tail <- .zplot_joint_marginal(case$object, samples, predictive,
        case$selection, 1.5, TRUE, control)
      expect_equal(as.numeric(tail$fitted), expected_tail, tolerance = 2e-6)
      expect_equal(tail$EDR, as.numeric(tail$extrapolated), tolerance = 0)
      if (sampling == "condition") {
        # Correlated retained sampling errors have the same scalar marginals.
        correlated <- make_case(estimate, sampling, rep("a", K), correlated = TRUE)
        actual <- .zplot_joint_marginal(correlated$object, samples, predictive,
          correlated$selection, z, FALSE, control)
        expect_equal(actual, first, tolerance = 1e-13)
      }
    }
  }
  # Dependence among integrated sources must retain the full-event route.
  dependent <- make_case("integrate", "integrate", rep("a", K), correlated = TRUE)
  expect_null(.zplot_diagonal_product_marginal(dependent$object, samples, predictive,
    dependent$selection, z, FALSE, control))
  dependent <- make_case("integrate", "integrate", rep("a", K), other = "integrate")
  expect_null(.zplot_diagonal_product_marginal(dependent$object, samples, predictive,
    dependent$selection, z, FALSE, control))
})
