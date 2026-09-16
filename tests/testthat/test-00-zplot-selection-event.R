.zplot_event_context <- function(rule, two_sided = FALSE, weights = c(1, .2), S = 1L,
                                 sei = c(.8, 1.2)) {

  cutoff <- if (two_sided) 1 else 0
  prior <- BayesTools::prior_weightfunction(
    if (two_sided) "two-sided" else "one-sided",
    if (two_sided) 2 * stats::pnorm(-cutoff) else .5,
    BayesTools::wf_fixed(weights),
    model = BayesTools::selection_model(weight_rule = rule, group = "paper")
  )
  context <- .selection_spec(list(outcome = list(bias = prior)), rep(0, length(sei)), sei, "positive")
  context$omega <- matrix(context$fixed_omega[1L, ], S, length(context$p_cuts) - 1L, byrow = TRUE)
  context$alpha <- numeric(S)
  context$phack_kind <- integer(S)
  context$kernel_mode <- rep(SELKERNEL_STEP, S)
  context$vector_rule <- rep(if (rule == "product") 0L else if (two_sided) 2L else 1L, S)
  context$use_normal <- rep(FALSE, S)
  context
}


test_that("full-event z projections retain product and best normalization", {

  sei <- c(.8, 1.2)
  mean <- c(.3, -.2)
  sd <- c(.9, 1.4)
  covariance <- matrix(c(sd[1]^2, 0, sd[2]^2), 1L)
  z <- c(-2, -.3, .4, 2)
  plan <- c(set_selection_likelihood_control(), list(designs = list()))
  for (two_sided in c(FALSE, TRUE)) for (rule in c("product", "best")) {
    for (weights in list(c(1, .2), c(1, 2), c(1, 0))) {
      context <- .zplot_event_context(rule, two_sided, weights)
      nonsignificant <- if (two_sided) {
        stats::pnorm(sei, mean, sd) - stats::pnorm(-sei, mean, sd)
      } else stats::pnorm(0, mean, sd)
      A <- weights[1] + diff(weights) * nonsignificant
      joint <- weights[1] + diff(weights) * prod(nonsignificant)
      expected <- vapply(z, function(value) {
        selected <- if (two_sided) abs(value) > 1 else value > 0
        w <- if (selected) weights[1] else weights[2]
        conditional <- if (selected) rep(weights[1], 2) else
          weights[1] + diff(weights) * rev(nonsignificant)
        mean(sei * stats::dnorm(value * sei, mean, sd) *
          if (rule == "product") w / A else conditional / joint)
      }, numeric(1L))
      got <- .zplot_full_event_projection(z, matrix(mean, 1L), covariance,
        sei, context, FALSE, plan)
      expect_equal(as.vector(got), expected, tolerance = 1e-12)
      q <- 2 * sei
      expected_tail <- mean((
        stats::pnorm(q, mean, sd, lower.tail = FALSE) * weights[1] +
          stats::pnorm(-q, mean, sd) * if (two_sided) weights[1] else
            if (rule == "product") weights[2] else
              weights[1] + diff(weights) * rev(nonsignificant)
      ) / if (rule == "product") A else joint)
      got_tail <- .zplot_full_event_projection(2, matrix(mean, 1L), covariance,
        sei, context, TRUE, plan)
      expect_equal(as.numeric(got_tail), expected_tail, tolerance = 1e-12)
      context$sign <- -1L
      reversed <- .zplot_full_event_projection(-z, matrix(-mean, 1L), covariance,
        sei, context, FALSE, plan)
      expect_equal(reversed, got, tolerance = 1e-12)
    }
  }
})


test_that("fresh retained contexts are integrated outside best normalization", {

  mean <- c(.3, -.2)
  sd <- c(.9, 1.4)
  sei <- c(.8, 1.2)
  loading <- c(.45, .2)
  z <- c(-1.5, .3, 2)
  context <- .zplot_event_context("best")
  control <- set_selection_likelihood_control(points_per_scramble = 256L,
    max_points_per_scramble = 4096L, relative_tolerance = .004)
  plan <- c(control, list(designs = list(), row_blocks = list(1:2)))
  expected <- vapply(z, function(value) {
    stats::integrate(function(h) {
      vapply(h, function(draw) {
        mu <- mean + loading * draw
        below <- stats::pnorm(0, mu, sd)
        conditional <- if (value > 0) c(1, 1) else 1 - .8 * rev(below)
        mean(sei * stats::dnorm(value * sei, mu, sd) * conditional) /
          (1 - .8 * prod(below)) * stats::dnorm(draw)
      }, numeric(1L))
    }, -Inf, Inf, rel.tol = 1e-9)$value
  }, numeric(1L))
  set.seed(167)
  before <- .Random.seed
  got <- .zplot_full_event_context_mixture(z, matrix(mean, 1L),
    array(diag(sd^2), c(1L, 2L, 2L)),
    array(tcrossprod(loading), c(1L, 2L, 2L)), sei, context, FALSE, control, plan)
  expect_equal(as.vector(got), expected, tolerance = .001)
  expect_identical(.Random.seed, before)
  singleton <- .zplot_full_event_context_mixture(.3, matrix(mean[1L]),
    array(sd[1L]^2, c(1L, 1L, 1L)), array(0, c(1L, 1L, 1L)),
    sei[1L], context, FALSE, control, c(control, list(row_blocks = list(1L))))
  expect_identical(dim(singleton), c(1L, 1L))
})


test_that("diagonal best bin reuse matches full events at physical boundaries", {

  plan <- c(set_selection_likelihood_control(), list(designs = list()))
  mean <- matrix(c(.2, -.1, .3, -.4, .1, .2), 3L, 2L)
  covariance <- matrix(c(.7, 0, 1.3), 1L)
  sei <- c(.8, 1.2)
  for (side in c(FALSE, TRUE)) for (direction in c(-1L, 1L)) {
    for (weights in list(c(1, .2), c(1, 0), c(1, 2))) {
      selection <- .zplot_event_context("best", side, weights, 3L)
      selection$kernel_mode[2L] <- SELKERNEL_NORMAL
      selection$use_normal[2L] <- TRUE
      selection$sign <- direction
      z <- sort(unique(c(-3, -.2, .3, 3,
        direction * selection$z_lower[is.finite(selection$z_lower)])))
      reference <- .zplot_full_event_projection(z, mean, covariance, sei,
        selection, FALSE, plan, diagonal_best = FALSE)
      actual <- .zplot_full_event_projection(z, mean, covariance, sei,
        selection, FALSE, plan)
      expect_equal(actual, reference, tolerance = 1e-12)
      expect_equal(.zplot_full_event_projection(1, mean, covariance, sei,
        selection, TRUE, plan), .zplot_full_event_projection(1, mean, covariance,
          sei, selection, TRUE, plan, diagonal_best = FALSE), tolerance = 1e-12)
    }
  }
  selection <- .zplot_event_context("best", weights = c(1, 0))
  rare <- .zplot_full_event_projection(-40, matrix(c(-40, -40), 1L),
    matrix(c(1, 0, 1), 1L), c(1, 1), selection, FALSE, plan)
  # Either independent far-tail event can select the pair; near the other
  # component's mode the projected density retains half its Gaussian mass.
  expect_equal(as.numeric(rare), stats::dnorm(0) / 2, tolerance = 1e-12)
})


test_that("event refinement preserves accepted rows and exact zero numerators", {

  plan <- list(points_per_scramble = 4L, max_points_per_scramble = 8L,
    relative_tolerance = .01, scrambles = 8L, seed = 371L,
    designs = list(`2` = "old dimension design", quadrature = "keep"))
  calls <- list()
  initial <- list(log_numerator = c(-Inf, -.3), relative_mcse = c(0, .02))
  result <- .selection_joint_checked_event(function(current, rows) {
    calls[[length(calls) + 1L]] <<- list(plan = current, rows = rows)
    if (is.null(rows)) initial else list(log_numerator = -.4, relative_mcse = .002)
  }, plan, "Test event", "Increase the integration budget.")
  expect_identical(result, list(log_numerator = c(-Inf, -.4), relative_mcse = c(0, .002)))
  expect_identical(calls[[2L]]$rows, 2L)
  expect_identical(calls[[2L]]$plan$designs, list(quadrature = "keep"))
  expect_identical(calls[[2L]]$plan[c("seed", "scrambles")], plan[c("seed", "scrambles")])
  expect_identical(.selection_joint_checked_event(function(plan, rows) initial,
    modifyList(plan, list(relative_tolerance = .03)), "Test", "Remedy."), initial)
  expect_error(.selection_joint_checked_event(function(plan, rows) {
    list(log_mass = 0, relative_mcse = .025)
  }, plan, "Test event", "Increase the integration budget."),
    "Test event was rejected by diagnostics: relative integration MCSE was 0.025. Increase the integration budget.",
    fixed = TRUE)
})


test_that("declared rank-one sampling events retain their scalar projections", {

  sampling <- list(rank = 1L, diagonal = c(0, 0), loading = matrix(c(1, -2), 2L, 1L))
  expect_identical(.selection_joint_declared_rank_one_loading(sampling, matrix(0, 2L, 2L)),
    c(1, -2))
  expect_null(.selection_joint_declared_rank_one_loading(sampling, diag(c(1e-12, 0))))
  changed <- sampling
  changed$diagonal <- c(1e-12, 0)
  expect_null(.selection_joint_declared_rank_one_loading(changed, matrix(0, 2L, 2L)))
  control <- set_selection_likelihood_control()
  plan <- c(control, list(row_blocks = list(1:2), designs = list()))
  context <- .zplot_event_context("best")
  z <- c(-2, -.3, .4, 2)
  actual <- .zplot_full_event_context_mixture(z, matrix(c(0, 0), 1L),
    array(tcrossprod(sampling$loading), c(1L, 2L, 2L)),
    array(0, c(1L, 2L, 2L)), c(1, 2), context, FALSE, control, plan,
    sampling_factor_blocks = list(sampling), random_covariance = array(0, c(1L, 2L, 2L)))
  # Y=(T,-2T), so one p-value is always in the preferred halfspace.
  expect_equal(as.vector(actual), stats::dnorm(z), tolerance = 1e-12)
  random_source <- .zplot_full_event_context_mixture(z, matrix(c(0, 0), 1L),
    array(tcrossprod(sampling$loading), c(1L, 2L, 2L)),
    array(0, c(1L, 2L, 2L)), c(1, 2), context, FALSE, control, plan,
    block_factors = list(list(rank = 1L, residual_sd = matrix(0, 1L, 2L),
      loading = matrix(c(1, -2), 1L))),
    random_covariance = array(tcrossprod(sampling$loading), c(1L, 2L, 2L)))
  expect_equal(random_source, actual, tolerance = 1e-12)
  expect_error(.zplot_full_event_projection(.3, matrix(c(0, 0), 1L), NULL,
    c(1, 2), context, FALSE, plan, rank_one_loading = matrix(c(0, 0), 1L)),
    "Selected partial-vector density is unavailable because the observed rank-one law has no Lebesgue density.",
    fixed = TRUE)
})


test_that("whole-sampling conditional zplots preserve scalar and deterministic marginal laws", {

  sei <- c(.2, .3)
  mu <- .2
  tau <- .3
  z <- c(-1, .5, 2)
  for (estimate in c("condition", "integrate")) {
    object <- bselmodel(yi = c(.2, -.1), sei = sei, measure = "GEN",
      prior_unit_information_sd = 1, only_priors = TRUE,
      prior_heterogeneity = BayesTools::prior("point", list(location = tau)),
      prior_bias = BayesTools::prior_weightfunction("one-sided", .5,
        BayesTools::wf_fixed(c(1, .4)), model = BayesTools::selection_model(
          estimate_random_effects = estimate, known_sampling_variance = "condition")))
    samples <- cbind(mu = mu, tau = tau)
    density <- .zplot_selection_marginal(object, samples, z, NULL, "marginal",
      set_selection_likelihood_control())
    tail <- .zplot_selection_marginal(object, samples, NULL, 2, "marginal",
      set_selection_likelihood_control())
    if (estimate == "condition") {
      expected_density <- vapply(z, function(value) {
        mean(sei * stats::dnorm(value * sei, mu, sqrt(sei^2 + tau^2)))
      }, numeric(1L))
      expected_tail <- mean(stats::pnorm(-2 * sei, mu, sqrt(sei^2 + tau^2)) +
        stats::pnorm(2 * sei, mu, sqrt(sei^2 + tau^2), lower.tail = FALSE))
    } else {
      # Independent one-dimensional integration over retained sampling errors.
      expected_density <- vapply(z, function(value) mean(vapply(sei, function(s) {
        q <- value * s
        variance <- s^2 + tau^2
        conditional <- stats::integrate(function(e) {
          stats::dnorm(e, s^2 / variance * (q - mu), sqrt(s^2 * tau^2 / variance)) /
            (.4 + .6 * stats::pnorm((mu + e) / tau))
        }, -Inf, Inf, rel.tol = 1e-10)$value
        s * stats::dnorm(q, mu, sqrt(variance)) * ifelse(value > 0, 1, .4) * conditional
      }, numeric(1L))), numeric(1L))
      expected_tail <- mean(vapply(sei, function(s) {
        stats::integrate(function(e) {
          stats::dnorm(e, 0, s) * (.4 * stats::pnorm(-2 * s, mu + e, tau) +
            stats::pnorm(2 * s, mu + e, tau, lower.tail = FALSE)) /
            (.4 + .6 * stats::pnorm((mu + e) / tau))
        }, -Inf, Inf, rel.tol = 1e-10)$value
      }, numeric(1L)))
    }
    expect_equal(as.numeric(density$fitted), expected_density, tolerance = 1e-6)
    expect_equal(as.numeric(tail$fitted), expected_tail, tolerance = 1e-6)
    expect_true(.zplot_vector_selection_target(object))
  }
})


test_that("rank-deficient factor projections retain scalar Gaussian densities", {

  control <- set_selection_likelihood_control()
  plan <- c(control, list(
    statistical_target = "whole_sampling_error_selection", row_blocks = list(1:3),
    quadrature = .selection_joint_cluster_quadrature_rules(SELNORM_CLUSTER_QUADRATURE_ORDERS),
    factor_quadrature = .selection_joint_factor_quadrature_rules()
  ))
  qmc <- BayesTools::selection_qmc_design(dimensions = 6L,
    points = control$max_points_per_scramble, scrambles = control$scrambles, seed = control$seed)
  loading <- rbind(c(1, 0), c(0, 1), c(1, 1))
  sei <- c(1, 1, sqrt(2))
  context <- .zplot_event_context("product", weights = c(1, .4), sei = sei)
  z <- c(-1, .5, 2)
  omega <- .4
  # The four quadrants of two independent standard normals give this exact
  # normalizer for W(Z1) W(Z2) W(Z1+Z2).
  normalizer <- (1 + omega + omega^2 + omega^3) / 4
  reference <- function(value) {
    first_mass <- if (value > 0) {
      .5 + omega * (.5 - stats::pnorm(-value)) + omega^2 * stats::pnorm(-value)
    } else {
      stats::pnorm(value) + omega * (.5 - stats::pnorm(value)) + .5 * omega^2
    }
    # Given Z1+Z2=q, Z1 is Normal(q/2,1/2), independently of its sum.
    third_mass <- if (value > 0) {
      omega + (1 - omega) * (2 * stats::pnorm(value) - 1)
    } else {
      omega + (omega^2 - omega) * (2 * stats::pnorm(-value) - 1)
    }
    stats::dnorm(value) * ifelse(value > 0, 1, omega) *
      (2 * first_mass + third_mass) / (3 * normalizer)
  }
  density <- .selection_factor_projection(matrix(0, 1L, 3L), matrix(0, 1L, 3L),
    matrix(loading, 1L), sei, context, plan, z, FALSE, rep(1L, 3L), qmc)
  expected <- vapply(z, reference, numeric(1L))
  expect_equal(as.numeric(density$density), expected, tolerance = .001)
  expect_equal(as.numeric(density$log_density), -log(normalizer), tolerance = .001)
  tail <- .selection_factor_projection(matrix(0, 1L, 3L), matrix(0, 1L, 3L),
    matrix(loading, 1L), sei, context, plan, 2, TRUE, rep(1L, 3L), qmc)
  expected_tail <- stats::integrate(function(x) vapply(x, reference, numeric(1L)),
    -Inf, -2, rel.tol = 1e-11)$value +
    stats::integrate(function(x) vapply(x, reference, numeric(1L)),
      2, Inf, rel.tol = 1e-11)$value
  expect_equal(as.numeric(tail$density), expected_tail, tolerance = .001)
  context$sign <- -1L
  reversed <- .selection_factor_projection(matrix(0, 1L, 3L), matrix(0, 1L, 3L),
    matrix(loading, 1L), sei, context, plan, -z, FALSE, rep(1L, 3L), qmc)
  reversed_tail <- .selection_factor_projection(matrix(0, 1L, 3L), matrix(0, 1L, 3L),
    matrix(loading, 1L), sei, context, plan, 2, TRUE, rep(1L, 3L), qmc)
  expect_equal(reversed$density, density$density, tolerance = 1e-10)
  expect_equal(reversed_tail$density, tail$density, tolerance = 1e-10)

  # A zero candidate coordinate is a retained context coordinate. Its marginal
  # remains Gaussian; native point-density projection must not approximate it.
  mean <- matrix(c(.2, -.1, .3), 1L)
  loading <- matrix(c(.3, 0, -.3), 3L, 1L)
  latent <- diag(c(.25, .36, .49))
  context <- .zplot_event_context("product", weights = c(1, 1), sei = rep(1, 3L))
  factors <- list(list(rank = 1L, diagonal = matrix(0, 1L, 3L),
    residual_sd = matrix(0, 1L, 3L), loading = matrix(loading, 1L)))
  actual <- .zplot_full_event_context_mixture(z, mean,
    array(tcrossprod(loading), c(1L, 3L, 3L)), array(latent, c(1L, 3L, 3L)),
    rep(1, 3L), context, FALSE, control, plan, block_factors = factors)
  expected <- vapply(z, function(value) {
    mean(stats::dnorm(value, as.numeric(mean), sqrt(diag(latent) + rowSums(loading^2))))
  }, numeric(1L))
  expect_equal(as.numeric(actual), expected, tolerance = .001)
})


test_that("factor projections report failed numerical diagnostics completely", {

  native_result <- list(log_density = 0, relative_mcse = .02, density = matrix(.2))
  project <- .selection_factor_projection
  environment(project) <- list2env(list(.Call = function(...) native_result),
    parent = environment(project))
  evaluate <- function() {
    project(matrix(0, 1L, 2L), matrix(1, 1L, 2L), matrix(0, 1L, 0L),
      c(.8, 1.2), .zplot_event_context("product"), set_selection_likelihood_control(),
      1, FALSE, c(1L, 1L), numeric())
  }
  error <- tryCatch(evaluate(), error = identity)
  expect_identical(conditionMessage(error), paste0(
    "Zplot selection projection was rejected by diagnostics: relative integration error was 0.02. ",
    "Increase 'max_points_per_scramble' or 'scrambles' in ",
    "'integration_control = set_selection_likelihood_control()'."))
  expect_null(conditionCall(error))
  native_result$relative_mcse <- 0
  native_result$density[,] <- NaN
  error <- tryCatch(evaluate(), error = identity)
  expect_identical(conditionMessage(error), "The selected factor projection returned invalid output.")
  expect_null(conditionCall(error))
})


test_that("vector zplot summaries expose the normalized reference without missing counts", {

  testthat::local_mocked_bindings(
    .zplot_fun.brma = function(...) list(EDR = c(.2, .3), weights = c(3, 4)),
    .package = "RoBMA"
  )
  yi <- rep(c(0, 3), 10L)
  vector <- bselmodel.mv(yi = yi, V = diag(20L), random = NULL,
    measure = "GEN", prior_unit_information_sd = 1, only_priors = TRUE)
  vector <- as_zplot(vector)
  expect_null(vector$zplot$estimates$weights)
  expect_identical(vector$zplot$target$reference, "pre_selection_gaussian")
  expect_identical(vector$zplot$target$missing_count, "unavailable")
  expect_identical(vector$zplot$target$publication_groups,
    .data_selection_model(vector$data)$groups)
  table <- as.data.frame(vector)
  expect_identical(table, data.frame(vector))
  expect_identical(table$parameter, c("EDR", "Soric FDR"))
  expect_identical(table$component, c("zplot", "zplot"))
  expect_true(any(grepl("^CI_", names(table))))
  printed <- paste(capture.output(print(summary(vector))), collapse = "\n")
  expect_match(printed,
    "Missing N is unavailable because relative selection weights do not identify an absolute publication probability.",
    fixed = TRUE)
  univariate <- bselmodel(yi = yi, sei = rep(1, 20L), measure = "GEN",
    prior_unit_information_sd = 1, only_priors = TRUE)
  univariate <- as_zplot(univariate)
  expect_identical(univariate$zplot$estimates$weights, c(3, 4))
  expect_true("Missing N" %in% as.data.frame(univariate)$parameter)
})

test_that("zplot reference scope follows contextual sources and resolved product models", {

  testthat::local_mocked_bindings(
    .zplot_fun.brma = function(...) list(EDR = c(.2, .3), weights = c(3, 4)),
    .package = "RoBMA"
  )
  dat <- data.frame(yi = rep(c(0, 3), 10L), sei = 1,
    study = rep(seq_len(10L), each = 2L), estimate = seq_len(20L))
  implicit <- bselmodel(yi = yi, sei = sei, cluster = study, data = dat,
    selection = BayesTools::selection_model(other_random_effects = "integrate", known_sampling_variance = "integrate"),
    measure = "GEN", prior_unit_information_sd = 1, only_priors = TRUE)
  explicit <- bselmodel(yi = yi, sei = sei, cluster = study, data = dat,
    prior_bias = prior_weightfunction("one-sided", steps = .025,
      model = selection_model(other_random_effects = "integrate",
        known_sampling_variance = "integrate", group = study)),
    measure = "GEN", prior_unit_information_sd = 1, only_priors = TRUE)
  expect_identical(.data_selection_model(implicit$data)$groups$group_index,
    .data_selection_model(explicit$data)$groups$group_index)
  expect_true(.data_selection_model(implicit$data)$applicability$other_random_effects)
  for (object in list(implicit, explicit)) {
    z <- as_zplot(object)
    expect_identical(z$zplot$target$reference, "pre_selection_gaussian")
    expect_identical(z$zplot$target$missing_count, "unavailable")
    expect_null(z$zplot$estimates$weights)
  }

  independent <- bselmodel(yi = yi, sei = sei, data = dat,
    selection = BayesTools::selection_model(other_random_effects = "integrate", known_sampling_variance = "integrate"),
    measure = "GEN", prior_unit_information_sd = 1, only_priors = TRUE)
  independent_explicit <- bselmodel(yi = yi, sei = sei, data = dat,
    prior_bias = prior_weightfunction("one-sided", steps = .025,
      model = selection_model(other_random_effects = "integrate",
        known_sampling_variance = "integrate", group = estimate)),
    measure = "GEN", prior_unit_information_sd = 1, only_priors = TRUE)
  conditional <- bselmodel(yi = yi, sei = sei, cluster = study, data = dat,
    measure = "GEN", prior_unit_information_sd = 1, only_priors = TRUE)
  conditional_explicit <- bselmodel(yi = yi, sei = sei, cluster = study, data = dat,
    prior_bias = prior_weightfunction("one-sided", steps = .025,
      model = selection_model(group = study)),
    measure = "GEN", prior_unit_information_sd = 1, only_priors = TRUE)
  expect_identical(.data_selection_model(independent$data)$applicability,
    list(estimate_random_effects = TRUE, other_random_effects = FALSE, known_sampling_variance = TRUE))
  expect_identical(.data_selection_model(independent$data)$groups$group_index,
    .data_selection_model(independent_explicit$data)$groups$group_index)
  expect_identical(.data_selection_model(conditional$data)$groups$group_index,
    .data_selection_model(conditional_explicit$data)$groups$group_index)
  for (object in list(independent, independent_explicit, conditional, conditional_explicit)) {
    z <- as_zplot(object)
    expect_identical(z$zplot$target$reference, "univariate_extrapolation")
    expect_identical(z$zplot$target$missing_count, "univariate_convention")
    expect_identical(z$zplot$estimates$weights, c(3, 4))
  }
})
