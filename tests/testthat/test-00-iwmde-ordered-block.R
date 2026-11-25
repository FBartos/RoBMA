.ordered_block_fixture <- function(prior = BayesTools::prior("normal", list(0, 1)),
                                   alpha = 3, beta = 1) {

  counts <- new.env(parent = emptyenv())
  for (name in c("inner_calls", "outer_nodes", "candidate_rows", "native_nodes", "maximum_chunk")) counts[[name]] <- 0L
  counts[["max_log_scale"]] <- 1
  share_prior <- BayesTools::prior("beta", list(alpha, beta))
  list(prior = prior, beta = share_prior, selected = FALSE, counts = counts,
    normal = .iwmde_qcmde_normal_prior_parameters(prior),
    normal_probability = .iwmde_ordered_probability_prior(prior),
    beta_probability = .iwmde_ordered_probability_prior(share_prior),
    gaussian = list(A = c(2, 0, 10), H = c(.3, .8, 0), C = c(-.7, 0, 0)),
    mu_a = 1, mu_b = 1, descriptor = list(k = -2.5), closed_total_numerator = TRUE,
    log_correction_lower = 0, log_correction_upper = 0)
}


test_that("direct Gaussian completion agrees with physical residuals near share boundaries", {

  shares <- c(2^-40, 1 - 2^-40, 1)
  for (i in seq_along(shares)) {
    share <- shares[[i]]
    node_reference <- .ordered_gaussian_reference(paste0("boundary-", i))
    expect_identical(share, node_reference[["share"]])
    block <- .ordered_block_fixture(BayesTools::prior("normal", list(0, 2^40)))
    block[["mu0"]] <- block[["mu_c"]] <- c(0, 0, 0)
    block[["mu_a"]] <- rep(1, 3L)
    block[["mu_b"]] <- rep(-1, 3L)
    block[["setup"]] <- list(yi = rep(1, 3L))
    block[["variance"]] <- rep(1e-8, 3L)
    block[["weights"]] <- c(1, 2, .5)
    information <- sum(block[["weights"]]) / 1e-8
    block[["gaussian"]] <- list(A = information * c(1, -2, 1), H = information * c(1, -1),
      C = sum(block[["weights"]] * stats::dnorm(1, 0, 1e-4, log = TRUE)))
    b <- 1 - share
    precision <- sum(block[["weights"]]) * b^2 / 1e-8 + 2^-80
    center <- sum(block[["weights"]]) * b / 1e-8 / precision
    sd <- 1 / sqrt(precision)
    original <- function(total) {
      vapply(total, function(t) sum(block[["weights"]] * stats::dnorm(1, b * t, 1e-4, log = TRUE)), numeric(1))
    }
    actual <- .iwmde_ordered_normal_completion(block, share)
    expected <- original(center) - .5 * (center / 2^40)^2 - .5 * log(2^80 * precision)
    floor <- 64 * .Machine$double.eps * 3 * max(1, abs(expected), abs(actual[["log_mass"]]))
    expect_lte(abs(actual[["log_mass"]] - expected), floor)
    expect_lte(abs(actual[["sd"]] - sd), 64 * .Machine$double.eps * sd)
    .ordered_gaussian_reference_inputs(node_reference, block)
    grid <- node_reference[["total"]]
    expect_length(grid, 3L)
    evaluated <- .iwmde_ordered_block_gaussian_log(block, grid, share)
    replayed <- original(grid)
    node_allowance <- node_reference[["allowance"]]
    expect_true(all(abs(evaluated - node_reference[["declared"]]) <= node_allowance))
    replay_mean <- outer(grid, rep(b, 3L))
    replay_sd <- matrix(1e-4, 3L, 3L)
    if (.ordered_gaussian_replay_matches(node_reference, replay_mean, replay_sd)) {
      expect_true(all(abs(replayed - node_reference[["rounded_replay"]]) <= node_allowance))
    }
    rounding_bound <- .ordered_gaussian_rounding_bound(node_reference, block[["weights"]],
      block[["variance"]], replay_mean, replay_sd)
    expect_true(all(node_reference[["forward_bound"]] <= rounding_bound))
    expect_true(all(abs(replayed - node_reference[["declared"]]) <= rounding_bound + node_allowance))
    if (share < 1) {
      anchor <- original(center) + stats::dnorm(center, 0, 2^40, log = TRUE)
      reference <- stats::integrate(function(z) exp(original(center + sd * z) +
        stats::dnorm(center + sd * z, 0, 2^40, log = TRUE) - anchor), -Inf, Inf, rel.tol = 1e-10, abs.tol = 0)
      expect_lte(abs(actual[["log_mass"]] - anchor - log(sd) - log(reference[["value"]])),
        reference[["abs.error"]] / reference[["value"]] + floor)
    }
    changed <- block
    changed[["setup"]][["yi"]] <- rep(2, 3L)
    expect_equal(.iwmde_ordered_block_gaussian_log(changed, 0, share),
      sum(block[["weights"]] * stats::dnorm(2, 0, 1e-4, log = TRUE)), tolerance = 0)
  }
})


test_that("positive unrepresentable Gaussian precision is not treated as constant", {

  block <- .ordered_block_fixture()
  block[["setup"]] <- list(yi = 1)
  block[["variance"]] <- block[["weights"]] <- 1
  block[["mu0"]] <- block[["mu_c"]] <- block[["mu_b"]] <- 0
  for (basis in c(1e-200, 1e200)) {
    block[["mu_a"]] <- basis
    expect_error(.iwmde_ordered_normal_completion(block, .5), class = "RoBMA_ordered_gaussian_completion_unrepresentable")
  }
})


test_that("block Gaussian completion includes residual and original truncation masses", {

  block <- .ordered_block_fixture(BayesTools::prior("normal", list(.3, 1.2), truncation = list(-.5, 1.4)))
  for (share in c(0, .03, .5, 1)) {
    actual <- .iwmde_ordered_normal_completion(block, share)
    independent <- stats::integrate(function(total) exp(.iwmde_ordered_block_gaussian_log(block, total, share)) *
      stats::dnorm(total, .3, 1.2) / (stats::pnorm(1.4, .3, 1.2) - stats::pnorm(-.5, .3, 1.2)),
      -.5, 1.4, rel.tol = 1e-12, abs.tol = 0)
    expect_lte(abs(exp(actual[["log_mass"]]) - independent[["value"]]),
      independent[["abs.error"]] + 128 * .Machine$double.eps * independent[["value"]])
    expect_lt(actual[["original_prior_interval_log_mass"]], 0)
  }
})


test_that("block total supplier and positive denominator agree with physical integrals", {

  block <- .ordered_block_fixture()
  denominator <- .iwmde_ordered_block_denominator(NULL, block)
  numerator <- .iwmde_ordered_block_numerator(NULL, block, "total")
  independent <- stats::integrate(function(share) vapply(share, function(s) {
    stats::integrate(function(total) stats::dnorm(total) *
      exp(-.7 + (.3 + .8 * s) * total - .5 * (2 + 10 * s^2) * total^2),
      -Inf, Inf, rel.tol = 1e-10, abs.tol = 0)[["value"]] * stats::dbeta(s, 3, 1)
  }, numeric(1)), 0, 1, rel.tol = 1e-9, abs.tol = 0)
  expect_lte(abs(exp(denominator[["log_value"]]) - independent[["value"]]),
    independent[["abs.error"]] + exp(denominator[["log_value"]]) * denominator[["relative_error"]])
  expect_equal(exp(numerator[["log_value"]]), stats::dnorm(0) * 1.5 * exp(-.7) / 2.5, tolerance = 1e-14)
  expect_identical(numerator[["relative_error"]], 0)
  expect_lte(denominator[["relative_error"]], .025 / 1.025)
})


test_that("proper CDF panels retain far-tail truncation and open quantile nodes", {

  for (prior in list(BayesTools::prior("normal", list(0, 1), truncation = list(8, 9)),
    BayesTools::prior("normal", list(0, 1), truncation = list(-9, -8)),
    BayesTools::prior("gamma", list(2, 1)), BayesTools::prior("beta", list(1, 2)))) {
    probability <- .iwmde_ordered_probability_prior(prior)
    cuts <- .iwmde_ordered_base_cuts(probability)
    panels <- .iwmde_ordered_probability_panels(cuts, probability)
    expect_equal(sum(exp(vapply(panels, `[[`, numeric(1), "log_mass"))), 1, tolerance = 1e-13)
    for (panel in panels) {
      nodes <- probability[["conditional"]](panel, .gauss_legendre_nodes(16)[["nodes"]])
      expect_true(all(nodes > panel[["lower"]] & nodes < panel[["upper"]]))
    }
  }
})


test_that("ordered root isolation preserves exact endpoint factors and close critical bands", {

  roots <- .iwmde_ordered_poly_roots(c(0, -.2, 1))
  expect_identical(as.numeric(roots[["roots"]][1, ]), c(0, 0))
  expect_true(roots[["roots"]][2, 1] <= .2 && roots[["roots"]][2, 2] >= .2)
  expect_error(.iwmde_ordered_poly_roots(c(1e-30, -2e-15, 1)),
    class = "RoBMA_ordered_root_isolation_unresolved")
})


test_that("Gamma share numerator uses the genuine inverse-weighted prior and signed Jacobian", {

  block <- .ordered_block_fixture(BayesTools::prior("gamma", list(2, 1)), alpha = 1, beta = 2)
  numerator <- .iwmde_ordered_block_numerator(NULL, block, "share")
  independent <- stats::integrate(function(total) stats::dgamma(total, 2, 1) / total *
    exp(-.7 + .3 * total - total^2) * 2 / 2.5, 0, Inf, rel.tol = 1e-10, abs.tol = 0)
  expect_lte(abs(exp(numerator[["log_value"]]) - independent[["value"]]),
    independent[["abs.error"]] + exp(numerator[["log_value"]]) * numerator[["relative_error"]])
  expect_identical(numerator[["method"]], "unshifted_gamma_inverse_supplier")
  block[["shift"]] <- .7
  shifted <- .iwmde_ordered_block_numerator(NULL, block, "share")
  independent_shift <- stats::integrate(function(total) stats::dgamma(total, 2, 1) / (total + .7) *
    exp(-.7 + .3 * total - total^2) * 2 / 2.5, 0, Inf, rel.tol = 1e-10, abs.tol = 0)
  expect_lte(abs(exp(shifted[["log_value"]]) - independent_shift[["value"]]),
    independent_shift[["abs.error"]] + exp(shifted[["log_value"]]) * shifted[["relative_error"]])
  expect_identical(shifted[["method"]], "true_affine_prior_inverse_supplier")
  expect_gt(abs(shifted[["log_value"]] - numerator[["log_value"]]), .1)
})


test_that("batched proper-CDF nodes exactly match scalar open-panel quantiles", {

  for (prior in list(BayesTools::prior("gamma", list(2, 1)),
    BayesTools::prior("normal", list(0, 1), truncation = list(8, 9)),
    BayesTools::prior("beta", list(1.5, 2)))) {
    probability <- .iwmde_ordered_probability_prior(prior)
    cuts <- .iwmde_ordered_base_cuts(probability)
    panels <- .iwmde_ordered_probability_panels(cuts, probability)
    fractions <- .gauss_legendre_nodes(16L)[["nodes"]]
    scalar <- unlist(lapply(panels, probability[["conditional"]], fraction = fractions), use.names = FALSE)
    expect_identical(probability[["conditional_many"]](panels, fractions), scalar)
    calls <- 0L
    result <- .iwmde_ordered_panel_rule(cuts, probability, 16L, function(nodes) {
      calls <<- calls + 1L
      list(log_value = rep(0, length(nodes)))
    })
    expect_identical(calls, 1L)
    expect_equal(result[["log_value"]], 0, tolerance = 1e-13)
  }
})


test_that("new source diagnostics preserve legacy table inspection", {

  template <- .iwmde_empty_public_density_diagnostics()
  added <- c("internal_method", "working_distribution", "normalization_basis", "ordered_error_status",
    "block_numerator_relative_error", "block_denominator_relative_error", "block_max_order", "block_max_panels", "block_native_nodes")
  legacy <- template[, setdiff(names(template), added), drop = FALSE]
  class(legacy) <- class(template)
  restored <- .density_diagnostics_validate(legacy)
  expect_identical(names(restored), names(template))
  expect_identical(vapply(restored, typeof, character(1)), vapply(template, typeof, character(1)))
})


test_that("Gamma completion agrees with independent full-support physical references", {

  # These exact binary inputs were integrated as positive physical kernels at
  # 80/105 digits, with different splits and an explicit infinite tail. They
  # include both shapes, rates .03/7, A=1e-8/1/1e8, K=12/240 and zero scores.
  references <- utils::read.csv(test_path("fixtures", "ordered-gamma-completion.csv"))
  expect_identical(nrow(references), 152L)
  for (i in seq_len(nrow(references))) {
    row <- references[i, ]
    block <- .ordered_block_fixture(BayesTools::prior("gamma", list(row$shape, row$rate)))
    block[["gaussian"]] <- list(A = row$A, H = row$H, C = row$C)
    context <- list(data = list(outcome = data.frame(yi = numeric(row$K))))
    result <- .iwmde_ordered_block_inner(context, block, 0)
    expect_identical(result[["method"]], "completed_gamma_source")
    completion <- result[["completion"]]
    expect_lte(abs(result[["log_value"]] - row$log_reference), completion[["prospective_floor"]])
    expect_identical(block[["counts"]][["max_log_scale"]], completion[["max_log_scale"]])
    expect_identical(block[["counts"]][["gamma_completion"]][["scale_updates"]], 1L)
    expect_lte(completion[["prospective_floor"]], log1p(.025 / 16))
  }
})


test_that("Gamma completion uses supplied effective priors and preserves constant provenance", {

  for (shape in c(2, 3)) for (rate in c(.03, 7)) {
    block <- .ordered_block_fixture(BayesTools::prior("gamma", list(shape, rate)), alpha = 1, beta = 2)
    result <- .iwmde_ordered_block_numerator(NULL, block, "share")
    reference <- stats::integrate(function(total) stats::dgamma(total, shape, rate = rate) / total *
      exp(-.7 + .3 * total - total^2) * 2 / 2.5, 0, Inf, rel.tol = 1e-11, abs.tol = 0)
    floor <- result[["completion"]][["prospective_floor"]]
    expect_lte(abs(exp(result[["log_value"]]) - reference[["value"]]),
      reference[["abs.error"]] + exp(result[["log_value"]]) * expm1(floor))
    expect_identical(result[["relative_error"]], 0)
    expect_identical(result[["method"]], "unshifted_gamma_inverse_supplier")
    expect_identical(result[["completion"]][["effective_prior"]][["parameters"]],
      BayesTools::prior("gamma", list(shape - 1, rate))[["parameters"]])
    expect_identical(block[["counts"]][["gamma_completion"]][["method"]], "completed_gamma_source")
  }
  block <- .ordered_block_fixture(BayesTools::prior("gamma", list(2, 7)))
  block[["gaussian"]] <- list(A = 0, H = 0, C = 1000)
  result <- .iwmde_ordered_block_inner(NULL, block, 0)
  expect_identical(result[["method"]], "declared_constant_source")
  expect_identical(result[["log_value"]], 1000)
  expect_identical(block[["counts"]][["max_log_scale"]], 1)
  expect_null(block[["counts"]][["gamma_completion"]])
  block[["source_state"]] <- list()
  expect_identical(.iwmde_ordered_block_inner(NULL, block, 0)[["method"]], "declared_constant_source")
  expect_null(block[["counts"]][["gamma_completion"]])
  block[["gaussian"]][["H"]] <- 1
  expect_error(.iwmde_ordered_block_inner(NULL, block, 0), class = "RoBMA_ordered_inner_geometry_unresolved")
})


test_that("Gamma completion scope and rejected scales preserve the original CDF dispatcher", {

  make_block <- function(A = 1, H = 1, prior = BayesTools::prior("gamma", list(1, 1))) {
    block <- .ordered_block_fixture(prior)
    block[["gaussian"]] <- list(A = A, H = H, C = 0)
    block
  }
  context <- list(data = list(outcome = data.frame(yi = numeric(12))))
  for (reason in c("source_state", "selected", "invalid_geometry", "dynamic_prior", "truncated_prior", "shape_outside_scope")) {
    block <- make_block()
    probability <- block[["normal_probability"]]
    if (reason == "source_state") block[["source_state"]] <- list()
    if (reason == "selected") block[["selected"]] <- TRUE
    if (reason == "invalid_geometry") block[["gaussian"]] <- NULL
    if (reason == "dynamic_prior") probability[["prior"]][["parameters"]][["rate"]] <- expression(exp(tau))
    if (reason == "truncated_prior") probability <- .iwmde_ordered_probability_prior(
      BayesTools::prior("gamma", list(1, 1), truncation = list(0, 10)))
    if (reason == "shape_outside_scope") probability <- .iwmde_ordered_probability_prior(BayesTools::prior("gamma", list(3, 1)))
    result <- .iwmde_ordered_gamma_completion(context, block, 0, probability)
    expect_false(result[["eligible"]])
    expect_identical(result[["reason"]], reason)
    expect_identical(block[["counts"]][["max_log_scale"]], 1)
  }
  dynamic <- BayesTools::prior("gamma", list(2, 1))
  dynamic[["parameters"]][["rate"]] <- expression(exp(tau))
  expect_error(.iwmde_ordered_probability_prior(dynamic), class = "RoBMA_ordered_proper_cdf_prior")
  rejected <- make_block(H = 100001)
  candidate <- .iwmde_ordered_gamma_completion(context, rejected, 0, rejected[["normal_probability"]])
  expect_identical(candidate[["reason"]], "prospective_floor")
  actual <- .iwmde_ordered_block_inner(context, rejected, 0)
  baseline <- make_block(H = 100001)
  original_completion <- .iwmde_ordered_gamma_completion
  local_mocked_bindings(.iwmde_ordered_gamma_completion = function(...) {
    out <- original_completion(...)
    out[["eligible"]] <- FALSE
    out[["reason"]] <- "forced_original_cdf"
    out
  }, .package = "RoBMA")
  original <- .iwmde_ordered_block_inner(context, baseline, 0)
  fields <- c("log_value", "log_error", "log_exterior", "relative_error", "components",
    "fine_panels", "refinements", "coalesced_cuts", "method")
  expect_identical(actual[fields], original[fields])
  expect_identical(as.numeric(actual[["cuts"]]), as.numeric(original[["cuts"]]))
  expect_identical(as.numeric(actual[["fine_cuts"]]), as.numeric(original[["fine_cuts"]]))
  expect_identical(actual[["guard"]][c("constant", "log_lower", "log_upper", "fraction", "center", "width", "radius", "trace")],
    original[["guard"]][c("constant", "log_lower", "log_upper", "fraction", "center", "width", "radius", "trace")])
  expect_identical(rejected[["counts"]][["max_log_scale"]], baseline[["counts"]][["max_log_scale"]])
  expect_lt(rejected[["counts"]][["max_log_scale"]], candidate[["max_log_scale"]])
  expect_identical(rejected[["counts"]][["gamma_completion"]][["fallback_reasons"]], c(prospective_floor = 1L))
  failure_block <- make_block(.01, 7 - .8, BayesTools::prior("gamma", list(2, 7)))
  expect_error(.iwmde_ordered_block_inner(context, failure_block, 0), class = "RoBMA_ordered_panel_budget_exceeded")
})


test_that("Gamma rows can mix exact completions and CDF nodes on one share partition", {

  block <- .ordered_block_fixture(BayesTools::prior("gamma", list(2, 7)))
  block[["gaussian"]] <- list(A = 1, H = c(-.9, 1), C = 0)
  context <- list(data = list(outcome = data.frame(yi = numeric(12))))
  shares <- c(0, .5, 1)
  results <- lapply(shares, function(share) .iwmde_ordered_block_inner(context, block, share))
  expect_identical(vapply(results, `[[`, character(1), "method"),
    c("proper_cdf_source", "completed_gamma_source", "completed_gamma_source"))
  for (i in seq_along(shares)) {
    H <- -.9 + shares[[i]]
    reference <- stats::integrate(function(total) stats::dgamma(total, 2, rate = 7) *
      exp(H * total - .5 * total^2), 0, Inf, rel.tol = 1e-11, abs.tol = 0)
    reported <- results[[i]][["relative_error"]]
    if (is.null(reported)) reported <- expm1(results[[i]][["completion"]][["prospective_floor"]])
    expect_lte(abs(exp(results[[i]][["log_value"]]) - reference[["value"]]),
      reference[["abs.error"]] + reported * exp(results[[i]][["log_value"]]))
  }
  expect_identical(block[["counts"]][["gamma_completion"]][["accepted"]], 2L)
  expect_identical(block[["counts"]][["gamma_completion"]][["fallback"]], 1L)
  near <- .ordered_block_fixture(BayesTools::prior("gamma", list(1, 1)))
  near[["gaussian"]] <- list(A = 1, H = 90001, C = 0)
  expect_identical(.iwmde_ordered_block_inner(context, near, 0)[["method"]], "completed_gamma_source")
  near[["gaussian"]][["H"]] <- 100001
  expect_identical(.iwmde_ordered_block_inner(context, near, 0)[["method"]], "proper_cdf_source")
})


test_that("huge negative Gamma scores retain the stable guard and original typed fallback", {

  block <- .ordered_block_fixture(BayesTools::prior("gamma", list(2, 3)))
  block[["gaussian"]] <- list(A = 1e-16, H = 0, C = 0)
  context <- list(data = list(outcome = data.frame(yi = numeric(12))))
  completion <- .iwmde_ordered_gamma_completion(context, block, 0, block[["normal_probability"]])
  expect_identical(completion[["z"]], -3e8)
  expect_identical(completion[["reason"]], "stable_negative_bound")
  expect_true(is.na(completion[["condition"]]))
  actual <- tryCatch(.iwmde_ordered_block_inner(context, block, 0), error = identity)
  expect_identical(block[["counts"]][["gamma_completion"]][["fallback_reasons"]], c(stable_negative_bound = 1L))
  expect_identical(block[["counts"]][["gamma_completion"]][["accepted"]], 0L)
  expect_identical(block[["counts"]][["gamma_completion"]][["z_range"]], c(Inf, -Inf))
  presented <- .iwmde_ordered_block_counts(block)[["gamma_completion"]]
  expect_identical(presented[["z_range"]], c(NA_real_, NA_real_))
  expect_identical(presented[["condition_range"]], c(NA_real_, NA_real_))
  baseline <- .ordered_block_fixture(BayesTools::prior("gamma", list(2, 3)))
  baseline[["gaussian"]] <- block[["gaussian"]]
  local_mocked_bindings(.iwmde_ordered_gamma_completion = function(...) {
    completion[["reason"]] <- "forced_original_cdf"
    completion
  }, .package = "RoBMA")
  original <- tryCatch(.iwmde_ordered_block_inner(context, baseline, 0), error = identity)
  expect_identical(actual[["method"]], "proper_cdf_source")
  expect_identical(actual[c("log_value", "log_error", "log_exterior", "relative_error", "components",
    "fine_panels", "refinements", "coalesced_cuts", "method")], original[c("log_value", "log_error",
    "log_exterior", "relative_error", "components", "fine_panels", "refinements", "coalesced_cuts", "method")])
  expect_identical(block[["counts"]][["max_log_scale"]], baseline[["counts"]][["max_log_scale"]])
})


test_that("proper-CDF block denominators retain general Gamma and truncated Normal support", {

  for (prior in list(BayesTools::prior("gamma", list(3, 1)),
    BayesTools::prior("normal", list(.3, 1.2), truncation = list(-.5, 1.4)))) {
    block <- .ordered_block_fixture(prior)
    actual <- .iwmde_ordered_block_denominator(NULL, block)
    support <- .iwmde_prior_support(prior)
    reference <- stats::integrate(function(share) vapply(share, function(s) {
      stats::integrate(function(total) exp(-.7 + (.3 + .8 * s) * total -
        .5 * (2 + 10 * s^2) * total^2 + BayesTools::lpdf(prior, total)),
        support[[1L]], support[[2L]], rel.tol = 1e-10, abs.tol = 0)[["value"]] * stats::dbeta(s, 3, 1)
    }, numeric(1)), 0, 1, rel.tol = 1e-9, abs.tol = 0)
    expect_lte(abs(exp(actual[["log_value"]]) - reference[["value"]]),
      reference[["abs.error"]] + actual[["relative_error"]] * exp(actual[["log_value"]]))
    if (prior[["distribution"]] == "gamma") {
      expect_identical(block[["counts"]][["gamma_completion"]][["accepted"]], 0L)
      expect_gt(block[["counts"]][["gamma_completion"]][["fallback_reasons"]][["shape_outside_scope"]], 0L)
    }
  }
})
