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
