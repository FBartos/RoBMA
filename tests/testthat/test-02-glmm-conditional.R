context("Conditional GLMM log tails")


test_that("conditional binomial matrices retain paired finite log tails", {

  ai <- c(0L, 2L, 3L)
  ci <- c(0L, 1L, 0L)
  n1 <- c(1L, 2L, 5L)
  n2 <- c(1L, 3L, 4L)
  mu <- matrix(c(80, -80, .3, -.8, .7, -.4), 2L)
  eta <- matrix(c(0, 0, -.2, .9, -1.1, .6), 2L)
  # Scalar analytic count formulas, retaining coefficients once per cell.
  count_term <- function(count, value) if (count == 0L) 0 else count * value
  expected <- vapply(seq_along(ai), function(k) {
    vapply(seq_len(nrow(mu)), function(s) {
      arm1 <- eta[s, k] + mu[s, k] / 2
      arm2 <- eta[s, k] - mu[s, k] / 2
      lp <- function(x) -max(0, -x) - log1p(exp(-abs(x)))
      lq <- function(x) -max(0, x) - log1p(exp(-abs(x)))
      lchoose(n1[k], ai[k]) + lchoose(n2[k], ci[k]) +
        count_term(ai[k], lp(arm1)) + count_term(n1[k] - ai[k], lq(arm1)) +
        count_term(ci[k], lp(arm2)) + count_term(n2[k] - ci[k], lq(arm2))
    }, numeric(1))
  }, numeric(nrow(mu)))
  value <- .outcome_pdf.binom_conditional(ai, ci, n1, n2, mu, eta)
  expect_identical(dim(value), c(2L, 3L))
  expect_equal(value, expected, tolerance = 1e-12)
  expect_equal(value[, 1L], rep(-40 - 2 * log1p(exp(-40)), 2L), tolerance = 1e-12)
  powers <- c(.5, 1.7, 2)
  expect_equal(.outcome_pdf_sum.binom_conditional(ai, ci, n1, n2, mu, eta, powers),
    rowSums(sweep(expected, 2L, powers, "*")), tolerance = 1e-12)
  for (effect in c(Inf, -Inf, NA_real_, NaN)) {
    for (counts in list(c(2L, 0L), c(0L, 3L))) {
      matrix_value <- .outcome_pdf.binom_conditional(counts[1L], counts[2L],
        2L, 3L, matrix(effect), matrix(0))
      native_value <- .outcome_pdf_sum.binom_conditional(counts[1L], counts[2L],
        2L, 3L, matrix(effect), matrix(0))
      expect_equal(as.numeric(matrix_value), native_value)
    }
  }
})


test_that("conditional Poisson matrices retain combined log rates and limits", {

  x1 <- c(1L, 0L, 3L)
  x2 <- c(0L, 0L, 2L)
  t1 <- c(1, 2, 7)
  t2 <- c(1, 3, 5)
  mu <- matrix(c(0, 0, .3, -.8, .7, -.4), 2L)
  phi <- matrix(c(-800, -800, -800, -800, -1.1, .6), 2L)
  expected <- vapply(seq_along(x1), function(k) {
    vapply(seq_len(nrow(mu)), function(s) {
      rate1 <- log(t1[k]) + phi[s, k] + mu[s, k] / 2
      rate2 <- log(t2[k]) + phi[s, k] - mu[s, k] / 2
      term <- function(x, rate) (if (x == 0L) 0 else x * rate) - exp(rate) - lgamma(x + 1)
      term(x1[k], rate1) + term(x2[k], rate2)
    }, numeric(1))
  }, numeric(nrow(mu)))
  value <- .outcome_pdf.pois_conditional(x1, x2, t1, t2, mu, phi)
  expect_identical(dim(value), c(2L, 3L))
  expect_equal(value, expected, tolerance = 1e-12)
  expect_identical(value[, 1L], c(-800, -800))
  expect_identical(value[, 2L], c(0, 0))
  powers <- c(.5, 1.7, 2)
  expect_equal(.outcome_pdf_sum.pois_conditional(x1, x2, t1, t2, mu, phi, powers),
    rowSums(sweep(expected, 2L, powers, "*")), tolerance = 1e-12)
  for (log_rate in c(800, Inf, -Inf, NA_real_, NaN)) {
    for (counts in list(c(0L, 0L), c(1L, 0L))) {
      value <- .outcome_pdf.pois_conditional(counts[1L], counts[2L], 1, 1,
        matrix(0), matrix(log_rate))
      native <- .outcome_pdf_sum.pois_conditional(counts[1L], counts[2L], 1, 1,
        matrix(0), matrix(log_rate))
      expect_equal(as.numeric(value), native)
    }
  }
  expect_identical(.poisson_log_kernel(0L, c(-Inf, -800)), c(0, 0))
  expect_identical(.poisson_log_kernel(1L, c(800, Inf)), c(-Inf, -Inf))
})


.glmm_conditional_callback_replay <- function(family, weighted = FALSE) {

  effect_prior <- BayesTools::prior("normal", list(mean = 0, sd = 100))
  tau_prior <- BayesTools::prior("normal", list(mean = 0, sd = 1),
    truncation = list(lower = 0, upper = Inf))
  nuisance <- if (family == "bin") BayesTools::prior("beta", list(1, 1)) else
    BayesTools::prior("normal", list(mean = 0, sd = 1000))
  n <- if (weighted) 2L else 1L
  args <- list(data = data.frame(row = seq_len(n)), mods = ~1,
    prior_effect = effect_prior, prior_heterogeneity = tau_prior,
    only_priors = TRUE, silent = TRUE)
  if (weighted) args$weights <- c(.5, 1.7)
  if (family == "bin") {
    args <- c(args, list(ai = if (weighted) c(1L, 2L) else 0L,
      ci = if (weighted) c(0L, 1L) else 0L,
      n1i = rep(3L, n), n2i = rep(2L, n), measure = "OR", prior_baserate = nuisance))
    if (!weighted) args$n1i <- args$n2i <- 1L
  } else {
    args <- c(args, list(x1i = if (weighted) c(1L, 2L) else 1L,
      x2i = if (weighted) c(0L, 1L) else 0L, t1i = rep(1, n), t2i = rep(1, n),
      measure = "IRR", prior_lograte = nuisance))
  }
  object <- do.call(RoBMA::brma.glmm, args)
  compiled <- do.call(BayesTools::JAGS_formula,
    .object_bayestools_formula_args(object, parameter = "mu", source = "mods"))
  parameters <- list(mu = if (weighted) .35 else if (family == "bin") 80 else 0,
    tau = .2, theta = if (weighted) c(.25, -.5) else 0)
  if (family == "bin") parameters$pi <- if (weighted) c(.2, .7) else .5 else
    parameters$phi <- if (weighted) c(-.7, .3) else -800
  predicates <- c(is_scale = ".is_data_scale", is_multilevel = ".is_data_multilevel",
    is_weights = ".is_data_weights", is_known_v = ".is_data_known_v", is_random = ".is_data_random")
  flags <- lapply(predicates, function(name) get(name, asNamespace("RoBMA"))(object$data))
  flags$is_PET <- .is_priors_PET(object$priors)
  flags$is_PEESE <- .is_priors_PEESE(object$priors)
  flags$is_weightfunction <- .is_priors_weightfunction(object$priors)
  value <- do.call(.log_posterior, c(list(parameters = parameters,
    data = .create_fit_data(object$data, object$priors), model_data = object$data,
    effect_direction = .data_effect_direction(object$data),
    outcome_type = .data_outcome_type(object$data)), flags))
  list(object = object, compiled = compiled, parameters = parameters, value = value,
    prior_logs = c(BayesTools::lpdf(effect_prior, parameters$mu),
      BayesTools::lpdf(tau_prior, parameters$tau),
      BayesTools::lpdf(nuisance, if (family == "bin") parameters$pi else parameters$phi)))
}


test_that("genuine compiler-backed GLMM callbacks keep tails and apply powers once", {

  for (family in c("bin", "pois")) {
    tail <- .glmm_conditional_callback_replay(family)
    expect_s3_class(tail$object, "only_priors.brma")
    expect_null(tail$object[["fit", exact = TRUE]])
    expect_identical(tail$compiled$data$N_mu, 1L)
    expect_true(all(is.finite(tail$prior_logs)))
    expect_equal(tail$value, if (family == "bin") -40 else -800, tolerance = 1e-12)
    replay <- .glmm_conditional_callback_replay(family, weighted = TRUE)
    p <- replay$parameters
    effect <- p$mu + p$tau * p$theta
    if (family == "bin") {
      expected <- stats::dbinom(c(1L, 2L), 3L,
        stats::plogis(stats::qlogis(p$pi) + effect / 2), log = TRUE) +
        stats::dbinom(c(0L, 1L), 2L,
          stats::plogis(stats::qlogis(p$pi) - effect / 2), log = TRUE)
    } else {
      expected <- stats::dpois(c(1L, 2L), exp(p$phi + effect / 2), log = TRUE) +
        stats::dpois(c(0L, 1L), exp(p$phi - effect / 2), log = TRUE)
    }
    expect_equal(replay$value, sum(c(.5, 1.7) * expected), tolerance = 1e-12)
  }
})
