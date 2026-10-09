test_that("two-level total random effects keep the estimate names", {

  object <- brma(
    yi = c(-.2, .1, .4), sei = c(.2, .2, .2), measure = "OR",
    prior_unit_information_sd = 1, only_priors = TRUE
  )
  samples <- cbind(mu = c(-.2, 0, .2, .4), tau = c(.1, .2, .3, .4))

  estimate <- ranef(object, component = "estimate", .posterior_samples = samples)
  total    <- ranef(object, component = "total", .posterior_samples = samples)
  expect_identical(colnames(total), colnames(estimate))
  expect_false(anyDuplicated(colnames(total)) > 0L)
  expect_equal(unname(as.matrix(total)), unname(as.matrix(estimate)), tolerance = 0)

  # Expanded totals are unchanged.
  expanded <- ranef(object, component = "total", expand = TRUE, .posterior_samples = samples)
  expect_identical(colnames(expanded), colnames(estimate))
})

.ranef_routing_mv <- function(scale = NULL, drop_unused = FALSE) {

  dat <- data.frame(
    yi = c(.1, .2, .3, .4), vi = c(.04, .05, .04, .05),
    id = factor(c("a", "a", "b", "b"), levels = c("a", "b", "unused")),
    x  = c(0, 1, 0, 1)
  )
  if (drop_unused) dat[["id"]] <- droplevels(dat[["id"]])
  args <- list(
    yi = quote(yi), V = quote(vi), data = dat, random = ~1 | id,
    known_v_parameterization = "block_mvn", measure = "GEN",
    prior_unit_information_sd = 1, only_priors = TRUE, silent = TRUE
  )
  if (!is.null(scale)) args[["scale"]] <- scale
  object <- do.call(brma.mv, args)
  design <- .fitted_formula_design(object, "mu", required = TRUE)
  term   <- design[["random_effects"]][[1L]]
  list(object = object, term = term)
}

test_that("unique random effects skip declared grouping levels without rows", {

  fixture <- .ranef_routing_mv()
  expect_identical(fixture$term[["group_levels"]], c("a", "b", "unused"))
  sd_name <- fixture$term[["sd_parameter_names"]][[1L]]
  samples <- matrix(c(.05, .2, .06, .3), 2L, byrow = TRUE,
                    dimnames = list(NULL, c("mu", sd_name)))

  unique_effects <- ranef(fixture$object, .posterior_samples = samples)
  expanded       <- ranef(fixture$object, expand = TRUE, .posterior_samples = samples)
  expect_identical(colnames(unique_effects), c("u_id[a]", "u_id[b]"))
  expect_equal(unname(as.matrix(unique_effects)), unname(as.matrix(expanded))[, c(1L, 3L)],
               tolerance = 1e-12)
})

test_that("observed grouping levels follow declared order with and without metadata", {

  term <- list(group_map = c(3L, 3L, 1L), group_levels = c("a", "unused", "c"))
  expect_identical(
    .ranef_observed_group_levels(term),
    list(group_map = c(2L, 2L, 1L), group_levels = c("a", "c"))
  )
  term[["group_observed_levels"]] <- c("c", "a")
  expect_identical(
    .ranef_observed_group_levels(term),
    list(group_map = c(2L, 2L, 1L), group_levels = c("a", "c"))
  )
})

.ranef_routing_scaled_mv <- function(predictor) {

  dat <- data.frame(
    yi = c(.1, .2, .3, .4), vi = c(.04, .05, .04, .05),
    id = factor(c("a", "a", "b", "b")), x = c(0, 1, 0, 1), g = c(0, 0, 1, 1)
  )
  scale <- stats::as.formula(paste("~", predictor))
  object <- brma.mv(
    yi = yi, V = vi, data = dat, random = ~1 | id, scale = list(id = scale),
    known_v_parameterization = "block_mvn", measure = "GEN",
    prior_unit_information_sd = 1, only_priors = TRUE, silent = TRUE
  )
  parameter <- paste0("log_tau_id_", predictor)
  samples <- matrix(
    c(.05, .2, log(2), .06, .3, log(2)), 2L, byrow = TRUE,
    dimnames = list(NULL, c("mu", "log_tau_id_intercept", parameter))
  )
  list(object = object, dat = dat, samples = samples)
}

test_that("unique random effects refuse a standard deviation that varies within a group", {

  fixture <- .ranef_routing_scaled_mv("x")
  expect_error(
    ranef(fixture$object, .posterior_samples = fixture$samples),
    class = "RoBMA_ranef_expand_required"
  )

  # The expanded output keeps one contribution per row: the Gaussian
  # conditional mean of a group-shared effect with row-specific SDs.
  expanded <- ranef(fixture$object, expand = TRUE, .posterior_samples = fixture$samples)
  dat      <- fixture$dat
  x_scaled <- (dat$x - mean(dat$x)) / stats::sd(dat$x)
  same     <- outer(as.integer(dat$id), as.integer(dat$id), "==")
  expected <- t(vapply(1:2, function(draw) {
    tau <- fixture$samples[draw, "log_tau_id_intercept"] *
      exp(fixture$samples[draw, "log_tau_id_x"] * x_scaled)
    covariance <- outer(tau, tau) * same
    drop(covariance %*% solve(diag(dat$vi) + covariance,
                              dat$yi - fixture$samples[draw, "mu"]))
  }, numeric(4L)))
  expect_equal(unname(as.matrix(expanded)), expected, tolerance = 1e-10)
})

test_that("unique random effects are kept for a scale design constant within groups", {

  fixture  <- .ranef_routing_scaled_mv("g")
  unique   <- ranef(fixture$object, .posterior_samples = fixture$samples)
  expanded <- ranef(fixture$object, expand = TRUE, .posterior_samples = fixture$samples)
  expect_identical(colnames(unique), c("u_id[a]", "u_id[b]"))
  expect_equal(unname(as.matrix(unique)), unname(as.matrix(expanded))[, c(1L, 3L)],
               tolerance = 1e-12)
})

test_that("clustered unique random effects refuse a scale design that varies within a cluster", {

  dat <- data.frame(
    yi = c(.1, .2, .3, .4), sei = rep(.2, 4L), cl = c("a", "a", "b", "b"),
    x = c(0, 1, 0, 1), g = c(0, 0, 1, 1)
  )
  fit <- function(predictor) brma(
    yi = yi, sei = sei, cluster = cl, scale = stats::as.formula(paste("~", predictor)),
    data = dat, measure = "GEN", prior_unit_information_sd = 1, only_priors = TRUE,
    silent = TRUE
  )
  samples <- function(predictor) {
    out <- rbind(c(.2, .8, .1, .5, .4, -.3), c(.3, .4, .2, .6, -.2, .5))
    colnames(out) <- c("log_tau_intercept", paste0("log_tau_", predictor), "mu", "rho",
                       "gamma[1]", "gamma[2]")
    out
  }

  varying <- fit("x")
  expect_error(
    ranef(varying, .posterior_samples = samples("x"), simplify = FALSE),
    class = "RoBMA_ranef_expand_required"
  )
  expanded <- ranef(varying, .posterior_samples = samples("x"), simplify = FALSE,
                    expand = TRUE)
  expect_identical(ncol(expanded[["cluster"]]), 4L)

  constant <- fit("g")
  unique   <- ranef(constant, .posterior_samples = samples("g"), simplify = FALSE)
  expanded <- ranef(constant, .posterior_samples = samples("g"), simplify = FALSE,
                    expand = TRUE)
  expect_equal(unname(as.matrix(unique[["cluster"]])),
               unname(as.matrix(expanded[["cluster"]]))[, c(1L, 3L)],
               tolerance = 1e-12)
})

test_that("a common scale is inspected per nested random-effect block", {

  dat <- data.frame(
    yi     = c(0.08, 0.13, 0.18, 0.20, 0.01, 0.05),
    study  = rep(c("s1", "s2", "s3"), each = 2L),
    effect = rep(c("a", "b"), 3L),
    x      = c(0, 1, 0, 1, 0, 1),
    g      = c(0, 0, 1, 1, 2, 2)
  )
  V <- matrix(0, 6L, 6L)
  for (g in 1:3) {
    index <- (2L * g - 1L):(2L * g)
    V[index, index] <- matrix(c(0.04, 0.018, 0.018, 0.05), 2L)
  }
  blocks <- function(predictor) {
    object <- brma.mv(
      yi = yi, V = V, data = dat, random = ~ 1 | study / effect,
      scale = stats::as.formula(paste("~", predictor)),
      known_v_parameterization = "block_mvn", measure = "GEN",
      prior_unit_information_sd = 1, only_priors = TRUE, silent = TRUE
    )
    terms <- .fitted_formula_design(object, "mu", required = TRUE)[["random_effects"]]
    stats::setNames(
      vapply(terms, function(term) {
        .ranef_block_sd_varies_within_groups(object, term, term[["group_map"]])
      }, logical(1L)),
      vapply(terms, `[[`, character(1L), "block_name")
    )
  }

  # Singleton groups (effect within study) cannot vary; the study block does
  # whenever the scale predictor changes between its rows.
  expect_identical(blocks("x"), c(effect_study = FALSE, study = TRUE))
  expect_identical(blocks("g"), c(effect_study = FALSE, study = FALSE))
})
