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

test_that("unique random effects are judged within the returned output cells", {

  # One group whose rows alternate between the two levels of an indicator
  # (compound-symmetric) term: the compiled group map is 1 1 1 1, the returned
  # cells are (group, level) combinations.
  dat <- data.frame(
    yi = c(.1, .2, .3, .4), vi = c(.04, .05, .04, .05), id = factor(rep("a", 4L)),
    g  = factor(c("u", "v", "u", "v")), x = c(0, 1, 0, 1), z = c(0, 0, 1, 1)
  )
  fixture <- function(predictor) {
    object <- brma.mv(
      yi = yi, V = vi, data = dat, random = ~ cs(g | id),
      scale = stats::as.formula(paste("~", predictor)),
      known_v_parameterization = "block_mvn", measure = "GEN",
      prior_unit_information_sd = 1, only_priors = TRUE, silent = TRUE
    )
    samples <- matrix(
      c(.05, .3, .2, log(2), .06, .5, .3, log(2)), 2L, byrow = TRUE,
      dimnames = list(NULL, c("mu_intercept", "mu__xREx__id_rho",
                              "log_tau_intercept", paste0("log_tau_", predictor)))
    )
    list(object = object, samples = samples)
  }

  # x is constant within each (group, level) cell: the default output works
  # and carries the expanded contributions of the first row of each cell.
  constant <- fixture("x")
  term     <- .fitted_formula_design(constant$object, "mu")[["random_effects"]][[1L]]
  cells    <- .ranef_unique_level_term(term, term[["block_name"]])[["group_map"]]
  expect_identical(as.integer(term[["group_map"]]), rep(1L, 4L))
  expect_identical(as.integer(cells), c(1L, 2L, 1L, 2L))
  expect_true(.ranef_scale_design_varies_within_cells(
    constant$object, "log_tau", term[["group_map"]]
  ))
  expect_false(.ranef_scale_design_varies_within_cells(
    constant$object, "log_tau", cells
  ))
  unique_effects <- ranef(constant$object, .posterior_samples = constant$samples)
  expanded       <- ranef(constant$object, expand = TRUE,
                          .posterior_samples = constant$samples)
  expect_identical(ncol(unique_effects), 2L)
  expect_equal(unname(as.matrix(unique_effects)), unname(as.matrix(expanded))[, 1:2],
               tolerance = 1e-12)
  expect_equal(unname(as.matrix(expanded))[, 1L], unname(as.matrix(expanded))[, 3L],
               tolerance = 1e-12)

  # z changes within the first cell (rows 1 and 3): refused.
  varying <- fixture("z")
  expect_error(
    ranef(varying$object, .posterior_samples = varying$samples),
    class = "RoBMA_ranef_expand_required"
  )
  expect_identical(
    ncol(ranef(varying$object, expand = TRUE, .posterior_samples = varying$samples)),
    4L
  )
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
        .ranef_block_sd_varies_within_cells(
          object, term,
          .ranef_unique_level_term(term, term[["block_name"]])[["group_map"]]
        )
      }, logical(1L)),
      vapply(terms, `[[`, character(1L), "block_name")
    )
  }

  # Singleton groups (effect within study) cannot vary; the study block does
  # whenever the scale predictor changes between its rows.
  expect_identical(blocks("x"), c(effect_study = FALSE, study = TRUE))
  expect_identical(blocks("g"), c(effect_study = FALSE, study = FALSE))
})

# A constructor-compiled BMA.mv object with a mixture location prior and a
# synthetic posterior, finalized through the BayesTools fit contract.
.ranef_conditional_fixture <- function() {

  dat <- data.frame(yi = c(.1, .2, .3, .4), obs = factor(1:4),
                    study = factor(c("a", "a", "b", "b")))
  V <- matrix(c(.04, .01, 0, 0, .01, .05, 0, 0, 0, 0, .04, .012, 0, 0, .012, .05), 4L)
  object <- BMA.mv(
    yi = yi, V = V, random = ~ 1 | study, data = dat, measure = "GEN",
    prior_unit_information_sd = 1, prior_effect = BayesTools::prior("normal", list(0, 1)),
    prior_effect_null = BayesTools::prior("spike", list(0)),
    known_v_parameterization = "block_mvn", only_priors = TRUE, silent = TRUE
  )
  args <- .create_jags_formula_args(object[["data"]], object[["priors"]])
  compiled <- lapply(names(args[["formula_list"]]), function(parameter) {
    BayesTools::JAGS_formula(
      formula = args[["formula_list"]][[parameter]], parameter = parameter,
      data = args[["formula_data_list"]][[parameter]],
      prior_list = args[["formula_prior_list"]][[parameter]],
      formula_scale = args[["formula_scale_list"]][[parameter]],
      prior_random = args[["formula_random_prior_list"]][[parameter]],
      random_effects_compile = args[["formula_random_effects_compile_list"]][[parameter]]
    )
  })
  names(compiled) <- names(args[["formula_list"]])
  design <- lapply(compiled, `[[`, "formula_design")
  allocation <- design[["mu"]][["random_allocations"]][[1L]]
  gate <- allocation[["inclusion"]][[1L]]
  samples <- cbind(
    mu_intercept = c(.1, .2, .3, .4), mu_intercept_indicator = c(1, 2, 1, 2),
    sd = c(.3, .4, .5, .6), prob = .5, ind = 1
  )
  colnames(samples) <- c("mu_intercept", "mu_intercept_indicator",
                         allocation[["source_node"]], gate[["prob_name"]],
                         gate[["indicator_name"]])

  priors <- list()
  for (item in compiled) priors <- c(priors, item[["prior_list"]])
  fit <- structure(
    list(mcmc = coda::mcmc.list(coda::mcmc(samples)), sample = nrow(samples)),
    class = c("runjags", "list")
  )
  attr(fit, "formula_design") <- design
  attr(fit, "formula_scale")  <- lapply(compiled, `[[`, "formula_scale")
  attr(fit, "prior_list")     <- c(priors, .create_fit_priors(object[["data"]],
                                                                object[["priors"]]))
  object[["formula_design"]] <- design
  object[["fit"]]            <- as_bayestools_fit(fit)
  list(object = object, samples = samples)
}

test_that("conditional multivariate ranef selects posterior rows after the decomposition", {

  fixture <- .ranef_conditional_fixture()
  object  <- fixture[["object"]]
  keep    <- .conditional_parameter_rows(
    object, .conditional_effect_parameters(object), fixture[["samples"]], "OR"
  )
  expect_identical(keep, c(FALSE, TRUE, FALSE, TRUE))

  unconditional <- ranef(object, simplify = FALSE, expand = TRUE,
                         .posterior_samples = fixture[["samples"]])
  conditional   <- ranef(object, simplify = FALSE, expand = TRUE, conditional = TRUE,
                         .posterior_samples = fixture[["samples"]])
  expect_named(conditional, names(unconditional))
  for (block in names(unconditional)) {
    expect_equal(
      unname(as.matrix(conditional[[block]])),
      unname(as.matrix(unconditional[[block]]))[keep, , drop = FALSE],
      tolerance = 1e-12
    )
  }

  # The unique-level and total routes select the same rows.
  unique_effects <- ranef(object, conditional = TRUE, .posterior_samples = fixture[["samples"]])
  expect_identical(nrow(as.matrix(unique_effects)), 2L)
  total <- ranef(object, component = "total", expand = TRUE, conditional = TRUE,
                 .posterior_samples = fixture[["samples"]])
  expect_equal(
    unname(as.matrix(total)),
    unname(Reduce(`+`, lapply(unconditional, as.matrix)))[keep, , drop = FALSE],
    tolerance = 1e-12
  )
  expect_error(
    ranef(object, conditional = NA, .posterior_samples = fixture[["samples"]]),
    "conditional"
  )
})
