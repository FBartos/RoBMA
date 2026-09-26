context("Hypothesis plans: one eligibility source for hypothesis() and hypothesis_quantities()")

# Small single-chain fits: every expectation is either an identity on the
# same posterior draws, an analytic prior ordinate, or a Monte Carlo
# comparison with its error stated.
.plan_fit_cache <- new.env(parent = emptyenv())

.plan_fits <- function() {

  if (!is.null(.plan_fit_cache[["fits"]])) {
    return(.plan_fit_cache[["fits"]])
  }
  fit <- function(data, mods, prior_mods = NULL, contrast = "treatment") {
    arguments <- list(
      yi = quote(yi), sei = quote(sei), mods = mods, data = data,
      measure = "SMD", set_contrast_factor_predictors = contrast,
      chains = 1, sample = 1000, burnin = 200, adapt = 100, seed = 1,
      silent = TRUE
    )
    if (!is.null(prior_mods)) {
      arguments[["prior_mods"]] <- prior_mods
    }
    suppressWarnings(do.call(brma, arguments))
  }
  set.seed(1)
  k <- 48L
  factors <- data.frame(
    g1  = factor(rep(c(5, 10, 20), length.out = k), levels = c(5, 10, 20)),
    g2  = factor(rep(1:4, length.out = k)[sample.int(k)], levels = 1:4),
    sei = stats::runif(k, 0.1, 0.3)
  )
  factors[["yi"]] <- stats::rnorm(
    k,
    c(0, 0.2, 0.4)[as.integer(factors[["g1"]])] +
      c(0, 0.1, 0.3, 0.5)[as.integer(factors[["g2"]])],
    factors[["sei"]]
  )
  fits <- lapply(
    c(treatment = "treatment", meandif = "meandif",
      orthonormal = "orthonormal", independent = "independent"),
    function(contrast) {
      fit(
        data     = factors,
        mods     = if (identical(contrast, "independent")) ~ g1 + g2 - 1 else ~ g1 + g2,
        contrast = contrast
      )
    }
  )
  for (n_levels in 2:4) {
    set.seed(3)
    labels <- c("a", "b", "c", "d")[seq_len(n_levels)]
    ordered <- data.frame(
      g   = factor(rep(labels, length.out = k), levels = labels),
      sei = stats::runif(k, 0.1, 0.3)
    )
    ordered[["yi"]] <- stats::rnorm(
      k, seq(0, 0.3, length.out = n_levels)[as.integer(ordered[["g"]])],
      ordered[["sei"]]
    )
    fits[[paste0("ordered_", n_levels)]] <- fit(
      data       = ordered,
      mods       = ~ g,
      prior_mods = list(g = BayesTools::prior_ordered(
        BayesTools::prior("normal", list(0, 1))
      ))
    )
  }
  set.seed(5)
  scaled <- data.frame(
    x   = stats::rnorm(k, 1, 2),
    g   = factor(rep(c("a", "b"), length.out = k)),
    sei = stats::runif(k, 0.1, 0.3)
  )
  scaled[["yi"]] <- stats::rnorm(k, 0.1 * scaled[["x"]], scaled[["sei"]])
  fits[["scaled_interaction"]] <- fit(data = scaled, mods = ~ x * g)
  fits[["model_averaged"]] <- suppressWarnings(BMA.norm(
    yi = yi, sei = sei, mods = ~ g1, data = factors, measure = "SMD",
    chains = 1, sample = 1000, burnin = 200, adapt = 100, seed = 1,
    silent = TRUE
  ))
  # A RoBMA ensemble with publication-bias components; its bias rows are not
  # hypothesis targets.
  fits[["robma_mixture"]] <- suppressWarnings(RoBMA(
    yi = yi, sei = sei, mods = ~ g1, data = factors, measure = "SMD",
    chains = 1, sample = 1000, burnin = 200, adapt = 100, seed = 1,
    silent = TRUE
  ))
  .plan_fit_cache[["fits"]] <- fits

  return(fits)
}


test_that("hypothesis_quantities() renders the plans that hypothesis() executes", {

  skip_on_cran()
  fits <- .plan_fits()
  for (name in names(fits)) {
    .expect_plans_consistent(fits[[name]], info = name)
  }
  fixtures <- list(
    gated     = gated_random_object(),
    ungated   = single_sd_random_object(BayesTools::prior("gamma", list(2, 2))),
    fixed     = single_sd_random_object(BayesTools::prior("spike", list(location = .2))),
    nested    = nested_allocation_random_object()
  )
  for (name in names(fixtures)) {
    .expect_plans_consistent(fixtures[[name]], info = name, run_precomputed = FALSE)
  }
})


# A small known-V fit with two scale formulas (one per random component).
.two_scale_fit_cache <- new.env(parent = emptyenv())

.two_scale_fit <- function() {

  if (!is.null(.two_scale_fit_cache[["fit"]])) {
    return(.two_scale_fit_cache[["fit"]])
  }
  data <- data.frame(
    yi     = c(0.08, 0.13, 0.18, 0.20, 0.01, 0.05),
    study  = rep(c("s1", "s2", "s3"), each = 2L),
    effect = rep(c("a", "b"), 3L),
    x      = c(0, 1, 0, 1, 0, 1)
  )
  V <- kronecker(diag(3L), matrix(c(0.04, 0.018, 0.018, 0.05), nrow = 2L))
  fit <- suppressWarnings(brma.mv(
    yi = yi, V = V, data = data, measure = "GEN",
    random = list(study = ~ 1 | study, effect = ~ 1 | study:effect),
    scale  = list(study = ~ x, effect = ~ x),
    prior_unit_information_sd = 1,
    chains = 1, sample = 500, burnin = 200, adapt = 200, seed = 1,
    silent = TRUE,
    convergence_checks = set_convergence_checks(max_Rhat = NULL, min_ESS = NULL)
  ))
  .two_scale_fit_cache[["fit"]] <- fit

  return(fit)
}


test_that("hypothesis_quantities() names quantities with shared aliases by their selector", {

  skip_on_cran()
  # Two scale formulas share the aliases 'intercept' and 'x' in the scale
  # component, where hypothesis() refuses them as referring to several
  # parameters; the location intercept keeps its alias 'intercept'.
  fit <- .two_scale_fit()

  expect_error(
    hypothesis(fit, "intercept = 0.5", component = "scale"),
    "Hypothesis references multiple model parameters",
    fixed = TRUE
  )
  # qCMDE/IWMDE are refused for every quantity (several scale formulas), and
  # every refused statement stops with its plan's refusal.
  quantities <- .expect_plans_consistent(
    fit, info = "two scale formulas", run_precomputed = TRUE
  )
  scale_rows <- quantities[["component"]] == "scale"
  expect_true(all(quantities[["point_test"]][scale_rows]))
  expect_true(all(quantities[["direction_test"]][scale_rows]))
  # The shared aliases are not listed for the scale quantities; every listed
  # alias names its own quantity (checked by .expect_plans_consistent()).
  expect_false(any(quantities[["alias"]][scale_rows] %in% c("intercept", "x")))
  expect_true("intercept" %in% quantities[["alias"]][
    quantities[["parameter"]] == "mu_intercept"
  ])
  expect_setequal(
    unique(quantities[["parameter"]][scale_rows]),
    c("log_tau_study_intercept", "log_tau_study_x",
      "log_tau_effect_intercept", "log_tau_effect_x")
  )

  metadata <- .brma_parameter_catalog_metadata(fit)
  entries  <- metadata[["entries"]]
  roots <- vapply(seq_len(nrow(entries)), function(i) {
    .hypothesis_quantities_reference_root(
      metadata,
      as.list(entries[i, setdiff(names(entries), "aliases"), drop = FALSE])
    )
  }, character(1))
  expected <- c(
    mu_intercept             = "intercept",
    log_tau_study_intercept  = "log_tau_study_intercept",
    log_tau_study_x          = "log_tau_study_x",
    log_tau_effect_intercept = "log_tau_effect_intercept",
    log_tau_effect_x         = "log_tau_effect_x"
  )
  expect_setequal(entries[["parameter"]], names(expected))
  expect_identical(
    roots[match(names(expected), entries[["parameter"]])],
    unname(expected)
  )
})


test_that("qCMDE/IWMDE point hypotheses are refused for several scale formulas", {

  skip_on_cran()
  fit    <- .two_scale_fit()
  reason <- paste0(
    "qCMDE/IWMDE density estimation is unavailable for models with ",
    "several scale formulas. Use density_method = 'KDE'."
  )
  expect_identical(
    .iwmde_capability(object = fit, density_method = "qCMDE"),
    list(available = FALSE, reason = reason)
  )

  # The location intercept and the scale slopes list KDE only, with the
  # reason; the scale intercepts were KDE-only before (exp(affine) targets).
  quantities <- hypothesis_quantities(fit)
  for (parameter in c("mu_intercept", "log_tau_study_x", "log_tau_effect_x")) {
    rows <- quantities[quantities[["parameter"]] == parameter, , drop = FALSE]
    expect_identical(unique(rows[["point_test_methods"]]), "KDE", info = parameter)
    expect_match(unique(rows[["reason"]]), reason, fixed = TRUE, info = parameter)
  }

  # hypothesis() with its default method stops with the method refusal; KDE
  # evaluates the same statement.
  expect_error(
    hypothesis(fit, "intercept = 0", component = "mods"),
    reason,
    fixed = TRUE,
    class = "RoBMA_hypothesis_method"
  )
  expect_error(
    hypothesis(fit, "log_tau_study_x = 0", component = "scale",
               density_method = "IWMDE"),
    reason,
    fixed = TRUE,
    class = "RoBMA_hypothesis_method"
  )
  kde <- suppressWarnings(hypothesis(
    fit, "intercept = 0", component = "mods", density_method = "KDE"
  ))
  expect_true(is.finite(attr(kde, "raw_BF")))
  # A context built without the capability check stops with the reason.
  expect_error(.iwmde_context(fit), reason, fixed = TRUE)
})


test_that("factor levels of every contrast have point tests and level contrasts", {

  skip_on_cran()
  fits <- .plan_fits()
  for (name in c("treatment", "meandif", "orthonormal", "independent")) {
    quantities <- hypothesis_quantities(fits[[name]])
    terms <- quantities[quantities[["term"]] %in% c("g1", "g2"), , drop = FALSE]
    expect_true(all(terms[["point_test"]]), info = name)
    expect_true(all(terms[["contrast_test"]]), info = name)
    expect_identical(unique(terms[["point_test_methods"]]), "KDE, qCMDE, IWMDE", info = name)
    expect_identical(unique(terms[["contrast_test_methods"]]), "KDE, qCMDE, IWMDE", info = name)
    expect_identical(unique(terms[["reason"]]), "", info = name)
  }
})


test_that("statements on every level of a term with a fixed level are refused as fixed", {

  skip_on_cran()
  fits <- .plan_fits()
  # A whole-term point or region event includes the treatment reference level,
  # which the contrast fixes at 0 (its region prior mass is 0).
  for (statement in c("g1 = 0", "g1 > 0")) {
    expect_error(
      suppressWarnings(hypothesis(fits[["treatment"]], statement, density_method = "KDE")),
      "The quantity 'g1[5]' is fixed by the fitted model; posterior hypothesis tests are undefined.",
      fixed = TRUE,
      class = "RoBMA_hypothesis_fixed",
      info  = statement
    )
  }
  # Comparisons with the other levels remain defined, and a term without a
  # fixed level keeps its whole-term region test.
  expect_s3_class(
    suppressWarnings(hypothesis(fits[["treatment"]], "g1[10] > g1[5]", density_method = "KDE")),
    "BayesTools_hypothesis_BF"
  )
  expect_s3_class(
    suppressWarnings(hypothesis(fits[["meandif"]], "g1 > 0", density_method = "KDE", seed = 1)),
    "BayesTools_hypothesis_BF"
  )
})


test_that("mean-difference level point hypotheses follow the exact Savage-Dickey computation", {

  skip_on_cran()
  fit  <- .plan_fits()[["meandif"]]
  mcmc <- as.matrix(fit[["fit"]][["mcmc"]])
  # The level 'g1[10]' is a linear combination of the two mean-difference
  # coordinates with the contrast row of level 10.
  weights <- BayesTools::contr.meandif(3)[2L, ]
  draws   <- as.numeric(mcmc[, c("mu_g1[1]", "mu_g1[2]")] %*% weights)

  # Prior ordinate: the canonical BayesTools prior density of the catalog
  # level, and analytically the normal with the level's variance under the
  # independent mNormal coordinates.
  catalog   <- BayesTools::parameter_catalog(fit[["fit"]])
  selection <- BayesTools::parameter_catalog_resolve(catalog, alias = "mu_g1[10]")
  prior_ordinate <- exp(BayesTools::prior_density_ordinate(
    BayesTools::parameter_prior_density(fit[["fit"]], selection), 0
  )[["log_density"]])
  prior_sd <- fit[["priors"]][["mods"]][["g1"]][["parameters"]][["sd"]]
  expect_equal(prior_ordinate, stats::dnorm(0, 0, prior_sd * sqrt(sum(weights^2))),
               tolerance = 1e-10)

  # The KDE route: the prior ordinate over the exact Gaussian kernel sum of
  # the level draws at the null (bandwidth bw.nrd0).
  bandwidth <- stats::bw.nrd0(draws)
  kernel    <- stats::dnorm(0, mean = draws, sd = bandwidth)
  kde <- suppressWarnings(hypothesis(fit, "g1[10] = 0", density_method = "KDE",
                                     columns = "all"))
  expect_equal(as.numeric(kde[["prior"]]), prior_ordinate, tolerance = 1e-10)
  expect_equal(attr(kde, "raw_BF"), prior_ordinate / mean(kernel), tolerance = 1e-10)

  # The normal approximation: the prior ordinate over the normal density of
  # the level draws.
  normal <- suppressWarnings(hypothesis(fit, "g1[10] = 0", density_method = "normal",
                                        columns = "all"))
  expect_equal(as.numeric(normal[["prior"]]), prior_ordinate, tolerance = 1e-10)
  expect_equal(
    as.numeric(normal[["posterior"]]),
    stats::dnorm(0, mean(draws), stats::sd(draws)),
    tolerance = 1e-6
  )
  # Long-run reference: under the normal approximation the kernel sum
  # estimates the normal density widened by the bandwidth; the kernel sum
  # agrees with it within four Monte Carlo standard errors of its terms.
  expected_kernel <- stats::dnorm(0, mean(draws), sqrt(stats::var(draws) + bandwidth^2))
  kernel_mcse     <- stats::sd(kernel) / sqrt(coda::effectiveSize(kernel))
  expect_lt(abs(mean(kernel) - expected_kernel), 4 * kernel_mcse)
})


test_that("levels of mean-difference multivariate t factors have exact point tests", {

  skip_on_cran()
  # The multivariate t prior mt(0, s^2 I, nu) of the mean-difference
  # coordinates b of a factor. A level a' b is univariate
  # t(0, s ||a||, nu) with the prior's degrees of freedom (BayesTools
  # 0.3.1.127): its point hypotheses have that exact prior ordinate.
  set.seed(1)
  k <- 48L
  data <- data.frame(
    g1  = factor(rep(c(5, 10, 20), length.out = k), levels = c(5, 10, 20)),
    sei = stats::runif(k, 0.1, 0.3)
  )
  data[["yi"]] <- stats::rnorm(
    k, c(0, 0.2, 0.4)[as.integer(data[["g1"]])], data[["sei"]]
  )
  fit <- suppressWarnings(brma(
    yi = yi, sei = sei, mods = ~ g1, data = data, measure = "SMD",
    set_contrast_factor_predictors = "meandif",
    prior_mods = list(g1 = BayesTools::prior_factor(
      "mt", list(location = 0, scale = 0.5, df = 3), contrast = "meandif"
    )),
    chains = 1, sample = 1000, burnin = 200, adapt = 100, seed = 1,
    silent = TRUE
  ))
  prior <- fit[["priors"]][["mods"]][["g1"]]
  expect_identical(prior[["distribution"]], "mt")

  quantities <- hypothesis_quantities(fit)
  levels <- quantities[quantities[["term"]] == "g1", , drop = FALSE]
  expect_true(all(levels[["point_test"]]))
  expect_true(all(levels[["contrast_test"]]))
  expect_identical(unique(levels[["point_test_methods"]]), "KDE, qCMDE, IWMDE")
  expect_identical(unique(levels[["reason"]]), "")

  # The level '10' has the contrast row a of contr.meandif(3).
  weights  <- BayesTools::contr.meandif(3)[2L, ]
  scale    <- prior[["parameters"]][["scale"]] * sqrt(sum(weights^2))
  df       <- prior[["parameters"]][["df"]]
  expected <- stats::dt(0 / scale, df = df) / scale
  catalog   <- BayesTools::parameter_catalog(fit[["fit"]])
  selection <- BayesTools::parameter_catalog_resolve(catalog, alias = "mu_g1[10]")
  ordinate  <- BayesTools::prior_density_ordinate(
    BayesTools::parameter_prior_density(fit[["fit"]], selection), 0
  )
  expect_true(ordinate[["exact"]])
  expect_equal(exp(ordinate[["log_density"]]), expected, tolerance = 1e-10)

  # The KDE Bayes factor: that ordinate over the exact Gaussian kernel sum of
  # the level draws at 0 (bandwidth bw.nrd0).
  mcmc      <- as.matrix(fit[["fit"]][["mcmc"]])
  draws     <- as.numeric(mcmc[, c("mu_g1[1]", "mu_g1[2]")] %*% weights)
  kernel    <- stats::dnorm(0, mean = draws, sd = stats::bw.nrd0(draws))
  kde <- suppressWarnings(hypothesis(fit, "g1[10] = 0", density_method = "KDE",
                                     columns = "all"))
  expect_equal(as.numeric(kde[["prior"]]), expected, tolerance = 1e-10)
  expect_equal(attr(kde, "raw_BF"), expected / mean(kernel), tolerance = 1e-10)
})


test_that("a two-level ordered factor tests its level and the level contrast alike", {

  skip_on_cran()
  fit <- .plan_fits()[["ordered_2"]]
  # The level 'g[b]' is the ordered total, N(0, 1), and 'g[a]' is fixed at 0.
  statements <- c("g[b] = 0", "g[b] = g[a]", "g[b] - g[a] = 0")
  for (method in c("KDE", "qCMDE")) {
    results <- lapply(statements, function(statement) {
      suppressWarnings(hypothesis(
        fit, statement, density_method = method, columns = "all", seed = 1,
        density_control = if (identical(method, "qCMDE")) list(n_points = 20, samples = 100)
      ))
    })
    for (result in results) {
      expect_equal(as.numeric(result[["prior"]]), stats::dnorm(0), tolerance = 1e-10, info = method)
      expect_equal(attr(result, "raw_BF"), attr(results[[1L]], "raw_BF"),
                   tolerance = 1e-8, info = method)
    }
  }
  expect_error(
    suppressWarnings(hypothesis(fit, "g[a] = 0", density_method = "KDE")),
    "The quantity 'g[a]' is fixed by the fitted model; posterior hypothesis tests are undefined.",
    fixed = TRUE,
    class = "RoBMA_hypothesis_fixed"
  )
})


test_that("ordered levels with a Beta share have an infinite prior ordinate at zero", {

  skip_on_cran()
  fit <- .plan_fits()[["ordered_3"]]
  # 'g[b] - g[a]' is the total T ~ N(0, 1) times the share w ~ Beta(1, 1):
  # f(x) = int_0^1 phi(x / w) / w dw = int_|x|^Inf phi(t) / t dt, which
  # diverges like -log|x| at 0.
  for (method in c("KDE", "qCMDE", "IWMDE")) {
    expect_error(
      suppressWarnings(hypothesis(fit, "g[b] = g[a]", density_method = method)),
      "Prior density at point hypothesis 'g[b] - g[a] = 0' is infinite, so the Savage-Dickey density ratio is undefined.",
      fixed = TRUE,
      class = "BayesTools_infinite_ordinate",
      info  = method
    )
  }
  density_01 <- stats::integrate(
    function(t) stats::dnorm(t) / t, lower = 0.1, upper = Inf, rel.tol = 1e-12
  )[["value"]]
  off_null <- suppressWarnings(hypothesis(fit, "g[b] - g[a] = 0.1",
                                          density_method = "KDE", columns = "all"))
  expect_equal(as.numeric(off_null[["prior"]]), density_01, tolerance = 1e-6)
  # 'g[c]' is the total itself: an exact normal ordinate.
  total <- suppressWarnings(hypothesis(fit, "g[c] = g[a]", density_method = "KDE",
                                       columns = "all"))
  expect_equal(as.numeric(total[["prior"]]), stats::dnorm(0), tolerance = 1e-10)
})


test_that("a random-effect variance and its standard deviation give the same Bayes factor", {

  object <- gated_random_object()
  for (value in c(0.3, 0.5)) {
    sd  <- hypothesis(object, paste0("`(mu) study: tau(intercept)` = ", value),
                      density_method = "KDE", columns = "all", seed = 1)
    var <- hypothesis(object, paste0("`(mu) study: tau2(intercept)` = ", value^2),
                      density_method = "KDE", columns = "all", seed = 1)
    expect_equal(attr(var, "raw_BF"), attr(sd, "raw_BF"), tolerance = 1e-12)
    # The variance densities are the SD densities over the derivative 2 * sd.
    expect_equal(as.numeric(var[["prior"]]), as.numeric(sd[["prior"]]) / (2 * value),
                 tolerance = 1e-12)
    expect_equal(as.numeric(var[["posterior"]]), as.numeric(sd[["posterior"]]) / (2 * value),
                 tolerance = 1e-12)
    expect_identical(var[["Null"]], paste0("(mu) study: tau2(intercept) = ", value^2))
  }
  # A statement comparing a variance point with a region has the Bayes factor
  # of the point against the encompassing model (evaluated through the SD, as
  # the point statement) over that of the region against the encompassing
  # model (on the variance draws, as the region statement); the inverse with
  # the region on the left.
  tau2 <- "`(mu) study: tau2(intercept)`"
  run <- function(statement) {
    hypothesis(object, statement, density_method = "KDE", columns = "all", seed = 1)
  }
  point    <- run(paste0(tau2, " = 0.09 vs ", tau2, " != 0.09"))
  region   <- run(paste0(tau2, " > 0.09 vs ", tau2, " >= 0"))
  mixed    <- run(paste0(tau2, " = 0.09 vs ", tau2, " > 0.09"))
  reversed <- run(paste0(tau2, " > 0.09 vs ", tau2, " = 0.09"))
  expect_equal(attr(mixed, "raw_BF"), attr(point, "raw_BF") / attr(region, "raw_BF"),
               tolerance = 1e-12)
  expect_equal(attr(reversed, "raw_BF"), attr(region, "raw_BF") / attr(point, "raw_BF"),
               tolerance = 1e-12)
  expect_identical(mixed[["Alternative"]], "(mu) study: tau2(intercept) = 0.09")
  expect_identical(mixed[["Null"]], "(mu) study: tau2(intercept) > 0.09")
  expect_identical(mixed[["method"]], "transitive Savage-Dickey")
  # The variance draws are the squared SD draws, so the statement equals the
  # same statement on the SD.
  sd_mixed <- run("`(mu) study: tau(intercept)` = 0.3 vs `(mu) study: tau(intercept)` > 0.3")
  expect_equal(attr(mixed, "raw_BF"), attr(sd_mixed, "raw_BF"), tolerance = 1e-12)
  # Mixed, point and region statements in one call keep their order.
  combined <- run(c(
    paste0(tau2, " > 0.09"),
    paste0(tau2, " = 0.09 vs ", tau2, " > 0.09"),
    paste0(tau2, " = 0.09")
  ))
  expect_equal(
    attr(combined, "raw_BF"),
    c(attr(run(paste0(tau2, " > 0.09")), "raw_BF"), attr(mixed, "raw_BF"),
      1 / attr(point, "raw_BF")),
    tolerance = 1e-12
  )
  # Draws conditional on the inclusion of the gated component: the same
  # identity on the conditional draws.
  run_conditional <- function(statement) {
    hypothesis(object, statement, density_method = "KDE", columns = "all",
               seed = 1, conditional = TRUE)
  }
  expect_equal(
    attr(run_conditional(paste0(tau2, " = 0.09 vs ", tau2, " > 0.09")), "raw_BF"),
    attr(run_conditional(paste0(tau2, " = 0.09 vs ", tau2, " != 0.09")), "raw_BF") /
      attr(run_conditional(paste0(tau2, " > 0.09 vs ", tau2, " >= 0")), "raw_BF"),
    tolerance = 1e-12
  )
})


test_that("point hypotheses with an inexact prior ordinate are refused with its class", {

  object <- nested_allocation_random_object()
  quantities <- hypothesis_quantities(object)
  nested <- quantities[quantities[["parameter"]] == "(mu) study: tau(x)", , drop = FALSE]
  expect_false(unique(nested[["point_test"]]))
  expect_true(unique(nested[["direction_test"]]))
  expect_match(unique(nested[["reason"]]), "nested variance allocation", fixed = TRUE)
  refusal <- tryCatch(
    hypothesis(object, "`(mu) study: tau(x)` = 0.4", density_method = "KDE"),
    error = function(error) error
  )
  expect_s3_class(refusal, "BayesTools_inexact_ordinate")
  expect_match(conditionMessage(refusal), "nested variance allocation", fixed = TRUE)
  # Its prior density is a plotting grid without deterministic provenance, so
  # region tests take the prior probability from prior draws.
  expect_s3_class(
    suppressWarnings(hypothesis(object, "`(mu) study: tau(x)` > 0.4", density_method = "KDE",
                                seed = 1, n_samples = 1000)),
    "BayesTools_hypothesis_BF"
  )
})


test_that("region prior probabilities come from prior draws only for densities without provenance", {

  # The provenance signal is BayesTools::prior_density_has_provenance(): the
  # plotting grid of a nested-allocation SD has none, so its region prior
  # probabilities come from prior draws; the ungated SD has its exact prior.
  nested <- .hypothesis_plans(
    nested_allocation_random_object(), "`(mu) study: tau(x)` > 0.4"
  )[[1L]]
  expect_false(BayesTools::prior_density_has_provenance(nested[["prior_density"]]))
  expect_true(nested[["prior_draws"]])
  ungated <- .hypothesis_plans(
    single_sd_random_object(BayesTools::prior("gamma", list(2, 2))),
    "`(mu) tau(intercept)` > 0.4"
  )[[1L]]
  expect_false(ungated[["prior_draws"]])
  # A combination without a structural ordinate route (three gamma terms) has
  # provenance: BayesTools evaluates its region probabilities, although its
  # ordinate method is "unsupported_provenance" like that of a grid without
  # provenance. The method field does not route regions.
  gammas <- BayesTools:::.prior_linear_combination_density(
    list(
      x = BayesTools::prior("gamma", list(2, 1)),
      y = BayesTools::prior("gamma", list(2, 1)),
      z = BayesTools::prior("gamma", list(2, 1))
    ),
    c(x = 1, y = 1, z = 1)
  )
  expect_identical(BayesTools::prior_density_ordinate(gammas, 6)[["method"]],
                   "unsupported_provenance")
  expect_true(.hypothesis_plan_density_has_provenance(gammas))
  expect_false(.hypothesis_plan_density_has_provenance(NULL))
})


test_that("level contrasts of model-averaged factors condition on the included term", {

  skip_on_cran()
  fit <- .plan_fits()[["model_averaged"]]
  # The levels share the null component of the averaged prior, so their
  # contrast has an atom at 0 unless the models without the term are left
  # out.
  refusal <- tryCatch(
    suppressWarnings(hypothesis(fit, "g1[10] = g1[5]", density_method = "KDE")),
    error = function(error) error
  )
  expect_s3_class(refusal, "BayesTools_linear_target_unavailable")
  expect_match(conditionMessage(refusal), "conditional = TRUE", fixed = TRUE)
  conditional <- suppressWarnings(hypothesis(
    fit, "g1[10] = g1[5]", density_method = "KDE", conditional = TRUE
  ))
  expect_true(is.finite(attr(conditional, "raw_BF")))
  # A single level off its atom keeps the continuous ordinate.
  level <- suppressWarnings(hypothesis(fit, "g1[10] = 0.1", density_method = "KDE"))
  expect_true(is.finite(attr(level, "raw_BF")))
})


test_that("levels of model-averaged mean-difference factors have exact point tests", {

  skip_on_cran()
  fits <- .plan_fits()
  # The prior of 'g1' mixes a null spike and mNormal(0, sd) with equal
  # weights. A mean-difference level is its contrast row times the
  # coordinates; the rows of contr.meandif(3) have unit norm, so every level
  # has the continuous prior ordinate (1/2) * dnorm(x, 0, sd) off its atom at
  # 0. The level '5' is a single coordinate (row (0, 1)), '10' and '20'
  # combine both. The KDE posterior ordinate is the continuous mass (the
  # nonzero level draws, from the models with the term) times the Gaussian
  # kernel sum of those draws.
  design <- BayesTools::contr.meandif(3)
  cases  <- data.frame(level = c("5", "10", "20"), row = 1:3, value = c(-0.2, 0.05, 0.15),
                       stringsAsFactors = FALSE)
  for (name in c("model_averaged", "robma_mixture")) {
    fit   <- fits[[name]]
    prior <- fit[["priors"]][["mods"]][["g1"]]
    expect_identical(attr(prior, "components"), c("null", "alternative"), info = name)
    expect_equal(attr(prior, "prior_weights"), c(1, 1), info = name)
    sd   <- prior[[2L]][["parameters"]][["sd"]]
    mcmc <- as.matrix(fit[["fit"]][["mcmc"]])
    quantities <- hypothesis_quantities(fit)
    rows <- quantities[quantities[["term"]] %in% "g1", , drop = FALSE]
    expect_true(all(rows[["point_test"]]), info = name)
    expect_identical(unique(rows[["point_test_methods"]]), "KDE, qCMDE, IWMDE", info = name)
    expect_false(any(rows[["contrast_test"]]), info = name)
    for (i in seq_len(nrow(cases))) {
      statement <- paste0("g1[", cases[["level"]][[i]], "] = ", cases[["value"]][[i]])
      info      <- paste(name, statement)
      draws    <- as.numeric(mcmc[, c("mu_g1[1]", "mu_g1[2]")] %*% design[cases[["row"]][[i]], ])
      included <- draws[draws != 0]
      prior_ordinate <- 0.5 * stats::dnorm(cases[["value"]][[i]], 0,
                                           sd * sqrt(sum(design[cases[["row"]][[i]], ]^2)))
      posterior_ordinate <- length(included) / length(draws) * mean(stats::dnorm(
        cases[["value"]][[i]], included, stats::bw.nrd0(included)
      ))
      kde <- suppressWarnings(hypothesis(fit, statement, density_method = "KDE",
                                         columns = "all"))
      expect_equal(as.numeric(kde[["prior"]]), prior_ordinate, tolerance = 1e-10, info = info)
      expect_equal(as.numeric(kde[["posterior"]]), posterior_ordinate, tolerance = 1e-10,
                   info = info)
      expect_equal(attr(kde, "raw_BF"), prior_ordinate / posterior_ordinate,
                   tolerance = 1e-10, info = info)
    }
    # qCMDE uses the same exact prior ordinate for the single-coordinate level.
    qcmde <- suppressWarnings(hypothesis(
      fit, "g1[5] = -0.2", density_method = "qCMDE", columns = "all", seed = 1,
      density_control = list(n_points = 20, samples = 50)
    ))
    expect_equal(as.numeric(qcmde[["prior"]]), 0.5 * stats::dnorm(-0.2, 0, sd),
                 tolerance = 1e-10, info = name)
    expect_true(is.finite(attr(qcmde, "raw_BF")), info = name)
  }
})
