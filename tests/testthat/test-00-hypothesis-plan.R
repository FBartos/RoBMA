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


# hypothesis() evaluates a statement exactly when its plan admits the method,
# and otherwise stops with the plan's refusal (its first class and message).
# hypothesis_quantities() renders the same plans. qCMDE/IWMDE statements the
# plans admit are evaluated for the first point and contrast statement of
# each object ('run_precomputed'), with a small density budget.
.expect_plans_consistent <- function(object, info, run_precomputed = TRUE) {

  metadata   <- .brma_parameter_catalog_metadata(object)
  quantities <- hypothesis_quantities(object)
  cache      <- .hypothesis_plan_cache()
  control    <- list(n_points = 20, samples = 50)
  run <- function(statement, component, method) {
    tryCatch(
      suppressWarnings(hypothesis(
        object, statement, component = component, density_method = method,
        density_control = if (method %in% c("qCMDE", "IWMDE")) control,
        seed = 1, n_samples = 500
      )),
      error = function(error) error
    )
  }
  ran_precomputed <- character()
  entries <- metadata[["entries"]]
  entries <- entries[entries[["component"]] != "bias", , drop = FALSE]
  for (i in seq_len(nrow(entries))) {
    entry <- as.list(entries[i, setdiff(names(entries), "aliases"), drop = FALSE])
    rows  <- quantities[quantities[["parameter"]] == entry[["parameter"]] &
                          quantities[["component"]] == entry[["component"]], ,
                        drop = FALSE]
    plans <- .hypothesis_quantities_plans(object, entry, metadata, cache)
    rendered <- .hypothesis_quantities_render_plans(
      plans   = plans,
      bracket = identical(entry[["role"]], "formula_coefficient_group")
    )
    row_info <- paste(info, entry[["parameter"]])
    expect_gt(nrow(rows), 0L)
    for (column in c("point_test", "direction_test", "contrast_test",
                     "point_test_methods", "contrast_test_methods", "reason")) {
      expect_identical(
        unique(rows[[column]]), rendered[[column]],
        info = paste(row_info, column)
      )
    }
    for (type in c("point", "region", "contrast")) {
      for (plan in plans[[type]]) {
        statement <- BayesTools::hypothesis_render(plan[["statement"]])
        for (method in .hypothesis_plan_advertised_methods()) {
          refusal <- .hypothesis_plan_status(plan, method)
          case    <- paste(row_info, statement, method)
          precomputed <- method %in% c("qCMDE", "IWMDE")
          if (is.null(refusal) && precomputed &&
              (!run_precomputed || paste(type, method) %in% ran_precomputed ||
                 identical(type, "region"))) {
            next
          }
          out <- run(statement, entry[["component"]], method)
          if (is.null(refusal)) {
            # An admitted statement runs; qCMDE/IWMDE ordinates may still be
            # rejected by their numerical diagnostics.
            expect_false(
              inherits(out, "error") &&
                !(precomputed && inherits(out, "RoBMA_density_ordinate_error")),
              info = paste(case, if (inherits(out, "error")) conditionMessage(out))
            )
            if (precomputed) {
              ran_precomputed <- c(ran_precomputed, paste(type, method))
            }
          } else {
            expect_s3_class(out, refusal[["class"]][[1L]])
            expect_identical(conditionMessage(out), refusal[["reason"]], info = case)
          }
        }
      }
    }
  }

  invisible(quantities)
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
  # A statement mixing point and region sides is evaluated as a whole.
  expect_error(
    hypothesis(object, "`(mu) study: tau2(intercept)` = 0.09 vs `(mu) study: tau2(intercept)` > 0.09",
               density_method = "KDE"),
    class = "RoBMA_hypothesis_statement"
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
