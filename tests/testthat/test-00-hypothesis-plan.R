context("Hypothesis plans: one eligibility source for hypothesis() and hypothesis_quantities()")

source(testthat::test_path("common-functions.R"))

# Small single-chain fits: every expectation is either an identity on the
# same posterior draws, an analytic prior ordinate, or a Monte Carlo
# comparison with its error stated.
.plan_fit_cache <- new.env(parent = emptyenv())

.plan_fits <- function() {

  if (!is.null(.plan_fit_cache[["fits"]])) {
    return(.plan_fit_cache[["fits"]])
  }
  fit <- function(data, mods, prior_mods = NULL, contrast = "treatment",
                  scale = NULL) {
    arguments <- list(
      yi = quote(yi), sei = quote(sei), mods = mods, data = data,
      measure = "SMD", set_contrast_factor_predictors = contrast,
      chains = 1, sample = 1000, burnin = 200, adapt = 100, seed = 1,
      silent = TRUE
    )
    if (!is.null(prior_mods)) {
      arguments[["prior_mods"]] <- prior_mods
    }
    if (!is.null(scale)) {
      arguments[["scale"]] <- scale
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
  # A factor term in the location and the scale formula: its alias 'g1'
  # names a quantity of each component.
  fits[["location_scale"]] <- fit(
    data = factors, mods = ~ g1, scale = ~ g1, contrast = "meandif"
  )
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


# hypothesis_quantities() of a plan fit, computed once per file: the matrix
# test below stores the tables it checks, and later tests reuse them.
.plan_quantities <- function(name) {

  quantities <- .plan_fit_cache[["quantities"]]
  if (is.null(quantities)) {
    quantities <- list()
  }
  if (is.null(quantities[[name]])) {
    quantities[[name]] <- hypothesis_quantities(.plan_fits()[[name]])
    .plan_fit_cache[["quantities"]] <- quantities
  }

  return(quantities[[name]])
}


test_that("hypothesis_quantities() renders the plans that hypothesis() executes", {

  skip_on_cran()
  fits       <- .plan_fits()
  quantities <- list()
  # Every refusal is checked in both profiles. The qCMDE/IWMDE ordinates of
  # admitted statements are numerical computations; they run in the
  # certification profile (case 'iwmde-qcmde').
  for (name in names(fits)) {
    quantities[[name]] <- .expect_plans_consistent(
      fits[[name]], info = name, run_precomputed = is_certification_profile()
    )
  }
  .plan_fit_cache[["quantities"]] <- quantities
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


# A small known-V fit with a factor moderator (for marginal means) and two
# scale formulas (one per random component).
.two_scale_fit_cache <- new.env(parent = emptyenv())

.two_scale_fit <- function() {

  if (!is.null(.two_scale_fit_cache[["fit"]])) {
    return(.two_scale_fit_cache[["fit"]])
  }
  data <- data.frame(
    yi     = c(0.08, 0.13, 0.18, 0.20, 0.01, 0.05),
    study  = rep(c("s1", "s2", "s3"), each = 2L),
    effect = rep(c("a", "b"), 3L),
    x      = c(0, 1, 0, 1, 0, 1),
    g      = factor(c("u", "v", "v", "u", "u", "v"))
  )
  V <- kronecker(diag(3L), matrix(c(0.04, 0.018, 0.018, 0.05), nrow = 2L))
  fit <- suppressWarnings(brma.mv(
    yi = yi, V = V, data = data, measure = "GEN", mods = ~ g,
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
    fixed = TRUE,
    class = "RoBMA_hypothesis_ambiguous"
  )
  # An ambiguous statement is a statement problem that 'component' resolves,
  # not an unavailable test.
  error <- tryCatch(
    hypothesis(fit, "intercept = 0.5", component = "scale"),
    error = identity
  )
  expect_identical(
    class(error),
    c("RoBMA_hypothesis_ambiguous", "RoBMA_hypothesis_statement", "error",
      "condition")
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
    mu_g                     = "mu_g",
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
  classes <- c("RoBMA_density_method_scale_components",
               "RoBMA_density_method_unavailable")
  expect_identical(
    .iwmde_capability(object = fit, density_method = "qCMDE"),
    list(available = FALSE, reason = reason, class = classes)
  )

  # The location intercept and the scale slopes list KDE only, with the
  # reason; the scale intercepts were KDE-only before (exp(affine) targets).
  quantities <- hypothesis_quantities(fit)
  for (parameter in c("mu_intercept", "mu_g", "log_tau_study_x", "log_tau_effect_x")) {
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
  # The hypothesis() refusal also carries the parent class of the
  # capability refusals.
  expect_error(
    hypothesis(fit, "intercept = 0", component = "mods"),
    class = "RoBMA_density_method_unavailable"
  )
  kde <- suppressWarnings(hypothesis(
    fit, "intercept = 0", component = "mods", density_method = "KDE"
  ))
  expect_true(is.finite(attr(kde, "raw_BF")))
  # A context built without the capability check stops with the reason.
  expect_error(.iwmde_context(fit), reason, fixed = TRUE)

  # Outside hypothesis(), every qCMDE/IWMDE entry point stops with the
  # classes of the capability refusal before any density is estimated.
  means   <- marginal_means(fit, density_method = "KDE")
  control <- list(n_points = 20, samples = 50)
  calls   <- list(
    "plot()"             = function(method) {
      plot(fit, "g", density_method = method, density_control = control)
    },
    "marginal_means()"   = function(method) {
      marginal_means(fit, density_method = method, density_control = control)
    },
    "marginal-means plot()" = function(method) {
      plot(means, "g", density_method = method, density_control = control)
    },
    ".iwmde_context()"   = function(method) .iwmde_context(fit)
  )
  for (method in c("qCMDE", "IWMDE")) {
    for (name in names(calls)) {
      for (class in classes) {
        expect_error(
          suppressMessages(calls[[name]](method)),
          class = class,
          info  = paste(name, method)
        )
      }
    }
  }
  expect_error(
    hypothesis(means, "g[u] = 0", density_method = "qCMDE",
               density_control = control),
    class = "RoBMA_density_method_unavailable"
  )
  # hypothesis() refusals, on the fit and on its marginal means, carry the
  # classes of the capability refusal after their own: one cause has the same
  # classes at every entry point.
  refusal_classes <- c(
    "RoBMA_hypothesis_method", "RoBMA_hypothesis_unavailable", classes,
    "error", "condition"
  )
  for (method in c("qCMDE", "IWMDE")) {
    refusals <- list(
      fit   = tryCatch(
        hypothesis(fit, "intercept = 0", component = "mods",
                   density_method = method),
        error = identity
      ),
      means = tryCatch(
        hypothesis(means, "g[u] = 0", density_method = method,
                   density_control = control),
        error = identity
      )
    )
    for (name in names(refusals)) {
      expect_identical(
        class(refusals[[name]]), refusal_classes,
        info = paste(name, method)
      )
    }
  }
})


test_that("factor levels of every contrast have point tests and level contrasts", {

  skip_on_cran()
  for (name in c("treatment", "meandif", "orthonormal", "independent")) {
    quantities <- .plan_quantities(name)
    terms <- quantities[quantities[["term"]] %in% c("g1", "g2"), , drop = FALSE]
    expect_true(all(terms[["point_test"]]), info = name)
    expect_true(all(terms[["contrast_test"]]), info = name)
    expect_identical(unique(terms[["point_test_methods"]]), "KDE, qCMDE, IWMDE", info = name)
    expect_identical(unique(terms[["contrast_test_methods"]]), "KDE, qCMDE, IWMDE", info = name)
    expect_identical(unique(terms[["reason"]]), "", info = name)
  }
})


test_that("levels of a factor alias shared by the location and scale formulas resolve with 'component'", {

  skip_on_cran()
  fit <- .plan_fits()[["location_scale"]]
  # 'g1' names the location term 'mu_g1' and the scale term 'log_tau_g1'.
  # With 'component', a level reference 'g1[<level>]' is the level of that
  # component's term: the statement gives the result of the same statement
  # with the parameter name.
  cases <- list(
    list(component = "mods",     parameter = "mu_g1"),
    list(component = "location", parameter = "mu_g1"),
    list(component = "scale",    parameter = "log_tau_g1")
  )
  for (case in cases) {
    for (statement in c("g1[10] = 0.1", "g1[10] = g1[5]", "g1[20] > 0.1")) {
      info    <- paste(case[["component"]], statement)
      aliased <- suppressWarnings(hypothesis(
        fit, statement, component = case[["component"]],
        density_method = "KDE", seed = 1
      ))
      named   <- suppressWarnings(hypothesis(
        fit, gsub("g1[", paste0(case[["parameter"]], "["), statement, fixed = TRUE),
        density_method = "KDE", seed = 1
      ))
      expect_identical(attr(aliased, "raw_BF"), attr(named, "raw_BF"), info = info)
      expect_true(all(is.finite(attr(aliased, "raw_BF"))), info = info)
    }
  }
  # Without 'component', the level reference names a level of both terms and
  # stays ambiguous (the same statement resolves with 'component' above).
  expect_error(
    suppressWarnings(hypothesis(fit, "g1[10] = 0.1", density_method = "KDE")),
    class = "RoBMA_hypothesis_ambiguous"
  )
  expect_error(
    .hypothesis_brma_select_parameter(fit, "g1[10] = 0.1", component = "auto"),
    class = "RoBMA_hypothesis_ambiguous"
  )
  # The contrast coefficient 'g1{1}' of both formulas: its quantities have
  # no RoBMA parameter entries, so BayesTools' ambiguity is re-raised with
  # its classes after the classes of the ambiguous references.
  error <- tryCatch(
    hypothesis(fit, "g1{1} = 0", density_method = "KDE"),
    error = identity
  )
  expect_identical(
    class(error),
    c("RoBMA_hypothesis_ambiguous", "RoBMA_hypothesis_statement",
      "BayesTools_parameter_ambiguous", "BayesTools_parameter_resolution_error",
      "error", "condition")
  )
  expect_identical(
    conditionMessage(error),
    paste0(
      "Parameter alias 'g1{1}' is ambiguous; use 'namespace' or 'component' ",
      "to select one quantity."
    )
  )
})


test_that("statements that do not match 'component' are component mismatches at every entry point", {

  skip_on_cran()
  fit <- .plan_fits()[["location_scale"]]
  mismatch <- c("RoBMA_hypothesis_statement", "RoBMA_component_mismatch",
                "error", "condition")
  # A reference of one parameter of another component, references of
  # parameters of other components, and a level alias of several parameters
  # none of which belongs to the component.
  cases <- list(
    list(statement = "mu_g1[10] = 0.1",            component = "scale",
         message   = paste0("The 'hypothesis' argument selects component = ",
                            "'mods' but 'component' was set to 'scale'.")),
    list(statement = "mu_g1[10] = log_tau_g1[10]", component = "random",
         message   = "The hypothesis does not resolve to component = 'random'."),
    list(statement = "g1[10] = 0.1",               component = "random",
         message   = "The hypothesis does not resolve to component = 'random'."),
    # References without a level are resolved within 'component' first;
    # names of another component's parameters (one, or several) are
    # mismatches too, not unknown names.
    list(statement = "mu_g1 > 0",                  component = "scale",
         message   = paste0("The 'hypothesis' argument selects component = ",
                            "'mods' but 'component' was set to 'scale'.")),
    list(statement = "log_tau_g1 > 0",             component = "mods",
         message   = paste0("The 'hypothesis' argument selects component = ",
                            "'scale' but 'component' was set to 'mods'.")),
    list(statement = "g1 > 0",                     component = "random",
         message   = "The hypothesis does not resolve to component = 'random'."),
    # A statement that also references a parameter of the component: the
    # name of the other component's parameter is a mismatch too, not an
    # unknown name.
    list(statement = "mu_intercept > log_tau_intercept", component = "mods",
         message   = paste0("The 'hypothesis' argument selects component = ",
                            "'scale' but 'component' was set to 'mods'.")),
    list(statement = "mu_intercept > log_tau_intercept", component = "scale",
         message   = paste0("The 'hypothesis' argument selects component = ",
                            "'mods' but 'component' was set to 'scale'.")),
    # Statements with levels are resolved without 'component'; a reference
    # of another component's parameter next to one of the component is a
    # mismatch too, in either order, and not a reference dropped before
    # evaluation.
    list(statement = "log_tau_g1[10] > mu_g1[10]", component = "mods",
         message   = paste0("The 'hypothesis' argument selects component = ",
                            "'scale' but 'component' was set to 'mods'.")),
    list(statement = "mu_g1[10] > log_tau_g1[10]", component = "mods",
         message   = paste0("The 'hypothesis' argument selects component = ",
                            "'scale' but 'component' was set to 'mods'.")),
    list(statement = "log_tau_g1[10] > mu_g1[10]", component = "scale",
         message   = paste0("The 'hypothesis' argument selects component = ",
                            "'mods' but 'component' was set to 'scale'.")),
    list(statement = "`(log_tau) g1[10]` > `(mu) g1[10]`", component = "mods",
         message   = paste0("The 'hypothesis' argument selects component = ",
                            "'scale' but 'component' was set to 'mods'.")),
    list(statement = "`(mu) g1[10]` > log_tau_intercept", component = "mods",
         message   = paste0("The 'hypothesis' argument selects component = ",
                            "'scale' but 'component' was set to 'mods'.")),
    list(statement = "`(mu) g1[10]` > log_tau_intercept", component = "scale",
         message   = paste0("The 'hypothesis' argument selects component = ",
                            "'mods' but 'component' was set to 'scale'."))
  )
  for (case in cases) {
    error <- tryCatch(
      hypothesis(fit, case[["statement"]], component = case[["component"]],
                 density_method = "KDE"),
      error = identity
    )
    expect_identical(class(error), mismatch, info = case[["statement"]])
    expect_identical(conditionMessage(error), case[["message"]], info = case[["statement"]])
  }
  # plot() and the prior functions raise the same cause without the
  # hypothesis class.
  calls <- list(
    "plot()"        = function() plot(fit, parameter_mods = "g1", component = "scale"),
    "plot_prior()"  = function() plot_prior(fit, parameter_mods = "g1", component = "scale"),
    "print_prior()" = function() print_prior(fit, parameter_scale = "g1", component = "mods")
  )
  for (name in names(calls)) {
    error <- tryCatch(calls[[name]](), error = identity)
    expect_identical(
      class(error), c("RoBMA_component_mismatch", "error", "condition"),
      info = name
    )
  }

  # The alias of a factor term of another component, which also names the
  # term's contrast coefficients, is a mismatch too.
  error <- tryCatch(
    hypothesis(.plan_fits()[["meandif"]], "g1 > 0", component = "scale",
               density_method = "KDE"),
    error = identity
  )
  expect_identical(class(error), mismatch)
  expect_identical(
    conditionMessage(error),
    paste0("The 'hypothesis' argument selects component = 'mods' but ",
           "'component' was set to 'scale'.")
  )
})


test_that("BayesTools resolution errors that escape planning are statement errors", {

  skip_on_cran()
  fit <- .plan_fits()[["location_scale"]]
  # The factor alias 'g1' of the location and the scale formula next to a
  # name that is unknown within an explicit 'component' (a name of the
  # other component, or a display label): the statement is resolved without
  # 'component', where the shared alias is ambiguous, so it stops with
  # BayesTools' refusal of its first unresolved reference (its message,
  # classes, and fields) after the statement class.
  ambiguous <- c("RoBMA_hypothesis_statement", "BayesTools_parameter_ambiguous",
                 "BayesTools_parameter_resolution_error", "error", "condition")
  cases <- list(
    list(statement = "g1 > log_tau_intercept",  component = "mods"),
    list(statement = "g1 > mu_intercept",       component = "scale"),
    list(statement = "g1 > `(mu) intercept`",   component = "mods")
  )
  for (case in cases) {
    error <- tryCatch(
      hypothesis(fit, case[["statement"]], component = case[["component"]],
                 density_method = "KDE"),
      error = identity
    )
    expect_identical(class(error), ambiguous, info = case[["statement"]])
    expect_identical(
      conditionMessage(error),
      paste0("Parameter alias 'g1' is ambiguous; use 'namespace' or ",
             "'component' to select one quantity."),
      info = case[["statement"]]
    )
    expect_identical(error[["alias"]], "g1", info = case[["statement"]])
    expect_true(is.data.frame(error[["candidates"]]), info = case[["statement"]])
  }
  # The other component's name first, or a level of the shared alias (which
  # resolves within the component's term): the name is the unknown name.
  for (statement in c("log_tau_intercept < g1", "g1[10] > log_tau_intercept")) {
    error <- tryCatch(
      hypothesis(fit, statement, component = "mods", density_method = "KDE"),
      error = identity
    )
    expect_identical(class(error), c(.hypothesis_not_found_class(), "error", "condition"),
                     info = statement)
    expect_identical(error[["alias"]], "log_tau_intercept", info = statement)
  }

  # A level of a factor term next to the whole term resolves to one
  # parameter; BayesTools refuses the whole term only when it evaluates the
  # statement on the levels, which is a statement error too.
  error <- tryCatch(
    hypothesis(fit, "g1[10] > mu_g1", component = "mods", density_method = "KDE"),
    error = identity
  )
  expect_identical(class(error), c(.hypothesis_not_found_class(), "error", "condition"))
  expect_identical(
    conditionMessage(error),
    "Hypothesis expression references unknown quantity 'mu_g1'."
  )
  expect_identical(error[["alias"]], "mu_g1")

  # Any BayesTools resolution error that escapes the planning or the
  # evaluation, on fitted objects and on marginal means, keeps its message,
  # classes, and fields after the statement class; a statement error keeps
  # its classes.
  means <- marginal_means(.plan_fits()[["treatment"]], density_method = "KDE")
  calls <- list(
    fit   = function() hypothesis(fit, "mu_intercept > 0", density_method = "KDE"),
    means = function() hypothesis(means, "g1[10] > 0")
  )
  unresolved <- structure(
    class = c("BayesTools_parameter_other", "BayesTools_parameter_resolution_error",
              "error", "condition"),
    list(message = "Unresolved reference.", call = NULL, alias = "z")
  )
  ambiguous_statement <- structure(
    class = c(.hypothesis_ambiguous_class(), "BayesTools_parameter_ambiguous",
              "BayesTools_parameter_resolution_error", "error", "condition"),
    list(message = "Ambiguous reference.", call = NULL)
  )
  check_routes <- function(stage) {
    for (route in names(calls)) {
      info <- paste(stage, route)
      refusal <<- unresolved
      error <- tryCatch(calls[[route]](), error = identity)
      expect_identical(
        class(error),
        c("RoBMA_hypothesis_statement", "BayesTools_parameter_other",
          "BayesTools_parameter_resolution_error", "error", "condition"),
        info = info
      )
      expect_identical(conditionMessage(error), "Unresolved reference.", info = info)
      expect_identical(error[["alias"]], "z", info = info)
      refusal <<- ambiguous_statement
      error <- tryCatch(calls[[route]](), error = identity)
      expect_identical(class(error), class(ambiguous_statement), info = info)
    }
  }
  refusal <- NULL
  local({
    testthat::local_mocked_bindings(
      .hypothesis_plan                = function(...) stop(refusal),
      .hypothesis_plan_marginal_means = function(...) stop(refusal),
      .package = "RoBMA"
    )
    check_routes("planning")
  })
  local({
    testthat::local_mocked_bindings(
      .hypothesis_plan_execute                = function(...) stop(refusal),
      .hypothesis_plan_execute_marginal_means = function(...) stop(refusal),
      .package = "RoBMA"
    )
    check_routes("evaluation")
  })
})


test_that("publication-bias parameters are refused as hypothesis targets", {

  skip_on_cran()
  fit   <- .plan_fits()[["robma_mixture"]]
  error <- tryCatch(
    hypothesis(fit, "PET = 0", density_method = "KDE"),
    error = identity
  )
  expect_identical(
    class(error),
    c("RoBMA_hypothesis_target", "RoBMA_hypothesis_unavailable", "error",
      "condition")
  )
  expect_identical(
    conditionMessage(error),
    "Hypothesis tests for publication-bias parameters are not supported."
  )
})


test_that("catalog quantities without a RoBMA parameter are refused as targets, not as stale fits", {

  skip_on_cran()
  fit    <- .plan_fits()[["robma_mixture"]]
  target <- c("RoBMA_hypothesis_target", "RoBMA_hypothesis_unavailable",
              "error", "condition")
  # Quantities of the publication-bias prior (weight-function coordinates by
  # their interval or index, and the mixture indicator) have the
  # publication-bias refusal.
  for (statement in c("`omega[0,0.025]` > 0.5", "omega[6] > 0.5",
                      "bias_indicator > 0.5")) {
    error <- tryCatch(
      hypothesis(fit, statement, density_method = "KDE"),
      error = identity
    )
    expect_identical(class(error), target, info = statement)
    expect_identical(
      conditionMessage(error),
      "Hypothesis tests for publication-bias parameters are not supported.",
      info = statement
    )
  }
  # Other quantities without a RoBMA parameter (the inclusion indicators of
  # model-averaged parameters) are refused as targets, named as the
  # statement references them.
  cases <- c(
    "mu_g1_indicator > 0.5"     = "mu_g1_indicator",
    "`(mu) g1_indicator` > 0.5" = "(mu) g1_indicator",
    "tau_indicator = 0.5"       = "tau_indicator"
  )
  for (statement in names(cases)) {
    error <- tryCatch(
      hypothesis(fit, statement, density_method = "KDE"),
      error = identity
    )
    expect_identical(class(error), target, info = statement)
    expect_identical(
      conditionMessage(error),
      paste0(
        "Hypothesis tests are unavailable for '", cases[[statement]], "'. ",
        "Use hypothesis_quantities() to list the quantities that ",
        "hypothesis() tests."
      ),
      info = statement
    )
  }
  # The same next to a parameter of an explicit 'component', to which such a
  # quantity does not belong.
  cases <- c(
    "mu_intercept > mu_g1_indicator" = "mu_g1_indicator",
    "mu_g1_indicator < mu_intercept" = "mu_g1_indicator"
  )
  for (statement in names(cases)) {
    error <- tryCatch(
      hypothesis(fit, statement, component = "mods", density_method = "KDE"),
      error = identity
    )
    expect_identical(class(error), target, info = statement)
    expect_identical(
      conditionMessage(error),
      paste0(
        "Hypothesis tests are unavailable for '", cases[[statement]], "'. ",
        "Use hypothesis_quantities() to list the quantities that ",
        "hypothesis() tests."
      ),
      info = statement
    )
  }
  # Only a parameter catalog without any RoBMA entries needs a refit.
  metadata <- .brma_parameter_catalog_metadata(fit)
  metadata[["entries"]] <- metadata[["entries"]][0L, , drop = FALSE]
  error <- tryCatch(
    .hypothesis_brma_select_parameter(fit, "mu_g1 > 0", "auto", metadata),
    error = identity
  )
  expect_identical(
    class(error),
    c("RoBMA_refit_required", "BayesTools_refit_required", "error", "condition")
  )
  expect_match(conditionMessage(error), "Refit the model", fixed = TRUE)
})


test_that("the selection of plot() refuses catalog quantities without a RoBMA parameter as hypothesis() does", {

  skip_on_cran()
  fit    <- .plan_fits()[["robma_mixture"]]
  target <- c("RoBMA_hypothesis_target", "RoBMA_hypothesis_unavailable",
              "error", "condition")
  unavailable <- function(argument, quantity) {
    paste0(
      "The specified ", argument, " '", quantity, "' is unavailable: it is a ",
      "quantity of the fitted model, not a model parameter."
    )
  }
  # The quantities that hypothesis() refuses as targets (publication-bias
  # quantities and inclusion indicators): the same classes at every entry
  # point of the parameter selection, with a message naming the quantity.
  for (quantity in c("omega[0,0.025]", "omega[6]", "bias_indicator",
                     "mu_g1_indicator", "(mu) g1_indicator", "tau_indicator")) {
    tested <- tryCatch(
      hypothesis(fit, paste0("`", quantity, "` > 0.5"), density_method = "KDE"),
      error = identity
    )
    expect_identical(class(tested), target, info = quantity)
    errors <- list(
      select      = tryCatch(.brma_parameter_select_entry(fit, quantity),
                             error = identity),
      plot        = tryCatch(plot(fit, quantity), error = identity),
      diagnostic  = tryCatch(plot_diagnostic(fit, quantity, type = "trace"),
                             error = identity),
      plot_prior  = tryCatch(plot_prior(fit, parameter = quantity), error = identity),
      print_prior = tryCatch(print_prior(fit, parameter = quantity), error = identity)
    )
    for (entry_point in names(errors)) {
      info <- paste(quantity, entry_point)
      expect_identical(class(errors[[entry_point]]), class(tested), info = info)
      expect_identical(conditionMessage(errors[[entry_point]]),
                       unavailable("parameter", quantity), info = info)
    }
  }
  # The message names the selecting argument.
  error <- tryCatch(
    .brma_parameter_select_entry(fit, "mu_g1_indicator", argument = "parameter_mods"),
    error = identity
  )
  expect_identical(class(error), target)
  expect_identical(conditionMessage(error),
                   unavailable("parameter_mods", "mu_g1_indicator"))
  # Only a parameter catalog without any RoBMA entries needs a refit.
  metadata <- .brma_parameter_catalog_metadata(fit)
  metadata[["entries"]] <- metadata[["entries"]][0L, , drop = FALSE]
  testthat::local_mocked_bindings(
    .brma_parameter_catalog_metadata = function(object) metadata,
    .package = "RoBMA"
  )
  for (quantity in c("mu_g1", "mu_g1_indicator")) {
    error <- tryCatch(.brma_parameter_select_entry(fit, quantity), error = identity)
    expect_identical(
      class(error),
      c("RoBMA_refit_required", "BayesTools_refit_required", "error", "condition"),
      info = quantity
    )
    expect_identical(
      conditionMessage(error),
      paste0(
        "Resolved parameter metadata are unavailable. Refit the model with ",
        "the current RoBMA/BayesTools build."
      ),
      info = quantity
    )
  }
})


test_that("display labels of parameters are evaluable with and without 'component'", {

  skip_on_cran()
  fit      <- .plan_fits()[["location_scale"]]
  draws    <- do.call(rbind, lapply(fit[["fit"]][["mcmc"]], as.matrix))
  metadata <- .brma_parameter_catalog_metadata(fit)
  # The posterior draws of a catalog quantity from its declared extraction
  # key: a fitted coordinate, or a factor level's linear weights on the
  # fitted coordinates.
  quantity_draws <- function(name) {
    quantities <- metadata[["catalog"]][["quantities"]]
    key <- quantities[["extraction_key"]][[match(name, quantities[["canonical_name"]])]]
    weights <- if (is.null(key[["weights"]])) 1 else key[["weights"]]
    drop(draws[, key[["dependencies"]], drop = FALSE] %*% weights)
  }
  # A display label (the location and the scale intercept, and a factor
  # level) tests the parameter it labels.
  cases <- list(
    list(label = "(mu) intercept",      parameter = "mu_intercept", component = "mods"),
    list(label = "(log_tau) intercept", parameter = "log_tau_intercept", component = "scale"),
    list(label = "(mu) g1[10]",         parameter = "mu_g1[10]", component = "mods")
  )
  for (case in cases) {
    values <- quantity_draws(case[["parameter"]])
    value  <- signif(unname(stats::quantile(values, 0.3)), 6L)
    inside <- mean(values > value)
    reference <- hypothesis(
      fit, paste0(case[["parameter"]], " > ", value),
      density_method = "KDE", columns = "all", seed = 1
    )
    for (component in c("auto", case[["component"]])) {
      info <- paste(case[["label"]], component)
      out  <- hypothesis(
        fit, paste0("`", case[["label"]], "` > ", value), component = component,
        density_method = "KDE", columns = "all", seed = 1
      )
      # The posterior odds of the region are those of the labelled draws.
      expect_equal(out[["posterior"]], inside / (1 - inside), tolerance = 1e-12,
                   info = info)
      for (column in c("BF", "prior", "posterior")) {
        expect_identical(unclass(out[[column]]), unclass(reference[[column]]),
                         info = paste(info, column))
      }
      expect_identical(out[["Alternative"]],
                       paste0(case[["label"]], " > ", value), info = info)
    }
  }
  # The labels of the scale factor term evaluate its levels too.
  out <- hypothesis(fit, "`(log_tau) g1[10]` > 0", component = "scale",
                    density_method = "KDE", columns = "all", seed = 1)
  inside <- mean(quantity_draws("log_tau_g1[10]") > 0)
  expect_equal(out[["posterior"]], inside / (1 - inside), tolerance = 1e-12)
  # A display label of the component next to a name of another component's
  # parameter is a component mismatch in either order, as a name known only
  # outside the component is.
  mismatch <- c("RoBMA_hypothesis_statement", "RoBMA_component_mismatch",
                "error", "condition")
  cases <- list(
    list(component = "mods", other = "scale", statements = c(
      "`(mu) intercept` > log_tau_intercept",
      "log_tau_intercept < `(mu) intercept`",
      "`(mu) intercept` > `(log_tau) intercept`",
      "`(log_tau) intercept` < `(mu) intercept`"
    )),
    list(component = "scale", other = "mods", statements = c(
      "`(log_tau) intercept` > mu_intercept",
      "mu_intercept < `(log_tau) intercept`"
    ))
  )
  for (case in cases) {
    for (statement in case[["statements"]]) {
      info  <- paste(case[["component"]], statement)
      error <- tryCatch(
        hypothesis(fit, statement, component = case[["component"]],
                   density_method = "KDE"),
        error = identity
      )
      expect_identical(class(error), mismatch, info = info)
      expect_identical(
        conditionMessage(error),
        paste0(
          "The 'hypothesis' argument selects component = '", case[["other"]],
          "' but 'component' was set to '", case[["component"]], "'."
        ),
        info = info
      )
    }
  }
})


test_that("statement refusals of the plans are statement errors, not unavailable tests", {

  skip_on_cran()
  fits      <- .plan_fits()
  gated     <- gated_random_object()
  statement <- c("RoBMA_hypothesis_statement", "error", "condition")
  # Each statement has to be restated; the tests themselves are available.
  cases <- list(
    list(
      object = fits[["treatment"]], component = "auto",
      hypothesis = "2 * mu_intercept = 0",
      message = paste0(
        "Point-null hypotheses require a direct parameter or level ",
        "reference; unsupported point expression in: '2 * mu_intercept = 0'."
      )
    ),
    list(
      object = fits[["meandif"]], component = "auto",
      hypothesis = "mu_g1[10] * mu_g1[20] = 0",
      message = paste0(
        "A linear target must be a linear combination of levels of 'mu_g1' ",
        "and numbers."
      )
    ),
    list(
      object = gated, component = "random",
      hypothesis = "2 * `(mu) esid: tau(intercept)` = 0.3",
      message = paste0(
        "Point-null tests for random-effect quantities require a direct ",
        "scalar parameter reference."
      )
    ),
    list(
      object = gated, component = "random",
      hypothesis = paste0(
        "`(mu) esid: tau2(intercept)` = 0.04 vs ",
        "`(mu) esid: tau2(intercept)` = 0.09"
      ),
      message = paste0(
        "Point hypotheses on random-effect variance 'esid: tau2' are ",
        "evaluated through its standard deviation and must compare one ",
        "point value per statement."
      )
    ),
    # One quantity named by several names in one statement.
    list(
      object = fits[["treatment"]], component = "auto",
      hypothesis = "mu_intercept > 0 & intercept < 1",
      message = paste0(
        "Hypothesis names one quantity by several names ('mu_intercept', ",
        "'intercept'). Use one of them throughout the statement."
      )
    ),
    list(
      object = fits[["treatment"]], component = "mods",
      hypothesis = "g1[10] > mu_g1[10]",
      message = paste0(
        "Hypothesis names one quantity by several names ('g1[10]', ",
        "'mu_g1[10]'). Use one of them throughout the statement."
      )
    )
  )
  for (case in cases) {
    error <- tryCatch(
      hypothesis(case[["object"]], case[["hypothesis"]],
                 component = case[["component"]], density_method = "KDE"),
      error = identity
    )
    expect_identical(class(error), statement, info = case[["hypothesis"]])
    expect_identical(conditionMessage(error), case[["message"]],
                     info = case[["hypothesis"]])
  }
  # The same on marginal means, whose intercept is also named 'mu'.
  means <- marginal_means(fits[["treatment"]], density_method = "KDE")
  error <- tryCatch(hypothesis(means, "mu > 0 & intercept < 1"), error = identity)
  expect_identical(class(error), statement)
  expect_identical(
    conditionMessage(error),
    paste0("Hypothesis names one quantity by several names ('mu', ",
           "'intercept'). Use one of them throughout the statement.")
  )
  # Different levels of one term, and one name repeated, are one statement.
  expect_s3_class(
    suppressWarnings(hypothesis(fits[["treatment"]], "g1[10] > mu_g1[5]",
                                density_method = "KDE", seed = 1)),
    "BayesTools_hypothesis_BF"
  )
})


test_that("unresolved references have the same classes on fits and marginal means", {

  skip_on_cran()
  fit   <- .two_scale_fit()
  means <- marginal_means(fit, density_method = "KDE")
  not_found <- c("RoBMA_hypothesis_statement", "BayesTools_parameter_not_found",
                 "BayesTools_parameter_resolution_error", "error", "condition")
  # An unknown name: the statement class followed by BayesTools' condition,
  # whose message and fields the fit path keeps.
  on_fit <- tryCatch(hypothesis(fit, "foo = 0", density_method = "KDE"), error = identity)
  on_means <- tryCatch(hypothesis(means, "foo > 0"), error = identity)
  expect_identical(class(on_fit), not_found)
  expect_identical(class(on_means), not_found)
  expect_identical(conditionMessage(on_fit), "No public parameter quantity matches 'foo'.")
  expect_identical(on_fit[["alias"]], "foo")
  expect_true("mu_intercept" %in% on_fit[["available"]])
  # The same with an explicit component in which the name is unknown too.
  expect_identical(
    class(tryCatch(
      hypothesis(fit, "foo = 0", component = "mods", density_method = "KDE"),
      error = identity
    )),
    not_found
  )
  # Next to a name of another component's parameter, the unknown name is the
  # one refused.
  beside <- tryCatch(
    hypothesis(fit, "mu_intercept > foo", component = "scale", density_method = "KDE"),
    error = identity
  )
  expect_identical(class(beside), not_found)
  expect_identical(conditionMessage(beside), "No public parameter quantity matches 'foo'.")
  expect_identical(beside[["alias"]], "foo")
  # An unknown level of a factor term.
  level_on_fit   <- tryCatch(hypothesis(fit, "g[w] = 0", density_method = "KDE"), error = identity)
  level_on_means <- tryCatch(hypothesis(means, "g[w] > 0"), error = identity)
  expect_identical(class(level_on_fit), not_found)
  expect_identical(class(level_on_means), not_found)
  expect_identical(
    conditionMessage(level_on_means),
    "Hypothesis references unknown level 'w' for parameter 'mu_g'."
  )
  # An unknown level in a linear combination of levels: on marginal means
  # BayesTools' refusal of the linear target, on the fit BayesTools' refusal
  # of the reference; both with the statement class, BayesTools' classes,
  # and BayesTools' fields.
  for (statement in c("g[u] - g[zz] = 0", "g[zz] = g[u]")) {
    combination <- list(
      fit   = tryCatch(hypothesis(fit, statement, density_method = "KDE"),
                       error = identity),
      means = tryCatch(hypothesis(means, statement), error = identity)
    )
    for (route in names(combination)) {
      error <- combination[[route]]
      info  <- paste(statement, route)
      expect_identical(class(error), not_found, info = info)
      expect_true(is.character(error[["alias"]]) && length(error[["alias"]]) > 0L,
                  info = info)
      expect_true(any(grepl("g", error[["alias"]], fixed = TRUE)), info = info)
      expect_true(is.character(error[["available"]]) &&
                    length(error[["available"]]) > 0L, info = info)
    }
    expect_identical(combination[["means"]][["alias"]], "mu_g[zz]", info = statement)
    expect_identical(
      conditionMessage(combination[["means"]]),
      "Hypothesis references unknown level 'zz' for parameter 'mu_g'.",
      info = statement
    )
  }
  # A statement without references: on the fit BayesTools' refusal (with the
  # classes BayesTools gives it) after the statement class; on marginal
  # means RoBMA's statement error with the same classes, also with
  # 'parameter' and next to a statement with references.
  metadata <- .brma_parameter_catalog_metadata(fit)
  refusal  <- tryCatch(
    BayesTools::hypothesis_resolve(
      BayesTools::hypothesis_parse("1 > 0"), metadata[["catalog"]]
    ),
    error = identity
  )
  expect_identical(
    conditionMessage(refusal),
    "The hypothesis contains no parameter symbols to resolve."
  )
  no_parameters <- unique(c("RoBMA_hypothesis_statement", class(refusal)))
  expect_identical(
    no_parameters,
    c("RoBMA_hypothesis_statement", "BayesTools_hypothesis_no_parameters",
      "BayesTools_parameter_resolution_error", "error", "condition")
  )
  on_fit <- list(
    auto   = tryCatch(hypothesis(fit, "1 > 0", density_method = "KDE"),
                      error = identity),
    scale  = tryCatch(hypothesis(fit, "1 > 0", component = "scale",
                                 density_method = "KDE"),
                      error = identity),
    beside = tryCatch(hypothesis(fit, c("g[u] > 0", "1 > 0"),
                                 density_method = "KDE"),
                      error = identity)
  )
  for (name in names(on_fit)) {
    expect_identical(class(on_fit[[name]]), no_parameters, info = name)
    expect_identical(conditionMessage(on_fit[[name]]), conditionMessage(refusal),
                     info = name)
  }
  on_means <- list(
    plain     = tryCatch(hypothesis(means, "1 > 0"), error = identity),
    parameter = tryCatch(hypothesis(means, "1 > 0", parameter = "g"),
                         error = identity),
    beside    = tryCatch(hypothesis(means, c("g[u] > 0", "1 > 0")),
                         error = identity)
  )
  for (name in names(on_means)) {
    expect_identical(class(on_means[[name]]), no_parameters, info = name)
    expect_identical(conditionMessage(on_means[[name]]),
                     "Hypothesis must reference a marginal-means parameter.",
                     info = name)
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
