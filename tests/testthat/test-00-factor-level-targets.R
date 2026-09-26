context("Factor-level hypothesis targets and plots")

# Level labels are selectors, never JAGS coordinate positions. The labels
# 5/10/20 differ from every coordinate name, and the labels 1:4 coincide with
# the coordinate names of other levels (treatment level 3 is coordinate 2).
# Small single-chain fits suffice: every expectation below is an identity on
# the same posterior draws.
.factor_level_target_cache <- new.env(parent = emptyenv())

.factor_level_target_fits <- function() {

  if (!is.null(.factor_level_target_cache[["fits"]])) {
    return(.factor_level_target_cache[["fits"]])
  }

  set.seed(1)
  k    <- 48L
  data <- data.frame(
    g1  = factor(rep(c(5, 10, 20), length.out = k), levels = c(5, 10, 20)),
    g2  = factor(rep(1:4, length.out = k)[sample.int(k)], levels = 1:4),
    sei = stats::runif(k, 0.1, 0.3)
  )
  data[["yi"]] <- stats::rnorm(
    k,
    c(0, 0.2, 0.4)[as.integer(data[["g1"]])] +
      c(0, 0.1, 0.3, 0.5)[as.integer(data[["g2"]])],
    data[["sei"]]
  )
  fits <- lapply(c(treatment = "treatment", meandif = "meandif"), function(contrast) {
    suppressWarnings(brma(
      yi = yi, sei = sei, mods = ~ g1 + g2, data = data, measure = "SMD",
      set_contrast_factor_predictors = contrast,
      chains = 1, sample = 1000, burnin = 200, adapt = 100,
      seed = 1, silent = TRUE
    ))
  })
  .factor_level_target_cache[["fits"]] <- fits

  return(fits)
}


# The weights on the fitted coordinates of the level target that hypothesis()
# plans for a point statement on one level.
.factor_level_target_weights <- function(fit, hypothesis) {

  plan <- .hypothesis_plans(fit, hypothesis)[[1L]]
  expect_identical(plan[["kind"]], "linear")
  expect_identical(plan[["route"]], "levels")

  plan[["targets"]][[1L]][["spec"]][["weights"]]
}


test_that("factor-level point hypotheses use the level's own fitted coordinate", {

  skip_on_cran()
  fit   <- .factor_level_target_fits()[["treatment"]]
  mcmc  <- as.matrix(fit[["fit"]][["mcmc"]])
  cases <- data.frame(
    hypothesis = c("g1[10] = 0", "g1[20] = 0", "g2[2] = 0", "g2[3] = 0", "g2[4] = 0"),
    term       = c("g1", "g1", "g2", "g2", "g2"),
    level      = c("10", "20", "2", "3", "4"),
    coordinate = c("mu_g1[1]", "mu_g1[2]", "mu_g2[1]", "mu_g2[2]", "mu_g2[3]"),
    stringsAsFactors = FALSE
  )
  compared <- character()

  for (i in seq_len(nrow(cases))) {
    # Savage-Dickey ratio of the treatment coordinate: its normal prior
    # density over the exact Gaussian kernel sum of its draws at the null
    # (bandwidth bw.nrd0). The comparison needs a null within three
    # bandwidths of the draws (levels 20 and 4 lie far from zero; their
    # ordinates are kernel tails).
    prior  <- fit[["priors"]][["mods"]][[cases[["term"]][i]]]
    draws  <- mcmc[, cases[["coordinate"]][i]]
    bw     <- stats::bw.nrd0(draws)
    if (0 >= min(draws) - 3 * bw && 0 <= max(draws) + 3 * bw) {
      height    <- mean(stats::dnorm(0, mean = draws, sd = bw))
      direct_BF <- stats::dnorm(
        0,
        prior[["parameters"]][["mean"]],
        prior[["parameters"]][["sd"]]
      ) / height
      result <- suppressWarnings(hypothesis(
        fit, cases[["hypothesis"]][i], density_method = "KDE"
      ))
      expect_equal(
        attr(result, "raw_BF"), direct_BF,
        tolerance = 1e-10, info = cases[["hypothesis"]][i]
      )
      compared <- c(compared, cases[["hypothesis"]][i])
    }
    expect_identical(
      .factor_level_target_weights(fit, cases[["hypothesis"]][i]),
      stats::setNames(1, cases[["coordinate"]][i])
    )
  }
  expect_identical(compared, c("g1[10] = 0", "g2[2] = 0", "g2[3] = 0"))
})


test_that("mean-difference level point hypotheses are linear targets of the contrast coordinates", {

  skip_on_cran()
  fit <- .factor_level_target_fits()[["meandif"]]

  # The first of four mean-difference levels has the design row (0, 1, 0) in
  # floating point: its extraction key is one coordinate with unit weight, the
  # coordinate of the contrast coefficient g2{2}. It is still a level: its
  # target is the design row, as for the other levels.
  quantities <- BayesTools::parameter_catalog(fit[["fit"]])[["quantities"]]
  unit_key   <- quantities[["extraction_key"]][[
    which(quantities[["canonical_name"]] == "mu_g2[1]")
  ]]
  coefficient_key <- quantities[["extraction_key"]][[
    which(quantities[["canonical_name"]] == "mu_g2{2}")
  ]]
  if (identical(unname(as.numeric(unit_key[["weights"]])), 1)) {
    expect_identical(unit_key[["dependencies"]], coefficient_key[["dependencies"]])
  }

  cases <- list(
    list(hypothesis = "g1[10] = 0", term = "g1", level = 2L, n_levels = 3L),
    list(hypothesis = "g2[2] = 0",  term = "g2", level = 2L, n_levels = 4L),
    list(hypothesis = "g2[1] = 0",  term = "g2", level = 1L, n_levels = 4L)
  )
  for (case in cases) {
    contrast <- BayesTools::contr.meandif(case[["n_levels"]])[case[["level"]], ]
    expected <- stats::setNames(
      contrast,
      paste0("mu_", case[["term"]], "[", seq_along(contrast), "]")
    )
    expect_equal(
      .factor_level_target_weights(fit, case[["hypothesis"]]),
      expected[expected != 0],
      tolerance = 1e-12,
      info = case[["hypothesis"]]
    )
    result <- suppressWarnings(hypothesis(fit, case[["hypothesis"]], density_method = "KDE"))
    expect_true(is.finite(attr(result, "raw_BF")), info = case[["hypothesis"]])
  }
  region   <- suppressWarnings(hypothesis(fit, "g1[10] > 0", density_method = "KDE"))
  contrast <- suppressWarnings(hypothesis(fit, "g1[10] = g1[5]", density_method = "KDE"))
  expect_true(is.finite(attr(region, "raw_BF")))
  expect_true(is.finite(attr(contrast, "raw_BF")))
})


test_that("contrast-coefficient selectors stop naming the level-label form", {

  skip_on_cran()
  fits <- .factor_level_target_fits()

  # Mean-difference coefficients are catalog quantities without RoBMA entries.
  # hypothesis() refuses them as a statement to restate on levels.
  for (selector in c("g1{1} = 0", "g1{1} > 0", "(mu) g1{1} = 0", "mu_g1{1} = 0")) {
    expect_error(
      suppressWarnings(hypothesis(fits[["meandif"]], selector, density_method = "KDE")),
      paste0(
        "Hypotheses on factor contrast coefficients such as 'g1{1}' are not ",
        "supported. State them on factor levels by their labels, such as ",
        "'g1[5]'."
      ),
      fixed = TRUE,
      info = selector
    )
    error <- tryCatch(
      suppressWarnings(hypothesis(fits[["meandif"]], selector, density_method = "KDE")),
      error = identity
    )
    expect_identical(
      class(error),
      c("RoBMA_hypothesis_statement", "error", "condition"),
      info = selector
    )
  }
  error <- tryCatch(
    plot(fits[["meandif"]], parameter = "g1{1}", plot_type = "ggplot"),
    error = identity
  )
  expect_identical(
    conditionMessage(error),
    paste0(
      "Factor contrast coefficients such as 'g1{1}' cannot be selected. ",
      "Select factor levels by their labels, such as 'g1[5]', or the whole ",
      "term 'g1'."
    )
  )
  expect_false(inherits(error, "RoBMA_hypothesis_statement"))

  # Treatment levels are the coefficients: BayesTools refuses the contrast
  # selector of a level coordinate and names its level form. hypothesis()
  # re-raises it as a statement to restate on that level, with its classes
  # and fields; plot() keeps BayesTools' condition.
  refusal <- tryCatch(
    suppressWarnings(hypothesis(fits[["treatment"]], "g1{1} = 0", density_method = "KDE")),
    BayesTools_selector_unavailable = function(condition) condition
  )
  expect_s3_class(refusal, "BayesTools_selector_unavailable")
  expect_identical(refusal[["selector"]], "g1{1}")
  expect_identical(refusal[["level"]], "g1[10]")
  bayestools <- tryCatch(
    BayesTools::hypothesis_parse(
      "g1{1} = 0",
      catalog        = .brma_parameter_catalog_metadata(fits[["treatment"]])[["catalog"]],
      simplify_names = TRUE
    ),
    error = identity
  )
  expect_s3_class(bayestools, "BayesTools_selector_unavailable")
  expect_identical(
    class(refusal),
    c("RoBMA_hypothesis_statement", class(bayestools))
  )
  expect_identical(conditionMessage(refusal), conditionMessage(bayestools))
  plotted <- tryCatch(
    plot(fits[["treatment"]], parameter = "g1{1}", plot_type = "ggplot"),
    error = identity
  )
  expect_s3_class(plotted, "BayesTools_selector_unavailable")
  expect_false(inherits(plotted, "RoBMA_hypothesis_statement"))

  # Random-slope quantities such as 'tau(g1{3})' are no fixed contrast
  # coefficients: a parse error, not the coefficient message.
  expect_error(
    suppressWarnings(hypothesis(fits[["treatment"]], "tau(g1{3}) > 0.05")),
    "Could not parse hypothesis expression 'tau(g1{3}) > 0.05'.",
    fixed = TRUE
  )
})


test_that("factor-cell plots draw the selected level for every contrast", {

  skip_on_cran()
  fits  <- .factor_level_target_fits()
  cases <- data.frame(
    parameter = c("g1[10]", "g1[20]", "g2[2]", "g2[3]", "g2[4]"),
    term      = c("g1", "g1", "g2", "g2", "g2"),
    level     = c(2L, 3L, 2L, 3L, 4L),
    stringsAsFactors = FALSE
  )

  for (contrast in names(fits)) {
    fit  <- fits[[contrast]]
    mcmc <- as.matrix(fit[["fit"]][["mcmc"]])
    for (i in seq_len(nrow(cases))) {
      term        <- cases[["term"]][i]
      n_levels    <- if (identical(term, "g1")) 3L else 4L
      coordinates <- paste0("mu_", term, "[", seq_len(n_levels - 1L), "]")
      # Level values from the fitted coordinates and the contrast matrix.
      contrasts <- if (identical(contrast, "treatment")) {
        stats::contr.treatment(n_levels)
      } else {
        BayesTools::contr.meandif(n_levels)
      }
      expected <- as.numeric(
        mcmc[, coordinates, drop = FALSE] %*% contrasts[cases[["level"]][i], ]
      )

      entry <- .brma_parameter_select_entry(
        fit, cases[["parameter"]][i], allow_factor_cells = TRUE
      )
      cell <- .plot_brma_factor_cell_samples(
        object                    = fit,
        entry                     = entry,
        standardized_coefficients = FALSE,
        conditional               = FALSE,
        precomputed               = FALSE
      )
      expect_equal(
        as.numeric(cell[["samples"]][[1L]]), expected,
        tolerance = 1e-10, info = paste(contrast, cases[["parameter"]][i])
      )
    }
  }

  plotted <- suppressWarnings(plot(
    fits[["meandif"]], parameter = "g1[10]", plot_type = "ggplot",
    density_method = "KDE"
  ))
  expect_s3_class(plotted, "ggplot")
})


test_that("repeated hypothesis rows are numbered by statement", {

  skip_on_cran()
  fit <- .factor_level_target_fits()[["treatment"]]

  # Rows of several statements on one quantity use BayesTools' "mu (1)"
  # scheme; "g11" would read like another parameter and "g1.1" collides with
  # the names `[.data.frame` gives duplicated rows.
  single <- suppressWarnings(hypothesis(
    fit, "g1[10] = g1[20]", density_method = "KDE"
  ))
  expect_identical(rownames(single), "g1")

  contrast <- suppressWarnings(hypothesis(
    fit, c("g1[10] = g1[20]", "g1[10] - g1[20] = 0.1"), density_method = "KDE"
  ))
  expect_identical(rownames(contrast), c("g1 (1)", "g1 (2)"))
  expect_identical(
    rownames(contrast[c(1L, 1L, 2L), ]),
    c("g1 (1)", "g1 (1).1", "g1 (2)")
  )

  # The boundary-valued region of level 20 keys its warning to its own row.
  regions <- suppressWarnings(hypothesis(
    fit, c("g1[10] > 0", "g1[20] > 0"), density_method = "KDE"
  ))
  expect_identical(rownames(regions), c("g1 (1)", "g1 (2)"))
  expect_identical(names(attr(regions, "warnings")), "g1 (2)")

  # Statements on another quantity in between keep the statement numbers of
  # the full hypothesis ("g1 (3)", not the second g1 row "g1 (2)").
  groups <- suppressWarnings(hypothesis(
    fit, c("g1[10] > 0", "g2[2] > 0", "g1[20] > 0"), density_method = "KDE"
  ))
  expect_identical(rownames(groups), c("g1 (1)", "g2", "g1 (3)"))
  expect_identical(names(attr(groups, "warnings")), "g1 (3)")
})


test_that("hypothesis_quantities reports point and contrast tests for the levels of every contrast", {

  skip_on_cran()
  fits <- .factor_level_target_fits()

  # Treatment levels are fitted coefficients (the reference level is fixed
  # and not tested itself); contrasts with it are defined.
  treatment <- hypothesis_quantities(fits[["treatment"]])
  treatment <- treatment[treatment[["term"]] %in% c("g1", "g2"), , drop = FALSE]
  expect_true(all(treatment[["point_test"]]))
  expect_true(all(treatment[["contrast_test"]]))
  expect_identical(unique(treatment[["point_test_methods"]]), "KDE, qCMDE, IWMDE")
  expect_false(any(nzchar(treatment[["reason"]])))

  # Mean-difference levels are linear combinations of the contrast
  # coefficients with exact prior ordinates.
  meandif <- hypothesis_quantities(fits[["meandif"]])
  g1 <- meandif[meandif[["alias"]] == "g1", , drop = FALSE]
  expect_true(g1[["point_test"]])
  expect_true(g1[["direction_test"]])
  expect_true(g1[["contrast_test"]])
  expect_identical(g1[["point_test_methods"]], "KDE, qCMDE, IWMDE")
  expect_identical(g1[["contrast_test_methods"]], "KDE, qCMDE, IWMDE")
  expect_identical(g1[["reason"]], "")
  for (level in c("5", "10", "20")) {
    result <- suppressWarnings(hypothesis(
      fits[["meandif"]], paste0("g1[", level, "] = 0"), density_method = "KDE"
    ))
    expect_true(is.finite(attr(result, "raw_BF")), info = level)
  }

  # Ordered coding: level 'mid' is the ordered total times a Dirichlet(1, 1)
  # share and 'hi' the total itself. The share has a log singularity at 0,
  # so the ordinate of 'mid' at 0 (and of the contrasts of adjacent levels)
  # is infinite; 'hi' has the normal ordinate of the total. Point tests of
  # 'mid' off 0 are exact.
  set.seed(3)
  k    <- 48L
  data <- data.frame(
    g   = factor(rep(c("lo", "mid", "hi"), length.out = k), levels = c("lo", "mid", "hi")),
    sei = stats::runif(k, 0.1, 0.3)
  )
  data[["yi"]] <- stats::rnorm(k, c(0, 0.2, 0.3)[as.integer(data[["g"]])], data[["sei"]])
  ordered <- suppressWarnings(brma(
    yi = yi, sei = sei, mods = ~ g, data = data, measure = "SMD",
    prior_mods = list(g = BayesTools::prior_ordered(BayesTools::prior("normal", list(0, 1)))),
    chains = 1, sample = 500, burnin = 100, adapt = 100, seed = 1, silent = TRUE
  ))
  quantities <- hypothesis_quantities(ordered)
  g <- quantities[quantities[["alias"]] == "g", , drop = FALSE]
  expect_true(g[["point_test"]])
  expect_identical(g[["point_test_methods"]], "KDE, qCMDE, IWMDE")
  expect_false(g[["contrast_test"]])
  expect_match(
    g[["reason"]],
    "Prior density at point hypothesis 'mu_g[lo] - mu_g[mid] = 0' is infinite",
    fixed = TRUE
  )
  expect_error(
    suppressWarnings(hypothesis(ordered, "g[mid] = 0", density_method = "KDE")),
    class = "BayesTools_infinite_ordinate"
  )
  for (method in c("KDE", "qCMDE")) {
    for (statement in c("g[mid] = 0.1", "g[hi] = 0")) {
      result <- suppressWarnings(
        hypothesis(ordered, statement, density_method = method)
      )
      expect_true(is.finite(attr(result, "raw_BF")), info = paste(statement, method))
    }
  }

  # With numeric labels, the first increment's coordinate 'mu_g[1]' reads as
  # the reference level '1': hypotheses name the level the hypothesis wrote.
  data[["g"]] <- factor(
    c("1", "2", "3")[as.integer(data[["g"]])],
    levels = c("1", "2", "3")
  )
  numeric_labels <- suppressWarnings(brma(
    yi = yi, sei = sei, mods = ~ g, data = data, measure = "SMD",
    prior_mods = list(g = BayesTools::prior_ordered(BayesTools::prior("normal", list(0, 1)))),
    chains = 1, sample = 500, burnin = 100, adapt = 100, seed = 1, silent = TRUE
  ))
  for (statement in c("g[2] = 0.1", "g[3] = 0.1")) {
    expect_s3_class(
      suppressWarnings(hypothesis(numeric_labels, statement, density_method = "KDE")),
      "data.frame"
    )
  }
  expect_error(
    suppressWarnings(hypothesis(numeric_labels, "g[1] = 0.1", density_method = "KDE")),
    "The quantity 'g[1]' is fixed by the fitted model; posterior hypothesis tests are undefined.",
    fixed = TRUE
  )
})


test_that("formula coefficient routes follow the fitted coefficient transform", {

  # The route of an original-scale coefficient is the map type and support
  # that BayesTools declares, and the coefficient's row of the transform
  # matrix on the fitted coordinates.
  transform <- structure(
    list(
      schema_version    = 2L,
      target_scale      = "original",
      target_names      = c("mu_x", "log_tau_intercept"),
      matrix            = rbind(
        mu_x              = c(mu_x = 0.5, log_tau_intercept = 0),
        log_tau_intercept = c(mu_x = 0, log_tau_intercept = 1)
      ),
      source_transforms = c(mu_x = "identity", log_tau_intercept = "log"),
      output_transforms = list(mu_x = "identity", log_tau_intercept = "exp"),
      targets           = formula_transform_targets(
        c(mu_x = "affine", log_tau_intercept = "exp_affine"),
        c(mu_x = "identity", log_tau_intercept = "exp")
      )
    ),
    class = "BayesTools_formula_coefficient_transform"
  )
  testthat::local_mocked_bindings(
    JAGS_formula_coefficient_transform = function(...) transform,
    .package = "BayesTools"
  )
  route <- function(parameter) {
    .brma_formula_coefficient_route(
      object   = list(fit = structure(list(), class = "BayesTools_fit")),
      selected = list(
        parameter = parameter,
        component = "mods",
        entry     = list(formula_parameter = "mu", role = "fixed_coefficient")
      )
    )
  }

  # The qCMDE/IWMDE parameter spec of a plotted original-scale coefficient:
  # a nonlinear map is refused with the original-scale density-method
  # classes; a structurally fixed coefficient is no density-method cause.
  plot_spec <- function(parameter) {
    tryCatch(
      .plot_brma_formula_parameter_spec(
        object          = list(fit = structure(list(), class = "BayesTools_fit")),
        parameter       = parameter,
        parameter_entry = list(
          formula_parameter = "mu",
          role              = "fixed_coefficient",
          component         = "mods"
        ),
        standardized_coefficients = FALSE
      ),
      error = identity
    )
  }
  original_scale <- c("RoBMA_density_method_original_scale",
                      "RoBMA_density_method_unavailable", "error", "condition")

  affine <- route("mu_x")
  expect_identical(affine[["type"]], "affine")
  expect_identical(affine[["weights"]], c(mu_x = 0.5))
  expect_identical(affine[["support"]], c(-Inf, Inf))
  expect_identical(plot_spec("mu_x"), list(type = "linear", weights = c(mu_x = 0.5)))
  exp_affine <- route("log_tau_intercept")
  expect_identical(exp_affine[["type"]], "exp_affine")
  expect_identical(exp_affine[["support"]], c(0, Inf))
  expect_identical(class(plot_spec("log_tau_intercept")), original_scale)

  transform[["targets"]][["map_type"]][[1L]] <- "unsupported"
  expect_identical(
    route("mu_x")[["reason"]],
    "The fitted nonlinear joint coefficient transform for 'mu_x' is not supported by hypothesis()."
  )
  error <- plot_spec("mu_x")
  expect_identical(class(error), original_scale)
  expect_identical(conditionMessage(error), route("mu_x")[["reason"]])
  transform[["matrix"]]["mu_x", "mu_x"] <- 0
  expect_identical(
    route("mu_x")[["reason"]],
    "The fitted coefficient 'mu_x' is structurally fixed and has no posterior hypothesis route."
  )
  error <- plot_spec("mu_x")
  expect_identical(conditionMessage(error), route("mu_x")[["reason"]])
  expect_false(inherits(error, "RoBMA_density_method_unavailable"))

  # Point values outside the open support of a transformed coefficient are
  # refused, naming the coefficient.
  refusal <- .hypothesis_plan_support_refusal(
    0, c(0, Inf), "transformed coefficient 'log_tau_intercept'"
  )
  expect_identical(
    refusal[["reason"]],
    paste0(
      "Point-null value 0 is outside or on the boundary of the open support ",
      "for transformed coefficient 'log_tau_intercept'."
    )
  )
  expect_identical(refusal[["class"]], c("RoBMA_hypothesis_target", "RoBMA_hypothesis_unavailable"))
  expect_null(.hypothesis_plan_support_refusal(
    0.5, c(0, Inf), "transformed coefficient 'log_tau_intercept'"
  ))
})


test_that("linear-target refusals of BayesTools keep their classes", {

  skip_on_cran()
  fit   <- .factor_level_target_fits()[["meandif"]]
  means <- suppressWarnings(marginal_means(fit, density_method = "KDE"))
  # An unclassed refusal of the statement's form is a statement to restate;
  # a classed BayesTools condition (here stale draw metadata) keeps its
  # classes instead of being taken for a statement error. An unknown level
  # (BayesTools' resolution error) is a statement error with BayesTools'
  # classes, as at the other entry points of hypothesis().
  refusals <- list(
    unclassed = simpleError("A linear target must be a linear combination."),
    classed   = structure(
      class = c("BayesTools_stale_metadata", "error", "condition"),
      list(message = "Draw metadata are stale.", call = NULL)
    ),
    not_found = structure(
      class = c("BayesTools_parameter_not_found",
                "BayesTools_parameter_resolution_error", "error", "condition"),
      list(message = "Hypothesis references unknown level '99' for parameter 'mu_g1'.",
           call = NULL, alias = "mu_g1[99]", available = "mu_g1[5]")
    )
  )
  expected <- list(
    unclassed = "RoBMA_hypothesis_statement",
    classed   = "BayesTools_stale_metadata",
    not_found = c("RoBMA_hypothesis_statement", "BayesTools_parameter_not_found",
                  "BayesTools_parameter_resolution_error")
  )
  for (name in names(refusals)) {
    testthat::local_mocked_bindings(
      hypothesis_linear_target = function(...) stop(refusals[[name]]),
      .package = "BayesTools"
    )
    plan <- .hypothesis_plans(fit, "g1[10] = g1[20]")[[1L]]
    expect_identical(plan[["route"]], "combination", info = name)
    expect_identical(plan[["refusal"]][["class"]], expected[[name]], info = name)
    expect_identical(
      plan[["refusal"]][["reason"]], conditionMessage(refusals[[name]]),
      info = name
    )
    # The marginal-means combination route refuses with the same classes.
    error <- tryCatch(hypothesis(means, "g1[10] = g1[20]"), error = identity)
    expect_identical(class(error), c(expected[[name]], "error", "condition"), info = name)
    expect_identical(conditionMessage(error), conditionMessage(refusals[[name]]), info = name)
  }
})


test_that("level refusals name the level by its selector", {

  skip_on_cran()
  fit <- .factor_level_target_fits()[["treatment"]]

  expect_error(
    suppressWarnings(hypothesis(fit, "g1[7] = 0", density_method = "KDE")),
    class = "BayesTools_parameter_not_found"
  )
  refusal <- tryCatch(
    suppressWarnings(hypothesis(fit, "g1[5] > 0", density_method = "KDE")),
    error = function(error) error
  )
  expect_s3_class(refusal, "RoBMA_hypothesis_fixed")
  expect_identical(
    conditionMessage(refusal),
    "The quantity 'g1[5]' is fixed by the fitted model; posterior hypothesis tests are undefined."
  )
})


test_that("point hypotheses on the fixed reference level give one reason for every method", {

  skip_on_cran()
  fit <- .factor_level_target_fits()[["treatment"]]

  # The treatment reference level g1[5] is fixed at 0 by the contrast: its
  # plan refuses every statement on it, whatever the density method.
  statements <- list(
    "g1[5] = 0",
    "g1[5] = 0.1",
    c("g1[5] = 0", "g1[10] = 0")
  )
  for (statement in statements) {
    info <- paste(statement, collapse = ", ")
    kde  <- tryCatch(
      suppressWarnings(hypothesis(fit, statement, density_method = "KDE")),
      error = conditionMessage
    )
    expect_identical(
      kde,
      "The quantity 'g1[5]' is fixed by the fitted model; posterior hypothesis tests are undefined.",
      info = info
    )
    expect_error(
      suppressWarnings(hypothesis(fit, statement)),
      kde,
      fixed = TRUE,
      info  = info
    )
    for (method in c("qCMDE", "IWMDE")) {
      expect_error(
        suppressWarnings(hypothesis(fit, statement, density_method = method)),
        kde,
        fixed = TRUE,
        info  = paste(info, method)
      )
    }
  }
  expect_error(
    suppressWarnings(hypothesis(fit, "g1[5] = 0", density_method = "KDE")),
    class = "RoBMA_hypothesis_fixed"
  )

  # The other levels keep their qCMDE ordinates, and contrasts with the
  # reference level are defined.
  expect_s3_class(
    suppressWarnings(hypothesis(fit, c("g1[10] = 0", "g1[20] = 0"), seed = 1)),
    "data.frame"
  )
  expect_s3_class(
    suppressWarnings(hypothesis(fit, "g1[10] = g1[5]", density_method = "KDE")),
    "data.frame"
  )
})


test_that("fixed reference-level statements keep one reason under qCMDE for model-averaged fits", {

  skip_on_cran()
  set.seed(1)
  k    <- 48L
  data <- data.frame(
    g1  = factor(rep(c(5, 10, 20), length.out = k), levels = c(5, 10, 20)),
    sei = stats::runif(k, 0.1, 0.3)
  )
  data[["yi"]] <- stats::rnorm(
    k, c(0, 0.2, 0.4)[as.integer(data[["g1"]])], data[["sei"]]
  )
  fit <- suppressWarnings(BMA.norm(
    yi = yi, sei = sei, mods = ~ g1, data = data, measure = "SMD",
    set_contrast_factor_predictors = "treatment",
    chains = 1, sample = 1000, burnin = 200, adapt = 100,
    seed = 1, silent = TRUE
  ))

  # A call mixing the fixed reference level with another level is refused by
  # the fixed level's plan before any qCMDE/IWMDE ordinate is computed.
  message_of <- function(expr) tryCatch(suppressWarnings(expr), error = conditionMessage)
  for (statement in list("g1[5] = 0", c("g1[5] = 0", "g1[10] = 0.1"),
                         c("g1[10] = 0.1", "g1[5] = 0"))) {
    info <- paste(statement, collapse = ", ")
    kde  <- message_of(hypothesis(fit, statement, density_method = "KDE"))
    expect_identical(
      kde,
      "The quantity 'g1[5]' is fixed by the fitted model; posterior hypothesis tests are undefined.",
      info = info
    )
    for (method in c("qCMDE", "IWMDE")) {
      expect_identical(
        message_of(hypothesis(fit, statement, density_method = method)),
        kde,
        info = paste(info, method)
      )
    }
  }
  # A sampled level at its null atom names the inclusion Bayes factor.
  expect_error(
    suppressWarnings(hypothesis(fit, "g1[10] = 0", density_method = "KDE")),
    "inclusion Bayes factor",
    fixed = TRUE,
    class = "BayesTools_point_mass_at_null"
  )
})
