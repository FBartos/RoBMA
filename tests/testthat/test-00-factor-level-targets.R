source(testthat::test_path("common-functions.R"))

context("Factor-level hypothesis targets and plots")

# Level labels are selectors, never JAGS coordinate positions. The labels
# 5/10/20 differ from every coordinate name, and the labels 1:4 coincide with
# the coordinate names of other levels (treatment level 3 is coordinate 2).
# Small single-chain fits suffice: every expectation below is an identity on
# the same posterior draws.
.factor_level_target_fits <- function() {

  list(
    treatment = load_fit("hypothesis_plan_treatment"),
    meandif = load_fit("hypothesis_plan_meandif")
  )
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
  ordered <- load_fit("lane_factor_ordered_labels")
  data <- load_info("lane_factor_ordered_labels")$data
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
  # Observe the actual public BF input without replacing its computation.
  captured <- NULL
  compute <- BayesTools::hypothesis_BF
  testthat::local_mocked_bindings(
    hypothesis_BF = function(posterior, ...) {
      captured <<- c(list(posterior = posterior), list(...))
      do.call(compute, captured)
    },
    .package = "BayesTools"
  )
  # A finite Bayes factor needs only a small qCMDE grid and row sample.
  control <- list(n_points = 20, samples = 50)
  for (method in c("KDE", "qCMDE")) {
    for (statement in c("g[mid] = 0.1", "g[hi] = 0")) {
      result <- suppressWarnings(hypothesis(
        ordered, statement, density_method = method,
        density_control = if (identical(method, "qCMDE")) control
      ))
      expect_true(is.finite(attr(result, "raw_BF")), info = paste(statement, method))
      if (identical(method, "qCMDE")) {
        level <- if (identical(statement, "g[mid] = 0.1")) "mid" else "hi"
        value <- if (identical(level, "mid")) 0.1 else 0
        .expect_iwmde_child_identity(captured[["posterior"]][[level]], "mu_g", level, value)
        condition <- BayesTools::posterior_metadata(captured[["posterior"]][[level]], "condition")
        expect_true(condition[["averaged"]])
        expect_identical(condition[["condition_key"]], "<averaged>")
        wrong <- captured
        ordinate <- BayesTools::posterior_metadata(wrong[["posterior"]][[level]], "posterior_ordinate")
        ordinate[["parameter"]] <- "mu_g[lo]"
        BayesTools::posterior_metadata(wrong[["posterior"]][[level]], "posterior_ordinate") <- ordinate
        expect_error(do.call(compute, wrong),
          "requires valid posterior ordinate or posterior density metadata", fixed = TRUE)
      }
    }
  }

  # The semantic marginal evaluates ordered totals and allocations before the
  # cell helper selects it. Preserve that value and its complete draw metadata,
  # including conditioning, instead of re-reading displayed contrast columns.
  marginal_posterior <- BayesTools::marginal_posterior
  selected_marginal <- NULL
  testthat::local_mocked_bindings(
    marginal_posterior = function(...) {
      selected_marginal <<- marginal_posterior(...)
      selected_marginal
    },
    .package = "BayesTools"
  )
  testthat::local_mocked_bindings(
    .brma_draws_with_values = function(...) {
      stop("Ordered cells must retain their semantic marginal values.")
    },
    .package = "RoBMA"
  )
  for (selector in c("g[lo]", "g[mid]", "g[hi]")) {
    entry <- .brma_parameter_select_entry(ordered, selector, allow_factor_cells = TRUE)
    expected <- as.numeric(as.matrix(BayesTools::parameter_draws(
      ordered[["fit"]], entry[["selection"]]
    )))
    for (standardized in c(FALSE, TRUE)) {
      cell <- .plot_brma_factor_cell_samples(
        ordered, entry, standardized, conditional = FALSE, precomputed = FALSE
      )
      value <- cell[["samples"]][[1L]]
      marginal_value <- selected_marginal[[attr(value, "level_name", exact = TRUE)]]
      expect_identical(as.numeric(value), as.numeric(marginal_value))
      expect_identical(attr(value, "bayestools_meta", exact = TRUE),
        attr(marginal_value, "bayestools_meta", exact = TRUE))
      for (field in c("atoms", "prior_density", "condition", "linear_weights")) {
        expect_identical(BayesTools::posterior_metadata(value, field),
          BayesTools::posterior_metadata(marginal_value, field), info = field)
      }
      expect_identical(as.numeric(value), expected)
      expect_identical(cell[["structural"]], identical(selector, "g[lo]"))
      if (cell[["structural"]]) {
        expect_identical(BayesTools::posterior_metadata(value, "atoms")[["mass"]], 1)
      }
    }
  }
  grDevices::pdf(NULL)
  tryCatch({
    for (selector in c("g", "g[mid]", "g[hi]")) {
      expect_null(plot(ordered, parameter = selector, prior = TRUE,
        density_method = "KDE", plot_type = "base"))
      plotted <- plot(ordered, parameter = selector, prior = TRUE,
        density_method = "KDE", plot_type = "ggplot")
      expect_s3_class(plotted, "ggplot")
      layers <- ggplot2::ggplot_build(plotted)[["data"]]
      if (selector %in% c("g[mid]", "g[hi]")) {
        expect_length(layers, 2L)
        expect_identical(unique(layers[[1L]][["colour"]]), "black")
        expect_identical(unique(layers[[1L]][["linetype"]]), 2)
        expect_identical(unique(layers[[2L]][["colour"]]), "black")
        expect_identical(unique(layers[[2L]][["linetype"]]), 1)
      }
    }
    styles <- list(
      list(posterior = list(col = "purple"), prior = NULL, color = "purple", lty = 2),
      list(posterior = list(color = "green4"), prior = NULL, color = "green4", lty = 2),
      list(posterior = list(colour = "blue"), prior = NULL, color = "blue", lty = 2),
      list(posterior = list(col = "purple", color = "green4", colour = "blue"),
        prior = list(col = "red", color = "blue", lty = 3), color = "red", lty = 3),
      list(posterior = list(), prior = list(color = "red", colour = "blue", linetype = "dotted"),
        color = "red", lty = "dotted"),
      list(posterior = list(), prior = list(colour = "red", lty = 3, linetype = "dotdash"),
        color = "red", lty = "dotdash")
    )
    for (style in styles) {
      plotted <- do.call(plot, c(list(x = ordered, parameter = "g[hi]",
        prior = TRUE, density_method = "KDE", plot_type = "ggplot",
        dots_prior = style[["prior"]]), style[["posterior"]]))
      layers <- ggplot2::ggplot_build(plotted)[["data"]]
      expect_identical(unique(layers[[1L]][["colour"]]), style[["color"]])
      expect_identical(unique(layers[[1L]][["linetype"]]), style[["lty"]])
      expect_identical(unique(layers[[2L]][["colour"]]),
        if (!length(style[["posterior"]])) "black" else if (is.null(style[["prior"]])) {
          style[["color"]]
        } else "purple")
    }
  }, finally = grDevices::dev.off())
  # BayesTools 0.3.0 warned and used all draws for a parameter that is not
  # conditional (restored in 0.3.1.258): the warning, not an error.
  expect_warning(.plot_brma_factor_cell_samples(ordered, entry, FALSE, TRUE, FALSE),
    "The parameter 'mu_g' is not a conditional parameter.", fixed = TRUE)
  entry <- .brma_parameter_select_entry(ordered, "g[mid]", allow_factor_cells = TRUE)
  expect_error(BayesTools::parameter_draws(ordered[["fit"]], entry[["selection"]],
    model_samples = as.matrix(ordered[["fit"]][["mcmc"]])[, "mu_g[1]", drop = FALSE]),
    class = "BayesTools_ordered_coordinates_unavailable")

  # With numeric labels, the first increment's coordinate 'mu_g[1]' reads as
  # the reference level '1': hypotheses name the level the hypothesis wrote.
  numeric_labels <- load_fit("lane_factor_numeric_labels")
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
  for (selector in c("g[1]", "g[2]", "g[3]")) {
    entry <- .brma_parameter_select_entry(numeric_labels, selector, allow_factor_cells = TRUE)
    cell <- .plot_brma_factor_cell_samples(numeric_labels, entry, FALSE, FALSE, FALSE)
    expect_identical(as.numeric(cell[["samples"]][[1L]]),
      as.numeric(as.matrix(BayesTools::parameter_draws(numeric_labels[["fit"]], entry[["selection"]]))))
  }
})


test_that("formula coefficient routes follow the fitted coefficient transform", {

  # Synthetic compiler/consumer rows, not posterior-inference evidence.
  location_compiled <- BayesTools::JAGS_formula(
    formula = ~ 1 + x, parameter = "mu", data = data.frame(x = c(2, 4, 6)),
    prior_list = list(
      intercept = BayesTools::prior("normal", list(0, 1)),
      x = BayesTools::prior("normal", list(0, 1))
    ), formula_scale = TRUE
  )
  location_fit <- coda::mcmc.list(coda::mcmc(cbind(
    mu_intercept = c(-.2, 0, .2), mu_x = c(-.1, 0, .1)
  )))
  attr(location_fit, "prior_list") <- location_compiled[["prior_list"]]
  attr(location_fit, "formula_design") <- list(mu = location_compiled[["formula_design"]])
  attr(location_fit, "formula_scale") <- list(mu = location_compiled[["formula_scale"]])
  location_fit <- as_bayestools_fit(location_fit)
  formula <- ~ 1 + x
  attr(formula, "log(intercept)") <- TRUE
  compiled <- BayesTools::JAGS_formula(
    formula = formula, parameter = "log_tau",
    data = data.frame(x = c(2, 4, 6)),
    prior_list = list(
      intercept = BayesTools::prior("gamma", list(2, 2)),
      x = BayesTools::prior("normal", list(0, 1))
    ), formula_scale = TRUE
  )
  fit <- coda::mcmc.list(coda::mcmc(cbind(
    log_tau_intercept = c(.2, .3, .4), log_tau_x = c(-.1, 0, .1)
  )))
  attr(fit, "prior_list") <- compiled[["prior_list"]]
  attr(fit, "formula_design") <- list(log_tau = compiled[["formula_design"]])
  attr(fit, "formula_scale") <- list(log_tau = compiled[["formula_scale"]])
  fit <- as_bayestools_fit(fit)

  scale_fit <- fit
  route <- function(parameter, fit = location_fit, formula_parameter = "mu",
                    component = "mods") {
    .brma_formula_coefficient_route(
      object = list(fit = fit),
      selected = list(parameter = parameter, component = component,
        entry = list(formula_parameter = formula_parameter, role = "fixed_coefficient"))
    )
  }
  plot_spec <- function(parameter, fit = location_fit, formula_parameter = "mu",
                        component = "mods") {
    tryCatch(
      .plot_brma_formula_parameter_spec(
        object = list(fit = fit), parameter = parameter,
        parameter_entry = list(formula_parameter = formula_parameter,
          role = "fixed_coefficient", component = component),
        standardized_coefficients = FALSE
      ), error = identity
    )
  }
  original_scale <- c("RoBMA_density_method_original_scale",
                      "RoBMA_density_method_unavailable", "error", "condition")

  affine <- route("mu_x")
  expect_identical(affine[["type"]], "affine")
  expect_identical(affine[["weights"]], c(mu_x = 0.5))
  expect_identical(affine[["support"]], c(-Inf, Inf))
  expect_identical(plot_spec("mu_x"), list(type = "linear", weights = c(mu_x = 0.5)))
  exp_affine <- route("log_tau_intercept", scale_fit, "log_tau", "scale")
  expect_identical(exp_affine[["type"]], "exp_affine")
  expect_identical(exp_affine[["support"]], c(0, Inf))
  expect_identical(class(plot_spec("log_tau_intercept", scale_fit, "log_tau", "scale")),
    original_scale)

  coefficient_transform <- BayesTools::JAGS_formula_coefficient_transform
  transform <- coefficient_transform(location_fit, "mu")
  # Artificial consumer-branch test: modify a copy of a genuine descriptor;
  # this unsupported map is not presented as a valid compiler output.
  transform[["targets"]][["map_type"]][transform[["targets"]][["target"]] == "mu_x"] <- "unsupported"
  testthat::local_mocked_bindings(
    JAGS_formula_coefficient_transform = function(...) transform,
    .package = "BayesTools"
  )
  expect_identical(
    route("mu_x")[["reason"]],
    "The fitted nonlinear joint coefficient transform for 'mu_x' is not supported by hypothesis()."
  )
  error <- plot_spec("mu_x")
  expect_identical(class(error), original_scale)
  expect_identical(conditionMessage(error), route("mu_x")[["reason"]])

  # Genuine all-point compiler metadata establish structural fixedness;
  # constant retained continuous draws do not establish it.
  fixed_compiled <- BayesTools::JAGS_formula(
    formula = ~ 1 + x, parameter = "mu", data = data.frame(x = c(2, 4, 6)),
    prior_list = list(intercept = BayesTools::prior("point", list(3)),
      x = BayesTools::prior("point", list(0))), formula_scale = TRUE
  )
  fixed_design <- fixed_compiled[["formula_design"]]
  fixed_scale <- BayesTools::JAGS_formula_scale_finalize(
    fixed_compiled[["formula_scale"]], fixed_design,
    fixed_compiled[["prior_list"]], owner_scope = "fit"
  )
  fixed_design[["formula_scale"]] <- fixed_scale
  fixed_fit <- structure(list(), class = "BayesTools_fit")
  attr(fixed_fit, "prior_list") <- fixed_compiled[["prior_list"]]
  attr(fixed_fit, "formula_design") <- list(mu = fixed_design)
  attr(fixed_fit, "formula_scale") <- list(mu = fixed_scale)
  attr(fixed_fit, "parameter_map") <- BayesTools:::.bt_build_parameter_map(
    columns = character(), prior_list = fixed_compiled[["prior_list"]],
    formula_design = list(mu = fixed_design), formula_scale = list(mu = fixed_scale)
  )
  fixed_fit <- BayesTools:::.bt_attach_fit_contract(fixed_fit)
  # Use the genuine producer saved before mocking the unsupported branch.
  transform <- coefficient_transform(fixed_fit, "mu")
  expect_identical(
    route("mu_x", fixed_fit)[["reason"]],
    "The fitted coefficient 'mu_x' is structurally fixed and has no posterior hypothesis route."
  )
  error <- plot_spec("mu_x", fixed_fit)
  expect_identical(conditionMessage(error), route("mu_x", fixed_fit)[["reason"]])
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
  # classes instead of being taken for a statement error. BayesTools'
  # resolution errors (an unknown level, a statement without parameter
  # symbols, a level of another component) are statement errors with
  # BayesTools' classes and fields, as at the other entry points of
  # hypothesis().
  resolution <- function(class, message, ...) {
    structure(
      class = c(class, "BayesTools_parameter_resolution_error", "error",
                "condition"),
      list(message = message, call = NULL, ...)
    )
  }
  refusals <- list(
    unclassed = simpleError("A linear target must be a linear combination."),
    classed   = structure(
      class = c("BayesTools_stale_metadata", "error", "condition"),
      list(message = "Draw metadata are stale.", call = NULL)
    ),
    not_found = resolution(
      "BayesTools_parameter_not_found",
      "Hypothesis references unknown level '99' for parameter 'mu_g1'.",
      alias = "mu_g1[99]", available = "mu_g1[5]"
    ),
    no_parameters = resolution(
      "BayesTools_hypothesis_no_parameters",
      "The hypothesis contains no parameter symbols to resolve."
    ),
    component_mismatch = resolution(
      "BayesTools_hypothesis_component_mismatch",
      "The level in hypothesis symbol 'g1[10]' does not match the requested catalog component 'mods'.",
      symbol = "g1[10]", component = "mods"
    )
  )
  statement <- function(condition) {
    c("RoBMA_hypothesis_statement",
      setdiff(class(condition), c("error", "condition")))
  }
  expected <- list(
    unclassed          = "RoBMA_hypothesis_statement",
    classed            = "BayesTools_stale_metadata",
    not_found          = statement(refusals[["not_found"]]),
    no_parameters      = statement(refusals[["no_parameters"]]),
    component_mismatch = statement(refusals[["component_mismatch"]])
  )
  fields <- function(condition) {
    unclass(condition)[setdiff(names(condition), c("message", "call"))]
  }
  # The plans are built with the mocked function, so none of them is kept.
  withr::defer(.hypothesis_plan_cache_clear())
  for (name in names(refusals)) {
    testthat::local_mocked_bindings(
      hypothesis_linear_target = function(...) stop(refusals[[name]]),
      .package = "BayesTools"
    )
    .hypothesis_plan_cache_clear()
    if (name == "classed") {
      for (read in list(
          function() .hypothesis_plans(fit, "g1[10] = g1[20]"),
          function() hypothesis(fit, "g1[10] = g1[20]", density_method = "KDE"),
          function() hypothesis(means, "g1[10] = g1[20]"))) {
        condition <- expect_error(read(), "Draw metadata are stale.", fixed = TRUE,
          class = "BayesTools_stale_metadata")
        expect_identical(class(condition), class(refusals[[name]]))
        expect_null(conditionCall(condition))
      }
      next
    }
    plan <- .hypothesis_plans(fit, "g1[10] = g1[20]")[[1L]]
    expect_identical(plan[["route"]], "combination", info = name)
    expect_identical(plan[["refusal"]][["class"]], expected[[name]], info = name)
    expect_identical(
      plan[["refusal"]][["reason"]], conditionMessage(refusals[[name]]),
      info = name
    )
    # Both combination routes refuse with the same classes, message, and
    # fields of the BayesTools condition.
    errors <- list(
      fit   = tryCatch(
        hypothesis(fit, "g1[10] = g1[20]", density_method = "KDE"),
        error = identity
      ),
      means = tryCatch(hypothesis(means, "g1[10] = g1[20]"), error = identity)
    )
    for (route in names(errors)) {
      info  <- paste(name, route)
      error <- errors[[route]]
      expect_identical(class(error), c(expected[[name]], "error", "condition"),
                       info = info)
      expect_identical(conditionMessage(error), conditionMessage(refusals[[name]]),
                       info = info)
      expect_identical(fields(error), fields(refusals[[name]]), info = info)
    }
  }

  # A refit request of BayesTools is no refusal of the statement: it stops
  # planning, unconverted, at both entry points.
  refit <- structure(
    class = c("BayesTools_refit_required", "error", "condition"),
    list(message = "Fitted metadata are unsupported. Refit the model.",
         call = NULL)
  )
  testthat::local_mocked_bindings(
    hypothesis_linear_target = function(...) stop(refit),
    .package = "BayesTools"
  )
  .hypothesis_plan_cache_clear()
  for (error in list(
    tryCatch(.hypothesis_plans(fit, "g1[10] = g1[20]"), error = identity),
    tryCatch(hypothesis(fit, "g1[10] = g1[20]", density_method = "KDE"),
             error = identity),
    tryCatch(hypothesis(means, "g1[10] = g1[20]"), error = identity)
  )) {
    expect_identical(error, refit)
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

  captured <- NULL
  compute <- BayesTools::hypothesis_BF
  testthat::local_mocked_bindings(
    hypothesis_BF = function(posterior, ...) {
      captured <<- c(list(posterior = posterior), list(...))
      do.call(compute, captured)
    },
    .package = "BayesTools"
  )
  # The other levels keep their qCMDE ordinates, and contrasts with the
  # reference level are defined.
  expect_s3_class(
    suppressWarnings(hypothesis(fit, c("g1[10] = 0", "g1[20] = 0"), seed = 1)),
    "data.frame"
  )
  for (level in c("10", "20")) {
    .expect_iwmde_child_identity(captured[["posterior"]][[level]], "mu_g1", level, 0)
    wrong <- captured
    ordinate <- BayesTools::posterior_metadata(wrong[["posterior"]][[level]], "posterior_ordinate")
    ordinate[["parameter"]] <- "mu_g1[5]"
    BayesTools::posterior_metadata(wrong[["posterior"]][[level]], "posterior_ordinate") <- ordinate
    expect_error(do.call(compute, wrong),
      "requires valid posterior ordinate or posterior density metadata", fixed = TRUE)
  }
  expect_s3_class(
    suppressWarnings(hypothesis(fit, "g1[10] = g1[5]", density_method = "KDE")),
    "data.frame"
  )
})


test_that("fixed reference-level statements keep one reason under qCMDE for model-averaged fits", {

  skip_on_cran()
  fit <- load_fit("lane_factor_averaged_reference")
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
