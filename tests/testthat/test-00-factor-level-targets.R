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


# The coefficient transformation target that hypothesis() pairs with a level,
# resolved through the same steps as hypothesis.brma().
.factor_level_prior_target <- function(fit, hypothesis) {

  metadata <- .brma_parameter_catalog_metadata(fit)
  parsed   <- BayesTools::hypothesis_parse(
    hypothesis,
    catalog        = metadata[["catalog"]],
    simplify_names = TRUE
  )
  selected <- .hypothesis_brma_select_parameter(
    object     = fit,
    hypothesis = parsed,
    component  = "auto",
    metadata   = metadata
  )
  rewritten <- .hypothesis_brma_rewrite(
    hypothesis = parsed,
    aliases    = selected[["aliases"]],
    parameter  = selected[["parameter"]]
  )
  point_refs <- .hypothesis_brma_point_refs(
    hypothesis     = rewritten,
    parameter      = selected[["parameter"]],
    require_direct = TRUE
  )
  targets <- .hypothesis_brma_formula_coefficient_level_targets(
    object     = fit,
    selected   = selected,
    point_refs = point_refs
  )

  return(vapply(targets, `[[`, character(1), "target"))
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
    # density over the unbounded Gaussian KDE of its draws at the null. The
    # comparison needs a null inside the KDE grid (levels 20 and 4 lie far
    # from zero; their ordinates are extrapolated).
    prior     <- fit[["priors"]][["mods"]][[cases[["term"]][i]]]
    posterior <- stats::density(mcmc[, cases[["coordinate"]][i]])
    height    <- stats::approx(posterior[["x"]], posterior[["y"]], xout = 0)[["y"]]
    if (!is.na(height)) {
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
      .factor_level_prior_target(fit, cases[["hypothesis"]][i]),
      stats::setNames(cases[["coordinate"]][i], cases[["level"]][i])
    )
  }
  expect_identical(compared, c("g1[10] = 0", "g2[2] = 0", "g2[3] = 0"))
})


test_that("mean-difference level point hypotheses stop without a single coordinate", {

  skip_on_cran()
  fit <- .factor_level_target_fits()[["meandif"]]

  # The first of four mean-difference levels has the design row (0, 1, 0) in
  # floating point: its extraction key is one coordinate with unit weight, the
  # coordinate of the contrast coefficient g2{2}. It is still not a level cell.
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
    c(hypothesis = "g1[10] = 0", level = "g1[10]", other = "g1[5]"),
    c(hypothesis = "g2[2] = 0",  level = "g2[2]",  other = "g2[1]"),
    c(hypothesis = "g2[1] = 0",  level = "g2[1]",  other = "g2[2]")
  )
  for (case in cases) {
    expect_error(
      suppressWarnings(hypothesis(fit, case[["hypothesis"]], density_method = "KDE")),
      paste0(
        "Point hypotheses on factor level '", case[["level"]], "' are not ",
        "supported: the level is a linear combination of the fitted contrast ",
        "coefficients (mean-difference, orthonormal, or ordered contrasts), ",
        "not a fitted coefficient itself, and point hypotheses on a single ",
        "level require a level fitted as its own coefficient (treatment or ",
        "independent contrasts). Test the level with a region hypothesis such ",
        "as '", case[["level"]], " > 0' or a level contrast such as '",
        case[["level"]], " = ", case[["other"]], "'."
      ),
      fixed = TRUE,
      info = case[["hypothesis"]]
    )
  }
  # The alternatives the message names are available.
  region   <- suppressWarnings(hypothesis(fit, "g1[10] > 0", density_method = "KDE"))
  contrast <- suppressWarnings(hypothesis(fit, "g1[10] = g1[5]", density_method = "KDE"))
  expect_true(is.finite(attr(region, "raw_BF")))
  expect_true(is.finite(attr(contrast, "raw_BF")))
})


test_that("contrast-coefficient selectors stop naming the level-label form", {

  skip_on_cran()
  fits <- .factor_level_target_fits()

  # Mean-difference coefficients are catalog quantities without RoBMA entries.
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
  }
  expect_error(
    plot(fits[["meandif"]], parameter = "g1{1}", plot_type = "ggplot"),
    paste0(
      "Factor contrast coefficients such as 'g1{1}' cannot be selected. ",
      "Select factor levels by their labels, such as 'g1[5]', or the whole ",
      "term 'g1'."
    ),
    fixed = TRUE
  )

  # Treatment levels are the coefficients: the selector has no quantity and
  # the example level skips the structural reference level 5.
  expect_error(
    suppressWarnings(hypothesis(fits[["treatment"]], "g1{1} = 0", density_method = "KDE")),
    paste0(
      "Hypotheses on factor contrast coefficients such as 'g1{1}' are not ",
      "supported. State them on factor levels by their labels, such as ",
      "'g1[10]'."
    ),
    fixed = TRUE
  )

  # Random-slope quantities such as 'tau(g1{3})' are no fixed contrast
  # coefficients, and '{j}' after a name that is no fixed factor term is a
  # parse error: neither gets the coefficient message, which names the fixed
  # coefficient selector wherever it appears.
  metadata <- .brma_parameter_catalog_metadata(fits[["treatment"]])
  expect_null(.hypothesis_brma_check_coefficient_selector("tau(g1{3}) > 0.05", metadata))
  expect_null(.hypothesis_brma_check_coefficient_selector("(mu) tau(g1{3}) > 0", metadata))
  expect_null(.hypothesis_brma_check_coefficient_selector("h{1} > 0", metadata))
  expect_error(
    .hypothesis_brma_check_coefficient_selector("tau(g2{1}) > 0 & g1{1} > 0", metadata),
    "Hypotheses on factor contrast coefficients such as 'g1{1}' are not supported.",
    fixed = TRUE
  )
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


test_that("hypothesis_quantities reports point tests only for fitted level coefficients", {

  skip_on_cran()
  fits <- .factor_level_target_fits()

  # Treatment levels are fitted coefficients (the reference level is fixed).
  treatment <- hypothesis_quantities(fits[["treatment"]])
  treatment <- treatment[treatment[["term"]] %in% c("g1", "g2"), , drop = FALSE]
  expect_true(all(treatment[["point_test"]]))
  expect_identical(unique(treatment[["point_test_methods"]]), "KDE, qCMDE, IWMDE")
  expect_false(any(nzchar(treatment[["reason"]])))

  # Mean-difference levels are linear combinations of the contrast
  # coefficients; hypothesis() stops their point hypotheses.
  meandif <- hypothesis_quantities(fits[["meandif"]])
  g1 <- meandif[meandif[["alias"]] == "g1", , drop = FALSE]
  expect_false(g1[["point_test"]])
  expect_true(g1[["direction_test"]])
  expect_identical(g1[["point_test_methods"]], "")
  expect_identical(g1[["reason"]], paste0(
    "Point hypotheses are not supported for levels 'g1[5]', 'g1[10]', ",
    "'g1[20]': each is a linear combination of the fitted contrast ",
    "coefficients (mean-difference, orthonormal, or ordered contrasts), not ",
    "a fitted coefficient itself. Region hypotheses and level contrasts are ",
    "available for all levels."
  ))
  expect_false(any(meandif[["point_test"]][meandif[["term"]] == "g2"]))
  for (level in c("5", "10", "20")) {
    expect_error(
      suppressWarnings(hypothesis(
        fits[["meandif"]], paste0("g1[", level, "] = 0"), density_method = "KDE"
      )),
      "Point hypotheses on factor level",
      fixed = TRUE,
      info = level
    )
  }

  # Ordered coding: the first increment is a fitted coefficient, but its
  # prior (the ordered total times its allocation) has no exact ordinate; the
  # later level is a sum of increments. No level supports point hypotheses.
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
  expect_false(g[["point_test"]])
  expect_identical(g[["point_test_methods"]], "")
  expect_match(
    g[["reason"]],
    "Point hypotheses are not supported for level 'g[hi]': it is a linear",
    fixed = TRUE
  )
  expect_match(
    g[["reason"]],
    paste0(
      "Point hypotheses are not supported for level 'g[mid]': the induced ",
      "prior of its fitted coefficient has no exact ordinate"
    ),
    fixed = TRUE
  )
  expect_error(
    suppressWarnings(hypothesis(ordered, "g[hi] = 0", density_method = "KDE")),
    "Point hypotheses on factor level 'g[hi]'",
    fixed = TRUE
  )
  # The stop names the level selector, not the backend coordinate 'mu_g[1]'.
  for (method in c("KDE", "qCMDE")) {
    expect_error(
      suppressWarnings(hypothesis(ordered, "g[mid] = 0.1", density_method = method)),
      paste0(
        "The induced prior ordinate for factor level 'g[mid]' is not exact ",
        "enough for a point-null Bayes factor."
      ),
      fixed = TRUE,
      info = method
    )
  }

  # With numeric labels, the first increment's coordinate 'mu_g[1]' reads as
  # the reference level '1': the stop names the level the hypothesis wrote.
  data[["g"]] <- factor(
    c("1", "2", "3")[as.integer(data[["g"]])],
    levels = c("1", "2", "3")
  )
  numeric_labels <- suppressWarnings(brma(
    yi = yi, sei = sei, mods = ~ g, data = data, measure = "SMD",
    prior_mods = list(g = BayesTools::prior_ordered(BayesTools::prior("normal", list(0, 1)))),
    chains = 1, sample = 500, burnin = 100, adapt = 100, seed = 1, silent = TRUE
  ))
  message <- tryCatch(
    suppressWarnings(hypothesis(numeric_labels, "g[2] = 0.1", density_method = "KDE")),
    error = conditionMessage
  )
  expect_identical(message, paste0(
    "The induced prior ordinate for factor level 'g[2]' is not exact enough ",
    "for a point-null Bayes factor."
  ))
})


test_that("formula coefficient messages name factor levels by their selector", {

  # A level target holds its backend coordinate ('mu_g[1]' for level '2');
  # messages name the level selector instead. Scalar targets keep their
  # parameter name.
  transform <- structure(
    list(
      schema_version    = 1L,
      target_scale      = "original",
      target_names      = "mu_g[1]",
      matrix            = matrix(
        1, 1L, 1L, dimnames = list("mu_g[1]", "mu_g[1]")
      ),
      source_transforms = c("mu_g[1]" = "log"),
      output_transforms = list("mu_g[1]" = "identity")
    ),
    class = "BayesTools_formula_coefficient_transform"
  )
  level  <- list(
    formula_parameter = "mu",
    target            = "mu_g[1]",
    target_i          = 1L,
    transform         = transform,
    level_selector    = "g[2]"
  )
  scalar <- level[setdiff(names(level), "level_selector")]

  expect_identical(
    .hypothesis_brma_formula_transform_route(level)[["reason"]],
    "The fitted nonlinear joint coefficient transform for 'g[2]' is not supported by hypothesis()."
  )
  expect_identical(
    .hypothesis_brma_formula_transform_route(scalar)[["reason"]],
    "The fitted nonlinear joint coefficient transform for 'mu_g[1]' is not supported by hypothesis()."
  )
  fixed <- level
  fixed[["transform"]][["matrix"]][1L, 1L] <- 0
  expect_identical(
    .hypothesis_brma_formula_transform_route(fixed)[["reason"]],
    "The fitted coefficient 'g[2]' is structurally fixed and has no posterior hypothesis route."
  )
  uncertified <- level
  uncertified[["transform"]][["target_scale"]] <- "fitted"
  expect_identical(
    .hypothesis_brma_formula_transform_route(uncertified)[["reason"]],
    paste0(
      "The fitted coefficient transform for 'g[2]' lacks the certified ",
      "structural metadata required for hypothesis testing."
    )
  )

  supported <- level
  supported[["route"]] <- list(type = "exp_affine", support = c(0, Inf))
  expect_error(
    .hypothesis_brma_check_formula_point_support(
      data.frame(value = 0), supported
    ),
    "open support for factor level 'g[2]'.",
    fixed = TRUE
  )
  scalar_supported <- supported[setdiff(names(supported), "level_selector")]
  expect_error(
    .hypothesis_brma_check_formula_point_support(
      data.frame(value = 0), scalar_supported
    ),
    "open support for transformed coefficient 'mu_g[1]'.",
    fixed = TRUE
  )

  # The prior-ordinate stop and the qCMDE/IWMDE transform reason.
  exact <- TRUE
  testthat::local_mocked_bindings(
    JAGS_formula_prior_density = function(...) list(),
    prior_density_ordinate     = function(...) list(exact = exact),
    .package = "BayesTools"
  )
  prior_target <- .hypothesis_brma_formula_prior_target(
    object       = list(fit = NULL),
    samples      = list(),
    hypothesis   = NULL,
    target_info  = supported,
    point_values = 0.1,
    force_linear = TRUE
  )
  expect_identical(
    prior_target[["parameter_spec"]][["reason"]],
    paste0(
      "qCMDE/IWMDE does not support the fitted nonlinear joint transform ",
      "for 'g[2]'. Use density_method = 'KDE' or ",
      "standardized_coefficients = TRUE."
    )
  )
  exact <- FALSE
  expect_error(
    .hypothesis_brma_formula_prior_target(
      object       = list(fit = NULL),
      samples      = list(),
      hypothesis   = NULL,
      target_info  = supported,
      point_values = 0.1,
      force_linear = TRUE
    ),
    paste0(
      "The induced prior ordinate for factor level 'g[2]' is not exact ",
      "enough for a point-null Bayes factor."
    ),
    fixed = TRUE
  )
})


test_that("an ambiguous level resolution names the level selector", {

  testthat::local_mocked_bindings(
    JAGS_formula_coefficient_transform = function(...) structure(
      list(schema_version = 1L),
      class = "BayesTools_formula_coefficient_transform"
    ),
    parameter_catalog = function(...) list(quantities = data.frame(
      quantity_id = character(), stringsAsFactors = FALSE
    )),
    .package = "BayesTools"
  )
  selected <- list(
    parameter  = "mu_g",
    component  = "mods",
    aliases    = list(g = "mu_g", mu_g = "mu_g"),
    entry      = list(
      role              = "formula_coefficient_group",
      formula_parameter = "mu"
    ),
    resolution = list(occurrences = data.frame(
      level          = c("2", "2"),
      canonical_name = c("mu_g[2]", "mu_h[2]"),
      quantity_id    = c("q1", "q2"),
      stringsAsFactors = FALSE
    ))
  )
  expect_error(
    .hypothesis_brma_formula_coefficient_level_targets(
      object     = list(fit = NULL),
      selected   = selected,
      point_refs = data.frame(level = "2", value = 0.1, stringsAsFactors = FALSE)
    ),
    "Factor level 'g[2]' is ambiguous in the fitted parameter catalog.",
    fixed = TRUE
  )
})


test_that("point hypotheses on the fixed reference level give the KDE reason for every method", {

  skip_on_cran()
  fit <- .factor_level_target_fits()[["treatment"]]

  # The treatment reference level g1[5] is fixed at 0 by the contrast: its
  # declared atom decides point statements on it, whatever the density
  # method, and no qCMDE/IWMDE ordinate enters.
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
    expect_type(kde, "character")
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
  expect_match(
    tryCatch(
      suppressWarnings(hypothesis(fit, "g1[5] = 0", density_method = "KDE")),
      error = conditionMessage
    ),
    "declared point mass at the exact null hypothesis value",
    fixed = TRUE
  )

  # The other levels keep their qCMDE ordinates.
  expect_s3_class(
    suppressWarnings(hypothesis(fit, c("g1[10] = 0", "g1[20] = 0"), seed = 1)),
    "data.frame"
  )
})


test_that("fixed reference-level statements keep the model-averaged KDE reason under qCMDE", {

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

  # For model-averaged fits the KDE reason names the inclusion Bayes factor.
  # A call mixing the fixed reference level with another level evaluates the
  # fixed-level statement before the qCMDE/IWMDE ordinates and must keep it.
  message_of <- function(expr) tryCatch(suppressWarnings(expr), error = conditionMessage)
  for (statement in list("g1[5] = 0", c("g1[5] = 0", "g1[10] = 0.1"),
                         c("g1[10] = 0.1", "g1[5] = 0"))) {
    info <- paste(statement, collapse = ", ")
    kde  <- message_of(hypothesis(fit, statement, density_method = "KDE"))
    expect_match(kde, "inclusion Bayes factor", fixed = TRUE, info = info)
    for (method in c("qCMDE", "IWMDE")) {
      expect_identical(
        message_of(hypothesis(fit, statement, density_method = method)),
        kde,
        info = paste(info, method)
      )
    }
  }
})
