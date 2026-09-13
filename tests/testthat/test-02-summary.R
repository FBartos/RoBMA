context("Summary")

source(testthat::test_path("common-functions.R"))

test_that("summary.brma does not print known-V backend metadata", {

  section_names <- c(
    "inclusion_components",
    "inclusion_mods",
    "inclusion_scale",
    "inclusion_random",
    "estimates",
    "estimates_conditional",
    "estimates_mods",
    "estimates_mods_conditional",
    "estimates_scale",
    "estimates_scale_conditional",
    "estimates_random",
    "estimates_random_conditional",
    "estimates_bias",
    "estimates_bias_conditional"
  )
  out <- as.list(stats::setNames(vector("list", length(section_names)), section_names))
  out[["name"]] <- "Bayesian Multivariate Meta-Analysis"
  out[["known_v_backend"]] <- list(
    known_v                             = TRUE,
    known_v_parameterization           = "block_mvn",
    known_v_parameterization_requested = "auto"
  )
  class(out) <- "summary.brma"

  output <- capture.output(print(out))

  expect_false(any(grepl("Known-V backend", output, fixed = TRUE)))
})

test_that("summary.brma coerces displayed sections to one data frame", {

  inclusion <- data.frame(
    prior_prob   = 0.5,
    post_prob    = 0.75,
    inclusion_BF = 3,
    row.names    = "Effect",
    check.names  = FALSE
  )
  location <- data.frame(
    Mean        = c(0.1, 0.2),
    SD          = c(0.01, 0.02),
    `0.025`     = c(0.08, 0.16),
    `0.975`     = c(0.12, 0.24),
    row.names   = c("sensitivity", "specificity"),
    check.names = FALSE
  )
  random <- data.frame(
    `Random name`      = "study",
    `Random grouping`  = "study",
    `Random structure` = "diag",
    Mean               = 0.3,
    SD                 = 0.03,
    `0.025`            = 0.24,
    `0.975`            = 0.36,
    row.names          = "tau(sensitivity)",
    check.names        = FALSE
  )
  sections <- c(
    "inclusion_components", "inclusion_mods", "inclusion_scale",
    "inclusion_random",
    "estimates", "estimates_conditional",
    "estimates_mods", "estimates_mods_conditional",
    "estimates_scale", "estimates_scale_conditional",
    "estimates_random", "estimates_random_conditional",
    "estimates_bias", "estimates_bias_conditional"
  )
  out <- c(
    list(
      name            = "Bayesian Multivariate Meta-Analysis",
      known_v_backend = list()
    ),
    as.list(stats::setNames(vector("list", length(sections)), sections))
  )
  out[["inclusion_components"]] <- inclusion
  out[["estimates_mods"]]       <- location
  out[["estimates_random"]]     <- random
  class(out) <- "summary.brma"

  frame <- as.data.frame(out)

  expect_identical(names(frame)[1:2], c("component", "parameter"))
  expect_identical(
    frame[["component"]],
    c("inclusion", "location", "location", "random")
  )
  expect_identical(
    frame[["parameter"]],
    c("Effect", "sensitivity", "specificity", "tau(sensitivity)")
  )
  expect_true(all(c("CI_0.025", "CI_0.975") %in% names(frame)))
  expect_false(any(c(
    "Random name", "Random grouping", "Random structure"
  ) %in% names(frame)))
  expect_true(is.na(frame[["Mean"]][1L]))
  expect_true(all(is.na(frame[["prior_prob"]][-1L])))
  expect_identical(data.frame(out), data.frame(frame))
})


test_that("random inclusion labels identify aggregate and component SDs", {

  labels <- c("(mu) inclusion", "(mu) inclusion(component)",
              "(mu) inclusion(study:esid)", "(mu) sd_total",
              "(mu) component: sd_total", "(mu) study:esid: sd",
              "(mu) study: sd")
  quantities <- data.frame(
    canonical_name = labels, display_label = labels,
    role = c(rep("random_inclusion", 3), rep("random_sd", 4)),
    quantity = c(rep("inclusion", 3), "sd_total", "sd_total", "sd", "sd"),
    owner_name = c("", "", "", "", "component", "study:esid", "study")
  )
  quantities[["extraction_key"]] <- I(list(
    list(source_parameter = "shared_gate"),
    list(source_parameter = "component_gate"),
    list(source_parameter = "study_gate"),
    list(dependencies = "shared_gate"),
    list(dependencies = "component_gate"),
    list(random_block = "study_esid", dependencies = c("component_gate", "study_gate")),
    list(random_block = "study", dependencies = "shared_gate")
  ))
  testthat::local_mocked_bindings(
    parameter_catalog = function(...) list(quantities = quantities),
    .package = "BayesTools"
  )
  testthat::local_mocked_bindings(
    .random_component_inclusion_map = function(...) {
      list(study_esid = "study_gate", study = "shared_gate")
    },
    .package = "RoBMA"
  )

  expect_identical(
    .summary_random_inclusion_labels(list(), labels[1:3], labels[1:3]),
    c("tau_total", "component: tau_total", "study:esid: tau")
  )
})

test_that("scalar multivariate summary labels respect declared SD priors", {

  args <- list(
    yi                        = c(0.1, 0.2, 0.3, 0.4),
    vi                        = rep(0.01, 4),
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE,
    silent                    = TRUE
  )
  priors <- list(
    NULL,
    prior("point", list(location = 0.3)),
    prior("normal", list(mean = 0, sd = 1), truncation = list(lower = 0))
  )
  model_types <- c("Fixed-Effect", "Random-Effects", "Random-Effects")
  footnote <- paste0(
    "exp(intercept) is the baseline heterogeneity SD (%s), already ",
    "exponentiated. Other coefficients are changes in log(SD); ",
    "exp(coefficient) is an SD multiplier."
  )

  for (i in seq_along(priors)) {
    arguments <- args
    if (!is.null(priors[[i]])) {
      arguments[["random"]] <- ~ 1 | estimate
      arguments[["data"]] <- data.frame(estimate = seq_along(args[["yi"]]))
      arguments[["selection"]] <- BayesTools::selection_model(group = "estimate")
      arguments[["prior_heterogeneity"]] <- priors[[i]]
    }
    object <- do.call(bselmodel.mv, arguments)
    expect_identical(
      .summary.brma_model_names(object),
      paste("Bayesian Multivariate", model_types[i], "Selection Model (k = 4)")
    )
    expect_identical(.summary_scale_footnotes(object), if (is.null(priors[[i]])) {
      sprintf(footnote, "tau")
    } else {
      paste0("exp(intercept) is the baseline heterogeneity SD (tau) of the ",
        "indicated target, already ",
        "exponentiated. Other coefficients are changes in log(SD); ",
        "exp(coefficient) is an SD multiplier.")
    })
  }

  for (cluster in list(NULL, c("a", "a", "b", "b"))) {
    object <- do.call(
      brma, c(args, list(prior_heterogeneity = NULL, cluster = cluster))
    )
    expect_identical(.summary_scale_footnotes(object), sprintf(footnote, "tau"))
    expect_identical(
      .summary.brma_model_names(object),
      if (is.null(cluster)) {
        "Bayesian Fixed-Effect Model (k = 4)"
      } else {
        "Bayesian Multilevel Fixed-Effect Model (k = 4, clusters = 2)"
      }
    )
  }
})

skip_if_no_fits()
fit_names <- list_fits()
fits      <- lazy_fits(fit_names, validate = FALSE)

summary_sections <- c(
  "name",
  "known_v_backend",
  "inclusion_components",
  "inclusion_mods",
  "inclusion_scale",
  "inclusion_random",
  "estimates",
  "estimates_conditional",
  "estimates_mods",
  "estimates_mods_conditional",
  "estimates_scale",
  "estimates_scale_conditional",
  "estimates_random",
  "estimates_random_conditional",
  "estimates_bias",
  "estimates_bias_conditional"
)

summary_table_sections <- setdiff(summary_sections, c("name", "known_v_backend"))

summary_common_parameters <- function(fit) {

  prior_names <- names(attr(fit[["fit"]], "prior_list"))
  return(intersect(c("mu", "tau", "rho"), prior_names))
}

expect_summary_subtable <- function(table, name, section) {

  if (length(table) == 0L) {
    return(invisible(NULL))
  }

  info <- paste0("summary section '", section, "' for '", name, "'")
  expect_true(is.matrix(table) || is.data.frame(table), info = info)
  expect_true(nrow(table) > 0L, info = info)
  expect_true(ncol(table) > 0L, info = info)
  expect_false(is.null(rownames(table)), info = info)
  expect_true(all(nzchar(rownames(table))), info = info)

  values <- unlist(
    as.data.frame(table)[vapply(as.data.frame(table), is.numeric, TRUE)],
    use.names = FALSE
  )
  if (length(values) > 0L) {
    expect_true(all(is.na(values) | !is.nan(values)), info = info)
  }

  title <- attr(table, "title")
  expect_true(
    is.character(title) && length(title) == 1L && nzchar(title),
    info = paste0(info, " title")
  )
}

expect_summary_sections <- function(summary_object, fit, name, conditional = FALSE) {

  common_parameters <- summary_common_parameters(fit)

  expect_equal(
    length(summary_object[["estimates"]]) > 0L,
    length(common_parameters) > 0L,
    info = paste0("common estimate section for '", name, "'")
  )
  expect_equal(
    length(summary_object[["estimates_mods"]]) > 0L,
    .is_mods(fit) || .is_random(fit),
    info = paste0("moderator estimate section for '", name, "'")
  )
  expect_equal(
    length(summary_object[["estimates_scale"]]) > 0L,
    .is_scale(fit),
    info = paste0("scale estimate section for '", name, "'")
  )
  expect_equal(
    length(summary_object[["estimates_bias"]]) > 0L,
    .is_bias(fit),
    info = paste0("bias estimate section for '", name, "'")
  )
  expect_equal(
    length(summary_object[["estimates_random"]]) > 0L,
    .summary_random_components_enabled(fit),
    info = paste0("random estimate section for '", name, "'")
  )
  expect_equal(
    length(summary_object[["inclusion_components"]]) > 0L,
    .is_RoBMA(fit),
    info = paste0("component inclusion section for '", name, "'")
  )

  conditional_sections <- c(
    "estimates_conditional",
    "estimates_mods_conditional",
    "estimates_scale_conditional",
    "estimates_random_conditional",
    "estimates_bias_conditional"
  )
  conditional_present <- vapply(
    summary_object[conditional_sections],
    function(x) length(x) > 0L,
    TRUE
  )

  if (conditional) {
    expect_true(
      any(conditional_present),
      info = paste0("conditional sections for '", name, "'")
    )
  } else {
    expect_false(
      any(conditional_present),
      info = paste0("conditional sections for '", name, "'")
    )
  }
}

expect_summary_contract <- function(summary_object, fit, name,
                                    conditional = FALSE) {

  expect_s3_class(summary_object, "summary.brma")
  expected_sections <- summary_sections
  sensitivity       <- fit[["selection_sensitivity_diagnostics"]]
  if (!is.null(sensitivity)) {
    expected_sections <- c(expected_sections, "selection_sensitivity_diagnostics")
    expect_identical(summary_object[["selection_sensitivity_diagnostics"]], sensitivity)
  }
  selection_model <- .data_selection_model(fit[["data"]])
  if (!is.null(selection_model)) {
    expected_sections <- c(expected_sections, "selection_model", "selection_sampling")
    expect_identical(summary_object[["selection_model"]], selection_model)
    expect_identical(
      summary_object[["selection_sampling"]],
      .selection_postfit_target_metadata(fit[["data"]])[["sampling_structure"]]
    )
  }
  expect_named(summary_object, expected_sections)
  expect_type(summary_object[["name"]], "character")
  expect_true(length(summary_object[["name"]]) == 1L)
  expect_true(nzchar(summary_object[["name"]]))
  expect_type(summary_object[["known_v_backend"]], "list")
  expect_identical(
    isTRUE(summary_object[["known_v_backend"]][["known_v"]]),
    inherits(fit, "brma.mv") && .is_data_known_v(fit[["data"]])
  )
  if (isTRUE(summary_object[["known_v_backend"]][["known_v"]])) {
    known_V <- .data_known_v_data(fit[["data"]])
    expect_equal(
      summary_object[["known_v_backend"]][["known_v_parameterization"]],
      known_V[["parameterization"]]
    )
    expect_equal(
      summary_object[["known_v_backend"]][["known_v_parameterization_requested"]],
      known_V[["parameterization_requested"]]
    )
  }

  expect_identical(attr(summary_object, "mods"), .is_mods(fit))
  expect_identical(attr(summary_object, "scale"), .is_scale(fit))
  expect_identical(attr(summary_object, "random"), .is_random(fit))
  expect_identical(attr(summary_object, "multilevel"), .is_multilevel(fit))
  expect_identical(attr(summary_object, "bias"), .is_bias(fit))
  expect_identical(attr(summary_object, "RoBMA"), .is_RoBMA(fit))
  expect_identical(attr(summary_object, "outcome_type"), .outcome_type(fit))

  for (section in summary_table_sections) {
    expect_summary_subtable(summary_object[[section]], name, section)
  }
  expect_summary_sections(summary_object, fit, name, conditional = conditional)
}

expect_printed_summary <- function(summary_object, name) {

  output <- capture.output(print(summary_object))
  expect_true(any(nzchar(output)), info = paste0("printed summary for '", name, "'"))
  expect_true(any(grepl("Bayesian", output, fixed = TRUE)),
              info = paste0("printed summary model name for '", name, "'"))
  expect_false(any(grepl("Known-V backend", output, fixed = TRUE)),
               info = paste0("printed known-V backend for '", name, "'"))
  expect_false(any(grepl("__xXx__", output, fixed = TRUE)),
               info = paste0("printed summary labels for '", name, "'"))
}


test_that("summary.brma returns a stable object contract", {

  for (name in names(fits)) {
    fit <- fits[[name]]
    out <- summary(fit)

    expect_summary_contract(out, fit, name)
    expect_printed_summary(out, name)
  }
})

test_that("scalar multivariate summaries retain RoBMA tau naming", {

  scalar_names <- c("brma.mv_block_mvn", "brma.mv_block_mvn_fixed_random_null")
  skip_if_missing_fits(scalar_names)

  for (name in scalar_names) {
    fit <- fits[[name]]
    map <- BayesTools::parameter_map(fit[["fit"]])
    reference <- .summary_estimates_table(
      object                   = fit,
      probs                    = c(.025, .50, .975),
      include_mcmc_diagnostics = FALSE,
      is_robma                 = FALSE,
      transform_factors        = TRUE,
      transform_scaled         = TRUE,
      keep_parameters          = c("mu", "tau"),
      random_effects_summary   = "none",
      title                    = "Estimates"
    )
    out <- summary(fit, include_mcmc_diagnostics = FALSE)

    expect_identical(rownames(out[["estimates"]]), c("mu", "tau"))
    expect_identical(attr(out[["estimates"]], "parameters"), c("mu", "tau"))
    restored <- out[["estimates"]]
    rownames(restored) <- rownames(reference)
    expect_identical(restored, reference)
    expect_identical(BayesTools::parameter_map(fit[["fit"]]), map)
    expect_identical(
      out[["name"]],
      paste(
        "Bayesian Multivariate",
        "Fixed-Effect",
        sprintf("Model (k = %i)", nrow(fit[["data"]][["outcome"]]))
      )
    )

    frame <- as.data.frame(out)
    expect_identical(frame[["parameter"]][frame[["component"]] == "common"],
                     c("mu", "tau"))
    expect_true(all(c("CI_0.025", "CI_0.975") %in% names(frame)))
    expect_identical(data.frame(out), frame)
  }
})

test_that("summary.brma prints known-R random-effect parameters", {

  name <- "brma.mv_block_mvn_known_R"
  skip_if_missing_fits(name)

  out    <- summary(fits[[name]], include_mcmc_diagnostics = FALSE)
  output <- capture.output(print(out))

  expect_summary_contract(out, fits[[name]], name)
  expect_true("tau" %in% rownames(out[["estimates_random"]]))
  expect_true(any(grepl("tau", output, fixed = TRUE)))
  expect_false(any(grepl("tau_mult", output, fixed = TRUE)))
  expect_false(any(grepl("group_covariance", output, fixed = TRUE)))
})


test_that("summary.brma omits a fixed-zero intercept with moderators", {

  name <- "brma.mv_v14_ishak2007_har"
  skip_if_missing_fits(name)

  out <- summary(fits[[name]], include_mcmc_diagnostics = FALSE)

  expect_false(any(grepl(
    "intercept",
    tolower(rownames(out[["estimates_mods"]])),
    fixed = TRUE
  )))
  expect_true(any(grepl("time_factor", rownames(out[["estimates_mods"]]),
                        fixed = TRUE)))
})

test_that("summary.brma options change table schema", {

  name <- "bcg_meta-analysis"
  skip_if_missing_fits(name)

  out <- summary(
    fits[[name]],
    probs                    = c(0.01, 0.99),
    include_mcmc_diagnostics = FALSE
  )
  expect_summary_contract(out, fits[[name]], name)

  cols <- colnames(out[["estimates"]])
  expect_true(all(c("Mean", "SD", "0.01", "0.99") %in% cols))
  expect_false(any(grepl("error\\(MCMC\\)|ESS|R-hat", cols)))
})

test_that("summary.brma validates standardized coefficient flag", {

  name <- "bcg_meta-analysis"
  skip_if_missing_fits(name)

  expect_error(
    summary(fits[[name]], standardized_coefficients = "no"),
    regexp = "standardized_coefficients"
  )
  expect_error(
    summary(fits[[name]], standardized_coefficients = c(TRUE, FALSE)),
    regexp = "standardized_coefficients"
  )
})

test_that("summary.brma controls BayesTools diagnostic columns", {

  name <- "dat.lehmann2018_RoBMA"
  skip_if_missing_fits(name)

  old_options <- options(
    BayesTools.JAGS_estimates_diagnostic_columns = c("ESS", "R_hat"),
    BayesTools.JAGS_BF_diagnostic_columns        = c("ESS", "MCMC_error", "BF_error_percent")
  )
  on.exit(options(old_options), add = TRUE)

  out <- summary(fits[[name]], include_mcmc_diagnostics = TRUE)
  estimate_diagnostics <- c("MCMC_error", "MCMC_SD_error", "ESS", "R_hat")
  BF_diagnostics       <- c("ESS", "MCMC_error", "BF_error_percent")

  expect_true(all(estimate_diagnostics %in% colnames(out[["estimates"]])))
  expect_true("BF_error_percent" %in% colnames(out[["inclusion_components"]]))
  expect_false(any(c("ESS", "MCMC_error") %in%
                     colnames(out[["inclusion_components"]])))

  out_none <- summary(fits[[name]], include_mcmc_diagnostics = FALSE)
  expect_false(any(estimate_diagnostics %in% colnames(out_none[["estimates"]])))
  expect_false(any(BF_diagnostics %in%
                     colnames(out_none[["inclusion_components"]])))
})

test_that("summary.brma inclusion summaries support BF direction and log scale", {

  name <- "dat.lehmann2018_RoBMA"
  skip_if_missing_fits(name)

  bf_column <- function(x) {

    return(as.data.frame(x[["inclusion_components"]])[["inclusion_BF"]])
  }

  out_default <- summary(fits[[name]], include_mcmc_diagnostics = FALSE)
  out_log     <- summary(fits[[name]], include_mcmc_diagnostics = FALSE, logBF = TRUE)
  out_BF01    <- summary(fits[[name]], include_mcmc_diagnostics = FALSE, BF01 = TRUE)
  out_both    <- summary(
    fits[[name]],
    include_mcmc_diagnostics = FALSE,
    logBF                    = TRUE,
    BF01                     = TRUE
  )

  expect_false(isTRUE(attr(bf_column(out_default), "logBF")))
  expect_false(isTRUE(attr(bf_column(out_default), "BF01")))
  expect_true(isTRUE(attr(bf_column(out_log), "logBF")))
  expect_false(isTRUE(attr(bf_column(out_log), "BF01")))
  expect_false(isTRUE(attr(bf_column(out_BF01), "logBF")))
  expect_true(isTRUE(attr(bf_column(out_BF01), "BF01")))
  expect_true(isTRUE(attr(bf_column(out_both), "logBF")))
  expect_true(isTRUE(attr(bf_column(out_both), "BF01")))

  default_effect <- as.numeric(bf_column(out_default)["Effect"])
  expect_equal(
    as.numeric(bf_column(out_log)["Effect"]),
    log(default_effect),
    tolerance = sqrt(.Machine$double.eps)
  )
  expect_equal(
    as.numeric(bf_column(out_BF01)["Effect"]),
    1 / default_effect,
    tolerance = sqrt(.Machine$double.eps)
  )
  expect_equal(
    as.numeric(bf_column(out_both)["Effect"]),
    log(1 / default_effect),
    tolerance = sqrt(.Machine$double.eps)
  )
})

test_that("summary.brma inclusion subtables preserve row-level BF bounds", {

  inclusion_BF <- BayesTools::format_BF(
    c(14999, 1.257, 0.588),
    inclusion = TRUE
  )
  attr(inclusion_BF, "bound_operator") <- c(">", NA_character_, NA_character_)
  class(inclusion_BF) <- unique(c("BayesTools_BF", class(inclusion_BF)))

  table <- data.frame(
    prior_prob   = c(0.5, 0.5, 0.5),
    post_prob    = c(1.0, 0.557, 0.370),
    inclusion_BF = inclusion_BF,
    row.names    = c("(mu) intercept", "(mu) tailor", "tau"),
    check.names  = FALSE
  )
  class(table)              <- c("BayesTools_table", class(table))
  attr(table, "type")       <- c("prior_prob", "post_prob", "inclusion_BF")
  attr(table, "parameters") <- c("mu_intercept", "mu_tailor", "tau")
  attr(table, "rownames")   <- TRUE

  mod_table <- .summary.inclusion_subtable(
    table      = table,
    indices    = 2L,
    row_labels = "tailor",
    title      = "Meta-Regression Inclusion"
  )

  expect_identical(
    attr(mod_table[["inclusion_BF"]], "bound_operator"),
    NA_character_
  )
  expect_false(any(grepl(">1.257", capture.output(print(mod_table)),
                         fixed = TRUE)))
})

test_that("summary.brma standardized coefficients use the standardized scale", {

  name <- "bangertdrowns2004_location-scale"
  skip_if_missing_fits(name)

  fit             <- fits[[name]]
  out_default     <- summary(fit)
  out_standardized <- summary(fit, standardized_coefficients = TRUE)

  expect_summary_contract(out_standardized, fit, name)
  expect_true("ni100" %in% rownames(out_default[["estimates_mods"]]))
  expect_true("ni100" %in% rownames(out_standardized[["estimates_mods"]]))
  expect_gt(
    abs(out_default[["estimates_mods"]]["ni100", "Mean"] -
          out_standardized[["estimates_mods"]]["ni100", "Mean"]),
    sqrt(.Machine$double.eps)
  )

  expect_equal(
    capture.output(print(fit)),
    capture.output(print(out_default)),
    info = "print.brma delegates to summary.brma"
  )
})

test_that("RoBMA conditional summaries expose conditional sections", {

  name <- "dat.lehmann2018_RoBMA_mods"
  skip_if_missing_fits(c("bcg_meta-analysis", name))

  out <- summary(fits[[name]], conditional = TRUE)
  expect_summary_contract(out, fits[[name]], name, conditional = TRUE)

  conditional_titles <- vapply(
    c(
      "estimates_conditional",
      "estimates_mods_conditional",
      "estimates_scale_conditional",
      "estimates_random_conditional",
      "estimates_bias_conditional"
    ),
    function(section) {
      title <- attr(out[[section]], "title")
      if (is.null(title)) "" else title
    },
    character(1)
  )
  expect_true(any(grepl("Conditional", conditional_titles, fixed = TRUE)))

  expect_error(
    summary(fits[["bcg_meta-analysis"]], conditional = TRUE),
    "RoBMA objects"
  )
})

test_that("RoBMA inclusion summaries use user-facing labels", {

  skip_if_missing_fits(c("dat.lehmann2018_RoBMA", "dat.lehmann2018_RoBMA_mods2"))

  out_simple <- summary(fits[["dat.lehmann2018_RoBMA"]])
  expect_true(all(c("Effect", "Heterogeneity", "Publication Bias") %in%
                    rownames(out_simple[["inclusion_components"]])))

  out_mods2 <- summary(fits[["dat.lehmann2018_RoBMA_mods2"]])
  expect_false(any(grepl("__xXx__", rownames(out_mods2[["inclusion_mods"]]),
                         fixed = TRUE)))
  expect_true("Preregistered:Gender" %in%
                rownames(out_mods2[["inclusion_mods"]]))
})
