context("Summary")

source(testthat::test_path("common-functions.R"))

test_that("summary.brma prints known-V backend metadata", {

  section_names <- c(
    "inclusion_components",
    "inclusion_mods",
    "inclusion_scale",
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

  expect_true(any(grepl(
    "Known-V backend: block_mvn (requested: auto)",
    output,
    fixed = TRUE
  )))
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
  expect_named(summary_object, summary_sections)
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
  known_v_label <- .brma_mv_known_v_backend_label(summary_object[["known_v_backend"]])
  expect_identical(
    any(grepl("Known-V backend", output, fixed = TRUE)),
    !is.null(known_v_label),
    info = paste0("printed known-V backend for '", name, "'")
  )
  if (!is.null(known_v_label)) {
    expect_true(
      any(grepl(known_v_label, output, fixed = TRUE)),
      info = paste0("printed known-V backend label for '", name, "'")
    )
  }
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

test_that("summary.brma prints known-R random-effect parameters", {

  name <- "brma.mv_block_mvn_known_R"
  skip_if_missing_fits(name)

  out    <- summary(fits[[name]], include_mcmc_diagnostics = FALSE)
  output <- capture.output(print(out))

  expect_summary_contract(out, fits[[name]], name)
  expect_true(any(grepl("sd_multiplier(", output, fixed = TRUE)))
  expect_false(any(grepl("group_covariance", output, fixed = TRUE)))
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
