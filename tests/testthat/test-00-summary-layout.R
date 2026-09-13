context("Multivariate summary layout")

.summary_layout_estimates <- function(parameters, title, values,
                                      diagnostics = NULL,
                                      random_metadata = FALSE) {

  out <- data.frame(
    Mean        = values,
    SD          = rep(0.0123456789, length(parameters)),
    `0.025`     = values - 0.1,
    `0.975`     = values + 0.1,
    row.names   = parameters,
    check.names = FALSE
  )
  types <- rep("estimate", ncol(out))
  if (!is.null(diagnostics)) {
    out[[diagnostics]] <- rep(123.456789, nrow(out))
    types <- c(types, diagnostics)
  }
  if (random_metadata) {
    out[["Random name"]] <- rep("study", nrow(out))
    out[["Random grouping"]] <- rep("study", nrow(out))
    out[["Random structure"]] <- rep("diag", nrow(out))
    types <- c(types, rep("string", 3L))
  }
  class(out) <- c("BayesTools_table", "BayesTools_runjags_summary", "data.frame")
  attr(out, "type") <- types
  attr(out, "parameters") <- paste0("parameter:", parameters)
  attr(out, "rownames") <- TRUE
  attr(out, "title") <- title
  attr(out, "footnotes") <- paste(title, "footnote.")
  attr(out, "warnings") <- paste(title, "warning.")
  out
}

.summary_layout_fixture <- function(mods = FALSE, scale = FALSE,
                                    conditional = TRUE) {

  sections <- c(
    "inclusion_components", "inclusion_mods", "inclusion_scale",
    "inclusion_random", "estimates", "estimates_conditional",
    "estimates_mods", "estimates_mods_conditional", "estimates_scale",
    "estimates_scale_conditional", "estimates_random",
    "estimates_random_conditional", "estimates_bias",
    "estimates_bias_conditional"
  )
  out <- stats::setNames(rep(list(list()), length(sections)), sections)
  out[["name"]] <- "Bayesian Multivariate Meta-Analysis"
  for (suffix in c("", if (conditional) "_conditional")) {
    prefix <- if (nzchar(suffix)) "Conditional " else ""
    out[[paste0("estimates", suffix)]] <- .summary_layout_estimates(
      "rho", paste0(prefix, if (mods || scale) "Common Estimates" else "Estimates"),
      0.123456789123,
      diagnostics = if (!nzchar(suffix)) "R_hat"
    )
    out[[paste0("estimates_mods", suffix)]] <- .summary_layout_estimates(
      if (mods) c("intercept", "x") else "mu",
      paste0(prefix, "Location"),
      if (mods) c(0.234567891234, 0.345678912345) else 0.234567891234
    )
    out[[paste0("estimates_random", suffix)]] <- .summary_layout_estimates(
      c("tau_total", "tau2_prop(study)"), paste0(prefix, "Random"),
      c(0.456789123456, 0.567891234567),
      diagnostics = if (!nzchar(suffix)) "ESS",
      random_metadata = TRUE
    )
    if (scale) {
      out[[paste0("estimates_scale", suffix)]] <- .summary_layout_estimates(
        c("exp(intercept)", "z"), paste0(prefix, "Scale"),
        c(0.678912345678, 0.789123456789)
      )
    }
    out[[paste0("estimates_bias", suffix)]] <- .summary_layout_estimates(
      "omega[0.05,1]", paste0(prefix, "Publication Bias"), 0.891234567891
    )
  }
  out[["inclusion_components"]] <- data.frame(
    prior_prob   = 0.5,
    post_prob    = 0.75,
    inclusion_BF = 3,
    row.names    = "Effect"
  )
  class(out[["inclusion_components"]]) <- c("BayesTools_table", "data.frame")
  attr(out[["inclusion_components"]], "type") <-
    c("prior_prob", "post_prob", "inclusion_BF")
  attr(out[["inclusion_components"]], "title") <- "Component Inclusion"
  class(out) <- "summary.brma"
  attr(out, "mods") <- mods
  attr(out, "scale") <- scale
  attr(out, "random") <- TRUE
  out
}

test_that("multivariate estimate layout follows moderator and scale sections", {

  for (mods in c(FALSE, TRUE)) {
    for (scale in c(FALSE, TRUE)) {
      out <- .summary_layout_fixture(mods, scale)
      original <- out
      prepared <- .summary_brma_prepare_print_sections(out)
      common_parameters <- c("rho", if (!mods) "mu", "tau_total", "tau2_prop(study)")
      common_title <- if (mods || scale) "Common Estimates" else "Estimates"

      for (suffix in c("", "_conditional")) {
        prefix <- if (nzchar(suffix)) "Conditional " else ""
        estimates <- prepared[[paste0("estimates", suffix)]]
        expect_identical(rownames(estimates), common_parameters)
        expect_identical(attr(estimates, "title"), paste0(prefix, common_title))
        expect_length(prepared[[paste0("estimates_random", suffix)]], 0L)
        if (mods) {
          expected_location <- out[[paste0("estimates_mods", suffix)]]
          attr(expected_location, "title") <- paste0(
            prefix, if (scale) "Location" else "Meta-Regression"
          )
          expect_identical(
            prepared[[paste0("estimates_mods", suffix)]],
            expected_location
          )
        } else {
          expect_length(prepared[[paste0("estimates_mods", suffix)]], 0L)
        }
        expect_identical(
          prepared[[paste0("estimates_scale", suffix)]],
          out[[paste0("estimates_scale", suffix)]]
        )
        expect_identical(
          prepared[[paste0("estimates_bias", suffix)]],
          out[[paste0("estimates_bias", suffix)]]
        )
      }
      expect_identical(prepared[["inclusion_components"]], out[["inclusion_components"]])

      output <- capture.output(print(out))
      expect_true(any(output == common_title))
      expect_true(any(output == paste0("Conditional ", common_title)))
      expect_false(any(output %in% c("Random", "Conditional Random")))
      expect_equal(any(output == "Meta-Regression"), mods && !scale)
      expect_equal(any(output == "Location"), mods && scale)
      expect_equal(any(output == "Scale"), scale)
      expect_true(any(output == "Publication Bias"))
      expect_true(any(output == "Component Inclusion"))

      frame <- as.data.frame(out)
      expect_identical(data.frame(out), frame)
      expect_identical(
        frame[["parameter"]][frame[["component"]] == "common"],
        common_parameters
      )
      expect_identical(
        frame[["parameter"]][frame[["component"]] == "conditional common"],
        common_parameters
      )
      expect_equal(any(frame[["component"]] == "location"), mods)
      expect_equal(any(frame[["component"]] == "scale"), scale)
      expect_false(any(frame[["component"]] %in% c("random", "conditional random")))
      expect_identical(out, original)
    }
  }
})

test_that("combined estimates retain formatting metadata and full precision", {

  out <- .summary_layout_fixture()
  prepared <- .summary_brma_prepare_print_sections(out)
  estimates <- prepared[["estimates"]]
  source_sections <- c("estimates", "estimates_mods", "estimates_random")
  source_tables <- out[source_sections]
  expected_parameters <- unlist(lapply(source_tables, attr, which = "parameters"),
                                use.names = FALSE)
  expect_s3_class(estimates, "BayesTools_table")
  expect_identical(attr(estimates, "parameters"), expected_parameters)
  expect_identical(attr(estimates, "rownames"), TRUE)
  expect_identical(attr(estimates, "type"), c(rep("estimate", 4L), "R_hat", "ESS"))
  for (attribute in c("footnotes", "warnings")) {
    expect_identical(
      attr(estimates, attribute),
      unlist(lapply(source_tables, attr, which = attribute), use.names = FALSE)
    )
  }
  expect_identical(names(estimates), c("Mean", "SD", "0.025", "0.975", "R_hat", "ESS"))
  expect_true(all(is.na(estimates[["ESS"]][1:2])))
  expect_true(all(is.na(estimates[["R_hat"]][-1L])))
  expect_true(all(vapply(estimates, is.numeric, TRUE)))

  conditional <- prepared[["estimates_conditional"]]
  expect_false(any(c("ESS", "R_hat") %in% names(conditional)))
  expect_identical(attr(conditional, "type"), rep("estimate", 4L))
  expect_identical(attr(conditional, "parameters"), expected_parameters)
  expect_identical(
    attr(conditional, "warnings"),
    unlist(lapply(out[paste0(source_sections, "_conditional")], attr,
                  which = "warnings"), use.names = FALSE)
  )

  frame <- as.data.frame(out)
  common <- frame[frame[["component"]] == "common", , drop = FALSE]
  expected_values <- unlist(lapply(source_tables, `[[`, "Mean"), use.names = FALSE)
  expect_identical(common[["Mean"]], expected_values)
  expect_identical(common[["CI_0.025"]], expected_values - 0.1)
  expect_identical(common[["CI_0.975"]], expected_values + 0.1)
  expect_false(any(c("Random name", "Random grouping", "Random structure") %in% names(frame)))
  expect_true(is.numeric(frame[["inclusion_BF"]]))
  expect_identical(frame[["inclusion_BF"]][frame[["component"]] == "inclusion"], 3)
})

test_that("empty random estimates still consolidate an intercept-only location", {

  out <- .summary_layout_fixture(conditional = FALSE)
  out[["estimates"]] <- list()
  out[["estimates_random"]] <- list()
  original <- out
  prepared <- .summary_brma_prepare_print_sections(out)

  expect_identical(rownames(prepared[["estimates"]]), "mu")
  expect_identical(attr(prepared[["estimates"]], "title"), "Estimates")
  expect_length(prepared[["estimates_mods"]], 0L)
  expect_length(prepared[["estimates_random"]], 0L)
  expect_length(prepared[["estimates_conditional"]], 0L)
  expect_false(any(as.data.frame(out)[["component"]] == "conditional common"))
  expect_identical(out, original)

  out[["estimates_mods"]] <- list()
  empty <- .summary_brma_prepare_print_sections(out)
  expect_length(empty[["estimates"]], 0L)
  expect_false(any(capture.output(print(out)) == "Estimates"))
})

test_that("empty estimate tables do not hide labels of populated sections", {

  out <- .summary_layout_fixture(conditional = FALSE)
  out[["estimates"]] <- BayesTools::runjags_estimates_empty_table(
    probs = c(0.025, 0.975), remove_diagnostics = TRUE, title = "Estimates"
  )
  prepared <- .summary_brma_prepare_print_sections(out)

  expect_identical(
    rownames(prepared[["estimates"]]),
    c("mu", "tau_total", "tau2_prop(study)")
  )
  expect_identical(attr(prepared[["estimates"]], "rownames"), TRUE)
  expect_true(any(grepl("tau_total", capture.output(print(out)), fixed = TRUE)))
})

test_that("models without multivariate random terms keep their existing layout", {

  out <- .summary_layout_fixture(mods = TRUE, scale = TRUE, conditional = FALSE)
  out[["estimates_random"]] <- list()
  attr(out, "random") <- FALSE
  expect_identical(.summary_brma_prepare_print_sections(out), out)

  attr(out, "random") <- NULL
  expect_identical(.summary_brma_prepare_print_sections(out), out)
})
