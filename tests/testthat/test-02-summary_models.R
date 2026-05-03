context("Summary models")

source(testthat::test_path("common-functions.R"))

skip_if_no_fits()
fit_names <- list_fits()
fits      <- lazy_fits(fit_names, validate = FALSE)

summary_model_names <- intersect(
  names(fits),
  catalog_fits(class = c("BMA.norm", "BMA.glmm", "RoBMA"))
)

expect_summary_models_marginal_table <- function(table, name, component) {

  info <- paste0("summary_models marginal '", component, "' for '", name, "'")

  expect_s3_class(table, "RoBMA_summary_models_marginal")
  expect_true(is.data.frame(table), info = info)
  expect_true(nrow(table) > 0L, info = info)
  expect_true(all(c("Hypothesis", "prior_prob", "post_prob", "inclusion_BF") %in%
                    colnames(table)), info = info)
  expect_true(all(table[["Hypothesis"]] %in% c("Null", "Alternative")),
              info = info)
  expect_true(all(is.finite(table[["prior_prob"]])), info = info)
  expect_true(all(is.finite(table[["post_prob"]])), info = info)
  expect_equal(sum(table[["prior_prob"]]), 1, tolerance = sqrt(.Machine$double.eps),
               info = info)
  expect_equal(sum(table[["post_prob"]]), 1, tolerance = sqrt(.Machine$double.eps),
               info = info)
  expect_equal(attr(table, "title"), component, info = info)
}

expect_summary_models_marginal <- function(out, name) {

  expect_s3_class(out, "summary_models.RoBMA")
  expect_equal(out[["type"]], "marginal")
  expect_type(out[["name"]], "character")
  expect_true(length(out[["name"]]) == 1L && nzchar(out[["name"]]))
  expect_type(out[["marginal"]], "list")
  expect_true(length(out[["marginal"]]) > 0L,
              info = paste0("marginal components for '", name, "'"))

  for (component in names(out[["marginal"]])) {
    expect_summary_models_marginal_table(
      out[["marginal"]][[component]],
      name      = name,
      component = component
    )
  }
}

expect_summary_models_individual <- function(out, marginal, name) {

  table <- out[["individual"]]
  info  <- paste0("summary_models individual for '", name, "'")

  expect_s3_class(out, "summary_models.RoBMA")
  expect_equal(out[["type"]], "individual")
  expect_s3_class(table, "RoBMA_summary_models_individual")
  expect_true(is.data.frame(table), info = info)
  expect_true(nrow(table) > 0L, info = info)
  expect_true(all(c("prior_prob", "post_prob", "inclusion_BF") %in%
                    colnames(table)), info = info)
  expect_true(all(names(marginal[["marginal"]]) %in% colnames(table)),
              info = info)
  expect_equal(
    nrow(table),
    prod(vapply(marginal[["marginal"]], nrow, integer(1))),
    info = info
  )
  expect_true(all(is.finite(table[["prior_prob"]])), info = info)
  expect_true(all(is.finite(table[["post_prob"]])), info = info)
  expect_equal(sum(table[["prior_prob"]]), 1, tolerance = sqrt(.Machine$double.eps),
               info = info)
  expect_equal(sum(table[["post_prob"]]), 1, tolerance = sqrt(.Machine$double.eps),
               info = info)
  expect_equal(attr(table, "title"), "Individual Models", info = info)
}

expect_printed_summary_models <- function(out, name) {

  output <- capture.output(print(out))
  expect_true(any(nzchar(output)),
              info = paste0("printed summary_models for '", name, "'"))
  expect_false(any(grepl("__xXx__", output, fixed = TRUE)),
               info = paste0("printed summary_models labels for '", name, "'"))
}


test_that("summary_models marginal summaries have stable structure", {

  for (name in summary_model_names) {
    out <- summary_models(fits[[name]], type = "marginal")

    expect_summary_models_marginal(out, name)
    expect_printed_summary_models(out, name)
  }
})

test_that("summary_models individual summaries have stable structure", {

  for (name in summary_model_names) {
    marginal   <- summary_models(fits[[name]], type = "marginal")
    individual <- summary_models(fits[[name]], type = "individual")

    expect_summary_models_individual(individual, marginal, name)
    expect_printed_summary_models(individual, name)
  }
})

test_that("summary_models individual weights marginalize to component weights", {

  for (name in summary_model_names) {
    marginal   <- summary_models(fits[[name]], type = "marginal")
    individual <- summary_models(fits[[name]], type = "individual")[["individual"]]

    for (component in names(marginal[["marginal"]])) {
      component_table <- marginal[["marginal"]][[component]]

      for (level in rownames(component_table)) {
        selected <- individual[[component]] == level

        expect_equal(
          sum(individual[["prior_prob"]][selected]),
          component_table[level, "prior_prob"],
          tolerance = sqrt(.Machine$double.eps),
          info      = paste(name, component, level, "prior marginalization")
        )
        expect_equal(
          sum(individual[["post_prob"]][selected]),
          component_table[level, "post_prob"],
          tolerance = sqrt(.Machine$double.eps),
          info      = paste(name, component, level, "posterior marginalization")
        )
      }
    }
  }
})

test_that("summary_models is RoBMA-only", {

  skip_if_missing_fits("bcg_meta-analysis")

  expect_error(
    summary_models(fits[["bcg_meta-analysis"]]),
    "RoBMA objects"
  )
})

test_that("summary_models decodes interaction labels", {

  skip_if_missing_fits("dat.lehmann2018_RoBMA_mods2")

  marginal_output <- capture.output(print(summary_models(
    fits[["dat.lehmann2018_RoBMA_mods2"]],
    type = "marginal"
  )))
  individual_output <- capture.output(print(summary_models(
    fits[["dat.lehmann2018_RoBMA_mods2"]],
    type = "individual"
  )))

  expect_false(any(grepl("__xXx__", marginal_output, fixed = TRUE)))
  expect_false(any(grepl("__xXx__", individual_output, fixed = TRUE)))
  expect_true(any(grepl("Preregistered:Gender", marginal_output, fixed = TRUE)))
  expect_true(any(grepl("Preregistered:Gender", individual_output, fixed = TRUE)))
})
