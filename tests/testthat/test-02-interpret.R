context("Interpret")

source(testthat::test_path("common-functions.R"))

REFERENCE_DIR <<- testthat::test_path("..", "results", "interpret")

skip_if_no_fits()

interpret_fit_names <- c(
  "bcg_meta-analysis",
  "bcg_BMA.glmm",
  "dat.lehmann2018_RoBMA",
  "dat.lehmann2018_RoBMA_mods2",
  "dat.lehmann2018_RoBMA_3lvl_mods_scale",
  "bangertdrowns2004_location-scale",
  "konstantopoulos2011_3lvl"
)
skip_if_missing_fits(interpret_fit_names)

fits <- lazy_fits(interpret_fit_names, validate = FALSE)

interpret_output_text <- function(object, ...) {

  return(paste(capture.output(print(interpret(object, ...))), collapse = "\n"))
}

expect_interpret_contract <- function(output, name, priors = FALSE) {

  expect_s3_class(output, "interpret.brma")
  expect_type(unclass(output), "character")
  expect_true(length(output) > 0L, info = paste0("non-empty output for ", name))
  expect_true(all(nzchar(output)), info = paste0("non-empty sections for ", name))
  expect_false(any(grepl("__xXx__", output, fixed = TRUE)),
               info = paste0("decoded labels for ", name))
  expect_true("model" %in% names(output), info = paste0("model section for ", name))

  printed <- capture.output(print(output))
  expect_true(any(grepl("Bayesian", printed, fixed = TRUE)),
              info = paste0("printed model name for ", name))
  expect_false(any(grepl("__xXx__", printed, fixed = TRUE)),
               info = paste0("printed decoded labels for ", name))
  if (priors) {
    expect_true(any(printed == ""), info = paste0("prior separator for ", name))
  } else {
    expect_false(any(printed == ""), info = paste0("no empty printed lines for ", name))
  }
  expect_false(any(grepl("prior =|posterior =|inverse =", printed)),
               info = paste0("compact BF reporting for ", name))

  records <- attr(output, "records")
  expect_s3_class(records, "BayesTools_interpret_records")
  expect_true(all(c(
    "record_id", "kind", "section", "item_id", "source", "row",
    "BF_canonical_value", "BF_canonical_bound_operator",
    "central_name", "central_value"
  ) %in% colnames(records)), info = paste0("records schema for ", name))
  expect_true("model.header.header" %in% records[["record_id"]],
              info = paste0("records header for ", name))
}

test_that("interpret.brma returns a stable object contract", {

  for (name in names(fits)) {
    output <- interpret(fits[[name]])
    expect_interpret_contract(output, name)
  }
})

test_that("interpret.brma handles options and errors", {

  output <- interpret(
    fits[["dat.lehmann2018_RoBMA_mods2"]],
    scope       = c("components", "moderators"),
    central     = "mean",
    probs       = c(.10, .90),
    digits      = 2
  )
  expect_interpret_contract(output, "dat.lehmann2018_RoBMA_mods2 options")
  expect_true(any(grepl("80% CrI", output, fixed = TRUE)))
  expect_true(any(grepl("mean", output, fixed = TRUE)))
  expect_true(any(grepl("Location inclusion", output, fixed = TRUE)))

  brma_output <- interpret(fits[["bcg_meta-analysis"]])
  expect_true(any(grepl("Pooled effect: mean", brma_output, fixed = TRUE)))
  expect_false(any(grepl("model-averaged", brma_output, fixed = TRUE)))

  robma_output <- interpret(fits[["dat.lehmann2018_RoBMA"]])
  expect_true(any(grepl("Pooled effect: model-averaged mean", robma_output,
                        fixed = TRUE)))
  expect_true(any(grepl("BF01 =", robma_output, fixed = TRUE)))
  robma_records <- attr(robma_output, "records")
  heterogeneity_record <- robma_records[
    robma_records[["record_id"]] == "components.heterogeneity.evidence",
    ,
    drop = FALSE
  ]
  expect_equal(heterogeneity_record[["BF_canonical_bound_operator"]], ">")
  expect_true(is.finite(heterogeneity_record[["BF_canonical_value"]]))

  exp_output <- interpret(fits[["bcg_BMA.glmm"]], transform = "EXP")
  expect_true(any(grepl("Pooled effect: model-averaged mode", exp_output,
                        fixed = TRUE)))
  expect_true(any(grepl("Pooled effect estimates are summarized", exp_output,
                        fixed = TRUE)))

  transformed_mods <- interpret(
    fits[["dat.lehmann2018_RoBMA_mods2"]],
    scope          = "moderators",
    output_measure = "COR"
  )
  expect_false(any(grepl("on the correlation scale", transformed_mods,
                         fixed = TRUE)))

  prior_output <- interpret(fits[["dat.lehmann2018_RoBMA"]], priors = TRUE)
  expect_interpret_contract(prior_output, "dat.lehmann2018_RoBMA priors",
                            priors = TRUE)
  prior_text <- attr(prior_output, "priors")
  expect_true(any(grepl("Prior distributions", prior_text, fixed = TRUE)))
  expect_true(any(grepl("mu:", prior_text, fixed = TRUE)))

  prior_printed <- capture.output(print(prior_output))
  expect_true(any(prior_printed == ""))
  expect_equal(sum(prior_printed == ""), 1L)
  expect_true(any(grepl("Prior distributions", prior_printed, fixed = TRUE)))

  expect_warning(
    legacy <- interpret(fits[["dat.lehmann2018_RoBMA"]], output_scale = "r"),
    "deprecated"
  )
  expect_true(any(grepl("correlation", legacy, fixed = TRUE)))

  expect_error(interpret(1), "brma objects")
  expect_error(
    interpret(fits[["bcg_meta-analysis"]], conditional = TRUE),
    "RoBMA objects"
  )
  expect_error(
    interpret(fits[["bcg_meta-analysis"]], unused = TRUE),
    "Unused argument"
  )
})


test_that("interpret.brma summarizes heterogeneity component lists separately", {

  total_samples <- .new_brma_samples(
    matrix(
      c(0.10, 0.12, 0.14, 0.16),
      ncol     = 1,
      dimnames = list(NULL, "tau")
    ),
    n_chains = 1,
    n_iter   = 4,
    title    = "total heterogeneity"
  )
  study_samples <- .new_brma_samples(
    matrix(
      c(0.20, 0.22, 0.24, 0.26),
      ncol     = 1,
      dimnames = list(NULL, "tau")
    ),
    n_chains = 1,
    n_iter   = 4,
    title    = "study heterogeneity"
  )
  sources <- .interpret_add_heterogeneity_sources(
    sources       = list(),
    heterogeneity = list(
      "total heterogeneity" = total_samples,
      study                 = study_samples
    ),
    conditional   = FALSE,
    probs         = c(.025, .975),
    central       = "mean"
  )
  plan <- .interpret_brma_plan(
    summary_object   = list(inclusion_components = NULL),
    effect_transform = NULL,
    scope            = "estimates",
    sources          = sources
  )
  records <- .interpret_records_standard(sources = sources, plan = plan)
  text    <- .interpret_records_text(records, digits = 3, averaged = FALSE)

  expect_true("estimates.heterogeneity.total.heterogeneity.estimate" %in%
    records[["record_id"]])
  expect_true("estimates.heterogeneity.study.estimate" %in%
    records[["record_id"]])
  expect_true(any(grepl("Pooled heterogeneity (total heterogeneity)",
                        text, fixed = TRUE)))
  expect_true(any(grepl("Pooled heterogeneity (study)", text, fixed = TRUE)))
})


test_that("interpret.brma printed output matches reference text", {

  cases <- list(
    list(
      name     = "bcg_meta-analysis",
      filename = "interpret-bcg_meta-analysis-core.txt",
      args     = list()
    ),
    list(
      name     = "bcg_BMA.glmm",
      filename = "interpret-bcg_BMA.glmm-exp.txt",
      args     = list(transform = "EXP")
    ),
    list(
      name     = "dat.lehmann2018_RoBMA",
      filename = "interpret-dat.lehmann2018_RoBMA-core.txt",
      args     = list()
    ),
    list(
      name     = "dat.lehmann2018_RoBMA",
      filename = "interpret-dat.lehmann2018_RoBMA-core-priors.txt",
      args     = list(priors = TRUE)
    ),
    list(
      name     = "dat.lehmann2018_RoBMA",
      filename = "interpret-dat.lehmann2018_RoBMA-bias.txt",
      args     = list(scope = "bias")
    ),
    list(
      name     = "dat.lehmann2018_RoBMA_mods2",
      filename = "interpret-dat.lehmann2018_RoBMA_mods2-all-conditional.txt",
      args     = list(scope = "all", conditional = TRUE)
    ),
    list(
      name     = "dat.lehmann2018_RoBMA_3lvl_mods_scale",
      filename = "interpret-dat.lehmann2018_RoBMA_3lvl_mods_scale-core.txt",
      args     = list()
    ),
    list(
      name     = "bangertdrowns2004_location-scale",
      filename = "interpret-bangertdrowns2004_location-scale-all.txt",
      args     = list(scope = "all")
    ),
    list(
      name     = "konstantopoulos2011_3lvl",
      filename = "interpret-konstantopoulos2011_3lvl-core.txt",
      args     = list()
    )
  )

  for (case in cases) {
    output <- do.call(
      interpret_output_text,
      c(list(object = fits[[case[["name"]]]]), case[["args"]])
    )

    test_reference_text(
      text     = output,
      filename = case[["filename"]],
      info_msg = paste0("Interpret reference mismatch for ", case[["name"]])
    )
  }
})
