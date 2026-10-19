test_that("fitted and ranef reject unsupported and method-controlled dots", {

  object <- .review_cleanup_plain_replay()
  expect_identical(fitted(object, quiet = TRUE), fitted(object, quiet = FALSE))
  for (argument in c("type", "unrecognized", ".posterior_samples")) {
    args <- stats::setNames(list(TRUE), argument)
    expect_error(do.call(fitted, c(list(object = object), args)),
      paste0("Unused argument in fitted(): '", argument, "'"), fixed = TRUE)
  }
  for (argument in c("newdata", "type", "quiet")) {
    for (value in list(TRUE, NULL)) {
      args <- stats::setNames(list(value), argument)
      expect_error(do.call(ranef, c(list(object = object), args)),
        paste0("Unused argument in ranef(): '", argument, "'"), fixed = TRUE)
    }
  }
  expect_s3_class(ranef(object), "brma_samples")
})

test_that("analogous unique-level refusals retain the existing expansion condition", {

  slope <- review_scaled_allocation_replay(standardized = FALSE)
  error <- tryCatch(ranef(slope$object), error = identity)
  expect_s3_class(error, "RoBMA_ranef_expand_required")
  expect_s3_class(error, "RoBMA_unsupported_request")
  expect_identical(conditionMessage(error), paste0(
    "Unique-level random effects are unavailable for random-slope block 'study'. ",
    "Use 'expand = TRUE'."))

  multiple <- review_heterogeneity_parent_replay("mean_variance")
  error <- tryCatch(ranef(multiple$object, component = "total"), error = identity)
  expect_s3_class(error, "RoBMA_ranef_expand_required")
  expect_s3_class(error, "RoBMA_unsupported_request")
  expect_identical(conditionMessage(error), paste0(
    "'component = \"total\"' across multiple random-effect blocks requires ",
    "'expand = TRUE'."))

  # Deliberately missing declared labels is a metadata guard control, not an
  # admitted fitted-model numerical failure.
  term <- review_heterogeneity_har_replay(direct = TRUE, n = 24L)$term
  term$column_names <- NULL
  error <- tryCatch(.ranef_unique_level_term(term, "study"), error = identity)
  expect_s3_class(error, "RoBMA_ranef_expand_required")
  expect_s3_class(error, "RoBMA_unsupported_request")
  expect_identical(conditionMessage(error), paste0(
    "Unique-level random effects are unavailable for block 'study' ",
    "because its semantic coefficient labels are missing. Use 'expand = TRUE'."))
})
