test_that("fitted and ranef warn about unsupported and method-controlled dots", {

  object <- .review_cleanup_plain_replay()
  expect_identical(fitted(object, quiet = TRUE), fitted(object, quiet = FALSE))
  for (argument in c("type", "unrecognized", ".posterior_samples")) {
    args <- stats::setNames(list(TRUE), argument)
    expect_warning(result <- do.call(fitted, c(list(object = object), args)),
      paste0("Unused argument in fitted(): '", argument, "'"), fixed = TRUE)
    expect_identical(result, fitted(object))
  }
  for (argument in c("newdata", "type", "quiet")) {
    for (value in list(TRUE, NULL)) {
      args <- stats::setNames(list(value), argument)
      expect_warning(do.call(ranef, c(list(object = object), args)),
        paste0("Unused argument in ranef(): '", argument, "'"), fixed = TRUE)
    }
  }
  expect_s3_class(ranef(object), "brma_samples")
})

test_that("blup, true_effects and ranef warn once under their own name and forward only predict arguments", {

  object <- .review_cleanup_plain_replay()
  for (fun in c("blup", "true_effects", "ranef")) {
    messages <- character()
    value <- withCallingHandlers(
      do.call(fun, list(object = object, unrecognized = TRUE)),
      warning = function(w) {
        messages <<- c(messages, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    # one warning, named after the method that was called, not after the
    # predict() calls inside it (ranef() makes several)
    expect_identical(messages, paste0("Unused argument in ", fun, "(): 'unrecognized'"), info = fun)
    expect_equal(value, do.call(fun, list(object = object)), info = fun)
    # predict() arguments that the method does not fix are forwarded silently
    expect_no_warning(do.call(fun, list(object = object, as_measure = TRUE)))
  }
  # an explicit NULL output_measure or transform of ranef() is silent, as at 4.0.0
  for (argument in c("output_measure", "transform")) {
    expect_no_warning(value <- do.call(ranef,
      stats::setNames(list(object, NULL), c("object", argument))))
    expect_equal(value, ranef(object))
  }
})

test_that("the arguments that a method fixes are R's error, only unknown arguments warn", {

  # As at 4.0.0, blup(), true_effects(), pooled_effect() and
  # pooled_heterogeneity() set `newdata`, `type` and `quiet` themselves, so one
  # supplied in `...` is an argument matched twice and R raises its own error;
  # a fitted-data result must not stand in for a requested new-data one.
  object <- .review_cleanup_plain_replay()
  for (fun in c("blup", "true_effects", "pooled_effect", "pooled_heterogeneity")) {
    for (argument in c("newdata", "type", "quiet")) {
      expect_error(
        do.call(fun, stats::setNames(list(object, "terms"), c("object", argument))),
        paste0("formal argument \"", argument, "\" matched by multiple actual arguments"),
        fixed = TRUE, info = paste(fun, argument)
      )
    }
  }
  expect_error(blup(object, newdata = data.frame(x = 1)),
    "formal argument \"newdata\" matched by multiple actual arguments", fixed = TRUE)
  expect_error(pooled_effect(object, type = "terms"),
    "formal argument \"type\" matched by multiple actual arguments", fixed = TRUE)

  # an unknown argument is unused and warns once under the name of the method
  for (fun in c("pooled_effect", "pooled_heterogeneity")) {
    expect_warning(
      do.call(fun, list(object = object, unrecognized = TRUE)),
      paste0("Unused argument in ", fun, ".brma(): 'unrecognized'"), fixed = TRUE
    )
  }
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
