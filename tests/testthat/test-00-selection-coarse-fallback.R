context("Fallback and per-fit record of the corrected selection sampler")
source(testthat::test_path("common-functions.R"))

# Each corrected sampler (one per parameter and chain, per model build) chooses
# its kernel once, at the end of adaptation: with the fallback enabled
# ('selection.coarse_min_acceptance' > 0) it switches to the exact slice kernel
# when (a) it made fewer than 50 adaptive proposals, (c) no surrogate it
# evaluated while adapting differed from the exact density, or (b) its adaptive
# acceptance rate is below the threshold, the first rule that holds being the
# recorded reason. Its post-adaptation work is counted per fit in every process
# that runs chains, and a low rate warns. The R-level record, warning and option
# logic run here; hand-written models and small fits run in the certification
# case 'numerical-kernels' (helpers in helper-selection-coarse-exactness.R).

.fallback_report_names <- c("samplers", "fallback_few_proposals", "fallback_surrogate_exact",
  "fallback_low_acceptance", "sampling_proposed", "sampling_accepted", "sampling_judged", "min_sampler_rate")

.fallback_record_names <- c("sampler", "min_acceptance", "samplers", "fallback_samplers",
  "fallback_few_proposals", "fallback_surrogate_exact", "fallback_low_acceptance",
  "sampling_proposed", "sampling_accepted", "sampling_judged", "acceptance_rate", "min_sampler_rate", "complete")

# One process's report as the native routine returns it; the judged
# acceptances of scalar samplers are their accepted proposals.
.fallback_report_value <- function(samplers = 0, few = 0, exact = 0, low = 0, proposed = 0,
                                   accepted = 0, rate = NA_real_, judged = accepted) {

  stats::setNames(list(samplers, few, exact, low, proposed, accepted, judged, rate),
                  .fallback_report_names)
}

.fallback_settings <- list(sampler = "slice.coarse", coarse_min_acceptance = .9)

.fallback_count_names_record <- c("samplers", "fallback_samplers", "fallback_few_proposals",
  "fallback_surrogate_exact", "fallback_low_acceptance", "sampling_proposed", "sampling_accepted",
  "sampling_judged")


test_that("the corrected sampler is the default and its threshold is a number between 0 and 1", {

  old_options <- RoBMA.options()
  on.exit(do.call(RoBMA.options, old_options), add = TRUE)
  defaults <- lapply(.RoBMA_option_schema[c("selection.sampler", "selection.coarse_min_acceptance")],
                     `[[`, "default")
  expect_identical(defaults, list(selection.sampler = "slice.coarse",
                                  selection.coarse_min_acceptance = 0.90))
  for (value in c(0, .5, 1)) {
    RoBMA.options(selection.coarse_min_acceptance = value)
    expect_identical(RoBMA.get_option("selection.coarse_min_acceptance"), value)
    expect_identical(.selection_runtime_settings()[["coarse_min_acceptance"]], value)
  }

  before <- RoBMA.options()
  name   <- "Option 'selection.coarse_min_acceptance' must be"
  cases  <- list(
    list(-.1, paste(name, ">= 0.")), list(1.1, paste(name, "<= 1.")),
    list(NA_real_, paste(name, "a finite number.")), list(Inf, paste(name, "a finite number.")),
    list("high", paste(name, "a finite number.")), list(c(.1, .2), paste(name, "a finite number."))
  )
  for (case in cases) {
    error <- tryCatch(RoBMA.options(selection.coarse_min_acceptance = case[[1L]]), error = identity)
    expect_s3_class(error, "error")
    expect_identical(conditionMessage(error), case[[2L]])
    expect_null(conditionCall(error))
    expect_identical(RoBMA.options(), before)
  }
})


test_that("the native settings take an optional threshold and the per-fit counters are read and cleared", {

  old_options <- RoBMA.options()
  on.exit(do.call(RoBMA.options, old_options), add = TRUE)
  control <- function(settings) {
    .Call("RoBMA_selnorm_sampler_control", NULL, settings, FALSE, PACKAGE = "RoBMA")
  }
  # settings frozen before the threshold existed have four values
  expect_named(control(c(.05, .05, .01, 3)), c("enabled", "proposed", "accepted", "correction_rejected"))
  expect_named(control(c(.05, .05, .01, 3, .9)), c("enabled", "proposed", "accepted", "correction_rejected"))
  for (threshold in c(-.1, 1.5, NA, NaN, Inf)) {
    expect_error(control(c(.05, .05, .01, 3, threshold)),
                 "The coarse acceptance threshold must be a number between 0 and 1.", fixed = TRUE)
  }
  expect_error(control(c(.05, .05, .01, 3, .9, 1)), paste0("'settings' must contain three grid steps, ",
    "a rule count and an optional acceptance threshold."), fixed = TRUE)

  .selection_sampler_report_read()
  report <- .Call("RoBMA_selnorm_sampler_report", FALSE, PACKAGE = "RoBMA")
  expect_identical(report, .fallback_report_value(rate = NA_real_))
  expect_error(.Call("RoBMA_selnorm_sampler_report", NA, PACKAGE = "RoBMA"),
               "'clear' must be one non-missing logical value.", fixed = TRUE)
})


test_that("the counters of earlier fits are discarded at the start of a fit and nowhere else", {

  old_options <- RoBMA.options()
  on.exit(do.call(RoBMA.options, old_options), add = TRUE)
  reads <- 0L
  testthat::local_mocked_bindings(.selection_sampler_report_read = function() {
    reads <<- reads + 1L
    .fallback_report_value()
  })
  settings <- .selection_runtime_settings()
  context  <- list(phase = "start", role = "local", chains = 2L, processes = 1L,
                   process_id = 1L, process_chains = 2L, parallel = FALSE)
  .selection_runtime_configure(settings, context)
  expect_identical(reads, 1L)
  .selection_runtime_configure(settings, utils::modifyList(context, list(phase = "finish")))
  .selection_runtime_configure(settings)
  RoBMA.options(selection.coarse_min_acceptance = .5)
  expect_identical(reads, 1L)

  # the stored setup closure keeps the body of RoBMA 4.1.1 and discards through
  # the namespace function, so fits saved by 4.1.1 replay it unchanged
  setup <- .selection_runtime_setup()
  expect_identical(deparse(body(setup), control = NULL), deparse(quote({
    utils::getFromNamespace(".selection_runtime_configure", "RoBMA")(settings, context)
    invisible(NULL)
  }), control = NULL))
  expect_invisible(setup(context))
  expect_identical(reads, 2L)
  setup(utils::modifyList(context, list(phase = "finish")))
  expect_identical(reads, 2L)

  # the report callback only reads (and clears) the counters of its process
  report <- .selection_sampler_report_callback()
  expect_identical(environment(report), baseenv())
  expect_identical(report(utils::modifyList(context, list(phase = "report"))), .fallback_report_value())
  expect_identical(reads, 3L)
})


test_that("the record sums the runs of every process and adds the input record once", {

  values <- list(
    list(.fallback_report_value(3, 1, 0, 0, 300, 290, .95), .fallback_report_value(3, 0, 1, 1, 100, 99, NA)),
    list(.fallback_report_value(0, 0, 0, 0, 50, 40, .8), .fallback_report_value(0, 0, 0, 0, 50, 50, NA))
  )
  record <- .selection_sampler_record(values, .fallback_settings)
  expect_identical(record, list(
    sampler = "slice.coarse", min_acceptance = .9, samplers = 6,
    fallback_samplers = 3,
    fallback_few_proposals = 1, fallback_surrogate_exact = 1, fallback_low_acceptance = 1,
    sampling_proposed = 500, sampling_accepted = 479, sampling_judged = 479, acceptance_rate = 479 / 500,
    min_sampler_rate = .8, complete = TRUE
  ))

  # an extension adds the record of its input once; without new runs (a failed
  # extension returns its input) the record is the input's
  extended <- .selection_sampler_record(values[1L], .fallback_settings, previous = record, extension = TRUE)
  expect_identical(unlist(extended[c("samplers", "fallback_samplers", "sampling_proposed", "sampling_accepted")]),
                   c(samplers = 12, fallback_samplers = 6, sampling_proposed = 900, sampling_accepted = 868))
  expect_identical(extended[["min_sampler_rate"]], .8)
  expect_true(extended[["complete"]])
  expect_identical(.selection_sampler_record(NULL, .fallback_settings, previous = record, extension = TRUE,
                                             new_draws = FALSE), record)

  # the worst-sampler rate is NA when no sampler of any run has enough proposals
  none <- .selection_sampler_record(list(list(.fallback_report_value(2, 0, 0, 0, 100, 100))), .fallback_settings)
  expect_identical(none[["min_sampler_rate"]], NA_real_)
  expect_true(none[["complete"]])

  # completeness: FALSE whenever draws exist that no report covers
  # fresh fit with reports: complete
  expect_true(record[["complete"]])
  for (empty in list(NULL, list())) {
    # fresh fit without any report
    fresh <- .selection_sampler_record(empty, .fallback_settings)
    expect_identical(fresh[["acceptance_rate"]], NA_real_)
    expect_false(fresh[["complete"]])
    # extension with new draws and no report (its runs lost every report)
    lost <- .selection_sampler_record(empty, .fallback_settings, previous = record, extension = TRUE,
                                      new_draws = TRUE)
    expect_false(lost[["complete"]])
    expect_identical(lost[c(.fallback_count_names_record)], record[c(.fallback_count_names_record)])
    # failed extension (its input, no new draws): the input's record and state
    expect_identical(.selection_sampler_record(empty, .fallback_settings, previous = record, extension = TRUE,
                                               new_draws = FALSE), record)
    expect_false(.selection_sampler_record(empty, .fallback_settings, extension = TRUE, new_draws = FALSE,
                                           previous = utils::modifyList(record, list(complete = FALSE)))[["complete"]])
  }
  # extension with new draws and reports: the input's state
  expect_true(.selection_sampler_record(values[1L], .fallback_settings, previous = record,
                                        extension = TRUE)[["complete"]])
  expect_false(.selection_sampler_record(values[1L], .fallback_settings, extension = TRUE,
                                         previous = utils::modifyList(record, list(complete = FALSE)))[["complete"]])

  # part of the chain without a count: an extension of a fit without a record,
  # a process whose report failed, or an input that was incomplete
  expect_false(.selection_sampler_record(values, .fallback_settings, previous = NULL, extension = TRUE)[["complete"]])
  expect_false(.selection_sampler_record(list(list(NULL, values[[1L]][[1L]])), .fallback_settings)[["complete"]])
  incomplete <- utils::modifyList(record, list(complete = FALSE))
  expect_false(.selection_sampler_record(values, .fallback_settings, previous = incomplete, extension = TRUE)[["complete"]])

  # settings frozen before the threshold existed (RoBMA 4.1.1) mean 0
  legacy <- .selection_sampler_record(values, list(sampler = "default"))
  expect_identical(legacy[c("sampler", "min_acceptance")], list(sampler = "default", min_acceptance = 0))
})


# A fitted object as .selection_sampler_warn() reads it.
.fallback_warning_fit <- function(record, warnings = NULL) {

  structure(list(), selection_sampler = record, warnings = warnings)
}

.fallback_warning_record <- function(proposed, accepted, rate = NA_real_, threshold = .9, low = 1) {

  .selection_sampler_record(list(list(.fallback_report_value(4, 0, 0, low, proposed, accepted, rate))),
                            list(sampler = "slice.coarse", coarse_min_acceptance = threshold))
}

.fallback_capture_warnings <- function(expr) {

  caught <- list()
  value  <- withCallingHandlers(expr, warning = function(w) {
    caught[[length(caught) + 1L]] <<- w
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = caught)
}


test_that("a low acceptance warns once with its classes and keeps one message with the fit", {

  low    <- .fallback_warning_record(proposed = 1000, accepted = 850, rate = .7)
  result <- .fallback_capture_warnings(.selection_sampler_warn(.fallback_warning_fit(low, "an earlier warning")))
  expect_length(result$warnings, 1L)
  condition <- result$warnings[[1L]]
  expect_identical(class(condition), c("RoBMA_selection_sampler_low_acceptance", "RoBMA_selection_sampler_warning",
                                       "simpleWarning", "warning", "condition"))
  expect_null(conditionCall(condition))
  message <- paste0(
    "The corrected selection sampler accepted 85% of its proposals per direction after adaptation (lowest rate of a ",
    "single sampler 70%), below 'selection.coarse_min_acceptance' = 0.9; 1 of 4 samplers switched to the ",
    "exact slice sampler at the end of adaptation. Sampling is exact but inefficient. Set ",
    "'RoBMA.options(selection.sampler = \"slice\")' before fitting to avoid the corrected sampler; an ",
    "existing fit keeps the settings it was fitted with."
  )
  expect_identical(conditionMessage(condition), message)
  expect_identical(result$value, c("an earlier warning", message))

  # an extension replaces the earlier message by the one of its cumulative record
  worst  <- .fallback_warning_record(proposed = 4000, accepted = 3800, rate = .5)
  result <- .fallback_capture_warnings(.selection_sampler_warn(.fallback_warning_fit(worst, result$value)))
  expect_length(result$warnings, 1L)
  expect_identical(result$value[1L], "an earlier warning")
  expect_length(result$value, 2L)
  expect_match(result$value[2L], "accepted 95% of its proposals per direction after adaptation (lowest rate of a single sampler 50%)",
               fixed = TRUE)

  # no warning, and no message kept: fewer than 1,000 counted proposals, rates at
  # or above the threshold, or a threshold of 0
  quiet <- list(
    .fallback_warning_record(proposed = 999, accepted = 0),
    .fallback_warning_record(proposed = 1000, accepted = 900, rate = .9),
    .fallback_warning_record(proposed = 1000, accepted = 1000),
    .fallback_warning_record(proposed = 5000, accepted = 0, rate = 0, threshold = 0)
  )
  for (record in quiet) {
    result <- .fallback_capture_warnings(.selection_sampler_warn(.fallback_warning_fit(record, c("other", message))))
    expect_length(result$warnings, 0L)
    expect_identical(result$value, "other")
  }
  expect_null(.selection_sampler_warn(.fallback_warning_fit(quiet[[1L]])))
})


test_that("the warning judges a block per direction, pooled with the scalar samplers by their proposals", {

  # a block of 2 members accepting 820 of 1,000 proposals: 0.906 per direction,
  # above the threshold 0.90, so no warning (raw, 0.82, it would warn)
  pair <- .selection_sampler_record(list(list(.fallback_report_value(
    1, 0, 0, 0, 1000, 820, rate = sqrt(.82), judged = 1000 * sqrt(.82)))), .fallback_settings)
  expect_identical(pair[["sampling_accepted"]], 820)
  expect_lt(abs(pair[["acceptance_rate"]] - sqrt(.82)), 1e-15)
  result <- .fallback_capture_warnings(.selection_sampler_warn(.fallback_warning_fit(pair)))
  expect_length(result$warnings, 0L)

  # a block of 5 members at a raw 0.60 (0.903 per direction) and three scalar
  # samplers at 0.95, 1,000 proposals each: the judged rate pools the
  # per-direction rates by proposals, 0.938; the raw rate would be 0.8625
  judged <- 1000 * .6^(1 / 5) + 3 * 950
  mixed  <- .selection_sampler_record(list(list(.fallback_report_value(
    4, 0, 0, 0, 4000, 600 + 3 * 950, rate = .6^(1 / 5), judged = judged))), .fallback_settings)
  expect_identical(mixed[["acceptance_rate"]], judged / 4000)
  expect_gt(mixed[["acceptance_rate"]], .93)
  result <- .fallback_capture_warnings(.selection_sampler_warn(.fallback_warning_fit(mixed)))
  expect_length(result$warnings, 0L)

  # a judged rate below the threshold warns with that rate
  slow <- .selection_sampler_record(list(list(.fallback_report_value(
    1, 0, 0, 0, 1000, 640, rate = .8, judged = 800))), .fallback_settings)
  result <- .fallback_capture_warnings(.selection_sampler_warn(.fallback_warning_fit(slow)))
  expect_length(result$warnings, 1L)
  expect_match(conditionMessage(result$warnings[[1L]]), "accepted 80% of its proposals per direction", fixed = TRUE)

  # the judged acceptances add up over runs, processes and extensions; a report
  # or an input record of an earlier build, without them, counts as scalar
  older <- .fallback_report_value(2, 0, 0, 0, 1000, 950)
  older[["sampling_judged"]] <- NULL
  pooled <- .selection_sampler_record(list(list(older), list(.fallback_report_value(
    1, 0, 0, 0, 1000, 640, judged = 800))), .fallback_settings)
  expect_identical(pooled[["sampling_judged"]], 1750)
  expect_identical(pooled[["acceptance_rate"]], 1750 / 2000)
  input <- pooled[setdiff(names(pooled), "sampling_judged")]
  extended <- .selection_sampler_record(list(list(.fallback_report_value(1, 0, 0, 0, 1000, 640, judged = 800))),
                                        .fallback_settings, previous = input, extension = TRUE)
  expect_identical(extended[["sampling_judged"]], 1590 + 800)
  fitted <- structure(list(fit = structure(list(), selection_sampler = input)), class = "brma")
  expect_identical(selection_sampler_info(fit = fitted)[["sampling_judged"]], 1590)
})


test_that("fallbacks alone do not warn, and the stored fallback message of development versions is dropped", {

  # every sampler of a model-averaged fit may fall back without a warning (the
  # allocation block is declared, so the corrected block sampler or BayesTools'
  # sampler updates it); only a low acceptance warns
  values <- list(list(.fallback_report_value(10, 1, 2, 3, 900, 890, .95)))
  record <- .selection_sampler_record(values, .fallback_settings)
  expect_false("allocation_block_withheld" %in% names(record))
  stored <- paste0(
    "6 of 10 corrected selection samplers switched to the exact slice sampler at the end of adaptation ",
    "(3 for low acceptance, 1 for too few adaptive proposals, 2 whose surrogate was exact). This ",
    "model-averaged fit does not declare the joint block of its heterogeneity allocation under ",
    "'selection.sampler = \"slice.coarse\"', so the allocation SD and weights may now be updated one at a ",
    "time and mix more slowly. Sampling is exact. Set 'RoBMA.options(selection.sampler = \"slice\")' before ",
    "fitting to keep the block; an existing fit keeps the settings it was fitted with."
  )
  # an extension of a fit of a development version refreshes its warnings: the
  # stored message goes, the others stay
  previous <- c(record, list(allocation_block_withheld = TRUE))
  extended <- .selection_sampler_record(values, .fallback_settings, previous = previous, extension = TRUE)
  expect_false("allocation_block_withheld" %in% names(extended))
  result <- .fallback_capture_warnings(.selection_sampler_warn(.fallback_warning_fit(extended,
                                                                                    c("an earlier warning", stored))))
  expect_length(result$warnings, 0L)
  expect_identical(result$value, "an earlier warning")
})


test_that("selection_sampler_info() reports a fit's record, or the process counters without one", {

  record <- .selection_sampler_record(list(list(.fallback_report_value(2, 0, 0, 1, 1200, 1100, .9))),
                                      .fallback_settings)
  fitted <- structure(list(fit = structure(list(), selection_sampler = record)), class = "brma")
  info   <- selection_sampler_info(fit = fitted)
  expect_identical(info, as.data.frame(record, stringsAsFactors = FALSE))
  expect_identical(names(info), .fallback_record_names)
  expect_identical(nrow(info), 1L)

  # a fit without a record: its frozen settings and NA counts, incomplete
  bare <- structure(list(fit = structure(list(), selection_runtime = list(sampler = "default"))), class = "brma")
  info <- selection_sampler_info(fit = bare)
  expect_identical(names(info), .fallback_record_names)
  expect_identical(info[["sampler"]], "default")
  expect_identical(info[["min_acceptance"]], 0)
  expect_true(all(is.na(info[setdiff(.fallback_record_names, c("sampler", "min_acceptance", "complete"))])))
  expect_false(info[["complete"]])

  # a record of a development version with the withdrawn allocation flag: the
  # flag is not reported
  dev_record <- c(record[c("sampler", "min_acceptance")], list(allocation_block_withheld = TRUE),
                  record[setdiff(names(record), c("sampler", "min_acceptance"))])
  dev_fit <- structure(list(fit = structure(list(), selection_sampler = dev_record)), class = "brma")
  expect_identical(selection_sampler_info(fit = dev_fit), as.data.frame(record, stringsAsFactors = FALSE))

  expect_error(selection_sampler_info(clear = TRUE, fit = fitted),
               "'clear' cannot be used with 'fit': a fit's record is not cleared.", fixed = TRUE)
  expect_error(selection_sampler_info(fit = list()), "'fit' must be a fitted 'brma' object or NULL.", fixed = TRUE)
  expect_named(selection_sampler_info(), c("component", "proposed", "accepted", "correction_rejected",
                                           "acceptance_rate"))
})


test_that("errors of the corrected sampler point to the ordinary samplers", {

  for (message in c("Error in node mu\nCorrected slice requires a finite current full density",
                    "Error in node mu\nSurrogate slice transition failed")) {
    hinted <- .sampler_error_hint(message)
    expect_identical(hinted, paste0(
      message, "\nThe error was raised by the corrected selection sampler that RoBMA uses by default. ",
      "Fit the model again after 'RoBMA.options(selection.sampler = \"slice\")' to sample these parameters ",
      "with the ordinary samplers, and please report the model to the RoBMA developers."
    ))
  }
  expect_identical(.sampler_error_hint("Error in node mu\nSlicer stuck"), "Error in node mu\nSlicer stuck")

  # a member whose two bounds are too wide apart for the block samplers' scale:
  # BayesTools' sampler would stop it under "slice" too, so the hint undeclares
  # the blocks
  for (message in c(paste0("Error in node p\nCorrected slice (block): the support bounds of this node are too wide ",
                           "for the sampler's scale (their difference exceeds the largest double)"),
                    paste0("Error in node p\nFactor slice sampler: the support bounds of this node are too wide ",
                           "for the sampler's scale (their difference exceeds the largest double); ",
                           "declare the block without this node"))) {
    hinted <- .sampler_error_hint(message)
    expect_identical(hinted, paste0(
      message, "\nThe error was raised by the sampler of a parameter block that RoBMA uses by default. ",
      "Fit the model again after 'RoBMA.options(sampler.blocks = FALSE)' to sample these parameters with ",
      "the ordinary samplers, and please report the model to the RoBMA developers."
    ))
  }
  # every other error of either block sampler gets the same option; a member on
  # a bound also gets the interior start
  for (message in c("Error in node p\nCorrected slice (block): current value lies on a bound of its support",
                    "Error in node p\nFactor slice sampler: current value lies on a bound of its support")) {
    hinted <- .sampler_error_hint(message)
    expect_true(grepl("RoBMA.options(sampler.blocks = FALSE)", hinted, fixed = TRUE), info = message)
    expect_false(grepl("selection.sampler", hinted, fixed = TRUE), info = message)
    expect_true(grepl("A start strictly inside the support of the parameter avoids this error", hinted, fixed = TRUE),
                info = message)
  }
  block_error <- "Error in node mu\nCorrected slice (block) requires a finite current full density"
  expect_identical(.sampler_error_hint(block_error), paste0(
    block_error, "\nThe error was raised by the sampler of a parameter block that RoBMA uses by default. ",
    "Fit the model again after 'RoBMA.options(sampler.blocks = FALSE)' to sample these parameters with ",
    "the ordinary samplers, and please report the model to the RoBMA developers."
  ))
})


# ---------------------------------------------------------------------------- #
# Hand-written models (certification case 'numerical-kernels').
# ---------------------------------------------------------------------------- #

# The moments of draws against a reference: mean and variance within 4 MCSE
# (plus twice the reference's error allowance) with at least `min_ess`
# effective draws; returns the MCSE.
.fallback_expect_moments <- function(values, reference, reference_error, label, min_ess = 400) {

  observables <- cbind(mean = values, variance = (values - mean(values))^2)
  ess      <- as.numeric(coda::effectiveSize(coda::mcmc(observables)))
  mcse     <- sqrt(apply(observables, 2L, stats::var) / ess)
  estimate <- c(mean(values), stats::var(values))
  expect_true(all(is.finite(mcse)) && all(mcse > 0), info = label)
  expect_true(all(ess >= min_ess), info = label)
  expect_true(all(abs(estimate - reference) <= 4 * mcse + 2 * reference_error), info = label)
  mcse
}

# A stress fixture: its exact and surrogate-only posteriors and two chains of
# the same model and seed, with the fallback off (the sampler stays corrected:
# its achieved acceptance rate, the identity of its counts and its exactness)
# and on at the default threshold (the sampler switches to the exact kernel at
# the end of adaptation, counts nothing after it, and the chain is exact).
.fallback_expect_stress_fixture <- function(fixture, below, draws = 20000L) {

  label <- fixture$label
  off   <- fixture$chain(0, draws)
  on    <- fixture$chain(.9, draws)
  burnin <- 500

  # fallback off: corrected throughout, at the asserted acceptance rate
  report <- off$report
  expect_identical(off$samplers, "RoBMA::CoarseCorrectedSlice", info = label)
  expect_identical(unlist(report[c("samplers", "fallback_few_proposals", "fallback_surrogate_exact",
                                   "fallback_low_acceptance", "sampling_proposed")]),
                   c(samplers = 1, fallback_few_proposals = 0, fallback_surrogate_exact = 0,
                     fallback_low_acceptance = 0, sampling_proposed = burnin + draws), info = label)
  rate <- report[["sampling_accepted"]] / report[["sampling_proposed"]]
  expect_lt(rate, below, label = paste(label, "acceptance rate"))
  expect_identical(off$counters[["proposed"]], as.numeric(draws), info = label)

  reference <- c(fixture$exact$mean, fixture$exact$variance)
  allowance <- c(fixture$exact$mean_error, fixture$exact$variance_error)
  mcse <- .fallback_expect_moments(off$draws, reference, allowance, paste(label, "corrected"))
  # power: the posterior of the surrogate alone is another distribution
  expect_gte(abs(fixture$exact$mean - fixture$surrogate_only$mean), 8 * mcse[1L],
             label = paste(label, "distance of the surrogate-only mean"))

  # fallback on: the sampler switched for its low adaptive acceptance
  expect_identical(on$samplers, "RoBMA::CoarseCorrectedSlice", info = label)
  expect_identical(on$report, .fallback_report_value(samplers = 1, low = 1), info = label)
  expect_identical(on$counters[["proposed"]], 0, info = label)
  .fallback_expect_moments(on$draws, reference, allowance, paste(label, "exact kernel"))

  invisible(list(rate = rate, mcse = mcse))
}

.fallback_cluster_fixture <- function(kind, omega = .05) {

  setup <- .coarse_exactness_setup(kind, omega = omega)
  exact_tracker <- new.env()
  exact_tracker$relative_error <- 0
  exact <- .coarse_exactness_moments(setup, function(mu) {
    .coarse_exactness_quadpack_log_normalizer(setup, mu, exact_tracker)
  }, exact_tracker)
  surrogate_tracker <- new.env()
  surrogate_tracker$relative_error <- 0
  surrogate_only <- .coarse_exactness_moments(setup, function(mu) {
    .coarse_exactness_rule_log_normalizer(setup, mu, setup$orders[3L])
  }, surrogate_tracker)
  list(label = paste(kind, "omega", omega), exact = exact, surrogate_only = surrogate_only,
       chain = function(min_acceptance, draws) {
         .coarse_exactness_chain(setup, draws, min_acceptance = min_acceptance)
       })
}

.fallback_dense_fixture <- function() {

  setup <- .coarse_dense_setup()
  tracker <- new.env()
  tracker$relative_error <- 0
  list(label = "dense", exact = .coarse_dense_moments(setup, surrogate = FALSE, tracker),
       surrogate_only = .coarse_dense_moments(setup, surrogate = TRUE, tracker),
       chain = function(min_acceptance, draws) {
         .coarse_hand_chain(.coarse_dense_hand_model(setup, min_acceptance = min_acceptance), "mu", draws)
       })
}


test_that("stress fixtures fall back at the default threshold and are exact with and without the fallback", {

  skip_if_not_certification("The stress fixtures (about 60 s) run in the numerical-kernels case.")
  skip_if_not_installed("rjags")
  skip_if_not_installed("coda")

  # crude rules of cluster and factor blocks, stronger selection of a cluster
  # block and a crude grid of a dense block: one rate below 0.9 at least, one
  # below 0.5. The corrected chain of the strong selection mixes slowly (about
  # 200 effective draws of the mean in 20,000), so it draws more.
  rates <- c(
    cluster = .fallback_expect_stress_fixture(.fallback_cluster_fixture("cluster"), below = .9)$rate,
    strong  = .fallback_expect_stress_fixture(.fallback_cluster_fixture("cluster", omega = .01), below = .5,
                                              draws = 80000L)$rate,
    factor  = .fallback_expect_stress_fixture(.fallback_cluster_fixture("factor"), below = .9)$rate,
    dense   = .fallback_expect_stress_fixture(.fallback_dense_fixture(), below = .9)$rate
  )
  expect_true(all(rates > 0))
})


# The kernel decision of one cluster sampler: a model with the crude schedule,
# `adapt` adaptive iterations, and the per-fit report after `draws` more.
.fallback_decision <- function(min_acceptance, adapt, draws = 200L, vector_rule = 0, mode = 1,
                               definitions = "") {

  setup <- .coarse_exactness_setup("cluster")
  hand  <- .coarse_hand_model("cluster", 4L, 1L, setup$rules,
    data = list(residual_sd = setup$residual_sd, loading = as.numeric(setup$loading)),
    definitions = paste0("  mu ~ dt(0, 1, 4)\n  for (j in 1:4) { means[j] <- mu }", definitions),
    vector_rule = vector_rule, mode = mode, inits = list(mu = -.25), min_acceptance = min_acceptance)
  .coarse_hand_chain(hand, "mu", draws, adapt = adapt)
}

# the report of one sampler that switched for `reason` (and counted nothing after it)
.fallback_switched <- function(reason) {

  report <- .fallback_report_value(samplers = 1)
  report[[reason]] <- 1
  report
}


test_that("a sampler decides its kernel once, at the end of adaptation, by the first rule that holds", {

  skip_if_not_certification("These hand-written models (about 5 s) run in the numerical-kernels case.")
  skip_if_not_installed("rjags")
  burnin <- 500

  # (a) fewer than 50 adaptive proposals; a threshold that no rate is below
  # isolates the rule
  tiny <- 1e-9
  for (adapt in c(0L, 49L)) {
    expect_identical(.fallback_decision(tiny, adapt)$report, .fallback_switched("fallback_few_proposals"),
                     info = adapt)
  }
  kept <- .fallback_decision(tiny, 50L)$report
  expect_identical(unlist(kept[c("samplers", "fallback_few_proposals", "fallback_surrogate_exact",
                                 "fallback_low_acceptance", "sampling_proposed")]),
                   c(samplers = 1, fallback_few_proposals = 0, fallback_surrogate_exact = 0,
                     fallback_low_acceptance = 0, sampling_proposed = burnin + 200))

  # (b) the adaptive acceptance rate against the threshold: the adaptation of a
  # model does not depend on the threshold, so the rate of its 500 adaptive
  # proposals (the process counters count adaptation) sets the boundary
  setup <- .coarse_exactness_setup("cluster")
  hand  <- .coarse_hand_model("cluster", 4L, 1L, setup$rules,
    data = list(residual_sd = setup$residual_sd, loading = as.numeric(setup$loading)),
    definitions = "  mu ~ dt(0, 1, 4)\n  for (j in 1:4) { means[j] <- mu }", inits = list(mu = -.25))
  RoBMA::selection_sampler_info(clear = TRUE)
  for (batch in seq_len(5L)) rjags::adapt(hand$model, 100L, end.adaptation = FALSE)
  adaptive <- RoBMA::selection_sampler_info()
  expect_identical(adaptive[["proposed"]], 500)
  rate <- adaptive[["accepted"]] / adaptive[["proposed"]]
  expect_true(rate > .1 && rate < .9)
  at <- .fallback_decision(rate, 500L)$report
  expect_identical(at[["fallback_low_acceptance"]], 0)
  expect_identical(at[["sampling_proposed"]], burnin + 200)
  expect_identical(.fallback_decision(rate * (1 + 1e-12), 500L)$report,
                   .fallback_switched("fallback_low_acceptance"))
  expect_identical(.fallback_decision(.9, 500L)$report, .fallback_switched("fallback_low_acceptance"))
  # a threshold of 0 never switches
  expect_identical(.fallback_decision(0, 0L)$report[["fallback_few_proposals"]], 0)

  # (c) no surrogate that differed from the exact density: a vector rule that is
  # a best-weight rule in every state of an indicator (the surrogate falls back to
  # the exact density), or a kernel mode that is the normal kernel in every state
  indicator <- "\n  ind ~ dcat(c(.5, .5))"
  expect_identical(.fallback_decision(.9, 500L, vector_rule = "vr", definitions = paste0(indicator, "\n  vr <- ind"))$report,
                   .fallback_switched("fallback_surrogate_exact"))
  expect_identical(.fallback_decision(.9, 500L, mode = "md", definitions = paste0(indicator, "\n  md <- 0 * ind"))$report,
                   .fallback_switched("fallback_surrogate_exact"))
  # a dense block under the normal kernel in every state: its anchor normalizer
  # is the exact one, so it falls back by the same rule
  dense <- .coarse_dense_hand_model(.coarse_dense_setup(), min_acceptance = .9, mode = "md",
                                    definitions = paste0(indicator, "\n  md <- 0 * ind"), inits = list(ind = 1))
  expect_identical(dense$samplers[dense$samplers != "base::Finite"], "RoBMA::CoarseCorrectedSlice")
  expect_identical(.coarse_hand_chain(dense, "mu", 200L)$report, .fallback_switched("fallback_surrogate_exact"))

  # the order of the rules: (a) before (c) and before (b); (c) and (b) cannot
  # both hold, since a surrogate equal to the exact density accepts every proposal
  expect_identical(.fallback_decision(.9, 49L, vector_rule = "vr", definitions = paste0(indicator, "\n  vr <- ind"))$report,
                   .fallback_switched("fallback_few_proposals"))
  expect_identical(.fallback_decision(.9, 49L)$report, .fallback_switched("fallback_few_proposals"))

  # a second adaptOff() changes nothing: the sampler that switched is counted
  # once and keeps the exact kernel
  setup <- .coarse_exactness_setup("cluster")
  hand  <- .coarse_hand_model("cluster", 4L, 1L, setup$rules,
    data = list(residual_sd = setup$residual_sd, loading = as.numeric(setup$loading)),
    definitions = "  mu ~ dt(0, 1, 4)\n  for (j in 1:4) { means[j] <- mu }", inits = list(mu = -.25),
    min_acceptance = .9)
  for (batch in seq_len(5L)) rjags::adapt(hand$model, 100L, end.adaptation = FALSE)
  expect_true(rjags::adapt(hand$model, 0L, end.adaptation = TRUE))
  expect_identical(.Call("RoBMA_selnorm_sampler_report", FALSE, PACKAGE = "RoBMA"),
                   .fallback_switched("fallback_low_acceptance"))
  .Call("adapt_off", hand$model$ptr(), PACKAGE = "rjags")
  stats::update(hand$model, n.iter = 100L, progress.bar = "none")
  expect_identical(.selection_sampler_report_read(), .fallback_switched("fallback_low_acceptance"))
})


test_that("the worst-sampler rate is the rate of one sampler over its whole run", {

  skip_if_not_certification("These hand-written models (about 3 s) run in the numerical-kernels case.")
  skip_if_not_installed("rjags")
  setup <- .coarse_exactness_setup("cluster", omega = .01)
  build <- function(definitions = "  mu ~ dt(0, 1, 4)\n  for (j in 1:4) { means[j] <- mu }",
                    data = list(residual_sd = setup$residual_sd), inits = list(mu = -.25), ...) {
    .coarse_hand_model("cluster", 4L, 1L, setup$rules,
      data = c(list(loading = as.numeric(setup$loading), omega = setup$omega), data),
      definitions = definitions, inits = inits, ...)
  }

  # one sampler (a strongly selected block, whose early proposals are no better
  # than its later ones): its rate over burn-in and draws, never a lower one
  report <- .coarse_hand_chain(build(), "mu", 2000L)$report
  expect_identical(report[["samplers"]], 1)
  expect_identical(report[["sampling_proposed"]], 500 + 2000)
  expect_identical(report[["min_sampler_rate"]], report[["sampling_accepted"]] / report[["sampling_proposed"]])
  expect_lt(report[["min_sampler_rate"]], 1)
  # two reads in a row: the second has nothing
  expect_identical(.selection_sampler_report_read(), .fallback_report_value(rate = NA_real_))

  # a sampler destroyed before the read still counts with its rate
  hand <- build()
  rjags::adapt(hand$model, 500L, end.adaptation = TRUE)
  stats::update(hand$model, n.iter = 1000L, progress.bar = "none")
  rm(hand)
  invisible(gc())
  report <- .selection_sampler_report_read()
  expect_identical(report[["sampling_proposed"]], 1000)
  expect_identical(report[["min_sampler_rate"]], report[["sampling_accepted"]] / report[["sampling_proposed"]])

  # fewer than 1,000 proposals in the run: no rate
  hand <- build()
  rjags::adapt(hand$model, 500L, end.adaptation = TRUE)
  stats::update(hand$model, n.iter = 999L, progress.bar = "none")
  report <- .selection_sampler_report_read()
  expect_identical(report[["sampling_proposed"]], 999)
  expect_identical(report[["min_sampler_rate"]], NA_real_)
  # the next run of the same sampler is its own window
  stats::update(hand$model, n.iter = 1000L, progress.bar = "none")
  report <- .selection_sampler_report_read()
  expect_identical(report[["min_sampler_rate"]], report[["sampling_accepted"]] / report[["sampling_proposed"]])

  # two samplers (the location and the residual SD): the lower of their rates,
  # at most the pooled rate
  report <- .coarse_hand_chain(build(
    definitions = "  mu ~ dt(0, 1, 4)\n  sdr ~ dunif(0.3, 1.5)\n  for (j in 1:4) { means[j] <- mu\n    residual_sd[j] <- sdr }",
    data = list(), inits = list(mu = -.25, sdr = .5)), c("mu", "sdr"), 2000L)$report
  expect_identical(report[["samplers"]], 2)
  expect_lte(report[["min_sampler_rate"]], report[["sampling_accepted"]] / report[["sampling_proposed"]])
  expect_gt(report[["min_sampler_rate"]], 0)
})


# ---------------------------------------------------------------------------- #
# Fits (certification case 'numerical-kernels').
# ---------------------------------------------------------------------------- #

# Muffles every warning except those of the corrected sampler.
.fallback_quiet <- function(expr) {

  withCallingHandlers(expr, warning = function(w) {
    if (!inherits(w, "RoBMA_selection_sampler_warning")) invokeRestart("muffleWarning")
  })
}

# A small dependent-row selection fit: three Assink et al. studies give two
# rank-one blocks and a singleton (as in test-00-selection-coarse-blocks.R).
.fallback_fit <- function(parallel = FALSE, chains = 2L, seed = 1L, sample = 100L, ...) {

  data("dat.assink2016", package = "metadat", envir = environment())
  dat <- dat.assink2016[dat.assink2016$study %in% c(1, 2, 4), ]
  dat[["study"]] <- as.integer(factor(dat[["study"]]))
  V <- metafor::vcalc(vi, cluster = study, obs = esid, rho = .6, data = dat)
  .fallback_quiet(RoBMA::bselmodel.mv(
    yi = yi, V = V, random = ~ 1 | study / esid, data = dat, measure = "SMD",
    chains = chains, adapt = 50L, burnin = 50L, sample = sample, seed = seed, parallel = parallel,
    silent = TRUE, convergence_checks = RoBMA::set_convergence_checks(max_Rhat = NULL, min_ESS = NULL),
    ...))
}

.fallback_record <- function(fit) attr(fit[["fit"]], "selection_sampler", exact = TRUE)

.fallback_values <- function(fit) attr(fit[["fit"]], "runtime_report_values", exact = TRUE)

# the counts of runtime report values, summed over runs and processes
.fallback_value_sums <- function(values) {

  processes <- unlist(values, recursive = FALSE)
  counts    <- .fallback_report_names[1:6]  # the integer counts
  sums <- vapply(counts, function(name) sum(vapply(processes, `[[`, numeric(1L), name)), numeric(1L))
  c(sums, fallback_samplers = sum(sums[c("fallback_few_proposals", "fallback_surrogate_exact",
                                          "fallback_low_acceptance")]))
}

# Every post-adaptation update of a sampler that stayed corrected is counted:
# burn-in, retained draws and runjags' two closing draws per backend run.
.fallback_expect_identity <- function(counts, iterations, label) {

  expect_identical(counts[["sampling_proposed"]],
                   iterations * (counts[["samplers"]] - counts[["fallback_samplers"]]), info = label)
}

.fallback_count_names <- c("samplers", "fallback_samplers", "fallback_few_proposals",
  "fallback_surrogate_exact", "fallback_low_acceptance", "sampling_proposed", "sampling_accepted")


test_that("serial and parallel fits record their corrected samplers, and update() adds once", {

  skip_on_cran()
  skip_if_not_certification("These fits (about 40 s) run in the numerical-kernels case.")
  skip_if_not_installed("rjags")
  skip_if_not_installed("metadat")
  skip_if_not_installed("metafor")
  coarse <- "RoBMA::CoarseCorrectedSlice"

  # 1,000 draws: each sampler makes at least 1,000 proposals, so it has a rate
  serial   <- .fallback_fit(parallel = FALSE, chains = 3L, sample = 1000L)
  parallel <- .fallback_fit(parallel = TRUE, chains = 3L, sample = 1000L)
  fits     <- list(serial = serial, parallel = parallel)
  for (label in names(fits)) {
    fit    <- fits[[label]]
    record <- .fallback_record(fit)
    values <- .fallback_values(fit)
    expect_identical(names(record), .fallback_record_names, info = label)
    expect_identical(record[c("sampler", "min_acceptance", "complete")],
                     list(sampler = "slice.coarse", min_acceptance = .9, complete = TRUE), info = label)
    # one run, with one report per process that ran chains, whose sum is the record
    expect_length(values, 1L)
    expect_length(values[[1L]], if (label == "serial") 1L else 3L)
    expect_identical(unlist(record[.fallback_count_names]), .fallback_value_sums(values)[.fallback_count_names],
                     info = label)
    # one sampler object per corrected node and chain, and one per chain for the
    # declared block of this single model (the corrected block sampler)
    expect_true(any(fit[["fit"]][["samplers"]][["Sampler"]] == "RoBMA::CoarseCorrectedBlock"), info = label)
    expect_identical(record[["samplers"]], 3 * (sum(fit[["fit"]][["samplers"]][["Sampler"]] == coarse) + 1),
                     info = label)
    expect_gt(record[["samplers"]], 0)
    .fallback_expect_identity(record, 50 + 1000 + 2, label)
    # the worst-sampler rate is the lowest of the processes' rates
    process_rates <- vapply(values[[1L]], `[[`, numeric(1L), "min_sampler_rate")
    expect_false(anyNA(process_rates), info = label)
    expect_identical(record[["min_sampler_rate"]], min(process_rates), info = label)
    expect_identical(selection_sampler_info(fit = fit), as.data.frame(record, stringsAsFactors = FALSE))
  }
  expect_identical(.fallback_record(serial)[["samplers"]], .fallback_record(parallel)[["samplers"]])

  # update() rebuilds the model, adapts again, and adds its runs to the record
  # of its input exactly once, with the same or another number of workers
  extensions <- list(
    serial   = .fallback_quiet(stats::update(serial, sample_extend = 60L)),
    parallel = .fallback_quiet(stats::update(parallel, sample_extend = 60L, parallel = TRUE, cores = 2L))
  )
  for (label in names(extensions)) {
    input    <- if (label == "serial") serial else parallel
    extended <- extensions[[label]]
    values   <- .fallback_values(extended)
    expect_length(values, 1L)
    expect_length(values[[1L]], if (label == "serial") 1L else 2L)
    added <- .fallback_value_sums(values)
    .fallback_expect_identity(added, 60 + 2, label)
    record <- .fallback_record(extended)
    expect_identical(unlist(record[.fallback_count_names]),
                     unlist(.fallback_record(input)[.fallback_count_names]) + added[.fallback_count_names],
                     info = label)
    expect_true(record[["complete"]], info = label)
  }
})


test_that("the JASP path, alternating fits and leftover counts keep each fit's record its own", {

  skip_on_cran()
  skip_if_not_certification("These fits (about 20 s) run in the numerical-kernels case.")
  skip_if_not_installed("rjags")
  skip_if_not_installed("metadat")
  skip_if_not_installed("metafor")
  old_options <- RoBMA.options()
  on.exit(do.call(RoBMA.options, old_options), add = TRUE)

  # the JASP path draws in several backend runs of one compiled model: adaptation
  # and burn-in with one draw, then five runs of a fifth of the draws, each run
  # with its two closing draws
  jasp   <- .fallback_fit(is_JASP = TRUE)
  values <- .fallback_values(jasp)
  record <- .fallback_record(jasp)
  expect_length(values, 6L)
  expect_identical(unlist(record[.fallback_count_names]), .fallback_value_sums(values)[.fallback_count_names])
  .fallback_expect_identity(record, 50 + (1 + 2) + 5 * (100 / 5 + 2), "JASP")

  # counts that a model left in the process (an earlier fit or a failed attempt)
  # do not enter the record of the next fit
  first <- .fallback_record(.fallback_fit())
  setup <- .coarse_exactness_setup("cluster")
  hand  <- .coarse_hand_model("cluster", 4L, 1L, setup$rules,
    data = list(residual_sd = setup$residual_sd, loading = as.numeric(setup$loading)),
    definitions = "  mu ~ dt(0, 1, 4)\n  for (j in 1:4) { means[j] <- mu }", inits = list(mu = -.25))
  rjags::adapt(hand$model, 100L, end.adaptation = TRUE)
  stats::update(hand$model, n.iter = 300L, progress.bar = "none")
  expect_identical(.fallback_record(.fallback_fit()), first)

  # the report callback reads and clears, as after a failed attempt
  stats::update(hand$model, n.iter = 300L, progress.bar = "none")
  report <- .selection_sampler_report_callback()
  expect_identical(report(list(phase = "report"))[["sampling_proposed"]], 300)
  expect_identical(report(list(phase = "report")), .fallback_report_value(rate = NA_real_))

  # a fit with the ordinary samplers, between two identical fits
  RoBMA.options(selection.sampler = "slice")
  ordinary <- .fallback_fit()
  RoBMA.options(selection.sampler = "slice.coarse")
  expect_identical(.fallback_record(ordinary)[c("sampler", "samplers", "sampling_proposed", "complete")],
                   list(sampler = "slice", samplers = 0, sampling_proposed = 0, complete = TRUE))
  expect_false(any(ordinary[["fit"]][["samplers"]][["Sampler"]] == "RoBMA::CoarseCorrectedSlice"))
  expect_identical(.fallback_record(.fallback_fit()), first)
})


test_that("a failed extension keeps the input record and a fit saved before the record counts from its extension", {

  skip_on_cran()
  skip_if_not_certification("These fits (about 30 s) run in the numerical-kernels case.")
  skip_if_not_installed("rjags")
  skip_if_not_installed("metadat")
  skip_if_not_installed("metafor")
  coarse <- "RoBMA::CoarseCorrectedSlice"
  fit    <- .fallback_fit()

  # a failed extension returns its input without report values
  local({
    testthat::local_mocked_bindings(JAGS_extend = function(fit, ...) {
      attr(fit, "runtime_report_values") <- NULL
      fit
    }, .package = "BayesTools")
    failed <- .fallback_quiet(stats::update(fit, sample_extend = 50L))
    expect_identical(.fallback_record(failed), .fallback_record(fit))
  })

  # an extension that drew new samples but lost every report: incomplete, with
  # the input's counts
  extend <- BayesTools::JAGS_extend
  local({
    testthat::local_mocked_bindings(JAGS_extend = function(fit, ...) {
      extended <- extend(fit, ...)
      attr(extended, "runtime_report_values") <- NULL
      extended
    }, .package = "BayesTools")
    unreported <- .fallback_quiet(stats::update(fit, sample_extend = 50L))
    expect_gt(.fit_draw_count(unreported[["fit"]]), .fit_draw_count(fit[["fit"]]))
    record <- .fallback_record(unreported)
    expect_false(record[["complete"]])
    expect_identical(record[.fallback_count_names], .fallback_record(fit)[.fallback_count_names])
  })

  # A fit saved by RoBMA 4.1.1: its setup closure (the body of 4.1.1) carries
  # settings without the threshold, and it stores no report callback and no
  # record. Its extension keeps its sampler with the fallback off and counts
  # from the extension on.
  legacy_fit <- function(sampler) {
    legacy   <- fit
    settings <- attr(legacy[["fit"]], "selection_runtime", exact = TRUE)
    settings[["coarse_min_acceptance"]] <- NULL
    settings[["sampler"]] <- sampler
    setup <- function(context) {
      utils::getFromNamespace(".selection_runtime_configure", "RoBMA")(
        settings, context
      )
      invisible(NULL)
    }
    environment(setup) <- list2env(list(settings = settings), parent = baseenv())
    attr(setup, "selection_runtime") <- settings
    attr(legacy[["fit"]], "runtime_setup")         <- setup
    attr(legacy[["fit"]], "selection_runtime")     <- settings
    attr(legacy[["fit"]], "runtime_report")        <- NULL
    attr(legacy[["fit"]], "runtime_report_values") <- NULL
    attr(legacy[["fit"]], "selection_sampler")     <- NULL
    legacy
  }
  for (parallel in c(FALSE, TRUE)) {
    label <- if (parallel) "parallel" else "serial"
    corrected <- .fallback_quiet(stats::update(legacy_fit("coarse_corrected"), sample_extend = 60L,
                                                parallel = parallel))
    record    <- .fallback_record(corrected)
    expect_identical(record[c("sampler", "min_acceptance", "fallback_samplers", "complete")],
                     list(sampler = "coarse_corrected", min_acceptance = 0, fallback_samplers = 0, complete = FALSE),
                     info = label)
    expect_true(any(corrected[["fit"]][["samplers"]][["Sampler"]] == coarse), info = label)
    .fallback_expect_identity(record, 60 + 2, label)

    ordinary <- .fallback_quiet(stats::update(legacy_fit("default"), sample_extend = 60L, parallel = parallel))
    expect_identical(.fallback_record(ordinary)[c("sampler", "min_acceptance", "samplers", "complete")],
                     list(sampler = "default", min_acceptance = 0, samplers = 0, complete = FALSE), info = label)
    expect_false(any(ordinary[["fit"]][["samplers"]][["Sampler"]] == coarse), info = label)
  }
})


test_that("a fit with a low acceptance warns and keeps one message through update()", {

  skip_on_cran()
  skip_if_not_certification("These fits (about 10 s) run in the numerical-kernels case.")
  skip_if_not_installed("rjags")
  skip_if_not_installed("metadat")
  skip_if_not_installed("metafor")
  old_options <- RoBMA.options()
  on.exit(do.call(RoBMA.options, old_options), add = TRUE)

  # a serial fit whose report gives a rate of 80% after adaptation
  testthat::local_mocked_bindings(.selection_sampler_report_read = function() {
    .fallback_report_value(samplers = 6, low = 1, proposed = 2000, accepted = 1600, rate = .75)
  })
  stored <- function(fit) {
    warnings <- attr(fit[["fit"]], "warnings", exact = TRUE)
    warnings[startsWith(warnings, .selection_sampler_warning_prefix)]
  }
  fit <- NULL
  expect_warning(fit <- .fallback_fit(), class = "RoBMA_selection_sampler_low_acceptance")
  expect_length(stored(fit), 1L)
  expect_match(stored(fit), "accepted 80% of its proposals after adaptation", fixed = TRUE)
  expect_true(stored(fit) %in% fit[["fit"]][["warnings"]])

  extended <- NULL
  expect_warning(extended <- .fallback_quiet(stats::update(fit, sample_extend = 50L)), class = "RoBMA_selection_sampler_warning")
  expect_length(stored(extended), 1L)
  expect_match(stored(extended), "2 of 12 samplers switched", fixed = TRUE)

  # a threshold of 0 never warns
  RoBMA.options(selection.coarse_min_acceptance = 0)
  quiet <- .fallback_capture_warnings(.fallback_fit())
  expect_false(any(vapply(quiet$warnings, inherits, logical(1L), "RoBMA_selection_sampler_warning")))
  expect_length(stored(quiet$value), 0L)
})
