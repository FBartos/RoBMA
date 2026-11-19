#' Inspect or clear the native selection cache
#'
#' @description Reports the bounded native selection-normalizer cache in the
#' current R process. Exact caching reuses identical numerical inputs and does
#' not change the likelihood or its integration diagnostics.
#' Post-fit requests using bounded CDF interpolation bypass these caches,
#' including when they fall back to the ordinary primitive. Fitting retains
#' its ordinary-primitive cache.
#'
#' @param clear clear the requested entries and reset their counters before
#' reporting them. Clearing \code{"all"} also releases the shared storage.
#' @param component report the \code{"exact"} cache, the \code{"coarse"}
#' proposal cache, or \code{"all"} caches.
#'
#' @details The default \code{selection.cache_max_bytes = "auto"} chooses the
#' smaller of 4 GiB and one quarter of currently available system RAM. Available
#' RAM is reported by \code{ps::ps_system_memory()}, excluding swap. A numeric
#' setting supplies an explicit total byte budget; \code{0} disables storage.
#' If the available-memory query fails, a warning is issued and caching stays
#' disabled until a budget can be resolved or supplied explicitly.
#' Loading the package does not query available memory. The first selection
#' computation or cache inspection initializes the current process's cache, including computations
#' on a saved fit. Explicit option changes and fitting setup also initialize it.
#' Exact and coarse entries share one byte budget for directly allocated key,
#' control, list-node and index blocks, with global least-recently-used eviction.
#' Storage is allocated lazily. The ordinary sampler does not use coarse entries.
#'
#' The budget is resolved once before fitting and divided equally among chains.
#' Shares are pooled for chains assigned to the same process: ten parallel
#' workers receive one tenth each, while ten sequential chains share the total.
#' Integer-byte rounding can leave a few bytes unused. The coordinator releases
#' its cache while workers run and restores its local budget when they stop.
#' Separate simultaneous fits are separate budgets; this is not a machine-wide
#' memory manager. Stored fit extensions retain the resolved budget and sampler
#' settings, redistributing the budget for their current worker topology.
#'
#' This function reports only the current process, not separate workers or saved
#' fitted-object caches. By default entries live in the native process. Set
#' \code{RoBMA.options(selection.cache_retain = TRUE)} to retain them with the
#' fit; see \code{\link{remove_selection_cache}}. Changing
#' the budget or clearing all entries releases the shared storage. Clearing one
#' component releases its keys, nodes and controls; shared index buckets can remain
#' allocated. Switching sampler options does not clear existing entries or replace
#' already compiled samplers.
#'
#' The exact cache accelerates successful deterministic covariance-envelope
#' normalizers for multivariate step selection. Unsupported numerical routes
#' continue through their ordinary evaluator. The coarse cache can also store its
#' prescribed unit-normalizer fallback. If cache storage cannot be allocated,
#' the same exact or surrogate evaluator runs without storage; memory capacity
#' never determines the returned normalizer or its numerical diagnostics.
#'
#' @return A data frame with the cache component, byte counts, entry count, and
#' lookup, eviction and allocation counters. The \code{"total"} row reports the
#' shared capacity and directly allocated blocks, including retained controls and
#' index buckets; component byte counts are \code{NA} because storage is shared.
#' Allocator metadata, the surrounding R process and JAGS model are excluded;
#' this is not a process-RSS limit. Eviction releases entry blocks, while controls
#' and index buckets can remain allocated. Oversize lookup keys increment misses
#' without incrementing allocation failures. Failed insertion, capacity or
#' allocation attempts can increment allocation failures; these counts do not
#' necessarily indicate exhaustion of system memory.
#' @export
selection_cache_info <- function(clear = FALSE, component = c("all", "exact", "coarse")) {

  BayesTools::check_bool(clear, "clear", allow_NA = FALSE)
  component <- match.arg(component)
  if (!.selection_runtime_available()) {
    stop("Selection cache information is unavailable because the native RoBMA routines are not loaded.",
         call. = FALSE)
  }
  .selection_runtime_ensure()
  mask <- if (clear) switch(component, exact = 1L, coarse = 2L, all = 3L) else 0L
  info <- .Call("RoBMA_selnorm_cache_control", NULL, mask, PACKAGE = "RoBMA")
  components <- if (component == "all") c("total", "exact", "coarse") else component
  memory <- c("capacity_bytes", "allocated_bytes", "peak_bytes")
  counts <- names(info[["exact"]])
  rows <- lapply(components, function(part) {

    bytes <- if (part == "total") info[memory] else
      stats::setNames(rep(list(NA_real_), length(memory)), memory)
    counters <- if (part == "total") {
      stats::setNames(lapply(counts, function(name) {
        info[["exact"]][[name]] + info[["coarse"]][[name]]
      }), counts)
    } else info[[part]]
    data.frame(component = part, as.data.frame(c(bytes, counters)),
      row.names = NULL, check.names = FALSE)
  })
  return(do.call(rbind, rows))
}

#' Remove a fitted model's selection cache
#'
#' @description Removes retained selection-normalizer cache entries from a fitted
#' object without changing its draws, inference, or ability to obtain more draws.
#' @param object a fitted \code{brma} object, including \code{bselmodel.mv()} and
#' \code{RoBMA()} fits.
#' @details Set \code{RoBMA.options(selection.cache_retain = TRUE)} before fitting
#' or extending to retain exact and coarse cache entries. The default is
#' \code{FALSE}. This option is checked for each RoBMA fit or \code{update()} call. Existing
#' retained entries are used during extension even when retention is turned off;
#' in that case the returned fit does not retain a new cache.
#'
#' Retained entries are ordinary serialized data and survive \code{saveRDS()}.
#' They are restored only with matching package code, native builds, R and JAGS.
#' Incompatible or damaged cache data are discarded with a warning; extension
#' proceeds with ordinary likelihood evaluation. Cache contents never replace
#' the model, numerical integration controls, or JAGS random-number state.
#' Reconstructed JAGS models can still require adaptation before collecting new
#' samples; retaining the normalizer cache does not retain sampler tuning.
#'
#' Each worker saves its own entries. Extension with the same worker count
#' restores the corresponding worker cache. With a different worker count,
#' previous worker caches are assigned round-robin without duplication, and
#' entries that do not fit the destination's share of the memory budget are
#' dropped. Increasing the worker count can leave some workers with an empty
#' cache initially. Cache reuse can improve speed but is not guaranteed.
#'
#' Successful capture releases the native cache. A retained object's compact
#' cache payload is bounded by the fitted total native budget, apart from small
#' R metadata overhead. The byte budget is not a limit on total R memory:
#' capture, serialization and extension can temporarily hold both a snapshot
#' and native storage, and other copies or saved fits can retain their own
#' snapshots. Save or load operations can therefore require substantial memory.
#'
#' Assign the result back to the object to release its reference to the retained
#' data. Memory remains in use if another object references the same snapshot.
#' This helper does not clear the current process's independent native cache;
#' use \code{selection_cache_info(clear = TRUE)} for that. Future retention is
#' controlled by \code{selection.cache_retain}.
#' @return The fitted object with its retained selection cache removed.
#' @examples \dontrun{
#' RoBMA.options(selection.cache_retain = TRUE)
#' fit <- bselmodel.mv(yi, V = V, random = ~ 1 | study / esid, data = dat)
#' saveRDS(fit, "fit.rds")
#' fit <- update(readRDS("fit.rds"), sample_extend = 1000)
#' fit <- remove_selection_cache(fit)
#' RoBMA.options(selection.cache_retain = FALSE)
#' }
#' @export
remove_selection_cache <- function(object) {

  if (!inherits(object, "brma")) {
    stop("'object' must be a 'brma' object.", call. = FALSE)
  }
  attr(object[["fit"]], "runtime_state") <- NULL
  return(object)
}

#' Inspect corrected selection-sampler counters
#'
#' @description Reports the work of the corrected selection sampler: the record
#' of a fitted model, or the counters of the current R process.
#' @param clear reset the process counters before reporting them. Cannot be
#' combined with \code{fit}.
#' @param fit a fitted selection model (for example from \code{bselmodel.mv()}
#' or \code{RoBMA.mv()}) whose record to report, or \code{NULL} for the
#' process counters.
#' @details The corrected sampler is the default,
#' \code{RoBMA.options(selection.sampler = "slice.coarse")}; set
#' \code{selection.sampler = "slice"} before fitting to use the other samplers
#' instead, or \code{"jags.default"} for only the samplers JAGS assigns.
#' Parameters that the options \code{sampler.blocks} and
#' \code{sampler.gates} declare to the factor slice and gated samplers of the
#' BayesTools JAGS module (see \code{\link{RoBMA.options}}) keep those samplers
#' under \code{"slice.coarse"}, so set the two options to \code{FALSE} for the
#' corrected sampler to update them as well.
#' Every proposal is corrected against the full original likelihood and prior.
#' Only the selection normalizer used to construct proposals is coarsened.
#' A reversible slice transition for that surrogate is followed by a
#' Metropolis correction using both original and surrogate densities at the
#' current and proposed states. The original likelihood keeps its existing
#' numerical integration controls and diagnostics; the correction does not
#' replace those numerical integrals by symbolic exact values.
#'
#' Dense multivariate step blocks. The surrogate rounds means and the
#' logarithms of positive weights to fixed grids. Each covariance diagonal is
#' rounded upward on a grid whose origin is that row's squared original
#' sampling SE; off-diagonal covariances, original SEs, selection cutoffs, and
#' zero weights are retained. These are proposal grids, not changes to fitted
#' variance components or the input covariance. The Gaussian numerator and
#' observed weights use the actual parameters. At most the first
#' \code{selection.coarse_max_rules} rules of the fitted quadrature schedule
#' are tried, with its existing integration criterion. If a valid anchor cannot
#' be constructed or its covariance-envelope calculation is unsupported or
#' unsuccessful, its prescribed normalizer is one; the proposal stage does not
#' fall back to QMC. This fallback is independent of cache history and memory
#' capacity and remains subject to the full-target correction. It does not
#' replace full-target integration failures.
#'
#' Rank-one cluster and factor step blocks. The surrogate is the exact
#' numerator (Gaussian log density and observed log weights at the actual
#' parameters) minus the log normalizer of one rule of the block's own
#' Gauss-Hermite schedule: rule number \code{selection.coarse_max_rules}, or the
#' last rule of a shorter schedule, evaluated at the actual parameters. No grid
#' rounding, change test, certificate, refinement, sparse grid, QMC or cache
#' enters, so the surrogate is a deterministic function of the sampled state;
#' the correction evaluates the exact block, which keeps its ordinary numerical
#' controls and diagnostics. Where that rule is unavailable at an update (a
#' best-weight rule in the current state, or a non-finite normalizer) the
#' surrogate is the exact density of the block. A block takes part unless its
#' fixed controls rule it out: a weight rule fixed at a best-weight rule, a
#' kernel mode fixed at the normal kernel, or, for a factor block, a rule beyond
#' the node budget at the declared rank (rank two to six at the default
#' setting). In model-averaged fits the weight rule and the kernel mode are
#' nodes of the bias indicator, so their blocks take part, with the exact
#' density in the states of the indicator that have a best-weight rule. Blocks
#' fixed at a best-weight rule, blocks of higher rank and blocks with a fixed
#' normal kernel keep the ordinary sampler.
#'
#' A scalar update is eligible when at least one stochastic child of the updated
#' node is a dense multivariate block or a coarse-capable cluster or factor
#' block; its other children enter the surrogate through their exact densities.
#' Discrete, vector-valued and unrelated updates, and updates whose selection
#' children are all without a coarse rule, retain ordinary JAGS samplers.
#' Whatever the surrogate, the corrected chain targets the exact full
#' conditional; a coarser rule changes the acceptance rate and the cost, not
#' the posterior.
#'
#' Fallback. Each corrected sampler (one per parameter and chain) chooses its
#' kernel once, at the end of adaptation, from its adaptation updates: with
#' \code{selection.coarse_min_acceptance} above 0 (default \code{0.90}) it
#' updates its parameter with the exact slice sampler for the rest of the fit
#' when it made fewer than 50 proposals during adaptation, when no surrogate it
#' evaluated during adaptation differed from the exact density (no block's
#' surrogate differs under the normal kernel, and a cluster or factor block's
#' does not where its rule is unavailable, as in a model average whose
#' selection branches all use the best-weight rule), or when its acceptance
#' rate during adaptation is below the threshold; the first of these that
#' holds is the recorded reason.
#' Both kernels leave the posterior unchanged and the choice is not revisited,
#' so the chain after adaptation is a single exact kernel. A threshold of
#' \code{0} disables the fallback and the warning. Every model build (a fit,
#' and each parallel worker's or extension's rebuild of the model) adapts and
#' decides again.
#'
#' Record. A fitted selection model keeps the record of its corrected samplers,
#' collected from every process that ran chains and summed over the fit and
#' its extensions (\code{update()}): \code{sampler} and \code{min_acceptance},
#' the settings the model was fitted with (\code{sampler} shows the value the
#' fit was made with, so fits of earlier versions show \code{"coarse_corrected"}
#' or \code{"default"}; \code{min_acceptance} is \code{0} for fits made before
#' the threshold existed); \code{samplers}, the corrected sampler objects
#' constructed, one per parameter and chain and counted again for every model
#' build; \code{fallback_samplers}, those that switched to the exact slice
#' sampler, split by reason into \code{fallback_few_proposals},
#' \code{fallback_surrogate_exact} and \code{fallback_low_acceptance};
#' \code{sampling_proposed} and \code{sampling_accepted}, the proposals made
#' after adaptation by samplers that stayed corrected and the proposals whose
#' correction accepted them (unchanged values included), counting burn-in,
#' retained and extension draws and the closing draws the sampler makes after
#' recording the chain end; \code{acceptance_rate}, their ratio
#' (\code{NA} without proposals); \code{min_sampler_rate}, the lowest
#' acceptance rate of one sampler over one backend run (of the fit, of an
#' automatic extension or of an \code{update()}), among samplers with at least
#' 1,000 proposals in that run (\code{NA} when none has). A fit made of several
#' backend runs (automatic or \code{update()} extensions, or fits that sample in
#' several runs, as JASP does) takes this minimum over more and shorter runs;
#' while a sampler still makes at least 1,000 proposals in a run, such a fit
#' reaches the warning more readily than one long run of the same model, and
#' with fewer proposals per run it has no worst-sampler rate. The
#' last element is \code{complete}, \code{FALSE} when part of the chain was drawn
#' without a record (an extension of a fit saved before the record existed, a
#' process whose report failed, or draws added without any report), in which
#' case the counts cover the recorded part only.
#' Adaptation updates and updates of samplers that fell back are not counted.
#'
#' Warning. A fit or an extension whose record has at least 1,000 counted
#' proposals and whose acceptance rate, or \code{min_sampler_rate}, is below
#' \code{min_acceptance} warns with the condition classes
#' \code{RoBMA_selection_sampler_low_acceptance} and
#' \code{RoBMA_selection_sampler_warning}: sampling is exact but inefficient.
#' The message is kept with the fit's warnings, replaced at every extension.
#'
#' Process counters. Without \code{fit}, the function reports the counters of
#' all corrected samplers in this process, including adaptation; clear them
#' between phases for a sampling-only comparison. They do not report separate
#' parallel workers. Inspect them between sampling calls.
#'
#' Already compiled samplers retain their grid, coarse rule-budget and
#' threshold settings. Coarse settings can affect sampling efficiency even when
#' correction acceptance is high. Compare effective samples per second for
#' parameters of interest, not only iteration times or this acceptance rate.
#' @return Without \code{fit}, a data frame with proposal, acceptance and
#' correction-rejection counts, and an acceptance rate (\code{NA} before any
#' proposals). With \code{fit}, a one-row data frame with the record described
#' in details; a fit without a record (not a selection model, or fitted before
#' the record existed) gives \code{NA} counts and \code{complete = FALSE}.
#' @export
selection_sampler_info <- function(clear = FALSE, fit = NULL) {

  BayesTools::check_bool(clear, "clear", allow_NA = FALSE)
  if (!is.null(fit)) {
    if (clear) {
      stop("'clear' cannot be used with 'fit': a fit's record is not cleared.", call. = FALSE)
    }
    if (!inherits(fit, "brma")) {
      stop("'fit' must be a fitted 'brma' object or NULL.", call. = FALSE)
    }
    return(.selection_sampler_record_frame(fit))
  }
  if (!.selection_runtime_available()) {
    stop("Selection sampler information is unavailable because the native RoBMA routines are not loaded.",
         call. = FALSE)
  }
  info <- .Call("RoBMA_selnorm_sampler_control", NULL, NULL, clear, PACKAGE = "RoBMA")
  return(data.frame(component = "coarse_corrected", proposed = info[["proposed"]],
    accepted = info[["accepted"]], correction_rejected = info[["correction_rejected"]],
    acceptance_rate = if (info[["proposed"]] > 0) info[["accepted"]] / info[["proposed"]] else NA_real_))
}

# The counts of the per-fit record, in the order of the record.
.selection_sampler_count_names <- c("samplers", "fallback_samplers", "fallback_few_proposals",
  "fallback_surrogate_exact", "fallback_low_acceptance", "sampling_proposed", "sampling_accepted")

# The record of a selection fit from the runtime report values of this call
# (a list of retained runs, each a list with one element per process) and,
# for an extension, the record of the input fit. 'new_draws' says whether this
# call added draws to the fit (a fresh fit always does; a failed extension
# returns its input and adds none).
.selection_sampler_record <- function(values, settings, previous = NULL, extension = FALSE,
                                      new_draws = TRUE) {

  reasons <- c("fallback_few_proposals", "fallback_surrogate_exact", "fallback_low_acceptance")
  counts  <- stats::setNames(numeric(length(.selection_sampler_count_names)),
                             .selection_sampler_count_names)
  rates   <- numeric()
  # Draws that no report covers make the record incomplete: new draws without
  # any report, or a process whose report is missing. A call without new draws
  # (a failed extension) changes nothing.
  complete <- !new_draws || length(values) > 0L
  for (run in values) {
    for (process in run) {
      if (is.null(process)) {
        complete <- FALSE
        next
      }
      process[["fallback_samplers"]] <- sum(unlist(process[reasons]))
      for (name in .selection_sampler_count_names) {
        counts[[name]] <- counts[[name]] + process[[name]]
      }
      rates <- c(rates, process[["min_sampler_rate"]])
    }
  }
  if (extension) {
    if (is.null(previous)) {
      complete <- FALSE
    } else {
      counts   <- counts + unlist(previous[.selection_sampler_count_names])
      rates    <- c(rates, previous[["min_sampler_rate"]])
      complete <- complete && isTRUE(previous[["complete"]])
    }
  }
  rates <- rates[!is.na(rates)]

  return(c(
    list(sampler = settings[["sampler"]],
         min_acceptance = .selection_runtime_min_acceptance(settings)),
    as.list(counts),
    list(acceptance_rate  = if (counts[["sampling_proposed"]] > 0)
           counts[["sampling_accepted"]] / counts[["sampling_proposed"]] else NA_real_,
         min_sampler_rate = if (length(rates)) min(rates) else NA_real_,
         complete         = complete)
  ))
}

.selection_sampler_record_frame <- function(object) {

  record <- attr(object[["fit"]], "selection_sampler", exact = TRUE)
  if (is.null(record)) {
    settings <- attr(object[["fit"]], "selection_runtime", exact = TRUE)
    record <- c(
      list(sampler = if (is.null(settings[["sampler"]])) NA_character_ else settings[["sampler"]],
           min_acceptance = if (is.null(settings)) NA_real_ else .selection_runtime_min_acceptance(settings)),
      as.list(stats::setNames(rep(NA_real_, length(.selection_sampler_count_names)),
                              .selection_sampler_count_names)),
      list(acceptance_rate = NA_real_, min_sampler_rate = NA_real_, complete = FALSE)
    )
  }
  return(as.data.frame(record, stringsAsFactors = FALSE))
}

# The beginning of the low-acceptance message, which identifies it among the
# fit's stored warnings.
.selection_sampler_warning_prefix <- "The corrected selection sampler accepted "

# Judge a fit's cumulative record: with at least 1,000 counted proposals, a
# pooled or worst-sampler acceptance rate below the threshold warns once and
# replaces an earlier message of this warning among the fit's warnings, which
# extensions carry forward. Returns the fit's warnings.
.selection_sampler_warn <- function(fit) {

  record   <- attr(fit, "selection_sampler", exact = TRUE)
  warnings <- as.character(attr(fit, "warnings", exact = TRUE))
  warnings <- warnings[!startsWith(warnings, .selection_sampler_warning_prefix)]
  threshold <- record[["min_acceptance"]]
  low <- threshold > 0 && record[["sampling_proposed"]] >= 1000 && (
    record[["acceptance_rate"]] < threshold ||
      (!is.na(record[["min_sampler_rate"]]) && record[["min_sampler_rate"]] < threshold))
  if (low) {
    message <- paste0(
      .selection_sampler_warning_prefix,
      format(100 * record[["acceptance_rate"]], digits = 3), "% of its proposals after adaptation",
      if (!is.na(record[["min_sampler_rate"]])) paste0(
        " (lowest rate of a single sampler ", format(100 * record[["min_sampler_rate"]], digits = 3), "%)"),
      ", below 'selection.coarse_min_acceptance' = ", threshold, "; ",
      record[["fallback_samplers"]], " of ", record[["samplers"]],
      " samplers switched to the exact slice sampler at the end of adaptation. ",
      "Sampling is exact but inefficient. Set 'RoBMA.options(selection.sampler = \"slice\")' ",
      "before fitting to avoid the corrected sampler; an existing fit keeps the settings it was fitted with."
    )
    warning(structure(
      list(message = message, call = NULL,
           acceptance_rate = record[["acceptance_rate"]],
           min_sampler_rate = record[["min_sampler_rate"]],
           threshold = threshold),
      class = c("RoBMA_selection_sampler_low_acceptance", "RoBMA_selection_sampler_warning",
                "simpleWarning", "warning", "condition")
    ))
    warnings <- c(warnings, message)
  }
  if (length(warnings) == 0L) {
    return(NULL)
  }
  return(warnings)
}

.selection_cache_configure <- function(capacity_bytes) {

  invisible(.Call("RoBMA_selnorm_cache_control", as.numeric(capacity_bytes), 0L,
                  PACKAGE = "RoBMA"))
}

.selection_runtime_available <- function() {

  all(vapply(c("RoBMA_selnorm_cache_control",
    "RoBMA_selnorm_sampler_control", "RoBMA_selnorm_sampler_report"),
    is.loaded, logical(1L), PACKAGE = "RoBMA"))
}

.selection_runtime_ensure <- function() {

  if (isTRUE(RoBMA.private[["selection_runtime_initialized"]])) {
    return(invisible(TRUE))
  }
  if (!.selection_runtime_available()) return(invisible(FALSE))
  .selection_runtime_configure(.selection_runtime_settings())
  invisible(TRUE)
}

.selection_runtime_settings <- function(options = .RoBMA_current_options(),
                                         capacity_bytes = NULL) {

  requested <- options[["selection.cache_max_bytes"]]
  available <- NA_real_
  capacity <- if (is.null(capacity_bytes)) requested else capacity_bytes
  if (is.null(capacity_bytes) && identical(requested, "auto")) {
    available <- tryCatch(ps::ps_system_memory()[["avail"]], error = function(e) NA_real_)
    if (!is.numeric(available) || length(available) != 1L ||
        !is.finite(available) || available < 0) {
      warning("Automatic selection cache memory budget is unavailable; caching is disabled. Set 'selection.cache_max_bytes' to a numeric byte budget.",
        call. = FALSE)
      available <- NA_real_
      capacity <- 0
    } else {
      capacity <- floor(min(4 * 1024^3, available / 4))
    }
  }
  list(cache_max_bytes = as.numeric(capacity),
    cache_budget_request = requested, available_memory_bytes = as.numeric(available),
    sampler = options[["selection.sampler"]],
    coarse_grid = options[["selection.coarse_grid"]],
    coarse_max_rules = options[["selection.coarse_max_rules"]],
    coarse_min_acceptance = options[["selection.coarse_min_acceptance"]])
}

# The fallback threshold of frozen settings; settings frozen before it existed
# (fits of RoBMA 4.1.1) mean 0, the corrected sampler without fallback.
.selection_runtime_min_acceptance <- function(settings) {

  threshold <- settings[["coarse_min_acceptance"]]
  if (is.null(threshold)) 0 else threshold
}

.selection_runtime_configure <- function(settings, context = NULL) {

  if (!.selection_runtime_available()) {
    stop("Selection runtime configuration is unavailable because the native RoBMA routines are not loaded.",
         call. = FALSE)
  }
  capacity <- settings[["cache_max_bytes"]]
  if (!is.null(context) && context[["phase"]] == "start") {
    if (context[["role"]] == "coordinator") {
      capacity <- 0
    } else if (isTRUE(context[["parallel"]])) {
      capacity <- floor(capacity / context[["chains"]]) * context[["process_chains"]]
    }
  }
  .selection_cache_configure(capacity)
  grid <- settings[["coarse_grid"]][c("mean", "variance", "log_weight")]
  .Call("RoBMA_selnorm_sampler_control", .selection_sampler_corrected(settings[["sampler"]]),
        as.numeric(c(grid, settings[["coarse_max_rules"]],
                     .selection_runtime_min_acceptance(settings))),
        FALSE, PACKAGE = "RoBMA")
  # A fit starts its sampler record from zero: counts that earlier fits in this
  # process left behind are discarded (the report of each run clears its own).
  if (!is.null(context) && context[["phase"]] == "start") {
    .selection_sampler_report_read()
  }
  RoBMA.private[["selection_runtime_initialized"]] <- TRUE
  invisible(NULL)
}

# The per-fit counters of the corrected sampler in this process, read and
# cleared in one native call.
.selection_sampler_report_read <- function() {

  .Call("RoBMA_selnorm_sampler_report", TRUE, PACKAGE = "RoBMA")
}

# The runtime report of selection fits (BayesTools::JAGS_fit(runtime_report)):
# called after every backend run in each process that ran chains, it returns
# that process's counters of the run. It only reads native counters; it does
# not configure the cache or the sampler factory.
.selection_sampler_report_callback <- function() {

  callback <- function(context) {
    utils::getFromNamespace(".selection_sampler_report_read", "RoBMA")()
  }
  environment(callback) <- baseenv()
  return(callback)
}

.selection_runtime_setup <- function() {

  # A small base-environment closure transports only frozen settings to workers.
  settings <- .selection_runtime_settings()
  # The closure keeps only the frozen settings, so the worker has to resolve the
  # configuration function from its own loaded namespace rather than carry one.
  callback <- function(context) {
    utils::getFromNamespace(".selection_runtime_configure", "RoBMA")(
      settings, context
    )
    invisible(NULL)
  }
  environment(callback) <- list2env(list(settings = settings), parent = baseenv())
  attr(callback, "selection_runtime") <- settings
  return(callback)
}

.selection_cache_build <- function() {

  list(packages = BayesTools::JAGS_package_builds(c("RoBMA", "BayesTools")),
    R = R.version.string, platform = R.version[["platform"]],
    JAGS = as.character(rjags::jags.version()))
}

.selection_cache_restore <- function(state) {

  # Begin with fit-owned storage even when this is the first retained fit in a
  # session. No cache pointer is held by the fitted object or its callbacks.
  selection_cache_info(clear = TRUE)
  state <- Filter(Negate(is.null), state)
  # Resolve the running build rather than trusting a serialized callback from
  # an earlier installation, including direct BayesTools::JAGS_extend() calls.
  build <- if (length(state)) .selection_cache_build() else NULL
  compatible <- vapply(state, function(shard) {

    is.list(shard) && identical(shard[["build"]], build) &&
      is.raw(shard[["payload"]])
  }, logical(1L))
  if (any(!compatible)) {
    warning("Retained selection cache is unavailable with mismatched or missing build metadata. Extending with an empty cache for those entries.",
      call. = FALSE)
  }
  snapshots <- lapply(state[compatible], `[[`, "payload")
  tryCatch(
    .Call("RoBMA_selnorm_cache_restore", snapshots, PACKAGE = "RoBMA"),
    error = function(error) {

      selection_cache_info(clear = TRUE)
      warning("Retained selection cache is unavailable: ", conditionMessage(error),
        " Extending with an empty cache.", call. = FALSE)
      NULL
    }
  )
}

.selection_cache_runtime <- function(fit = NULL) {

  retain <- RoBMA.get_option("selection.cache_retain")
  if (!retain && is.null(attr(fit, "runtime_state", exact = TRUE))) {
    return(NULL)
  }
  callback <- function(context, state = NULL) {

    if (context[["phase"]] == "restore") {
      utils::getFromNamespace(".selection_cache_restore", "RoBMA")(state)
    } else if (context[["phase"]] == "capture" && retain) {
      payload <- .Call("RoBMA_selnorm_cache_snapshot", PACKAGE = "RoBMA")
      build <- if (length(payload)) {
        utils::getFromNamespace(".selection_cache_build", "RoBMA")()
      } else NULL
      RoBMA::selection_cache_info(clear = TRUE)
      if (length(payload) > 0L) return(list(build = build, payload = payload))
    }
    NULL
  }
  environment(callback) <- list2env(list(retain = retain),
    parent = baseenv())
  return(callback)
}
