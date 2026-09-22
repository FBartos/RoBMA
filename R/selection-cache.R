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
#' computation initializes the current process's cache, including computations
#' on a saved fit. Explicit option changes and fitting setup also initialize it.
#' Exact and coarse entries share one lazily allocated pool with global
#' least-recently-used eviction. The ordinary sampler does not use coarse entries.
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
#' the budget or clearing all entries releases the pool. Clearing one component
#' leaves the shared pool available for reuse. Switching sampler options does
#' not clear existing entries or replace already compiled samplers.
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
#' shared capacity and allocated pool/index storage; component byte counts are
#' \code{NA} because pool blocks are shared. They exclude the surrounding R
#' process and JAGS model. The allocated size can remain unchanged after eviction.
#' Allocation-failure counts include attempts refused by the configured cap;
#' they do not necessarily indicate exhaustion of system memory.
#' @export
selection_cache_info <- function(clear = FALSE, component = c("all", "exact", "coarse")) {

  BayesTools::check_bool(clear, "clear", allow_NA = FALSE)
  component <- match.arg(component)
  if (!.selection_runtime_available()) {
    stop("Selection cache information is unavailable because the native RoBMA routines are not loaded.",
         call. = FALSE)
  }
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
#' @description Counts completed coarse slice proposals and their full-target
#' acceptance decisions in the current R process.
#' @param clear reset counters before reporting them.
#' @details Enable the alternative for future fits with
#' \code{RoBMA.options(selection.sampler = "coarse_corrected")}.
#' Every proposal is corrected against the full original likelihood and prior.
#' Only the selection normalizer used to construct proposals is coarsened.
#' A reversible slice transition for that surrogate is followed by a
#' Metropolis correction using both original and surrogate densities at the
#' current and proposed states. The original likelihood keeps its existing
#' numerical integration controls and diagnostics; the correction does not
#' replace those numerical integrals by symbolic exact values.
#'
#' The surrogate rounds means and the logarithms of positive weights to fixed
#' grids. Each covariance diagonal is rounded upward on a grid whose origin
#' is that row's squared original sampling SE; off-diagonal covariances,
#' original SEs, selection cutoffs, and zero weights are retained. These are
#' proposal grids, not changes to fitted variance components or the input
#' covariance. The Gaussian numerator and observed weights use the actual
#' parameters.
#'
#' At most the first \code{selection.coarse_max_rules} rules of the fitted
#' quadrature schedule are tried, with its existing integration criterion.
#' If a valid anchor cannot be constructed or its covariance-envelope calculation
#' is unsupported or unsuccessful, its prescribed normalizer is one; the
#' proposal stage does not fall back to QMC. This fallback is independent of
#' cache history and memory capacity and remains subject to the full-target
#' correction. It does not replace full-target integration failures.
#'
#' Already compiled samplers retain their grid and coarse rule-budget settings.
#' Discrete, vector-valued, and unrelated updates retain ordinary JAGS samplers.
#' Coarse settings can affect sampling efficiency even when correction
#' acceptance is high. Compare effective samples per second for parameters of
#' interest, not only iteration times or this acceptance rate.
#' Counters combine all corrected samplers in this process, including warmup;
#' clear them between phases for a sampling-only comparison. They do not report
#' separate parallel workers. Inspect them between sampling calls.
#' @return A data frame with proposal, acceptance and correction-rejection
#' counts, and an acceptance rate (\code{NA} before any proposals).
#' @export
selection_sampler_info <- function(clear = FALSE) {

  BayesTools::check_bool(clear, "clear", allow_NA = FALSE)
  if (!.selection_runtime_available()) {
    stop("Selection sampler information is unavailable because the native RoBMA routines are not loaded.",
         call. = FALSE)
  }
  info <- .Call("RoBMA_selnorm_sampler_control", NULL, NULL, clear, PACKAGE = "RoBMA")
  return(data.frame(component = "coarse_corrected", proposed = info[["proposed"]],
    accepted = info[["accepted"]], correction_rejected = info[["correction_rejected"]],
    acceptance_rate = if (info[["proposed"]] > 0) info[["accepted"]] / info[["proposed"]] else NA_real_))
}

.selection_cache_configure <- function(capacity_bytes) {

  invisible(.Call("RoBMA_selnorm_cache_control", as.numeric(capacity_bytes), 0L,
                  PACKAGE = "RoBMA"))
}

.selection_runtime_available <- function() {

  all(vapply(c("RoBMA_selnorm_cache_control",
    "RoBMA_selnorm_sampler_control"), is.loaded, logical(1L), PACKAGE = "RoBMA"))
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
    coarse_max_rules = options[["selection.coarse_max_rules"]])
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
  .Call("RoBMA_selnorm_sampler_control", settings[["sampler"]] == "coarse_corrected",
        as.numeric(c(grid, settings[["coarse_max_rules"]])), FALSE, PACKAGE = "RoBMA")
  RoBMA.private[["selection_runtime_initialized"]] <- TRUE
  invisible(NULL)
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

  list(packages = BayesTools:::.JAGS_package_builds(c("RoBMA", "BayesTools")),
    R = R.version.string, platform = R.version[["platform"]],
    JAGS = as.character(rjags::jags.version()))
}

.selection_cache_restore <- function(state) {

  # Begin with a fit-owned pool even when this is the first retained fit in a
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
