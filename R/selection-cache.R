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
#' The selection-weight parameters of a branch of a publication-bias mixture
#' that \code{sampler.gates} declares are updated together by the corrected
#' gated block sampler when one of the branch's selection blocks takes part
#' (below); other gated parameters (PET and PEESE slopes, and weights whose
#' selection blocks do not take part) keep the gated sampler of the BayesTools
#' JAGS module. A block that \code{sampler.blocks} declares
#' (see \code{\link{RoBMA.options}}) is updated by the corrected block sampler
#' when one of its selection blocks takes part (below), and by the factor slice
#' sampler of the BayesTools JAGS module otherwise. A scalar indicator (the bias
#' indicator, the heterogeneity inclusion indicator and effect and moderator
#' indicators with fixed prior probabilities) with such a selection block among
#' its children is updated by the screened indicator sampler, unless a joint
#' sampler updates it or \code{sampler.probabilities} records its conditional
#' probabilities. Independent step-selection rows (the estimates of
#' \code{RoBMA()} and the single-estimate clusters of \code{RoBMA.mv()}) take
#' part in the gated weights and the bias indicator through a row cache
#' (below).
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
#'
#' Step-selection rows. The weights of a branch whose estimates enter as
#' independent step-selection rows, and an indicator that sets the weights or
#' the kernel of such rows (the bias indicator), take part through a row
#' cache: while only those weights or only the indicator move, the
#' probabilities of each row's selection intervals are fixed, so they are
#' computed once per update and each surrogate evaluation is a weighted sum
#' and a logarithm per row. This is the rows' exact density up to rounding; a
#' row whose mean moves with the indicator (the PET and PEESE branches), or
#' whose normalizer is too small for the cached sum, enters with its exact
#' density. A proposal is evidence for the cache (below) only when its
#' surrogate and exact changes differ by more than a relative \eqn{10^{-9}}.
#' Indicators that change only the rows' means or variances (heterogeneity,
#' effect and moderator indicators) keep their samplers.
#'
#' Blocks. The corrected block sampler makes one factor slice sweep on the
#' surrogate over all parameters of a declared block, in forward or backward
#' direction order chosen at random at each update, and accepts or rejects the
#' whole sweep with the same correction. The random order makes the proposal
#' reversible, so the corrected chain targets the exact posterior of the block;
#' the directions and widths adapt as in the factor slice sampler. Gamma
#' parameters that enter the model only through their proportions (the
#' Dirichlet weights of a heterogeneity allocation, the weights of a cumulative
#' weight function) are updated through their proportions, and their sum is
#' drawn from its gamma distribution at each update.
#'
#' Gated weights. While its branch is inactive, the weights of a branch are
#' drawn from their prior, as in the gated sampler; while it is active, they
#' are updated by the corrected block sweep, whose directions adapt only on
#' those updates.
#'
#' Indicators. The screened indicator sampler scores every level of the
#' indicator with the surrogates of its selection children, each level's score
#' computed from that level's model state alone. At each update it proposes a
#' level from these scores by one of two moves, chosen at random independently
#' of the state: with probability 3/4 a Metropolized Gibbs step (a level other
#' than the current one, accepted on the scores), with probability 1/4 a draw
#' from the scores with the current level included (which keeps a thinned
#' chain from alternating between two equally likely levels). A change of
#' level is corrected with one Metropolis step against the exact scores of
#' both levels, so the chain targets the exact posterior, and an update that
#' keeps the current level costs no exact evaluation.
#' Vector-valued and unrelated updates, and updates whose selection
#' children are all without a coarse rule (other than the step-selection rows
#' above), retain ordinary JAGS samplers.
#' For any surrogate computed from the model state it scores, the corrected
#' chain targets the exact full conditional; a coarser rule changes the
#' acceptance rate and the cost, not the posterior.
#'
#' Fallback. Each corrected sampler (one per parameter, block or indicator and
#' chain) chooses its kernel once, at the end of adaptation, from its evidence
#' proposals during adaptation: proposals for which at least one selection
#' block evaluated a surrogate that differs from its exact density (none does
#' under the normal kernel of a model-averaged fit's normal branches, or where
#' its rule is unavailable); for an indicator, changes of level whose
#' surrogate differed at either level. With
#' \code{selection.coarse_min_acceptance} above 0 (default \code{0.90}) it
#' updates its parameter with the exact slice sampler (a block with the factor
#' slice sweep of the exact posterior, an indicator by enumerating its levels)
#' for the rest of the fit when it made at least 50 evidence proposals during
#' adaptation and their acceptance rate is below the threshold. With fewer it
#' stays corrected and the warning below judges it after adaptation. A sampler
#' whose selection blocks can never differ from their exact densities (decided
#' from the fixed settings of the blocks) is not created. A block of d
#' directions is judged per direction: its acceptance rate to the power 1/d,
#' here, for \code{min_sampler_rate} and for the warning below.
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
#' constructed, one per parameter, block or indicator and chain and counted
#' again for every model build; \code{fallback_samplers}, those that switched
#' to the exact kernel, and \code{fallback_low_acceptance}, the same count
#' under the one remaining reason (records of development versions also split
#' out two retired reasons, which are dropped);
#' \code{sampling_proposed} and \code{sampling_accepted}, the evidence proposals
#' made after adaptation by samplers that stayed corrected and the proposals
#' whose correction accepted them (unchanged values included), counting burn-in
#' and all sampling iterations through the end of each chain or extension,
#' including iterations between retained draws when thinning;
#' \code{sampling_judged}, the sum over the samplers
#' of their proposals times their acceptance rate per direction (rate^(1/d)
#' for a block of d parameters, the accepted proposals of a scalar sampler,
#' over each backend run); \code{acceptance_rate}, the judged, per-direction
#' rate \code{sampling_judged / sampling_proposed}, which equals the raw rate
#' \code{sampling_accepted / sampling_proposed} for fits without blocks
#' (\code{NA} without proposals); \code{min_sampler_rate}, the lowest
#' acceptance rate of one sampler over one backend run (of the fit, of an
#' automatic extension or of an \code{update()}), among samplers with at least
#' 1,000 proposals in that run (\code{NA} when none has). A fit made of several
#' backend runs (automatic or \code{update()} extensions) takes this minimum
#' over more and shorter runs. While a sampler still makes at least 1,000
#' proposals in a run, such a fit reaches the warning more readily than one
#' long run of the same model, and
#' with fewer proposals per run it has no worst-sampler rate. In JASP, each
#' sampling run contributes one report per process across all its progress ticks.
#' The last element is \code{complete}, \code{FALSE} when part of the chain was drawn
#' without a record (an extension of a fit saved before the record existed, a
#' process whose report failed, or draws added without any report), in which
#' case the counts cover the recorded part only.
#' Adaptation updates, proposals that are not evidence (for example in the
#' normal branches of a model average) and updates of samplers that fell back
#' are not counted.
#'
#' Warning. A fit or an extension whose record has at least 1,000 counted
#' proposals and whose judged \code{acceptance_rate}, or \code{min_sampler_rate}, is below
#' \code{min_acceptance} warns with the condition classes
#' \code{RoBMA_selection_sampler_low_acceptance} and
#' \code{RoBMA_selection_sampler_warning}: sampling is exact but inefficient.
#' The message is kept with the fit's warnings, replaced at every extension.
#'
#' Process counters. Without \code{fit}, the function reports the counters of
#' all corrected samplers in this process (every proposal with a correction,
#' evidence or not, and for an indicator every change of level), including
#' adaptation; clear them
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
.selection_sampler_count_names <- c("samplers", "fallback_samplers", "fallback_low_acceptance",
  "sampling_proposed", "sampling_accepted", "sampling_judged")

# Fallback reasons of development versions (the sampler kept its exact kernel
# after fewer than 50 adaptive proposals, or without a surrogate that differed
# from the exact density); a stored record's elements are dropped.
.selection_sampler_retired_reasons <- c("fallback_few_proposals", "fallback_surrogate_exact")

# A report or record without judged acceptances (made by an earlier build) is
# read as one of scalar samplers, whose judged acceptances are their accepted
# proposals.
.selection_sampler_judged <- function(counts) {

  if (!is.null(counts) && is.null(counts[["sampling_judged"]])) {
    counts[["sampling_judged"]] <- counts[["sampling_accepted"]]
  }
  return(counts)
}

# The record of a selection fit from the runtime report values of this call
# (a list of retained runs, each a list with one element per process) and,
# for an extension, the record of the input fit. 'new_draws' says whether this
# call added draws to the fit (a fresh fit always does; a failed extension
# returns its input and adds none).
.selection_sampler_record <- function(values, settings, previous = NULL, extension = FALSE,
                                      new_draws = TRUE) {

  reasons <- "fallback_low_acceptance"
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
      process <- .selection_sampler_judged(process)
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
      counts   <- counts + unlist(.selection_sampler_judged(previous)[.selection_sampler_count_names])
      rates    <- c(rates, previous[["min_sampler_rate"]])
      complete <- complete && isTRUE(previous[["complete"]])
    }
  }
  rates <- rates[!is.na(rates)]

  return(c(
    list(sampler = settings[["sampler"]],
         min_acceptance = .selection_runtime_min_acceptance(settings)),
    as.list(counts),
    # the judged rate: per direction for blocks, raw for scalar samplers
    list(acceptance_rate  = if (counts[["sampling_proposed"]] > 0)
           counts[["sampling_judged"]] / counts[["sampling_proposed"]] else NA_real_,
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
  # records of development versions carried a flag of the withheld allocation
  # block and two retired fallback reasons, and earlier records no judged
  # acceptances
  record[["allocation_block_withheld"]] <- NULL
  record[.selection_sampler_retired_reasons] <- NULL
  if (is.null(record[["sampling_judged"]])) {
    names  <- names(record)
    at     <- match("sampling_accepted", names)
    record <- c(record[seq_len(at)], list(sampling_judged = record[["sampling_accepted"]]),
                record[-seq_len(at)])
  }
  return(as.data.frame(record, stringsAsFactors = FALSE))
}

# The beginning of the low-acceptance message, which identifies it among the
# fit's stored warnings.
.selection_sampler_warning_prefix <- "The corrected selection sampler accepted "

# The beginning of the fallback message that development versions stored with
# model-averaged fits (the warning is withdrawn); such a stored message is
# dropped when the fit's warnings are refreshed.
.selection_sampler_fallback_pattern <- paste0(
  "^[0-9]+ of [0-9]+ corrected selection samplers switched to the exact slice sampler ",
  "at the end of adaptation \\(")

# Judge a fit's cumulative record: with at least 1,000 counted proposals, a
# pooled judged (per-direction) or worst-sampler acceptance rate below the
# threshold warns once. The message replaces an earlier one among the fit's
# warnings, which extensions carry forward. Returns the fit's warnings.
.selection_sampler_warn <- function(fit) {

  record   <- attr(fit, "selection_sampler", exact = TRUE)
  warnings <- as.character(attr(fit, "warnings", exact = TRUE))
  warnings <- warnings[!startsWith(warnings, .selection_sampler_warning_prefix) &
                         !grepl(.selection_sampler_fallback_pattern, warnings)]
  threshold <- record[["min_acceptance"]]
  low <- threshold > 0 && record[["sampling_proposed"]] >= 1000 && (
    record[["acceptance_rate"]] < threshold ||
      (!is.na(record[["min_sampler_rate"]]) && record[["min_sampler_rate"]] < threshold))
  if (low) {
    message <- paste0(
      .selection_sampler_warning_prefix,
      format(100 * record[["acceptance_rate"]], digits = 3), "% of its proposals per direction after adaptation",
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
    .location_block_report_read()
    .location_block_evidence_read()
    .location_bias_evidence_read()
    .location_scale_report_read()
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
