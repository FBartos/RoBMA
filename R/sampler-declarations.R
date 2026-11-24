# ============================================================================ #
# sampler-declarations.R
# ============================================================================ #
#
# RoBMA passes declarations to BayesTools::JAGS_fit() for sampler ownership,
# probability recording and the selection-weight representation. Joint
# declarations use the bound priors and data that the backend has compiled.
#
# ============================================================================ #


# The declarations of a fresh fit, from the data, priors and fit priors of .fit().
.fit_sampler_declarations <- function(data, priors, fit_priors) {

  module_samplers <- .sampler_module_samplers()
  requested <- list(
    joint          = RoBMA.get_option("sampler.joint"),
    nested_weights = RoBMA.get_option("sampler.nested.weights"),
    probabilities  = RoBMA.get_option("sampler.probabilities"),
    specialized    = module_samplers
  )
  nested_weights  <- module_samplers && requested[["nested_weights"]] &&
    inherits(fit_priors[["bias"]], "prior.bias_mixture") &&
    .is_priors_weightfunction(priors) && !.selection_all_sources_conditioned(data)
  joint <- .sampler_joint_declared(data, fit_priors, requested[["joint"]])
  requested[["joint_sampler"]] <- attr(joint, "sampler", exact = TRUE)

  return(list(
    sampler_gates           = .sampler_gates_declared(fit_priors),
    sampler_blocks          = .sampler_blocks_declared(data, priors),
    sampler_joint           = joint,
    sampler_nested_weights = nested_weights,
    sampler_probabilities   = if (module_samplers) requested[["probabilities"]] else character(),
    requested               = requested
  ))
}


# Keep the original request and refresh ownership from the last successful run.
# The backend keeps the markers, nested layout and probability monitor plans.
.sampler_request_record <- function(fit, requested = NULL) {

  declarations <- attr(fit, "sampler_declarations", exact = TRUE)
  if (is.null(requested)) requested <- declarations[["robma"]][["request"]]
  if (is.null(requested)) return(fit)

  record <- list()
  workers      <- fit[["samplers_by_worker"]]
  assignments <- if (is.null(workers)) list(fit[["samplers"]]) else lapply(workers, `[[`, "samplers")
  disabled     <- !requested[["specialized"]]
  disabled_reason <- "The 'jags.default' sampler policy disables specialized samplers."
  policy       <- requested[["joint"]]
  expected     <- switch(policy, auto = requested[["joint_sampler"]], bias = "RoBMA::IndependentBiasJoint",
                         gaussian = "RoBMA::GaussianJoint", location = "RoBMA::LocationBlock",
                         clustered = "RoBMA::JointCluster", clustered.scale = "RoBMA::JointCluster")
  if (!is.null(expected)) {
    claimed <- lapply(assignments, function(table) table[["Node"]][table[["Sampler"]] == expected])
    owned   <- Reduce(intersect, claimed)
    mixed   <- any(!vapply(claimed, setequal, logical(1), y = owned))
    budget  <- requested[["joint_budget"]]
    status <- if (disabled) "disabled" else if (mixed) "mixed" else if (length(owned) > 0L) "active" else
      if (!is.null(budget)) "budget_bypassed" else
      if (is.null(declarations[["joint"]])) "graph_declined" else "native_declined"
    reason <- switch(status,
      disabled = disabled_reason,
      mixed = "Workers assigned different nodes to the joint sampler. See the per-worker sampler assignments.",
      budget_bypassed = sprintf(
        "The fit has %d retained estimates and exceeds the automatic Gaussian joint sampler limit of %d.",
        budget[["rows"]], budget[["limit"]]),
      graph_declined = switch(expected,
        "RoBMA::IndependentBiasJoint" = paste(
          "This update requires independent step-selection rows and eligible effect and bias mixtures.",
          "Moderated fits require a supported additive design and ordinary selection weights."),
        "RoBMA::GaussianJoint" = if (policy == "auto")
          "This update requires an unweighted independent Gaussian fit, one standardized numeric moderator and eligible normal mixtures." else
          "This update requires an unweighted independent Gaussian fit, a single continuous moderator and eligible normal mixtures.",
        "RoBMA::LocationBlock" = paste(
          "This update requires a fit of a Gaussian outcome (multilevel, random-effect, or with a gated moderator) with a location term",
          "whose prior mixes an exact zero with a normal slab, and no exact probability readout of its indicator."),
        "This update requires an eligible clustered fit with one factor moderator and supported prior laws."),
      native_declined = if (expected == "RoBMA::GaussianJoint")
        "The compiled graph or Gaussian proposal coverage did not qualify. See the fit's sampler assignments." else
        "The native sampler did not claim the declared nodes. See the fit's sampler assignments."
    )
    record[["joint"]] <- list(policy = policy, status = status, reason = reason, sampler = expected, nodes = owned)
    if (expected == "RoBMA::LocationBlock") {
      record[["joint"]][["counters"]] <- .location_block_record(attr(fit, "runtime_report_values", exact = TRUE))
      record[["joint"]][["scale"]] <- .location_scale_record(attr(fit, "runtime_report_values", exact = TRUE))
    }
  }

  if (requested[["nested_weights"]]) {
    nested <- declarations[["nested_weights"]]
    if (is.null(nested)) nested <- list(realized = FALSE, branches = integer(),
      reason = if (disabled) disabled_reason else "The fit has no applicable selection-weight mixture.")
    record[["nested_weights"]] <- nested
  }

  if (length(requested[["probabilities"]]) > 0L) {
    requested_nodes <- requested[["probabilities"]]
    nodes <- declarations[["probabilities"]][["node"]]
    if (is.null(nodes)) nodes <- requested_nodes
    end_state <- fit[["end_state"]]
    chains <- if (is.null(workers)) list(seq_along(end_state)) else lapply(workers, `[[`, "chains")
    owners <- vapply(nodes, function(node) {
      node_owners <- unique(vapply(seq_along(assignments), function(worker) {
        .sampler_probability_owner(node, assignments[[worker]], end_state[chains[[worker]]])
      }, character(1)))
      if (length(node_owners) == 1L) node_owners[[1L]] else "mixed"
    }, character(1), USE.NAMES = FALSE)
    recorded <- vapply(nodes, function(node) {
      inherits(fit[["sampler_probabilities"]][[node]], "mcmc.list")
    }, logical(1))
    # The backend resolves fitted aliases when it binds the monitor. A recorded
    # trace identifies the actual finite sampler even when JAGS reports another
    # spelling for its stochastic node in the assignment table.
    owners[recorded] <- "BayesTools::Finite"
    record[["probabilities"]] <- data.frame(
      node = nodes, requested_node = requested_nodes, sampler = owners, status = if (disabled) "disabled" else
        ifelse(recorded, "active", "unavailable"),
      reason = if (disabled) disabled_reason else
        ifelse(recorded, NA_character_, "The fitted sampler supplies no conditional branch probabilities.")
    )
  }
  evidence <- .sampler_evidence_record(fit, declarations, assignments, workers)
  if (!is.null(evidence)) record[["evidence"]] <- evidence
  if (length(record) > 0L) {
    record[["request"]] <- requested
    declarations[["robma"]] <- record
    attr(fit, "sampler_declarations") <- declarations
  }
  return(fit)
}


# The evidence record of indicators that do not switch, from the runtime
# reports of this call's retained runs (updates after adaptation): per
# location gate and chain the flips proposed and accepted in each direction,
# the proposals the free screen removed and, for those that reached stage 1,
# the count, sum and maximum of the log collapsed odds of the other state
# against the current one (exact in Gaussian labels, Laplace in selection
# labels); per screened indicator and chain the updates, changes of level and
# the mean and maximum full-conditional probability of leaving the current
# level; for a bias label owned by the location block (`bias`, per chain) the
# proposals and acceptances of its label update, Gaussian-family move, slope
# update and cross-family swap, and per swap direction the count, sum and
# maximum of the stage-1 log ratio of every attempt and of the swaps the exact
# correction kept, with the sweeps ending in each label (`bias_labels`).
# Counts and sums are added and maxima taken over runs and processes;
# NULL without any record. The native screened record numbers its samplers in
# creation order, which is the order of their nodes in the sampler table; a
# mismatch of the counts leaves the names missing.
.sampler_evidence_record <- function(fit, declarations, assignments, workers) {

  values <- attr(fit, "runtime_report_values", exact = TRUE)
  if (length(values) == 0L) return(NULL)
  gates <- character()
  for (marker in declarations[["joint"]][["markers"]]) {
    arguments <- marker[["arguments"]]
    if (identical(marker[["tag"]], "robma_location_group") && !isTRUE(arguments[[3L]] == 0)) {
      gates[as.character(arguments[[5L]])] <- arguments[[1L]]
    }
  }
  global_of <- function(process) {
    chains <- if (is.null(workers)) NULL else workers[[process]][["chains"]]
    function(local) if (is.null(chains)) as.integer(local) + 1L else chains[as.integer(local) + 1L]
  }
  # the screened samplers' nodes, matched in the last run (whose final state the
  # fit keeps) and reused for the same sampler in earlier runs of the process
  known <- list()
  last  <- values[[length(values)]]
  for (process in seq_along(last)) {
    rows <- last[[process]][["screened"]]
    if (!is.matrix(rows) || nrow(rows) == 0L) next
    table  <- assignments[[min(process, length(assignments))]]
    owners <- table[["Node"]][table[["Sampler"]] == "RoBMA::ScreenedFinite"]
    known[[process]] <- stats::setNames(.sampler_screened_names(rows, owners, global_of(process), fit),
                                        rows[, "serial"])
  }
  location <- list()
  screened <- list()
  bias <- list()
  bias_labels <- list()
  indicator <- "bias_indicator"
  for (marker in declarations[["joint"]][["markers"]]) {
    if (identical(marker[["tag"]], "robma_location_bias")) indicator <- marker[["arguments"]][[1L]]
  }
  for (run in values) {
    for (process in seq_along(run)) {
      report <- run[[process]]
      if (is.null(report)) next
      global <- global_of(process)
      rows <- report[["location_evidence"]]
      if (is.matrix(rows) && nrow(rows) > 0L) {
        location[[length(location) + 1L]] <- data.frame(
          node = unname(gates[as.character(rows[, "slot"])]), chain = global(rows[, "chain"]),
          rows[, setdiff(colnames(rows), c("slot", "chain")), drop = FALSE], stringsAsFactors = FALSE)
      }
      rows <- report[["location_bias"]]
      if (is.list(rows) && is.matrix(rows[["chains"]]) && nrow(rows[["chains"]]) > 0L) {
        counts <- rows[["chains"]]
        bias[[length(bias) + 1L]] <- data.frame(node = indicator, chain = global(counts[, "chain"]),
          counts[, setdiff(colnames(counts), "chain"), drop = FALSE], stringsAsFactors = FALSE)
      }
      if (is.list(rows) && is.matrix(rows[["labels"]]) && nrow(rows[["labels"]]) > 0L) {
        labels <- rows[["labels"]]
        bias_labels[[length(bias_labels) + 1L]] <- data.frame(chain = global(labels[, "chain"]),
          label = as.integer(labels[, "label"]), sweeps = labels[, "sweeps"])
      }
      rows <- report[["screened"]]
      if (is.matrix(rows) && nrow(rows) > 0L) {
        names <- if (process <= length(known)) known[[process]] else NULL
        node  <- if (is.null(names)) rep(NA_character_, nrow(rows)) else unname(names[as.character(rows[, "serial"])])
        screened[[length(screened) + 1L]] <- data.frame(
          node = node, chain = global(rows[, "chain"]),
          rows[, c("levels", "updates", "changes", "leave_sum", "leave_max"), drop = FALSE], stringsAsFactors = FALSE)
      }
    }
  }
  combine <- function(parts, maxima) {
    if (length(parts) == 0L) return(NULL)
    all <- do.call(rbind, parts)
    key <- paste(all[["node"]], all[["chain"]])
    out <- do.call(rbind, lapply(split(all, factor(key, levels = unique(key))), function(rows) {
      row <- rows[1L, , drop = FALSE]
      for (column in setdiff(names(rows), c("node", "chain", "levels"))) {
        row[[column]] <- if (column %in% maxima) max(rows[[column]]) else sum(rows[[column]])
      }
      row
    }))
    out <- out[order(out[["node"]], out[["chain"]]), , drop = FALSE]
    rownames(out) <- NULL
    out
  }
  bias <- combine(bias, c("leave_selection_log_ratio_max", "leave_gaussian_log_ratio_max",
                          "kept_leave_selection_log_ratio_max", "kept_leave_gaussian_log_ratio_max"))
  if (!is.null(bias)) {
    for (column in grep("_max$", names(bias), value = TRUE)) bias[[column]][!is.finite(bias[[column]])] <- NA_real_
  }
  if (length(bias_labels) > 0L) {
    bias_labels <- do.call(rbind, bias_labels)
    bias_labels <- stats::aggregate(sweeps ~ chain + label, data = bias_labels, FUN = sum)
    bias_labels <- bias_labels[order(bias_labels[["chain"]], bias_labels[["label"]]), , drop = FALSE]
    rownames(bias_labels) <- NULL
  } else {
    bias_labels <- NULL
  }
  location <- combine(location, "log_odds_max")
  if (!is.null(location)) location[["log_odds_max"]][!is.finite(location[["log_odds_max"]])] <- NA_real_
  screened <- combine(screened, "leave_max")
  if (!is.null(screened)) {
    screened[["leave_mean"]] <- screened[["leave_sum"]] / screened[["updates"]]
    screened <- screened[, c("node", "chain", "levels", "updates", "changes", "leave_mean", "leave_max")]
  }
  if (is.null(location) && is.null(screened) && is.null(bias)) return(NULL)
  out <- list(location = location, screened = screened)
  if (!is.null(bias)) out[c("bias", "bias_labels")] <- list(bias, bias_labels)
  return(out)
}


# The node of each screened sampler's rows (one process): native samplers know
# no node names, and JAGS orders its sampler table differently from the order
# the factory creates them, so a serial is matched to the screened node whose
# final value equals the sampler's last value in every chain of the process and
# whose monitored draws lie in the sampler's support. A serial without exactly
# one such node is left unnamed.
.sampler_screened_names <- function(rows, owners, global, fit) {

  end_state <- fit[["end_state"]]
  draws     <- tryCatch(as.matrix(fit[["mcmc"]]), error = function(e) NULL)
  value_of  <- function(node, chain) {
    state <- end_state[[chain]]
    base  <- sub("\\[.*$", "", node)
    value <- state[[base]]
    if (is.null(value)) return(NA_real_)
    if (base == node) return(as.numeric(value)[1L])
    index <- as.integer(strsplit(sub("^.*\\[(.*)\\]$", "\\1", node), ",", fixed = TRUE)[[1L]])
    if (length(index) == 1L) return(as.numeric(value)[index])
    as.numeric(value[matrix(index, nrow = 1L)])
  }
  names <- vapply(unique(rows[, "serial"]), function(serial) {
    own    <- rows[rows[, "serial"] == serial, , drop = FALSE]
    chains <- global(own[, "chain"])
    lower  <- own[1L, "lower"]
    upper  <- lower + own[1L, "levels"] - 1
    fits <- vapply(owners, function(node) {
      if (!isTRUE(all(vapply(chains, function(chain) value_of(node, chain), numeric(1)) == own[, "last"]))) return(FALSE)
      if (is.null(draws) || !node %in% colnames(draws)) return(TRUE)
      all(draws[, node] >= lower & draws[, node] <= upper)
    }, logical(1))
    if (sum(fits) == 1L) owners[fits] else NA_character_
  }, character(1))
  unname(names[match(rows[, "serial"], unique(rows[, "serial"]))])
}


# JAGS can report the bare or indexed spelling of a singleton. Its final
# state establishes the shape even for coordinates omitted from monitoring.
.sampler_probability_owner <- function(node, table, states) {

  index <- match(node, table[["Node"]])
  if (!is.na(index)) return(table[["Sampler"]][[index]])

  base <- sub("\\[1\\]$", "", node)
  if (grepl("[", base, fixed = TRUE) || length(states) == 0L ||
      !all(vapply(states, function(state) {
        value <- state[[base]]
        is.numeric(value) && length(value) == 1L && length(dim(value)) <= 1L
      }, logical(1)))) return(NA_character_)

  alias <- if (node == base) paste0(base, "[1]") else base
  index <- match(alias, table[["Node"]])
  if (is.na(index)) return(NA_character_)
  return(table[["Sampler"]][[index]])
}


# A fitting error raised by a BayesTools sampler or by RoBMA's corrected
# selection sampler keeps its message and gains a hint: the option that returns
# the parameters to the stock JAGS samplers, and a request to report the model
# (the density could not be computed, or a declaration is wrong).
.sampler_error_hint <- function(message) {

  # Errors of either block sampler (RoBMA's corrected block sampler and
  # BayesTools' factor slice sampler, which takes the same block under "slice"
  # and stops the same way): only undeclaring the blocks helps. A member on a
  # bound of its support is avoided by a start strictly inside it.
  bound <- if (grepl("current value lies on a bound of its support", message, fixed = TRUE)) {
    "A start strictly inside the support of the parameter avoids this error; otherwise fit"
  } else {
    "Fit"
  }
  if (grepl("Corrected slice (block", message, fixed = TRUE) ||
      grepl("the support bounds of this node are too wide", message, fixed = TRUE)) {
    return(paste0(
      message, "\nThe error was raised by the sampler of a parameter block that RoBMA uses by default. ",
      bound, " the model again after 'RoBMA.options(sampler.blocks = FALSE)' to sample these parameters with ",
      "the ordinary samplers, and please report the model to the RoBMA developers."
    ))
  }
  if (grepl("Factor slice sampler:", message, fixed = TRUE)) {
    return(paste0(
      message, "\nThe error was raised by a sampler of the BayesTools JAGS module that RoBMA uses by default. ",
      bound, " the model again after 'RoBMA.options(sampler.blocks = FALSE)' to sample these parameters with ",
      "the ordinary JAGS samplers, and please report the model to the RoBMA developers."
    ))
  }

  # The corrected scalar and gated block samplers and the screened indicator
  # sampler: the 'slice' tier returns their parameters to the other samplers.
  if (grepl("Corrected slice", message, fixed = TRUE) ||
      grepl("Corrected finite sampler", message, fixed = TRUE) ||
      grepl("Surrogate slice transition failed", message, fixed = TRUE)) {
    return(paste0(
      message, "
The error was raised by the corrected selection sampler that RoBMA uses by default. ",
      "Fit the model again after 'RoBMA.options(selection.sampler = \"slice\")' to sample these parameters ",
      "with the ordinary samplers, and please report the model to the RoBMA developers."
    ))
  }

  if (!grepl("gated node changes the likelihood", message, fixed = TRUE)) {
    return(message)
  }

  return(paste0(
    message, "\nThe error was raised by a sampler of the BayesTools JAGS module that RoBMA uses by default. ",
    "Fit the model again after 'RoBMA.options(sampler.gates = FALSE)' to sample these parameters with the ordinary JAGS ",
    "samplers, and please report the model to the RoBMA developers."
  ))
}


# Whether the samplers of the JAGS modules may be declared at all:
# 'selection.sampler = "jags.default"' keeps only the samplers JAGS assigns.
.sampler_module_samplers <- function() {

  return(!isTRUE(RoBMA.get_option("selection.sampler") %in% "jags.default"))
}


# Gates: a bias-mixture prior, with or without a weight function branch. Under
# 'selection.sampler = "slice.coarse"' the weights of a branch go to RoBMA's
# corrected gated block when its selection children take part (selection
# blocks with a coarse rule, or independent step-selection rows through the
# row cache); the other gated parameters keep BayesTools' gated slice sampler.
.sampler_gates_declared <- function(fit_priors) {

  return(.sampler_module_samplers() && isTRUE(RoBMA.get_option("sampler.gates")) &&
           inherits(fit_priors[["bias"]], "prior.bias_mixture"))
}


# Blocks: NULL unless a selection likelihood is present; otherwise a function
# of the combined prior list (formula priors under their final names and with
# their random-effect attributes), which BayesTools::JAGS_fit() calls once.
# Under 'selection.sampler = "slice.coarse"' RoBMA's corrected block sampler
# takes a declared block with a coarse-capable selection child; BayesTools'
# factor slice sampler takes the others.
.sampler_blocks_declared <- function(data, priors) {

  if (!.sampler_module_samplers() || !isTRUE(RoBMA.get_option("sampler.blocks")) ||
      !.is_priors_weightfunction(priors) || .selection_all_sources_conditioned(data)) {
    return(NULL)
  }

  scale_parameters <- vapply(
    .data_scale_component_specs(data), function(spec) spec[["parameter"]], character(1)
  )

  return(.sampler_block_function(scale_parameters))
}


# A closure over the names of the scale formulas only (it keeps no data).
.sampler_block_function <- function(scale_parameters) {

  force(scale_parameters)
  return(function(prior_list) .sampler_block_members(prior_list, scale_parameters))
}


# The block of a combined prior list: list(<prior names>) or NULL.
.sampler_block_members <- function(prior_list, scale_parameters = character()) {

  if (!.sampler_heterogeneity_estimated(prior_list, scale_parameters)) {
    return(NULL)
  }

  # BayesTools resolves stochastic leaves. Keep supported empty priors in the
  # declaration so it retains the model-level grouping of the supplied priors.
  sources      <- BayesTools::JAGS_sampler_sources(prior_list)
  eligible     <- unique(sources[["prior"]][sources[["status"]] %in% c("eligible", "empty")])
  unit_latents <- c("theta", "gamma", "sampling_z", "pi", "phi")
  members      <- setdiff(eligible, unit_latents)

  if (.sampler_is_ensemble(prior_list)) {
    # Model-averaged fits: the plain allocation SD and the Dirichlet allocation
    # weights of a random-effect variance allocation only.
    members <- members[vapply(prior_list[members], .sampler_is_allocation_block_prior, logical(1))]
  }

  if (length(members) < 2L) {
    return(NULL)
  }

  return(list(members))
}


.sampler_prior_scalar <- function(prior) {

  return(BayesTools::is.prior.simple(prior) &&
           !BayesTools::is.prior.factor(prior) &&
           !BayesTools::is.prior.vector(prior))
}

# An ensemble only groups the SD and Dirichlet weights of variance allocations.
.sampler_is_allocation_block_prior <- function(prior) {

  if (isTRUE(attr(prior, "random_allocation_sd", exact = TRUE))) {
    return(.sampler_prior_scalar(prior))
  }
  return(!is.null(attr(prior, "random_allocation", exact = TRUE)) &&
           identical(prior[["distribution"]], "dirichlet"))
}


# A model-averaged fit: a mixture, a spike-and-slab or an inclusion prior.
.sampler_is_ensemble <- function(prior_list) {

  return(any(vapply(prior_list, function(prior) {
    BayesTools::is.prior.mixture(prior) || BayesTools::is.prior.spike_and_slab(prior) ||
      !is.null(attr(prior, "random_allocation_inclusion", exact = TRUE))
  }, logical(1))))
}


# A heterogeneity parameter is estimated: a non-point prior of 'tau', 'rho', a
# scale formula (location-scale models), a random-effect SD or an allocation SD.
.sampler_heterogeneity_estimated <- function(prior_list, scale_parameters = character()) {

  heterogeneity <- vapply(names(prior_list), function(name) {
    prior <- prior_list[[name]]
    name %in% c("tau", "rho") ||
      isTRUE(attr(prior, "parameter", exact = TRUE) %in% scale_parameters) ||
      isTRUE(attr(prior, "random_sd", exact = TRUE)) ||
      isTRUE(attr(prior, "random_allocation_sd", exact = TRUE))
  }, logical(1))
  point <- vapply(prior_list, BayesTools::is.prior.point, logical(1))

  return(any(heterogeneity & !point))
}
