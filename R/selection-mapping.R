# ============================================================================ #
# Selection-Model Mapping Helpers
# ============================================================================ #
#
# R-side weighted-normal helpers must mirror JAGS dwnorm_mix: each posterior
# sample's bias branch may use only a subset of the global cutpoint grid.
#
# ============================================================================ #


# ---------------------------------------------------------------------------- #
# .selection_has_mapping
# ---------------------------------------------------------------------------- #
#
# @return TRUE when branch-specific selection mapping can be used.
#
# ---------------------------------------------------------------------------- #
.selection_has_mapping <- function(bias_indicator, crit_yi_mapping, crit_yi_mapping_max) {

  return(
    !is.null(bias_indicator) &&
      !is.null(crit_yi_mapping) &&
      !is.null(crit_yi_mapping_max)
  )
}


# ---------------------------------------------------------------------------- #
# .selection_branch_apply
# ---------------------------------------------------------------------------- #
#
# Apply normal or mapped weighted callbacks branch-by-branch.
#
# @param S                    number of posterior samples.
# @param bias_indicator       integer vector of branch IDs, length S.
# @param crit_yi_mapping      cutpoint mapping matrix, cuts x branches.
# @param crit_yi_mapping_max  active cutpoint count per branch.
# @param normal_fun           function(rows).
# @param weighted_fun         function(rows, map_idx).
#
# @return numeric vector of length S.
#
# ---------------------------------------------------------------------------- #
.selection_branch_apply <- function(S, bias_indicator, crit_yi_mapping,
                                    crit_yi_mapping_max, normal_fun,
                                    weighted_fun) {

  bias_indicator <- as.integer(bias_indicator)

  if (length(bias_indicator) != S) {
    stop("'bias_indicator' length does not match posterior samples.", call. = FALSE)
  }

  out <- rep(NA_real_, S)

  for (branch in sort(unique(bias_indicator))) {
    rows        <- which(bias_indicator == branch)
    mapping_max <- crit_yi_mapping_max[branch]

    if (is.na(mapping_max) || mapping_max == 0) {
      out[rows] <- normal_fun(rows)
    } else {
      map_idx <- as.integer(crit_yi_mapping[seq_len(mapping_max), branch])
      out[rows] <- weighted_fun(rows, map_idx)
    }
  }

  return(out)
}


# ---------------------------------------------------------------------------- #
# .selection_active_omega
# ---------------------------------------------------------------------------- #
#
# Mirror DWNMIX.cc omega mapping: first global omega plus weights after each
# active threshold.
#
# ---------------------------------------------------------------------------- #
.selection_active_omega <- function(omega, rows, map_idx) {

  return(cbind(omega[rows, 1], omega[rows, map_idx + 1, drop = FALSE]))
}
