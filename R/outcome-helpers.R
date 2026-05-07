# ============================================================================ #
# outcome-helpers.R
# ============================================================================ #
#
# Shared outcome-data extraction helpers used by residuals, diagnostics, plots,
# likelihood evaluation, and zplot code.
#
# ============================================================================ #


# ---------------------------------------------------------------------------- #
# .outcome_data_yi
# ---------------------------------------------------------------------------- #
#
# Extract or compute observed effect sizes from brma object.
#
# For normal models, yi is directly available in the data.
# For GLMM models, yi is computed from raw frequency data using formulas
# equivalent to metafor::escalc with default zero-cell adjustment (add = 0.5).
#
# @param object brma object
#
# @return numeric vector of observed effect sizes
#
# ---------------------------------------------------------------------------- #
.outcome_data_yi <- function(object) {

  outcome_type <- .outcome_type(object)
  outcome_data <- object[["data"]][["outcome"]]

  if (outcome_type == "norm") {
    # normal models: yi is directly available
    return(outcome_data[["yi"]])
  } else if (outcome_type == "bin") {
    # binomial models: compute log odds ratio from 2x2 table
    # using metafor::escalc formula with default zero-cell adjustment (add = 0.5)
    ai <- outcome_data[["ai"]]
    bi <- outcome_data[["n1i"]] - outcome_data[["ai"]] # bi = n1i - ai
    ci <- outcome_data[["ci"]]
    di <- outcome_data[["n2i"]] - outcome_data[["ci"]] # di = n2i - ci
    add <- 0.5

    # vectorized computation with conditional zero-cell adjustment
    needs_adj <- (ai == 0 | bi == 0 | ci == 0 | di == 0)
    ai_adj <- ai + add * needs_adj
    bi_adj <- bi + add * needs_adj
    ci_adj <- ci + add * needs_adj
    di_adj <- di + add * needs_adj

    return(log((ai_adj * di_adj) / (bi_adj * ci_adj)))
  } else if (outcome_type == "pois") {
    # Poisson models: compute log incidence rate ratio
    x1i <- outcome_data[["x1i"]]
    t1i <- outcome_data[["t1i"]]
    x2i <- outcome_data[["x2i"]]
    t2i <- outcome_data[["t2i"]]
    add <- 0.5

    # vectorized computation with conditional zero-cell adjustment
    needs_adj <- (x1i == 0 | x2i == 0)
    x1i_adj <- x1i + add * needs_adj
    x2i_adj <- x2i + add * needs_adj

    return(log((x1i_adj / t1i) / (x2i_adj / t2i)))
  }
}


# ---------------------------------------------------------------------------- #
# .outcome_data_sei
# ---------------------------------------------------------------------------- #
#
# Extract or compute sampling standard errors from brma object.
#
# For normal models, sei is directly available in the data.
# For GLMM models, sei is computed from raw frequency data using formulas
# equivalent to metafor::escalc with default zero-cell adjustment (add = 0.5).
#
# @param object brma object
#
# @return numeric vector of sampling standard errors
#
# ---------------------------------------------------------------------------- #
.outcome_data_sei <- function(object) {

  outcome_type <- .outcome_type(object)
  outcome_data <- object[["data"]][["outcome"]]

  if (outcome_type == "norm") {
    # normal models: sei is directly available
    return(outcome_data[["sei"]])
  } else if (outcome_type == "bin") {
    # binomial models: compute approximate SE for log odds ratio
    # SE(logOR) = sqrt(1/ai + 1/bi + 1/ci + 1/di) with zero-cell adjustment
    ai <- outcome_data[["ai"]]
    bi <- outcome_data[["n1i"]] - outcome_data[["ai"]]
    ci <- outcome_data[["ci"]]
    di <- outcome_data[["n2i"]] - outcome_data[["ci"]]
    add <- 0.5

    # vectorized computation with conditional zero-cell adjustment
    needs_adj <- (ai == 0 | bi == 0 | ci == 0 | di == 0)
    ai_adj <- ai + add * needs_adj
    bi_adj <- bi + add * needs_adj
    ci_adj <- ci + add * needs_adj
    di_adj <- di + add * needs_adj

    return(sqrt(1 / ai_adj + 1 / bi_adj + 1 / ci_adj + 1 / di_adj))
  } else if (outcome_type == "pois") {
    # Poisson models: compute approximate SE for log IRR
    # SE(logIRR) = sqrt(1/x1i + 1/x2i)
    x1i <- outcome_data[["x1i"]]
    x2i <- outcome_data[["x2i"]]
    add <- 0.5

    # vectorized computation with conditional zero-cell adjustment
    needs_adj <- (x1i == 0 | x2i == 0)
    x1i_adj <- x1i + add * needs_adj
    x2i_adj <- x2i + add * needs_adj

    return(sqrt(1 / x1i_adj + 1 / x2i_adj))
  }
}


# ---------------------------------------------------------------------------- #
# .outcome_data_vi
# ---------------------------------------------------------------------------- #
#
# Extract or compute sampling variances from brma object.
#
# Convenience wrapper that returns sei^2 from .outcome_data_sei().
#
# @param object brma object
#
# @return numeric vector of sampling variances
#
# ---------------------------------------------------------------------------- #
.outcome_data_vi <- function(object) {

  return(.outcome_data_sei(object)^2)
}


# ---------------------------------------------------------------------------- #
# .outcome_data_weights
# ---------------------------------------------------------------------------- #
#
# Extract likelihood weights from the outcome data, defaulting to unit weights.
#
# ---------------------------------------------------------------------------- #
.outcome_data_weights <- function(object) {

  outcome_data <- object[["data"]][["outcome"]]

  if (.is_weights(object)) {
    return(as.numeric(outcome_data[["weights"]]))
  }

  return(rep(1, nrow(outcome_data)))
}
