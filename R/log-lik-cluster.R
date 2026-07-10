# ============================================================================ #
# Cluster-Unit Log-Likelihood Dispatch
# ============================================================================ #

.log_lik_cluster_from_setup <- function(setup) {

  outcome_type      <- setup[["outcome_type"]]
  is_weightfunction <- setup[["is_weightfunction"]]
  yi                <- setup[["yi"]]
  sei               <- setup[["sei"]]
  cluster_setup     <- setup

  if (outcome_type == "norm" && !is_weightfunction &&
      is.null(setup[["weights"]])) {
    if (setup[["effect_direction"]] == "negative") {
      cluster_setup[["mu"]] <- -cluster_setup[["mu"]]
      yi                    <- -yi
    }

    return(.log_lik_cluster_norm_analytic(
      setup = cluster_setup,
      yi    = yi,
      vi    = sei^2
    ))
  }

  if (outcome_type == "norm") {
    if (setup[["effect_direction"]] == "negative" && !is_weightfunction) {
      cluster_setup[["mu"]] <- -cluster_setup[["mu"]]
      yi                    <- -yi
    }

    selection_context <- if (is_weightfunction) {
      .selection_context_from_parts(
        fit                  = setup[["fit"]],
        data                 = setup[["data"]],
        priors               = setup[["priors"]],
        posterior_samples    = setup[["posterior_samples"]],
        effect_direction     = setup[["effect_direction"]]
      )
    } else {
      NULL
    }

    return(.log_lik_cluster_norm_quadrature(
      setup             = cluster_setup,
      yi                = yi,
      sei               = sei,
      is_weightfunction = is_weightfunction,
      selection_context = selection_context
    ))
  }

  if (outcome_type %in% c("bin", "pois")) {
    return(.log_lik_cluster_glmm(
      setup        = cluster_setup,
      data         = setup[["data"]],
      priors       = setup[["priors"]],
      outcome_type = outcome_type
    ))
  }

  stop("Unsupported outcome type for cluster-unit log-likelihood.",
       call. = FALSE)
}


.log_lik_cluster_sum_from_setup <- function(setup) {

  outcome_type      <- setup[["outcome_type"]]
  is_weightfunction <- setup[["is_weightfunction"]]
  yi                <- setup[["yi"]]
  sei               <- setup[["sei"]]
  cluster_setup     <- setup

  if (outcome_type == "norm" && !is_weightfunction &&
      is.null(setup[["weights"]])) {
    if (setup[["effect_direction"]] == "negative") {
      cluster_setup[["mu"]] <- -cluster_setup[["mu"]]
      yi                    <- -yi
    }

    return(rowSums(.log_lik_cluster_norm_analytic(
      setup = cluster_setup,
      yi    = yi,
      vi    = sei^2
    )))
  }

  if (outcome_type == "norm") {
    if (setup[["effect_direction"]] == "negative" && !is_weightfunction) {
      cluster_setup[["mu"]] <- -cluster_setup[["mu"]]
      yi                    <- -yi
    }

    selection_context <- if (is_weightfunction) {
      .selection_context_from_parts(
        fit                  = setup[["fit"]],
        data                 = setup[["data"]],
        priors               = setup[["priors"]],
        posterior_samples    = setup[["posterior_samples"]],
        effect_direction     = setup[["effect_direction"]]
      )
    } else {
      NULL
    }

    return(.log_lik_cluster_norm_quadrature_sum(
      setup             = cluster_setup,
      yi                = yi,
      sei               = sei,
      is_weightfunction = is_weightfunction,
      selection_context = selection_context
    ))
  }

  if (outcome_type %in% c("bin", "pois")) {
    return(.log_lik_cluster_glmm_sum(
      setup        = cluster_setup,
      data         = setup[["data"]],
      priors       = setup[["priors"]],
      outcome_type = outcome_type
    ))
  }

  stop("Unsupported outcome type for cluster-unit log-likelihood.",
       call. = FALSE)
}



# ---------------------------------------------------------------------------- #
# .log_lik_cluster.brma
# ---------------------------------------------------------------------------- #
#
# Cluster-unit likelihood for multilevel models.
#
# Each column is the joint held-out-cluster log-likelihood contribution.
#
# @param object brma object.
#
# @return S x G matrix of log-likelihood values.
#
# ---------------------------------------------------------------------------- #
.log_lik_cluster.brma <- function(object) {

  return(.log_lik_from_posterior_samples(
    fit                  = object[["fit"]],
    posterior_samples    = .get_posterior_samples(object[["fit"]]),
    data                 = object[["data"]],
    priors               = object[["priors"]],
    unit                 = "cluster",
    add_metadata         = TRUE,
    data_hash            = .get_outcome_hash(object)
  ))
}



# ---------------------------------------------------------------------------- #
# .add_cluster_log_lik_metadata
# ---------------------------------------------------------------------------- #
#
# @param log_lik         S x G log-likelihood matrix.
# @param cluster_indices named list of cluster index vectors.
# @param data_hash       character; hash of the outcome target.
#
# @return log-likelihood matrix with names and metadata.
#
# ---------------------------------------------------------------------------- #
.add_cluster_log_lik_metadata <- function(log_lik, cluster_indices, data_hash) {

  cluster_labels    <- names(cluster_indices)
  colnames(log_lik) <- paste0("log_lik_cluster[", cluster_labels, "]")
  attr(log_lik, "RoBMA_target") <- list(
    unit               = "cluster",
    conditioning_depth = "cluster",
    target             = "cluster_joint",
    n                  = length(cluster_indices),
    targets            = cluster_labels,
    data_hash          = data_hash
  )

  return(log_lik)
}
