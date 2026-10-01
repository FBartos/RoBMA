# ============================================================================ #
# GLMM Cluster-Unit Log-Likelihoods
# ============================================================================ #

.log_lik_cluster_glmm <- function(setup, data, priors, outcome_type) {

  stop(
    "Cluster-unit GLMM log-likelihood is unavailable: certified nested ",
    "adaptive quadrature for the cluster and estimate-level latent effects ",
    "has not been implemented. Use unit = 'estimate'.",
    call. = FALSE
  )
}


.log_lik_cluster_glmm_sum <- function(setup, data, priors, outcome_type) {

  .log_lik_cluster_glmm(
    setup        = setup,
    data         = data,
    priors       = priors,
    outcome_type = outcome_type
  )
}
