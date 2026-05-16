# ============================================================================ #
# GLMM Cluster-Unit Log-Likelihoods
# ============================================================================ #

# ---------------------------------------------------------------------------- #
# .log_lik_cluster_glmm.brma
# ---------------------------------------------------------------------------- #
#
# Gamma-quadrature cluster-unit likelihood for binomial and Poisson multilevel
# models. Conditional on gamma, existing estimate-level GLMM quadrature helpers
# are reused.
#
# @param object brma object.
#
# @return S x G cluster-unit log-likelihood matrix.
#
# ---------------------------------------------------------------------------- #
.log_lik_cluster_glmm.brma <- function(object) {

  setup       <- .log_lik_cluster_setup.brma(object)
  data        <- setup[["data"]]
  priors      <- setup[["priors"]]
  outcome_type <- .outcome_type(object)

  log_lik <- .log_lik_cluster_glmm(
    setup        = setup,
    data         = data,
    priors       = priors,
    outcome_type = outcome_type
  )

  return(.add_cluster_log_lik_metadata(log_lik, setup[["cluster"]], setup[["data_hash"]]))
}



.log_lik_cluster_glmm <- function(setup, data, priors, outcome_type) {

  if (.has_native_glmm_cluster()) {
    return(.log_lik_cluster_glmm_native(
      setup        = setup,
      data         = data,
      priors       = priors,
      outcome_type = outcome_type
    ))
  }

  return(.log_lik_cluster_glmm_r(
    setup        = setup,
    data         = data,
    priors       = priors,
    outcome_type = outcome_type
  ))
}


.log_lik_cluster_glmm_sum <- function(setup, data, priors, outcome_type) {

  if (.has_native_glmm_row_sum(outcome_type, cluster = TRUE)) {
    return(.log_lik_cluster_glmm_native_sum(
      setup        = setup,
      data         = data,
      priors       = priors,
      outcome_type = outcome_type
    ))
  }

  return(rowSums(.log_lik_cluster_glmm(
    setup        = setup,
    data         = data,
    priors       = priors,
    outcome_type = outcome_type
  )))
}



# ---------------------------------------------------------------------------- #
# .log_lik_cluster_glmm_r
# ---------------------------------------------------------------------------- #
#
# R-composed reference implementation for GLMM cluster likelihoods.
#
# ---------------------------------------------------------------------------- #
.log_lik_cluster_glmm_r <- function(setup, data, priors, outcome_type,
                                    n_theta = 15,
                                    n_gamma = .get_cluster_likelihood_n_gamma(),
                                    n_pi = 30, n_phi = 30) {

  log_lik_fun <- switch(
    outcome_type,
    "bin" = function(idx, mu_node) {
      .outcome_pdf.binom(
        ai         = data[["outcome"]][["ai"]][idx],
        ci         = data[["outcome"]][["ci"]][idx],
        n1i        = data[["outcome"]][["n1i"]][idx],
        n2i        = data[["outcome"]][["n2i"]][idx],
        mu_samples = mu_node,
        tau_within = setup[["tau_within"]][, idx, drop = FALSE],
        prior_pi   = priors[["outcome"]][["pi"]],
        weights    = setup[["weights"]][idx],
        n_theta    = n_theta,
        n_pi       = n_pi
      )
    },
    "pois" = function(idx, mu_node) {
      .outcome_pdf.pois(
        x1i        = data[["outcome"]][["x1i"]][idx],
        x2i        = data[["outcome"]][["x2i"]][idx],
        t1i        = data[["outcome"]][["t1i"]][idx],
        t2i        = data[["outcome"]][["t2i"]][idx],
        mu_samples = mu_node,
        tau_within = setup[["tau_within"]][, idx, drop = FALSE],
        prior_phi  = priors[["outcome"]][["phi"]],
        weights    = setup[["weights"]][idx],
        n_theta    = n_theta,
        n_phi      = n_phi
      )
    }
  )

  log_lik <- .log_lik_cluster_gamma_quadrature(
    cluster_indices     = setup[["cluster"]],
    mu_samples          = setup[["mu"]],
    tau_between_samples = setup[["tau_between"]],
    log_lik_fun         = log_lik_fun,
    n_gamma             = n_gamma
  )

  return(log_lik)
}



# ---------------------------------------------------------------------------- #
# .log_lik_cluster_glmm_native
# ---------------------------------------------------------------------------- #
#
# Native batched implementation for GLMM cluster likelihoods.
#
# ---------------------------------------------------------------------------- #
.log_lik_cluster_glmm_native <- function(setup, data, priors, outcome_type,
                                         n_theta = 15,
                                         n_gamma = .get_cluster_likelihood_n_gamma(),
                                         n_pi = 30, n_phi = 30) {

  gh_theta <- .gauss_hermite_nodes(n_theta)
  gh_gamma <- .gauss_hermite_nodes(n_gamma)
  cluster  <- .cluster_indices_flatten(setup[["cluster"]])
  weights  <- if (is.null(setup[["weights"]])) NULL else .native_numeric_vector(setup[["weights"]])

  if (outcome_type == "bin") {
    pi_grid <- .glmm_binom_logit_pi_grid(
      ai       = data[["outcome"]][["ai"]],
      ci       = data[["outcome"]][["ci"]],
      n1i      = data[["outcome"]][["n1i"]],
      n2i      = data[["outcome"]][["n2i"]],
      prior_pi = priors[["outcome"]][["pi"]],
      n_pi     = n_pi
    )

    return(.Call(
      "RoBMA_glmm_binom_cluster_loglik",
      .native_integer_vector(data[["outcome"]][["ai"]]),
      .native_integer_vector(data[["outcome"]][["ci"]]),
      .native_integer_vector(data[["outcome"]][["n1i"]]),
      .native_integer_vector(data[["outcome"]][["n2i"]]),
      .native_numeric_matrix(setup[["mu"]]),
      .native_numeric_matrix(setup[["tau_within"]]),
      .native_numeric_matrix(setup[["tau_between"]]),
      cluster[["index"]],
      cluster[["size"]],
      weights,
      .native_numeric_vector(gh_theta[["nodes"]]),
      .native_numeric_vector(gh_theta[["log_weights"]]),
      .native_numeric_matrix(pi_grid[["grid"]]),
      .native_numeric_matrix(pi_grid[["log_weights"]]),
      .native_numeric_vector(gh_gamma[["nodes"]]),
      .native_numeric_vector(gh_gamma[["log_weights"]]),
      PACKAGE = "RoBMA"
    ))
  }

  if (outcome_type == "pois") {
    phi_grid <- .glmm_pois_log_phi_grid(
      x1i       = data[["outcome"]][["x1i"]],
      x2i       = data[["outcome"]][["x2i"]],
      t1i       = data[["outcome"]][["t1i"]],
      t2i       = data[["outcome"]][["t2i"]],
      prior_phi = priors[["outcome"]][["phi"]],
      n_phi     = n_phi
    )

    return(.Call(
      "RoBMA_glmm_pois_cluster_loglik",
      .native_integer_vector(data[["outcome"]][["x1i"]]),
      .native_integer_vector(data[["outcome"]][["x2i"]]),
      .native_numeric_vector(data[["outcome"]][["t1i"]]),
      .native_numeric_vector(data[["outcome"]][["t2i"]]),
      .native_numeric_matrix(setup[["mu"]]),
      .native_numeric_matrix(setup[["tau_within"]]),
      .native_numeric_matrix(setup[["tau_between"]]),
      cluster[["index"]],
      cluster[["size"]],
      weights,
      .native_numeric_vector(gh_theta[["nodes"]]),
      .native_numeric_vector(gh_theta[["log_weights"]]),
      .native_numeric_matrix(phi_grid[["grid"]]),
      .native_numeric_matrix(phi_grid[["log_weights"]]),
      .native_numeric_vector(gh_gamma[["nodes"]]),
      .native_numeric_vector(gh_gamma[["log_weights"]]),
      PACKAGE = "RoBMA"
    ))
  }

  stop("Unsupported outcome type for GLMM cluster-unit log-likelihood.",
       call. = FALSE)
}


.log_lik_cluster_glmm_native_sum <- function(setup, data, priors,
                                             outcome_type, n_theta = 15,
                                             n_gamma = .get_cluster_likelihood_n_gamma(),
                                             n_pi = 30, n_phi = 30) {

  gh_theta <- .gauss_hermite_nodes(n_theta)
  gh_gamma <- .gauss_hermite_nodes(n_gamma)
  cluster  <- .cluster_indices_flatten(setup[["cluster"]])
  weights  <- if (is.null(setup[["weights"]])) NULL else .native_numeric_vector(setup[["weights"]])

  if (outcome_type == "bin") {
    pi_grid <- .glmm_binom_logit_pi_grid(
      ai       = data[["outcome"]][["ai"]],
      ci       = data[["outcome"]][["ci"]],
      n1i      = data[["outcome"]][["n1i"]],
      n2i      = data[["outcome"]][["n2i"]],
      prior_pi = priors[["outcome"]][["pi"]],
      n_pi     = n_pi
    )

    return(.Call(
      "RoBMA_glmm_binom_cluster_loglik_row_sum",
      .native_integer_vector(data[["outcome"]][["ai"]]),
      .native_integer_vector(data[["outcome"]][["ci"]]),
      .native_integer_vector(data[["outcome"]][["n1i"]]),
      .native_integer_vector(data[["outcome"]][["n2i"]]),
      .native_numeric_matrix(setup[["mu"]]),
      .native_numeric_matrix(setup[["tau_within"]]),
      .native_numeric_matrix(setup[["tau_between"]]),
      cluster[["index"]],
      cluster[["size"]],
      weights,
      .native_numeric_vector(gh_theta[["nodes"]]),
      .native_numeric_vector(gh_theta[["log_weights"]]),
      .native_numeric_matrix(pi_grid[["grid"]]),
      .native_numeric_matrix(pi_grid[["log_weights"]]),
      .native_numeric_vector(gh_gamma[["nodes"]]),
      .native_numeric_vector(gh_gamma[["log_weights"]]),
      PACKAGE = "RoBMA"
    ))
  }

  if (outcome_type == "pois") {
    phi_grid <- .glmm_pois_log_phi_grid(
      x1i       = data[["outcome"]][["x1i"]],
      x2i       = data[["outcome"]][["x2i"]],
      t1i       = data[["outcome"]][["t1i"]],
      t2i       = data[["outcome"]][["t2i"]],
      prior_phi = priors[["outcome"]][["phi"]],
      n_phi     = n_phi
    )

    return(.Call(
      "RoBMA_glmm_pois_cluster_loglik_row_sum",
      .native_integer_vector(data[["outcome"]][["x1i"]]),
      .native_integer_vector(data[["outcome"]][["x2i"]]),
      .native_numeric_vector(data[["outcome"]][["t1i"]]),
      .native_numeric_vector(data[["outcome"]][["t2i"]]),
      .native_numeric_matrix(setup[["mu"]]),
      .native_numeric_matrix(setup[["tau_within"]]),
      .native_numeric_matrix(setup[["tau_between"]]),
      cluster[["index"]],
      cluster[["size"]],
      weights,
      .native_numeric_vector(gh_theta[["nodes"]]),
      .native_numeric_vector(gh_theta[["log_weights"]]),
      .native_numeric_matrix(phi_grid[["grid"]]),
      .native_numeric_matrix(phi_grid[["log_weights"]]),
      .native_numeric_vector(gh_gamma[["nodes"]]),
      .native_numeric_vector(gh_gamma[["log_weights"]]),
      PACKAGE = "RoBMA"
    ))
  }

  stop("Unsupported outcome type for GLMM cluster-unit log-likelihood.",
       call. = FALSE)
}
