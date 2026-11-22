# Describe the legacy clustered Gaussian graph using the compiled moderator
# design and prior sources. Multivariate and scale-regression graphs have
# different latent coordinates and covariance structure.
.sampler_joint_clustered <- function(data, prior_list, formula_design, fit_data, scale = FALSE) {

  if (.data_outcome_type(data) != "norm" || !.is_data_multilevel(data) ||
      .is_data_known_v(data) || .is_data_random(data) || .is_data_scale(data) ||
      .is_data_weights(data) || .data_effect_direction(data) != "positive" ||
      any(data[["outcome"]][["bias_exempt"]]) ||
      !identical(names(formula_design), "mu")) return(NULL)

  design <- formula_design[["mu"]]
  terms  <- design[["model_terms_type"]]
  if (length(terms) != 2L || !identical(unname(terms), c("continuous", "factor")) ||
      names(terms)[1L] != "intercept" || length(design[["random_effects"]]) > 0L ||
      length(design[["expressions"]]) > 0L) return(NULL)

  name_map     <- design[["name_map"]]
  fixed        <- name_map[name_map[["kind"]] == "fixed", , drop = FALSE]
  intercept    <- fixed[["jags_name"]][fixed[["term"]] == "intercept"]
  moderator    <- fixed[["jags_name"]][fixed[["term"]] != "intercept"]
  matrix_name  <- design[["jags_data_names"]][[names(terms)[2L]]]
  contrasts    <- fit_data[[matrix_name]]
  if (length(intercept) != 1L || length(moderator) != 1L ||
      !is.matrix(contrasts) || ncol(contrasts) < 1L ||
      !all(is.finite(contrasts)) || nrow(contrasts) != fit_data[["K"]]) return(NULL)

  mu    <- .sampler_joint_normal_mixture(prior_list, intercept)
  group <- .sampler_cluster_moderator(prior_list, moderator, ncol(contrasts))
  tau   <- .sampler_joint_normal_mixture(prior_list, "tau", tau = TRUE)
  rho   <- .sampler_cluster_rho_source(prior_list)
  if (is.null(mu) || is.null(group) || is.null(tau) || is.null(rho) ||
      !all(names(prior_list) %in% c(intercept, moderator, "tau", "rho", "gamma", "bias"))) return(NULL)

  study <- as.integer(fit_data[["cluster"]])
  n_studies <- length(unique(study))
  if (length(study) != fit_data[["K"]] || n_studies < 2L ||
      !identical(sort(unique(study)), seq_len(n_studies))) return(NULL)
  gamma <- BayesTools::JAGS_sampler_sources(prior_list, "gamma")
  gamma <- gamma[["node"]][gamma[["status"]] == "eligible"]
  if (length(gamma) != n_studies) return(NULL)

  rows <- .sampler_cluster_rows(data, prior_list, fit_data)
  if (is.null(rows)) return(NULL)

  member <- function(node, role, slot) {
    list(name = paste0("RoBMA_cluster_member_", role, "_", slot),
         distribution = "drobma_cluster_member", arguments = list(node, role, slot))
  }
  markers <- c(list(member(group[["gate"]], 1L, 1L),
                    member(mu[["source"]], 2L, 1L),
                    member(group[["source"]], 3L, 1L)),
               lapply(seq_along(gamma), function(i) member(gamma[i], 4L, i)))
  if (scale) markers <- c(markers, list(member(tau[["source"]], 5L, 1L), member(rho, 6L, 1L)))

  # Store the actual fitted contrast values in the same row order as the
  # likelihood. The sampler uses only fixed geometry and live nuisance nodes.
  cluster_data <- list(
    RoBMA_cluster_X       = as.numeric(contrasts[rows[["order"]], , drop = FALSE]),
    RoBMA_cluster_study   = study[rows[["order"]]],
    RoBMA_cluster_y       = rows[["response"]],
    RoBMA_cluster_se      = rows[["se"]],
    RoBMA_cluster_var     = rows[["variance"]],
    RoBMA_cluster_bin     = rows[["bin"]]
  )
  markers[[length(markers) + 1L]] <- list(
    name = "RoBMA_cluster_setup", distribution = "drobma_cluster_setup",
    arguments = list(group[["gate"]], "RoBMA_cluster_X", "RoBMA_cluster_study",
                     "RoBMA_cluster_y", "RoBMA_cluster_se", "RoBMA_cluster_var",
                     "RoBMA_cluster_bin", rows[["lower"]], rows[["upper"]],
                     rows[["omega"]], rows[["mode"]], rows[["telescope"]],
                     mu[["gate"]], mu[["label"]], group[["label"]], "tau", "rho",
                     tau[["gate"]], tau[["label"]], rows[["PET"]], rows[["PEESE"]],
                     n_studies, 5L)
  )
  list(markers = markers, data = cluster_data)
}


.sampler_cluster_moderator <- function(prior_list, name, width) {

  prior <- prior_list[[name]]
  if (!BayesTools::is.prior.mixture(prior) || BayesTools::is.prior.spike_and_slab(prior) ||
      length(prior) != 2L || !is.null(attr(prior, "multiply_by", exact = TRUE))) return(NULL)
  mixture <- BayesTools::JAGS_mixture_sources(prior_list, name)
  branches <- mixture[["branches"]]
  families <- vapply(branches, `[[`, character(1), "family")
  if (!setequal(families, c("point", "vector")) ||
      !all(is.finite(mixture[["probabilities"]]) & mixture[["probabilities"]] > 0)) return(NULL)
  slab_index <- which(families == "vector")
  slab <- branches[[slab_index]]
  spike <- branches[[which(families == "point")]]
  p <- slab[["prior"]]
  if (!length(spike[["location"]]) || !all(spike[["location"]] == 0) ||
      length(slab[["coordinates"]]) != width || !isTRUE(slab[["stochastic"]]) ||
      !is.null(attr(p, "multiply_by", exact = TRUE)) ||
      !isTRUE(p[["parameters"]][["sd"]] > 0) ||
      !isTRUE(p[["truncation"]][["lower"]] == -Inf) ||
      !isTRUE(p[["truncation"]][["upper"]] == Inf)) return(NULL)
  list(gate = mixture[["gate"]], source = slab[["source"]], label = slab_index)
}


.sampler_cluster_rho_source <- function(prior_list) {

  prior <- prior_list[["rho"]]
  source <- "rho"
  if (BayesTools::is.prior.mixture(prior) && !BayesTools::is.prior.spike_and_slab(prior) &&
      length(prior) == 1L) {
    branch <- BayesTools::JAGS_mixture_sources(prior_list, "rho")[["branches"]][[1L]]
    prior <- branch[["prior"]]
    source <- branch[["source"]]
  }
  if (!.sampler_joint_scalar_prior(prior) || !identical(prior[["distribution"]], "beta") ||
      !isTRUE(prior[["truncation"]][["lower"]] == 0) ||
      !isTRUE(prior[["truncation"]][["upper"]] == 1)) return(NULL)
  source
}


.sampler_cluster_rows <- function(data, prior_list, fit_data) {

  k <- fit_data[["K"]]
  bias <- prior_list[["bias"]]
  if (is.null(bias)) {
    if (.is_data_joint_selection(data) || length(fit_data[["yi"]]) != k ||
        length(fit_data[["sei"]]) != k || !all(fit_data[["sei"]] > 0)) return(NULL)
    return(list(order = seq_len(k), response = fit_data[["yi"]], se = fit_data[["sei"]],
                variance = fit_data[["sei"]]^2, bin = rep(1L, k), lower = -1, upper = 1,
                omega = 1, mode = 0, telescope = 0, PET = 0, PEESE = 0))
  }

  if (!inherits(bias, "prior.bias_mixture") || !.sampler_joint_independent_selection(data) ||
      !.selection_retains_other_random(data) || !isTRUE(fit_data[["sel_sign"]] == 1) ||
      !isTRUE(fit_data[["sel_joint_singleton_n"]] == k)) return(NULL)
  spec     <- .data_selection_execution_plan(data)[["selection_spec"]]
  branches <- BayesTools::JAGS_mixture_sources(prior_list, "bias")[["branches"]]
  families <- vapply(branches, `[[`, character(1), "family")
  if (!all(families %in% c("none", "PET", "PEESE", "selection")) ||
      !all(fit_data[["sel_joint_singleton_weights"]] == 1) ||
      !all(fit_data[["sel_joint_singleton_sei"]] > 0) ||
      !all(fit_data[["sel_joint_singleton_obs_bin"]] > 0)) return(NULL)

  list(order = fit_data[["sel_joint_singleton_row"]], response = fit_data[["sel_joint_singleton_y"]],
       se = fit_data[["sel_joint_singleton_sei"]], variance = fit_data[["sel_joint_singleton_sampling_variance"]],
       bin = fit_data[["sel_joint_singleton_obs_bin"]], lower = "sel_z_lower", upper = "sel_z_upper",
       omega = spec[["jags_omega"]], mode = if (isTRUE(spec[["jags_use_step_switch"]])) spec[["jags_kernel_mode"]] else 1,
       telescope = "sel_telescope_probabilities", PET = if ("PET" %in% families) "PET" else 0,
       PEESE = if ("PEESE" %in% families) "PEESE" else 0)
}
