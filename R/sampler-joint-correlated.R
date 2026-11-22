# The correlated proposal conditions on the retained study contribution and
# fixed moderator coefficients. The fitted selection graph is left intact.
.sampler_joint_correlated <- function(data, prior_list, formula_design, fit_data) {

  if (.data_outcome_type(data) != "norm" || !.is_data_known_v(data) ||
      !.is_data_random(data) || .is_data_scale(data) || .is_data_weights(data) ||
      !.is_data_joint_selection(data) || .selection_retains_sampling(data) ||
      .data_effect_direction(data) != "positive" ||
      any(data[["outcome"]][["bias_exempt"]]) ||
      !inherits(prior_list[["bias"]], "prior.bias_mixture") ||
      !identical(names(formula_design), "mu")) return(NULL)

  model <- .data_selection_model(data)
  plan  <- .data_selection_execution_plan(data)
  if (!identical(plan[["selection_spec"]][["mode"]], "step") ||
      !length(plan[["dependent_blocks"]]) ||
      !all(vapply(model[["branches"]][model[["active_branches"]]], function(branch) {
        identical(branch[["weight_rule"]], "product")
      }, logical(1)))) return(NULL)

  k <- fit_data[["K"]]
  design <- formula_design[["mu"]]
  random <- .sampler_correlated_random_sources(design, plan, k)
  fixed  <- .sampler_correlated_fixed_sources(design, prior_list, fit_data, k)
  if (is.null(random) || is.null(fixed)) return(NULL)

  bias <- BayesTools::JAGS_mixture_sources(prior_list, "bias")
  branches <- bias[["branches"]]
  family_codes <- c(none = 0L, selection = 1L, PET = 2L, PEESE = 3L)
  codes <- unname(family_codes[vapply(branches, `[[`, character(1), "family")])
  if (anyNA(codes) || !any(codes == 1L) || !any(codes != 1L) ||
      !all(is.finite(bias[["probabilities"]]) & bias[["probabilities"]] > 0)) return(NULL)

  n_bins <- length(fit_data[["sel_z_lower"]])
  if (n_bins < 2L || length(fit_data[["sel_z_upper"]]) != n_bins ||
      !all(lengths(lapply(branches[codes == 1L], `[[`, "omega")) == n_bins)) return(NULL)
  slope_branches <- which(codes %in% c(2L, 3L))
  if (!all(vapply(branches[slope_branches], function(branch) {
    .sampler_joint_scalar_prior(branch[["prior"]]) && !is.na(branch[["source"]])
  }, logical(1)))) return(NULL)

  zero   <- "RoBMA_correlated_zero"
  one    <- "RoBMA_correlated_one"
  slopes <- rep(zero, length(branches))
  slopes[slope_branches] <- vapply(branches[slope_branches], `[[`, character(1), "source")
  omega <- lapply(branches, function(branch) {
    if (branch[["family"]] == "selection") branch[["omega"]] else rep(one, n_bins)
  })

  V <- .known_v_covariance_matrix(.data_known_v_data(data))
  if (!identical(dim(V), c(k, k)) || !any(V[lower.tri(V)] != 0)) return(NULL)
  covariance <- eigen(V, symmetric = TRUE)
  if (!all(is.finite(covariance[["values"]]) & covariance[["values"]] > 0)) return(NULL)

  payload <- list(RoBMA_correlated_zero = 0, RoBMA_correlated_one = 1,
    RoBMA_correlated_rows = c(data[["outcome"]][["yi"]], data[["outcome"]][["sei"]], diag(V)),
    RoBMA_correlated_eigenvalues = covariance[["values"]],
    RoBMA_correlated_basis = as.double(covariance[["vectors"]]),
    RoBMA_correlated_design = as.double(fixed[["design"]]))

  mu <- fixed[["mu"]]
  return(list(markers = list(list(
    name         = "RoBMA_joint_correlated",
    distribution = "drobma_joint_correlated",
    arguments    = list(bias[["gate"]], mu[["source"]], mu[["gate"]], mu[["label"]],
      random[["sd"]], random[["retained"]], fixed[["moderator"]], "RoBMA_correlated_design",
      codes, slopes, unlist(omega, use.names = FALSE), "RoBMA_correlated_rows",
      "RoBMA_correlated_eigenvalues", "RoBMA_correlated_basis", "sel_z_lower", "sel_z_upper")
  )), data = payload))
}


.sampler_correlated_random_sources <- function(design, plan, k) {

  random <- design[["random_effects"]]
  blocks <- vapply(random, .random_effect_term_block_name, character(1))
  integrated <- plan[["random_covariance"]][["term_names"]]
  if (length(random) != 2L || length(integrated) != 1L || !integrated %in% blocks ||
      anyDuplicated(blocks)) return(NULL)

  estimate <- random[[match(integrated, blocks)]]
  retained <- random[[which(blocks != integrated)]]
  unit_intercept <- function(term) {
    X <- term[["model_matrix"]]
    identical(dim(X), c(k, 1L)) && all(X == 1) &&
      identical(term[["model_terms"]], "intercept") &&
      !.random_effect_term_has_known_group_covariance(term)
  }
  if (!unit_intercept(estimate) || !unit_intercept(retained) ||
      !identical(estimate[["compile_mode"]], "marginalized") ||
      !identical(retained[["compile_mode"]], "sampled") ||
      length(unique(estimate[["group_map"]])) != k ||
      !identical(estimate[["structure"]], "us") ||
      length(estimate[["sd_parameter_names"]]) != 1L ||
      length(retained[["parameter_stem"]]) != 1L) return(NULL)

  return(list(sd = estimate[["sd_parameter_names"]], retained = retained[["parameter_stem"]]))
}


.sampler_correlated_fixed_sources <- function(design, prior_list, fit_data, k) {

  if (isTRUE(design[["log_intercept"]]) || length(design[["expression_specs"]])) return(NULL)
  map <- design[["name_map"]]
  fixed <- map[map[["kind"]] == "fixed" & map[["role"]] == "coefficient", , drop = FALSE]
  intercept <- fixed[["jags_name"]][fixed[["term"]] == "intercept"]
  moderator <- fixed[fixed[["term"]] != "intercept", , drop = FALSE]
  if (length(intercept) != 1L || nrow(moderator) > 1L) return(NULL)
  mu <- .sampler_joint_normal_mixture(prior_list, intercept)
  if (is.null(mu)) return(NULL)
  if (!nrow(moderator)) {
    return(list(mu = mu, moderator = "RoBMA_correlated_zero", design = numeric(k)))
  }

  term <- moderator[["term"]]
  matrix_name <- design[["jags_data_names"]][[term]]
  if (length(matrix_name) != 1L) return(NULL)
  X <- fit_data[[matrix_name]]
  if (is.numeric(X) && is.null(dim(X))) X <- matrix(X, ncol = 1L)
  if (!is.matrix(X) || !is.numeric(X) || nrow(X) != k || ncol(X) < 1L ||
      !all(is.finite(X)) || !design[["model_terms_type"]][[term]] %in%
        c("continuous", "factor")) return(NULL)

  return(list(mu = mu, moderator = moderator[["jags_name"]], design = X))
}
