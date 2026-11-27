# Qualify data using the same descriptors that choose the fitted likelihood.
# Prior and design eligibility is checked after BayesTools binds the formulas.
.sampler_joint_declared <- function(data, fit_priors, policy) {

  if (!.sampler_module_samplers()) return(NULL)
  if (policy == "off") return(NULL)
  if (policy %in% c("clustered", "clustered.scale")) {
    checked_data <- data
    return(function(prior_list, formula_design, data) {
      .sampler_joint_clustered(checked_data, prior_list, formula_design, data,
                               scale = policy == "clustered.scale")
    })
  }
  location <- function(moderator = FALSE) {
    probabilities <- RoBMA.get_option("sampler.probabilities")
    structure(function(prior_list, formula_design, data) {
      .sampler_joint_location(prior_list, formula_design, data, probabilities, moderator)
    }, sampler = "RoBMA::LocationBlock")
  }
  if (policy %in% c("auto", "location") && .sampler_joint_location_data(data, policy)) return(location())
  if (policy == "location") return(NULL)
  independent <- .data_outcome_type(data) == "norm" && !.is_data_multilevel(data) && !.is_data_random(data)
  if (!independent) return(NULL)

  if (!.is_data_scale(data) && !.is_data_weights(data)) {
    sign <- if (.data_effect_direction(data) == "negative") -1 else 1
    if (policy %in% c("auto", "gaussian") && !.is_data_known_v(data) &&
        .is_data_mods(data) && is.null(fit_priors[["bias"]])) {
      auto_dispatch <- if (policy == "auto") new.env(parent = emptyenv()) else NULL
      # "auto": a Gaussian candidate that its bound design declines falls back
      # to the location block (resolved once the design is bound; the
      # dispatch records the sampler taken). The row budget is a deliberate
      # exception: a fit over it keeps the ordinary samplers.
      fallback <- if (policy == "auto") location(TRUE) else NULL
      return(structure(function(prior_list, formula_design, data) {
        markers <- .sampler_joint_gaussian(prior_list, formula_design, data, sign, auto_dispatch)
        if (is.null(markers) && !is.null(fallback) && is.null(auto_dispatch[["budget"]])) {
          markers <- fallback(prior_list, formula_design, data)
          if (!is.null(markers)) auto_dispatch[["sampler"]] <- attr(fallback, "sampler", exact = TRUE)
        }
        markers
      }, sampler = "RoBMA::GaussianJoint", auto_dispatch = auto_dispatch,
      fallback = attr(fallback, "sampler", exact = TRUE)))
    }
    if (policy %in% c("auto", "bias") &&
        (!.is_data_mods(data) || (policy == "bias" && !.is_data_known_v(data))) &&
        inherits(fit_priors[["bias"]], "prior.bias_mixture") &&
        .sampler_joint_independent_selection(data)) {
      return(structure(function(prior_list, formula_design, data) {
        .sampler_joint_bias(prior_list, formula_design, data, sign)
      }, sampler = "RoBMA::IndependentBiasJoint"))
    }
  }
  # Independent fits that no other joint sampler takes: the location block
  # (without latents) when the fit has a gated moderator (decided by
  # measurement on the Assink meta-regression).
  if (policy == "auto" && .is_data_mods(data)) return(location(TRUE))
  return(NULL)
}


.sampler_joint_independent_selection <- function(data) {

  if (!.is_data_joint_selection(data) || .selection_retains_sampling(data) ||
      .selection_retains_estimate(data) || .data_known_v_correlated(data) ||
      .data_known_v_rank(data) > 0L) return(FALSE)
  model <- .data_selection_model(data)
  plan  <- .data_selection_execution_plan(data)
  rows  <- seq_len(nrow(data[["outcome"]]))

  return(length(plan[["dependent_blocks"]]) == 0L &&
           setequal(plan[["singleton_rows"]], rows) &&
           is.null(plan[["random_covariance"]]) &&
           identical(plan[["selection_spec"]][["mode"]], "step") &&
           all(vapply(model[["branches"]][model[["active_branches"]]], function(branch) {
             identical(branch[["weight_rule"]], "product")
           }, logical(1))))
}


# Only literal independent scalar laws can supply the proposal coordinates.
.sampler_joint_scalar_prior <- function(prior) {

  return(.sampler_prior_scalar(prior) &&
           !BayesTools::is.prior.R2(prior) && !BayesTools::is.prior.regularized(prior) &&
           !BayesTools::is.prior.regularized_term(prior) &&
           is.null(attr(prior, "multiply_by", exact = TRUE)) &&
           all(vapply(prior[["parameters"]], function(value) {
             is.numeric(value) && length(value) == 1L && is.finite(value)
           }, logical(1))))
}


# Resolve a two-label mixture containing an exact zero and a normal slab.
# The physical spike/slab ordering is the one emitted by the compiler.
.sampler_joint_normal_mixture <- function(prior_list, name, tau = FALSE) {

  prior <- prior_list[[name]]
  if (!BayesTools::is.prior.mixture(prior) || BayesTools::is.prior.spike_and_slab(prior) ||
      length(prior) != 2L || !is.null(attr(prior, "multiply_by", exact = TRUE))) return(NULL)

  mixture <- BayesTools::JAGS_mixture_sources(prior_list, name)
  branches <- mixture[["branches"]]
  family <- vapply(branches, `[[`, character(1), "family")
  if (!setequal(family, c("point", "scalar")) ||
      !all(is.finite(mixture[["probabilities"]]) & mixture[["probabilities"]] > 0)) return(NULL)

  spike_index <- which(family == "point")
  slab_index  <- which(family == "scalar")
  spike       <- branches[[spike_index]]
  slab        <- branches[[slab_index]]
  slab_prior  <- slab[["prior"]]
  if (!is.numeric(spike[["location"]]) || !isTRUE(spike[["location"]] == 0) ||
      !.sampler_joint_scalar_prior(spike[["prior"]]) ||
      !isTRUE(slab[["stochastic"]]) || !.sampler_joint_scalar_prior(slab_prior) ||
      !identical(slab_prior[["distribution"]], "normal")) return(NULL)

  parameters <- slab_prior[["parameters"]]
  truncation <- slab_prior[["truncation"]]
  sd <- parameters[["sd"]]
  precision <- 1 / sd^2
  lower <- if (tau) 0 else -Inf
  if (!isTRUE(sd > 0 && is.finite(precision) && precision > 0 &&
              truncation[["upper"]] == Inf) ||
      !isTRUE(truncation[["lower"]] == lower) ||
      (tau && !isTRUE(parameters[["mean"]] == 0))) return(NULL)

  return(list(gate = mixture[["gate"]], source = slab[["source"]], label = slab_index))
}


.sampler_joint_bias <- function(prior_list, formula_design, data, sign) {

  location <- "mu"
  intercept <- "mu"
  moderators <- NULL
  if (length(formula_design) > 0L) {
    moderators <- .sampler_joint_bias_moderators(prior_list, formula_design, data)
    if (is.null(moderators)) return(NULL)
    location  <- moderators[["prior_names"]]
    intercept <- moderators[["intercept"]]
  }
  prior_names <- names(prior_list)
  if (!all(c(location, "bias") %in% prior_names) ||
      !all(prior_names %in% c(location, "bias", "tau")) || length(prior_list[["bias"]]) < 2L) return(NULL)
  mu <- .sampler_joint_normal_mixture(prior_list, intercept)
  if (is.null(mu)) return(NULL)

  # Heterogeneity is left in place. Its source law must be independent of
  # the joint coordinates so the row variances are fixed during the proposal.
  tau <- prior_list[["tau"]]
  if (!is.null(tau)) {
    tau_components <- if (BayesTools::is.prior.mixture(tau) &&
                          !BayesTools::is.prior.spike_and_slab(tau)) tau else list(tau)
    if (!all(vapply(tau_components, .sampler_joint_scalar_prior, logical(1)))) return(NULL)
  }

  bias <- BayesTools::JAGS_mixture_sources(prior_list, "bias")
  if (!is.null(moderators) && isTRUE(bias[["nested_weights"]])) return(NULL)
  branches <- bias[["branches"]]
  family_codes <- c(none = 0, selection = 1, PET = 2, PEESE = 3)
  codes <- unname(family_codes[vapply(branches, `[[`, character(1), "family")])
  if (anyNA(codes) || !any(codes == 1) || !any(codes != 1) ||
      !all(is.finite(bias[["probabilities"]]) & bias[["probabilities"]] > 0)) return(NULL)

  n_bins <- length(data[["sel_z_lower"]])
  singleton_rows <- data[["sel_joint_singleton_row"]]
  singleton_se   <- data[["sel_joint_singleton_sei"]]
  n <- data[["K"]]
  if (n_bins < 2L || !isTRUE(data[["sel_joint_singleton_n"]] == n) ||
      length(data[["sel_joint_singleton_y"]]) != n || length(singleton_se) != n ||
      !all(singleton_se > 0) ||
      !all(data[["sel_joint_singleton_weights"]] == 1)) return(NULL)
  if (!all(lengths(lapply(branches[codes == 1], `[[`, "omega")) == n_bins)) return(NULL)

  slope_branches <- which(codes %in% c(2, 3))
  if (length(slope_branches) > 0L) {
    if (!identical(unname(data[["sei"]][singleton_rows]), unname(singleton_se))) return(NULL)
    exempt <- data[["bias_exempt"]]
    if (!is.null(exempt) && !identical(data[["sel_joint_singleton_obs_bin"]] == 0L,
                                      exempt[singleton_rows] == 1L)) return(NULL)
  }
  if (!all(vapply(branches[slope_branches], function(branch) {
    .sampler_joint_scalar_prior(branch[["prior"]]) && !is.na(branch[["source"]])
  }, logical(1)))) return(NULL)

  zero   <- "RoBMA_joint_zero"
  one    <- "RoBMA_joint_one"
  slopes <- rep(zero, length(branches))
  slopes[slope_branches] <- vapply(branches[slope_branches], `[[`, character(1), "source")
  omega  <- lapply(branches, function(branch) {
    if (branch[["family"]] == "selection") branch[["omega"]] else rep(one, n_bins)
  })

  tag <- "robma_joint_bias"
  arguments <- list(bias[["gate"]], mu[["gate"]], mu[["source"]], mu[["label"]],
                    sign, codes, slopes, unlist(omega, use.names = FALSE))
  marker_data <- stats::setNames(list(0, 1), c(zero, one))
  if (!is.null(moderators)) {
    tag <- "robma_joint_bias_moderated"
    arguments <- c(arguments, list("sel_joint_singleton_y", moderators[["coefficients"]],
                                   "RoBMA_joint_moderator_design"))
    marker_data[["RoBMA_joint_moderator_design"]] <- moderators[["design"]]
  }
  return(list(markers = list(list(name = "RoBMA_joint_bias", tag = tag,
                                 arguments = arguments)), data = marker_data))
}


.sampler_joint_gaussian_auto_max_rows <- 800L

.sampler_joint_gaussian <- function(prior_list, formula_design, data, sign, auto_dispatch = NULL) {

  if (!identical(names(formula_design), "mu")) return(NULL)
  design <- formula_design[["mu"]]
  terms  <- design[["model_terms_type"]]
  slope  <- setdiff(names(terms), "intercept")
  if (length(terms) != 2L || !all(terms == "continuous") || length(slope) != 1L ||
      isTRUE(design[["log_intercept"]]) || length(design[["expression_specs"]]) > 0L ||
      length(design[["random_effects"]]) > 0L || ncol(design[["model_matrix"]]) != 2L) return(NULL)

  map <- design[["name_map"]]
  coefficients <- map[map[["kind"]] == "fixed" & map[["role"]] == "coefficient", , drop = FALSE]
  intercept_prior <- coefficients[["jags_name"]][coefficients[["term"]] == "intercept"]
  slope_prior <- coefficients[["jags_name"]][coefficients[["term"]] != "intercept"]
  slope_data  <- design[["jags_data_names"]][[slope]]
  if (length(intercept_prior) != 1L || length(slope_prior) != 1L ||
      length(slope_data) != 1L || !setequal(names(prior_list), c(intercept_prior, slope_prior, "tau")) ||
      !is.numeric(data[[slope_data]]) || !is.null(dim(data[[slope_data]])) ||
      length(data[[slope_data]]) != data[["K"]] || !all(data[["sei"]] > 0)) return(NULL)
  automatic <- !is.null(auto_dispatch)
  if (automatic && !slope_prior %in% names(design[["formula_scale"]])) return(NULL)

  intercept <- .sampler_joint_normal_mixture(prior_list, intercept_prior)
  moderator <- .sampler_joint_normal_mixture(prior_list, slope_prior)
  tau       <- .sampler_joint_normal_mixture(prior_list, "tau", tau = TRUE)
  if (is.null(intercept) || is.null(moderator) || is.null(tau)) return(NULL)

  max_rows <- .sampler_joint_gaussian_auto_max_rows
  if (automatic && data[["K"]] > max_rows) {
    auto_dispatch[["budget"]] <- list(rows = data[["K"]], limit = max_rows)
    return(NULL)
  }

  return(list(markers = list(list(
    name         = "RoBMA_joint_gaussian",
    tag          = "robma_gaussian_joint",
    arguments    = list(intercept[["gate"]], intercept[["source"]],
                        moderator[["gate"]], moderator[["source"]], tau[["gate"]], tau[["source"]],
                        c(intercept[["label"]], moderator[["label"]], tau[["label"]]),
                        slope_data, "yi", "sei", sign)
  ))))
}
