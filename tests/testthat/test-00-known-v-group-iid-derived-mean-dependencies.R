test_that("allocated SD updates carry declared fixed-mean dependencies", {

  r <- function(name) get(name, envir = asNamespace("RoBMA"), inherits = FALSE)
  bt <- function(name) get(name, envir = asNamespace("BayesTools"), inherits = FALSE)
  helpers <- list(review_compile_object = review_compile_object,
    review_prior_columns = review_prior_columns,
    review_coefficient_name = review_coefficient_name,
    review_attach_replay = review_attach_replay)
  dat <- data.frame(yi = rep(0, 4L), x = c(1, 2, 1, 2),
    study = factor(c("a", "a", "b", "b")), estimate = factor(seq_len(4L)))
  gamma <- BayesTools::prior("gamma", list(shape = 2, rate = 1))
  constructor <- function(x_prior) RoBMA::brma.mv(yi = yi, V = diag(4), mods = ~0 + x,
    random = list(study = ~id(1 | study), estimate = ~id(1 | estimate)),
    data = dat, measure = "GEN", prior_unit_information_sd = 1,
    prior_mods = list(x = x_prior),
    prior_effect = BayesTools::prior("point", list(location = 0)),
    prior_heterogeneity = BayesTools::prior_random(allocation = list(
      BayesTools::random_variance_allocation(name = "split",
        terms = c(study = "study", estimate = "estimate"), target = "block",
        scale = "total_variance", sd = gamma,
        weights = BayesTools::prior("dirichlet", list(alpha = c(1, 1)))))),
    standardize_continuous_predictors = FALSE, only_priors = TRUE, silent = TRUE)
  discovery <- helpers$review_compile_object(constructor(gamma))
  discovery_terms <- discovery$compiled$mu$formula_design$random_effects
  grouped <- which(vapply(discovery_terms, function(term) {
    identical(term$block_name, "study") && anyDuplicated(term$group_map) > 0L
  }, logical(1)))
  stopifnot(length(grouped) == 1L)
  grouped_sd <- unique(discovery_terms[[grouped]]$sd_parameter_names)
  stopifnot(is.character(grouped_sd), length(grouped_sd) == 1L, !is.na(grouped_sd))

  for (dependent in c(TRUE, FALSE)) {
    x_prior <- gamma
    attr(x_prior, "multiply_by") <- if (dependent) grouped_sd else sqrt(.5)
    setup <- helpers$review_compile_object(constructor(x_prior))
    design <- setup$compiled$mu$formula_design
    terms <- design$random_effects
    stopifnot(length(terms) == 2L,
      setequal(vapply(terms, `[[`, character(1), "block_name"), c("study", "estimate")))
    term <- terms[[which(vapply(terms, function(x) {
      identical(x$block_name, "study")
    }, logical(1)))]]
    allocation <- term$sd_binding$allocations[[1L]]
    source <- allocation$source$name
    weight <- allocation$weight_name
    stopifnot(identical(allocation$source$kind, "prior"), isTRUE(allocation$source$owned),
      identical(allocation$source$shape, "scalar"), identical(allocation$scale, "total_variance"),
      identical(allocation$n_targets, 2L), identical(unique(term$sd_parameter_names), grouped_sd))
    priors <- c(setup$compiled$mu$prior_list,
      r(".create_fit_priors")(setup$object$data, setup$object$priors))
    stopifnot(!anyDuplicated(names(priors)))
    x_name <- helpers$review_coefficient_name(setup, "x")
    intercept <- helpers$review_coefficient_name(setup, "intercept")
    stopifnot(identical(priors[[source]]$distribution, "gamma"),
      priors[[source]]$parameters$shape == 2, priors[[source]]$parameters$rate == 1,
      identical(priors[[x_name]]$distribution, "gamma"),
      priors[[x_name]]$parameters$shape == 2, priors[[x_name]]$parameters$rate == 1,
      identical(priors[[weight]]$distribution, "dirichlet"),
      identical(as.numeric(priors[[weight]]$parameters$alpha), c(1, 1)))

    # Use actual compile modes/name maps: only the study block emits latents.
    name_map <- design$name_map
    additional <- unique(setup$compiled$mu$add_parameters)
    latent_bases <- character()
    latent_columns <- character()
    shapes <- list()
    for (random_term in terms) {
      compile_mode <- r(".random_effect_term_compile_mode")(random_term)
      base <- name_map$jags_name[name_map$kind == "random" &
        name_map$term == random_term$block_name & name_map$role == "_xRE_Zx"]
      stopifnot(compile_mode %in% c("sampled", "marginalized"),
        is.null(random_term$latent_layout), isTRUE(random_term$monitor$latent),
        random_term$n_columns == 1L, random_term$n_groups >= 1L,
        random_term$n_groups == as.integer(random_term$n_groups))
      columns <- character()
      if (identical(compile_mode, "sampled")) {
        stopifnot(length(base) == 1L, base %in% additional)
        index <- expand.grid(group = seq_len(random_term$n_groups),
          column = seq_len(random_term$n_columns))
        columns <- paste0(base, "[", index$group, ",", index$column, "]")
      } else {
        stopifnot(length(base) == 0L,
          !paste0(random_term$parameter_stem, "_xRE_Zx") %in% additional)
      }
      shapes[[random_term$block_name]] <- list(compile_mode = compile_mode,
        n_groups = random_term$n_groups, n_columns = random_term$n_columns)
      latent_bases <- c(latent_bases, base)
      latent_columns <- c(latent_columns, columns)
    }
    stopifnot(identical(shapes$study$compile_mode, "sampled"),
      identical(shapes$estimate$compile_mode, "marginalized"),
      setequal(latent_bases, grep("_xRE_Zx$", additional, value = TRUE)),
      !anyDuplicated(latent_columns),
      length(latent_columns) == sum(vapply(shapes, function(shape) {
        if (identical(shape$compile_mode, "sampled")) {
          shape$n_groups * shape$n_columns
        } else 0L
      }, integer(1))))
    nodes <- bt(".bt_deterministic_nodes")(prior_list = priors,
      formula_design = lapply(setup$compiled, `[[`, "formula_design"))
    sd_nodes <- Filter(function(node) identical(node$family, "random_sd"), nodes)
    sd_names <- vapply(sd_nodes, `[[`, character(1), "node")
    stopifnot(length(sd_nodes) == 2L, grouped_sd %in% sd_names,
      setequal(sd_names, unlist(lapply(terms, `[[`, "sd_parameter_names"), use.names = FALSE)))
    prior_columns <- helpers$review_prior_columns(priors)
    prior_bases <- unique(sub("\\[.*$", "", prior_columns))
    stopifnot(length(setdiff(additional, c(latent_bases, sd_names, prior_bases))) == 0L)
    weights <- paste0(weight, "[", seq_len(allocation$n_targets), "]")
    auxiliary <- r(".iwmde_simplex_auxiliary_columns")(weight, allocation$n_targets)
    columns <- unique(c(prior_columns, auxiliary, latent_columns, sd_names))
    samples <- matrix(.6, 24L, length(columns), dimnames = list(NULL, columns))
    samples[, intercept] <- 0
    samples[, x_name] <- c(1, seq(.65, 1.75, length.out = 23L))
    samples[, source] <- c(1, seq(.7, 1.8, length.out = 23L))
    samples[, weights] <- .5
    samples[, auxiliary] <- 1
    lookup <- bt(".bt_deterministic_lookup")(samples, prior_list = priors)
    for (node in sd_nodes) {
      derived <- bt(".bt_deterministic_node_evaluate")(node, lookup)
      stopifnot(is.matrix(derived), nrow(derived) == 24L,
        identical(colnames(derived), node$coordinates),
        all(is.finite(derived)), all(derived >= 0))
      samples[, node$coordinates] <- derived
    }
    required <- unique(c(BayesTools::JAGS_to_monitor(priors), additional))
    supplied <- unique(sub("\\[.*$", "", colnames(samples)))
    stopifnot(all(required %in% supplied), all(is.finite(samples)),
      all(samples[, source] > 0), all(samples[, x_name] > 0),
      all(samples[, weights] > 0 & samples[, weights] < 1),
      all(rowSums(samples[, weights, drop = FALSE]) == 1),
      all(samples[, latent_columns, drop = FALSE] == .6))
    object <- helpers$review_attach_replay(setup, samples)
    BayesTools::JAGS_validate_fit_contract(object$fit, requires = c("name_encoding",
      "formula_name_map", "formula_design", "parameter_map", "draw_geometry"))
    stopifnot(identical(class(object), class(setup$object)))
    selection <- BayesTools::parameter_catalog_resolve(
      BayesTools::parameter_catalog(object$fit), "(mu) split: sd_total", "mu")
    stopifnot(nrow(selection$quantities) == 1L, !selection$quantities$internal)
    selector <- r(".brma_random_parameter_io_labels")(selection$quantities, "selector")
    target <- r(".brma_random_parameter_density_target")(object, selector)
    stopifnot(is.null(target$reason), identical(target$parameter, source))
    context <- r(".iwmde_context")(object)
    spec <- r(".iwmde_parameter_spec")(context, target$parameter, target$parameter_spec)
    stopifnot(identical(spec$status, "ok"), identical(spec$type, "primitive"))
    states <- r(".iwmde_row_states")(context, 1L,
      parameter = target$parameter, parameter_spec = spec)
    replacement <- r(".iwmde_replacement_spec")(context, target$parameter, spec)
    active <- states[[1L]]$active_setup
    plan <- r(".iwmde_known_v_random_group_iid_plan")(context, target$parameter,
      replacement, priors = active$priors)
    stopifnot(length(states) == 1L, identical(states[[1L]]$state_scope, "global"),
      identical(states[[1L]]$likelihood_mode, "marginal"),
      identical(replacement$type, "scalar"), identical(replacement$name, source),
      identical(plan$target_mode, "source_sd"), identical(plan$source_parameter, source))

    values <- c(1, 2)
    args <- list(context = context, parameter = target$parameter, values = values,
      row_states = states, replacement = replacement)
    candidates <- do.call(r(".iwmde_build_replacement_samples"), args)
    stopifnot(all(candidates$valid), nrow(candidates$samples) == 2L,
      identical(as.numeric(candidates$samples[, source]), values))
    expect_setequal(candidates$changed_parameters, c(source, sd_names))
    # The first finite candidate equals the source row, but declared names remain.
    for (value in values) {
      candidate <- r(".iwmde_build_replacement_samples")(
        context, target$parameter, value, states, replacement)
      expect_setequal(candidate$changed_parameters, c(source, sd_names))
      dependencies <- BayesTools::JAGS_formula_coordinate_dependencies(
        object$fit, candidate$changed_parameters)
      # The public query includes the conditional random predictor separately
      # from the fixed-mean multiplier on the same changed SD coordinate.
      expect_true(any(dependencies$formula_parameter == "mu" &
        dependencies$dependency_type == "random_effect"))
      expect_identical(any(dependencies$formula_parameter == "mu" &
        dependencies$dependency_type == "multiplier"), dependent)
      expect_identical(r(".iwmde_q_grid_fixed_mu_is_invariant")(
        context, candidate), !dependent)
    }
    expect_true(all(is.finite(candidates$samples[, sd_names, drop = FALSE])))
    expect_true(max(abs(unname(candidates$samples[, sd_names, drop = FALSE]) -
      outer(values, rep(sqrt(.5), length(sd_names))))) <= 1e-12)
    state_priors <- states[[1L]]$prior_list
    expect_setequal(names(state_priors), c(source, x_name, intercept, weight))
    nuisance_columns <- setdiff(colnames(candidates$samples), c(source, sd_names))
    for (column in nuisance_columns) {
      expect_identical(as.numeric(candidates$samples[, column]),
        rep(as.numeric(samples[1L, column]), 2L))
    }
    nuisance_priors <- state_priors[setdiff(names(state_priors), source)]
    nuisance_logs <- vapply(names(nuisance_priors), function(owner) {
      r(".iwmde_log_prior_rows")(candidates$samples,
        stats::setNames(list(nuisance_priors[[owner]]), owner))
    }, numeric(2L))
    expect_true(all(is.finite(nuisance_logs)))
    expect_true(all(nuisance_logs[1L, ] == nuisance_logs[2L, ]))
    full_prior <- r(".iwmde_log_prior_rows")(candidates$samples, state_priors)
    independent_prior_delta <- dgamma(values, 2, 1, log = TRUE) -
      dgamma(values[[1L]], 2, 1, log = TRUE)
    expect_true(all(is.finite(full_prior)))
    expect_true(max(abs(full_prior - full_prior[[1L]] - independent_prior_delta)) <= 1e-12)

    beta <- samples[1L, x_name]
    p <- samples[1L, weights[[allocation$index]]]
    stopifnot(beta == 1, samples[1L, source] == 1, p == .5)
    Z <- cbind(a = as.integer(dat$study == "a"), b = as.integer(dat$study == "b"))
    oracle_mu <- t(vapply(values, function(value) {
      dat$x * beta * if (dependent) value * sqrt(p) else sqrt(.5)
    }, numeric(4L)))
    actual_mu <- r(".iwmde_predictor_evaluate_fixed_mu")(context, active, candidates$samples)
    expect_identical(as.numeric(actual_mu), as.numeric(oracle_mu))
    oracle_log_lik <- vapply(seq_along(values), function(i) {
      covariance <- diag(4) + values[[i]]^2 * (p * tcrossprod(Z) + (1 - p) * diag(4))
      factor <- chol(covariance)
      solved <- forwardsolve(t(factor), dat$yi - oracle_mu[i, ])
      -.5 * (4 * log(2 * pi) + 2 * sum(log(diag(factor))) + sum(solved^2))
    }, numeric(1))
    oracle_delta <- oracle_log_lik - oracle_log_lik[[1L]] + independent_prior_delta
    routes <- list(
      dispatcher = do.call(r(".iwmde_log_q_grid"), args),
      generic = do.call(r(".iwmde_log_q_grid_batch"), args),
      scalar = do.call(r(".iwmde_log_q_grid_scalar"), args))
    direct <- do.call(r(".iwmde_log_q_grid_known_v_random_group_iid"),
      c(args, list(active_setup = active)))
    if (dependent) {
      expect_null(direct)
    } else {
      expect_true(is.matrix(direct))
      routes$direct <- direct
    }
    # Four well-conditioned observations, with no integration or Monte Carlo
    # error. Preserve the independently checked 190 log-difference allowance.
    for (route in routes) {
      expect_true(is.matrix(route))
      expect_identical(dim(route), c(2L, 1L))
      expect_true(all(is.finite(route)))
      delta <- as.numeric(route) - as.numeric(route)[[1L]]
      expect_true(max(abs(delta - oracle_delta)) <= 1e-12)
    }

    if (dependent) {
      # A finite negative source produces invalid derived SDs and no SD write.
      # Names describe the selected declaration even when no row can be used.
      invalid <- r(".iwmde_build_replacement_samples")(
        context, target$parameter, -1, states, replacement)
      expect_identical(invalid$valid, FALSE)
      expect_setequal(invalid$changed_parameters, c(source, sd_names))
      expect_identical(as.numeric(invalid$samples[, sd_names]),
        as.numeric(samples[1L, sd_names]))
    }
  }
})
