# Synthetic replay of constructor objects for review regressions: the genuine
# formula designs and priors of an 'only_priors' object are compiled by
# BayesTools and finalized as a BayesTools fit around supplied posterior
# draws, so that post-fit methods run without sampling.

review_heterogeneity_har_replay <- function(scaled = FALSE, n = 120L,
    weights = c(.5, .5), rho = 0, time = c(1, 2, 1, 2), direct = FALSE,
    component_name = NULL) {

  dat <- data.frame(yi = seq(.1, .4, length.out = length(time)),
    study = factor(rep(c("s1", "s2"), each = length(time) / 2L)), time = time)
  args <- list(yi = quote(yi), V = diag(.04, nrow(dat)),
    random = ~ har(time | study), data = dat, measure = "GEN",
    prior_unit_information_sd = 1, only_priors = TRUE, silent = TRUE,
    standardize_continuous_predictors = FALSE,
    prior_effect = BayesTools::prior("point", list(location = 0)),
    prior_heterogeneity = BayesTools::prior("normal", list(mean = 0, sd = 1),
      truncation = list(lower = 0, upper = Inf)))
  if (scaled) args[["scale"]] <- ~ time
  if (!is.null(component_name)) {
    args[["random"]] <- stats::setNames(list(~ har(time | study)), component_name)
  }
  if (direct) args[["prior_heterogeneity"]] <- BayesTools::prior_random(
    sd = BayesTools::prior("normal", list(mean = 0, sd = 1),
      truncation = list(lower = 0, upper = Inf)))
  setup <- review_compile_object(do.call(brma.mv, args))
  term <- setup$compiled$mu$formula_design$random_effects[[1L]]
  allocation <- if (direct) NULL else term$sd_binding$allocations[[1L]]
  priors <- list()
  for (compiled in setup$compiled) priors <- c(priors, compiled$prior_list)
  priors <- c(priors, .create_fit_priors(setup$object$data, setup$object$priors))
  weight_columns <- if (direct) character() else paste0(allocation$weight_name, "[", 1:2, "]")
  auxiliary <- if (direct) character() else .iwmde_simplex_auxiliary_columns(allocation$weight_name, 2L)
  additional <- unique(unlist(lapply(setup$compiled, `[[`, "add_parameters")))
  latent_base <- paste0(term$parameter_stem, "_xRE_Zx")
  index <- expand.grid(group = seq_len(term$n_groups), column = seq_len(term$n_columns))
  latent <- paste0(latent_base, "[", index$group, ",", index$column, "]")
  columns <- unique(c(review_prior_columns(priors), auxiliary,
    setdiff(additional, latent_base), latent))
  samples <- matrix(0, n, length(columns), dimnames = list(NULL, columns))
  if (!direct) {
    samples[, weight_columns] <- rep(weights, each = n)
    samples[, auxiliary] <- rep(2 * weights, each = n)
  }
  samples[, term$correlation$sample_name] <- if (identical(term$correlation$rho_scale, "fisher_z")) atanh(rho) else rho
  if (scaled) {
    samples[, review_coefficient_name(setup, "intercept", "log_tau")] <- 2
    samples[, review_coefficient_name(setup, "time", "log_tau")] <- log(2)
    output <- setup$compiled$log_tau$formula_design$name_map
    output <- output$jags_name[output$kind == "formula_output"]
    generated <- matrix(rep(c(log(2) + time * log(2), 2 * 2^time), each = n), nrow = n,
      dimnames = list(NULL, c(paste0(output, "[", seq_along(time), "]"),
        paste0(allocation$source$name, "[", seq_along(time), "]"))))
    samples <- cbind(samples, generated)
  } else {
    if (direct) {
      samples[, term$sd_parameter_names] <- rep(2 * sqrt(2 * weights), each = n)
    } else {
      samples[, allocation$source$name] <- 2
      samples[, term$sd_binding$sd_component_names] <- rep(2 * sqrt(2 * weights), each = n)
    }
  }
  object <- review_attach_replay(setup, samples)
  BayesTools::JAGS_validate_fit_contract(object$fit,
    requires = c("name_encoding", "formula_name_map", "formula_design", "parameter_map", "draw_geometry"))
  list(object = object, samples = samples, setup = setup, term = term, allocation = allocation)
}

review_heterogeneity_parent_replay <- function(role) {

  dat <- data.frame(yi = c(.1, .2, .3, .4), time = c(1, 2, 1, 2),
    study = factor(c("s1", "s1", "s2", "s2")), other = factor(1:4))
  args <- list(yi = quote(yi), V = diag(.04, 4L),
    random = list(combined = ~har(time | study) + id(1 | other)),
    scale = list(combined = ~time), data = dat, prior_unit_information_sd = 1,
    prior_effect = BayesTools::prior("point", list(location = 0)), measure = "GEN",
    standardize_continuous_predictors = FALSE, only_priors = TRUE, silent = TRUE)
  discovery <- do.call(brma.mv, args)
  blocks <- vapply(discovery$formula_design$mu$random_effects, `[[`, character(1), "block_name")
  source <- .data_scale_component_specs(discovery$data)[[1L]]$source
  args$prior_heterogeneity <- BayesTools::prior_random(allocation = list(
    BayesTools::random_variance_allocation(name = "root",
      terms = c(study = blocks[1L], other = blocks[2L]),
      sd_source = BayesTools::random_sd_source(source, shape = "row"),
      weights = BayesTools::prior("dirichlet", list(alpha = c(1, 1))),
      inclusion = list(study = BayesTools::prior("spike", list(location = .5)))),
    BayesTools::random_variance_allocation(name = "child", terms = blocks[1L],
      target = "sd_component", parent = BayesTools::allocation_ref("root", "study"),
      scale = role, weights = BayesTools::prior("dirichlet", list(alpha = c(1, 1))))))
  setup <- review_compile_object(do.call(brma.mv, args))
  design <- setup$compiled$mu$formula_design
  priors <- list()
  for (compiled in setup$compiled) priors <- c(priors, compiled$prior_list)
  priors <- c(priors, .create_fit_priors(setup$object$data, setup$object$priors))
  additional <- unique(unlist(lapply(setup$compiled, `[[`, "add_parameters")))
  latent <- unlist(lapply(design$random_effects, function(term) {
    base <- paste0(term$parameter_stem, "_xRE_Zx")
    index <- expand.grid(group = seq_len(term$n_groups), column = seq_len(term$n_columns))
    paste0(base, "[", index$group, ",", index$column, "]")
  }))
  bases <- vapply(design$random_effects, function(term) paste0(term$parameter_stem, "_xRE_Zx"), character(1))
  auxiliary <- unlist(lapply(design$random_allocations, function(allocation) {
    .iwmde_simplex_auxiliary_columns(allocation$weight_name, 2L)
  }))
  columns <- unique(c(review_prior_columns(priors), setdiff(additional, bases), auxiliary, latent))
  samples <- matrix(0, 4L, length(columns), dimnames = list(NULL, columns))
  root <- design$random_allocations$root
  child <- design$random_allocations$child
  for (allocation in list(root, child)) {
    weights <- if (allocation$label == "root") c(.25, .75) else c(.2, .8)
    samples[, paste0(allocation$weight_name, "[", 1:2, "]")] <- rep(weights, each = 4L)
    samples[, .iwmde_simplex_auxiliary_columns(allocation$weight_name, 2L)] <- rep(2 * weights, each = 4L)
  }
  gate <- root$inclusion$study$indicator_name
  samples[, gate] <- c(0, 1, 0, 1)
  parameter <- .data_scale_component_specs(setup$object$data)[[1L]]$parameter
  samples[, review_coefficient_name(setup, "intercept", parameter)] <- 2
  samples[, review_coefficient_name(setup, "time", parameter)] <- log(2)
  output <- setup$compiled[[parameter]]$formula_design$name_map
  output <- output$jags_name[output$kind == "formula_output"]
  generated <- matrix(rep(c(log(2) + dat$time * log(2), 2 * 2^dat$time), each = 4L), nrow = 4L,
    dimnames = list(NULL, c(paste0(output, "[", seq_len(4L), "]"),
      paste0(source, "[", seq_len(4L), "]"))))
  samples <- cbind(samples, generated)
  object <- review_attach_replay(setup, samples)
  list(object = object, setup = setup, samples = samples, root = root, child = child,
    gate = gate, blocks = blocks)
}

review_scaled_allocation_replay <- function(shift = 0, standardized = TRUE) {

  dat <- data.frame(yi = seq(-.1, .2, length.out = 8L),
    study = rep(letters[1:4], each = 2L),
    x = rep(c(-sqrt(3.5), sqrt(3.5)), 4L) + shift)
  setup <- review_compile_object(brma.mv(
    yi = yi, V = diag(.1, 8L), mods = ~ x, random = ~ diag(1 + x | study),
    data = dat, measure = "GEN", prior_unit_information_sd = 1,
    standardize_continuous_predictors = standardized, only_priors = TRUE, silent = TRUE
  ))
  compiled <- setup$compiled$mu
  term <- compiled$formula_design$random_effects[[1L]]
  allocation <- term$sd_binding$allocations[[1L]]
  priors <- c(compiled$prior_list, .create_fit_priors(setup$object$data, setup$object$priors))
  weights <- paste0(allocation$weight_name, "[", seq_len(allocation$n_targets), "]")
  auxiliary <- .iwmde_simplex_auxiliary_columns(allocation$weight_name, allocation$n_targets)
  map <- compiled$formula_design$name_map
  latent <- map$jags_name[map$kind == "random" & map$role == "_xRE_Zx"]
  latent <- unlist(lapply(seq_len(term$n_columns), function(j) {
    paste0(latent, "[", seq_len(term$n_groups), ",", j, "]")
  }), use.names = FALSE)
  columns <- unique(c(review_prior_columns(priors), auxiliary, latent))
  samples <- matrix(0, 24L, length(columns), dimnames = list(NULL, columns))
  samples[, allocation$source$name] <- seq(.8, 1.2, length.out = 24L)
  samples[, weights] <- .5
  samples[, auxiliary] <- 1
  seed <- review_attach_replay(setup, samples)
  registry <- BayesTools::JAGS_deterministic_nodes(seed$fit)
  nodes <- registry$node[registry$family == "random_sd"]
  samples <- cbind(samples, BayesTools::JAGS_evaluate_deterministic(seed$fit, draws = samples, nodes = nodes))
  object <- review_attach_replay(setup, samples)
  BayesTools::JAGS_validate_fit_contract(object$fit,
    requires = c("name_encoding", "formula_name_map", "formula_design", "parameter_map", "draw_geometry"))
  list(object = object, samples = samples)
}

review_random_mean_dependency_replay <- function(factor = FALSE, dependent = TRUE) {

  dat <- data.frame(yi = c(0, 0), x = c(1, 2), study = factor(c("a", "b")))
  prior <- BayesTools::prior("gamma", list(shape = 2, rate = 1))
  random <- if (factor) ~ us(1 + x | study) else ~ id(0 + x | study)
  constructor <- function(x_prior) brma.mv(
    yi = yi, V = diag(2), mods = ~ 0 + x, random = random, data = dat,
    measure = "GEN", prior_unit_information_sd = 1, prior_mods = list(x = x_prior),
    prior_heterogeneity = BayesTools::prior_random(sd = prior),
    standardize_continuous_predictors = FALSE, only_priors = TRUE, silent = TRUE)
  discovery <- review_compile_object(constructor(prior))
  source <- discovery$compiled$mu$formula_design$random_effects[[1L]]$sd_parameter_names[[1L]]
  x_prior <- prior
  if (dependent) attr(x_prior, "multiply_by") <- source
  setup <- review_compile_object(constructor(x_prior))
  term <- setup$compiled$mu$formula_design$random_effects[[1L]]
  correlation <- if (factor) term$correlation$primitive_names else character()
  columns <- unique(c(review_prior_columns(setup$compiled$mu$prior_list), correlation))
  samples <- matrix(.6, 24L, length(columns), dimnames = list(NULL, columns))
  samples[, "mu_intercept"] <- 0
  samples[, "mu_x"] <- c(1, seq(.65, 1.75, length.out = 23L))
  samples[, source] <- c(1, seq(.7, 1.8, length.out = 23L))
  if (factor) samples[, correlation] <- .2
  object <- review_attach_replay(setup, samples)
  quantities <- BayesTools::parameter_catalog(object$fit)$quantities
  selected <- quantities[quantities$role == "random_sd", , drop = FALSE]
  selected <- selected[vapply(selected$extraction_key, function(key) identical(key$source_parameter, source), logical(1)), , drop = FALSE]
  stopifnot(nrow(selected) == 1L)
  selector <- .brma_random_parameter_io_labels(selected, "selector")
  target <- .brma_random_parameter_density_target(object, selector)
  context <- .iwmde_context(object)
  spec <- .iwmde_parameter_spec(context, target$parameter, target$parameter_spec)
  states <- .iwmde_row_states(context, 1L, parameter = target$parameter, parameter_spec = spec)
  list(object = object, context = context, states = states, parameter = target$parameter,
    replacement = .iwmde_replacement_spec(context, target$parameter, spec))
}

# The compiled formula designs of an 'only_priors' object.
review_compile_object <- function(object) {

  args <- .create_jags_formula_args(object[["data"]], object[["priors"]])
  compiled <- lapply(names(args[["formula_list"]]), function(parameter) {
    BayesTools::JAGS_formula(
      formula       = args[["formula_list"]][[parameter]],
      parameter     = parameter,
      data          = args[["formula_data_list"]][[parameter]],
      prior_list    = args[["formula_prior_list"]][[parameter]],
      formula_scale = args[["formula_scale_list"]][[parameter]],
      prior_random  = args[["formula_random_prior_list"]][[parameter]],
      random_effects_compile = args[["formula_random_effects_compile_list"]][[parameter]]
    )
  })
  names(compiled) <- names(args[["formula_list"]])
  object[["formula_design"]] <- lapply(compiled, `[[`, "formula_design")

  list(object = object, compiled = compiled)
}

# The JAGS name of the coefficient of a term of a compiled formula parameter.
review_coefficient_name <- function(setup, term, parameter = "mu") {

  map <- setup[["object"]][["formula_design"]][[parameter]][["name_map"]]
  out <- map[["jags_name"]][map[["kind"]] == "fixed" & map[["term"]] == term]
  stopifnot(length(out) == 1L)

  out
}

# The object of 'setup' with 'samples' (draws x JAGS names; a list of such
# matrices gives one chain each) as its fitted posterior, finalized through the
# BayesTools fit contract.
review_attach_replay <- function(setup, samples) {

  object <- setup[["object"]]
  chains <- if (is.list(samples)) samples else list(samples)
  fit <- structure(
    list(mcmc = coda::mcmc.list(lapply(chains, coda::mcmc)), sample = nrow(chains[[1L]])),
    class = c("BayesTools_fit", "runjags", "list")
  )
  attr(fit, "formula_design") <- object[["formula_design"]]
  scale_info <- lapply(setup[["compiled"]], `[[`, "formula_scale")
  if (any(lengths(scale_info) > 0L)) {
    attr(fit, "formula_scale") <- scale_info
  }
  compiled_priors <- list()
  for (compiled in setup[["compiled"]]) {
    compiled_priors <- c(compiled_priors, compiled[["prior_list"]])
  }
  attr(fit, "prior_list") <- c(
    compiled_priors, .create_fit_priors(object[["data"]], object[["priors"]])
  )
  object[["fit"]] <- as_bayestools_fit(fit)

  object
}

# A replay of a model with one factor 'g' (6 estimates, 'levels' repeated in
# equal blocks) on the location, a heterogeneity fixed at .1 and 'n' draws.
# With 'variable_intercept = FALSE' the intercept is a fixed zero. The slope
# draws are independent standard normal ('independent') or deterministic
# normal quantiles.
review_factor_replay <- function(levels, contrast, variable_intercept = TRUE,
                                 independent = FALSE, distribution = "mnormal",
                                 n = 100L, measure = "GEN") {

  dat <- data.frame(
    yi  = seq(-.1, .2, length.out = 6),
    sei = .1,
    g   = factor(rep(levels, each = 6L / length(levels)), levels = levels)
  )
  object <- brma(
    yi = yi, sei = sei, mods = ~ g, data = dat, measure = measure,
    prior_unit_information_sd = 1,
    prior_effect = if (variable_intercept) {
      BayesTools::prior("normal", list(mean = 0, sd = 1))
    } else {
      BayesTools::prior("point", list(location = 0))
    },
    prior_mods = list(g = BayesTools::prior_factor(
      distribution, list(mean = 0, sd = 1), contrast = contrast
    )),
    prior_heterogeneity = BayesTools::prior("point", list(location = .1)),
    only_priors = TRUE, silent = TRUE
  )
  setup  <- review_compile_object(object)
  g_name <- review_coefficient_name(setup, "g")
  p      <- length(levels) - 1L
  if (independent) {
    withr::local_seed(6123)
    slopes    <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
    intercept <- if (variable_intercept) stats::rnorm(n) else rep(0, n)
  } else {
    slopes <- matrix(stats::qnorm((seq_len(n) - .5) / n), nrow = n, ncol = p)
    if (p == 2L) slopes[, 2L] <- rev(slopes[, 1L]) / 2 + .1
    intercept <- if (variable_intercept) {
      rev(stats::qnorm((seq_len(n) - .5) / n)) / 3 + .2
    } else {
      rep(0, n)
    }
  }
  colnames(slopes) <- BayesTools:::.JAGS_prior_factor_names(
    g_name, setup[["compiled"]][["mu"]][["prior_list"]][[g_name]]
  )
  samples <- cbind(slopes, intercept, tau = rep(.1, n))
  colnames(samples)[p + 1L] <- review_coefficient_name(setup, "intercept")

  review_attach_replay(setup, samples)
}

# The coefficient columns that the sampler reports for a list of priors.
review_prior_columns <- function(priors) {

  unique(unlist(lapply(names(priors), function(parameter) {
    BayesTools:::.prior_linear_prior_columns(parameter, priors[[parameter]])
  })))
}

# Replays an allocated random-component source through the public
# constructor/compiler with supplied-row declarations.
review_allocated_source_replay <- function(kind, constant_weights = FALSE,
    source_prior = NULL, source_values = NULL, source_indicators = NULL,
    scale = "total_variance") {
  if (is.null(source_prior)) source_prior <- switch(kind,
    mixture = BayesTools::prior_mixture(list(BayesTools::prior('point', list(location = 1)),
      BayesTools::prior('gamma', list(shape = 2, rate = 2))), is_null = c(FALSE, FALSE)),
    point = BayesTools::prior('point', list(location = 1)),
    zero = BayesTools::prior('point', list(location = 0)),
    uniform = BayesTools::prior('uniform', list(a = 0, b = 1)))
  dat <- data.frame(yi = rep(0, 4), study = factor(c('a', 'a', 'b', 'b')),
    estimate = factor(seq_len(4)))
  heterogeneity <- BayesTools::prior_random(allocation = list(
    BayesTools::random_variance_allocation(name = 'split',
      terms = c(study = 'study', estimate = 'estimate'), target = 'block',
      scale = scale, sd = source_prior,
      weights = BayesTools::prior('dirichlet', list(alpha = c(1, 1))))))
  setup <- review_compile_object(RoBMA::brma.mv(yi = yi, V = diag(4), mods = ~1,
    random = list(study = ~id(1 | study), estimate = ~id(1 | estimate)),
    data = dat, measure = 'GEN', prior_unit_information_sd = 1,
    prior_effect = BayesTools::prior('point', list(location = 0)),
    prior_heterogeneity = heterogeneity, standardize_continuous_predictors = FALSE,
    only_priors = TRUE, silent = TRUE))
  allocation <- setup$compiled$mu$formula_design$random_effects[[1L]]$sd_binding$allocations[[1L]]
  source_parameter <- allocation$source$name
  weights <- paste0(allocation$weight_name, '[', 1:2, ']')
  auxiliaries <- .iwmde_simplex_auxiliary_columns(allocation$weight_name, 2L)
  indicator <- if (BayesTools::is.prior.mixture(source_prior)) paste0(source_parameter, '_indicator') else character()
  columns <- unique(c(review_prior_columns(setup$compiled$mu$prior_list), source_parameter,
    weights, auxiliaries, indicator))
  n <- if (kind == 'uniform') 32L else 24L
  samples <- matrix(.5, n, length(columns), dimnames = list(NULL, columns))
  samples[, 'mu_intercept'] <- 0
  samples[, source_parameter] <- switch(kind, mixture = c(rep(1, 2), seq(.4, 1.5, length.out = 22)),
    point = rep(1, n), zero = rep(0, n), uniform = seq_len(n) / 33)
  if (!is.null(source_values)) samples[, source_parameter] <- source_values
  samples[, weights[1L]] <- if (kind == 'uniform') c(.25, .75, ((seq_len(30) * 11) %% 31 + .5) / 31) else
    rep(c(.25, .75, seq(.15, .85, length.out = 10)), 2)
  samples[, weights[2L]] <- 1 - samples[, weights[1L]]
  if (constant_weights) samples[, weights] <- .5
  samples[, auxiliaries] <- 2 * samples[, weights]
  if (length(indicator)) samples[, indicator] <- if (!is.null(source_indicators)) {
    source_indicators
  } else if (kind == "mixture") c(rep(1, 2), rep(2, 22)) else rep(1:2, length.out = n)
  stopifnot(all(is.finite(samples)), all(samples[, weights] > 0), all(samples[, weights] < 1),
    max(abs(rowSums(samples[, weights]) - 1)) < 1e-15,
    max(abs(samples[, auxiliaries] / rowSums(samples[, auxiliaries]) - samples[, weights])) < 1e-15)
  object <- review_attach_replay(setup, samples)
  BayesTools::JAGS_validate_fit_contract(object$fit, requires = c("name_encoding", "formula_name_map", "formula_design", "parameter_map", "draw_geometry"))
  stopifnot(identical(class(object), class(setup$object)))
  list(object = object, samples = samples, allocation = allocation,
    source_parameter = source_parameter, weights = weights)
}

review_bma_gate_replay <- function(study_label = "study") {
  dat <- data.frame(yi = c(.1, .2, .3, .4),
    study = factor(c("s1", "s1", "s2", "s2")), other = factor(1:4))
  setup <- review_compile_object(BMA.mv(yi = yi, V = diag(.04, 4L),
    random = stats::setNames(list(~ 1 | study, ~ 1 | other), c(study_label, "other")), data = dat,
    prior_effect = BayesTools::prior("point", list(location = 0)), prior_effect_null = NULL,
    prior_heterogeneity = BayesTools::prior("gamma", list(shape = 2, rate = 2)),
    measure = "GEN", prior_unit_information_sd = 1, only_priors = TRUE, silent = TRUE))
  design <- setup$compiled$mu$formula_design
  allocation <- design$random_allocations[[1L]]
  priors <- c(setup$compiled$mu$prior_list, .create_fit_priors(setup$object$data, setup$object$priors))
  additional <- unique(unlist(lapply(setup$compiled, `[[`, "add_parameters")))
  latent_bases <- vapply(design$random_effects, function(term) paste0(term$parameter_stem, "_xRE_Zx"), character(1))
  latent <- unlist(lapply(design$random_effects, function(term) {
    base <- paste0(term$parameter_stem, "_xRE_Zx")
    if (!base %in% additional) return(character())
    index <- expand.grid(group = seq_len(term$n_groups), column = seq_len(term$n_columns))
    paste0(base, "[", index$group, ",", index$column, "]")
  }), use.names = FALSE)
  weight_columns <- paste0(allocation$weight_name, "[", 1:2, "]")
  auxiliary <- .iwmde_simplex_auxiliary_columns(allocation$weight_name, 2L)
  columns <- unique(c(review_prior_columns(priors), setdiff(additional, latent_bases), latent, auxiliary))
  monitors <- BayesTools::JAGS_to_monitor(priors)
  required <- unlist(lapply(monitors, function(base) {
    expanded <- columns[startsWith(columns, paste0(base, "["))]
    if (length(expanded)) expanded else base
  }), use.names = FALSE)
  columns <- unique(c(columns, required))
  samples <- matrix(0, 4L, length(columns), dimnames = list(NULL, columns))
  for (name in names(priors)) {
    prior <- priors[[name]]
    if (BayesTools::is.prior.point(prior)) samples[, name] <- prior$parameters$location
    if (BayesTools::is.prior.mixture(prior)) samples[, paste0(name, "_indicator")] <- 1
  }
  samples[, allocation$source$name] <- 1
  samples[, weight_columns] <- rep(c(.25, .75), each = 4L)
  samples[, auxiliary] <- 2 * samples[, weight_columns]
  gates <- vapply(allocation$inclusion, `[[`, character(1), "indicator_name")
  component_index <- match(c(study_label, "other"), allocation$component_names)
  gate_index <- match(component_index,
    vapply(allocation$inclusion, `[[`, integer(1), "index"))
  stopifnot(!anyNA(component_index), !anyNA(gate_index))
  samples[, gates[[gate_index[[1L]]]]] <- c(0, 1, 0, 1)
  samples[, gates[[gate_index[[2L]]]]] <- c(0, 0, 1, 1)
  preliminary <- review_attach_replay(setup, samples)
  registry <- BayesTools::JAGS_deterministic_nodes(preliminary$fit)
  nodes <- registry$node[registry$family == "random_sd"]
  generated <- BayesTools::JAGS_evaluate_deterministic(preliminary$fit, draws = samples, nodes = nodes)
  existing <- intersect(nodes, colnames(samples))
  if (length(existing)) samples[, existing] <- generated[, existing, drop = FALSE]
  extra <- setdiff(nodes, colnames(samples))
  if (length(extra)) samples <- cbind(samples, generated[, extra, drop = FALSE])
  stopifnot(all(required %in% colnames(samples)), all(is.finite(samples)))
  object <- review_attach_replay(setup, samples)
  BayesTools::JAGS_validate_fit_contract(object$fit, requires = c(
    "name_encoding", "formula_name_map", "formula_design", "parameter_map", "draw_geometry"))
  list(object = object, setup = setup, allocation = allocation, samples = samples,
    required = required, gates = gates, nodes = nodes)
}

review_fixed_zero_location_replay <- function() {
  dat <- data.frame(yi = rep(c(.1, 10), each = 3), study = factor(rep(1:3, 2)),
    g = factor(rep(c("A", "B"), each = 3)))
  setup <- review_compile_object(brma.mv(yi = yi, V = diag(.04, 6),
    mods = ~ 0 + g, random = ~ 1 | study, data = dat, measure = "GEN",
    prior_effect = BayesTools::prior("point", list(0)),
    prior_mods = list(g = BayesTools::prior_factor("normal", list(0, 1), contrast = "independent")),
    prior_heterogeneity = BayesTools::prior("point", list(.1)),
    prior_unit_information_sd = 1, only_priors = TRUE, silent = TRUE))
  priors <- c(setup$compiled$mu$prior_list, .create_fit_priors(setup$object$data, setup$object$priors))
  design <- setup$compiled$mu$formula_design
  additional <- setup$compiled$mu$add_parameters
  latent_bases <- vapply(design$random_effects, function(term) paste0(term$parameter_stem, "_xRE_Zx"), character(1))
  latent <- unlist(lapply(design$random_effects, function(term) {
    base <- paste0(term$parameter_stem, "_xRE_Zx")
    if (!base %in% additional) return(character())
    index <- expand.grid(group = seq_len(term$n_groups), column = seq_len(term$n_columns))
    paste0(base, "[", index$group, ",", index$column, "]")
  }), use.names = FALSE)
  columns <- unique(c(review_prior_columns(priors), setdiff(additional, latent_bases), latent))
  n <- 240L
  samples <- matrix(0, n, length(columns), dimnames = list(NULL, columns))
  for (name in names(priors)) if (BayesTools::is.prior.point(priors[[name]])) {
    samples[, name] <- priors[[name]]$parameters$location
  }
  name <- review_coefficient_name(setup, "g")
  factor_columns <- BayesTools:::.JAGS_prior_factor_names(name, priors[[name]])
  quantiles <- stats::qnorm((seq_len(n) - .5) / n)
  samples[, factor_columns] <- cbind(.1 + .12 * quantiles, 9.84 + .12 * rev(quantiles))
  object <- review_attach_replay(setup, samples)
  BayesTools::JAGS_validate_fit_contract(object$fit, requires = c("name_encoding",
    "formula_name_map", "formula_design", "parameter_map", "draw_geometry"))
  object
}


.plot_range_point_replay <- function(object) {

  setup <- review_compile_object(object)
  priors <- c(do.call(c, unname(lapply(setup$compiled, `[[`, "prior_list"))),
    .create_fit_priors(object$data, object$priors))
  monitors <- BayesTools::JAGS_to_monitor(priors)
  columns <- setdiff(review_prior_columns(priors), "omega")
  primitive <- matrix(0, 24L, length(columns), dimnames = list(NULL, columns))
  for (name in names(priors)) {
    if (BayesTools::is.prior.point(priors[[name]])) {
      primitive[, name] <- priors[[name]]$parameters$location
    }
  }
  preliminary <- review_attach_replay(setup, primitive)
  generated <- BayesTools::JAGS_evaluate_deterministic(preliminary$fit, draws = primitive)
  complete <- cbind(primitive, generated[, setdiff(colnames(generated), columns), drop = FALSE])
  required <- unlist(lapply(monitors, function(name) {
    if (name == "omega") return(paste0("omega[", seq_len(nrow(priors$omega$bins)), "]"))
    name
  }))
  stopifnot(!anyDuplicated(required), all(required %in% colnames(complete)),
    all(is.finite(complete)))
  samples <- complete[, required, drop = FALSE]
  object <- review_attach_replay(setup, samples)
  BayesTools::JAGS_validate_fit_contract(object$fit, requires = c(
    "name_encoding", "formula_name_map", "formula_design", "parameter_map", "draw_geometry"))
  list(object = object, samples = samples)
}


.review_cleanup_plain_replay <- function(known_v = FALSE, constant_mu = FALSE,
    fixed_tau = TRUE, singular = FALSE) {

  dat <- data.frame(yi = seq(.1, .6, .1), sei = .2, x = rep(c(-1, 2), 3))
  args <- list(yi = dat$yi, data = dat, measure = "GEN",
    prior_unit_information_sd = 1, only_priors = TRUE, silent = TRUE,
    prior_effect = prior("normal", list(0, 1)),
    prior_heterogeneity = if (fixed_tau) prior("point", list(0)) else
      prior("normal", list(0, 1), truncation = list(0, Inf)))
  if (known_v) {
    V <- diag(.04, 6L)
    V[1L, 2L] <- V[2L, 1L] <- .01
    args$V <- V
    args$prior_heterogeneity <- NULL
    args$mods <- if (singular) ~ x else ~ 1
    args$standardize_continuous_predictors <- FALSE
  } else args$sei <- dat$sei
  setup <- review_compile_object(do.call(if (known_v) brma.mv else brma, args))
  n <- if (known_v) 120L else 24L
  mu <- if (constant_mu) rep(.2, n) else if (known_v)
    .35 + .04 * stats::qnorm((seq_len(n) - .5) / n) else rep(c(0, 0, 1, 2), 6L)
  samples <- matrix(mu, ncol = 1L, dimnames = list(NULL, "mu"))
  if (!known_v) samples <- cbind(samples, tau = if (fixed_tau) 0 else .3)
  if (known_v) colnames(samples)[1L] <- review_coefficient_name(setup, "intercept")
  if (singular) {
    samples <- cbind(samples, mu)
    colnames(samples)[ncol(samples)] <- review_coefficient_name(setup, "x")
  }
  object <- review_attach_replay(setup, samples)
  BayesTools::JAGS_validate_fit_contract(object$fit, requires = c("name_encoding",
    "formula_name_map", "formula_design", "parameter_map", "draw_geometry"))
  object
}
