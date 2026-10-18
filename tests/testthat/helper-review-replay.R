# Synthetic replay of constructor objects for review regressions: the genuine
# formula designs and priors of an 'only_priors' object are compiled by
# BayesTools and finalized as a BayesTools fit around supplied posterior
# draws, so that post-fit methods run without sampling.

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
                                 n = 100L) {

  dat <- data.frame(
    yi  = seq(-.1, .2, length.out = 6),
    sei = .1,
    g   = factor(rep(levels, each = 6L / length(levels)), levels = levels)
  )
  object <- brma(
    yi = yi, sei = sei, mods = ~ g, data = dat, measure = "GEN",
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
