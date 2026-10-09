# Synthetic replay of constructor objects for review regressions: the genuine
# formula designs and priors of an 'only_priors' object are compiled by
# BayesTools and finalized as a BayesTools fit around supplied posterior
# draws, so that post-fit methods run without sampling.

# The compiled formula designs of an 'only_priors' object.
review_compile_object <- function(object) {

  args <- .create_jags_formula_args(object[["data"]], object[["priors"]])
  compiled <- lapply(names(args[["formula_list"]]), function(parameter) {
    BayesTools::JAGS_formula(
      formula       = args[["formula_list"]][[parameter]],
      parameter     = parameter,
      data          = args[["formula_data_list"]][[parameter]],
      prior_list    = args[["formula_prior_list"]][[parameter]],
      formula_scale = args[["formula_scale_list"]][[parameter]]
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

# The object of 'setup' with 'samples' (draws x JAGS names) as its fitted
# posterior, finalized through the BayesTools fit contract.
review_attach_replay <- function(setup, samples) {

  object <- setup[["object"]]
  fit <- structure(
    list(mcmc = coda::mcmc.list(coda::mcmc(samples)), sample = nrow(samples)),
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
