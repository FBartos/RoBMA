.assign_prior.heterogeneity_allocation <- function(prior) {

  # use default settings if prior distribution is not specified
  if (missing(prior) || is.null(prior)) {
    prior <- BayesTools::prior("beta", parameters = list("alpha" = 1, "beta" = 1))
    return(prior)
  }

  # check the user specified prior distribution
  prior <- .check_prior.restricted_01(prior, prior_name = "heterogeneity_allocation")

  return(prior)
}

.assign_prior.heterogeneity_allocation_mixture <- function(prior, prior_null) {

  ### check the alternative prior input
  # prior can be either
  # - unspecified (set default as in brma)
  # - NULL/FALSE (omit prior distribution)
  # - single prior distribution
  # - list of prior distributions
  # each of these prior distributions needs to satisfy the same properties as in brma
  # if only a single distribution is specified, place it in a list for a simpler treatment
  if (missing(prior) || BayesTools::is.prior(prior)) {
    priors_alt <- list(.assign_prior.heterogeneity_allocation(prior = prior))
  } else if (is.null(prior) || isFALSE(prior)) {
    priors_alt <- list()
  } else if (is.list(prior) && all(sapply(prior, BayesTools::is.prior))) {
    priors_alt <- prior
    for (i in seq_along(priors_alt)) {
      priors_alt[[i]] <- .assign_prior.heterogeneity_allocation(prior = priors_alt[[i]])
    }
  } else {
    stop("The 'prior_heterogeneity_allocation' must be either prior distribution or a list of prior distributions.", call. = FALSE)
  }

  ### check the null prior input
  # prior can be either
  # - empty (omit prior as default)
  # - NULL/FALSE (omit prior distribution)
  # - single prior distribution
  # - list of prior distributions
  # each of these prior distributions needs to satisfy the same properties as in brma
  # if only a single distribution is specified, place it in a list for a simpler treatment
  if (missing(prior_null) || is.null(prior_null) || isFALSE(prior_null)) {
    priors_null <- list()
  } else if (BayesTools::is.prior(prior_null)) {
    priors_null <- list(.assign_prior.heterogeneity_allocation(prior = prior_null))
  } else if (is.list(prior_null) && all(sapply(prior_null, BayesTools::is.prior))) {
    priors_null <- prior_null
    for (i in seq_along(priors_null)) {
      priors_null[[i]] <- .assign_prior.heterogeneity_allocation(prior = priors_null[[i]])
    }
  } else {
    stop("The 'prior_heterogeneity_allocation_null' must be either prior distribution or a list of prior distributions.", call. = FALSE)
  }

  if (length(priors_null) + length(priors_alt) == 0) {
    stop("At least one prior distribution needs to be defined for 'heterogeneity_allocation'.", call. = FALSE)
  }

  ### specify the corresponding mixture prior
  prior <- BayesTools::prior_mixture(
    prior_list = c(priors_null, priors_alt),
    is_null    = c(rep(TRUE, length(priors_null)), rep(FALSE, length(priors_alt)))
  )

  return(prior)
}
