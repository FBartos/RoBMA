# ============================================================================ #
# random-marginal-vcov.R
# ============================================================================ #
#
# Shared adapter from brma.mv fitted objects to the BayesTools random-effect
# marginal covariance backend.
#
# ============================================================================ #


.brma_mv_random_effects_marginal_vcov <- function(
    object, posterior_samples, blocks = NULL, diagonal_only = FALSE) {

  formula_fit <- .posterior_formula_fit(
    fit               = object[["fit"]],
    posterior_samples = posterior_samples,
    formula_design    = TRUE
  )
  formula_design <- if (.is_scale(object)) {
    .predict_known_v_formula_design_with_row_source_values(
      object = object,
      data   = object[["data"]]
    )
  } else {
    .fitted_formula_design(object, "mu", required = TRUE)
  }
  attr(formula_fit, "formula_design") <- list(mu = formula_design)

  location_priors <- attr(object[["fit"]], "prior_list")
  if (is.null(location_priors)) {
    location_priors <- formula_design[["prior_list"]]
  }
  if (is.null(location_priors)) {
    location_priors <- object[["priors"]][["location"]]
  }

  return(BayesTools::random_effects_marginal_vcov(
    fit               = formula_fit,
    parameter         = "mu",
    data              = object[["data"]][["location"]],
    posterior_samples = posterior_samples,
    prior_list        = location_priors,
    blocks            = blocks,
    diagonal_only     = diagonal_only
  ))
}
