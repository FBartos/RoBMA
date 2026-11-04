.review_scalar_source_wide_replay <- function() {

  dat <- data.frame(yi = c(.1, .2), x = c(1, 0), y = c(0, 1), study = c("a", "a"))
  setup <- review_compile_object(brma.mv(yi = yi, V = diag(.04, 2L),
    random = ~ diag(1 + x + y | study), data = dat, measure = "GEN",
    prior_unit_information_sd = 1,
    prior_effect = BayesTools::prior("point", list(location = 0)),
    prior_heterogeneity = BayesTools::prior("point", list(location = 1)),
    standardize_continuous_predictors = FALSE, only_priors = TRUE, silent = TRUE))
  term <- setup$compiled$mu$formula_design$random_effects[[1L]]
  allocation <- term$sd_binding$allocations[[1L]]
  stopifnot(term$n_columns == 3L, nrow(term$model_matrix) == 2L,
    identical(allocation$source$shape, "scalar"), identical(allocation$scale, "mean_variance"))
  priors <- c(setup$compiled$mu$prior_list, .create_fit_priors(setup$object$data, setup$object$priors))
  map <- setup$compiled$mu$formula_design$name_map
  latent_base <- map$jags_name[map$kind == "random" & map$role == "_xRE_Zx"]
  stopifnot(length(latent_base) == 1L, latent_base %in% setup$compiled$mu$add_parameters)
  index <- expand.grid(group = seq_len(term$n_groups), column = seq_len(term$n_columns))
  latent <- paste0(latent_base, "[", index$group, ",", index$column, "]")
  auxiliary <- .iwmde_simplex_auxiliary_columns(allocation$weight_name, allocation$n_targets)
  columns <- unique(c(review_prior_columns(priors), auxiliary, latent))
  primitive <- matrix(0, 24L, length(columns), dimnames = list(NULL, columns))
  for (name in names(priors)) if (BayesTools::is.prior.point(priors[[name]])) {
    primitive[, name] <- priors[[name]]$parameters$location
  }
  weights <- c(.2, .3, .5)
  weight_columns <- paste0(allocation$weight_name, "[", 1:3, "]")
  primitive[, weight_columns] <- rep(weights, each = 24L)
  primitive[, auxiliary] <- rep(3 * weights, each = 24L)
  initial <- review_attach_replay(setup, primitive)
  generated <- BayesTools::JAGS_evaluate_deterministic(initial$fit, draws = primitive)
  complete <- cbind(primitive, generated[, setdiff(colnames(generated), columns), drop = FALSE])
  monitors <- unique(c(BayesTools::JAGS_to_monitor(priors), setup$compiled$mu$add_parameters))
  required <- unique(unlist(lapply(monitors, function(base) {
    expanded <- colnames(complete)[startsWith(colnames(complete), paste0(base, "["))]
    if (length(expanded)) expanded else base
  }), use.names = FALSE))
  stopifnot(all(required %in% colnames(complete)), all(is.finite(complete)))
  samples <- complete[, required, drop = FALSE]
  object <- review_attach_replay(setup, samples)
  BayesTools::JAGS_validate_fit_contract(object$fit, requires = c(
    "name_encoding", "formula_name_map", "formula_design", "parameter_map", "draw_geometry"))
  stopifnot(identical(class(object), class(setup$object)))
  list(object = object, setup = setup, samples = samples, allocation = allocation,
    weights = weights, required = required)
}

test_that("scalar SD sources never index observation rows by design columns", {

  replay <- .review_scalar_source_wide_replay()
  expected <- sqrt(3 * replay$weights)
  labels <- .brma_mv_sd_component_allocation_labels(replay$allocation)
  rows <- paste0("tau(", labels, ")")
  for (object in list(replay$object, replay$setup$object)) {
    result <- tryCatch(summary_heterogeneity(object, component = "study",
      .posterior_samples = replay$samples), error = identity)
    expect_false(inherits(result, "error"))
    if (inherits(result, "error")) next
    expect_true(all(rows %in% rownames(result$estimates)))
    expect_equal(as.numeric(result$estimates[rows, "Mean"]), expected, tolerance = 1e-13)
  }
})
