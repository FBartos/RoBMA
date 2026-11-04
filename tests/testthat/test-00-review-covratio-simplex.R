.review_ordered_influence <- function(kind = "mods", fixed = FALSE) {

  dat <- data.frame(yi = seq(.1, .6, .1), sei = 1, x = rep(c(-1, 0, 1), 2),
    f = ordered(rep(c("lo", "mid", "hi"), 2), levels = c("lo", "mid", "hi")))
  normal <- prior("normal", list(0, 1))
  point <- prior("point", list(0))
  ordered <- BayesTools::prior_ordered(if (fixed) prior("point", list(2)) else normal,
    id = if (kind == "shared") "shared" else NULL)
  args <- list(yi = dat$yi, sei = dat$sei, data = dat, prior_effect = normal,
    only_priors = TRUE, silent = TRUE, measure = "GEN", prior_unit_information_sd = 1)
  args <- c(args, if (kind == "scale") list(scale = ~ f, prior_scale = list(f = ordered),
    prior_heterogeneity = prior("normal", list(0, 1), truncation = list(0, Inf))) else list(mods = if (kind == "shared") ~ f * x else ~ f,
      prior_mods = if (kind == "shared") list(f = ordered, x = normal, `f:x` = ordered) else list(f = ordered),
      prior_heterogeneity = point))
  setup <- review_compile_object(do.call(brma, args))
  priors <- c(do.call(c, unname(lapply(setup$compiled, `[[`, "prior_list"))),
    .create_fit_priors(setup$object$data, setup$object$priors))
  grid <- expand.grid(intercept = if (kind == "scale") c(.5, 1.5) else c(-1, 1),
    total = c(1, 2), pair = 1:3, replicate = 1:2)
  primitive <- matrix(0, 24L, length(review_prior_columns(priors)),
    dimnames = list(NULL, review_prior_columns(priors)))
  seed <- review_attach_replay(setup, primitive)
  specs <- BayesTools::JAGS_ordered_parameter_spec(seed$fit)
  intercept <- review_coefficient_name(setup, "intercept", if (kind == "scale") "log_tau" else "mu")
  primitive[, intercept] <- grid$intercept
  append <- function(values, names) {
    values <- matrix(values, nrow = 24L, ncol = length(names), dimnames = list(NULL, names))
    primitive <<- cbind(primitive[, !colnames(primitive) %in% names, drop = FALSE], values)
  }
  for (j in seq_along(specs)) {
    spec <- specs[[j]]
    if (length(spec$total_names)) append(if (fixed) 2 else if (j == 1L) grid$total else grid$replicate + 2,
      spec$total_names)
    for (allocation in spec$allocations) if (!all(allocation$gamma_coordinates %in% colnames(primitive))) {
      append(rbind(c(1, 3), c(2, 2), c(3, 1))[grid$pair, ], allocation$gamma_coordinates)
    }
  }
  seed <- review_attach_replay(setup, primitive)
  generated <- BayesTools::JAGS_evaluate_deterministic(seed$fit, draws = primitive)
  pool <- cbind(primitive, generated[, setdiff(colnames(generated), colnames(primitive)), drop = FALSE])
  monitors <- BayesTools::JAGS_to_monitor(priors)
  names <- unique(unlist(lapply(monitors, function(name) {
    if (name %in% colnames(pool)) name else colnames(pool)[startsWith(colnames(pool), paste0(name, "["))]
  })))
  object <- review_attach_replay(setup, pool[, names, drop = FALSE])
  share <- c(.25, .5, .75)[grid$pair]
  first <- cbind(grid$intercept, if (!fixed) grid$total, share,
    if (kind == "shared" && !fixed) grid$replicate + 2)
  last <- first
  last[, if (fixed) 2L else 3L] <- 1 - share
  weights <- 1 + .2 * grid$intercept + .1 * (grid$total - 1) + .4 * share
  list(object = object, first = first, last = last, weights = weights / sum(weights))
}

test_that("COVRATIO uses declared simplex dimension and shared source identities", {

  ratio <- function(x, weights) det(stats::cov.wt(x, wt = weights, method = "ML")$cov) /
    det(stats::cov.wt(x, wt = rep(1/nrow(x), nrow(x)), method = "ML")$cov)
  for (kind in c("mods", "scale", "shared")) for (fixed in c(FALSE, TRUE)) {
    fixture <- .review_ordered_influence(kind, fixed)
    object <- fixture$object
    component <- if (kind == "scale") "scale" else "mods"
    expected <- ratio(fixture$first, fixture$weights)
    expect_equal(expected, ratio(fixture$last, fixture$weights), tolerance = 1e-12)
    expect_equal(as.numeric(covratio(object, component = component,
      .weights = matrix(1/24, 24L, 6L))), rep(1, 6L), tolerance = 1e-12)
    expect_equal(as.numeric(covratio(object, component = component,
      .weights = matrix(fixture$weights, 24L, 6L))), rep(expected, 6L), tolerance = 1e-12)
    # Individual shares remain available in DFBETAS, independently of the
    # determinant coordinate basis used only by COVRATIO.
    specs <- BayesTools::JAGS_ordered_parameter_spec(object$fit, draws = .get_posterior_samples(object$fit))
    shares <- unlist(lapply(specs, function(spec) BayesTools::parameter_labels(
      spec$sampled_parameters$label_parts[spec$sampled_parameters$allocation_names],
      style = "table", formula_prefix = TRUE)), use.names = FALSE)
    expect_true(all(shares %in% colnames(dfbetas(object, component = component,
      standardized_coefficients = TRUE, transform_factors = FALSE,
      .weights = matrix(1/24, 24L, 6L)))))
  }
})
