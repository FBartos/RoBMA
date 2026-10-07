context("Ordered fitted-source draw exports")


# Formula metadata are real; the two chains contain independently specified
# primitives and deliberately unrelated raw increments. No fitting is needed
# for these extraction identities; live monitor evidence belongs to verification.
.ordered_draws_fixture <- function(formula = ~ f, priors, columns) {

  data <- data.frame(
    f = ordered(rep(c("base", "middle [nested] level", "high"), 8),
      levels = c("base", "middle [nested] level", "high")),
    h = ordered(rep(c("base", "middle [nested] level", "high"), each = 8),
      levels = c("base", "middle [nested] level", "high")),
    g = factor(rep(c("A", "B"), each = 12)), x = seq_len(24)
  )
  normal <- BayesTools::prior("normal", list(0, 1))
  compiled <- BayesTools::JAGS_formula(formula, "mu", data,
    c(list(intercept = normal), priors))
  values <- cbind(mu_intercept = seq_len(12) / 10,
    do.call(cbind, lapply(columns, function(x) rep(x, length.out = 12))))
  colnames(values) <- c("mu_intercept", names(columns))
  chains <- lapply(c(1, 2), function(chain) {
    draw_values <- values
    draw_values[, "mu_intercept"] <- draw_values[, "mu_intercept"] * chain
    coda::mcmc(draw_values, start = 11, thin = 2)
  })
  fit <- structure(list(mcmc = coda::mcmc.list(chains), sample = 12L),
    class = c("runjags", "BayesTools_fit", "list"))
  # Promote this declared compiler recipe to the fitted owner, as the current
  # source-complete test builders do; no primitive declaration is reconstructed.
  design <- compiled[["formula_design"]]
  scale <- BayesTools:::.bt_formula_scale_finalize(design$formula_scale, design,
    prior_list = compiled[["prior_list"]], owner_scope = "fit")
  design <- BayesTools:::.bt_formula_point_finalize(design,
    compiled[["prior_list"]], design$source_data, parameter_names = colnames(values))
  design$formula_scale <- scale
  attr(fit, "prior_list") <- compiled[["prior_list"]]
  attr(fit, "formula_design") <- list(mu = design)
  attr(fit, "formula_scale") <- list(mu = scale)
  fit <- BayesTools:::.bt_attach_parameter_map(fit)
  fit <- BayesTools:::.bt_attach_draw_geometry(fit)
  fit <- BayesTools:::.bt_attach_fit_contract(fit)
  structure(list(fit = fit), class = "brma")
}

.ordered_draws_columns <- function(total = c(0, 2, 3)) {

  list(mu_f_ordered_total = total,
    `prior_par_eta_mu_f_ordered_alloc_f_1[1]` = c(1, 2, 3),
    `prior_par_eta_mu_f_ordered_alloc_f_1[2]` = c(3, 2, 1),
    `mu_f[1]` = 99, `mu_f[2]` = 98)
}


test_that("ordered defaults export totals and Gamma shares even on zero-total rows", {

  object <- .ordered_draws_fixture(priors = list(f = BayesTools::prior_ordered(
    BayesTools::prior("normal", list(0, 1)))), columns = .ordered_draws_columns())
  extracted <- .brma_to_mcmc.list(object)
  expected_names <- c("mu_intercept", "mu_f_ordered_total",
    "mu_f_ordered_allocation[middle [nested] level]", "mu_f_ordered_allocation[high]")
  expect_identical(colnames(extracted[[1L]]), expected_names)
  for (chain in extracted) {
    expect_identical(coda::mcpar(chain), c(11, 33, 2))
    expect_identical(as.numeric(chain[, expected_names[[2L]]]), rep(c(0, 2, 3), 4))
    expect_identical(as.numeric(chain[, expected_names[[3L]]]), rep(c(.25, .5, .75), 4))
    expect_identical(as.numeric(chain[, expected_names[[4L]]]), rep(c(.75, .5, .25), 4))
  }
  auxiliary <- .brma_to_mcmc.list(object, include_auxiliary = TRUE)
  expect_true(all(c("mu_f[1]", "mu_f[2]", expected_names) %in% colnames(auxiliary[[1L]])))
  expect_identical(as.numeric(auxiliary[[1L]][, "mu_f[1]"]), rep(99, 12))
  expect_false(any(startsWith(colnames(auxiliary[[1L]]), "prior_par_eta_")))

  converters <- list(as_draws, as_draws_array, as_draws_df,
    as_draws_list, as_draws_matrix, as_draws_rvars)
  expected <- do.call(rbind, lapply(extracted, as.matrix))
  for (converter in converters) {
    draws <- converter(object)
    restored <- posterior::as_draws_array(draws)
    expect_identical(posterior::variables(restored), expected_names)
    expect_identical(posterior::nchains(restored), 2L)
    expect_identical(posterior::niterations(restored), 12L)
    expect_identical(as.numeric(posterior::as_draws_matrix(draws)), as.numeric(expected))
  }
})


test_that("point totals materialize every slice and fixed allocations add no shares", {

  point <- BayesTools::prior_ordered(BayesTools::prior("point", list(2.5)),
    allocation = c(.25, .75))
  object <- .ordered_draws_fixture(~ f * g,
    priors = list(f = point,
      g = BayesTools::prior_factor("normal", list(0, 1), contrast = "independent"),
      `f:g` = point),
    columns = list(`mu_f[1]` = 99, `mu_f[2]` = 98,
      `mu_g[1]` = 1, `mu_g[2]` = 2,
      `mu_f__xXx__g[1]` = 97, `mu_f__xXx__g[2]` = 96,
      `mu_f__xXx__g[3]` = 95, `mu_f__xXx__g[4]` = 94))
  extracted <- .brma_to_mcmc.list(object)
  totals <- c("mu_f_ordered_total", "mu_f__xXx__g_ordered_total[1]",
    "mu_f__xXx__g_ordered_total[2]")
  expect_true(all(totals %in% colnames(extracted[[1L]])))
  expect_identical(as.numeric(as.matrix(extracted[[1L]])[, totals]), rep(2.5, 36))
  expect_length(attr(extracted, "ordered_allocation_names", exact = TRUE), 0L)
  rvars <- as_draws_rvars(object)
  expect_identical(dim(rvars[["mu_f__xXx__g_ordered_total"]]), 2L)
  expect_identical(dim(rvars[["mu_g"]]), 2L)
  expect_identical(as.numeric(posterior::as_draws_matrix(rvars)[, totals]), rep(2.5, 72))

  public_calls <- 0L
  empty <- coda::mcmc.list(lapply(seq_len(2), function(i) {
    coda::mcmc(matrix(numeric(), 12L, 0L), start = 11, thin = 2)
  }))
  local_mocked_bindings(JAGS_materialize_draws = function(fit, parameters = NULL,
      include_internal = FALSE) {
    expect_null(parameters)
    expect_false(include_internal)
    public_calls <<- public_calls + 1L
    empty
  }, .package = "BayesTools")
  fixed <- .brma_to_mcmc.list(object)
  expect_identical(public_calls, 1L)
  expect_identical(colnames(fixed[[1L]]), totals)
  expect_identical(as.numeric(fixed[[1L]]), rep(2.5, 36))
  expect_identical(coda::mcpar(fixed[[1L]]), c(11, 33, 2))
  expect_false("RoBMAorderedRowCount" %in% colnames(fixed[[1L]]))
})


test_that("missing stochastic sources never use the fixed row-count fallback", {

  for (ordered in list(
      BayesTools::prior_ordered(BayesTools::prior("normal", list(0, 1))),
      BayesTools::prior_ordered(BayesTools::prior("point", list(2.5))),
      BayesTools::prior_ordered(BayesTools::prior("normal", list(0, 1)), allocation = c(.25, .75)),
      BayesTools::prior_ordered(BayesTools::prior("bernoulli", list(.5)), allocation = c(.25, .75)),
      BayesTools::prior_ordered(BayesTools::prior("point", list(expression(mu_intercept))),
        allocation = c(.25, .75)),
      BayesTools::prior_ordered(BayesTools::prior_mixture(list(
        BayesTools::prior("point", list(0)), BayesTools::prior("point", list(2.5)))),
        allocation = c(.25, .75)),
      BayesTools::prior_ordered(BayesTools::prior_spike_and_slab(
        BayesTools::prior("point", list(2.5)), BayesTools::prior("point", list(.5))),
        allocation = c(.25, .75)))) {
    object <- .ordered_draws_fixture(priors = list(f = ordered),
      columns = list(`mu_f[1]` = 99, `mu_f[2]` = 98))
    error <- tryCatch(.brma_to_mcmc.list(object), error = identity)
    expect_identical(class(error), c("RoBMA_ordered_draws_unavailable",
      "BayesTools_ordered_coordinates_unavailable", "BayesTools_ordered_unavailable",
      "RoBMA_refit_required", "BayesTools_refit_required", "error", "condition"))
    expect_identical(conditionMessage(error), paste0("Ordered source draws are unavailable. ",
      "Refit the model with the current RoBMA/BayesTools build."))
  }
})


test_that("shared allocations appear under each consuming term", {

  ordered <- BayesTools::prior_ordered(BayesTools::prior("normal", list(0, 1)), id = "shared")
  columns <- c(.ordered_draws_columns(), list(mu_x = 7,
    mu_f__xXx__x_ordered_total = c(4, 0, 6),
    `mu_f__xXx__x[1]` = 97, `mu_f__xXx__x[2]` = 96))
  names(columns)[2:3] <- c("prior_par_eta_ordered_alloc_shared_f[1]",
    "prior_par_eta_ordered_alloc_shared_f[2]")
  object <- .ordered_draws_fixture(~ f * x,
    priors = list(f = ordered, x = BayesTools::prior("normal", list(0, 1)), `f:x` = ordered),
    columns = columns)
  extracted <- .brma_to_mcmc.list(object)
  names <- attr(extracted, "ordered_allocation_names", exact = TRUE)
  expect_length(names, 4L)
  expect_identical(as.numeric(extracted[[1L]][, names[[1L]]]),
    as.numeric(extracted[[1L]][, names[[3L]]]))
  expect_identical(posterior::variables(posterior::as_draws_array(as_draws_rvars(object))),
    colnames(extracted[[1L]]))
})


test_that("nested spike total auxiliaries retain their declared event and inclusion", {

  ordered <- BayesTools::prior_ordered(BayesTools::prior_spike_and_slab(
    BayesTools::prior("normal", list(0, 1)), BayesTools::prior("beta", list(2, 3))))
  columns <- .ordered_draws_columns(total = c(0, 2, 0))
  columns <- c(columns, list(mu_f_ordered_total_variable = c(1, 2, 3),
    mu_f_ordered_total_indicator = c(0, 1, 0), mu_f_ordered_total_inclusion = c(.2, .4, .6)))
  object <- .ordered_draws_fixture(priors = list(f = ordered), columns = columns)
  extracted <- .brma_to_mcmc.list(object)
  expect_false("mu_f_ordered_total_variable" %in% colnames(extracted[[1L]]))
  expect_true(all(c("mu_f_ordered_total_indicator", "mu_f_ordered_total_inclusion") %in%
    colnames(extracted[[1L]])))
  expect_identical(as.numeric(extracted[[1L]][, "mu_f_ordered_total_indicator"]), rep(c(0, 1, 0), 4))
  expect_identical(as.numeric(extracted[[1L]][, "mu_f_ordered_total_inclusion"]), rep(c(.2, .4, .6), 4))
  expect_true("mu_f_ordered_total_variable" %in%
    colnames(.brma_to_mcmc.list(object, include_auxiliary = TRUE)[[1L]]))
  # A cell fixture must also have coherent semantic increments; the export
  # fixture above deliberately makes them unrelated to its source primitives.
  cell_columns <- columns
  cell_columns[["mu_f[1]"]] <- c(0, 1, 0)
  cell_columns[["mu_f[2]"]] <- c(0, 1, 0)
  cell_object <- .ordered_draws_fixture(priors = list(f = ordered), columns = cell_columns)
  catalog <- BayesTools::parameter_catalog(cell_object[["fit"]])
  quantities <- catalog[["quantities"]]
  middle <- which(vapply(quantities[["extraction_key"]], function(key) {
    identical(key[["type"]], "factor_level") &&
      identical(key[["dependencies"]], "mu_f[1]") &&
      identical(as.numeric(key[["weights"]]), 1)
  }, logical(1)))
  expect_length(middle, 1L)
  parameter <- quantities[["canonical_name"]][middle]
  selection <- BayesTools::parameter_catalog_resolve(catalog, parameter)
  entry <- list(parent_parameter = "mu_f", parameter = parameter,
    formula_parameter = "mu", selection = selection)
  marginal_posterior <- BayesTools::marginal_posterior
  selected_marginal <- NULL
  display_transform <- FALSE
  marginal_calls <- 0L
  testthat::local_mocked_bindings(marginal_posterior = function(...) {
    marginal_calls <<- marginal_calls + 1L
    selected_marginal <<- marginal_posterior(...)
    if (display_transform && marginal_calls == 2L) {
      selected_marginal <<- lapply(selected_marginal, BayesTools::posterior_transform,
        transformation = "exp")
    }
    selected_marginal
  }, .package = "BayesTools")
  for (conditional in c(FALSE, TRUE)) {
    cell <- .plot_brma_factor_cell_samples(cell_object, entry, FALSE, conditional, FALSE)
    value <- cell[["samples"]][[1L]]
    expected <- selected_marginal[[attr(value, "level_name", exact = TRUE)]]
    expect_identical(as.numeric(value), as.numeric(expected))
    expect_identical(attr(value, "bayestools_meta", exact = TRUE),
      attr(expected, "bayestools_meta", exact = TRUE))
    expect_length(value, if (conditional) 8L else 24L)
    expect_identical(BayesTools::posterior_metadata(value, "condition")[["averaged"]],
      !conditional)
    expect_identical(BayesTools::posterior_metadata(value, "atoms"),
      BayesTools::posterior_metadata(expected, "atoms"))
    if (conditional) expect_identical(as.numeric(value), rep(1, 8L))
  }
  # A controlled display supplier changes the semantic values after the raw
  # catalog-weight match. The helper must retain that transformed, conditioned
  # marginal even though its source contrast columns still have fitted values.
  display_transform <- TRUE
  marginal_calls <- 0L
  displayed <- .plot_brma_factor_cell_samples(cell_object, entry, FALSE, TRUE, FALSE)
  value <- displayed[["samples"]][[1L]]
  expected <- selected_marginal[[attr(value, "level_name", exact = TRUE)]]
  expect_identical(as.numeric(value), rep(exp(1), 8L))
  expect_identical(attr(value, "bayestools_meta", exact = TRUE),
    attr(expected, "bayestools_meta", exact = TRUE))
  expect_identical(BayesTools::posterior_metadata(value, "output_transformations"),
    BayesTools::posterior_metadata(expected, "output_transformations"))
  expect_false(BayesTools::posterior_metadata(value, "condition")[["averaged"]])
  display_transform <- FALSE
  columns[["mu_f_ordered_total_variable"]] <- NULL
  snapshot <- .ordered_draws_fixture(priors = list(f = ordered), columns = columns)
  expect_identical(as.numeric(.brma_to_mcmc.list(snapshot)[[1L]][, "mu_f_ordered_total"]),
    rep(c(0, 2, 0), 4))
  columns[["prior_par_eta_mu_f_ordered_alloc_f_1[1]"]] <- NULL
  missing <- .ordered_draws_fixture(priors = list(f = ordered), columns = columns)
  expect_error(.plot_brma_factor_cell_samples(missing, entry, FALSE, FALSE, FALSE),
    class = "BayesTools_ordered_coordinates_unavailable")
  error <- tryCatch(.brma_to_mcmc.list(missing), error = identity)
  expect_identical(class(error), c("RoBMA_ordered_draws_unavailable",
    "BayesTools_ordered_coordinates_unavailable", "BayesTools_ordered_unavailable",
    "RoBMA_refit_required", "BayesTools_refit_required", "error", "condition"))
  expect_identical(conditionMessage(error), paste0("Ordered source draws are unavailable. ",
    "Refit the model with the current RoBMA/BayesTools build."))
})


test_that("ordered mixture component pseudo-priors are auxiliary and their indicator is public", {

  ordered <- BayesTools::prior_ordered(BayesTools::prior_mixture(list(
    BayesTools::prior("normal", list(0, 1)), BayesTools::prior("normal", list(1, 2)))))
  columns <- c(.ordered_draws_columns(), list(
    mu_f_ordered_total_component_1 = c(0, 2, 3),
    mu_f_ordered_total_component_2 = c(4, 5, 6),
    mu_f_ordered_total_indicator = 1))
  object <- .ordered_draws_fixture(priors = list(f = ordered), columns = columns)
  extracted <- .brma_to_mcmc.list(object)
  components <- c("mu_f_ordered_total_component_1", "mu_f_ordered_total_component_2")
  expect_false(any(components %in% colnames(extracted[[1L]])))
  expect_identical(as.numeric(extracted[[1L]][, "mu_f_ordered_total_indicator"]), rep(1, 12))
  expect_true(all(components %in% colnames(
    .brma_to_mcmc.list(object, include_auxiliary = TRUE)[[1L]])))
})


test_that("sampled multi-slice totals and allocations retain fitted slices in every format", {

  ordered <- BayesTools::prior_ordered(BayesTools::prior("normal", list(0, 1)))
  columns <- c(.ordered_draws_columns(), list(`mu_g[1]` = 1, `mu_g[2]` = 2,
    `mu_f__xXx__g_ordered_total[1]` = c(4, 0, 6),
    `mu_f__xXx__g_ordered_total[2]` = c(8, 9, 0),
    `prior_par_eta_mu_f__xXx__g_ordered_alloc_f_1[1]` = 1,
    `prior_par_eta_mu_f__xXx__g_ordered_alloc_f_1[2]` = 4,
    `prior_par_eta_mu_f__xXx__g_ordered_alloc_f_2[1]` = 2,
    `prior_par_eta_mu_f__xXx__g_ordered_alloc_f_2[2]` = 3,
    `mu_f__xXx__g[1]` = 97, `mu_f__xXx__g[2]` = 96,
    `mu_f__xXx__g[3]` = 95, `mu_f__xXx__g[4]` = 94))
  object <- .ordered_draws_fixture(~ f * g, priors = list(f = ordered,
    g = BayesTools::prior_factor("normal", list(0, 1), contrast = "independent"),
    `f:g` = ordered), columns = columns)
  extracted <- .brma_to_mcmc.list(object)
  allocations <- attr(extracted, "ordered_allocation_names", exact = TRUE)
  expect_length(allocations, 6L)
  expect_identical(as.numeric(extracted[[1L]][, allocations[[3L]]]), rep(.2, 12))
  expect_identical(as.numeric(extracted[[1L]][, allocations[[5L]]]), rep(.4, 12))
  expected <- do.call(rbind, lapply(extracted, as.matrix))
  for (converter in list(as_draws, as_draws_array, as_draws_df,
      as_draws_list, as_draws_matrix, as_draws_rvars)) {
    restored <- posterior::as_draws_array(converter(object))
    expect_identical(posterior::variables(restored), colnames(expected))
    expect_identical(as.numeric(posterior::as_draws_matrix(restored)), as.numeric(expected))
    expect_identical(posterior::nchains(restored), 2L)
    expect_identical(posterior::niterations(restored), 12L)
  }
})


test_that("two ordered factors export each declared sampled allocation", {

  ordered <- BayesTools::prior_ordered(BayesTools::prior("normal", list(0, 1)), id = "shared")
  columns <- .ordered_draws_columns()
  names(columns)[2:3] <- c("prior_par_eta_ordered_alloc_shared_f[1]",
    "prior_par_eta_ordered_alloc_shared_f[2]")
  columns <- c(columns, list(mu_h_ordered_total = c(4, 0, 6),
    `prior_par_eta_ordered_alloc_shared_h[1]` = 1,
    `prior_par_eta_ordered_alloc_shared_h[2]` = 4,
    `mu_h[1]` = 97, `mu_h[2]` = 96,
    mu_f__xXx__h_ordered_total = c(8, 9, 0),
    `mu_f__xXx__h[1]` = 95, `mu_f__xXx__h[2]` = 94,
    `mu_f__xXx__h[3]` = 93, `mu_f__xXx__h[4]` = 92))
  object <- .ordered_draws_fixture(~ f * h,
    priors = list(f = ordered, h = ordered, `f:h` = BayesTools::prior_ordered(
      BayesTools::prior("normal", list(0, 1)), id = c(f = "shared", h = "shared"))),
    columns = columns)
  extracted <- .brma_to_mcmc.list(object)
  allocations <- attr(extracted, "ordered_allocation_names", exact = TRUE)
  expect_length(allocations, 8L)
  expect_identical(as.numeric(extracted[[1L]][, allocations[[1L]]]),
    as.numeric(extracted[[1L]][, allocations[[5L]]]))
  expect_identical(as.numeric(extracted[[1L]][, allocations[[3L]]]),
    as.numeric(extracted[[1L]][, allocations[[7L]]]))
  expect_identical(posterior::variables(posterior::as_draws_array(as_draws_rvars(object))),
    colnames(extracted[[1L]]))
})


test_that("rvars restore scalar allocation labels while preserving genuine vector indices", {

  names <- c("RoBMAorderedAllocation1", "RoBMAorderedAllocation2[1]",
    "RoBMAorderedAllocation2[2]", "ordinary[1]", "ordinary[2]",
    "ordered_total[1]", "ordered_total[2]", "allocation[1][middle [nested] level]",
    "allocation[2][middle [nested] level]", "other_allocation[middle [nested] level]")
  values <- matrix(seq_len(12 * length(names)), 12, dimnames = list(NULL, names))
  extracted <- coda::mcmc.list(coda::mcmc(values, start = 11, thin = 2),
    coda::mcmc(values + 1000, start = 11, thin = 2))
  attr(extracted, "ordered_allocation_names") <- tail(names, 3L)
  local_mocked_bindings(.brma_to_mcmc.list = function(...) extracted, .package = "RoBMA")
  draws <- as_draws_rvars(structure(list(), class = "brma"))
  expect_identical(dim(draws[["ordinary"]]), 2L)
  expect_identical(dim(draws[["RoBMAorderedAllocation2"]]), 2L)
  expect_identical(dim(draws[["ordered_total"]]), 2L)
  expect_true(all(vapply(draws[tail(names, 3L)], function(x) identical(dim(x), 1L), logical(1))))
  restored <- posterior::as_draws_array(draws)
  expect_identical(posterior::variables(restored), names)
  expect_identical(posterior::nchains(restored), 2L)
  expect_identical(posterior::niterations(restored), 12L)
  expect_identical(as.numeric(posterior::as_draws_matrix(draws)), as.numeric(rbind(values, values + 1000)))
})
