context("Influence cleanup preserves declared targets and draw populations")

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

test_that("known-V public influence diagnostics use canonical accurate target metadata", {

  object <- .review_cleanup_plain_replay(known_v = TRUE)
  # This computes actual scores and PSIS on the constructor/compiler replay.
  object <- add_loo(object)
  results <- list(dfbetas = dfbetas(object), covratio = covratio(object),
    dffits = dffits(object), cooks = cooks.distance(object))
  methods <- c(dfbetas = "dfbetas()", covratio = "covratio()",
    dffits = "dffits()", cooks = "cooks.distance()")
  for (name in names(results)) {
    target <- attr(results[[name]], "RoBMA_target")
    expect_identical(target$method, methods[[name]])
    expect_identical(target$deletion_target, "estimate")
    expect_identical(target$reported_target, if (name %in% c("dffits", "cooks"))
      "fixed_location_fitted_value" else "estimate-unit PSIS influence")
    expect_null(attr(results[[name]], "brma_mv_target"))
  }
  expect_identical(attr(dfbetas(object, return_loo_estimates = TRUE),
    "RoBMA_target"), attr(results$dfbetas, "RoBMA_target"))
})

test_that("known-V undefined COVRATIO return paths retain their PSIS target", {

  weights <- matrix(1 / 120, 120L, 6L)
  constant <- covratio(.review_cleanup_plain_replay(known_v = TRUE,
    constant_mu = TRUE), .weights = weights)
  singular <- covratio(.review_cleanup_plain_replay(known_v = TRUE,
    singular = TRUE), .weights = weights)
  for (result in list(constant, singular)) {
    expect_s3_class(result, "covratio.brma")
    expect_true(all(is.nan(result)))
    expect_type(attr(result, "note"), "character")
    expect_identical(attr(result, "RoBMA_target")$method, "covratio()")
    expect_identical(attr(result, "RoBMA_target")$deletion_target, "estimate")
    expect_null(attr(result, "brma_mv_target"))
  }
})

test_that("joint and covariance registry rows have no false estimate deletion", {

  for (method in c("add_marglik()/bridge_sampler()", "hatvalues()",
      "rstandard(conditioning_depth = 'marginal')", "vif()", "cluster/block LOO")) {
    expect_null(.brma_mv_target_metadata(method)$deletion_target, info = method)
  }
  for (method in c("log_lik(unit = 'estimate')", "add_loo()/loo()", "dfbetas()",
      "covratio()", "dffits()", "cooks.distance()")) {
    expect_identical(.brma_mv_target_metadata(method)$deletion_target,
      "estimate", info = method)
  }
})

# The historical31 genuine gated BMA.mv replay, adapted to current public
# deterministic and fit-finalization helpers. Component witnesses remain
# separate from the exact monitor expansion stored in the final object.
.review_cleanup_gated_correlation <- function() {

  dat <- data.frame(yi = seq(.1, .6, .1),
    study = rep(c("s1", "s2", "s3"), each = 2L), x = rep(c(-1, 2), 3L))
  setup <- review_compile_object(BMA.mv(yi = yi, V = diag(.0225, 6L),
    random = ~ us(1 + x | study), data = dat, measure = "GEN",
    prior_unit_information_sd = 1, only_priors = TRUE, silent = TRUE))
  priors <- c(do.call(c, unname(lapply(setup$compiled, `[[`, "prior_list"))),
    .create_fit_priors(setup$object$data, setup$object$priors))
  primitive <- cbind(mu_intercept_indicator = rep(c(1, 2), 12L),
    mu_intercept_component_1 = 0,
    mu_intercept_component_2 = seq(.1, .6, length.out = 24L),
    rep(c(.2, .4, .6), 8L), .5, rep(c(0, 1), each = 12L),
    rep(c(1, 2, 3), 8L), rep(c(3, 2, 1), 8L), rep(c(.15, .35, .65, .85), 6L))
  colnames(primitive)[4:9] <- c("mu__xRE_ALLOCx_heterogeneity__allocation_sd",
    "mu__xRE_ALLOCx_heterogeneity__include_component_1_prob",
    "mu__xRE_ALLOCx_heterogeneity__include_component_1_indicator",
    "prior_par_eta_mu__xRE_ALLOCx_study_components__weight[1]",
    "prior_par_eta_mu__xRE_ALLOCx_study_components__weight[2]",
    "mu__xREx__study_xRE_CORx_lkj_u[1]")
  allocation <- primitive[, 7:8, drop = FALSE]
  allocation <- allocation / rowSums(allocation)
  colnames(allocation) <- paste0("mu__xRE_ALLOCx_study_components__weight[", 1:2, "]")
  latent <- matrix(0, 24L, 6L, dimnames = list(NULL,
    as.vector(outer(1:3, 1:2, function(i, j)
      paste0("mu__xREx__study_xRE_Zx[", i, ",", j, "]")))))
  witness <- cbind(primitive, allocation, latent)
  seed <- review_attach_replay(setup, witness)
  generated <- BayesTools::JAGS_evaluate_deterministic(seed$fit, draws = witness)
  pool <- cbind(witness, generated[, setdiff(colnames(generated), colnames(witness)), drop = FALSE])
  declared <- unique(c(BayesTools::JAGS_to_monitor(priors),
    unlist(lapply(setup$compiled, `[[`, "add_parameters"), use.names = FALSE)))
  bases <- sub("\\[.*$", "", colnames(pool))
  missing <- declared[!vapply(declared, function(name)
    any(colnames(pool) == name | bases == name), logical(1))]
  stopifnot(length(missing) == 0L, all(is.finite(pool)),
    !anyDuplicated(colnames(pool)), all(allocation > 0), all(allocation < 1),
    all(rowSums(allocation) == 1))
  samples <- pool[, colnames(pool) %in% declared | bases %in% declared, drop = FALSE]
  object <- review_attach_replay(setup, samples)
  BayesTools::JAGS_validate_fit_contract(object$fit, requires = c("name_encoding",
    "formula_name_map", "formula_design", "parameter_map", "draw_geometry"))
  stopifnot(identical(class(object), class(setup$object)))
  list(object = object, samples = samples, witness = witness,
    gate = primitive[, 6L])
}

test_that("declared undefined correlations refuse the complete influence request clearly", {

  fixture <- .review_cleanup_gated_correlation()
  object <- fixture$object
  expect_identical(dim(fixture$samples), c(24L, 26L))
  expect_true(all(is.finite(fixture$samples)))
  owned <- BayesTools::JAGS_evaluate_deterministic(object$fit)
  shared <- intersect(colnames(owned), colnames(fixture$samples))
  expect_equal(unname(owned[, shared, drop = FALSE]),
    unname(fixture$samples[, shared, drop = FALSE]), tolerance = 0)
  catalog <- BayesTools::parameter_catalog(object$fit)
  quantity <- catalog$quantities
  quantity <- quantity[quantity$quantity == "cor" & !quantity$internal, , drop = FALSE]
  selection <- BayesTools::parameter_catalog_resolve(catalog,
    alias = quantity$canonical_name, namespace = quantity$namespace)
  public <- as.matrix(BayesTools::parameter_draws(object$fit, selection))
  expect_identical(which(is.na(public[, 1L])), 1:12)
  bundle <- .brma_random_parameter_bundle(object)
  parameter <- bundle$specs$parameter[bundle$specs$quantity == "cor"]
  selected <- .diagnostic_random_parameter_samples(object, parameter)
  expect_identical(BayesTools::posterior_metadata(as.matrix(selected), "undefined_draws"),
    stats::setNames("correlation", colnames(selected)))
  expect_identical(which(is.na(selected[, 1L])), 1:12)
  message <- function(diagnostic) paste0(diagnostic,
    " is unavailable for random-effect quantity '", parameter,
    "' (12 of 24 posterior draws defined) because the correlation is undefined ",
    "when either SD is zero. Use 'component = \"random\"' and 'parameter' to ",
    "explicitly select an everywhere-defined quantity.")
  weights <- matrix(1 / 24, 24L, 6L)
  expect_error(dfbetas(object, component = "random", .weights = weights),
    message("DFBETAS"), fixed = TRUE)
  expect_error(dfbetas(object, component = "random", parameter = parameter,
    .weights = weights), message("DFBETAS"), fixed = TRUE)
  expect_error(covratio(object, component = "random", parameter = parameter,
    .weights = weights), message("COVRATIO"), fixed = TRUE)

  on <- object
  on$fit <- BayesTools::JAGS_with_draws(object$fit,
    coda::mcmc.list(coda::mcmc(fixture$samples[fixture$gate == 1, , drop = FALSE])))
  BayesTools::JAGS_validate_fit_contract(on$fit, requires = c("name_encoding",
    "formula_name_map", "formula_design", "parameter_map", "draw_geometry"))
  expect_equal(as.numeric(as.matrix(dfbetas(on, component = "random",
    parameter = parameter, .weights = matrix(1 / 12, 12L, 6L)))),
    rep(0, 6L), tolerance = 1e-12)
  expect_equal(as.numeric(covratio(on, component = "random", parameter = parameter,
    .weights = matrix(1 / 12, 12L, 6L))), rep(1, 6L), tolerance = 1e-12)
  everywhere <- bundle$specs$parameter[bundle$specs$quantity == "sd"][1L]
  expect_true(all(is.finite(as.matrix(dfbetas(object, component = "random",
    parameter = everywhere, .weights = weights)))))
})

test_that("undeclared missing and infinite influence draws remain invalid", {

  weights <- matrix(1 / 24, 24L, 6L)
  for (value in c(NA_real_, Inf)) {
    samples <- matrix(seq_len(24), 24L, 1L)
    samples[1L, 1L] <- value
    expect_error(.dfbetas_internal(samples, weights),
      "'samples' must be a non-empty matrix of finite values.", fixed = TRUE)
    expect_error(.covratio_internal(samples, weights),
      "'samples' must be a non-empty matrix of finite values.", fixed = TRUE)
  }
})

test_that("GLMM LOO-PIT availability is checked before posterior setup", {

  object <- brma.glmm(ai = c(1L, 2L), bi = c(9L, 8L),
    ci = c(3L, 4L), di = c(7L, 6L), measure = "OR", only_data = TRUE)
  expect_error(rstudent(object), paste0("LOO-PIT residuals are not available ",
    "for binomial or Poisson GLMMs because a discrete PIT convention has not ",
    "been defined."), fixed = TRUE)
})

test_that("a declared empty fixed-scale set returns the existing undefined COVRATIO result", {

  object <- .review_cleanup_plain_replay()
  weights <- matrix(1 / 24, 24L, 6L)
  expect_identical(.fixed_tau_prior_value(object$priors), 0)
  table <- BayesTools::JAGS_estimates_table(object$fit, keep_parameters = "tau",
    remove_inclusion = TRUE, random_effects_summary = "none",
    remove_diagnostics = TRUE, return_samples = TRUE)
  expect_identical(dim(as.matrix(table)), c(24L, 0L))
  result <- covratio(object, type = "scale", .weights = weights)
  expect_true(all(is.nan(result)))
  expect_identical(names(result), .diagnostic_study_labels(object))
  expect_s3_class(result, "covratio.brma")
  expect_identical(attr(result, "note"), paste0("COVRATIO could not be computed ",
    "because no parameters with non-zero posterior variance remain; values ",
    "are reported as NaN."))
  expect_identical(as.data.frame(result), data.frame(result))
  expect_identical(names(as.data.frame(result)), c("component", "parameter", "value"))

  constant <- covratio(.review_cleanup_plain_replay(fixed_tau = FALSE),
    type = "scale", .weights = weights)
  expect_true(all(is.nan(constant)))
  expect_identical(names(constant), names(result))
  expect_identical(class(constant), class(result))
  expect_true(endsWith(attr(constant, "note"), attr(result, "note")))
  expect_error(covratio(object, component = "random", parameter = "not a quantity",
    .weights = weights))
  expect_error(.covratio_internal(matrix(numeric(), 0L, 0L), weights),
    "'samples' must be a non-empty matrix of finite values.", fixed = TRUE)
})
