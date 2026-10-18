test_that("RoBMA formula groups honor supplied source rows and current recipes", {

  object <- review_factor_replay(c("a", "b", "c"), "treatment",
    distribution = "normal", n = 24L)
  catalog <- .brma_parameter_catalog_metadata(object)[["catalog"]]
  group <- catalog[["quantities"]]
  group <- group[group[["provider"]] == "RoBMA" &
    group[["role"]] == "formula_coefficient_group", , drop = FALSE]
  expect_equal(nrow(group), 1L)
  selection <- BayesTools::parameter_catalog_resolve(
    catalog, alias = group[["canonical_name"]], component = "mods"
  )
  dependencies <- group[["extraction_key"]][[1L]][["dependencies"]]
  supplied <- matrix(seq_len(14L) / 10, nrow = 7L,
    dimnames = list(NULL, dependencies))
  original <- serialize(object[["fit"]], NULL)
  for (source in list(supplied, coda::mcmc(supplied, start = 101, thin = 3))) {
    draws <- BayesTools::parameter_draws(object, selection, model_samples = source)
    expect_identical(as.matrix(draws), supplied)
    expect_identical(attr(draws[[1L]], "mcpar"),
      if (inherits(source, "mcmc")) c(101, 119, 3) else c(1, 7, 1))
  }
  reordered <- cbind(extra = 1:7, supplied[7:1, rev(dependencies), drop = FALSE])
  draws <- BayesTools::parameter_draws(object, selection,
    model_samples = coda::mcmc(reordered, start = 201, thin = 2))
  expect_identical(as.matrix(draws), supplied[7:1, , drop = FALSE])
  expect_identical(attr(draws[[1L]], "mcpar"), c(201, 213, 2))
  expect_equal(nrow(as.matrix(BayesTools::parameter_draws(object, selection))), 24L)
  expect_error(BayesTools::parameter_draws(object, selection, unexpected = TRUE),
    "Unused argument in parameter_draws.brma(): 'unexpected'", fixed = TRUE)
  for (source in list(as.data.frame(supplied),
      `colnames<-`(supplied, rep(dependencies[1L], 2L)))) {
    expect_error(BayesTools::parameter_draws(object, selection, model_samples = source),
      "'model_samples' must be a numeric matrix with unique column names.", fixed = TRUE)
  }
  expect_error(BayesTools::parameter_draws(object, selection,
    model_samples = supplied[, 1L, drop = FALSE]),
    paste0("The supplied 'model_samples' are missing declared source coordinates: ",
      dependencies[2L], "."), fixed = TRUE)
  stale <- selection
  stale[["quantities"]][["quantity_id"]] <- "RoBMA::not_current"
  forged <- selection
  forged[["quantities"]][["extraction_key"]][[1L]][["dependencies"]] <- "tau"
  for (invalid in list(stale, forged)) {
    expect_error(BayesTools::parameter_draws(object, invalid, model_samples = supplied),
      "'selection' does not identify a current RoBMA parameter quantity.", fixed = TRUE)
  }
  scalar <- BayesTools::parameter_catalog_resolve(
    catalog, alias = "g[b]", namespace = "mu"
  )
  draws <- BayesTools::parameter_draws(object, scalar,
    model_samples = coda::mcmc(supplied, start = 101, thin = 3))
  expect_identical(as.numeric(draws[[1L]][, 1L]), supplied[, 1L])
  expect_identical(attr(draws[[1L]], "mcpar"), c(101, 119, 3))
  expect_identical(serialize(object[["fit"]], NULL), original)
})

test_that("pooled totals and prediction intervals use full allocation covariance", {

  x <- review_heterogeneity_har_replay()
  anchor <- pooled_heterogeneity(x$object)
  total <- pooled_heterogeneity(x$object, component = "total")
  expect_identical(colnames(anchor), "tau_common")
  expect_equal(as.numeric(anchor), rep(2, 120L), tolerance = 0)
  expect_identical(colnames(total), "tau_total")
  expect_equal(as.numeric(total), rep(sqrt(2), 120L), tolerance = 1e-14)
  withr::local_seed(2401)
  expected <- rnorm(120L, 0, sqrt(2))
  set.seed(2401)
  seed <- .Random.seed
  kind <- RNGkind()
  effect <- pooled_effect(x$object)
  expect_equal(as.numeric(attr(effect, "prediction_samples")), expected, tolerance = 1e-14)
  expect_identical(.Random.seed, seed)
  expect_identical(RNGkind(), kind)
  unequal <- review_heterogeneity_har_replay(weights = c(.2, .8), rho = .4,
    time = c(1, 1, 1, 2))
  qbar <- c(.75, .25)
  sd <- 2 * sqrt(2 * c(.2, .8))
  G <- diag(sd) %*% matrix(c(1, .4, .4, 1), 2L) %*% diag(sd)
  expected <- sqrt(drop(t(qbar) %*% G %*% qbar))
  expect_equal(as.numeric(pooled_heterogeneity(unequal$object, component = "total")),
    rep(expected, 120L), tolerance = 1e-13)
  direct <- review_heterogeneity_har_replay(weights = c(.2, .8), rho = .4,
    time = c(1, 1, 1, 2), direct = TRUE)
  expect_equal(as.numeric(pooled_heterogeneity(direct$object, component = "total")),
    rep(expected, 120L), tolerance = 1e-13)
  row <- review_heterogeneity_har_replay(scaled = TRUE)
  expect_equal(as.numeric(pooled_heterogeneity(row$object, component = "total")),
    rep(4, 120L), tolerance = 1e-13)
})
