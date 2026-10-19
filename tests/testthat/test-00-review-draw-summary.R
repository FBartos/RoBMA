test_that("finite component SDs reject an overflowing total variance", {

  component <- matrix(sqrt(.75 * .Machine$double.xmax), 2L, 1L)
  expect_true(all(is.finite(component)))
  expect_error(.total_brma_mv_heterogeneity_samples(list(component, component)),
    "Total heterogeneity requires finite nonnegative variance.", fixed = TRUE)
  expect_identical(.total_brma_mv_heterogeneity_samples(list(matrix(0, 2L), matrix(0, 2L))), matrix(0, 2L))
  expect_equal(.total_brma_mv_heterogeneity_samples(list(matrix(3, 2L), matrix(4, 2L))), matrix(5, 2L))
})

test_that("completed allocations resolve to the authoritative prepared binding", {

  x <- review_heterogeneity_har_replay(scaled = TRUE)
  original <- .fitted_formula_design(x$object, "mu", required = TRUE)
  complete <- .brma_mv_resolve_sd_component_allocation(
    original$random_allocations[[1L]], original)
  expect_gt(length(complete$leaf_terms), 1L)
  prepared <- .brma_mv_pooled_formula_design(x$object)
  allocation <- .brma_mv_resolve_sd_component_allocation(complete, prepared)
  expect_true(is.function(allocation$source$source$values))
  expect_identical(allocation$source,
    prepared$random_effects[[1L]]$sd_binding$allocations[[1L]]$source)
})

test_that("RoBMA formula groups honor supplied source rows and current recipes", {

  object <- review_factor_replay(c("a", "b", "c"), "treatment",
    distribution = "normal", n = 24L)
  catalog <- .brma_parameter_catalog_metadata(object)[["catalog"]]
  group <- catalog[["quantities"]]
  group <- group[group[["provider"]] == "RoBMA" &
    group[["role"]] == "formula_coefficient_group", , drop = FALSE]
  expect_equal(nrow(group), 1L)
  selection <- BayesTools::parameter_catalog_resolve(
    catalog, alias = group[["canonical_name"]], namespace = group[["namespace"]],
    component = group[["component"]]
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
  normalized <- selection
  rownames(normalized[["quantities"]]) <- "unimportant_row_name"
  expect_identical(as.matrix(BayesTools::parameter_draws(object, normalized,
    model_samples = supplied)), supplied)
  corruptions <- list()
  for (field in c("schema_version", "parameter_map_version", "quantity_id")) {
    wrong <- selection
    wrong[[field]] <- if (identical(field, "quantity_id")) "RoBMA::foreign" else 999L
    corruptions[[paste0(field, "_wrong")]] <- wrong
    for (value in list(NULL, logical(), NA, list(1L), c(1L, 2L))) {
      malformed <- selection
      malformed[field] <- list(value)
      if (is.null(value)) malformed[[field]] <- NULL
      corruptions[[length(corruptions) + 1L]] <- malformed
    }
  }
  extra <- selection
  extra[["extra"]] <- TRUE
  corruptions[["extra_field"]] <- extra
  corruptions[["reordered_fields"]] <- selection[rev(names(selection))]
  for (invalid in corruptions) {
    expect_error(BayesTools::parameter_draws(object, invalid, model_samples = supplied),
      "'selection' does not identify a current RoBMA parameter quantity.", fixed = TRUE)
  }
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

test_that("row-source allocation anchors use prepared pooled source values", {

  x <- review_heterogeneity_har_replay(scaled = TRUE)
  for (object in list(x$object, x$setup$object)) {
    for (component in c("all", "study")) {
      anchor <- pooled_heterogeneity(object, component = component,
        .posterior_samples = x$samples)
      expect_identical(colnames(anchor), "tau_common")
      expect_equal(as.numeric(anchor), rep(sqrt(32), 120L), tolerance = 1e-13)
    }
    expect_equal(as.numeric(pooled_heterogeneity(object, component = "total",
      .posterior_samples = x$samples)), rep(4, 120L), tolerance = 1e-13)
    missing <- x$samples[, colnames(x$samples) != "log_tau_time", drop = FALSE]
    expect_error(pooled_heterogeneity(object, .posterior_samples = missing),
      "Pooled SD-component allocation is unavailable; missing source coordinates: log_tau_time.",
      fixed = TRUE)
  }
  observed <- summary_heterogeneity(x$object, component = "total")
  expect_equal(observed$estimates["tau", "Mean"], sqrt(40), tolerance = 1e-13)
})

test_that("row-source anchors retain admitted parent weights and gates once", {

  for (role in c("mean_variance", "total_variance")) {
    x <- review_heterogeneity_parent_replay(role)
    expected <- sqrt(32) * sqrt(.25) * c(0, 1, 0, 1)
    for (object in list(x$object, x$setup$object)) {
      anchor <- pooled_heterogeneity(object, component = x$blocks[1L],
        .posterior_samples = x$samples)
      expect_equal(as.numeric(anchor), expected, tolerance = 1e-13)
      expect_identical(colnames(anchor), if (role == "mean_variance") "tau_common" else "tau_total")
      missing <- x$samples[, colnames(x$samples) != x$gate, drop = FALSE]
      expect_error(pooled_heterogeneity(object, component = x$blocks[1L], .posterior_samples = missing),
        paste0("Cannot evaluate marginalized random-effect allocation; missing posterior column: ", x$gate), fixed = TRUE)
      weight <- paste0(x$root$weight_name, "[1]")
      missing <- x$samples[, colnames(x$samples) != weight, drop = FALSE]
      expect_error(pooled_heterogeneity(object, component = x$blocks[1L], .posterior_samples = missing),
        paste0("Cannot evaluate marginalized random-effect allocation; missing posterior column: ", weight), fixed = TRUE)
      multiplier <- if (role == "mean_variance") 2 else 1
      study_variance <- 32 * .25 * multiplier * sum(c(.5, .5)^2 * c(.2, .8)) * c(0, 1, 0, 1)
      total <- sqrt(study_variance + 32 * .75)
      expect_equal(as.numeric(pooled_heterogeneity(object, component = "total", .posterior_samples = x$samples)),
        total, tolerance = 1e-13)
    }
  }
})

test_that("fitted row-source summaries merge absolute quantities with catalog rows", {

  for (time in list(c(1, 2, 1, 2), c(2, 1, 1, 1))) {
    x <- review_heterogeneity_har_replay(scaled = TRUE, time = time)
    expected <- c(tau_common = sqrt(mean((2 * 2^time)^2)),
      tau2_common = mean((2 * 2^time)^2),
      `tau(time[1])` = 4, `tau(time[2])` = 8,
      `tau2(time[1])` = 16, `tau2(time[2])` = 64, rho = 0,
      `tau_mult(time[1])` = 1, `tau_mult(time[2])` = 1,
      `tau2_mult(time[1])` = 1, `tau2_mult(time[2])` = 1)
    for (component in c("all", "study")) {
      result <- summary_heterogeneity(x$object, component = component)
      expect_identical(rownames(result$estimates), names(expected))
      expect_equal(as.numeric(result$estimates[, "Mean"]), unname(expected), tolerance = 1e-13)
      frame <- as.data.frame(result)
      expect_identical(frame, data.frame(result))
      expect_type(frame$Mean, "double")
      expect_equal(frame$Mean, unname(expected), tolerance = 1e-13)
      expect_true(all(frame$component == "study"))
    }
    total <- summary_heterogeneity(x$object, component = "total")
    expect_equal(total$estimates["tau", "Mean"], expected[["tau_common"]], tolerance = 1e-13)
  }
  scalar <- review_heterogeneity_har_replay()
  expect_equal(nrow(summary_heterogeneity(scalar$object)$estimates), 11L)
  for (role in c("mean_variance", "total_variance")) {
    x <- review_heterogeneity_parent_replay(role)
    samples <- .brma_mv_allocation_sample_lists(x$object, x$samples)[[x$blocks[1L]]]
    gate <- c(0, 1, 0, 1)
    multiplier <- if (role == "mean_variance") 2 else 1
    for (i in 1:2) {
      sd <- c(4, 8)[i] * sqrt(.25) * sqrt(multiplier * c(.2, .8)[i]) * gate
      sd_name <- paste0(x$blocks[1L], ": tau(time[", i, "])")
      var_name <- paste0(x$blocks[1L], ": tau2(time[", i, "])")
      expect_equal(samples[[sd_name]], sd, tolerance = 1e-13)
      expect_equal(samples[[var_name]], sd^2, tolerance = 1e-13)
    }
    aggregate <- paste0("child: ", if (role == "mean_variance") "tau_common" else "tau_total")
    expect_equal(samples[[aggregate]], sqrt(40) * sqrt(.25) * gate, tolerance = 1e-13)
    expect_equal(sum(grepl("tau_mult(time[1])", names(samples), fixed = TRUE)), 1L)
    expect_equal(sum(grepl("rho", names(samples), fixed = TRUE)), 1L)
    for (index in list(seq_len(4L), which(gate == 0), which(gate == 1))) {
      out <- summary_heterogeneity(x$object, component = x$blocks[1L],
        .posterior_samples = x$samples[index, , drop = FALSE])
      expect_identical(as.data.frame(out), data.frame(out))
      expect_equal(out$estimates[aggregate, "Mean"],
        mean(sqrt(40) * sqrt(.25) * gate[index]), tolerance = 1e-13)
    }
  }
})
