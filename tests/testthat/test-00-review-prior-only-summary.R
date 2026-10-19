test_that("single named prior-only allocation summaries retain declared owners", {

  for (name in list(NULL, "named")) {
    replay <- review_heterogeneity_har_replay(n = 24L, component_name = name)
    fitted <- summary_heterogeneity(replay$object, .posterior_samples = replay$samples)
    prior_only <- summary_heterogeneity(replay$setup$object, .posterior_samples = replay$samples)
    expected <- rownames(fitted$estimates)
    # Correlation availability is a separate contract. Compare the common
    # declared SD/variance/allocation family only.
    expected <- expected[!grepl("rho", expected, fixed = TRUE)]
    expect_identical(rownames(prior_only$estimates), expected)
    expect_equal(unname(as.matrix(prior_only$estimates)),
      unname(as.matrix(fitted$estimates[expected, , drop = FALSE])), tolerance = 1e-13)
    if (!is.null(name)) {
      expect_true(isTRUE(replay$term$component_visible))
      expect_false(isTRUE(replay$term$has_explicit_name))
      expect_true(all(startsWith(expected, "named: ")))
    }
  }
})

test_that("prior-only block summaries use realized declared allocation activity", {

  replay <- review_bma_gate_replay()
  expected <- list(tau_total = c(0, .5, sqrt(.75), 1),
    tau2_total = c(0, .25, .75, 1),
    `tau2_prop(study)` = c(NA, 1, 0, .25),
    `tau2_prop(other)` = c(NA, 0, 1, .75))
  values <- lapply(list(fitted = replay$object, prior_only = replay$setup$object), function(object) {
    lists <- .brma_mv_allocation_sample_lists(object, replay$samples)
    expect_length(lists, 1L)
    samples <- lists[[1L]]
    for (name in names(expected)) {
      expect_equal(as.numeric(samples[[name]]), expected[[name]], tolerance = 1e-13)
    }
    for (name in c("tau2_prop(study)", "tau2_prop(other)")) {
      expect_identical(BayesTools::posterior_metadata(samples[[name]], "undefined_draws"),
        "allocation_active")
    }
    .summary_heterogeneity_brma_mv_allocations(object, replay$samples,
      probs = c(.025, .975))[[1L]]$estimates
  })
  expect_equal(unname(as.matrix(values$prior_only[names(expected), , drop = FALSE])),
    unname(as.matrix(values$fitted[names(expected), , drop = FALSE])), tolerance = 1e-13)
  expect_identical(unname(attr(values$prior_only, "footnotes")),
    unname(attr(values$fitted, "footnotes")))
  expect_gt(length(attr(values$fitted, "footnotes")), 0L)

  # Source zero does not turn declared active gates off. Recompute the actual
  # compiler deterministic monitors after changing these primitive states.
  samples <- replay$samples
  samples[, replay$allocation$source$name] <- 0
  samples[, replay$gates] <- 1
  generated <- BayesTools::JAGS_evaluate_deterministic(replay$object$fit,
    draws = samples, nodes = replay$nodes)
  samples[, colnames(generated)] <- generated
  for (object in list(replay$object, replay$setup$object)) {
    actual <- .brma_mv_allocation_sample_lists(object, samples)[[1L]]
    expect_equal(as.numeric(actual[["tau_total"]]), rep(0, 4L), tolerance = 0)
    expect_equal(as.numeric(actual[["tau2_prop(study)"]]), rep(.25, 4L), tolerance = 0)
    expect_equal(as.numeric(actual[["tau2_prop(other)"]]), rep(.75, 4L), tolerance = 0)
  }
})

test_that("prior-only child summaries apply declared parent factors once", {

  for (role in c("mean_variance", "total_variance")) {
    replay <- review_heterogeneity_parent_replay(role)
    BayesTools::JAGS_validate_fit_contract(replay$object$fit, requires = c(
      "name_encoding", "formula_name_map", "formula_design", "parameter_map", "draw_geometry"))
    component <- replay$blocks[1L]
    fitted <- .brma_mv_allocation_sample_lists(replay$object, replay$samples)[[component]]
    prior_only <- .brma_mv_allocation_sample_lists(replay$setup$object, replay$samples)[[component]]
    gate <- c(0, 1, 0, 1)
    multiplier <- if (role == "mean_variance") 2 else 1
    for (i in 1:2) {
      name <- paste0(component, ": tau(time[", i, "])")
      expected <- c(4, 8)[i] * sqrt(.25) * sqrt(multiplier * c(.2, .8)[i]) * gate
      expect_equal(as.numeric(fitted[[name]]), expected, tolerance = 1e-13)
      expect_equal(as.numeric(prior_only[[name]]), expected, tolerance = 1e-13)
    }
  }
})
