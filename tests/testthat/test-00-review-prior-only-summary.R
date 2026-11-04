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

# Two random blocks over the same nested groups: the summaries of their
# allocations share the name "study/esid", so the complete result names the
# second one "study/esid_1".
.review_duplicate_allocation_names_replay <- function() {

  dat <- data.frame(yi = seq(.1, .8, length.out = 8L), study = factor(rep(1:2, each = 4L)),
    esid = factor(rep(1:2, 4L)))
  setup <- review_compile_object(brma.mv(yi = yi, V = diag(.1, 8L),
    random = list(a = ~ 1 | study/esid, b = ~ 1 | study/esid), data = dat,
    measure = "GEN", prior_unit_information_sd = 1, only_priors = TRUE, silent = TRUE))
  columns <- review_prior_columns(setup$compiled$mu$prior_list)
  samples <- matrix(1, 24L, length(columns), dimnames = list(NULL, columns))
  for (allocation in setup$compiled$mu$formula_design$random_allocations) {
    weights <- paste0(allocation$weight_name, "[", seq_len(allocation$n_targets), "]")
    samples[, weights] <- 1 / allocation$n_targets
  }

  list(setup = setup, samples = samples, fitted = review_attach_replay(setup, samples))
}

test_that("heterogeneity selectors accept the unique names of the complete result", {

  replay <- .review_duplicate_allocation_names_replay()
  for (object in list(replay$fitted, replay$setup$object)) {
    complete <- summary_heterogeneity(object, .posterior_samples = replay$samples)
    allocations <- c("a+b", "study/esid", "study/esid_1")
    expect_identical(names(complete)[1:3], allocations)
    all_lists <- .brma_mv_allocation_sample_lists(object, replay$samples)
    expect_identical(names(all_lists), allocations)
    expect_true(all(allocations %in%
      attr(all_lists, "available_components", exact = TRUE)))
    for (name in allocations) {
      selected <- summary_heterogeneity(object, component = name,
        .posterior_samples = replay$samples)
      expect_identical(selected[["component"]], name)
      expect_identical(selected[["estimates"]], complete[[name]][["estimates"]])
      lists <- .brma_mv_allocation_sample_lists(object, replay$samples,
        component = name)
      expect_identical(names(lists), name)
      expect_identical(lists[[name]], all_lists[[name]])
    }
  }
})

test_that("component selection evaluates only the selected allocation draws", {

  replay <- .review_duplicate_allocation_names_replay()
  calls <- 0L
  testthat::local_mocked_bindings(
    parameter_draws = local({
      original <- BayesTools::parameter_draws
      function(...) {
        calls <<- calls + 1L
        original(...)
      }
    }),
    .package = "BayesTools"
  )
  for (object in list(replay$fitted, replay$setup$object)) {
    calls <- 0L
    complete <- .brma_mv_allocation_sample_lists(object, replay$samples)
    complete_calls <- calls
    calls <- 0L
    selected <- .brma_mv_allocation_sample_lists(object, replay$samples,
      component = "study/esid_1")
    expect_identical(names(selected), "study/esid_1")
    expect_identical(selected[["study/esid_1"]], complete[["study/esid_1"]])
    if (!is.null(object[["fit"]])) {
      expect_gt(calls, 0L)
      expect_lt(calls, complete_calls)
    }
  }
})

# A nested allocation over scalar SD sources: the root splits an SD of 2 with
# variance shares (.25, .75), and the child allocation of the first share
# divides it over the leaves of the study block with shares (.2, .8).
.review_nested_scalar_source_replay <- function() {

  dat <- data.frame(yi = c(.1, .2, .3, .4), time = c(1, 2, 1, 2),
    study = factor(c("s1", "s1", "s2", "s2")), other = factor(1:4))
  args <- list(yi = quote(yi), V = diag(.04, 4L),
    random = list(combined = ~ har(time | study) + id(1 | other)), data = dat,
    measure = "GEN", prior_unit_information_sd = 1,
    prior_effect = BayesTools::prior("point", list(location = 0)),
    standardize_continuous_predictors = FALSE, only_priors = TRUE, silent = TRUE)
  blocks <- vapply(do.call(brma.mv, args)$formula_design$mu$random_effects,
    `[[`, character(1), "block_name")
  args$prior_heterogeneity <- BayesTools::prior_random(allocation = list(
    BayesTools::random_variance_allocation(name = "root",
      terms = c(study = blocks[1L], other = blocks[2L]),
      sd = BayesTools::prior("gamma", list(2, 2)),
      weights = BayesTools::prior("dirichlet", list(c(1, 1)))),
    BayesTools::random_variance_allocation(name = "child", terms = blocks[1L],
      target = "sd_component", parent = BayesTools::allocation_ref("root", "study"),
      scale = "mean_variance",
      weights = BayesTools::prior("dirichlet", list(c(1, 1))))))
  setup <- review_compile_object(do.call(brma.mv, args))
  allocations <- setup$compiled$mu$formula_design$random_allocations
  columns <- review_prior_columns(setup$compiled$mu$prior_list)
  samples <- matrix(.2, 24L, length(columns), dimnames = list(NULL, columns))
  samples[, allocations$root$source$name] <- 2
  for (allocation in allocations) {
    weights <- if (allocation$label == "root") c(.25, .75) else c(.2, .8)
    samples[, paste0(allocation$weight_name, "[", 1:2, "]")] <- rep(weights, each = 24L)
  }

  list(setup = setup, samples = samples, blocks = blocks, allocations = allocations,
    fitted = review_attach_replay(setup, samples))
}

test_that("prior-only pooled and summary SDs apply the declared parent factors", {

  replay <- .review_nested_scalar_source_replay()
  block <- replay$blocks[1L]
  # The root SD 2 scaled by the study variance share .25 is 1.
  fitted <- pooled_heterogeneity(replay$fitted, .posterior_samples = replay$samples)
  prior_only <- pooled_heterogeneity(replay$setup$object, .posterior_samples = replay$samples)
  expect_equal(as.numeric(fitted[[block]]), rep(1, 24L), tolerance = 1e-13)
  expect_equal(as.numeric(prior_only[[block]]), as.numeric(fitted[[block]]),
    tolerance = 1e-13)
  for (object in list(replay$fitted, replay$setup$object)) {
    summary <- summary_heterogeneity(object, component = block,
      .posterior_samples = replay$samples)
    row <- grep("tau_common$", rownames(summary$estimates))
    expect_length(row, 1L)
    expect_equal(unname(summary$estimates[row, "Mean"]), 1, tolerance = 1e-13)
  }
})
