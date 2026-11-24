context("Row cache of singleton step-selection rows in the corrected samplers")
source(testthat::test_path("common-functions.R"))
skip_on_cran()

# The row-cache surrogate of singleton step-selection rows
# (src/samplers/StepRowCache.h) in the gated weight block and the screened
# indicator: which nodes take it, its value at frozen states of weight and
# label moves against the rows' graph densities (and its mutants), the nested
# weight layout of a RoBMA() fit, corrections accepted up to rounding, and
# unchanged draws without step-selection rows. Fixtures in
# helper-step-row-cache.R; the exactness of the chains against quadrature
# references runs in test-00-selection-corrected-weights-exactness.R.

.rows_build <- function(...) {

  invisible(gc())
  .rows_model(...)
}


test_that("step-selection rows give their weights to the gated block and their label to the screened indicator", {

  skip_if_not_installed("rjags")
  setup <- .rows_setup()
  nodes <- c("b", "r", "eta[1]", "eta[2]", "eta[3]", "pet")
  # PET reaches the rows' means: it keeps the gated slice sampler
  expect_identical(unname(.rows_build(setup)$samplers[nodes]),
                   c("RoBMA::ScreenedFinite", rep("RoBMA::GatedCorrectedBlock", 4L), "BayesTools::GatedSlice"))
  # the eligibility before the cache, and the 'slice' tier
  stock <- c("base::Finite", rep("BayesTools::GatedSlice", 5L))
  expect_identical(unname(.with_row_flags(list(disable = TRUE), .rows_build(setup))$samplers[nodes]), stock)
  expect_identical(unname(.rows_build(setup, sampler = "slice")$samplers[nodes]), stock)
  # an indicator that reaches the rows only through their SDs gains nothing
  # from the cache and keeps its sampler
  het <- .rows_build(setup, heterogeneity = TRUE)
  expect_identical(unname(het$samplers[c("b", "h")]), c("RoBMA::ScreenedFinite", "base::Finite"))
})


test_that("the row cache equals the rows' graph densities at frozen states of weight and label moves", {

  skip_if_not_installed("rjags")
  for (direction in c(1, -1)) {
    for (telescope in c(1, 0)) {
      setup <- .rows_setup(direction = direction, telescope = telescope)
      for (b in 1:5) {
        .rows_build(setup, b = b, r = .35, eta = c(.7, 1.6, .4), pet = .6)
        gaps <- .rows_frozen_gaps()
        label <- sprintf("direction %d, telescope %d, label %d", direction, telescope, b)
        expect_lt(gaps[["label"]], 1e-10, label = paste(label, "label scores"))
        expect_lt(gaps[["weights"]], 1e-10, label = paste(label, "weight moves"))
      }
    }
  }

  # under label 5 row 8's normalizer (weight 1e-300 on its bin) is below the
  # linear threshold: the cache evaluates it, and only it, by its node
  setup <- .rows_setup()
  .rows_build(setup, b = 5L)
  .rows_test_hook(clear = TRUE)
  probe <- .rows_probe()
  counters <- .rows_test_hook(clear = TRUE)$counters
  expect_gt(nrow(probe), 0L)
  expect_identical(counters[["direct"]], as.numeric(nrow(probe)))
  expect_identical(counters[["cached"]], 7 * nrow(probe))
  expect_lt(.rows_relative_gap(probe[, "cache"], probe[, "graph"]), 1e-10)
})


test_that("the row cache's mutants differ from the rows' graph densities", {

  skip_if_not_installed("rjags")
  setup <- .rows_setup()
  gaps <- function(mutant) {
    .with_row_flags(list(mutant = mutant), do.call(rbind, lapply(2:4, function(b) {
      .rows_build(setup, b = b, r = .35, eta = c(.7, 1.6, .4), pet = .6)
      .rows_frozen_gaps()
    })))
  }
  # the power omitted (rows 2, 4 and 6), the observed-bin weight omitted
  expect_gt(max(gaps(1L)), 1e-3)
  expect_gt(max(gaps(2L)), 1e-3)
  # a stale cache: the PET label moves the means
  expect_gt(max(gaps(3L)[, "label"]), 1e-3)
})


test_that("the nested and ordinary weight layouts of a RoBMA() fit keep the cache identities", {

  skip_if_not_installed("rjags")
  skip_if_not_installed("metadat")
  bias_side <- c("bias_indicator", "PET_1", "PEESE_1")
  # the default policy: the location block owns the bias label and the slopes
  # (frozen-state identities of its label scores in
  # test-00-location-bias-native.R), the gated block the weights; without
  # joint samplers the screened indicator owns the label and the gated slices
  # the slopes
  owners <- list(auto = rep("RoBMA::LocationBlock", 3L),
                 off = c("RoBMA::ScreenedFinite", "BayesTools::GatedSlice", "BayesTools::GatedSlice"))
  for (joint in names(owners)) for (nested in c(FALSE, TRUE)) {
    fit <- .rows_robma_fit(nested, joint = joint)
    table <- fit[["samplers"]]
    expect_identical(table[["Sampler"]][match(bias_side, table[["Node"]])], owners[[joint]])
    weights <- table[["Node"]][grepl("^(eta_component|omega_ratio|BayesTools_nested)", table[["Node"]])]
    expect_true(all(table[["Sampler"]][match(weights, table[["Node"]])] == "RoBMA::GatedCorrectedBlock"))
    for (label in 2:9) {
      context <- sprintf("sampler.joint %s, nested %s, label %d", joint, nested, label)
      hand <- .rows_recompile(fit, label)
      expect_identical(unname(hand$samplers[bias_side]), owners[[joint]], label = context)
      gaps <- .rows_frozen_gaps()
      if (joint == "off") {
        expect_lt(gaps[["label"]], 1e-10, label = paste(context, "label scores"))
      } else {
        expect_true(is.na(gaps[["label"]]), label = paste(context, "no screened label"))
      }
      expect_lt(gaps[["weights"]], 1e-10, label = paste(context, "weight moves"))
    }
  }
})


test_that("the corrections of cached rows are accepted up to rounding and are no evidence", {

  skip_if_not_installed("rjags")
  # the rows of the weights fixture (helper-selection-corrected-weights.R),
  # whose three labels mix
  setup <- .rows_weights_setup(1:3, x = rbind(c(2.4, 1.3, .6, -.3), c(2.1, 1.7, .9, -.8), c(2.9, 1.1, .2, -.5),
                                              c(2.2, 1.5, .4, -.2)))
  run <- function(flags = list()) {
    invisible(gc())
    .corrected_test_hook(clear = TRUE)
    chain <- .with_row_flags(flags, .block_chain(.rows_weights_model(setup), c("b", "r", "eta"), 3000L,
                                                 adapt = 200L, burnin = 200L))
    list(draws = chain$draws, kinds = .corrected_test_hook(clear = TRUE)$kinds)
  }
  exact <- run()
  # every label is visited
  expect_true(all(1:3 %in% exact$draws[, "b"]))
  for (kind in c("gated", "finite")) {
    expect_gt(exact$kinds[kind, "proposed"], 100)
    expect_identical(exact$kinds[kind, "accepted"], exact$kinds[kind, "proposed"], label = paste(kind, "accepted"))
    expect_identical(exact$kinds[kind, "evidence_proposed"], 0, label = paste(kind, "evidence"))
  }

  # a perturbed cache is corrected: its proposals are evidence and some are rejected
  perturbed <- run(list(perturbation = 1))
  for (kind in c("gated", "finite")) {
    expect_gt(perturbed$kinds[kind, "evidence_proposed"], 0)
    expect_lt(perturbed$kinds[kind, "accepted"], perturbed$kinds[kind, "proposed"])
  }
})


test_that("models without step-selection rows draw as before", {

  skip_if_not_installed("rjags")
  # the weights fixture of helper-selection-corrected-weights.R: cluster
  # children only
  setup <- .weights_setup(1:3)
  run <- function(disable) {
    invisible(gc())
    .with_row_flags(list(disable = disable),
                    .block_chain(.weights_model(setup), c("b", "r", "eta"), 2000L, adapt = 200L, burnin = 200L))
  }
  on <- run(FALSE)
  off <- run(TRUE)
  expect_identical(on$samplers, off$samplers)
  expect_identical(on$draws, off$draws)
})
