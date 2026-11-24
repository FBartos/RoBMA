context("Corrected samplers of selection weights and indicators")
source(testthat::test_path("common-functions.R"))
skip_on_cran()

# The gated weight block (RoBMA::GatedCorrectedBlock), the screened indicator
# (RoBMA::ScreenedFinite) and the gamma-sum coordinates of the corrected block
# samplers (src/samplers/GatedCorrectedBlock.cc, ScreenedFinite.cc,
# GammaGroup.cc, CorrectedBlockKernel.cc): which nodes they claim and decline,
# their scores at frozen states, and draws that the two-entry node memo and
# the samplers' absence leave unchanged. Fixtures in
# helper-selection-corrected-weights.R; the exactness of the chains against
# quadrature references runs in test-00-selection-corrected-weights-exactness.R
# and the fits in the certification case 'numerical-kernels' below.

.weights_build <- function(...) {

  invisible(gc())
  .weights_model(...)
}


test_that("selection weights get the gated block sampler and indicators the screened one", {

  skip_if_not_installed("rjags")
  setup <- .weights_setup(1:3)
  hand <- .weights_build(setup)
  expect_identical(unname(hand$samplers[c("b", "r", "eta[1]", "eta[2]", "eta[3]")]),
                   c("RoBMA::ScreenedFinite", rep("RoBMA::GatedCorrectedBlock", 4L)))
  # eta is held by two log-ratio directions (its sum is drawn), r by one; both
  # start at 2.5 prior standard deviations of their coordinates (logit of a
  # Beta(1, 1), log ratio of two Gamma(1, 1)): 2.5 sqrt(2 trigamma(1))
  widths <- .weights_block_widths()
  expect_identical(sort(vapply(widths, `[[`, numeric(1L), "d")), c(1, 2))
  expect_lt(max(abs(unlist(lapply(widths, `[[`, "widths")) - 2.5 * sqrt(2 * trigamma(1)))), 1e-12)

  # the heterogeneity fixture: the indicator is screened and the declared block
  # holds the allocation SD and one log-ratio of its Dirichlet weights
  het <- .heterogeneity_model(.heterogeneity_setup())
  expect_identical(unname(het$samplers[c("h", "s", "e[1]", "e[2]")]),
                   c("RoBMA::ScreenedFinite", rep("RoBMA::CoarseCorrectedBlock", 3L)))
  expect_identical(.weights_block_widths()[[length(.weights_block_widths())]][["d"]], 2)
})


test_that("the samplers decline what they cannot claim exactly or usefully", {

  skip_if_not_installed("rjags")
  setup <- .weights_setup(1:3)
  stock <- c("base::Finite", rep("BayesTools::GatedSlice", 4L))
  nodes <- c("b", "r", "eta[1]", "eta[2]", "eta[3]")

  # the test hook's switch, and the 'slice' tier
  hand <- .with_corrected_flags(list(disable_new = TRUE), .weights_build(setup))
  expect_identical(unname(hand$samplers[nodes]), stock)
  expect_identical(unname(.weights_build(setup, sampler = "slice")$samplers[nodes]), stock)

  # a fixed normal kernel: no child can use a surrogate that differs from its
  # exact density, so none of the samplers is created
  normal <- .weights_build(setup, mode = "0")
  expect_identical(unname(normal$samplers[nodes]), stock)

  # a gated parameter that reaches the rows' mean (a PET-like slope, gated by a
  # label of its own) keeps the gated sampler; the weights stay with the gated
  # block (a group of one gate node is claimed or declined whole)
  slope <- .weights_build(setup, means = "mean + slope * equals(b, 1)",
    extra = "  slope ~ dnorm(0, 1)T(0,)\n  gate_slope ~ dbt_gate(slope, equals(b, 1))",
    extra_data = list(mean = setup$mean, gate_slope = 0), inits = list(b = 1L, r = .5, eta = c(1, 1, 1), slope = .1))
  expect_identical(unname(slope$samplers[c("slope", "r")]), c("BayesTools::GatedSlice", "RoBMA::GatedCorrectedBlock"))

  # gamma members that reach a node other than their normalization, or whose
  # rates differ, keep their own (log) coordinates: three directions
  for (variant in list(list(extra = "  twice <- 2 * eta[1]"), list(rates = c("1", "1", "2")))) {
    hand <- do.call(.weights_build, c(list(setup), variant))
    expect_identical(unname(hand$samplers["eta[1]"]), "RoBMA::GatedCorrectedBlock")
    expect_true(3 %in% vapply(.weights_block_widths(), `[[`, numeric(1L), "d"))
  }

  # an indicator with a finite readout tag (sampler.probabilities) keeps the
  # enumerating sampler, whose recorded probabilities are the exact conditionals
  tagged <- .weights_build(setup, extra = "  tag <- bt_finite_tag(b)")
  expect_false(identical(unname(tagged$samplers["b"]), "RoBMA::ScreenedFinite"))
})


test_that("the screened scores are the surrogate and exact densities of every level at frozen states", {

  skip_if_not_installed("rjags")
  setup <- .weights_setup(1:3)
  states <- list(list(b = 1L, r = .3, eta = c(.5, 1, 2)), list(b = 2L, r = .8, eta = c(2, .2, .7)),
                 list(b = 3L, r = .05, eta = c(1, 3, .1)))
  for (state in states) {
    hand <- .weights_build(setup, inits = state)
    scores <- .weights_screened_scores()
    scores <- scores[[length(scores)]]
    p <- state$eta / sum(state$eta)
    omegas <- list(c(1, 1, 1), c(1, state$r, state$r), c(1, p[2] + p[3], p[3]))
    modes <- c(0, 1, 1)
    for (level in 1:3) {
      surrogate <- log(setup$prior[level]) + .weights_children_log(setup, omegas[[level]], modes[level], "surrogate")
      exact <- log(setup$prior[level]) + .weights_children_log(setup, omegas[[level]], modes[level], "exact")
      expect_lt(abs(scores$surrogate[level] - surrogate), 1e-9 * max(1, abs(surrogate)))
      expect_lt(abs(scores$exact[level] - exact), 1e-9 * max(1, abs(exact)))
    }
    # the selection levels' surrogates differ from their exact densities, the
    # normal level's does not
    expect_gt(min(abs(scores$surrogate[2:3] - scores$exact[2:3])), 1e-6)
    expect_lt(abs(scores$surrogate[1] - scores$exact[1]), 1e-9 * max(1, abs(scores$exact[1])))
    rm(hand)
  }
})


test_that("the two-entry memo and absent samplers leave the draws unchanged", {

  skip_if_not_installed("rjags")
  setup <- .weights_setup(1:3)
  chain <- function(flags, ...) {
    .with_corrected_flags(flags, .block_chain(.weights_build(setup, ...), c("b", "r", "eta"), 300L,
                                              adapt = 100L, burnin = 100L))
  }
  counts <- function() .corrected_test_hook(clear = TRUE)$memo
  invisible(counts())
  one <- chain(list(memo_slots = 1))
  one_counts <- counts()
  two <- chain(list(memo_slots = 2))
  two_counts <- counts()
  none <- chain(list(memo_slots = 0))
  expect_identical(one$draws, two$draws)
  expect_identical(none$draws, two$draws)
  expect_gt(two_counts[["cluster_hits"]], one_counts[["cluster_hits"]])
  expect_identical(two_counts[["cluster_lookups"]], one_counts[["cluster_lookups"]])

  # without a coarse-capable child (the normal kernel throughout) the samplers
  # are absent and the draws are those of the stock samplers
  normal_on  <- chain(list(), mode = "0")
  normal_off <- chain(list(disable_new = TRUE), mode = "0")
  expect_identical(normal_on$draws, normal_off$draws)
})


test_that("the screened proposals and the corrected sum direction keep the indicator and the total on their laws", {

  skip_if_not_installed("rjags")
  # a quick check of the per-kind counters: every update of the screened
  # indicator either keeps its level without an exact score or makes one
  # corrected proposal with two; the gated blocks only propose while open
  # (the plain draw, Liu's step and their default mixture)
  setup <- .weights_setup(1:3)
  invisible(.corrected_test_hook(clear = TRUE))
  for (proposal in c(0, 2, 3)) {
    chain <- .with_corrected_flags(list(finite_proposal = proposal), .block_chain(
      .weights_build(setup), c("b", "eta"), 1000L, adapt = 100L, burnin = 100L))
    kinds <- .corrected_test_hook(clear = TRUE)$kinds
    expect_identical(kinds["finite", "exact_scores"], 2 * kinds["finite", "proposed"], info = paste(proposal))
    expect_lt(kinds["gated", "proposed"], kinds["gated", "updates"])
    expect_gt(kinds["finite", "accepted"], 0)
    # the eta totals keep their Gamma(3, 1) law: their mean within 5 standard errors
    total <- rowSums(chain$draws[, paste0("eta[", 1:3, "]")])
    expect_lt(abs(mean(total) - 3), 5 * sqrt(3 / coda::effectiveSize(total)), label = paste("total", proposal))
  }
})


test_that("the screened indicator's default kernel is aperiodic between two tied labels", {

  skip_if_not_installed("rjags")
  skip_if_not_installed("BayesTools")
  # The tie fixture (helper-step-row-cache.R): both labels score the same, so
  # Liu's step alone swaps them at every update and a chain thinned by two
  # keeps its initial label. The default kernel takes Liu's step (moves with
  # probability one) with probability 3/4 and the plain draw (moves with
  # probability 1/2) with 1/4: the label changes with probability q = 7/8 per
  # update and 2 q (1 - q) = 7/32 per two. A two-state chain that switches
  # with probability q per retained draw keeps each label with frequency 1/2
  # and variance (1 - q) / (4 n q) over n draws; the bound is 4 of its
  # standard deviations.
  fit <- .rows_tie_fit()
  syntax <- paste(as.character(fit[["model"]]), collapse = "\n")
  n <- 4000L
  frequency <- function(entry, thin, flags = list()) {
    invisible(gc())
    hand <- .with_corrected_flags(flags, .block_hand_model(syntax, fit[["data"]],
                                                           inits = list(bias_indicator = entry)))
    expect_identical(unname(hand$samplers["bias_indicator"]), "RoBMA::ScreenedFinite")
    rjags::adapt(hand$model, 0L, end.adaptation = TRUE)
    draws <- rjags::coda.samples(hand$model, "bias_indicator", n.iter = thin * n, thin = thin,
                                 progress.bar = "none")
    mean(as.matrix(draws[[1L]])[, 1L] == 2)
  }
  for (thin in c(2L, 1L)) {
    q <- if (thin == 1L) 7 / 8 else 2 * (7 / 8) * (1 / 8)
    bound <- 4 * sqrt((1 - q) / (4 * n * q))
    for (entry in 1:2) {
      label <- sprintf("thin %d, initial label %d", thin, entry)
      expect_lt(abs(frequency(entry, thin) - .5), bound, label = label)
    }
  }
  # Liu's step alone retains one label at even thinning, the one its phase
  # (the initial label) gives
  liu <- vapply(1:2, function(entry) frequency(entry, 2L, list(finite_proposal = 2)), numeric(1L))
  expect_identical(sort(liu), c(0, 1))
})


test_that("the screened evidence records the move probability of the kernel that ran", {

  skip_if_not_installed("rjags")
  # two levels with equal exact and surrogate scores (normal kernel, equal
  # means): the plain draw and the exact enumeration change the level with
  # probability 1/2, Liu's step with probability one, and the correction
  # accepts every proposal, so the recorded probabilities predict the changes
  syntax <- paste0("model {\n  b ~ dcat(prior[1:2])\n  m <- means[b]\n  mode <- modes[b]\n",
                   "  for (j in 1:2) { w[j] <- weights[b, j] }\n",
                   "  y ~ dselnorm_step_switch(m, 1, 1, 1, w[1:2], lo[1:2], hi[1:2], 1, 1, mode, 1)\n}")
  data <- list(prior = c(.5, .5), means = c(0, 0), modes = c(0, 0), weights = rbind(c(1, .2), c(1, .8)), y = .25,
               lo = c(0, -1e300), hi = c(1e300, 0))
  updates <- 400
  # (force_fallback is read at the end of adaptation)
  record <- function(flags) {
    invisible(gc())
    .with_corrected_flags(flags, {
      hand <- .block_hand_model(syntax, data, inits = list(b = 1L))
      expect_identical(unname(hand$samplers["b"]), "RoBMA::ScreenedFinite")
      rjags::adapt(hand$model, 0L, end.adaptation = TRUE)
    })
    invisible(.selection_sampler_report_read())
    stats::update(hand$model, n.iter = updates, progress.bar = "none")
    rows <- .selection_sampler_report_read()$screened
    expect_identical(nrow(rows), 1L)
    rows[1L, ]
  }
  liu <- record(list(finite_proposal = 2))
  expect_identical(unname(liu[c("updates", "changes", "leave_sum", "leave_max")]), c(updates, updates, updates, 1))
  for (flags in list(list(finite_proposal = 0), list(force_fallback = TRUE))) {
    half <- record(flags)
    expect_identical(unname(half[c("updates", "leave_sum", "leave_max")]), c(updates, updates / 2, .5))
    # the observed changes are their own count: binomial(400, 1/2) within 4 SD
    expect_lt(abs(half[["changes"]] - updates / 2), 4 * sqrt(updates / 4))
  }
  # the default mixture: Liu's step (1) or the plain draw (1/2) per update
  mixture <- record(list())
  expect_identical(mixture[["leave_max"]], 1)
  expect_gt(mixture[["leave_sum"]], updates / 2)
  expect_lt(mixture[["leave_sum"]], updates)
  expect_true(mixture[["leave_sum"]] * 2 == round(mixture[["leave_sum"]] * 2))
  expect_lt(abs(mixture[["changes"]] - mixture[["leave_sum"]]), 4 * sqrt(updates / 4))
})


# ---------------------------------------------------------------------------- #
# Fits (certification case 'numerical-kernels').
# ---------------------------------------------------------------------------- #

# A model-averaged fit of three Assink et al. studies with a full sampling
# covariance (rank-one cluster blocks), short.
.weights_fit <- function(parallel = FALSE, chains = 2L, seed = 1L, sample = 200L, V_rho = .6, vi = FALSE, ...) {

  data("dat.assink2016", package = "metadat", envir = environment())
  dat <- dat.assink2016[dat.assink2016$study %in% c(1, 2, 4, 9), ]
  dat[["study"]] <- as.integer(factor(dat[["study"]]))
  V <- if (vi) dat$vi else metafor::vcalc(vi, cluster = study, obs = esid, rho = V_rho, data = dat)
  suppressWarnings(RoBMA::RoBMA.mv(
    yi = yi, V = V, random = ~ 1 | study / esid, data = dat, measure = "SMD",
    chains = chains, adapt = 100L, burnin = 100L, sample = sample, seed = seed, parallel = parallel,
    silent = TRUE, convergence_checks = RoBMA::set_convergence_checks(max_Rhat = NULL, min_ESS = NULL),
    ...))
}

.weights_fit_samplers <- function(fit) unique(fit[["fit"]][["samplers"]][["Sampler"]])


test_that("model-averaged fits claim their weights and indicators, also in parallel and through update()", {

  skip_on_cran()
  skip_if_not_certification("These fits (about 60 s) run in the numerical-kernels case.")
  skip_if_not_installed("rjags")
  skip_if_not_installed("metadat")
  skip_if_not_installed("metafor")

  sampler_of <- function(fit, node) {
    table <- fit[["fit"]][["samplers"]]
    table[["Sampler"]][match(node, table[["Node"]])]
  }
  # the record counts every corrected sampler object, one per claimed node,
  # block or indicator and chain: the location block, and the samplers its
  # heterogeneity scale runs embedded (named in the scale record, not in the
  # sampler table)
  objects <- function(fit) {
    table <- fit[["fit"]][["samplers"]]
    scale <- attr(fit[["fit"]], "sampler_declarations")[["robma"]][["joint"]][["scale"]]
    embedded <- if (is.null(scale)) character(0L) else
      unlist(strsplit(scale[["decisions"]][["embedded"]], ",", fixed = TRUE))
    length(unique(table[["Index.No"]][grepl("^RoBMA::(Coarse|Gated|Screened|LocationBlock)", table[["Sampler"]])])) +
      sum(grepl("^RoBMA::(Coarse|Gated|Screened)", embedded))
  }
  bias_side <- c("bias_indicator", "PET_1", "PEESE_1")
  allocation_sd <- "mu__xRE_ALLOCx_heterogeneity__allocation_sd"

  # the default policy: the location block owns the bias label, the PET and
  # PEESE slopes and (through its heterogeneity scale) the allocation SD; the
  # weights keep the gated block and the heterogeneity indicator the screened
  # indicator
  serial <- .weights_fit()
  expect_identical(sampler_of(serial, c(bias_side, allocation_sd)), rep("RoBMA::LocationBlock", 4L))
  expect_identical(sampler_of(serial, "eta_component_3[1]"), "RoBMA::GatedCorrectedBlock")
  expect_identical(sampler_of(serial, "mu__xRE_ALLOCx_heterogeneity__include_component_1_indicator"),
                   "RoBMA::ScreenedFinite")
  record <- selection_sampler_info(fit = serial)
  expect_identical(record[["samplers"]], 2 * objects(serial))
  expect_true(record[["complete"]])

  extended <- suppressWarnings(stats::update(serial, sample_extend = 50L))
  expect_true(all(c("RoBMA::ScreenedFinite", "RoBMA::GatedCorrectedBlock", "RoBMA::LocationBlock") %in%
                    .weights_fit_samplers(extended)))
  expect_identical(selection_sampler_info(fit = extended)[["samplers"]], 4 * objects(serial))

  parallel <- .weights_fit(parallel = TRUE)
  for (worker in parallel[["fit"]][["samplers_by_worker"]]) {
    expect_true(all(c("RoBMA::ScreenedFinite", "RoBMA::GatedCorrectedBlock", "RoBMA::LocationBlock") %in%
                      worker[["samplers"]][["Sampler"]]))
  }
  expect_identical(selection_sampler_info(fit = parallel)[["samplers"]], 2 * objects(serial))

  # without joint samplers the screened indicator owns the bias label, the
  # gated slices the slopes and the corrected block the allocation SD
  options <- RoBMA.options()["sampler.joint"]
  on.exit(do.call(RoBMA.options, options), add = TRUE)
  RoBMA.options(sampler.joint = "off")
  ordinary <- .weights_fit()
  expect_identical(sampler_of(ordinary, c(bias_side, allocation_sd)),
                   c("RoBMA::ScreenedFinite", "BayesTools::GatedSlice", "BayesTools::GatedSlice",
                     "RoBMA::CoarseCorrectedBlock"))
  expect_identical(sampler_of(ordinary, "eta_component_3[1]"), "RoBMA::GatedCorrectedBlock")
  record <- selection_sampler_info(fit = ordinary)
  expect_identical(record[["samplers"]], 2 * objects(ordinary))
  expect_true(record[["complete"]])
  extended <- suppressWarnings(stats::update(ordinary, sample_extend = 50L))
  expect_identical(sampler_of(extended, bias_side), sampler_of(ordinary, bias_side))
  expect_identical(selection_sampler_info(fit = extended)[["samplers"]], 4 * objects(ordinary))
})


test_that("fits without coarse-capable children, and the memo alone, draw as before", {

  skip_on_cran()
  skip_if_not_certification("These fits (about 60 s) run in the numerical-kernels case.")
  skip_if_not_installed("rjags")
  skip_if_not_installed("metadat")
  skip_if_not_installed("metafor")
  draws <- function(fit) coda::as.mcmc.list(fit[["fit"]])

  # V = vi: every row is a singleton switch node and no child has a coarse
  # rung. The row cache gives the weights to the gated block
  # (test-00-selection-step-row-cache.R); without it neither the samplers nor
  # the gamma coordinates enter
  vi_rows <- .weights_fit(vi = TRUE)
  vi_new  <- .with_row_flags(list(disable = TRUE), .weights_fit(vi = TRUE))
  vi_old  <- .with_corrected_flags(list(disable_new = TRUE), .weights_fit(vi = TRUE))
  expect_true("RoBMA::GatedCorrectedBlock" %in% .weights_fit_samplers(vi_rows))
  expect_false(any(grepl("^RoBMA::(Gated|Screened)", .weights_fit_samplers(vi_new))))
  expect_identical(draws(vi_new), draws(vi_old))

  # a full sampling covariance with the samplers of this round off: one or two
  # memo entries draw the same
  one <- .with_corrected_flags(list(disable_new = TRUE, memo_slots = 1), .weights_fit())
  two <- .with_corrected_flags(list(disable_new = TRUE, memo_slots = 2), .weights_fit())
  expect_identical(draws(one), draws(two))
})
