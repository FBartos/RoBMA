context("Corrected block sampler: assignment, kernel decision, rollback and record")
source(testthat::test_path("common-functions.R"))
skip_on_cran()

# A block of parameters declared by BayesTools' dbt_block markers with a
# coarse-capable selection child gets RoBMA::CoarseCorrectedBlock: a factor slice
# sweep on the surrogate in a random direction order, corrected by one
# Metropolis step against the full conditional. Hand-written models here (the
# helpers are in helper-selection-corrected-block.R); small fits run in the
# certification case 'numerical-kernels', the exactness fixtures in
# test-00-selection-corrected-block-exactness.R.

.block_sampler <- "RoBMA::CoarseCorrectedBlock"

# A block of mu and nu on one crude cluster child of mu + nu, with further
# `markers` and model `extra` lines.
.block_pair_model <- function(markers = c("dbt_block(mu, 1)", "dbt_block(nu, 1)"), extra = "",
                              data = list(), inits = list(), vector_rule = 0, ...) {

  marker_names <- sprintf("marker_%d", seq_along(markers))
  syntax <- paste0("model {\n  mu ~ dt(0, 1, 4)\n  nu ~ dt(0, 1, 4)\n  for (j in 1:4) { means[j] <- mu + nu }\n",
    .block_cluster_child("y[1:4]", "means[1:4]", 4L, vector_rule = vector_rule), "\n", extra, "\n",
    if (length(markers)) paste0("  ", marker_names, " ~ ", markers, collapse = "\n"), "\n}")
  data <- c(.block_cluster_data(), list(y = rep(0, 4)),
            stats::setNames(as.list(rep(0, length(markers))), marker_names), data)
  .block_hand_model(syntax, data, inits = c(list(mu = -.25, nu = .1), inits), ...)
}

# The samplers of nodes, by node name.
.block_assigned <- function(hand, nodes) unname(hand$samplers[nodes])


test_that("a declared block with a coarse-capable selection child gets the corrected block sampler", {

  skip_if_not_installed("rjags")
  coarse <- "RoBMA::CoarseCorrectedSlice"
  factor <- "BayesTools::FactorSlice"

  # the RoBMA factories are asked before BayesTools' (RoBMA's module loads later)
  factories <- rjags::list.factories("sampler")[["factory"]]
  expect_identical(factories[1:3], c(.block_sampler, coarse, factor))

  expect_identical(.block_assigned(.block_pair_model(), c("mu", "nu")), rep(.block_sampler, 2L))
  # the block factory follows the corrected sampler's configuration
  for (sampler in c("slice", "jags.default")) {
    expect_identical(.block_assigned(.block_pair_model(sampler = sampler), c("mu", "nu")), rep(factor, 2L),
                     info = sampler)
  }
  # a block whose children have no surrogate evaluator stays with BayesTools
  plain <- .block_hand_model(paste0("model {\n  x1 ~ dnorm(0, 1)\n  x2 ~ dnorm(x1, 1)\n  z ~ dnorm(x1 + x2, 1)\n",
                                    .block_markers(c("x1", "x2"))$syntax, "\n}"),
                             c(list(z = .3), .block_markers(c("x1", "x2"))$data))
  expect_identical(.block_assigned(plain, c("x1", "x2")), rep(factor, 2L))
  # a fixed best-weight vector rule leaves the cluster child without a coarse rule
  expect_identical(.block_assigned(.block_pair_model(vector_rule = 1), c("mu", "nu")), rep(factor, 2L))

  # several ids: each group with a selection child is its own block
  two <- .block_pair_model(
    markers = c("dbt_block(mu, 1)", "dbt_block(nu, 1)", "dbt_block(la, 7)", "dbt_block(lb, 7)"),
    extra = paste0("  la ~ dt(0, 1, 4)\n  lb ~ dt(0, 1, 4)\n  for (j in 1:4) { means2[j] <- la + lb }\n",
                   .block_cluster_child("y2[1:4]", "means2[1:4]", 4L)),
    data = list(y2 = rep(0, 4)), inits = list(la = .1, lb = -.1))
  expect_identical(.block_assigned(two, c("mu", "nu", "la", "lb")), rep(.block_sampler, 4L))
  expect_identical(sum(names(rjags::list.samplers(two$model)) == .block_sampler), 2L)

  # a block mixing a member without a selection child is claimed whole
  mixed <- .block_pair_model(markers = c("dbt_block(mu, 1)", "dbt_block(nu, 1)", "dbt_block(x, 1)"),
                             extra = "  x ~ dnorm(0, 1)\n  z ~ dnorm(x, 1)", data = list(z = .4), inits = list(x = 0))
  expect_identical(.block_assigned(mixed, c("mu", "nu", "x")), rep(.block_sampler, 3L))
  expect_identical(sum(names(rjags::list.samplers(mixed$model)) == .block_sampler), 1L)

  # a gated member leaves the group (the gate takes precedence); the stranded
  # singleton is declined by both block factories and, being declared, by the
  # scalar corrected factory: it keeps the stock sampler
  gated <- .block_pair_model(markers = c("dbt_block(mu, 1)", "dbt_block(nu, 1)", "dbt_gate(nu, 1)"))
  expect_identical(.block_assigned(gated, c("mu", "nu")), c("base::RealSlicer", "BayesTools::GatedSlice"))
  single <- .block_pair_model(markers = "dbt_block(mu, 3)")
  expect_identical(.block_assigned(single, c("mu", "nu")), c("base::RealSlicer", coarse))
})


test_that("with BayesTools' module loaded after RoBMA's, BayesTools' factor slice sampler takes the block", {

  skip_if_not_installed("rjags")
  BayesTools::BayesTools_load_JAGS_module(quiet = TRUE)
  gc()
  on.exit({
    gc()
    rjags::unload.module("RoBMA", quiet = TRUE)
    RoBMA:::.load_RoBMA_module()
  }, add = TRUE)
  rjags::unload.module("BayesTools", quiet = TRUE)
  BayesTools::BayesTools_load_JAGS_module(quiet = TRUE)
  factories <- rjags::list.factories("sampler")[["factory"]]
  expect_lt(match("BayesTools::FactorSlice", factories), match(.block_sampler, factories))
  expect_identical(.block_assigned(.block_pair_model(), c("mu", "nu")), rep("BayesTools::FactorSlice", 2L))
})


test_that("the module order is restored for the remaining tests", {

  skip_if_not_installed("rjags")
  factories <- rjags::list.factories("sampler")[["factory"]]
  expect_lt(match(.block_sampler, factories), match("BayesTools::FactorSlice", factories))
})


test_that("after a fallback the block runs BayesTools' factor slice update, draw for draw", {

  skip_if_not_installed("rjags")
  # no adaptive update and the fallback on: rule (a) switches the block to the
  # exact kernel at once, the same calls in the same order as
  # BayesTools::FactorSlice with unadapted axes
  fallback <- .block_chain(.block_pair_model(min_acceptance = .9), c("mu", "nu"), 300L, adapt = 0L, burnin = 0L)
  factor   <- .block_chain(.block_pair_model(sampler = "slice"), c("mu", "nu"), 300L, adapt = 0L, burnin = 0L)
  expect_identical(unname(fallback$samplers[c("mu", "nu")]), rep(.block_sampler, 2L))
  expect_identical(unname(factor$samplers[c("mu", "nu")]), rep("BayesTools::FactorSlice", 2L))
  expect_identical(RoBMA:::.selection_sampler_report_read()[["fallback_few_proposals"]], 0)
  expect_identical(fallback$report[["sampling_proposed"]], 0)
  expect_identical(fallback$counters[["proposed"]], 0)
  expect_gt(length(unique(fallback$draws[, "mu"])), 100L)
  expect_identical(fallback$draws, factor$draws)
})


# The kernel decision of one block of two members: the report after `adapt`
# adaptive iterations, burn-in and draws.
.block_decision <- function(threshold, adapt, ...) {

  .block_chain(.block_pair_model(min_acceptance = threshold, ...), "mu", 100L, adapt = adapt, burnin = 100L)$report
}

.block_switched <- function(reason) {

  report <- list(samplers = 1, fallback_few_proposals = 0, fallback_surrogate_exact = 0,
                 fallback_low_acceptance = 0, sampling_proposed = 0, sampling_accepted = 0,
                 sampling_judged = 0, min_sampler_rate = NA_real_)
  report[[reason]] <- 1
  report
}


test_that("a block decides its kernel once, at the end of adaptation, judging rate^(1/d)", {

  skip_if_not_installed("rjags")
  # (a) fewer than 50 adaptive proposals
  tiny <- 1e-9
  for (adapt in c(0L, 49L)) {
    expect_identical(.block_decision(tiny, adapt), .block_switched("fallback_few_proposals"), info = adapt)
  }
  kept <- .block_decision(tiny, 50L)
  expect_identical(unlist(kept[c("samplers", "fallback_few_proposals", "fallback_surrogate_exact",
                                 "fallback_low_acceptance", "sampling_proposed")]),
                   c(samplers = 1, fallback_few_proposals = 0, fallback_surrogate_exact = 0,
                     fallback_low_acceptance = 0, sampling_proposed = 200))

  # (b) the adaptive rate of the block (the process counters count adaptation;
  # the adaptation does not depend on the threshold) judged as rate^(1/2)
  hand <- .block_pair_model()
  RoBMA::selection_sampler_info(clear = TRUE)
  rjags::adapt(hand$model, 300L, end.adaptation = FALSE)
  adaptive <- RoBMA::selection_sampler_info()
  expect_identical(adaptive[["proposed"]], 300)
  rate <- adaptive[["accepted"]] / adaptive[["proposed"]]
  expect_true(rate > .1 && rate < .9)
  per_direction <- sqrt(rate)
  # between the block rate and its per-direction equivalent: kept
  expect_identical(.block_decision((rate + per_direction) / 2, 300L)[["fallback_low_acceptance"]], 0)
  expect_identical(.block_decision(per_direction * (1 - 1e-12), 300L)[["fallback_low_acceptance"]], 0)
  expect_identical(.block_decision(per_direction * (1 + 1e-12), 300L), .block_switched("fallback_low_acceptance"))

  # (c) a surrogate that is the exact density in every state: a vector rule that
  # is a best-weight rule for both values of an indicator
  exact <- .block_decision(.9, 300L, vector_rule = "vr", extra = "  ind ~ dcat(c(.5, .5))\n  vr <- ind",
                           inits = list(ind = 1))
  expect_identical(exact, .block_switched("fallback_surrogate_exact"))
  # (a) comes before (c)
  expect_identical(.block_decision(.9, 49L, vector_rule = "vr", extra = "  ind ~ dcat(c(.5, .5))\n  vr <- ind",
                                   inits = list(ind = 1)), .block_switched("fallback_few_proposals"))

  # a later adaptOff() changes nothing: the block that switched counts once
  hand <- .block_pair_model(min_acceptance = .9)
  rjags::adapt(hand$model, 300L, end.adaptation = FALSE)
  expect_true(rjags::adapt(hand$model, 0L, end.adaptation = TRUE))
  expect_identical(.Call("RoBMA_selnorm_sampler_report", FALSE, PACKAGE = "RoBMA"),
                   .block_switched("fallback_low_acceptance"))
  .Call("adapt_off", hand$model$ptr(), PACKAGE = "rjags")
  stats::update(hand$model, n.iter = 100L, progress.bar = "none")
  expect_identical(RoBMA:::.selection_sampler_report_read(), .block_switched("fallback_low_acceptance"))
})


test_that("a block started where the surrogate is poor learns its axes and samples the posterior", {

  skip_if_not_installed("rjags")
  skip_if_not_installed("coda")
  # While adapting, a rejected correction is followed by one exact sweep, so the
  # chain moves from a start where the corrected block rejects most proposals
  # and the axes learn from moving states: the widths it learns are those of an
  # exact sampler at stationarity within a factor of 3, with the fallback off
  # (the corrected kernel after adaptation) and on (its adaptive rate is below
  # 0.90, so the exact sweep after adaptation), and the draws agree with the
  # reference.
  for (min_acceptance in c(0, .9)) {
    label  <- paste("min_acceptance", min_acceptance)
    result <- .block_poor_start(min_acceptance)
    expect_identical(unname(result$samplers[c("a", "b")]), rep(.block_sampler, 2L), info = label)
    expect_lt(result$adaptive_rate, .9, label = label)
    expect_identical(result$report[["fallback_low_acceptance"]], as.numeric(min_acceptance > 0), info = label)
    # the chain moves: the block left its entry state in at least 90% of the
    # last 200 adaptive updates (a chain frozen at the start moved in a few)
    expect_identical(result$adaptive_updates, 500L, info = label)
    expect_gte(result$moved_late, .9, label = label)
    expect_length(result$widths, 2L)
    ratio <- result$widths / result$reference_widths
    expect_true(all(ratio > 1 / 3 & ratio < 3), info = paste(label, paste(round(ratio, 3), collapse = " ")))
    comparison <- result$comparison[result$comparison$observable %in% c("a", "b", "var_a", "var_b", "ab"), ]
    expect_true(all(comparison$ess >= 150), info = paste(label, paste(round(comparison$ess), collapse = " ")))
    expect_true(all(abs(comparison$estimate - comparison$reference) <= 4.5 * comparison$mcse + 2 * comparison$allowance),
                info = paste(label, paste(round(comparison$z, 2), collapse = " ")))
  }
})


test_that("a rejection restores the block-entry values and an exhausted interval the direction-entry values", {

  skip_if_not_installed("rjags")
  # rejections: every rejected proposal repeats the previous draw exactly, and
  # no accepted one does (no interval was exhausted)
  chain <- .block_chain(.block_pair_model(), c("mu", "nu"), 400L, adapt = 0L, burnin = 0L)
  repeated <- sum(vapply(2:400, function(i) identical(chain$draws[i, ], chain$draws[i - 1L, ]), logical(1L)))
  expect_identical(chain$exhausted, 0)
  expect_gt(chain$counters[["correction_rejected"]], 50)
  # the first draw follows the initial state, which is not monitored
  expect_true(repeated %in% (chain$counters[["correction_rejected"]] - 0:1))

  # exhaustion: a constant term of -5e17 makes every point of every slice equal
  # to its level in floating point, so every interval is exhausted and each
  # direction must restore its entry values; g's transform log(g) does not return
  # g exactly through exp(), so a restore through the inverse transform would move it
  candidates <- seq(.501, 3, by = .001)
  g0 <- candidates[which(exp(log(candidates)) != candidates)[1L]]
  hand <- .block_pair_model(markers = c("dbt_block(mu, 1)", "dbt_block(nu, 1)", "dbt_block(g, 1)"),
                            extra = "  g ~ dgamma(2, 1)\n  big ~ dnorm(0 * mu + 0 * nu + 0 * g, 1)",
                            data = list(big = 1e9), inits = list(g = g0))
  expect_identical(.block_assigned(hand, c("mu", "nu", "g")), rep(.block_sampler, 3L))
  chain <- .block_chain(hand, c("mu", "nu", "g"), 20L, adapt = 0L, burnin = 0L)
  expect_identical(chain$exhausted, 60)
  expect_identical(unname(chain$draws[, "g"]), rep(g0, 20L))
  expect_identical(unname(chain$draws[, "mu"]), rep(-.25, 20L))
  expect_identical(unname(chain$draws[, "nu"]), rep(.1, 20L))
  expect_identical(chain$counters[["proposed"]], 20)
})


test_that("a member whose two bounds are too wide apart stops the fit, a wide one-sided bound does not", {

  skip_if_not_installed("rjags")
  members <- c("dbt_block(mu, 1)", "dbt_block(nu, 1)", "dbt_block(p, 1)")
  # p ~ N(0, 1) above -1e6: u = log(p + 1e6) holds p to far finer than its
  # posterior scale
  wide <- .block_pair_model(markers = members, extra = "  p ~ dnorm(0, 1) T(-1e6, )", inits = list(p = .3))
  expect_identical(.block_assigned(wide, c("mu", "nu", "p")), rep(.block_sampler, 3L))
  chain <- .block_chain(wide, c("mu", "nu", "p"), 200L, adapt = 300L, burnin = 0L)
  expect_gt(length(unique(chain$draws[, "p"])), 50L)
  expect_true(all(abs(chain$draws[, "p"]) < 6))
  # on (-1e308, 1e308) the span upper - lower is not finite: the block sampler
  # refuses the member when it is built, and BayesTools' sampler under "slice"
  span <- "  p ~ dnorm(0, 1) T(-1e308, 1e308)"
  for (sampler in c("slice.coarse", "slice")) {
    message <- tryCatch({
      .block_pair_model(markers = members, extra = span, inits = list(p = .3), sampler = sampler)
      NA_character_
    }, error = function(error) conditionMessage(error))
    expected <- if (sampler == "slice") "Factor slice sampler" else "Corrected slice (block)"
    expect_true(grepl(paste0(expected, ": the support bounds of this node are too wide for the sampler's scale"),
                      message, fixed = TRUE), info = sampler)
    # the hint undeclares the blocks
    expect_true(grepl("RoBMA.options(sampler.blocks = FALSE)", .sampler_error_hint(message), fixed = TRUE),
                info = sampler)
  }
})


test_that("the record counts a block once per chain and model build, with the scalar samplers", {

  skip_if_not_installed("rjags")
  for (chains in 1:2) {
    .block_pair_model(chains = chains)
    expect_identical(RoBMA:::.selection_sampler_report_read()[["samplers"]], as.numeric(chains), info = chains)
  }
  # a block and a scalar corrected node
  .block_pair_model(chains = 2L, extra = paste0("  kappa ~ dt(0, 1, 4)\n  for (j in 1:4) { means3[j] <- kappa }\n",
                                                .block_cluster_child("y3[1:4]", "means3[1:4]", 4L)),
                    data = list(y3 = rep(0, 4)), inits = list(kappa = 0))
  expect_identical(RoBMA:::.selection_sampler_report_read()[["samplers"]], 4)

  # the judged acceptances of one block of two over its window: its proposals
  # times its rate per direction, rate^(1/2)
  chain  <- .block_chain(.block_pair_model(), "mu", 600L, adapt = 0L, burnin = 0L)
  report <- chain$report
  expect_identical(report[["sampling_proposed"]], 600)
  expect_lt(report[["sampling_accepted"]], 600)
  expected <- 600 * sqrt(report[["sampling_accepted"]] / 600)
  expect_lt(abs(report[["sampling_judged"]] - expected), 1e-9 * expected)
})


test_that("errors of the block sampler say 'Corrected slice' and their hint recovers the fit", {

  skip_if_not_installed("rjags")
  # a member on a bound of its support
  bounded <- function(markers, p0, ...) {
    .block_pair_model(markers = markers, extra = "  p ~ dunif(0, 1)", inits = list(p = p0), ...)
  }
  declared <- c("dbt_block(mu, 1)", "dbt_block(p, 1)")
  hand <- bounded(declared, 0)
  expect_identical(.block_assigned(hand, c("mu", "p")), rep(.block_sampler, 2L))
  message <- tryCatch({
    stats::update(hand$model, n.iter = 1L, progress.bar = "none")
    NA_character_
  }, error = function(error) conditionMessage(error))
  expect_true(grepl("Corrected slice (block): current value lies on a bound of its support", message, fixed = TRUE))
  hinted <- .sampler_error_hint(message)
  expect_true(grepl("A start strictly inside the support of the parameter avoids this error", hinted, fixed = TRUE))
  expect_true(grepl("RoBMA.options(sampler.blocks = FALSE)", hinted, fixed = TRUE))
  # under "slice" BayesTools' sampler takes the block and stops the same way:
  # the hint is the same
  slice <- bounded(declared, 0, sampler = "slice")
  message <- tryCatch({
    stats::update(slice$model, n.iter = 1L, progress.bar = "none")
    NA_character_
  }, error = function(error) conditionMessage(error))
  expect_true(grepl("Factor slice sampler: current value lies on a bound of its support", message, fixed = TRUE))
  expect_true(grepl("RoBMA.options(sampler.blocks = FALSE)", .sampler_error_hint(message), fixed = TRUE))
  # both remedies recover the fit: the model without the block declaration (what
  # sampler.blocks = FALSE fits), and a start strictly inside the support
  for (recovered in list(bounded(character(), 0), bounded(declared, .5))) {
    stats::update(recovered$model, n.iter = 20L, progress.bar = "none")
    values <- rjags::coda.samples(recovered$model, "p", n.iter = 20L, progress.bar = "none")
    expect_true(all(as.numeric(values[[1L]]) > 0 & as.numeric(values[[1L]]) < 1))
  }
})


# ---------------------------------------------------------------------------- #
# Fits (certification case 'numerical-kernels').
# ---------------------------------------------------------------------------- #

.block_fit_data <- function(studies = c(1, 2, 4, 12, 13)) {

  data("dat.assink2016", package = "metadat", envir = environment())
  dat <- dat.assink2016[dat.assink2016$study %in% studies, ]
  dat[["study"]] <- as.integer(factor(dat[["study"]]))
  V <- metafor::vcalc(vi, cluster = study, type = deltype, obs = esid, rho = c(.7, .5), data = dat)
  list(data = dat, V = V)
}

.block_fit <- function(options = list(), chains = 1L, parallel = FALSE) {

  previous <- RoBMA::RoBMA.options()[c("sampler.blocks", "sampler.gates", "selection.sampler")]
  on.exit(do.call(RoBMA::RoBMA.options, previous), add = TRUE)
  if (length(options)) do.call(RoBMA::RoBMA.options, options)
  setup <- .block_fit_data()
  suppressWarnings(RoBMA::bselmodel.mv(yi = yi, V = setup$V, random = ~ 1 | study / esid, data = setup$data,
    measure = "SMD", chains = chains, adapt = 100L, burnin = 100L, sample = 100L, seed = 1L, parallel = parallel,
    silent = TRUE, convergence_checks = RoBMA::set_convergence_checks(max_Rhat = NULL, min_ESS = NULL)))
}


test_that("a selection model with a full sampling covariance samples its declared block with the block sampler", {

  skip_if_not_certification("These fits (about 20 s) run in the numerical-kernels case.")
  skip_if_not_installed("rjags")
  skip_if_not_installed("metadat")
  skip_if_not_installed("metafor")

  default  <- .block_fit()
  samplers <- default[["fit"]][["samplers"]]
  blocks   <- attr(default[["fit"]], "sampler_declarations", exact = TRUE)[["blocks"]]
  expect_null(default[["fit"]][["errors"]])
  expect_identical(sum(samplers[["Sampler"]] == .block_sampler), nrow(blocks))
  expect_false(any(samplers[["Sampler"]] == "BayesTools::FactorSlice"))
  record <- RoBMA::selection_sampler_info(fit = default)
  expect_identical(record[["samplers"]], 1 + sum(samplers[["Sampler"]] == "RoBMA::CoarseCorrectedSlice"))

  # the other options keep their samplers
  slice <- .block_fit(list(selection.sampler = "slice"))
  expect_identical(sum(slice[["fit"]][["samplers"]][["Sampler"]] == "BayesTools::FactorSlice"), nrow(blocks))
  expect_identical(slice[["fit"]][["samplers"]][["Node"]], samplers[["Node"]])
  for (options in list(list(selection.sampler = "jags.default"), list(sampler.blocks = FALSE))) {
    other <- .block_fit(options)
    expect_false(any(other[["fit"]][["samplers"]][["Sampler"]] %in% c(.block_sampler, "BayesTools::FactorSlice")),
                 info = names(options))
  }

  # chains built in parallel workers (which load BayesTools' module, then RoBMA's)
  parallel <- .block_fit(chains = 2L, parallel = TRUE)
  expect_null(parallel[["fit"]][["errors"]])
  expect_identical(sum(parallel[["fit"]][["samplers"]][["Sampler"]] == .block_sampler), nrow(blocks))
  expect_identical(RoBMA::selection_sampler_info(fit = parallel)[["samplers"]], 2 * record[["samplers"]])
})
