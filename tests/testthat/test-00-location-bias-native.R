context("Location block sampler: the bias label")
source(testthat::test_path("common-functions.R"))
skip_on_cran()

# With a bias mixture of none, step-selection, PET and PEESE labels the
# location block also owns the bias indicator and the slope sources
# (src/samplers/LocationBlockBias.h): dormant slopes are drawn from their
# priors (D), the label is updated by Liu's Metropolized Gibbs step on its
# surrogate full conditional (L0), and the sweep keeps one exact correction.
# The fixture and its quadrature reference are in helper-location-bias.R.


test_that("the location block owns the bias label, its slopes and the nested node's label and declines what it cannot verify", {

  skip_if_not_installed("rjags")
  members <- c("bias_indicator", "PET_1", "PEESE_1", "g[1]", "g[2]", "slab[1]", "slab[2]", paste0("Z[", 1:3, "]"))
  for (nested in c(FALSE, TRUE)) {
    setup <- .bias_setup(nested = nested)
    model <- .bias_model(setup)
    chain <- .location_hand_chain(model, c("bias_indicator", "PET", "PEESE"), draws = 2000L, burnin = 200L,
                                  adapt = 200L)
    expect_identical(chain$decision, "claimed with the bias label", info = paste("nested", nested))
    expect_true(all(chain$samplers[members] == "RoBMA::LocationBlock"), info = paste("nested", nested))
    # the nested weight node is a child, not a member
    if (nested) expect_false(identical(unname(chain$samplers["a"]), "RoBMA::LocationBlock"))
    counters <- chain$counters
    expect_gt(counters[["bias_label_accepted"]], 100)
    expect_gt(counters[["bias_dormant_proposed"]], 1000)
    expect_identical(counters[["bias_dormant_accepted"]], counters[["bias_dormant_proposed"]])
    # the swap L2 is attempted in a quarter of the sweeps (a coin independent of the state)
    expect_lt(abs(counters[["bias_label_swap_proposed"]] / counters[["sweeps"]] - 0.25), 0.05)
    # every label is visited and the physical slopes are zero outside their labels
    expect_setequal(unique(chain$draws[, "bias_indicator"]), 1:4)
    expect_true(all(chain$draws[chain$draws[, "bias_indicator"] != 3, "PET"] == 0))
    expect_true(all(chain$draws[chain$draws[, "bias_indicator"] != 4, "PEESE"] == 0))
    # the exact correction is one up to rounding whenever both ends are Gaussian
    expect_lt(counters[["max_correction_gaussian"]], 1e-10)
  }

  # a two-sided slope prior keeps the bias nodes with their ordinary samplers
  # and the location members with the block
  setup <- .bias_setup()
  model <- .bias_model(setup)
  model$syntax <- sub("PET_1 ~ dt(0, pet_precision, 1) T(0, )", "PET_1 ~ dnorm(0, pet_precision)", model$syntax,
                      fixed = TRUE)
  chain <- .location_hand_chain(model, "bias_indicator", draws = 200L, burnin = 50L, adapt = 50L)
  expect_identical(chain$decision, "claimed; bias label declined: bias: slope source law or support")
  expect_true(all(chain$samplers[c("g[1]", "g[2]", "Z[1]")] == "RoBMA::LocationBlock"))
  expect_false(any(chain$samplers[c("bias_indicator", "PET_1", "PEESE_1")] == "RoBMA::LocationBlock"))

  # a slope entering a mean through anything but its chart declines the label
  model <- .bias_model(setup)
  model$syntax <- sub("+ PET * pet_col[1]", "+ PET * pet_col[1] + 0.01 * PET_1", model$syntax, fixed = TRUE)
  chain <- .location_hand_chain(model, "bias_indicator", draws = 200L, burnin = 50L, adapt = 50L)
  expect_identical(chain$decision, "claimed; bias label declined: structure: slope source outside its chart")
  expect_false(any(chain$samplers[c("bias_indicator", "PET_1")] == "RoBMA::LocationBlock"))
})


test_that("public fits give the bias label to the location block unless its probabilities are read", {

  skip_if_not_installed("rjags")
  options <- RoBMA.options()[c("sampler.joint", "sampler.probabilities", "sampler.nested.weights")]
  on.exit(do.call(RoBMA.options, options), add = TRUE)
  dat <- data.frame(
    study = rep(1:8, each = 3L),
    x     = rep(c(-1, 0, 1), 8L) + rep(seq(-0.2, 0.2, length.out = 8L), each = 3L),
    yi    = 0.15 + 0.3 * (rep(c(-1, 0, 1), 8L) + rep(seq(-0.2, 0.2, length.out = 8L), each = 3L)) +
      rep(c(0.03, -0.02, 0.05, -0.04, 0.01, 0.02, -0.03, 0.04), 3L),
    vi    = rep(c(0.02, 0.03, 0.025), 8L)
  )
  fit <- function() {
    .location_one_chain(RoBMA.mv(yi = yi, V = vi, mods = ~ x, random = ~ 1 | study, data = dat, measure = "SMD",
                                 adapt = 100L, burnin = 50L, sample = 100L, chains = 1L, seed = 4L, autofit = FALSE,
                                 silent = TRUE))
  }
  owned <- function(fit) {
    table <- fit[["fit"]][["samplers"]]
    table[["Node"]][table[["Sampler"]] == "RoBMA::LocationBlock"]
  }
  bias <- c("bias_indicator", "PET_1", "PEESE_1")
  RoBMA.options(sampler.joint = "auto", sampler.probabilities = character())
  for (nested in c(TRUE, FALSE)) {
    RoBMA.options(sampler.nested.weights = nested)
    fitted <- fit()
    expect_true(all(bias %in% owned(fitted)), info = paste("nested", nested))
    expect_true(all(c("mu_intercept_indicator", "mu_x_indicator") %in% owned(fitted)), info = paste("nested", nested))
    record <- attr(fitted[["fit"]], "sampler_declarations", exact = TRUE)[["robma"]][["joint"]]
    expect_identical(record[["status"]], "active")
    expect_gt(record[["counters"]][["bias_label_proposed"]], 0)
  }
  # an exact probability readout of the indicator keeps its enumerating sampler
  RoBMA.options(sampler.nested.weights = TRUE, sampler.probabilities = "bias_indicator")
  readout <- fit()
  expect_false(any(bias %in% owned(readout)))
  expect_true("mu_intercept_indicator" %in% owned(readout))
  expect_identical(attr(readout[["fit"]], "sampler_declarations", exact = TRUE)[["robma"]][["probabilities"]][["status"]],
                   "active")
})


test_that("frozen-state identities of the Gaussian-family moves hold for regular and prior cells", {

  skip_if_not_installed("rjags")
  # the L1 probe records the analytic log ratio F_c'(theta') / F_c(theta)
  # (stage1) and the ratio computed from the inner target and the proposal
  # densities at both states (full); the S probe the slice move's acceptance
  # ratio, one when the collapsed Gaussian algebra is exact
  identity <- function(setup, draws = 1500L) {
    chain <- .location_hand_chain(.bias_model(setup), "bias_indicator", draws = draws, burnin = 100L, adapt = 100L,
                                  probe = TRUE)
    expect_true(grepl("with the bias label", chain$decision), info = chain$decision)
    probe <- chain$probe
    family <- probe[probe$type == "family" & is.finite(probe$stage1), ]
    slope  <- probe[probe$type == "slope", ]
    expect_gt(nrow(family), 50L)
    expect_lt(max(abs(family$full - family$stage1) / pmax(1, abs(family$stage1))), 1e-10)
    if (nrow(slope) > 0L) expect_lt(max(abs(slope$full)), 1e-10)
    list(family = family, slope = slope)
  }
  # the fixture: gated intercept with slab mean 0.2, half-Cauchy PET of scale
  # 0.5, Gamma(2, 2) PEESE
  base <- identity(.bias_setup())
  expect_gt(nrow(base$slope), 10L)
  # near-equal standard errors (the PET slope nearly collinear with the
  # intercept) and a PEESE slope that enters no mean (not identified: prior
  # cells): moves between regular and prior cells in both directions
  flat <- .bias_setup()
  flat$se <- .21 + 1e-4 * c(1, -1, 2, 0, -2, 1, 0, -1, 2, 1)
  flat$residual_sd <- flat$se * sqrt(1 - .85)
  flat$loading <- flat$se * sqrt(.85)
  flat$peese_column <- rep(0, 10L)
  cells <- identity(flat)$family
  expect_gt(sum(cells$selection_part %in% c(1, 3)), 0L)   # from a prior cell
  expect_gt(sum(cells$selection_part %in% c(2, 3)), 0L)   # to a prior cell
  expect_gt(sum(cells$selection_part == 0), 0L)
  # an ungated intercept, no PEESE label and a Gamma(2, 1) PET prior (zero
  # density at zero) with data that pull the slope below zero
  ungated <- .bias_setup()
  ungated$gated_intercept <- FALSE
  ungated$peese <- FALSE
  ungated$pet_prior <- "dgamma(2, 1)"
  ungated$y <- 0.6 - 2 * ungated$se + c(.02, -.01, .03, -.02, .01, .02, -.03, .01, .02, -.01)
  identity(ungated)
})


test_that("frozen-state identities of the cross-family swap and of the label scores hold", {

  skip_if_not_installed("rjags")
  # the probe dry-runs L2 at its slot of every sweep (in reversed sweeps after
  # the location components moved the state): from the proposal it rebuilds
  # the shears, E_0, every label's evidences and the label proposals, the
  # reverse move's stage-1 and full ratios and the accepted probability
  # flows of both directions; at the base point it compares the fitted
  # quadratic of the Gaussian plane with the children's own densities; and at
  # L0's slot it compares every label's exact score from the block's buffers
  # with the graph's densities of the working state
  for (variant in list(list(), list(nested = TRUE), list(second_selection = TRUE))) {
    setup <- do.call(.bias_setup, variant)
    chain <- .location_hand_chain(.bias_model(setup), "bias_indicator", draws = 1200L, burnin = 100L,
                                  adapt = 100L, probe = TRUE)
    label <- paste(names(variant), collapse = " ")
    probe <- chain$probe
    swap  <- probe[probe$type == "label_swap", ]
    reciprocal <- swap[swap$kind == 0, ]
    rebuilt    <- swap[swap$kind == 1, ]
    plane      <- swap[swap$kind == 2, ]
    shared     <- swap[swap$kind == 3, ]
    base_point <- swap[swap$kind == 4, ]
    scores     <- probe[probe$type == "label", ]
    expect_gt(nrow(reciprocal), 500L, label = label)
    # stage 1 of the reverse move is the inverse; so is the full ratio up to
    # the rounding of the Laplace laws refitted at the other end: the base
    # point's functions rebuilt from either end agree to rounding (below), and
    # the damped Newton fits take every step whose change of the objective is
    # at the rounding level, so both ends take the same steps. Measured maxima
    # over chain seeds 101-130 of each variant (40,500 records, Windows build):
    # 5.4e-12 relative for stage 1 and 9.4e-12 for the full ratio. A missing
    # Hastings term is at least 1e-3, so the bounds are 1e-9, which leave room
    # for the rounding of other platforms and compilers. From the forward
    # move's base quantities the reverse is the inverse to rounding
    expect_lt(max(abs(reciprocal$full) / pmax(1, abs(reciprocal$stage1))), 1e-9, label = label)
    expect_lt(max(abs(reciprocal$log_odds) / pmax(1, abs(reciprocal$selection_part))), 1e-9, label = label)
    expect_lt(max(abs(shared$stage1) / pmax(1, abs(reciprocal$stage1))), 1e-12, label = label)
    expect_lt(max(abs(shared$full) / pmax(1, abs(reciprocal$selection_part))), 1e-12, label = label)
    expect_lt(max(base_point$stage1, base_point$full), 1e-12, label = label)
    # the skip rule's M_0 (the largest sum of absolute child terms of a label
    # at the base point) is the same from both ends (measured at most 1e-15
    # relative), and no sweep of these fixtures skips L2
    expect_gt(min(shared$selection_part), 0, label = label)
    expect_lt(max(abs(shared$selection_part - shared$log_odds) / shared$selection_part), 1e-12, label = label)
    expect_identical(chain$counters[["swap_unrepresentable"]], 0, label = label)
    # nor any dry run at its base point (kind 8: M_0, Q, the working means'
    # representability criterion): Q stays far below its tolerance 1e-8 on
    # these fixtures (measured at most 4.6e-13 on the default fixture over
    # chain seeds 1-40)
    expect_false(any(swap$kind == 8), label = label)
    expect_gt(chain$counters[["max_swap_working"]], 0, label = label)
    expect_lt(chain$counters[["max_swap_working"]], 1e-10, label = label)
    # E_0, the shears, the evidences and the label proposals agree; the
    # accepted flows of the two directions are equal up to the same rounding
    # of the refitted laws (measured maximum 3.5e-11 over the same seeds; the
    # bound is 1e-9, as above)
    expect_lt(max(rebuilt$stage1), 1e-12, label = label)
    expect_lt(max(rebuilt$full), 1e-9, label = label)
    expect_lt(max(rebuilt$selection_part), 1e-12, label = label)
    expect_lt(max(abs(rebuilt$log_odds)), 1e-9, label = label)
    expect_lt(max(plane$stage1 / plane$full), 1e-10, label = label)
    expect_gt(nrow(scores), 500L, label = label)
    expect_lt(max(scores$stage1 / scores$full), 1e-12, label = label)
    expect_gt(chain$counters[["bias_label_swap_accepted"]], 100, label = label)
    # the live swap rebuilds its proposal end as an attempt that starts there
    # computes it, from that state alone: after every accepted swap an
    # attempt started at the new state, forced back, reproduces the rebuilt
    # geometry, base point (Q among its quantities), label proposals, stage 1
    # and proposal densities bitwise (kind 6); the skip rule's residuals stay
    # below its tolerance in every dry run (kind 5), and no live swap is
    # skipped
    reproduced <- swap[swap$kind == 6, ]
    rule       <- swap[swap$kind == 5, ]
    expect_gte(nrow(reproduced), chain$counters[["bias_label_swap_accepted"]], label = label)
    expect_true(all(reproduced$stage1 == 0 & reproduced$full == 0 & reproduced$selection_part == 0), label = label)
    expect_true(all(rule$stage1 == 1), label = label)
    expect_lt(max(rule$full, rule$selection_part), 1e-9, label = label)
    expect_identical(chain$counters[["swap_irreproducible"]], 0, label = label)
  }
  # Q is a function of the end's state (its rows' sizes before the move to
  # the base point), so the two builds of one end from the same working
  # state give it bitwise equal; mutant 31 takes the rebuilt end's sizes
  # from the attempt's own end, which both ends' Q still admit on this
  # fixture, and the reproduction records the difference
  chain <- .location_hand_chain(.bias_model(.bias_setup()), "bias_indicator", draws = 300L, burnin = 50L,
                                adapt = 50L, probe = TRUE, mutant = 31L)
  reproduced <- chain$probe[chain$probe$type == "label_swap" & chain$probe$kind == 6, ]
  expect_gt(nrow(reproduced), 10L)
  expect_true(all(is.finite(reproduced$stage1)))
  expect_gt(max(reproduced$stage1), 0)
})


test_that("the cross-family swap makes no proposal while its base point is not representable", {

  skip_if_not_installed("rjags")
  # L2 is the identity when the absolute child terms of one label at its base
  # point sum to more than 1e9 (M_0; swap_unrepresentable in
  # LocationBlockBiasMoves.h). The bias fixture with its first row at
  # 1000 / 1e-6 is fitted only through its cluster latent: from the chain's
  # start the latent is far from that row, so while adapting (the precision
  # guard is silent then) the base point's terms sum to about 5e17 and L2
  # makes no proposal; once the latent absorbs the row, M_0 is representable
  # (about 1.2e7) but the working means are not for the swap: Q, the size of
  # their rounding in units of the rows' SDs weighted by the standardized
  # residuals (swap_working_criterion), is at least 4e-6 against the
  # tolerance 1e-8, and L2 still makes no proposal. L2 is scheduled in every
  # sweep (switch swap = 1) and every attempt is traced.
  extreme <- function(y1, se1) {
    setup <- .bias_setup()
    setup$y[1]  <- y1
    setup$se[1] <- se1
    setup$residual_sd <- setup$se * sqrt(1 - .85)
    setup$loading     <- setup$se * sqrt(.85)
    setup
  }
  trace <- function() .Call("RoBMA_location_block_test_trace", PACKAGE = "RoBMA")
  invisible(trace())
  RoBMA:::.location_block_report_read()
  jags <- .location_compile(.bias_model(extreme(1000, 1e-6)),
                            config = c(refresh = -1L, fit_rule = -1L, timing = -1L, swap = 1L, trace = 1L, coin = -1L))
  rjags::adapt(jags, 100L, end.adaptation = FALSE)
  adapting <- RoBMA:::.location_block_report_read()
  attempts <- trace()
  expect_gt(adapting[["swap_unrepresentable"]], 0)
  # a skipped sweep is neither a proposal nor an attempt: the sweeps split
  # into skips and proposals, and only proposals are traced
  expect_identical(adapting[["factorization_skips"]], 0)
  expect_identical(adapting[["swap_unrepresentable"]] + adapting[["bias_label_swap_proposed"]], adapting[["sweeps"]])
  expect_lte(NROW(attempts), adapting[["bias_label_swap_proposed"]])
  expect_lte(adapting[["bias_label_swap_accepted"]], NROW(attempts))
  # the other moves run (L0 proposes once the other labels' masses are
  # positive, which the first sweeps' scores of about -5e17 leave at zero)
  expect_gt(adapting[["bias_label_proposed"]], 0)
  expect_gt(adapting[["bias_label_accepted"]], 0)
  expect_gt(adapting[["bias_gaussian_family_proposed"]], 0)
  expect_gt(adapting[["corrections_accepted"]], 0)
  # after adaptation (the counters and the evidence record read after the
  # adapting sweeps were cleared): every attempt is skipped at its base point
  # by Q, so nothing is proposed, traced or accepted, no admitted Q is kept,
  # and the evidence record counts the skips. These chains start where the
  # latent has absorbed the row (the state the default build reaches with
  # this seed after adaptation): from the far start above, a chain can stay
  # away from the row through adaptation (1 seed of 31 on a build whose
  # RoBMA fuses multiply-adds), and the precision guard then stops the
  # sampling, as it should; the expectations below are about the absorbed state
  absorbed <- list(Z = c(1420.46, 1.97, -949.72), slab = c(786.94, 0), g = c(2L, 1L))
  sampled_after <- function(mutant, setup = extreme(1000, 1e-6), draws = 200L, inits = list()) {
    invisible(trace())
    RoBMA:::.location_block_report_read()
    jags <- .location_compile(.bias_model(setup), mutant = mutant,
                              inits = c(list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = 4127L), inits),
                              config = c(refresh = -1L, fit_rule = -1L, timing = -1L, swap = 1L, trace = 1L, coin = -1L))
    rjags::adapt(jags, 100L, end.adaptation = TRUE)
    RoBMA:::.location_block_report_read()
    RoBMA:::.location_bias_evidence_read()
    invisible(trace())
    invisible(rjags::coda.samples(jags, "bias_indicator", n.iter = draws, progress.bar = "none"))
    list(counters = RoBMA:::.location_block_report_read(),
         record = RoBMA:::.location_bias_evidence_read()[["chains"]][1L, ], attempts = trace())
  }
  skipped <- sampled_after(0L, inits = absorbed)
  expect_identical(skipped$counters[["factorization_skips"]], 0)
  expect_identical(skipped$counters[["swap_unrepresentable"]], skipped$counters[["sweeps"]])
  expect_identical(skipped$counters[["bias_label_swap_proposed"]], 0)
  expect_identical(skipped$counters[["bias_label_swap_accepted"]], 0)
  expect_identical(skipped$counters[["max_swap_working"]], 0)
  expect_identical(NROW(skipped$attempts), 0L)
  expect_gt(unname(skipped$record[["swap_unrepresentable"]]), 0)
  expect_identical(unname(skipped$record[["swap_unrepresentable"]]), skipped$counters[["swap_unrepresentable"]])
  expect_identical(unname(skipped$record[["swap_proposed"]]), 0)
  # mutant 24 leaves Q out of the rule: swaps run and are accepted at this
  # row, with admitted Q above the tolerance (so the expectations above
  # fail). The proposal end's Laplace laws, refitted from working means whose
  # rounding the row's term amplifies (M_0 about 1.2e7), then reproduce the
  # forward end's to 1e-4, not to the tolerance 1e-9, in many attempts: those
  # swaps are the identity after stage 1 (swap_irreproducible, in the
  # counters and the evidence record; the trace holds them as passed
  # attempts without a stage-2 ratio), and every swap that went on to stage
  # 2 had both residuals within the tolerance. Whether an attempt reproduces
  # within the tolerance and is accepted follows the last bits of that
  # rounding: over 1000 draws at least 6 swaps were accepted for each of the
  # seeds 1-30 and 4127, with and without fused multiply-adds in RoBMA
  sampling <- sampled_after(24L, draws = 1000L, inits = absorbed)
  record   <- sampling$record
  expect_identical(sampling$counters[["swap_unrepresentable"]], 0)
  expect_identical(unname(record[["swap_unrepresentable"]]), 0)
  expect_identical(unname(record[["swap_proposed"]]), sampling$counters[["sweeps"]])
  expect_gt(sampling$counters[["bias_label_swap_accepted"]], 0)
  expect_gt(sampling$counters[["max_swap_working"]], 1e-8)
  expect_gt(sampling$counters[["swap_irreproducible"]], 0)
  expect_identical(unname(record[["swap_irreproducible"]]), sampling$counters[["swap_irreproducible"]])
  expect_lte(sampling$counters[["max_swap_stage1_residual"]], 1e-9)
  expect_lte(sampling$counters[["max_swap_full_residual"]], 1e-9)
  passed <- sampling$attempts[, "stage1_passed"] == 1
  expect_identical(as.numeric(sum(passed & is.nan(sampling$attempts[, "stage2"]))),
                   sampling$counters[["swap_irreproducible"]])
  expect_true(all(is.finite(sampling$attempts[sampling$attempts[, "accepted"] == 1, "stage2"])))
  # the rebuild's rule with Q in place: at 20 / 1e-4 (Q up to about 4e-9) a
  # few attempts are the identity after stage 1; mutant 19 scores the
  # reverse move with the forward end's base point, as the swap did before:
  # its residuals vanish, nothing is skipped after stage 1, and the
  # expectation of such skips fails
  moderate <- sampled_after(0L, extreme(20, 1e-4), draws = 600L)
  expect_identical(moderate$counters[["swap_unrepresentable"]], 0)
  expect_gt(moderate$counters[["swap_irreproducible"]], 0)
  expect_lte(moderate$counters[["max_swap_working"]], 1e-8)
  expect_identical(sampled_after(19L, extreme(20, 1e-4), draws = 600L)$counters[["swap_irreproducible"]], 0)
  # the skip after the target label's draw (an entered selection label's
  # terms join M_0 only there): mutant 25 takes every such base point as not
  # representable, so on the default fixture every attempt from a Gaussian
  # label (none, PET, PEESE) is skipped after its draw, counted as
  # swap_unrepresentable after adaptation in the counters and the evidence
  # record, and neither proposed nor traced, while swaps from the selection
  # label run; without the mutant no attempt is skipped
  after_draw <- sampled_after(25L, .bias_setup())
  counters   <- after_draw$counters
  expect_gt(counters[["swap_unrepresentable"]], 0)
  expect_identical(unname(after_draw$record[["swap_unrepresentable"]]), counters[["swap_unrepresentable"]])
  expect_identical(counters[["swap_unrepresentable"]] + counters[["bias_label_swap_proposed"]], counters[["sweeps"]])
  expect_identical(NROW(after_draw$attempts), as.integer(counters[["bias_label_swap_proposed"]]))
  expect_false(any(after_draw$attempts[, "from"] %in% c(1, 3, 4)))
  expect_gt(counters[["bias_label_swap_accepted"]], 0)
  ordinary <- sampled_after(0L, .bias_setup())
  expect_identical(ordinary$counters[["swap_unrepresentable"]], 0)
  expect_true(any(ordinary$attempts[, "from"] %in% c(1, 3, 4)))
  # data located far from zero relative to their SEs: every estimate and the
  # intercept's prior mean shifted by 30,000, the SEs times 100 (|y| / SE
  # from 1,200 to 1,670). The working means were read and carried where the
  # intercept holds the rows' size, which the move to the base point removes
  # from the means but not from their rounding; Q takes the rows' sizes at
  # the end's state (about 2.2e-8 here; at the base point it was 1.5e-12)
  # and skips every attempt after adaptation, counted in the counters and the
  # evidence record, while the other label moves run. Without the Q rule
  # (mutant 24) swaps are accepted
  relocated <- function(shift, scale) {
    setup <- .bias_setup()
    setup$y  <- setup$y + shift
    setup$se <- setup$se * scale
    setup$residual_sd <- setup$se * sqrt(1 - .85)
    setup$loading     <- setup$se * sqrt(.85)
    setup$prior_mean[1]   <- setup$prior_mean[1] + shift
    setup$gated_intercept <- FALSE
    setup
  }
  far <- sampled_after(0L, relocated(30000, 100))
  expect_identical(far$counters[["factorization_skips"]], 0)
  expect_identical(far$counters[["swap_unrepresentable"]], far$counters[["sweeps"]])
  expect_identical(far$counters[["bias_label_swap_proposed"]], 0)
  expect_identical(far$counters[["bias_label_swap_accepted"]], 0)
  expect_identical(unname(far$record[["swap_unrepresentable"]]), far$counters[["swap_unrepresentable"]])
  expect_gt(far$counters[["bias_label_accepted"]], 0)
  expect_gt(sampled_after(24L, relocated(30000, 100))$counters[["bias_label_swap_accepted"]], 0)
  # with the probe: every dry run is skipped at its base point (kind 8: M_0,
  # Q, after the draw), wherever M_0 is representable by Q, which stays at
  # least 100 times above the tolerance (measured at least 4.1e-6); without
  # the Q rule (mutant 24) the dry runs' rebuild rule skips most attempts,
  # the residuals being far above the tolerance, and an attempt started at
  # an accepted swap's state reproduces its rebuilt end bitwise here too
  chain <- .location_hand_chain(.bias_model(extreme(1000, 1e-6)), "bias_indicator", draws = 600L,
                                burnin = 100L, adapt = 100L, probe = TRUE, inits = absorbed)
  swap <- chain$probe[chain$probe$type == "label_swap", ]
  base_skips <- swap[swap$kind == 8, ]
  expect_gt(nrow(base_skips), 500L)
  expect_false(any(swap$kind %in% c(5, 6, 7)))
  expect_gt(min(base_skips$full[base_skips$stage1 <= 1e9]), 1e-6)
  expect_identical(chain$counters[["bias_label_swap_accepted"]], 0)
  # (4000 draws: the reproduction is checked after each accepted swap, at
  # least 4 for each of the seeds 1-30 and 4127 with and without fused
  # multiply-adds in RoBMA; with 600 draws some seeds accepted none)
  chain <- .location_hand_chain(.bias_model(extreme(1000, 1e-6)), "bias_indicator", draws = 4000L,
                                burnin = 100L, adapt = 100L, probe = TRUE, mutant = 24L, inits = absorbed)
  swap <- chain$probe[chain$probe$type == "label_swap", ]
  rule <- swap[swap$kind == 5 & swap$stage1 == 1, ]
  expect_gt(mean(rule$full > 1e-9 | rule$selection_part > 1e-9), 0.5)
  expect_gt(max(rule$full, rule$selection_part), 1e-6)
  reproduced <- swap[swap$kind == 6, ]
  expect_gt(nrow(reproduced), 0L)
  expect_true(all(reproduced$stage1 == 0 & reproduced$full == 0 & reproduced$selection_part == 0))
  expect_gt(chain$counters[["swap_irreproducible"]], 0)
  expect_lte(chain$counters[["max_swap_full_residual"]], 1e-9)

  # moderately extreme rows: nothing is unrepresentable, swaps run, and L2's
  # reciprocity keeps the bound of the ordinary fixtures (1e-9) at 2 / 1e-3
  # (M_0 about 31; measured 1.2e-10 for the full ratio, 4.5e-11 for stage 1).
  # At 20 / 1e-4 (M_0 about 4.6e3) the residual, from the Laplace laws
  # refitted at the other end, reaches 5e-9 (the full ratio's) and 2e-9
  # (stage 1's) with this seed: those attempts are skipped, and the swaps
  # that went on to stage 2 keep the bound. Neither row skips an attempt by
  # Q, the working means' representability criterion (a dry run skipped at
  # its base point with M_0 within its limit, kind 8): Q stays below its
  # tolerance 1e-8 (measured at most 4.7e-11 at 2 / 1e-3 and 4.1e-9 at
  # 20 / 1e-4 over chain seeds 1-40, L2 in every sweep).
  for (row in list(c(2, 1e-3), c(20, 1e-4))) {
    label <- paste(row, collapse = " / ")
    chain <- .location_hand_chain(.bias_model(extreme(row[1], row[2])), "bias_indicator", draws = 600L,
                                  burnin = 100L, adapt = 100L, probe = TRUE)
    swap <- chain$probe[chain$probe$type == "label_swap", ]
    reciprocal <- swap[swap$kind == 0, ]
    shared     <- swap[swap$kind == 3, ]
    expect_identical(chain$counters[["swap_unrepresentable"]], 0, label = label)
    expect_false(any(swap$kind == 8 & swap$stage1 <= 1e9), label = label)
    expect_lte(chain$counters[["max_swap_working"]], 1e-8, label = label)
    expect_gt(chain$counters[["bias_label_swap_accepted"]], 50, label = label)
    expect_gt(nrow(reciprocal), 500L, label = label)
    expect_lt(max(shared$selection_part, shared$log_odds), 1e4, label = label)
    expect_lt(max(abs(shared$selection_part - shared$log_odds) / shared$selection_part), 1e-12, label = label)
    expect_lte(max(chain$counters[c("max_swap_stage1_residual", "max_swap_full_residual")]), 1e-9, label = label)
    if (row[1] == 2) {
      expect_identical(chain$counters[["swap_irreproducible"]], 0, label = label)
      expect_lt(max(abs(reciprocal$full) / pmax(1, abs(reciprocal$stage1))), 1e-9, label = label)
      expect_lt(max(abs(reciprocal$log_odds) / pmax(1, abs(reciprocal$selection_part))), 1e-9, label = label)
    }
  }
})


test_that("a Q refusal at the swap's rebuilt end alone makes the attempt the identity", {

  skip_if_not_installed("rjags")
  # Q is checked at the attempt's own end (swap_start) and at its rebuilt
  # proposal end (swap_back, counted as swap_irreproducible). Mutant 29 takes
  # the rebuilt end's Q as 10 times its tolerance and changes nothing else;
  # on the default fixture the own ends pass (Q about 1e-13), M_0 is far
  # within its limit and both reciprocity residuals are about 1e-13, so only
  # that check refuses: every attempt that passes stage 1 is the identity
  # and none reaches stage 2. Mutant 30 is mutant 29 without the rebuilt
  # end's Q check: the same attempts reach stage 2 and are accepted.
  run <- function(mutant) {
    chain <- .location_hand_chain(.bias_model(.bias_setup()), "bias_indicator", draws = 300L, burnin = 50L,
                                  adapt = 50L, probe = TRUE, mutant = mutant,
                                  config = c(refresh = -1L, fit_rule = -1L, timing = -1L, swap = 1L, trace = -1L,
                                             coin = -1L))
    list(counters = chain$counters, swap = chain$probe[chain$probe$type == "label_swap", ])
  }
  forced <- run(29L)
  counters <- forced$counters
  expect_identical(counters[["swap_unrepresentable"]], 0)
  expect_identical(counters[["bias_label_swap_proposed"]], counters[["sweeps"]])
  expect_gt(counters[["bias_label_swap_stage1"]], 50)
  expect_identical(counters[["swap_irreproducible"]], counters[["bias_label_swap_stage1"]])
  expect_identical(counters[["bias_label_swap_accepted"]], 0)
  # the largest admitted Q is the own ends'
  expect_lt(counters[["max_swap_working"]], 1e-10)
  # no dry run is skipped at its own end (kind 8) or rebuilt (kind 5 with 1),
  # and no live attempt reaches stage 2 (kind 10)
  expect_false(any(forced$swap$kind %in% c(8, 10)))
  expect_gt(sum(forced$swap$kind == 5), 100L)
  expect_true(all(forced$swap$stage1[forced$swap$kind == 5] == 0))
  dropped <- run(30L)
  counters <- dropped$counters
  swap     <- dropped$swap
  expect_identical(counters[["swap_unrepresentable"]], 0)
  expect_identical(counters[["swap_irreproducible"]], 0)
  expect_gt(counters[["bias_label_swap_accepted"]], 50)
  expect_gt(sum(swap$kind == 10), 50L)
  # the rebuilt ends' Q is the forced one, M_0 of both ends and both
  # residuals pass
  expect_identical(counters[["max_swap_working"]], 10 * 1e-8)
  expect_lte(max(counters[c("max_swap_stage1_residual", "max_swap_full_residual")]), 1e-9)
  rule   <- swap[swap$kind == 5, ]
  shared <- swap[swap$kind == 3, ]
  expect_gt(nrow(rule), 100L)
  expect_true(all(rule$stage1 == 1))
  expect_lt(max(rule$full, rule$selection_part), 1e-9)
  expect_lt(max(shared$selection_part, shared$log_odds), 1e9)
})


test_that("the swap's live stage 2 takes the exact ratio of its two evaluated ends", {

  skip_if_not_installed("rjags")
  # Stage 2 of a live L2 attempt takes each proposal density from the base
  # point of the end that draws it: r = delta + log q'(x | x') - log q(x' | x),
  # q' at the rebuilt proposal end. With the probe, every live attempt that
  # the rule admits records the acceptance's input (r - stage 1), stage 1,
  # and from the attempt's quantities the exact ratio and the forward end's
  # own ratio delta + log q(x | x') - log q(x' | x) (kind 10). With a row at
  # 20 / 1e-4 the two ends' proposal densities differ by rounding that the
  # Laplace fits amplify, so the own-end ratio differs from the exact one on
  # admitted pairs; the input must be the exact ratio minus stage 1. Mutant
  # 26 takes the own-end ratio (both proposal densities from the forward
  # end's base point) and fails that.
  setup <- .bias_setup()
  setup$y[1]  <- 20
  setup$se[1] <- 1e-4
  setup$residual_sd <- setup$se * sqrt(1 - .85)
  setup$loading     <- setup$se * sqrt(.85)
  stage2 <- function(mutant) {
    chain <- .location_hand_chain(.bias_model(setup), "bias_indicator", draws = 300L, burnin = 100L, adapt = 100L,
                                  probe = TRUE, mutant = mutant,
                                  config = c(refresh = -1L, fit_rule = -1L, timing = -1L, swap = 1L, trace = -1L,
                                             coin = -1L))
    live <- chain$probe[chain$probe$type == "label_swap" & chain$probe$kind == 10, ]
    live$scale <- pmax(1, abs(live$selection_part), abs(live$full))
    live
  }
  # columns: stage1 the acceptance's input, full stage 1, selection_part the
  # exact ratio, log_odds the own-end ratio
  live <- stage2(0L)
  expect_gt(nrow(live), 100L)
  expect_true(all(is.finite(live$stage1)))
  expect_lte(max(abs(live$stage1 - (live$selection_part - live$full)) / live$scale), 4 * .Machine$double.eps)
  expect_gt(sum(abs(live$log_odds - live$selection_part) / live$scale > 4e3 * .Machine$double.eps), 10L)
  own <- stage2(26L)
  expect_gt(max(abs(own$stage1 - (own$selection_part - own$full)) / own$scale), 4e3 * .Machine$double.eps)
})


test_that("a discarded swap leaves the sweep's evidence flag as it was", {

  skip_if_not_installed("rjags")
  # The evidence flag about the surrogate (_not_exact) marks a sweep that
  # evaluated a surrogate that is not exact; L2's snapshot restores it, so
  # the swap's evaluations at its base point, at its rebuilt end and at a
  # proposal it discards do not mark the sweep. Mutant 27 runs only the
  # dormant draws D and L2 in every sweep and rejects every proposal of L2 at
  # stage 2. The chain starts and stays at the selection label, whose step
  # kernel the fixture's surrogate evaluates with a coarse normalizer, which
  # is not exact: every sweep has a selection endpoint and starts with the
  # flag false; D redraws the dormant PET and PEESE slopes from their priors
  # without evaluating a child, so the state changes and the exact
  # correction runs (an unchanged state is accepted before any evidence is
  # counted); L2 evaluates such surrogates at its base point (the selection
  # label's Laplace evidence), at the current state before the move and, for
  # attempts that pass stage 1, at the proposal's rebuilt end, then discards
  # the proposal. No sweep is evidence. Mutant 28 is mutant 27 with the
  # snapshot leaving the flag as the evaluations set it, as the build before
  # the snapshot held it did: every sweep is evidence.
  run <- function(mutant) {
    RoBMA:::.location_block_report_read()
    jags <- .location_compile(.bias_model(.bias_setup()), mutant = mutant,
                              inits = list(bias_indicator = 2L, .RNG.name = "base::Mersenne-Twister",
                                           .RNG.seed = 4127L),
                              config = c(refresh = -1L, fit_rule = -1L, timing = -1L, swap = 1L, trace = -1L,
                                         coin = -1L))
    rjags::adapt(jags, 50L, end.adaptation = TRUE)
    RoBMA:::.location_block_report_read()
    draws <- rjags::coda.samples(jags, "bias_indicator", n.iter = 100L, progress.bar = "none")
    list(draws = as.matrix(draws[[1L]]), counters = RoBMA:::.location_block_report_read())
  }
  kept <- run(27L)
  counters <- kept$counters
  expect_true(all(kept$draws[, "bias_indicator"] == 2))
  expect_identical(counters[["sweeps"]], 100)
  expect_identical(counters[["sweeps_selection"]], counters[["sweeps"]])
  expect_identical(counters[["bias_label_proposed"]], 0)
  expect_identical(counters[["bias_dormant_accepted"]], 2 * counters[["sweeps"]])
  expect_identical(counters[["corrections_accepted"]], counters[["sweeps"]])
  expect_identical(counters[["bias_label_swap_proposed"]], counters[["sweeps"]])
  expect_gt(counters[["bias_label_swap_stage1"]], 10)
  expect_identical(counters[["bias_label_swap_accepted"]], 0)
  expect_gt(counters[["surrogate_evaluations"]], 0)
  expect_identical(counters[["sweeps_evidence"]], 0)
  marked <- run(28L)$counters
  expect_identical(marked[["sweeps_selection"]], marked[["sweeps"]])
  expect_identical(marked[["sweeps_evidence"]], marked[["sweeps"]])
})


test_that("the swap's gate probabilities come from one centered pair of evidences", {

  # The gate of L2's target cell: log P(on) = -log(1 + e^-d) and log P(off) =
  # -log(1 + e^d) from d = log e_on - log e_off, for the draw and the score
  # alike. At evidences of size 1e9 (the reviewer's pair) the difference d =
  # 17 is exact, so the reference values below are exact up to the rounding
  # of log1p and exp; a total log(e_on + e_off) would round to the larger
  # evidence, and the draw from on - total would never take the smaller cell
  # while a score from off - total gave it 4.1e-8.
  gate <- function(off, on) .Call("RoBMA_location_block_test_numeric", "swap_gate", c(off, on), numeric(),
                                  PACKAGE = "RoBMA")
  on  <- -999000000
  off <- -999000017
  reference <- c(on = -log1p(exp(-17)), off = -17 - log1p(exp(-17)))
  value <- gate(off, on)
  expect_lt(max(abs(value - reference)), 1e-15 * 17)
  expect_gt(exp(value[2L]), 4.1e-8)
  expect_lt(abs(sum(exp(value)) - 1), 1e-15)
  # the reversed pair, a pair of size one, and cells without mass on a side
  expect_lt(max(abs(gate(on, off) - rev(reference))), 1e-15 * 17)
  expect_lt(max(abs(gate(log(.3), log(.7)) - log(c(.7, .3)))), 1e-15)
  expect_identical(gate(-Inf, -5), c(0, -Inf))
  expect_identical(gate(-5, -Inf), c(-Inf, 0))
  expect_identical(gate(-Inf, -Inf), c(-Inf, -Inf))
  # the normalization by a total, as the swap drew the gate before
  total <- max(on, off) + log(exp(on - max(on, off)) + exp(off - max(on, off)))
  expect_identical(exp(on - total), 1)
})


test_that("the bias label's evidence record shows a family that is never visited", {

  skip_if_not_installed("rjags")
  options <- RoBMA.options()[c("sampler.joint", "sampler.probabilities", "sampler.nested.weights")]
  on.exit(do.call(RoBMA.options, options), add = TRUE)
  RoBMA.options(sampler.joint = "auto", sampler.probabilities = character(), sampler.nested.weights = TRUE)
  dat <- data.frame(
    study = rep(1:8, each = 3L),
    yi    = 0.25 + rep(c(0.03, -0.02, 0.05, -0.04, 0.01, 0.02, -0.03, 0.04), 3L),
    vi    = rep(c(0.02, 0.03, 0.025), 8L)
  )
  # the Gaussian family (a PET label) has a negligible prior weight: never
  # visited; every swap out of the selection family is very unfavourable
  fit <- .location_one_chain(RoBMA.mv(
    yi = yi, V = vi, random = ~ 1 | study, data = dat, measure = "SMD",
    prior_bias = list(prior_weightfunction("one-sided", c(.05)), prior_PET("cauchy", list(0, 1), list(0, Inf),
                                                                             prior_weights = 1e-12)),
    prior_bias_null = FALSE, adapt = 100L, burnin = 50L, sample = 150L, chains = 1L, seed = 3L, autofit = FALSE,
    silent = TRUE))
  table <- fit[["fit"]][["samplers"]]
  expect_true("bias_indicator" %in% table[["Node"]][table[["Sampler"]] == "RoBMA::LocationBlock"])
  evidence <- attr(fit[["fit"]], "sampler_declarations", exact = TRUE)[["robma"]][["evidence"]]
  bias   <- evidence[["bias"]]
  labels <- evidence[["bias_labels"]]
  expect_identical(bias[["node"]], "bias_indicator")
  expect_identical(sum(labels[["sweeps"]]), 200)
  pet <- labels[["label"]] == 2L
  expect_identical(sum(labels[["sweeps"]][pet]), 0)
  expect_identical(bias[["leave_gaussian_attempts"]], 0)
  expect_gt(bias[["leave_selection_attempts"]], 0)
  expect_lt(bias[["leave_selection_log_ratio_max"]], -10)
  expect_identical(bias[["kept_leave_selection"]], 0)
  expect_identical(bias[["swap_proposed"]], bias[["leave_selection_attempts"]])
  expect_gt(bias[["label_proposed"]], 0)
  expect_lt(bias[["leave_selection_pass_probability_sum"]], 1e-4)
  # the label update gives the Gaussian family a negligible probability in
  # every sweep and never moves there; nothing leaves the Gaussian family
  expect_identical(bias[["leave_selection_label_updates"]], 200)
  expect_lt(bias[["leave_selection_label_probability_sum"]], 1e-6)
  expect_identical(bias[["leave_selection_label_crossings"]], 0)
  expect_lt(bias[["leave_selection_label_log_ratio_max"]], -10)
  expect_identical(bias[["leave_gaussian_label_updates"]], 0)
  expect_true(is.na(bias[["leave_gaussian_label_log_ratio_max"]]))
  # the retained draws: one censored run in the selection family
  sojourns <- evidence[["bias_sojourns"]]
  expect_identical(sojourns[["first_family"]], "selection")
  expect_identical(sojourns[["draws"]], 150L)
  expect_identical(sojourns[["initial_run"]], 150L)
  expect_true(sojourns[["initial_censored"]])
  expect_identical(c(sojourns[["leave_selection_crossings"]], sojourns[["leave_gaussian_crossings"]]), c(0L, 0L))
  expect_identical(c(sojourns[["longest_selection_run"]], sojourns[["longest_gaussian_run"]]), c(150L, 0L))
})


test_that("the bias label's family record agrees with the moves of the label update and the swap", {

  skip_if_not_installed("rjags")
  # A fixture whose label mixes between the families. The family changes that
  # L0 makes (and the stage-1 passes of L2) are Bernoulli draws with the
  # recorded probabilities given the past, so count minus probability sum is a
  # martingale with variance sum p (1 - p) <= sum p: |count - sum| < 4
  # sqrt(sum) fails with probability below 1e-4 for sums of this size
  RoBMA:::.location_bias_evidence_read()
  chain <- .location_hand_chain(.bias_model(.bias_setup()), "bias_indicator", draws = 4000L, burnin = 200L,
                                adapt = 200L)
  record <- as.list(RoBMA:::.location_bias_evidence_read()[["chains"]][1L, ])
  for (family in c("selection", "gaussian")) {
    column <- function(name) record[[paste0("leave_", family, "_label_", name)]]
    expect_gt(column("probability_sum"), 25, label = family)
    expect_lt(abs(column("crossings") - column("probability_sum")), 4 * sqrt(column("probability_sum")),
              label = paste(family, "label crossings"))
    expect_lt(column("log_ratio_max"), Inf, label = family)
  }
  expect_identical(record[["leave_selection_label_updates"]] + record[["leave_gaussian_label_updates"]],
                   record[["label_proposed"]])
  passes <- record[["leave_selection_pass_probability_sum"]] + record[["leave_gaussian_pass_probability_sum"]]
  expect_gt(passes, 25)
  expect_lt(abs(record[["swap_stage1"]] - passes), 4 * sqrt(passes))
  expect_lte(record[["leave_selection_pass_probability_sum"]], record[["leave_selection_attempts"]])
  expect_lte(record[["leave_gaussian_pass_probability_sum"]], record[["leave_gaussian_attempts"]])

  # L2's first stage is a plain Metropolis step, so its recorded pass
  # probabilities are min(1, exp(stage 1)) of the kernel that runs: the sums
  # equal those over the traced attempts after adaptation (the last ones)
  trace <- function() .Call("RoBMA_location_block_test_trace", PACKAGE = "RoBMA")
  invisible(trace())
  RoBMA:::.location_bias_evidence_read()
  .location_hand_chain(.bias_model(.bias_setup()), "bias_indicator", draws = 1600L, adapt = 100L, burnin = 0L,
                       config = c(refresh = -1L, fit_rule = -1L, timing = -1L, swap = -1L, trace = 1L))
  attempts <- trace()
  record <- as.list(RoBMA:::.location_bias_evidence_read()[["chains"]][1L, ])
  n <- record[["leave_selection_attempts"]] + record[["leave_gaussian_attempts"]]
  expect_gt(n, 100)
  rows <- attempts[utils::tail(seq_len(nrow(attempts)), n), , drop = FALSE]
  selection <- rows[, "from"] == 2
  expect_lt(abs(sum(exp(pmin(rows[selection, "stage1"], 0))) - record[["leave_selection_pass_probability_sum"]]),
            1e-9 * n)
  expect_lt(abs(sum(exp(pmin(rows[!selection, "stage1"], 0))) - record[["leave_gaussian_pass_probability_sum"]]),
            1e-9 * n)
})


test_that("the swap's first stage passes with the recorded Metropolis probabilities", {

  skip_if_not_installed("rjags")
  # L2 in every sweep and the coin of the two-state moves forced to Barker's
  # (which L2 does not use): the stage-1 passes after adaptation are
  # Bernoulli draws with the recorded probabilities min(1, exp(stage 1))
  # given the past, so |passes - their sum| < 4 sqrt(sum) as above. Mutant
  # 20 runs Barker's first stage, exp(s) / (1 + exp(s)), well below the
  # recorded probability for the favourable attempts that make up most of
  # this fixture's: it fails the bound.
  passes <- function(mutant) {
    RoBMA:::.location_bias_evidence_read()
    .location_hand_chain(.bias_model(.bias_setup()), "bias_indicator", draws = 2000L, burnin = 0L, adapt = 100L,
                         mutant = mutant,
                         config = c(refresh = -1L, fit_rule = -1L, timing = -1L, swap = 1L, trace = -1L, coin = 1L))
    record <- as.list(RoBMA:::.location_bias_evidence_read()[["chains"]][1L, ])
    c(attempts = record[["leave_selection_attempts"]] + record[["leave_gaussian_attempts"]],
      passed = record[["swap_stage1"]],
      expected = record[["leave_selection_pass_probability_sum"]] + record[["leave_gaussian_pass_probability_sum"]])
  }
  exact <- passes(0L)
  expect_gt(exact[["attempts"]], 1500)
  expect_lt(abs(exact[["passed"]] - exact[["expected"]]), 4 * sqrt(exact[["expected"]]))
  barker <- passes(20L)
  expect_gt(abs(barker[["passed"]] - barker[["expected"]]), 4 * sqrt(barker[["expected"]]))
})


test_that("the label update's record holds the move probability of the kernel that ran", {

  skip_if_not_installed("rjags")
  # At each L0 update the probe records the frozen state's label masses
  # exp(score - max), the kernel that ran and the probability of moving to
  # the other family that the evidence record adds. Recomputed here from the
  # masses: the Gibbs draw's sum of p_l over the other family, and Liu's
  # step's sum of p_l / max(1 - p_k, 1 - p_l). With the coin forced each way
  # the record follows the forced kernel, under the default coin the kernel
  # that ran; the coin's mixture (1/8 Gibbs, 7/8 Liu) differs from both.
  # Without adaptation every update enters the evidence record, whose sums
  # equal the probe's.
  updates <- function(coin) {
    RoBMA:::.location_bias_evidence_read()
    chain <- .location_hand_chain(.bias_model(.bias_setup()), "bias_indicator", draws = 300L, burnin = 0L,
                                  adapt = 0L, probe = TRUE,
                                  config = c(refresh = -1L, fit_rule = -1L, timing = -1L, swap = -1L, trace = -1L,
                                             coin = coin))
    record <- as.list(RoBMA:::.location_bias_evidence_read()[["chains"]][1L, ])
    rows  <- chain$probe[chain$probe$type == "label_record", ]
    heads <- which(rows$kind == 1)
    out <- do.call(rbind, lapply(heads, function(h) {
      k      <- rows$selection_part[h]
      labels <- rows[h + seq_len(rows$log_odds[h]), ]
      m      <- labels$full
      gaussian <- labels$selection_part == 1
      other  <- which(gaussian != gaussian[k])
      others <- function(j) sum(m[-j])
      data.frame(gibbs = rows$stage1[h] == 1, recorded = rows$full[h], from_gaussian = gaussian[k],
                 proposes = others(k) > 0, gibbs_cross = sum(m[other]) / sum(m),
                 liu_cross = sum(vapply(other, function(l) if (m[l] > 0) m[l] / max(others(k), others(l)) else 0, 0)))
    }))
    expect_identical(nrow(out), 300L)
    for (family in c("gaussian", "selection")) {
      own <- out$from_gaussian == (family == "gaussian")
      expect_identical(sum(own), as.integer(record[[paste0("leave_", family, "_label_updates")]]), label = family)
      expect_lt(abs(sum(out$recorded[own]) - record[[paste0("leave_", family, "_label_probability_sum")]]),
                1e-12 * max(1, sum(out$recorded[own])), label = family)
    }
    out
  }
  relative <- function(a, b) max(abs(a - b) / pmax(abs(b), 1e-300))
  gibbs <- updates(1L)
  expect_true(all(gibbs$gibbs[gibbs$proposes]))
  expect_lt(relative(gibbs$recorded[gibbs$proposes], gibbs$gibbs_cross[gibbs$proposes]), 1e-12)
  liu <- updates(2L)
  expect_false(any(liu$gibbs))
  expect_lt(relative(liu$recorded, liu$liu_cross), 1e-12)
  mixed <- updates(0L)
  ran <- mixed[mixed$proposes, ]
  expect_gt(sum(ran$gibbs), 10)
  expect_gt(sum(!ran$gibbs), 100)
  expect_lt(relative(ran$recorded, ifelse(ran$gibbs, ran$gibbs_cross, ran$liu_cross)), 1e-12)
  expect_gt(max(abs(ran$recorded - (ran$gibbs_cross / 8 + 7 * ran$liu_cross / 8))), 1e-3)
})


test_that("the bias label's sojourns match a direct count of the family runs", {

  # an independent count over the draws, one at a time
  reference <- function(selection) {
    n <- length(selection)
    runs <- list(selection = 0L, gaussian = 0L)
    crossings <- c(selection = 0L, gaussian = 0L)
    current <- 1L
    initial <- NA_integer_
    for (i in seq_len(n)) {
      if (i < n && selection[i + 1L] == selection[i]) {
        current <- current + 1L
        next
      }
      family <- if (selection[i]) "selection" else "gaussian"
      if (is.na(initial)) initial <- current
      runs[[family]] <- max(runs[[family]], current)
      if (i < n) crossings[[family]] <- crossings[[family]] + 1L
      current <- 1L
    }
    list(first_family = if (selection[1L]) "selection" else "gaussian", initial_run = initial,
         initial_censored = initial == n, leave_selection_crossings = crossings[["selection"]],
         leave_gaussian_crossings = crossings[["gaussian"]], longest_selection_run = runs[["selection"]],
         longest_gaussian_run = runs[["gaussian"]])
  }
  # labels 1 none, 2-3 step selection, 4 PET, 5 PEESE
  codes <- c(0, 1, 1, 2, 3)
  set.seed(7L)
  sequences <- list(c(4, 4, 1, 2, 3, 2, 5, 1, 1, 3, 3, 3, 4), rep(2L, 9L), rep(5L, 6L), 3L, c(1, 2),
                    sample(1:5, 500L, replace = TRUE, prob = c(.1, .4, .4, .05, .05)))
  for (labels in sequences) {
    summary <- RoBMA:::.sampler_bias_sojourns(codes[labels] == 1)
    expected <- reference(codes[labels] == 1)
    expect_identical(summary[["draws"]], length(labels))
    for (name in names(expected)) {
      expect_identical(summary[[name]], expected[[name]], label = paste(name, paste(head(labels, 6L), collapse = " ")))
    }
  }

  # the record maps labels to families by the bias declaration's codes, per chain
  fit <- list(data = list(RoBMA_location_bias_codes = codes),
              mcmc = list(cbind(mu = 0, bias_indicator = sequences[[1L]]), cbind(mu = 0, bias_indicator = 2)))
  declarations <- list(joint = list(markers = list(list(tag = "robma_location_bias"))))
  record <- RoBMA:::.sampler_bias_sojourn_record(fit, declarations, "bias_indicator")
  expect_identical(record[["chain"]], 1:2)
  expect_identical(record[["first_family"]], c("gaussian", "selection"))
  expect_identical(record[["initial_run"]], c(3L, 1L))
  expect_identical(record[["leave_gaussian_crossings"]], c(2L, 0L))
  expect_identical(record[["initial_censored"]], c(FALSE, TRUE))
  expect_null(RoBMA:::.sampler_bias_sojourn_record(fit, declarations, "other_indicator"))
  expect_null(RoBMA:::.sampler_bias_sojourn_record(fit, list(joint = list(markers = list())), "bias_indicator"))
})


test_that("fixed component slots and reversed composites keep a small discrete chain stationary", {

  # The sweep logic on a four-state chain (two labels x two values of a
  # coordinate): components reversible for a surrogate target, applied as a
  # palindromic sweep and followed by one exact correction against the target,
  # leave the target invariant only when the slot list does not depend on the
  # entry state and composite components reverse their internal order (sweep
  # conditions 1 and 3 of the location block).
  target    <- c(0.1, 0.2, 0.3, 0.4)     # states 1 = (S, 0), 2 = (S, 1), 3 = (G, 0), 4 = (G, 1)
  surrogate <- c(0.18, 0.12, 0.3, 0.4)   # equal to the target in the Gaussian family
  # a label move that does not commute with the value move
  label_move <- matrix(c(0, 0, .7, .3, 0, 0, .3, .7, .7, .3, 0, 0, .3, .7, 0, 0), 4, byrow = TRUE)
  value_move <- matrix(c(0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0), 4, byrow = TRUE)
  metropolis <- function(proposal) {
    K <- proposal * pmin(1, outer(1 / surrogate, surrogate))
    diag(K) <- 0
    diag(K) <- 1 - rowSums(K)
    K
  }
  L <- metropolis(label_move)
  V <- metropolis(value_move)
  # a component acting only in the Gaussian family (the joint refresh)
  J <- diag(4)
  J[3:4, 3:4] <- matrix(surrogate[3:4] / sum(surrogate[3:4]), 2, 2, byrow = TRUE)
  corrected <- function(Q) {
    w <- target / surrogate
    P <- Q * pmin(1, outer(1 / w, w))
    diag(P) <- 0
    diag(P) <- 1 - rowSums(P)
    P
  }
  stationary <- function(Q) max(abs(target %*% corrected(Q) - target))
  # fixed slots: the list (L, V, J) forward or reversed with probability 1/2
  fixed <- 0.5 * (L %*% V %*% J) + 0.5 * (J %*% V %*% L)
  expect_lt(stationary(fixed), 1e-15)
  # slots from the entry state: no J slot in sweeps entered in the selection family
  entry <- rbind((0.5 * (L %*% V) + 0.5 * (V %*% L))[1:2, ], fixed[3:4, ])
  expect_gt(stationary(entry), 1e-3)
  # a composite (V then L) that keeps its order in the reversed sweep
  forward <- 0.5 * (V %*% L %*% J) + 0.5 * (J %*% V %*% L)
  expect_gt(stationary(forward), 1e-3)
})


test_that("the label update sums the other labels' masses directly at both endpoints", {

  skip_if_not_installed("rjags")
  # Labels with a prior mass of 1e-30 against a dominant Gaussian label: the
  # sum of the others' masses subtracted from the total was zero, so L0 never
  # proposed although its record reported a positive crossing probability.
  # Summed directly, every update proposes (and essentially never moves);
  # labels without prior mass are never proposed.
  run <- function(weights) {
    setup <- .bias_setup()
    setup$bias_weights <- weights
    RoBMA:::.location_bias_evidence_read()
    chain <- .location_hand_chain(.bias_model(setup), "bias_indicator", draws = 100L, burnin = 0L, adapt = 1L,
                                  config = c(refresh = -1L, fit_rule = -1L, timing = -1L, swap = 0L, trace = -1L))
    list(chain = chain, evidence = RoBMA:::.location_bias_evidence_read()[["chains"]][1L, ])
  }
  tiny <- run(c(1, 1e-30, 1e-30, 1e-30))
  expect_identical(tiny$chain$decision, "claimed with the bias label")
  expect_identical(unname(tiny$evidence[["leave_gaussian_label_updates"]]), 100)
  expect_gt(tiny$evidence[["leave_gaussian_label_probability_sum"]], 0)
  expect_identical(unname(tiny$evidence[["label_proposed"]]), 100)
  expect_true(all(tiny$chain$draws[, "bias_indicator"] == 1))
  # Zero support of the other family: the factory refuses labels without
  # prior mass, so here the selection label has no likelihood (its step
  # weight of the bin observed for the first estimate is zero). From the
  # Gaussian labels L0 never proposes it and records no crossing probability.
  setup <- .bias_setup()
  setup$y[1L] <- -0.25
  setup$omega <- 0
  jags <- .location_compile(.bias_model(setup), inits = list(bias_indicator = 1L, .RNG.name = "base::Mersenne-Twister",
                                                             .RNG.seed = 4127L))
  rjags::adapt(jags, 1L, end.adaptation = TRUE)
  RoBMA:::.location_bias_evidence_read()
  draws <- as.matrix(rjags::coda.samples(jags, "bias_indicator", n.iter = 100L, progress.bar = "none")[[1L]])
  zero <- RoBMA:::.location_bias_evidence_read()[["chains"]][1L, ]
  expect_true(all(draws[, "bias_indicator"] != 2))
  expect_identical(unname(zero[["leave_gaussian_label_updates"]]), 100)
  expect_identical(unname(zero[["leave_gaussian_label_probability_sum"]]), 0)
  expect_identical(unname(zero[["leave_gaussian_label_crossings"]]), 0)
  expect_gt(zero[["label_proposed"]], 0)
})


test_that("the bias swap's and the label update's identities hold with the coin forced each way", {

  skip_if_not_installed("rjags")
  # coin 1: the Gibbs draw of L0 in every update; coin 2: never (L2's first
  # stage is Metropolis's under either). The reverse ratios and flows keep
  # their bounds (the test above), and L0's buffered scores agree with the
  # graph's.
  for (coin in 1:2) {
    config <- c(refresh = -1L, fit_rule = -1L, timing = -1L, swap = -1L, trace = -1L, coin = coin)
    chain <- .location_hand_chain(.bias_model(.bias_setup()), "bias_indicator", draws = 1200L, burnin = 100L,
                                  adapt = 100L, probe = TRUE, config = config)
    probe <- chain$probe
    swap <- probe[probe$type == "label_swap", ]
    reciprocal <- swap[swap$kind == 0, ]
    rebuilt <- swap[swap$kind == 1, ]
    scores <- probe[probe$type == "label", ]
    label <- paste("coin", coin)
    expect_gt(nrow(reciprocal), 500L, label = label)
    expect_lt(max(abs(reciprocal$full) / pmax(1, abs(reciprocal$stage1))), 1e-9, label = label)
    expect_lt(max(abs(reciprocal$log_odds) / pmax(1, abs(reciprocal$selection_part))), 1e-9, label = label)
    expect_lt(max(abs(rebuilt$log_odds)), 1e-9, label = label)
    expect_lt(max(scores$stage1 / scores$full), 1e-12, label = label)
    expect_gt(chain$counters[["bias_label_accepted"]], 0, label = label)
  }
})


test_that("a forced rejection after a label change restores every member and descendant", {

  skip_if_not_installed("rjags")
  # Mutant 15 makes the exact correction reject every sweep that changed the
  # bias label. With the probe, each such sweep reads the graph's endpoints
  # independently of the block's buffers, at the entry and after the
  # restoring write: the members (gates, slabs, latents, the indicator and
  # the slope sources), every parent value of every stochastic child of the
  # view (means, weights, kernel modes, the nested node's parameters) and the
  # densities (JAGS's own); every difference is zero.
  for (variant in list(list(), list(nested = TRUE))) {
    chain <- .location_hand_chain(.bias_model(do.call(.bias_setup, variant)), "bias_indicator", draws = 400L,
                                  burnin = 100L, adapt = 100L, mutant = 15L, probe = TRUE)
    label <- paste(c("variant", names(variant)), collapse = " ")
    records <- chain$probe[chain$probe$type == "restoration", ]
    expect_gt(nrow(records), 20L, label = label)
    expect_true(all(records$stage1 == 0 & records$full == 0 & records$selection_part == 0 &
                      records$log_odds == 0), label = label)
    expect_length(unique(chain$draws[, "bias_indicator"]), 1L)
    expect_lt(chain$counters[["corrections_accepted"]], chain$counters[["sweeps"]], label = label)
  }
  # the probe's comparison of two vectors of graph values: a value that is not
  # finite against a different value, on either side, and vectors of different
  # lengths count as infinitely different; equal values and two NaN as equal
  difference <- function(a, b) .Call("RoBMA_location_block_test_numeric", "endpoint_difference", as.double(a),
                                     as.double(b), PACKAGE = "RoBMA")
  expect_identical(difference(c(1, 2), c(1, 2)), 0)
  expect_identical(difference(c(1, 3), c(1, 2.5)), 0.5)
  expect_identical(difference(c(1, 2), c(1, NaN)), Inf)
  expect_identical(difference(c(NaN, 2), c(1, 2)), Inf)
  expect_identical(difference(c(1, NaN, 2), c(1, NaN, 2)), 0)
  expect_identical(difference(c(1, 2), c(1, Inf)), Inf)
  expect_identical(difference(c(-Inf, 2), c(Inf, 2)), Inf)
  expect_identical(difference(c(Inf, 2), c(Inf, 2)), 0)
  expect_identical(difference(c(1, 2), 1), Inf)
  # a NaN pair after a finite difference keeps it; a mismatch after it is infinite
  expect_identical(difference(c(1, NaN), c(1.5, NaN)), 0.5)
  expect_identical(difference(c(1, 2), c(1.5, NaN)), Inf)
})
