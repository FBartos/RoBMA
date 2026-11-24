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
    # point's functions rebuilt from either end agree to rounding (below), but
    # the damped Newton fit can take a different last step (up to 6e-8 of the
    # law's SD). Measured maximum 8.1e-8 relative over 13 chain seeds of each
    # variant; from the forward move's base quantities the reverse is the
    # inverse to rounding
    expect_lt(max(abs(reciprocal$full) / pmax(1, abs(reciprocal$stage1))), 1e-10, label = label)
    expect_lt(max(abs(reciprocal$log_odds) / pmax(1, abs(reciprocal$selection_part))), 1e-7, label = label)
    expect_lt(max(abs(shared$stage1) / pmax(1, abs(reciprocal$stage1))), 1e-12, label = label)
    expect_lt(max(abs(shared$full) / pmax(1, abs(reciprocal$selection_part))), 1e-12, label = label)
    expect_lt(max(base_point$stage1, base_point$full), 1e-12, label = label)
    # E_0, the shears, the evidences and the label proposals agree; the
    # accepted flows of the two directions are equal up to the same rounding
    # of the refitted laws (measured maximum 2.9e-8)
    expect_lt(max(rebuilt$stage1), 1e-12, label = label)
    expect_lt(max(rebuilt$full), 1e-9, label = label)
    expect_lt(max(rebuilt$selection_part), 1e-12, label = label)
    expect_lt(max(abs(rebuilt$log_odds)), 1e-7, label = label)
    expect_lt(max(plane$stage1 / plane$full), 1e-10, label = label)
    expect_gt(nrow(scores), 500L, label = label)
    expect_lt(max(scores$stage1 / scores$full), 1e-12, label = label)
    expect_gt(chain$counters[["bias_label_swap_accepted"]], 100, label = label)
  }
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
