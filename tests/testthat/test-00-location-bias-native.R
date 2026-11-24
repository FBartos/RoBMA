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
