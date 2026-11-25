context("Location block sampler of model-averaged .mv fits")
source(testthat::test_path("common-functions.R"))
skip_on_cran()

# RoBMA::LocationBlock (src/samplers/LocationBlock.cc) owns the location gates,
# their normal slabs and the standardized study effects. Every component move
# is reversible for the inner target and one exact correction per update
# compares it with the graph's own densities (helper-location-block.R).


test_that("Gaussian children: closed-form gate probabilities, exact stages, no selection evaluation", {

  skip_if_not_installed("rjags")
  setup <- .location_gaussian_setup()
  model <- .location_gaussian_model(setup)
  reference <- .location_gaussian_reference(setup)

  chain <- .location_hand_chain(model, c("g", "b", "b0"), draws = 40000L, probe = TRUE)
  expect_identical(chain$decision, "claimed")
  members <- c("b0", paste0("g[", 1:3, "]"), paste0("slab[", 1:3, "]"), paste0("Z[", 1:4, "]"))
  expect_true(all(chain$samplers[members] == "RoBMA::LocationBlock"))
  comparison <- .location_compare(.location_gaussian_observables(chain$draws), reference)
  expect_true(attr(comparison, "passes"), info = paste(round(comparison$z, 2), collapse = " "))

  # every stage-2 ratio and every correction is 1 up to rounding, and nothing
  # but Gaussian terms is evaluated
  counters <- chain$counters
  expect_equal(counters[["corrections_accepted"]], counters[["sweeps"]])
  expect_lt(counters[["max_stage2_gaussian"]], 1e-10)
  expect_lt(counters[["max_correction_gaussian"]], 1e-10)
  expect_identical(unname(counters[c("sweeps_selection", "surrogate_evaluations", "exact_evaluations")]),
                   c(0, 0, 0))
  expect_gt(counters[["flip_on_accepted"]], 1000)
  expect_identical(counters[["factorization_skips"]], 0)
  # frozen-state identity: the collapsed odds and line laws equal the direct
  # Gaussian evaluation of priors, children and proposals
  expect_lt(max(abs(chain$probe$full - chain$probe$stage1)), 1e-10)
  expect_setequal(unique(chain$probe$type), c("flip_on", "flip_off", "refresh", "study", "joint", "ellipse"))

  # power: a shear that differs between a flip and its reverse is not exact
  # and breaks the identity; so does a flip without its reverse proposal
  # density
  shear <- .location_hand_chain(model, c("g", "b", "b0"), draws = 40000L, mutant = 1L, probe = TRUE)
  expect_false(attr(.location_compare(.location_gaussian_observables(shear$draws), reference), "passes"))
  expect_gt(max(abs(shear$probe$full - shear$probe$stage1)), 1)
  dormant <- .location_hand_chain(model, c("g", "b", "b0"), draws = 2000L, mutant = 2L, probe = TRUE)
  expect_gt(max(abs(dormant$probe$full - dormant$probe$stage1)), 1)
})


test_that("the exact correction pairs each node's terms, so a large shared term keeps the others' differences", {

  skip_if_not_installed("rjags")
  # Three independent estimates, a gated N(0, 1) moderator and no
  # heterogeneity. The third row's log density (about -5e25) is shared by both
  # endpoints of every sweep, while the moderator moves its mean by 1e-30 per
  # unit: a correction that subtracts the endpoints' totals loses the first two
  # rows' differences (it gave a conditional mean of -0.11 against 0.53).
  options <- RoBMA.options()["sampler.joint"]
  on.exit(do.call(RoBMA.options, options), add = TRUE)
  RoBMA.options(sampler.joint = "location")
  data <- data.frame(y = c(.8, .8, 1e8), se = c(1, 1, 1e-5), x = c(1, 1, 1e-30))
  fit <- .location_one_chain(BMA(
    yi = y, sei = se, mods = ~ 0 + x, data = data, measure = "GEN", prior_mods = prior("normal", list(0, 1)),
    prior_heterogeneity = prior("point", list(0)), prior_heterogeneity_null = FALSE, prior_unit_information_sd = 1,
    standardize_continuous_predictors = FALSE, adapt = 100L, burnin = 100L, sample = 4000L, chains = 1L,
    seed = 516L, autofit = FALSE, silent = TRUE))
  request <- attr(fit[["fit"]], "sampler_declarations", exact = TRUE)[["robma"]][["joint"]]
  expect_identical(request[["status"]], "active")
  # the conditional law given the gate on is N(score / precision, 1 / precision)
  precision <- 1 + sum(data$x^2 / data$se^2)
  score     <- sum(data$x * data$y / data$se^2)
  p_on      <- stats::plogis(-log(precision) / 2 + score^2 / (2 * precision))
  mean_on   <- score / precision
  draws <- as.matrix(fit[["fit"]][["mcmc"]])
  on    <- draws[, "mu_x_indicator"] == 2
  comparison <- .location_compare(
    cbind(on = on * 1, b = draws[, "mu_x"], b2 = draws[, "mu_x"]^2),
    c(on = p_on, b = p_on * mean_on, b2 = p_on * (1 / precision + mean_on^2)))
  expect_true(attr(comparison, "passes"), info = paste(round(comparison$z, 2), collapse = " "))
  expect_lt(abs(mean(draws[on, "mu_x"]) - mean_on), 0.1)
  # every correction of the Gaussian label is one up to rounding
  expect_equal(request[["counters"]][["corrections_accepted"]], request[["counters"]][["sweeps"]])
  expect_lt(request[["counters"]][["max_correction_gaussian"]], 1e-9)
})


test_that("the factory verifies each design column against its baseline, and an untracked mean change is an error", {

  skip_if_not_installed("rjags")
  BayesTools::BayesTools_load_JAGS_module(quiet = TRUE)
  # mu = offset + beta with a gated N(0, 1) slab and offset ~ N(0, 1): from a
  # starting offset of 2^55 the probe's unit perturbation of the slab rounds
  # away and the column would pass as zero (the gate then froze)
  syntax <- paste(c(
    "model {", "  g ~ dcat(p[])", "  slab ~ dnorm(0, 1)", "  offset ~ dnorm(0, 1)",
    "  beta <- zero * (g == 1) + slab * (g == 2)", "  mu <- offset + beta", "  y ~ dnorm(mu, 1)",
    "  t1 <- robma_location_group(g, slab, 2, beta, 1)", "  t2 <- robma_location_setup(g, 1, 0)", "}"),
    collapse = "\n")
  reference <- stats::dnorm(.8, 0, sqrt(c(2, 3)))
  reference <- c(on = reference[[2L]] / sum(reference))
  run <- function(offset, gate) {
    connection <- textConnection(syntax)
    on.exit(close(connection))
    model <- rjags::jags.model(connection, data = list(p = c(.5, .5), zero = 0, y = .8), n.chains = 1L,
      n.adapt = 0L, quiet = TRUE, inits = list(g = gate, slab = 0, offset = offset,
                                               .RNG.name = "base::Mersenne-Twister", .RNG.seed = 376L))
    decision <- .location_hook()$decision
    rjags::adapt(model, 50L, end.adaptation = TRUE)
    draws <- as.matrix(rjags::coda.samples(model, "g", n.iter = 4000L, progress.bar = "none")[[1L]])
    list(decision = decision, owners = names(rjags::list.samplers(model)), on = (draws[, "g"] == 2) * 1)
  }
  for (gate in 1:2) {
    near <- run(0, gate)
    expect_identical(near$decision, "claimed without latents")
    far <- run(2^55, gate)
    expect_identical(far$decision,
                     "design probe: a coordinate's perturbation is not representable against the baseline means")
    expect_false("RoBMA::LocationBlock" %in% far$owners)
    # declined or correct, never frozen
    for (fit in list(near, far)) {
      comparison <- .location_compare(cbind(on = fit$on), reference)
      expect_true(attr(comparison, "passes"), info = paste(gate, round(comparison$z, 2)))
    }
  }

  # mutant 12 drops the design column of the last coordinate (here the slab)
  # after the probe: the graph then changes a mean that no component of the
  # sweep moved, an error naming the node that leaves the update's entry state
  # in place
  .location_hook(mutant = 12L)
  on.exit(.location_hook(mutant = 0L), add = TRUE)
  connection <- textConnection(syntax)
  jags <- rjags::jags.model(connection, data = list(p = c(.5, .5), zero = 0, y = .8), n.chains = 1L,
    n.adapt = 0L, quiet = TRUE, inits = list(g = 2L, slab = 0.3, offset = 0,
                                             .RNG.name = "base::Mersenne-Twister", .RNG.seed = 376L))
  close(connection)
  .location_hook(mutant = 0L)
  expect_true("RoBMA::LocationBlock" %in% names(rjags::list.samplers(jags)))
  rjags::adapt(jags, 0L, end.adaptation = TRUE)
  error <- NULL
  for (iteration in 1:50) {
    jags$sync()
    entry <- jags$state()[[1L]]
    error <- tryCatch({
      stats::update(jags, n.iter = 1L, progress.bar = "none")
      NULL
    }, error = conditionMessage)
    if (!is.null(error)) break
  }
  expect_match(error, "RoBMA_location_untracked_mean", fixed = TRUE)
  expect_match(error, "node y", fixed = TRUE)
  jags$sync()
  expect_identical(jags$state()[[1L]], entry)
})


test_that("each sweep runs the fit's fixed slot list, exactly reversed in a reversed sweep", {

  skip_if_not_installed("rjags")
  trace <- c(refresh = -1L, fit_rule = -1L, timing = -1L, swap = -1L, trace = 2L)
  run <- function(model, mutant = 0L, config = trace, sweeps = 300L) {
    .location_sequences()
    jags <- .location_compile(model, mutant = mutant, config = config)
    RoBMA:::.location_block_report_read()
    rjags::adapt(jags, 100L, end.adaptation = TRUE)
    counters <- RoBMA:::.location_block_report_read()
    stats::update(jags, n.iter = sweeps, progress.bar = "none")
    list(sequences = .location_sequences(), counters = counters)
  }
  # the slot lists of the forward sweeps, and those of the reversed sweeps
  # mirrored, by entry label
  slot_lists <- function(sequences) {
    unique(lapply(sequences, function(x) {
      slots <- if (x$reversed) x$slots[rev(seq_len(nrow(x$slots))), , drop = FALSE] else x$slots
      slots
    }))
  }

  # The bias fixture: its label changes inside sweeps, so sweeps are entered
  # in selection and in Gaussian labels (sweep conditions 1 and 3).
  model <- .bias_model(.bias_setup())
  exact <- run(model)$sequences
  expect_gt(length(exact), 300L)
  expect_true(any(vapply(exact, `[[`, logical(1), "selection")) && any(!vapply(exact, `[[`, logical(1), "selection")))
  expect_true(any(vapply(exact, `[[`, logical(1), "reversed")) && any(!vapply(exact, `[[`, logical(1), "reversed")))
  expect_length(slot_lists(exact), 1L)
  # mutant 11 takes the joint refresh's slot from the entry label
  expect_gt(length(slot_lists(run(model, mutant = 11L)$sequences)), 1L)

  # The two-coordinate slab with every group's refresh the slice update
  # (refresh variant 5): the composite's coordinates run in reverse order in
  # reversed sweeps; mutant 13 keeps their forward order.
  forced <- trace
  forced[["refresh"]] <- 5L
  model <- .location_factor_model(.location_selection_setup())
  coordinate_orders <- function(sequences) {
    out <- lapply(sequences, function(x) {
      own <- x$coordinates[x$coordinates[, 1L] == 1, 2L]
      if (length(own) == 0L) return(NULL)
      list(reversed = x$reversed, order = own)
    })
    Filter(Negate(is.null), out)
  }
  reversing <- function(orders) {
    all(vapply(orders, function(x) identical(x$order, if (x$reversed) c(1, 0) else c(0, 1)), logical(1)))
  }
  slice <- run(model, config = forced)
  orders <- coordinate_orders(slice$sequences)
  expect_true(any(vapply(orders, `[[`, logical(1), "reversed")) && any(!vapply(orders, `[[`, logical(1), "reversed")))
  expect_true(reversing(orders))
  expect_length(slot_lists(slice$sequences), 1L)
  # the record counts every group whose refresh became the slice update
  expect_identical(unname(slice$counters[c("slice_refresh_groups", "slice_refresh_multi")]), c(2, 1))
  expect_false(reversing(coordinate_orders(run(model, mutant = 13L, config = forced)$sequences)))
})


test_that("an error after the correction's graph write restores every member and descendant", {

  skip_if_not_installed("rjags")
  # Mutant 14 throws once after the exact correction writes the final state
  # (bias fixture: label nodes, slopes, gates, slabs and latents written).
  # The members are restored, and a model compiled from the restored state
  # (with the RNG state) takes the same steps, so no descendant is stale.
  model <- .bias_model(.bias_setup())
  state <- function(jags, internal = FALSE) {
    jags$sync()
    jags$state(internal = internal)[[1L]]
  }
  step <- function(jags) {
    tryCatch({
      stats::update(jags, n.iter = 1L, progress.bar = "none")
      NULL
    }, error = conditionMessage)
  }
  live <- .location_compile(model, mutant = 14L)
  rjags::adapt(live, 0L, end.adaptation = TRUE)
  error <- NULL
  for (iteration in 1:20) {
    entry <- state(live)
    error <- step(live)
    if (!is.null(error)) break
  }
  expect_match(error, "injected after the final graph write", fixed = TRUE)
  expect_identical(state(live), entry)
  replay <- .location_compile(model, inits = state(live, internal = TRUE))
  rjags::adapt(replay, 0L, end.adaptation = TRUE)
  for (iteration in 1:5) {
    expect_null(step(live))
    expect_null(step(replay))
    expect_identical(state(live, internal = TRUE), state(replay, internal = TRUE))
  }
})


test_that("the location block falls back to exact mode at the end of adaptation below its threshold", {

  skip_if_not_installed("rjags")
  # The selection fixture's crude surrogate: about 90% of the evidence sweeps
  # kept while adapting. Below a threshold of 0.99 the block falls back at the
  # first adaptOff() and evaluates no surrogate evidence afterwards; at 0 it
  # stays corrected.
  model <- .location_selection_model(.location_selection_setup())
  capture <- function(min_acceptance) {
    jags <- .location_compile(model, min_acceptance = min_acceptance)
    RoBMA:::.location_block_report_read()
    rjags::adapt(jags, 200L, end.adaptation = TRUE)
    warm <- RoBMA:::.location_block_report_read()
    stats::update(jags, n.iter = 100L, progress.bar = "none")
    list(warm = warm, sample = RoBMA:::.location_block_report_read())
  }
  strict <- capture(0.99)
  expect_gte(strict$warm[["sweeps_evidence"]], 50)
  expect_lt(strict$warm[["corrections_evidence_accepted"]], 0.99 * strict$warm[["sweeps_evidence"]])
  expect_identical(strict$warm[["fallback_exact"]], 1)
  expect_identical(strict$sample[["fallback_exact"]], 0)
  expect_identical(strict$sample[["sweeps_evidence"]], 0)
  expect_lt(strict$sample[["max_correction"]], 1e-10)
  lenient <- capture(0)
  expect_identical(lenient$warm[["fallback_exact"]], 0)
  expect_gt(lenient$sample[["sweeps_evidence"]], 0)
})


test_that("the closed-form tilt of a step-selection row has the derivatives of its density", {

  # The analytic Newton path of the location block (StepRowNormalizer.h): the
  # row's log density and its first two derivatives in the mean in closed form,
  # against Richardson-extrapolated central differences of the distribution's
  # own density (DSELNORMSTEPSWITCH::logDensity).
  probe <- function(y, mean, row, omega, lower, upper, h) {
    at <- mean + c(-h, h, -h / 2, h / 2, 0)
    out <- .Call("RoBMA_step_row_derivative_probe", y, mean, as.double(row), as.double(omega),
                 as.double(lower), as.double(upper), at, PACKAGE = "RoBMA")
    f <- out$exact
    first  <- function(k) (f[2 * k] - f[2 * k - 1]) / (2 * h / k)
    second <- function(k) (f[2 * k] - 2 * f[5] + f[2 * k - 1]) / (h / k)^2
    c(value = out$value, exact = f[5], first = out$first, first_fd = (4 * first(2) - first(1)) / 3,
      second = out$second, second_fd = (4 * second(2) - second(1)) / 3, ok = out$ok,
      # the block's direct row density is the distribution's, bit for bit
      direct = identical(out$direct, out$exact))
  }
  two_sided <- list(omega = c(1, .4), lower = c(1.96, -1e300), upper = c(1e300, 1.96))
  three_bins <- list(omega = c(1, .5, .1), lower = c(1.96, 0, -1e300), upper = c(1e300, 1.96, 0))
  one_sided <- list(omega = c(1, 0), lower = c(1.96, -1e300), upper = c(1e300, 1.96))
  rows <- list(
    # an ordinary row (sd, se, power, bin, sign, telescoping)
    ordinary  = c(list(y = .3, mean = .1, row = c(.25, .2, 1, 2, 1, 1)), two_sided),
    # the observation at a cut point
    cut_point = c(list(y = 1.96 * .2, mean = .2, row = c(.25, .2, 1, 1, 1, 1)), two_sided),
    # a negative effect direction, three bins and an observation weight
    negative  = c(list(y = -.5, mean = -.3, row = c(.3, .2, .5, 1, -1, 1)), three_bins),
    # far in the lower tail: P(observed) = Q(41.6) underflows outside the log scale
    lower_tail = c(list(y = .5, mean = -10, row = c(.25, .2, 1, 1, 1, 1)), one_sided),
    # far in the upper tail of a two-sided selection
    upper_tail = c(list(y = .5, mean = 9, row = c(.25, .2, 1, 1, 1, 1)), two_sided)
  )
  for (name in names(rows)) {
    r <- rows[[name]]
    result <- probe(r$y, r$mean, r$row, r$omega, r$lower, r$upper, h = .005)
    expect_true(result[["ok"]] == 1, info = name)
    expect_true(result[["direct"]] == 1, info = name)
    expect_true(all(is.finite(result)), info = name)
    expect_lte(abs(result[["value"]] - result[["exact"]]), 1e-10 * max(1, abs(result[["exact"]])), label = name)
    expect_lte(abs(result[["first"]] - result[["first_fd"]]), 1e-6 * max(1, abs(result[["first_fd"]])), label = name)
    expect_lte(abs(result[["second"]] - result[["second_fd"]]), 1e-6 * max(1, abs(result[["second_fd"]])),
               label = name)
  }
})


.location_test_data <- function() {

  # 8 studies of 3 estimates with independent sampling errors (several
  # likelihood blocks share each study effect), a three-level factor, and a
  # numeric moderator with overwhelming evidence; no heterogeneity
  data.frame(
    study = rep(1:8, each = 3L),
    group = factor(rep(c("a", "b", "c"), 8L)),
    x     = rep(c(-1, 0, 1), 8L) + rep(seq(-0.2, 0.2, length.out = 8L), each = 3L),
    yi    = 0.15 + 0.8 * (rep(c(-1, 0, 1), 8L) + rep(seq(-0.2, 0.2, length.out = 8L), each = 3L)) +
      rep(c(0.03, -0.02, 0.05, -0.04, 0.01, 0.02, -0.03, 0.04), 3L),
    vi    = rep(c(0.002, 0.003, 0.0025), 8L)
  )
}


test_that("frozen-state identities of the line moves hold in selection and Gaussian labels", {

  skip_if_not_installed("rjags")
  dat <- .location_test_data()
  on.exit(.location_hook(probe = FALSE), add = TRUE)
  RoBMA:::.location_block_report_read()

  .location_hook(probe = TRUE)
  selection <- .location_one_chain(RoBMA.mv(
    yi = yi, V = vi, mods = ~ group + x, random = ~ 1 | study, data = dat,
    measure = "SMD", prior_effect = prior("normal", list(mean = 0.2, sd = 0.5)),
    adapt = 100L, burnin = 50L, sample = 150L, chains = 1L, seed = 11L, autofit = FALSE, silent = TRUE))
  records <- .location_probe_records(.location_hook()$records)
  # the elliptical slice step (type "ellipse"): log T as the block computes it
  # equals the separately evaluated change of the inner target minus that of
  # the Gaussian working model
  ellipse <- records[records$type == "ellipse", ]
  expect_true(any(ellipse$selection) && any(!ellipse$selection))
  expect_lt(max(abs(ellipse$stage1 - ellipse$full) /
                  pmax(1, abs(ellipse$selection_part), abs(ellipse$log_odds))), 1e-9)
  expect_gt(max(abs(ellipse$stage1[ellipse$selection])), 0)
  expect_identical(max(abs(ellipse$stage1[!ellipse$selection])), 0)
  # the location moves (the bias-label records, LocationBlockReport.h, are
  # tested in test-00-location-bias-native.R)
  records <- records[records$type %in% c("flip_on", "flip_off", "refresh", "study", "swap", "joint"), ]
  request <- attr(selection[["fit"]], "sampler_declarations", exact = TRUE)[["robma"]][["joint"]]
  expect_identical(request[["status"]], "active")
  # both labels, overwhelming evidence and the heterogeneity gate off occur
  expect_true(any(records$selection) && any(!records$selection))
  expect_gt(max(abs(records$log_odds)), 1000)
  heterogeneity <- as.matrix(selection[["fit"]][["mcmc"]])
  expect_true(any(heterogeneity[, grep("include_component_1_indicator$", colnames(heterogeneity))] == 0))
  scale <- pmax(1, abs(records$log_odds), abs(records$stage1))
  expect_lt(max(abs(records$full - records$stage1) / scale), 1e-10)
  # the selection part is log A~(old) - log A~(new): zero in Gaussian labels
  expect_identical(max(abs(records$selection_part[!records$selection])), 0)
  expect_gt(max(abs(records$selection_part[records$selection])), 0)

  gaussian <- .location_one_chain(BMA.mv(
    yi = yi, V = vi, mods = ~ group + x, random = ~ 1 | study, data = dat,
    measure = "SMD", prior_effect = prior("normal", list(mean = 0.2, sd = 0.5)),
    adapt = 100L, burnin = 50L, sample = 150L, chains = 1L, seed = 11L, autofit = FALSE, silent = TRUE))
  records <- .location_probe_records(.location_hook(probe = FALSE)$records)
  expect_false(any(records$selection))
  ellipse <- records[records$type == "ellipse", ]
  expect_gt(nrow(ellipse), 0L)
  expect_lt(max(abs(ellipse$stage1 - ellipse$full) /
                  pmax(1, abs(ellipse$selection_part), abs(ellipse$log_odds))), 1e-9)
  records <- records[records$type != "ellipse", ]
  expect_lt(max(abs(records$full - records$stage1) / pmax(1, abs(records$log_odds), abs(records$stage1))), 1e-10)
  counters <- attr(gaussian[["fit"]], "sampler_declarations", exact = TRUE)[["robma"]][["joint"]][["counters"]]
  expect_lt(counters[["max_correction_gaussian"]], 1e-10)
  expect_identical(unname(counters[c("surrogate_evaluations", "exact_evaluations")]), c(0, 0))
})


test_that("the location block owns eligible fits, keeps its sampler on update() and in parallel, and declines otherwise", {

  skip_if_not_installed("rjags")
  options <- RoBMA.options()[c("sampler.joint", "selection.sampler", "sampler.probabilities")]
  on.exit(do.call(RoBMA.options, options), add = TRUE)
  dat <- .location_test_data()
  fit_bma <- function(..., mods = ~ group, chains = 1L, parallel = FALSE) {
    .location_one_chain(BMA.mv(yi = yi, V = vi, mods = mods, random = ~ 1 | study, data = dat, measure = "SMD",
           adapt = 60L, burnin = 50L, sample = 100L, chains = chains, parallel = parallel, seed = 5L,
           autofit = FALSE, silent = TRUE, ...))
  }
  owned <- function(fit) {
    table <- fit[["fit"]][["samplers"]]
    table[["Node"]][table[["Sampler"]] == "RoBMA::LocationBlock"]
  }
  request <- function(fit) attr(fit[["fit"]], "sampler_declarations", exact = TRUE)[["robma"]][["joint"]]
  # with the study effects, the heterogeneity scale's source (LocationScale.cc)
  members <- c("mu_intercept_indicator", "mu_intercept_component_2", "mu_group_indicator",
               "mu_group_component_2", paste0("mu__xREx__study_xRE_Zx[", 1:8, ",1]"),
               "mu__xRE_ALLOCx_heterogeneity__allocation_sd")

  eligible <- fit_bma()
  expect_setequal(owned(eligible), members)
  expect_identical(request(eligible)[["status"]], "active")
  expect_identical(request(eligible)[["sampler"]], "RoBMA::LocationBlock")
  counters <- request(eligible)[["counters"]]
  expect_equal(counters[["corrections_accepted"]], counters[["sweeps"]])

  extended <- .location_one_chain(update(eligible, sample_extend = 40L))
  expect_setequal(owned(extended), members)
  expect_identical(nrow(as.matrix(extended[["fit"]][["mcmc"]])), 140L)

  parallel <- fit_bma(chains = 2L, parallel = TRUE)
  workers <- parallel[["fit"]][["samplers_by_worker"]]
  if (is.null(workers)) workers <- list(list(samplers = parallel[["fit"]][["samplers"]]))
  for (worker in workers) {
    expect_setequal(worker[["samplers"]][["Node"]][worker[["samplers"]][["Sampler"]] == "RoBMA::LocationBlock"],
                    members)
  }
  expect_identical(request(parallel)[["status"]], "active")

  RoBMA.options(sampler.joint = "off")
  expect_length(owned(fit_bma()), 0L)
  RoBMA.options(sampler.joint = "auto", selection.sampler = "jags.default")
  expect_length(owned(fit_bma()), 0L)
  RoBMA.options(selection.sampler = options[["selection.sampler"]])

  # a non-normal slab is not a member; without a normal gated term the block
  # is declined
  heavy <- fit_bma(mods = ~ x, prior_effect = prior("t", list(0, 0.5, 3)), prior_mods = prior("t", list(0, 0.25, 3)))
  expect_length(owned(heavy), 0L)
  expect_identical(request(heavy)[["status"]], "graph_declined")
  partial <- fit_bma(mods = ~ x, prior_mods = prior("t", list(0, 0.25, 3)))
  expect_setequal(owned(partial), c("mu_intercept_indicator", "mu_intercept_component_2",
                                    paste0("mu__xREx__study_xRE_Zx[", 1:8, ",1]"),
                                    "mu__xRE_ALLOCx_heterogeneity__allocation_sd"))

  # an exact probability readout of a location gate keeps the enumerating sampler
  RoBMA.options(sampler.probabilities = "mu_intercept_indicator")
  readout <- fit_bma()
  expect_length(owned(readout), 0L)
  expect_identical(request(readout)[["status"]], "graph_declined")
  expect_identical(attr(readout[["fit"]], "sampler_declarations", exact = TRUE)[["robma"]][["probabilities"]][["status"]],
                   "active")
})


test_that("the evidence record tells a gate decided by the data from a mixing one", {

  skip_if_not_installed("rjags")
  dat <- .location_test_data()
  # a moderator without any association with the effects: its gate mixes
  dat$noise <- rep(c(-1, 1, 0.5, -0.5), 6L)
  fit <- .location_one_chain(BMA.mv(
    yi = yi, V = vi, mods = ~ x + noise, random = ~ 1 | study, data = dat, measure = "SMD",
    adapt = 100L, burnin = 100L, sample = 1000L, chains = 2L, seed = 3L, autofit = FALSE, silent = TRUE))
  record <- attr(fit[["fit"]], "sampler_declarations", exact = TRUE)[["robma"]][["evidence"]][["location"]]
  expect_setequal(unique(record$node), c("mu_intercept_indicator", "mu_x_indicator", "mu_noise_indicator"))
  expect_setequal(unique(record$chain), 1:2)
  # updates after adaptation: one inner flip proposal per gate and sweep,
  # every one reaching stage 1 in a Gaussian label
  expect_true(all(record$updates == 1100))
  expect_true(all(record$inner_flips_on_proposed + record$inner_flips_off_proposed == 1100))
  expect_true(all(record$inner_odds_evaluated == 1100))
  expect_true(all(record$inner_screened_out == 0))
  decided <- record[record$node == "mu_x_indicator", ]
  expect_true(all(decided$switched_on + decided$switched_off == 0))
  expect_true(all(decided$inner_flips_on_accepted + decided$inner_flips_off_accepted == 0))
  expect_true(all(decided$inner_log_odds_max < -50))
  mixing <- record[record$node == "mu_noise_indicator", ]
  expect_true(all(mixing$switched_on > 10 & mixing$switched_off > 10))
  expect_true(all(abs(mixing$switched_on - mixing$switched_off) <= 1))
  expect_true(all(mixing$inner_log_odds_max > -5))
  draws <- as.matrix(fit[["fit"]][["mcmc"]])
  expect_true(all(draws[, "mu_x_indicator"] == draws[1L, "mu_x_indicator"]))
})


test_that("the gate transitions of the evidence record are those the chain retained", {

  skip_if_not_installed("rjags")
  # The selection fixture: sweeps that the exact correction rejects (after
  # inner flips) and swaps that change two gates at once. The transitions of
  # each gate in the record are exactly those of its retained draws.
  BayesTools::BayesTools_load_JAGS_module(quiet = TRUE)
  old <- RoBMA::RoBMA.options()[c("selection.sampler", "selection.coarse_max_rules", "selection.coarse_min_acceptance")]
  on.exit(do.call(RoBMA::RoBMA.options, old), add = TRUE)
  RoBMA::RoBMA.options(selection.sampler = "slice.coarse", selection.coarse_max_rules = 3L,
                       selection.coarse_min_acceptance = 0)
  model <- .location_selection_model(.location_selection_setup())
  connection <- textConnection(model$syntax)
  jags <- rjags::jags.model(connection, data = model$data, n.chains = 1L, n.adapt = 0L, quiet = TRUE,
                            inits = list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = 4127L))
  close(connection)
  rjags::adapt(jags, 100L, end.adaptation = TRUE)
  RoBMA:::.location_block_report_read()
  RoBMA:::.location_block_evidence_read()
  jags$sync()
  entry <- jags$state()[[1L]]$g
  draws <- as.matrix(rjags::coda.samples(jags, "g", n.iter = 600L, progress.bar = "none")[[1L]])
  counters <- RoBMA:::.location_block_report_read()
  record <- RoBMA:::.location_block_evidence_read()
  path <- rbind(entry, draws)
  for (j in 1:3) {
    row <- record[record[, "slot"] == j, ]
    steps <- diff(path[, j])
    expect_identical(unname(row[["updates"]]), 600)
    expect_identical(unname(row[["switched_on"]]), as.numeric(sum(steps == 1)))
    expect_identical(unname(row[["switched_off"]]), as.numeric(sum(steps == -1)))
  }
  # both cases occur: rejected sweeps after accepted inner flips, and swaps
  expect_lt(counters[["corrections_accepted"]], counters[["sweeps"]])
  expect_gt(sum(record[, "inner_flips_on_accepted"]), sum(record[, "switched_on"]))
  expect_gt(counters[["swap_accepted"]], 0)
})


test_that("the evidence record names the screened indicators of a selection fit", {

  skip_if_not_installed("rjags")
  dat <- .location_test_data()
  # correlated sampling errors within study: cluster selection blocks, whose
  # surrogate gives the bias and heterogeneity indicators the screened sampler
  V <- diag(dat$vi)
  for (s in unique(dat$study)) {
    rows <- which(dat$study == s)
    V[rows, rows] <- 0.5 * sqrt(outer(dat$vi[rows], dat$vi[rows]))
    diag(V)[rows] <- dat$vi[rows]
  }
  fit <- .location_one_chain(RoBMA.mv(
    yi = dat$yi, V = V, mods = ~ x, random = ~ 1 | study, data = dat, measure = "SMD",
    adapt = 100L, burnin = 50L, sample = 100L, chains = 1L, seed = 4L, autofit = FALSE, silent = TRUE))
  table <- fit[["fit"]][["samplers"]]
  screened_nodes <- table[["Node"]][table[["Sampler"]] == "RoBMA::ScreenedFinite"]
  record <- attr(fit[["fit"]], "sampler_declarations", exact = TRUE)[["robma"]][["evidence"]][["screened"]]
  expect_gt(length(screened_nodes), 0L)
  expect_setequal(record$node, screened_nodes)
  expect_true(all(record$updates == 150))
  expect_true(all(record$leave_mean >= 0 & record$leave_mean <= record$leave_max & record$leave_max <= 1))
  # the names follow the samplers: the bias indicator has a level per bias branch
  if ("bias_indicator" %in% record$node) expect_gt(record$levels[record$node == "bias_indicator"], 2)
  heterogeneity <- grep("include_component_1_indicator$", record$node)
  if (length(heterogeneity)) expect_identical(record$levels[heterogeneity], 2)
})


test_that("the automatic policy gives independent fits with a gated moderator the location block", {

  skip_if_not_installed("rjags")
  dat <- .location_test_data()
  owner <- function(fit) {
    table <- fit[["fit"]][["samplers"]]
    unique(table[["Sampler"]][grepl("Joint|Location", table[["Sampler"]])])
  }
  record <- function(fit) attr(fit[["fit"]], "sampler_declarations", exact = TRUE)[["robma"]][["joint"]]
  independent <- .location_one_chain(RoBMA(
    yi = dat$yi, sei = sqrt(dat$vi), mods = ~ x, data = dat, measure = "SMD",
    adapt = 100L, burnin = 50L, sample = 100L, chains = 1L, seed = 6L, autofit = FALSE, silent = TRUE))
  expect_identical(owner(independent), "RoBMA::LocationBlock")
  expect_identical(record(independent)[["status"]], "active")
  # an unweighted Gaussian meta-regression whose unstandardized moderator the
  # automatic Gaussian candidate declines once its design is bound falls back
  # to the location block (the row budget, test-00-joint-gaussian-policy.R,
  # keeps the ordinary samplers)
  x <- seq(-1, 1, length.out = 24L)
  unstandardized <- .location_one_chain(BMA(
    yi = 0.05 + 0.20 * x + 0.03 * sin(seq_along(x)), sei = rep(0.2, length(x)), mods = ~ x, data = data.frame(x),
    measure = "SMD", standardize_continuous_predictors = FALSE, adapt = 100L, burnin = 50L, sample = 100L,
    chains = 1L, seed = 112L, autofit = FALSE, silent = TRUE))
  expect_identical(owner(unstandardized), "RoBMA::LocationBlock")
  expect_identical(record(unstandardized)[["status"]], "active")
  expect_identical(record(unstandardized)[["sampler"]], "RoBMA::LocationBlock")
  expect_setequal(record(unstandardized)[["nodes"]], c("mu_intercept_indicator", "mu_intercept_component_2",
                                                       "mu_x_indicator", "mu_x_component_2"))
  expect_gt(record(unstandardized)[["counters"]][["sweeps"]], 0)
  # without moderators the independent joint bias sampler keeps the fit
  unmoderated <- .location_one_chain(RoBMA(
    yi = dat$yi, sei = sqrt(dat$vi), data = dat, measure = "SMD",
    adapt = 100L, burnin = 50L, sample = 100L, chains = 1L, seed = 6L, autofit = FALSE, silent = TRUE))
  expect_identical(owner(unmoderated), "RoBMA::IndependentBiasJoint")
})


test_that("a gate whose collapsed odds are one does not alternate deterministically", {

  skip_if_not_installed("rjags")
  # No data on the gated coefficient and even prior odds: each flip's first
  # stage has ratio one. A Metropolis step then flips every sweep, so a chain
  # thinned by two keeps one state; Barker's acceptance after a coin with
  # probability 1/4 per attempt breaks the period. Under thin = 2 the gate-on
  # frequency of 4,000 draws is within 0.06 (about four of its standard
  # errors, lag-two autocorrelation 0.56) of 0.5 from either start.
  BayesTools::BayesTools_load_JAGS_module(quiet = TRUE)
  model <- list(syntax = paste(c(
    "model {", "  b0 ~ dnorm(0, 1)", "  g ~ dcat(p[])", "  slab ~ dnorm(0, 1)",
    "  beta <- zero * (g == 1) + slab * (g == 2)",
    "  for (i in 1:3) { mu[i] <- b0 + beta * x[i]", "    y[i] ~ dnorm(mu[i], 1) }",
    "  t0 <- robma_location_group(b0, b0, 0, b0, 1)", "  t1 <- robma_location_group(g, slab, 2, beta, 2)",
    "  t2 <- robma_location_setup(b0, 2, 0)", "}"), collapse = "\n"),
    data = list(p = c(.5, .5), zero = 0, x = c(0, 0, 0), y = c(.2, .5, -.1)))
  frequency <- function(gate, coin) {
    config <- c(refresh = -1L, fit_rule = -1L, timing = -1L, swap = -1L, trace = -1L, coin = coin)
    jags <- .location_compile(model, inits = list(g = gate, slab = 0.1, b0 = 0, .RNG.name = "base::Mersenne-Twister",
                                                  .RNG.seed = 91L), config = config)
    rjags::adapt(jags, 100L, end.adaptation = TRUE)
    draws <- as.matrix(rjags::coda.samples(jags, "g", n.iter = 8000L, thin = 2L, progress.bar = "none")[[1L]])
    mean(draws[, "g"] == 2)
  }
  for (gate in 1:2) {
    expect_lt(abs(frequency(gate, 0L) - 0.5), 0.06, label = paste("start", gate))
    # without the coin the thinned chain keeps its first state
    expect_true(frequency(gate, 2L) %in% c(0, 1), label = paste("start", gate, "without the coin"))
  }
})


test_that("the frozen-state identities hold with the coin forced each way", {

  skip_if_not_installed("rjags")
  # Barker's first stage (coin 1) and the Metropolis one (coin 2) leave every
  # dry-run identity of the Gaussian moves intact; the forced coin changes
  # which flips are accepted
  setup <- .location_gaussian_setup()
  model <- .location_gaussian_model(setup)
  accepted <- numeric()
  for (coin in 1:2) {
    config <- c(refresh = -1L, fit_rule = -1L, timing = -1L, swap = -1L, trace = -1L, coin = coin)
    chain <- .location_hand_chain(model, c("g", "b", "b0"), draws = 2000L, probe = TRUE, config = config)
    expect_lt(max(abs(chain$probe$full - chain$probe$stage1)), 1e-10, label = paste("coin", coin))
    expect_lt(chain$counters[["max_correction_gaussian"]], 1e-10, label = paste("coin", coin))
    accepted[coin] <- chain$counters[["flip_on_accepted"]] + chain$counters[["flip_off_accepted"]]
  }
  expect_lt(accepted[[1L]], accepted[[2L]])
})
