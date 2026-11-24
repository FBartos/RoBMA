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
  # updates after adaptation: one flip proposal per gate and sweep, every one
  # reaching stage 1 in a Gaussian label
  expect_true(all(record$flips_on_proposed + record$flips_off_proposed == 1100))
  expect_true(all(record$odds_evaluated == 1100))
  expect_true(all(record$screened_out == 0))
  decided <- record[record$node == "mu_x_indicator", ]
  expect_true(all(decided$flips_on_accepted + decided$flips_off_accepted == 0))
  expect_true(all(decided$log_odds_max < -50))
  mixing <- record[record$node == "mu_noise_indicator", ]
  expect_true(all(mixing$flips_on_accepted > 10 & mixing$flips_off_accepted > 10))
  expect_true(all(mixing$log_odds_max > -5))
  draws <- as.matrix(fit[["fit"]][["mcmc"]])
  expect_true(all(draws[, "mu_x_indicator"] == draws[1L, "mu_x_indicator"]))
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
  # without moderators the independent joint bias sampler keeps the fit
  unmoderated <- .location_one_chain(RoBMA(
    yi = dat$yi, sei = sqrt(dat$vi), data = dat, measure = "SMD",
    adapt = 100L, burnin = 50L, sample = 100L, chains = 1L, seed = 6L, autofit = FALSE, silent = TRUE))
  expect_identical(owner(unmoderated), "RoBMA::IndependentBiasJoint")
})
