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
  records <- records[records$type != "ellipse", ]
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
  members <- c("mu_intercept_indicator", "mu_intercept_component_2", "mu_group_indicator",
               "mu_group_component_2", paste0("mu__xREx__study_xRE_Zx[", 1:8, ",1]"))

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
                                    paste0("mu__xREx__study_xRE_Zx[", 1:8, ",1]")))

  # an exact probability readout of a location gate keeps the enumerating sampler
  RoBMA.options(sampler.probabilities = "mu_intercept_indicator")
  readout <- fit_bma()
  expect_length(owned(readout), 0L)
  expect_identical(request(readout)[["status"]], "graph_declined")
  expect_identical(attr(readout[["fit"]], "sampler_declarations", exact = TRUE)[["robma"]][["probabilities"]][["status"]],
                   "active")
})
