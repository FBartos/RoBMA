context("Sampler declarations in fitted models")

# Fits with and without the sampler declarations RoBMA passes to
# BayesTools::JAGS_fit() (options 'sampler.blocks' and 'sampler.gates'). The
# declared nodes get the BayesTools samplers and every other node keeps the
# sampler it has without declarations. Certification case 'numerical-kernels'.

.declaration_fits_data <- function(studies = c(1, 2, 4, 12, 13)) {

  data("dat.assink2016", package = "metadat", envir = environment())
  dat <- dat.assink2016[dat.assink2016$study %in% studies, ]
  dat[["study"]] <- as.integer(factor(dat[["study"]]))
  # an ordered moderator of three levels (for the guard of an ordered total)
  dat[["o"]]     <- ordered(c("lo", "mid", "hi")[dat[["esid"]] %% 3L + 1L], levels = c("lo", "mid", "hi"))
  V <- metafor::vcalc(vi, cluster = study, type = deltype, obs = esid, rho = c(.7, .5), data = dat)

  return(list(data = dat, V = V))
}

# A small fit with the declarations on or off (both options), restoring them.
.declaration_fit <- function(constructor, args, declarations, sample = 100L, chains = 1L) {

  previous <- RoBMA::RoBMA.options()[c("sampler.blocks", "sampler.gates")]
  on.exit(do.call(RoBMA::RoBMA.options, previous), add = TRUE)
  RoBMA::RoBMA.options(sampler.blocks = declarations, sampler.gates = declarations)
  args <- utils::modifyList(list(
    measure = "SMD", chains = chains, adapt = 50L, burnin = 50L, sample = sample, seed = 1L,
    parallel = FALSE, silent = TRUE,
    convergence_checks = RoBMA::set_convergence_checks(max_Rhat = NULL, min_ESS = NULL)
  ), args)

  return(suppressWarnings(do.call(constructor, args)))
}

.declaration_markers <- function(fit) {

  lines <- strsplit(paste(as.character(fit[["fit"]][["model"]]), collapse = "\n"), "\n", fixed = TRUE)[[1L]]
  return(sum(grepl("BayesTools_(gate|block)_[0-9]+ ~ dbt_(gate|block)\\(", lines)))
}

# Declared nodes have the BayesTools samplers (one per declaration), and the
# fit without declarations has none of them and no marker; every node without a
# BayesTools sampler has the same sampler in both fits.
.expect_declared_assignment <- function(on, off, block_priors, n_gates, label) {

  declarations <- attr(on[["fit"]], "sampler_declarations", exact = TRUE)
  on_samplers  <- on[["fit"]][["samplers"]]
  off_samplers <- off[["fit"]][["samplers"]]
  bayestools   <- startsWith(on_samplers[["Sampler"]], "BayesTools::")

  expect_identical(unique(declarations[["blocks"]][["prior"]]), block_priors, info = label)
  expect_identical(nrow(declarations[["gates"]]), as.integer(n_gates), info = label)
  expect_identical(nrow(declarations[["declined"]]), 0L, info = label)
  expect_identical(sum(on_samplers[["Sampler"]] == "BayesTools::FactorSlice"),
                   nrow(declarations[["blocks"]]), info = label)
  expect_identical(sum(on_samplers[["Sampler"]] == "BayesTools::GatedSlice"),
                   nrow(declarations[["gates"]]), info = label)
  expect_identical(.declaration_markers(on),
                   nrow(declarations[["blocks"]]) + nrow(declarations[["gates"]]), info = label)

  expect_null(attr(off[["fit"]], "sampler_declarations", exact = TRUE))
  expect_false(any(startsWith(off_samplers[["Sampler"]], "BayesTools::")), info = label)
  expect_identical(.declaration_markers(off), 0L, info = label)
  expect_identical(nrow(on_samplers), nrow(off_samplers), info = label)
  kept <- on_samplers[!bayestools, c("Node", "Sampler")]
  expect_identical(
    off_samplers[["Sampler"]][match(kept[["Node"]], off_samplers[["Node"]])],
    kept[["Sampler"]], info = label
  )
}


test_that("declared nodes get the BayesTools samplers and the others keep theirs", {

  skip_on_cran()
  skip_if_not_certification("These fits (about 33 s) run in the numerical-kernels case.")
  skip_if_not_installed("rjags")
  skip_if_not_installed("metadat")
  skip_if_not_installed("metafor")

  # the JAGS modules of an ordinary RoBMA fit: no glm module
  expect_false("glm" %in% rjags::list.modules())
  setup <- .declaration_fits_data()
  dat   <- setup[["data"]]
  V     <- setup[["V"]]
  allocation <- c("mu__xRE_ALLOCx_heterogeneity__allocation_sd", "mu__xRE_ALLOCx_heterogeneity__weight")
  models <- list(
    bselmodel    = list(constructor = RoBMA::bselmodel, n_gates = 0L,
                        args = list(yi = quote(yi), vi = quote(vi), data = dat),
                        blocks = c("mu", "tau", "omega")),
    bselmodel.mv = list(constructor = RoBMA::bselmodel.mv, n_gates = 0L,
                        args = list(yi = quote(yi), V = V, random = ~ 1 | study / esid, data = dat),
                        blocks = c("mu_intercept", allocation, "omega")),
    RoBMA        = list(constructor = RoBMA::RoBMA, n_gates = 17L,
                        args = list(yi = quote(yi), vi = quote(vi), data = dat),
                        blocks = character()),
    RoBMA.mv     = list(constructor = RoBMA::RoBMA.mv, n_gates = 17L,
                        args = list(yi = quote(yi), V = V, random = ~ 1 | study / esid, data = dat),
                        blocks = c("mu__xRE_ALLOCx_heterogeneity__allocation_sd",
                                   "mu__xRE_ALLOCx_component_1_split__weight"))
  )
  for (name in names(models)) {
    model <- models[[name]]
    on    <- .declaration_fit(model[["constructor"]], model[["args"]], declarations = TRUE)
    off   <- .declaration_fit(model[["constructor"]], model[["args"]], declarations = FALSE)
    expect_null(on[["fit"]][["errors"]], info = name)
    .expect_declared_assignment(on, off, model[["blocks"]], model[["n_gates"]], label = name)
    # gates belong to the bias mixture only
    gates <- attr(on[["fit"]], "sampler_declarations", exact = TRUE)[["gates"]]
    expect_true(all(gates[["prior"]] == "bias"), info = name)
    expect_true(all(grepl("^equals\\(bias_indicator, [0-9]+\\)$", gates[["active"]])), info = name)
  }
})


test_that("guarded priors keep the stock samplers", {

  skip_on_cran()
  skip_if_not_certification("These fits (about 13 s) run in the numerical-kernels case.")
  skip_if_not_installed("rjags")
  skip_if_not_installed("metadat")
  skip_if_not_installed("metafor")

  setup <- .declaration_fits_data()
  dat   <- setup[["data"]]
  V     <- setup[["V"]]

  # an ordered moderator whose total is spike-and-slab: its nodes are not
  # declared; the other global parameters are blocked
  ordered <- list(yi = quote(yi), vi = quote(vi), mods = ~ o, data = dat,
    prior_mods = list(o = BayesTools::prior_ordered(
      BayesTools::prior_spike_and_slab(BayesTools::prior("gamma", list(2, 2))), allocation = c(.25, .75))))
  on  <- .declaration_fit(RoBMA::bselmodel, ordered, declarations = TRUE)
  off <- .declaration_fit(RoBMA::bselmodel, ordered, declarations = FALSE)
  .expect_declared_assignment(on, off, c("mu_intercept", "tau", "omega"), 0L, label = "ordered total")
  ordered_nodes <- grepl("^(mu_o|prior_par_eta_mu_o)", on[["fit"]][["samplers"]][["Node"]])
  expect_true(any(ordered_nodes))
  expect_false(any(startsWith(on[["fit"]][["samplers"]][["Sampler"]][ordered_nodes], "BayesTools::")))

  # a spike-and-slab random-effect SD: an ensemble without allocation
  # parameters, so nothing is declared and every sampler is unchanged
  spike <- list(yi = quote(yi), V = V, random = ~ 1 | study, data = dat,
    prior_heterogeneity = BayesTools::prior_random(sd = BayesTools::prior_spike_and_slab(
      BayesTools::prior("normal", list(0, 1), list(0, Inf)))))
  on  <- .declaration_fit(RoBMA::bselmodel.mv, spike, declarations = TRUE)
  off <- .declaration_fit(RoBMA::bselmodel.mv, spike, declarations = FALSE)
  .expect_declared_assignment(on, off, character(), 0L, label = "spike-and-slab SD")
  expect_identical(on[["fit"]][["samplers"]], off[["fit"]][["samplers"]])
})


test_that("declarations hold through parallel chains, autofit and update()", {

  skip_on_cran()
  skip_if_not_certification("These fits (about 33 s) run in the numerical-kernels case.")
  skip_if_not_installed("rjags")
  skip_if_not_installed("metadat")
  skip_if_not_installed("metafor")

  setup <- .declaration_fits_data()
  dat   <- setup[["data"]]
  V     <- setup[["V"]]
  args  <- list(yi = quote(yi), V = V, random = ~ 1 | study / esid, data = dat)
  n_markers <- 5L

  # parallel chains: workers load the BayesTools module for the markers
  parallel <- .declaration_fit(RoBMA::bselmodel.mv, c(args, list(parallel = TRUE)),
                               declarations = TRUE, chains = 2L)
  expect_null(parallel[["fit"]][["errors"]])
  expect_identical(.declaration_markers(parallel), n_markers)
  expect_identical(nrow(attr(parallel[["fit"]], "sampler_declarations", exact = TRUE)[["blocks"]]), n_markers)
  expect_identical(dim(parallel[["fit"]][["mcmc"]][[2L]])[1L], 100L)

  # automatic fitting extends the declared model
  autofit <- .declaration_fit(RoBMA::RoBMA.mv, c(args, list(autofit = TRUE,
    autofit_control = RoBMA::set_autofit_control(max_Rhat = 1.0001, min_ESS = 1e6,
                                                  sample_extend = 50L, max_extend = 1L))),
    declarations = TRUE)
  expect_null(autofit[["fit"]][["errors"]])
  expect_identical(nrow(autofit[["fit"]][["mcmc"]][[1L]]), 150L)
  expect_identical(.declaration_markers(autofit), 20L)
  expect_identical(sum(autofit[["fit"]][["samplers"]][["Sampler"]] == "BayesTools::FactorSlice"), 3L)

  # update() replays the stored model with its markers, serially and in parallel
  for (cores in c(1L, 2L)) {
    updated <- suppressWarnings(update(parallel, sample_extend = 50L, parallel = cores > 1L,
                                       cores = cores, silent = TRUE))
    expect_null(updated[["fit"]][["errors"]], info = cores)
    expect_identical(nrow(updated[["fit"]][["mcmc"]][[1L]]), 150L, info = cores)
    expect_identical(.declaration_markers(updated), n_markers, info = cores)
    expect_identical(attr(updated[["fit"]], "sampler_declarations", exact = TRUE),
                     attr(parallel[["fit"]], "sampler_declarations", exact = TRUE), info = cores)
  }
})


test_that("a serial update() of a selection model adapts its factor slice block", {

  skip_on_cran()
  skip_if_not_certification("These fits (about 7 s) run in the numerical-kernels case.")
  skip_if_not_installed("rjags")
  skip_if_not_installed("metadat")

  data("dat.assink2016", package = "metadat", envir = environment())
  caught  <- character()
  collect <- function(condition) {
    caught <<- c(caught, conditionMessage(condition))
    tryInvokeRestart(if (inherits(condition, "warning")) "muffleWarning" else "muffleMessage")
  }
  fit <- withCallingHandlers(RoBMA::bselmodel(
    yi = yi, vi = vi, data = dat.assink2016, measure = "SMD", chains = 2L, adapt = 500L,
    burnin = 500L, sample = 2000L, seed = 1L, parallel = FALSE, silent = TRUE,
    convergence_checks = RoBMA::set_convergence_checks(max_Rhat = NULL, min_ESS = NULL)
  ), warning = collect, message = collect)
  expect_identical(sum(fit[["fit"]][["samplers"]][["Sampler"]] == "BayesTools::FactorSlice"), 3L)
  updated <- withCallingHandlers(
    update(fit, sample_extend = 2000L, parallel = FALSE, silent = TRUE),
    warning = collect, message = collect)

  # no adaptation warning, in the conditions raised or stored with the fit
  expect_false(any(grepl("adapt", c(caught, updated[["fit"]][["warnings"]],
                                    attr(updated[["fit"]], "warnings", exact = TRUE)),
                         ignore.case = TRUE)))
  # the extension's effective sample size per draw of the block members is
  # comparable to the fit's; a block sampled unadapted after the recompile
  # (coordinate axes, unit widths) would fall far below it
  members <- c("mu", "tau", "omega[2]")
  draws   <- updated[["fit"]][["mcmc"]]
  expect_identical(nrow(draws[[1L]]), 4000L)
  ess <- function(rows) {
    coda::effectiveSize(coda::mcmc.list(lapply(draws, function(chain) {
      coda::mcmc(chain[rows, members, drop = FALSE])
    })))
  }
  ratio <- ess(2001:4000) / ess(1:2000)
  expect_true(all(ratio > .5), info = paste(names(ratio), round(ratio, 2), collapse = ", "))
})
