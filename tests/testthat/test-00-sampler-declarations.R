context("Sampler declarations passed to BayesTools")

# .fit() passes 'sampler_gates' and 'sampler_blocks' to BayesTools::JAGS_fit():
# gates for bias mixtures, and for selection models with estimated
# heterogeneity a function of the combined prior list that returns the block.
# These tests build objects without fitting and call the block function on the
# combined prior list JAGS_fit() would build (formula priors under their final
# names, then the fit priors).

.declaration_data <- function() {

  data <- data.frame(
    yi    = c(.10, .25, -.05, .32, .18, .40, .05, .22),
    vi    = c(.020, .030, .025, .040, .020, .035, .030, .028),
    study = c(1, 1, 2, 2, 3, 3, 4, 4),
    esid  = 1:8,
    x     = c(-1.2, -.4, .3, .9, -.7, .1, 1.4, .5),
    g     = factor(c("a", "b", "c", "a", "b", "c", "a", "b")),
    o     = ordered(c("lo", "mid", "hi", "lo", "mid", "hi", "lo", "mid"), levels = c("lo", "mid", "hi"))
  )
  V <- diag(data[["vi"]])
  V[1L, 2L] <- V[2L, 1L] <- .01

  return(list(data = data, V = V))
}

.declaration_prior_list <- function(object) {

  data         <- object[["data"]]
  priors       <- object[["priors"]]
  formula_args <- .create_jags_formula_args(data = data, priors = priors)
  formula_priors <- list()
  for (parameter in names(formula_args[["formula_list"]])) {
    output <- BayesTools::JAGS_formula(
      formula                = formula_args[["formula_list"]][[parameter]],
      parameter              = parameter,
      data                   = formula_args[["formula_data_list"]][[parameter]],
      prior_list             = formula_args[["formula_prior_list"]][[parameter]],
      formula_scale          = formula_args[["formula_scale_list"]][[parameter]],
      prior_random           = formula_args[["formula_random_prior_list"]][[parameter]],
      random_effects_compile = formula_args[["formula_random_effects_compile_list"]][[parameter]]
    )
    formula_priors <- c(formula_priors, output[["prior_list"]])
  }

  return(c(formula_priors, .create_fit_priors(data = data, priors = priors)))
}

# The declarations of an object: the gates flag, and the blocks as the block
# function returns them for the combined prior list ("none" when .fit() passes
# no function at all).
.declarations_of <- function(object) {

  fit_priors   <- .create_fit_priors(data = object[["data"]], priors = object[["priors"]])
  declarations <- .fit_sampler_declarations(
    data = object[["data"]], priors = object[["priors"]], fit_priors = fit_priors
  )
  blocks <- declarations[["sampler_blocks"]]
  if (!is.null(blocks)) {
    expect_true(is.function(blocks))
    blocks <- list(value = blocks(.declaration_prior_list(object)))
  } else {
    blocks <- "none"
  }

  return(list(gates = declarations[["sampler_gates"]], blocks = blocks))
}

.expect_declarations <- function(object, gates, blocks, label) {

  declarations <- .declarations_of(object)
  expect_identical(declarations[["gates"]], gates, info = label)
  if (identical(blocks, "none")) {
    expect_identical(declarations[["blocks"]], "none", info = label)
  } else {
    expect_identical(declarations[["blocks"]][["value"]], blocks, info = label)
  }
}

allocation_sd <- "mu__xRE_ALLOCx_heterogeneity__allocation_sd"
allocation_w  <- "mu__xRE_ALLOCx_heterogeneity__weight"
split_w       <- "mu__xRE_ALLOCx_component_1_split__weight"


test_that("the sampler options are logical, default on, and outside the selection prefix", {

  old_options <- RoBMA.options()
  on.exit(do.call(RoBMA.options, old_options), add = TRUE)

  for (name in c("sampler.gates", "sampler.blocks")) {
    expect_true(old_options[[name]], info = name)
    expect_false(startsWith(name, "selection."), info = name)
    for (value in list(NA, 1, "TRUE", c(TRUE, FALSE), NULL)) {
      expect_error(do.call(RoBMA.options, stats::setNames(list(value), name)),
                   "must be TRUE or FALSE", info = name)
    }
    updated <- do.call(RoBMA.options, stats::setNames(list(FALSE), name))
    expect_false(updated[[name]], info = name)
    expect_false(RoBMA.get_option(name), info = name)
  }
})


test_that("normal-likelihood, PET and PEESE models declare nothing", {

  setup <- .declaration_data()
  dat   <- setup[["data"]]
  V     <- setup[["V"]]
  glmm  <- list(ai = c(1, 2, 3), bi = c(9, 8, 7), ci = c(2, 1, 2), di = c(8, 9, 8),
                measure = "OR", only_priors = TRUE)

  objects <- list(
    brma      = brma(yi = yi, vi = vi, data = dat, measure = "SMD", only_priors = TRUE),
    bPET      = bPET(yi = yi, vi = vi, data = dat, measure = "SMD", only_priors = TRUE),
    bPEESE    = bPEESE(yi = yi, vi = vi, data = dat, measure = "SMD", only_priors = TRUE),
    BMA       = BMA(yi = yi, vi = vi, data = dat, measure = "SMD", only_priors = TRUE),
    brma.mv   = brma.mv(yi = yi, V = V, random = ~ 1 | study / esid, data = dat,
                        measure = "SMD", only_priors = TRUE),
    BMA.mv    = BMA.mv(yi = yi, V = V, random = ~ 1 | study / esid, data = dat,
                       measure = "SMD", only_priors = TRUE),
    bPET.mv   = bPET.mv(yi = yi, V = V, random = ~ 1 | study / esid, data = dat,
                        measure = "SMD", only_priors = TRUE),
    brma.glmm = do.call(brma.glmm, glmm),
    BMA.glmm  = do.call(BMA.glmm, glmm)
  )
  for (name in names(objects)) {
    .expect_declarations(objects[[name]], gates = FALSE, blocks = "none", label = name)
  }
})


test_that("single selection models block their global parameters", {

  setup <- .declaration_data()
  dat   <- setup[["data"]]
  V     <- setup[["V"]]

  # univariate: effect, heterogeneity and the weight function; the cluster
  # effects 'gamma' and the estimate effects 'theta' are unit-level latents
  .expect_declarations(
    bselmodel(yi = yi, vi = vi, data = dat, measure = "SMD", only_priors = TRUE),
    gates = FALSE, blocks = list(c("mu", "tau", "omega")), label = "bselmodel")
  .expect_declarations(
    bselmodel(yi = yi, vi = vi, cluster = study, data = dat, measure = "SMD", only_priors = TRUE),
    gates = FALSE, blocks = list(c("mu", "tau", "rho", "omega")), label = "bselmodel cluster")
  retained <- bselmodel(yi = yi, vi = vi, data = dat, measure = "SMD", only_priors = TRUE,
                        selection = selection_model(estimate_random_effects = "condition"))
  expect_true("theta" %in% names(.declaration_prior_list(retained)))
  .expect_declarations(retained, gates = FALSE, blocks = list(c("mu", "tau", "omega")),
                       label = "bselmodel estimate retained")

  # meta-regression: continuous and treatment-factor moderators join, a
  # mean-difference factor (a vector prior) does not
  .expect_declarations(
    bselmodel(yi = yi, vi = vi, mods = ~ x + g, data = dat, measure = "SMD", only_priors = TRUE),
    gates = FALSE, blocks = list(c("mu_intercept", "mu_x", "mu_g", "tau", "omega")),
    label = "bselmodel treatment")
  .expect_declarations(
    bselmodel(yi = yi, vi = vi, mods = ~ x + g, data = dat, measure = "SMD",
              set_contrast_factor_predictors = "meandif", only_priors = TRUE),
    gates = FALSE, blocks = list(c("mu_intercept", "mu_x", "tau", "omega")),
    label = "bselmodel meandif")

  # location-scale: the scale formula is the heterogeneity
  .expect_declarations(
    bselmodel(yi = yi, vi = vi, scale = ~ x, data = dat, measure = "SMD", only_priors = TRUE),
    gates = FALSE, blocks = list(c("log_tau_intercept", "log_tau_x", "mu", "omega")),
    label = "bselmodel location-scale")

  # multilevel allocation, with sampling variances and with a known V: the
  # allocation SD and the Dirichlet allocation weights join the block
  allocation_block <- list(c("mu_intercept", allocation_sd, allocation_w, "omega"))
  .expect_declarations(
    bselmodel.mv(yi = yi, vi = vi, random = ~ 1 | study / esid, data = dat,
                 measure = "SMD", only_priors = TRUE),
    gates = FALSE, blocks = allocation_block, label = "bselmodel.mv vi")
  .expect_declarations(
    bselmodel.mv(yi = yi, V = V, random = ~ 1 | study / esid, data = dat,
                 measure = "SMD", only_priors = TRUE),
    gates = FALSE, blocks = allocation_block, label = "bselmodel.mv V")
  .expect_declarations(
    bselmodel.mv(yi = yi, V = V, mods = ~ x + g, random = ~ 1 | study / esid, data = dat,
                 measure = "SMD", only_priors = TRUE),
    gates = FALSE, blocks = list(c("mu_intercept", "mu_x", "mu_g", allocation_sd, allocation_w, "omega")),
    label = "bselmodel.mv V regression")
  # a retained sampling source adds the unit-level 'sampling_z', never a member
  sampling <- bselmodel.mv(yi = yi, V = V, random = ~ 1 | study / esid, data = dat, measure = "SMD",
                           selection = selection_model(known_sampling_variance = "condition"),
                           only_priors = TRUE)
  expect_true("sampling_z" %in% names(.declaration_prior_list(sampling)))
  .expect_declarations(sampling, gates = FALSE, blocks = allocation_block,
                       label = "bselmodel.mv V sampling retained")

  # location-scale .mv: the scale formula replaces the allocation SD, so the
  # block holds the effect, the allocation weights, the scale coefficients and
  # the weight function
  .expect_declarations(
    bselmodel.mv(yi = yi, V = V, random = ~ 1 | study / esid, scale = ~ x, data = dat,
                 measure = "SMD", only_priors = TRUE),
    gates = FALSE, blocks = list(c("mu_intercept", allocation_w, "log_tau_intercept", "log_tau_x", "omega")),
    label = "bselmodel.mv V location-scale")

  # a direct random-effect SD (prior_random(sd = )) is not a member
  .expect_declarations(
    bselmodel.mv(yi = yi, V = V, random = ~ 1 | study, data = dat, measure = "SMD",
                 prior_heterogeneity = BayesTools::prior_random(
                   sd = BayesTools::prior("normal", list(0, 1), list(0, Inf))),
                 only_priors = TRUE),
    gates = FALSE, blocks = list(c("mu_intercept", "omega")), label = "bselmodel.mv direct SD")
})


test_that("model-averaged fits gate their bias branches and block only allocation parameters", {

  setup <- .declaration_data()
  dat   <- setup[["data"]]
  V     <- setup[["V"]]

  # a univariate ensemble has no allocation, so no block
  .expect_declarations(
    RoBMA(yi = yi, vi = vi, data = dat, measure = "SMD", only_priors = TRUE),
    gates = TRUE, blocks = NULL, label = "RoBMA")
  .expect_declarations(
    RoBMA(yi = yi, vi = vi, cluster = study, data = dat, measure = "SMD", only_priors = TRUE),
    gates = TRUE, blocks = NULL, label = "RoBMA cluster")

  # multilevel allocation: the allocation SD and the split weights; the
  # inclusion probability and the mixture of the effect are not members
  for (variances in c("vi", "V")) {
    object <- if (variances == "vi") {
      RoBMA.mv(yi = yi, V = vi, random = ~ 1 | study / esid, data = dat, measure = "SMD",
               only_priors = TRUE)
    } else {
      RoBMA.mv(yi = yi, V = V, random = ~ 1 | study / esid, data = dat, measure = "SMD",
               only_priors = TRUE)
    }
    .expect_declarations(object, gates = TRUE, blocks = list(c(allocation_sd, split_w)),
                         label = paste("RoBMA.mv", variances))
  }

  # without heterogeneity, and with the scale formula in place of the
  # allocation SD (the weights alone are one prior), no block
  .expect_declarations(
    RoBMA.mv(yi = yi, V = V, data = dat, measure = "SMD", only_priors = TRUE),
    gates = TRUE, blocks = NULL, label = "RoBMA.mv fixed")
  .expect_declarations(
    RoBMA.mv(yi = yi, V = V, random = ~ 1 | study / esid, scale = ~ x, data = dat,
             measure = "SMD", only_priors = TRUE),
    gates = TRUE, blocks = NULL, label = "RoBMA.mv location-scale")
})


test_that("guarded priors are never declared and inapplicable models declare no block", {

  setup <- .declaration_data()
  dat   <- setup[["data"]]
  V     <- setup[["V"]]
  half_normal <- BayesTools::prior("normal", list(0, 1), list(0, Inf))

  # an ordered moderator with a spike-and-slab total is not a member
  ordered <- bselmodel(yi = yi, vi = vi, mods = ~ o, data = dat, measure = "SMD", only_priors = TRUE,
    prior_mods = list(o = BayesTools::prior_ordered(
      BayesTools::prior_spike_and_slab(BayesTools::prior("gamma", list(2, 2))), allocation = c(.25, .75))))
  expect_true("mu_o" %in% names(.declaration_prior_list(ordered)))
  .expect_declarations(ordered, gates = FALSE, blocks = list(c("mu_intercept", "tau", "omega")),
                       label = "ordered spike-and-slab total")

  # a spike-and-slab random-effect SD makes the model an ensemble without
  # allocation parameters: no block
  .expect_declarations(
    bselmodel.mv(yi = yi, V = V, random = ~ 1 | study, data = dat, measure = "SMD",
                 prior_heterogeneity = BayesTools::prior_random(
                   sd = BayesTools::prior_spike_and_slab(half_normal)),
                 only_priors = TRUE),
    gates = FALSE, blocks = NULL, label = "prior_random spike-and-slab SD")

  # a spike-and-slab or mixture allocation SD is not a member, and the weights
  # alone are one prior
  discovery <- bselmodel.mv(yi = yi, V = V, random = ~ 1 | study / esid, data = dat,
                            measure = "SMD", only_priors = TRUE)
  blocks <- vapply(discovery[["formula_design"]][["mu"]][["random_effects"]], `[[`,
                   character(1), "block_name")
  allocation <- function(sd) BayesTools::prior_random(allocation = list(
    BayesTools::random_variance_allocation(
      name = "heterogeneity", terms = stats::setNames(blocks, blocks), sd = sd,
      weights = BayesTools::prior("dirichlet", list(alpha = c(1, 1))))))
  for (sd in list(BayesTools::prior_spike_and_slab(half_normal),
                  BayesTools::prior_mixture(list(half_normal, BayesTools::prior("normal", list(0, 2), list(0, Inf)))))) {
    object <- bselmodel.mv(yi = yi, V = V, random = ~ 1 | study / esid, data = dat, measure = "SMD",
                           prior_heterogeneity = allocation(sd), only_priors = TRUE)
    expect_true(allocation_sd %in% names(.declaration_prior_list(object)))
    .expect_declarations(object, gates = FALSE, blocks = NULL,
                         label = paste("allocation SD", class(sd)[2L]))
  }
  # the same allocation with a plain SD is blocked
  .expect_declarations(
    bselmodel.mv(yi = yi, V = V, random = ~ 1 | study / esid, data = dat, measure = "SMD",
                 prior_heterogeneity = allocation(half_normal), only_priors = TRUE),
    gates = FALSE, blocks = list(c("mu_intercept", allocation_sd, allocation_w, "omega")),
    label = "allocation plain SD")

  # point heterogeneity: univariate and fixed-effect .mv
  .expect_declarations(
    bselmodel(yi = yi, vi = vi, prior_heterogeneity = NULL, data = dat, measure = "SMD",
              only_priors = TRUE),
    gates = FALSE, blocks = NULL, label = "bselmodel fixed")
  .expect_declarations(
    bselmodel.mv(yi = yi, V = V, data = dat, measure = "SMD", only_priors = TRUE),
    gates = FALSE, blocks = NULL, label = "bselmodel.mv fixed")

  # all variation sources conditioned: the ordinary Gaussian likelihood
  conditioned <- bselmodel.mv(yi = yi, V = V, random = ~ 1 | study / esid, data = dat,
    measure = "SMD", only_priors = TRUE, selection = selection_model(
      estimate_random_effects = "condition", other_random_effects = "condition",
      known_sampling_variance = "condition"))
  expect_true(.selection_all_sources_conditioned(conditioned[["data"]]))
  .expect_declarations(conditioned, gates = FALSE, blocks = "none", label = "all sources conditioned")
})


test_that("unit-level latent parameters are never block members", {

  # independent factor priors under the names of RoBMA's unit-level latents are
  # eligible by type, and still excluded by name
  latent <- BayesTools::prior_factor_levels(
    BayesTools::prior_factor("normal", parameters = list(mean = 0, sd = 1), contrast = "independent"), 3L)
  prior_list <- list(
    mu  = BayesTools::prior("normal", list(0, 1)),
    tau = BayesTools::prior("normal", list(0, 1), list(0, Inf))
  )
  for (name in c("theta", "gamma", "sampling_z", "pi", "phi")) {
    prior_list[[name]] <- latent
  }
  prior_list[["eta"]] <- latent
  expect_identical(.sampler_block_members(prior_list), list(c("mu", "tau", "eta")))
  expect_null(.sampler_block_members(prior_list[c("tau", "theta", "gamma")]))
})


test_that("the options switch the declarations off", {

  old_options <- RoBMA.options()
  on.exit(do.call(RoBMA.options, old_options), add = TRUE)
  setup <- .declaration_data()
  dat   <- setup[["data"]]
  V     <- setup[["V"]]
  robma <- RoBMA.mv(yi = yi, V = V, random = ~ 1 | study / esid, data = dat, measure = "SMD",
                    only_priors = TRUE)
  selection <- bselmodel(yi = yi, vi = vi, data = dat, measure = "SMD", only_priors = TRUE)

  RoBMA.options(sampler.gates = FALSE)
  .expect_declarations(robma, gates = FALSE, blocks = list(c(allocation_sd, split_w)),
                       label = "gates off")
  RoBMA.options(sampler.gates = TRUE, sampler.blocks = FALSE)
  .expect_declarations(robma, gates = TRUE, blocks = "none", label = "blocks off")
  .expect_declarations(selection, gates = FALSE, blocks = "none", label = "blocks off")
})


test_that(".fit() passes the declarations to JAGS_fit() for fresh fits only", {

  old_options <- RoBMA.options()
  on.exit(do.call(RoBMA.options, old_options), add = TRUE)
  setup  <- .declaration_data()
  object <- bselmodel(yi = yi, vi = vi, data = setup[["data"]], measure = "SMD", only_priors = TRUE)
  forwarded <- NULL
  testthat::local_mocked_bindings(
    JAGS_fit = function(...) {
      forwarded <<- list(...)
      simpleError("Declaration forwarding probe")
    },
    JAGS_extend = function(...) {
      forwarded <<- list(...)
      simpleError("Declaration forwarding probe")
    },
    .package = "BayesTools")

  .fit(object)
  expect_false(forwarded[["sampler_gates"]])
  expect_true(is.function(forwarded[["sampler_blocks"]]))

  RoBMA.options(sampler.gates = FALSE, sampler.blocks = FALSE)
  .fit(object)
  expect_true("sampler_blocks" %in% names(forwarded))
  expect_false(forwarded[["sampler_gates"]])
  expect_null(forwarded[["sampler_blocks"]])

  # an extension replays the stored model and passes no declaration
  object[["fit"]] <- list("extension probe")
  .fit(object, extend = TRUE)
  expect_false(any(c("sampler_gates", "sampler_blocks") %in% names(forwarded)))
})
