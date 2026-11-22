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
    prior_list <- .declaration_prior_list(object)
    blocks <- list(value = blocks(prior_list), prior_list = prior_list)
  } else {
    blocks <- "none"
  }

  return(list(gates = declarations[["sampler_gates"]], blocks = blocks))
}

# The declared blocks also pass BayesTools' block declaration check, so these
# tests fail against a BayesTools that refuses a member RoBMA declares.
.expect_declarations <- function(object, gates, blocks, label) {

  declarations <- .declarations_of(object)
  expect_identical(declarations[["gates"]], gates, info = label)
  if (identical(blocks, "none")) {
    expect_identical(declarations[["blocks"]], "none", info = label)
  } else {
    expect_identical(declarations[["blocks"]][["value"]], blocks, info = label)
    refusal <- tryCatch({
      BayesTools:::.bt_sampler_block_declarations(declarations[["blocks"]][["value"]],
                                                  declarations[["blocks"]][["prior_list"]])
      NULL
    }, error = conditionMessage)
    expect_null(refusal, info = label)
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

  # a direct random-effect SD (prior_random(sd = )) is a member
  .expect_declarations(
    bselmodel.mv(yi = yi, V = V, random = ~ 1 | study, data = dat, measure = "SMD",
                 prior_heterogeneity = BayesTools::prior_random(
                   sd = BayesTools::prior("normal", list(0, 1), list(0, Inf))),
                 only_priors = TRUE),
    gates = FALSE, blocks = list(c("mu_intercept", "mu__xREx__study_intercept", "omega")),
    label = "bselmodel.mv direct SD")
})


test_that("direct random-effect SDs and scalar structured correlations join the block of a single model only", {

  setup <- .declaration_data()
  dat   <- setup[["data"]]
  V     <- setup[["V"]]
  direct <- BayesTools::prior_random(sd = BayesTools::prior("normal", list(0, 1), list(0, Inf)))

  # single models: the direct SDs of nested intercepts, the SD and correlation
  # of cs(), and the correlation of hcs() next to its default allocation
  .expect_declarations(
    bselmodel.mv(yi = yi, V = V, random = ~ 1 | study / esid, data = dat, measure = "SMD",
                 prior_heterogeneity = direct, only_priors = TRUE),
    gates = FALSE,
    blocks = list(c("mu_intercept", "mu__xREx__esid_study_intercept", "mu__xREx__study_intercept", "omega")),
    label = "bselmodel.mv direct nested SDs")
  .expect_declarations(
    bselmodel.mv(yi = yi, V = V, random = ~ cs(g | study), data = dat, measure = "SMD", only_priors = TRUE),
    gates = FALSE, blocks = list(c("mu_intercept", "mu__xREx__study_sd", "mu__xREx__study_rho", "omega")),
    label = "bselmodel.mv cs")
  .expect_declarations(
    bselmodel.mv(yi = yi, V = V, random = ~ hcs(g | study), data = dat, measure = "SMD", only_priors = TRUE),
    gates = FALSE, blocks = list(c("mu_intercept", allocation_sd, allocation_w, "mu__xREx__study_rho", "omega")),
    label = "bselmodel.mv hcs")
  # ar1() and car() have a single SD and a correlation; har() has the
  # correlation next to its default allocation
  .expect_declarations(
    bselmodel.mv(yi = yi, V = V, random = ~ ar1(g | study), data = dat, measure = "SMD", only_priors = TRUE),
    gates = FALSE, blocks = list(c("mu_intercept", "mu__xREx__study_sd", "mu__xREx__study_rho", "omega")),
    label = "bselmodel.mv ar1")
  .expect_declarations(
    bselmodel.mv(yi = yi, V = V, random = ~ har(g | study), data = dat, measure = "SMD", only_priors = TRUE),
    gates = FALSE, blocks = list(c("mu_intercept", allocation_sd, allocation_w, "mu__xREx__study_rho", "omega")),
    label = "bselmodel.mv har")
  dat[["time"]] <- rep(1:2, 4)
  .expect_declarations(
    bselmodel.mv(yi = yi, V = V, random = ~ car(time | study), data = dat, measure = "SMD", only_priors = TRUE),
    gates = FALSE, blocks = list(c("mu_intercept", "mu__xREx__study_sd", "mu__xREx__study_rho", "omega")),
    label = "bselmodel.mv car")

  # an unstructured correlation is an LKJ primitive, not a prior-list entry:
  # the block is unchanged
  .expect_declarations(
    bselmodel.mv(yi = yi, V = V, random = ~ g | study, data = dat, measure = "SMD", only_priors = TRUE),
    gates = FALSE, blocks = list(c("mu_intercept", allocation_sd, allocation_w, "omega")),
    label = "bselmodel.mv LKJ")

  # ensembles with the same random structure keep their allocation members
  .expect_declarations(
    RoBMA.mv(yi = yi, V = V, random = ~ cs(g | study), data = dat, measure = "SMD", only_priors = TRUE),
    gates = TRUE, blocks = NULL, label = "RoBMA.mv cs")
  .expect_declarations(
    RoBMA.mv(yi = yi, V = V, random = ~ hcs(g | study), data = dat, measure = "SMD", only_priors = TRUE),
    gates = TRUE, blocks = list(c(allocation_sd, "mu__xRE_ALLOCx_study_components__weight")),
    label = "RoBMA.mv hcs")

  # a random-effect scalar that is a summary, part of an allocation or neither
  # an SD nor a correlation is not eligible
  half_normal <- BayesTools::prior("normal", list(0, 1), list(0, Inf))
  as_random   <- function(...) do.call(structure, c(list(half_normal), list(...)))
  declared <- function(prior) {
    "candidate" %in% .sampler_block_members(list(mu = half_normal, tau = half_normal, candidate = prior))[[1L]]
  }
  expect_true(declared(as_random(random_sd = TRUE)))
  expect_true(declared(as_random(random_correlation = TRUE)))
  expect_false(declared(as_random(random_sd = TRUE, random_summary = "sd")))
  expect_false(declared(as_random(random_sd = TRUE, random_allocation = "heterogeneity")))
  expect_false(declared(as_random(random_summary = "sd")))
})


test_that("regularized moderators join the block by coefficient, centered shrinkage factors do not", {

  dat    <- .declaration_data()[["data"]]
  global <- prior("cauchy", list(0, 0.2), list(0, Inf))

  # a shrinkage slope and a treatment factor join as their coefficients; the
  # group's global, slab and local scales are not prior-list entries
  .expect_declarations(
    bselmodel(yi = yi, vi = vi, mods = ~ x + g, data = dat, measure = "SMD",
              prior_mods = list(pool = prior_regularized("horseshoe", scale = global, slab_scale = 1)),
              only_priors = TRUE),
    gates = FALSE, blocks = list(c("mu_intercept", "mu_x", "mu_g", "tau", "omega")),
    label = "bselmodel regularized treatment")
  # a centered factor shrunk as one exchangeable unit is not a member
  .expect_declarations(
    bselmodel(yi = yi, vi = vi, mods = ~ x + g, data = dat, measure = "SMD",
              prior_mods = list(pool = prior_regularized("lasso", scale = global)),
              set_contrast_factor_predictors = "meandif", only_priors = TRUE),
    gates = FALSE, blocks = list(c("mu_intercept", "mu_x", "tau", "omega")),
    label = "bselmodel regularized meandif")
})

test_that("block sources decline bound R-squared slopes and retain fixed weight priors", {

  compiled <- BayesTools::JAGS_formula(
    ~ x, parameter = "mu", data = data.frame(x = c(-1, 0, 2)),
    prior_list = list(
      intercept = BayesTools::prior("normal", list(0, 1)),
      signal = BayesTools::prior_R2(terms = "x", model = BayesTools::R2_model(sigma = 1))
    )
  )
  expect_true(BayesTools::is.prior.R2(compiled[["prior_list"]][["mu_x"]]))
  priors <- c(compiled[["prior_list"]], list(
    tau = BayesTools::prior("normal", list(0, 1), list(0, Inf)),
    omega = BayesTools::prior_weightfunction("one-sided", steps = .025,
                                             weights = BayesTools::wf_fixed(c(1, .5)))
  ))
  members <- .sampler_block_members(priors)[[1L]]
  expect_false("mu_x" %in% members)
  expect_true(all(c("mu_intercept", "tau", "omega") %in% members))
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
  # inclusion probability and the mixture of the effect are not members. The
  # block is left out only where the corrected selection sampler (the default)
  # can take the allocation parameters: with a dependent selection block and a
  # product-rule branch. On sampling variances alone every selection block is a
  # single estimate, and best-weight branches on rank-one or factor blocks give
  # no coarse rule, so the block stays there under both settings.
  old_options <- RoBMA.options()
  on.exit(do.call(RoBMA.options, old_options), add = TRUE)
  best   <- selection_model(weight_rule = "best", group = "study")
  # a constant within-study correlation: rank-one blocks only (no dense block)
  V_rank_one <- metafor::vcalc(vi, cluster = study, obs = esid, rho = .6, data = dat)
  models <- list(
    vi   = RoBMA.mv(yi = yi, V = vi, random = ~ 1 | study / esid, data = dat, measure = "SMD",
                    only_priors = TRUE),
    V    = RoBMA.mv(yi = yi, V = V, random = ~ 1 | study / esid, data = dat, measure = "SMD",
                    only_priors = TRUE),
    best = RoBMA.mv(yi = yi, V = V_rank_one, random = ~ 1 | study / esid, data = dat, measure = "SMD",
                    prior_bias = list(
                      BayesTools::prior_weightfunction("one-sided", steps = .025, model = best),
                      BayesTools::prior_weightfunction("one-sided", steps = c(.025, .05), model = best)),
                    only_priors = TRUE)
  )
  expect_true(all(.data_selection_execution_plan(models[["vi"]][["data"]])[["block_methods"]] == "singleton"))
  expect_true(any(.data_selection_execution_plan(models[["V"]][["data"]])[["block_methods"]] != "singleton"))
  best_methods <- .data_selection_execution_plan(models[["best"]][["data"]])[["block_methods"]]
  expect_true(any(best_methods %in% c("rank_one", "factor")))
  expect_false(any(best_methods == "dense"))
  for (sampler in c("slice", "slice.coarse")) {
    RoBMA.options(selection.sampler = sampler)
    for (variances in names(models)) {
      # the allocation block is declared under both settings; the corrected
      # block sampler takes it where a selection block takes part
      .expect_declarations(models[[variances]], gates = TRUE, blocks = list(c(allocation_sd, split_w)),
                           label = paste("RoBMA.mv", variances, sampler))
    }
  }
  # single selection models keep their block under both settings
  selection <- bselmodel.mv(yi = yi, V = V, random = ~ 1 | study / esid, data = dat, measure = "SMD",
                            only_priors = TRUE)
  for (sampler in c("slice", "slice.coarse")) {
    RoBMA.options(selection.sampler = sampler)
    .expect_declarations(selection, gates = FALSE,
                         blocks = list(c("mu_intercept", allocation_sd, allocation_w, "omega")),
                         label = paste("bselmodel.mv", sampler))
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


test_that("an ensemble's block is its allocation SD and weights, through the block function too", {

  setup  <- .declaration_data()
  object <- RoBMA.mv(yi = yi, V = setup[["V"]], random = ~ 1 | study / esid, data = setup[["data"]],
                     measure = "SMD", only_priors = TRUE)
  prior_list <- .declaration_prior_list(object)
  expect_identical(.sampler_block_members(prior_list), list(c(allocation_sd, split_w)))
  expect_identical(.sampler_block_function(character())(prior_list), list(c(allocation_sd, split_w)))
})


test_that("the native factor rule budget is the corrected sampler's", {

  # the sampler's rule: order^rank nodes within the budget of 2^18
  expect_identical(.Call("RoBMA_selnorm_factor_rule_affordable", c(7, 9, 7, 63, 95), c(6, 6, 7, 3, 3),
                         PACKAGE = "RoBMA"), c(TRUE, FALSE, FALSE, TRUE, FALSE))
  expect_error(.Call("RoBMA_selnorm_factor_rule_affordable", 7.5, 6, PACKAGE = "RoBMA"))
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

  # the ensemble's allocation block exists with the corrected selection sampler off
  RoBMA.options(sampler.gates = FALSE, selection.sampler = "slice")
  .expect_declarations(robma, gates = FALSE, blocks = list(c(allocation_sd, split_w)),
                       label = "gates off")
  RoBMA.options(sampler.gates = TRUE, sampler.blocks = FALSE)
  .expect_declarations(robma, gates = TRUE, blocks = "none", label = "blocks off")
  .expect_declarations(selection, gates = FALSE, blocks = "none", label = "blocks off")

  # "jags.default" declares nothing whatever the two options say, also in a
  # model-averaged fit without a selection likelihood
  pet_peese <- RoBMA(yi = yi, vi = vi, data = dat, measure = "SMD", only_priors = TRUE,
                     prior_bias = list(BayesTools::prior_PET("cauchy", list(0, 1), list(0, Inf)),
                                       BayesTools::prior_PEESE("cauchy", list(0, 5), list(0, Inf))))
  RoBMA.options(sampler.gates = TRUE, sampler.blocks = TRUE, selection.sampler = "slice")
  .expect_declarations(robma, gates = TRUE, blocks = list(c(allocation_sd, split_w)), label = "slice")
  .expect_declarations(selection, gates = FALSE, blocks = list(c("mu", "tau", "omega")), label = "slice")
  .expect_declarations(pet_peese, gates = TRUE, blocks = "none", label = "slice PET-PEESE")
  RoBMA.options(selection.sampler = "jags.default")
  expect_true(RoBMA.get_option("sampler.gates"))
  expect_true(RoBMA.get_option("sampler.blocks"))
  .expect_declarations(robma, gates = FALSE, blocks = "none", label = "jags.default")
  .expect_declarations(selection, gates = FALSE, blocks = "none", label = "jags.default")
  .expect_declarations(pet_peese, gates = FALSE, blocks = "none", label = "jags.default PET-PEESE")

  # a named value is stored unnamed and behaves as the unnamed one
  RoBMA.options(selection.sampler = "slice")
  RoBMA.options(selection.sampler = c(stock = "jags.default"))
  expect_identical(RoBMA.get_option("selection.sampler"), "jags.default")
  expect_true(RoBMA.get_option("sampler.gates"))
  expect_true(RoBMA.get_option("sampler.blocks"))
  .expect_declarations(robma, gates = FALSE, blocks = "none", label = "named jags.default")
  .expect_declarations(selection, gates = FALSE, blocks = "none", label = "named jags.default")
  .expect_declarations(pet_peese, gates = FALSE, blocks = "none", label = "named jags.default PET-PEESE")
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
  expect_null(forwarded[["sampler_joint"]])
  # a selection fit collects the corrected sampler's record from every process
  expect_true(is.function(forwarded[["runtime_report"]]))

  RoBMA.options(sampler.gates = FALSE, sampler.blocks = FALSE)
  .fit(object)
  expect_true("sampler_blocks" %in% names(forwarded))
  expect_false(forwarded[["sampler_gates"]])
  expect_null(forwarded[["sampler_blocks"]])

  joint <- RoBMA(yi = yi, vi = vi, data = setup[["data"]], measure = "SMD", only_priors = TRUE)
  RoBMA.options(sampler.joint = "auto", sampler.nested.weights = FALSE,
                sampler.probabilities = "tau_indicator")
  .fit(joint)
  expect_true(is.function(forwarded[["sampler_joint"]]))
  expect_identical(forwarded[["jags_modules"]], "RoBMA")
  expect_identical(forwarded[["sampler_probabilities"]], "tau_indicator")
  expect_false(forwarded[["sampler_nested_weights"]])

  RoBMA.options(sampler.nested.weights = TRUE)
  .fit(joint)
  expect_true(forwarded[["sampler_nested_weights"]])
  expect_null(forwarded[["sampler_joint"]])

  # an extension replays the stored model and passes no declaration; it passes
  # the report explicitly (fits saved before the record store none)
  object[["fit"]] <- list("extension probe")
  .fit(object, extend = TRUE)
  expect_false(any(c("sampler_gates", "sampler_blocks", "sampler_joint", "sampler_nested_weights",
                     "sampler_probabilities") %in% names(forwarded)))
  expect_true(is.function(forwarded[["runtime_report"]]))

  # other fits pass no report
  ordinary <- brma(yi = yi, vi = vi, data = setup[["data"]], measure = "SMD", only_priors = TRUE)
  .fit(ordinary)
  expect_true("runtime_report" %in% names(forwarded))
  expect_null(forwarded[["runtime_report"]])
})


test_that(".fit() declares an ensemble's allocation block and records fallbacks without a warning", {

  old_options <- RoBMA.options()
  on.exit(do.call(RoBMA.options, old_options), add = TRUE)
  RoBMA.options(selection.sampler = "slice.coarse", sampler.blocks = TRUE)
  setup  <- .declaration_data()
  object <- RoBMA.mv(yi = yi, V = setup[["V"]], random = ~ 1 | study / esid, data = setup[["data"]],
                     measure = "SMD", only_priors = TRUE)
  prior_list <- .declaration_prior_list(object)
  # one process's report: one of four corrected samplers fell back
  report <- list(samplers = 4, fallback_few_proposals = 0, fallback_surrogate_exact = 0,
                 fallback_low_acceptance = 1, sampling_proposed = 500, sampling_accepted = 500,
                 min_sampler_rate = NA_real_)
  fake_fit <- function() structure(list(), runtime_report_values = list(list(report)))
  declared <- NULL
  testthat::local_mocked_bindings(
    JAGS_fit = function(...) {
      # JAGS_fit() calls the block function once on the combined prior list
      declared <<- list(...)[["sampler_blocks"]](prior_list)
      fake_fit()
    },
    # BayesTools carries the input fit's stored warnings into the extension
    JAGS_extend = function(fit, ...) structure(fake_fit(), warnings = attr(fit, "warnings", exact = TRUE)),
    JAGS_check_convergence = function(...) TRUE,
    .package = "BayesTools")

  fitted <- NULL
  expect_no_warning(fitted <- .fit(object))
  expect_identical(declared, list(c(allocation_sd, split_w)))
  record <- attr(fitted, "selection_sampler")
  expect_false("allocation_block_withheld" %in% names(record))
  expect_identical(record[["fallback_samplers"]], 1)
  expect_null(fitted[["warnings"]])

  # an extension of a fit of a development version drops its stored fallback message
  stored <- paste0("1 of 4 corrected selection samplers switched to the exact slice sampler at the end of ",
                   "adaptation (1 for low acceptance). This model-averaged fit does not declare ...")
  attr(fitted, "warnings") <- c("an earlier warning", stored)
  object[["fit"]] <- fitted
  extended <- NULL
  expect_no_warning(extended <- .fit(object, extend = TRUE))
  expect_identical(attr(extended, "warnings"), "an earlier warning")
  expect_identical(extended[["warnings"]], "an earlier warning")
})


test_that("an error of a BayesTools sampler names the option that turns it off", {

  setup  <- .declaration_data()
  object <- bselmodel(yi = yi, vi = vi, data = setup[["data"]], measure = "SMD", only_priors = TRUE)
  fit_error <- function(message) {
    testthat::local_mocked_bindings(
      JAGS_fit = function(...) simpleError(message),
      .package = "BayesTools")
    .fit(object)[["errors"]]
  }

  messages <- c(
    blocks = "Error in node mu\nFactor slice sampler: failure to calculate the log density of the block at a proposal",
    gates  = "Error in node PET_1\ngated node changes the likelihood while its gate is closed"
  )
  for (option in names(messages)) {
    errors <- fit_error(messages[[option]])
    expect_length(errors, 1L)
    # the original message stays, followed by the hint
    expect_true(startsWith(errors, messages[[option]]), info = option)
    expect_true(grepl(paste0("RoBMA.options(sampler.", option, " = FALSE)"), errors, fixed = TRUE),
                info = option)
    expect_true(grepl("please report the model", errors, fixed = TRUE), info = option)
  }

  # other errors are passed on unchanged
  expect_identical(fit_error("Error in node y\nFailure to calculate log density"),
                   "Error in node y\nFailure to calculate log density")
})
