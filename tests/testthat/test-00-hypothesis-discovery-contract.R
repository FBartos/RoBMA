test_that("hypothesis plans refuse methods by the runtime qCMDE/IWMDE capability", {

  # A point plan whose targets are eligible: every method is available up to
  # the capability of the fitted model.
  plan <- list(point = TRUE, refusal = NULL, targets = list(), method_refusals = list())
  glmm <- structure(list(), class = c("brma.glmm", "brma"))
  # The engine checks a plan against the density backend and the model's
  # capability on every call: the provider names the model (a stand-in here)
  # and the density methods it offers.
  provider <- list(
    object            = glmm,
    capability_object = glmm,
    density_methods   = c("KDE", "normal", "qCMDE", "IWMDE")
  )
  status <- function(plan, method) engine_plan_status(plan, glmm, method, provider)
  check <- function(plan, method) {
    BayesTools:::.hypothesis_model_plan_check(
      plan, provider, BayesTools:::.hypothesis_model_backend(method)
    )
  }
  plan <- BayesTools:::.hypothesis_model_plan_finish(plan)
  capability <- .iwmde_capability(object = glmm, density_method = "IWMDE")

  expect_null(status(plan, "KDE"))
  expect_null(status(plan, "qCMDE"))
  expect_null(status(plan, "normal"))
  expect_identical(
    capability[["class"]],
    c("RoBMA_density_method_glmm", "RoBMA_density_method_unavailable")
  )
  # A capability refusal keeps the hypothesis classes and also carries the
  # classes of the capability refusal: the cause and the parent class.
  expect_identical(
    status(plan, "IWMDE"),
    list(
      reason = capability[["reason"]],
      class  = c("BayesTools_hypothesis_method", "BayesTools_hypothesis_unavailable",
                 "RoBMA_density_method_glmm",
                 "RoBMA_density_method_unavailable")
    )
  )
  expect_error(
    check(plan, "IWMDE"),
    capability[["reason"]],
    fixed = TRUE,
    class = "BayesTools_hypothesis_method"
  )
  for (class in c("RoBMA_density_method_glmm",
                  "RoBMA_density_method_unavailable")) {
    expect_error(
      check(plan, "IWMDE"),
      class = class
    )
  }
  # Outside hypothesis(), the capability check stops with the classes of the
  # capability refusal.
  expect_error(
    .iwmde_check_density_method_supported(glmm, "IWMDE"),
    class = "RoBMA_density_method_glmm"
  )
  expect_error(
    .iwmde_check_density_method_supported(glmm, "IWMDE"),
    class = "RoBMA_density_method_unavailable"
  )
  expect_invisible(.iwmde_check_density_method_supported(glmm, "qCMDE"))
  rendered <- BayesTools:::.hypothesis_model_quantities_render_plans(
    list(point = list(plan), region = list(plan), contrast = NULL),
    bracket = FALSE, provider = provider
  )
  expect_identical(rendered[["point_test_methods"]], "KDE, qCMDE")
  expect_identical(rendered[["reason"]], capability[["reason"]])
  expect_true(is.na(rendered[["contrast_test"]]))

  # The capability of the fitted model comes before the quantity's own
  # refusal of the method, as in plot(); a method the model supports keeps
  # the quantity's refusal.
  own <- BayesTools:::.hypothesis_model_refusal("The quantity refuses the method.", "method")
  plan_own <- BayesTools:::.hypothesis_model_plan_finish(
    list(point = TRUE, refusal = NULL, targets = list(),
         method_refusals = list(qCMDE = own, IWMDE = own))
  )
  expect_identical(status(plan_own, "qCMDE"), own)
  expect_identical(
    status(plan_own, "IWMDE"),
    status(plan, "IWMDE")
  )
  expect_null(status(plan_own, "KDE"))
  # An object without the fitted model the method needs (a marginal-means
  # object without its source fit) is refused before the capability.
  missing_source <- BayesTools:::.hypothesis_model_refusal("The source fit is missing.", "target")
  plan_source <- BayesTools:::.hypothesis_model_plan_finish(
    list(point = TRUE, refusal = NULL, targets = list(),
         method_refusals = list(qCMDE = own, IWMDE = own),
         source_refusals = list(ordinate = missing_source))
  )
  for (method in c("qCMDE", "IWMDE")) {
    expect_identical(status(plan_source, method), missing_source,
                     info = method)
  }

  # A refusal of the whole plan (a fixed quantity) applies to point and
  # region statements and every method.
  fixed <- BayesTools:::.hypothesis_model_plan_finish(
    list(point = FALSE, refusal = BayesTools:::.hypothesis_model_plan_fixed_refusal("mu"),
         targets = list(), method_refusals = list())
  )
  rendered <- BayesTools:::.hypothesis_model_quantities_render_plans(
    list(point = list(fixed), region = list(fixed), contrast = NULL),
    bracket = FALSE, provider = provider
  )
  expect_false(rendered[["point_test"]])
  expect_false(rendered[["direction_test"]])
  expect_identical(rendered[["point_test_methods"]], "")
  expect_match(rendered[["reason"]], "fixed by the fitted model", fixed = TRUE)
})


test_that("hypothesis discovery shares the runtime qCMDE/IWMDE capability", {

  object <- single_sd_random_object(BayesTools::prior("gamma", list(2, 2)))
  capability <- .iwmde_capability(
    object         = object,
    density_method = "qCMDE"
  )
  out <- hypothesis_quantities(object)
  intercept <- out[out[["parameter"]] == "mu_intercept", , drop = FALSE]

  expect_false(capability[["available"]])
  expect_identical(unique(intercept[["point_test_methods"]]), "KDE")
  expect_identical(unique(intercept[["reason"]]), capability[["reason"]])
  expect_error(
    hypothesis(object, "mu_intercept = 0", density_method = "qCMDE"),
    capability[["reason"]],
    fixed = TRUE,
    class = "BayesTools_hypothesis_method"
  )
  expect_error(
    hypothesis(object, "mu_intercept = 0", density_method = "qCMDE"),
    class = "RoBMA_density_method_unavailable"
  )
  classes <- c("RoBMA_density_method_random_unknown_v",
               "RoBMA_density_method_unavailable")
  expect_identical(capability[["class"]], classes)
  expect_error(
    .check_iwmde_available(object, "qCMDE/IWMDE hypothesis()"),
    capability[["reason"]],
    fixed = TRUE
  )
  for (class in classes) {
    expect_error(
      hypothesis(object, "mu_intercept = 0", density_method = "qCMDE"),
      class = class
    )
    expect_error(
      .check_iwmde_available(object, "qCMDE/IWMDE hypothesis()"),
      class = class
    )
    expect_error(
      .iwmde_check_density_method_supported(object, "qCMDE"),
      class = class
    )
    expect_error(
      plot(object, "mu", density_method = "qCMDE"),
      class = class
    )
  }

  known_v_object <- object
  attr(known_v_object[["data"]], "known_V") <- TRUE
  expect_true(.iwmde_capability(
    object         = known_v_object,
    density_method = "qCMDE"
  )[["available"]])
  expect_invisible(
    .check_iwmde_available(known_v_object, "qCMDE/IWMDE hypothesis()")
  )

  # Objects without a fit cannot be tested.
  unfitted <- hypothesis_quantities(structure(list(), class = "brma"))
  expect_false(any(unfitted[["point_test"]]))
  expect_identical(unique(unfitted[["reason"]]), "'hypothesis' requires a fitted brma object.")
})


test_that("qCMDE/IWMDE are refused for objects that are no RoBMA models", {

  # With RoBMA loaded, a BayesTools fit is described by the default provider of
  # BayesTools, which offers the registered ordinate backends of RoBMA.
  model <- single_sd_random_object(BayesTools::prior("gamma", list(2, 2)))
  fit   <- model[["fit"]]
  expect_s3_class(fit, "BayesTools_fit")
  expect_false(inherits(fit, "brma"))
  expect_true(all(
    c("KDE", "qCMDE", "IWMDE") %in%
      BayesTools:::.hypothesis_model_quantities_methods(engine_provider(fit))
  ))

  # Both backends refuse the object as a method, naming the supported objects.
  plan <- BayesTools:::.hypothesis_model_plan_finish(
    list(point = TRUE, refusal = NULL, targets = list(), method_refusals = list())
  )
  for (method in c("qCMDE", "IWMDE")) {
    backend <- BayesTools:::.hypothesis_model_backend(method)
    for (object in list(fit, list(), NULL)) {
      refusal <- backend[["available"]](object, plan, NULL)
      expect_identical(refusal[["type"]], "method", info = method)
      expect_match(
        refusal[["reason"]],
        "qCMDE/IWMDE evaluate ordinates of RoBMA models (fitted 'brma' objects and their marginal means)",
        fixed = TRUE, info = method
      )
      expect_match(refusal[["reason"]], "Use density_method = 'KDE'.", fixed = TRUE, info = method)
      expect_identical(
        refusal[["class"]],
        c("RoBMA_density_method_model_class", "RoBMA_density_method_unavailable"),
        info = method
      )
    }
    expect_match(backend[["available"]](fit, plan, NULL)[["reason"]],
                 "class 'BayesTools_fit'", fixed = TRUE, info = method)
    # The engine stops with that refusal before any statement is evaluated.
    expect_identical(
      engine_plan_status(plan, fit, method)[["reason"]],
      backend[["available"]](fit, plan, NULL)[["reason"]],
      info = method
    )
    expect_error(
      BayesTools::hypothesis_BF(fit, "mu_intercept = 0", density_method = method, seed = 1),
      "evaluate ordinates of RoBMA models", fixed = TRUE,
      class = "BayesTools_hypothesis_method", info = method
    )
    expect_error(
      BayesTools::hypothesis_BF(fit, "mu_intercept = 0", density_method = method, seed = 1),
      class = "RoBMA_density_method_unavailable", info = method
    )
  }

  # Discovery lists the methods that evaluate the fit, none of qCMDE/IWMDE.
  quantities <- hypothesis_quantities(fit)
  expect_gt(nrow(quantities), 0L)
  expect_true(all(quantities[["point_test"]]))
  expect_false(any(grepl("qCMDE|IWMDE", quantities[["point_test_methods"]])))
  expect_true(all(grepl("^KDE", quantities[["point_test_methods"]])))
  expect_match(unique(quantities[["reason"]]), "evaluate ordinates of RoBMA models", fixed = TRUE)
  expect_s3_class(
    BayesTools::hypothesis_BF(fit, "mu_intercept = 0", density_method = "KDE", seed = 1),
    "BayesTools_hypothesis_BF"
  )

  # The RoBMA models keep their methods: the fitted brma objects (the capability
  # objects of both providers) are available up to their own capability.
  known_v_model <- model
  attr(known_v_model[["data"]], "known_V") <- TRUE
  glmm <- structure(list(), class = c("brma.glmm", "brma"))
  for (method in c("qCMDE", "IWMDE")) {
    backend <- BayesTools:::.hypothesis_model_backend(method)
    expect_null(backend[["available"]](known_v_model, plan, NULL), info = method)
    expect_identical(
      backend[["available"]](model, plan, NULL)[["class"]],
      .iwmde_capability(model, density_method = method)[["class"]],
      info = method
    )
  }
  expect_null(BayesTools:::.hypothesis_model_backend("qCMDE")[["available"]](glmm, plan, NULL))
  expect_identical(
    BayesTools:::.hypothesis_model_backend("IWMDE")[["available"]](glmm, plan, NULL)[["class"]],
    c("RoBMA_density_method_glmm", "RoBMA_density_method_unavailable")
  )
  expect_match(
    unique(hypothesis_quantities(known_v_model)[["point_test_methods"]])[[1L]],
    "KDE, qCMDE, IWMDE", fixed = TRUE
  )
})


test_that("qCMDE/IWMDE are unavailable for component-specific scale formulas", {

  # The IWMDE context evaluates one scale formula ('log_tau'): a named
  # 'scale' list gives one formula ('log_tau_<component>') per random
  # component, with one or several entries.
  data <- data.frame(
    yi     = c(0.08, 0.13, 0.18, 0.20, 0.01, 0.05),
    study  = rep(c("s1", "s2", "s3"), each = 2L),
    effect = rep(c("a", "b"), 3L),
    x      = c(0, 1, 0, 1, 0, 1)
  )
  V <- kronecker(diag(3L), matrix(c(0.04, 0.018, 0.018, 0.05), nrow = 2L))
  prior_object <- function(random, scale) {
    suppressWarnings(brma.mv(
      yi = yi, V = V, data = data, measure = "GEN", random = random,
      scale = scale, prior_unit_information_sd = 1, only_priors = TRUE
    ))
  }
  objects <- list(
    plain   = prior_object(list(study = ~ 1 | study), ~ x),
    one     = prior_object(list(study = ~ 1 | study), list(study = ~ x)),
    several = prior_object(
      list(study = ~ 1 | study, effect = ~ 1 | study:effect),
      list(study = ~ x, effect = ~ x)
    )
  )
  reasons <- c(
    one     = paste0(
      "qCMDE/IWMDE density estimation is unavailable for models with a ",
      "component-specific scale formula. Use density_method = 'KDE'."
    ),
    several = paste0(
      "qCMDE/IWMDE density estimation is unavailable for models with ",
      "several scale formulas. Use density_method = 'KDE'."
    )
  )

  for (method in c("qCMDE", "IWMDE")) {
    expect_true(.iwmde_capability(
      object = objects[["plain"]], density_method = method
    )[["available"]], info = method)
    for (name in names(reasons)) {
      expect_identical(
        .iwmde_capability(object = objects[[name]], density_method = method),
        list(
          available = FALSE,
          reason    = reasons[[name]],
          class     = c("RoBMA_density_method_scale_components",
                        "RoBMA_density_method_unavailable")
        ),
        info = paste(name, method)
      )
    }
  }
  # Every caller building the IWMDE context stops with the same reason, not
  # with the error of evaluating a missing scale formula.
  plain_inputs <- .iwmde_formula_inputs(
    objects[["plain"]][["data"]], objects[["plain"]][["priors"]]
  )
  expect_s3_class(plain_inputs[["scale"]][["formula"]], "formula")
  for (name in names(reasons)) {
    expect_error(
      .iwmde_formula_inputs(objects[[name]][["data"]], objects[[name]][["priors"]]),
      reasons[[name]],
      fixed = TRUE,
      info  = name
    )
    for (class in c("RoBMA_density_method_scale_components",
                    "RoBMA_density_method_unavailable")) {
      expect_error(
        .iwmde_formula_inputs(objects[[name]][["data"]], objects[[name]][["priors"]]),
        class = class,
        info  = name
      )
    }
  }
})


test_that("random discovery uses the authoritative likelihood-aware target", {

  object <- single_sd_random_object(BayesTools::prior("gamma", list(2, 2)))
  attr(object[["data"]], "known_V") <- TRUE
  parameter <- "(mu) tau(intercept)"
  # The plans are built with the mocked density target, so none is kept.
  plan_of <- function() {
    BayesTools::hypothesis_cache_clear()
    engine_plans(object, paste0("`", parameter, "` = 0.3"))[[1L]]
  }
  withr::defer(BayesTools::hypothesis_cache_clear())

  # A quantity without a scalar random-component coordinate: the refusal
  # of the density target (as plot() raises it) follows the hypothesis
  # classes.
  testthat::local_mocked_bindings(
    .brma_random_parameter_density_target = function(object, parameter, ...) {
      .iwmde_unavailable("unsupported", "random_target")
    },
    .package = "RoBMA"
  )
  unsupported <- plan_of()
  expect_null(engine_plan_status(unsupported, object, "KDE"))
  for (method in c("qCMDE", "IWMDE")) {
    expect_identical(
      engine_plan_status(unsupported, object, method),
      list(reason = "unsupported",
           class  = c("BayesTools_hypothesis_target", "BayesTools_hypothesis_unavailable",
                      "RoBMA_density_method_random_target",
                      "RoBMA_density_method_unavailable"),
           cause  = "random_target")
    )
  }

  testthat::local_mocked_bindings(
    .brma_random_parameter_density_target = function(object, parameter, ...) {
      list(parameter = "rho_raw", parameter_spec = list(type = "primitive"))
    },
    .package = "RoBMA"
  )
  supported <- plan_of()
  for (method in c("KDE", "qCMDE", "IWMDE")) {
    expect_null(engine_plan_status(supported, object, method), info = method)
  }
})


# A marginal mean of a fixture: draws with an exact normal prior density and
# declared atoms (none, or a point mass carrying all mass for a fixed mean).
.discovery_marginal_mean <- function(values, fixed = FALSE) {

  with_draw_metadata(
    values,
    class         = c("marginal_posterior.simple", "marginal_posterior"),
    prior_density = BayesTools::prior("normal", list(0, 1)),
    atoms         = if (fixed) {
      BayesTools::posterior_atom_attribute(point_masses = data.frame(x = 0, mass = 1))
    } else {
      BayesTools::posterior_atom_attribute()
    }
  )
}


.discovery_marginal_intercept <- function() {

  object <- review_factor_replay(c("a", "b"), "treatment", distribution = "normal")
  means <- marginal_means(object, bf = FALSE, n_samples = 128L)
  means[["inference"]][["averaged"]][["mu_intercept"]]
}


.discovery_marginal_object <- function(term_map, means, source_object) {

  structure(list(
    term_map      = term_map,
    inference     = list(averaged = means, conditional = means),
    source_object = source_object
  ), class = "marginal_means.brma")
}


test_that("marginal discovery shares its source-model IWMDE capability", {

  source_object <- structure(
    list(data = structure(list(), random = TRUE)),
    class = c("brma.mv", "brma")
  )
  object <- .discovery_marginal_object(
    term_map      = data.frame(term = "intercept", parameter = "mu_intercept"),
    means         = list(mu_intercept = .discovery_marginal_intercept()),
    source_object = source_object
  )

  out <- hypothesis_quantities(object)

  # A source object without a fit cannot compute qCMDE/IWMDE ordinates; that
  # refusal precedes the capability of the source model (a random-formula
  # model without known V), as in the marginal-means plot().
  expect_identical(out[["point_test_methods"]], "KDE")
  expect_match(out[["reason"]], "does not contain the source fitted brma", fixed = TRUE)

  source_object[["fit"]] <- list(fitted = TRUE)
  object[["source_object"]] <- source_object
  out <- hypothesis_quantities(object)
  capability <- .iwmde_capability(object = source_object, density_method = "IWMDE")
  expect_identical(out[["point_test_methods"]], "KDE")
  expect_identical(out[["reason"]], capability[["reason"]])
})


test_that("marginal-means discovery reports selectable levels", {

  object <- .discovery_marginal_object(
    term_map = data.frame(
      term      = c("intercept", "group"),
      parameter = c("mu_intercept", "mu_group")
    ),
    means = list(
      mu_intercept = .discovery_marginal_intercept(),
      mu_group     = list(
        a = .discovery_marginal_mean(c(-1, 0)),
        b = .discovery_marginal_mean(c(0, 1))
      )
    ),
    source_object = structure(list(fit = list(fitted = TRUE)), class = c("brma.glmm", "brma"))
  )

  out <- hypothesis_quantities(object)
  scalar <- out[out[["parameter"]] == "mu_intercept", , drop = FALSE]
  grouped <- out[out[["parameter"]] == "mu_group", , drop = FALSE]

  expect_identical(scalar[["bracket"]], "mu_intercept[intercept]")
  expect_identical(grouped[["bracket"]], c("mu_group[a]", "mu_group[b]"))
  expect_identical(unique(out[["point_test_methods"]]), "KDE, qCMDE")
  expect_true(is.na(scalar[["contrast_test"]]))
  expect_identical(grouped[["contrast_test"]], c(FALSE, FALSE))
})


test_that("marginal-means discovery rejects a fixed no-intercept scalar", {

  object <- .discovery_marginal_object(
    term_map = data.frame(term = "intercept", parameter = "mu_intercept"),
    means    = list(mu_intercept = list(
      intercept = .discovery_marginal_mean(rep(0, 20L), fixed = TRUE)
    )),
    source_object = structure(list(fit = list(fitted = TRUE)), class = "brma")
  )

  out <- hypothesis_quantities(object)

  expect_identical(out[["bracket"]], "mu_intercept[intercept]")
  expect_false(out[["point_test"]])
  expect_false(out[["direction_test"]])
  expect_identical(out[["point_test_methods"]], "")
  expect_match(out[["reason"]], "fixed by the fitted model")
  expect_error(
    hypothesis(object, "mu_intercept[intercept] > 0"),
    class = "BayesTools_hypothesis_fixed"
  )
})


test_that("marginal-means discovery keeps sampled siblings of a fixed level", {

  # The fixed level is declared by its atoms, not by constant draws: the
  # sampled levels' constant draws do not fix them.
  object <- .discovery_marginal_object(
    term_map = data.frame(term = "x", parameter = "mu_x"),
    means    = list(mu_x = list(
      `-1SD` = .discovery_marginal_mean(c(-1, -2, -3)),
      `0SD`  = .discovery_marginal_mean(c(0, 0, 0), fixed = TRUE),
      `1SD`  = .discovery_marginal_mean(c(1, 1, 1))
    )),
    source_object = structure(list(fit = list(fitted = TRUE)), class = "brma")
  )

  out <- hypothesis_quantities(object)
  fixed <- out[["bracket"]] == "mu_x[0SD]"

  expect_identical(
    out[["bracket"]],
    c("mu_x[-1SD]", "mu_x[0SD]", "mu_x[1SD]")
  )
  expect_identical(out[["point_test"]], c(TRUE, FALSE, TRUE))
  expect_identical(out[["direction_test"]], c(TRUE, FALSE, TRUE))
  expect_identical(
    out[["point_test_methods"]][!fixed],
    rep("KDE, qCMDE, IWMDE", 2L)
  )
  expect_identical(out[["point_test_methods"]][fixed], "")
  expect_match(out[["reason"]][fixed], "fixed by the fitted model")
})


