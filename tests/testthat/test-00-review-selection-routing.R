test_that("F3-4 sole bias branches use complete backend monitors without an indicator", {

  wf <- BayesTools::prior_weightfunction("one-sided", .025, BayesTools::wf_fixed(c(1, .5)))
  for (kind in c("null", "false", "bsel", "two")) {
    args <- list(yi = c(.1, .2, .3, .4), sei = rep(.2, 4), measure = "GEN",
      prior_unit_information_sd = 1, prior_bias = wf, only_priors = TRUE, silent = TRUE)
    if (kind == "null") args["prior_bias_null"] <- list(NULL)
    if (kind == "false") args$prior_bias_null <- FALSE
    setup <- review_compile_object(do.call(if (kind == "bsel") bselmodel else RoBMA, args))
    priors <- .create_fit_priors(setup$object$data, setup$object$priors)
    primitive <- if (kind == "bsel") cbind(mu = rep(.1, 24), tau = .15) else
      cbind(mu_indicator = rep(1:2, 12), mu_component_2 = .1,
        tau_indicator = rep(c(1, 1, 2, 2), 6), tau_component_2 = .15)
    if (kind == "two") primitive <- cbind(primitive, bias_indicator = rep(rep(1:2, each = 4), 3))
    preliminary <- review_attach_replay(setup, primitive)
    generated <- BayesTools::JAGS_evaluate_deterministic(preliminary$fit, draws = primitive)
    samples <- cbind(primitive, generated)
    monitor <- unlist(lapply(BayesTools::JAGS_to_monitor(priors), function(name) {
      if (name == "omega") c("omega[1]", "omega[2]") else name
    }))
    expect_true(all(monitor %in% colnames(samples)))
    samples <- samples[, monitor, drop = FALSE]
    object <- review_attach_replay(setup, samples)
    expect_silent(BayesTools::JAGS_validate_fit_contract(object$fit,
      requires = c("name_encoding", "formula_name_map", "formula_design", "parameter_map", "draw_geometry")))
    expect_identical(class(object), class(setup$object))
    route <- .selection_row_routing(object$priors, samples)
    expect_identical(route$bias_indicator, if (kind == "two") as.integer(samples[, "bias_indicator"]) else rep(1L, 24))
    withr::local_seed(26010)
    actual <- stats::predict(object, type = "response", quiet = TRUE)
    expect_identical(dim(actual), c(24L, 4L))
    expect_true(all(is.finite(actual)))
    set.seed(26010)
    expect_identical(as.vector(actual), as.vector(stats::predict(object, type = "response", quiet = TRUE)))
    if (kind %in% c("null", "false")) {
      expect_false("bias_indicator" %in% monitor)
      supplied <- cbind(samples, bias_indicator = 1)
      expect_identical(.selection_row_routing(object$priors, supplied)$bias_indicator, rep(1L, 24))
      for (value in c(.5, NaN, Inf, 2)) {
        supplied[1L, "bias_indicator"] <- value
        expect_error(.selection_row_routing(object$priors, supplied), "posterior model indicator")
      }
      for (mu_branch in 1:2) for (tau_branch in 1:2) {
        row <- which(samples[, "mu_indicator"] == mu_branch & samples[, "tau_indicator"] == tau_branch)[1L]
        matched <- review_attach_replay(setup, samples[rep(row, 24L), , drop = FALSE])
        control_args <- args
        control_args$prior_bias_null <- NULL
        control_args$prior_effect <- priors$mu[[mu_branch]]
        control_args$prior_heterogeneity <- priors$tau[[tau_branch]]
        control_setup <- review_compile_object(do.call(bselmodel, control_args))
        control_samples <- cbind(mu = rep(samples[row, "mu"], 24L), tau = samples[row, "tau"])
        initial <- review_attach_replay(control_setup, control_samples)
        generated <- BayesTools::JAGS_evaluate_deterministic(initial$fit, draws = control_samples)
        control_samples <- cbind(control_samples, generated[, setdiff(colnames(generated), colnames(control_samples)), drop = FALSE])
        control <- review_attach_replay(control_setup, control_samples)
        expect_identical(unname(control_samples[, "mu"]), unname(samples[rep(row, 24L), "mu"]))
        expect_identical(unname(control_samples[, "tau"]), unname(samples[rep(row, 24L), "tau"]))
        expect_true(all(c("mu", "tau", "omega[1]", "omega[2]") %in% colnames(control_samples)))
        set.seed(26010)
        target_prediction <- stats::predict(matched, type = "response", quiet = TRUE)
        set.seed(26010)
        control_prediction <- stats::predict(control, type = "response", quiet = TRUE)
        expect_lte(max(abs(as.vector(target_prediction) - as.vector(control_prediction))), 1e-12)
      }
    }
    if (kind == "two") {
      expect_error(.selection_row_routing(object$priors,
        samples[, colnames(samples) != "bias_indicator", drop = FALSE]),
        "Missing posterior model indicator: 'bias_indicator'.", fixed = TRUE)
      expect_identical(route$use_normal, samples[, "bias_indicator"] == 1)
    }
  }
  normal <- list(outcome = list(bias = BayesTools::prior_none()))
  expect_true(all(.selection_row_routing(normal, matrix(0, 3L, 0L))$use_normal))
})
