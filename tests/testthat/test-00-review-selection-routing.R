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
