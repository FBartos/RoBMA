test_that("a random inclusion gate belongs only to its declared block", {

  result <- BayesTools::JAGS_formula(
    formula = ~ 1 + random(1 | study, name = "study", covariance = "diag") +
      random(1 | esid, name = "esid", covariance = "diag"),
    parameter = "mu",
    data = data.frame(study = factor(c("a", "a", "b", "b")), esid = factor(1:4)),
    prior_list = list(intercept = BayesTools::prior("normal", list(0, 1))),
    prior_random = BayesTools::prior_random(
      sd = BayesTools::prior("gamma", list(2, 2)),
      allocation = list(BayesTools::random_variance_allocation(
        name = "gated", terms = c(study = "study"),
        sd = BayesTools::prior("gamma", list(2, 2)),
        inclusion = list(study = BayesTools::prior("spike", list(location = .5)))
      ))
    )
  )
  design <- result[["formula_design"]]
  allocation <- design[["random_allocations"]][["gated"]]
  source <- allocation[["source_node"]]
  gate <- allocation[["inclusion"]][["study"]][["indicator_name"]]
  sd_name <- design[["random_effects"]][[2L]][["sd_parameter_names"]]
  samples <- cbind(
    mu_intercept = c(-.1, .1, .2, .3), source = c(.2, .3, .4, .5),
    gate = c(0, 1, 0, 1), sd = c(.6, .7, .8, .9)
  )
  colnames(samples)[2:4] <- c(source, gate, sd_name)
  fit <- coda::mcmc.list(coda::mcmc(samples))
  class(fit) <- c("BayesTools_fit", class(fit))
  attr(fit, "prior_list") <- result[["prior_list"]]
  attr(fit, "formula_design") <- list(mu = design)
  fit <- BayesTools:::.bt_attach_parameter_map(fit)
  fit <- BayesTools:::.bt_attach_draw_geometry(fit)
  fit <- BayesTools:::.bt_attach_fit_contract(fit)
  object <- structure(list(
    fit = fit, data = structure(list(), random = TRUE)
  ), class = c("RoBMA", "brma.mv", "brma"))

  for (quantity in c("tau", "tau2")) {
    study_parameter <- paste0("study: ", quantity)
    esid_parameter <- paste0("esid: ", quantity)

    study_posterior <- .brma_random_parameter_mixed_posterior(object, study_parameter)[[1L]]
    esid_posterior <- .brma_random_parameter_mixed_posterior(object, esid_parameter)[[1L]]
    study_atoms <- BayesTools::posterior_metadata(study_posterior, "atoms")
    esid_atoms <- BayesTools::posterior_metadata(esid_posterior, "atoms")
    expect_equal(as.numeric(study_atoms[["locations"]]), 0)
    expect_equal(study_atoms[["mass"]], .5)
    expect_length(esid_atoms[["mass"]], 0L)
    expect_equal(as.numeric(esid_posterior), if (quantity == "tau") {
      samples[, sd_name]
    } else {
      samples[, sd_name]^2
    })

    # Conditioning keeps the draws with the study gate on; the ungated esid
    # block has no inclusion event.
    study_conditional <- .brma_random_parameter_mixed_posterior(
      object, study_parameter, conditional = TRUE
    )[[1L]]
    expect_length(study_conditional, 2L)
    expect_true(BayesTools::posterior_atoms_free(study_conditional))
    expect_error(
      .brma_random_parameter_mixed_posterior(object, esid_parameter, conditional = TRUE),
      "the quantity has no inclusion gate"
    )
  }
})


test_that("variance proportions require a possible positive parent allocation", {

  # Parent (root) gate (0, 0, 1, 1) and child study gate (0, 1, 0, 1): the
  # study share is defined only where the parent allocation is included,
  # also where the child gate is on.
  object  <- shared_gate_random_object(
    root_gate  = c(0, 0, 1, 1),
    study_gate = c(0, 1, 0, 1)
  )
  selected <- .brma_random_parameter_select(object, "(mu) split: tau2_prop(study)")
  state <- .iwmde_gate_states(
    list(object = object, posterior_samples = as.matrix(object[["fit"]][[1L]])),
    list(gate_selection = .brma_random_parameter_gate_selection(object, selected))
  )
  expect_identical(state[["defined"]], c(FALSE, FALSE, TRUE, TRUE))
  expect_identical(state[["point_zero"]], c(FALSE, FALSE, TRUE, FALSE))
  expect_identical(state[["continuous"]], c(FALSE, FALSE, FALSE, TRUE))
})


test_that("conditional random-effect plots refuse qCMDE/IWMDE with the density-method classes", {

  # A known-V stand-in passes the qCMDE/IWMDE capability check of the model;
  # its conditional random-effect plots are KDE-only.
  object <- shared_gate_random_object()
  attr(object[["data"]], "known_V") <- TRUE
  attr(object[["data"]], "measure") <- "GEN"
  expect_true(.iwmde_capability(
    object = object, density_method = "IWMDE"
  )[["available"]])
  for (method in c("qCMDE", "IWMDE")) {
    for (class in c("RoBMA_density_method_conditional_random",
                    "RoBMA_density_method_unavailable")) {
      expect_error(
        plot(object, "study: tau", component = "random", conditional = TRUE,
             density_method = method, plot_type = "ggplot"),
        class = class,
        info  = method
      )
    }
    # hypothesis() refuses the same cause with its own classes followed by
    # the density-method classes.
    error <- tryCatch(
      hypothesis(object, "`(mu) study: tau(intercept)` = 0.3",
                 component = "random", conditional = TRUE,
                 density_method = method, n_samples = 1000L, seed = 1),
      error = identity
    )
    expect_identical(
      class(error),
      c("RoBMA_hypothesis_method", "RoBMA_hypothesis_unavailable",
        "RoBMA_density_method_conditional_random",
        "RoBMA_density_method_unavailable", "error", "condition"),
      info = method
    )
  }
})


test_that("random-effect quantities without a scalar coordinate refuse qCMDE/IWMDE plots with the density-method classes", {

  # The variance of a random component has no scalar random-component
  # coordinate for a qCMDE/IWMDE density curve.
  object <- shared_gate_random_object()
  attr(object[["data"]], "known_V") <- TRUE
  attr(object[["data"]], "measure") <- "GEN"
  for (method in c("qCMDE", "IWMDE")) {
    error <- tryCatch(
      plot(object, "(mu) study: tau2(intercept)", component = "random",
           density_method = method, plot_type = "ggplot"),
      error = identity
    )
    expect_identical(
      class(error),
      c("RoBMA_density_method_random_target",
        "RoBMA_density_method_unavailable", "error", "condition"),
      info = method
    )
    expect_match(
      conditionMessage(error),
      "because it has no supported scalar random-component coordinate",
      fixed = TRUE,
      info  = method
    )
  }
})
