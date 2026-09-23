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
    study <- .brma_random_parameter_select(object, study_parameter)
    esid <- .brma_random_parameter_select(object, esid_parameter)
    expect_identical(.brma_random_parameter_inclusion_indicator(object, study), gate)
    expect_null(.brma_random_parameter_inclusion_indicator(object, esid))

    study_posterior <- .brma_random_parameter_mixed_posterior(object, study_parameter)[[1L]]
    esid_posterior <- .brma_random_parameter_mixed_posterior(object, esid_parameter)[[1L]]
    study_atoms <- attr(study_posterior, "posterior_atoms", exact = TRUE)
    esid_atoms <- attr(esid_posterior, "posterior_atoms", exact = TRUE)
    expect_equal(as.numeric(study_atoms[["locations"]]), 0)
    expect_equal(study_atoms[["mass"]], .5)
    expect_length(esid_atoms[["mass"]], 0L)
    expect_equal(as.numeric(esid_posterior), if (quantity == "tau") {
      samples[, sd_name]
    } else {
      samples[, sd_name]^2
    })
    expect_error(
      .brma_random_parameter_mixed_posterior(object, esid_parameter, conditional = TRUE),
      "owned by one independently gated component"
    )
  }
})


test_that("variance proportions require a possible positive parent allocation", {

  gate_prior <- function(indicator, probability) {
    prior <- BayesTools::prior("spike", list(location = probability))
    attr(prior, "random_allocation_indicator") <- indicator
    prior
  }
  fit <- structure(list(), prior_list = list(
    parent = gate_prior("parent", 0),
    child = gate_prior("child", .5)
  ))
  metadata <- list(
    quantity = "var_prop", index = 1L,
    component_indicators = c("child", NA_character_),
    parent_indicators = "parent"
  )
  prior <- .brma_random_parameter_allocation_gate_prior(list(fit = fit), metadata)
  expect_equal(prior[["continuous_mass"]], 0)
  expect_equal(nrow(prior[["points"]]), 0L)

  state <- .brma_random_parameter_allocation_gate_state(
    metadata, cbind(parent = c(0, 0), child = c(0, 1))
  )
  expect_false(any(state[["defined"]]))
})
