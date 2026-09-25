# A fitted-object stand-in with a gated total-variance allocation of a
# half-normal SD over 'study' (inclusion gate with prior probability 0.5) and
# 'esid' (ungated), Dirichlet(1, 1) shares, and synthetic draws.
.gated_random_object <- function(n = 200L) {

  result <- BayesTools::JAGS_formula(
    formula = ~ 1 + random(1 | study, name = "study", covariance = "diag") +
      random(1 | esid, name = "esid", covariance = "diag"),
    parameter = "mu",
    data = data.frame(study = factor(c("a", "a", "b", "b")), esid = factor(1:4)),
    prior_list = list(intercept = BayesTools::prior("normal", list(0, 1))),
    prior_random = BayesTools::prior_random(
      sd = BayesTools::prior("gamma", list(2, 2)),
      allocation = list(BayesTools::random_variance_allocation(
        name = "split", terms = c(study = "study", esid = "esid"),
        sd = BayesTools::prior("normal", list(0, 1), list(0, Inf)),
        inclusion = list(study = BayesTools::prior("spike", list(location = .5)))
      ))
    )
  )
  share <- seq(0.02, 0.98, length.out = n)
  samples <- cbind(
    mu_intercept = rep(c(-0.1, 0.1), length.out = n),
    mu__xRE_ALLOCx_split__allocation_sd = 0.3 + 0.4 * share,
    "mu__xRE_ALLOCx_split__weight[1]" = share,
    "mu__xRE_ALLOCx_split__weight[2]" = 1 - share,
    mu__xRE_ALLOCx_split__include_study_indicator = rep(c(0, 1), length.out = n)
  )
  fit <- coda::mcmc.list(coda::mcmc(samples))
  class(fit) <- c("BayesTools_fit", class(fit))
  attr(fit, "prior_list") <- result[["prior_list"]]
  attr(fit, "formula_design") <- list(mu = result[["formula_design"]])
  fit <- BayesTools:::.bt_attach_parameter_map(fit)
  fit <- BayesTools:::.bt_attach_draw_geometry(fit)
  fit <- BayesTools:::.bt_attach_fit_contract(fit)

  structure(
    list(fit = fit, data = structure(list(), random = TRUE)),
    class = c("RoBMA", "brma.mv", "brma")
  )
}

.gated_random_hypothesis <- function(object, parameter, statement,
                                     density_method = "KDE") {

  .hypothesis_brma_random(
    object                    = object,
    parameter                 = parameter,
    hypothesis                = BayesTools::hypothesis_rewrite(
      BayesTools::hypothesis_parse(statement),
      c(theta = parameter)
    ),
    standardized_coefficients = FALSE,
    conditional               = FALSE,
    logBF                     = FALSE,
    BF01                      = FALSE,
    seed                      = 1,
    density_method            = density_method,
    n_samples                 = 1000L,
    columns                   = "all"
  )
}


test_that("gated random point hypotheses use the exact BayesTools prior densities", {

  object <- .gated_random_object()
  # The study SD is T * gate * sqrt(w) with T ~ half-normal(0, 1),
  # P(gate = 1) = 0.5 and w ~ Beta(1, 1): its continuous part at y is
  # 0.5 * int_0^1 f_T(y / sqrt(s)) / sqrt(s) ds, and the variance density
  # at y^2 is the SD density at y divided by 2 y.
  sd_density <- function(y) {
    0.5 * stats::integrate(
      function(s) 2 * stats::dnorm(y / sqrt(s)) / sqrt(s),
      lower = 0, upper = 1, rel.tol = 1e-12
    )$value
  }
  cases <- list(
    list(parameter = "(mu) study: tau(intercept)", value = 0.3,
         prior = sd_density(0.3)),
    list(parameter = "(mu) study: tau2(intercept)", value = 0.09,
         prior = sd_density(0.3) / (2 * 0.3))
  )
  for (case in cases) {
    out <- .gated_random_hypothesis(
      object, case[["parameter"]], paste("theta =", case[["value"]])
    )
    expect_equal(as.numeric(out[["prior"]]), case[["prior"]], tolerance = 1e-8,
                 info = case[["parameter"]])

    # The route is the BayesTools mixed posterior of the catalog quantity.
    selection <- .brma_parameter_select_entry(
      object, case[["parameter"]], component = "random"
    )[["selection"]]
    samples <- list(BayesTools::parameter_mixed_posterior(object[["fit"]], selection))
    names(samples) <- "theta"
    class(samples) <- c("as_mixed_posteriors", "mixed_posteriors", "list")
    attr(samples, "prior_list") <- list(theta = BayesTools::prior_none())
    direct <- BayesTools::hypothesis_BF(
      posterior  = BayesTools::marginal_posterior(
        samples, "theta", prior_samples = TRUE, use_formula = FALSE,
        n_samples = 1000L
      ),
      hypothesis = paste("theta =", case[["value"]]),
      parameter  = "theta",
      seed       = 1
    )
    expect_equal(attr(out, "raw_BF"), attr(direct, "raw_BF"), tolerance = 1e-12,
                 info = case[["parameter"]])
  }
})


test_that("random point hypotheses refuse prior point masses and nonregular values", {

  object <- .gated_random_object()
  for (case in list(
    list(parameter = "(mu) split: tau2_prop(study)", value = 0),
    list(parameter = "(mu) split: tau2_prop(esid)", value = 1),
    list(parameter = "(mu) study: tau2(intercept)", value = 0)
  )) {
    expect_error(
      .gated_random_hypothesis(object, case[["parameter"]],
                               paste("theta =", case[["value"]])),
      paste0("because its prior has a point mass at ", case[["value"]], "."),
      fixed = TRUE,
      info = case[["parameter"]]
    )
  }
  # Interior proportions are regular points of the continuous part.
  expect_s3_class(
    .gated_random_hypothesis(object, "(mu) split: tau2_prop(study)", "theta = 0.3"),
    "BayesTools_hypothesis_BF"
  )
  # The ungated variance has an infinite prior ordinate at 0.
  expect_error(
    .gated_random_hypothesis(object, "(mu) esid: tau2(intercept)", "theta = 0"),
    class = "BayesTools_infinite_ordinate"
  )
})


test_that("hypothesis discovery lists point tests of gated random quantities", {

  object <- .gated_random_object()
  quantities <- hypothesis_quantities(object)
  random <- quantities[quantities[["component"]] == "random", , drop = FALSE]

  expect_gt(nrow(random), 0L)
  expect_true(all(random[["point_test"]]))
  expect_true(all(grepl("KDE", random[["point_test_methods"]], fixed = TRUE)))
  expect_true(all(random[["reason"]] == ""))
})


test_that("realized allocation gates define aggregate and proportion states", {

  allocation <- list(
    scale          = "total_variance",
    n_targets      = 2L,
    inclusion      = list(
      list(index = 1L, indicator_name = "gate_study"),
      list(index = 2L, indicator_name = "gate_drug")
    ),
    parent_factors = list()
  )
  raw_samples <- cbind(
    gate_study = c(0, 1, 0, 1),
    gate_drug  = c(0, 0, 1, 1)
  )

  total_metadata <- .brma_random_parameter_allocation_gate_metadata(list(
    spec = list(quantity = "sd_total", allocation_index = NA_integer_),
    allocation_definition = allocation
  ))
  total_state <- .brma_random_parameter_allocation_gate_state(
    total_metadata,
    raw_samples
  )
  expect_identical(total_state[["point_zero"]], c(TRUE, FALSE, FALSE, FALSE))
  expect_identical(total_state[["continuous"]], c(FALSE, TRUE, TRUE, TRUE))

  proportion_selected <- list(
    spec = list(quantity = "var_prop", allocation_index = 1L),
    allocation_definition = allocation
  )
  proportion_state <- .brma_random_parameter_allocation_gate_state(
    .brma_random_parameter_allocation_gate_metadata(proportion_selected),
    raw_samples
  )
  expect_identical(proportion_state[["defined"]], c(FALSE, TRUE, TRUE, TRUE))
  expect_identical(proportion_state[["point_zero"]], c(FALSE, FALSE, TRUE, FALSE))
  expect_identical(proportion_state[["point_one"]], c(FALSE, TRUE, FALSE, FALSE))
  expect_identical(proportion_state[["continuous"]], c(FALSE, FALSE, FALSE, TRUE))

  inherited_allocation <- allocation
  inherited_allocation[["inclusion"]] <- list()
  inherited_allocation[["parent_factors"]] <- list(
    list(inclusion_name = "gate_study")
  )
  inherited_selected <- proportion_selected
  inherited_selected[["allocation_definition"]] <- inherited_allocation
  inherited_state <- .brma_random_parameter_allocation_gate_state(
    .brma_random_parameter_allocation_gate_metadata(inherited_selected),
    raw_samples
  )
  expect_identical(inherited_state[["defined"]], c(FALSE, TRUE, FALSE, TRUE))
  expect_identical(
    inherited_state[["continuous"]],
    inherited_state[["defined"]]
  )
})
