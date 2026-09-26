# Fixture: gated_random_object() (helper-contracts.R), a gated total-variance
# allocation of a half-normal SD over 'study' (inclusion gate with prior
# probability 0.5) and 'esid' (ungated), Dirichlet(1, 1) shares, and
# synthetic draws.
.gated_random_hypothesis <- function(object, parameter, value,
                                     density_method = "KDE") {

  hypothesis(
    object,
    paste0("`", parameter, "` = ", value),
    density_method = density_method,
    seed           = 1,
    n_samples      = 1000L,
    columns        = "all"
  )
}


test_that("gated random point hypotheses use the exact BayesTools prior densities", {

  object <- gated_random_object()
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
  # The route is the BayesTools mixed posterior of the catalog quantity of
  # the standard deviation (the variance is tested through it).
  selection <- .brma_parameter_select_entry(
    object, "(mu) study: tau(intercept)", component = "random"
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
    hypothesis = "theta = 0.3",
    parameter  = "theta",
    seed       = 1
  )
  for (case in cases) {
    out <- .gated_random_hypothesis(object, case[["parameter"]], case[["value"]])
    expect_equal(as.numeric(out[["prior"]]), case[["prior"]], tolerance = 1e-8,
                 info = case[["parameter"]])
    expect_equal(attr(out, "raw_BF"), attr(direct, "raw_BF"), tolerance = 1e-12,
                 info = case[["parameter"]])
  }
})


test_that("random point hypotheses refuse prior point masses and nonregular values", {

  object <- gated_random_object()
  for (case in list(
    list(parameter = "(mu) split: tau2_prop(study)", value = 0),
    list(parameter = "(mu) split: tau2_prop(esid)", value = 1),
    list(parameter = "(mu) study: tau2(intercept)", value = 0)
  )) {
    expect_error(
      .gated_random_hypothesis(object, case[["parameter"]], case[["value"]]),
      paste0("because its prior has a point mass at ", case[["value"]], "."),
      fixed = TRUE,
      class = "BayesTools_point_mass_at_null",
      info  = case[["parameter"]]
    )
  }
  # Interior proportions are regular points of the continuous part.
  expect_s3_class(
    .gated_random_hypothesis(object, "(mu) split: tau2_prop(study)", 0.3),
    "BayesTools_hypothesis_BF"
  )
  # The ungated variance has an infinite prior ordinate at 0.
  expect_error(
    .gated_random_hypothesis(object, "(mu) esid: tau2(intercept)", 0),
    class = "BayesTools_infinite_ordinate"
  )
})


test_that("hypothesis discovery lists point tests of gated random quantities", {

  object <- gated_random_object()
  quantities <- hypothesis_quantities(object)
  random <- quantities[quantities[["component"]] == "random", , drop = FALSE]

  # Point tests are available at values off the gate atoms; qCMDE/IWMDE are
  # unavailable for random-formula models without known V.
  capability <- .iwmde_capability(object = object, density_method = "qCMDE")
  expect_gt(nrow(random), 0L)
  expect_true(all(random[["point_test"]]))
  expect_true(all(random[["point_test_methods"]] == "KDE"))
  expect_true(all(grepl(capability[["reason"]], random[["reason"]], fixed = TRUE) |
                    grepl("no supported scalar random-component coordinate",
                          random[["reason"]], fixed = TRUE)))

  # A quantity without a scalar random-component coordinate is refused with
  # the hypothesis classes followed by the density-method classes of the
  # same refusal by plot().
  for (method in c("qCMDE", "IWMDE")) {
    error <- tryCatch(
      .gated_random_hypothesis(object, "(mu) split: tau2_prop(study)", 0.3,
                               density_method = method),
      error = identity
    )
    expect_identical(
      class(error),
      c("RoBMA_hypothesis_target", "RoBMA_hypothesis_unavailable",
        "RoBMA_density_method_random_target",
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


test_that("realized allocation gates define aggregate and proportion states", {

  # The study share is gated (gate 0, 1, 0, 1); the esid share is not.
  object  <- gated_random_object(n = 4L)
  samples <- as.matrix(object[["fit"]][[1L]])
  gate    <- samples[, "mu__xRE_ALLOCx_split__include_study_indicator"] == 1
  context <- list(object = object, posterior_samples = samples)
  states  <- function(parameter) {
    selected  <- .brma_random_parameter_select(object, parameter)
    selection <- .brma_random_parameter_gate_selection(object, selected)
    expect_false(is.null(selection), info = parameter)
    .iwmde_gate_states(context, list(gate_selection = selection))
  }

  # The ungated esid share keeps the total positive in every draw.
  total <- states("(mu) split: tau_total")
  expect_identical(total[["defined"]], rep(TRUE, 4L))
  expect_identical(total[["continuous"]], rep(TRUE, 4L))
  expect_identical(total[["point_zero"]], rep(FALSE, 4L))

  study <- states("(mu) split: tau2_prop(study)")
  expect_identical(study[["defined"]], rep(TRUE, 4L))
  expect_identical(study[["point_zero"]], !gate)
  expect_identical(study[["point_one"]], rep(FALSE, 4L))
  expect_identical(study[["continuous"]], gate)

  esid <- states("(mu) split: tau2_prop(esid)")
  expect_identical(esid[["point_zero"]], rep(FALSE, 4L))
  expect_identical(esid[["point_one"]], !gate)
  expect_identical(esid[["continuous"]], gate)

  # The allocation chain of the esid SD has no gate: it is continuous in
  # every draw.
  esid_sd <- states("(mu) esid: tau(intercept)")
  expect_identical(esid_sd[["continuous"]], rep(TRUE, 4L))
  expect_identical(esid_sd[["point_zero"]], rep(FALSE, 4L))
})
