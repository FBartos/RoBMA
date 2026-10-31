context("qCMDE/IWMDE classify inclusion indicators by their declared role")

# A normal meta-regression on one continuous predictor 'predictor' with 400
# synthetic draws; the heterogeneity has the prior 'tau_prior'.
.review_predictor_replay <- function(predictor, tau_prior) {

  dat <- data.frame(yi = c(.1, .25, .3, .4, .2, .15, .05, .35), sei = .2)
  dat[[predictor]] <- c(-1, 0, 1, 2, 3, 1, -2, 0)
  object <- brma(
    yi = dat[["yi"]], sei = dat[["sei"]],
    mods = stats::as.formula(paste("~", predictor)), data = dat, measure = "GEN",
    prior_unit_information_sd = 1,
    prior_effect = BayesTools::prior("normal", list(0, 1)),
    prior_mods = BayesTools::prior("normal", list(0, 1)),
    prior_heterogeneity = tau_prior, only_priors = TRUE, silent = TRUE
  )
  setup <- review_compile_object(object)
  withr::local_seed(11)
  n <- 400L
  samples <- cbind(
    stats::rnorm(n, .1, .1), stats::rnorm(n, .05, .1),
    abs(stats::rnorm(n, .15, .08)) + .02
  )
  colnames(samples) <- c(
    review_coefficient_name(setup, "intercept"),
    review_coefficient_name(setup, predictor), "tau"
  )

  review_attach_replay(setup, samples)
}

.review_indicator_bf <- function(object, statement) {

  hypothesis(
    object, statement, columns = "all", density_method = "qCMDE",
    density_control = list(n_points = 20, samples = 100), seed = 1
  )
}

test_that("qCMDE hypotheses run for a continuous predictor named like an indicator", {

  tau_prior <- BayesTools::prior("normal", list(0, 1), truncation = list(0, Inf))
  object  <- .review_predictor_replay("treatment_indicator", tau_prior)
  control <- .review_predictor_replay("treatment_flag", tau_prior)

  # point, linear and non-linear targets run like those of the renamed predictor
  for (statement in c("treatment_indicator = 0", "intercept = 0", "tau = 0.15")) {
    renamed  <- sub("treatment_indicator", "treatment_flag", statement)
    expected <- .review_indicator_bf(control, renamed)
    observed <- tryCatch(.review_indicator_bf(object, statement), error = identity)
    expect_s3_class(observed, "data.frame")
    if (is.data.frame(observed)) {
      expect_true(is.finite(observed[["BF"]][[1L]]), info = statement)
      expect_equal(observed[["BF"]], expected[["BF"]], tolerance = 1e-8, info = statement)
    }
  }
})

test_that("a declared coefficient named like an indicator is no indicator", {

  tau_prior <- BayesTools::prior("normal", list(0, 1), truncation = list(0, Inf))
  context <- .iwmde_context(.review_predictor_replay("treatment_indicator", tau_prior))

  expect_identical(context[["indicator_names"]], character())
  expect_false(.iwmde_parameter_is_indicator(
    "mu_treatment_indicator", context[["indicator_names"]]
  ))
  expect_identical(
    .iwmde_parameter_spec(context, "mu_treatment_indicator")[["status"]], "ok"
  )
})

test_that("inclusion indicators are those of the declared mixture priors", {

  dat <- data.frame(
    yi = c(.1, .2, .3, .4, .2, .1), sei = .2,
    x = c(-1, 0, 1, 2, 0, 1), treatment_indicator = c(-1, 0, 1, 2, 3, 1)
  )
  object <- RoBMA(
    yi = yi, sei = sei, mods = ~ x + treatment_indicator, data = dat,
    measure = "GEN", prior_unit_information_sd = 1, only_priors = TRUE, silent = TRUE
  )
  setup <- review_compile_object(object)
  priors <- list()
  for (compiled in setup[["compiled"]]) priors <- c(priors, compiled[["prior_list"]])
  priors <- c(priors, .create_fit_priors(object[["data"]], object[["priors"]]))
  indicators <- c(
    "mu_x_indicator", "mu_treatment_indicator_indicator", "tau_indicator",
    "bias_indicator"
  )
  columns <- c(review_prior_columns(priors), indicators)
  samples <- matrix(.5, 20L, length(columns), dimnames = list(NULL, columns))
  samples[, indicators] <- rep(c(1, 2), 10)
  context <- .iwmde_context(review_attach_replay(setup, samples))

  # the coefficient 'mu_treatment_indicator' ends like an indicator, its
  # inclusion indicator 'mu_treatment_indicator_indicator' is one
  expect_setequal(context[["indicator_names"]], indicators)
  for (name in indicators) {
    expect_true(.iwmde_parameter_is_indicator(name, context[["indicator_names"]]),
                info = name)
  }
  expect_false(.iwmde_parameter_is_indicator(
    "mu_treatment_indicator", context[["indicator_names"]]
  ))
})

test_that("indicator names follow the declared priors and coordinate roles", {

  mixture <- BayesTools::prior_mixture(
    list(BayesTools::prior("spike", list(0)), BayesTools::prior("normal", list(0, 1))),
    is_null = c(TRUE, FALSE)
  )
  simple <- BayesTools::prior("normal", list(0, 1))
  samples <- cbind(
    theta_indicator = c(1, 2), gate_indicator = c(1, 2), sd_indicator = c(1, 2),
    other_indicator = c(1, 2), tau = c(1, 1)
  )

  # only declared mixtures have an indicator; an undeclared column that ends
  # like an indicator and a context without priors have none
  expect_identical(
    .iwmde_indicator_names(samples, list(theta = mixture, tau = simple)),
    "theta_indicator"
  )
  expect_identical(.iwmde_indicator_names(samples), character())
  # inclusion gates of random-effect variance allocations are no model indicators
  gated <- structure(simple, random_allocation_indicator = "gate_indicator")
  expect_identical(
    .iwmde_indicator_names(samples, list(theta = mixture, gate = mixture, g = gated)),
    "theta_indicator"
  )
  # declared coordinates: random-effect inclusion indicators count and a
  # formula coefficient never does
  coordinates <- data.frame(
    coordinate_name = c("sd_indicator", "theta_indicator"),
    role = c("random_inclusion_indicator", "fixed_coefficient"),
    stringsAsFactors = FALSE
  )
  expect_identical(
    .iwmde_indicator_names(samples, list(theta = mixture, tau = simple), coordinates),
    "sd_indicator"
  )
  expect_null(.iwmde_parameter_coordinates(list()))
})

test_that("inclusion gates of variance allocations stay discrete targets", {

  # the gate is no model indicator, but its states are discrete all the same
  context <- .iwmde_context(shared_gate_random_object(root_gate = c(0, 1, 0, 1)))
  gate <- grep("include_component_indicator$", colnames(context[["posterior_samples"]]),
               value = TRUE)

  expect_length(gate, 1L)
  expect_false(gate %in% context[["indicator_names"]])
  expect_identical(.iwmde_parameter_spec(context, gate)[["status"]], "unsupported")
})
