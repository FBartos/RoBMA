probe_state("test-00-bias-exempt-input.R")
context("Row-level publication bias exemptions")

test_that("native selection bins require a logical exemption mask", {

  context <- list(obs_bin = 1:3, bias_exempt = c(TRUE, FALSE, TRUE))
  expect_identical(.selection_kernel_obs_bin(context), c(0L, 2L, 0L))
  expect_error(.selection_kernel_obs_bin(context, c(1L, 0L, 1L)),
               "'bias_exempt' must contain one logical TRUE or FALSE per observed bin.", fixed = TRUE)
})

test_that("bias_exempt is a strict logical input aligned to retained outcome rows", {

  dat <- data.frame(
    yi = c(.2, -.1, NA, .3, .4), sei = rep(.1, 5),
    registered = c(NA, TRUE, NA, FALSE, TRUE)
  )
  for (constructor in list(bPET, bPEESE, bselmodel, RoBMA,
                           bPET.mv, bPEESE.mv, bselmodel.mv, RoBMA.mv)) {
    # It follows cluster (random for .mv) and, like other data inputs, has no default.
    arguments <- formals(constructor)
    anchor <- if ("cluster" %in% names(arguments)) "cluster" else "random"
    expect_identical(match("bias_exempt", names(arguments)), match(anchor, names(arguments)) + 1L)
    expect_identical(deparse(arguments[["bias_exempt"]]), "")
    expect_warning(
      fit <- constructor(yi = yi, sei = sei, data = dat, subset = 2:5,
                         bias_exempt = registered, only_data = TRUE),
      "1 observation\\(s\\) removed due to missing values"
    )
    expect_identical(fit[["data"]][["outcome"]][["bias_exempt"]], c(TRUE, FALSE, TRUE))
    expect_equal(fit[["data"]][["outcome"]][["yi"]], c(-.1, .3, .4))
  }
  dat <- dat[c(2, 4, 5), ]
  fitted_data <- function(...) {
    bPET(yi = yi, sei = sei, data = dat, only_data = TRUE, ...)[["data"]]
  }
  inverted <- bPET(yi = yi, sei = sei, data = dat, bias_exempt = !registered,
                   only_data = TRUE)
  expect_identical(inverted[["data"]][["outcome"]][["bias_exempt"]], c(FALSE, TRUE, FALSE))
  expect_identical(fitted_data(bias_exempt = rep(FALSE, 3)), fitted_data())
  expect_error(fitted_data(bias_exempt = TRUE),
               "The 'bias_exempt' argument must have length '3'.", fixed = TRUE)
  expect_error(fitted_data(bias_exempt = c(TRUE, NA, FALSE)),
               "The 'bias_exempt' argument must be TRUE or FALSE on every retained row.",
               fixed = TRUE)
  for (invalid in list(1:3, c("yes", "no", "yes"), factor(1:3), matrix(TRUE, 3, 1))) {
    expect_error(fitted_data(bias_exempt = invalid),
                 "The 'bias_exempt' argument must be a logical vector.", fixed = TRUE)
  }
})

test_that("correlated PET and PEESE bridge densities retain exempt outcomes", {

  yi <- c(-.24, .18, .30)
  V <- matrix(c(.16, .04, -.02, .04, .25, .03, -.02, .03, .36), 3, 3)
  exempt <- c(TRUE, FALSE, TRUE)
  for (parameter in c("PET", "PEESE")) {
    constructor <- if (parameter == "PET") bPET.mv else bPEESE.mv
    for (direction in c("positive", "negative")) {
      object <- constructor(yi = yi, V = V, bias_exempt = exempt,
                            measure = "GEN", prior_unit_information_sd = 1,
                            effect_direction = direction, only_priors = TRUE)
      data <- object[["data"]]
      fit_data <- .create_fit_data(data, object[["priors"]])
      expect_equal(fit_data[["sei"]], sqrt(diag(V)))
      expect_identical(fit_data[["bias_exempt"]], as.integer(exempt))
      fit_data <- .marglik_add_random_covariance_bridge_data(
        fit_data, data, marginalizing = TRUE, sampling_latent_marginalized = TRUE
      )
      parameters <- list(mu = -.05, tau = 0)
      parameters[[parameter]] <- .6
      actual <- .log_posterior(
        parameters = parameters, data = fit_data, is_scale = FALSE,
        is_multilevel = FALSE, is_weights = FALSE, is_known_v = TRUE,
        is_PET = parameter == "PET", is_PEESE = parameter == "PEESE",
        is_weightfunction = FALSE, effect_direction = direction, outcome_type = "norm",
        model_data = data, sampling_latent_marginalized = TRUE
      )
      slope <- if (direction == "positive") .6 else -.6
      predictor <- if (parameter == "PET") sqrt(diag(V)) else diag(V)
      mean <- -.05 + slope * predictor * !exempt
      expect_equal(actual, mvtnorm::dmvnorm(yi, mean, V, log = TRUE), tolerance = 1e-12)
    }
  }
})

test_that("ordered selection ordinates retain only eligible observed weight factors", {

  yi <- c(.1, .11, .6, .12)
  sei <- rep(.2, 4)
  weights <- c(3, 2, 4, .5)
  exempt <- c(TRUE, FALSE, TRUE, FALSE)
  prior <- BayesTools::prior_weightfunction("one-sided", .05, BayesTools::wf_fixed(c(1, .25)))
  object <- bselmodel(yi = yi, sei = sei, weights = weights, bias_exempt = exempt,
                      prior_bias = prior, measure = "SMD", effect_direction = "positive",
                      only_priors = TRUE)
  data <- object[["data"]]
  priors <- object[["priors"]]
  draws <- matrix(c(.1, .01), 1L, dimnames = list(NULL, c("mu", "tau")))
  selection <- .selection_context_from_parts(object[["fit"]], data, priors, draws, "positive")
  state <- list(active_setup = list(priors = priors), ordered_chart = list(source = "total"))
  setup <- list(yi = yi, sei = sei, selection_sei = sei, weights = weights,
                mu = matrix(.1, 1L, 4L), tau_within = matrix(.01, 1L, 4L))
  expect_equal(.iwmde_ordered_selection_observed_constant(list(data = data), state, setup, selection),
               (2 + .5) * log(.25))
  selection[["bias_exempt"]] <- rev(exempt)
  expect_error(.iwmde_ordered_selection_observed_constant(list(data = data), state, setup, selection),
               class = "RoBMA_ordered_selection_route")
})

test_that("prediction exemption flags are logical newdata columns, not fitted predictors", {

  dat <- data.frame(yi = c(.1, -.2, .4), sei = c(.5, .7, .9), bias_exempt = c(-2, 0, 3),
                    registered = c(FALSE, TRUE, FALSE))
  common <- list(sei = quote(sei), data = dat, measure = "GEN",
                 prior_unit_information_sd = 1, only_priors = TRUE, silent = TRUE)
  # A fitted moderator named 'bias_exempt' keeps its values and flags nothing.
  newdata <- data.frame(sei = c(.4, .6), bias_exempt = c(1.5, -1))
  plain <- do.call(bPET, c(list(yi = yi ~ bias_exempt), common))
  prepared <- .prepare_newdata(plain, newdata, type = "terms")
  expect_identical(prepared[["mods"]][["bias_exempt"]], newdata$bias_exempt)
  expect_null(prepared[["outcome"]][["bias_exempt"]])
  # With exempt fitted rows, that column cannot also supply new rows' flags.
  masked <- do.call(bPET, c(list(yi = yi ~ bias_exempt, bias_exempt = quote(registered)), common))
  expect_error(.prepare_newdata(masked, newdata, type = "terms"),
    paste0("Publication-bias predictions from this fit need exemption flags for new rows, ",
           "but 'bias_exempt' names a variable of the fitted model, so 'newdata' cannot ",
           "supply them."), fixed = TRUE)
  expect_identical(.prepare_newdata(masked, newdata, type = "terms", bias_adjusted = TRUE)[["mods"]][["bias_exempt"]],
                   newdata$bias_exempt)
  # Callers that construct rows, such as regplot() grids, pass their flags.
  grid <- .prepare_newdata(masked, newdata, type = "terms", bias_exempt = c(FALSE, TRUE))
  expect_identical(grid[["mods"]][["bias_exempt"]], newdata$bias_exempt)
  expect_identical(grid[["outcome"]][["bias_exempt"]], c(FALSE, TRUE))
  # Flags enter the parser by value and never stand in for a missing predictor.
  hidden <- dat
  names(hidden)[names(hidden) == "bias_exempt"] <- ".RoBMA_newdata_parser_bias_exempt"
  hidden_fit <- do.call(bPET, c(list(yi = yi ~ .RoBMA_newdata_parser_bias_exempt,
    bias_exempt = quote(registered)), utils::modifyList(common, list(data = hidden))))
  expect_error(.prepare_newdata(hidden_fit, data.frame(sei = .4, bias_exempt = TRUE), type = "terms"),
    ".RoBMA_newdata_parser_bias_exempt", fixed = TRUE)
  # A best-rule publication-group column of that name is not a flag either.
  grouped <- data.frame(yi = dat$yi, bias_exempt = c("a", "a", "b"), registered = dat$registered)
  group_model <- list(yi = quote(yi), V = diag(dat$sei^2), data = grouped,
    measure = "GEN", prior_unit_information_sd = 1, only_priors = TRUE, silent = TRUE,
    prior_bias = BayesTools::prior_weightfunction("one-sided", .5, BayesTools::wf_fixed(c(1, .2)),
      model = selection_model(weight_rule = "best", group = "bias_exempt")))
  flags <- data.frame(sei = .4, bias_exempt = TRUE)
  expect_null(.prepare_newdata_bias_exempt(do.call(bselmodel.mv, group_model)$data, flags))
  group_model$bias_exempt <- quote(registered)
  expect_error(.prepare_newdata_bias_exempt(do.call(bselmodel.mv, group_model)$data, flags),
    "names a variable of the fitted model", fixed = TRUE)
  # A character flag is an error, never the name of another column.
  exempt <- do.call(bPET, c(list(yi = quote(yi), bias_exempt = quote(registered)), common))
  flags <- data.frame(sei = .4, bias_exempt = "registered", registered = TRUE)
  expect_error(.prepare_newdata(exempt, flags, type = "terms"),
    "The 'bias_exempt' column in 'newdata' must be logical.", fixed = TRUE)
  flags$bias_exempt <- TRUE
  expect_identical(.prepare_newdata(exempt, flags, type = "terms")[["outcome"]][["bias_exempt"]], TRUE)
})
