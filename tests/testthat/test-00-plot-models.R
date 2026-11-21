.plot_models_replay <- function(contrast = "treatment", interaction = FALSE) {

  data <- data.frame(
    yi = seq(-.2, .2, length.out = 9), sei = .2,
    x = rep(1:3, each = 3), g = factor(rep(c("a", "b", "c"), 3))
  )
  setup <- review_compile_object(RoBMA(
    yi = yi, sei = sei, data = data, measure = "GEN",
    mods = if (interaction) ~ x * g else ~ x + g, scale = ~ x,
    set_contrast_factor_predictors = contrast,
    standardize_continuous_predictors = interaction,
    prior_heterogeneity_null = NULL,
    prior_effect_null = list(None = BayesTools::prior("point", list(0))),
    prior_bias = list(PET = BayesTools::prior_PET("normal", list(0, 1))),
    prior_unit_information_sd = 1, only_priors = TRUE
  ))
  priors <- c(setup[["compiled"]][["mu"]][["prior_list"]],
    setup[["compiled"]][["log_tau"]][["prior_list"]],
    .create_fit_priors(setup[["object"]][["data"]], setup[["object"]][["priors"]]))
  indicators <- c(
    mu_intercept = "221222", mu_x = "222222", mu_g = "222211",
    log_tau_intercept = "111111", log_tau_x = "222222", bias = "121211"
  )
  if (interaction) {
    indicators <- c(indicators, stats::setNames("222211", review_coefficient_name(setup, "x:g")))
  }
  columns <- unique(c(review_prior_columns(priors), paste0(names(indicators), "_indicator"), "PET"))
  samples <- matrix(0, 6L, length(columns), dimnames = list(NULL, columns))
  for (parameter in names(indicators)) {
    samples[, paste0(parameter, "_indicator")] <-
      as.integer(strsplit(indicators[[parameter]], "", fixed = TRUE)[[1L]])
  }
  samples[, "mu_intercept"] <- c(-1, .4, 0, .6, .4, .6)
  samples[, "mu_x"] <- 1
  samples[, "log_tau_intercept"] <- .2
  samples[, "log_tau_x"] <- -.1
  samples[, "PET"] <- .1
  factor_columns <- review_prior_columns(priors["mu_g"])
  samples[, factor_columns] <- cbind(c(1, 1, 3, -1, 0, 0), c(3, -1, 7, 1, 0, 0))
  if (interaction) {
    interaction_columns <- review_prior_columns(priors[review_coefficient_name(setup, "x:g")])
    samples[, interaction_columns] <- cbind(c(5, 5, 5, 5, 0, 0), c(5, 5, 5, 5, 0, 0))
  }
  list(fit = review_attach_replay(setup, list(samples[1:3, ], samples[4:6, ])),
       samples = samples, priors = priors)
}


test_that("posterior model plot retains the full individual-model weights", {

  fit <- .plot_models_replay()[["fit"]]
  plot_data <- .plot_models_data(fit)
  expected <- summary_models(fit, type = "individual", include_mcmc_diagnostics = FALSE)
  posterior <- expected[["individual"]][["post_prob"]]
  expect_equal(unname(plot_data[["probabilities"]]), posterior[posterior > 0])
  expect_identical(unique(plot_data[["cells"]][["label"]]), c(
    "Effect", "Heterogeneity", "Publication Bias", "(mu) x",
    "(mu) g[b]", "(mu) g[c]", "(tau) x"
  ))
  expect_error(plot_models(fit, "x"), "ambiguous")
  subset <- .plot_models_data(fit, c("(log_tau) x", "g[c]", "mu"),
    labels = c("mu_g[c]" = "Delivery C"))
  expect_identical(subset[["probabilities"]], plot_data[["probabilities"]])
  expect_identical(unique(subset[["cells"]][["label"]]), c("(tau) x", "Delivery C", "Effect"))
  expect_error(.plot_models_data(fit, "g[a]"), "unavailable")
  expect_identical(unique(.plot_models_data(fit, "x", component = "location")[["cells"]][["row"]]),
                   "mu_x")
  expect_setequal(plot_data[["cells"]][["state"]][plot_data[["cells"]][["row"]] == "mu"],
                  c("positive", "negative", "excluded"))
  expect_setequal(plot_data[["cells"]][["state"]][plot_data[["cells"]][["row"]] == "mu_g[b]"],
                  c("positive", "included", "excluded"))
})


test_that("centered factor levels apply the fitted contrast to draws before taking medians", {

  for (contrast in c("meandif", "orthonormal")) {
    replay <- .plot_models_replay(contrast)
    fit <- replay[["fit"]]
    data <- .plot_models_data(fit, c("g", "mu"),
      factor_view = c(g = "levels"), group_labels = TRUE)
    expect_identical(unique(data[["cells"]][["label"]]),
                     c("g[dif: a]", "g[dif: b]", "g[dif: c]", "Effect"))
    components <- .summary_models_components(fit)
    joint <- .plot_models_joint(components, replay[["samples"]])
    prior <- replay[["priors"]][["mu_g"]]
    design <- attr(prior, "factor_design", exact = TRUE)
    if (is.list(design)) design <- design[["design"]]
    estimates <- replay[["samples"]][, review_prior_columns(replay[["priors"]]["mu_g"])] %*%
      t(design)
    states <- vapply(joint[["draws"]], function(draws) {
      if (replay[["samples"]][draws[[1L]], "mu_g_indicator"] == 1L) {
        return(rep("excluded", 3L))
      }
      signed <- sign(apply(estimates[draws, , drop = FALSE], 2L, stats::median))
      c("negative", "included", "positive")[signed + 2L]
    }, character(3L))
    cells <- data[["cells"]]
    expect_equal(unname(as.vector(t(states))), cells[["state"]][grepl("^mu_g", cells[["row"]])])
    terms <- .plot_models_data(fit, "g", factor_view = c(g = "terms"))
    expect_identical(unique(terms[["cells"]][["row"]]), "mu_g")
    expect_setequal(terms[["cells"]][["state"]], c("included", "excluded"))
  }
})


test_that("mixture rows use declared members and publication-bias colors are scoped", {

  fit <- .plot_models_replay()[["fit"]]
  blocks <- .plot_models_data(fit, c("mu", "bias"),
    mixture_view = c(bias = "individual"),
    col_bias = c(no_bias = "gray80", "bias[PET]" = "darkgreen"))
  expect_identical(unique(blocks[["cells"]][["row"]]),
                   c("mu", "bias{None}", "bias{PET}"))
  bias <- blocks[["cells"]][blocks[["cells"]][["group"]] == "Components" &
    blocks[["cells"]][["row"]] != "mu", ]
  expect_true(all(table(bias[["model"]][bias[["state"]] != "excluded"]) == 1L))
  expect_identical(unname(blocks[["colors"]]), c("gray80", "darkgreen"))
  categories <- .plot_models_data(fit, c("mu", "bias"), mixture_view = "model",
    col_bias = c(None = "gray80", "bias[PET]" = "darkgreen"))
  expect_identical(unique(categories[["cells"]][["row"]]), c("mu", "bias"))
  expect_false(any(categories[["colors"]][startsWith(names(categories[["colors"]]), "mu ")] %in%
                     c("gray80", "darkgreen")))
  expect_s3_class(plot_models(fit, "bias", legend = TRUE, plot_type = "ggplot"), "ggplot")
})


test_that("directional colors retain the fitted predictor scale with centered interactions", {

  replay <- .plot_models_replay(interaction = TRUE)
  fit <- replay[["fit"]]
  data <- .plot_models_data(fit, c("g[b]", "mu"))
  group_b <- data[["cells"]][data[["cells"]][["row"]] == "mu_g[b]", ]
  # Uncentering would subtract the same positive interaction contribution
  # from both treatment draws and a positive slope contribution from both
  # intercept draws. The plot uses the fitted gate and its parameterization.
  joint <- .plot_models_joint(.summary_models_components(fit), replay[["samples"]])
  selected <- which(vapply(joint[["draws"]], function(rows) 2L %in% rows, logical(1)))
  x <- rep(1:3, each = 3)
  offset <- mean(x) / stats::sd(x)
  expect_lt(stats::median(c(1, -1) - 5 * offset), 0)
  expect_identical(group_b[["state"]][[selected]], "included")
  effect <- data[["cells"]][data[["cells"]][["row"]] == "mu", ]
  expect_lt(stats::median(c(.4, .6) - offset), 0)
  expect_identical(effect[["state"]][[selected]], "positive")
  excluded <- which(vapply(joint[["draws"]], function(rows) 3L %in% rows, logical(1)))
  expect_identical(effect[["state"]][[excluded]], "excluded")
})


test_that("declared, always-included coefficients can accompany model components", {

  data <- data.frame(yi = c(.1, .2, .15, .25), sei = .2, x = c(-1, 1, -1, 1))
  setup <- review_compile_object(BMA(yi = yi, sei = sei, mods = ~ x, data = data,
    measure = "GEN", prior_mods = list(intercept = BayesTools::prior("normal", list(.4, .2))),
    prior_heterogeneity_null = NULL, prior_unit_information_sd = 1,
    standardize_continuous_predictors = FALSE, only_priors = TRUE))
  samples <- cbind(mu_intercept = c(-.4, .4), mu_x = c(0, .3),
    mu_x_indicator = c(1L, 2L), tau = .2, tau_indicator = 1L)
  fit <- review_attach_replay(setup, samples)
  data <- .plot_models_data(fit, "mu", mixture_view = "individual")
  expect_identical(data[["cells"]][["row"]], c("mu", "mu"))
  expect_identical(data[["cells"]][["state"]], c("negative", "positive"))
  expect_identical(unname(data[["probabilities"]]), c(.5, .5))
})
