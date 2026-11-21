.plot_models_replay <- function(contrast = "treatment", interaction = FALSE,
                                bias_names = "PET", null_bias_name = "None") {

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
    prior_bias = stats::setNames(lapply(bias_names, function(name) {
      BayesTools::prior_PET("normal", list(0, 1))
    }), bias_names),
    prior_bias_null = stats::setNames(list(BayesTools::prior_none()), null_bias_name),
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
  expect_error(plot_models(fit, c("bias", "")), "'parameter' selectors must be non-empty")
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
})


test_that("plot_models methods dispatch to both renderers with the displayed row count", {

  expect_error(plot_models(NULL), "requires a fitted model-averaged RoBMA object")
  fit <- .plot_models_replay()[["fit"]]
  expect_error(plot_models(fit, unused_setting = TRUE), "Unused argument")
  options <- list(x = fit, parameter = c("mu", "g", "bias"), group_labels = TRUE,
                  mixture_view = c(bias = "individual"), legend = TRUE)
  ggplot <- do.call(plot_models, c(options, list(plot_type = "ggplot")))
  expect_s3_class(ggplot, "ggplot")
  expect_identical(attr(ggplot, "n_rows"), 5L)

  withr::local_pdf(NULL)
  base <- withVisible(do.call(plot_models, options))
  expect_identical(base, list(value = NULL, visible = FALSE))
})


test_that("bias palette keys identify one purpose even when member names collide", {

  fit <- .plot_models_replay(bias_names = c("no_bias", "PET", "bias[PET]"),
                            null_bias_name = "bias")[["fit"]]
  binary <- c(bias = "darkgreen", no_bias = "gray80")
  for (view in c("model", "individual")) {
    data <- .plot_models_data(fit, "bias", mixture_view = view, col_bias = binary)
    colors <- stats::setNames(data[["colors"]], data[["state_labels"]])
    expect_identical(colors[["Publication Bias: bias"]], "gray80")
    expect_identical(colors[["Publication Bias: no_bias"]], "darkgreen")

    data <- .plot_models_data(fit, "bias", mixture_view = view,
      col_bias = c(binary, "bias[bias]" = "black", "bias[no_bias]" = "skyblue",
                   "bias[PET]" = "tomato"))
    colors <- stats::setNames(data[["colors"]], data[["state_labels"]])
    expect_identical(colors[["Publication Bias: bias"]], "black")
    expect_identical(colors[["Publication Bias: no_bias"]], "skyblue")
    expect_identical(colors[["Publication Bias: PET"]], "tomato")
    expect_identical(colors[["Publication Bias: bias[PET]"]], "darkgreen")
  }
})


test_that("color keys validate against plotted components and mixture legends track labels", {

  fit <- .plot_models_replay()[["fit"]]
  expect_error(plot_models(fit, col = c(postive = "green")), "Unknown 'col' key: 'postive'")
  expect_error(plot_models(fit, col = c(bias = "green")), "Unknown 'col' key: 'bias'")
  expect_error(plot_models(fit, "bias", col_bias = c(PETT = "green")),
               "Unknown 'col_bias' key: 'PETT'")
  expect_error(plot_models(fit, "mu", col_bias = c(PET = "green")),
               "'col_bias' requires a displayed publication-bias row")
  category <- .plot_models_data(fit, "bias", labels = c(bias = "Bias mechanism"),
                                mixture_view = "model")
  expect_identical(unname(category[["state_labels"]]),
                   c("Bias mechanism: None", "Bias mechanism: PET"))
  member <- .plot_models_data(fit, "bias", labels = c("bias[PET]" = "PET fitted"),
                              mixture_view = "individual")
  expect_identical(unname(member[["state_labels"]]),
                   c("Publication Bias: None", "PET fitted"))
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


test_that("heterogeneity allocation remains a full joint dimension when heterogeneity is off", {

  data <- data.frame(yi = seq(-.3, .3, length.out = 9), sei = .15, cl = rep(1:3, each = 3))
  setup <- review_compile_object(BMA(yi = yi, sei = sei, cluster = cl, data = data,
    measure = "GEN", prior_unit_information_sd = 1, only_priors = TRUE,
    prior_heterogeneity_allocation = list(Variable = BayesTools::prior("beta", list(1, 1))),
    prior_heterogeneity_allocation_null = list(Equal = BayesTools::prior("point", list(.5)))))
  joint <- rbind(c(1, 2, 1), c(1, 2, 2), c(2, 1, 1), c(2, 1, 2), c(2, 2, 1), c(2, 2, 2))
  counts <- c(1L, 2L, 3L, 4L, 5L, 6L)
  observed <- joint[rep(seq_len(nrow(joint)), counts), , drop = FALSE]
  colnames(observed) <- paste0(c("mu", "tau", "rho"), "_indicator")
  samples <- cbind(observed, mu = .3 * (observed[, 1L] == 2),
    tau = .3 * (observed[, 2L] == 2), rho = ifelse(observed[, 3L] == 2, .4, .5),
    `gamma[1]` = 0, `gamma[2]` = 0, `gamma[3]` = 0)
  fit <- review_attach_replay(setup, samples)
  data <- .plot_models_data(fit)
  chosen <- .plot_models_data(fit, "heterogeneity_allocation", component = "scale")
  expect_identical(unique(data[["cells"]][["row"]]), c("mu", "tau", "rho"))
  expect_identical(chosen[["probabilities"]], data[["probabilities"]])
  expect_identical(unique(chosen[["cells"]][["label"]]), "Heterogeneity Allocation")
  expect_error(.plot_models_data(fit, "rho", component = "mods"), "unavailable")

  expected <- stats::setNames(counts / sum(counts), apply(joint, 1L, paste, collapse = ":"))
  expect_equal(data[["probabilities"]][names(expected)], expected)
  model_summary <- summary_models(fit, type = "individual", include_mcmc_diagnostics = FALSE)
  weights <- model_summary[["individual"]][["post_prob"]]
  expect_equal(unname(data[["probabilities"]]), weights[weights > 0])
  expect_identical(chosen[["cells"]][["state"]][
    match(c("2:1:1", "2:1:2"), chosen[["cells"]][["model"]])], c("excluded", "included"))

  member <- .plot_models_data(fit, "rho", mixture_view = "individual")
  expect_identical(unique(member[["cells"]][["row"]]), c("rho{Equal}", "rho{Variable}"))
  expect_true(all(table(member[["cells"]][["model"]][
    member[["cells"]][["state"]] != "excluded"]) == 1L))
})


test_that("the effect row follows the presence of the fitted intercept", {

  data <- data.frame(yi = seq(-.2, .2, length.out = 6), sei = .2,
    g = factor(rep(c("a", "b"), 3)))
  for (intercept in c(FALSE, TRUE)) {
    setup <- review_compile_object(BMA(yi = yi, sei = sei, data = data, measure = "GEN",
      mods = if (intercept) ~ g else ~ 0 + g,
      prior_mods = if (intercept) list(intercept = BayesTools::prior("point", list(0))) else NULL,
      set_contrast_factor_predictors = if (intercept) "treatment" else "independent",
      prior_unit_information_sd = 1, only_priors = TRUE))
    prior <- setup[["compiled"]][["mu"]][["prior_list"]][["mu_g"]]
    columns <- review_prior_columns(list(mu_g = prior))
    samples <- cbind(mu_intercept = c(0, 0), tau = .2, tau_indicator = 2L,
      mu_g_indicator = c(1L, 2L),
      matrix(c(0, .3), nrow = 2L, ncol = length(columns), dimnames = list(NULL, columns)))
    out <- .plot_models_data(review_attach_replay(setup, samples))
    expect_identical("mu" %in% out[["cells"]][["row"]], intercept)
  }
})


.plot_models_ordered_replay <- function(total, on, spike = FALSE) {

  data <- data.frame(yi = seq(.1, .6, by = .1), sei = .2,
    f = ordered(rep(c("a", "b", "c"), 2)))
  setup <- review_compile_object(BMA(yi = yi, sei = sei, data = data,
    mods = ~ f, measure = "GEN", prior_mods_null = list(f = NULL),
    prior_mods = list(f = BayesTools::prior_ordered(total, allocation = c(.25, .75))),
    prior_unit_information_sd = 1, only_priors = TRUE))
  sums <- 2 * as.numeric(on)
  samples <- cbind(mu_intercept = rep(.2, length(on)), mu_intercept_indicator = 2L,
    tau = .2, tau_indicator = 2L, mu_f_ordered_total = sums,
    `mu_f[1]` = .25 * sums, `mu_f[2]` = sums)
  if (BayesTools::is.prior.mixture(total)) {
    samples <- cbind(samples, mu_f_ordered_total_indicator = as.integer(on) + !spike)
    if (spike) {
      samples <- cbind(samples, mu_f_ordered_total_variable = 2,
        mu_f_ordered_total_inclusion = .4)
    } else {
      samples <- cbind(samples, mu_f_ordered_total_component_2 = 2)
    }
  }
  review_attach_replay(setup, samples)
}


test_that("ordered totals contribute their actual gate to the full joint model", {

  zero <- BayesTools::prior("point", list(0))
  slab <- BayesTools::prior("normal", list(0, 1))
  mixture <- BayesTools::prior_mixture(list(Zero = zero, Diff = slab), is_null = c(TRUE, FALSE))
  spike <- BayesTools::prior_spike_and_slab(slab, BayesTools::prior("beta", list(2, 3)))
  on <- c(FALSE, TRUE, FALSE, TRUE, FALSE, FALSE)

  for (total in list(mixture, spike)) {
    fit <- .plot_models_ordered_replay(total, on, spike = identical(total, spike))
    marginal <- summary_models(fit, include_mcmc_diagnostics = FALSE)[["marginal"]][["Location: f"]]
    included <- marginal[["Hypothesis"]] == "Alternative"
    expect_equal(marginal[["post_prob"]][included], mean(on))
    expect_equal(marginal[["prior_prob"]][included], if (identical(total, spike)) .4 else .5)
    individual <- summary_models(fit, type = "individual", include_mcmc_diagnostics = FALSE)
    data <- .plot_models_data(fit, "f")
    expected <- individual[["individual"]][["post_prob"]]
    expect_equal(unname(data[["probabilities"]]), expected[expected > 0])
    expect_identical(unique(data[["cells"]][["row"]]), c("mu_f[b]", "mu_f[c]"))
    expect_setequal(data[["cells"]][["state"]], c("positive", "excluded"))
    displayed <- data[["cells"]][["row"]] == "mu_f[b]" & data[["cells"]][["state"]] == "positive"
    expect_equal(unname(data[["probabilities"]][data[["cells"]][["model"]][displayed]]), mean(on))
    hidden <- .plot_models_data(fit, "mu")
    expect_identical(hidden[["probabilities"]], data[["probabilities"]])
  }

  fit <- .plot_models_ordered_replay(zero, rep(FALSE, 6L))
  fixed <- .plot_models_data(fit, "f")
  expect_identical(unname(fixed[["probabilities"]]), 1)
  expect_identical(unique(fixed[["cells"]][["state"]]), "excluded")
})
