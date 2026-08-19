if (file.exists("helper-scenarios.R")) source("helper-scenarios.R") else source("tests/scenarios/helper-scenarios.R")
scenario_start("hoogeveen2023")
# testthat::test_file("tests/scenarios/test-hoogeveen2023.R")

### Description
# Compare the exact rank-one sampling-covariance analysis of the Many-Analysts
# Religion Project with independent analysis effects and with the known quality
# covariance R.

testthat::test_that("Hoogeveen rank-one sampling covariance and known quality R", {

  data("Hoogeveen2023", package = "RoBMA")

  dat <- Hoogeveen2023[trimws(Hoogeveen2023[["type"]]) == "beta", ]
  dat$median_sei  <- median(dat$sei)
  dat$quality     <- dat$team_knowledge / 5
  dat$analysis    <- factor(seq_len(nrow(dat)))

  V_median  <- tcrossprod(dat$median_sei)
  R_quality <- diag(1 / dat$quality)
  dimnames(R_quality) <- list(levels(dat$analysis), levels(dat$analysis))

  R_correlated <- diag(1 / ((dat$yi - 0.02)/max(dat$yi - 0.02)))
  dimnames(R_correlated) <- list(levels(dat$analysis), levels(dat$analysis))

  uisd <- estimate_unit_information_sd(sei = dat$median_sei, ni = dat$ni)

  ### metafor reference analyses ----
  fit_metafor_mv <- metafor::rma.mv(yi = yi, V = V_median, random = ~ 1 | analysis,
                                    data = dat, method = "REML", test = "t")

  fit_metafor_mv_quality <- metafor::rma.mv(yi = yi, V = V_median, random = ~ 1 | analysis,
                                            R = list(analysis = R_quality), Rscale = "none",
                                            data = dat, method = "REML", test = "t")
  # fit with an altered rcor matrix to induce high dependency on R
  fit_metafor_mv_rcor <- metafor::rma.mv(yi = yi, V = V_median, random = ~ 1 | analysis,
                                            R = list(analysis = R_correlated), Rscale = "none",
                                            data = dat, method = "REML", test = "t")


  ### RoBMA analyses without known R ----
  fit_brma_mv <- scenario_fit("fit_brma_mv", {
    temp_fit <- brma.mv(yi = yi, V = V_median, random = ~ 1 | analysis, data = dat, measure = "GEN",
                        prior_unit_information_sd = uisd, seed = 1, silent = TRUE)
    temp_fit <- add_loo(temp_fit)
    temp_fit <- add_marglik(temp_fit)
    return(temp_fit)
  })
  fit_brma_mv_null <- scenario_fit("fit_brma_mv_null", {
    temp_fit <- brma.mv(yi = yi, V = V_median, random = ~ 1 | analysis, data = dat, measure = "GEN",
                        prior_effect = prior("spike", list(0)), prior_unit_information_sd = uisd,
                        seed = 1, silent = TRUE)
    temp_fit <- add_loo(temp_fit)
    temp_fit <- add_marglik(temp_fit)
    return(temp_fit)
  })

  ### RoBMA analyses with known quality R ----
  fit_brma_mv_quality <- scenario_fit("fit_brma_mv_quality", {
    temp_fit <- brma.mv(yi = yi, V = V_median, random = ~ 1 | analysis,
                        R = list(analysis = R_quality), Rscale = "none", data = dat, measure = "GEN",
                        prior_unit_information_sd = uisd, seed = 1, silent = TRUE)
    temp_fit <- add_loo(temp_fit)
    temp_fit <- add_marglik(temp_fit)
    return(temp_fit)
  })
  fit_brma_mv_quality_null <- scenario_fit("fit_brma_mv_quality_null", {
    temp_fit <- brma.mv(yi = yi, V = V_median, random = ~ 1 | analysis,
                        R = list(analysis = R_quality), Rscale = "none", data = dat, measure = "GEN",
                        prior_effect = prior("spike", list(0)), prior_unit_information_sd = uisd,
                        seed = 1, silent = TRUE)
    temp_fit <- add_loo(temp_fit)
    temp_fit <- add_marglik(temp_fit)
    return(temp_fit)
  })

  # fits with an altered rcor matrix to induce high dependency on R
  fit_brma_mv_rcor <- scenario_fit("fit_brma_mv_rcor", {
    temp_fit <- brma.mv(yi = yi, V = V_median, random = ~ 1 | analysis,
                        R = list(analysis = R_correlated), Rscale = "none", data = dat, measure = "GEN",
                        prior_unit_information_sd = uisd, seed = 1, silent = TRUE)
    temp_fit <- add_loo(temp_fit)
    temp_fit <- add_marglik(temp_fit)
    return(temp_fit)
  })
  fit_brma_mv_rcor_null <- scenario_fit("fit_brma_mv_rcor_null", {
    temp_fit <- brma.mv(yi = yi, V = V_median, random = ~ 1 | analysis,
                        R = list(analysis = R_correlated), Rscale = "none", data = dat, measure = "GEN",
                        prior_effect = prior("spike", list(0)), prior_unit_information_sd = uisd,
                        seed = 1, silent = TRUE)
    temp_fit <- add_loo(temp_fit)
    temp_fit <- add_marglik(temp_fit)
    return(temp_fit)
  })

  ### model summaries ----
  # TODO: the summarized parameters prbly need some improvements
  scenario_text("summary-no-known-R", summary(fit_brma_mv, include_mcmc_diagnostics = FALSE))
  scenario_text("summary-known-R",    summary(fit_brma_mv_quality, include_mcmc_diagnostics = FALSE))
  scenario_text("summary-cor-R",      summary(fit_brma_mv_rcor, include_mcmc_diagnostics = FALSE))

  scenario_text("summary_het-no-known-R", summary_heterogeneity(fit_brma_mv))
  scenario_text("summary_het-known-R",    summary_heterogeneity(fit_brma_mv_quality))
  scenario_text("summary_het-cor-R",      summary_heterogeneity(fit_brma_mv_rcor))


  summarize_metafor <- function(fit) {
    c(mu = as.numeric(stats::coef(fit)[[1L]]), se = fit[["se"]][[1L]],
      lower = fit[["ci.lb"]][[1L]], upper = fit[["ci.ub"]][[1L]], tau = sqrt(fit[["sigma2"]][[1L]]))
  }
  summarize_brma    <- function(fit) {
    mu  <- as.numeric(pooled_effect(fit))
    tau <- as.numeric(pooled_heterogeneity(fit))
    c(mu = mean(mu), se = sd(mu), lower = unname(quantile(mu, .025)),
      upper = unname(quantile(mu, .975)), tau = median(tau))
  }

  # compare RoBMA and metafor models
  scenario_text("model-comparison", {
    data.frame(
      model = rep(c("without known R", "with known R", "with corr R"), each = 2L),
      implementation = rep(c("metafor", "RoBMA"), 3L),
      round(rbind(
        summarize_metafor(fit_metafor_mv),         summarize_brma(fit_brma_mv),
        summarize_metafor(fit_metafor_mv_quality), summarize_brma(fit_brma_mv_quality),
        summarize_metafor(fit_metafor_mv_rcor),    summarize_brma(fit_brma_mv_rcor)
      ),4),row.names = NULL
    )
  })


  ### bridge-sampling and density Bayes factors ----
  effect_bayes_factors <- function(fit, null_fit) {
    data.frame(
      bridge  = bf(fit, null_fit)[["bf"]],
      KDE     = hypothesis(fit, "mu = 0")[["BF"]][[1L]],
      qCMDE   = hypothesis(fit, "mu = 0", density_method = "qCMDE", density_control = list(samples = 1000L))[["BF"]][[1L]],
      normal  = hypothesis(fit, "mu = 0", density_method = "normal")[["BF"]][[1L]]
    )
  }

  scenario_text("effect-bayes-factors", rbind(
    "without known R" = effect_bayes_factors(fit_brma_mv,         fit_brma_mv_null),
    "with known R"    = effect_bayes_factors(fit_brma_mv_quality, fit_brma_mv_quality_null),
    "with corr R "    = effect_bayes_factors(fit_brma_mv_rcor,    fit_brma_mv_rcor_null)
  ))


  # Both models use conditional estimate deletion for the exact known V.
  scenario_text("model-fit-loo", loo_model_weights(
    fit_brma_mv, fit_brma_mv_null,
    fit_brma_mv_quality, fit_brma_mv_quality_null,
    fit_brma_mv_rcor, fit_brma_mv_rcor_null))

  ### posterior and random effects ----
  scenario_plot("posterior-mu", {
    plot(fit_brma_mv, "mu", prior = TRUE, xlim = c(0, .30), col = "blue", ylim = c(0, 40))
    lines(fit_brma_mv, "mu", lty = 2, col = "blue", density_method = "qCMDE")

    lines(fit_brma_mv_quality, "mu", col = "green")
    lines(fit_brma_mv_quality, "mu", lty = 2, col = "green", density_method = "qCMDE")

    lines(fit_brma_mv_rcor, "mu", col = "red")
    lines(fit_brma_mv_rcor, "mu", lty = 2, col = "red", density_method = "qCMDE")
  })

  scenario_plot("posterior-tau", {
    plot(fit_brma_mv, "sd(intercept)", prior = TRUE, xlim = c(0, .10), col = "blue")
    # TODO: examine this takes very long:
    lines(fit_brma_mv, "sd(intercept)", lty = 2, col = "blue", density_method = "qCMDE")

    # TODO: does not export sd(intercept)
#    lines(fit_brma_mv_quality, "sd(intercept)", col = "green")
#    lines(fit_brma_mv_quality, "sd(intercept)", lty = 2, col = "green", density_method = "qCMDE")
#
#    lines(fit_brma_mv_rcor, "sd(intercept)", col = "red")
#    lines(fit_brma_mv_rcor, "sd(intercept)", lty = 2, col = "red", density_method = "qCMDE")
  })

  scenario_plot("random-effects", {
    par(mfrow = c(1, 3), mar = c(4, 4, 2, 1))
    scenario_agreement_plot(
      metafor::ranef(fit_metafor_mv)[["analysis"]][["intrcpt"]],
      as.data.frame(ranef(fit_brma_mv))[["Mean"]],
      main = "No known R"
    )
    scenario_agreement_plot(
      metafor::ranef(fit_metafor_mv_quality)[["analysis"]][["intrcpt"]],
      as.data.frame(ranef(fit_brma_mv_quality))[["Mean"]],
      main = "Known quality R"
    )
    scenario_agreement_plot(
      metafor::ranef(fit_metafor_mv_rcor)[["analysis"]][["intrcpt"]],
      as.data.frame(ranef(fit_brma_mv_rcor))[["Mean"]],
      main = "Known quality R"
    )
  })

  ### diagnostics ----
  plot_marginal_diagnostics <- function(fit_metafor, fit_brma) {
    metafor_values <- list(
      "Residuals"      = as.numeric(stats::residuals(fit_metafor)),
      "Rstandard"      = stats::rstandard(fit_metafor)[["z"]],
      "Hat values"     = as.numeric(stats::hatvalues(fit_metafor)),
      "Cooks distance" = stats::cooks.distance(fit_metafor),
      "DFBETAS"        = unlist(stats::dfbetas(fit_metafor))
    )
    brma_values <- list(
      "Residuals"      = residuals(fit_brma, type = "outcome", conditioning_depth = "marginal"),
      "Rstandard"      = rstandard(fit_brma, conditioning_depth = "marginal")[["z"]],
      "Hat values"     = hatvalues(fit_brma),
      "Cooks distance" = cooks.distance(fit_brma),
      "DFBETAS"        = unlist(dfbetas(fit_brma))
    )

    par(mfrow = c(3, 2), mar = c(4, 4, 2, 1))
    for (diagnostic in names(metafor_values)) {
      scenario_agreement_plot(metafor_values[[diagnostic]], brma_values[[diagnostic]], diagnostic)
    }
    invisible(NULL)
  }

  scenario_plot("diagnostics-no-known-R", plot_marginal_diagnostics(fit_metafor_mv, fit_brma_mv))
  scenario_plot("diagnostics-known-R",    plot_marginal_diagnostics(fit_metafor_mv_quality, fit_brma_mv_quality))
  scenario_plot("diagnostics-cor-R",      plot_marginal_diagnostics(fit_metafor_mv_rcor,    fit_brma_mv_rcor))

  ### diagnostic plots ----
  scenario_plot("funnel", {
    par(mfrow = c(2, 2))
    funnel(fit_brma_mv,         main = "Funnel: no known R")
    bfunnel(fit_brma_mv,        main = "Bayesian funnel: no known R")
    funnel(fit_brma_mv_quality, main = "Funnel: known R")
    funnel(fit_brma_mv_rcor,    main = "Funnel: cor R")
  })
  scenario_plot("qqnorm", {
    par(mfrow = c(1, 3))
    qqnorm(fit_brma_mv,         main = "No known R")
    qqnorm(fit_brma_mv_quality, main = "Known quality R")
    qqnorm(fit_brma_mv_rcor,    main = "Cor R")
  })
  scenario_plot("zplot", {
    zplot(fit_brma_mv_quality, from = 0, to = 20, step = 1, main = "Known quality R")
  })
})
