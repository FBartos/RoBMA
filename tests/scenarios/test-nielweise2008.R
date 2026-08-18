if (file.exists("helper-scenarios.R")) source("helper-scenarios.R") else source("tests/scenarios/helper-scenarios.R")
scenario_start("nielweise2008")
# testthat::test_file("tests/scenarios/test-nielweise2008.R")

### Description
# compare Poisson incidence-rate models fitted with metafor, brma.glmm(), and
# BMA.glmm(), including posterior prediction of the observed event counts

testthat::test_that("Niel-Weise Poisson incidence-rate models", {

  set.seed(1)
  data("dat.nielweise2008", package = "metadat")
  dat <- dat.nielweise2008
  dat$label <- paste0(dat$study, ". ", dat$authors)

  ### Model fits ----
  fit_metafor <- scenario_fit("fit_metafor", {
    metafor::rma.glmm(measure = "IRR", x1i = x1i, t1i = t1i, x2i = x2i, t2i = t2i, slab = label, data = dat, model = "UM.FS")
  })
  fit_brma <- scenario_fit("fit_brma", {
    temp_fit <- brma.glmm(x1i = x1i, t1i = t1i, x2i = x2i, t2i = t2i, slab = label, data = dat, measure = "IRR", seed = 1, silent = TRUE)
    temp_fit <- add_loo(temp_fit)
    temp_fit <- add_marglik(temp_fit)
    return(temp_fit)
  })
  fit_brma_null <- scenario_fit("fit_brma_null", {
    temp_fit <- brma.glmm(x1i = x1i, t1i = t1i, x2i = x2i, t2i = t2i, slab = label,
                          data = dat, prior_effect = NULL, measure = "IRR", seed = 1, silent = TRUE)
    temp_fit <- add_loo(temp_fit)
    temp_fit <- add_marglik(temp_fit)
    return(temp_fit)
  })
  fit_BMA <- scenario_fit("fit_BMA", {
    temp_fit <- BMA.glmm(x1i = x1i, t1i = t1i, x2i = x2i, t2i = t2i, slab = label, data = dat, measure = "IRR", seed = 1, silent = TRUE)
    temp_fit <- add_loo(temp_fit)
    return(temp_fit)
  })

  ### model summaries ----
  fit_metafor
  scenario_text("summary-fit_brma",       summary(fit_brma))
  scenario_text("summary-fit_BMA",        summary(fit_BMA, conditional = TRUE))
  scenario_text("summary-models-fit_BMA", summary_models(fit_BMA))
  scenario_text("summary-heterogenity",   summary_heterogeneity(fit_brma))

  ### basic fit plots ----
  set.seed(1)
  scenario_plot("fit_posterior_mu", {
    plot(fit_brma, "mu", ylim = c(0, 2.5), prior = TRUE)
    lines(fit_brma, "mu", density_method = "qCMDE", lty = 2)

    lines(fit_BMA, "mu", conditional = TRUE, col = "blue")
    lines(fit_BMA, "mu", conditional = TRUE, density_method = "qCMDE", lty = 2, col = "blue")
  })

  scenario_plot("fit_posterior_tau", {
    plot(fit_brma, "tau", ylim = c(0, 3), prior = TRUE)
    lines(fit_brma, "tau", density_method = "qCMDE", lty = 2)

    lines(fit_BMA, "tau", col = "blue")
    lines(fit_BMA, "tau", density_method = "qCMDE", lty = 2, col = "blue")
  })

  # metafor::forest(fit_metafor)
  scenario_plot("fit_forest", metafor::forest(as_metafor_forest(fit_brma)))
  scenario_plot("fit_BMA_forest", metafor::forest(as_metafor_forest(fit_BMA)))

  ### hypothesis comparisons ----
  BF_mu_1 <- hypothesis(fit_brma, "mu = 0")[["BF"]]
  BF_mu_2 <- hypothesis(fit_brma, "mu = 0", density_method = "qCMDE")[["BF"]]
  BF_mu_3 <- hypothesis(fit_BMA, "mu = 0", density_method = "qCMDE", conditional = TRUE)[["BF"]]
  BF_mu   <- bf(fit_brma, fit_brma_null)[["bf"]]

  scenario_text("bf-effect", c("marglik" = BF_mu, "KDE" = BF_mu_1, "qCMDE" = BF_mu_2, "qCMDE (BMA)" = BF_mu_3))

  ### plot diagnostics ----
  scenario_plot("fit_funnel", funnel(fit_brma))
  scenario_plot("residuals",  {scenario_agreement_plot(residuals(fit_metafor), residuals(fit_brma))})

  ### predicted effects ----
  compare_preds <- function(fit_metafor, fit_RoBMA) {
    cbind.data.frame(
      "metafor"  = t(data.frame(predict(fit_metafor))[,-2]),
      "brma"     = t(unname(data.frame(pooled_effect(fit_RoBMA))[,-2]))
    )
  }

  set.seed(1)
  scenario_text("pooled-effect",  compare_preds(fit_metafor, fit_brma))
  scenario_text("pooled-effect-EXP",   pooled_effect(fit_brma, transform = "EXP"))
  scenario_text("pooled-heterogenity", pooled_heterogeneity(fit_brma))
})
