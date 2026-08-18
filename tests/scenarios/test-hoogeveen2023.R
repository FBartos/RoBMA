if (file.exists("helper-scenarios.R")) source("helper-scenarios.R") else source("tests/scenarios/helper-scenarios.R")
scenario_start("hoogeveen2023")
# testthat::test_file("tests/scenarios/test-hoogeveen2023.R")

### Description
# Compare the exact rank-one sampling-covariance analysis of the Many-Analysts
# Religion Project with independent analysis effects and with the known quality
# covariance R.

testthat::test_that("Hoogeveen rank-one sampling covariance and known quality R", {

  skip("not fully implemented")
  data("Hoogeveen2023", package = "RoBMA")

  dat <- Hoogeveen2023[trimws(Hoogeveen2023[["type"]]) == "beta", ]
  dat[["sei"]]      <- rep(stats::median(dat[["sei"]]), nrow(dat))
  dat[["quality"]]  <- dat[["team_knowledge"]] / 5
  dat[["analysis"]] <- factor(seq_len(nrow(dat)))

  v_dependent <- tcrossprod(dat[["sei"]])
  r_quality   <- diag(1 / dat[["quality"]])
  dimnames(r_quality) <- list(levels(dat[["analysis"]]), levels(dat[["analysis"]]))

  uisd <- estimate_unit_information_sd(sei = dat[["sei"]], ni = dat[["ni"]])

  ### metafor reference analyses ----
  fit_metafor_mv <- metafor::rma.mv(yi = yi, V = v_dependent, random = ~ 1 | analysis,
                                    data = dat, method = "REML", test = "t")

  fit_metafor_mv_quality <- metafor::rma.mv(yi = yi, V = v_dependent, random = ~ 1 | analysis,
                                            R = list(analysis = r_quality), Rscale = "none",
                                            data = dat, method = "REML", test = "t")

  ### RoBMA analyses without known R ----
  fit_brma_mv <- scenario_fit("fit_brma_mv", {
    temp_fit <- brma.mv(yi = yi, V = v_dependent, random = ~ 1 | analysis, data = dat, measure = "GEN",
                        prior_unit_information_sd = uisd, seed = 1, silent = TRUE)
    temp_fit <- add_loo(temp_fit)
    temp_fit <- add_marglik(temp_fit)
    return(temp_fit)
  })
  fit_brma_mv_null <- scenario_fit("fit_brma_mv_null", {
    temp_fit <- brma.mv(yi = yi, V = v_dependent, random = ~ 1 | analysis, data = dat, measure = "GEN",
                        prior_effect = prior("spike", list(0)), prior_unit_information_sd = uisd,
                        seed = 1, silent = TRUE)
    temp_fit <- add_marglik(temp_fit)
    return(temp_fit)
  })

  ### RoBMA analyses with known quality R ----
  fit_brma_mv_quality <- scenario_fit("fit_brma_mv_quality", {
    temp_fit <- brma.mv(yi = yi, V = v_dependent, random = ~ 1 | analysis,
                        R = list(analysis = r_quality), Rscale = "none", data = dat, measure = "GEN",
                        prior_unit_information_sd = uisd, seed = 1, silent = TRUE)
    temp_fit <- add_loo(temp_fit)
    temp_fit <- add_marglik(temp_fit)
    return(temp_fit)
  })
  fit_brma_mv_quality_null <- scenario_fit("fit_brma_mv_quality_null", {
    temp_fit <- brma.mv(yi = yi, V = v_dependent, random = ~ 1 | analysis,
                        R = list(analysis = r_quality), Rscale = "none", data = dat, measure = "GEN",
                        prior_effect = prior("spike", list(0)), prior_unit_information_sd = uisd,
                        seed = 1, silent = TRUE)
    temp_fit <- add_marglik(temp_fit)
    return(temp_fit)
  })

  ### model summaries ----
  # The two implementations target the same exact-V model; remaining
  # differences reflect classical point estimation versus posterior inference.
  scenario_text("summary-no-known-R", {
    print(fit_metafor_mv)
    summary(fit_brma_mv, include_mcmc_diagnostics = FALSE)
  })
  scenario_text("summary-known-R", {
    print(fit_metafor_mv_quality)
    summary(fit_brma_mv_quality, include_mcmc_diagnostics = FALSE)
  })

  summarize_metafor <- function(fit) {
    c(mu = as.numeric(stats::coef(fit)[[1L]]), se = fit[["se"]][[1L]],
      lower = fit[["ci.lb"]][[1L]], upper = fit[["ci.ub"]][[1L]],
      tau   = sqrt(fit[["sigma2"]][[1L]])
    )
  }
  summarize_brma <- function(fit) {
    mu  <- as.numeric(pooled_effect(fit))
    tau <- as.numeric(pooled_heterogeneity(fit))
    c(mu    = mean(mu),
      se    = stats::sd(mu),
      lower = unname(stats::quantile(mu, .025)),
      upper = unname(stats::quantile(mu, .975)),
      tau   = mean(tau)
    )
  }

  # The known-R model is expected to differ because quality changes the
  # row-specific heterogeneity.
  scenario_text("model-comparison", {
    data.frame(
      model = rep(c("without known R", "with known R"), each = 2L),
      implementation = rep(c("metafor", "RoBMA"), 2L),
      rbind(
        summarize_metafor(fit_metafor_mv),         summarize_brma(fit_brma_mv),
        summarize_metafor(fit_metafor_mv_quality), summarize_brma(fit_brma_mv_quality)
      ),
      row.names = NULL
    )
  })

  scenario_text("prior-specification", {
    cat("Exact V without known R:\n")
    print_prior(fit_brma_mv)
    cat("\nExact V with known quality R:\n")
    print_prior(fit_brma_mv_quality)
    invisible(NULL)
  })

  ### bridge-sampling and density Bayes factors ----
  effect_bayes_factors <- function(fit, null_fit) {
    data.frame(
      bridge = bf(fit, null_fit)[["bf"]],
      KDE     = hypothesis(fit, "mu = 0")[["BF"]][[1L]],
      qCMDE   = hypothesis(fit, "mu = 0", density_method = "qCMDE", density_control = list(samples = 2000L))[["BF"]][[1L]],
      normal  = hypothesis(fit, "mu = 0", density_method = "normal")[["BF"]][[1L]]
    )
  }

  # All columns estimate BF10 for the same fitted model. This extreme tail can
  # expose density-estimator limitations; non-finite or discrepant values are
  # intentionally left visible for maintainer review.
  scenario_text("effect-bayes-factors", rbind(
    `without known R` = effect_bayes_factors(fit_brma_mv,         fit_brma_mv_null),
    `with known R`    = effect_bayes_factors(fit_brma_mv_quality, fit_brma_mv_quality_null)
  ))

  loo_performance <- function(fit) {
    fit_loo  <- loo(fit)
    pareto_k <- fit_loo[["diagnostics"]][["pareto_k"]]
    data.frame(
      logml              = logml(fit),
      looic              = fit_loo[["estimates"]]["looic", 1L],
      maximum_pareto_k   = max(pareto_k),
      pareto_k_above_0.7 = sum(pareto_k > .7),
      pareto_k_above_1   = sum(pareto_k > 1)
    )
  }

  # Both models use conditional estimate deletion for the exact known V.
  scenario_text("model-fit", rbind(
    `without known R` = loo_performance(fit_brma_mv),
    `with known R`    = loo_performance(fit_brma_mv_quality)
  ))

  ### pooled effects and random effects ----
  scenario_plot("posterior-mu", {
    plot(fit_brma_mv, "mu", prior = TRUE, xlim = c(0, .18), lwd = 2, col = "blue")
    lines(fit_brma_mv_quality, "mu", lwd = 2, col = "darkgreen")
    legend("topleft", legend = c("without known R", "with known R"),
           col = c("blue", "darkgreen"), lty = 1, lwd = 2, bty = "n")
  })

  # BLUPs should show the same pattern, but posterior integration means that
  # pointwise identity with metafor is not expected.
  scenario_plot("random-effects", {
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    scenario_agreement_plot(
      metafor::ranef(fit_metafor_mv, expand = TRUE)[["analysis"]][["intrcpt"]],
      as.data.frame(ranef(fit_brma_mv, component = "total", expand = TRUE))[["Mean"]],
      main = "No known R"
    )
    scenario_agreement_plot(
      metafor::ranef(fit_metafor_mv_quality, expand = TRUE)[["analysis"]][["intrcpt"]],
      as.data.frame(ranef(fit_brma_mv_quality, component = "total", expand = TRUE))[["Mean"]],
      main = "Known quality R"
    )
  })

  ### diagnostics ----
  # These are visual comparisons rather than identities: RoBMA influence
  # diagnostics integrate posterior uncertainty and use PSIS where applicable.
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

  ### user-facing diagnostic plots ----
  scenario_plot("funnel", {
    par(mfrow = c(2, 2))
    funnel(fit_brma_mv,         main = "Funnel: no known R")
    bfunnel(fit_brma_mv,        main = "Bayesian funnel: no known R")
    funnel(fit_brma_mv_quality, main = "Funnel: known R")
    bfunnel(fit_brma_mv_quality, main = "Bayesian funnel: known R")
  })
  scenario_plot("qqnorm", {
    par(mfrow = c(1, 2))
    qqnorm(fit_brma_mv,         main = "No known R")
    qqnorm(fit_brma_mv_quality, main = "Known quality R")
  })
  scenario_plot("zplot", {
    par(mfrow = c(1, 2))
    zplot(fit_brma_mv,         to = 14, main = "No known R")
    zplot(fit_brma_mv_quality, to = 14, main = "Known quality R")
  })
})
