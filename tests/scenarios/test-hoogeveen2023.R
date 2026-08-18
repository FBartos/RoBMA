if (file.exists("helper-scenarios.R")) source("helper-scenarios.R") else source("tests/scenarios/helper-scenarios.R")
scenario_start("hoogeveen2023")
# testthat::test_file("tests/scenarios/test-hoogeveen2023.R")

### Description
# Compare the exact rank-one sampling-covariance analysis from the Many-
# Analysts Religion Project with its two-stage scalar representation. Repeat
# both analyses with the known diagonal quality covariance R. The scenario
# exercises singular known-V fitting, known-R row multipliers, bridge and
# density Bayes factors, prediction, random effects, and diagnostics.

testthat::test_that("Hoogeveen rank-one sampling covariance and known quality R", {

  set.seed(1)
  data("Hoogeveen2023", package = "RoBMA")

  dat <- Hoogeveen2023[trimws(Hoogeveen2023[["type"]]) == "beta", ]
  dat[["sei"]]      <- rep(stats::median(dat[["sei"]]), nrow(dat))
  dat[["quality"]]  <- dat[["team_knowledge"]] / 5
  dat[["analysis"]] <- factor(seq_len(nrow(dat)))

  v_dependent <- tcrossprod(dat[["sei"]])
  r_quality   <- diag(1 / dat[["quality"]])
  dimnames(r_quality) <- list(
    levels(dat[["analysis"]]),
    levels(dat[["analysis"]])
  )
  uisd <- estimate_unit_information_sd(
    sei = dat[["sei"]],
    ni  = dat[["ni"]]
  )

  # Both packages warn because V is intentionally positive semidefinite. The
  # fitted one-to-one random-intercept variance regularizes its null space.
  muffle_rank_one_warning <- function(code) {

    withCallingHandlers(
      code,
      warning = function(warning) {

        message <- conditionMessage(warning)
        expected <- startsWith(
          message,
          "The 'V' argument is positive semidefinite, not positive definite"
        ) || identical(message, "'V' appears to be not positive definite.")
        if (expected) {
          invokeRestart("muffleWarning")
        }
      }
    )
  }
  muffle_known_v_loo_warning <- function(code) {

    withCallingHandlers(
      code,
      warning = function(warning) {

        if (startsWith(
          conditionMessage(warning),
          "Estimate-unit LOO for brma.mv() known-V models uses conditional"
        )) {
          invokeRestart("muffleWarning")
        }
      }
    )
  }

  scenario_text("input-structure", data.frame(
    estimates                 = nrow(dat),
    V_rank                    = qr(v_dependent)[["rank"]],
    V_dimension               = nrow(v_dependent),
    common_se                 = unique(dat[["sei"]]),
    minimum_quality           = min(dat[["quality"]]),
    maximum_quality           = max(dat[["quality"]]),
    minimum_R_diagonal        = min(diag(r_quality)),
    maximum_R_diagonal        = max(diag(r_quality)),
    unit_information_sd       = uisd,
    stringsAsFactors          = FALSE
  ))

  ### metafor reference analyses ----
  fit_metafor_mv <- muffle_rank_one_warning(metafor::rma.mv(
    yi     = yi,
    V      = v_dependent,
    random = ~ 1 | analysis,
    data   = dat,
    method = "REML",
    test   = "t"
  ))
  fit_metafor_stage1 <- metafor::rma(
    yi     = yi,
    vi     = sei^2,
    data   = dat,
    method = "REML"
  )
  fit_metafor_scalar <- metafor::rma(
    yi     = yi,
    vi     = sei^2 * nrow(dat),
    tau2   = fit_metafor_stage1[["tau2"]],
    data   = dat,
    method = "REML"
  )

  fit_metafor_mv_quality <- muffle_rank_one_warning(metafor::rma.mv(
    yi     = yi,
    V      = v_dependent,
    random = ~ 1 | analysis,
    R      = list(analysis = r_quality),
    Rscale = "none",
    data   = dat,
    method = "REML",
    test   = "t"
  ))
  fit_metafor_quality_stage1 <- metafor::rma(
    yi     = yi,
    vi     = sei^2 / quality,
    data   = dat,
    method = "REML"
  )
  fit_metafor_scalar_quality <- metafor::rma(
    yi     = yi,
    vi     = sei^2 * nrow(dat) / quality,
    tau2   = fit_metafor_quality_stage1[["tau2"]],
    data   = dat,
    method = "REML"
  )

  ### RoBMA analyses without known R ----
  fit_brma_mv <- scenario_fit("fit_brma_mv", {
    temp_fit <- muffle_rank_one_warning(brma.mv(
      yi                        = yi,
      V                         = v_dependent,
      random                    = ~ 1 | analysis,
      data                      = dat,
      measure                   = "GEN",
      prior_unit_information_sd = uisd,
      seed                      = 1,
      silent                    = TRUE
    ))
    temp_fit <- muffle_known_v_loo_warning(add_loo(temp_fit))
    temp_fit <- add_marglik(temp_fit)
    return(temp_fit)
  })
  fit_brma_mv_null <- scenario_fit("fit_brma_mv_null", {
    temp_fit <- muffle_rank_one_warning(brma.mv(
      yi                        = yi,
      V                         = v_dependent,
      random                    = ~ 1 | analysis,
      data                      = dat,
      measure                   = "GEN",
      prior_effect              = prior("spike", list(0)),
      prior_unit_information_sd = uisd,
      seed                      = 1,
      silent                    = TRUE
    ))
    temp_fit <- add_marglik(temp_fit)
    return(temp_fit)
  })
  fit_brma_stage1 <- scenario_fit("fit_brma_stage1", {
    brma(
      yi                        = yi,
      vi                        = sei^2,
      data                      = dat,
      measure                   = "GEN",
      prior_unit_information_sd = uisd,
      seed                      = 1,
      silent                    = TRUE
    )
  })
  tau_no_r <- as.numeric(stats::coef(fit_brma_stage1)[["tau"]])

  fit_brma_scalar <- scenario_fit("fit_brma_scalar", {
    temp_fit <- brma(
      yi                        = yi,
      vi                        = sei^2 * nrow(dat),
      data                      = dat,
      measure                   = "GEN",
      prior_heterogeneity       = prior("spike", list(tau_no_r)),
      prior_unit_information_sd = uisd,
      seed                      = 1,
      silent                    = TRUE
    )
    temp_fit <- add_loo(temp_fit)
    temp_fit <- add_marglik(temp_fit)
    return(temp_fit)
  })
  fit_brma_scalar_null <- scenario_fit("fit_brma_scalar_null", {
    temp_fit <- brma(
      yi                        = yi,
      vi                        = sei^2 * nrow(dat),
      data                      = dat,
      measure                   = "GEN",
      prior_effect              = prior("spike", list(0)),
      prior_heterogeneity       = prior("spike", list(tau_no_r)),
      prior_unit_information_sd = uisd,
      seed                      = 1,
      silent                    = TRUE
    )
    temp_fit <- add_marglik(temp_fit)
    return(temp_fit)
  })

  ### RoBMA analyses with known quality R ----
  fit_brma_mv_quality <- scenario_fit("fit_brma_mv_quality", {
    temp_fit <- muffle_rank_one_warning(brma.mv(
      yi                        = yi,
      V                         = v_dependent,
      random                    = ~ 1 | analysis,
      R                         = list(analysis = r_quality),
      Rscale                    = "none",
      data                      = dat,
      measure                   = "GEN",
      prior_unit_information_sd = uisd,
      seed                      = 1,
      silent                    = TRUE
    ))
    temp_fit <- muffle_known_v_loo_warning(add_loo(temp_fit))
    temp_fit <- add_marglik(temp_fit)
    return(temp_fit)
  })
  fit_brma_mv_quality_null <- scenario_fit("fit_brma_mv_quality_null", {
    temp_fit <- muffle_rank_one_warning(brma.mv(
      yi                        = yi,
      V                         = v_dependent,
      random                    = ~ 1 | analysis,
      R                         = list(analysis = r_quality),
      Rscale                    = "none",
      data                      = dat,
      measure                   = "GEN",
      prior_effect              = prior("spike", list(0)),
      prior_unit_information_sd = uisd,
      seed                      = 1,
      silent                    = TRUE
    ))
    temp_fit <- add_marglik(temp_fit)
    return(temp_fit)
  })
  fit_brma_quality_stage1 <- scenario_fit("fit_brma_quality_stage1", {
    brma(
      yi                        = yi,
      vi                        = sei^2 / quality,
      data                      = dat,
      measure                   = "GEN",
      prior_unit_information_sd = uisd,
      seed                      = 1,
      silent                    = TRUE
    )
  })
  tau_quality <- as.numeric(stats::coef(fit_brma_quality_stage1)[["tau"]])

  fit_brma_scalar_quality <- scenario_fit("fit_brma_scalar_quality", {
    temp_fit <- brma(
      yi                        = yi,
      vi                        = sei^2 * nrow(dat) / quality,
      data                      = dat,
      measure                   = "GEN",
      prior_heterogeneity       = prior("spike", list(tau_quality)),
      prior_unit_information_sd = uisd,
      seed                      = 1,
      silent                    = TRUE
    )
    temp_fit <- add_loo(temp_fit)
    temp_fit <- add_marglik(temp_fit)
    return(temp_fit)
  })
  fit_brma_scalar_quality_null <- scenario_fit("fit_brma_scalar_quality_null", {
    temp_fit <- brma(
      yi                        = yi,
      vi                        = sei^2 * nrow(dat) / quality,
      data                      = dat,
      measure                   = "GEN",
      prior_effect              = prior("spike", list(0)),
      prior_heterogeneity       = prior("spike", list(tau_quality)),
      prior_unit_information_sd = uisd,
      seed                      = 1,
      silent                    = TRUE
    )
    temp_fit <- add_marglik(temp_fit)
    return(temp_fit)
  })

  ### summaries and reference agreement ----
  summarize_metafor <- function(fit) {

    tau <- if (inherits(fit, "rma.mv")) {
      sqrt(fit[["sigma2"]][[1L]])
    } else {
      sqrt(fit[["tau2"]])
    }
    return(data.frame(
      mu    = as.numeric(stats::coef(fit)[[1L]]),
      se    = fit[["se"]][[1L]],
      lower = fit[["ci.lb"]][[1L]],
      upper = fit[["ci.ub"]][[1L]],
      p     = fit[["pval"]][[1L]],
      tau   = tau
    ))
  }
  scenario_text("metafor-model-comparison", {
    output <- rbind(
      exact_V              = summarize_metafor(fit_metafor_mv),
      scalar               = summarize_metafor(fit_metafor_scalar),
      exact_V_known_R      = summarize_metafor(fit_metafor_mv_quality),
      scalar_quality       = summarize_metafor(fit_metafor_scalar_quality)
    )
    output[["model"]] <- rownames(output)
    rownames(output)  <- NULL
    output[c("model", "mu", "se", "lower", "upper", "p", "tau")]
  })

  summarize_known_v <- function(fit) {

    known_v <- attr(fit[["data"]], "known_V_data")
    return(data.frame(
      requested_parameterization = known_v[["parameterization_requested"]],
      selected_parameterization  = known_v[["parameterization"]],
      effective_backend          = known_v[["effective_backend"]],
      singular                   = known_v[["singular"]],
      rank                       = known_v[["rank"]],
      dependency_blocks          = length(known_v[["block_indices"]]),
      stringsAsFactors           = FALSE
    ))
  }
  scenario_text("known-V-backend", rbind(
    no_known_R = summarize_known_v(fit_brma_mv),
    known_R    = summarize_known_v(fit_brma_mv_quality)
  ))

  scenario_text("summary-fit_brma_mv", {
    summary(fit_brma_mv, include_mcmc_diagnostics = FALSE)
  })
  scenario_text("summary-fit_brma_scalar", {
    summary(fit_brma_scalar, include_mcmc_diagnostics = FALSE)
  })
  scenario_text("summary-fit_brma_mv_quality", {
    summary(fit_brma_mv_quality, include_mcmc_diagnostics = FALSE)
  })
  scenario_text("summary-fit_brma_scalar_quality", {
    summary(fit_brma_scalar_quality, include_mcmc_diagnostics = FALSE)
  })

  bind_sample_summaries <- function(objects) {

    tables <- lapply(objects, function(object) {
      table           <- as.data.frame(object)
      rownames(table) <- NULL
      return(table)
    })
    output            <- do.call(rbind, tables)
    rows_per_model    <- vapply(tables, nrow, integer(1))
    output[["model"]] <- rep(names(objects), rows_per_model)
    output[c("model", setdiff(names(output), "model"))]
  }
  scenario_text("pooled-effects", bind_sample_summaries(list(
    exact_V         = pooled_effect(fit_brma_mv),
    scalar          = pooled_effect(fit_brma_scalar),
    exact_V_known_R = pooled_effect(fit_brma_mv_quality),
    scalar_quality  = pooled_effect(fit_brma_scalar_quality)
  )))
  scenario_text("pooled-heterogeneity", bind_sample_summaries(list(
    exact_V         = pooled_heterogeneity(fit_brma_mv),
    scalar          = pooled_heterogeneity(fit_brma_scalar),
    exact_V_known_R = pooled_heterogeneity(fit_brma_mv_quality),
    scalar_quality  = pooled_heterogeneity(fit_brma_scalar_quality)
  )))
  scenario_text("heterogeneity-exact-V", list(
    identity_R = summary_heterogeneity(fit_brma_mv),
    known_R    = summary_heterogeneity(fit_brma_mv_quality)
  ))

  scenario_text("prior-specification", {
    cat("Exact V without known R:\n")
    print_prior(fit_brma_mv)
    cat("\nExact V with known quality R:\n")
    print_prior(fit_brma_mv_quality)
    cat("\nScalar stage-one model:\n")
    print_prior(fit_brma_stage1)
    invisible(NULL)
  })

  ### bridge-sampling and density Bayes factors ----
  effect_bayes_factors <- function(fit, null_fit) {

    set.seed(1)
    bf_kde <- hypothesis(fit, "mu = 0")[["BF"]][[1L]]
    bf_qcmde <- hypothesis(
      fit,
      "mu = 0",
      density_method  = "qCMDE",
      density_control = list(samples = 2000L)
    )[["BF"]][[1L]]
    bf_normal <- hypothesis(
      fit,
      "mu = 0",
      density_method = "normal"
    )[["BF"]][[1L]]
    return(data.frame(
      bridge = bf(fit, null_fit)[["bf"]],
      KDE     = bf_kde,
      qCMDE   = bf_qcmde,
      normal  = bf_normal
    ))
  }
  scenario_text("effect-bayes-factors", rbind(
    exact_V = effect_bayes_factors(
      fit_brma_mv,
      fit_brma_mv_null
    ),
    scalar = effect_bayes_factors(
      fit_brma_scalar,
      fit_brma_scalar_null
    ),
    exact_V_known_R = effect_bayes_factors(
      fit_brma_mv_quality,
      fit_brma_mv_quality_null
    ),
    scalar_quality = effect_bayes_factors(
      fit_brma_scalar_quality,
      fit_brma_scalar_quality_null
    )
  ))

  loo_performance <- function(fit, target) {

    fit_loo  <- loo(fit)
    pareto_k <- fit_loo[["diagnostics"]][["pareto_k"]]
    return(data.frame(
      loo_target              = target,
      log_marginal_likelihood = logml(fit),
      elpd_loo                = fit_loo[["estimates"]]["elpd_loo", 1L],
      looic                   = fit_loo[["estimates"]]["looic", 1L],
      maximum_pareto_k        = max(pareto_k),
      pareto_k_above_0.7      = sum(pareto_k > 0.7),
      pareto_k_above_1        = sum(pareto_k > 1)
    ))
  }
  scenario_text("model-fit", rbind(
    exact_V = loo_performance(
      fit_brma_mv,
      "conditional estimate deletion"
    ),
    scalar = loo_performance(
      fit_brma_scalar,
      "independent estimate deletion"
    ),
    exact_V_known_R = loo_performance(
      fit_brma_mv_quality,
      "conditional estimate deletion"
    ),
    scalar_quality = loo_performance(
      fit_brma_scalar_quality,
      "independent estimate deletion"
    )
  ))

  ### posterior and random-effect comparisons ----
  scenario_plot("posterior-mu", {
    plot(
      fit_brma_mv,
      "mu",
      prior = TRUE,
      xlim  = c(0, 0.18),
      lwd   = 2,
      col   = "blue"
    )
    lines(fit_brma_scalar, "mu", lwd = 2, col = "darkgreen")
    lines(fit_brma_mv_quality, "mu", lwd = 2, col = "blue", lty = 2)
    lines(
      fit_brma_scalar_quality,
      "mu",
      lwd = 2,
      col = "darkgreen",
      lty = 2
    )
    legend(
      "topleft",
      legend = c("exact V", "scalar", "exact V + known R", "scalar quality"),
      col    = c("blue", "darkgreen", "blue", "darkgreen"),
      lty    = c(1, 1, 2, 2),
      lwd    = 2,
      bty    = "n"
    )
  })

  random_effect_means <- function(fit) {

    as.data.frame(ranef(fit, component = "total", expand = TRUE))[["Mean"]]
  }
  scenario_plot("random-effects", {
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    scenario_agreement_plot(
      metafor::ranef(fit_metafor_mv, expand = TRUE)[["analysis"]][["intrcpt"]],
      random_effect_means(fit_brma_mv),
      main = "No known R"
    )
    scenario_agreement_plot(
      metafor::ranef(fit_metafor_mv_quality, expand = TRUE)[["analysis"]][["intrcpt"]],
      random_effect_means(fit_brma_mv_quality),
      main = "Known quality R"
    )
  })

  ### numerical diagnostics ----
  plot_marginal_diagnostics <- function(fit_metafor, fit_robma) {

    metafor_values <- list(
      "Residuals"      = as.numeric(stats::residuals(fit_metafor)),
      "Rstandard"      = stats::rstandard(fit_metafor)[["z"]],
      "Hat values"     = as.numeric(stats::hatvalues(fit_metafor)),
      "Cooks distance" = stats::cooks.distance(fit_metafor),
      "DFBETAS"        = unlist(stats::dfbetas(fit_metafor))
    )
    robma_values <- list(
      "Residuals" = residuals(
        fit_robma,
        type               = "outcome",
        conditioning_depth = "marginal"
      ),
      "Rstandard" = suppressWarnings(rstandard(
        fit_robma,
        conditioning_depth = "marginal"
      ))[["z"]],
      "Hat values"     = suppressWarnings(hatvalues(fit_robma)),
      "Cooks distance" = suppressWarnings(cooks.distance(fit_robma)),
      "DFBETAS"        = unlist(suppressWarnings(dfbetas(fit_robma)))
    )

    par(mfrow = c(3, 2), mar = c(4, 4, 2, 1))
    for (diagnostic in names(metafor_values)) {
      scenario_agreement_plot(
        metafor_values[[diagnostic]],
        robma_values[[diagnostic]],
        diagnostic
      )
    }
    return(invisible(NULL))
  }
  scenario_plot("diagnostics-no-known-R", {
    plot_marginal_diagnostics(fit_metafor_mv, fit_brma_mv)
  })
  scenario_plot("diagnostics-known-R", {
    plot_marginal_diagnostics(fit_metafor_mv_quality, fit_brma_mv_quality)
  })

  ### graphical diagnostics ----
  set.seed(1)
  scenario_plot("funnel", {
    par(mfrow = c(2, 2))
    funnel(fit_brma_mv, main = "Funnel: no known R")
    bfunnel(fit_brma_mv, main = "Bayesian funnel: no known R")
    funnel(fit_brma_mv_quality, main = "Funnel: known R")
    bfunnel(fit_brma_mv_quality, main = "Bayesian funnel: known R")
  })
  scenario_plot("qqnorm", {
    par(mfrow = c(1, 2))
    qqnorm(fit_brma_mv, main = "No known R")
    qqnorm(fit_brma_mv_quality, main = "Known quality R")
  })
  scenario_plot("zplot", {
    par(mfrow = c(1, 2))
    zplot(fit_brma_mv, to = 14, main = "No known R")
    zplot(fit_brma_mv_quality, to = 14, main = "Known quality R")
  })
})
