if (file.exists("helper-scenarios.R")) source("helper-scenarios.R") else source("tests/scenarios/helper-scenarios.R")
REGENERATE_SCENARIO_FILES <- FALSE
SHOW_SCENARIO_OUTPUT      <- FALSE
scenario_start("bangertdrowns2004")
# testthat::test_file("tests/scenarios/test-bangertdrowns2004.R")

### Description
# use the bangertdrowns2004 dataset to test location-scale related features

testthat::test_that("Bangertdrowns location-scale models", {
  set.seed(1)

  data(dat.bangertdrowns2004, package = "metadat")
  dat.bangertdrowns2004$ni100 <- dat.bangertdrowns2004$ni / 100

  # fit a standard random-effects model
  metafor_simple   <- metafor::rma(yi, vi, data = dat.bangertdrowns2004)
  metafor_l        <- metafor::rma(yi, vi, mods = ~ ni100, data = dat.bangertdrowns2004)
  metafor_s        <- metafor::rma(yi, vi, scale = ~ ni100, data = dat.bangertdrowns2004)
  metafor_ls       <- metafor::rma(yi, vi, mods = ~ ni100, scale = ~ ni100, data = dat.bangertdrowns2004)

  fit_simple <- scenario_fit("fit_simple", {
    tmp <- brma(yi = yi, vi = vi, data = dat.bangertdrowns2004, measure = "SMD", seed = 1)
    tmp <- add_loo(tmp)
    tmp <- add_marglik(tmp)
    return(tmp)
  })
  fit_l <- scenario_fit("fit_l", {
    tmp <- brma(yi = yi, vi = vi, mods = ~ ni100, data = dat.bangertdrowns2004, measure = "SMD", seed = 1)
    tmp <- add_loo(tmp)
    tmp <- add_marglik(tmp)
    return(tmp)
  })
  fit_s <- scenario_fit("fit_s", {
    tmp <- brma(yi = yi, vi = vi, scale = ~ ni100, data = dat.bangertdrowns2004, measure = "SMD", seed = 1)
    tmp <- add_loo(tmp)
    tmp <- add_marglik(tmp)
    return(tmp)
  })
  fit_ls <- scenario_fit("fit_ls", {
    tmp <- brma(yi = yi, vi = vi, mods = ~ ni100, scale = ~ ni100, data = dat.bangertdrowns2004, measure = "SMD", seed = 1)
    tmp <- add_loo(tmp)
    tmp <- add_marglik(tmp)
    return(tmp)
  })

  ### simple summary ----
  metafor_simple
  scenario_text("fit_simple_summary", {summary(fit_simple)})
  metafor_l
  scenario_text("fit_l_summary",  {summary(fit_l)})
  metafor_s # exp(-1.9612/2) -> to get log intercept
  scenario_text("fit_s_summary",  {summary(fit_s)})
  metafor_ls # exp(-1.9209/2) -> to get log intercept
  scenario_text("fit_ls_summary", {summary(fit_ls)})

  ### pooled effect summary
  set.seed(1)
  predict(metafor_simple)
  scenario_text("fit_simple_effect", {pooled_effect(fit_simple)})
  apply(as.data.frame(predict(metafor_l)), 2 , mean)
  scenario_text("fit_l_effect",      {pooled_effect(fit_l)})
  apply(as.data.frame(predict(metafor_s)), 2 , mean)
  scenario_text("fit_s_effect",      {pooled_effect(fit_s)})
  apply(as.data.frame(predict(metafor_ls)), 2 , mean)
  scenario_text("fit_ls_effect",     {pooled_effect(fit_ls)})

  set.seed(1)
  confint(metafor_simple)
  scenario_text("fit_simple_heterogeneity", {pooled_heterogeneity(fit_simple)})
  confint(metafor_l)
  scenario_text("fit_l_heterogeneity",      {pooled_heterogeneity(fit_l)})
  apply(as.data.frame(predict(metafor_s, newscale = TRUE)), 2 , function(x) sqrt(mean(exp(x))))
  scenario_text("fit_s_heterogeneity",      {pooled_heterogeneity(fit_s)})
  apply(as.data.frame(predict(metafor_ls, newscale = TRUE)), 2 , function(x) sqrt(mean(exp(x))))
  scenario_text("fit_ls_heterogeneity",     {pooled_heterogeneity(fit_ls)})


  ### simple model comparison ----
  scenario_text("fit_bf_comparison",  {post_prob(fit_simple, fit_l, fit_s, fit_ls)})
  scenario_text("fit_loo_comparison", {loo_model_weights(fit_simple, fit_l, fit_s, fit_ls)})

  ### basic fit plots ----
  set.seed(1)
  scenario_plot("fit_posterior_mu", {
    plot(fit_simple, "mu", ylim = c(0, 10), xlim = c(-0.2, 0.6), prior = TRUE)
    lines(fit_simple, "mu", density_method = "IWMDE", lty = 2)

    lines(fit_l, "mu", col = "blue")
    lines(fit_l, "mu", density_method = "IWMDE", lty = 2, col = "blue")

    lines(fit_s, "mu", col = "red")
    lines(fit_s, "mu", density_method = "qCMDE", lty = 2, col = "red")

    lines(fit_ls, "mu", col = "green")
    lines(fit_ls, "mu", density_method = "qCMDE", lty = 2, col = "green")
  })

  set.seed(1)
  scenario_plot("fit_posterior_tau", {
    plot(fit_simple, "tau", ylim = c(0, 10), xlim = c(0.0, 0.6), prior = TRUE)
    lines(fit_simple, "tau", density_method = "qCMDE", lty = 2)

    lines(fit_l, "tau", col = "blue")
    lines(fit_l, "tau", density_method = "qCMDE", lty = 2, col = "blue")

    # need to use `standardized_coefficients = TRUE`, otherwise its the intercept estimate at ni100 = 0
    lines(fit_s, "tau", col = "red", standardized_coefficients = TRUE)
    lines(fit_s, "tau", density_method = "qCMDE", lty = 2, col = "red", standardized_coefficients = TRUE, density_control = list(samples = 2500))

    lines(fit_ls, "tau", col = "green", standardized_coefficients = TRUE)
    lines(fit_ls, "tau", density_method = "qCMDE", lty = 2, col = "green", standardized_coefficients = TRUE, density_control = list(samples = 2500))
  })

  set.seed(1)
  scenario_plot("fit_posterior_mods", {
    plot(fit_l, "ni100", ylim = c(0, 20), xlim = c(-0.5, 0.5), prior = TRUE)
    lines(fit_l, "ni100", density_method = "IWMDE", lty = 2)

    lines(fit_ls, "ni100", col = "blue", component = "mods")
    lines(fit_ls, "ni100", density_method = "IWMDE", lty = 2, col = "blue", component = "mods")
  })

  set.seed(1)
  scenario_plot("fit_posterior_scale", {
    par(mfrow = c(1, 2))
    plot(fit_s, "ni100", ylim = c(0, 2), xlim = c(-1, 1), prior = TRUE)
    lines(fit_s, "ni100", density_method = "IWMDE", lty = 2)

    lines(fit_ls, "ni100", col = "blue", component = "scale")
    lines(fit_ls, "ni100", density_method = "IWMDE", lty = 2, col = "blue", component = "scale")

    plot(fit_s, "ni100", ylim = c(0, 3), xlim = c(0, 3), prior = TRUE, transform = "EXP")
    lines(fit_s, "ni100", density_method = "IWMDE", lty = 2, transform = "EXP")

    lines(fit_ls, "ni100", col = "blue", component = "scale", transform = "EXP")
    lines(fit_ls, "ni100", density_method = "IWMDE", lty = 2, col = "blue", component = "scale", transform = "EXP")
  })

  ### regression plots ----
  set.seed(1)
  scenario_plot("fit_l",  {
    par(mfrow = c(1, 2))
    regplot(fit_l, mod = "ni100", pi = TRUE, si = TRUE, ylim = c(-1, 1.5))
    metafor::regplot(metafor_l, mod = "ni100", pi = TRUE, ylim = c(-1, 1.5))
  })
  scenario_plot("fit_s",  {
    par(mfcol = c(1, 2))
    regplot(fit_s, mod = "ni100", pi = TRUE, si = TRUE, ylim = c(-1, 1.5))
    # metafor does not allow this
    # metafor::regplot(metafor_s, mod = "ni100")
  })
  scenario_plot("fit_ls",  {
    par(mfcol = c(1, 2))
    regplot(fit_ls, mod = "ni100", pi = TRUE, si = TRUE, ylim = c(-1, 1.5))
    metafor::regplot(metafor_ls, mod = "ni100", ylim = c(-1, 1.5)) # Cannot draw prediction interval for the given model.
  })

})
