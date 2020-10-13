context("(4) fitting and updating functions")
skip_on_cran()


# test objects
saved_fits   <- readRDS(file = "../results/saved_fits.RDS")
saved_fits2  <- readRDS(file = "../results/saved_fits2.RDS")
updated_fits <- readRDS(file = "../results/updated_fits.RDS")
remove_time  <- function(fit){
  for(m in 1:length(fit$models)){
    if(is.null(fit$models[[m]]$fit))next
    fit$models[[m]]$fit$timetaken       <- NULL
    fit$models[[m]]$fit$runjags.version <- NULL
  }
  return(fit)
}

# create mock data
k  <- 3
n  <- rep(20, 3)
r  <- c(.1, .2, .3)
t  <- psych::r2t(r, n - 2)
d  <- psych::t2d(t, n)
se <- c(psych::d.ci(d, n)[,2] - psych::d.ci(d, n)[,1]) / 1.96

# fit default priors model configuration
fit_default <- RoBMA(t = t, n = n, chains = 2, burnin = 1000, iter = 4000, control = list(silent = TRUE), seed = 666)

test_that("Default model fit works", {
  fit_default <- remove_time(fit_default)
  expect_equal(saved_fits[[1]], fit_default)
})


# fit custom fits
fit_custom1 <- RoBMA(r = r, n = n, chains = 2, burnin = 1000, iter = 4000, control = list(autofit = FALSE, silent = TRUE), seed = 666, save = "min",
                     priors_mu = NULL, priors_tau = NULL, priors_omega = NULL)
fit_custom2 <- RoBMA(r = r, n = n, chains = 2, burnin = 1000, iter = 4000, control = list(autofit = FALSE, silent = TRUE), seed = 666, save = "min",
                     priors_mu_null = NULL, priors_tau_null = NULL, priors_omega_null = NULL)
fit_custom3 <- RoBMA(r = r, n = n, chains = 2, burnin = 1000, iter = 4000, control = list(autofit = FALSE, silent = TRUE), seed = 666, save = "min",
                     priors_mu = NULL, priors_tau = NULL, priors_omega_null = NULL,
                     priors_omega = list(
                       prior("two.sided", parameters = list(steps = c(.20),           alpha = c(1,10))),
                       prior("one.sided", parameters = list(steps = c(.10, .30),      alpha = c(1,5,1))),
                       prior("one.sided", parameters = list(steps = c(.15, .25, .75), alpha1 = c(1,1,1), alpha2 = c(1,1)))
                     ))
fit_custom4 <- RoBMA(r = r, n = n, chains = 2, burnin = 1000, iter = 4000, control = list(autofit = FALSE, silent = TRUE), seed = 666, save = "min",
                     priors_mu = NULL, priors_tau = NULL, priors_omega = NULL,
                     priors_mu_null = list(
                       prior("point",     parameters = list(location = 1)),
                       prior("normal",    parameters = list(mean = 0, sd = 1)),
                       prior("normal",    parameters = list(mean = 1, sd = 1),                truncation = list(lower = 0, upper = 2), prior_odds = 2),
                       prior("cauchy",    parameters = list(location = 1, scale = 1),         truncation = list(lower = 0, upper = 2)),
                       prior("t",         parameters = list(location = 1, scale = 1, df = 1), truncation = list(lower = 0, upper = 2)),
                       prior("t",         parameters = list(location = 1, scale = 1, df = 5), truncation = list(lower = -2, upper = 2)),
                       prior("gamma",     parameters = list(shape = 1, rate  = 2),            truncation = list(lower = 1, upper = Inf)),
                       prior("gamma",     parameters = list(shape = 1, scale = 1/2),          truncation = list(lower = 1, upper = Inf)),
                       prior("invgamma",  parameters = list(shape = 1, scale = .15),          truncation = list(lower = 0, upper = Inf)),
                       prior("uniform",   parameters = list(a = 2, b = 3))
                     ))
fit_custom5 <- RoBMA(r = r, n = n, chains = 2, burnin = 1000, iter = 4000, control = list(autofit = FALSE, silent = TRUE), seed = 666, save = "min",
                     priors_mu = NULL, priors_tau_null = NULL, priors_omega = NULL,
                     priors_tau = list(
                       prior("point",     parameters = list(location = 0)),
                       prior("normal",    parameters = list(mean = 1, sd = 1),                truncation = list(lower = 0, upper = 2), prior_odds = 2),
                       prior("cauchy",    parameters = list(location = 1, scale = 1),         truncation = list(lower = 0, upper = 2)),
                       prior("t",         parameters = list(location = 1, scale = 1, df = 1), truncation = list(lower = 0, upper = 2)),
                       prior("t",         parameters = list(location = 1, scale = 1, df = 5), truncation = list(lower = .5, upper = 2)),
                       prior("gamma",     parameters = list(shape = 1, rate  = 2),            truncation = list(lower = 1, upper = Inf)),
                       prior("gamma",     parameters = list(shape = 1, scale = 1/2),          truncation = list(lower = 1, upper = Inf)),
                       prior("invgamma",  parameters = list(shape = 1, scale = .15),          truncation = list(lower = 0, upper = Inf)),
                       prior("uniform",   parameters = list(a = .5, b = 1.5))
                     ))

test_that("Custom models works", {

  fit_custom1 <- remove_time(fit_custom1)
  fit_custom2 <- remove_time(fit_custom2)
  fit_custom3 <- remove_time(fit_custom3)
  fit_custom4 <- remove_time(fit_custom4)
  fit_custom5 <- remove_time(fit_custom5)

  expect_equal(saved_fits[[2]], fit_custom1)
  expect_equal(saved_fits[[3]], fit_custom2)
  expect_equal(saved_fits[[4]], fit_custom3)
  expect_equal(saved_fits[[5]], fit_custom4)
  expect_equal(saved_fits[[6]], fit_custom5)
})


# fit failling models
fit_fail1  <- suppressWarnings(RoBMA(t = t, n = n, chains = 2, burnin = 1000, iter = 4000, control = list(autofit = FALSE, silent = TRUE, bridge_max_iter = 2), seed = 666, save = "min",
                    priors_tau = NULL, priors_omega = NULL,
                    priors_mu  = prior("uniform", parameters = list(a = 10, b = 11))))

fit_fail2  <- suppressWarnings(RoBMA(t = t, n = n, chains = 2, burnin = 1000, iter = 4000, control = list(autofit = FALSE, silent = TRUE, allow_min_ESS = 2000), seed = 666, save = "min",
                    priors_tau = NULL, priors_omega = NULL,
                    priors_mu  = prior("uniform", parameters = list(a = 10, b = 11))))


test_that("Error handling works", {

  fit_fail1 <- remove_time(fit_fail1)
  fit_fail2 <- remove_time(fit_fail2)

  expect_equal(saved_fits[[7]], fit_fail1)
  expect_equal(saved_fits[[8]], fit_fail2)

})


# test additional settings
fit_settings  <- RoBMA(t = t, n = n, chains = 3, burnin = 1001, iter = 8002, thin = 2,
                       control = list(autofit = FALSE, silent = TRUE, adapt = 101), seed = 666,
                       priors_tau = NULL, priors_omega = NULL, priors_mu_null = NULL,
                       priors_mu  = prior("uniform", parameters = list(a = -.5, b = .5)))


test_that("Main settings work", {

  expect_equal(fit_settings$models[[1]]$fit$thin,   2)
  expect_equal(fit_settings$models[[1]]$fit$burnin, 1102)
  expect_equal(fit_settings$models[[1]]$fit$sample, 8002)

  expect_equal(length(fit_settings$models[[1]]$fit$mcmc),   3)
  expect_equal(dim(fit_settings$models[[1]]$fit$mcmc[[1]]), c(8002, 1))

})


# test updating model fits
fit_update1 <- RoBMA(t = t, n = n, chains = 2, burnin = 1000, iter = 4000, control = list(autofit = FALSE, silent = TRUE), seed = 666,
                     priors_mu = NULL, priors_tau = NULL, priors_omega = NULL)

fit_update1a<- update(fit_update1,
                      prior_mu  = prior("normal",        parameters = list(mean = 1, sd = 1)),
                      prior_tau = prior("uniform",       parameters = list(a = 0, b = .5), prior_odds = 2),
                      prior_omega_null = prior("point",  parameters = list(location = 1)))

fit_update1b<- update(fit_update1,
                      prior_mu  = prior("normal",        parameters = list(mean = 1, sd = 1)),
                      prior_tau = prior("uniform",       parameters = list(a = 0, b = .5)),
                      prior_omega_null = prior("point",  parameters = list(location = 1)),
                      prior_odds = 2)

fit_update1c<- update(fit_update1,
                      prior_mu  = prior("normal",        parameters = list(mean = 1, sd = 1)),
                      prior_tau = prior("uniform",       parameters = list(a = 0, b = .5)),
                      prior_omega_null = prior("point",  parameters = list(location = 1)))
fit_update1c<- update(fit_update1c, prior_odds = c(1, 2))


test_that("Adding models and changing prior odds works", {

  fit_update1a <- remove_time(fit_update1a)
  fit_update1b <- remove_time(fit_update1b)
  fit_update1c <- remove_time(fit_update1c)

  expect_equal(updated_fits[[1]]$RoBMA, fit_update1a$RoBMA)
  expect_equal(updated_fits[[1]]$RoBMA, fit_update1b$RoBMA)
  expect_equal(updated_fits[[1]]$RoBMA, fit_update1c$RoBMA)

})



# test refitting failed models (important that it updates only the failed model)
fit_update2  <- suppressWarnings(RoBMA(t = t, n = n, chains = 2, burnin = 1000, iter = 4000, control = list(autofit = FALSE, silent = TRUE, allow_min_ESS = 2000), seed = 666,
                                       priors_tau = NULL, priors_omega = NULL,
                                       priors_mu  = list(
                                         prior("uniform", parameters = list(a = 1, b = 2)),
                                         prior("uniform", parameters = list(a = 0, b = 1)))))
fit_update2 <- suppressWarnings(update(fit_update2, iter = 8000))

test_that("Updating failed models works", {

  fit_update2 <- remove_time(fit_update2)
  expect_equal(updated_fits[[2]], fit_update2)

})

# only changing the settings
fit_update3  <- suppressWarnings(RoBMA(t = t, n = n, chains = 2, burnin = 1000, iter = 4000, control = list(autofit = FALSE, silent = TRUE, allow_min_ESS = 6000), seed = 666,
                                       priors_tau = NULL, priors_omega = NULL,
                                       priors_mu  = prior("uniform", parameters = list(a = 0, b = 1))))
fit_update3 <- suppressWarnings(update(fit_update3, refit_failed = FALSE, control = list(allow_min_ESS = NULL)))

test_that("Updating failed models works", {

  fit_update3 <- remove_time(fit_update3)
  expect_equal(updated_fits[[3]], fit_update3)

})


# test model preview
test_that("Model preview works", {
  expect_equal(
    capture.output(check_setup(models = T)),
    c("Robust Bayesian Meta-Analysis (Set-Up)",
      "              Models Prior prob.",
      "Effect          6/12       0.500",
      "Heterogeneity   6/12       0.500",
      "Pub. bias       8/12       0.500",
      ""                                ,
      "Models Overview"                 ,
      "                  Prior mu                 Prior tau                       Prior omega Prior prob.",
      "1                 Spike(0)                  Spike(0)                          Spike(1)       0.125",
      "2                 Spike(0)                  Spike(0)         Two-sided((0.05), (1, 1))       0.062",
      "3                 Spike(0)                  Spike(0) Two-sided((0.1, 0.05), (1, 1, 1))       0.062",
      "4                 Spike(0) InvGamma(1, 0.15)[0, Inf]                          Spike(1)       0.125",
      "5                 Spike(0) InvGamma(1, 0.15)[0, Inf]         Two-sided((0.05), (1, 1))       0.062",
      "6                 Spike(0) InvGamma(1, 0.15)[0, Inf] Two-sided((0.1, 0.05), (1, 1, 1))       0.062",
      "7  Normal(0, 1)[-Inf, Inf]                  Spike(0)                          Spike(1)       0.125",
      "8  Normal(0, 1)[-Inf, Inf]                  Spike(0)         Two-sided((0.05), (1, 1))       0.062",
      "9  Normal(0, 1)[-Inf, Inf]                  Spike(0) Two-sided((0.1, 0.05), (1, 1, 1))       0.062",
      "10 Normal(0, 1)[-Inf, Inf] InvGamma(1, 0.15)[0, Inf]                          Spike(1)       0.125",
      "11 Normal(0, 1)[-Inf, Inf] InvGamma(1, 0.15)[0, Inf]         Two-sided((0.05), (1, 1))       0.062",
      "12 Normal(0, 1)[-Inf, Inf] InvGamma(1, 0.15)[0, Inf] Two-sided((0.1, 0.05), (1, 1, 1))       0.062")
  )
})


# additional model fits (not used for plotting and other tests)
fit_new1  <- suppressWarnings(RoBMA(OR = c(1.1, 1.05, 1.15), lCI = c(1, 0.95, 1.05), uCI = c(1.20, 1.15, 1.25),
                   chains = 2, burnin = 1000, iter = 4000, control = list(silent = TRUE), seed = 666,
                   priors_mu_null = NULL, priors_tau_null = NULL, priors_omega_null = NULL,
                   priors_omega = prior("two.sided", parameters = list(steps = c(.20), alpha = c(1,10)))))

fit_new2a <- RoBMA(y = d,  se = se, chains = 2, burnin = 1000, iter = 4000, control = list(silent = TRUE), seed = 666,
                   priors_mu_null = NULL, priors_tau_null = NULL, priors_omega_null = NULL,
                   priors_mu    = prior("normal",    parameters = list(mean = 1, sd = 1)),
                   priors_omega = prior("one.sided", parameters = list(steps = c(.20), alpha = c(1,1))))

fit_new2b <- RoBMA(y = -d, se = se, chains = 2, burnin = 1000, iter = 4000, control = list(silent = TRUE), seed = 666,
                   effect_direction = "negative", priors_mu_null = NULL, priors_tau_null = NULL, priors_omega_null = NULL,
                   priors_mu    = prior("normal",    parameters = list(mean = -1, sd = 1)),
                   priors_omega = prior("one.sided", parameters = list(steps = c(.20), alpha = c(1,1))))

test_that("OR model fit works", {
  fit_new1 <- remove_time(fit_new1)
  expect_equal(saved_fits2[[1]], fit_new1)
})

test_that("Direction change model fit works", {
  fit_new2a <- remove_time(fit_new2a)
  fit_new2b <- remove_time(fit_new2b)
  expect_equal(saved_fits2[[2]], fit_new2a)
  expect_equal(saved_fits2[[3]], fit_new2b)
})




#### creating / updating the test settings ####
if(FALSE){
  saved_fits <- list(fit_default, fit_custom1, fit_custom2, fit_custom3, fit_custom4, fit_custom5, fit_fail1, fit_fail2)
  for(i in 1:length(saved_fits)){
    saved_fits[[i]] <- remove_time(saved_fits[[i]])
  }
  saveRDS(saved_fits, file = "tests/results/saved_fits.RDS", compress  = "xz")

  updated_fits <- list(fit_update1a, fit_update2, fit_update3)
  for(i in 1:3){
    updated_fits[[i]] <- remove_time(updated_fits[[i]])
  }
  saveRDS(updated_fits, file = "tests/results/updated_fits.RDS", compress  = "xz")

  saved_fits2 <- list(fit_new1, fit_new2a, fit_new2b)
  for(i in 1:length(saved_fits2)){
    saved_fits2[[i]] <- remove_time(saved_fits2[[i]])
  }
  saveRDS(saved_fits2, file = "tests/results/saved_fits2.RDS", compress  = "xz")

}
