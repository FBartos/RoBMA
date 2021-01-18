## ----setup, include = FALSE---------------------------------------------------
is_check <- ("CheckExEnv" %in% search()) || any(c("_R_CHECK_TIMINGS_",
             "_R_CHECK_LICENSE_") %in% names(Sys.getenv())) || !file.exists("../prefitted/Bem_update1.RDS") 
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = !is_check
)

## -----------------------------------------------------------------------------
Bem2011 <- data.frame(
  study = c( "1",  "2",  "3",  "4",  "5",  "6",  "7",  "8",  "9"),
  t     = c(2.51, 2.39, 2.55, 2.03, 2.23, 2.41, 1.31, 1.92, 2.96),
  N     = c( 100,  150,   97,   99,  100,  150,  200,  100,   50)
)

## -----------------------------------------------------------------------------
library(RoBMA)

fit <- RoBMA(t = Bem2011$t, n = Bem2011$N, study_names = Bem2011$study,
             priors_mu = NULL, priors_tau = NULL, priors_omega = NULL,
             priors_mu_null    = prior("spike", parameters = list(location = 0)),
             priors_tau_null   = prior("spike", parameters = list(location = 0)),
             priors_omega_null = prior("spike", parameters = list(location = 1)),
             control = list(silent = TRUE), seed = 666)

## ----fig.height = 3.25, fig.width = 4, fig.align = "center"-------------------
plot(prior("normal", parameters = list(mean = .15, sd = .10)))

## ----include = FALSE----------------------------------------------------------
# these fits are relatively fast, but we reduce the knitting time considerably
fit <- readRDS(file = "../prefitted/Bem_update1.RDS")

## -----------------------------------------------------------------------------
summary(fit, type = "models")

summary(fit, type = "individual")

## ----fig.height = 3.25, fig.width = 4, fig.align = "center"-------------------
plot(prior("one.sided", parameters = list(steps = c(0.05, .10), alpha = c(1,1,1))))

## ----include = FALSE----------------------------------------------------------
fit <- readRDS(file = "../prefitted/Bem_update2.RDS")

## -----------------------------------------------------------------------------
summary(fit, type = "models")

## -----------------------------------------------------------------------------
summary(fit)

## ----fig.height = 3.25, fig.width = 4, fig.align = "center"-------------------
plot(fit, parameter = "mu", prior = TRUE)

## ----fig.height = 3.25, fig.width = 4, fig.align = "center"-------------------
plot(fit, parameter = "omega", prior = TRUE)

