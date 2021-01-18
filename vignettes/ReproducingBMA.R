## ----setup, include = FALSE---------------------------------------------------
is_check <- ("CheckExEnv" %in% search()) || any(c("_R_CHECK_TIMINGS_",
             "_R_CHECK_LICENSE_") %in% names(Sys.getenv())) || !file.exists("../prefitted/BMA_PowerPoseTest.RDS") 
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = !is_check
)

## ----include = FALSE----------------------------------------------------------
library(RoBMA)
# we pre-load the RoBMA models, the fitting time is around 2-5 minutes
fit_BMA_test   <- readRDS(file = "../prefitted/BMA_PowerPoseTest.RDS")
fit_BMA_est    <- readRDS(file = "../prefitted/BMA_PowerPoseEst.RDS")
fit_RoBMA_test <- readRDS(file = "../prefitted/PowerPoseTest.RDS")
fit_RoBMA_est  <- readRDS(file = "../prefitted/PowerPoseEst.RDS")

## -----------------------------------------------------------------------------
data("power_pose", package = "metaBMA")
power_pose[,c("study", "effectSize", "SE")]


## -----------------------------------------------------------------------------
fit_BMA_test$inclusion

round(fit_BMA_est$estimates,2)

## -----------------------------------------------------------------------------
summary(fit_RoBMA_test)

summary(fit_RoBMA_est, conditional = TRUE)

## ----fig.height = 3.25, fig.width = 4, fig.align = "center"-------------------
plot(fit_RoBMA_est, parameter = "mu", prior = TRUE)

## ----fig.height = 3.25, fig.width = 4, fig.align = "center"-------------------
plot(fit_RoBMA_est, parameter = "mu", prior = TRUE, type = "conditional", plot_type = "ggplot")

## ----fig.height = 5, fig.width = 7.5, fig.align = "center"--------------------
plot(fit_RoBMA_est, parameter = "mu", type = "models")

## ----fig.height = 4.5, fig.width = 5, fig.align = "center"--------------------
plot(fit_RoBMA_est, parameter = "forest")

