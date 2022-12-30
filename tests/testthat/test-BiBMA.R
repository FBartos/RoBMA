library(RoBMA)
x1 <- c(0, 5, 10,  5, 3)
x2 <- c(1, 4, 15, 10, 8)
n1 <- c(8, 15,30, 30,10)
n2 <- c(7, 14,32, 30,10)


fit <- BiBMA(x1 = rep(x1, 10), x2 = rep(x2, 10), n1 = rep(n1, 10), n2 = rep(n2, 10), silent = FALSE, seed = 1)

summary(fit)
summary(fit, "i")

summary(fit$models[[4]]$fit)
fit$models[[4]]$fit_summary

debug(RoBMA:::.fit_BiBMA_model)
debug(BiBMA)

xxx$models[[1]]$fit$model


metafor::rma.glmm(measure = "OR",
                  ci = rep(x1, 10), n2i = rep(n1, 10),
                  ai = rep(x2, 10), n1i = rep(n2, 10), model = "UM.FS")
