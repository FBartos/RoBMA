library(RoBMA)
x1 <- c(0, 5, 10,  5, 3)
x2 <- c(1, 4, 15, 10, 8)
n1 <- c(8, 15,30, 30,10)
n2 <- c(7, 14,32, 30,10)

nrep <- 5
fit <- BiBMA(x1 = rep(x1, nrep), x2 = rep(x2, nrep), n1 = rep(n1, nrep), n2 = rep(n2, nrep), silent = FALSE, seed = 1, parallel = TRUE)

# Robust Bayesian meta-analysis
# Components summary:
#   Models Prior prob. Post. prob. Inclusion BF
# Effect           2/4       0.500       1.000 4.928403e+09
# Heterogeneity    2/4       0.500       0.311 4.520000e-01
# Bias             0/4       0.000       0.000 0.000000e+00
#
# Model-averaged estimates:
#   Mean Median 0.025 0.975
# mu  0.769  0.767 0.563 0.981
# tau 0.040  0.000 0.000 0.249

summary(fit)
summary(fit, "i")
summary(fit, "m")

# Call:
#   BiBMA(x1 = rep(x1, nrep), x2 = rep(x2, nrep), n1 = rep(n1, nrep),
#         n2 = rep(n2, nrep), seed = 1, silent = FALSE)
#
# Robust Bayesian meta-analysis
# Models overview:
#   Model Prior Effect Prior Heterogeneity Prior Bias Prior prob. log(marglik) Post. prob. Inclusion BF
# 1     Spike(0)            Spike(0)                  0.250      -275.55       0.000        0.000
# 2     Spike(0)   InvGamma(1, 0.15)                  0.250      -271.98       0.000        0.000
# 3 Normal(0, 1)            Spike(0)                  0.250      -250.00       0.689        6.643
# 4 Normal(0, 1)   InvGamma(1, 0.15)                  0.250      -250.80       0.311        1.355

summary(fit$models[[4]]$fit)
fit$models[[4]]$fit$model
fit$models[[4]]$fit_summary
fit$models[[4]]$fit_summaries$logOR

debug(RoBMA:::.fit_BiBMA_model)
debug(BiBMA)

xxx$models[[1]]$fit$model


metafor::rma.glmm(measure = "OR",
                  ci = rep(x1, nrep), n2i = rep(n1, nrep),
                  ai = rep(x2, nrep), n1i = rep(n2, nrep), model = "UM.FS")
