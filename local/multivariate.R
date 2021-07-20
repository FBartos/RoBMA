devtools::install_github("josue-rodriguez/bayeslincom")
devtools::install_github("donaldRwilliams/blsmeta")

# matrix multiplication
# inprod(X_location[i, ], beta)

library(psymetadata)
library(blsmeta)
fit <-  blsmeta(yi = yi,
                vi = vi,
                es_id = es_id,
                study_id = study_id,
                data = gnambs2020)

fit
fit$model_code
cat(fit$model_code)

data <- list(
  y         = gnambs2020$yi,
  v         = gnambs2020$vi,
  study_id  = gnambs2020$study_id,
  K         = nrow(gnambs2020),
  J         = length(unique(gnambs2020$study_id))
)

model1 <- "
model{

tau2 ~ dunif(0, 1)
tau3 ~ dunif(0, 1)
mu   ~ dnorm(0, 1)

for(k in 1:K){
  re_2[k] ~ dnorm(0, 1/pow(tau2, 2))
}

for(j in 1:J){
  re_3[j] ~ dnorm(0, 1/pow(tau3, 2))
}

for(k in 1:K){
  y[k] ~ dnorm(mu + re_2[k] + re_3[study_id[k]], 1 / v[k])
}

}
"

model2 <- "
model{

mu   ~ dnorm(0, 1)
tau2 ~ dunif(0, 1)
tau3 ~ dunif(0, 1)

for(j in 1:J){
  re_3[j] ~ dnorm(0, 1/pow(tau3, 2))
}

for(k in 1:K){
  y[k] ~ dnorm(mu + re_3[study_id[k]], 1 / (v[k] + pow(tau2, 2) ) )
}

}
"

model3 <- "
model{

mu   ~ dnorm(0, 1)
tau2 ~ dunif(0, 1)
tau3 ~ dunif(0, 1)

for(j in 1:J){
  tau3_z[j] ~ dchisqr(4)
}

for(k in 1:K){
  y[k] ~ dnorm(mu, 1 / (v[k] + pow(tau2, 2) + tau3_z[study_id[k]] * pow(tau3, 2) / 4 ) )
}

}
"

mu    <- .3
tau3  <- .5
tau2  <- .7
K     <- 100
n_rep <- 4

studies  <- rnorm(K, mu, tau3)
data     <- lapply(1:K, function(k) data.frame(study = k, true = rnorm(n_rep, studies[k], tau2)))
data     <- do.call(rbind, data)
data$se  <- runif(nrow(effects), 0.10, .30)
data$y   <- rnorm(nrow(data), data$true, data$se)
data$id  <- 1:nrow(data)

data_l <- list(
  y         = data$y,
  v         = data$se^2,
  study_id  = data$study,
  K         = nrow(data),
  J         = length(unique(data$study))
)

fit  <- blsmeta(yi = y, sei = se, es_id = id, study_id = study, data = data)
fit2 <- runjags::run.jags(model2, c("mu", "tau2", "tau3"), data_l, sample = 10000)
fit3 <- runjags::run.jags(model3, c("mu", "tau2", "tau3"), data_l, sample = 10000)



k_ID     <- do.call(c, lapply(1:K, function(k) rep(k, n_rep)))
obs_se   <- runif(length(k_ID), 0.10, .30)
obs_y    <- rnorm(K, mu, tau3)[k_ID] + rnorm(length(k_ID), 0, tau2) + rnorm(length(k_ID), 0, obs_se)
data_a   <- data.frame(
  y      = obs_y,
  se     = obs_se,
  study  = k_ID,
  id     = 1:length(k_ID)
)
fit_a  <- blsmeta(yi = y, sei = se, es_id = id, study_id = study, data = data_a)


k_ID     <- do.call(c, lapply(1:K, function(k) rep(k, sample(1:10, 1))))
obs_se   <- runif(length(k_ID), 0.10, .30)
obs_y    <- rnorm(K, mu, tau3)[k_ID] + rnorm(length(k_ID), 0, sqrt(tau2^2 + obs_se^2))
data_b   <- data.frame(
  y      = obs_y,
  se     = obs_se,
  study  = k_ID,
  id     = 1:length(k_ID)
)
fit_b  <- blsmeta(yi = y, sei = se, es_id = id, study_id = study, data = data_b)


k_ID     <- do.call(c, lapply(1:K, function(k) rep(k, sample(1:10, 1))))
obs_se   <- runif(length(k_ID), 0.10, .30)
obs_y    <- rnorm(K, mu, tau3)[k_ID] + rnorm(length(k_ID), 0, sqrt(tau2^2 + obs_se^2))
data_c   <- data.frame(
  y      = obs_y,
  se     = obs_se,
  study  = k_ID,
  id     = 1:length(k_ID)
)
fit_c  <- blsmeta(yi = y, sei = se, es_id = id, study_id = study, data = data_c)









omega1      <- 0.15
mu          <- 0.10
K           <- 10000
k_ID        <- do.call(c, lapply(1:K, function(k) rep(k, sample(1:10, 1))))
obs_se      <- runif(length(k_ID), 0.10, .30)
true_study  <- mu + rnorm(K, 0, tau3)[k_ID]
true_effect <- true_study  + rnorm(length(k_ID), 0, tau2)
obs_effect  <- true_effect + rnorm(length(k_ID), 0, obs_se)


sel_effect  <- rbinom(length(k_ID), 1, ifelse(obs_effect / obs_se > 1.96, 1, omega1))

hist(true_study)
hist(true_study[as.logical(sel_effect)])


plot(density(true_study), col = "blue", lwd = 2)
lines(density(true_study[as.logical(sel_effect)]), col = "red", lwd = 2)

plot(density(true_effect), col = "blue", lwd = 2)
lines(density(true_effect[as.logical(sel_effect)]), col = "red", lwd = 2)

plot(density(obs_effect), col = "blue", lwd = 2)
lines(density(obs_effect[as.logical(sel_effect)]), col = "red", lwd = 2)



k_ID        <- do.call(c, lapply(1:K, function(k) rep(k, sample(1:10, 1))))
obs_se      <- runif(length(k_ID), 0.10, .30)
dev_study   <- rnorm(K, 0, tau3)[k_ID]
dev_effect  <- rnorm(length(k_ID), 0, tau2)
dev_error   <- rnorm(length(k_ID), 0, obs_se)

obs_effect  <- mu + dev_study + dev_effect + dev_error

