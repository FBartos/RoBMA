library(RoBMA)
module_location <- RoBMA:::RoBMA.private$module_location
load_all()
RoBMA:::RoBMA.private$module_location <- module_location
rjags::load.module("RoBMA", path = module_location)

undebug(BayesTools::JAGS_bridgesampling)

m19 <- .fit_RoBMA_model(object, 19)
m55 <- .fit_RoBMA_model(object, 55)


m19$fit_summary
m55$fit_summary

m19$marglik
m55$marglik

library(RoBMA)
data <- cbind.data.frame(
  d        = c(0.25, 0.15, 0.20, 0.40, 0.42, 0.45, 0.31, 0.32, 0.30),
  se       = runif(9, 0.10, 0.20),
  "x_cont" = rnorm(9),
  "x_fac3" = c("A", "A", "A", "B", "B", "B", "C", "C", "C")
)
data <- rbind(data, data, data)
fit <- RoBMA.reg(~ x_cont + x_fac3, data = data, test_predictors = "x_fac3",
                 prior_factors = prior_factor("mnormal", parameters = list(mean = 0, sd = 0.25), contrast = "orthonormal"),
                 chains = 2, sample = 2000, parallel = TRUE)
summary(fit)

summary(fit, "m")
summary(fit, "i")
saveRDS(fit,file = "tests/testthat/fit.RDS")

fit <- readRDS(file = "tests/testthat/fit.RDS")

### summary settings ----
object <- fit
type = "ensemble"; conditional = FALSE;
output_scale = NULL; probs = c(.025, .975); logBF = FALSE; BF01 = FALSE;
short_name = FALSE; remove_spike_0 = FALSE;

summary(fit, "m", short_name = T)


### RoBMA.reg settings ----
formula = ~ x_cont + x_fac3
data    = data
test_predictors = "x_fac3"
study_names = data$study
study_ids = NULL
transformation     = "fishers_z"
prior_scale        = "cohens_d"
standardize_predictors = TRUE
effect_direction       = "positive"

# prior specification
priors       = NULL
model_type   = NULL

priors_effect         = prior(distribution = "normal",    parameters = list(mean  = 0, sd = 1))
priors_heterogeneity  = prior(distribution = "invgamma",  parameters = list(shape = 1, scale = .15))
priors_bias           = list(
  prior_weightfunction(distribution = "two.sided", parameters = list(alpha = c(1, 1),       steps = c(0.05)),             prior_weights = 1/12),
  prior_weightfunction(distribution = "two.sided", parameters = list(alpha = c(1, 1, 1),    steps = c(0.05, 0.10)),       prior_weights = 1/12),
  prior_weightfunction(distribution = "one.sided", parameters = list(alpha = c(1, 1),       steps = c(0.05)),             prior_weights = 1/12),
  prior_weightfunction(distribution = "one.sided", parameters = list(alpha = c(1, 1, 1),    steps = c(0.025, 0.05)),      prior_weights = 1/12),
  prior_weightfunction(distribution = "one.sided", parameters = list(alpha = c(1, 1, 1),    steps = c(0.05, 0.5)),        prior_weights = 1/12),
  prior_weightfunction(distribution = "one.sided", parameters = list(alpha = c(1, 1, 1, 1), steps = c(0.025, 0.05, 0.5)), prior_weights = 1/12),
  prior_PET(distribution   = "Cauchy", parameters = list(0,1), truncation = list(0, Inf),  prior_weights = 1/4),
  prior_PEESE(distribution = "Cauchy", parameters = list(0,5), truncation = list(0, Inf),  prior_weights = 1/4)
)
priors_effect_null         = prior(distribution = "point", parameters = list(location = 0))
priors_heterogeneity_null  = prior(distribution = "point", parameters = list(location = 0))
priors_bias_null           = prior_none()
priors_rho                 = prior("beta", parameters = list(alpha = 1, beta = 1))
priors_rho_null            = NULL

prior_covariates       = prior("normal", parameters = list(mean = 0, sd = 0.5))
prior_covariates_null  = prior("spike",  parameters = list(location = 0))
prior_factors          = prior_factor("mnormal", parameters = list(mean = 0, sd = 0.50), contrast = "orthonormal")
prior_factors_null     = prior("spike",  parameters = list(location = 0))

chains = 1; sample = 2000; burnin = 1000; adapt = 500; thin = 1; parallel = FALSE;
autofit = FALSE; autofit_control = set_autofit_control(); convergence_checks = set_convergence_checks(max_Rhat = 999, min_ESS = 1, max_error = 999, max_SD_error = 1);

# additional settings
save = "all"; seed = NULL; silent = FALSE

dots         <- .RoBMA_collect_dots()
