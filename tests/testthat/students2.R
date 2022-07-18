library(RoBMA)
data("Kroupova2021", package = "RoBMA")
Kroupova2021    <- Kroupova2021[!is.na(Kroupova2021$se),]

fit_fixed <- RoBMA.reg(
  # specify the model formula and data input
  formula   = ~ 1,
  data      = Kroupova2021,

  # specify slightly informative prior for the effect size parameter (on Cohens'd scale)
  priors_effect = prior("normal", parameters = list(mean = 0, sd = 1)),
  prior_scale   = "cohens_d",

  # remove the remaining model components
  priors_bias          = NULL,
  priors_heterogeneity = NULL,

  # some additional settings
  parallel = TRUE, seed = 1
)
fit_fixed_weighted <- RoBMA.reg(
  # specify the model formula and data input
  formula   = ~ 1,
  data      = Kroupova2021,
  study_ids = Kroupova2021$study,

  # specify slightly informative prior for the effect size parameter (on Cohens'd scale)
  priors_effect = prior("normal", parameters = list(mean = 0, sd = 1)),
  prior_scale   = "cohens_d",

  # remove the remaining model components
  priors_bias          = NULL,
  priors_heterogeneity = NULL,

  # some additional settings
  parallel = TRUE, seed = 1, weighted = TRUE
)
fit_random_weighted <- RoBMA.reg(
  # specify the model formula and data input
  formula   = ~ 1,
  data      = Kroupova2021,
  study_ids = Kroupova2021$study,

  # specify slightly informative prior for the effect size parameter (on Cohens'd scale)
  priors_effect = prior("normal", parameters = list(mean = 0, sd = 1)),
  priors_heterogeneity = prior(distribution = "invgamma", parameters = list(shape = 1, scale = 0.15)),
  prior_scale   = "cohens_d",

  # remove the remaining model components
  priors_bias               = NULL,
  priors_heterogeneity_null = NULL,

  # some additional settings
  parallel = TRUE, seed = 1, weighted = TRUE
)
Kroupova2021$z_se  <- se_r2se_z(se_r = Kroupova2021$se, r = Kroupova2021$r)
fit_BMA_weighted <- RoBMA.reg(
  # specify the model formula and data input
  formula   = ~ 1,
  data      = Kroupova2021,
  study_ids = Kroupova2021$study,

  # specify slightly informative prior for the effect size parameter (on Cohens'd scale)
  priors_effect = prior("normal", parameters = list(mean = 0, sd = 1)),
  priors_heterogeneity = prior(distribution = "invgamma", parameters = list(shape = 1, scale = 0.15)),
  prior_scale   = "cohens_d",

  # remove the remaining model components
  priors_bias = NULL,

  # some additional settings
  parallel = TRUE, seed = 1, weighted = TRUE
)
fit_BMA_regression_weighted <- RoBMA.reg(
  # specify the model formula and data input
  formula   = ~ 1 + z_se,
  data      = Kroupova2021,
  study_ids = Kroupova2021$study,
  test_predictors        = "z_se",
  standardize_predictors = FALSE,

  # specify slightly informative prior for the effect size parameter (on Cohens'd scale)
  priors = list("z_se" = prior("normal", parameters = list(mean = 0, sd = 1))),
  priors_effect = prior("normal", parameters = list(mean = 0, sd = 1)),
  priors_heterogeneity = prior(distribution = "invgamma", parameters = list(shape = 1, scale = 0.15)),
  prior_scale   = "cohens_d",

  # remove the remaining model components
  priors_bias = NULL,

  # some additional settings
  parallel = TRUE, seed = 1, weighted = TRUE
)
fit_RoBMA <- RoBMA.reg(
  # specify the model formula and data input
  formula    = ~ 1,
  data       = Kroupova2021,
  model_type       = "PSMA",
  effect_direction = "negative",

  # specify slightly informative prior for the effect size parameter (on Cohens'd scale)
  # some additional settings
  parallel = TRUE, seed = 1
)
fit_RoBMA_weighted <- RoBMA.reg(
  # specify the model formula and data input
  formula    = ~ 1,
  data       = Kroupova2021,
  study_ids  = Kroupova2021$study,
  model_type       = "PSMA",
  effect_direction = "negative",

  # specify slightly informative prior for the effect size parameter (on Cohens'd scale)
  # some additional settings
  parallel = TRUE, seed = 1, weighted = TRUE
)
fit_RoBMA_regression_full_weighted <- RoBMA.reg(
  # specify the model formula and data input
  formula    = ~ 1 + education_outcome + students_gender + location + design + endogenity_control,
  data       = Kroupova2021,
  study_ids  = Kroupova2021$study,
  test_predictors = "",

  model_type       = "PSMA",
  effect_direction = "negative",
  priors     = list(
    education_outcome  = prior_factor("mnormal", list(mean = 0, sd = 0.5), contrast = "orthonormal"),
    students_gender    = prior_factor("mnormal", list(mean = 0, sd = 0.5), contrast = "orthonormal"),
    location           = prior_factor("mnormal", list(mean = 0, sd = 0.5), contrast = "orthonormal"),
    design             = prior_factor("mnormal", list(mean = 0, sd = 0.5), contrast = "orthonormal"),
    endogenity_control = prior_factor("mnormal", list(mean = 0, sd = 0.5), contrast = "orthonormal")
  ),

  # some additional settings
  parallel = TRUE, seed = 1, weighted = TRUE
)

debug(RoBMA:::.fit_RoBMA_model)

apply(Kroupova2021, 2, anyNA)



saveRDS(fit_fixed,                 file = "../models/MetaRegression/fit_fixed.RDS", compress = "xz")
saveRDS(fit_fixed_weighted,        file = "../models/MetaRegression/fit_fixed_weighted.RDS", compress = "xz")
saveRDS(fit_random_weighted,                  file = "../models/MetaRegression/fit_random_weighted.RDS", compress = "xz")
saveRDS(fit_BMA_weighted,                  file = "../models/MetaRegression/fit_BMA_weighted.RDS", compress = "xz")
saveRDS(fit_BMA_regression_weighted,                  file = "../models/MetaRegression/fit_BMA_regression_weighted.RDS", compress = "xz")
saveRDS(fit_RoBMA,                  file = "../models/MetaRegression/fit_RoBMA.RDS", compress = "xz")
saveRDS(fit_RoBMA_weighted,                  file = "../models/MetaRegression/fit_RoBMA_weighted.RDS", compress = "xz")



library(RoBMA)
data("Kroupova2021", package = "RoBMA")
Kroupova2021    <- Kroupova2021[!is.na(Kroupova2021$se),]

job::job({


  missing <- list.files("../models/MetaRegression/fit_RoBMA_regression_full_weighted/")
  missing <- gsub("m_", "", missing)
  missing <- gsub(".RDS", "", missing)

  fit_RoBMA_regression_full_weighted <- RoBMA.reg(
    # specify the model formula and data input
    formula    = ~ 1 + education_outcome + students_gender + location + design + endogenity_control,
    data       = Kroupova2021,
    study_ids  = Kroupova2021$study,
    test_predictors = "",

    model_type       = "PSMA",
    effect_direction = "negative",
    priors     = list(
      education_outcome  = prior_factor("mnormal", list(mean = 0, sd = 0.5), contrast = "orthonormal"),
      students_gender    = prior_factor("mnormal", list(mean = 0, sd = 0.5), contrast = "orthonormal"),
      location           = prior_factor("mnormal", list(mean = 0, sd = 0.5), contrast = "orthonormal"),
      design             = prior_factor("mnormal", list(mean = 0, sd = 0.5), contrast = "orthonormal"),
      endogenity_control = prior_factor("mnormal", list(mean = 0, sd = 0.5), contrast = "orthonormal")
    ),

    # some additional settings
    parallel = FALSE, seed = 1, weighted = TRUE, do_not_fit = TRUE
  )

  missing <- seq_along(fit_RoBMA_regression_full_weighted$models)[!seq_along(fit_RoBMA_regression_full_weighted$models) %in% as.numeric(missing)]

  for(i in missing){
    temp_model <- RoBMA:::.fit_RoBMA_model(fit_RoBMA_regression_full_weighted, i)
    saveRDS(temp_model, file = paste0("../models/MetaRegression/fit_RoBMA_regression_full_weighted/", "m_", i, ".RDS"), compress = "xz")
  }

})



fit_RoBMA_regression_weighted <- RoBMA.reg(
  # specify the model formula and data input
  formula    = ~ 1 + education_outcome + students_gender + location + design + endogenity_control,
  data       = Kroupova2021,
  study_ids  = Kroupova2021$study,
  test_predictors = c("education_outcome", "students_gender", "location", "design", "endogenity_control"),

  model_type       = "PSMA",
  effect_direction = "negative",
  priors     = list(
    education_outcome  = prior_factor("mnormal", list(mean = 0, sd = 0.5), contrast = "orthonormal"),
    students_gender    = prior_factor("mnormal", list(mean = 0, sd = 0.5), contrast = "orthonormal"),
    location           = prior_factor("mnormal", list(mean = 0, sd = 0.5), contrast = "orthonormal"),
    design             = prior_factor("mnormal", list(mean = 0, sd = 0.5), contrast = "orthonormal"),
    endogenity_control = prior_factor("mnormal", list(mean = 0, sd = 0.5), contrast = "orthonormal")
  ),

  # some additional settings
  parallel = FALSE, seed = 1, weighted = TRUE, do_not_fit = TRUE
)

missing <- list.files("../models/MetaRegression/fit_RoBMA_regression_weighted/")
missing <- gsub("m_", "", missing)
missing <- gsub(".RDS", "", missing)
missing <- seq_along(fit_RoBMA_regression_weighted$models)[!seq_along(fit_RoBMA_regression_weighted$models) %in% as.numeric(missing)]

cl <- parallel::makeCluster(23)
parallel::clusterExport(cl, c("fit_RoBMA_regression_weighted"))
parallel::parSapplyLB(cl, missing, function(i){

  library(RoBMA)
  temp_model <- RoBMA:::.fit_RoBMA_model(fit_RoBMA_regression_weighted, i)
  saveRDS(temp_model, file = paste0("../models/MetaRegression/fit_RoBMA_regression_weighted/", "m_", i, ".RDS"), compress = "xz")

})

parallel::stopCluster(cl)











for(i in seq_along(fit_RoBMA_regression_full_weighted$models)){
  fit_RoBMA_regression_full_weighted$models[[i]] <- readRDS(paste0("../models/MetaRegression/fit_RoBMA_regression_full_weighted/", "m_", i, ".RDS"))
}
object <- fit_RoBMA_regression_full_weighted

sum(!RoBMA:::.get_model_convergence(object))

object$models        <- BayesTools::models_inference(object[["models"]])
object$RoBMA         <- RoBMA:::.ensemble_inference(object)
object$coefficients  <- RoBMA:::.compute_coeficients(object[["RoBMA"]])

### collect and print errors and warnings
object$add_info[["errors"]]   <- c(object$add_info[["errors"]],   RoBMA:::.get_model_errors(object))
object$add_info[["warnings"]] <- c(object$add_info[["warnings"]], RoBMA:::.get_model_warnings(object))
RoBMA:::.print_errors_and_warnings(object)

class(object) <- c("RoBMA", "RoBMA.reg")

fit_RoBMA_regression_full_weighted  <- object
saveRDS(fit_RoBMA_regression_full_weighted , file = "../models/MetaRegression/fit_RoBMA_regression_full_weighted .RDS", compress = "xz")


object <- .remove_model_posteriors(object)
object <- .remove_model_margliks(object)
saveRDS(object, file = "../models/MetaRegression/fit_RoBMA_regression_full_weighted .RDS", compress = "xz")




# fit_RoBMA_regression_weighted
for(i in seq_along(fit_RoBMA_regression_weighted$models)){
  fit_RoBMA_regression_weighted$models[[i]] <- readRDS(paste0("../models/MetaRegression/fit_RoBMA_regression_weighted/", "m_", i, ".RDS"))
}
object <- fit_RoBMA_regression_weighted

sum(!RoBMA:::.get_model_convergence(object))

object$models        <- BayesTools::models_inference(object[["models"]])
object$RoBMA         <- RoBMA:::.ensemble_inference(object)
object$coefficients  <- RoBMA:::.compute_coeficients(object[["RoBMA"]])

### collect and print errors and warnings
object$add_info[["errors"]]   <- c(object$add_info[["errors"]],   RoBMA:::.get_model_errors(object))
object$add_info[["warnings"]] <- c(object$add_info[["warnings"]], RoBMA:::.get_model_warnings(object))
RoBMA:::.print_errors_and_warnings(object)

class(object) <- c("RoBMA", "RoBMA.reg")

fit_RoBMA_regression_weighted <- object
saveRDS(fit_RoBMA_regression_weighted, file = "../models/MetaRegression/fit_RoBMA_regression_weighted.RDS", compress = "xz")


object <- .remove_model_posteriors(object)
object <- .remove_model_margliks(object)
saveRDS(object, file = "../models/MetaRegression/fit_RoBMA_regression_weighted.RDS", compress = "xz")
