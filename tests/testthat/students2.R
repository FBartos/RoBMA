library(RoBMA)
data("Kroupova2021", package = "RoBMA")
Kroupova2021    <- Kroupova2021[!is.na(Kroupova2021$se),]


fit_BFE <- RoBMA.reg(
  # specify the model formula and data input
  formula   = ~ 1,
  data      = Kroupova2021,

  # specify slightly informative prior for the effect size parameter (on Fisher's z scale)
  priors_effect = prior("normal", parameters = list(mean = 0, sd = 1)),
  prior_scale   = "fishers_z",

  # remove the remaining model components
  priors_bias          = NULL,
  priors_heterogeneity = NULL,
  priors_effect_null   = NULL,

  # some additional settings
  parallel = TRUE, seed = 1
)
fit_wBFE <- RoBMA.reg(
  # specify the model formula and data input
  formula   = ~ 1,
  data      = Kroupova2021,
  study_ids = Kroupova2021$study,
  weighted  = TRUE,

  # specify slightly informative prior for the effect size parameter (on Fisher's z scale)
  priors_effect = prior("normal", parameters = list(mean = 0, sd = 1)),
  prior_scale   = "fishers_z",

  # remove the remaining model components
  priors_bias          = NULL,
  priors_heterogeneity = NULL,
  priors_effect_null   = NULL,

  # some additional settings
  parallel = TRUE, seed = 1
)
fit_BFE10 <- RoBMA.reg(
  # specify the model formula and data input
  formula   = ~ 1,
  data      = Kroupova2021,

  # specify informative prior for the effect size parameter under the alternative hypothesis
  # and a specify a null hypothesis of no effect
  priors_effect       = prior("normal", parameters = list(mean = 0, sd = 0.25)),
  priors_effect_null  = prior("spike",  parameters = list(location = 0)),
  prior_scale         = "fishers_z",

  # remove the remaining model components
  priors_bias          = NULL,
  priors_heterogeneity = NULL,

  # some additional settings
  parallel = TRUE, seed = 1
)
fit_wBFE10 <- RoBMA.reg(
  # specify the model formula and data input
  formula   = ~ 1,
  data      = Kroupova2021,
  study_ids = Kroupova2021$study,
  weighted  = TRUE,

  # specify informative prior for the effect size parameter under the alternative hypothesis
  # and a specify a null hypothesis of no effect
  priors_effect       = prior("normal", parameters = list(mean = 0, sd = 0.25)),
  priors_effect_null  = prior("spike",  parameters = list(location = 0)),
  prior_scale         = "fishers_z",

  # remove the remaining model components
  priors_bias          = NULL,
  priors_heterogeneity = NULL,

  # some additional settings
  parallel = TRUE, seed = 1
)
fit_wPSMA <- RoBMA.reg(
  # specify the model formula and data input
  formula          = ~ 1,
  data             = Kroupova2021,
  effect_direction = "negative",
  study_ids        = Kroupova2021$study,
  weighted         = TRUE,

  # specify informative prior for the effect size parameter under the alternative hypothesis
  # and a specify a null hypothesis of no effect
  priors_effect       = prior("normal", parameters = list(mean = 0, sd = 0.25)),
  priors_effect_null  = prior("spike",  parameters = list(location = 0)),
  prior_scale         = "fishers_z",

  # some additional settings
  parallel = TRUE, seed = 1
)
fit_wBFE_reg <- RoBMA.reg(
  # specify the model formula and data input
  formula         = ~ 1 + location,
  data            = Kroupova2021,
  study_ids       = Kroupova2021$study,
  test_predictors = "",
  weighted        = TRUE,

  # specify slightly informative prior for the effect size parameter (on Fisher's z scale)
  priors       = list(
    location = prior_factor("mnormal", list(mean = 0, sd = 0.50), contrast = "orthonormal")
  ),
  priors_effect = prior("normal", parameters = list(mean = 0, sd = 1)),
  prior_scale   = "fishers_z",

  # remove the remaining model components
  priors_bias          = NULL,
  priors_heterogeneity = NULL,
  priors_effect_null   = NULL,

  # some additional settings
  parallel = TRUE, seed = 1
)
fit_wBFE_reg10 <- RoBMA.reg(
  # specify the model formula and data input
  formula         = ~ 1 + location,
  data            = Kroupova2021,
  study_ids       = Kroupova2021$study,
  test_predictors = "location",
  weighted        = TRUE,

  # specify slightly informative prior for the effect size parameter (on Fisher's z scale)
  priors       = list(
    location = prior_factor("mnormal", list(mean = 0, sd = 0.10), contrast = "orthonormal")
  ),
  priors_effect = prior("normal", parameters = list(mean = 0, sd = 0.25)),
  prior_scale   = "fishers_z",

  # remove the remaining model components
  priors_bias          = NULL,
  priors_heterogeneity = NULL,
  priors_effect_null   = NULL,

  # some additional settings
  parallel = TRUE, seed = 1
)
fit_wPSMA_reg <- RoBMA.reg(
  # specify the model formula and data input
  formula          = ~ 1 + education_outcome + students_gender + location + design + endogenity_control,
  data             = Kroupova2021,
  effect_direction = "negative",
  study_ids        = Kroupova2021$study,
  weighted         = TRUE,
  test_predictors  = "",

  # specify informative prior for the effect size parameter under the alternative hypothesis
  # and a specify a null hypothesis of no effect
  priors             = list(
    education_outcome  = prior_factor("mnormal", list(mean = 0, sd = 0.10), contrast = "orthonormal"),
    students_gender    = prior_factor("mnormal", list(mean = 0, sd = 0.10), contrast = "orthonormal"),
    location           = prior_factor("mnormal", list(mean = 0, sd = 0.10), contrast = "orthonormal"),
    design             = prior_factor("mnormal", list(mean = 0, sd = 0.10), contrast = "orthonormal"),
    endogenity_control = prior_factor("mnormal", list(mean = 0, sd = 0.10), contrast = "orthonormal")
  ),
  priors_effect       = prior("normal", parameters = list(mean = 0, sd = 0.25)),
  priors_effect_null  = prior("spike",  parameters = list(location = 0)),
  prior_scale         = "fishers_z",

  # some additional settings
  parallel = TRUE, seed = 1
)

fit_wPSMA_reg10 <- RoBMA.reg(
  # specify the model formula and data input
  formula          = ~ 1 + education_outcome + students_gender + location + design + endogenity_control,
  data             = Kroupova2021,
  effect_direction = "negative",
  study_ids        = Kroupova2021$study,
  weighted         = TRUE,
  test_predictors  = c("education_outcome", "students_gender", "location", "design", "endogenity_control"),

  # specify informative prior for the effect size parameter under the alternative hypothesis
  # and a specify a null hypothesis of no effect
  priors             = list(
    education_outcome  = prior_factor("mnormal", list(mean = 0, sd = 0.10), contrast = "orthonormal"),
    students_gender    = prior_factor("mnormal", list(mean = 0, sd = 0.10), contrast = "orthonormal"),
    location           = prior_factor("mnormal", list(mean = 0, sd = 0.10), contrast = "orthonormal"),
    design             = prior_factor("mnormal", list(mean = 0, sd = 0.10), contrast = "orthonormal"),
    endogenity_control = prior_factor("mnormal", list(mean = 0, sd = 0.10), contrast = "orthonormal")
  ),
  priors_effect       = prior("normal", parameters = list(mean = 0, sd = 0.25)),
  priors_effect_null  = prior("spike",  parameters = list(location = 0)),
  prior_scale         = "fishers_z",

  # some additional settings
  parallel = FALSE, seed = 1, do_not_fit = TRUE
)

missing <- list.files("../models/MetaRegression/fit_wPSMA_reg10/")
missing <- gsub("m_", "", missing)
missing <- gsub(".RDS", "", missing)
missing <- seq_along(fit_wPSMA_reg10$models)[!seq_along(fit_wPSMA_reg10$models) %in% as.numeric(missing)]

cl <- parallel::makeCluster(23)
parallel::clusterExport(cl, c("fit_wPSMA_reg10"))
parallel::parSapplyLB(cl, missing, function(i){

  library(RoBMA)
  temp_model <- RoBMA:::.fit_RoBMA_model(fit_wPSMA_reg10, i)
  saveRDS(temp_model, file = paste0("../models/MetaRegression/fit_wPSMA_reg10/", "m_", i, ".RDS"), compress = "xz")

})

parallel::stopCluster(cl)

fit_wPSMA_reg10 <- RoBMA.reg(
  # specify the model formula and data input
  formula          = ~ 1 + education_outcome + students_gender + location + design + endogenity_control,
  data             = Kroupova2021,
  effect_direction = "negative",
  study_ids        = Kroupova2021$study,
  weighted         = TRUE,
  test_predictors  = c("education_outcome", "students_gender", "location", "design", "endogenity_control"),

  # specify informative prior for the effect size parameter under the alternative hypothesis
  # and a specify a null hypothesis of no effect
  priors             = list(
    education_outcome  = prior_factor("mnormal", list(mean = 0, sd = 0.10), contrast = "orthonormal"),
    students_gender    = prior_factor("mnormal", list(mean = 0, sd = 0.10), contrast = "orthonormal"),
    location           = prior_factor("mnormal", list(mean = 0, sd = 0.10), contrast = "orthonormal"),
    design             = prior_factor("mnormal", list(mean = 0, sd = 0.10), contrast = "orthonormal"),
    endogenity_control = prior_factor("mnormal", list(mean = 0, sd = 0.10), contrast = "orthonormal")
  ),
  priors_effect       = prior("normal", parameters = list(mean = 0, sd = 0.25)),
  priors_effect_null  = prior("spike",  parameters = list(location = 0)),
  prior_scale         = "fishers_z",

  # some additional settings
  parallel = FALSE, seed = 1, do_not_fit = TRUE
)

for(i in seq_along(fit_PSMA_reg10$models)){
  fit_PSMA_reg10$models[[i]] <- readRDS(paste0("../models/MetaRegression/fit_PSMA_reg10/", "m_", i, ".RDS"))
}

fit_PSMA_reg10$models        <- BayesTools::models_inference(fit_PSMA_reg10[["models"]])
fit_PSMA_reg10$RoBMA         <- RoBMA:::.ensemble_inference(fit_PSMA_reg10)
fit_PSMA_reg10$coefficients  <- RoBMA:::.compute_coeficients(fit_PSMA_reg10[["RoBMA"]])

fit_PSMA_reg10$add_info[["errors"]]   <- c(fit_PSMA_reg10$add_info[["errors"]],   RoBMA:::.get_model_errors(fit_PSMA_reg10))
fit_PSMA_reg10$add_info[["warnings"]] <- c(fit_PSMA_reg10$add_info[["warnings"]], RoBMA:::.get_model_warnings(fit_PSMA_reg10))

fit_PSMA_reg10 <- RoBMA:::.remove_model_posteriors(fit_PSMA_reg10)
fit_PSMA_reg10 <- RoBMA:::.remove_model_margliks(fit_PSMA_reg10)

class(fit_PSMA_reg10) <- c("RoBMA", "RoBMA.reg")

fit_wPSMA <- RoBMA:::.remove_model_posteriors(fit_wPSMA)
fit_wPSMA <- RoBMA:::.remove_model_margliks(fit_wPSMA)

fit_wPSMA_reg <- RoBMA:::.remove_model_posteriors(fit_wPSMA_reg)
fit_wPSMA_reg <- RoBMA:::.remove_model_margliks(fit_wPSMA_reg)

saveRDS(fit_BFE,  file = "../models/MetaRegression/fit_BFE.RDS",  compress = "xz")
saveRDS(fit_wBFE, file = "../models/MetaRegression/fit_wBFE.RDS", compress = "xz")
saveRDS(fit_BFE10,  file = "../models/MetaRegression/fit_BFE10.RDS",  compress = "xz")
saveRDS(fit_wBFE10, file = "../models/MetaRegression/fit_wBFE10.RDS", compress = "xz")
saveRDS(fit_wPSMA,       file = "../models/MetaRegression/fit_wPSMA.RDS",       compress = "xz")
saveRDS(fit_wBFE_reg ,   file = "../models/MetaRegression/fit_wBFE_reg.RDS",    compress = "xz")
saveRDS(fit_wBFE_reg10 , file = "../models/MetaRegression/fit_wBFE_reg10.RDS",  compress = "xz")
saveRDS(fit_wPSMA_reg,   file = "../models/MetaRegression/fit_wPSMA_reg.RDS",   compress = "xz")
saveRDS(fit_PSMA_reg10, file = "../models/MetaRegression/fit_PSMA_reg10.RDS", compress = "xz")




fit_PSMA_reg10 <- RoBMA.reg(
  # specify the model formula and data input
  formula          = ~ 1 + education_outcome + students_gender + location + design + endogenity_control,
  data             = Kroupova2021,
  effect_direction = "negative",
  test_predictors  = c("education_outcome", "students_gender", "location", "design", "endogenity_control"),

  # specify informative prior for the effect size parameter under the alternative hypothesis
  # and a specify a null hypothesis of no effect
  priors             = list(
    education_outcome  = prior_factor("mnormal", list(mean = 0, sd = 0.10), contrast = "orthonormal"),
    students_gender    = prior_factor("mnormal", list(mean = 0, sd = 0.10), contrast = "orthonormal"),
    location           = prior_factor("mnormal", list(mean = 0, sd = 0.10), contrast = "orthonormal"),
    design             = prior_factor("mnormal", list(mean = 0, sd = 0.10), contrast = "orthonormal"),
    endogenity_control = prior_factor("mnormal", list(mean = 0, sd = 0.10), contrast = "orthonormal")
  ),
  priors_effect       = prior("normal", parameters = list(mean = 0, sd = 0.25)),
  priors_effect_null  = prior("spike",  parameters = list(location = 0)),
  prior_scale         = "fishers_z",

  # some additional settings
  parallel = FALSE, seed = 1, do_not_fit = TRUE
)

missing <- list.files("../models/MetaRegression/fit_PSMA_reg10/")
missing <- gsub("m_", "", missing)
missing <- gsub(".RDS", "", missing)
missing <- seq_along(fit_PSMA_reg10$models)[!seq_along(fit_PSMA_reg10$models) %in% as.numeric(missing)]

cl <- parallel::makeCluster(23)
parallel::clusterExport(cl, c("fit_PSMA_reg10"))
parallel::parSapplyLB(cl, missing, function(i){

  library(RoBMA)
  temp_model <- RoBMA:::.fit_RoBMA_model(fit_PSMA_reg10, i)
  saveRDS(temp_model, file = paste0("../models/MetaRegression/fit_PSMA_reg10/", "m_", i, ".RDS"), compress = "xz")

})

parallel::stopCluster(cl)

fit_PSMA_reg10 <- RoBMA.reg(
  # specify the model formula and data input
  formula          = ~ 1 + education_outcome + students_gender + location + design + endogenity_control,
  data             = Kroupova2021,
  effect_direction = "negative",
  test_predictors  = c("education_outcome", "students_gender", "location", "design", "endogenity_control"),

  # specify informative prior for the effect size parameter under the alternative hypothesis
  # and a specify a null hypothesis of no effect
  priors             = list(
    education_outcome  = prior_factor("mnormal", list(mean = 0, sd = 0.10), contrast = "orthonormal"),
    students_gender    = prior_factor("mnormal", list(mean = 0, sd = 0.10), contrast = "orthonormal"),
    location           = prior_factor("mnormal", list(mean = 0, sd = 0.10), contrast = "orthonormal"),
    design             = prior_factor("mnormal", list(mean = 0, sd = 0.10), contrast = "orthonormal"),
    endogenity_control = prior_factor("mnormal", list(mean = 0, sd = 0.10), contrast = "orthonormal")
  ),
  priors_effect       = prior("normal", parameters = list(mean = 0, sd = 0.25)),
  priors_effect_null  = prior("spike",  parameters = list(location = 0)),
  prior_scale         = "fishers_z",

  # some additional settings
  parallel = FALSE, seed = 1, do_not_fit = TRUE
)

for(i in seq_along(fit_PSMA_reg10$models)){
  fit_PSMA_reg10$models[[i]] <- readRDS(paste0("../models/MetaRegression/fit_PSMA_reg10/", "m_", i, ".RDS"))
}

fit_PSMA_reg10$models        <- BayesTools::models_inference(fit_PSMA_reg10[["models"]])
fit_PSMA_reg10$RoBMA         <- RoBMA:::.ensemble_inference(fit_PSMA_reg10)
fit_PSMA_reg10$coefficients  <- RoBMA:::.compute_coeficients(fit_PSMA_reg10[["RoBMA"]])

fit_PSMA_reg10$add_info[["errors"]]   <- c(fit_PSMA_reg10$add_info[["errors"]],   RoBMA:::.get_model_errors(fit_PSMA_reg10))
fit_PSMA_reg10$add_info[["warnings"]] <- c(fit_PSMA_reg10$add_info[["warnings"]], RoBMA:::.get_model_warnings(fit_PSMA_reg10))

fit_PSMA_reg10 <- RoBMA:::.remove_model_posteriors(fit_PSMA_reg10)
fit_PSMA_reg10 <- RoBMA:::.remove_model_margliks(fit_PSMA_reg10)

class(fit_PSMA_reg10) <- c("RoBMA", "RoBMA.reg")

saveRDS(fit_PSMA_reg10, file = "../models/MetaRegression/fit_PSMA_reg10.RDS", compress = "xz")
summary(fit_PSMA_reg10)







fit_3PP_reg10 <- RoBMA.reg(
  # specify the model formula and data input
  formula          = ~ 1 + education_outcome + students_gender + location + design + endogenity_control,
  data             = Kroupova2021,
  study_ids        = Kroupova2021$study,
  effect_direction = "negative",
  test_predictors  = c("education_outcome", "students_gender", "location", "design", "endogenity_control"),

  # specify informative prior for the effect size parameter under the alternative hypothesis
  # and a specify a null hypothesis of no effect
  priors             = list(
    education_outcome  = prior_factor("mnormal", list(mean = 0, sd = 0.10), contrast = "orthonormal"),
    students_gender    = prior_factor("mnormal", list(mean = 0, sd = 0.10), contrast = "orthonormal"),
    location           = prior_factor("mnormal", list(mean = 0, sd = 0.10), contrast = "orthonormal"),
    design             = prior_factor("mnormal", list(mean = 0, sd = 0.10), contrast = "orthonormal"),
    endogenity_control = prior_factor("mnormal", list(mean = 0, sd = 0.10), contrast = "orthonormal")
  ),
  priors_effect       = prior("normal", parameters = list(mean = 0, sd = 0.25)),
  priors_effect_null  = prior("spike",  parameters = list(location = 0)),
  priors_bias         = list(
    prior_PET(distribution = "Cauchy", parameters = list(0, 1), truncation = list(0, Inf), prior_weights = 1/2),
    prior_PEESE(distribution = "Cauchy", parameters = list(0, 5), truncation = list(0, Inf), prior_weights = 1/2)
  ),
  prior_scale         = "fishers_z",

  # some additional settings
  parallel = FALSE, silent = FALSE, seed = 1, do_not_fit = TRUE
)

missing <- list.files("../models/MetaRegression/fit_3PP_reg10/")
missing <- gsub("m_", "", missing)
missing <- gsub(".RDS", "", missing)
missing <- seq_along(fit_3PP_reg10$models)[!seq_along(fit_3PP_reg10$models) %in% as.numeric(missing)]

for(i in missing){
  temp_model <- RoBMA:::.fit_RoBMA_model(fit_3PP_reg10, i)
  saveRDS(temp_model, file = paste0("../models/MetaRegression/fit_3PP_reg10/", "m_", i, ".RDS"), compress = "xz")
}

cl <- parallel::makeCluster(23)
parallel::clusterExport(cl, c("fit_3PP_reg10"))
parallel::parSapplyLB(cl, missing, function(i){

  library(RoBMA)
  temp_model <- RoBMA:::.fit_RoBMA_model(fit_3PP_reg10, i)
  saveRDS(temp_model, file = paste0("../models/MetaRegression/fit_3PP_reg10/", "m_", i, ".RDS"), compress = "xz")

})

parallel::stopCluster(cl)

fit_3PP_reg10 <- RoBMA.reg(
  # specify the model formula and data input
  formula          = ~ 1 + education_outcome + students_gender + location + design + endogenity_control,
  data             = Kroupova2021,
  study_ids        = Kroupova2021$study,
  effect_direction = "negative",
  test_predictors  = c("education_outcome", "students_gender", "location", "design", "endogenity_control"),

  # specify informative prior for the effect size parameter under the alternative hypothesis
  # and a specify a null hypothesis of no effect
  priors             = list(
    education_outcome  = prior_factor("mnormal", list(mean = 0, sd = 0.10), contrast = "orthonormal"),
    students_gender    = prior_factor("mnormal", list(mean = 0, sd = 0.10), contrast = "orthonormal"),
    location           = prior_factor("mnormal", list(mean = 0, sd = 0.10), contrast = "orthonormal"),
    design             = prior_factor("mnormal", list(mean = 0, sd = 0.10), contrast = "orthonormal"),
    endogenity_control = prior_factor("mnormal", list(mean = 0, sd = 0.10), contrast = "orthonormal")
  ),
  priors_effect       = prior("normal", parameters = list(mean = 0, sd = 0.25)),
  priors_effect_null  = prior("spike",  parameters = list(location = 0)),
  priors_bias         = list(
    prior_PET(distribution = "Cauchy", parameters = list(0, 1), truncation = list(0, Inf), prior_weights = 1/2),
    prior_PEESE(distribution = "Cauchy", parameters = list(0, 5), truncation = list(0, Inf), prior_weights = 1/2)
  ),
  prior_scale         = "fishers_z",

  # some additional settings
  parallel = FALSE, seed = 1, do_not_fit = TRUE
)

for(i in seq_along(fit_3PP_reg10$models)){
  fit_3PP_reg10$models[[i]] <- readRDS(paste0("../models/MetaRegression/fit_3PP_reg10/", "m_", i, ".RDS"))
}

fit_3PP_reg10$models        <- BayesTools::models_inference(fit_3PP_reg10[["models"]])
fit_3PP_reg10$RoBMA         <- RoBMA:::.ensemble_inference(fit_3PP_reg10)
fit_3PP_reg10$coefficients  <- RoBMA:::.compute_coeficients(fit_3PP_reg10[["RoBMA"]])

fit_3PP_reg10$add_info[["errors"]]   <- c(fit_3PP_reg10$add_info[["errors"]],   RoBMA:::.get_model_errors(fit_3PP_reg10))
fit_3PP_reg10$add_info[["warnings"]] <- c(fit_3PP_reg10$add_info[["warnings"]], RoBMA:::.get_model_warnings(fit_3PP_reg10))

fit_3PP_reg10 <- RoBMA:::.remove_model_posteriors(fit_3PP_reg10)
fit_3PP_reg10 <- RoBMA:::.remove_model_margliks(fit_3PP_reg10)

class(fit_3PP_reg10) <- c("RoBMA", "RoBMA.reg")

saveRDS(fit_3PP_reg10, file = "../models/MetaRegression/fit_3PP_reg10.RDS", compress = "xz")
summary(fit_3PP_reg10)
