context("(4) Fitting and updating functions")
skip_on_cran()
skip_on_covr()

# test objects
saved_files <- paste0("fit_", 1:17, ".RDS")
saved_fits  <- list()
for(i in seq_along(saved_files)){
  saved_fits[[i]] <- readRDS(file = file.path("../results/fits", saved_files[i]))
}

# functions simplifying the comparison
remove_time  <- function(fit){
  for(m in seq_along(fit[["models"]])){
    if(is.null(fit$models[[m]]$fit))next
    fit$models[[m]]$fit$timetaken       <- NULL
    fit$models[[m]]$fit$runjags.version <- NULL
  }
  if(!is.null(fit[["model"]])){
    fit$model$fit$timetaken       <- NULL
    fit$model$fit$runjags.version <- NULL
  }

  return(fit)
}
clean_all    <- function(fit, only_samples = TRUE, remove_call = FALSE){
  if(only_samples){
    fit$data     <- NULL
    fit$add_info <- NULL
    fit$control  <- NULL
    fit$models   <- NULL
    fit$model    <- NULL
  }
  if(remove_call){
    fit$call <- NULL
  }
  return(fit)
}
try_parallel <- function(x, rep = 3){
  temp_fit <- NULL
  i        <- 0
  while(is.null(temp_fit) & i < rep){
    temp_fit <- tryCatch(eval(x), error = function(e) NULL)
    i        <- i + 1
  }
  return(temp_fit)
}


# create mock data
k    <- 3
n    <- c(15, 20, 50)
r    <- c(.1, .2, .3)
d    <- r2d(r)
d_se <- se_d(d, n)


test_that("Default model (RoBMA-PSMA) works", {

  fit1 <- try_parallel(RoBMA(d = d, se = d_se, seed = 1, parallel = TRUE))
  fit1 <- remove_time(fit1)
  expect_equal(clean_all(saved_fits[[1]]), clean_all(fit1))

  fit4 <- try_parallel(RoBMA(r = r, n = n, seed = 1, model_type = "PSMA", parallel = TRUE,
                             sample = 2500, burnin = 1000, adapt = 500, chains = 2, autofit = FALSE,
                             convergence_checks = set_convergence_checks(max_Rhat = 2, min_ESS = 10, max_error = 1, max_SD_error = 1), algorithm = "ss"))
  fit4 <- remove_time(fit4)
  expect_equal(clean_all(saved_fits[[4]]), clean_all(fit4))

  fit5 <- try_parallel(RoBMA(d = d, se = d_se, seed = 1, model_type = "PSMA", transformation = "logOR", parallel = TRUE,
                             sample = 2500, burnin = 1000, adapt = 500, chains = 2, autofit = FALSE,
                             convergence_checks = set_convergence_checks(max_Rhat = 2, min_ESS = 10, max_error = 1, max_SD_error = 1), algorithm = "ss"))
  fit5 <- remove_time(fit5)
  expect_equal(clean_all(saved_fits[[5]]), clean_all(fit5))

  fit6 <- try_parallel(RoBMA(d = -d, se = d_se, seed = 1, model_type = "PSMA", effect_direction = "negative", parallel = TRUE,
                             sample = 2500, burnin = 1000, adapt = 500, chains = 2, autofit = FALSE,
                             convergence_checks = set_convergence_checks(max_Rhat = 2, min_ESS = 10, max_error = 1, max_SD_error = 1), algorithm = "ss"))
  fit6 <- remove_time(fit6)
  expect_equal(clean_all(saved_fits[[6]]), clean_all(fit6))

  # verify that the transformations and etc holds, up to MCMC error
  expect_equal(coef(fit1)[1:2], coef(fit4)[1:2], 0.02)
  expect_equal(coef(fit1)[1:2], coef(fit5)[1:2], 0.02)
  # the effect size is in the opposite direction for fit6
  expect_equal(coef(fit1)[1],  -coef(fit6)[1], 0.01)
  expect_equal(coef(fit1)[2],   coef(fit6)[2], 0.01)

  # lower precision for weights
  expect_equal(coef(fit1)[3:8], coef(fit4)[3:8], 0.03)
  expect_equal(coef(fit1)[3:8], coef(fit5)[3:8], 0.03)
  expect_equal(coef(fit1)[3:8], coef(fit6)[3:8], 0.03)

  # PET is also stable (and in this case PEESE as well, since it is low)
  expect_equal(coef(fit1)[9:10], coef(fit4)[9:10], 0.02)
  expect_equal(coef(fit1)[9:10], coef(fit6)[9:10], 0.02)
  expect_equal(coef(fit1)[9:10], coef(fit6)[9:10], 0.02)
})

test_that("RoBMA-2w works", {
  fit2 <- try_parallel(RoBMA(d = d, se = d_se, seed = 1, parallel = TRUE, model_type = "2w",
                             sample = 2500, burnin = 1000, adapt = 500, chains = 2, autofit = FALSE, algorithm = "ss",
                             convergence_checks = set_convergence_checks(max_Rhat = 2, min_ESS = 10, max_error = 1, max_SD_error = 1)))
  fit2 <- remove_time(fit2)
  expect_equal(clean_all(saved_fits[[2]]), clean_all(fit2))
})

test_that("RoBMA-PP works", {
  fit3 <- try_parallel(RoBMA(d = d, se = d_se, seed = 1, parallel = TRUE, model_type = "PP",
                             sample = 500, burnin = 250, adapt = 100, chains = 2, autofit = FALSE,
                             convergence_checks = set_convergence_checks(max_Rhat = 2, min_ESS = 10, max_error = 1, max_SD_error = 1)))
  fit3 <- remove_time(fit3)
  expect_equal(clean_all(saved_fits[[3]]), clean_all(fit3))
})

test_that("Custom models - only alternative", {
  fit7 <- try_parallel(RoBMA(d = d, se = d_se, seed = 1,
                priors_bias = list(
                  prior_weightfunction("one-sided", list(c(0.10), c(1, 1))),
                  prior_PET("normal", list(0, 1))
                ),
                priors_effect_null = NULL, priors_heterogeneity_null = NULL, priors_bias_null = NULL, parallel = TRUE,
                sample = 2500, burnin = 1000, adapt = 500, chains = 2, autofit = FALSE, algorithm = "ss",
                convergence_checks = set_convergence_checks(max_Rhat = 2, min_ESS = 10, max_error = 1, max_SD_error = 1)))
  fit7 <- remove_time(fit7)
  expect_equal(clean_all(saved_fits[[7]]), clean_all(fit7))
})

test_that("Custom models - only null", {
  fit8 <- try_parallel(RoBMA(d = d, se = d_se, seed = 1,
                priors_effect = NULL, priors_heterogeneity = NULL, priors_bias = NULL, parallel = TRUE,
                sample = 500, burnin = 250, adapt = 100, chains = 2, autofit = FALSE,
                convergence_checks = set_convergence_checks(max_Rhat = 2, min_ESS = 10, max_error = 1, max_SD_error = 1)))
  fit8 <- remove_time(fit8)
  expect_equal(clean_all(saved_fits[[8]]), clean_all(fit8))
})

test_that("Custom models - only null (non-point)", {
  fit9 <- try_parallel(RoBMA(d = d, se = d_se, seed = 1,
                priors_effect_null = prior("normal", list(0, 1)),
                priors_heterogeneity_null = prior("invgamma", list(1, 0.15)),
                priors_bias_null = list(
                  prior_weightfunction("one-sided", list(c(0.10), c(1, 1))),
                  prior_PET("normal", list(0, 1))
                ),
                priors_effect = NULL, priors_heterogeneity = NULL, priors_bias = NULL, parallel = TRUE,
                sample = 2500, burnin = 1000, adapt = 500, chains = 2, autofit = FALSE, algorithm = "ss",
                convergence_checks = set_convergence_checks(max_Rhat = 2, min_ESS = 10, max_error = 1, max_SD_error = 1)))
  fit9 <- remove_time(fit9)
  expect_equal(clean_all(saved_fits[[9]]), clean_all(fit9))
})

test_that("Custom models - fixed weightfunctions", {
  fit10 <- try_parallel(RoBMA(d = d, se = d_se, seed = 1,
                priors_bias = prior_weightfunction("one-sided.fixed", list(c(0.10), c(1, .5))),
                priors_effect_null = NULL, priors_heterogeneity_null = NULL, priors_bias_null = NULL, parallel = TRUE,
                sample = 500, burnin = 250, adapt = 100, chains = 2, autofit = FALSE,
                convergence_checks = set_convergence_checks(max_Rhat = 2, min_ESS = 10, max_error = 1, max_SD_error = 1)))
  fit10 <- remove_time(fit10)
  expect_equal(clean_all(saved_fits[[10]]), clean_all(fit10))
})

test_that("Custom models - unknown effect size", {
  fit11 <- try_parallel(RoBMA(y = d, se = d_se, seed = 1,
                priors_bias = list(
                  prior_weightfunction("two-sided", list(c(0.10), c(1, 1))),
                  prior_PET("normal", list(0, 1))
                ), parallel = TRUE,
                sample = 2500, burnin = 1000, adapt = 500, chains = 2, autofit = FALSE, algorithm = "ss",
                convergence_checks = set_convergence_checks(max_Rhat = 2, min_ESS = 10, max_error = 1, max_SD_error = 1)))
  fit11 <- remove_time(fit11)
  expect_equal(clean_all(saved_fits[[11]]), clean_all(fit11))
})

test_that("Custom models - updating models", {
  fit12 <- try_parallel(RoBMA(d = d, se = d_se, seed = 1,
                priors_effect = NULL, priors_heterogeneity = NULL, priors_bias = NULL, parallel = TRUE,
                sample = 500, burnin = 250, adapt = 100, chains = 2, autofit = FALSE,
                convergence_checks = set_convergence_checks(max_Rhat = 2, min_ESS = 10, max_error = 1, max_SD_error = 1)))
  fit12 <- update(fit12, prior_effect = prior("normal", list(0, 1), list(0, 1)), prior_heterogeneity = prior_none(), prior_bias = prior_none())
  fit12 <- update(fit12, prior_weights = c(1, 2))
  fit12 <- remove_time(fit12)
  expect_equal(clean_all(saved_fits[[12]]), clean_all(fit12))
})

test_that("Main settings work", {
  fit_settings <- try_parallel(RoBMA(d = d, se = d_se, seed = 1,
                                     thin = 2, sample = 300, burnin = 250, adapt = 100, chains = 1, autofit = FALSE,
                                     convergence_checks = set_convergence_checks(max_Rhat = 2, min_ESS = 10, max_error = 1, max_SD_error = 1),
                        priors_effect_null = NULL, priors_heterogeneity = NULL, priors_bias = NULL, parallel = TRUE))

  expect_equal(fit_settings$models[[1]]$fit$thin,   2)
  expect_equal(fit_settings$models[[1]]$fit$burnin, 350)
  expect_equal(fit_settings$models[[1]]$fit$sample, 300)

  expect_equal(length(fit_settings$models[[1]]$fit$mcmc),   1)
  expect_equal(dim(fit_settings$models[[1]]$fit$mcmc[[1]]), c(300, 2))
})

test_that("Convergence warnings work", {
  fit_warnings <- suppressWarnings(
    try_parallel(RoBMA(d = d, se = d_se, seed = 1, autofit = FALSE,
                       priors_effect_null = NULL, priors_heterogeneity = NULL, priors_bias = NULL, parallel = TRUE,
                       thin = 2, sample = 500, burnin = 500, adapt = 100, chains = 2,
          convergence_checks = set_convergence_checks(max_Rhat = 1.01, min_ESS = 1000, max_error = 0.001, max_SD_error = 0.002)))
  )

  expect_equal(
    suppressWarnings(RoBMA::check_RoBMA(fit_warnings)),
    c(
      "Model (1): R-hat 1.026 is larger than the set target (1.01).",
      "Model (1): ESS 829 is lower than the set target (1000).",
      "Model (1): MCMC error 0.00733 is larger than the set target (0.001).",
      "Model (1): MCMC SD error 0.035 is larger than the set target (0.002)."
    )
  )
})

test_that("3-level models work", {

  fit13 <- suppressWarnings(try_parallel(RoBMA(d = d, se = d_se, study_ids = c(1,1,2), seed = 1, parallel = TRUE,
                              autofit = FALSE, thin = 2, sample = 500, burnin = 250, adapt = 100, chains = 1,
                              convergence_checks = set_convergence_checks(max_Rhat = 2, min_ESS = 10, max_error = 1, max_SD_error = 1))))
  fit13 <- remove_time(fit13)
  expect_equal(clean_all(saved_fits[[13]]), clean_all(fit13))

})

test_that("weighted models work", {

  temp_data <- combine_data(
    d         = c(d[1], rep(d[2], 2), rep(d[3], 3)),
    se        = c(d_se[1], rep(d_se[2], 2), rep(d_se[3], 3)),
    weight    = c(1, rep(1/2, 2), rep(1/3, 3))
  )

  fit1w <- try_parallel(suppressWarnings(RoBMA(data = temp_data, seed = 1, parallel = TRUE,
                                               sample = 2500, burnin = 1000, adapt = 500, chains = 2, autofit = FALSE, algorithm = "ss",
                                               convergence_checks = set_convergence_checks(max_Rhat = 2, min_ESS = 10, max_error = 1, max_SD_error = 1))))
  fit1w <- remove_time(fit1w)
  # changed after reducing the number of samples
  # expect_equal(clean_all(saved_fits[[1]], remove_call = TRUE), clean_all(fit1w, remove_call = TRUE))
})

test_that("BMA regression work", {

  df_reg <- data.frame(
    d       = c(rep(-1, 5), rep(0, 5), rep(1, 5)),
    se      = rep(0.1, 15),
    mod_cat = c(rep("A", 5), rep("B", 5), rep("C", 5)),
    mod_con = c((1:15)/15)
  )

  set.seed(1)
  fit_14 <- try_parallel(suppressWarnings(RoBMA.reg(~ mod_cat + mod_con, data = df_reg, priors_bias = NULL, seed = 1, parallel = TRUE,
                                                    sample = 500, burnin = 250, adapt = 100, chains = 2, autofit = FALSE,
                                                    convergence_checks = set_convergence_checks(max_Rhat = 2, min_ESS = 10, max_error = 1, max_SD_error = 1))))
  fit_14 <- remove_time(fit_14)
  expect_equal(clean_all(saved_fits[[14]]), clean_all(fit_14))
})

test_that("RoBMA (simplified) regression with custom priors work", {

  df_reg <- data.frame(
    d       = scale(c((1:15)/15)),
    se      = rep(0.1, 15),
    mod_con = scale(c((1:15)/15))
  )

  set.seed(1)
  fit_15 <- try_parallel(suppressWarnings(RoBMA.reg(~ mod_con, data = df_reg,
                                                    priors = list(
                                                      "mod_con" = list(
                                                        "null" = prior("normal", list(0,    0.05)),
                                                        "alt"  = prior("normal", list(0.30, 0.15))
                                                      )
                                                    ),
                                                    priors_effect_null   = NULL,
                                                    priors_heterogeneity = NULL,
                                                    priors_bias          = list(
                                                      prior_weightfunction(distribution = "two.sided", parameters = list(alpha = c(1, 1), steps = c(0.05)), prior_weights = 1/2),
                                                      prior_PET(distribution = "Cauchy", parameters = list(0, 1), truncation = list(0, Inf), prior_weights = 1/2)
                                                    ),
                                                    seed = 1, parallel = TRUE,
                                                    sample = 2500, burnin = 1000, adapt = 500, chains = 2, autofit = FALSE, algorithm = "ss",
                                                    convergence_checks = set_convergence_checks(max_Rhat = 2, min_ESS = 10, max_error = 1, max_SD_error = 1))))
  fit_15 <- remove_time(fit_15)
  expect_equal(clean_all(saved_fits[[15]]), clean_all(fit_15))
})

test_that("BiBMA works", {

  fit_16 <- try_parallel(BiBMA(x1 = 0:4, x2 = 2:6, n1 = rep(20,5), n2 = rep(20, 5), seed = 1, parallel = TRUE,
                               sample = 500, burnin = 250, adapt = 100, chains = 2, autofit = FALSE,
                               convergence_checks = set_convergence_checks(max_Rhat = 2, min_ESS = 10, max_error = 1, max_SD_error = 1)))

  # check update
  fit_16a <- update(fit_16, extend_all = TRUE, autofit_control = set_autofit_control(sample_extend = 500))
  expect_equal(nrow(fit_16$models[[1]]$fit$mcmc[[1]]) + 500, nrow(fit_16a$models[[1]]$fit$mcmc[[1]]))

  fit_16 <- remove_time(fit_16)
  expect_equal(clean_all(saved_fits[[16]]), clean_all(fit_16))
})

test_that("BiBMA.reg works", {

  df_reg <- data.frame(
    x1 = c(5, 6, 4, 5, 6, 12, 11, 10, 13, 12),
    x2 = c(6, 5, 6, 4, 5,  6,  5,  6,  4,  5),
    n1 = rep(20, 10),
    n2 = rep(20, 10),
    mod_cat = c(rep("A", 5), rep("B", 5))
  )

  fit_17 <- try_parallel(BiBMA.reg( ~ mod_cat, data = df_reg, seed = 1, parallel = TRUE, rescale_priors = 2,
                               sample = 1500, burnin = 750, adapt = 750, chains = 2, autofit = FALSE, algorithm = "ss"))

  fit_17 <- remove_time(fit_17)
  expect_equal(clean_all(saved_fits[[17]]), clean_all(fit_17))
})


#### creating / updating the test settings ####
if(FALSE){
  saved_fits <- list(fit1, fit2, fit3, fit4, fit5, fit6, fit7, fit8, fit9, fit10, fit11, fit12, fit13, fit_14, fit_15, fit_16)

  for(i in 1:length(saved_fits)){
    saved_fits[[i]] <- remove_time(saved_fits[[i]])
  }

  for(i in 1:length(saved_fits)){
    saveRDS(saved_fits[[i]], file = file.path("tests/results/fits/", paste0("fit_",i,".RDS")), compress  = "xz")
  }

  # package version update
  # test objects
  saved_files <- paste0("fit_", 1:16, ".RDS")
  saved_fits  <- list()
  for(i in seq_along(saved_files)){
    temp_fit <- readRDS(file = file.path("tests/results/fits", saved_files[i]))
  #  temp_fit <- RoBMA:::.update_object(temp_fit)
    temp_fit <- remove_time(temp_fit)
    saveRDS(temp_fit, file = file.path("tests/results/fits/", paste0("fit_",i,".RDS")), compress  = "xz")
  }

  # single model update
  fit_14 <- remove_time(fit_14)
  saveRDS(fit_14, file = file.path("tests/results/fits/", paste0("fit_",14,".RDS")), compress  = "xz")
}
