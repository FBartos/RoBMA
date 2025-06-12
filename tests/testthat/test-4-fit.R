context("(4) Fitting and updating functions")
skip_on_cran()
skip_on_covr()

# Create temporary directory for fitted models that will be used by subsequent tests
temp_fits_dir <- file.path(tempdir(), "RoBMA_test_fits")
dir.create(temp_fits_dir, showWarnings = FALSE, recursive = TRUE)

# Store the temp directory path for other tests to use
Sys.setenv(ROBMA_TEST_FITS_DIR = temp_fits_dir)



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

  fit1 <- try_parallel(RoBMA(d = d, se = d_se, seed = 1, parallel = TRUE,
                             sample = 1000, burnin = 250, adapt = 250, chains = 2, autofit = FALSE,
                             convergence_checks = set_convergence_checks(max_Rhat = 1.25, min_ESS = 100, max_error = 1, max_SD_error = 1)))
  fit1 <- remove_time(fit1)
  expect_true(is.RoBMA(fit1))
  saveRDS(fit1, file = file.path(temp_fits_dir, "fit_1.RDS"), compress = "xz")

  fit4 <- try_parallel(RoBMA(r = r, n = n, seed = 1, model_type = "PSMA", parallel = TRUE,
                             sample = 2500, burnin = 1000, adapt = 500, chains = 2, autofit = FALSE,
                             convergence_checks = set_convergence_checks(max_Rhat = 2, min_ESS = 10, max_error = 1, max_SD_error = 1), algorithm = "ss"))
  fit4 <- remove_time(fit4)
  expect_true(is.RoBMA(fit4))
  saveRDS(fit4, file = file.path(temp_fits_dir, "fit_4.RDS"), compress = "xz")

  fit5 <- try_parallel(RoBMA(d = d, se = d_se, seed = 1, model_type = "PSMA", transformation = "logOR", parallel = TRUE,
                             sample = 2500, burnin = 1000, adapt = 500, chains = 2, autofit = FALSE,
                             convergence_checks = set_convergence_checks(max_Rhat = 2, min_ESS = 10, max_error = 1, max_SD_error = 1), algorithm = "ss"))
  fit5 <- remove_time(fit5)
  expect_true(is.RoBMA(fit5))
  saveRDS(fit5, file = file.path(temp_fits_dir, "fit_5.RDS"), compress = "xz")

  fit6 <- try_parallel(RoBMA(d = -d, se = d_se, seed = 1, model_type = "PSMA", effect_direction = "negative", parallel = TRUE,
                             sample = 2500, burnin = 1000, adapt = 500, chains = 2, autofit = FALSE,
                             convergence_checks = set_convergence_checks(max_Rhat = 2, min_ESS = 10, max_error = 1, max_SD_error = 1), algorithm = "ss"))
  fit6 <- remove_time(fit6)
  expect_true(is.RoBMA(fit6))
  saveRDS(fit6, file = file.path(temp_fits_dir, "fit_6.RDS"), compress = "xz")

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
  expect_true(is.RoBMA(fit2))
  saveRDS(fit2, file = file.path(temp_fits_dir, "fit_2.RDS"), compress = "xz")
})

test_that("RoBMA-PP works", {
  fit3 <- try_parallel(RoBMA(d = d, se = d_se, seed = 1, parallel = TRUE, model_type = "PP",
                             sample = 500, burnin = 250, adapt = 100, chains = 2, autofit = FALSE,
                             convergence_checks = set_convergence_checks(max_Rhat = 2, min_ESS = 10, max_error = 1, max_SD_error = 1)))
  fit3 <- remove_time(fit3)
  expect_true(is.RoBMA(fit3))
  saveRDS(fit3, file = file.path(temp_fits_dir, "fit_3.RDS"), compress = "xz")
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
  expect_true(is.RoBMA(fit7))
  saveRDS(fit7, file = file.path(temp_fits_dir, "fit_7.RDS"), compress = "xz")
})

test_that("Custom models - only null", {
  fit8 <- try_parallel(RoBMA(d = d, se = d_se, seed = 1,
                priors_effect = NULL, priors_heterogeneity = NULL, priors_bias = NULL, parallel = TRUE,
                sample = 500, burnin = 250, adapt = 100, chains = 2, autofit = FALSE,
                convergence_checks = set_convergence_checks(max_Rhat = 2, min_ESS = 10, max_error = 1, max_SD_error = 1)))
  fit8 <- remove_time(fit8)
  expect_true(is.RoBMA(fit8))
  saveRDS(fit8, file = file.path(temp_fits_dir, "fit_8.RDS"), compress = "xz")
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
  expect_true(is.RoBMA(fit9))
  saveRDS(fit9, file = file.path(temp_fits_dir, "fit_9.RDS"), compress = "xz")
})

test_that("Custom models - fixed weightfunctions", {
  fit10 <- try_parallel(RoBMA(d = d, se = d_se, seed = 1,
                priors_bias = prior_weightfunction("one-sided.fixed", list(c(0.10), c(1, .5))),
                priors_effect_null = NULL, priors_heterogeneity_null = NULL, priors_bias_null = NULL, parallel = TRUE,
                sample = 500, burnin = 250, adapt = 100, chains = 2, autofit = FALSE,
                convergence_checks = set_convergence_checks(max_Rhat = 2, min_ESS = 10, max_error = 1, max_SD_error = 1)))
  fit10 <- remove_time(fit10)
  expect_true(is.RoBMA(fit10))
  saveRDS(fit10, file = file.path(temp_fits_dir, "fit_10.RDS"), compress = "xz")
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
  expect_true(is.RoBMA(fit11))
  saveRDS(fit11, file = file.path(temp_fits_dir, "fit_11.RDS"), compress = "xz")
})

test_that("Custom models - updating models", {
  fit12 <- try_parallel(RoBMA(d = d, se = d_se, seed = 1,
                priors_effect = NULL, priors_heterogeneity = NULL, priors_bias = NULL, parallel = TRUE,
                sample = 500, burnin = 250, adapt = 100, chains = 2, autofit = FALSE,
                convergence_checks = set_convergence_checks(max_Rhat = 2, min_ESS = 10, max_error = 1, max_SD_error = 1)))
  fit12 <- update(fit12, prior_effect = prior("normal", list(0, 1), list(0, 1)), prior_heterogeneity = prior_none(), prior_bias = prior_none())
  fit12 <- update(fit12, prior_weights = c(1, 2))
  fit12 <- remove_time(fit12)
  expect_true(is.RoBMA(fit12))
  saveRDS(fit12, file = file.path(temp_fits_dir, "fit_12.RDS"), compress = "xz")
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
  expect_true(is.RoBMA(fit13))
  saveRDS(fit13, file = file.path(temp_fits_dir, "fit_13.RDS"), compress = "xz")

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
  expect_true(is.RoBMA(fit1w))
  # changed after reducing the number of samples
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
  expect_true(is.RoBMA.reg(fit_14))
  saveRDS(fit_14, file = file.path(temp_fits_dir, "fit_14.RDS"), compress = "xz")
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
  expect_true(is.RoBMA.reg(fit_15))
  saveRDS(fit_15, file = file.path(temp_fits_dir, "fit_15.RDS"), compress = "xz")
})

test_that("BiBMA works", {

  fit_16 <- try_parallel(BiBMA(x1 = 0:4, x2 = 2:6, n1 = rep(20,5), n2 = rep(20, 5), seed = 1, parallel = TRUE,
                               sample = 500, burnin = 250, adapt = 100, chains = 2, autofit = FALSE,
                               convergence_checks = set_convergence_checks(max_Rhat = 2, min_ESS = 10, max_error = 1, max_SD_error = 1)))

  # check update
  fit_16a <- update(fit_16, extend_all = TRUE, autofit_control = set_autofit_control(sample_extend = 500))
  expect_equal(nrow(fit_16$models[[1]]$fit$mcmc[[1]]) + 500, nrow(fit_16a$models[[1]]$fit$mcmc[[1]]))

  fit_16 <- remove_time(fit_16)
  expect_true(is.BiBMA(fit_16))
  saveRDS(fit_16, file = file.path(temp_fits_dir, "fit_16.RDS"), compress = "xz")
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
  expect_true(is.BiBMA(fit_17))
  saveRDS(fit_17, file = file.path(temp_fits_dir, "fit_17.RDS"), compress = "xz")
})



