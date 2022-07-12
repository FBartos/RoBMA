context("(4) Fitting and updating functions")
skip_on_cran()
skip_on_covr()

# test objects
saved_files <- paste0("fit_", 1:13, ".RDS")
saved_fits  <- list()
for(i in seq_along(saved_files)){
  saved_fits[[i]] <- readRDS(file = file.path("../results/fits", saved_files[i]))
}

# functions simplifying the comparison
remove_time  <- function(fit){
  for(m in 1:length(fit$models)){
    if(is.null(fit$models[[m]]$fit))next
    fit$models[[m]]$fit$timetaken       <- NULL
    fit$models[[m]]$fit$runjags.version <- NULL
  }
  return(fit)
}
clean_all  <- function(fit, only_samples = TRUE, remove_call = FALSE){
  if(only_samples){
    fit$data     <- NULL
    fit$add_info <- NULL
    fit$control  <- NULL
    fit$models   <- NULL
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

  fit4 <- try_parallel(RoBMA(r = r, n = n, seed = 1, model_type = "PSMA", parallel = TRUE))
  fit4 <- remove_time(fit4)
  expect_equal(clean_all(saved_fits[[4]]), clean_all(fit4))

  fit5 <- try_parallel(RoBMA(d = d, se = d_se, seed = 1, model_type = "PSMA", transformation = "logOR", parallel = TRUE))
  fit5 <- remove_time(fit5)
  expect_equal(clean_all(saved_fits[[5]]), clean_all(fit5))

  fit6 <- try_parallel(RoBMA(d = -d, se = d_se, seed = 1, model_type = "PSMA", effect_direction = "negative", parallel = TRUE))
  fit6 <- remove_time(fit6)
  expect_equal(clean_all(saved_fits[[6]]), clean_all(fit6))

  # verify that the transformations and etc holds, up to MCMC error
  expect_equal(coef(fit1)[1:2], coef(fit4)[1:2], 0.005)
  expect_equal(coef(fit1)[1:2], coef(fit5)[1:2], 0.005)
  # the effect size is in the opposite direction for fit6
  expect_equal(coef(fit1)[1],  -coef(fit6)[1], 0.005)
  expect_equal(coef(fit1)[2],   coef(fit6)[2], 0.005)

  # lower precision for weights
  expect_equal(coef(fit1)[3:8], coef(fit4)[3:8], 0.02)
  expect_equal(coef(fit1)[3:8], coef(fit5)[3:8], 0.02)
  expect_equal(coef(fit1)[3:8], coef(fit6)[3:8], 0.02)

  # PET is also stable (and in this case PEESE as well, since it is low)
  expect_equal(coef(fit1)[9:10], coef(fit4)[9:10], 0.005)
  expect_equal(coef(fit1)[9:10], coef(fit6)[9:10], 0.005)
  expect_equal(coef(fit1)[9:10], coef(fit6)[9:10], 0.005)
})

test_that("RoBMA-2w works", {
  fit2 <- try_parallel(RoBMA(d = d, se = d_se, seed = 1, parallel = TRUE, model_type = "2w"))
  fit2 <- remove_time(fit2)
  expect_equal(clean_all(saved_fits[[2]]), clean_all(fit2))
})

test_that("RoBMA-PP works", {
  fit3 <- try_parallel(RoBMA(d = d, se = d_se, seed = 1, parallel = TRUE, model_type = "PP"))
  fit3 <- remove_time(fit3)
  expect_equal(clean_all(saved_fits[[3]]), clean_all(fit3))
})

test_that("Custom models - only alternative", {
  fit7 <- try_parallel(RoBMA(d = d, se = d_se, seed = 1,
                priors_bias = list(
                  prior_weightfunction("one-sided", list(c(0.10), c(1, 1))),
                  prior_PET("normal", list(0, 1))
                ),
                priors_effect_null = NULL, priors_heterogeneity_null = NULL, priors_bias_null = NULL, parallel = TRUE))
  fit7 <- remove_time(fit7)
  expect_equal(clean_all(saved_fits[[7]]), clean_all(fit7))
})

test_that("Custom models - only null", {
  fit8 <- try_parallel(RoBMA(d = d, se = d_se, seed = 1,
                priors_effect = NULL, priors_heterogeneity = NULL, priors_bias = NULL, parallel = TRUE))
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
                priors_effect = NULL, priors_heterogeneity = NULL, priors_bias = NULL, parallel = TRUE))
  fit9 <- remove_time(fit9)
  expect_equal(clean_all(saved_fits[[9]]), clean_all(fit9))
})

test_that("Custom models - fixed weightfunctions", {
  fit10 <- try_parallel(RoBMA(d = d, se = d_se, seed = 1,
                priors_bias = prior_weightfunction("one-sided.fixed", list(c(0.10), c(1, .5))),
                priors_effect_null = NULL, priors_heterogeneity_null = NULL, priors_bias_null = NULL, parallel = TRUE))
  fit10 <- remove_time(fit10)
  expect_equal(clean_all(saved_fits[[10]]), clean_all(fit10))
})

test_that("Custom models - unknown effect size", {
  fit11 <- try_parallel(RoBMA(y = d, se = d_se, seed = 1,
                priors_bias = list(
                  prior_weightfunction("two-sided", list(c(0.10), c(1, 1))),
                  prior_PET("normal", list(0, 1))
                ), parallel = TRUE))
  fit11 <- remove_time(fit11)
  expect_equal(clean_all(saved_fits[[11]]), clean_all(fit11))
})

test_that("Custom models - updating models", {
  fit12 <- try_parallel(RoBMA(d = d, se = d_se, seed = 1,
                priors_effect = NULL, priors_heterogeneity = NULL, priors_bias = NULL, parallel = TRUE))
  fit12 <- update(fit12, prior_effect = prior("normal", list(0, 1), list(0, 1)), prior_heterogeneity = prior_none(), prior_bias = prior_none())
  fit12 <- update(fit12, prior_weights = c(1, 2))
  fit12 <- remove_time(fit12)
  expect_equal(clean_all(saved_fits[[12]]), clean_all(fit12))
})

test_that("Main settings work", {
  fit_settings <- try_parallel(RoBMA(d = d, se = d_se, seed = 1, autofit = FALSE, thin = 2, sample = 1000, burnin = 500, adapt = 100, chains = 1,
                        priors_effect_null = NULL, priors_heterogeneity = NULL, priors_bias = NULL, parallel = TRUE))

  expect_equal(fit_settings$models[[1]]$fit$thin,   2)
  expect_equal(fit_settings$models[[1]]$fit$burnin, 500 + 100)
  expect_equal(fit_settings$models[[1]]$fit$sample, 1000)

  expect_equal(length(fit_settings$models[[1]]$fit$mcmc),   1)
  expect_equal(dim(fit_settings$models[[1]]$fit$mcmc[[1]]), c(1000, 2))
})

test_that("Convergence warnings work", {
  fit_warnings <- suppressWarnings(
    try_parallel(RoBMA(d = d, se = d_se, seed = 1, autofit = FALSE, thin = 2, sample = 500, burnin = 500, adapt = 100, chains = 2,
          convergence_checks = set_convergence_checks(max_Rhat = 1.01, min_ESS = 1000, max_error = 0.001, max_SD_error = 0.002),
          priors_effect_null = NULL, priors_heterogeneity = NULL, priors_bias = NULL, parallel = TRUE))
  )

  expect_equal(
    suppressWarnings(RoBMA::check_RoBMA(fit_warnings)),
    c(
      "Model (1): ESS 829 is lower than the set target (1000).",
      "Model (1): MCMC error 0.00733 is larger than the set target (0.001).",
      "Model (1): MCMC SD error 0.035 is larger than the set target (0.002)."
    )
  )
})

test_that("3-level models work", {

  fit13 <- suppressWarnings(try_parallel(RoBMA(d = d, se = d_se, study_ids = c(1,1,2), seed = 1, parallel = TRUE,
                              autofit = FALSE, thin = 2, sample = 1000, burnin = 500, adapt = 100, chains = 1)))
  fit13 <- remove_time(fit13)
  expect_equal(clean_all(saved_fits[[13]]), clean_all(fit13))

})

test_that("weighted models work", {

  # all weights == 1 should correspond to fit1
  temp_data <- combine_data(
    d         = d,
    se        = d_se,
    study_ids = seq_along(d)
  )
  attr(temp_data, "all_independent") <- FALSE

  fit1w <- try_parallel(suppressWarnings(RoBMA(data = temp_data, weighted = TRUE, seed = 1, parallel = TRUE)))
  fit1w <- remove_time(fit1w)
  expect_equal(clean_all(saved_fits[[1]], remove_call = TRUE), clean_all(fit1w, remove_call = TRUE))

  # check that the models are actually weighted
  expect_true(all(grepl("dwwn", sapply(c(2:7,  11:16, 20:25, 29:34), function(i) as.character(fit1w$models[[i]]$fit$model)))))
  expect_true(all(grepl("dwn",  sapply(c(8:10, 17:19, 26:28, 35:36), function(i) as.character(fit1w$models[[i]]$fit$model)))))
})

#### creating / updating the test settings ####
if(FALSE){
  saved_fits <- list(fit1, fit2, fit3, fit4, fit5, fit6, fit7, fit8, fit9, fit10, fit11, fit12, fit13)

  for(i in 1:length(saved_fits)){
    saved_fits[[i]] <- remove_time(saved_fits[[i]])
  }

  for(i in 1:length(saved_fits)){
    saveRDS(saved_fits[[i]], file = file.path("tests/results/fits/", paste0("fit_",i,".RDS")),   compress  = "xz")
  }
}
