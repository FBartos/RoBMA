context("(3) priors related functions")
skip_on_cran()


# test objects
saved_priors <- readRDS(file = "saved_priors.RDS")
print_priors <- read.table("print_priors.txt", sep = ";", stringsAsFactors = FALSE)

# create priors
p_point        <- prior("point",     parameters = list(location = 1))
p_normal       <- prior("normal",    parameters = list(mean = 0, sd = 1))
p_normal_trunc <- prior("normal",    parameters = list(mean = 1, sd = 1),                truncation = list(lower = 0, upper = 2), prior_odds = 2)
p_cauchy       <- prior("cauchy",    parameters = list(location = 1, scale = 1),         truncation = list(lower = 0, upper = 2))
p_t_cauchy     <- prior("t",         parameters = list(location = 1, scale = 1, df = 1), truncation = list(lower = 0, upper = 2))
p_t            <- prior("t",         parameters = list(location = 1, scale = 1, df = 5), truncation = list(lower = -2, upper = 2))
p_gamma1       <- prior("gamma",     parameters = list(shape = 1, rate  = 2),            truncation = list(lower = 1, upper = Inf))
p_gamma2       <- prior("gamma",     parameters = list(shape = 1, scale = 1/2),          truncation = list(lower = 1, upper = Inf))
p_invgamma     <- prior("invgamma",  parameters = list(shape = 1, scale = .15),          truncation = list(lower = 0, upper = Inf))
p_uniform      <- prior("uniform",   parameters = list(a = 2, b = 3))
p_two.sided1   <- prior("two.sided", parameters = list(steps = c(.05),           alpha = c(1,1)))
p_two.sided2   <- prior("two.sided", parameters = list(steps = c(.05, .10),      alpha = c(1,1,1)))
p_two.sided3   <- prior("two.sided", parameters = list(steps = c(.05, .10, .80), alpha = c(1,3,5,10)))
p_one.sided1   <- prior("one.sided", parameters = list(steps = c(.05, .10),      alpha = c(1,1,1)))
p_one.sided2   <- prior("one.sided", parameters = list(steps = c(.05, .10, .80), alpha1 = c(1,1,1), alpha2 = c(1,1)))

fitted_priors  <- list(p_point, p_normal, p_normal_trunc, p_cauchy, p_t_cauchy, p_t, p_gamma1, p_gamma2, p_invgamma, p_uniform, p_two.sided1, p_two.sided2, p_two.sided3, p_one.sided1, p_one.sided2)


test_that("Priors structure matches", {
  expect_equal(length(saved_priors), length(fitted_priors))
  for(i in 1:length(saved_priors)){
    expect_equal(saved_priors[i], fitted_priors[i])
  }

  # special cases
  expect_equal(p_cauchy, p_t_cauchy)
  expect_equal(p_gamma1, p_gamma2)
})

test_that("Prior print works", {

  for(i in 1:length(saved_priors)){
    expect_equal(print_priors[i,1], print(fitted_priors[[i]], silent = TRUE))
  }


})

test_that("Priors plots work", {
  set.seed(666)
  for(i in 1:length(fitted_priors)){
    expect_doppelganger(paste0("prior_plot_",i),plot(fitted_priors[[i]], plot_type = "ggplot", samples = 1000, points = 100))
  }

  expect_doppelganger("prior_plot_transformed_1", plot(fitted_priors[[2]],  plot_type = "ggplot", samples = 10000, points = 100, mu_transform = "cohens_d"))
  expect_doppelganger("prior_plot_transformed_2", plot(fitted_priors[[3]],  plot_type = "ggplot", samples = 10000, points = 100, mu_transform = "cohens_d"))
  expect_doppelganger("prior_plot_transformed_3", plot(fitted_priors[[2]],  plot_type = "ggplot", samples = 10000, points = 100, mu_transform = "fishers_z"))
  expect_doppelganger("prior_plot_transformed_4", plot(fitted_priors[[3]],  plot_type = "ggplot", samples = 10000, points = 100, mu_transform = "fishers_z"))

  plot_weights1 <- plot(fitted_priors[[14]], plot_type = "ggplot", samples = 1000, points = 100, weights = T)
  plot_weights2 <- plot(fitted_priors[[15]], plot_type = "ggplot", samples = 1000, points = 100, weights = T)

  expect_equal(length(plot_weights1), 2)
  expect_equal(length(plot_weights2), 3)

  for(i in 1:2){
    expect_doppelganger(paste0("prior_plot_weights1_",i),plot_weights1[[i]])
  }

  for(i in 1:3){
    expect_doppelganger(paste0("prior_plot_weights2_",i),plot_weights2[[i]])
  }
})


#### creating / updating the test settings ####
# saved_priors <- list(p_point, p_normal, p_normal_trunc, p_cauchy, p_t_cauchy, p_t, p_gamma1, p_gamma2, p_invgamma, p_uniform, p_two.sided1, p_two.sided2, p_two.sided3, p_one.sided1, p_one.sided2)
# saveRDS(saved_priors, file = "tests/testthat/saved_priors.RDS")
# saved_priors <- readRDS(file = "tests/testthat/saved_priors.RDS")
#
# sink(file = "tests/testthat/print_priors.txt", append = TRUE)
# for(i in 1:length(fitted_priors)){
#   print(fitted_priors[[i]])
#   cat("\n")
# }
# sink(file = NULL)
