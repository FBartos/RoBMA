context("(8) Effect size transformations")
skip_on_cran()

# Based on examples from: Borenstein, M., Hedges, L. V., Higgins, J. P., & Rothstein, H. R. (2011). Introduction to meta-analysis. John Wiley & Sons.
test_that("standard error of Cohen's d  (p. 28)", {
  expect_equal(
    round(se_d(0.5970, 100), 4),
    0.2044
  )
})
test_that("correlation to Fisher's z (p. 42)", {
  expect_equal(
    round(r2z(0.5), 4),
    0.5493
  )
})
test_that("standard error of Fisher's z (p. 42)", {
  expect_equal(
    round(se_z(100), 4),
    0.1015
  )
})
test_that("logOR to Cohen's d (p. 47)", {
  expect_equal(
    round(logOR2d(0.9069), 4),
    0.5000
  )
})
test_that("variance of logOR to Cohen's d (p. 47)", {
  expect_equal(
    round(se_logOR2se_d(sqrt(0.0676)), 3),
    round(sqrt(0.0205), 3)
  )
})
test_that("Cohen's d to logOR (p. 47)", {
  expect_equal(
    round(d2logOR(0.5), 4),
    0.9069
  )
})
test_that("variance Cohen's d to logOR (p. 48)", {
  expect_equal(
    round(se_d2se_logOR(sqrt(0.0205)), 3),
    round(sqrt(0.0676), 3)
  )
})
test_that("r to Cohen's d (p. 48)", {
  expect_equal(
    round(r2d(0.50), 4),
    1.1547
  )
})
test_that("variance of r to Cohen's d (p. 48)", {
  expect_equal(
    round(se_r2se_d(sqrt(0.0058), 0.50), 3),
    round(sqrt(0.0550), 3)
  )
})
test_that("Cohen's d to r (p. 49)", {
  expect_equal(
    round(d2r(1.1547), 4),
    0.5
  )
})
test_that("variance of Cohen's d to r (p. 49)", {
  expect_equal(
    round(se_d2se_r(sqrt(0.0550), 1.1547), 3),
    round(sqrt(0.0058), 3)
  )
})

# sanity checks for back and forth transformations
test_that("back and forth transformations for Cohen's d work", {

  d     <- 0.5
  se_d  <- 0.2

  expect_equal(se_d, se_d(d, n_d(d, se_d)))

  expect_equal(d, r2d(d2r(d)))
  expect_equal(d, z2d(d2z(d)))
  expect_equal(d, logOR2d(d2logOR(d)))

  expect_equal(se_d, se_r2se_d(se_d2se_r(se_d, d), d2r(d)))
  expect_equal(se_d, se_z2se_d(se_d2se_z(se_d, d), d2z(d)))
  expect_equal(se_d, se_logOR2se_d(se_d2se_logOR(se_d)))

})
test_that("back and forth transformations for correlation work", {

  r     <- 0.15
  se_r  <- 0.10

  expect_equal(se_r, se_r(r, n_r(r, se_r)))

  expect_equal(r, d2r(r2d(r)))
  expect_equal(r, z2r(r2z(r)))
  expect_equal(r, logOR2r(r2logOR(r)))

  expect_equal(se_r, se_d2se_r(se_r2se_d(se_r, r), r2d(r)))
  expect_equal(se_r, se_z2se_r(se_r2se_z(se_r, r), r2z(r)))
  expect_equal(se_r, se_logOR2se_r(se_r2se_logOR(se_r, r), r2logOR(r)))

})
test_that("back and forth transformations for Fisher's z work", {

  z     <- 0.66
  se_z  <- 0.36

  expect_equal(se_z, se_z(n_z(se_z)))

  expect_equal(z, d2z(z2d(z)))
  expect_equal(z, r2z(z2r(z)))
  expect_equal(z, logOR2z(z2logOR(z)))

  expect_equal(se_z, se_d2se_z(se_z2se_d(se_z, z), z2d(z)))
  expect_equal(se_z, se_r2se_z(se_z2se_r(se_z, z), z2r(z)))
  expect_equal(se_z, se_logOR2se_z(se_z2se_logOR(se_z, z), z2logOR(z))) #

})
test_that("back and forth transformations for logOR works", {

  logOR     <- -0.50
  se_logOR  <-  0.16

  expect_equal(logOR, d2logOR(logOR2d(logOR)))
  expect_equal(logOR, r2logOR(logOR2r(logOR)))
  expect_equal(logOR, z2logOR(logOR2z(logOR)))

  expect_equal(se_logOR, se_d2se_logOR(se_logOR2se_d(se_logOR)))
  expect_equal(se_logOR, se_r2se_logOR(se_logOR2se_r(se_logOR, logOR), logOR2r(logOR))) #
  expect_equal(se_logOR, se_z2se_logOR(se_logOR2se_z(se_logOR, logOR), logOR2z(logOR))) #

})

# additional tests for the data combination and transformation functions
test_that("data combination works", {

  # data set preparation function
  orig_data <- data.frame(
    z  = seq(-3, 3, .10),
    se = se_z(seq(20, 20 + 60*10, length.out = 61))
  )

  orig_data_z <- data.frame(
    z   = orig_data$z,
    se  = orig_data$se,
    v   = orig_data$se^2,
    lCI = orig_data$z - orig_data$se * qnorm(0.975),
    uCI = orig_data$z + orig_data$se * qnorm(0.975),
    n   = n_z(orig_data$se),
    t   = orig_data$z / orig_data$se
  )
  orig_data_d <- data.frame(
    d   = z2d(orig_data_z$z),
    se  = se_z2se_d(orig_data_z$se, orig_data_z$z),
    v   = se_z2se_d(orig_data_z$se, orig_data_z$z)^2,
    lCI = z2d(orig_data_z$lCI),
    uCI = z2d(orig_data_z$uCI),
    n   = n_d(z2d(orig_data_z$z), se_z2se_d(orig_data_z$se, orig_data_z$z))
  )
  orig_data_d$t <- RoBMA:::.t_dn(orig_data_d$d, orig_data_d$n)

  orig_data_r <- data.frame(
    r   = z2r(orig_data_z$z),
    se  = se_z2se_r(orig_data_z$se, orig_data_z$z),
    v   = se_z2se_r(orig_data_z$se, orig_data_z$z)^2,
    lCI = z2r(orig_data_z$lCI),
    uCI = z2r(orig_data_z$uCI),
    n   = n_r(z2r(orig_data_z$z), se_z2se_r(orig_data_z$se, orig_data_z$z))
  )
  orig_data_r$t <- RoBMA:::.t_rn(orig_data_r$r, orig_data_r$n)

  orig_data_logOR <- data.frame(
    logOR = z2logOR(orig_data_z$z),
    se    = se_z2se_logOR(orig_data_z$se, orig_data_z$z),
    v     = se_z2se_logOR(orig_data_z$se, orig_data_z$z)^2,
    lCI   = z2logOR(orig_data_z$z) - se_z2se_logOR(orig_data_z$se, orig_data_z$z) * qnorm(0.975),
    uCI   = z2logOR(orig_data_z$z) + se_z2se_logOR(orig_data_z$se, orig_data_z$z) * qnorm(0.975)
  )
  orig_data_logOR$t <- orig_data_logOR$logOR / orig_data_logOR$se


  test_data   <- matrix(NA, ncol = 0, nrow = nrow(orig_data))
  for(var in c("d", "r", "z", "logOR", "se", "v", "n", "t", "lCI", "uCI")){
    test_data <- cbind(test_data, NA)
    colnames(test_data)[ncol(test_data)] <- var
  }
  test_data <- data.frame(test_data)

  set.seed(666)
  use_effect_measure      <- sample(c("z", "d", "r", "logOR"),   nrow(test_data), TRUE)
  use_variability_measure <- sample(c("n", "t", "se", "v", "CI"), nrow(test_data), TRUE)
  use_variability_measure[use_effect_measure == "logOR" & use_variability_measure == "n"] <- "se"

  # add scrambled records
  for(i in 1:nrow(test_data)){

    temp_df <- switch(
      use_effect_measure[i],
      "z"     = orig_data_z,
      "r"     = orig_data_r,
      "d"     = orig_data_d,
      "logOR" = orig_data_logOR)

    test_data[i, use_effect_measure[i]] <- temp_df[i, use_effect_measure[i]]
    if(use_variability_measure[i] == "CI"){
      test_data[i, c("lCI", "uCI")] <- temp_df[i, c("lCI", "uCI")]
    }else{
      test_data[i, use_variability_measure[i]] <- temp_df[i, use_variability_measure[i]]
    }
  }

  parsed_z <- combine_data(data = test_data)

  colnames(parsed_z)[1] <- "z"
  expect_equal(as.matrix(parsed_z[, c("z", "se")]), as.matrix(orig_data_z[, c("z", "se")]))


  parsed_d <- combine_data(data = test_data, transformation = "cohens_d")

  colnames(parsed_d)[1] <- "d"
  expect_equal(as.matrix(parsed_d[, c("d", "se")]), as.matrix(orig_data_d[, c("d", "se")]))


  parsed_r <- combine_data(data = test_data, transformation = "r")

  colnames(parsed_r)[1] <- "r"
  expect_equal(as.matrix(parsed_r[, c("r", "se")]), as.matrix(orig_data_r[, c("r", "se")]))


  parsed_logOR <- combine_data(data = test_data, transformation = "logOR")

  colnames(parsed_logOR)[1] <- "logOR"
  expect_equal(as.matrix(parsed_logOR[, c("logOR", "se")]), as.matrix(orig_data_logOR[, c("logOR", "se")]))

  # mixing OR and logOR
  c_logOR     <- c(0.5, 0.1, 0.3)
  c_logOR.lCI <- c_logOR - 0.1
  c_logOR.uCI <- c_logOR + 0.1

  comb_logOR <- combine_data(logOR = c_logOR, lCI = c_logOR.lCI, uCI = c_logOR.uCI, transformation = "fishers_z")
  comb_OR    <- combine_data(OR = exp(c_logOR), lCI = exp(c_logOR.lCI), uCI = exp(c_logOR.uCI), transformation = "fishers_z")
  expect_equivalent(comb_logOR, comb_OR)

  # using combine_data.bi
  com.RoBMA   <- .combine_data.bi(x1 = c(0, 1), x2 = c(1, 2), n1  = c(5, 6), n2  = c(6, 5), return_all = TRUE, transformation = "logOR")
  com.metafor <- metafor::escalc( ai = c(0, 1), ci = c(1, 2), n1i = c(5, 6), n2i = c(6, 5), measure = "OR")

  com.RoBMA <- data.frame(
    yi = com.RoBMA$logOR,
    vi = com.RoBMA$se^2
  )
  expect_equivalent(com.RoBMA, com.metafor)

})

# test that input checks work
test_that("input check works", {

  expect_error(combine_data(d = 0, lCI = 0.1, uCI = 0.2), "All effect sizes must be within the CI intervals.")
  expect_error(combine_data(r = 0, lCI = 0.1, uCI = 0.2), "All effect sizes must be within the CI intervals.")
  expect_error(combine_data(y = 0, lCI = 0.1, uCI = 0.2), "All effect sizes must be within the CI intervals.")
  expect_error(combine_data(z = 0, lCI = 0.1, uCI = 0.2), "All effect sizes must be within the CI intervals.")
  expect_error(combine_data(d = 0, lCI = 0.1, uCI = 0.0), "'lCI' must be lower than 'uCI'.")
  expect_error(combine_data(d = 0, se = 0), "The 'se' must be higher than 0.")
  expect_error(combine_data(r = 1, se = 1), "The 'r' must be lower than 1.")
  expect_error(combine_data(r =-1, se = 1), "The 'r' must be higher than -1.")
  expect_error(combine_data(r = 0, d  = 1), "The data do not contain any variability measure.")
  expect_error(combine_data(se = 1, v  = 1), "The data do not contain any effect size measure.")
  expect_error(combine_data(d = c(0,1), se = c(NA, 1)), "At least one variability measure per study must be supplied.")
  expect_error(combine_data(d = c(NA,1), se = c(1, 1)), "At least one effect size measure per study must be supplied.")
  expect_error(combine_data(d = c(NA,1), r = c(0, 0), se = c(1, 1)), "Only one effect size measure per study can be supplied.")
  expect_error(combine_data(d = c(0, 1), v = c(2, 1), se = c(1, 1)), "Only one variability measure per study can be supplied.")
  expect_error(combine_data(d = c(0, 1), lCI = c(NA, 1), se = c(1, NA)), "Either none or both CI bounds must be supplied.")
  expect_error(combine_data(d = c(NA,1), y = c(1, NA), se = c(1, 1)), "Standardized and general effect sizes cannot be combined.")

})

# density transformations work (utilizing Jacobian)
test_that("density transformations works", {

  # visually compares histogram of transformed samples to a transformed density of untransformed samples
  test_transformation <- function(samples, from, to){

    hist(do.call(eval(parse(text = paste0(".", from, "2", to, "$fun"))), args = list(samples)), freq = FALSE, breaks = 50,
         main = paste0(from, "->", to), xlab = to)

    if(from == "OR"){
      den <- density(samples, from = 0)
    }else{
      den <- density(samples)
    }

    den_x_trans <- do.call(eval(parse(text = paste0(".", from, "2", to, "$fun"))), args = list(den$x))
    lines(
      den_x_trans,
      den$y * do.call(eval(parse(text = paste0(".", from, "2", to, "$jac"))), args = list(den_x_trans))
    )
  }

  set.seed(1)
  x_samples <- rnorm(100000, .3, .15)

  for(from in c("d", "r", "z", "logOR", "OR")){
    for(to in c("d", "r", "z", "logOR", "OR")){
      if(from == to){
        next
      }else if(from == "OR"){
        vdiffr::expect_doppelganger(paste0("transformation_",from,"2",to), function()test_transformation(x_samples[x_samples > 0], from, to))
      }else{
        vdiffr::expect_doppelganger(paste0("transformation_",from,"2",to), function()test_transformation(x_samples, from, to))
      }
    }
  }

})

# density scaling work (utilizing Jacobian)
test_that("density scaling works", {

  # visually compares histogram of transformed samples to a transformed density of untransformed samples
  test_transformation <- function(samples, from, to){

    hist(do.call(eval(parse(text = paste0(".", from, "2", to, "$fun"))), args = list(samples)), freq = FALSE, breaks = 50,
         main = paste0(from, "->", to), xlab = to)

    if(from == "OR"){
      den <- density(samples, from = 0)
    }else{
      den <- density(samples)
    }

    den_x_trans <- do.call(eval(parse(text = paste0(".", from, "2", to, "$fun"))), args = list(den$x))
    lines(
      den_x_trans,
      den$y * do.call(eval(parse(text = paste0(".", from, "2", to, "$jac"))), args = list(den_x_trans))
    )
  }

  set.seed(1)
  x_samples <- rnorm(100000, .3, .15)

  for(from in c("d", "r", "z", "logOR")){
    for(to in c("d", "r", "z", "logOR")){
      if(from == to){
        next
      }else{
        vdiffr::expect_doppelganger(paste0("scaling_",from,"2",to), function()test_transformation(x_samples, from, to))
      }
    }
  }

})

# JAGS implementation works
# density scaling work (utilizing jakobian)
test_that("density scaling works", {

  # re-load the module
  RoBMA:::.load_RoBMA_module()


  # scaling
  var_names <- expand.grid(from = c("d", "r", "z", "logOR"), to = c("d", "r", "z", "logOR"))
  var_names <- var_names[var_names$from != var_names$to, ]

  model_syntax <-
    "model
    {
      x ~ dnorm(-.45, pow(0.10, -2))

      d2r_scaled     = scale_d2r(x)
      d2z_scaled     = scale_d2z(x)
      d2logOR_scaled = scale_d2logOR(x)

      r2d_scaled     = scale_r2d(x)
      r2z_scaled     = scale_r2z(x)
      r2logOR_scaled = scale_r2logOR(x)

      z2d_scaled     = scale_z2d(x)
      z2r_scaled     = scale_z2r(x)
      z2logOR_scaled = scale_z2logOR(x)

      logOR2d_scaled = scale_logOR2d(x)
      logOR2r_scaled = scale_logOR2r(x)
      logOR2z_scaled = scale_logOR2z(x)

    }"

  model <- rjags::jags.model(file = textConnection(model_syntax),quiet = TRUE)
  fit   <- rjags::coda.samples(model = model, variable.names = c("x", paste0(var_names$from, "2", var_names$to, "_scaled")),
                               n.iter = 100, quiet = TRUE, progress.bar = "none")

  for(from in c("d", "r", "z", "logOR")){
    for(to in c("d", "r", "z", "logOR")){
      if(from == to){
        next
      }else{
        expect_equal(RoBMA:::.scale(fit[,"x"][[1]], from, to), fit[, paste0(from, "2", to, "_scaled")][[1]])
      }
    }
  }


  # effect sizes
  model_syntax <-
    "model
    {
      x ~ dnorm(-.45, pow(0.10, -2))

      d2r_transformed     = d2r(x)
      d2z_transformed     = d2z(x)
      d2logOR_transformed = d2logOR(x)

      r2d_transformed     = r2d(x)
      r2z_transformed     = r2z(x)
      r2logOR_transformed = r2logOR(x)

      z2d_transformed     = z2d(x)
      z2r_transformed     = z2r(x)
      z2logOR_transformed = z2logOR(x)

      logOR2d_transformed = logOR2d(x)
      logOR2r_transformed = logOR2r(x)
      logOR2z_transformed = logOR2z(x)

    }"

  model <- rjags::jags.model(file = textConnection(model_syntax),quiet = TRUE)
  fit   <- rjags::coda.samples(model = model, variable.names = c("x", paste0(var_names$from, "2", var_names$to, "_transformed")),
                               n.iter = 100, quiet = TRUE, progress.bar = "none")

  for(from in c("d", "r", "z", "logOR")){
    for(to in c("d", "r", "z", "logOR")){
      if(from == to){
        next
      }else{
        expect_equal(do.call(RoBMA:::.get_transformation(from, to), arg = list(fit[,"x"][[1]])),
                     fit[, paste0(from, "2", to, "_transformed")][[1]])
      }
    }
  }

  # standard errors
  model_syntax <-
    "model
    {
      se  ~ dunif(0.10, .30)
      x   ~ dunif(-.55, -.35)

      d2r_se_transformed     = se_d2se_r(se, x)
      d2z_se_transformed     = se_d2se_z(se, x)
      d2logOR_se_transformed = se_d2se_logOR(se, x)

      r2d_se_transformed     = se_r2se_d(se, x)
      r2z_se_transformed     = se_r2se_z(se, x)
      r2logOR_se_transformed = se_r2se_logOR(se, x)

      z2d_se_transformed     = se_z2se_d(se, x)
      z2r_se_transformed     = se_z2se_r(se, x)
      z2logOR_se_transformed = se_z2se_logOR(se, x)

      logOR2d_se_transformed = se_logOR2se_d(se, x)
      logOR2r_se_transformed = se_logOR2se_r(se, x)
      logOR2z_se_transformed = se_logOR2se_z(se, x)
    }"

  model <- rjags::jags.model(file = textConnection(model_syntax),quiet = TRUE)
  fit   <- rjags::coda.samples(model = model, variable.names = c("se", "x", paste0(var_names$from, "2", var_names$to, "_se_transformed")),
                               n.iter = 100, quiet = TRUE, progress.bar = "none")

  for(from in c("d", "r", "z", "logOR")){
    for(to in c("d", "r", "z", "logOR")){
      if(from == to){
        next
      }else{
        expect_equal(do.call(RoBMA:::.get_transformation_se(from, to), arg = list(fit[,"se"][[1]], fit[,"x"][[1]])),
                     fit[, paste0(from, "2", to, "_se_transformed")][[1]])
      }
    }
  }


})
