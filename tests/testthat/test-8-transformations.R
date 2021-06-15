context("(8) Effect size transformations")

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
    t   = orig_data$z / orig_data$se,
    p   = pt(orig_data$z / orig_data$se, n_z(orig_data$se) - 2, lower.tail = FALSE)
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
  orig_data_d$p <- pt(orig_data_d$t, orig_data_d$n - 2, lower.tail = FALSE)

  orig_data_r <- data.frame(
    r   = z2r(orig_data_z$z),
    se  = se_z2se_r(orig_data_z$se, orig_data_z$z),
    v   = se_z2se_r(orig_data_z$se, orig_data_z$z)^2,
    lCI = z2r(orig_data_z$lCI),
    uCI = z2r(orig_data_z$uCI),
    n   = n_r(z2r(orig_data_z$z), se_z2se_r(orig_data_z$se, orig_data_z$z))
  )
  orig_data_r$t <- RoBMA:::.t_rn(orig_data_r$r, orig_data_r$n)
  orig_data_r$p <- pt(orig_data_r$t, orig_data_r$n - 2, lower.tail = FALSE)

  orig_data_logOR <- data.frame(
    logOR = z2logOR(orig_data_z$z),
    se    = se_z2se_logOR(orig_data_z$se, orig_data_z$z),
    v     = se_z2se_logOR(orig_data_z$se, orig_data_z$z)^2,
    lCI   = z2logOR(orig_data_z$z) - se_z2se_logOR(orig_data_z$se, orig_data_z$z) * qnorm(0.975),
    uCI   = z2logOR(orig_data_z$z) + se_z2se_logOR(orig_data_z$se, orig_data_z$z) * qnorm(0.975)
  )
  orig_data_logOR$t <- orig_data_logOR$logOR / orig_data_logOR$se
  orig_data_logOR$p <- pnorm(orig_data_logOR$t, lower.tail = FALSE)


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

    # don't check p-values - they are not consistent across transformation
    test_data[i, "p"] <- temp_df[i, "p"]

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

})
