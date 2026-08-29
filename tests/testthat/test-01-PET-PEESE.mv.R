source(testthat::test_path("common-functions.R"))
skip_on_cran()


.fit_pet_peese_mv_input <- function() {

  data <- data.frame(
    yi    = c(0.05, 0.16, 0.22, 0.34, 0.08, 0.19),
    study = factor(rep(c("s1", "s2", "s3"), each = 2L)),
    obs   = factor(seq_len(6L)),
    x     = rep(c(0, 1), 3L)
  )
  V <- kronecker(
    diag(3L),
    matrix(c(0.010, 0.004, 0.004, 0.014), nrow = 2L)
  )

  list(data = data, V = V)
}


test_that("bPET.mv fits known-V PET regression with random effects", {

  skip_refit_if_cached("bPET.mv")
  input <- .fit_pet_peese_mv_input()
  fit <- bPET.mv(
    yi = yi, V = input[["V"]], mods = ~ x, random = ~ 1 | study,
    data = input[["data"]], measure = "GEN",
    prior_unit_information_sd = 1,
    known_v_parameterization = "whitened",
    chains = 2, sample = 300, burnin = 150, adapt = 500,
    seed = 201, silent = TRUE,
    convergence_checks = set_convergence_checks(
      max_Rhat = NULL,
      min_ESS  = NULL
    )
  )
  fit <- suppressWarnings(add_loo(fit, unit = "estimate"))
  fit <- add_marglik(
    fit,
    repetitions = 1,
    maxiter     = 2000,
    silent      = TRUE
  )
  save_fit("bPET.mv_random", fit, info = input)

  expect_identical(
    class(fit),
    c("bPET.mv", "bPET", "brma.mv", "brma.norm", "brma")
  )
  expect_s3_class(fit[["loo"]][["estimate"]], "loo")
  expect_true(is.finite(as.numeric(logml(fit))))
})


test_that("bPEESE.mv fits known-V PEESE regression with random effects", {

  skip_refit_if_cached("bPEESE.mv")
  input <- .fit_pet_peese_mv_input()
  fit <- bPEESE.mv(
    yi = yi, V = input[["V"]], mods = ~ x,
    random = list(study = ~ 1 | study, estimate = ~ 1 | obs),
    data = input[["data"]], measure = "GEN",
    prior_unit_information_sd = 1,
    known_v_parameterization = "block_mvn",
    chains = 2, sample = 300, burnin = 150, adapt = 500,
    seed = 202, silent = TRUE,
    convergence_checks = set_convergence_checks(
      max_Rhat = NULL,
      min_ESS  = NULL
    )
  )
  fit <- suppressWarnings(add_loo(fit, unit = "estimate"))
  fit <- add_marglik(
    fit,
    repetitions = 1,
    maxiter     = 2000,
    silent      = TRUE
  )
  save_fit("bPEESE.mv_random", fit, info = input)

  expect_identical(
    class(fit),
    c("bPEESE.mv", "bPEESE", "brma.mv", "brma.norm", "brma")
  )
  expect_s3_class(fit[["loo"]][["estimate"]], "loo")
  expect_true(is.finite(as.numeric(logml(fit))))
})
