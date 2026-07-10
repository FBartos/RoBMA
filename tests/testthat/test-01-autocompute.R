context("Autocompute options")

source(testthat::test_path("common-functions.R"))
skip_on_cran()
skip_if_not_installed("metadat")
skip_if_not_installed("metafor")

test_that("autocompute options are applied for brma()", {

  old_options <- list(
    autocompute.loo     = RoBMA.get_option("autocompute.loo"),
    autocompute.waic    = RoBMA.get_option("autocompute.waic"),
    autocompute.marglik = RoBMA.get_option("autocompute.marglik")
  )
  on.exit(do.call(RoBMA.options, old_options), add = TRUE)

  # Using a very simple model to keep the test-01 fit cheap.
  data(dat.bcg, package = "metadat")
  dat <- metafor::escalc(
    measure = "RR",
    ai      = tpos,
    bi      = tneg,
    ci      = cpos,
    di      = cneg,
    data    = dat.bcg
  )
  dat <- dat[1:3, ]

  RoBMA.options(
    autocompute.loo     = FALSE,
    autocompute.waic    = FALSE,
    autocompute.marglik = FALSE
  )
  fit_default <- suppressWarnings(brma(
    yi      = yi,
    vi      = vi,
    data    = dat,
    measure = "RR",
    seed    = 1,
    silent  = TRUE,
    chains  = 2,
    sample  = 250,
    burnin  = 50,
    adapt   = 100
  ))

  expect_null(fit_default$loo)
  expect_null(fit_default$waic)
  expect_null(fit_default$marglik)

  RoBMA.options(
    autocompute.loo     = TRUE,
    autocompute.waic    = TRUE,
    autocompute.marglik = TRUE
  )
  fit_auto <- suppressWarnings(brma(
    yi      = yi,
    vi      = vi,
    data    = dat,
    measure = "RR",
    seed    = 1,
    silent  = TRUE,
    chains  = 2,
    sample  = 250,
    burnin  = 50,
    adapt   = 100
  ))

  expect_s3_class(fit_auto$loo$estimate, "loo")
  expect_s3_class(fit_auto$waic$estimate, "waic")
  expect_s3_class(fit_auto$marglik, "bridge")
})

test_that("autocompute options are applied for brma.mv() known-V storage", {

  old_options <- list(
    autocompute.loo     = RoBMA.get_option("autocompute.loo"),
    autocompute.waic    = RoBMA.get_option("autocompute.waic"),
    autocompute.marglik = RoBMA.get_option("autocompute.marglik")
  )
  on.exit(do.call(RoBMA.options, old_options), add = TRUE)

  yi <- c(0.10, 0.20, 0.15)
  fit_args <- list(
    yi                        = yi,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    seed                      = 1,
    silent                    = TRUE,
    chains                    = 2,
    sample                    = 250,
    burnin                    = 50,
    adapt                     = 500,
    convergence_checks        = set_convergence_checks(
      max_Rhat = NULL,
      min_ESS  = NULL
    )
  )

  RoBMA.options(
    autocompute.loo     = FALSE,
    autocompute.waic    = FALSE,
    autocompute.marglik = FALSE
  )
  fit_default <- suppressWarnings(do.call(
    brma.mv,
    c(fit_args, list(V = c(0.04, 0.05, 0.06)))
  ))
  expect_null(fit_default[["loo"]])
  expect_null(fit_default[["waic"]])
  expect_null(fit_default[["marglik"]])

  representations <- list(
    diagonal = list(
      V                        = c(0.04, 0.05, 0.06),
      known_v_parameterization = "latent",
      storage                  = "diagonal",
      backend                  = "diagonal"
    ),
    blocks = list(
      V = list(
        matrix(c(0.04, 0.01, 0.01, 0.05), nrow = 2L),
        matrix(0.06, nrow = 1L)
      ),
      known_v_parameterization = "whitened",
      storage                  = "blocks",
      backend                  = "whitened"
    ),
    dense = list(
      V = matrix(
        c(
          0.04, 0.01, 0.005,
          0.01, 0.05, 0.008,
          0.005, 0.008, 0.06
        ),
        nrow = 3L,
        byrow = TRUE
      ),
      known_v_parameterization = "block_mvn",
      storage                  = "dense",
      backend                  = "block_mvn"
    )
  )

  RoBMA.options(
    autocompute.loo     = TRUE,
    autocompute.waic    = TRUE,
    autocompute.marglik = TRUE
  )
  for (name in names(representations)) {
    representation <- representations[[name]]
    fit_auto <- suppressWarnings(do.call(
      brma.mv,
      c(
        fit_args,
        representation[c("V", "known_v_parameterization")]
      )
    ))
    known_V <- .data_known_v_data(fit_auto[["data"]])
    target  <- attr(fit_auto[["marglik"]], "RoBMA_target", exact = TRUE)

    expect_s3_class(fit_auto[["loo"]][["estimate"]], "loo")
    expect_s3_class(fit_auto[["waic"]][["estimate"]], "waic")
    expect_s3_class(fit_auto[["marglik"]], "bridge")
    expect_true(is.finite(logml(fit_auto)), info = name)
    expect_equal(known_V[["storage"]], representation[["storage"]], info = name)
    expect_equal(
      target[["known_v_effective_backend"]],
      representation[["backend"]],
      info = name
    )
  }
})

test_that("product-space constructors autocompute LOO and WAIC only", {

  old_options <- list(
    autocompute.loo     = RoBMA.get_option("autocompute.loo"),
    autocompute.waic    = RoBMA.get_option("autocompute.waic"),
    autocompute.marglik = RoBMA.get_option("autocompute.marglik")
  )
  on.exit(do.call(RoBMA.options, old_options), add = TRUE)

  data(dat.bcg, package = "metadat")
  dat <- metafor::escalc(
    measure = "RR",
    ai      = tpos,
    bi      = tneg,
    ci      = cpos,
    di      = cneg,
    data    = dat.bcg
  )
  dat <- dat[1:3, ]

  RoBMA.options(
    autocompute.loo     = TRUE,
    autocompute.waic    = TRUE,
    autocompute.marglik = TRUE
  )
  fit_auto <- suppressWarnings(BMA.norm(
    yi      = yi,
    vi      = vi,
    data    = dat,
    measure = "RR",
    seed    = 1,
    silent  = TRUE,
    chains  = 2,
    sample  = 250,
    burnin  = 50,
    adapt   = 100
  ))

  expect_s3_class(fit_auto$loo$estimate, "loo")
  expect_s3_class(fit_auto$waic$estimate, "waic")
  expect_null(fit_auto$marglik)
})
