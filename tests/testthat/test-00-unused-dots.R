test_that("brma.norm is exported alias for brma", {

  expect_true("brma.norm" %in% getNamespaceExports("RoBMA"))
  expect_identical(brma.norm, brma)
})

test_that("brma.mv is exported", {

  expect_true("brma.mv" %in% getNamespaceExports("RoBMA"))
})

test_that("top-level fitting constructors reject stale algorithm argument", {

  norm_args <- list(
    yi         = c(0.10, 0.20, 0.15),
    sei        = c(0.05, 0.06, 0.07),
    only_data  = TRUE,
    algorithm  = "ss"
  )
  glmm_args <- list(
    ai         = c(5, 8, 6),
    bi         = c(45, 42, 44),
    ci         = c(4, 7, 5),
    di         = c(46, 43, 45),
    only_data  = TRUE,
    algorithm  = "ss"
  )

  norm_constructors <- list(
    brma      = brma,
    brma.norm = brma.norm,
    brma.mv   = brma.mv,
    RoBMA     = RoBMA,
    BMA       = BMA,
    BMA.norm  = BMA.norm,
    bPET      = bPET,
    bPEESE    = bPEESE,
    bselmodel = bselmodel
  )
  glmm_constructors <- list(
    brma.glmm = brma.glmm,
    BMA.glmm  = BMA.glmm
  )

  for (constructor in norm_constructors) {
    args <- norm_args
    if (identical(constructor, brma.mv)) {
      args[["V"]]  <- diag(args[["sei"]]^2)
      args[["sei"]] <- NULL
    }
    expect_error(
      do.call(constructor, args),
      "Unused argument.*algorithm"
    )
  }
  for (constructor in glmm_constructors) {
    expect_error(
      do.call(constructor, glmm_args),
      "Unused argument.*algorithm"
    )
  }
})

test_that("top-level fitting constructors reject other named dots", {

  expect_error(
    RoBMA(
      yi               = c(0.10, 0.20, 0.15),
      sei              = c(0.05, 0.06, 0.07),
      legacy_algorithm = "ss",
      only_data        = TRUE
    ),
    "Unused argument.*legacy_algorithm"
  )
})

test_that("brma.mv rejects cluster through unused dots", {

  expect_error(
    brma.mv(
      yi        = c(0.10, 0.20, 0.15),
      V         = diag(c(0.05, 0.06, 0.07)^2),
      cluster   = c(1, 1, 2),
      measure   = "GEN",
      only_data = TRUE
    ),
    "'cluster' is not supported in brma.mv()"
  )
})

test_that("unused-dot warning helper reports ignored arguments", {

  expect_warning(
    .warn_unused_dots(
      dots    = list(legacy_argument = 1),
      allowed = character(),
      caller  = "test()"
    ),
    "Unused argument.*legacy_argument"
  )
  expect_silent(.warn_unused_dots(
    dots    = list(lwd = 1),
    allowed = "lwd",
    caller  = "test()"
  ))
})

test_that("predict.brma rejects unused dots before prediction setup", {

  object <- structure(list(), class = "brma")

  expect_error(
    predict.brma(object, legacy_argument = TRUE),
    "Unused argument.*legacy_argument"
  )
})

test_that("top-level fitting constructors validate internal dot types", {

  expect_error(
    RoBMA(
      yi        = c(0.10, 0.20, 0.15),
      sei       = c(0.05, 0.06, 0.07),
      measure   = "SMD",
      only_data = 1
    ),
    "only_data"
  )

  expect_error(
    brma.glmm(
      ai             = c(5, 8, 6),
      bi             = c(45, 42, 44),
      ci             = c(4, 7, 5),
      di             = c(46, 43, 45),
      only_data      = TRUE,
      is_JASP_prefix = TRUE
    ),
    "is_JASP_prefix"
  )
})
