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
    expect_error(
      do.call(constructor, norm_args),
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
