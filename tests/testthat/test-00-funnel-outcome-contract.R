test_that("outcome-mode funnels reject study-specific fitted quantities", {

  make_object <- function(mods = FALSE, scale = FALSE,
                          classes = "brma") {

    data <- list(outcome = data.frame(yi = c(.1, .2), sei = c(.2, .3)))
    attr(data, "mods")  <- mods
    attr(data, "scale") <- scale
    structure(list(data = data), class = classes)
  }

  expected <- paste0(
    "Outcome-mode funnel plots are not supported for models with location ",
    "or scale predictors.*Use 'residual = TRUE'"
  )

  expect_error(
    funnel(make_object(mods = TRUE), residual = FALSE, as_data = TRUE),
    expected
  )
  expect_error(
    funnel(make_object(scale = TRUE), residual = FALSE, as_data = TRUE),
    expected
  )
})


test_that("intercept-only outcome funnels and automatic residual routing remain", {

  make_object <- function(mods = FALSE, scale = FALSE,
                          classes = "brma") {

    data <- list(outcome = data.frame(yi = c(.1, .2), sei = c(.2, .3)))
    attr(data, "mods")  <- mods
    attr(data, "scale") <- scale
    structure(list(data = data), class = classes)
  }

  testthat::local_mocked_bindings(
    .funnel_data_outcome = function(...) list(mode = "outcome"),
    .funnel_data_residual = function(...) list(mode = "residual"),
    .package = "RoBMA"
  )

  intercept_only <- list(
    normal  = make_object(),
    known_V = make_object(classes = c("brma.mv", "brma")),
    GLMM    = make_object(classes = c("brma.glmm", "brma"))
  )
  for (object in intercept_only) {
    expect_identical(
      funnel(object, residual = FALSE, as_data = TRUE),
      list(mode = "outcome")
    )
  }

  expect_identical(
    funnel(make_object(mods = TRUE), type = "outcome", as_data = TRUE),
    list(mode = "residual")
  )
  expect_identical(
    funnel(make_object(scale = TRUE), type = "outcome", as_data = TRUE),
    list(mode = "residual")
  )
})
