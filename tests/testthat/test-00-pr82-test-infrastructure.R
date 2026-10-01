test_that("fit source fingerprints retain changed assertions", {

  source(testthat::test_path("common-functions.R"), local = TRUE)
  first <- quote(test_that("fit", {
    fit <- constructor()
    expect_equal(coef(fit), 1)
  }))
  changed <- quote(test_that("fit", {
    fit <- constructor()
    expect_equal(coef(fit), 2)
  }))
  removed <- quote(test_that("fit", {
    fit <- constructor()
  }))
  fingerprints <- lapply(list(first, changed, removed), .source_hash_normalize_expr)
  expect_false(identical(fingerprints[[1]], fingerprints[[2]]))
  expect_false(identical(fingerprints[[1]], fingerprints[[3]]))
})

test_that("hypothesis aliases forward the object and arguments", {

  object <- structure(list(marker = "model"), class = "alias_fixture")
  captured <- NULL
  local_mocked_bindings(
    hypothesis = function(object, ...) {
      captured <<- list(object = object, arguments = list(...))
      "hypothesis result"
    },
    .package = "RoBMA"
  )
  for (alias in list(bf_hypothesis, BF_hypothesis)) {
    expect_identical(
      alias(object, "mu > 0", density_method = "KDE"),
      "hypothesis result"
    )
    expect_identical(captured[["object"]], object)
    expect_identical(captured[["arguments"]], list("mu > 0", density_method = "KDE"))
  }
})

test_that("scalar multivariate coefficients use the public mu label", {

  data <- structure(list(), mods = FALSE)
  for (intercept in c("intercept", "(mu) intercept")) {
    object <- structure(
      list(data = data, coefficients = stats::setNames(.25, intercept)),
      class = c("brma.mv", "brma")
    )
    expect_identical(coef(object), c(mu = .25))
    expect_identical(names(object[["coefficients"]]), intercept)
  }

  moderated <- structure(
    list(
      data = structure(list(), mods = TRUE),
      coefficients = c("(mu) intercept" = .1, ablat = .2, "alloc[systematic]" = .3)
    ),
    class = c("brma.mv", "brma")
  )
  expect_identical(
    coef(moderated),
    c("(mu) intercept" = .1, ablat = .2, "alloc[systematic]" = .3)
  )
})
