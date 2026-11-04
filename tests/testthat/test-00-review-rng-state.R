context("Preserved R random-number state")

# Every test restores its caller's public RNG state. No test promises to save
# Box-Muller's hidden normal cache after reseeding or consuming normal draws.
.review_rng_cleanup <- function() {

  saved <- .rng_state()
  withr::defer(.rng_state_restore(saved), envir = parent.frame())
}


test_that("Inversion scopes restore public state and subsequent normal draws", {

  .review_rng_cleanup()
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  for (fail in c(FALSE, TRUE)) for (change_kind in c(FALSE, TRUE)) {
    set.seed(203L)
    expected <- stats::rnorm(12L)
    set.seed(203L)
    saved <- .rng_state()
    run_scope <- function() .with_preserved_rng({
      if (change_kind) RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
      stats::runif(7L)
      stats::rnorm(24L)
      if (fail) stop("scope failure", call. = FALSE)
      NULL
    })
    if (fail) {
      expect_error(run_scope(), "scope failure", fixed = TRUE)
    } else {
      expect_null(run_scope())
    }
    expect_identical(RNGkind(), saved[["kind"]])
    expect_identical(.Random.seed, saved[["seed"]])
    expect_identical(stats::rnorm(12L), expected)
  }
})


test_that("no-op and uniform scopes keep a pending Box-Muller normal", {

  .review_rng_cleanup()
  RNGkind("Mersenne-Twister", "Box-Muller", "Rejection")
  for (uniform in c(FALSE, TRUE)) for (fail in c(FALSE, TRUE)) {
    set.seed(203L)
    stats::rnorm(1L)
    expected <- stats::rnorm(5L)
    set.seed(203L)
    stats::rnorm(1L)
    saved <- .rng_state()
    run_scope <- function() .with_preserved_rng({
      if (uniform) stats::runif(7L)
      if (fail) stop("scope failure", call. = FALSE)
      NULL
    })
    if (fail) {
      expect_error(run_scope(), "scope failure", fixed = TRUE)
    } else {
      expect_null(run_scope())
    }
    expect_identical(RNGkind(), saved[["kind"]])
    expect_identical(.Random.seed, saved[["seed"]])
    expect_identical(stats::rnorm(5L), expected)
  }
})


test_that("preserved scopes remove a seed that was originally absent", {

  .review_rng_cleanup()
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  for (initialize in c(FALSE, TRUE)) for (fail in c(FALSE, TRUE)) {
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
    saved_kind <- RNGkind()
    run_scope <- function() .with_preserved_rng({
      stats::runif(7L)
      if (fail) stop("scope failure", call. = FALSE)
      NULL
    }, initialize = initialize)
    if (fail) {
      expect_error(run_scope(), "scope failure", fixed = TRUE)
    } else {
      expect_null(run_scope())
    }
    expect_identical(RNGkind(), saved_kind)
    expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
  }
})


test_that("unchanged legacy Rounding does not emit a restore warning", {

  .review_rng_cleanup()
  expect_warning(
    RNGkind("Mersenne-Twister", "Inversion", "Rounding"),
    "non-uniform 'Rounding' sampler used", fixed = TRUE
  )
  set.seed(203L)
  saved <- .rng_state()
  expect_silent(.with_preserved_rng(stats::runif(7L)))
  expect_identical(RNGkind(), saved[["kind"]])
  expect_identical(.Random.seed, saved[["seed"]])
})
