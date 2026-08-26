# ============================================================================ #
# Current loo comparison contract
# ============================================================================ #

source(testthat::test_path("common-functions.R"))

test_that("outcome hashes use a cross-version canonical encoding", {

  normal_object <- brma(
    yi                        = c(0.1, -0.2),
    sei                       = c(0.3, 0.4),
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  binomial_object <- brma.glmm(
    ai        = c(1L, 2L),
    bi        = c(9L, 8L),
    ci        = c(3L, 4L),
    di        = c(7L, 6L),
    measure   = "OR",
    only_data = TRUE
  )
  poisson_object <- brma.glmm(
    x1i       = c(3L, 5L),
    x2i       = c(2L, 4L),
    t1i       = c(10, 12),
    t2i       = c(11, 13),
    measure   = "IRR",
    only_data = TRUE
  )
  clustered_object <- brma(
    yi        = c(0.1, -0.2, 0.3),
    sei       = c(0.3, 0.4, 0.5),
    weights   = c(1, 2, 3),
    cluster   = c("b", "a", "b"),
    measure   = "GEN",
    only_data = TRUE
  )
  unclustered_object <- brma(
    yi        = c(0.1, -0.2, 0.3),
    sei       = c(0.3, 0.4, 0.5),
    weights   = c(1, 2, 3),
    measure   = "GEN",
    only_data = TRUE
  )
  known_v_object <- brma.mv(
    yi        = c(0.1, 0.2, 0.3),
    V         = matrix(c(
      0.04, 0.01, 0.00,
      0.01, 0.09, 0.00,
      0.00, 0.00, 0.16
    ), nrow = 3, byrow = TRUE),
    measure   = "GEN",
    only_data = TRUE
  )

  expect_identical(.get_outcome_hash(normal_object),
                   "v1:297cae5202260996")
  expect_identical(.get_outcome_hash(binomial_object),
                   "v1:3081885a00045b6b")
  expect_identical(.get_outcome_hash(poisson_object),
                   "v1:3a9cb7d444143c10")
  expect_identical(.get_outcome_hash(clustered_object),
                   "v1:1cb5b0cd7ad933fb")
  expect_identical(.get_outcome_hash(clustered_object),
                   .get_outcome_hash(unclustered_object))
  expect_identical(.get_outcome_hash(known_v_object),
                   "v1:0fb2f0ab7e994106")
  expect_identical(
    .outcome_hash_bytes(list(value = -0)),
    .outcome_hash_bytes(list(value = 0))
  )
  expect_false(identical(
    .outcome_hash_bytes(list(value = 1L)),
    .outcome_hash_bytes(list(value = 1))
  ))
  expect_false(identical(
    .outcome_hash_bytes(list(a = 1, b = 2)),
    .outcome_hash_bytes(list(b = 2, a = 1))
  ))
  expect_false(identical(
    .outcome_hash_bytes(list(1)),
    .outcome_hash_bytes(structure(list(1), names = ""))
  ))
  expect_false(identical(
    .outcome_hash_bytes(list(value = 0)),
    .outcome_hash_bytes(list(value = NA_real_))
  ))
  expect_false(identical(
    .outcome_hash_bytes(list(value = NA_real_)),
    .outcome_hash_bytes(list(value = NaN))
  ))
  expect_false(identical(
    .outcome_hash_bytes(list(value = Inf)),
    .outcome_hash_bytes(list(value = -Inf))
  ))
  expect_false(identical(
    .outcome_hash_bytes(list(value = NA_integer_)),
    .outcome_hash_bytes(list(value = NA_real_))
  ))

  utf8_value   <- "\u00e9"
  latin1_value <- iconv(utf8_value, from = "UTF-8", to = "latin1")
  Encoding(utf8_value) <- "UTF-8"
  expect_false(is.na(latin1_value))
  expect_identical(
    .outcome_hash_bytes(list(value = utf8_value)),
    .outcome_hash_bytes(list(value = latin1_value))
  )
})


test_that("loo_compare preserves the released numeric comparison table", {

  set.seed(11)
  log_lik_a <- matrix(stats::rnorm(2000, mean = -2, sd = 0.3), ncol = 10)
  log_lik_b <- log_lik_a + rep(seq(-0.02, 0.02, length.out = 10), each = 200)

  loo_a <- suppressWarnings(loo::loo(log_lik_a))
  loo_b <- suppressWarnings(loo::loo(log_lik_b))
  target <- list(
    unit             = "estimate",
    retained_context = "remaining_data",
    data_hash        = "loo-contract",
    target           = "estimate_log_score"
  )
  attr(loo_a, "RoBMA_target") <- target
  attr(loo_b, "RoBMA_target") <- target

  out      <- loo_compare(loo_a, loo_b)
  upstream <- do.call(
    get("loo_compare.default", envir = asNamespace("loo"), inherits = FALSE),
    list(loo_a, loo_b)
  )
  legacy_columns <- c(
    "elpd_diff", "se_diff", "elpd_loo", "se_elpd_loo",
    "p_loo", "se_p_loo", "looic", "se_looic"
  )

  expect_identical(class(out), c("compare.loo", "matrix", "array"))
  expect_true(is.numeric(out))
  expect_identical(colnames(out), legacy_columns)
  expect_identical(rownames(out), as.character(upstream[["model"]]))
  expect_equal(
    as.vector(out),
    as.vector(as.matrix(upstream[, legacy_columns, drop = FALSE]))
  )
  expect_true(all(is.finite(out)))
  expect_true("candidate" %in% rownames(loo_compare(loo_a, candidate = loo_b)))

  simplified     <- capture.output(print(out))
  full           <- capture.output(print(out, simplify = FALSE))
  upstream_print <- capture.output(print(upstream))
  expect_match(paste(simplified, collapse = "\n"), "elpd_diff")
  expect_false(any(grepl("elpd_loo", simplified, fixed = TRUE)))
  expect_true(any(grepl("elpd_loo", full, fixed = TRUE)))
  expect_true(any(grepl("model", upstream_print, fixed = TRUE)))
})


test_that("loo_compare preserves the released numeric WAIC table", {

  set.seed(12)
  log_lik_a <- matrix(stats::rnorm(2000, mean = -2, sd = 0.3), ncol = 10)
  log_lik_b <- log_lik_a + rep(seq(-0.02, 0.02, length.out = 10), each = 200)

  waic_a <- suppressWarnings(loo::waic(log_lik_a))
  waic_b <- suppressWarnings(loo::waic(log_lik_b))
  target <- list(
    unit             = "estimate",
    retained_context = "remaining_data",
    data_hash        = "waic-contract",
    target           = "estimate_log_score"
  )
  attr(waic_a, "RoBMA_target") <- target
  attr(waic_b, "RoBMA_target") <- target

  out      <- loo_compare(waic_a, waic_b)
  upstream <- do.call(
    get("loo_compare.default", envir = asNamespace("loo"), inherits = FALSE),
    list(waic_a, waic_b)
  )
  legacy_columns <- c(
    "elpd_diff", "se_diff", "elpd_waic", "se_elpd_waic",
    "p_waic", "se_p_waic", "waic", "se_waic"
  )

  expect_identical(class(out), c("compare.loo", "matrix", "array"))
  expect_true(is.numeric(out))
  expect_identical(colnames(out), legacy_columns)
  expect_identical(rownames(out), as.character(upstream[["model"]]))
  expect_equal(
    as.vector(out),
    as.vector(as.matrix(upstream[, legacy_columns, drop = FALSE]))
  )
  expect_true(all(is.finite(out)))

  simplified <- capture.output(print(out))
  full       <- capture.output(print(out, simplify = FALSE))
  expect_false(any(grepl("elpd_waic", simplified, fixed = TRUE)))
  expect_true(any(grepl("elpd_waic", full, fixed = TRUE)))
})
