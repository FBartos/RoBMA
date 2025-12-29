context("Fit helper functions")


# ============================================================================
# Tests for .replace_intercept_with_log1
# ============================================================================

test_that(".replace_intercept_with_log1 handles formula with -1 at end", {

  formula <- ~ x + a - 1
  result  <- RoBMA:::.replace_intercept_with_log1(formula)

  expect_s3_class(result, "formula")
  expect_equal(paste0(deparse(result), collapse = " "), "~x + a + log(1)")
})

test_that(".replace_intercept_with_log1 handles formula with 0 at start", {

  formula <- ~ 0 + x + a
  result  <- RoBMA:::.replace_intercept_with_log1(formula)

  expect_s3_class(result, "formula")
  expect_equal(paste0(deparse(result), collapse = " "), "~x + a + log(1)")
})

test_that(".replace_intercept_with_log1 handles formula with only 0", {

  formula <- ~ 0
  result  <- RoBMA:::.replace_intercept_with_log1(formula)

  expect_s3_class(result, "formula")
  expect_equal(paste0(deparse(result), collapse = " "), "~log(1)")
})

test_that(".replace_intercept_with_log1 handles formula with only -1", {

  formula <- ~ -1
  result  <- RoBMA:::.replace_intercept_with_log1(formula)

  expect_s3_class(result, "formula")
  expect_equal(paste0(deparse(result), collapse = " "), "~log(1)")
})

test_that(".replace_intercept_with_log1 handles single predictor with -1", {

  formula <- ~ x - 1
  result  <- RoBMA:::.replace_intercept_with_log1(formula)

  expect_s3_class(result, "formula")
  expect_equal(paste0(deparse(result), collapse = " "), "~x + log(1)")
})

test_that(".replace_intercept_with_log1 handles single predictor with 0", {

  formula <- ~ 0 + x
  result  <- RoBMA:::.replace_intercept_with_log1(formula)

  expect_s3_class(result, "formula")
  expect_equal(paste0(deparse(result), collapse = " "), "~x + log(1)")
})

test_that(".replace_intercept_with_log1 handles formula with interaction terms", {

  formula <- ~ x * a - 1
  result  <- RoBMA:::.replace_intercept_with_log1(formula)

  expect_s3_class(result, "formula")
  expect_equal(paste0(deparse(result), collapse = " "), "~x * a + log(1)")
})

test_that(".replace_intercept_with_log1 handles formula with nested terms", {

  formula <- ~ x + I(x^2) - 1
  result  <- RoBMA:::.replace_intercept_with_log1(formula)

  expect_s3_class(result, "formula")
  expect_equal(paste0(deparse(result), collapse = " "), "~x + I(x^2) + log(1)")
})

test_that(".replace_intercept_with_log1 handles formula with multiple predictors and 0", {

  formula <- ~ 0 + a + b + c
  result  <- RoBMA:::.replace_intercept_with_log1(formula)

  expect_s3_class(result, "formula")
  expect_equal(paste0(deparse(result), collapse = " "), "~a + b + c + log(1)")
})

test_that(".replace_intercept_with_log1 handles formula with + - 1 pattern", {

  formula <- ~ x + y + - 1
  result  <- RoBMA:::.replace_intercept_with_log1(formula)

  expect_s3_class(result, "formula")
  expect_equal(paste0(deparse(result), collapse = " "), "~x + y + log(1)")
})

test_that(".replace_intercept_with_log1 preserves response variable", {

  formula <- y ~ x + a - 1
  result  <- RoBMA:::.replace_intercept_with_log1(formula)

  expect_s3_class(result, "formula")
  expect_equal(paste0(deparse(result), collapse = " "), "y ~ x + a + log(1)")
})

test_that(".replace_intercept_with_log1 handles formula with spaces around -1", {

  formula <- as.formula("~ x +   a   -   1")
  result  <- RoBMA:::.replace_intercept_with_log1(formula)

  expect_s3_class(result, "formula")
  expect_equal(paste0(deparse(result), collapse = " "), "~x + a + log(1)")
})

test_that(".replace_intercept_with_log1 handles complex formula", {

  formula <- ~ poly(x, 2) + factor(group) + offset(log(n)) - 1
  result  <- RoBMA:::.replace_intercept_with_log1(formula)

  expect_s3_class(result, "formula")
  expect_equal(paste0(deparse(result), collapse = " "), "~poly(x, 2) + factor(group) + offset(log(n)) + log(1)")
})
