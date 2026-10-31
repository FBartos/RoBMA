test_that("normal inputs complete missing standard errors from variances", {

  input <- data.frame(yi = c(0.1, 0.2, 0.3),
                      vi = c(0.04, NA_real_, NA_real_),
                      sei = c(NA_real_, 0.3, NA_real_))
  outcome <- .check_and_list_data.outcome.norm(
    quote(brma(yi = yi, vi = vi, sei = sei)), input, environment(),
    effect_direction = "positive"
  )
  expect_equal(outcome[["data_outcome"]][["sei"]], c(0.2, 0.3, NA_real_),
               tolerance = 0)
  expect_error(.check_and_list_data.outcome.norm(
    quote(brma(yi = 0.1, vi = 1, sei = 2)), NULL, environment(),
    effect_direction = "positive"
  ), "^The provided 'vi' and 'sei' values are inconsistent\\.$")
})

test_that("multivariate diagonal inputs complete complementary missingness", {

  inputs <- .check_and_list_data.mv_known_v_input(
    V = NULL, vi = c(0.04, NA_real_, NA_real_),
    sei = c(NA_real_, 0.3, NA_real_), k = 3L
  )
  expect_identical(inputs[["missing_for_na"]], c(FALSE, FALSE, TRUE))
  expect_equal(inputs[["V"]], c(0.04, 0.09, NA_real_), tolerance = 0)
  expect_true(.check_and_list_data.mv_validate_hidden_inputs(
    inputs[["hidden"]], !inputs[["missing_for_na"]]
  ))
})

test_that("function names are not random-effect data variables", {

  expect_error(.check_and_list_data.random_variable("t", NULL, baseenv()),
               "^Cannot find the random-effect variable \\('t'\\)\\.$")
  expect_identical(.check_and_list_data.random_variable(
    "t", data.frame(t = 1:3), baseenv()), 1:3)
  values <- list2env(list(t = c("a", "a", "b")), parent = baseenv())
  expect_identical(.check_and_list_data.random_variable("t", NULL, values),
                   c("a", "a", "b"))
})
