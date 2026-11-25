probe_state("test-00-review-input-missing.R")
test_that("normal inputs complete missing standard errors from variances", {

  input <- data.frame(yi = c(0.1, 0.2, 0.3),
                      vi = c(0.04, NA_real_, NA_real_),
                      sei = c(NA_real_, 0.3, NA_real_))
  expect_warning(object <- brma(yi = yi, vi = vi, sei = sei, data = input, only_data = TRUE),
                 "1 observation(s) removed", fixed = TRUE)
  expect_equal(object$data$outcome$sei, c(0.2, 0.3), tolerance = 0)
  expect_error(brma(yi = 0.1, vi = 1, sei = 2, only_data = TRUE),
               "^The provided 'vi' and 'sei' values are inconsistent\\.$")
})

test_that("multivariate diagonal inputs complete complementary missingness", {

  inputs <- .check_and_list_data.mv_known_v_input(
    V = NULL, vi = c(0.04, NA_real_, NA_real_),
    sei = c(NA_real_, 0.3, NA_real_), k = 3L
  )
  expect_identical(inputs[["missing_for_na"]], c(FALSE, FALSE, TRUE))
  expect_equal(inputs[["V"]], c(0.04, 0.09, NA_real_), tolerance = 0)
  expect_identical(inputs$sampling_inputs$vi, c(0.04, NA_real_, NA_real_))
  expect_identical(inputs$sampling_inputs$sei, c(NA_real_, 0.3, NA_real_))
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
