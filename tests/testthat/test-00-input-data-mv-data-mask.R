test_that("brma.mv evaluates diagonal sampling inputs in data", {

  data <- data.frame(
    y = c(0.10, 0.20, -0.05),
    v = c(0.04, 0.09, 0.16),
    s = c(0.20, 0.30, 0.40)
  )

  from_vi <- brma.mv(
    yi                        = y,
    vi                        = v,
    data                      = data,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_data                 = TRUE
  )
  from_sei <- brma.mv(
    yi                        = y,
    sei                       = s,
    data                      = data,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_data                 = TRUE
  )
  from_expression <- brma.mv(
    yi                        = y,
    vi                        = s^2,
    data                      = data,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_data                 = TRUE
  )

  for (object in list(from_vi, from_sei, from_expression)) {
    expect_equal(
      .known_v_materialize(.data_known_v_data(object[["data"]])),
      diag(data[["v"]])
    )
  }
})


test_that("brma.mv diagonal sampling inputs are named-only formals", {

  formal_names <- names(formals(brma.mv))

  expect_gt(match("vi", formal_names), match("...", formal_names))
  expect_gt(match("sei", formal_names), match("...", formal_names))
})
