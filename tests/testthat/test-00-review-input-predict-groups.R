test_that("marginalized new-effect draws preserve supplied group dependence", {

  withr::local_seed(82)
  scales <- matrix(rep(c(1, 2, 3), each = 1000L), 1000L, 3L)
  draws <- .predict_known_v_marginalized_random_term_draws(
    list(group_label = "study"),
    list(location = data.frame(study = c("new_a", "new_a", "new_b"))),
    scales
  )
  expect_identical(draws[, 2L], 2 * draws[, 1L])
  expect_false(identical(draws[, 3L], 3 * draws[, 1L]))
  expect_true(all(is.finite(draws)))
})

test_that("missing marginalized grouping metadata fails before drawing", {

  withr::local_seed(82)
  before <- .Random.seed
  locations <- list(NULL, data.frame(study = "a"),
                    data.frame(other = c("a", "a")),
                    data.frame(study = c("a", NA_character_)))
  for (location in locations) {
    expect_error(.predict_known_v_marginalized_random_term_draws(
      list(group_label = "study"), list(location = location), matrix(1, 4L, 2L)
    ), paste0("^Marginalized random-effect prediction requires complete grouping ",
               "variables in 'newdata'\\.$"))
  }
  expect_error(.predict_known_v_marginalized_random_term_draws(
    list(group_label = "study:subgroup"),
    list(location = data.frame(study = c("a", "a"), subgroup = c("x", NA))),
    matrix(1, 4L, 2L)
  ), "requires complete grouping")
  expect_identical(.Random.seed, before)
})
