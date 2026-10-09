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

# S05: the grouping variables' values form the key, not their display strings.
.review_colliding_tuple_object <- function() {

  data <- data.frame(a = c("a:b", "a"), b = c("c", "b:c"))
  point <- function(value) BayesTools::prior("point", list(value))
  brma.mv(
    yi = c(0, 0), V = c(.1, .1), random = ~ 1 | a:b, data = data,
    measure = "GEN", prior_unit_information_sd = 1,
    prior_effect = point(0), prior_heterogeneity = point(1),
    only_priors = TRUE, silent = TRUE
  )
}

test_that("marginalized group keys keep component-value tuples distinct", {

  location <- data.frame(
    a = c("a:b", "a", "a:b", "a:", "a"),
    b = c("c", "b:c", "c", "b", ":b")
  )
  keys <- .predict_known_v_marginalized_random_group_keys(
    list(group_label = "a:b"), list(location = location), nrow(location)
  )

  # ("a:b", "c") and ("a", "b:c") share the display string "a:b:c".
  expect_identical(keys[[1L]], keys[[3L]])
  expect_length(unique(keys), 4L)
  expect_false(identical(keys[[1L]], keys[[2L]]))
  # ("a:", "b") and ("a", ":b") share the display string "a::b".
  expect_false(identical(keys[[4L]], keys[[5L]]))
  expect_false(any(keys[4:5] %in% keys[1:3]))
})

test_that("marginalized draws are independent for distinct colliding tuples", {

  object <- .review_colliding_tuple_object()
  term   <- .data_marginalized_random_effects(object[["data"]])[[1L]]
  S      <- 4000L
  posterior_samples <- matrix(
    1, nrow = S, ncol = 1L, dimnames = list(NULL, term[["sd_parameter_names"]][[1L]])
  )
  new_data <- object[["data"]]
  new_data[["outcome"]] <- data.frame(
    yi = c(0, 0, 0), sei = c(0, 0, 0), slab = paste0("new", 1:3)
  )
  new_data[["location"]] <- data.frame(
    a = c("a:b", "a", "a:b"), b = c("c", "b:c", "c")
  )

  withr::local_seed(82)
  draws <- .predict_known_v_marginalized_random_draws(
    object = object, data = new_data, posterior_samples = posterior_samples
  )

  # A repeated tuple shares one draw; distinct tuples do not.
  expect_identical(draws[, 1L], draws[, 3L])
  expect_false(any(draws[, 1L] == draws[, 2L]))
  expect_lt(abs(stats::cor(draws[, 1L], draws[, 2L])), 4 / sqrt(S))
})

test_that("public marginal true-effect draws are independent for colliding tuples", {

  object  <- .review_colliding_tuple_object()
  term    <- .data_marginalized_random_effects(object[["data"]])[[1L]]
  S       <- 4000L
  samples <- cbind(mu = rep(0, S), rep(1, S))
  colnames(samples)[[2L]] <- term[["sd_parameter_names"]][[1L]]

  withr::local_seed(42)
  actual <- as.matrix(predict(
    object, type = "estimate", quiet = TRUE, .posterior_samples = samples
  ))

  expect_false(any(actual[, 1L] == actual[, 2L]))
  expect_lt(abs(stats::cor(actual[, 1L], actual[, 2L])), 4 / sqrt(S))
})
