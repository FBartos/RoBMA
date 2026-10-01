test_that("conditioning preserves the row dimension for one draw", {

  samples <- matrix(c(1, 0), 1L, dimnames = list(NULL, c("a", "b")))
  local_mocked_bindings(
    .conditional_parameter_rows_single = function(parameter, object, posterior_samples) {
      posterior_samples[, parameter] == 1
    }, .package = "RoBMA"
  )
  expect_identical(.conditional_parameter_rows(list(), c("a", "b"), samples), TRUE)
  expect_warning(
    expect_identical(.conditional_parameter_rows(list(), c("a", "b"), samples, "AND"), FALSE),
    "No samples left after conditioning.", fixed = TRUE
  )
  expect_identical(.conditional_parameter_rows(list(), "a", samples), TRUE)
})

test_that("random inclusion draws reject missing and infinite values clearly", {

  local_mocked_bindings(
    .random_inclusion_sd_names = function(...) c(gate = "tau"),
    .random_slab_sd_names = function(...) character(), .package = "RoBMA"
  )
  value <- NA_real_
  local_mocked_bindings(
    JAGS_materialize_draws = function(...) list(matrix(value, 1L,
      dimnames = list(NULL, "gate"))), .package = "BayesTools"
  )
  for (value in c(NA_real_, NaN, Inf, 0.5)) {
    expect_error(.brma_random_inclusion_indicator_chains(list()),
                 "Random-effect indicator draws are not whole numbers.", fixed = TRUE)
  }
})

test_that("collapsed total random effects retain grouping labels", {

  samples <- matrix(1:8, 4L, dimnames = list(NULL, c("u_study[group a]", "u_study[group b]")))
  result <- .select_ranef_components(list(study = samples), component = "total",
    labels = paste0("estimate", 1:6), n_chains = 1L, n_iter = 4L,
    probs = c(0.025, 0.975), data = list())
  expect_identical(colnames(result), c("u[group a]", "u[group b]"))
  expect_equal(unname(as.matrix(result)), unname(samples), ignore_attr = TRUE)
})
