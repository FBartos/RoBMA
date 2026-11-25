test_that("empty allocation priors give their validation error without warnings", {

  expect_warning(expect_error(
    .assign_prior.heterogeneity_allocation_mixture(list(), list()),
    "^At least one prior distribution needs to be defined for 'heterogeneity_allocation'\\.$"
  ), NA)
})
