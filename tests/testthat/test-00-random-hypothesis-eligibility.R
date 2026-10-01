test_that("random hypothesis eligibility uses structural rather than sampled constancy", {

  for (fixed in c(FALSE, TRUE)) {
    prior <- if (fixed) {
      BayesTools::prior("spike", list(location = .2))
    } else {
      BayesTools::prior("gamma", list(shape = 2, rate = 2))
    }
    object <- single_sd_random_object(prior)
    quantities <- hypothesis_quantities(object)
    random <- quantities[quantities[["component"]] == "random", , drop = FALSE]

    expect_gt(nrow(random), 0L)
    expect_identical(unique(random[["direction_test"]]), !fixed)
    expect_identical(unique(random[["point_test"]]), !fixed)
    expect_identical(any(grepl("fixed by the fitted model", random[["reason"]])), fixed)
  }
})


test_that("sampled random quantities take region prior odds from the exact prior", {

  # The posterior draws of the study SD are constant (0.2), yet the quantity
  # is sampled: its region hypotheses use the prior odds of its canonical
  # prior density, Gamma(2, 2), not prior draws.
  object <- single_sd_random_object(
    BayesTools::prior("gamma", list(shape = 2, rate = 2))
  )
  out <- suppressWarnings(hypothesis(
    object, "`(mu) tau(intercept)` > .5 vs `(mu) tau(intercept)` <= .5",
    density_method = "KDE", seed = 1, n_samples = 20L, columns = "all"
  ))
  expect_equal(
    as.numeric(out[["prior"]]),
    stats::pgamma(.5, 2, 2, lower.tail = FALSE) / stats::pgamma(.5, 2, 2),
    tolerance = 1e-10
  )
})
