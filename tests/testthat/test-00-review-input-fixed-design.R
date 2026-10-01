test_that("point-fixed formulas map columns without an intercept", {

  point <- function(value) BayesTools::prior("point", list(location = value))
  predictors <- data.frame(a = c(-1, 0, 2), b = c(4, 2, -1))
  for (with_intercept in c(FALSE, TRUE)) {
    X <- stats::model.matrix(if (with_intercept) ~ a + b else ~ 0 + a + b,
                            predictors)
    design <- list(
      model_matrix = X, assign = attr(X, "assign"), parameter = "mu",
      model_terms = c(if (with_intercept) "intercept", "a", "b"),
      prior_list = c(if (with_intercept) list(mu_intercept = point(5)),
                     list(mu_a = point(2), mu_b = point(3)))
    )
    expect_equal(.formula_design_fixed_values(design),
                 2 * predictors$a + 3 * predictors$b +
                   if (with_intercept) 5 else 0,
                 tolerance = 0)
  }
})
