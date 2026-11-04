context("Marginal means of factors with a fixed-zero intercept")

# The draws of the fitted slope coordinates of the replay.
.review_slope_draws <- function(object) {

  draws <- as.matrix(object[["fit"]][["mcmc"]][[1L]])

  draws[, startsWith(colnames(draws), "mu_g"), drop = FALSE]
}

.review_level_table <- function(object) {

  `class<-`(summary(marginal_means(object)), "data.frame")
}

test_that("treatment-coded levels keep the reference at the fixed zero", {

  object <- review_factor_replay(
    c("a", "b"), "treatment", variable_intercept = FALSE,
    independent = TRUE, distribution = "normal"
  )
  slope <- .review_slope_draws(object)[, 1L]
  table <- .review_level_table(object)

  expect_identical(rownames(table), c("g[a]", "g[b]"))
  expect_equal(unname(unlist(table["g[a]", c("Mean", "SD")])), c(0, 0))
  expect_equal(table["g[b]", "Mean"], mean(slope), tolerance = 1e-12)
  expect_equal(table["g[b]", "SD"], stats::sd(slope), tolerance = 1e-12)
  expect_equal(
    unname(unlist(table["g[b]", c("0.025", "0.5", "0.975")])),
    unname(stats::quantile(slope, c(.025, .5, .975))),
    tolerance = 1e-12
  )
})

test_that("mean-difference levels are the contrast rows times the draws", {

  object <- review_factor_replay(
    c("A", "B", "C"), "meandif", variable_intercept = FALSE, independent = TRUE
  )
  slopes <- .review_slope_draws(object)
  contrasts <- BayesTools::contr.meandif(3)
  table <- .review_level_table(object)

  expect_identical(rownames(table), c("g[A]", "g[B]", "g[C]"))
  for (i in 1:3) {
    level <- as.numeric(slopes %*% contrasts[i, ])
    expect_equal(table[i, "Mean"], mean(level), tolerance = 1e-12, info = i)
    expect_equal(
      unname(unlist(table[i, c("0.025", "0.5", "0.975")])),
      unname(stats::quantile(level, c(.025, .5, .975))),
      tolerance = 1e-12, info = i
    )
  }
})
