test_that("nested sample exports retain heterogeneous prediction intervals", {

  draws <- matrix(seq_len(40L), 20L, dimnames = list(NULL, c("a", "b")))
  plain <- .new_brma_samples(draws, 1L, 20L, "Plain")
  predictive <- .new_brma_samples(draws + .25, 1L, 20L, "Predictive",
    prediction_samples = draws * 2)
  plain_table <- as.data.frame(summary(plain))
  predictive_table <- as.data.frame(summary(predictive))
  columns <- c("component", "parameter", "Mean", "Median", "CI_0.025",
    "CI_0.975", "PI_0.025", "PI_0.975")
  plain_table[["PI_0.025"]] <- rep(NA_real_, 2L)
  plain_table[["PI_0.975"]] <- rep(NA_real_, 2L)
  plain_table[["component"]] <- "outer/plain"
  predictive_table[["component"]] <- "outer/predictive"
  for (pi_first in c(TRUE, FALSE)) {
    leaves <- if (pi_first) {
      list(predictive = predictive, plain = plain)
    } else list(plain = plain, predictive = predictive)
    result <- .new_brma_samples_list(list(outer = leaves))
    expected <- if (pi_first) {
      rbind(predictive_table[columns], plain_table[columns])
    } else rbind(plain_table[columns], predictive_table[columns])
    rownames(expected) <- NULL
    for (coerce in list(as.data.frame, data.frame)) {
      actual <- coerce(result)
      expect_identical(actual, expected)
      expect_identical(names(actual), columns)
      expect_identical(actual[["component"]],
        rep(paste0("outer/", names(leaves)), each = 2L))
      expect_identical(actual[["parameter"]], rep(c("a", "b"), 2L))
      expect_true(all(vapply(actual[-(1:2)], is.double, logical(1L))))
      expect_true(all(is.na(actual[actual$component == "outer/plain",
        c("PI_0.025", "PI_0.975")])))
      expect_equal(as.numeric(actual[actual$component == "outer/predictive",
        "PI_0.025"]), unname(apply(draws * 2, 2L, stats::quantile, probs = .025)))
    }
  }
})
