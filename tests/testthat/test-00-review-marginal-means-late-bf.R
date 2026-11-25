context("Late Bayes factor requests on qCMDE/IWMDE marginal means")

test_that("marginal means built without Bayes factors name the missing ordinates", {

  object <- review_factor_replay(c("A", "B", "C"), "meandif")

  for (method in c("qCMDE", "IWMDE")) {
    means <- marginal_means(object, bf = FALSE, density_method = method)
    expect_false(means[["bf"]])

    result  <- summary(means, bf = TRUE)
    reasons <- unlist(attr(result, "warnings", exact = TRUE))

    # the Bayes factors stay unavailable, with the reason that the ordinates
    # were never computed and the supported call, not a failed diagnostic
    expect_true(all(is.na(result[["inclusion_BF"]])), info = method)
    expect_gt(length(reasons), 0L)
    expect_true(all(grepl("built with 'bf = FALSE'", reasons, fixed = TRUE)), info = method)
    expect_true(all(grepl(
      paste0("marginal_means(object, bf = TRUE, density_method = \"", method, "\")"),
      reasons, fixed = TRUE
    )), info = method)
    expect_false(any(grepl("failed diagnostics", reasons, fixed = TRUE)), info = method)
  }
})
