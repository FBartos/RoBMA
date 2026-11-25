test_that("cluster comparison keys preserve ordered row membership", {

  make_score <- function(blocks) {
    score <- .add_cluster_log_lik_metadata(matrix(0, 5L, 2L), blocks, "same-data")
    structure(list(), class = "loo", RoBMA_target = attr(score, "RoBMA_target"))
  }
  first <- make_score(list(a = 1:2, b = 3:4))
  renamed <- make_score(list(x = 1:2, y = 3:4))
  changed <- make_score(list(a = c(1L, 3L), b = c(2L, 4L)))
  reversed <- make_score(list(a = 3:4, b = 1:2))
  expect_silent(.check_loo_compare_targets(list(first, renamed)))
  expect_error(.check_loo_compare_targets(list(first, changed)),
               "different cluster partitions", fixed = TRUE)
  expect_error(.check_loo_compare_targets(list(first, reversed)),
               "different cluster partitions", fixed = TRUE)
  attr(changed, "RoBMA_target")$cluster_partition <- NULL
  expect_error(.check_loo_compare_targets(list(first, changed)),
               "without cluster partitions", fixed = TRUE)
})

test_that("cached cluster scores reject changed partitions without changing outcome identity", {

  object <- brma(yi = c(0.1, 0.2, 0.3, 0.4), sei = rep(0.2, 4),
    cluster = c("a", "a", "b", "b"), measure = "SMD",
    only_priors = TRUE, silent = TRUE)
  key <- .current_predictive_target_key(object, "cluster")
  expect_silent(.check_cached_predictive_target(object, key, "cluster", "LOO"))
  object$data$outcome$cluster <- c(1L, 2L, 1L, 2L)
  expect_identical(.get_outcome_hash(object), key$data_hash)
  expect_error(.check_cached_predictive_target(object, key, "cluster", "LOO"),
               "current cluster partition", fixed = TRUE)
})
