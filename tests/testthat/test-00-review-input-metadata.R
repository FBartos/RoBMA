test_that("synthetic clusters do not overwrite a moderator column", {

  original <- data.frame(yi = c(.1, .2, .3), sei = rep(.1, 3),
                         group = c("a", "a", "b"), .RoBMA_cluster = 1:3)
  object <- brma(yi = yi, sei = sei, cluster = group, mods = ~ .RoBMA_cluster,
                  data = original, measure = "SMD", only_priors = TRUE,
                  silent = TRUE)
  prepared <- .prepare_newdata(object,
    data.frame(.RoBMA_cluster = c(99, 101)), type = "terms")
  expect_identical(prepared[["mods"]][[".RoBMA_cluster"]], c(99, 101))
  expect_equal(length(unique(prepared[["outcome"]][["cluster"]])), 2L)
})

