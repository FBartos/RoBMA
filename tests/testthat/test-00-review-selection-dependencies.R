test_that("F2-2 dependency components preserve sampling and random coefficient support", {

  data <- list(outcome = list(sei = rep(1, 4L)))
  plan <- function(C = diag(4L), loading = matrix(0, 4L, 0L), dense = FALSE, blocks = list()) {
    if (dense) list(representation = "dense", covariance = C) else
      list(representation = "diagonal_factor", diagonal = diag(C), loading = loading, supplied = list(blocks = blocks))
  }
  term <- list(block_name = "study", structure = "id", group_map = c(1L, 1L, 2L, 2L),
    n_columns = 2L, model_matrix = rbind(c(1, 0), c(0, 1), c(1, 0), c(0, 1)))
  known <- term; known$group_covariance <- list(type = "known", kernel = matrix(c(1, .3, .3, 1), 2L))
  overlap <- term; overlap$model_matrix[2L, 1L] <- 1
  diagonal <- known; diagonal$group_covariance$kernel <- diag(2L)
  crossed <- term; crossed$block_name <- "crossed"; crossed$group_map <- c(1L, 2L, 1L, 2L); crossed$structure <- "un"
  L <- cbind(c(1, 1, 0, 0), c(0, 1, 1, 0))
  C <- diag(4L); C[1L, 3L] <- C[3L, 1L] <- .2
  cases <- list(
    list(plan(), list(), diag(TRUE, 4L)),
    list(plan(loading = L), list(), tcrossprod(L) != 0 | diag(TRUE, 4L)),
    list(plan(C, dense = TRUE), list(), C != 0),
    list(plan(blocks = list(list(index = 1:4, covariance = C))), list(), C != 0),
    list(plan(), list(term), NULL), list(plan(), list(overlap), NULL),
    list(plan(), list(known), NULL), list(plan(), list(diagonal), NULL),
    list(plan(), list(term, crossed), NULL))
  for (case in cases) {
    adjacency <- case[[3L]]
    if (is.null(adjacency)) adjacency <- BayesTools::random_effects_dependency_matrix(case[[2L]], 4L)
    expected <- .known_v_block_indices(1 * adjacency)
    expect_identical(.selection_joint_dependency_blocks(data, case[[2L]], case[[1L]]), expected)
  }
  expect_identical(.selection_joint_dependency_blocks(data, list(known), plan()), list(c(1L, 3L), c(2L, 4L)))
  expect_identical(.selection_joint_dependency_blocks(data, list(term), plan()), as.list(1:4))
})

test_that("F2-2 independent preparation avoids the dense block finder", {

  testthat::local_mocked_bindings(.known_v_block_indices = function(...) stop("dense graph requested"), .package = "RoBMA")
  expect_identical(.selection_joint_dependency_blocks(list(outcome = list(sei = rep(1, 100L)))), as.list(seq_len(100L)))
})
