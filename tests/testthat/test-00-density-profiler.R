context("Density profiler controls")


test_that("density profiler controls resolve with their unchanged ordinate budgets", {

  path <- testthat::test_path("..", "..", "tools", "profile-iwmde-qcmde.R")
  skip_if_not(file.exists(path), "Developer profiler requires a source checkout.")
  expressions <- parse(path, keep.source = FALSE)
  definition <- Filter(function(expr) {
    is.call(expr) && identical(expr[[1L]], as.name("<-")) &&
      identical(expr[[2L]], as.name(".profile_iwmde_estimator"))
  }, as.list(expressions))[[1L]][[3L]]
  find_estimator <- function(expr) {
    if (!is.call(expr)) return(NULL)
    if (identical(expr[[1L]], as.name(".iwmde_estimate"))) return(expr)
    for (part in as.list(expr)) {
      found <- find_estimator(part)
      if (!is.null(found)) return(found)
    }
    return(NULL)
  }
  estimator <- find_estimator(definition)
  controls <- eval(estimator[["density_control"]],
                   list(n_points = 40L, max_samples = 40L,
                        normalization_points = 101L))
  resolved <- tryCatch(RoBMA:::.iwmde_density_control_resolve(
    controls, density_method = "qCMDE", purpose = "ordinate"
  ), error = identity)
  expect_false(inherits(resolved, "error"),
               info = if (inherits(resolved, "error")) conditionMessage(resolved))
  if (!inherits(resolved, "error")) {
    expect_identical(resolved[["samples"]], 40L)
    expect_identical(resolved[["n_points"]], 40L)
    expect_identical(resolved[["normalization_points"]], 101L)
    expect_identical(resolved[["normalization_prob"]], .9999)
    expect_identical(resolved[["display_grid"]], "ordinate")
  }
})
