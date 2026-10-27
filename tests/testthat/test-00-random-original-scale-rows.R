context("Random-effect rows without an original-scale value")

# A scaled predictor in factor-specific slopes of an unstructured block,
# us(1 + f:x | study): centering x induces level-specific intercepts that are
# not columns of the block, so the original-scale correlations of the intercept
# are not defined. BayesTools reports them as NA with a note; the RoBMA tables
# relay the rows and the note in RoBMA's own terms. Small fit: the expectations
# are properties of the tables, not of the posterior.
.random_na_rows_env <- new.env(parent = emptyenv())

.random_na_rows_fit <- function() {

  if (is.null(.random_na_rows_env[["fit"]])) {
    set.seed(11)
    k_study <- 12L
    per     <- 4L
    dat <- data.frame(
      study = factor(rep(sprintf("s%02d", seq_len(k_study)), each = per)),
      x     = rnorm(k_study * per, mean = 10, sd = 2),
      f     = factor(rep(c("a", "b"), length.out = k_study * per))
    )
    dat[["vi"]] <- 0.02 + rexp(nrow(dat), 50)
    dat[["yi"]] <- 0.2 + 0.03 * (dat[["x"]] - 10) +
      rnorm(nrow(dat), 0, 0.15) + rnorm(nrow(dat), 0, sqrt(dat[["vi"]]))

    .random_na_rows_env[["fit"]] <- suppressWarnings(brma.mv(
      yi = yi, V = vi, mods = ~ x, random = ~ us(1 + f:x | study),
      data = dat, measure = "SMD", chains = 2, adapt = 100, burnin = 150,
      sample = 200, seed = 1, silent = TRUE, parallel = FALSE
    ))
  }

  return(.random_na_rows_env[["fit"]])
}


test_that("summary() translates the BayesTools scale advice to standardized_coefficients", {

  skip_on_cran()
  fit      <- .random_na_rows_fit()
  original <- summary(fit)
  rows     <- original[["estimates_random"]]
  na_rows  <- rownames(rows)[is.na(rows[, "Mean"])]

  expect_identical(na_rows, c("rho(intercept,f[a]:x)", "rho(intercept,f[b]:x)"))

  notes <- attr(rows, "footnotes", exact = TRUE)
  expect_identical(names(notes), na_rows)
  expect_true(all(startsWith(notes, paste0(na_rows, ": NA: "))))
  expect_true(all(endsWith(notes, " Use 'standardized_coefficients = TRUE' for fitted-scale values.")))
  expect_false(any(grepl("transform_scaled", notes, fixed = TRUE)))

  printed <- utils::capture.output(print(fit))
  expect_false(any(grepl("transform_scaled", printed, fixed = TRUE)))
  expect_equal(sum(grepl("Use 'standardized_coefficients = TRUE'", printed, fixed = TRUE)), 2L)

  # the stored object summary is not printed, but it carries the same note
  stored <- attr(fit[["summary"]], "footnotes", exact = TRUE)
  expect_length(stored, 2L)
  expect_false(any(grepl("transform_scaled", stored, fixed = TRUE)))
  expect_true(all(grepl("standardized_coefficients = TRUE", stored, fixed = TRUE)))

  # the advice works: the fitted-scale table has every value and no note
  standardized <- summary(fit, standardized_coefficients = TRUE)[["estimates_random"]]
  expect_false(anyNA(standardized[, "Mean"]))
  expect_true(all(na_rows %in% rownames(standardized)))
  expect_length(attr(standardized, "footnotes", exact = TRUE), 0L)
})
