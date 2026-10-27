context("Random-effect rows without an original-scale value")

# A scaled predictor in factor-specific slopes of an unstructured block,
# us(1 + f:x | study): centering x induces level-specific intercepts that are
# not columns of the block, so the original-scale correlations of the intercept
# are not defined. BayesTools reports them as NA with a note; the RoBMA tables
# relay the rows and the note in RoBMA's own terms. Small fit: the expectations
# are properties of the tables, not of the posterior.
.random_na_rows_env <- new.env(parent = emptyenv())

.random_na_rows_data <- function() {

  set.seed(11)
  k_study <- 16L
  per     <- 4L
  dat <- data.frame(
    study = factor(rep(sprintf("s%02d", seq_len(k_study)), each = per)),
    x     = rnorm(k_study * per, mean = 10, sd = 2),
    z     = rnorm(k_study * per, mean = 3, sd = 1.5),
    f     = factor(rep(c("a", "b"), length.out = k_study * per))
  )
  dat[["vi"]] <- 0.02 + rexp(nrow(dat), 50)
  dat[["yi"]] <- 0.2 + 0.03 * (dat[["x"]] - 10) + 0.02 * dat[["z"]] +
    rnorm(nrow(dat), 0, 0.15) + rnorm(nrow(dat), 0, sqrt(dat[["vi"]]))

  return(dat)
}

# The fits by key: the unstructured block of the header, an exact block (no row
# without an original-scale value), and a diagonal block whose lower-order slope
# z is shifted by the factor-specific slopes of x:z.
.random_na_rows_fit <- function(key = "us_fx") {

  specs <- list(
    us_fx     = list(mods = ~ x,     random = ~ us(1 + f:x | study)),
    us_exact  = list(mods = ~ x + z, random = ~ us(1 + x + z | study)),
    diag_zfxz = list(mods = ~ x + z, random = ~ diag(1 + z + f:x:z | study))
  )
  if (is.null(.random_na_rows_env[[key]])) {
    spec <- specs[[key]]
    .random_na_rows_env[[key]] <- suppressWarnings(brma.mv(
      yi = yi, V = vi, mods = spec[["mods"]], random = spec[["random"]],
      data = .random_na_rows_data(), measure = "SMD", chains = 2, adapt = 100,
      burnin = 150, sample = 200, seed = 1, silent = TRUE, parallel = FALSE
    ))
  }

  return(.random_na_rows_env[[key]])
}

# The unavailable rows as the table of all draws gives them (the rule before the
# set was read from a few draws): catalog name -> row label.
.random_na_rows_reference <- function(object) {

  table <- .summary_estimates_table(
    object = object, probs = c(.025, .975), include_mcmc_diagnostics = FALSE,
    is_robma = .is_robust_RoBMA(object), keep_parameters = "random",
    transform_scaled = TRUE, random_effects_summary = "full", formula_prefix = FALSE
  )
  rows <- rownames(table)[is.na(table[, "Mean"]) &
    rownames(table) %in% names(attr(table, "footnotes", exact = TRUE))]

  return(stats::setNames(rows, attr(table, "parameters", exact = TRUE)[match(rows, rownames(table))]))
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


test_that("summary_heterogeneity() reports the rows that summary() shows as NA as NA with the same note", {

  skip_on_cran()
  fit      <- .random_na_rows_fit()
  random   <- summary(fit)[["estimates_random"]]
  result   <- summary_heterogeneity(fit)
  estimates <- result[["estimates"]]
  na_rows  <- rownames(estimates)[is.na(estimates[, "Mean"])]

  # the intercept of the block and its correlations have no original-scale
  # value; nothing else is lost
  expect_identical(na_rows, c(
    "tau(intercept)", "tau2(intercept)",
    "rho(intercept,f[a]:x)", "rho(intercept,f[b]:x)"
  ))
  expect_true(all(is.finite(as.matrix(estimates[setdiff(rownames(estimates), na_rows), ]))))
  expect_true(all(c("tau(f[a]:x)", "tau(f[b]:x)", "rho(f[a]:x,f[b]:x)") %in% rownames(estimates)))

  # every row that summary() reports as NA is NA here, with the identical note
  summary_na <- rownames(random)[is.na(random[, "Mean"])]
  expect_true(all(summary_na %in% na_rows))
  notes <- attr(estimates, "footnotes", exact = TRUE)
  expect_identical(names(notes), na_rows)
  expect_true(all(startsWith(notes, paste0(na_rows, ": NA: "))))
  summary_notes <- attr(random, "footnotes", exact = TRUE)
  expect_identical(notes[summary_na], summary_notes[summary_na])
  expect_true(all(endsWith(notes, " Use 'standardized_coefficients = TRUE' for fitted-scale values.")))

  printed <- utils::capture.output(print(result))
  expect_equal(sum(grepl(": NA: the original-scale value is unavailable", printed, fixed = TRUE)), 4L)

  # the data-frame export keeps the NA rows instead of values at the predictor mean
  exported <- as.data.frame(result)
  expect_true(all(is.na(exported[match(na_rows, exported[["parameter"]]), "Mean"])))
})


test_that("the unavailable rows are read from a few draws and equal those of the table of all draws", {

  skip_on_cran()
  expected <- list(
    us_fx     = c("(mu) sd(intercept)", "(mu) var(intercept)",
                  "(mu) cor(intercept,f[a]:x)", "(mu) cor(intercept,f[b]:x)"),
    us_exact  = character(),
    diag_zfxz = c("(mu) sd(intercept)", "(mu) var(intercept)", "(mu) sd(z)", "(mu) var(z)")
  )
  for (key in names(expected)) {
    fit       <- .random_na_rows_fit(key)
    unavailable <- .brma_mv_original_scale_unavailable(fit)
    reference <- .random_na_rows_reference(fit)
    expect_identical(sort(as.character(names(unavailable))), sort(expected[[key]]), info = key)
    expect_identical(sort(as.character(names(unavailable))),
                     sort(as.character(names(reference))), info = key)
    # one draw reads the same design
    expect_identical(sort(as.character(names(.brma_mv_original_scale_unavailable(fit, n_draws = 1L)))),
                     sort(expected[[key]]), info = key)
  }

  # the fit with an exact block (nothing unavailable) has no NA row in either table
  exact <- summary_heterogeneity(.random_na_rows_fit("us_exact"))[["estimates"]]
  expect_false(anyNA(as.matrix(exact)))
  expect_length(attr(exact, "footnotes", exact = TRUE), 0L)

  # the lower-order slope z of the diagonal block is shifted by f:x:z: its SD
  # is NA with the intercept's
  shifted <- summary_heterogeneity(.random_na_rows_fit("diag_zfxz"))[["estimates"]]
  expect_identical(rownames(shifted)[is.na(shifted[, "Mean"])],
                   c("tau(intercept)", "tau(z)", "tau2(intercept)", "tau2(z)"))
})
