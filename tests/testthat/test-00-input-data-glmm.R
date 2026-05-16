context("Input handling for brma.glmm")

skip_on_cran()
source(testthat::test_path("helper-contracts.R"))

# Test data for GLMM input specification tests (2x2 contingency tables)
test_data_glmm <- data.frame(
  # Cell counts (ai = events in treatment, bi = non-events in treatment,
  #              ci = events in control, di = non-events in control)
  ai  = c(10L, 15L, 12L,  8L, 20L),
  bi  = c(40L, 35L, 38L, 42L, 30L),
  ci  = c( 5L, 10L,  8L,  4L, 12L),
  di  = c(45L, 40L, 42L, 46L, 38L),
  # Marginal totals
  n1i = c(50L, 50L, 50L, 50L, 50L),  # ai + bi
  n2i = c(50L, 50L, 50L, 50L, 50L),  # ci + di
  # Optional variables
  wgt     = c(1.0, 1.5, 1.2, 0.8, 1.3),
  study   = c("Study A", "Study B", "Study C", "Study D", "Study E"),
  cluster = c("g1", "g1", "g2", "g2", "g3"),
  stringsAsFactors = FALSE
)

expect_glmm_outcome <- function(result, n, ai, ci, n1i = NULL, n2i = NULL) {

  expect_type(result, "list")
  expect_true("outcome" %in% names(result))
  expect_equal(nrow(result$outcome), n)
  expect_equal(result$outcome$ai, ai)
  expect_equal(result$outcome$ci, ci)

  if (!is.null(n1i)) {
    expect_equal(result$outcome$n1i, n1i)
  }
  if (!is.null(n2i)) {
    expect_equal(result$outcome$n2i, n2i)
  }
}


# ============================================================================
# Basic input specification tests
# ============================================================================

test_that("GLMM input rejects unsupported effect-size measures directly", {

  expect_error(
    brma.glmm(
      ai        = c(10L, 15L),
      bi        = c(40L, 35L),
      ci        = c(5L, 10L),
      di        = c(45L, 40L),
      measure   = "SMD",
      only_data = TRUE
    ),
    "OR.*IRR"
  )

  expect_error(
    BMA.glmm(
      ai        = c(10L, 15L),
      bi        = c(40L, 35L),
      ci        = c(5L, 10L),
      di        = c(45L, 40L),
      measure   = "RR",
      only_data = TRUE
    ),
    "OR.*IRR"
  )
})

test_that("GLMM accepts supported 2x2 input specifications", {

  input_cases <- list(
    list(
      label = "direct cell vectors",
      expr  = quote(brma.glmm(
        ai        = c(10L, 15L, 12L),
        bi        = c(40L, 35L, 38L),
        ci        = c(5L, 10L, 8L),
        di        = c(45L, 40L, 42L),
        only_data = TRUE
      )[["data"]]),
      n    = 3,
      ai   = c(10L, 15L, 12L),
      ci   = c(5L, 10L, 8L),
      n1i  = c(50L, 50L, 50L),
      n2i  = c(50L, 50L, 50L)
    ),
    list(
      label = "events with arm totals",
      expr  = quote(brma.glmm(
        ai        = c(10L, 15L, 12L),
        ci        = c(5L, 10L, 8L),
        n1i       = c(50L, 50L, 50L),
        n2i       = c(50L, 50L, 50L),
        only_data = TRUE
      )[["data"]]),
      n    = 3,
      ai   = c(10L, 15L, 12L),
      ci   = c(5L, 10L, 8L),
      n1i  = c(50L, 50L, 50L),
      n2i  = c(50L, 50L, 50L)
    ),
    list(
      label = "unquoted data columns",
      expr  = quote(brma.glmm(
        ai        = ai,
        bi        = bi,
        ci        = ci,
        di        = di,
        data      = test_data_glmm,
        only_data = TRUE
      )[["data"]]),
      n    = 5,
      ai   = test_data_glmm$ai,
      ci   = test_data_glmm$ci,
      n1i  = test_data_glmm$n1i,
      n2i  = test_data_glmm$n2i
    ),
    list(
      label = "quoted data columns",
      expr  = quote(brma.glmm(
        ai        = "ai",
        bi        = "bi",
        ci        = "ci",
        di        = "di",
        data      = test_data_glmm,
        only_data = TRUE
      )[["data"]]),
      n    = 5,
      ai   = test_data_glmm$ai,
      ci   = test_data_glmm$ci,
      n1i  = test_data_glmm$n1i,
      n2i  = test_data_glmm$n2i
    ),
    list(
      label = "mixed quoted and unquoted columns",
      expr  = quote(brma.glmm(
        ai        = ai,
        bi        = "bi",
        ci        = ci,
        di        = di,
        data      = test_data_glmm,
        only_data = TRUE
      )[["data"]]),
      n    = 5,
      ai   = test_data_glmm$ai,
      ci   = test_data_glmm$ci,
      n1i  = test_data_glmm$n1i,
      n2i  = test_data_glmm$n2i
    )
  )

  for (case in input_cases) {
    expect_glmm_outcome(
      eval(case$expr),
      n    = case$n,
      ai   = case$ai,
      ci   = case$ci,
      n1i  = case$n1i,
      n2i  = case$n2i
    )
  }
})


# ============================================================================
# Optional arguments tests
# ============================================================================

test_that("GLMM handles slab and cluster metadata", {

  skip_on_cran()

  default_result <- brma.glmm(
    ai        = ai,
    bi        = bi,
    ci        = ci,
    di        = di,
    data      = test_data_glmm,
    only_data = TRUE
  )[["data"]]

  result <- brma.glmm(
    ai        = ai,
    bi        = bi,
    ci        = ci,
    di        = di,
    cluster   = cluster,
    slab      = study,
    data      = test_data_glmm,
    only_data = TRUE
  )[["data"]]

  expect_equal(default_result$outcome$slab, paste0("Study ", 1:5))
  expect_equal(result$outcome$slab, test_data_glmm$study)
  expect_true(is.numeric(result$outcome$cluster))
  expect_equal(result$outcome$cluster, c(1, 1, 2, 2, 3))
})


# ============================================================================
# Input validation tests
# ============================================================================

test_that("GLMM rejects missing required inputs", {

  skip_on_cran()

  expect_error_cases(list(
    list(
      label  = "missing control arm",
      expr   = quote(brma.glmm(
        ai        = c(10L, 15L),
        bi        = c(40L, 35L),
        only_data = TRUE
      )),
      regexp = "provide either"
    ),
    list(
      label  = "missing ci",
      expr   = quote(brma.glmm(
        ai        = c(10L, 15L),
        bi        = c(40L, 35L),
        di        = c(45L, 40L),
        only_data = TRUE
      )),
      regexp = "provide either"
    )
  ))
})


test_that("GLMM validates cell counts are non-negative integers", {

  skip_on_cran()

  expect_error_cases(list(
    list(
      label  = "negative ai",
      expr   = quote(brma.glmm(
        ai        = c(10L, -5L),
        bi        = c(40L, 35L),
        ci        = c(5L, 10L),
        di        = c(45L, 40L),
        only_data = TRUE
      )),
      regexp = "ai"
    ),
    list(
      label  = "fractional ai",
      expr   = quote(brma.glmm(
        ai        = c(10, 15.5),
        bi        = c(40L, 35L),
        ci        = c(5L, 10L),
        di        = c(45L, 40L),
        only_data = TRUE
      )),
      regexp = "ai"
    )
  ))
})


test_that("GLMM validates redundant 2x2 totals against supplied cells", {

  skip_on_cran()

  result <- brma.glmm(
    ai        = c(10L, 15L),
    bi        = c(40L, 35L),
    ci        = c(5L, 10L),
    di        = c(45L, 40L),
    n1i       = c(50L, 50L),
    n2i       = c(50L, 50L),
    only_data = TRUE
  )[["data"]]

  expect_equal(result$outcome$n1i, c(50L, 50L))
  expect_equal(result$outcome$n2i, c(50L, 50L))

  expect_error_cases(list(
    list(
      label  = "bad supplied n1i",
      expr   = quote(brma.glmm(
        ai        = c(10L, 15L),
        bi        = c(40L, 35L),
        ci        = c(5L, 10L),
        di        = c(45L, 40L),
        n1i       = c(50L, 100L),
        n2i       = c(50L, 50L),
        only_data = TRUE
      )),
      regexp = "n1i.*ai [+] bi"
    ),
    list(
      label  = "bad supplied n2i",
      expr   = quote(brma.glmm(
        ai        = c(10L, 15L),
        bi        = c(40L, 35L),
        ci        = c(5L, 10L),
        di        = c(45L, 40L),
        n1i       = c(50L, 50L),
        n2i       = c(50L, 100L),
        only_data = TRUE
      )),
      regexp = "n2i.*ci [+] di"
    ),
    list(
      label  = "bad inferred bi",
      expr   = quote(brma.glmm(
        ai        = c(10L, 15L),
        bi        = c(40L, 99L),
        ci        = c(5L, 10L),
        n1i       = c(50L, 50L),
        n2i       = c(50L, 50L),
        only_data = TRUE
      )),
      regexp = "n1i.*ai [+] bi"
    ),
    list(
      label  = "bad inferred di",
      expr   = quote(brma.glmm(
        ai        = c(10L, 15L),
        ci        = c(5L, 10L),
        di        = c(45L, 99L),
        n1i       = c(50L, 50L),
        n2i       = c(50L, 50L),
        only_data = TRUE
      )),
      regexp = "n2i.*ci [+] di"
    )
  ))
})


test_that("GLMM input preserves likelihood weights", {

  skip_on_cran()

  bin_result <- brma.glmm(
    ai        = ai,
    bi        = bi,
    ci        = ci,
    di        = di,
    weights   = wgt,
    data      = test_data_glmm,
    only_data = TRUE
  )[["data"]]

  expect_true(attr(bin_result, "weights"))
  expect_equal(bin_result$outcome$weights, test_data_glmm$wgt)

  pois_result <- brma.glmm(
    x1i       = c(3L, 5L),
    x2i       = c(2L, 4L),
    t1i       = c(10, 12),
    t2i       = c(11, 13),
    weights   = c(1, 2),
    measure   = "IRR",
    only_data = TRUE
  )[["data"]]

  expect_true(attr(pois_result, "weights"))
  expect_equal(pois_result$outcome$weights, c(1, 2))

  expect_error(
    brma.glmm(
      ai        = ai,
      bi        = bi,
      ci        = ci,
      di        = di,
      weights   = c(1, 0, 1, 1, 1),
      data      = test_data_glmm,
      only_data = TRUE
    ),
    regexp = "weights"
  )
})


test_that("GLMM validates inferred cells and argument lengths", {

  skip_on_cran()

  expect_error_cases(list(
    list(
      label  = "computed bi is negative",
      expr   = quote(brma.glmm(
        ai        = c(10L, 60L),
        ci        = c(5L, 10L),
        n1i       = c(50L, 50L),
        n2i       = c(50L, 50L),
        only_data = TRUE
      )),
      regexp = "bi.*negative"
    ),
    list(
      label  = "computed di is negative",
      expr   = quote(brma.glmm(
        ai        = c(10L, 15L),
        ci        = c(5L, 60L),
        n1i       = c(50L, 50L),
        n2i       = c(50L, 50L),
        only_data = TRUE
      )),
      regexp = "di.*negative"
    ),
    list(
      label  = "cell-count length mismatch",
      expr   = quote(brma.glmm(
        ai        = c(10L, 15L, 12L),
        bi        = c(40L, 35L),
        ci        = c(5L, 10L, 8L),
        di        = c(45L, 40L, 42L),
        only_data = TRUE
      )),
      regexp = "bi"
    ),
    list(
      label  = "redundant-total length mismatch",
      expr   = quote(brma.glmm(
        ai        = c(10L, 15L, 12L),
        bi        = c(40L, 35L),
        ci        = c(5L, 10L, 8L),
        n1i       = c(50L, 50L, 50L),
        n2i       = c(50L, 50L, 50L),
        only_data = TRUE
      )),
      regexp = "bi"
    ),
    list(
      label  = "zero treatment arm total",
      expr   = quote(brma.glmm(
        ai        = c(0L, 10L),
        ci        = c(1L, 10L),
        n1i       = c(0L, 50L),
        n2i       = c(50L, 50L),
        only_data = TRUE
      )),
      regexp = "n1i.*positive"
    ),
    list(
      label  = "zero control arm total",
      expr   = quote(brma.glmm(
        ai        = c(1L, 10L),
        ci        = c(0L, 10L),
        n1i       = c(50L, 50L),
        n2i       = c(0L, 50L),
        only_data = TRUE
      )),
      regexp = "n2i.*positive"
    ),
    list(
      label  = "zero treatment arm total from cells",
      expr   = quote(brma.glmm(
        ai        = c(0L, 10L),
        bi        = c(0L, 40L),
        ci        = c(1L, 10L),
        di        = c(49L, 40L),
        only_data = TRUE
      )),
      regexp = "n1i.*positive"
    )
  ))
})


# ============================================================================
# Subset tests
# ============================================================================

test_that("GLMM applies and validates subsets", {

  skip_on_cran()

  subset_cases <- list(
    list(subset = c(TRUE, TRUE, FALSE, TRUE, FALSE), rows = c(1, 2, 4)),
    list(subset = c(1, 3, 5), rows = c(1, 3, 5))
  )

  for (case in subset_cases) {
    result <- brma.glmm(
      ai        = ai,
      bi        = bi,
      ci        = ci,
      di        = di,
      subset    = case$subset,
      data      = test_data_glmm,
      only_data = TRUE
    )[["data"]]

    expect_equal(nrow(result$outcome), length(case$rows))
    expect_equal(result$outcome$ai, test_data_glmm$ai[case$rows])
    expect_equal(result$outcome$ci, test_data_glmm$ci[case$rows])
  }

  expect_error(
    brma.glmm(
      ai     = ai,
      bi     = bi,
      ci     = ci,
      di     = di,
      subset = c(TRUE, NA, TRUE, TRUE, TRUE),
      data   = test_data_glmm,
      only_data = TRUE
    ),
    regexp = "subset.*NA"
  )
})


# ============================================================================
# NA handling tests
# ============================================================================

test_that("GLMM drops NA cell rows and refreshes derived metadata", {

  skip_on_cran()

  test_data_na <- data.frame(
    ai = c(10L, 15L, 12L, 8L),
    bi = c(40L, 35L, 38L, 42L),
    ci = c(5L, 10L, 8L, 4L),
    di = c(45L, 40L, 42L, 46L)
  )
  test_data_na[["ai"]][2] <- NA_integer_

  expect_warning(
    result <- brma.glmm(
      ai        = ai,
      bi        = bi,
      ci        = ci,
      di        = di,
      data      = test_data_na,
      only_data = TRUE
    )[["data"]],
    regexp = "1 observation.*removed"
  )

  expect_equal(nrow(result$outcome), 3)
  expect_equal(result$outcome$ai, c(10L, 12L, 8L))

  test_data_multi_na <- data.frame(
    ai = c(10L, 15L, 12L, 8L),
    bi = c(40L, 35L, 38L, 42L),
    ci = c(5L, 10L, 8L, 4L),
    di = c(45L, 40L, 42L, 46L)
  )
  test_data_multi_na[["bi"]][2] <- NA_integer_
  test_data_multi_na[["ci"]][3] <- NA_integer_
  test_data_multi_na[["di"]][4] <- NA_integer_

  expect_warning(
    result <- brma.glmm(
      ai        = ai,
      bi        = bi,
      ci        = ci,
      di        = di,
      data      = test_data_multi_na,
      only_data = TRUE
    )[["data"]],
    regexp = "3 observation.*removed"
  )

  expect_equal(nrow(result$outcome), 1)
  expect_equal(result$outcome$ai, 10L)

  test_data_na <- data.frame(
    ai  = c(10L, 15L, 12L),
    bi  = c(40L, NA, 38L),
    ci  = c(5L, 10L, 8L),
    di  = c(45L, 40L, 42L),
    n1i = c(50L, 50L, 50L),
    n2i = c(50L, 50L, 50L)
  )

  expect_warning(
    result <- brma.glmm(
      ai        = ai,
      bi        = bi,
      ci        = ci,
      di        = di,
      n1i       = n1i,
      n2i       = n2i,
      data      = test_data_na,
      only_data = TRUE
    )[["data"]],
    regexp = "1 observation.*removed"
  )

  expect_equal(nrow(result$outcome), 2)
  expect_equal(result$outcome$ai, c(10L, 12L))
  expect_equal(result$outcome$n1i, c(50L, 50L))

  test_data_na <- data.frame(
    ai = c(10L, 15L, 12L),
    bi = c(40L, 35L, 38L),
    ci = c(5L, 10L, 8L),
    di = c(45L, 40L, 42L)
  )
  test_data_na[["ai"]][2] <- NA_integer_

  expect_warning(
    result <- brma.glmm(
      ai   = ai,
      bi   = bi,
      ci   = ci,
      di   = di,
      data = test_data_na,
      only_data = TRUE
    )[["data"]],
    regexp = "removed"
  )

  expect_equal(result$outcome$slab, c("Study 1", "Study 2"))

  test_data_na <- data.frame(
    ai      = c(10L, 15L, 12L, 8L),
    bi      = c(40L, 35L, 38L, 42L),
    ci      = c(5L, 10L, 8L, 4L),
    di      = c(45L, 40L, 42L, 46L),
    cluster = c("g1", "g1", "g2", "g2")
  )
  test_data_na[["ai"]][2] <- NA_integer_

  expect_warning(
    result <- brma.glmm(
      ai        = ai,
      bi        = bi,
      ci        = ci,
      di        = di,
      cluster   = cluster,
      data      = test_data_na,
      only_data = TRUE
    )[["data"]],
    regexp = "removed"
  )

  expect_equal(result$outcome$cluster, c(1, 2, 2))
})


test_that("GLMM rejects NA cluster, but preserves NA slab", {

  skip_on_cran()

  test_data_na <- data.frame(
    ai      = c(10L, 15L, 12L, 8L),
    bi      = c(40L, 35L, 38L, 42L),
    ci      = c(5L, 10L, 8L, 4L),
    di      = c(45L, 40L, 42L, 46L),
    cluster = c("g1", "g1", "g2", "g3"),
    slab    = c("A", NA, "C", "D")
  )

  test_data_na[["cluster"]][2] <- NA_character_

  expect_error(
    brma.glmm(
      ai        = ai,
      bi        = bi,
      ci        = ci,
      di        = di,
      cluster   = cluster,
      data      = test_data_na,
      only_data = TRUE
    ),
    regexp = "argument must not contain missing values"
  )

  test_data_na[["cluster"]][2] <- "g1"

  result <- brma.glmm(
    ai        = ai,
    bi        = bi,
    ci        = ci,
    di        = di,
    cluster   = cluster,
    slab      = slab,
    data      = test_data_na,
    only_data = TRUE
  )[["data"]]

  expect_equal(nrow(result$outcome), 4)
  expect_true(is.na(result$outcome$slab[2]))
})


test_that("GLMM rejects all-NA observations", {

  skip_on_cran()

  test_data_all_na <- data.frame(
    ai = c(1L, 2L),
    bi = c(40L, 35L),
    ci = c(5L, 10L),
    di = c(45L, 40L)
  )
  test_data_all_na[["ai"]] <- c(NA_integer_, NA_integer_)

  expect_error(
    suppressWarnings(
      brma.glmm(
        ai   = ai,
        bi   = bi,
        ci   = ci,
        di   = di,
        data = test_data_all_na,
        only_data = TRUE
      )
    ),
    regexp = "No observations remaining"
  )
})


# ============================================================================
# Moderator (mods) tests
# ============================================================================

test_data_glmm_mods <- data.frame(
  ai         = c(10L, 15L, 12L, 8L, 20L),
  bi         = c(40L, 35L, 38L, 42L, 30L),
  ci         = c(5L, 10L, 8L, 4L, 12L),
  di         = c(45L, 40L, 42L, 46L, 38L),
  mod_cont   = c(1.5, 2.3, 1.8, 3.1, 0.9),
  mod_factor = factor(c("A", "B", "A", "B", "A")),
  scale_var  = c(0.5, 1.0, 0.8, 1.2, 0.6),
  stringsAsFactors = FALSE
)


test_that("GLMM parses moderator and scale formulas", {

  skip_on_cran()

  formula_cases <- list(
    list(
      expr        = quote(brma.glmm(
        ai        = ai,
        bi        = bi,
        ci        = ci,
        di        = di,
        mods      = ~ mod_cont + mod_factor,
        data      = test_data_glmm_mods,
        only_data = TRUE
      )[["data"]]),
      mods_names  = c("mod_cont", "mod_factor"),
      scale_names = character()
    ),
    list(
      expr        = quote(brma.glmm(
        ai        = ai,
        bi        = bi,
        ci        = ci,
        di        = di,
        scale     = ~ scale_var,
        data      = test_data_glmm_mods,
        only_data = TRUE
      )[["data"]]),
      mods_names  = character(),
      scale_names = "scale_var"
    ),
    list(
      expr        = quote(brma.glmm(
        ai        = ai,
        bi        = bi,
        ci        = ci,
        di        = di,
        mods      = ~ mod_cont + mod_factor,
        scale     = ~ scale_var,
        data      = test_data_glmm_mods,
        only_data = TRUE
      )[["data"]]),
      mods_names  = c("mod_cont", "mod_factor"),
      scale_names = "scale_var"
    )
  )

  for (case in formula_cases) {
    result <- eval(case$expr)

    if (length(case$mods_names)) {
      expect_equal(nrow(result$mods), 5)
      expect_true(all(case$mods_names %in% names(result$mods)))
    } else {
      expect_null(result$mods)
    }

    if (length(case$scale_names)) {
      expect_equal(nrow(result$scale), 5)
      expect_true(all(case$scale_names %in% names(result$scale)))
    } else {
      expect_null(result$scale)
    }
  }
})


test_that("GLMM drops rows with NA in mods", {

  skip_on_cran()

  test_data_na <- data.frame(
    ai       = c(10L, 15L, 12L, 8L),
    bi       = c(40L, 35L, 38L, 42L),
    ci       = c(5L, 10L, 8L, 4L),
    di       = c(45L, 40L, 42L, 46L),
    mod_cont = c(1.5, NA, 1.8, 3.1)
  )

  expect_warning(
    result <- brma.glmm(
      ai   = ai,
      bi   = bi,
      ci   = ci,
      di   = di,
      mods = ~ mod_cont,
      data = test_data_na,
      only_data = TRUE
    )[["data"]],
    regexp = "1 observation.*removed"
  )

  expect_equal(nrow(result$outcome), 3)
  expect_equal(nrow(result$mods), 3)
})


# ============================================================================
# Predictor validation tests
# ============================================================================

test_that("GLMM rejects degenerate moderator design matrices", {

  skip_on_cran()

  test_data_const <- data.frame(
    ai        = c(10L, 15L, 12L),
    bi        = c(40L, 35L, 38L),
    ci        = c(5L, 10L, 8L),
    di        = c(45L, 40L, 42L),
    const_var = c(5, 5, 5)
  )

  test_data_single <- data.frame(
    ai         = c(10L, 15L, 12L),
    bi         = c(40L, 35L, 38L),
    ci         = c(5L, 10L, 8L),
    di         = c(45L, 40L, 42L),
    single_fac = factor(c("A", "A", "A"))
  )

  expect_error_cases(list(
    list(
      label  = "constant numeric moderator",
      expr   = quote(brma.glmm(
        ai        = ai,
        bi        = bi,
        ci        = ci,
        di        = di,
        mods      = ~ const_var,
        data      = test_data_const,
        only_data = TRUE
      )),
      regexp = "zero variance"
    ),
    list(
      label  = "single-level factor moderator",
      expr   = quote(brma.glmm(
        ai        = ai,
        bi        = bi,
        ci        = ci,
        di        = di,
        mods      = ~ single_fac,
        data      = test_data_single,
        only_data = TRUE
      )),
      regexp = "only one level"
    ),
    list(
      label  = "subset drops all but one factor level",
      expr   = quote(brma.glmm(
        ai        = ai,
        bi        = bi,
        ci        = ci,
        di        = di,
        mods      = ~ mod_factor,
        subset    = c(TRUE, FALSE, TRUE, FALSE, TRUE),
        data      = test_data_glmm_mods,
        only_data = TRUE
      )),
      regexp = "only one level"
    )
  ))
})


# ============================================================================
# RoBMA_data class and attributes tests
# ============================================================================

test_that("GLMM returns RoBMA_data class with correct attributes", {

  skip_on_cran()

  result <- brma.glmm(
    ai        = ai,
    bi        = bi,
    ci        = ci,
    di        = di,
    cluster   = cluster,
    slab      = study,
    data      = test_data_glmm,
    only_data = TRUE
  )[["data"]]

  expect_s3_class(result, "RoBMA_data")
  expect_equal(attr(result, "k_final"), 5)
  expect_equal(attr(result, "n_dropped"), 0)
  expect_false(attr(result, "mods"))
  expect_false(attr(result, "scale"))
  expect_false(attr(result, "weights"))
  expect_true(attr(result, "slab"))
  expect_true(attr(result, "cluster"))
})


test_that("GLMM attributes reflect NA dropping", {

  skip_on_cran()

  test_data_na <- data.frame(
    ai = c(10L, 15L, 12L, 8L),
    bi = c(40L, 35L, 38L, 42L),
    ci = c(5L, 10L, 8L, 4L),
    di = c(45L, 40L, 42L, 46L)
  )
  test_data_na[["ai"]][2] <- NA_integer_

  expect_warning(
    result <- brma.glmm(
      ai   = ai,
      bi   = bi,
      ci   = ci,
      di   = di,
      data = test_data_na,
      only_data = TRUE
    )[["data"]],
    regexp = "removed"
  )

  expect_equal(attr(result, "k_final"), 3)
  expect_equal(attr(result, "n_dropped"), 1)
})
