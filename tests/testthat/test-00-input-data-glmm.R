context("Input handling for brma.glmm")

skip_on_cran()

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


# ============================================================================
# Basic input specification tests
# ============================================================================

test_that("GLMM input works with direct vectors (ai, bi, ci, di)", {

  result <- brma.glmm(
    ai = c(10L, 15L, 12L),
    bi = c(40L, 35L, 38L),
    ci = c(5L, 10L, 8L),
    di = c(45L, 40L, 42L),
    only_data = TRUE
  )[["data"]]

  expect_type(result, "list")
  expect_true("outcome" %in% names(result))
  expect_equal(nrow(result$outcome), 3)
  expect_equal(result$outcome$ai, c(10L, 15L, 12L))
  expect_equal(result$outcome$ci, c(5L, 10L, 8L))
  # n1i and n2i should be computed
  expect_equal(result$outcome$n1i, c(50L, 50L, 50L))
  expect_equal(result$outcome$n2i, c(50L, 50L, 50L))
})


test_that("GLMM input works with alternative format (ai, ci, n1i, n2i)", {

  skip_on_cran()

  result <- brma.glmm(
    ai  = c(10L, 15L, 12L),
    ci  = c(5L, 10L, 8L),
    n1i = c(50L, 50L, 50L),
    n2i = c(50L, 50L, 50L),
    only_data = TRUE
  )[["data"]]

  expect_type(result, "list")
  expect_equal(nrow(result$outcome), 3)
  expect_equal(result$outcome$ai, c(10L, 15L, 12L))
  # bi should be computed as n1i - ai
  expect_equal(result$outcome$ci, c(5L, 10L, 8))
  expect_equal(result$outcome$n1i, c(50L, 50L, 50L))
  expect_equal(result$outcome$n2i, c(50L, 50L, 50L))
})


test_that("GLMM input works with unquoted column names from data", {

  skip_on_cran()

  result <- brma.glmm(
    ai   = ai,
    bi   = bi,
    ci   = ci,
    di   = di,
    data = test_data_glmm,
    only_data = TRUE
  )[["data"]]

  expect_type(result, "list")
  expect_equal(nrow(result$outcome), 5)
  expect_equal(result$outcome$ai, test_data_glmm$ai)
  expect_equal(result$outcome$ci, test_data_glmm$ci)
})


test_that("GLMM input works with quoted string column names from data", {

  skip_on_cran()

  result <- brma.glmm(
    ai   = "ai",
    bi   = "bi",
    ci   = "ci",
    di   = "di",
    data = test_data_glmm,
    only_data = TRUE
  )[["data"]]

  expect_type(result, "list")
  expect_equal(nrow(result$outcome), 5)
  expect_equal(result$outcome$ai, test_data_glmm$ai)
  expect_equal(result$outcome$ci, test_data_glmm$ci)
})


test_that("GLMM input works with mixed specification styles", {

  skip_on_cran()

  external_weights <- c(2, 2, 2, 2, 2)

  result <- brma.glmm(
    ai      = ai,           # unquoted column name
    bi      = "bi",         # quoted string column name
    ci      = ci,
    di      = di,
    weights = external_weights,  # direct vector from environment
    data    = test_data_glmm,
    only_data = TRUE
  )[["data"]]

  expect_type(result, "list")
  expect_equal(nrow(result$outcome), 5)
  expect_equal(result$outcome$ai, test_data_glmm$ai)
  expect_equal(result$outcome$weights, external_weights)
})


# ============================================================================
# Optional arguments tests
# ============================================================================

test_that("GLMM input works with optional arguments", {

  skip_on_cran()

  result <- brma.glmm(
    ai        = ai,
    bi        = bi,
    ci        = ci,
    di        = di,
    weights   = wgt,
    study_ids = cluster,
    slab      = study,
    data      = test_data_glmm,
    only_data = TRUE
  )[["data"]]

  expect_equal(result$outcome$weights, test_data_glmm$wgt)
  expect_equal(result$outcome$slab, test_data_glmm$study)
  # study_ids should be converted to numeric factor
  expect_true(is.numeric(result$outcome$study_ids))
})


test_that("GLMM generates default slab when not provided", {

  skip_on_cran()

  result <- brma.glmm(
    ai   = ai,
    bi   = bi,
    ci   = ci,
    di   = di,
    data = test_data_glmm,
    only_data = TRUE
  )[["data"]]

  expect_equal(result$outcome$slab, paste0("Study ", 1:5))
})


test_that("GLMM correctly handles study_ids", {

  skip_on_cran()

  result <- brma.glmm(
    ai        = ai,
    bi        = bi,
    ci        = ci,
    di        = di,
    study_ids = cluster,
    data      = test_data_glmm,
    only_data = TRUE
  )[["data"]]

  # study_ids should be converted to numeric indices
  expect_true(is.numeric(result$outcome$study_ids))
  # g1, g1, g2, g2, g3 should become 1, 1, 2, 2, 3
  expect_equal(result$outcome$study_ids, c(1, 1, 2, 2, 3))
})


# ============================================================================
# Input validation tests
# ============================================================================

test_that("GLMM throws error when required inputs are missing", {

  skip_on_cran()

  # Missing ai, bi, ci, di and n1i, n2i
  expect_error(
    brma.glmm(
      ai = c(10L, 15L),
      bi = c(40L, 35L),
      only_data = TRUE
    ),
    regexp = "provide either"
  )

  # Missing ci when using ai, bi
  expect_error(
    brma.glmm(
      ai = c(10L, 15L),
      bi = c(40L, 35L),
      di = c(45L, 40L),
      only_data = TRUE
    ),
    regexp = "provide either"
  )
})


test_that("GLMM validates cell counts are non-negative integers", {

  skip_on_cran()

  # Negative values should error
  expect_error(
    brma.glmm(
      ai = c(10L, -5L),
      bi = c(40L, 35L),
      ci = c(5L, 10L),
      di = c(45L, 40L),
      only_data = TRUE
    ),
    regexp = "ai"
  )
})


test_that("GLMM validates bi computed from n1i - ai is non-negative", {

  skip_on_cran()

  # ai > n1i should error
  expect_error(
    brma.glmm(
      ai  = c(10L, 60L),  # 60 > 50
      ci  = c(5L, 10L),
      n1i = c(50L, 50L),
      n2i = c(50L, 50L),
      only_data = TRUE
    ),
    regexp = "bi.*negative"
  )
})


test_that("GLMM validates di computed from n2i - ci is non-negative", {

  skip_on_cran()

  # ci > n2i should error
  expect_error(
    brma.glmm(
      ai  = c(10L, 15L),
      ci  = c(5L, 60L),  # 60 > 50
      n1i = c(50L, 50L),
      n2i = c(50L, 50L),
      only_data = TRUE
    ),
    regexp = "di.*negative"
  )
})


test_that("GLMM validates consistent lengths", {

  skip_on_cran()

  expect_error(
    brma.glmm(
      ai = c(10L, 15L, 12L),
      bi = c(40L, 35L),  # Wrong length
      ci = c(5L, 10L, 8L),
      di = c(45L, 40L, 42L),
      only_data = TRUE
    ),
    regexp = "bi"
  )
})


# ============================================================================
# Subset tests
# ============================================================================

test_that("GLMM applies logical subset correctly", {

  skip_on_cran()

  result <- brma.glmm(
    ai     = ai,
    bi     = bi,
    ci     = ci,
    di     = di,
    subset = c(TRUE, TRUE, FALSE, TRUE, FALSE),
    data   = test_data_glmm,
    only_data = TRUE
  )[["data"]]

  expect_equal(nrow(result$outcome), 3)
  expect_equal(result$outcome$ai, test_data_glmm$ai[c(1, 2, 4)])
  expect_equal(result$outcome$ci, test_data_glmm$ci[c(1, 2, 4)])
})


test_that("GLMM applies numeric subset correctly", {

  skip_on_cran()

  result <- brma.glmm(
    ai     = ai,
    bi     = bi,
    ci     = ci,
    di     = di,
    subset = c(1, 3, 5),
    data   = test_data_glmm,
    only_data = TRUE
  )[["data"]]

  expect_equal(nrow(result$outcome), 3)
  expect_equal(result$outcome$ai, test_data_glmm$ai[c(1, 3, 5)])
})


test_that("GLMM subset rejects NA values", {

  skip_on_cran()

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

test_that("GLMM drops rows with NA in cell counts", {

  skip_on_cran()
  # length mismatch when NA values are present
  # Create data frame with integer NA

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
    regexp = "1 observation.*removed"
  )

  expect_equal(nrow(result$outcome), 3)
  expect_equal(result$outcome$ai, c(10L, 12L, 8L))
})


test_that("GLMM drops rows with NA in any cell count column", {

  skip_on_cran()

  # Create data frame with integer NAs in different columns
  test_data_na <- data.frame(
    ai = c(10L, 15L, 12L, 8L),
    bi = c(40L, 35L, 38L, 42L),
    ci = c(5L, 10L, 8L, 4L),
    di = c(45L, 40L, 42L, 46L)
  )
  test_data_na[["bi"]][2] <- NA_integer_
  test_data_na[["ci"]][3] <- NA_integer_
  test_data_na[["di"]][4] <- NA_integer_

  expect_warning(
    result <- brma.glmm(
      ai   = ai,
      bi   = bi,
      ci   = ci,
      di   = di,
      data = test_data_na,
      only_data = TRUE
    )[["data"]],
    regexp = "3 observation.*removed"
  )

  expect_equal(nrow(result$outcome), 1)
  expect_equal(result$outcome$ai, 10L)
})


test_that("GLMM regenerates slab after NA dropping", {

  skip_on_cran()

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

  # slab should be regenerated with correct numbering
  expect_equal(result$outcome$slab, c("Study 1", "Study 2"))
})


test_that("GLMM renumbers study_ids after NA dropping", {

  skip_on_cran()

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
      study_ids = cluster,
      data      = test_data_na,
      only_data = TRUE
    )[["data"]],
    regexp = "removed"
  )

  # study_ids should be renumbered: g1, g2, g2 -> 1, 2, 2
  expect_equal(result$outcome$study_ids, c(1, 2, 2))
})


test_that("GLMM does not drop rows with NA in weights/study_ids/slab", {

  skip_on_cran()

  test_data_na <- data.frame(
    ai      = c(10L, 15L, 12L, 8L),
    bi      = c(40L, 35L, 38L, 42L),
    ci      = c(5L, 10L, 8L, 4L),
    di      = c(45L, 40L, 42L, 46L),
    wgt     = c(1.0, NA, 1.2, 0.8)
  )

  # Should not trigger warning about NA removal
  result <- brma.glmm(
    ai      = ai,
    bi      = bi,
    ci      = ci,
    di      = di,
    weights = wgt,
    data    = test_data_na,
    only_data = TRUE
  )[["data"]]

  expect_equal(nrow(result$outcome), 4)
  expect_true(is.na(result$outcome$weights[2]))
})


test_that("GLMM throws error when all observations are NA", {

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


test_that("GLMM mods works with formula", {

  skip_on_cran()

  result <- brma.glmm(
    ai   = ai,
    bi   = bi,
    ci   = ci,
    di   = di,
    mods = ~ mod_cont + mod_factor,
    data = test_data_glmm_mods,
    only_data = TRUE
  )[["data"]]

  expect_true(!is.null(result$mods))
  expect_equal(nrow(result$mods), 5)
  expect_true("mod_cont" %in% names(result$mods))
  expect_true("mod_factor" %in% names(result$mods))
})


test_that("GLMM scale works with formula", {

  skip_on_cran()

  result <- brma.glmm(
    ai    = ai,
    bi    = bi,
    ci    = ci,
    di    = di,
    scale = ~ scale_var,
    data  = test_data_glmm_mods,
    only_data = TRUE
  )[["data"]]

  expect_true(!is.null(result$scale))
  expect_equal(nrow(result$scale), 5)
  expect_true("scale_var" %in% names(result$scale))
})


test_that("GLMM mods and scale work together", {

  skip_on_cran()

  result <- brma.glmm(
    ai    = ai,
    bi    = bi,
    ci    = ci,
    di    = di,
    mods  = ~ mod_cont + mod_factor,
    scale = ~ scale_var,
    data  = test_data_glmm_mods,
    only_data = TRUE
  )[["data"]]

  expect_true(!is.null(result$mods))
  expect_true(!is.null(result$scale))
  expect_equal(nrow(result$mods), 5)
  expect_equal(nrow(result$scale), 5)
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

test_that("GLMM constant mods variable throws error", {

  skip_on_cran()

  test_data_const <- data.frame(
    ai        = c(10L, 15L, 12L),
    bi        = c(40L, 35L, 38L),
    ci        = c(5L, 10L, 8L),
    di        = c(45L, 40L, 42L),
    const_var = c(5, 5, 5)
  )

  expect_error(
    brma.glmm(
      ai   = ai,
      bi   = bi,
      ci   = ci,
      di   = di,
      mods = ~ const_var,
      data = test_data_const,
      only_data = TRUE
    ),
    regexp = "zero variance"
  )
})


test_that("GLMM single-level factor in mods throws error", {

  skip_on_cran()

  test_data_single <- data.frame(
    ai         = c(10L, 15L, 12L),
    bi         = c(40L, 35L, 38L),
    ci         = c(5L, 10L, 8L),
    di         = c(45L, 40L, 42L),
    single_fac = factor(c("A", "A", "A"))
  )

  expect_error(
    brma.glmm(
      ai   = ai,
      bi   = bi,
      ci   = ci,
      di   = di,
      mods = ~ single_fac,
      data = test_data_single,
      only_data = TRUE
    ),
    regexp = "only one level"
  )
})


test_that("GLMM subset that reduces factor to single level throws error", {

  skip_on_cran()

  expect_error(
    brma.glmm(
      ai     = ai,
      bi     = bi,
      ci     = ci,
      di     = di,
      mods   = ~ mod_factor,
      subset = c(TRUE, FALSE, TRUE, FALSE, TRUE),  # Only level "A" remains
      data   = test_data_glmm_mods,
      only_data = TRUE
    ),
    regexp = "only one level"
  )
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
    weights   = wgt,
    study_ids = cluster,
    slab      = study,
    data      = test_data_glmm,
    only_data = TRUE
  )[["data"]]

  expect_s3_class(result, "RoBMA_data")
  expect_equal(attr(result, "k_final"), 5)
  expect_equal(attr(result, "n_dropped"), 0)
  expect_false(attr(result, "mods"))
  expect_false(attr(result, "scale"))
  expect_true(attr(result, "weights"))
  expect_true(attr(result, "slab"))
  expect_true(attr(result, "study_ids"))
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
