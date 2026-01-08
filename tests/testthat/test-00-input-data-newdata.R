context("Input handling for .prepare_newdata")

# ============================================================================
# Test data setup
# ============================================================================

# Test data for normal likelihood models
test_data_norm <- data.frame(
  yi        = c(0.2, 0.5, -0.1, 0.3, 0.4),
  sei       = c(0.1, 0.15, 0.12, 0.08, 0.11),
  mod_cont  = c(1.5, 2.3, 1.8, 3.1, 0.9),
  mod_fac   = factor(c("A", "B", "A", "B", "A")),
  scale_var = c(0.5, 1.0, 0.8, 1.2, 0.6),
  stringsAsFactors = FALSE
)

# Test data with vi instead of sei
test_data_norm_vi <- data.frame(
  yi        = c(0.2, 0.5, -0.1, 0.3, 0.4),
  vi        = c(0.01, 0.0225, 0.0144, 0.0064, 0.0121),
  mod_cont  = c(1.5, 2.3, 1.8, 3.1, 0.9),
  mod_fac   = factor(c("A", "B", "A", "B", "A")),
  stringsAsFactors = FALSE
)

# Test data for GLMM binomial models
test_data_glmm <- data.frame(
  ai        = c(10L, 15L, 12L, 8L, 20L),
  ci        = c(5L, 10L, 8L, 4L, 12L),
  n1i       = c(50L, 50L, 50L, 50L, 50L),
  n2i       = c(50L, 50L, 50L, 50L, 50L),
  mod_cont  = c(1.5, 2.3, 1.8, 3.1, 0.9),
  mod_fac   = factor(c("A", "B", "A", "B", "A")),
  stringsAsFactors = FALSE
)


# ============================================================================
# Helper function to compare data lists (ignoring slab differences)
# ============================================================================

compare_data_lists <- function(data1, data2, check_slab = FALSE, structure_only = FALSE) {
  # Check class
  expect_s3_class(data1, "RoBMA_data")
  expect_s3_class(data2, "RoBMA_data")

  # Check outcome dimensions
  if (!structure_only) {
    expect_equal(nrow(data1$outcome), nrow(data2$outcome))
  }
  expect_equal(ncol(data1$outcome), ncol(data2$outcome))
  expect_equal(names(data1$outcome), names(data2$outcome))

  # Check outcome column values (excluding slab which may differ)
  if (!structure_only) {
    outcome_cols <- setdiff(names(data1$outcome), if (check_slab) character(0) else "slab")
    for (col in outcome_cols) {
      expect_equal(
        data1$outcome[[col]], data2$outcome[[col]],
        info = paste("Mismatch in outcome column:", col)
      )
    }
  }

  # Check mods if present
  if (!is.null(data1$mods) && !is.null(data2$mods)) {
    if (!structure_only) {
      expect_equal(nrow(data1$mods), nrow(data2$mods))
    }
    expect_equal(names(data1$mods), names(data2$mods))
    if (!structure_only) {
      for (col in names(data1$mods)) {
        expect_equal(
          data1$mods[[col]], data2$mods[[col]],
          info = paste("Mismatch in mods column:", col)
        )
      }
    }
    # Check formula attribute
    expect_equal(
      attr(data1$mods, "formula"), attr(data2$mods, "formula"),
      info = "Mismatch in mods formula attribute"
    )
  } else {
    expect_equal(is.null(data1$mods), is.null(data2$mods))
  }

  # Check scale if present
  if (!is.null(data1$scale) && !is.null(data2$scale)) {
    if (!structure_only) {
      expect_equal(nrow(data1$scale), nrow(data2$scale))
    }
    expect_equal(names(data1$scale), names(data2$scale))
    if (!structure_only) {
      for (col in names(data1$scale)) {
        expect_equal(
          data1$scale[[col]], data2$scale[[col]],
          info = paste("Mismatch in scale column:", col)
        )
      }
    }
    expect_equal(
      attr(data1$scale, "formula"), attr(data2$scale, "formula"),
      info = "Mismatch in scale formula attribute"
    )
  } else {
    expect_equal(is.null(data1$scale), is.null(data2$scale))
  }

  # Check key attributes
  expect_equal(attr(data1, "outcome_type"), attr(data2, "outcome_type"))
  if (!structure_only) {
    expect_equal(attr(data1, "k_final"), attr(data2, "k_final"))
  }
  expect_equal(attr(data1, "mods"), attr(data2, "mods"))
  expect_equal(attr(data1, "scale"), attr(data2, "scale"))
  expect_equal(
    attr(data1, "standardize_continuous_predictors"),
    attr(data2, "standardize_continuous_predictors")
  )
  expect_equal(
    attr(data1, "set_contrast_factor_predictors"),
    attr(data2, "set_contrast_factor_predictors")
  )
  expect_equal(attr(data1, "effect_direction"), attr(data2, "effect_direction"))
}


# ============================================================================
# Basic brma.norm tests - outcome variables
# ============================================================================

test_that(".prepare_newdata works with brma.norm using sei", {

  skip_on_cran()

  # Create a fitted object (data only)
  fit <- brma.norm(
    yi   = yi,
    sei  = sei,
    data = test_data_norm,
    only_data = TRUE
  )

  # Prepare newdata using the same data
  newdata_result <- RoBMA:::.prepare_newdata(
    object                       = fit,
    newdata                      = test_data_norm,
    type                         = "terms",
    incorporate_publication_bias = TRUE
  )

  # Compare with the original data
  compare_data_lists(fit[["data"]], newdata_result)
})


test_that(".prepare_newdata works with brma.norm using vi", {

  skip_on_cran()

  # Create a fitted object using vi
  fit <- brma.norm(
    yi   = test_data_norm_vi$yi,
    sei  = sqrt(test_data_norm_vi$vi),
    only_data = TRUE
  )

  # Prepare newdata using the same data
  newdata_result <- RoBMA:::.prepare_newdata(
    object                       = fit,
    newdata                      = test_data_norm_vi,
    type                         = "terms",
    incorporate_publication_bias = TRUE
  )

  # Compare with the original data
  compare_data_lists(fit[["data"]], newdata_result)
})


test_that(".prepare_newdata works with different newdata (same structure)", {

  skip_on_cran()

  # Create a fitted object
  fit <- brma.norm(
    yi   = yi,
    sei  = sei,
    data = test_data_norm,
    only_data = TRUE
  )

  # Create new data with different values but same structure
  new_df <- data.frame(
    yi  = c(0.1, 0.6, 0.0),
    sei = c(0.2, 0.1, 0.15)
  )

  # Prepare using .prepare_newdata
  newdata_result <- RoBMA:::.prepare_newdata(
    object                       = fit,
    newdata                      = new_df,
    type                         = "terms",
    incorporate_publication_bias = TRUE
  )

  # Compare
  compare_data_lists(fit[["data"]], newdata_result, structure_only = TRUE)
})


# ============================================================================
# brma.norm with moderators (mods)
# ============================================================================

test_that(".prepare_newdata works with brma.norm and mods formula", {

  skip_on_cran()

  # Create a fitted object with mods
  fit <- brma.norm(
    yi   = yi,
    sei  = sei,
    mods = ~ mod_cont + mod_fac,
    data = test_data_norm,
    only_data = TRUE
  )

  # Prepare newdata
  newdata_result <- RoBMA:::.prepare_newdata(
    object                       = fit,
    newdata                      = test_data_norm,
    type                         = "terms",
    incorporate_publication_bias = TRUE
  )

  # Compare with the original data
  compare_data_lists(fit[["data"]], newdata_result)
})


test_that(".prepare_newdata with mods and different newdata", {

  skip_on_cran()

  # Create a fitted object with mods
  fit <- brma.norm(
    yi   = yi,
    sei  = sei,
    mods = ~ mod_cont + mod_fac,
    data = test_data_norm,
    only_data = TRUE
  )

  # Create new data with different values
  new_df <- data.frame(
    yi       = c(0.1, 0.3),
    sei      = c(0.1, 0.2),
    mod_cont = c(2.0, 1.0),
    mod_fac  = factor(c("A", "B"), levels = c("A", "B"))
  )

  # Prepare using .prepare_newdata
  newdata_result <- RoBMA:::.prepare_newdata(
    object                       = fit,
    newdata                      = new_df,
    type                         = "terms",
    incorporate_publication_bias = TRUE
  )

  # Compare
  compare_data_lists(fit[["data"]], newdata_result, structure_only = TRUE)
})


# ============================================================================
# brma.norm with scale predictors
# ============================================================================

test_that(".prepare_newdata works with brma.norm and scale formula", {

  skip_on_cran()

  # Create a fitted object with scale
  fit <- brma.norm(
    yi    = yi,
    sei   = sei,
    scale = ~ scale_var,
    data  = test_data_norm,
    only_data = TRUE
  )

  # Prepare newdata
  newdata_result <- RoBMA:::.prepare_newdata(
    object                       = fit,
    newdata                      = test_data_norm,
    type                         = "terms",
    incorporate_publication_bias = TRUE
  )

  # Compare with the original data
  compare_data_lists(fit[["data"]], newdata_result)
})


test_that(".prepare_newdata works with both mods and scale", {

  skip_on_cran()

  # Create a fitted object with both mods and scale
  fit <-  brma.norm(
    yi    = yi,
    sei   = sei,
    mods  = ~ mod_cont,
    scale = ~ scale_var,
    data  = test_data_norm,
    only_data = TRUE
  )

  # Prepare newdata
  newdata_result <- RoBMA:::.prepare_newdata(
    object                       = fit,
    newdata                      = test_data_norm,
    type                         = "terms",
    incorporate_publication_bias = TRUE
  )

  # Compare with the original data
  compare_data_lists(fit[["data"]], newdata_result)
})


# ============================================================================
# brma.glmm tests (binomial)
# ============================================================================

test_that(".prepare_newdata works with brma.glmm (binomial)", {

  skip_on_cran()

  # Create a fitted object
  fit <-  brma.glmm(
    ai   = ai,
    ci   = ci,
    n1i  = n1i,
    n2i  = n2i,
    data = test_data_glmm,
    only_data = TRUE
  )

  # Prepare newdata
  newdata_result <- RoBMA:::.prepare_newdata(
    object                       = fit,
    newdata                      = test_data_glmm,
    type                         = "terms",
    incorporate_publication_bias = TRUE
  )

  # Compare with the original data
  compare_data_lists(fit[["data"]], newdata_result)
})


test_that(".prepare_newdata works with brma.glmm and mods", {

  skip_on_cran()

  # Create a fitted object with mods
  fit <-  brma.glmm(
    ai   = ai,
    ci   = ci,
    n1i  = n1i,
    n2i  = n2i,
    mods = ~ mod_cont + mod_fac,
    data = test_data_glmm,
    only_data = TRUE
  )

  # Prepare newdata
  newdata_result <- RoBMA:::.prepare_newdata(
    object                       = fit,
    newdata                      = test_data_glmm,
    type                         = "terms",
    incorporate_publication_bias = TRUE
  )

  # Compare with the original data
  compare_data_lists(fit[["data"]], newdata_result)
})


# ============================================================================
# Input validation tests - missing columns
# ============================================================================

test_that(".prepare_newdata throws error when yi is missing", {

  skip_on_cran()

  fit <-  brma.norm(
    yi   = yi,
    sei  = sei,
    data = test_data_norm,
    only_data = TRUE
  )

  # newdata without yi
  bad_df <- data.frame(sei = c(0.1, 0.2))

  expect_error(
    RoBMA:::.prepare_newdata(
      object                       = fit,
      newdata                      = bad_df,
      type                         = "terms",
      incorporate_publication_bias = TRUE
    ),
    regexp = "yi"
  )
})


test_that(".prepare_newdata throws error when sei/vi is missing", {

  skip_on_cran()

  fit <-  brma.norm(
    yi   = yi,
    sei  = sei,
    data = test_data_norm,
    only_data = TRUE
  )

  # newdata without sei or vi
  bad_df <- data.frame(yi = c(0.1, 0.2))

  expect_error(
    RoBMA:::.prepare_newdata(
      object                       = fit,
      newdata                      = bad_df,
      type                         = "terms",
      incorporate_publication_bias = TRUE
    ),
    regexp = "sei.*vi"
  )
})


test_that(".prepare_newdata throws error when mods variables are missing", {

  skip_on_cran()

  fit <-  brma.norm(
    yi   = yi,
    sei  = sei,
    mods = ~ mod_cont + mod_fac,
    data = test_data_norm,
    only_data = TRUE
  )

  # newdata without mod_cont
  bad_df <- data.frame(
    yi      = c(0.1, 0.2),
    sei     = c(0.1, 0.2),
    mod_fac = factor(c("A", "B"))
  )

  expect_error(
    RoBMA:::.prepare_newdata(
      object                       = fit,
      newdata                      = bad_df,
      type                         = "terms",
      incorporate_publication_bias = TRUE
    ),
    regexp = "mod_cont"
  )
})


test_that(".prepare_newdata throws error when scale variables are missing", {

  skip_on_cran()

  fit <-  brma.norm(
    yi    = yi,
    sei   = sei,
    scale = ~ scale_var,
    data  = test_data_norm,
    only_data = TRUE
  )

  # newdata without scale_var
  bad_df <- data.frame(
    yi  = c(0.1, 0.2),
    sei = c(0.1, 0.2)
  )

  expect_error(
    RoBMA:::.prepare_newdata(
      object                       = fit,
      newdata                      = bad_df,
      type                         = "terms",
      incorporate_publication_bias = TRUE
    ),
    regexp = "scale_var"
  )
})


test_that(".prepare_newdata throws error when GLMM columns are missing", {

  skip_on_cran()

  fit <-  brma.glmm(
    ai   = ai,
    ci   = ci,
    n1i  = n1i,
    n2i  = n2i,
    data = test_data_glmm,
    only_data = TRUE
  )

  # newdata missing n2i
  bad_df <- data.frame(
    ai  = c(10L, 15L),
    ci  = c(5L, 10L),
    n1i = c(50L, 50L)
  )

  expect_error(
    RoBMA:::.prepare_newdata(
      object                       = fit,
      newdata                      = bad_df,
      type                         = "terms",
      incorporate_publication_bias = TRUE
    ),
    regexp = "n2i"
  )
})


# ============================================================================
# Settings preservation tests
# ============================================================================

test_that(".prepare_newdata preserves standardize_continuous_predictors setting", {

  skip_on_cran()

  # Create fit with standardize = TRUE
  fit_std <- brma.norm(
    yi   = yi,
    sei  = sei,
    mods = ~ mod_cont,
    data = test_data_norm,
    standardize_continuous_predictors = TRUE,
    only_data = TRUE
  )

  # Create fit with standardize = FALSE
  fit_no_std <- brma.norm(
    yi   = yi,
    sei  = sei,
    mods = ~ mod_cont,
    data = test_data_norm,
    standardize_continuous_predictors = FALSE,
    only_data = TRUE
  )

  # Prepare newdata for both
  result_std <- RoBMA:::.prepare_newdata(
    object                       = fit_std,
    newdata                      = test_data_norm,
    type                         = "terms",
    incorporate_publication_bias = TRUE
  )

  result_no_std <- RoBMA:::.prepare_newdata(
    object                       = fit_no_std,
    newdata                      = test_data_norm,
    type                         = "terms",
    incorporate_publication_bias = TRUE
  )

  # Check attributes are preserved

  expect_true(attr(result_std, "standardize_continuous_predictors"))
  expect_false(attr(result_no_std, "standardize_continuous_predictors"))
})


test_that(".prepare_newdata preserves set_contrast_factor_predictors setting", {

  skip_on_cran()

  # Create fit with treatment contrasts
  fit_treatment <- brma.norm(
    yi   = yi,
    sei  = sei,
    mods = ~ mod_fac,
    data = test_data_norm,
    set_contrast_factor_predictors = "treatment",
    only_data = TRUE
  )

  # Create fit with meandif contrasts
  fit_meandif <- brma.norm(
    yi   = yi,
    sei  = sei,
    mods = ~ mod_fac,
    data = test_data_norm,
    set_contrast_factor_predictors = "meandif",
    only_data = TRUE
  )

  # Prepare newdata for both
  result_treatment <- RoBMA:::.prepare_newdata(
    object                       = fit_treatment,
    newdata                      = test_data_norm,
    type                         = "terms",
    incorporate_publication_bias = TRUE
  )

  result_meandif <- RoBMA:::.prepare_newdata(
    object                       = fit_meandif,
    newdata                      = test_data_norm,
    type                         = "terms",
    incorporate_publication_bias = TRUE
  )

  # Check attributes are preserved
  expect_equal(attr(result_treatment, "set_contrast_factor_predictors"), "treatment")
  expect_equal(attr(result_meandif, "set_contrast_factor_predictors"), "meandif")
})


# ============================================================================
# NA handling tests
# ============================================================================

test_that(".prepare_newdata handles NA values in newdata", {

  skip_on_cran()

  fit <-  brma.norm(
    yi   = yi,
    sei  = sei,
    data = test_data_norm,
    only_data = TRUE
  )

  # newdata with NA
  new_df <- data.frame(
    yi  = c(0.1, NA, 0.3),
    sei = c(0.1, 0.2, 0.15)
  )

  # Should drop the NA row with a warning
  expect_warning(
    result <- RoBMA:::.prepare_newdata(
      object                       = fit,
      newdata                      = new_df,
      type                         = "terms",
      incorporate_publication_bias = TRUE
    ),
    regexp = "removed"
  )

  expect_equal(nrow(result$outcome), 2)
})


test_that(".prepare_newdata handles NA in mods", {

  skip_on_cran()

  fit <-  brma.norm(
    yi   = yi,
    sei  = sei,
    mods = ~ mod_cont,
    data = test_data_norm,
    only_data = TRUE
  )

  # newdata with NA in moderator
  new_df <- data.frame(
    yi       = c(0.1, 0.2, 0.3),
    sei      = c(0.1, 0.2, 0.15),
    mod_cont = c(1.5, NA, 2.0)
  )

  # Should drop the NA row
  expect_warning(
    result <- RoBMA:::.prepare_newdata(
      object                       = fit,
      newdata                      = new_df,
      type                         = "terms",
      incorporate_publication_bias = TRUE
    ),
    regexp = "removed"
  )

  expect_equal(nrow(result$outcome), 2)
  expect_equal(nrow(result$mods), 2)
})


# ============================================================================
# Edge cases
# ============================================================================

test_that(".prepare_newdata works with single row newdata", {

  skip_on_cran()

  fit <-  brma.norm(
    yi   = yi,
    sei  = sei,
    data = test_data_norm,
    only_data = TRUE
  )

  # Single row newdata
  new_df <- data.frame(yi = 0.5, sei = 0.1)

  result <- RoBMA:::.prepare_newdata(
    object                       = fit,
    newdata                      = new_df,
    type                         = "terms",
    incorporate_publication_bias = TRUE
  )

  expect_equal(nrow(result$outcome), 1)
  expect_equal(result$outcome$yi, 0.5)
  expect_equal(result$outcome$sei, 0.1)
})


test_that(".prepare_newdata works with newdata having extra columns", {

  skip_on_cran()

  fit <-  brma.norm(
    yi   = yi,
    sei  = sei,
    data = test_data_norm,
    only_data = TRUE
  )

  # newdata with extra columns
  new_df <- data.frame(
    yi        = c(0.1, 0.2),
    sei       = c(0.1, 0.2),
    extra_col = c("a", "b"),
    another   = c(1, 2)
  )

  # Should work, ignoring extra columns
  result <- RoBMA:::.prepare_newdata(
    object                       = fit,
    newdata                      = new_df,
    type                         = "terms",
    incorporate_publication_bias = TRUE
  )

  expect_equal(nrow(result$outcome), 2)
})
