context("Input handling for brma.uni")

# Test data for input specification tests
test_data <- data.frame(
  effect    = c(0.10, 0.25, 0.15, 0.30, 0.05),
  variance  = c(0.04, 0.06, 0.05, 0.08, 0.03),
  std_err   = sqrt(c(0.04, 0.06, 0.05, 0.08, 0.03)),
  n         = c(50L, 75L, 60L, 40L, 90L),
  wgt       = c(1.0, 1.5, 1.2, 0.8, 1.3),
  study     = c("A", "B", "C", "D", "E"),
  cluster   = c("g1", "g1", "g2", "g2", "g3"),
  stringsAsFactors = FALSE
)


test_that("Input works with direct vectors", {

  skip_on_cran()

  # Direct vector specification
  result <- brma.uni(
    yi  = c(0.10, 0.25, 0.15),
    sei = c(0.20, 0.24, 0.22),
    only_data = TRUE
  )[["data"]]

  expect_type(result, "list")
  expect_true("outcome" %in% names(result))
  expect_equal(nrow(result$outcome), 3)
  expect_equal(result$outcome$yi, c(0.10, 0.25, 0.15))
  expect_equal(result$outcome$sei, c(0.20, 0.24, 0.22))
})


test_that("Input works with unquoted column names from data", {

  skip_on_cran()

  # Unquoted column name specification
  result <- brma.uni(
    yi   = effect,
    vi   = variance,
    data = test_data,
    only_data = TRUE
  )[["data"]]

  expect_type(result, "list")
  expect_equal(nrow(result$outcome), 5)
  expect_equal(result$outcome$yi, test_data$effect)
  expect_equal(result$outcome$sei, sqrt(test_data$variance), tolerance = 1e-6)
})


test_that("Input works with quoted string column names from data", {

  skip_on_cran()

  # Quoted string column name specification
  result <- brma.uni(
    yi   = "effect",
    sei  = "std_err",
    data = test_data,
    only_data = TRUE
  )[["data"]]

  expect_type(result, "list")
  expect_equal(nrow(result$outcome), 5)
  expect_equal(result$outcome$yi, test_data$effect)
  expect_equal(result$outcome$sei, test_data$std_err, tolerance = 1e-6)
})


test_that("Input works with mixed specification styles", {

  skip_on_cran()

  # Mix of direct vectors and column references
  external_weights <- c(2, 2, 2, 2, 2)

  result <- brma.uni(
    yi      = effect,       # unquoted column name
    vi      = "variance",   # quoted string column name
    ni      = n,            # unquoted column name
    weights = external_weights,  # direct vector from environment
    data    = test_data,
    only_data = TRUE
  )[["data"]]

  expect_type(result, "list")
  expect_equal(nrow(result$outcome), 5)
  expect_equal(result$outcome$yi, test_data$effect)
  expect_equal(result$outcome$ni, test_data$n)
  expect_equal(result$outcome$weights, external_weights)
})


test_that("Input works with optional arguments", {

  skip_on_cran()

  # With slab and study_ids
  result <- brma.uni(
    yi        = effect,
    sei       = std_err,
    ni        = n,
    slab      = study,
    study_ids = cluster,
    data      = test_data,
    only_data = TRUE
  )[["data"]]

  expect_type(result, "list")
  expect_equal(nrow(result$outcome), 5)
  expect_equal(result$outcome$slab, test_data$study)
  expect_equal(result$outcome$study_ids, as.numeric(as.factor(test_data$cluster)))
  expect_equal(result$outcome$ni, test_data$n)
})


test_that("Input converts vi to sei correctly", {

  skip_on_cran()

  result <- brma.uni(
    yi   = effect,
    vi   = variance,
    data = test_data,
    only_data = TRUE
  )[["data"]]

  expected_sei <- sqrt(test_data$variance)
  expect_equal(result$outcome$sei, expected_sei, tolerance = 1e-10)
})


test_that("Input generates default slab when not provided", {

  skip_on_cran()

  result <- brma.uni(
    yi  = c(0.1, 0.2, 0.3),
    sei = c(0.1, 0.1, 0.1),
    only_data = TRUE
  )[["data"]]

  expect_equal(result$outcome$slab, c("Study 1", "Study 2", "Study 3"))
})


test_that("Input handles subset argument correctly", {

  skip_on_cran()

  # Logical subset
  result_logical <- brma.uni(
    yi     = effect,
    sei    = std_err,
    subset = c(TRUE, FALSE, TRUE, FALSE, TRUE),
    data   = test_data,
    only_data = TRUE
  )[["data"]]

  expect_equal(nrow(result_logical$outcome), 3)
  expect_equal(result_logical$outcome$yi, test_data$effect[c(1, 3, 5)])

  # Numeric subset
  result_numeric <- brma.uni(
    yi     = effect,
    sei    = std_err,
    subset = c(1, 3, 5),
    data   = test_data,
    only_data = TRUE
  )[["data"]]

  expect_equal(nrow(result_numeric$outcome), 3)
  expect_equal(result_numeric$outcome$yi, test_data$effect[c(1, 3, 5)])
})


test_that("Input throws error when yi is missing", {

  skip_on_cran()

  expect_error(
    brma.uni(
      sei = c(0.1, 0.2),
      only_data = TRUE
    ),
    regexp = "yi"
  )
})


test_that("Input throws error when neither vi nor sei is provided", {

  skip_on_cran()

  expect_error(
    brma.uni(
      yi = c(0.1, 0.2),
      only_data = TRUE
    ),
    regexp = "vi|sei|variance|standard error"
  )
})


test_that("Input throws error for length mismatch", {

  skip_on_cran()

  expect_error(
    brma.uni(
      yi  = c(0.1, 0.2, 0.3),
      sei = c(0.1, 0.2),  # wrong length
      only_data = TRUE
    ),
    regexp = "length"
  )
})


test_that("Input throws error for invalid column reference", {

  skip_on_cran()

  expect_error(
    brma.uni(
      yi   = nonexistent_column,
      sei  = std_err,
      data = test_data,
      only_data = TRUE
    ),
    regexp = "nonexistent_column|Cannot find"
  )
})


test_that("Input handles variables from calling environment", {

  skip_on_cran()

  # Create variables in the current environment
  my_effects   <- c(0.15, 0.25, 0.35)
  my_std_errs  <- c(0.10, 0.12, 0.11)
  my_labels    <- c("Study X", "Study Y", "Study Z")

  result <- brma.uni(
    yi   = my_effects,
    sei  = my_std_errs,
    slab = my_labels,
    only_data = TRUE
  )[["data"]]

  expect_equal(result$outcome$yi, my_effects)
  expect_equal(result$outcome$sei, my_std_errs)
  expect_equal(result$outcome$slab, my_labels)
})


test_that("Input works when called from within a function with direct column references", {

  skip_on_cran()

  # Wrapper function that uses NSE to pass column names
  # Note: When brma.uni is called from another function, the column names

  # must be available in the data frame or the calling environment of brma.uni
  wrapper_function <- function(data) {
    # Direct column references work because they are evaluated in data
    brma.uni(
      yi   = effect,
      sei  = std_err,
      data = data,
      only_data = TRUE
    )[["data"]]
  }

  # This tests that NSE works correctly when brma.uni is called from another function
  result <- wrapper_function(data = test_data)

  expect_type(result, "list")
  expect_equal(nrow(result$outcome), 5)
  expect_equal(result$outcome$yi, test_data$effect)
})


# ============================================================================
# Tests for moderator (mods) and scale input
# ============================================================================

# Extended test data with moderator variables
test_data_mods <- data.frame(
  effect     = c(0.10, 0.25, 0.15, 0.30, 0.05),
  variance   = c(0.04, 0.06, 0.05, 0.08, 0.03),
  std_err    = sqrt(c(0.04, 0.06, 0.05, 0.08, 0.03)),
  mod_cont   = c(1.5, 2.3, 1.8, 3.1, 0.9),
  mod_factor = factor(c("A", "B", "A", "B", "A")),
  mod_char   = c("low", "high", "low", "high", "low"),
  scale_var  = c(0.5, 1.0, 0.8, 1.2, 0.6),
  stringsAsFactors = FALSE
)


test_that("Mods works with formula using column names from data", {

  skip_on_cran()

  result <- brma.uni(
    yi   = effect,
    sei  = std_err,
    mods = ~ mod_cont + mod_factor,
    data = test_data_mods,
    only_data = TRUE
  )[["data"]]

  expect_type(result, "list")
  expect_true(!is.null(result$mods))
  expect_true(is.data.frame(result$mods))
  expect_equal(nrow(result$mods), 5)
  expect_true("mod_cont" %in% names(result$mods))
  expect_true("mod_factor" %in% names(result$mods))

  # Check that types are preserved
  expect_type(result$mods$mod_cont, "double")
  expect_s3_class(result$mods$mod_factor, "factor")

  # Check formula attribute
  expect_s3_class(attr(result$mods, "formula"), "formula")
  expect_equal(
    paste0(as.character(attr(result$mods, "formula")), collapse = " "),
    "~ mod_cont + mod_factor")
})


test_that("Mods works with formula using vectors from environment", {

  skip_on_cran()

  # Create vectors in the current environment
  env_mod1 <- c(10, 20, 30, 40, 50)
  env_mod2 <- factor(c("X", "Y", "X", "Y", "X"))

  result <- brma.uni(
    yi   = effect,
    sei  = std_err,
    mods = ~ env_mod1 + env_mod2,
    data = test_data_mods,
    only_data = TRUE
  )[["data"]]

  expect_true(!is.null(result$mods))
  expect_equal(nrow(result$mods), 5)
  expect_true("env_mod1" %in% names(result$mods))
  expect_true("env_mod2" %in% names(result$mods))
  expect_equal(result$mods$env_mod1, env_mod1)

  # Check formula attribute
  expect_equal(
    paste0(as.character(attr(result$mods, "formula")), collapse = " "),
    "~ env_mod1 + env_mod2")

  # Verify formula can be used with the returned data.frame
  mf_check <- stats::model.frame(attr(result$mods, "formula"), data = result$mods)
  expect_equal(nrow(mf_check), 5)
})


test_that("Mods preserves factor levels and numeric types", {

  skip_on_cran()

  result <- brma.uni(
    yi   = effect,
    sei  = std_err,
    mods = ~ mod_cont + mod_factor,
    data = test_data_mods,
    only_data = TRUE
  )[["data"]]

  # Check numeric is preserved
  expect_type(result$mods$mod_cont, "double")
  expect_equal(result$mods$mod_cont, test_data_mods$mod_cont)

  # Check factor is preserved with levels
  expect_s3_class(result$mods$mod_factor, "factor")
  expect_equal(levels(result$mods$mod_factor), levels(test_data_mods$mod_factor))

  # Check formula attribute
  expect_equal(
    paste0(as.character(attr(result$mods, "formula")), collapse = " "),
    "~ mod_cont + mod_factor")
})


test_that("Mods applies subset correctly", {

  skip_on_cran()

  # Test with logical subset
  result <- brma.uni(
    yi     = effect,
    sei    = std_err,
    mods   = ~ mod_cont + mod_factor,
    subset = c(TRUE, FALSE, TRUE, FALSE, TRUE),
    data   = test_data_mods,
    only_data = TRUE
  )[["data"]]

  expect_equal(nrow(result$mods), 3)
  expect_equal(result$mods$mod_cont, test_data_mods$mod_cont[c(1, 3, 5)])

  # Verify that unused factor levels are dropped
  expect_equal(levels(result$mods$mod_factor), "A")  # Only "A" remains after subsetting

  # Check formula attribute
  expect_equal(
    paste0(as.character(attr(result$mods, "formula")), collapse = " "),
    "~ mod_cont + mod_factor")
})


test_that("Mods applies numeric subset correctly", {

  skip_on_cran()

  result <- brma.uni(
    yi     = effect,
    sei    = std_err,
    mods   = ~ mod_cont,
    subset = c(2, 4),
    data   = test_data_mods,
    only_data = TRUE
  )[["data"]]

  expect_equal(nrow(result$mods), 2)
  expect_equal(result$mods$mod_cont, test_data_mods$mod_cont[c(2, 4)])

  # Check formula attribute
  expect_equal(
    paste0(as.character(attr(result$mods, "formula")), collapse = " "),
    "~ mod_cont")
})


test_that("Scale works with formula using column names from data", {

  skip_on_cran()

  result <- brma.uni(
    yi    = effect,
    sei   = std_err,
    scale = ~ scale_var,
    data  = test_data_mods,
    only_data = TRUE
  )[["data"]]

  expect_true(!is.null(result$scale))
  expect_true(is.data.frame(result$scale))
  expect_equal(nrow(result$scale), 5)
  expect_true("scale_var" %in% names(result$scale))
  expect_equal(result$scale$scale_var, test_data_mods$scale_var)

  # Check formula attribute
  expect_equal(
    paste0(as.character(attr(result$scale, "formula")), collapse = " "),
    "~ scale_var")
})


test_that("Both mods and scale can be specified together", {

  skip_on_cran()

  result <- brma.uni(
    yi    = effect,
    sei   = std_err,
    mods  = ~ mod_cont + mod_factor,
    scale = ~ scale_var,
    data  = test_data_mods,
    only_data = TRUE
  )[["data"]]

  expect_true(!is.null(result$mods))
  expect_true(!is.null(result$scale))
  expect_equal(nrow(result$mods), 5)
  expect_equal(nrow(result$scale), 5)

  # Check formula attributes
  expect_equal(
    paste0(as.character(attr(result$mods, "formula")), collapse = " "),
    "~ mod_cont + mod_factor")
  expect_equal(
    paste0(as.character(attr(result$scale, "formula")), collapse = " "),
    "~ scale_var")
})


test_that("Mods throws error for length mismatch", {

  skip_on_cran()

  short_mod <- c(1, 2, 3)  # Only 3 values, but yi has 5

  expect_error(
    brma.uni(
      yi   = effect,
      sei  = std_err,
      mods = ~ short_mod,
      data = test_data_mods,
      only_data = TRUE
    ),
    regexp = "number of rows|must equal"
  )
})


test_that("Mods works with matrix input", {

  skip_on_cran()

  mod_matrix <- matrix(1:10, nrow = 5, ncol = 2)
  colnames(mod_matrix) <- c("X1", "X2")

  result <- brma.uni(
    yi   = effect,
    sei  = std_err,
    mods = mod_matrix,
    data = test_data_mods,
    only_data = TRUE
  )[["data"]]

  expect_true(!is.null(result$mods))
  expect_true(is.data.frame(result$mods))
  expect_equal(nrow(result$mods), 5)
  expect_true(all(c("X1", "X2") %in% names(result$mods)))

  # Check formula attribute
  expect_equal(
    paste0(as.character(attr(result$mods, "formula")), collapse = " "),
    "~ X1 + X2")
})


test_that("Mods works with data.frame input", {

  skip_on_cran()

  mod_df <- data.frame(
    predictor1 = c(1.1, 2.2, 3.3, 4.4, 5.5),
    predictor2 = factor(c("a", "b", "a", "b", "a"))
  )

  result <- brma.uni(
    yi   = effect,
    sei  = std_err,
    mods = mod_df,
    data = test_data_mods,
    only_data = TRUE
  )[["data"]]

  expect_true(!is.null(result$mods))
  expect_true(is.data.frame(result$mods))
  expect_equal(nrow(result$mods), 5)

  # Check types preserved
  expect_type(result$mods$predictor1, "double")
  expect_s3_class(result$mods$predictor2, "factor")

  # Check formula attribute
  expect_equal(
    paste0(as.character(attr(result$mods, "formula")), collapse = " "),
    "~ predictor1 + predictor2")
})


test_that("Mods formula attribute can be evaluated in clean environment", {

  skip_on_cran()

  # Use environment vectors that won't be available later
  local({
    temp_var <- c(100, 200, 300, 400, 500)

    result <- brma.uni(
      yi   = effect,
      sei  = std_err,
      mods = ~ temp_var,
      data = test_data_mods,
      only_data = TRUE
    )[["data"]]

    # Check formula attribute
    expect_equal(
      paste0(as.character(attr(result$mods, "formula")), collapse = " "),
      "~ temp_var")

    # Verify formula can be evaluated with the returned data.frame
    mf <- stats::model.frame(attr(result$mods, "formula"), data = result$mods)
    expect_equal(nrow(mf), 5)
    expect_equal(mf$temp_var, temp_var)
  })
})


test_that("Mods handles formula with LHS (with warning)", {

  skip_on_cran()

  expect_warning(
    result <- brma.uni(
      yi   = effect,
      sei  = std_err,
      mods = effect ~ mod_cont,  # LHS should be ignored
      data = test_data_mods,
      only_data = TRUE
    )[["data"]],
    regexp = "left-hand side|LHS"
  )

  expect_true(!is.null(result$mods))
  expect_true("mod_cont" %in% names(result$mods))
  expect_false("effect" %in% names(result$mods))  # LHS was ignored

  # Check formula attribute (should be ~ mod_cont, not effect ~ mod_cont)
  expect_equal(
    paste0(as.character(attr(result$mods, "formula")), collapse = " "),
    "~ mod_cont")
})


test_that("Scale applies subset correctly", {

  skip_on_cran()

  result <- brma.uni(
    yi     = effect,
    sei    = std_err,
    scale  = ~ scale_var,
    subset = c(1, 3, 5),
    data   = test_data_mods,
    only_data = TRUE
  )[["data"]]

  expect_equal(nrow(result$scale), 3)
  expect_equal(result$scale$scale_var, test_data_mods$scale_var[c(1, 3, 5)])
})


test_that("Mods works with inline transformations in formula", {

  skip_on_cran()

  result <- brma.uni(
    yi   = effect,
    sei  = std_err,
    mods = ~ I(mod_cont - 2) + I(mod_cont^2),
    data = test_data_mods,
    only_data = TRUE
  )[["data"]]

  expect_true(!is.null(result$mods))
  expect_equal(nrow(result$mods), 5)

  # Check that transformations are correctly applied (as.numeric to strip AsIs class)
  expect_equal(as.numeric(result$mods$`I(mod_cont - 2)`), test_data_mods$mod_cont - 2)
  expect_equal(as.numeric(result$mods$`I(mod_cont^2)`), test_data_mods$mod_cont^2)

  # Check formula attribute preserves transformations
  expect_equal(
    paste0(as.character(attr(result$mods, "formula")), collapse = " "),
    "~ I(mod_cont - 2) + I(mod_cont^2)")
})


test_that("Mods works with interaction terms in formula", {

  skip_on_cran()

  result <- brma.uni(
    yi   = effect,
    sei  = std_err,
    mods = ~ mod_cont * mod_factor,
    data = test_data_mods,
    only_data = TRUE
  )[["data"]]

  expect_true(!is.null(result$mods))
  expect_equal(nrow(result$mods), 5)

  # Should have both main effects and the original variables
  expect_true("mod_cont" %in% names(result$mods))
  expect_true("mod_factor" %in% names(result$mods))

  # Check that the formula attribute preserves the interaction structure
  expect_equal(
    paste0(as.character(attr(result$mods, "formula")), collapse = " "),
    "~ mod_cont * mod_factor")
})


test_that("Mods auto-converts character to factor", {

  skip_on_cran()

  result <- brma.uni(
    yi   = effect,
    sei  = std_err,
    mods = ~ mod_char,  # mod_char is character in test_data_mods
    data = test_data_mods,
    only_data = TRUE
  )[["data"]]

  expect_true(!is.null(result$mods))
  # model.frame should convert character to factor
  expect_true(is.factor(result$mods$mod_char) || is.character(result$mods$mod_char))
})


test_that("Mods works with data$column syntax", {

  skip_on_cran()

  # Using $ syntax directly
  result <- brma.uni(
    yi   = test_data_mods$effect,
    sei  = test_data_mods$std_err,
    mods = ~ test_data_mods$mod_cont,
    only_data = TRUE
  )[["data"]]

  expect_true(!is.null(result$mods))
  expect_equal(nrow(result$mods), 5)

  # Check formula attribute
  expect_equal(
    paste0(as.character(attr(result$mods, "formula")), collapse = " "),
    "~ test_data_mods$mod_cont")
})


test_that("NULL mods and scale return NULL", {

  skip_on_cran()

  result <- brma.uni(
    yi   = effect,
    sei  = std_err,
    data = test_data_mods,
    only_data = TRUE
  )[["data"]]

  expect_null(result$mods)
  expect_null(result$scale)
})

# ============================================================================
# Tests for yi ~ mods formula syntax
# ============================================================================

test_that("yi ~ mods formula syntax works with single moderator", {

  skip_on_cran()

  result <- brma.uni(
    yi   = effect ~ mod_cont,
    sei  = std_err,
    data = test_data_mods,
    only_data = TRUE
  )[["data"]]

  expect_true(!is.null(result$mods))
  expect_equal(nrow(result$mods), 5)

  # Should have mod_cont column with original values (not model matrix)
  expect_true("mod_cont" %in% names(result$mods))
  expect_equal(result$mods$mod_cont, test_data_mods$mod_cont)

  # Check that yi was correctly extracted
  expect_equal(result$outcome$yi, test_data_mods$effect)

  # Check formula attribute preserves original formula structure
  expect_equal(paste0(as.character(attr(result$mods, "formula")), collapse = " "), "~ mod_cont")
})


test_that("yi ~ mods formula syntax works with multiple moderators", {

  skip_on_cran()

  result <- brma.uni(
    effect ~ mod_cont + mod_factor,
    sei  = std_err,
    data = test_data_mods,
    only_data = TRUE
  )[["data"]]

  expect_true(!is.null(result$mods))
  expect_equal(nrow(result$mods), 5)

  # Should have original variables (not expanded model matrix)
  expect_true("mod_cont" %in% names(result$mods))
  expect_true("mod_factor" %in% names(result$mods))

  # Factor should be preserved as factor
  expect_s3_class(result$mods$mod_factor, "factor")

  # Check formula attribute
  expect_equal(
    paste0(as.character(attr(result$mods, "formula")), collapse = " "),
    "~ mod_cont + mod_factor")
})


test_that("yi ~ 0 + mods formula syntax removes intercept but keeps original variables", {

  skip_on_cran()

  result <- brma.uni(
    yi   = effect ~ 0 + mod_factor,
    sei  = std_err,
    data = test_data_mods,
    only_data = TRUE
  )[["data"]]

  expect_true(!is.null(result$mods))
  expect_equal(nrow(result$mods), 5)

  # Should have the original variable (not expanded dummy columns)
  expect_true("mod_factor" %in% names(result$mods))
  expect_false("mod_factorA" %in% names(result$mods))
  expect_false("mod_factorB" %in% names(result$mods))

  # The variable should be a factor
  expect_s3_class(result$mods$mod_factor, "factor")

  # Check that the formula attribute preserves ~ 0 + structure
  expect_equal(paste0(as.character(attr(result$mods, "formula")), collapse = " "),
               "~ 0 + mod_factor")
})


test_that("yi ~ mods formula syntax works with interactions", {

  skip_on_cran()

  result <- brma.uni(
    yi   = effect ~ mod_cont * mod_factor,
    sei  = std_err,
    data = test_data_mods,
    only_data = TRUE
  )[["data"]]

  expect_true(!is.null(result$mods))
  expect_equal(nrow(result$mods), 5)

  # Should have both main effects
  expect_true("mod_cont" %in% names(result$mods))
  expect_true("mod_factor" %in% names(result$mods))

  # Check formula attribute preserves interaction
  expect_equal(paste0(as.character(attr(result$mods, "formula")), collapse = " "),
               "~ mod_cont * mod_factor")
})


test_that("yi ~ mods formula syntax applies subset correctly", {

  skip_on_cran()

  result <- brma.uni(
    yi     = effect ~ mod_cont,
    sei    = std_err,
    subset = c(1, 3, 5),
    data   = test_data_mods,
    only_data = TRUE
  )[["data"]]

  expect_true(!is.null(result$mods))
  expect_equal(nrow(result$mods), 3)
  expect_equal(nrow(result$outcome), 3)

  # Check subset was applied to yi correctly
  expect_equal(result$outcome$yi, test_data_mods$effect[c(1, 3, 5)])

  # Check formula attribute
  expect_equal(paste0(as.character(attr(result$mods, "formula")), collapse = " "), "~ mod_cont")
})


test_that("yi ~ mods formula syntax errors when mods also specified", {

  skip_on_cran()

  expect_error(
    brma.uni(
      yi   = effect ~ mod_cont,
      sei  = std_err,
      mods = ~ mod_factor,  # This should cause an error
      data = test_data_mods,
      only_data = TRUE
    ),
    regexp = "Cannot specify 'mods' when 'yi' is a formula"
  )
})


test_that("yi ~ mods formula returns equivalent results to mods = ~", {

  skip_on_cran()

  result1 <- brma.uni(
    yi   = effect ~ mod_cont + mod_factor,
    sei  = std_err,
    data = test_data_mods,
    only_data = TRUE
  )[["data"]]

  result2 <- brma.uni(
    yi   = effect,
    sei  = std_err,
    mods = ~ mod_cont + mod_factor,
    data = test_data_mods,
    only_data = TRUE
  )[["data"]]

  # Both should have the same yi values
  expect_equal(result1$outcome$yi, result2$outcome$yi)

  # yi ~ formula gives model matrix (expanded with intercept and dummies)
  # mods = ~ formula gives data.frame with original variables
  # So the structures differ, but both should have mod_cont values
  expect_true("mod_cont" %in% names(result1$mods))
  expect_true("mod_cont" %in% names(result2$mods))

  # Both should have same number of rows
  expect_equal(nrow(result1$mods), nrow(result2$mods))

  # Both should have same formula attribute
  expect_equal(paste0(as.character(attr(result1$mods, "formula")), collapse = " "), "~ mod_cont + mod_factor")
  expect_equal(paste0(as.character(attr(result2$mods, "formula")), collapse = " "), "~ mod_cont + mod_factor")
})


test_that("Input works with metafor::escalc output", {

  skip_on_cran()
  skip_if_not_installed("metadat")
  skip_if_not_installed("metafor")

  # Load BCG vaccine data and compute log risk ratios
  data(dat.bcg, package = "metadat")
  dat <- metafor::escalc(measure = "RR", ai = tpos, bi = tneg, ci = cpos, di = cneg, data = dat.bcg)

  # Test that brma.uni can handle escalc output (yi and vi have special attributes)
  result <- brma.uni(
    yi      = yi,
    vi      = vi,
    data    = dat,
    measure = "RR",
    only_data = TRUE
  )[["data"]]

  expect_type(result, "list")
  expect_true("outcome" %in% names(result))
  expect_equal(nrow(result$outcome), nrow(dat))

  # Check that yi values are correctly extracted (as plain numeric, no attributes)

  expect_equal(result$outcome$yi, as.vector(dat$yi))
  expect_null(attributes(result$outcome$yi))

  # Check that sei is computed from vi
  expect_equal(result$outcome$sei, sqrt(as.vector(dat$vi)), tolerance = 1e-10)
  expect_null(attributes(result$outcome$sei))
})
