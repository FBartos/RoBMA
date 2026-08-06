context("Input handling for brma.norm")

source(testthat::test_path("helper-contracts.R"))

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

expect_norm_outcome <- function(result, n, yi, sei = NULL, ni = NULL,
                                weights = NULL, slab = NULL, cluster = NULL) {

  expect_type(result, "list")
  expect_true("outcome" %in% names(result))
  expect_equal(nrow(result$outcome), n)
  expect_equal(result$outcome$yi, yi)

  if (!is.null(sei)) {
    expect_equal(result$outcome$sei, sei, tolerance = 1e-10)
  }
  if (!is.null(ni)) {
    expect_equal(result$outcome$ni, ni)
  }
  if (!is.null(weights)) {
    expect_equal(result$outcome$weights, weights)
  }
  if (!is.null(slab)) {
    expect_equal(result$outcome$slab, slab)
  }
  if (!is.null(cluster)) {
    expect_equal(result$outcome$cluster, cluster)
  }
}

formula_text <- function(x) {

  return(paste0(as.character(attr(x, "formula")), collapse = " "))
}


test_that("Input accepts supported outcome specifications", {

  skip_on_cran()

  external_weights <- c(2, 2, 2, 2, 2)

  input_cases <- list(
    list(
      expr = quote(brma.norm(
        yi        = c(0.10, 0.25, 0.15),
        sei       = c(0.20, 0.24, 0.22),
        only_data = TRUE
      )[["data"]]),
      n    = 3,
      yi   = c(0.10, 0.25, 0.15),
      sei  = c(0.20, 0.24, 0.22)
    ),
    list(
      expr = quote(brma.norm(
        yi        = effect,
        vi        = variance,
        data      = test_data,
        only_data = TRUE
      )[["data"]]),
      n    = 5,
      yi   = test_data$effect,
      sei  = sqrt(test_data$variance)
    ),
    list(
      expr = quote(brma.norm(
        yi        = "effect",
        sei       = "std_err",
        data      = test_data,
        only_data = TRUE
      )[["data"]]),
      n    = 5,
      yi   = test_data$effect,
      sei  = test_data$std_err
    ),
    list(
      expr = quote(brma.norm(
        yi        = effect,
        vi        = "variance",
        ni        = n,
        weights   = external_weights,
        data      = test_data,
        only_data = TRUE
      )[["data"]]),
      n       = 5,
      yi      = test_data$effect,
      sei     = sqrt(test_data$variance),
      ni      = test_data$n,
      weights = external_weights
    ),
    list(
      expr = quote(brma.norm(
        yi        = effect,
        sei       = std_err,
        ni        = n,
        slab      = study,
        cluster   = cluster,
        data      = test_data,
        only_data = TRUE
      )[["data"]]),
      n       = 5,
      yi      = test_data$effect,
      sei     = test_data$std_err,
      ni      = test_data$n,
      slab    = test_data$study,
      cluster = as.numeric(as.factor(test_data$cluster))
    )
  )

  for (case in input_cases) {
    expect_norm_outcome(
      eval(case$expr),
      n       = case$n,
      yi      = case$yi,
      sei     = case$sei,
      ni      = case$ni,
      weights = case$weights,
      slab    = case$slab,
      cluster = case$cluster
    )
  }
})


test_that("Input generates default slab when not provided", {

  skip_on_cran()

  result <- brma.norm(
    yi  = c(0.1, 0.2, 0.3),
    sei = c(0.1, 0.1, 0.1),
    only_data = TRUE
  )[["data"]]

  expect_equal(result$outcome$slab, c("Study 1", "Study 2", "Study 3"))
})


test_that("normal input requires representable and consistent sampling variances", {

  for (sei in c(1e-200, 1e200)) {
    expect_error(
      brma.norm(yi = 0, sei = sei, only_data = TRUE),
      "positive finite squared sampling variances"
    )
  }

  expect_silent(brma.norm(
    yi        = c(-1, 0, 1),
    vi        = c(1e-300, .04, 1e300),
    sei       = sqrt(c(1e-300, .04, 1e300)),
    only_data = TRUE
  ))
  expect_error(
    brma.norm(
      yi        = 0,
      vi        = 1e-300,
      sei       = 1e-150 * (1 + 1e-10),
      only_data = TRUE
    ),
    "inconsistent"
  )
})


test_that("Input handles subset argument", {

  skip_on_cran()

  # Logical subset
  result_logical <- brma.norm(
    yi     = effect,
    sei    = std_err,
    subset = c(TRUE, FALSE, TRUE, FALSE, TRUE),
    data   = test_data,
    only_data = TRUE
  )[["data"]]

  expect_equal(nrow(result_logical$outcome), 3)
  expect_equal(result_logical$outcome$yi, test_data$effect[c(1, 3, 5)])

  # Numeric subset
  result_numeric <- brma.norm(
    yi     = effect,
    sei    = std_err,
    subset = c(1, 3, 5),
    data   = test_data,
    only_data = TRUE
  )[["data"]]

  expect_equal(nrow(result_numeric$outcome), 3)
  expect_equal(result_numeric$outcome$yi, test_data$effect[c(1, 3, 5)])
})


test_that("Effect direction is passed and detected after final row selection", {

  skip_on_cran()

  result_explicit <- RoBMA(
    yi               = c(-0.10, -0.20, -0.30),
    sei              = c(0.10, 0.10, 0.10),
    effect_direction = "negative",
    only_data        = TRUE
  )[["data"]]

  expect_equal(attr(result_explicit, "effect_direction"), "negative")

  result_robma_detect <- RoBMA(
    yi               = c(-0.10, -0.20, -0.30),
    sei              = c(0.10, 0.10, 0.10),
    effect_direction = "detect",
    only_data        = TRUE
  )[["data"]]

  expect_equal(attr(result_robma_detect, "effect_direction"), "negative")

  result_subset <- bselmodel(
    yi               = c(100, 90, -1, -2),
    sei              = c(0.10, 0.10, 0.10, 0.10),
    subset           = 3:4,
    effect_direction = "detect",
    only_data        = TRUE
  )[["data"]]

  expect_equal(attr(result_subset, "effect_direction"), "negative")

  data_na_direction <- data.frame(
    yi  = c(100, 90, 80, -1, -2),
    sei = c(0.10, 0.10, 0.10, 0.10, 0.10),
    mod = c(NA, NA, NA, 1, 2)
  )

  expect_warning(
    result_na <- bselmodel(
      yi               = yi,
      sei              = sei,
      mods             = ~ mod,
      data             = data_na_direction,
      effect_direction = "detect",
      only_data        = TRUE
    )[["data"]],
    regexp = "3 observation.*removed due to missing"
  )

  expect_equal(attr(result_na, "effect_direction"), "negative")
})


test_that("Input rejects invalid outcome specifications", {

  skip_on_cran()

  expect_error_cases(list(
    list(
      label  = "missing yi",
      expr   = quote(brma.norm(sei = c(0.1, 0.2), only_data = TRUE)),
      regexp = "yi"
    ),
    list(
      label  = "missing variance",
      expr   = quote(brma.norm(yi = c(0.1, 0.2), only_data = TRUE)),
      regexp = "vi|sei|variance|standard error"
    ),
    list(
      label  = "length mismatch",
      expr   = quote(brma.norm(
        yi        = c(0.1, 0.2, 0.3),
        sei       = c(0.1, 0.2),
        only_data = TRUE
      )),
      regexp = "length"
    ),
    list(
      label  = "invalid column reference",
      expr   = quote(brma.norm(
        yi        = nonexistent_column,
        sei       = std_err,
        data      = test_data,
        only_data = TRUE
      )),
      regexp = "nonexistent_column|Cannot find"
    )
  ))
})


test_that("Input rejects unresolved supplied optional variables", {

  skip_on_cran()

  expect_error_cases(list(
    list(
      label  = "unresolved weights",
      expr   = quote(brma.norm(
        yi        = effect,
        sei       = std_err,
        weights   = wights,
        data      = test_data,
        only_data = TRUE
      )),
      regexp = "Cannot find.*weights"
    ),
    list(
      label  = "unresolved cluster",
      expr   = quote(brma.norm(
        yi        = effect,
        sei       = std_err,
        cluster   = cluser,
        data      = test_data,
        only_data = TRUE
      )),
      regexp = "Cannot find.*cluster"
    ),
    list(
      label  = "unresolved slab",
      expr   = quote(brma.norm(
        yi        = effect,
        sei       = std_err,
        slab      = lables,
        data      = test_data,
        only_data = TRUE
      )),
      regexp = "Cannot find.*slab"
    ),
    list(
      label  = "unresolved subset",
      expr   = quote(brma.norm(
        yi        = effect,
        sei       = std_err,
        subset    = subest,
        data      = test_data,
        only_data = TRUE
      )),
      regexp = "Cannot find.*subset"
    ),
    list(
      label  = "unresolved sample size",
      expr   = quote(brma.norm(
        yi        = effect,
        sei       = std_err,
        ni        = sample_size,
        data      = test_data,
        only_data = TRUE
      )),
      regexp = "Cannot find.*ni"
    )
  ))

  result <- brma.norm(
    yi        = effect,
    sei       = std_err,
    ni        = NULL,
    weights   = NULL,
    cluster   = NULL,
    slab      = NULL,
    subset    = NULL,
    data      = test_data,
    only_data = TRUE
  )[["data"]]

  expect_equal(nrow(result$outcome), nrow(test_data))
  expect_false(attr(result, "weights"))
  expect_false(attr(result, "cluster"))
  expect_false(attr(result, "slab"))
})


test_that("Input handles variables from calling environment", {

  skip_on_cran()

  # Create variables in the current environment
  my_effects   <- c(0.15, 0.25, 0.35)
  my_std_errs  <- c(0.10, 0.12, 0.11)
  my_labels    <- c("Study X", "Study Y", "Study Z")

  result <- brma.norm(
    yi   = my_effects,
    sei  = my_std_errs,
    slab = my_labels,
    only_data = TRUE
  )[["data"]]

  expect_equal(result$outcome$yi, my_effects)
  expect_equal(result$outcome$sei, my_std_errs)
  expect_equal(result$outcome$slab, my_labels)
})


test_that("Input resolves direct column references inside functions", {

  skip_on_cran()

  # Wrapper function that uses NSE to pass column names
  # Note: When brma.norm is called from another function, the column names

  # must be available in the data frame or the calling environment of brma.norm
  wrapper_function <- function(data) {
    # Direct column references work because they are evaluated in data
    brma.norm(
      yi   = effect,
      sei  = std_err,
      data = data,
      only_data = TRUE
    )[["data"]]
  }

  # This tests that NSE works correctly when brma.norm is called from another function
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


test_that("Moderator formulas preserve data columns, types, and formula metadata", {

  skip_on_cran()

  result <- brma.norm(
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
  expect_true(all(c("mod_cont", "mod_factor") %in% names(result$mods)))
  expect_equal(result$mods$mod_cont, test_data_mods$mod_cont)
  expect_type(result$mods$mod_cont, "double")
  expect_s3_class(result$mods$mod_factor, "factor")
  expect_equal(levels(result$mods$mod_factor), levels(test_data_mods$mod_factor))
  expect_s3_class(attr(result$mods, "formula"), "formula")
  expect_equal(formula_text(result$mods), "~ mod_cont + mod_factor")
})


test_that("Mods accepts formula using vectors from environment", {

  skip_on_cran()

  env_mod1 <- c(10, 20, 30, 40, 50)
  env_mod2 <- factor(c("X", "Y", "X", "Y", "X"))

  result <- brma.norm(
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
  expect_equal(formula_text(result$mods), "~ env_mod1 + env_mod2")

  mf_check <- stats::model.frame(attr(result$mods, "formula"), data = result$mods)
  expect_equal(nrow(mf_check), 5)
})


test_that("Mods and scale apply logical and numeric subsets", {

  skip_on_cran()

  result_logical <- brma.norm(
    yi     = effect,
    sei    = std_err,
    mods   = ~ mod_cont + mod_factor,
    subset = c(TRUE, TRUE, FALSE, TRUE, FALSE),
    data   = test_data_mods,
    only_data = TRUE
  )[["data"]]

  result_numeric <- brma.norm(
    yi     = effect,
    sei    = std_err,
    mods   = ~ mod_cont,
    subset = c(2, 4),
    data   = test_data_mods,
    only_data = TRUE
  )[["data"]]

  result_scale <- brma.norm(
    yi        = effect,
    sei       = std_err,
    scale     = ~ scale_var,
    subset    = c(1, 3, 5),
    data      = test_data_mods,
    only_data = TRUE
  )[["data"]]

  expect_equal(nrow(result_logical$mods), 3)
  expect_equal(result_logical$mods$mod_cont, test_data_mods$mod_cont[c(1, 2, 4)])
  expect_equal(levels(result_logical$mods$mod_factor), c("A", "B"))
  expect_equal(formula_text(result_logical$mods), "~ mod_cont + mod_factor")

  expect_equal(nrow(result_numeric$mods), 2)
  expect_equal(result_numeric$mods$mod_cont, test_data_mods$mod_cont[c(2, 4)])
  expect_equal(formula_text(result_numeric$mods), "~ mod_cont")

  expect_equal(nrow(result_scale$scale), 3)
  expect_equal(result_scale$scale$scale_var, test_data_mods$scale_var[c(1, 3, 5)])
})


test_that("Scale formulas are accepted alone and with moderators", {

  skip_on_cran()

  scale_result <- brma.norm(
    yi        = effect,
    sei       = std_err,
    scale     = ~ scale_var,
    data      = test_data_mods,
    only_data = TRUE
  )[["data"]]

  joint_result <- brma.norm(
    yi    = effect,
    sei   = std_err,
    mods  = ~ mod_cont + mod_factor,
    scale = ~ scale_var,
    data  = test_data_mods,
    only_data = TRUE
  )[["data"]]

  expect_true(is.data.frame(scale_result$scale))
  expect_equal(nrow(scale_result$scale), 5)
  expect_true("scale_var" %in% names(scale_result$scale))
  expect_equal(scale_result$scale$scale_var, test_data_mods$scale_var)
  expect_equal(formula_text(scale_result$scale), "~ scale_var")

  expect_true(!is.null(joint_result$mods))
  expect_true(!is.null(joint_result$scale))
  expect_equal(nrow(joint_result$mods), 5)
  expect_equal(nrow(joint_result$scale), 5)
  expect_equal(formula_text(joint_result$mods), "~ mod_cont + mod_factor")
  expect_equal(formula_text(joint_result$scale), "~ scale_var")
})


test_that("Mods rejects length mismatch", {

  skip_on_cran()

  short_mod <- c(1, 2, 3)  # Only 3 values, but yi has 5

  expect_error(
    brma.norm(
      yi   = effect,
      sei  = std_err,
      mods = ~ short_mod,
      data = test_data_mods,
      only_data = TRUE
    ),
    regexp = "number of rows|must equal"
  )
})


test_that("Mods accepts matrix input", {

  skip_on_cran()

  mod_matrix <- matrix(1:10, nrow = 5, ncol = 2)
  colnames(mod_matrix) <- c("X1", "X2")

  result <- brma.norm(
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


test_that("Mods accepts data.frame input", {

  skip_on_cran()

  mod_df <- data.frame(
    predictor1 = c(1.1, 2.2, 3.3, 4.4, 5.5),
    predictor2 = factor(c("a", "b", "a", "b", "a"))
  )

  result <- brma.norm(
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

    result <- brma.norm(
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


test_that("Persisted formula metadata does not retain caller frames", {

  make_object <- function(payload_size) {

    local({
      unrelated_payload <- raw(payload_size)
      temp_var <- c(100, 200, 300, 400, 500)

      brma.norm(
        yi        = effect,
        sei       = std_err,
        mods      = ~ temp_var,
        data      = test_data_mods,
        only_data = TRUE
      )
    })
  }

  baseline <- make_object(0L)
  payload  <- make_object(5e6L)

  expect_lt(
    abs(length(serialize(payload, NULL)) - length(serialize(baseline, NULL))),
    10000L
  )
  expect_identical(
    environment(attr(payload[["data"]][["mods"]], "formula")),
    baseenv()
  )
  expect_identical(
    attr(attr(payload[["data"]][["mods"]], "terms"), ".Environment"),
    baseenv()
  )

  model_frame <- stats::model.frame(
    attr(payload[["data"]][["mods"]], "formula"),
    data = payload[["data"]][["mods"]]
  )
  expect_equal(model_frame[["temp_var"]], c(100, 200, 300, 400, 500))
})


test_that("Mods handles formula with LHS (with warning)", {

  skip_on_cran()

  expect_warning(
    result <- brma.norm(
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


test_that("Mods accepts inline transformations in formula", {

  skip_on_cran()

  result <- brma.norm(
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


test_that("Mods accepts interaction terms in formula", {

  skip_on_cran()

  result <- brma.norm(
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

  result <- brma.norm(
    yi   = effect,
    sei  = std_err,
    mods = ~ mod_char,  # mod_char is character in test_data_mods
    data = test_data_mods,
    only_data = TRUE
  )[["data"]]

  expect_true(!is.null(result$mods))
  expect_true(is.factor(result$mods$mod_char))
  expect_equal(as.character(result$mods$mod_char), test_data_mods$mod_char)
})


test_that("yi formula auto-converts character moderators to factor", {

  skip_on_cran()

  result <- brma.norm(
    yi   = effect ~ mod_char,
    sei  = std_err,
    data = test_data_mods,
    only_data = TRUE
  )[["data"]]

  expect_true(!is.null(result$mods))
  expect_true(is.factor(result$mods$mod_char))
  expect_equal(as.character(result$mods$mod_char), test_data_mods$mod_char)
  expect_equal(paste0(as.character(attr(result$mods, "formula")), collapse = " "), "~ mod_char")
})


test_that("Mods accepts data$column syntax", {

  skip_on_cran()

  # Using $ syntax directly
  result <- brma.norm(
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

  result <- brma.norm(
    yi   = effect,
    sei  = std_err,
    data = test_data_mods,
    only_data = TRUE
  )[["data"]]

  expect_null(result$mods)
  expect_null(result$scale)
})


# ============================================================================
# Tests for predictor validation (constant variables and single-level factors)
# ============================================================================

test_that("Degenerate moderator and scale designs are rejected", {

  skip_on_cran()

  test_data_const <- data.frame(
    effect    = c(0.10, 0.25, 0.15),
    std_err   = c(0.2, 0.25, 0.22),
    const_var = c(5, 5, 5)
  )

  test_data_single <- data.frame(
    effect     = c(0.10, 0.25, 0.15),
    std_err    = c(0.2, 0.25, 0.22),
    single_fac = factor(c("A", "A", "A"))
  )

  test_data_const_scale <- data.frame(
    effect     = c(0.10, 0.25, 0.15),
    std_err    = c(0.2, 0.25, 0.22),
    mod_cont   = c(1.5, 2.3, 1.8),
    const_scale = c(1, 1, 1)
  )

  expect_error_cases(list(
    list(
      label  = "constant numeric moderator",
      expr   = quote(brma.norm(
        yi        = effect,
        sei       = std_err,
        mods      = ~ const_var,
        data      = test_data_const,
        only_data = TRUE
      )),
      regexp = "zero variance"
    ),
    list(
      label  = "single-level factor moderator",
      expr   = quote(brma.norm(
        yi        = effect,
        sei       = std_err,
        mods      = ~ single_fac,
        data      = test_data_single,
        only_data = TRUE
      )),
      regexp = "only one level"
    ),
    list(
      label  = "subset drops all but one factor level",
      expr   = quote(brma.norm(
        yi        = effect,
        sei       = std_err,
        mods      = ~ mod_factor,
        subset    = c(TRUE, FALSE, TRUE, FALSE, TRUE),
        data      = test_data_mods,
        only_data = TRUE
      )),
      regexp = "only one level"
    ),
    list(
      label  = "constant scale predictor",
      expr   = quote(brma.norm(
        yi        = effect,
        sei       = std_err,
        mods      = ~ mod_cont,
        scale     = ~ const_scale,
        data      = test_data_const_scale,
        only_data = TRUE
      )),
      regexp = "zero variance"
    )
  ))
})


test_that("Valid mods with variation passes validation", {

  skip_on_cran()

  test_data_valid <- data.frame(
    effect   = c(0.10, 0.25, 0.15, 0.20),
    std_err  = c(0.2, 0.25, 0.22, 0.18),
    mod_cont = c(1.5, 2.3, 1.8, 3.0),
    mod_fac  = factor(c("A", "B", "A", "B"))
  )

  result <- brma.norm(
    yi   = effect,
    sei  = std_err,
    mods = ~ mod_cont + mod_fac,
    data = test_data_valid,
    only_data = TRUE
  )[["data"]]

  expect_true(!is.null(result$mods))
  expect_equal(nrow(result$mods), 4)
})


# ============================================================================
# Tests for yi ~ mods formula syntax
# ============================================================================

test_that("yi ~ mods formula syntax preserves variables and formula metadata", {

  skip_on_cran()

  formula_cases <- list(
    list(
      expr      = quote(brma.norm(
        yi        = effect ~ mod_cont,
        sei       = std_err,
        data      = test_data_mods,
        only_data = TRUE
      )[["data"]]),
      rows      = seq_len(5),
      variables = "mod_cont",
      formula   = "~ mod_cont"
    ),
    list(
      expr      = quote(brma.norm(
        effect ~ mod_cont + mod_factor,
        sei       = std_err,
        data      = test_data_mods,
        only_data = TRUE
      )[["data"]]),
      rows      = seq_len(5),
      variables = c("mod_cont", "mod_factor"),
      formula   = "~ mod_cont + mod_factor"
    ),
    list(
      expr      = quote(brma.norm(
        yi        = effect ~ 0 + mod_factor,
        sei       = std_err,
        data      = test_data_mods,
        only_data = TRUE
      )[["data"]]),
      rows      = seq_len(5),
      variables = "mod_factor",
      absent    = c("mod_factorA", "mod_factorB"),
      formula   = "~ 0 + mod_factor"
    ),
    list(
      expr      = quote(brma.norm(
        yi        = effect ~ mod_cont * mod_factor,
        sei       = std_err,
        data      = test_data_mods,
        only_data = TRUE
      )[["data"]]),
      rows      = seq_len(5),
      variables = c("mod_cont", "mod_factor"),
      formula   = "~ mod_cont * mod_factor"
    ),
    list(
      expr      = quote(brma.norm(
        yi        = effect ~ mod_cont,
        sei       = std_err,
        subset    = c(1, 3, 5),
        data      = test_data_mods,
        only_data = TRUE
      )[["data"]]),
      rows      = c(1, 3, 5),
      variables = "mod_cont",
      formula   = "~ mod_cont"
    )
  )

  for (case in formula_cases) {
    result <- eval(case$expr)

    expect_true(!is.null(result$mods))
    expect_equal(nrow(result$mods), length(case$rows))
    expect_equal(nrow(result$outcome), length(case$rows))
    expect_equal(result$outcome$yi, test_data_mods$effect[case$rows])
    expect_true(all(case$variables %in% names(result$mods)))
    expect_equal(formula_text(result$mods), case$formula)

    if ("mod_cont" %in% case$variables) {
      expect_equal(result$mods$mod_cont, test_data_mods$mod_cont[case$rows])
    }
    if ("mod_factor" %in% case$variables) {
      expect_s3_class(result$mods$mod_factor, "factor")
    }
    if (!is.null(case$absent)) {
      expect_false(any(case$absent %in% names(result$mods)))
    }
  }
})


test_that("yi ~ mods formula syntax rejects duplicate mods argument", {

  skip_on_cran()

  expect_error(
    brma.norm(
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

  result1 <- brma.norm(
    yi   = effect ~ mod_cont + mod_factor,
    sei  = std_err,
    data = test_data_mods,
    only_data = TRUE
  )[["data"]]

  result2 <- brma.norm(
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


test_that("Input accepts metafor::escalc output", {

  skip_on_cran()
  skip_if_not_installed("metadat")
  skip_if_not_installed("metafor")

  # Load BCG vaccine data and compute log risk ratios
  data(dat.bcg, package = "metadat")
  dat <- metafor::escalc(measure = "RR", ai = tpos, bi = tneg, ci = cpos, di = cneg, data = dat.bcg)

  # Test that brma.norm can handle escalc output (yi and vi have special attributes)
  result <- brma.norm(
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


# ----- NA Handling Tests -----

# Test data with NAs for NA handling tests
test_data_na <- data.frame(
  effect     = c(0.10, NA, 0.15, 0.30, 0.05, 0.20),
  variance   = c(0.04, 0.06, NA, 0.08, 0.03, 0.05),
  std_err    = sqrt(c(0.04, 0.06, 0.05, 0.08, 0.03, 0.05)),
  n          = c(50L, 75L, 60L, 40L, 90L, 55L),
  study      = c("A", "B", "C", "D", "E", "F"),
  cluster    = c("g1", "g1", "g2", "g2", "g3", "g3"),
  mod_cont   = c(1.5, 2.0, NA, 3.0, 1.0, 2.5),
  mod_factor = factor(c("low", "high", "low", NA, "medium", "high")),
  scale_var  = c(0.5, NA, 0.7, 0.8, 0.6, 0.9),
  stringsAsFactors = FALSE
)


test_that("Outcome NAs drop rows with warning", {

  skip_on_cran()

  na_cases <- list(
    list(
      yi          = c(0.10, NA, 0.15, 0.30, 0.05),
      sei         = c(0.20, 0.24, 0.22, 0.28, 0.17),
      regexp      = "1 observation.*removed due to missing",
      expected_yi = c(0.10, 0.15, 0.30, 0.05),
      expected_se = c(0.20, 0.22, 0.28, 0.17)
    ),
    list(
      yi          = c(0.10, 0.25, 0.15, 0.30, 0.05),
      sei         = c(0.20, 0.24, NA, 0.28, 0.17),
      regexp      = "1 observation.*removed due to missing",
      expected_yi = c(0.10, 0.25, 0.30, 0.05),
      expected_se = c(0.20, 0.24, 0.28, 0.17)
    ),
    list(
      yi          = c(0.10, NA, 0.15, 0.30, 0.05),
      sei         = c(0.20, 0.24, 0.22, NA, 0.17),
      regexp      = "2 observation.*removed due to missing",
      expected_yi = c(0.10, 0.15, 0.05),
      expected_se = c(0.20, 0.22, 0.17)
    )
  )

  for (case in na_cases) {
    expect_warning(
      result <- brma.norm(
        yi        = case$yi,
        sei       = case$sei,
        only_data = TRUE
      )[["data"]],
      regexp = case$regexp
    )

    expect_norm_outcome(
      result,
      n   = length(case$expected_yi),
      yi  = case$expected_yi,
      sei = case$expected_se
    )
  }
})


test_that("NA in moderators drops rows with warning", {

  skip_on_cran()

  # Data with NA in moderator (row 3)
  test_mods <- data.frame(
    effect   = c(0.10, 0.25, 0.15, 0.30, 0.05),
    std_err  = c(0.20, 0.24, 0.22, 0.28, 0.17),
    mod_cont = c(1.5, 2.0, NA, 3.0, 1.0)
  )

  expect_warning(
    result <- brma.norm(
      yi   = effect,
      sei  = std_err,
      mods = ~ mod_cont,
      data = test_mods,
      only_data = TRUE
    )[["data"]],
    regexp = "1 observation.*removed due to missing"
  )

  expect_equal(nrow(result$outcome), 4)
  expect_equal(nrow(result$mods), 4)
  expect_equal(result$outcome$yi, c(0.10, 0.25, 0.30, 0.05))
  expect_equal(result$mods$mod_cont, c(1.5, 2.0, 3.0, 1.0))
})


test_that("NA in factor moderator drops rows with warning", {

  skip_on_cran()

  # Data with NA in factor moderator (row 4)
  test_mods <- data.frame(
    effect     = c(0.10, 0.25, 0.15, 0.30, 0.05),
    std_err    = c(0.20, 0.24, 0.22, 0.28, 0.17),
    mod_factor = factor(c("low", "high", "low", NA, "medium"))
  )

  expect_warning(
    result <- brma.norm(
      yi   = effect,
      sei  = std_err,
      mods = ~ mod_factor,
      data = test_mods,
      only_data = TRUE
    )[["data"]],
    regexp = "1 observation.*removed due to missing"
  )

  expect_equal(nrow(result$outcome), 4)
  expect_equal(nrow(result$mods), 4)

  # Check factor levels are properly dropped
  expect_equal(levels(result$mods$mod_factor), c("high", "low", "medium"))
})


test_that("NA in scale moderator drops rows with warning", {

  skip_on_cran()

  # Data with NA in scale variable (row 2)
  test_scale <- data.frame(
    effect    = c(0.10, 0.25, 0.15, 0.30, 0.05),
    std_err   = c(0.20, 0.24, 0.22, 0.28, 0.17),
    scale_var = c(0.5, NA, 0.7, 0.8, 0.6)
  )

  expect_warning(
    result <- brma.norm(
      yi    = effect,
      sei   = std_err,
      scale = ~ scale_var,
      data  = test_scale,
      only_data = TRUE
    )[["data"]],
    regexp = "1 observation.*removed due to missing"
  )

  expect_equal(nrow(result$outcome), 4)
  expect_equal(nrow(result$scale), 4)
  expect_equal(result$outcome$yi, c(0.10, 0.15, 0.30, 0.05))
  expect_equal(result$scale$scale_var, c(0.5, 0.7, 0.8, 0.6))
})


test_that("NAs across outcome, mods, and scale are all handled", {

  skip_on_cran()

  # Data with NAs in multiple places
  # yi NA in row 1, mod NA in row 3, scale NA in row 5
  test_combined <- data.frame(
    effect    = c(NA, 0.25, 0.15, 0.30, 0.05, 0.20),
    std_err   = c(0.20, 0.24, 0.22, 0.28, 0.17, 0.19),
    mod_cont  = c(1.5, 2.0, NA, 3.0, 1.0, 2.5),
    scale_var = c(0.5, 0.6, 0.7, 0.8, NA, 0.9)
  )

  expect_warning(
    result <- brma.norm(
      yi    = effect,
      sei   = std_err,
      mods  = ~ mod_cont,
      scale = ~ scale_var,
      data  = test_combined,
      only_data = TRUE
    )[["data"]],
    regexp = "3 observation.*removed due to missing"
  )

  # Should have rows 2, 4, 6 remaining
  expect_equal(nrow(result$outcome), 3)
  expect_equal(nrow(result$mods), 3)
  expect_equal(nrow(result$scale), 3)
  expect_equal(result$outcome$yi, c(0.25, 0.30, 0.20))
  expect_equal(result$mods$mod_cont, c(2.0, 3.0, 2.5))
  expect_equal(result$scale$scale_var, c(0.6, 0.8, 0.9))
})


test_that("NA in weights or cluster is rejected instead of dropped", {

  skip_on_cran()

  # Data with NA in weights
  expect_error(
    result <- brma.norm(
      yi        = c(0.10, 0.25, 0.15, 0.30, 0.05),
      sei       = c(0.20, 0.24, 0.22, 0.28, 0.17),
      weights   = c(1.0, 1.5, NA, 0.8, 1.3),
      cluster   = c("g1", "g1", "g2", "g2", "g3"),
      only_data = TRUE
    )[["data"]], regexp = "argument must not contain missing values"
  )

  # Data with NA in cluster
  expect_error(
    result <- brma.norm(
      yi        = c(0.10, 0.25, 0.15, 0.30, 0.05),
      sei       = c(0.20, 0.24, 0.22, 0.28, 0.17),
      weights   = c(1.0, 1.5, 1, 0.8, 1.3),
      cluster   = c("g1", NA, "g2", "g2", "g3"),
      only_data = TRUE
    )[["data"]], regexp = "argument must not contain missing values"
  )
})


test_that("Slab and cluster metadata are refreshed after NA dropping", {

  skip_on_cran()

  expect_warning(
    default_slab <- brma.norm(
      yi        = c(0.10, NA, 0.15, 0.30, 0.05),
      sei       = c(0.20, 0.24, 0.22, 0.28, 0.17),
      only_data = TRUE
    )[["data"]],
    regexp = "1 observation.*removed due to missing"
  )

  expect_warning(
    cluster_result <- brma.norm(
      yi        = c(NA, NA, 0.15, 0.30, 0.05),
      sei       = c(0.20, 0.24, 0.22, 0.28, 0.17),
      cluster   = c("g1", "g1", "g2", "g2", "g3"),
      only_data = TRUE
    )[["data"]],
    regexp = "2 observation.*removed due to missing"
  )

  expect_warning(
    user_slab <- brma.norm(
      yi        = c(0.10, NA, 0.15, 0.30, 0.05),
      sei       = c(0.20, 0.24, 0.22, 0.28, 0.17),
      slab      = c("A", "B", "C", "D", "E"),
      only_data = TRUE
    )[["data"]],
    regexp = "1 observation.*removed due to missing"
  )

  expect_equal(default_slab$outcome$slab, c("Study 1", "Study 2", "Study 3", "Study 4"))
  expect_equal(cluster_result$outcome$cluster, c(1, 1, 2))
  expect_equal(user_slab$outcome$slab, c("A", "C", "D", "E"))
})


test_that("All-NA yi is rejected", {

  skip_on_cran()

  # All NAs in yi should error (caught during validation by BayesTools::check_real)
  expect_error(
    brma.norm(
      yi   = c(NA, NA, NA),
      sei  = c(0.20, 0.24, 0.22),
      only_data = TRUE
    ),
    regexp = "numeric vector"
  )
})


test_that("NA handling rejects empty post-drop data", {

  skip_on_cran()

  # Each row has at least one NA
  expect_error(
    suppressWarnings(brma.norm(
      yi   = c(NA, 0.25, 0.15),
      sei  = c(0.20, NA, 0.22),
      mods = ~ mod_cont,
      data = data.frame(
        mod_cont = c(1.5, 2.0, NA)
      ),
      only_data = TRUE
    )),
    regexp = "No observations remaining"
  )
})


test_that("NA handling supports yi ~ mods formula syntax", {

  skip_on_cran()

  # Data with NA in formula moderator
  test_formula <- data.frame(
    effect   = c(0.10, 0.25, 0.15, 0.30, 0.05),
    std_err  = c(0.20, 0.24, 0.22, 0.28, 0.17),
    mod_cont = c(1.5, NA, 2.5, 3.0, 1.0)
  )

  expect_warning(
    result <- brma.norm(
      yi   = effect ~ mod_cont,
      sei  = std_err,
      data = test_formula,
      only_data = TRUE
    )[["data"]],
    regexp = "1 observation.*removed due to missing"
  )

  expect_equal(nrow(result$outcome), 4)
  expect_equal(nrow(result$mods), 4)
  expect_equal(result$outcome$yi, c(0.10, 0.15, 0.30, 0.05))
  expect_equal(result$mods$mod_cont, c(1.5, 2.5, 3.0, 1.0))
})


test_that("NA handling supports subset argument", {

  skip_on_cran()

  # Apply subset first, then NA dropping
  # Original data has 5 rows, subset keeps rows 1-4, NA in row 2 gets dropped
  expect_warning(
    result <- brma.norm(
      yi     = c(0.10, NA, 0.15, 0.30, 0.05),
      sei    = c(0.20, 0.24, 0.22, 0.28, 0.17),
      subset = 1:4,
      only_data = TRUE
    )[["data"]],
    regexp = "1 observation.*removed due to missing"
  )

  # After subset (4 rows) then NA drop (1 row), should have 3 rows
  expect_equal(nrow(result$outcome), 3)
  expect_equal(result$outcome$yi, c(0.10, 0.15, 0.30))
})

test_that("Numeric subset must use integer indices", {

  skip_on_cran()

  expect_error(
    brma.norm(
      yi        = c(0.10, 0.25, 0.15),
      sei       = c(0.20, 0.24, 0.22),
      subset    = c(1, 2.5),
      only_data = TRUE
    ),
    regexp = "integer numeric indices"
  )
})


test_that("Formula attributes are preserved after NA dropping", {

  skip_on_cran()

  test_mods <- data.frame(
    effect   = c(0.10, 0.25, 0.15, NA, 0.05),
    std_err  = c(0.20, 0.24, 0.22, 0.28, 0.17),
    mod_cont = c(1.5, 2.0, 2.5, 3.0, 1.0)
  )

  expect_warning(
    result <- brma.norm(
      yi   = effect,
      sei  = std_err,
      mods = ~ mod_cont,
      data = test_mods,
      only_data = TRUE
    )[["data"]],
    regexp = "1 observation.*removed due to missing"
  )

  # Formula attribute should be preserved
  expect_equal(
    paste0(as.character(attr(result$mods, "formula")), collapse = " "),
    "~ mod_cont"
  )
})


test_that("NA in subset argument is rejected", {

  skip_on_cran()

  # Logical subset with NA
  expect_error(
    brma.norm(
      yi     = c(0.10, 0.25, 0.15, 0.30, 0.05),
      sei    = c(0.20, 0.24, 0.22, 0.28, 0.17),
      subset = c(TRUE, NA, TRUE, TRUE, FALSE),
      only_data = TRUE
    ),
    regexp = "must not contain NA"
  )

  # Numeric subset with NA
  expect_error(
    brma.norm(
      yi     = c(0.10, 0.25, 0.15, 0.30, 0.05),
      sei    = c(0.20, 0.24, 0.22, 0.28, 0.17),
      subset = c(1, NA, 3),
      only_data = TRUE
    ),
    regexp = "must not contain NA"
  )
})
