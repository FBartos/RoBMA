context("Input handling for predict")

source(testthat::test_path("helper-contracts.R"))

skip_on_cran()

test_data_norm <- data.frame(
  yi        = c(0.2, 0.5, -0.1, 0.3, 0.4),
  sei       = c(0.1, 0.15, 0.12, 0.08, 0.11),
  mod_cont  = c(1.5, 2.3, 1.8, 3.1, 0.9),
  mod_fac   = factor(c("A", "B", "A", "B", "A")),
  scale_var = c(0.5, 1.0, 0.8, 1.2, 0.6),
  stringsAsFactors = FALSE
)

test_data_norm_vi <- data.frame(
  yi       = c(0.2, 0.5, -0.1, 0.3, 0.4),
  vi       = c(0.01, 0.0225, 0.0144, 0.0064, 0.0121),
  mod_cont = c(1.5, 2.3, 1.8, 3.1, 0.9),
  mod_fac  = factor(c("A", "B", "A", "B", "A")),
  stringsAsFactors = FALSE
)

test_data_glmm <- data.frame(
  ai       = c(10L, 15L, 12L, 8L, 20L),
  ci       = c(5L, 10L, 8L, 4L, 12L),
  n1i      = c(50L, 50L, 50L, 50L, 50L),
  n2i      = c(50L, 50L, 50L, 50L, 50L),
  mod_cont = c(1.5, 2.3, 1.8, 3.1, 0.9),
  mod_fac  = factor(c("A", "B", "A", "B", "A")),
  stringsAsFactors = FALSE
)

compare_data_lists <- function(data1, data2, check_slab = FALSE,
                               structure_only = FALSE, info = NULL) {

  expect_s3_class(data1, "RoBMA_data")
  expect_s3_class(data2, "RoBMA_data")

  if (!structure_only) {
    expect_equal(nrow(data1$outcome), nrow(data2$outcome), info = info)
  }
  expect_equal(ncol(data1$outcome), ncol(data2$outcome), info = info)
  expect_equal(names(data1$outcome), names(data2$outcome), info = info)

  if (!structure_only) {
    outcome_cols <- setdiff(names(data1$outcome), if (check_slab) character() else "slab")
    for (col in outcome_cols) {
      expect_equal(data1$outcome[[col]], data2$outcome[[col]],
                   info = paste(info, "outcome", col))
    }
  }

  for (slot in c("mods", "scale")) {
    expect_equal(is.null(data1[[slot]]), is.null(data2[[slot]]), info = info)
    if (!is.null(data1[[slot]])) {
      if (!structure_only) {
        expect_equal(nrow(data1[[slot]]), nrow(data2[[slot]]), info = info)
      }
      expect_equal(names(data1[[slot]]), names(data2[[slot]]), info = info)
      if (!structure_only) {
        for (col in names(data1[[slot]])) {
          expect_equal(data1[[slot]][[col]], data2[[slot]][[col]],
                       info = paste(info, slot, col))
        }
      }
      expect_equal(attr(data1[[slot]], "formula"), attr(data2[[slot]], "formula"),
                   info = paste(info, slot, "formula"))
    }
  }

  attrs <- c(
    "outcome_type", "mods", "scale",
    "standardize_continuous_predictors",
    "set_contrast_factor_predictors", "effect_direction"
  )
  for (attr_name in attrs) {
    expect_equal(attr(data1, attr_name), attr(data2, attr_name),
                 info = paste(info, "attribute", attr_name))
  }
  if (!structure_only) {
    expect_equal(attr(data1, "k_final"), attr(data2, "k_final"),
                 info = paste(info, "attribute k_final"))
  }
}

prepare_newdata_cases <- list(
  list(
    label = "normal same data with sei",
    fit = quote(brma.norm(
      yi = yi, sei = sei, data = test_data_norm, only_data = TRUE
    )),
    newdata = quote(test_data_norm)
  ),
  list(
    label = "normal same data with vi",
    fit = quote(brma.norm(
      yi = yi, vi = vi, data = test_data_norm_vi, only_data = TRUE
    )),
    newdata = quote(test_data_norm_vi)
  ),
  list(
    label = "normal different outcome rows",
    fit = quote(brma.norm(
      yi = yi, sei = sei, data = test_data_norm, only_data = TRUE
    )),
    newdata = quote(data.frame(yi = c(0.1, 0.6, 0.0), sei = c(0.2, 0.1, 0.15))),
    structure_only = TRUE
  ),
  list(
    label = "normal moderators same data",
    fit = quote(brma.norm(
      yi = yi, sei = sei, mods = ~ mod_cont + mod_fac,
      data = test_data_norm, only_data = TRUE
    )),
    newdata = quote(test_data_norm)
  ),
  list(
    label = "normal moderators different rows",
    fit = quote(brma.norm(
      yi = yi, sei = sei, mods = ~ mod_cont + mod_fac,
      data = test_data_norm, only_data = TRUE
    )),
    newdata = quote(data.frame(
      yi       = c(0.1, 0.3),
      sei      = c(0.1, 0.2),
      mod_cont = c(2.0, 1.0),
      mod_fac  = factor(c("A", "B"), levels = c("A", "B"))
    )),
    structure_only = TRUE
  ),
  list(
    label = "normal scale same data",
    fit = quote(brma.norm(
      yi = yi, sei = sei, scale = ~ scale_var,
      data = test_data_norm, only_data = TRUE
    )),
    newdata = quote(test_data_norm)
  ),
  list(
    label = "normal moderators and scale",
    fit = quote(brma.norm(
      yi = yi, sei = sei, mods = ~ mod_cont, scale = ~ scale_var,
      data = test_data_norm, only_data = TRUE
    )),
    newdata = quote(test_data_norm)
  ),
  list(
    label = "GLMM same data",
    fit = quote(brma.glmm(
      ai = ai, ci = ci, n1i = n1i, n2i = n2i,
      data = test_data_glmm, only_data = TRUE
    )),
    newdata = quote(test_data_glmm)
  ),
  list(
    label = "GLMM moderators same data",
    fit = quote(brma.glmm(
      ai = ai, ci = ci, n1i = n1i, n2i = n2i, mods = ~ mod_cont + mod_fac,
      data = test_data_glmm, only_data = TRUE
    )),
    newdata = quote(test_data_glmm)
  )
)

test_that(".prepare_newdata reconstructs response, moderator, and scale data", {

  for (case in prepare_newdata_cases) {
    fit <- eval(case[["fit"]])
    result <- RoBMA:::.prepare_newdata(
      object  = fit,
      newdata = eval(case[["newdata"]]),
      type    = "terms"
    )

    compare_data_lists(
      fit[["data"]],
      result,
      structure_only = isTRUE(case[["structure_only"]]),
      info           = case[["label"]]
    )
  }
})

test_that(".prepare_newdata inserts dummy outcomes only when the response is unused", {

  fit <- brma.norm(
    yi = yi, sei = sei, mods = ~ mod_cont,
    data = test_data_norm, only_data = TRUE
  )

  result <- RoBMA:::.prepare_newdata(
    object  = fit,
    newdata = data.frame(mod_cont = c(1.2, 2.4)),
    type    = "terms"
  )

  expect_equal(nrow(result[["outcome"]]), 2L)
  expect_equal(result[["outcome"]][["yi"]], c(0, 0))
  expect_equal(result[["outcome"]][["sei"]], c(0, 0))

  response_fit <- brma.norm(
    yi = yi, sei = sei, data = test_data_norm, only_data = TRUE
  )
  response <- RoBMA:::.prepare_newdata(
    object  = response_fit,
    newdata = data.frame(sei = c(0.1, 0.2)),
    type    = "response"
  )
  expect_equal(response[["outcome"]][["yi"]], c(0, 0))
  expect_equal(response[["outcome"]][["sei"]], c(0.1, 0.2))
})

test_that(".prepare_newdata rejects missing required variables", {

  fit_norm <- brma.norm(
    yi = yi, sei = sei, data = test_data_norm, only_data = TRUE
  )
  fit_mods <- brma.norm(
    yi = yi, sei = sei, mods = ~ mod_cont + mod_fac,
    data = test_data_norm, only_data = TRUE
  )
  fit_scale <- brma.norm(
    yi = yi, sei = sei, scale = ~ scale_var,
    data = test_data_norm, only_data = TRUE
  )
  fit_glmm <- brma.glmm(
    ai = ai, ci = ci, n1i = n1i, n2i = n2i,
    data = test_data_glmm, only_data = TRUE
  )

  expect_error_cases(list(
    list(
      label  = "normal response missing sei/vi",
      expr   = quote(RoBMA:::.prepare_newdata(
        object = fit_norm, newdata = data.frame(yi = c(0.1, 0.2)),
        type = "response"
      )),
      regexp = "sei.*vi"
    ),
    list(
      label  = "missing moderator",
      expr   = quote(RoBMA:::.prepare_newdata(
        object = fit_mods,
        newdata = data.frame(
          yi = c(0.1, 0.2), sei = c(0.1, 0.2),
          mod_fac = factor(c("A", "B"))
        ),
        type = "terms"
      )),
      regexp = "mod_cont"
    ),
    list(
      label  = "missing scale predictor",
      expr   = quote(RoBMA:::.prepare_newdata(
        object = fit_scale,
        newdata = data.frame(yi = c(0.1, 0.2), sei = c(0.1, 0.2)),
        type = "terms"
      )),
      regexp = "scale_var"
    ),
    list(
      label  = "missing GLMM response size",
      expr   = quote(RoBMA:::.prepare_newdata(
        object = fit_glmm,
        newdata = data.frame(ai = c(10L, 15L), ci = c(5L, 10L), n1i = c(50L, 50L)),
        type = "response"
      )),
      regexp = "n2i"
    )
  ))
})

test_that(".prepare_newdata permits bias-adjusted PET/PEESE terms without new standard errors", {

  fit <- bPET(
    yi = yi, sei = sei, data = test_data_norm, only_data = TRUE
  )
  fit[["priors"]] <- list(outcome = list(
    bias = BayesTools::prior_PET("normal", list(mean = 0, sd = 1))
  ))
  new_df <- data.frame(row = 1:2)

  expect_error(
    RoBMA:::.prepare_newdata(
      object = fit, newdata = new_df, type = "terms"
    ),
    regexp = "sei.*vi"
  )
  expect_no_error(RoBMA:::.prepare_newdata(
    object = fit, newdata = new_df, type = "terms", bias_adjusted = TRUE
  ))
})

test_that(".prepare_newdata preserves predictor transformation settings", {

  fit_std <- brma.norm(
    yi = yi, sei = sei, mods = ~ mod_cont,
    data = test_data_norm,
    standardize_continuous_predictors = TRUE,
    only_data = TRUE
  )
  fit_no_std <- brma.norm(
    yi = yi, sei = sei, mods = ~ mod_cont,
    data = test_data_norm,
    standardize_continuous_predictors = FALSE,
    only_data = TRUE
  )
  result_std <- RoBMA:::.prepare_newdata(fit_std, test_data_norm, type = "terms")
  result_no_std <- RoBMA:::.prepare_newdata(fit_no_std, test_data_norm, type = "terms")

  expect_true(attr(result_std, "standardize_continuous_predictors"))
  expect_false(attr(result_no_std, "standardize_continuous_predictors"))

  fit_treatment <- brma.norm(
    yi = yi, sei = sei, mods = ~ mod_fac,
    data = test_data_norm,
    set_contrast_factor_predictors = "treatment",
    only_data = TRUE
  )
  fit_meandif <- brma.norm(
    yi = yi, sei = sei, mods = ~ mod_fac,
    data = test_data_norm,
    set_contrast_factor_predictors = "meandif",
    only_data = TRUE
  )

  result_treatment <- RoBMA:::.prepare_newdata(fit_treatment, test_data_norm, type = "terms")
  result_meandif <- RoBMA:::.prepare_newdata(fit_meandif, test_data_norm, type = "terms")

  expect_equal(attr(result_treatment, "set_contrast_factor_predictors"), "treatment")
  expect_equal(attr(result_meandif, "set_contrast_factor_predictors"), "meandif")
})

test_that(".prepare_newdata drops rows with missing outcome or moderator values", {

  fit_norm <- brma.norm(
    yi = yi, sei = sei, data = test_data_norm, only_data = TRUE
  )
  expect_warning(
    result_norm <- RoBMA:::.prepare_newdata(
      object = fit_norm,
      newdata = data.frame(yi = c(0.1, NA, 0.3), sei = c(0.1, 0.2, 0.15)),
      type = "terms"
    ),
    regexp = "removed"
  )
  expect_equal(nrow(result_norm[["outcome"]]), 2L)

  fit_mods <- brma.norm(
    yi = yi, sei = sei, mods = ~ mod_cont,
    data = test_data_norm, only_data = TRUE
  )
  expect_warning(
    result_mods <- RoBMA:::.prepare_newdata(
      object = fit_mods,
      newdata = data.frame(
        yi = c(0.1, 0.2, 0.3), sei = c(0.1, 0.2, 0.15),
        mod_cont = c(1.5, NA, 2.0)
      ),
      type = "terms"
    ),
    regexp = "removed"
  )
  expect_equal(nrow(result_mods[["outcome"]]), 2L)
  expect_equal(nrow(result_mods[["mods"]]), 2L)
})

test_that(".prepare_newdata accepts prediction-only edge cases", {

  edge_cases <- list(
    list(
      label = "single row",
      fit = quote(brma.norm(
        yi = yi, sei = sei, data = test_data_norm, only_data = TRUE
      )),
      newdata = quote(data.frame(yi = 0.5, sei = 0.1)),
      check = function(result) {
        expect_equal(nrow(result[["outcome"]]), 1L)
        expect_equal(result[["outcome"]][["yi"]], 0.5)
      }
    ),
    list(
      label = "extra columns",
      fit = quote(brma.norm(
        yi = yi, sei = sei, data = test_data_norm, only_data = TRUE
      )),
      newdata = quote(data.frame(
        yi = c(0.1, 0.2), sei = c(0.1, 0.2),
        extra_col = c("a", "b"), another = c(1, 2)
      )),
      check = function(result) expect_equal(nrow(result[["outcome"]]), 2L)
    ),
    list(
      label = "single-level factor",
      fit = quote(brma.norm(
        yi = yi, sei = sei, mods = ~ mod_fac,
        data = test_data_norm, only_data = TRUE
      )),
      newdata = quote(data.frame(
        yi = c(0.1, 0.2, 0.3), sei = c(0.1, 0.1, 0.1),
        mod_fac = factor(c("A", "A", "A"), levels = c("A", "B"))
      )),
      check = function(result) {
        expect_equal(nrow(result[["mods"]]), 3L)
        expect_equal(as.character(result[["mods"]][["mod_fac"]]), rep("A", 3))
      }
    ),
    list(
      label = "zero-variance moderator",
      fit = quote(brma.norm(
        yi = yi, sei = sei, mods = ~ mod_cont,
        data = test_data_norm, only_data = TRUE
      )),
      newdata = quote(data.frame(
        yi = c(0.1, 0.2, 0.3), sei = c(0.1, 0.1, 0.1),
        mod_cont = c(2, 2, 2)
      )),
      check = function(result) expect_equal(result[["mods"]][["mod_cont"]], c(2, 2, 2))
    ),
    list(
      label = "single observation with mixed moderators",
      fit = quote(brma.norm(
        yi = yi, sei = sei, mods = ~ mod_cont + mod_fac,
        data = test_data_norm, only_data = TRUE
      )),
      newdata = quote(data.frame(
        yi = 0.25, sei = 0.1, mod_cont = 2.5,
        mod_fac = factor("B", levels = c("A", "B"))
      )),
      check = function(result) {
        expect_equal(nrow(result[["mods"]]), 1L)
        expect_equal(as.character(result[["mods"]][["mod_fac"]]), "B")
      }
    ),
    list(
      label = "zero-variance scale predictor",
      fit = quote(brma.norm(
        yi = yi, sei = sei, scale = ~ scale_var,
        data = test_data_norm, only_data = TRUE
      )),
      newdata = quote(data.frame(
        yi = c(0.1, 0.2), sei = c(0.1, 0.1), scale_var = c(1, 1)
      )),
      check = function(result) expect_equal(result[["scale"]][["scale_var"]], c(1, 1))
    ),
    list(
      label = "zero sei",
      fit = quote(brma.norm(
        yi = yi, sei = sei, data = test_data_norm, only_data = TRUE
      )),
      newdata = quote(data.frame(yi = c(0.1, 0.2, 0.3), sei = c(0, 0.1, 0))),
      check = function(result) expect_equal(result[["outcome"]][["sei"]], c(0, 0.1, 0))
    ),
    list(
      label = "zero vi",
      fit = quote(brma.norm(
        yi = yi, vi = vi, data = test_data_norm_vi, only_data = TRUE
      )),
      newdata = quote(data.frame(yi = c(0.1, 0.2), vi = c(0, 0))),
      check = function(result) expect_equal(result[["outcome"]][["sei"]], c(0, 0))
    ),
    list(
      label = "optional ni omitted",
      fit = quote(brma.norm(
        yi = yi, sei = sei, ni = ni,
        data = data.frame(yi = c(0.2, 0.5, -0.1),
                          sei = c(0.1, 0.15, 0.12),
                          ni = c(50, 100, 75)),
        only_data = TRUE
      )),
      newdata = quote(data.frame(yi = c(0.1, 0.2), sei = c(0, 0.1))),
      check = function(result) expect_equal(result[["outcome"]][["sei"]], c(0, 0.1))
    )
  )

  for (case in edge_cases) {
    fit <- eval(case[["fit"]])
    result <- RoBMA:::.prepare_newdata(
      object  = fit,
      newdata = eval(case[["newdata"]]),
      type    = "terms"
    )
    case[["check"]](result)
  }
})
