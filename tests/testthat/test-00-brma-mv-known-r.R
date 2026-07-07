known_r_data <- function() {

  data.frame(
    yi  = c(0.10, 0.20, 0.30, 0.40),
    id  = factor(c("b", "a", "c", "b"), levels = c("b", "a", "c")),
    lab = factor(c("l1", "l1", "l2", "l2"), levels = c("l1", "l2")),
    x   = c(0, 1, 0, 1)
  )
}

known_r_kernel <- function(extra = FALSE) {

  levels <- c("a", "b", "c")
  K <- matrix(
    c(
      4.0, 1.0, 0.5,
      1.0, 9.0, 2.0,
      0.5, 2.0, 16.0
    ),
    nrow  = 3,
    byrow = TRUE,
    dimnames = list(levels, levels)
  )

  if (!extra) {
    return(K)
  }

  K_extra <- diag(4)
  dimnames(K_extra) <- list(c(levels, "extra"), c(levels, "extra"))
  K_extra[levels, levels] <- K
  K_extra
}


test_that("brma.mv normalizes metafor Rscale aliases", {

  expect_equal(.brma_mv_normalize_Rscale(TRUE, "id"), c(id = "cor"))
  expect_equal(.brma_mv_normalize_Rscale(FALSE, "id"), c(id = "none"))
  expect_equal(.brma_mv_normalize_Rscale(0, "id"), c(id = "none"))
  expect_equal(.brma_mv_normalize_Rscale(1, "id"), c(id = "cor"))
  expect_equal(.brma_mv_normalize_Rscale(2, "id"), c(id = "cor0"))
  expect_equal(.brma_mv_normalize_Rscale(3, "id"), c(id = "cov0"))
  expect_equal(.brma_mv_normalize_Rscale("none", c("id", "lab")),
               c(id = "none", lab = "none"))

  unnamed_Rscale <- c("none")
  names(unnamed_Rscale) <- NA_character_
  expect_equal(.brma_mv_normalize_Rscale(unnamed_Rscale, "id"), c(id = "none"))

  partially_named_Rscale <- c("cor", "none")
  names(partially_named_Rscale) <- c("id", NA_character_)
  expect_error(
    .brma_mv_normalize_Rscale(partially_named_Rscale, c("id", "lab")),
    "must be named"
  )

  expect_equal(
    .brma_mv_normalize_Rscale(
      list(id = "none", lab = "cor"),
      c("id", "lab")
    ),
    c(id = "none", lab = "cor")
  )

  expect_error(
    .brma_mv_normalize_Rscale("bad", "id"),
    "Rscale"
  )
  expect_error(
    .brma_mv_normalize_Rscale(c("cor", "none"), c("id", "lab")),
    "must be named"
  )
  expect_error(
    .brma_mv_normalize_Rscale(list("cor", "none"), c("id", "lab")),
    "must be named"
  )
})


test_that("brma.mv stores known R metadata on random-effect formulas", {

  dat <- known_r_data()
  K   <- known_r_kernel()

  object <- brma.mv(
    yi                        = yi,
    V                         = diag(0.01, nrow(dat)),
    random                    = ~ 1 | id,
    R                         = K,
    Rscale                    = "none",
    data                      = dat,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_data                 = TRUE
  )

  random_effects <- attr(object[["data"]][["location"]], "random_effects")
  term           <- random_effects[["terms"]][[1L]]
  formula_terms  <- attr(
    attr(object[["data"]][["location"]], "formula"),
    "random_terms"
  )

  expect_s3_class(random_effects, "BayesTools_random_effects")
  expect_s3_class(term[["group_covariance"]], "random_group_covariance")
  expect_equal(term[["group_covariance"]][["scale"]], "none")
  expect_s3_class(formula_terms[[1L]][["group_covariance"]],
                  "random_group_covariance")
})


test_that("brma.mv passes known R into BayesTools formula designs", {

  dat <- known_r_data()
  K   <- known_r_kernel(extra = TRUE)

  object <- brma.mv(
    yi                        = yi,
    V                         = diag(0.01, nrow(dat)),
    random                    = list(id_component = ~ 1 | id,
                                     lab_component = ~ 1 | lab),
    R                         = list(id = K),
    Rscale                    = list(id = "none"),
    data                      = dat,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )

  design <- .fitted_formula_design(object, "mu", required = TRUE)
  terms  <- design[["random_effects"]]
  names(terms) <- vapply(terms, `[[`, character(1), "group_label")
  id_group_covariance  <- terms[["id"]][["group_covariance"]]
  lab_group_covariance <- terms[["lab"]][["group_covariance"]]

  expect_s3_class(id_group_covariance, "random_group_covariance_kernel")
  expect_equal(id_group_covariance[["scale"]], "none")
  expect_equal(id_group_covariance[["levels"]], levels(dat[["id"]]))
  expect_equal(id_group_covariance[["dropped_levels"]], "extra")
  expect_equal(id_group_covariance[["kernel"]],
               K[levels(dat[["id"]]), levels(dat[["id"]])])
  expect_null(lab_group_covariance)
})


test_that("brma.mv known R appears in JAGS formula syntax and data", {

  dat <- known_r_data()
  K   <- known_r_kernel()

  object <- brma.mv(
    yi                        = yi,
    V                         = diag(0.01, nrow(dat)),
    random                    = ~ 1 | id,
    R                         = K,
    Rscale                    = "none",
    data                      = dat,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )

  result <- BayesTools::JAGS_formula(
    formula       = .create_fit_formula_list(object[["data"]], "location"),
    parameter     = "mu",
    data          = object[["data"]][["location"]],
    prior_list    = object[["priors"]][["location"]],
    formula_scale = .data_standardize_continuous_predictors(object[["data"]]),
    prior_random  = object[["priors"]][["random"]]
  )

  expect_match(result[["formula_syntax"]], "_xRE_GROUP_Zx", fixed = TRUE)
  expect_true(any(grepl("_xRE_GROUP_PRECx", names(result[["data"]]),
                        fixed = TRUE)))
  expect_false(any(grepl("_xRE_PRECx", result[["formula_syntax"]],
                         fixed = TRUE)))
})


test_that("brma.mv Rscale none preserves marginal covariance scale", {

  dat <- known_r_data()
  K   <- known_r_kernel()

  object <- brma.mv(
    yi                        = yi,
    V                         = diag(0.01, nrow(dat)),
    random                    = ~ 1 | id,
    R                         = K,
    Rscale                    = "none",
    data                      = dat,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  result <- BayesTools::JAGS_formula(
    formula       = .create_fit_formula_list(object[["data"]], "location"),
    parameter     = "mu",
    data          = object[["data"]][["location"]],
    prior_list    = object[["priors"]][["location"]],
    formula_scale = .data_standardize_continuous_predictors(object[["data"]]),
    prior_random  = object[["priors"]][["random"]]
  )
  random_term <- result[["formula_design"]][["random_effects"]][[1L]]
  posterior <- matrix(
    0.5,
    nrow = 1,
    ncol = length(random_term[["sd_parameter_names"]]),
    dimnames = list(NULL, random_term[["sd_parameter_names"]])
  )
  out <- BayesTools::random_effects_marginal_vcov(
    result[["formula_design"]],
    posterior_samples = posterior,
    prior_list        = result[["prior_list"]]
  )

  expected <- unname(0.5^2 * K[as.character(dat[["id"]]),
                               as.character(dat[["id"]])])
  expect_equal(unname(out[["samples"]][1L, , ]), expected, tolerance = 1e-12)
  expect_equal(unname(diag(out[["samples"]][1L, , ])),
               0.5^2 * c(9, 4, 16, 9))
})


test_that("brma.mv rejects unsupported known R inputs", {

  dat <- known_r_data()
  K   <- known_r_kernel()

  expect_error(
    brma.mv(
      yi                        = yi,
      V                         = diag(0.01, nrow(dat)),
      R                         = K,
      data                      = dat,
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_data                 = TRUE
    ),
    "requires a 'random' formula"
  )
  expect_error(
    brma.mv(
      yi                        = yi,
      V                         = diag(0.01, nrow(dat)),
      random                    = ~ 1 | id,
      R                         = K,
      Rscale                    = "bad",
      data                      = dat,
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_data                 = TRUE
    ),
    "Rscale"
  )
  expect_error(
    brma.mv(
      yi                        = yi,
      V                         = diag(0.01, nrow(dat)),
      random                    = ~ diag(1 + x | id),
      R                         = K,
      data                      = dat,
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_data                 = TRUE
    ),
    "random slopes"
  )
  expect_error(
    brma.mv(
      yi                        = yi,
      V                         = diag(0.01, nrow(dat)),
      random                    = ~ cs(1 | id),
      R                         = K,
      data                      = dat,
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_data                 = TRUE
    ),
    "structure 'cs'"
  )
  expect_error(
    brma.mv(
      yi                        = yi,
      V                         = diag(0.01, nrow(dat)),
      random                    = ~ 1 | id,
      R                         = list(K, K),
      data                      = dat,
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_data                 = TRUE
    ),
    "fully named"
  )
})


test_that("brma.mv delegates known R matching and level validation to BayesTools", {

  dat       <- known_r_data()
  K         <- known_r_kernel()
  K_missing <- K[c("a", "b"), c("a", "b")]
  unmatched_message <- tryCatch(
    brma.mv(
      yi                        = yi,
      V                         = diag(0.01, nrow(dat)),
      random                    = ~ 1 | id,
      R                         = list(missing = K),
      data                      = dat,
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_data                 = TRUE
    ),
    error = function(e) conditionMessage(e)
  )

  expect_error(
    brma.mv(
      yi                        = yi,
      V                         = diag(0.01, nrow(dat)),
      random                    = list(one = ~ 1 | id, two = ~ 1 | id),
      R                         = list(id = K),
      data                      = dat,
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_data                 = TRUE
    ),
    "ambiguous"
  )
  expect_match(unmatched_message, "'R'")
  expect_false(grepl("group_covariance", unmatched_message, fixed = TRUE))
  expect_error(
    brma.mv(
      yi                        = yi,
      V                         = diag(0.01, nrow(dat)),
      random                    = ~ 1 | id,
      R                         = list(id = K_missing),
      data                      = dat,
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_priors               = TRUE
    ),
    "missing fitted level"
  )
})


test_that("brma.mv marginalizes supported one-to-one known R blocks", {

  dat <- data.frame(
    yi       = c(0.10, 0.20, 0.30),
    estimate = factor(c("e1", "e2", "e3"), levels = c("e1", "e2", "e3"))
  )
  K <- diag(c(4, 9, 16))
  dimnames(K) <- list(levels(dat[["estimate"]]), levels(dat[["estimate"]]))

  known_r <- brma.mv(
    yi                        = yi,
    V                         = diag(0.01, nrow(dat)),
    random                    = ~ 1 | estimate,
    R                         = K,
    Rscale                    = "none",
    data                      = dat,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  known_r_term <- .fitted_formula_design(known_r, "mu")[["random_effects"]][[1L]]
  known_r_terms <- .data_marginalized_random_effects(known_r[["data"]])
  fit_data      <- .create_fit_data(known_r[["data"]], known_r[["priors"]])
  syntax        <- .create_model_syntax(known_r[["data"]], known_r[["priors"]])

  result <- BayesTools::JAGS_formula(
    formula                = .create_fit_formula_list(known_r[["data"]], "location"),
    parameter              = "mu",
    data                   = known_r[["data"]][["location"]],
    prior_list             = known_r[["priors"]][["location"]],
    formula_scale          = .data_standardize_continuous_predictors(known_r[["data"]]),
    prior_random           = known_r[["priors"]][["random"]],
    random_effects_compile = .data_random_effects_compile(known_r[["data"]])
  )

  ordinary <- brma.mv(
    yi                        = yi,
    V                         = diag(0.01, nrow(dat)),
    random                    = ~ 1 | estimate,
    data                      = dat,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  ordinary_term <- .fitted_formula_design(ordinary, "mu")[["random_effects"]][[1L]]

  expect_equal(known_r_term[["compile_mode"]], "marginalized")
  expect_true(.data_has_marginalized_random_effects(known_r[["data"]]))
  expect_true(.data_has_marginalized_known_group_covariance(known_r[["data"]]))
  expect_length(known_r_terms, 1L)
  expect_equal(unname(known_r_terms[[1L]][["row_multiplier"]]), c(4, 9, 16))
  expect_equal(
    fit_data[[known_r_terms[[1L]][["row_multiplier_name"]]]],
    known_r_terms[[1L]][["row_multiplier"]]
  )
  expect_match(syntax, "mu__xRE_MVARx_estimate\\[i\\]")
  expect_false(grepl("_xRE_GROUP_Zx", result[["formula_syntax"]], fixed = TRUE))
  expect_false(any(grepl("_xRE_GROUP_PRECx", names(result[["data"]]),
                         fixed = TRUE)))
  expect_error(
    predict(known_r, newdata = dat, type = "estimate"),
    "new or reordered known-R rows"
  )
  expect_equal(ordinary_term[["compile_mode"]], "marginalized")
  expect_true(.data_has_marginalized_random_effects(ordinary[["data"]]))
})


test_that("brma.mv sampled known R path is preserved when requested", {

  dat <- data.frame(
    yi       = c(0.10, 0.20, 0.30),
    estimate = factor(c("e1", "e2", "e3"), levels = c("e1", "e2", "e3"))
  )
  K <- diag(c(4, 9, 16))
  dimnames(K) <- list(levels(dat[["estimate"]]), levels(dat[["estimate"]]))

  object <- brma.mv(
    yi                        = yi,
    V                         = diag(0.01, nrow(dat)),
    random                    = ~ 1 | estimate,
    R                         = K,
    Rscale                    = "none",
    data                      = dat,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    marginalize_estimate_level  = FALSE,
    only_priors               = TRUE
  )
  term <- .fitted_formula_design(object, "mu")[["random_effects"]][[1L]]
  result <- BayesTools::JAGS_formula(
    formula       = .create_fit_formula_list(object[["data"]], "location"),
    parameter     = "mu",
    data          = object[["data"]][["location"]],
    prior_list    = object[["priors"]][["location"]],
    formula_scale = .data_standardize_continuous_predictors(object[["data"]]),
    prior_random  = object[["priors"]][["random"]]
  )

  expect_equal(term[["compile_mode"]], "sampled")
  expect_false(.data_has_marginalized_random_effects(object[["data"]]))
  expect_match(result[["formula_syntax"]], "_xRE_GROUP_Zx", fixed = TRUE)
  expect_true(any(grepl("_xRE_GROUP_PRECx", names(result[["data"]]),
                        fixed = TRUE)))
})


test_that("brma.mv known R row multipliers respect Rscale", {

  dat <- data.frame(
    yi       = c(0.10, 0.20, 0.30),
    estimate = factor(c("e1", "e2", "e3"), levels = c("e1", "e2", "e3"))
  )
  K <- diag(c(4, 9, 16))
  dimnames(K) <- list(levels(dat[["estimate"]]), levels(dat[["estimate"]]))

  expected <- list(
    none = c(4, 9, 16),
    cor  = c(1, 1, 1),
    cor0 = c(1, 1, 1),
    cov0 = c(4, 9, 16)
  )

  for (Rscale in names(expected)) {
    object <- brma.mv(
      yi                        = yi,
      V                         = diag(0.01, nrow(dat)),
      random                    = ~ 1 | estimate,
      R                         = K,
      Rscale                    = Rscale,
      data                      = dat,
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_priors               = TRUE
    )
    term <- .data_marginalized_random_effects(object[["data"]])[[1L]]

    expect_equal(unname(term[["row_multiplier"]]), expected[[Rscale]],
                 info = Rscale)
  }
})


test_that("brma.mv keeps unsupported known R marginalization designs sampled", {

  repeated <- known_r_data()
  repeated_R <- known_r_kernel()

  repeated_object <- brma.mv(
    yi                        = yi,
    V                         = diag(0.01, nrow(repeated)),
    random                    = ~ 1 | id,
    R                         = repeated_R,
    Rscale                    = "none",
    data                      = repeated,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  repeated_term <- .fitted_formula_design(
    repeated_object,
    "mu"
  )[["random_effects"]][[1L]]

  offdiag <- data.frame(
    yi       = c(0.10, 0.20, 0.30),
    estimate = factor(c("e1", "e2", "e3"), levels = c("e1", "e2", "e3"))
  )
  offdiag_R <- matrix(
    c(
      4, 1, 0,
      1, 9, 0,
      0, 0, 16
    ),
    nrow  = 3,
    byrow = TRUE,
    dimnames = list(levels(offdiag[["estimate"]]),
                    levels(offdiag[["estimate"]]))
  )
  offdiag_object <- brma.mv(
    yi                        = yi,
    V                         = diag(0.01, nrow(offdiag)),
    random                    = ~ 1 | estimate,
    R                         = offdiag_R,
    Rscale                    = "none",
    data                      = offdiag,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  offdiag_term <- .fitted_formula_design(
    offdiag_object,
    "mu"
  )[["random_effects"]][[1L]]

  expect_equal(repeated_term[["compile_mode"]], "sampled")
  expect_false(.data_has_marginalized_random_effects(repeated_object[["data"]]))
  expect_equal(offdiag_term[["compile_mode"]], "sampled")
  expect_false(.data_has_marginalized_random_effects(offdiag_object[["data"]]))
})
