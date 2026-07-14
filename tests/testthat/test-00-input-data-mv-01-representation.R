context("brma.mv known-V representation and validation")

test_that("known-V consumers require the current representation", {

  data <- list(outcome = data.frame(yi = 0, sei = 1))
  attr(data, "known_V")      <- TRUE
  attr(data, "known_V_data") <- list(V = matrix(1, nrow = 1L, ncol = 1L))

  expect_error(
    .data_known_v_data(data),
    "known-V representation is incomplete"
  )

  attr(data, "known_V_data") <- .known_v_canonicalize(1)
  expect_error(
    .data_known_v_data(data),
    "prepared known-V metadata are incomplete"
  )
})

test_that("brma.mv stores and decomposes known V", {

  V <- matrix(
    c(
      0.04, 0.01, 0.00,
      0.01, 0.09, 0.00,
      0.00, 0.00, 0.16
    ),
    nrow = 3,
    byrow = TRUE
  )

  object <- brma.mv(
    yi                        = c(0.10, 0.20, 0.30),
    V                         = V,
    known_v_parameterization  = "latent",
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_data                 = TRUE
  )
  data <- object[["data"]]

  expect_true(attr(data, "known_V"))
  expect_equal(data[["outcome"]][["sei"]], sqrt(diag(V)))

  known_V <- attr(data, "known_V_data")
  expect_equal(.known_v_materialize(known_V), V)
  reconstructed <- diag(known_V[["residual_variance"]], nrow = 3)
  for (block in known_V[["latent_blocks"]]) {
    index <- block[["index"]]
    reconstructed[index, index] <- reconstructed[index, index] +
      tcrossprod(block[["B"]])
  }
  expect_equal(
    reconstructed,
    V,
    tolerance = 1e-10
  )
  expect_equal(
    known_V[["rank"]],
    sum(vapply(known_V[["latent_blocks"]], `[[`, integer(1), "rank"))
  )
  expect_null(known_V[["B"]])
  expect_null(known_V[["sampling_factor"]])
})


test_that("brma.mv supports list V input", {

  V1 <- matrix(c(0.04, 0.01, 0.01, 0.09), nrow = 2)
  V2 <- matrix(0.16, nrow = 1)

  object <- brma.mv(
    yi                        = c(0.10, 0.20, 0.30),
    V                         = list(V1, V2),
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_data                 = TRUE
  )

  known_V <- .data_known_v_data(object[["data"]])
  expect_equal(
    .known_v_materialize(known_V),
    rbind(cbind(V1, matrix(0, 2, 1)), cbind(matrix(0, 1, 2), V2))
  )
  expect_equal(known_V[["storage"]], "blocks")
  expect_false(any(c("V", "B", "whitening_matrix", "sampling_factor") %in%
                   names(known_V)))
})


test_that("connected dense known V is retained without duplicate block storage", {

  V <- matrix(
    c(
      0.04, 0.01, 0.005,
      0.01, 0.09, 0.015,
      0.005, 0.015, 0.16
    ),
    nrow = 3L,
    byrow = TRUE
  )
  object <- brma.mv(
    yi                        = c(0.10, 0.20, 0.30),
    V                         = V,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_data                 = TRUE
  )
  known_V <- .data_known_v_data(object[["data"]])

  expect_equal(known_V[["storage"]], "dense")
  expect_equal(known_V[["V"]], V)
  expect_null(known_V[["blocks"]])
  expect_equal(known_V[["block_indices"]], list(seq_len(nrow(V))))
  expect_identical(
    .known_v_correlated_blocks(known_V)[[1L]][["covariance"]],
    known_V[["V"]]
  )
})


test_that("known V list input rejects empty blocks", {

  expect_error(
    brma.mv(
      yi                        = c(0.10, 0.20),
      V                         = list(matrix(numeric(0), nrow = 0L, ncol = 0L)),
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_data                 = TRUE
    ),
    "non-empty square"
  )
})


test_that("known V validates dense structure before row filtering", {

  expect_error(
    brma.mv(
      yi        = c(0.10, 0.20),
      V         = matrix(letters[1:4], nrow = 2L),
      measure   = "GEN",
      only_data = TRUE
    ),
    "'V' argument must be numeric",
    fixed = TRUE
  )
  expect_error(
    brma.mv(
      yi        = c(0.10, 0.20, 0.30),
      V         = matrix(seq_len(6), nrow = 3L),
      measure   = "GEN",
      only_data = TRUE
    ),
    "'V' argument must be a non-empty square matrix",
    fixed = TRUE
  )
})


test_that("known V symmetry tolerance scales with covariance magnitude", {

  V <- matrix(
    c(1e8, 1e4, 1e4 + 1e-4, 2e8),
    nrow = 2L
  )
  object <- brma.mv(
    yi                        = c(0.10, 0.20),
    V                         = V,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_data                 = TRUE
  )

  expect_equal(
    .known_v_materialize(.data_known_v_data(object[["data"]])),
    (V + t(V)) / 2
  )
})


test_that("V_new list input rejects empty blocks", {

  expect_error(
    .known_v_newdata_prepare(
      V_new = list(matrix(numeric(0), nrow = 0L, ncol = 0L)),
      k     = 0L
    ),
    "V_new.*non-empty square"
  )
})


test_that("V_new structural errors name V_new", {

  expect_error(
    .known_v_newdata_prepare(V_new = "invalid", k = 1L),
    "'V_new' argument must be a variance vector",
    fixed = TRUE
  )
  expect_error(
    .known_v_newdata_prepare(
      V_new = matrix(letters[1:4], nrow = 2L),
      k     = 2L
    ),
    "'V_new' argument must be numeric",
    fixed = TRUE
  )
})


test_that("V_new symmetry tolerance scales with covariance magnitude", {

  V_new <- matrix(
    c(1e8, 1e4, 1e4 + 1e-4, 2e8),
    nrow = 2L
  )
  out <- .known_v_newdata_prepare(V_new, k = 2L)

  expect_equal(.known_v_materialize(out), (V_new + t(V_new)) / 2)
})


test_that("brma.mv supports variance vector V input", {

  vi <- c(0.04, 0.09, 0.16)
  object <- brma.mv(
    yi                        = c(0.10, 0.20, 0.30),
    V                         = vi,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_data                 = TRUE
  )

  known_V <- .data_known_v_data(object[["data"]])
  expect_equal(.known_v_materialize(known_V), diag(vi))
  expect_equal(object[["data"]][["outcome"]][["sei"]], sqrt(vi))
  expect_equal(known_V[["storage"]], "diagonal")
  expect_false(any(c("V", "B", "whitening_matrix", "sampling_factor") %in%
                   names(known_V)))
})


test_that("known V base covariance preserves singleton draw dimensions", {

  known_V <- .known_v_prepare(
    V                         = 0.04,
    keep_rows                 = TRUE,
    known_v_parameterization  = "auto",
    known_v_residual_fraction = NULL
  )
  covariance_samples <- array(c(0.01, 0.02), dim = c(2L, 1L, 1L))

  out <- .known_v_add_base_covariance(
    base_covariance    = known_V,
    covariance_samples = covariance_samples
  )

  expect_equal(dim(out), c(2L, 1L, 1L))
  expect_equal(out[, 1L, 1L], c(0.05, 0.06))
})


test_that("brma.mv supports hidden vi and sei diagonal known-V input", {

  vi <- c(0.04, 0.09, 0.16)

  vi_object <- brma.mv(
    yi                        = c(0.10, 0.20, -0.05),
    vi                        = vi,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_data                 = TRUE
  )
  sei_object <- brma.mv(
    yi                        = c(0.10, 0.20, -0.05),
    sei                       = sqrt(vi),
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_data                 = TRUE
  )
  both_object <- brma.mv(
    yi                        = c(0.10, 0.20, -0.05),
    vi                        = vi,
    sei                       = sqrt(vi),
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_data                 = TRUE
  )

  expect_equal(
    .known_v_materialize(.data_known_v_data(vi_object[["data"]])),
    diag(vi)
  )
  expect_equal(
    .known_v_materialize(.data_known_v_data(sei_object[["data"]])),
    diag(vi)
  )
  expect_equal(
    .known_v_materialize(.data_known_v_data(both_object[["data"]])),
    diag(vi)
  )
  expect_error(
    brma.mv(
      yi        = c(0.10, 0.20, -0.05),
      V         = vi,
      vi        = vi,
      measure   = "GEN",
      only_data = TRUE
    ),
    "Use only one of 'V' and hidden 'vi'/'sei'"
  )
  expect_error(
    brma.mv(
      yi        = c(0.10, 0.20, -0.05),
      vi        = vi,
      sei       = sqrt(vi) + 0.01,
      measure   = "GEN",
      only_data = TRUE
    ),
    "must be consistent"
  )
})


test_that("large diagonal known V retains linear storage", {

  K  <- 10000L
  vi <- seq(0.01, 0.02, length.out = K)

  object <- brma.mv(
    yi                        = rep(0, K),
    V                         = vi,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_data                 = TRUE
  )
  known_V <- .data_known_v_data(object[["data"]])

  expect_equal(.known_v_nrow(known_V), K)
  expect_equal(.known_v_diagonal(known_V), vi)
  expect_equal(known_V[["storage"]], "diagonal")
  expect_length(.known_v_correlated_blocks(known_V), 0L)
  expect_false(any(c("V", "B", "whitening_matrix", "sampling_factor") %in%
                   names(known_V)))
  expect_lt(as.numeric(utils::object.size(known_V)), 2 * 1024^2)
})


test_that("many small known-V blocks retain local backend artifacts", {

  n_blocks <- 1000L
  block    <- matrix(c(0.04, 0.01, 0.01, 0.09), nrow = 2L)
  V        <- rep(list(block), n_blocks)

  object <- brma.mv(
    yi                        = rep(0, 2L * n_blocks),
    V                         = V,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_data                 = TRUE
  )
  known_V <- .data_known_v_data(object[["data"]])

  expect_equal(.known_v_nrow(known_V), 2L * n_blocks)
  expect_equal(known_V[["storage"]], "blocks")
  expect_length(.known_v_correlated_blocks(known_V), n_blocks)
  expect_length(known_V[["whitening_blocks"]], n_blocks)
  expect_true(all(vapply(
    known_V[["whitening_blocks"]],
    function(x) identical(dim(x[["rotation"]]), c(2L, 2L)),
    logical(1)
  )))
  expect_false(any(c("V", "B", "whitening_matrix", "sampling_factor") %in%
                   names(known_V)))
  expect_lt(as.numeric(utils::object.size(known_V)), 10 * 1024^2)
})


test_that("known-V hashes are invariant to equivalent input representation", {

  yi       <- c(0.10, 0.20, -0.05)
  diagonal <- c(0.04, 0.09, 0.16)
  block    <- matrix(c(0.04, 0.01, 0.01, 0.09), nrow = 2L)

  make_object <- function(V) {
    brma.mv(
      yi                        = yi,
      V                         = V,
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_data                 = TRUE
    )
  }

  diagonal_objects <- list(
    make_object(diagonal),
    make_object(diag(diagonal)),
    make_object(lapply(diagonal, function(x) matrix(x, 1L, 1L)))
  )
  diagonal_hashes <- vapply(diagonal_objects, .get_outcome_hash, character(1))
  expect_length(unique(diagonal_hashes), 1L)

  block_dense <- rbind(
    cbind(block, matrix(0, 2L, 1L)),
    cbind(matrix(0, 1L, 2L), matrix(diagonal[[3L]], 1L, 1L))
  )
  block_objects <- list(
    make_object(block_dense),
    make_object(list(block, matrix(diagonal[[3L]], 1L, 1L)))
  )
  block_hashes <- vapply(block_objects, .get_outcome_hash, character(1))
  expect_length(unique(block_hashes), 1L)
})


test_that("internal mv data input defaults NULL known-V parameterization to auto", {

  yi <- c(0.10, 0.20)
  V  <- matrix(c(0.04, 0.01, 0.01, 0.09), nrow = 2)

  data <- .check_and_list_data(
    .call                              = quote(brma.mv(yi = yi, V = V)),
    .envir                             = environment(),
    class                              = "mv",
    measure                            = "GEN",
    set_contrast_factor_predictors    = "treatment",
    standardize_continuous_predictors = FALSE,
    random_group_covariance            = NULL,
    known_v_parameterization            = NULL,
    known_v_residual_fraction           = NULL,
    known_v_residual_fraction_specified = FALSE
  )
  known_V <- .data_known_v_data(data)

  expect_equal(known_V[["parameterization_requested"]], "auto")
  expect_equal(known_V[["parameterization"]], "whitened")
})


test_that("brma.mv warns and accepts rank-one all-correlated known V", {

  sei <- c(0.20, 0.30, 0.40)
  V   <- tcrossprod(sei)

  expect_warning(
    object <- brma.mv(
      yi                        = c(0.10, 0.20, 0.30),
      V                         = V,
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_data                 = TRUE
    ),
    "positive semidefinite"
  )

  known_V <- attr(object[["data"]], "known_V_data")

  expect_equal(.known_v_materialize(known_V), V)
  expect_equal(known_V[["parameterization"]], "whitened")
  expect_null(known_V[["sampling_factor"]])
  expect_lte(min(abs(known_V[["diagnostics"]][["min_whitening_variance"]])), 1e-10)

  expect_error(
    suppressWarnings(brma.mv(
      yi                        = c(0.10, 0.20, 0.30),
      V                         = V,
      known_v_parameterization  = "latent",
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_data                 = TRUE
    )),
    "cannot use known_v_parameterization = 'latent'"
  )

  dat <- data.frame(
    yi = c(0.10, 0.20, 0.30),
    x  = c(0, 1, 2)
  )
  old_max_block <- getOption("RoBMA.known_v_block_mvn_max_block_size", NULL)
  on.exit({
    options(RoBMA.known_v_block_mvn_max_block_size = old_max_block)
  }, add = TRUE)
  options(RoBMA.known_v_block_mvn_max_block_size = 2L)

  expect_warning(
    scale_object <- brma.mv(
      yi                        = yi,
      V                         = V,
      scale                     = ~ x,
      data                      = dat,
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_data                 = TRUE
    ),
    "positive semidefinite"
  )
  expect_equal(
    attr(scale_object[["data"]], "known_V_data")[["parameterization"]],
    "block_mvn"
  )
})


test_that("brma.mv aligns V after subset and missing rows", {

  V <- matrix(
    c(
      0.04, 0.01, 0.00, 0.00,
      0.01, 0.09, 0.02, 0.00,
      0.00, 0.02, 0.16, 0.03,
      0.00, 0.00, 0.03, 0.25
    ),
    nrow = 4,
    byrow = TRUE
  )
  dat <- data.frame(
    yi    = c(0.10, NA, 0.30, 0.40),
    x     = c(0, 1, 1, 0),
    study = c("a", "b", "c", "d"),
    label = paste0("Study ", 1:4)
  )

  expect_warning(
    object <- brma.mv(
      yi                        = yi,
      V                         = V,
      mods                      = ~ x,
      random                    = ~ 1 | study,
      data                      = dat,
      slab                      = label,
      subset                    = c(TRUE, TRUE, TRUE, TRUE),
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_data                 = TRUE
    ),
    "removed due to missing values"
  )

  expected_V <- V[c(1, 3, 4), c(1, 3, 4), drop = FALSE]
  expect_equal(
    .known_v_materialize(.data_known_v_data(object[["data"]])),
    expected_V
  )
  expect_equal(object[["data"]][["outcome"]][["yi"]], c(0.10, 0.30, 0.40))
  expect_equal(object[["data"]][["outcome"]][["slab"]],
               c("Study 1", "Study 3", "Study 4"))
  expect_true(attr(object[["data"]], "slab"))
  expect_equal(as.character(object[["data"]][["location"]][["study"]]), c("a", "c", "d"))

  expect_silent(
    subset_object <- brma.mv(
      yi                        = yi,
      V                         = V,
      mods                      = ~ x,
      random                    = ~ 1 | study,
      data                      = dat,
      slab                      = label,
      subset                    = c(TRUE, FALSE, TRUE, FALSE),
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_data                 = TRUE
    )
  )
  expect_equal(
    .known_v_materialize(.data_known_v_data(subset_object[["data"]])),
    V[c(1, 3), c(1, 3), drop = FALSE]
  )
  expect_equal(subset_object[["data"]][["outcome"]][["yi"]], c(0.10, 0.30))
  expect_equal(subset_object[["data"]][["outcome"]][["slab"]],
               c("Study 1", "Study 3"))
  expect_equal(as.character(subset_object[["data"]][["location"]][["study"]]),
               c("a", "c"))
})


test_that("brma.mv validates only retained known-V blocks", {

  valid_block <- matrix(
    c(0.04, 0.01, 0.01, 0.09),
    nrow = 2
  )
  indefinite_block <- matrix(
    c(0.04, 0.08, 0.08, 0.04),
    nrow = 2
  )
  nonfinite_block <- matrix(
    c(0.04, Inf, Inf, 0.09),
    nrow = 2
  )
  keep_valid <- c(TRUE, TRUE, FALSE, FALSE)
  keep_bad   <- !keep_valid
  yi         <- c(0.10, 0.20, 0.30, 0.40)

  expect_silent(
    indefinite_excluded <- brma.mv(
      yi                        = yi,
      V                         = list(valid_block, indefinite_block),
      subset                    = keep_valid,
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_data                 = TRUE
    )
  )
  expect_equal(
    .known_v_materialize(.data_known_v_data(indefinite_excluded[["data"]])),
    valid_block
  )

  expect_silent(
    nonfinite_excluded <- brma.mv(
      yi                        = yi,
      V                         = list(valid_block, nonfinite_block),
      subset                    = keep_valid,
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_data                 = TRUE
    )
  )
  expect_equal(
    .known_v_materialize(.data_known_v_data(nonfinite_excluded[["data"]])),
    valid_block
  )

  expect_error(
    brma.mv(
      yi                        = yi,
      V                         = list(valid_block, indefinite_block),
      subset                    = keep_bad,
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_data                 = TRUE
    ),
    "positive definite"
  )
  expect_error(
    brma.mv(
      yi                        = yi,
      V                         = list(valid_block, nonfinite_block),
      subset                    = keep_bad,
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_data                 = TRUE
    ),
    "finite non-missing"
  )
})


test_that("brma.mv drops diagonal-NA rows before known-V validation", {

  yi <- c(0.10, 0.20, 0.30)
  V  <- diag(c(0.04, NA_real_, 0.16))

  expect_warning(
    object <- brma.mv(
      yi                        = yi,
      V                         = V,
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_data                 = TRUE
    ),
    "removed due to missing values"
  )

  expect_equal(object[["data"]][["outcome"]][["yi"]], yi[c(1, 3)])
  expect_equal(
    .known_v_materialize(.data_known_v_data(object[["data"]])),
    V[c(1, 3), c(1, 3), drop = FALSE]
  )
})


test_that("brma.mv validates hidden vi and sei after row selection", {

  yi     <- c(0.10, 0.20, 0.30)
  subset <- c(TRUE, FALSE, TRUE)

  expect_warning(
    missing_vi <- brma.mv(
      yi                        = yi,
      vi                        = c(0.04, NA_real_, 0.16),
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_data                 = TRUE
    ),
    "removed due to missing values"
  )
  expect_equal(
    .known_v_materialize(.data_known_v_data(missing_vi[["data"]])),
    diag(c(0.04, 0.16))
  )

  expect_silent(
    excluded_signs <- brma.mv(
      yi                        = yi,
      vi                        = c(0.04, -0.09, 0.16),
      sei                       = c(0.20, -0.30, 0.40),
      subset                    = subset,
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_data                 = TRUE
    )
  )
  expect_equal(
    .known_v_materialize(.data_known_v_data(excluded_signs[["data"]])),
    diag(c(0.04, 0.16))
  )

  expect_silent(
    excluded_inconsistency <- brma.mv(
      yi                        = yi,
      vi                        = c(0.04, 0.09, 0.16),
      sei                       = c(0.20, 0.50, 0.40),
      subset                    = subset,
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_data                 = TRUE
    )
  )
  expect_equal(
    .known_v_materialize(.data_known_v_data(excluded_inconsistency[["data"]])),
    diag(c(0.04, 0.16))
  )

  expect_silent(
    excluded_nonfinite <- brma.mv(
      yi                        = yi,
      vi                        = c(0.04, Inf, 0.16),
      sei                       = c(0.20, Inf, 0.40),
      subset                    = subset,
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_data                 = TRUE
    )
  )
  expect_equal(
    .known_v_materialize(.data_known_v_data(excluded_nonfinite[["data"]])),
    diag(c(0.04, 0.16))
  )

  expect_error(
    brma.mv(
      yi                        = yi,
      vi                        = c(0.04, -0.09, 0.16),
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_data                 = TRUE
    ),
    "'vi' must contain positive finite values"
  )
  expect_error(
    brma.mv(
      yi                        = yi,
      sei                       = c(0.20, -0.30, 0.40),
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_data                 = TRUE
    ),
    "'sei' must contain positive finite values"
  )
  expect_error(
    brma.mv(
      yi                        = yi,
      vi                        = c(0.04, 0.09, 0.16),
      sei                       = c(0.20, 0.50, 0.40),
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_data                 = TRUE
    ),
    "must be consistent"
  )
})


test_that("brma.mv rejects invalid and unsupported inputs", {

  expect_error(
    brma.mv(
      yi                        = c(0.10, 0.20),
      V                         = matrix(c(1, 2, 2, 1), nrow = 2),
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_data                 = TRUE
    ),
    "positive definite"
  )

  sei <- c(0.20, 0.30, 0.40)
  R   <- matrix(
    c(
      1,  0,  0,
      0,  1, -1,
      0, -1,  1
    ),
    nrow  = 3,
    byrow = TRUE
  )
  V_singular_not_ones <- diag(sei) %*% R %*% diag(sei)
  expect_warning(
    singular_data <- brma.mv(
      yi                        = c(0.10, 0.20, 0.30),
      V                         = V_singular_not_ones,
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_data                 = TRUE
    ),
    "provisional"
  )
  expect_s3_class(singular_data, "brma.mv")
  expect_true(.data_known_v_data(singular_data[["data"]])[["singular"]])

  expect_error(
    brma.mv(
      yi                        = c(0.10, 0.20),
      V                         = diag(2),
      weights                   = c(1, 1),
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_data                 = TRUE
    ),
    "'weights' are not supported"
  )

  dat_unsupported <- data.frame(
    yi    = c(0.10, 0.20),
    study = c("s1", "s2")
  )
  expect_error(
    brma.mv(
      yi                        = yi,
      V                         = diag(2),
      cluster                   = study,
      data                      = dat_unsupported,
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_data                 = TRUE
    ),
    "'cluster' is not supported"
  )

  expect_error(
    brma.mv(
      yi                        = c(0.10, 0.20),
      V                         = diag(2),
      prior_bias                = TRUE,
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_data                 = TRUE
    ),
    "Selection/publication-bias priors are not supported"
  )

  expect_error(
    brma.mv(
      yi                         = c(0.10, 0.20),
      V                          = diag(2),
      known_v_parameterization   = "block_mvn",
      known_v_residual_fraction  = NA_real_,
      measure                    = "GEN",
      prior_unit_information_sd  = 1,
      only_data                  = TRUE
    ),
    "known_v_residual_fraction"
  )
})
