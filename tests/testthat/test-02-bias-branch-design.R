context("Bias-member design classes of cached ensembles")

# Product-space ensembles hold null, weight-function, PET, and PEESE members.
# Each posterior draw uses the design of its member: the null and
# weight-function members share the location-only design. The references are
# written in base R from the declared prior list, the posterior columns, and
# the union design columns.

source(testthat::test_path("common-functions.R"))

skip_if_no_fits()

fit_names <- list_fits()
fits      <- lazy_fits(fit_names, validate = FALSE)


# design kind (none, PET, PEESE) of every posterior draw
.bbc_draw_kinds <- function(object, rows = NULL) {

  samples <- .get_posterior_samples(object[["fit"]])
  if (!is.null(rows)) {
    samples <- samples[rows, , drop = FALSE]
  }
  members <- object[["priors"]][["outcome"]][["bias"]]
  kinds   <- vapply(members, function(prior) {
    if (BayesTools::is.prior.PET(prior)) {
      return("PET")
    }
    if (BayesTools::is.prior.PEESE(prior)) {
      return("PEESE")
    }
    return("none")
  }, character(1L))

  return(unname(kinds[samples[, "bias_indicator"]]))
}

# the location columns plus the bias column of the kind
.bbc_design <- function(union, kind) {

  bias <- c("PET", "PEESE")
  keep <- !(colnames(union) %in% bias) | colnames(union) == kind

  return(union[, keep, drop = FALSE])
}

# GVIF of every non-intercept column (all terms have one df)
.bbc_gvif <- function(vcov) {

  vcov <- vcov[-1L, -1L, drop = FALSE]
  if (ncol(vcov) == 1L) {
    return(1)
  }
  R <- stats::cov2cor(vcov)

  return(vapply(seq_len(ncol(R)), function(j) {
    det(R[-j, -j, drop = FALSE]) / det(R)
  }, numeric(1L)))
}

# heterogeneity SDs of the posterior draws, as the diagnostics evaluate them
.bbc_tau <- function(object) {

  K             <- nrow(.get_model_matrix(object))
  is_scale      <- .is_scale(object)
  is_multilevel <- .is_multilevel(object)

  .evaluate.brma.tau(
    fit               = object[["fit"]],
    scale_data        = object[["data"]][["scale"]],
    scale_formula     = if (is_scale) .create_fit_formula_list(
      data = object[["data"]],
      "scale"
    ) else NULL,
    scale_priors      = object[["priors"]][["scale"]],
    is_scale          = is_scale,
    is_multilevel     = is_multilevel,
    K                 = K,
    fixed_tau         = .fixed_tau_prior_value(object[["priors"]]),
    fixed_rho         = .fixed_rho_prior_value(object[["priors"]])
  )
}

# posterior mean of solve(X' M_s^-1 X) over the draws `rows`
.bbc_vcov <- function(object, X, tau, rows) {

  vi            <- .outcome_data_vi(object)
  weights       <- .outcome_data_weights(object)
  block_indices <- if (.is_multilevel(object)) {
    split(seq_len(nrow(X)), object[["data"]][["outcome"]][["cluster"]])
  } else {
    list()
  }

  Reduce(`+`, lapply(rows, function(s) {
    M <- diag((vi + tau[["tau_within"]][s, ]^2) / weights, nrow = length(vi))
    for (index in block_indices) {
      M[index, index] <- M[index, index] + tcrossprod(tau[["tau_between"]][s, index])
    }
    solve(crossprod(X, solve(M, X)))
  })) / length(rows)
}

.bbc_vif_reference <- function(object) {

  union <- .get_model_matrix(object)
  kinds <- .bbc_draw_kinds(object)
  tau   <- .bbc_tau(object)
  gvif  <- list()
  count <- c(none = 0, PET = 0, PEESE = 0)
  for (kind in names(count)) {
    rows        <- which(kinds == kind)
    count[kind] <- length(rows)
    X           <- .bbc_design(union, kind)
    gvif[[kind]] <- .bbc_gvif(.bbc_vcov(object, X, tau, rows))
  }

  list(
    gvif  = c(
      sum(count * vapply(gvif, `[`, numeric(1L), 1L)) / sum(count),
      gvif[["PET"]][2L],
      gvif[["PEESE"]][2L]
    ),
    count = count
  )
}


test_that("product-space VIF uses each member's design", {

  names <- c(
    "dat.lehmann2018_RoBMA_mods",
    "dat.lehmann2018_RoBMA_3lvl_mods_scale"
  )
  skip_if_missing_fits(names)

  for (name in names) {
    object    <- fits[[name]]
    reference <- .bbc_vif_reference(object)
    table     <- vif(object, posterior_correlation = FALSE)[["vif"]]

    expect_identical(table[["term"]], c("Preregistered", "PET", "PEESE"), info = name)
    expect_equal(table[["GVIF"]], reference[["gvif"]], tolerance = 1e-10, info = name)

    # the weight-function members belong to the location-only class
    members   <- object[["priors"]][["outcome"]][["bias"]]
    indicator <- .get_posterior_samples(object[["fit"]])[, "bias_indicator"]
    n_selection <- sum(vapply(members, .prior_has_selection, logical(1L))[indicator])
    n_null      <- sum(vapply(members, BayesTools::is.prior.none, logical(1L))[indicator])
    expect_gt(n_selection, 0)
    expect_equal(reference[["count"]][["none"]], n_selection + n_null, info = name)

    # the union design with both bias columns inflates their GVIF
    union <- .get_model_matrix(object)
    union_table <- .vif_table_from_vcov(
      vcov   = .vif_vcov_brma(object, union),
      assign = attr(union, "assign"),
      object = object,
      X      = union
    )
    expect_gt(union_table[["GVIF"]][2L], table[["GVIF"]][2L], label = name)
    expect_gt(union_table[["GVIF"]][3L], table[["GVIF"]][3L], label = name)
  }

  # the PSMA meta-regression of report 58: PET and PEESE are not collinear
  table <- vif(fits[["dat.lehmann2018_RoBMA_mods"]], posterior_correlation = FALSE)[["vif"]]
  expect_true(all(table[["GVIF"]][2:3] < 1.1))
})

test_that("product-space heterogeneity summaries use each member's design", {

  names <- c(
    "dat.lehmann2018_RoBMA_mods",
    "dat.lehmann2018_RoBMA_3lvl_mods_scale"
  )
  skip_if_missing_fits(names)

  for (name in names) {
    object <- fits[[name]]
    union  <- .get_model_matrix(object)
    kinds  <- .bbc_draw_kinds(object)
    tau    <- .bbc_tau(object)
    vi     <- .outcome_data_vi(object)
    K      <- length(vi)

    typical <- vapply(c("none", "PET", "PEESE"), function(kind) {
      X <- .bbc_design(union, kind)
      W <- diag(1 / vi)
      P <- W - W %*% X %*% solve(crossprod(X, W %*% X), t(X) %*% W)
      (K - ncol(X)) / sum(diag(P))
    }, numeric(1L))
    v_tilde <- unname(typical[kinds])

    within  <- tau[["tau_within"]]^2
    between <- if (is.null(tau[["tau_between"]])) 0 * within else tau[["tau_between"]]^2
    total   <- within + between
    expected <- c(
      tau2 = mean(rowMeans(total)),
      I2   = mean(rowMeans(100 * total / (total + v_tilde))),
      H2   = mean(rowMeans((total + v_tilde) / v_tilde))
    )
    if (.is_multilevel(object)) {
      expected[["I2 [within]"]]  <- mean(rowMeans(100 * within / (total + v_tilde)))
      expected[["I2 [between]"]] <- mean(rowMeans(100 * between / (total + v_tilde)))
    }

    summary <- as.data.frame(summary_heterogeneity(object))
    actual  <- summary[["Mean"]][match(names(expected), summary[["parameter"]])]

    expect_equal(actual, unname(expected), tolerance = 1e-10, info = name)
    expect_gt(length(unique(typical)), 1L)
  }
})


test_that("hat-based diagnostics still refuse ensembles with weight-function members", {

  name <- "dat.lehmann2018_RoBMA_mods"
  skip_if_missing_fits(name)

  expect_error(hatvalues(fits[[name]]), "not available for selection models")
  expect_error(rstandard(fits[[name]]), "not available for selection models")
})

# ---------------------------------------------------------------------------- #
# known-V product space: weight-function, null, PET, and PEESE members
# ---------------------------------------------------------------------------- #

test_that("known-V product-space hat values run on the design of each draw's member", {

  name <- "RoBMA.mv_marg_product_space"
  skip_if_missing_fits(name)

  object   <- fits[[name]]
  union    <- .get_model_matrix(object)
  samples  <- .get_posterior_samples(object[["fit"]])
  selected <- .thin_sample_rows(nrow(samples), 40)
  kinds    <- .bbc_draw_kinds(object, selected)
  expect_true(all(c("none", "PET", "PEESE") %in% kinds))

  covariance <- .known_v_marginal_covariance_samples(object, max_samples = 40)
  hat_draw <- function(X, M) {
    pivot <- qr(X)
    X     <- X[, pivot[["pivot"]][seq_len(pivot[["rank"]])], drop = FALSE]
    W     <- solve(M)
    diag(X %*% solve(crossprod(X, W %*% X), crossprod(X, W)))
  }
  expected <- t(vapply(seq_along(selected), function(i) {
    hat_draw(.bbc_design(union, kinds[i]), covariance[i, , ])
  }, numeric(nrow(union))))

  captured    <- NULL
  design_plan <- .hat_design_plan
  testthat::local_mocked_bindings(
    .hat_design_plan = function(...) {
      captured <<- design_plan(...)
      captured
    },
    .package = "RoBMA"
  )
  for (depth in c("marginal", "estimate")) {
    actual <- suppressWarnings(.compute_hat_matrix_samples(
      object, conditioning_depth = depth, max_samples = 40
    ))[["H_diag"]]

    expect_equal(actual, expected, tolerance = 1e-10, check.attributes = FALSE,
                 info = depth)
    # the classes are those of the selected rows, not of the first rows
    member_class <- .bias_design_member_classes(object[["priors"]])
    expect_identical(
      captured[["draw_class"]],
      unname(member_class[samples[selected, "bias_indicator"]]),
      info = depth
    )
    expect_false(identical(
      captured[["draw_class"]],
      unname(member_class[samples[seq_along(selected), "bias_indicator"]])
    ))
  }
})

test_that("known-V product-space VIF averages the draws of each design class", {

  name <- "RoBMA.mv_marg_product_space"
  skip_if_missing_fits(name)

  object      <- fits[[name]]
  union       <- .get_model_matrix(object)
  samples     <- .get_posterior_samples(object[["fit"]])
  selected    <- .thin_sample_rows(nrow(samples), 40)
  kinds       <- .bbc_draw_kinds(object, selected)
  covariance  <- .known_v_marginal_covariance_samples(object, max_samples = 40)

  # the location-only design has full rank: its covariance is unambiguous
  X        <- .bbc_design(union, "none")
  rows     <- which(kinds == "none")
  expected <- Reduce(`+`, lapply(rows, function(i) {
    solve(crossprod(X, solve(covariance[i, , ], X)))
  })) / length(rows)
  dimnames(expected) <- list(colnames(X), colnames(X))
  sample_info <- suppressWarnings(
    .vif_known_v_sample_info(object, max_samples = 40)
  )

  expect_equal(
    .vif_vcov_brma_mv(object, X, sample_info = sample_info, rows = rows),
    expected,
    tolerance = 1e-10,
    check.attributes = FALSE
  )

  table <- suppressWarnings(
    vif(object, posterior_correlation = FALSE, max_samples = 40)
  )[["vif"]]
  expect_identical(table[["term"]], c("x", "PET", "PEESE"))
  expect_true(all(is.finite(table[["GVIF"]])))
})
