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
    gvif[[kind]] <- if (length(rows)) {
      rowMeans(vapply(rows, function(s) .bbc_gvif(.bbc_vcov(object, X, tau, s)),
                      numeric(ncol(X) - 1L)))
    } else {
      rep(NA_real_, ncol(X) - 1L)
    }
  }

  list(
    gvif  = c(
      sum(count[count > 0L] * vapply(gvif[count > 0L], `[`, numeric(1L), 1L)) / sum(count),
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


# ---------------------------------------------------------------------------- #
# fitted ensembles (prepared by test-01-bias-branch-design.R)
# ---------------------------------------------------------------------------- #
#
# A PP ensemble with a moderator (null, PET and PEESE members, with LOO), the
# same without the moderator, a bPET fit, and a known-V RoBMA.mv() ensemble
# with PET and PEESE members, with the data they were fitted to. Each is read
# on first use; a test skips when the cache lacks the fit.

delayedAssign("bbd_fit", load_fit("bias_branch_PP_mods"))
delayedAssign("bbd_data", load_info("bias_branch_PP_mods")[["data"]])
delayedAssign("bbd_pet", load_fit("bias_branch_bPET_mods"))
delayedAssign("bbd_intercept_fit", load_fit("bias_branch_PP_intercept"))
delayedAssign("bbd_kv_fit", load_fit("bias_branch_known_v_ensemble"))
delayedAssign("bbd_kv_input", load_info("bias_branch_known_v_ensemble"))

# member kind (none, PET, PEESE) and heterogeneity SD of every posterior draw
.bbd_draws <- function(fit) {

  samples <- .get_posterior_samples(fit[["fit"]])
  bias    <- fit[["priors"]][["outcome"]][["bias"]]
  members <- if (BayesTools::is.prior.mixture(bias)) bias else list(bias)
  kinds   <- vapply(members, function(prior) {
    if (BayesTools::is.prior.PET(prior)) {
      return("PET")
    }
    if (BayesTools::is.prior.PEESE(prior)) {
      return("PEESE")
    }
    return("none")
  }, character(1L))
  indicator <- if ("bias_indicator" %in% colnames(samples)) {
    samples[, "bias_indicator"]
  } else {
    rep(1L, nrow(samples))
  }

  return(list(
    kind = unname(kinds[indicator]),
    tau  = unname(samples[, "tau"])
  ))
}

.bbd_design <- function(d, kind) {

  X <- cbind("(Intercept)" = 1, x = d[["x"]])
  if (kind == "PET") {
    X <- cbind(X, PET = d[["sei"]])
  } else if (kind == "PEESE") {
    X <- cbind(X, PEESE = d[["sei"]]^2)
  }

  return(X)
}

# weighted least squares pieces of one draw with marginal covariance M
.bbd_hat_draw <- function(X, M) {

  W <- solve(M)
  H <- X %*% solve(crossprod(X, W %*% X), crossprod(X, W))

  return(H)
}

# per-draw hat diagonals, a draws x K matrix
.bbd_hat <- function(d, draws, rows = seq_along(draws[["tau"]])) {

  vi <- d[["sei"]]^2
  t(vapply(rows, function(s) {
    X <- .bbd_design(d, draws[["kind"]][s])
    diag(.bbd_hat_draw(X, diag(vi + draws[["tau"]][s]^2)))
  }, numeric(nrow(d))))
}

# GVIF of every non-intercept column (all terms have one df)
.bbd_gvif <- function(vcov) {

  vcov <- vcov[-1L, -1L, drop = FALSE]
  p    <- ncol(vcov)
  if (p == 1L) {
    return(c(1))
  }
  R <- stats::cov2cor(vcov)

  return(vapply(seq_len(p), function(j) {
    det(R[-j, -j, drop = FALSE]) / det(R)
  }, numeric(1L)))
}

# reference VIF table of the ensemble: location term x, then PET, PEESE
.bbd_vif <- function(d, draws, rows = seq_along(draws[["tau"]]), adjusted = FALSE) {

  vi    <- d[["sei"]]^2
  gvif  <- list()
  count <- c(none = 0, PET = 0, PEESE = 0)
  for (kind in names(count)) {
    selected <- rows[draws[["kind"]][rows] == kind]
    count[kind] <- length(selected)
    if (length(selected) == 0L) {
      next
    }
    X <- .bbd_design(d, kind)
    values <- vapply(selected, function(s) {
      covariance <- solve(crossprod(X, X / (vi + draws[["tau"]][s]^2)))
      .bbd_gvif(covariance)
    }, numeric(ncol(X) - 1L))
    if (adjusted) values <- sqrt(values)
    gvif[[kind]] <- rowMeans(values)
  }
  present <- names(count)[count > 0]
  x_gvif  <- sum(vapply(present, function(kind) {
    count[[kind]] * gvif[[kind]][1L]
  }, numeric(1L))) / sum(count)

  return(c(
    x     = x_gvif,
    PET   = if (count[["PET"]] > 0) gvif[["PET"]][2L] else NA_real_,
    PEESE = if (count[["PEESE"]] > 0) gvif[["PEESE"]][2L] else NA_real_
  ))
}

# draw-wise GLS residuals, residual SDs and their ratio (marginal depth)
.bbd_residuals <- function(d, draws, V = NULL, rows = seq_along(draws[["tau"]])) {

  vi <- d[["sei"]]^2
  K  <- nrow(d)
  per_draw <- lapply(rows, function(s) {
    X <- .bbd_design(d, draws[["kind"]][s])
    M <- if (is.null(V)) diag(vi + draws[["tau"]][s]^2) else V + diag(draws[["tau"]][s]^2, K)
    A <- diag(K) - .bbd_hat_draw(X, M)
    residual <- as.vector(A %*% d[["yi"]])
    se       <- sqrt(diag(A %*% M %*% t(A)))
    list(resid = residual, se = se, z = residual / se)
  })

  return(data.frame(
    resid = rowMeans(vapply(per_draw, `[[`, numeric(K), "resid")),
    se    = rowMeans(vapply(per_draw, `[[`, numeric(K), "se")),
    z     = rowMeans(vapply(per_draw, `[[`, numeric(K), "z"))
  ))
}

# Run the hat computation with one design for every draw (the design of
# `class` from the union matrix), as a model with a single design would.
.bbd_single_design <- function(object, class, ...) {

  design <- .bias_class_design(.get_model_matrix(object), TRUE, TRUE, class)
  testthat::local_mocked_bindings(
    .get_model_matrix    = function(object) design,
    .bias_design_classes = function(...) NULL,
    .package             = "RoBMA"
  )

  return(.compute_hat_matrix_samples(object, ...))
}

.bbd_expect_class_rows <- function(object, kinds, ..., tolerance = 1e-10) {

  multi <- suppressWarnings(.compute_hat_matrix_samples(object, ...))
  for (class in 1:3) {
    rows   <- which(kinds == c("none", "PET", "PEESE")[class])
    single <- suppressWarnings(.bbd_single_design(object, class, ...))
    for (name in c("H_diag", "se", "resid")) {
      expect_equal(multi[[name]][rows, ], single[[name]][rows, ],
                   tolerance = tolerance, info = paste(name, "class", class))
    }
    expect_equal(multi[["H"]][rows, , ], single[["H"]][rows, , ],
                 tolerance = tolerance, info = paste("H class", class))
    expect_equal(multi[["M_diag"]][rows, ], single[["M_diag"]][rows, ],
                 tolerance = tolerance, info = paste("M_diag class", class))
  }
}

# per-draw hat diagonals of the known-V ensemble, marginal GLS
.bbd_kv_hat <- function(rows) {

  draws <- .bbd_draws(bbd_kv_fit)
  K     <- nrow(bbd_kv_input[["data"]])
  t(vapply(rows, function(s) {
    X <- .bbd_design(bbd_kv_input[["data"]], draws[["kind"]][s])
    M <- bbd_kv_input[["V"]] + diag(draws[["tau"]][s]^2, K)
    diag(.bbd_hat_draw(X, M))
  }, numeric(K)))
}

test_that("the ensemble has draws of every design class", {

  draws   <- .bbd_draws(bbd_fit)
  classes <- .bias_design_classes(bbd_fit)

  expect_true(all(c("none", "PET", "PEESE") %in% draws[["kind"]]))
  expect_identical(
    .bias_design_member_classes(bbd_fit[["priors"]]),
    c(1L, 2L, 3L)
  )
  expect_identical(
    classes[["class"]],
    unname(c(none = 1L, PET = 2L, PEESE = 3L)[draws[["kind"]]])
  )
  expect_identical(
    classes[["n_draws"]],
    as.integer(table(factor(draws[["kind"]], c("none", "PET", "PEESE"))))
  )
})

test_that("hat values average the hat matrix of each draw's member design", {

  draws    <- .bbd_draws(bbd_fit)
  expected <- .bbd_hat(bbd_data, draws)

  per_draw <- .compute_hat_matrix_samples(bbd_fit)[["H_diag"]]
  actual   <- hatvalues(bbd_fit)

  expect_equal(per_draw, expected, tolerance = 1e-10, check.attributes = FALSE)
  expect_equal(unname(actual), colMeans(expected), tolerance = 1e-10)

  # the trace is the mean number of design columns; the union design with both
  # bias columns would give four
  union <- vapply(seq_along(draws[["tau"]]), function(s) {
    X <- cbind(.bbd_design(bbd_data, "PET"), PEESE = bbd_data[["sei"]]^2)
    diag(.bbd_hat_draw(X, diag(bbd_data[["sei"]]^2 + draws[["tau"]][s]^2)))
  }, numeric(nrow(bbd_data)))
  expect_equal(sum(rowMeans(union)), 4, tolerance = 1e-10)
  expect_equal(
    sum(actual),
    mean(c(none = 2, PET = 3, PEESE = 3)[draws[["kind"]]]),
    tolerance = 1e-10
  )
  expect_lt(sum(actual), 3)

  # the full hat matrix of a draw uses the design of the draw's member
  full  <- .compute_hat_matrix_samples(bbd_fit, return_full_H = TRUE)[["H"]]
  first <- vapply(c("none", "PET", "PEESE"), function(kind) {
    which(draws[["kind"]] == kind)[1L]
  }, integer(1L))
  for (s in first) {
    X <- .bbd_design(bbd_data, draws[["kind"]][s])
    expect_equal(
      full[s, , ],
      .bbd_hat_draw(X, diag(bbd_data[["sei"]]^2 + draws[["tau"]][s]^2)),
      tolerance = 1e-10,
      info = draws[["kind"]][s]
    )
  }
})

test_that("a fitted single-bias model is one design", {

  expect_null(.bias_design_classes(bbd_pet))
})

test_that("a single bias term keeps its own hat values", {

  draws <- .bbd_draws(bbd_pet)
  draws[["kind"]] <- rep("PET", length(draws[["tau"]]))

  expect_equal(
    unname(hatvalues(bbd_pet)),
    colMeans(.bbd_hat(bbd_data, draws)),
    tolerance = 1e-10
  )
  expect_equal(sum(hatvalues(bbd_pet)), 3, tolerance = 1e-10)
})

test_that("internally standardized residuals use each draw's member design", {

  draws    <- .bbd_draws(bbd_fit)
  expected <- .bbd_residuals(bbd_data, draws)
  actual   <- rstandard(bbd_fit)

  expect_equal(unname(actual[["resid"]]), expected[["resid"]], tolerance = 1e-10)
  expect_equal(unname(actual[["se"]]), expected[["se"]], tolerance = 1e-10)
  expect_equal(unname(actual[["z"]]), expected[["z"]], tolerance = 1e-10)
})

test_that("the hat column of influence() uses each draw's member design", {

  draws <- .bbd_draws(bbd_fit)
  hat   <- suppressWarnings(influence(bbd_fit))[["inf"]][["hat"]]

  expect_equal(unname(hat), colMeans(.bbd_hat(bbd_data, draws)), tolerance = 1e-10)
})

test_that("VIF reports each member design's GVIF", {

  draws    <- .bbd_draws(bbd_fit)
  expected <- .bbd_vif(bbd_data, draws)
  table    <- vif(bbd_fit, posterior_correlation = FALSE)[["vif"]]

  expect_identical(table[["term"]], c("x", "PET", "PEESE"))
  expect_identical(table[["df"]], c(1L, 1L, 1L))
  expect_equal(table[["GVIF"]], unname(expected), tolerance = 1e-10)
  expect_equal(table[["GVIF^(1/(2*df))"]],
               unname(.bbd_vif(bbd_data, draws, adjusted = TRUE)), tolerance = 1e-10)
  expect_true(all(table[["GVIF"]] < 1.5))

  # the union design with both columns is nearly collinear
  union <- .vif_table_from_vcov(
    vcov   = .vif_vcov_brma(bbd_fit, .get_model_matrix(bbd_fit)),
    assign = attr(.get_model_matrix(bbd_fit), "assign"),
    object = bbd_fit,
    X      = .get_model_matrix(bbd_fit)
  )
  expect_gt(union[["GVIF"]][2L], 10)

  # a single bias term is one design
  single   <- vif(bbd_pet, posterior_correlation = FALSE)[["vif"]]
  draws    <- .bbd_draws(bbd_pet)
  draws[["kind"]] <- rep("PET", length(draws[["tau"]]))
  expected <- .bbd_vif(bbd_data, draws)
  expect_identical(single[["term"]], c("x", "PET"))
  expect_equal(single[["GVIF"]], unname(expected[c("x", "PET")]), tolerance = 1e-10)
})

test_that("VIF of a bias term without draws is unavailable", {

  samples   <- .get_posterior_samples(bbd_fit[["fit"]])
  indicator <- as.integer(samples[, "bias_indicator"])
  kinds     <- .bbd_draws(bbd_fit)[["kind"]]

  # no PEESE draws: the members of the PEESE class are relabelled as null
  # members; no null draws: the null members are relabelled as PET members
  scenarios <- list(
    list(
      relabel = c(1L, 2L, 1L), empty = "PEESE",
      kind    = ifelse(kinds == "PEESE", "none", kinds)
    ),
    list(
      relabel = c(2L, 2L, 3L), empty = "none",
      kind    = ifelse(kinds == "none", "PET", kinds)
    )
  )
  for (scenario in scenarios) {
    local({
      relabelled <- scenario[["relabel"]][indicator]
      testthat::local_mocked_bindings(
        .extract_bias_indicator = function(object, posterior_samples = NULL) {
          relabelled
        },
        .package = "RoBMA"
      )
      draws <- .bbd_draws(bbd_fit)
      draws[["kind"]] <- scenario[["kind"]]
      table <- vif(bbd_fit, posterior_correlation = FALSE)[["vif"]]
      expected <- .bbd_vif(bbd_data, draws)

      expect_identical(table[["term"]], c("x", "PET", "PEESE"))
      expect_equal(table[["GVIF"]], unname(expected), tolerance = 1e-10)
      if (scenario[["empty"]] == "PEESE") {
        expect_true(is.na(table[["GVIF"]][3L]))
        expect_true(is.na(table[["GVIF^(1/(2*df))"]][3L]))
      } else {
        expect_true(all(is.finite(table[["GVIF"]])))
      }
    })
  }
})

test_that("heterogeneity summaries use each member design's typical variance", {

  draws <- .bbd_draws(bbd_fit)
  vi    <- bbd_data[["sei"]]^2
  K     <- length(vi)
  typical <- function(X) {
    W <- diag(1 / vi)
    P <- W - W %*% X %*% solve(crossprod(X, W %*% X), t(X) %*% W)
    (K - ncol(X)) / sum(diag(P))
  }
  v_tilde <- vapply(c("none", "PET", "PEESE"), function(kind) {
    typical(.bbd_design(bbd_data, kind))
  }, numeric(1L))
  v <- unname(v_tilde[draws[["kind"]]])
  tau2 <- draws[["tau"]]^2

  summary <- as.data.frame(summary_heterogeneity(bbd_fit))
  mean_of <- function(parameter) summary[["Mean"]][summary[["parameter"]] == parameter]

  expect_equal(mean_of("I2"), mean(100 * tau2 / (tau2 + v)), tolerance = 1e-10)
  expect_equal(mean_of("H2"), mean((tau2 + v) / v), tolerance = 1e-10)
  expect_equal(mean_of("tau2"), mean(tau2), tolerance = 1e-10)

  # the union design gives a different typical variance
  union_v <- typical(cbind(.bbd_design(bbd_data, "PET"), PEESE = vi))
  expect_gt(abs(mean(100 * tau2 / (tau2 + union_v)) - mean_of("I2")), 1e-4)
})

test_that("every per-draw output equals the single-design output on its class rows", {

  kinds <- .bbd_draws(bbd_fit)[["kind"]]
  for (depth in c("marginal", "cluster", "estimate")) {
    .bbd_expect_class_rows(
      bbd_fit, kinds, conditioning_depth = depth, return_full_H = TRUE,
      return_se = TRUE, return_resid = TRUE
    )
  }
})

test_that("known-V hat values use each draw's member design", {

  draws <- .bbd_draws(bbd_kv_fit)
  S     <- length(draws[["tau"]])
  expect_true(all(c("none", "PET", "PEESE") %in% draws[["kind"]]))

  expected <- .bbd_kv_hat(seq_len(S))
  for (depth in c("marginal", "estimate")) {
    actual <- .compute_hat_matrix_samples(
      bbd_kv_fit, conditioning_depth = depth
    )[["H_diag"]]
    expect_equal(actual, expected, tolerance = 1e-8, check.attributes = FALSE,
                 info = depth)
  }
  expect_equal(unname(hatvalues(bbd_kv_fit)), colMeans(expected), tolerance = 1e-8)
})

test_that("known-V hat values of thinned draws use the designs of the selected draws", {

  draws    <- .bbd_draws(bbd_kv_fit)
  selected <- .thin_sample_rows(length(draws[["tau"]]), 50)
  expect_length(selected, 50L)

  actual <- expect_warning(
    .compute_hat_matrix_samples(bbd_kv_fit, max_samples = 50),
    "deterministically thinned"
  )[["H_diag"]]

  expect_equal(actual, .bbd_kv_hat(selected), tolerance = 1e-8,
               check.attributes = FALSE)
  expect_false(isTRUE(all.equal(
    actual, .bbd_kv_hat(1:50), tolerance = 1e-8, check.attributes = FALSE
  )))
})

test_that("known-V design classes follow the rows of the thinned posterior", {

  draws    <- .bbd_draws(bbd_kv_fit)
  selected <- .thin_sample_rows(length(draws[["tau"]]), 50)
  captured <- NULL
  design_plan <- .hat_design_plan
  testthat::local_mocked_bindings(
    .hat_design_plan = function(...) {
      captured <<- design_plan(...)
      captured
    },
    .package = "RoBMA"
  )

  suppressWarnings(.compute_hat_matrix_samples(bbd_kv_fit, max_samples = 50))

  classes <- c(none = 1L, PET = 2L, PEESE = 3L)
  expect_identical(
    captured[["draw_class"]],
    unname(classes[draws[["kind"]][selected]])
  )
  expect_false(identical(
    captured[["draw_class"]],
    unname(classes[draws[["kind"]][1:50]])
  ))
})

test_that("known-V standardized residuals use each draw's member design", {

  draws    <- .bbd_draws(bbd_kv_fit)
  expected <- .bbd_residuals(
    bbd_kv_input[["data"]], draws, V = bbd_kv_input[["V"]]
  )
  actual <- rstandard(bbd_kv_fit)

  expect_equal(unname(actual[["resid"]]), expected[["resid"]], tolerance = 1e-8)
  expect_equal(unname(actual[["se"]]), expected[["se"]], tolerance = 1e-8)
  expect_equal(unname(actual[["z"]]), expected[["z"]], tolerance = 1e-8)
})

test_that("known-V VIF reports each member design's GVIF", {

  d     <- bbd_kv_input[["data"]]
  V     <- bbd_kv_input[["V"]]
  draws <- .bbd_draws(bbd_kv_fit)
  K     <- nrow(d)
  reference <- function(rows) {
    count <- gvif <- list()
    for (kind in c("none", "PET", "PEESE")) {
      selected <- rows[draws[["kind"]][rows] == kind]
      count[[kind]] <- length(selected)
      X <- .bbd_design(d, kind)
      vcov <- Reduce(`+`, lapply(selected, function(s) {
        solve(crossprod(X, solve(V + diag(draws[["tau"]][s]^2, K), X)))
      })) / length(selected)
      gvif[[kind]] <- .bbd_gvif(vcov)
    }
    c(
      x     = sum(vapply(names(count), function(kind) {
        count[[kind]] * gvif[[kind]][1L]
      }, numeric(1L))) / sum(unlist(count)),
      PET   = gvif[["PET"]][2L],
      PEESE = gvif[["PEESE"]][2L]
    )
  }

  actual <- vif(bbd_kv_fit, posterior_correlation = FALSE)[["vif"]]
  expect_identical(actual[["term"]], c("x", "PET", "PEESE"))
  expect_equal(actual[["GVIF"]], unname(reference(seq_along(draws[["tau"]]))),
               tolerance = 1e-8)

  selected <- .thin_sample_rows(length(draws[["tau"]]), 60)
  thinned <- expect_warning(
    vif(bbd_kv_fit, posterior_correlation = FALSE, max_samples = 60),
    "deterministically thinned"
  )
  expect_equal(thinned[["vif"]][["GVIF"]], unname(reference(selected)),
               tolerance = 1e-8)
  expect_true(attr(thinned, "known_v_diagnostic")[["thinned"]])
})

test_that("every known-V per-draw output equals the single-design output on its class rows", {

  kinds <- .bbd_draws(bbd_kv_fit)[["kind"]]
  for (depth in c("marginal", "estimate")) {
    .bbd_expect_class_rows(
      bbd_kv_fit, kinds, conditioning_depth = depth, return_full_H = TRUE,
      return_se = TRUE, return_resid = TRUE, tolerance = 1e-8
    )
  }
})

test_that("intercept-only ensembles average the member designs' hat values", {

  fit <- bbd_intercept_fit

  draws <- .bbd_draws(fit)
  expect_true(all(c("none", "PET", "PEESE") %in% draws[["kind"]]))
  vi <- bbd_data[["sei"]]^2
  expected <- t(vapply(seq_along(draws[["tau"]]), function(s) {
    X <- cbind(rep(1, nrow(bbd_data)), switch(
      draws[["kind"]][s],
      PET   = bbd_data[["sei"]],
      PEESE = vi,
      NULL
    ))
    diag(.bbd_hat_draw(X, diag(vi + draws[["tau"]][s]^2)))
  }, numeric(nrow(bbd_data))))

  expect_equal(unname(hatvalues(fit)), colMeans(expected), tolerance = 1e-10)
  expect_equal(
    sum(hatvalues(fit)),
    mean(c(none = 1, PET = 2, PEESE = 2)[draws[["kind"]]]),
    tolerance = 1e-10
  )
})
