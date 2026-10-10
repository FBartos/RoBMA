context("Bias-member design classes")

# A RoBMA() ensemble mixes members with different bias terms: no member holds
# both the PET and the PEESE column. The design-based diagnostics use, for each
# posterior draw, the design of the draw's member. The references below are
# written in base R from the declared prior list, the posterior columns, and
# the data, not from the package's design helpers.

skip_on_cran()


.bbd_data <- function() {

  withr::local_seed(2)
  k <- 20
  d <- data.frame(sei = stats::runif(k, 0.08, 0.30), x = stats::rnorm(k))
  d[["yi"]] <- 0.25 + 0.1 * d[["x"]] + 0.5 * d[["sei"]] +
    stats::rnorm(k, 0, sqrt(0.02^2 + d[["sei"]]^2))
  d[["ni"]] <- 50

  return(d)
}

.bbd_fit_settings <- function() {
  list(
    chains = 1, burnin = 500, adapt = 300, sample = 1000, seed = 1,
    silent = TRUE, parallel = FALSE,
    convergence_checks = set_convergence_checks(max_Rhat = NULL, min_ESS = NULL)
  )
}

bbd_data <- .bbd_data()
bbd_args <- c(
  list(
    yi = quote(yi), sei = quote(sei), ni = quote(ni), mods = ~ x,
    data = bbd_data, measure = "GEN", standardize_continuous_predictors = FALSE
  ),
  .bbd_fit_settings()
)
bbd_fit <- suppressWarnings(do.call(RoBMA, c(bbd_args, list(model_type = "PP"))))
bbd_fit <- suppressWarnings(add_loo(bbd_fit))
bbd_pet <- suppressWarnings(do.call(bPET, bbd_args))


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
.bbd_vif <- function(d, draws, rows = seq_along(draws[["tau"]])) {

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
    vcov <- Reduce(`+`, lapply(selected, function(s) {
      solve(crossprod(X, X / (vi + draws[["tau"]][s]^2)))
    })) / length(selected)
    gvif[[kind]] <- .bbd_gvif(vcov)
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


.bbd_only_priors <- function(mods = TRUE) {

  keep <- c("yi", "sei", "ni", if (mods) "mods", "data", "measure",
            "standardize_continuous_predictors")

  return(c(bbd_args[keep], list(only_priors = TRUE, silent = TRUE)))
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

test_that("class designs equal the designs of single-bias models", {

  for (mods in c(TRUE, FALSE)) {
    args     <- .bbd_only_priors(mods)
    ensemble <- suppressWarnings(do.call(RoBMA, c(args, list(model_type = "PP"))))
    single   <- list(
      suppressWarnings(do.call(brma, args)),
      suppressWarnings(do.call(bPET, args)),
      suppressWarnings(do.call(bPEESE, args))
    )

    union <- .get_model_matrix(ensemble)
    expect_identical(
      colnames(union),
      c("(Intercept)", if (mods) "x", "PET", "PEESE")
    )
    for (class in 1:3) {
      expect_equal(
        .bias_class_design(union, TRUE, TRUE, class),
        .get_model_matrix(single[[class]]),
        info = paste("moderators:", mods, "class:", class)
      )
    }
  }

  # an ensemble of the null and PET members holds one bias column
  args    <- .bbd_only_priors()
  pet_null <- suppressWarnings(do.call(RoBMA, c(
    args,
    list(prior_bias = BayesTools::prior_PET("cauchy", list(0, 1), list(0, Inf)))
  )))
  union <- .get_model_matrix(pet_null)
  expect_identical(colnames(union), c("(Intercept)", "x", "PET"))
  expect_identical(.bias_design_member_classes(pet_null[["priors"]]), c(1L, 2L))
  expect_equal(
    .bias_class_design(union, TRUE, FALSE, 2L),
    .get_model_matrix(suppressWarnings(do.call(bPET, args)))
  )
  expect_equal(
    .bias_class_design(union, TRUE, FALSE, 1L),
    .get_model_matrix(suppressWarnings(do.call(brma, args)))
  )
})

test_that("models with one design use the union design", {

  args <- .bbd_only_priors()
  for (constructor in list(brma, bPET, bPEESE)) {
    expect_null(.bias_design_classes(suppressWarnings(do.call(constructor, args))))
  }
  expect_null(.bias_design_classes(bbd_pet))

  # members of one class form one design
  members <- list(
    BayesTools::prior_PET("cauchy", list(0, 1), list(0, Inf)),
    BayesTools::prior_PET("normal", list(0, 1), list(0, Inf))
  )
  two_pet <- suppressWarnings(do.call(RoBMA, c(
    args, list(prior_bias = members, prior_bias_null = NULL)
  )))
  expect_identical(.bias_design_member_classes(two_pet[["priors"]]), c(2L, 2L))
  expect_null(.bias_design_classes(two_pet))

  # the weight-function members share the location-only design with the null
  selection <- suppressWarnings(do.call(bselmodel, c(
    args[c("yi", "sei", "ni", "data", "measure", "only_priors", "silent")]
  )))
  expect_true(all(.bias_design_member_classes(selection[["priors"]]) == 1L))
  expect_null(.bias_design_classes(selection))

  expect_null(.bias_design_classes(
    structure(list(priors = list(outcome = list(bias = NULL))), class = "brma")
  ))
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
  expect_equal(table[["GVIF^(1/(2*df))"]], unname(sqrt(expected)), tolerance = 1e-10)
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

test_that("the ensemble table weights class GVIFs by their draw counts", {

  tables <- list(
    data.frame(term = "a", df = 1L, GVIF = 1, "GVIF^(1/(2*df))" = 1, check.names = FALSE),
    data.frame(term = c("a", "PET"), df = c(1L, 1L), GVIF = c(2, 3),
               "GVIF^(1/(2*df))" = sqrt(c(2, 3)), check.names = FALSE),
    NULL
  )
  designs <- list(matrix(0), matrix(0), matrix(0))

  full <- .vif_ensemble_table(tables, n_draws = c(2L, 6L, 0L), designs = designs)

  expect_identical(full[["term"]], c("a", "PET", "PEESE"))
  expect_equal(full[["GVIF"]], c(0.25 * 1 + 0.75 * 2, 3, NA_real_))
  expect_equal(full[["GVIF^(1/(2*df))"]], sqrt(full[["GVIF"]]))
  expect_identical(names(full), c("term", "df", "GVIF", "GVIF^(1/(2*df))"))
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

test_that("every per-draw output equals the single-design output on its class rows", {

  kinds <- .bbd_draws(bbd_fit)[["kind"]]
  for (depth in c("marginal", "cluster", "estimate")) {
    .bbd_expect_class_rows(
      bbd_fit, kinds, conditioning_depth = depth, return_full_H = TRUE,
      return_se = TRUE, return_resid = TRUE
    )
  }
})

test_that("exactly fitted rows are zeroed per design class", {

  design_plan <- list(
    draw_class         = c(1L, 2L, 1L, 2L),
    zero_residual_rows = list(2L, c(1L, 3L), integer())
  )
  x <- matrix(1, nrow = 4, ncol = 3)

  zeroed <- .hat_zero_class_rows(x, design_plan)

  expected <- x
  expected[c(1L, 3L), 2L] <- 0
  expected[c(2L, 4L), c(1L, 3L)] <- 0
  expect_identical(zeroed, expected)
})


# ---------------------------------------------------------------------------- #
# known-V ensemble: marginal GLS with the member design of each draw
# ---------------------------------------------------------------------------- #

.bbd_known_v_input <- function() {

  withr::local_seed(3)
  k_study <- 6L
  sei     <- stats::runif(2L * k_study, 0.08, 0.30)
  x       <- stats::rnorm(2L * k_study)
  d <- data.frame(
    sei = sei, x = x,
    yi  = 0.2 + 0.1 * x + 0.5 * sei + stats::rnorm(2L * k_study, 0, sei)
  )
  V <- matrix(0, 2L * k_study, 2L * k_study)
  for (s in seq_len(k_study)) {
    index <- 2L * s - 1:0
    V[index, index] <- 0.4 * tcrossprod(sei[index])
    diag(V[index, index]) <- sei[index]^2
  }

  return(list(data = d, V = V))
}

bbd_kv_input <- .bbd_known_v_input()
bbd_kv_fit <- suppressWarnings(RoBMA.mv(
  yi = yi, V = bbd_kv_input[["V"]], mods = ~ x, data = bbd_kv_input[["data"]],
  measure = "GEN",
  prior_bias = list(
    BayesTools::prior_PET("cauchy", list(0, 1), list(0, Inf)),
    BayesTools::prior_PEESE("cauchy", list(0, 1), list(0, Inf))
  ),
  prior_unit_information_sd = 1, standardize_continuous_predictors = FALSE,
  chains = 1, sample = 400, burnin = 200, adapt = 200, seed = 1, silent = TRUE,
  parallel = FALSE,
  convergence_checks = set_convergence_checks(max_Rhat = NULL, min_ESS = NULL)
))

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

  args <- c(
    bbd_args[c("yi", "sei", "ni", "data", "measure")],
    .bbd_fit_settings()
  )
  fit <- suppressWarnings(do.call(RoBMA, c(args, list(model_type = "PP"))))

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


# ---------------------------------------------------------------------------- #
# known-V design classes: failures name the draws of the request
# ---------------------------------------------------------------------------- #

# A replay of a known-V ensemble (null, PET and PEESE members, a random
# intercept) with 19 draws whose design classes interleave. The known
# covariance is ill-conditioned only in the draw `ill_row`, where the random
# SD is 0 and V alone holds a pair of estimates with correlation 1 - 1e-12.
# The odd rows are the draws that 'max_samples = 10' keeps, with the classes
# 1, 2, 1, 2, 3, 1, 2, 3, 1, 2.
.bbd_ill_conditioned_replay <- function(ill_row) {

  dat <- data.frame(
    yi  = c(.1, .3, .2, .4, .15, .2),
    x   = c(-1, 0, 1, 2, .3, -.4),
    sei = c(.2, .3, .4, .25, .15, .22),
    study = rep(1:3, each = 2)
  )
  V <- diag(dat[["sei"]]^2)
  V[1L, 2L] <- V[2L, 1L] <- (1 - 1e-12) * dat[["sei"]][1L] * dat[["sei"]][2L]
  object <- RoBMA.mv(
    yi = yi, V = V, mods = ~ x, random = ~ 1 | study, data = dat,
    measure = "GEN", prior_unit_information_sd = 1, only_priors = TRUE,
    known_v_parameterization = "block_mvn",
    prior_heterogeneity = BayesTools::prior("normal", list(0, 1), list(0, Inf)),
    prior_heterogeneity_null = NULL, standardize_continuous_predictors = FALSE,
    silent = TRUE,
    prior_bias = list(
      BayesTools::prior_PET("normal", list(0, 1), list(0, Inf)),
      BayesTools::prior_PEESE("normal", list(0, 1), list(0, Inf))
    )
  )
  setup  <- review_compile_object(object)
  priors <- c(
    setup[["compiled"]][["mu"]][["prior_list"]],
    .create_fit_priors(object[["data"]], object[["priors"]])
  )
  gates   <- BayesTools:::.bt_random_variance_allocation_inclusion_indicator_names(
    setup[["object"]][["formula_design"]]
  )
  columns <- unique(c(
    unlist(lapply(names(priors), function(name) {
      BayesTools:::.prior_linear_prior_columns(name, priors[[name]])
    })),
    "bias_indicator", "PET", "PEESE", "tau", "tau_indicator", gates
  ))

  n_draws <- 19L
  samples <- matrix(.1, n_draws, length(columns), dimnames = list(NULL, columns))
  for (name in intersect(names(priors), columns)) {
    location <- priors[[name]][["parameters"]][["location"]]
    if (!is.null(location)) {
      samples[, name] <- location
    }
  }
  samples[, grepl("indicator", columns)] <- 2
  if (length(gates)) {
    samples[, gates] <- 1
  }
  samples[, "bias_indicator"] <- rep(c(3, 3, 2, 1, 1, 2, 3, 1, 2), length.out = n_draws)
  samples[seq(1L, n_draws, by = 2L), "bias_indicator"] <- c(1, 2, 1, 2, 3, 1, 2, 3, 1, 2)
  tau <- rep(.2, n_draws)
  tau[ill_row] <- 0
  samples[, "tau"] <- tau
  samples[, "tau_indicator"] <- ifelse(tau > 0, 2, 1)
  for (column in columns[grepl("xREx.*intercept$|__allocation_sd$", columns)]) {
    samples[, column] <- tau
  }

  replay <- review_attach_replay(setup, samples)
  class(replay) <- setdiff(class(replay), "only_priors.brma")

  return(replay)
}

# the hat computation without the monitored random coefficients that the
# likelihood setup of a replay lacks (they matter at the estimate depth only)
.bbd_replay_hat_setup <- function(object, bias_adjusted = FALSE,
                                  posterior_samples = NULL,
                                  condition_local_effects = TRUE) {

  samples <- .get_posterior_samples(object[["fit"]], posterior_samples)

  return(list(
    posterior_samples = samples,
    data              = object[["data"]],
    yi                = .outcome_data_yi(object),
    S                 = nrow(samples),
    K                 = length(.outcome_data_yi(object))
  ))
}

test_that("known-V design classes report precision failures at the draws of the request", {

  skip_on_cran()
  unavailable <- function(expr) {
    tryCatch(expr, RoBMA_known_v_precision_unavailable = identity)
  }
  # request position of each ill-conditioned row, with and without thinning,
  # and the class of the draw
  cases <- list(
    list(row = 11L, thinned = 6L, class = 1L),
    list(row = 13L, thinned = 7L, class = 2L),
    list(row = 15L, thinned = 8L, class = 3L)
  )
  for (case in cases) {
    replay  <- .bbd_ill_conditioned_replay(case[["row"]])
    samples <- .get_posterior_samples(replay[["fit"]])
    expect_identical(
      .bias_design_classes(replay)[["class"]][case[["row"]]], case[["class"]]
    )

    # the draw of the class that fails, as the class batch itself reports it
    class_rows <- which(.bias_design_classes(replay)[["class"]] == case[["class"]])
    design     <- .bias_design_classes(replay)[["designs"]][[case[["class"]]]]
    local      <- unavailable(.vif_vcov_from_factor_plan(
      replay, design, samples[class_rows, , drop = FALSE],
      .data_known_v_data(replay[["data"]])
    ))
    expect_s3_class(local, "RoBMA_known_v_precision_unavailable")
    expect_identical(local[["rows"]], match(case[["row"]], class_rows))
    expect_lt(local[["reciprocal_condition"]], 1e-9)

    testthat::local_mocked_bindings(
      .estimate_likelihood_setup.brma = .bbd_replay_hat_setup,
      .package = "RoBMA"
    )
    for (request in c("all", "thinned")) {
      expected <- if (request == "all") case[["row"]] else case[["thinned"]]
      max_samples <- if (request == "all") Inf else 10
      calls <- list(
        vif = function() vif(
          replay, posterior_correlation = FALSE, max_samples = max_samples
        ),
        hat = function() .compute_hat_matrix_samples(
          replay, conditioning_depth = "marginal", max_samples = max_samples
        ),
        resid = function() .compute_hat_matrix_samples(
          replay, conditioning_depth = "marginal", max_samples = max_samples,
          return_se = TRUE, return_resid = TRUE
        )
      )
      for (name in names(calls)) {
        condition <- unavailable(suppressWarnings(calls[[name]]()))
        info <- paste(name, "row", case[["row"]], request)

        expect_s3_class(condition, "RoBMA_known_v_precision_unavailable", exact = FALSE)
        expect_identical(
          class(condition),
          c("RoBMA_known_v_precision_unavailable", "RoBMA_numerical_unavailable",
            "error", "condition")
        )
        expect_identical(condition[["rows"]], expected, info = info)
        expect_identical(
          condition[["reciprocal_condition"]], local[["reciprocal_condition"]],
          info = info
        )
        expect_match(
          conditionMessage(condition),
          paste0("at draw ", expected, " among the posterior draws it was applied to"),
          fixed = TRUE, info = info
        )
      }
    }
  }
})
