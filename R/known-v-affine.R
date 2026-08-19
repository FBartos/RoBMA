# ============================================================================ #
# Exact affine covariance families
# ============================================================================ #

.known_v_affine_spectral_plan <- function(base_covariance,
                                          update_covariance,
                                          reference_coefficient = 0) {

  if (!is.matrix(base_covariance) || !is.numeric(base_covariance) ||
      nrow(base_covariance) == 0L ||
      nrow(base_covariance) != ncol(base_covariance) ||
      !is.matrix(update_covariance) || !is.numeric(update_covariance) ||
      !identical(dim(update_covariance), dim(base_covariance)) ||
      any(!is.finite(base_covariance)) ||
      any(!is.finite(update_covariance)) ||
      !is.numeric(reference_coefficient) ||
      length(reference_coefficient) != 1L ||
      !is.finite(reference_coefficient)) {
    stop("Internal error: invalid affine covariance family.", call. = FALSE)
  }
  if (any(base_covariance != t(base_covariance)) ||
      any(update_covariance != t(update_covariance))) {
    stop(
      "Internal error: affine covariance matrices must be symmetric (maximum asymmetry: ",
      format(max(
        abs(base_covariance - t(base_covariance)),
        abs(update_covariance - t(update_covariance))
      ), scientific = TRUE),
      ").",
      call. = FALSE
    )
  }

  root <- tryCatch(chol(base_covariance), error = function(e) NULL)
  if (is.null(root)) {
    return(NULL)
  }
  transformed <- forwardsolve(t(root), update_covariance)
  transformed <- t(forwardsolve(t(root), t(transformed)))
  # The source is structurally symmetric; avoid separately rounded triangles.
  transformed[lower.tri(transformed)] <-
    t(transformed)[lower.tri(transformed)]
  decomposition <- tryCatch(
    eigen(transformed, symmetric = TRUE),
    error = function(e) NULL
  )
  if (is.null(decomposition) ||
      any(!is.finite(decomposition[["values"]])) ||
      any(!is.finite(decomposition[["vectors"]]))) {
    return(NULL)
  }

  list(
    root                  = root,
    eigenvalues           = unname(decomposition[["values"]]),
    eigenvectors          = unname(decomposition[["vectors"]]),
    log_determinant_base  = 2 * sum(log(diag(root))),
    reference_coefficient = reference_coefficient,
    base_diagonal         = diag(base_covariance),
    update_diagonal       = diag(update_covariance),
    dimension             = nrow(base_covariance)
  )
}


.known_v_affine_spectral_denominators <- function(plan, coefficients) {

  delta <- coefficients - plan[["reference_coefficient"]]
  1 + outer(delta, plan[["eigenvalues"]])
}


.known_v_affine_log_likelihood <- function(plan, residual, coefficients) {

  K <- plan[["dimension"]]
  if (!is.numeric(residual) || length(residual) != K ||
      any(!is.finite(residual)) || !is.numeric(coefficients) ||
      any(!is.finite(coefficients))) {
    stop("Internal error: invalid affine Gaussian likelihood inputs.",
         call. = FALSE)
  }

  denominators <- .known_v_affine_spectral_denominators(
    plan,
    coefficients
  )
  if (any(!is.finite(denominators)) || any(denominators <= 0)) {
    return(NULL)
  }
  whitened <- forwardsolve(t(plan[["root"]]), residual)
  spectral <- as.vector(crossprod(plan[["eigenvectors"]], whitened))
  log_determinant <- plan[["log_determinant_base"]] +
    rowSums(log(denominators))
  quadratic <- as.vector((1 / denominators) %*% spectral^2)

  -0.5 * (K * log(2 * pi) + log_determinant + quadratic)
}


.known_v_affine_gls_projection_setup <- function(plan, X, y) {

  K <- plan[["dimension"]]
  if (!is.matrix(X) || !is.numeric(X) || nrow(X) != K ||
      any(!is.finite(X)) || !is.numeric(y) || length(y) != K ||
      any(!is.finite(y))) {
    stop("Internal error: invalid affine GLS inputs.", call. = FALSE)
  }
  rhs <- cbind(y, X)
  whitened <- forwardsolve(t(plan[["root"]]), rhs)

  list(
    X              = X,
    y              = y,
    spectral_rhs   = crossprod(plan[["eigenvectors"]], whitened),
    solve_transform = backsolve(plan[["root"]], plan[["eigenvectors"]])
  )
}


.known_v_affine_gls_projection <- function(plan, X, y, coefficient,
                                           return_full_H = FALSE,
                                           setup = NULL) {

  K <- plan[["dimension"]]
  if (!is.matrix(X) || !is.numeric(X) || nrow(X) != K ||
      any(!is.finite(X)) || !is.numeric(y) || length(y) != K ||
      any(!is.finite(y)) || !is.numeric(coefficient) ||
      length(coefficient) != 1L || !is.finite(coefficient) ||
      !is.logical(return_full_H) || length(return_full_H) != 1L ||
      is.na(return_full_H)) {
    stop("Internal error: invalid affine GLS inputs.", call. = FALSE)
  }
  if (is.null(setup)) {
    setup <- .known_v_affine_gls_projection_setup(plan, X, y)
  } else if (!is.list(setup) || !identical(setup[["X"]], X) ||
             !identical(setup[["y"]], y)) {
    stop("Internal error: affine GLS setup does not match its inputs.",
         call. = FALSE)
  }

  denominators <- .known_v_affine_spectral_denominators(plan, coefficient)
  if (any(!is.finite(denominators)) || any(denominators <= 0)) {
    return(NULL)
  }
  weighted <- setup[["spectral_rhs"]] / as.vector(denominators)
  solved   <- setup[["solve_transform"]] %*% weighted
  Wy       <- solved[, 1L]
  WX       <- solved[, -1L, drop = FALSE]

  XtWX     <- crossprod(X, WX)
  XtWX_inv <- .hat_solve_crossprod(XtWX)
  beta_hat <- as.vector(XtWX_inv %*% crossprod(X, Wy))
  residual <- y - as.vector(X %*% beta_hat)
  XC       <- X %*% XtWX_inv
  H_diag   <- rowSums(XC * WX)
  coefficient_delta <- coefficient - plan[["reference_coefficient"]]
  covariance_diagonal <- plan[["base_diagonal"]] +
    coefficient_delta * plan[["update_diagonal"]]
  residual_variance <-
    covariance_diagonal - rowSums(XC * X)

  list(
    WX                = WX,
    XtWX_inv          = XtWX_inv,
    H                  = if (return_full_H) XC %*% t(WX) else NULL,
    H_diag             = H_diag,
    beta_hat           = beta_hat,
    residual           = residual,
    residual_variance  = residual_variance,
    covariance_diagonal = covariance_diagonal
  )
}


.known_v_random_affine_diagnostic_plan <- function(object, posterior_samples,
                                                    known_V) {

  if (!inherits(object, "brma.mv") || !is.matrix(posterior_samples) ||
      nrow(posterior_samples) == 0L) {
    return(NULL)
  }
  result <- tryCatch({
    catalog <- BayesTools::parameter_catalog(object[["fit"]])
    quantities <- catalog[["quantities"]]
    candidate_rows <- which(
      quantities[["formula_parameter"]] == "mu" &
        quantities[["role"]] %in% c("random_sd", "random_sd_ratio")
    )
    updates <- lapply(candidate_rows, function(index) {
      selection <- BayesTools::parameter_catalog_resolve(
        catalog,
        quantities[["canonical_name"]][[index]],
        "mu"
      )
      BayesTools::random_effects_marginal_update_plan(
        object[["fit"]],
        selection
      )
    })
    eligible <- vapply(updates, function(update) {
      source <- update[["source_parameter"]]
      identical(update[["family"]], "affine") &&
        identical(update[["update"]], "scale") &&
        identical(update[["coefficient_input"]], "source") &&
        !is.null(update[["invariant_covariance"]]) &&
        !is.null(update[["coefficient_transform"]]) &&
        is.character(source) && length(source) == 1L &&
        !is.na(source) && nzchar(source) &&
        source %in% colnames(posterior_samples)
    }, logical(1))
    updates <- updates[eligible]
    sources <- vapply(updates, `[[`, character(1), "source_parameter")
    if (length(updates) == 0L || length(unique(sources)) != 1L) {
      return(NULL)
    }
    if (length(updates) > 1L && !all(vapply(
      updates[-1L],
      function(update) identical(
        update[["invariant_covariance"]],
        updates[[1L]][["invariant_covariance"]]
      ),
      logical(1)
    ))) {
      return(NULL)
    }
    update <- updates[[1L]]
    invariant <- update[["invariant_covariance"]]
    source <- update[["source_parameter"]]
    source_values <- posterior_samples[, source]
    if (any(!is.finite(source_values)) || any(source_values < 0)) {
      return(NULL)
    }
    coefficients <- BayesTools::parameter_transform_forward(
      source_values,
      update[["coefficient_transform"]]
    )
    if (any(!is.finite(coefficients))) {
      return(NULL)
    }
    reference <- coefficients[[1L]]
    sampling_covariance <- .marglik_known_v_covariance_matrix(known_V)
    base_covariance <- sampling_covariance +
      invariant[["base_covariance"]] +
      reference * invariant[["update_covariance"]]
    plan <- .known_v_affine_spectral_plan(
      base_covariance       = base_covariance,
      update_covariance     = invariant[["update_covariance"]],
      reference_coefficient = reference
    )
    if (is.null(plan)) {
      return(NULL)
    }

    list(
      plan         = plan,
      coefficients = coefficients,
      source       = source
    )
  }, error = function(error) NULL)

  result
}
