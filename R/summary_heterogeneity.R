# ============================================================================ #
# brma.summary_heterogeneity.R
# ============================================================================ #
#
# This file contains the summary_heterogeneity method for brma model fits.
# It computes absolute and relative heterogeneity measures:
# - tau, tau^2: Absolute heterogeneity (SD and variance)
# - I^2: Percentage of total variance due to heterogeneity
# - H^2: Ratio of total variance to sampling variance
#
# For multilevel (3-level) models, heterogeneity is partitioned into:
# - tau [within]: Estimate-level heterogeneity
# - tau [between]: Cluster-level heterogeneity
# - I^2 [within], I^2 [between]: Partitioned relative heterogeneity
#
# Formulas follow metafor documentation and Higgins & Thompson (2002).
#
# ============================================================================ #


# ---------------------------------------------------------------------------- #
# summary_heterogeneity generic
# ---------------------------------------------------------------------------- #

#' @title Summary of Heterogeneity
#'
#' @description Computes method-specific absolute and relative heterogeneity
#' summaries for a fitted model.
#'
#' @param object a fitted model object
#' @param ... additional arguments passed to methods
#'
#' @return Method-specific return value containing heterogeneity estimates.
#'
#' @seealso [pooled_heterogeneity()], [summary.brma()]
#' @export
summary_heterogeneity <- function(object, ...) {

  UseMethod("summary_heterogeneity")
}


# ---------------------------------------------------------------------------- #
# summary_heterogeneity.brma
# ---------------------------------------------------------------------------- #

#' @title Summary of Heterogeneity for brma Objects
#'
#' @description Computes the absolute heterogeneity (tau, tau^2) and
#' relative measures of heterogeneity (I^2, H^2) for a fitted brma object.
#'
#' @param object a fitted brma object
#' @param probs quantiles of the posterior distribution to be displayed.
#' Defaults to \code{c(.025, .975)} for 95% credible intervals.
#' @param component heterogeneity component to return for \code{brma.mv()}
#' models. Defaults to \code{"all"}. Use \code{"total"} for the
#' variance-additive total heterogeneity. Random-formula allocation nodes can
#' be selected by their displayed name, such as \code{"study/esid"}.
#' @param ... additional arguments (currently ignored)
#'
#' @details
#' For standard (2-level) random-effects models, the function reports:
#' \itemize{
#'   \item \code{tau}: Between-study standard deviation
#'   \item \code{tau2}: Between-study variance
#'   \item \code{I2}: Percentage of total variance due to heterogeneity
#'   \item \code{H2}: Ratio of total to sampling variance
#' }
#'
#' For multilevel (3-level) models with nested effects, the function additionally
#' partitions heterogeneity into estimate-level and cluster-level components:
#' \itemize{
#'   \item \code{rho}: Proportion of heterogeneity variance allocated to clusters
#'   \item \code{tau [within]}: Estimate-level standard deviation
#'   \item \code{tau [between]}: Cluster-level standard deviation
#'   \item \code{tau2 [within]}: Estimate-level variance
#'   \item \code{tau2 [between]}: Cluster-level variance
#'   \item \code{I2 [within]}: Percentage of variance due to estimate-level heterogeneity
#'   \item \code{I2 [between]}: Percentage of variance due to cluster-level heterogeneity
#' }
#'
#' For location-scale models, \code{tau2} aggregates the observation-specific
#' heterogeneity variances \eqn{\tau_i^2}; the corresponding \code{tau} summary
#' is the square root of this aggregate variance. The relative \eqn{I^2} and
#' \eqn{H^2} measures average the observation-specific indices.
#'
#' The I^2 and H^2 statistics are computed following the metafor package
#' implementation, using the "typical" sampling variance formula from
#' \insertCite{higgins2002quantifying;textual}{RoBMA}. For multilevel models,
#' the partitioned I^2 follows the approach described in the metafor documentation.
#'
#' For \code{brma.mv()} models, the method reports \code{sd} and \code{var}
#' for each selected random-effect component. A genuine variance-additive
#' aggregate is reported as \code{sd_total} and \code{var_total}; a
#' mean-variance allocation scale is reported as \code{sd_common} and
#' \code{var_common}. This vocabulary also applies when no explicit
#' \code{random} formula is supplied because \code{brma.mv()} adds an
#' estimate-level random effect by default.
#' When a known group covariance \eqn{R} is supplied, \code{sd} is the fitted
#' multiplier of that covariance kernel. It need not equal every row's
#' marginal standard deviation when \eqn{\mathrm{diag}(R)} is not one.
#' Row-specific marginal standard deviations remain available from
#' \code{predict(type = "terms.scale")} and in covariance-based prediction and
#' diagnostics.
#' Relative \eqn{I^2} and \eqn{H^2} summaries are not reported for general
#' known-V covariance structures. For \code{component = "total"}, independent
#' component variances are summed before reporting \code{var_total};
#' \code{sd_total} is the square root of this variance draw.
#' When a random-formula model uses a shared total-SD plus variance-allocation
#' node, \code{component = "all"} also includes an allocation-node table with
#' the appropriate aggregate SD and variance plus
#' \code{var_prop(<block>)} rows.
#' For nested formulas
#' such as \code{random = ~ 1 | study / esid}, this table is displayed under
#' the user-facing path \code{"study/esid"}.
#' For heterogeneous structured random effects, such as
#' \code{random = ~ har(time | study)}, the allocation table reports
#' the exhaustive semantic allocation family: aggregate SD and variance,
#' level-specific component SDs and variances, correlations, SD multipliers, and
#' variance multipliers. A redundant block owner is omitted for a bare formula or
#' unnamed one-entry list and retained for explicitly named one-entry lists and
#' multiple components.
#'
#' @return A list of class \code{summary_heterogeneity.brma} containing:
#' \itemize{
#'   \item \code{estimates}: A \code{BayesTools_table} with heterogeneity statistics
#'   \item \code{component}: The public heterogeneity-component name
#' }
#' For decomposed \code{brma.mv()} models with multiple selected components,
#' a named list of such summary objects is returned.
#'
#' @examples \dontrun{
#' if (requireNamespace("metadat", quietly = TRUE)) {
#'   data(dat.lehmann2018, package = "metadat")
#'   fit <- brma(
#'     yi      = yi,
#'     vi      = vi,
#'     data    = dat.lehmann2018,
#'     measure = "SMD",
#'     seed    = 1,
#'     silent  = TRUE
#'   )
#'
#'   summary_heterogeneity(fit)
#' }
#' }
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso [pooled_heterogeneity()], [summary.brma()]
#' @export
summary_heterogeneity.brma <- function(object, probs = c(.025, .975),
                                       component = "all", ...) {

  # input validation
  BayesTools::check_real(probs, "probs", allow_NULL = TRUE, check_length = 0)
  if (inherits(object, "brma.mv")) {
    return(
      .summary_heterogeneity_brma_mv(
        object    = object,
        probs     = probs,
        component = component,
        ...
      )
    )
  }

  .check_univariate_heterogeneity_component(component)

  # extract model characteristics
  is_multilevel <- .is_multilevel(object)
  is_scale      <- .is_scale(object)

  # extract sampling variances using helper
  vi <- .outcome_data_vi(object)
  K  <- length(vi)

  # compute typical sampling variance (v_tilde)
  # following metafor formula (Higgins & Thompson, 2002, eq. 9)
  # using the generalized formula with projection matrix for meta-regression
  X <- .get_model_matrix(object)
  p <- ncol(X)
  W <- diag(1/vi, nrow = K)
  P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W

  v_tilde <- (K - p) / sum(diag(P))

  posterior_samples <- .get_posterior_samples(object[["fit"]])

  # extract tau samples using the evaluate helper
  tau_result <- .evaluate.brma.tau(
    fit               = object[["fit"]],
    scale_data        = object[["data"]][["scale"]],
    scale_formula     = if (is_scale) .create_fit_formula_list(data = object[["data"]], "scale") else NULL,
    scale_priors      = object[["priors"]][["scale"]],
    is_scale          = is_scale,
    is_multilevel     = is_multilevel,
    K                 = K,
    posterior_samples = posterior_samples,
    fixed_tau         = .fixed_tau_prior_value(object[["priors"]]),
    fixed_rho         = .fixed_rho_prior_value(object[["priors"]])
  )

  samples_list <- .summary_heterogeneity_samples(
    tau_within_samples  = tau_result[["tau_within"]],
    tau_between_samples = tau_result[["tau_between"]],
    rho_samples         = if (is_multilevel) tau_result[["rho"]] else NULL,
    v_tilde             = v_tilde,
    is_multilevel       = is_multilevel
  )

  # generate summary table
  estimates <- BayesTools::ensemble_estimates_table(
    samples    = samples_list,
    parameters = names(samples_list),
    probs      = probs,
    title      = "Heterogeneity Estimates:"
  )

  # create output object
  output <- list(
    estimates = estimates,
    component = "location"
  )

  class(output) <- "summary_heterogeneity.brma"

  return(output)
}


# ---------------------------------------------------------------------------- #
# .summary_heterogeneity_samples
# ---------------------------------------------------------------------------- #
#
# Compute heterogeneity summaries from observation-level tau samples.
#
.summary_heterogeneity_samples <- function(tau_within_samples, tau_between_samples,
                                           v_tilde, is_multilevel,
                                           rho_samples = NULL) {

  if (is.null(tau_between_samples)) {
    tau_between_samples <- matrix(0, nrow = nrow(tau_within_samples),
                                  ncol = ncol(tau_within_samples))
  }

  if (!is.matrix(tau_within_samples) || !is.matrix(tau_between_samples)) {
    stop("'tau_within_samples' and 'tau_between_samples' must be matrices.",
         call. = FALSE)
  }
  if (!identical(dim(tau_within_samples), dim(tau_between_samples))) {
    stop("'tau_within_samples' and 'tau_between_samples' must have matching dimensions.",
         call. = FALSE)
  }
  if (!is.numeric(v_tilde) || length(v_tilde) != 1 || !is.finite(v_tilde) || v_tilde <= 0) {
    stop("'v_tilde' must be a positive finite number.", call. = FALSE)
  }
  sigma2_within_matrix  <- tau_within_samples^2
  sigma2_between_matrix <- tau_between_samples^2
  sigma2_total_matrix   <- sigma2_within_matrix + sigma2_between_matrix

  sigma2_within  <- rowMeans(sigma2_within_matrix)
  sigma2_between <- rowMeans(sigma2_between_matrix)
  sigma2_total   <- sigma2_within + sigma2_between

  denominator_matrix <- sigma2_total_matrix + v_tilde

  I2_total <- rowMeans(100 * sigma2_total_matrix / denominator_matrix)
  H2       <- rowMeans(denominator_matrix / v_tilde)

  if (is_multilevel) {

    I2_within  <- rowMeans(100 * sigma2_within_matrix / denominator_matrix)
    I2_between <- rowMeans(100 * sigma2_between_matrix / denominator_matrix)
    if (is.null(rho_samples)) {
      rho_samples <- rep(NA_real_, length(sigma2_total))
      identified  <- sigma2_total > 0
      rho_samples[identified] <- sigma2_between[identified] /
        sigma2_total[identified]
    } else {
      rho_samples <- .resolve_heterogeneity_allocation(
        rho       = rho_samples,
        n_samples = nrow(tau_within_samples),
        context   = "Heterogeneity summary"
      )
    }

    return(list(
      "tau"             = sqrt(sigma2_total),
      "rho"             = rho_samples,
      "tau [within]"    = sqrt(sigma2_within),
      "tau [between]"   = sqrt(sigma2_between),
      "tau2"            = sigma2_total,
      "tau2 [within]"   = sigma2_within,
      "tau2 [between]"  = sigma2_between,
      "I2"              = I2_total,
      "I2 [within]"     = I2_within,
      "I2 [between]"    = I2_between,
      "H2"              = H2
    ))
  }

  return(list(
    "tau"  = sqrt(sigma2_total),
    "tau2" = sigma2_total,
    "I2"   = I2_total,
    "H2"   = H2
  ))
}


# ---------------------------------------------------------------------------- #
# print method for summary.brma_heterogeneity
# ---------------------------------------------------------------------------- #

#' @title Print Summary of Heterogeneity
#'
#' @description Prints the heterogeneity summary table.
#'
#' @param x a \code{summary_heterogeneity.brma} object
#' @param ... additional arguments (currently ignored)
#'
#' @return Returns \code{x} invisibly.
#'
#' @export
print.summary_heterogeneity.brma <- function(x, ...) {

  cat("\n")
  print(x$estimates)
  cat("\n")

  return(invisible(x))
}


#' @rdname summary_heterogeneity.brma
#' @param x a \code{summary_heterogeneity.brma_list} object
#' @export
print.summary_heterogeneity.brma_list <- function(x, ...) {

  cat("\n")
  for (i in seq_along(x)) {
    print(x[[i]][["estimates"]])
    cat("\n")
  }

  return(invisible(x))
}


#' @title Convert Heterogeneity Summaries to Data Frames
#'
#' @description Converts heterogeneity summaries to plain long data frames with
#' leading \code{component} and \code{parameter} columns. A multi-component
#' result can instead be returned as a named list with
#' \code{format = "list"}.
#'
#' @param x a \code{summary_heterogeneity.brma} or
#' \code{summary_heterogeneity.brma_list} object.
#' @param row.names \code{NULL} or a character vector giving the row names for
#' the result. Custom row names are unsupported when \code{format = "list"}.
#' @param optional logical; passed to the final data-frame coercion.
#' @param format for multi-component results, whether to return one
#' \code{"long"} data frame or a named \code{"list"} of data frames.
#' @param stringsAsFactors accepted for compatibility with \code{data.frame()}.
#' @param ... unused additional arguments.
#'
#' @return A plain \code{data.frame}, or a named list of data frames when
#' \code{format = "list"}.
#'
#' @export
as.data.frame.summary_heterogeneity.brma <- function(
    x, row.names = NULL, optional = FALSE, stringsAsFactors = FALSE, ...) {

  component <- x[["component"]]
  if (!is.character(component) || length(component) != 1L ||
      is.na(component) || !nzchar(component)) {
    component <- "location"
  }
  output <- .output_table_as_long_data_frame(
    table            = x[["estimates"]],
    component        = component,
    row.names        = row.names,
    optional         = optional,
    stringsAsFactors = stringsAsFactors
  )

  return(output)
}


#' @rdname as.data.frame.summary_heterogeneity.brma
#' @export
as.data.frame.summary_heterogeneity.brma_list <- function(
    x, row.names = NULL, optional = FALSE, format = c("long", "list"),
    stringsAsFactors = FALSE, ...) {

  format <- match.arg(format)
  tables <- lapply(names(x), function(component) {
    value <- x[[component]]
    value[["component"]] <- component
    as.data.frame(value, stringsAsFactors = stringsAsFactors)
  })
  names(tables) <- names(x)
  if (identical(format, "list")) {
    if (!is.null(row.names)) {
      stop("'row.names' is unsupported when format = 'list'.", call. = FALSE)
    }
    return(tables)
  }

  output <- .output_bind_long_data_frames(
    tables    = unname(tables),
    row.names = row.names,
    optional  = optional
  )

  return(output)
}
