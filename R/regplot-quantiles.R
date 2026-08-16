# ============================================================================ #
# regplot-quantiles.R
# ============================================================================ #

.plot_mixture_quantiles_native <- function(mean_samples, sd_samples, probs,
                                            weights, se = NULL,
                                            selected_rows = NULL,
                                            selection_context = NULL,
                                            caller) {

  mean_samples <- as.matrix(mean_samples)
  sd_samples   <- as.matrix(sd_samples)
  weights      <- as.numeric(weights)

  if (is.null(selection_context) || !any(selected_rows)) {
    return(.Call(
      "RoBMA_plot_normal_mixture_quantiles",
      .native_numeric_matrix(mean_samples),
      .native_numeric_matrix(sd_samples),
      .native_numeric_vector(probs),
      .native_numeric_vector(weights),
      PACKAGE = "RoBMA"
    ))
  }

  .selection_require_step_evaluable(selection_context, caller)
  native_static <- BayesTools::selection_native_static_args(selection_context)

  return(.Call(
    "RoBMA_plot_selnorm_mixture_quantiles",
    .native_numeric_matrix(mean_samples),
    .native_numeric_matrix(sd_samples),
    .native_numeric_vector(se),
    .native_numeric_vector(probs),
    .native_numeric_vector(weights),
    .native_integer_vector(selected_rows),
    .native_numeric_matrix(selection_context[["omega"]]),
    .native_numeric_vector(selection_context[["alpha"]]),
    .native_integer_vector(selection_context[["phack_kind"]]),
    .native_integer_vector(selection_context[["kernel_mode"]]),
    native_static[["z_lower"]],
    native_static[["z_upper"]],
    native_static[["sign"]],
    native_static[["phack_q"]],
    native_static[["phack_z_source"]],
    native_static[["phack_z_dest"]],
    native_static[["segment_bounds"]],
    native_static[["segment_step_bin"]],
    native_static[["segment_phack_region"]],
    native_static[["telescope_probabilities"]],
    PACKAGE = "RoBMA"
  ))
}


.regplot_mixture_interval_quantiles <- function(mean_samples, sd_samples,
                                                probs) {

  mean_samples <- as.matrix(mean_samples)
  sd_samples   <- as.matrix(sd_samples)
  quantiles <- .plot_mixture_quantiles_native(
    mean_samples = mean_samples,
    sd_samples   = sd_samples,
    probs        = probs,
    weights      = rep(1, nrow(mean_samples)),
    caller       = ".regplot_mixture_interval_quantiles()"
  )

  return(list(lower = quantiles[, 1L], upper = quantiles[, 2L]))
}


.regplot_selection_mixture_interval_quantiles <- function(
    x, mean_samples, sd_samples, se, probs, posterior_samples = NULL) {

  mean_samples <- as.matrix(mean_samples)
  sd_samples   <- as.matrix(sd_samples)
  setup        <- .regplot_selection_setup(
    x                 = x,
    posterior_samples = posterior_samples
  )
  if (!is.numeric(se) || !length(se) %in% c(1L, ncol(mean_samples)) ||
      any(!is.finite(se)) || any(se <= 0)) {
    stop(
      "Selection-model regression-plot standard errors must be positive and ",
      "have length one or one value per prediction.",
      call. = FALSE
    )
  }

  quantiles <- .plot_mixture_quantiles_native(
    mean_samples      = mean_samples,
    sd_samples        = sd_samples,
    se                = se,
    probs             = probs,
    weights           = rep(1, nrow(mean_samples)),
    selected_rows     = setup[["is_weightfunction"]],
    selection_context = setup[["selection"]],
    caller            = ".regplot_selection_mixture_interval_quantiles()"
  )

  return(list(lower = quantiles[, 1L], upper = quantiles[, 2L]))
}


.regplot_selection_setup <- function(x, posterior_samples = NULL) {

  if (is.null(posterior_samples)) {
    posterior_samples <- .get_posterior_samples(x[["fit"]])
  }
  S              <- nrow(posterior_samples)
  bias_indicator <- .extract_bias_indicator(
    x,
    posterior_samples = posterior_samples
  )
  selection      <- .selection_context(
    object            = x,
    posterior_samples = posterior_samples
  )
  use_normal     <- if (is.null(selection)) {
    rep(TRUE, S)
  } else {
    selection[["use_normal"]]
  }

  return(list(
    bias_indicator    = bias_indicator,
    is_weightfunction = !use_normal,
    selection         = selection
  ))
}
