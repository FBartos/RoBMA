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
  # The constructor has checked positive acceptance for retained contexts.
  # All-conditioned selection leaves the marginal Gaussian law unchanged.
  selection      <- if (.selection_all_sources_conditioned(x[["data"]])) {
    NULL
  } else .selection_context(
    object            = x,
    posterior_samples = posterior_samples
  )
  .plot_check_scalar_selection_target(x, selection, "regplot")
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


# Scalar selected quantiles have no retained-context or joint-event integral.
# Keep this availability check separate from choosing a full-event contour grid.
.plot_check_scalar_selection_target <- function(x, selection, family) {

  data <- x[["data"]]
  if (is.null(selection) || all(selection[["use_normal"]]) ||
      !.is_data_joint_selection(data) || .selection_all_sources_conditioned(data)) {
    return(invisible(NULL))
  }
  omega <- selection[["omega"]]
  constant <- omega[, 1L] > 0 & rowSums(omega != omega[, 1L]) == 0L
  selected <- which(!selection[["use_normal"]] & !constant)
  if (!length(selected)) return(invisible(NULL))

  model <- .data_selection_model(data)
  plan <- .data_selection_execution_plan(data)
  rules <- rep_len(selection[["vector_rule"]], nrow(omega))[selected]
  product_rule <- all(rules == 0L)
  # Conditioned sources no longer block a funnel contour: the band is the
  # mixture over their population law, which `.funnel_conditioned_expansion()`
  # evaluates. A funnel contour is the law of a new, independent estimate at
  # each hypothetical standard error; under the product rule that estimate is
  # its own dependency block, so the fitted rows' integration blocks do not
  # enter its scalar law. A best rule weighs an estimate against the others of
  # its publication and still needs that design, and regression-plot intervals
  # keep the stricter scalar rule over the fitted rows themselves.
  supports_conditioned <- identical(family, "funnel")
  available <- (supports_conditioned || !.selection_retains_sampling(data)) &&
    ((supports_conditioned && product_rule) ||
       all(lengths(plan[["row_blocks"]]) == 1L))
  if (!product_rule && any(lengths(model[["groups"]][["row_blocks"]]) != 1L)) {
    available <- FALSE
  }
  # A changed regression design can activate shared coefficient supports
  # that were disjoint in the fitted rows. No new-event certificate is passed.
  requires_zero <- if (supports_conditioned) {
    list()
  } else {
    Filter(function(source) isTRUE(source[["retained"]]) ||
      (identical(family, "regplot") && identical(source[["role"]], "other")),
      model[["sources"]][["random"]])
  }
  if (available && length(requires_zero)) {
    if (.is_data_random(data)) {
      design <- .fitted_formula_design(x, "mu", required = TRUE)
      terms <- design[["random_effects"]]
      term_names <- vapply(terms, .random_effect_term_block_name, character(1L))
      indices <- match(vapply(requires_zero, `[[`, character(1L), "name"), term_names)
      available <- !anyNA(indices) && all(vapply(terms[indices],
        .marglik_random_effect_fixed_zero, logical(1L), data = data,
        prior_list = design[["prior_list"]], K = nrow(data[["outcome"]])))
    } else {
      tau <- if (.is_data_scale(data)) NULL else .fixed_tau_prior_value(x[["priors"]])
      rho <- .fixed_rho_prior_value(x[["priors"]])
      available <- all(vapply(requires_zero, function(source) {
        isTRUE(tau == 0) || (.is_data_multilevel(data) &&
          ((identical(source[["role"]], "estimate") && isTRUE(rho == 1)) ||
           (identical(source[["role"]], "other") && isTRUE(rho == 0))))
      }, logical(1L)))
    }
  }
  if (!available) {
    if (identical(family, "funnel")) {
      stop("Selected funnel contours are unavailable for this joint selection configuration. ",
        "Set 'sampling_bias = FALSE', or use 'zplot()' to view its marginal selected distribution.",
        call. = FALSE)
    }
    stop("Selected regression-plot sampling intervals are unavailable for this joint selection configuration. ",
      "Set 'sampling_bias = FALSE' to draw bias-adjusted sampling intervals.",
      call. = FALSE)
  }
  invisible(NULL)
}
