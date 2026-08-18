.hypothesis_brma_level_contrast_BF <- function(
    object, posterior, hypothesis, parameter, density_method, density_control,
    logBF, BF01, seed, columns, standardized_coefficients = FALSE) {

  if (.is_RoBMA(object)) {
    stop(
      "Cross-level point hypotheses are currently available only for a ",
      "single fitted model; model-averaged atom conditioning is not defined.",
      call. = FALSE
    )
  }
  target <- BayesTools::hypothesis_level_contrast(
    posterior  = posterior,
    hypothesis = hypothesis,
    parameter  = parameter
  )
  precomputed <- .density_method_uses_precomputed(
    density_method,
    allow_normal = TRUE
  )
  if (precomputed) {
    if (!inherits(object, "brma") || is.null(object[["fit"]])) {
      stop("A fitted single-model brma object is required for this target.",
           call. = FALSE)
    }
    if (isTRUE(standardized_coefficients)) {
      stop(
        "Cross-level qCMDE/IWMDE hypotheses require coefficients on the ",
        "fitted scale; use standardized_coefficients = FALSE.",
        call. = FALSE
      )
    }
    .check_iwmde_available(object, "qCMDE/IWMDE cross-level hypothesis()")
    .iwmde_check_point_ordinate_supported(object, density_method)
    if (is.null(density_control[["normalization_points"]])) {
      density_control[["normalization_points"]] <- max(
        50L,
        density_control[["n_points"]]
      )
    }

    point_refs <- BayesTools::hypothesis_parse_point_reference(
      hypothesis     = BayesTools::hypothesis_render(target[["hypothesis"]]),
      allow_compound = FALSE
    )
    values <- unique(point_refs[["value"]])
    if (length(values) == 0L) {
      stop("Internal error: the level contrast has no point ordinate.",
           call. = FALSE)
    }

    target_posterior <- target[["posterior"]]
    parameter_spec <- list(
      type          = "linear",
      weights       = target[["weights"]],
      prior_density = attr(
        target_posterior,
        "prior_density",
        exact = TRUE
      )
    )
    estimate_cache <- .iwmde_estimate_cache()
    conditional <- .iwmde_first_nonempty_attr(
      target_posterior,
      c("effective_conditional", "conditional")
    )
    for (value in values) {
      target_posterior <- .hypothesis_brma_attach_iwmde_scalar(
        posterior             = target_posterior,
        raw_posterior         = target_posterior,
        context               = .iwmde_context(object),
        estimate_cache        = estimate_cache,
        parameter             = target[["parameter"]],
        parameter_label       = parameter,
        value                 = value,
        conditional           = conditional,
        n_points              = density_control[["n_points"]],
        samples               = density_control[["samples"]],
        target_relative_mcse  = density_control[["target_relative_mcse"]],
        normalization_points  = density_control[["normalization_points"]],
        normalization_prob    = density_control[["normalization_prob"]],
        density_method        = density_method,
        parameter_spec        = parameter_spec
      )
    }
    target[["posterior"]] <- target_posterior
  }

  out <- BayesTools::hypothesis_BF(
    posterior      = target[["posterior"]],
    hypothesis     = target[["hypothesis"]],
    parameter      = target[["parameter"]],
    logBF          = logBF,
    BF01           = BF01,
    seed           = seed,
    columns        = columns,
    density_method = if (precomputed) "precomputed" else density_method
  )
  if (precomputed) {
    out <- .hypothesis_brma_append_iwmde_warnings(
      table     = out,
      posterior = target[["posterior"]]
    )
  }
  rownames(out) <- rep(parameter, nrow(out))

  return(out)
}
