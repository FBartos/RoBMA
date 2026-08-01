# ---------------------------------------------------------------------------- #
# .regplot_get_moderator
# ---------------------------------------------------------------------------- #
#
# Identify and validate the moderator variable for plotting.
#
# @param x     brma object
# @param mod   index or name of moderator (NULL for auto-detect)
#
# @return list with name, type, and data for the moderator
#
# ---------------------------------------------------------------------------- #
.regplot_get_moderator <- function(x, mod) {

  mods_data <- x[["data"]][["mods"]]
  design    <- .fitted_formula_design(x, "mu", required = TRUE)

  if (is.null(mods_data)) {
    stop("No moderator data found in the model.", call. = FALSE)
  }

  mod_names <- design[["predictors"]]

  if (is.null(mod)) {
    # auto-detect: use single moderator or error if multiple
    if (length(mod_names) == 1) {
      mod_name <- mod_names[1]
    } else {
      stop("Multiple moderators present. Please specify 'mod' argument. ",
           "Available moderators: ", paste(mod_names, collapse = ", "),
           call. = FALSE)
    }
  } else if (is.numeric(mod)) {
    # index-based selection
    if (mod < 1 || mod > length(mod_names)) {
      stop("Invalid moderator index. Must be between 1 and ",
           length(mod_names), ".", call. = FALSE)
    }
    mod_name <- mod_names[mod]
  } else if (is.character(mod)) {
    # name-based selection
    if (!mod %in% mod_names) {
      stop("Moderator '", mod, "' not found. ",
           "Available moderators: ", paste(mod_names, collapse = ", "),
           call. = FALSE)
    }
    mod_name <- mod
  } else {
    stop("'mod' must be NULL, a character name, or a numeric index.",
         call. = FALSE)
  }

  # extract moderator data and determine type
  mod_values <- mods_data[[mod_name]]
  predictor_type <- design[["predictor_types"]][[mod_name]]

  if (identical(predictor_type, "factor") ||
      is.factor(mod_values) || is.character(mod_values)) {
    mod_type <- "categorical"
    mod_values <- factor(
      mod_values,
      levels = if (!is.null(design[["xlevels"]][[mod_name]])) {
        design[["xlevels"]][[mod_name]]
      } else {
        levels(factor(mod_values))
      }
    )
  } else {
    mod_type <- "continuous"
  }

  return(list(
    name = mod_name,
    type = mod_type,
    data = mod_values
  ))
}


# ---------------------------------------------------------------------------- #
# .regplot_get_by
# ---------------------------------------------------------------------------- #
#
# Identify and validate a grouping moderator for interaction displays.
#
# ---------------------------------------------------------------------------- #
.regplot_get_by <- function(x, by, mod_name) {

  mods_data <- x[["data"]][["mods"]]
  design    <- .fitted_formula_design(x, "mu", required = TRUE)

  if (is.null(by)) {
    by <- .regplot_detect_interaction_by(x, mod_name)
  }

  if (is.null(by)) {
    return(NULL)
  }

  mod_names <- design[["predictors"]]

  if (!is.character(by) || length(by) != 1L) {
    stop("'by' must be NULL or a single moderator name.", call. = FALSE)
  }
  if (by == mod_name) {
    stop("'by' must be different from 'mod'.", call. = FALSE)
  }
  if (!by %in% mod_names) {
    stop("Moderator '", by, "' not found. ",
         "Available moderators: ", paste(mod_names, collapse = ", "),
         call. = FALSE)
  }

  by_values <- mods_data[[by]]
  predictor_type <- design[["predictor_types"]][[by]]

  if (identical(predictor_type, "factor") ||
      is.factor(by_values) || is.character(by_values)) {
    by_type <- "categorical"
    by_values <- factor(
      by_values,
      levels = if (!is.null(design[["xlevels"]][[by]])) {
        design[["xlevels"]][[by]]
      } else {
        levels(factor(by_values))
      }
    )
  } else {
    by_type <- "continuous"
  }

  return(list(
    name = by,
    type = by_type,
    data = by_values
  ))
}


# ---------------------------------------------------------------------------- #
# .regplot_detect_interaction_by
# ---------------------------------------------------------------------------- #
#
# Auto-select the second variable from a unique two-way interaction containing
# the plotted moderator.
#
# ---------------------------------------------------------------------------- #
.regplot_detect_interaction_by <- function(x, mod_name) {

  mods_data <- x[["data"]][["mods"]]

  if (is.null(mods_data)) {
    return(NULL)
  }

  term_labels <- .fitted_formula_terms(
    object            = x,
    parameter         = "mu",
    include_intercept = FALSE,
    display           = TRUE
  )
  term_labels <- term_labels[grepl(":", term_labels, fixed = TRUE)]

  if (length(term_labels) == 0L) {
    return(NULL)
  }

  candidates <- character()

  for (term in term_labels) {
    variables <- strsplit(term, ":", fixed = TRUE)[[1]]
    variables <- trimws(variables)

    if (length(variables) == 2L && mod_name %in% variables) {
      candidates <- c(candidates, setdiff(variables, mod_name))
    }
  }

  candidates <- unique(candidates[candidates %in% names(mods_data)])

  if (length(candidates) == 0L) {
    return(NULL)
  }
  if (length(candidates) > 1L) {
    stop("The selected moderator enters multiple interactions. ",
         "Please specify 'by'. Candidate moderators: ",
         paste(candidates, collapse = ", "), ".", call. = FALSE)
  }

  return(candidates)
}


# ---------------------------------------------------------------------------- #
# .regplot_band_data_continuous
# ---------------------------------------------------------------------------- #
#
# Build continuous interval-band data with one row per polygon vertex while
# preserving the underlying interval bounds for programmatic access.
#
# @param xpred numeric vector of prediction x values
# @param lower numeric vector of lower interval bounds
# @param upper numeric vector of upper interval bounds
#
# @return data.frame with consistent polygon and interval metadata
#
# ---------------------------------------------------------------------------- #
.regplot_band_data_continuous <- function(xpred, lower, upper, group = "All",
                                          group_id = 1L) {

  band_x     <- c(xpred, rev(xpred))
  band_lower <- c(lower, rev(lower))
  band_upper <- c(upper, rev(upper))

  data.frame(
    x     = unname(band_x),
    y     = unname(c(lower, rev(upper))),
    lower = unname(band_lower),
    upper = unname(band_upper),
    xpred = unname(band_x),
    group = rep(group, length(band_x)),
    group_id = rep(group_id, length(band_x))
  )
}


# ---------------------------------------------------------------------------- #
# .regplot_prediction_grid
# ---------------------------------------------------------------------------- #
#
# Build moderator values for predictions. The selected moderator forms the
# x-axis; an optional `by` moderator expands the grid into interaction lines.
#
# ---------------------------------------------------------------------------- #
.regplot_prediction_grid <- function(x, mod_name, mod_type, mod_data,
                                     by_info, at, digits) {

  mods_data <- x[["data"]][["mods"]]
  design    <- .fitted_formula_design(x, "mu", required = TRUE)

  if (mod_type == "continuous") {
    if (is.null(at)) {
      x_range <- range(mod_data)
      at_pred <- seq(x_range[1], x_range[2], length.out = 101)
    } else {
      at_pred <- at
    }
    x_values  <- at_pred
    x_numeric <- at_pred
    x_levels  <- NULL
  } else {
    x_levels  <- levels(mod_data)
    at_pred   <- x_levels
    x_values  <- factor(at_pred, levels = x_levels)
    x_numeric <- seq_along(at_pred)
  }

  if (is.null(by_info)) {
    by_values <- list(
      values = NULL,
      labels = "All",
      levels = NULL
    )
  } else {
    by_values <- .regplot_by_values(by_info, digits)
  }

  grid <- expand.grid(
    x_id     = seq_along(at_pred),
    group_id = seq_along(by_values[["labels"]]),
    KEEP.OUT.ATTRS = FALSE
  )
  grid[["group"]] <- by_values[["labels"]][grid[["group_id"]]]
  grid[["x"]]     <- x_numeric[grid[["x_id"]]]
  if (mod_type == "categorical") {
    grid[["level"]] <- at_pred[grid[["x_id"]]]
  }

  n_pred       <- nrow(grid)
  newdata_list <- list()

  for (nm in design[["predictors"]]) {
    if (nm == mod_name) {
      if (mod_type == "continuous") {
        newdata_list[[nm]] <- x_values[grid[["x_id"]]]
      } else {
        newdata_list[[nm]] <- factor(at_pred[grid[["x_id"]]], levels = x_levels)
      }
    } else if (!is.null(by_info) && nm == by_info[["name"]]) {
      if (by_info[["type"]] == "continuous") {
        newdata_list[[nm]] <- by_values[["values"]][grid[["group_id"]]]
      } else {
        newdata_list[[nm]] <- factor(
          by_values[["values"]][grid[["group_id"]]],
          levels = by_values[["levels"]]
        )
      }
    } else {
      other_vals <- mods_data[[nm]]
      if (is.factor(other_vals) || is.character(other_vals)) {
        other_levels <- if (!is.null(design[["xlevels"]][[nm]])) {
          design[["xlevels"]][[nm]]
        } else {
          levels(factor(other_vals))
        }
        newdata_list[[nm]] <- factor(
          rep(other_levels[1], n_pred),
          levels = other_levels
        )
      } else {
        newdata_list[[nm]] <- rep(mean(other_vals), n_pred)
      }
    }
  }

  return(list(
    newdata  = as.data.frame(newdata_list),
    grid     = grid,
    at_pred  = at_pred,
    groups   = by_values[["labels"]],
    x_levels = x_levels
  ))
}


# ---------------------------------------------------------------------------- #
# .regplot_by_values
# ---------------------------------------------------------------------------- #
#
# Values used for the grouping moderator.
#
# ---------------------------------------------------------------------------- #
.regplot_by_values <- function(by_info, digits) {

  by_data <- by_info[["data"]]

  if (by_info[["type"]] == "continuous") {
    by_mean <- mean(by_data)
    by_sd   <- stats::sd(by_data)
    values  <- c(by_mean - by_sd, by_mean, by_mean + by_sd)
    labels  <- c("Mean - 1 SD", "Mean", "Mean + 1 SD")

    if (!is.finite(by_sd) || by_sd == 0) {
      values <- by_mean
      labels <- format(round(by_mean, digits = digits), trim = TRUE)
    }

    labels <- paste0(by_info[["name"]], ": ", labels)
    return(list(values = values, labels = labels, levels = NULL))
  }

  by_factor <- factor(by_data)

  return(list(
    values = levels(by_factor),
    labels = levels(by_factor),
    levels = levels(by_factor)
  ))
}


# ---------------------------------------------------------------------------- #
# .regplot_observed_groups
# ---------------------------------------------------------------------------- #
#
# Assign observed studies to plotting groups.
#
# ---------------------------------------------------------------------------- #
.regplot_observed_groups <- function(by_info, groups, n) {

  if (is.null(by_info)) {
    return(list(
      group    = rep("All", n),
      group_id = rep(1L, n)
    ))
  }

  by_data <- by_info[["data"]]

  if (by_info[["type"]] == "categorical") {
    group <- as.character(factor(by_data, levels = groups))
  } else {
    by_mean <- mean(by_data)
    by_sd   <- stats::sd(by_data)

    if (!is.finite(by_sd) || by_sd == 0 || length(groups) == 1L) {
      group <- rep(groups[1], length(by_data))
    } else {
      by_values <- c(by_mean - by_sd, by_mean, by_mean + by_sd)
      breaks    <- c(-Inf, mean(by_values[1:2]), mean(by_values[2:3]), Inf)
      group     <- as.character(cut(by_data, breaks = breaks, labels = groups))
    }
  }

  return(list(
    group    = group,
    group_id = match(group, groups)
  ))
}


# ---------------------------------------------------------------------------- #
# .regplot_add_dummy_outcome
# ---------------------------------------------------------------------------- #
#
# Add outcome columns for prediction paths that use sampling bias inputs.
#
# ---------------------------------------------------------------------------- #
.regplot_add_dummy_outcome <- function(x, newdata, n_pred, sei) {

  outcome_type <- .outcome_type(x)

  if (outcome_type == "norm") {
    newdata[["yi"]]  <- rep(0, n_pred)
    newdata[["sei"]] <- rep(sei, n_pred)
  } else if (outcome_type == "bin") {
    newdata[["ai"]]  <- rep(0, n_pred)
    newdata[["ci"]]  <- rep(0, n_pred)
    newdata[["n1i"]] <- rep(0, n_pred)
    newdata[["n2i"]] <- rep(0, n_pred)
  } else if (outcome_type == "pois") {
    newdata[["x1i"]] <- rep(0, n_pred)
    newdata[["x2i"]] <- rep(0, n_pred)
    newdata[["t1i"]] <- rep(0, n_pred)
    newdata[["t2i"]] <- rep(0, n_pred)
  }

  return(newdata)
}


# ---------------------------------------------------------------------------- #
# .regplot_add_random_prediction_columns
# ---------------------------------------------------------------------------- #
#
# Add missing random-formula variables to the plotting grid. Missing grouping
# variables get synthetic row-unique levels so pointwise bands do not imply a
# joint curve covariance. Missing random-slope variables fall back to the same
# representative-value convention used for fixed moderators.
#
# ---------------------------------------------------------------------------- #
.regplot_add_random_prediction_columns <- function(x, newdata) {

  if (!inherits(x, "brma.mv") || !.is_random(x)) {
    return(newdata)
  }

  design <- .fitted_formula_design(x, "mu", required = TRUE)
  terms  <- design[["random_effects"]]
  if (length(terms) == 0L) {
    return(newdata)
  }

  original_location <- x[["data"]][["location"]]
  n_pred            <- nrow(newdata)

  for (term in terms) {
    group_vars <- all.vars(term[["group_expr"]])
    slope_vars <- all.vars(term[["term_formula"]])
    term_vars  <- unique(c(group_vars, slope_vars))

    for (var in term_vars) {
      if (var %in% names(newdata)) {
        next
      }
      if (var %in% group_vars) {
        newdata[[var]] <- paste0(".RoBMA_regplot_", var, "_", seq_len(n_pred))
      } else {
        newdata[[var]] <- .regplot_representative_variable(
          original_location[[var]],
          n_pred      = n_pred,
          name        = var,
          random_term = term
        )
      }
    }
  }

  return(newdata)
}


.regplot_representative_variable <- function(x, n_pred, name, random_term = NULL) {

  if (is.null(x)) {
    stop(
      "regplot() cannot construct random-effect prediction data because ",
      "variable '", name, "' is missing from the fitted data.",
      call. = FALSE
    )
  }

  term_levels <- NULL
  term_type   <- NULL
  if (!is.null(random_term)) {
    term_levels <- random_term[["xlevels"]][[name]]
    term_type   <- random_term[["model_terms_type"]][[name]]
  }
  if (!is.null(term_levels) || identical(term_type, "factor")) {
    if (is.null(term_levels)) {
      term_levels <- levels(factor(x))
    }
    return(factor(rep(term_levels[1L], n_pred), levels = term_levels))
  }

  if (is.factor(x)) {
    return(factor(rep(levels(x)[1L], n_pred), levels = levels(x)))
  }
  if (is.character(x)) {
    values <- unique(x)
    return(rep(values[1L], n_pred))
  }
  if (is.logical(x)) {
    return(rep(x[which(!is.na(x))[1L]], n_pred))
  }
  if (is.numeric(x) || is.integer(x)) {
    return(rep(mean(x, na.rm = TRUE), n_pred))
  }

  stop(
    "regplot() cannot construct a representative value for random-effect ",
    "variable '", name, "'.",
    call. = FALSE
  )
}


.regplot_mv_random_sd_samples <- function(x, newdata, posterior_samples) {

  if (!inherits(x, "brma.mv") || !.is_random(x)) {
    return(NULL)
  }

  location_data   <- .regplot_add_random_prediction_columns(x, newdata)
  location_priors <- attr(x[["fit"]], "prior_list")
  if (is.null(location_priors)) {
    location_priors <- x[["priors"]][["location"]]
  }

  design         <- .fitted_formula_design(x, "mu", required = TRUE)
  sampled_blocks <- .formula_design_sampled_random_effect_blocks(design)
  variance       <- matrix(0, nrow = nrow(posterior_samples), ncol = nrow(newdata))

  if (length(sampled_blocks) > 0L) {
    prediction <- tryCatch(
      BayesTools::JAGS_predict_formula(
        fit             = .posterior_formula_fit(
          fit               = x[["fit"]],
          posterior_samples = posterior_samples,
          formula_design    = TRUE
        ),
        formula         = .create_fit_formula_list(data = x[["data"]], parameter = "location"),
        parameter       = "mu",
        data            = location_data,
        prior_list      = location_priors,
        formula_target  = "marginal",
        blocks          = sampled_blocks,
        marginal_method = "covariance",
        new_levels      = "sample"
      ),
      error = function(e) {
        if (grepl("row-indexed external SD source", conditionMessage(e), fixed = TRUE) &&
            .regplot_mv_random_intercept_only(x)) {
          return(NULL)
        }
        stop(e)
      }
    )

    if (is.null(prediction)) {
      return(.regplot_mv_random_pointwise_total_sd_samples(
        x                 = x,
        newdata           = newdata,
        posterior_samples = posterior_samples
      ))
    } else {
      vcov_samples <- prediction[["vcov"]][["samples"]]
      if (length(dim(vcov_samples)) != 3L ||
          dim(vcov_samples)[1L] != nrow(posterior_samples) ||
          dim(vcov_samples)[2L] != nrow(newdata) ||
          dim(vcov_samples)[3L] != nrow(newdata)) {
        stop("Random-effect covariance samples have inconsistent dimensions.",
             call. = FALSE)
      }

      sampled_variance <- t(vapply(
        seq_len(dim(vcov_samples)[1L]),
        function(s) diag(vcov_samples[s, , ]),
        numeric(nrow(newdata))
      ))
      sampled_variance <- .regplot_validate_variance(
        sampled_variance,
        "Sampled random-effect variance"
      )
      variance <- variance + sampled_variance
    }
  }

  marginalized_variance <- .regplot_mv_marginalized_variance_samples(
    x                 = x,
    newdata           = newdata,
    posterior_samples = posterior_samples
  )
  marginalized_variance <- .regplot_validate_variance(
    marginalized_variance,
    "Marginalized random-effect variance"
  )
  variance <- variance + marginalized_variance

  return(.regplot_variance_sd(
    variance,
    "Combined random-effect variance"
  ))
}


.regplot_variance_sd <- function(variance, context) {

  return(sqrt(.regplot_validate_variance(variance, context)))
}


.regplot_validate_variance <- function(variance, context) {

  if (!is.numeric(variance) || any(!is.finite(variance)) ||
      any(variance < 0)) {
    stop(context, " must be finite and non-negative.", call. = FALSE)
  }

  return(variance)
}


.regplot_mv_random_intercept_only <- function(x) {

  design <- .fitted_formula_design(x, "mu", required = TRUE)
  terms  <- design[["random_effects"]]
  if (length(terms) == 0L) {
    return(FALSE)
  }

  all(vapply(terms, function(term) {
    model_matrix <- term[["model_matrix"]]
    identical(as.integer(term[["n_columns"]]), 1L) &&
      is.matrix(model_matrix) &&
      ncol(model_matrix) == 1L &&
      all(is.finite(model_matrix)) &&
      all(model_matrix[, 1L] == 1)
  }, logical(1)))
}


.regplot_mv_random_pointwise_total_sd_samples <- function(x, newdata,
                                                         posterior_samples) {

  if (!.is_scale(x)) {
    stop(
      "regplot() cannot compute newdata random-effect bands for this ",
      "row-indexed random-formula model because it has no scale formula to ",
      "replay on the prediction grid.",
      call. = FALSE
    )
  }

  new_data_prepared <- .prepare_newdata(
    object  = x,
    newdata = newdata,
    type    = "estimate"
  )
  tau_result <- .evaluate.brma.tau(
    fit               = x[["fit"]],
    scale_data        = new_data_prepared[["scale"]],
    scale_formula     = .create_fit_formula_list(data = new_data_prepared, "scale"),
    scale_priors      = x[["priors"]][["scale"]],
    is_scale          = TRUE,
    is_multilevel     = FALSE,
    K                 = nrow(newdata),
    posterior_samples = posterior_samples
  )

  return(tau_result[["tau_total"]])
}


.regplot_mv_marginalized_variance_samples <- function(x, newdata,
                                                      posterior_samples) {

  S        <- nrow(posterior_samples)
  K        <- nrow(newdata)
  variance <- matrix(0, nrow = S, ncol = K)

  if (!.data_has_marginalized_random_effects(x[["data"]])) {
    return(variance)
  }

  if (.is_scale(x)) {
    new_data_prepared <- .prepare_newdata(
      object  = x,
      newdata = newdata,
      type    = "estimate"
    )
    tau_result <- .evaluate.brma.tau(
      fit               = x[["fit"]],
      scale_data        = new_data_prepared[["scale"]],
      scale_formula     = .create_fit_formula_list(data = new_data_prepared, "scale"),
      scale_priors      = x[["priors"]][["scale"]],
      is_scale          = TRUE,
      is_multilevel     = FALSE,
      K                 = K,
      posterior_samples = posterior_samples
    )
    return(tau_result[["tau_total"]]^2)
  }

  return(.evaluate_marginalized_random_variance(
    data              = x[["data"]],
    posterior_samples = posterior_samples,
    K                 = K
  ))
}


# ---------------------------------------------------------------------------- #
# .regplot_band_data_categorical
# ---------------------------------------------------------------------------- #
#
# Categorical intervals are rendered as box-style summaries.
#
# ---------------------------------------------------------------------------- #
.regplot_band_data_categorical <- function(grid, middle, lower, upper) {

  out <- data.frame(
    x        = grid[["x"]],
    y        = middle,
    middle   = middle,
    lower    = lower,
    upper    = upper,
    group    = grid[["group"]],
    group_id = grid[["group_id"]],
    level    = grid[["level"]],
    stringsAsFactors = FALSE
  )

  return(out)
}


.regplot_normalized_precision <- function(sei) {

  if (!is.numeric(sei) || length(sei) == 0L || any(!is.finite(sei)) ||
      any(sei <= 0)) {
    stop("Regression-plot standard errors must be finite and positive.",
         call. = FALSE)
  }

  relative_sd <- min(sei) / sei
  minimum     <- min(relative_sd)
  if (minimum == 1) {
    return(numeric(length(sei)))
  }

  return(
    ((relative_sd - minimum) * (relative_sd + minimum)) /
      ((1 - minimum) * (1 + minimum))
  )
}


# ---------------------------------------------------------------------------- #
# .regplot_data
# ---------------------------------------------------------------------------- #
#
# Generate data for regression plot.
#
# @param x        brma object
# @param mod_name name of moderator variable
# @param mod_type "continuous" or "categorical"
# @param mod_data moderator data values
# @param by_info  list with grouping variable info (or NULL)
# @param pred     logical; show prediction line
# @param ci       logical; show CI bands
# @param pi       logical; show PI bands
# @param si       logical; show SI bands
# @param level    confidence level (0-100)
# @param at       values at which to evaluate predictions
# @param psize    point sizes (or NULL for auto)
# @param plim     point size range
# @param transf   transformation function
# @param xlim     x-axis limits
# @param ylim     y-axis limits
# @param xlab     x-axis label
# @param ylab     y-axis label
# @param refline       reference line position
# @param sampling_bias logical; incorporate bias into predictions
# @param max_samples   maximum posterior samples for plot summaries
# @param reference_sei numeric; reference standard error for sampling paths
# @param dots          graphical parameters
#
# @return list with plot data components
#
# ---------------------------------------------------------------------------- #
.regplot_data <- function(x, mod_name, mod_type, mod_data, by_info,
                          pred, ci, pi, si, level, at, digits, psize, plim,
                          transf, xlim, ylim, xlab, ylab, refline,
                          sampling_bias, max_samples, reference_sei, dots) {

  yi      <- .outcome_data_yi(x)
  sei_obs <- .outcome_data_sei(x)
  K       <- length(yi)
  se_rep  <- if (is.null(reference_sei)) stats::median(sei_obs) else reference_sei

  alpha <- (100 - level) / 100
  probs <- c(alpha / 2, 1 - alpha / 2)

  if (is.null(psize)) {
    precision <- .regplot_normalized_precision(sei_obs)
    psize     <- plim[1] + precision * (plim[2] - plim[1])
  } else if (length(psize) == 1) {
    psize <- rep(psize, K)
  }

  yi_plot <- yi
  if (!is.null(transf)) {
    yi_plot <- transf(yi)
  }

  prediction_grid <- .regplot_prediction_grid(
    x         = x,
    mod_name  = mod_name,
    mod_type  = mod_type,
    mod_data  = mod_data,
    by_info   = by_info,
    at        = at,
    digits    = digits
  )
  grid   <- prediction_grid[["grid"]]
  groups <- prediction_grid[["groups"]]
  n_pred <- nrow(grid)

  prediction_se <- if (sampling_bias) se_rep else 0
  newdata <- .regplot_add_dummy_outcome(
    x       = x,
    newdata = prediction_grid[["newdata"]],
    n_pred  = n_pred,
    sei     = prediction_se
  )

  posterior_samples <- .get_posterior_samples(x[["fit"]])
  selected_rows     <- .thin_sample_rows(nrow(posterior_samples), max_samples)
  if (!is.null(selected_rows)) {
    posterior_samples <- posterior_samples[selected_rows, , drop = FALSE]
  }

  pred_samples <- predict.brma(
    object             = x,
    newdata            = newdata,
    type               = "terms",
    bias_adjusted      = !sampling_bias,
    quiet              = TRUE,
    .posterior_samples = posterior_samples
  )
  pred_samples <- as.matrix(pred_samples)

  if (nrow(pred_samples) != nrow(posterior_samples)) {
    stop("Posterior sample count mismatch in regplot().", call. = FALSE)
  }

  pred_mean  <- colMeans(pred_samples)
  pred_lower <- apply(pred_samples, 2, stats::quantile, probs = probs[1], names = FALSE)
  pred_upper <- apply(pred_samples, 2, stats::quantile, probs = probs[2], names = FALSE)

  if (pi || si) {
    if (inherits(x, "brma.mv") && .is_random(x)) {
      tau_total <- .regplot_mv_random_sd_samples(
        x                 = x,
        newdata           = newdata,
        posterior_samples = posterior_samples
      )
    } else {
      is_scale      <- .is_scale(x)
      is_multilevel <- .is_multilevel(x)

      scale_data    <- NULL
      scale_formula <- NULL
      if (is_scale) {
        new_data_prepared <- .prepare_newdata(object = x, newdata = newdata, type = "estimate")
        scale_data        <- new_data_prepared[["scale"]]
        scale_formula     <- .create_fit_formula_list(data = new_data_prepared, "scale")
      }

      tau_result <- .evaluate.brma.tau(
        fit               = x[["fit"]],
        scale_data        = scale_data,
        scale_formula     = scale_formula,
        scale_priors      = x[["priors"]][["scale"]],
        is_scale          = is_scale,
        is_multilevel     = is_multilevel,
        K                 = n_pred,
        posterior_samples = posterior_samples,
        fixed_tau         = .fixed_tau_prior_value(x[["priors"]]),
        fixed_rho         = .fixed_rho_prior_value(x[["priors"]])
      )
      tau_total <- .root_sum_squares(
        tau_result[["tau_within"]],
        tau_result[["tau_between"]]
      )
    }
  }

  if (pi) {
    pi_bounds <- .regplot_mixture_interval_quantiles(
      mean_samples = pred_samples,
      sd_samples   = tau_total,
      probs        = probs
    )
    pi_lower <- pi_bounds[["lower"]]
    pi_upper <- pi_bounds[["upper"]]
  } else {
    pi_lower <- NULL
    pi_upper <- NULL
  }

  if (si) {
    sd_si <- .root_sum_squares(tau_total, se_rep)

    if (sampling_bias && .is_weightfunction(x)) {
      si_bounds <- .regplot_selection_mixture_interval_quantiles(
        x                 = x,
        mean_samples      = pred_samples,
        sd_samples        = sd_si,
        se                = se_rep,
        probs             = probs,
        posterior_samples = posterior_samples
      )
    } else {
      si_bounds <- .regplot_mixture_interval_quantiles(
        mean_samples = pred_samples,
        sd_samples   = sd_si,
        probs        = probs
      )
    }

    si_lower <- si_bounds[["lower"]]
    si_upper <- si_bounds[["upper"]]
  } else {
    si_lower <- NULL
    si_upper <- NULL
  }

  if (!is.null(transf)) {
    pred_mean  <- transf(pred_mean)
    ci_lo      <- transf(pred_lower)
    ci_hi      <- transf(pred_upper)
    pred_lower <- pmin(ci_lo, ci_hi)
    pred_upper <- pmax(ci_lo, ci_hi)

    if (pi) {
      pi_lo    <- transf(pi_lower)
      pi_hi    <- transf(pi_upper)
      pi_lower <- pmin(pi_lo, pi_hi)
      pi_upper <- pmax(pi_lo, pi_hi)
    }
    if (si) {
      si_lo    <- transf(si_lower)
      si_hi    <- transf(si_upper)
      si_lower <- pmin(si_lo, si_hi)
      si_upper <- pmax(si_lo, si_hi)
    }
  }

  if (is.null(xlim)) {
    if (mod_type == "continuous") {
      xlim <- range(pretty(range(mod_data)))
    } else {
      xlim <- c(0.5, length(levels(mod_data)) + 0.5)
    }
  }

  if (is.null(ylim)) {
    all_y <- c(yi_plot, pred_mean)
    if (ci) all_y <- c(all_y, pred_lower, pred_upper)
    if (pi) all_y <- c(all_y, pi_lower, pi_upper)
    if (si) all_y <- c(all_y, si_lower, si_upper)
    ylim <- range(pretty(range(all_y, na.rm = TRUE)))
  }

  if (is.null(xlab)) xlab <- mod_name
  if (is.null(ylab)) ylab <- "Observed Effect Size"

  observed_groups <- .regplot_observed_groups(
    by_info = by_info,
    groups  = groups,
    n       = K
  )

  df_points <- data.frame(
    x        = if (mod_type == "continuous") mod_data else as.numeric(mod_data),
    y        = yi_plot,
    size     = psize,
    group    = observed_groups[["group"]],
    group_id = observed_groups[["group_id"]],
    stringsAsFactors = FALSE
  )

  if (mod_type == "categorical") {
    df_points$level <- mod_data
  }

  df_pred <- NULL
  if (pred) {
    df_pred <- data.frame(
      x        = grid[["x"]],
      y        = pred_mean,
      group    = grid[["group"]],
      group_id = grid[["group_id"]],
      stringsAsFactors = FALSE
    )
    if (mod_type == "categorical") {
      df_pred[["level"]] <- grid[["level"]]
    }
  }

  df_ci <- NULL
  if (ci) {
    if (mod_type == "continuous") {
      df_ci <- do.call(rbind, lapply(seq_along(groups), function(i) {
        rows <- grid[["group_id"]] == i
        .regplot_band_data_continuous(
          xpred    = grid[["x"]][rows],
          lower    = pred_lower[rows],
          upper    = pred_upper[rows],
          group    = groups[i],
          group_id = i
        )
      }))
      rownames(df_ci) <- NULL
    } else {
      df_ci <- .regplot_band_data_categorical(grid, pred_mean, pred_lower, pred_upper)
    }
  }

  df_pi <- NULL
  if (pi) {
    if (mod_type == "continuous") {
      df_pi <- do.call(rbind, lapply(seq_along(groups), function(i) {
        rows <- grid[["group_id"]] == i
        .regplot_band_data_continuous(
          xpred    = grid[["x"]][rows],
          lower    = pi_lower[rows],
          upper    = pi_upper[rows],
          group    = groups[i],
          group_id = i
        )
      }))
      rownames(df_pi) <- NULL
    } else {
      df_pi <- .regplot_band_data_categorical(grid, pred_mean, pi_lower, pi_upper)
    }
  }

  df_si <- NULL
  if (si) {
    if (mod_type == "continuous") {
      df_si <- do.call(rbind, lapply(seq_along(groups), function(i) {
        rows <- grid[["group_id"]] == i
        .regplot_band_data_continuous(
          xpred    = grid[["x"]][rows],
          lower    = si_lower[rows],
          upper    = si_upper[rows],
          group    = groups[i],
          group_id = i
        )
      }))
      rownames(df_si) <- NULL
    } else {
      df_si <- .regplot_band_data_categorical(grid, pred_mean, si_lower, si_upper)
    }
  }

  df_refline <- NULL
  if (!is.null(refline)) {
    df_refline <- data.frame(y = refline)
  }

  return(list(
    points   = df_points,
    pred     = df_pred,
    ci       = df_ci,
    pi       = df_pi,
    si       = df_si,
    refline  = df_refline,
    xlim     = xlim,
    ylim     = ylim,
    xlab     = xlab,
    ylab     = ylab,
    mod_type = mod_type,
    mod_name = mod_name,
    by_name  = if (!is.null(by_info)) by_info[["name"]] else NULL,
    by_type  = if (!is.null(by_info)) by_info[["type"]] else NULL,
    groups   = groups,
    levels   = if (mod_type == "categorical") levels(mod_data) else NULL,
    sei      = se_rep
  ))
}
