#' @title Input Data Specification
#' @name data_input
#' @description
#' Shared data-input arguments used by the RoBMA fitting functions.
#'
#' Normal models use approximate effect-size estimates supplied through
#' `yi` with either `vi` or `sei`. GLMM models use the raw two-arm count
#' arguments for binomial (`measure = "OR"`) or Poisson (`measure = "IRR"`)
#' outcomes.
#'
#' @param yi a vector of effect sizes, or a formula with the effect size on the
#' left-hand side and location moderators on the right-hand side (for example
#' `yi ~ x1 + x2`). If a formula is supplied, `mods` must not be specified.
#' @param vi a vector of sampling variances. Either `vi` or `sei` must be
#' supplied for normal models.
#' @param sei a vector of standard errors. Either `vi` or `sei` must be
#' supplied for normal models.
#' @param V a known working variance-covariance matrix, or a list of block
#' variance-covariance matrices, used by `brma.mv()`.
#' @param known_v_parameterization known-`V` backend used by `brma.mv()`.
#' `"auto"` chooses an exact backend when feasible; `"latent"` uses a latent
#' `D + BB'` decomposition; `"whitened"` uses an eigen-rotated normal
#' likelihood; `"block_mvn"` uses an exact native block multivariate-normal
#' likelihood.
#' @param known_v_residual_fraction proportion of the diagonal of `V` left as
#' conditional independent residual sampling variance in the latent `D + BB'`
#' representation. Defaults to `0.10`; values are validated for all backends
#' but affect only the latent backend. When `"auto"` selects `"whitened"` or
#' `"block_mvn"`, explicitly supplied values are disregarded by the likelihood.
#' Direct `"whitened"` or `"block_mvn"` requests warn when an explicit value is
#' supplied.
#' @param weights an optional vector of positive likelihood weights. For
#' normal/effect-size models, each weight powers the estimate likelihood. For
#' constructors with GLMM raw-count input, each weight powers the paired
#' two-arm likelihood for one study.
#' @param ni an optional vector of sample sizes. Used for `measure = "GEN"` or
#' when estimating the unit information standard deviation.
#' @param mods an optional matrix, data.frame, or formula specifying
#' location moderators (meta-regressors). Formula input is evaluated in `data`
#' and supports explicit columns and standard formula operators. Precompute
#' transformed predictors as columns in `data`; inline calls such as `I()`,
#' `poly()`, spline constructors, and user functions are not supported.
#' @param scale an optional matrix, data.frame, or formula specifying
#' scale predictors for location-scale models. Formula input is evaluated in
#' `data` and follows the same explicit-column formula grammar as `mods`.
#' @param random an optional formula or list of formulas specifying
#' BayesTools random-effect terms for `brma.mv()`. Use
#' [random-effect formula structure tags][random_effect_formula_tags] such as
#' `diag()`, `us()`/`un()`, `cs()`, `hcs()`, `ar1()`/`ar()`, `har()`, or
#' `car()` inside the formula. Plain `(expr | group)` syntax uses an
#' unstructured random-coefficient block; `expr || group` uses a diagonal
#' block. The coefficient-formula and structured-index tag families have
#' different left-side semantics; see the linked documentation. A bare formula
#' or unnamed one-entry list has no redundant component prefix. An explicitly
#' named one-entry list retains its name; lists with two or more entries
#' generate missing names as `component 1`, `component 2`, and so on.
#' @param cluster an optional vector of cluster identifiers for multilevel
#' meta-analysis.
#' @param data an optional data frame containing the variables.
#' @param slab an optional vector of study labels.
#' @param subset an optional logical or numeric vector specifying a subset of
#' data to be used.
#' @param measure a character string specifying the effect size measure.
#' Normal/effect-size constructors require an explicit value and support
#' `"SMD"`, `"ZCOR"`, `"RR"`, `"OR"`, `"HR"`, `"RD"`, `"IRR"`, and `"GEN"`.
#' Use `"GEN"` only for general effect sizes without a known unit information
#' standard deviation. GLMM raw-count constructors support only `"OR"` and
#' `"IRR"` and default to `"OR"`.
#' @param effect_direction direction used by publication-bias adjustments.
#' `"positive"` assumes statistically significant positive estimates are more
#' likely to be selected; `"negative"` mirrors the selection direction;
#' `"detect"` infers the direction from the fitted data.
#' @param ai a vector of the number of events in the treatment or experimental
#' group for binomial GLMM models.
#' @param bi a vector of the number of non-events in the treatment or
#' experimental group for binomial GLMM models.
#' @param ci a vector of the number of events in the control group for binomial
#' GLMM models.
#' @param di a vector of the number of non-events in the control group for
#' binomial GLMM models.
#' @param n1i a vector of the sample size in the treatment or experimental
#' group. If omitted for binomial GLMMs, it is computed as `ai + bi`.
#' @param n2i a vector of the sample size in the control group. If omitted for
#' binomial GLMMs, it is computed as `ci + di`.
#' @param x1i a vector of the number of events in the treatment/experimental group
#' (for Poisson data).
#' @param x2i a vector of the number of events in the control group (for Poisson data).
#' @param t1i a vector of the person-time in the treatment/experimental group.
#' @param t2i a vector of the person-time in the control group.
#'
NULL



# Internal helper function to extract a variable from a data frame or environment
# Supports three input formats:
#   1. Direct vector: yi = c(0.1, 0.2, 0.3)
#   2. Variable name (unquoted): yi = effect_size (where effect_size is a column in data)
#   3. String name (quoted): yi = "effect_size" (where "effect_size" is a column name in data)
#
# @param mf The matched call from the parent function
# @param data The data frame to search in (can be NULL)
# @param enclos The enclosing environment for evaluation
# @param name The name of the argument (for error messages)
# @param allow_NULL Logical; if TRUE, NULL values are allowed
# @return The extracted vector or NULL
.get_variable <- function(mf, data, enclos, name, allow_NULL = TRUE) {

  # Check if argument is in the call
  arg_index <- match(name, names(mf))

  # If argument was not specified in the call, return NULL
  if (is.na(arg_index)) {
    return(NULL)
  }

  # Get the unevaluated expression from the matched call
  mf_x <- mf[[arg_index]]

  # If the argument value itself is NULL, return NULL
  if (is.null(mf_x)) {
    return(NULL)
  }

  # First, try to evaluate the expression directly in the data frame (if provided)
  if (!is.null(data) && is.data.frame(data)) {
    out <- try(eval(mf_x, data, enclos), silent = TRUE)
  } else {
    out <- try(eval(mf_x, enclos), silent = TRUE)
  }

  # If evaluation succeeded and result is a valid vector, return it
  if (!inherits(out, "try-error") && !is.function(out)) {

    # If result is a single string matching a column name in data, extract that column
    if (is.character(out) && length(out) == 1 && !is.null(data) && is.data.frame(data) && out %in% names(data)) {
      out <- data[[out]]
    }

    # Strip attributes (e.g., from metafor::escalc output) from atomic vectors only
    # Do not strip from formulas, data.frames, lists, or other complex objects
    if (is.atomic(out) && is.null(dim(out)) && !inherits(out, "formula")) {
      out <- as.vector(out)
    }

    return(out)
  }

  # If direct evaluation failed, check if it's a string referring to a column
  if (is.character(mf_x) && length(mf_x) == 1 && !is.null(data) && is.data.frame(data) && mf_x %in% names(data)) {
    return(as.vector(data[[mf_x]]))
  }

  # If still not found, report error. Missing arguments and explicit NULL
  # values were handled above.
  stop(paste0("Cannot find the object/variable ('", deparse(mf_x),
              "') specified for the '", name, "' argument."), call. = FALSE)
}


# Internal function to extract common optional variables shared by all outcome handlers
# Handles weights, cluster, and slab extraction and validation
#
# @param .call Matched call from the calling function
# @param data The data frame (can be NULL)
# @param .envir The enclosing environment
# @param k The expected number of observations
# @param primary_var The name of the primary variable (for error messages, e.g., "yi" or "ai")
#
# @return A list with:
#   - weights: numeric vector or NULL
#   - cluster: vector or NULL
#   - slab: character vector or NULL
#   - slab_provided: logical
#   - cluster_provided: logical
.check_and_list_data.optional_vars <- function(.call, data, .envir, k, primary_var) {

  # Extract optional variables
  weights   <- .get_variable(.call, data, .envir, "weights",   allow_NULL = TRUE)
  cluster   <- .get_variable(.call, data, .envir, "cluster",   allow_NULL = TRUE)
  slab      <- .get_variable(.call, data, .envir, "slab",      allow_NULL = TRUE)

  # Track which optional fields were provided
  weights_provided   <- !is.null(weights)
  slab_provided      <- !is.null(slab)
  cluster_provided   <- !is.null(cluster)

  # Validate weights
  if (!is.null(weights) && anyNA(weights))
    stop("The 'weights' argument must not contain missing values.", call. = FALSE)
  if (!is.null(weights))
    BayesTools::check_real(weights, "weights", check_length = k, allow_NULL = TRUE, allow_NA = FALSE, lower = 0, allow_bound = FALSE)

  # Validate cluster
  if (!is.null(cluster) && length(cluster) != k)
    stop(paste0("The 'cluster' argument must have length ", k, " (same as '", primary_var, "')."), call. = FALSE)
  if (!is.null(cluster) && anyNA(cluster))
    stop("The 'cluster' argument must not contain missing values.", call. = FALSE)

  # Validate slab (study labels)
  if (!is.null(slab) && length(slab) != k)
    stop(paste0("The 'slab' argument must have length ", k, " (same as '", primary_var, "')."), call. = FALSE)

  return(list(
    weights            = weights,
    cluster            = cluster,
    slab               = slab,
    weights_provided   = weights_provided,
    slab_provided      = slab_provided,
    cluster_provided   = cluster_provided
  ))
}


# Internal function to check and list input data for RoBMA
# This is the main entry point for processing input data. It dispatches to
# the appropriate outcome handler based on the class argument.
#
# @param .call Matched call from the calling function (required for NSE)
# @param .envir The environment where the calling function was invoked (required for NSE)
# @param class The model class: "norm" for normal likelihood (yi, sei) or
#              "glmm" for generalized linear mixed models (ai, bi, ci, di, etc.)
# @param measure The effect size measure (used for "glmm" class to determine
#                outcome type: "OR" for binomial, "IRR" for Poisson)
# @param standardize_continuous_predictors Whether to standardize continuous predictors
# @param skip_validation Whether to skip predictor validation checks (variability/levels).
#        Should be TRUE when processing newdata for prediction. (validation done in BayesTools::JAGS_formula())
#
# @return A list with components:
#   - outcome: data.frame with outcome variables (columns depend on class)
#   - mods: moderator information (data.frame or NULL)
#   - scale: scale information (data.frame or NULL)
.check_and_list_data <- function(.call, .envir, class = "norm", measure,
                                 set_contrast_factor_predictors, standardize_continuous_predictors,
                                 effect_direction = "positive", skip_validation = FALSE,
                                 allow_na_drop = TRUE,
                                 random_effects_metadata = NULL,
                                 random_group_covariance = NULL,
                                 known_v_parameterization = "auto",
                                 known_v_residual_fraction = NULL,
                                 known_v_residual_fraction_specified = FALSE) {

  # check additional input
  .check_measure(measure, class = class)
  BayesTools::check_bool(standardize_continuous_predictors, "standardize_continuous_predictors", allow_NA = FALSE)
  BayesTools::check_char(set_contrast_factor_predictors, "set_contrast_factor_predictors", allow_values = c("treatment", "meandif", "orthonormal", "independent"), allow_NA = FALSE)
  if (missing(effect_direction)) {
    effect_direction <- "positive"
  }
  BayesTools::check_char(effect_direction, "effect_direction", allow_values = c("positive", "negative", "detect"))
  BayesTools::check_bool(skip_validation, "skip_validation")
  BayesTools::check_bool(allow_na_drop, "allow_na_drop")
  BayesTools::check_bool(known_v_residual_fraction_specified, "known_v_residual_fraction_specified")
  if (is.null(known_v_parameterization)) {
    known_v_parameterization <- "auto"
  }
  BayesTools::check_char(
    known_v_parameterization,
    "known_v_parameterization",
    allow_values = c("auto", "latent", "whitened", "block_mvn")
  )

  ### Extract the data argument first - other variables may reference columns within it
  data <- .get_variable(.call, NULL, .envir, "data", allow_NULL = TRUE)

  ### Step 1: Extract and validate outcome variables (dispatch based on class)
  outcome_result <- switch(
    class,
    "norm" = .check_and_list_data.outcome.norm(.call, data, .envir, effect_direction, skip_validation),
    "mv"   = .check_and_list_data.outcome.mv(.call, data, .envir, effect_direction, skip_validation),
    "glmm" = switch(
      measure,
      "OR"  = .check_and_list_data.outcome.bin(.call, data, .envir, skip_validation),
      "IRR" = .check_and_list_data.outcome.pois(.call, data, .envir, skip_validation)
    )
  )

  data_outcome       <- outcome_result$data_outcome
  k                  <- outcome_result$k
  mods_from_yi       <- outcome_result$mods_from_yi
  formula_yi         <- outcome_result$formula_yi
  weights_provided   <- outcome_result$weights_provided
  slab_provided      <- outcome_result$slab_provided
  cluster_provided   <- outcome_result$cluster_provided
  na_check_cols      <- outcome_result$na_check_cols
  outcome_type       <- outcome_result$outcome_type
  effect_direction   <- outcome_result$effect_direction
  known_V_input      <- outcome_result$known_V_input
  known_V_hidden     <- outcome_result$known_V_hidden

  ### Step 2: Extract moderator variables (mods and scale)
  if (!is.null(mods_from_yi)) {
    # Mods already extracted from yi formula
    data_mods <- as.data.frame(mods_from_yi)
    attr(data_mods, "formula")    <- formula_yi[-2]  # RHS only
    attr(data_mods, "formula_yi") <- formula_yi
  } else {
    data_mods <- .check_and_list_data.predictors(
      .call  = .call,
      data   = data,
      .envir = .envir,
      name   = "mods",
      k      = k
    )
  }

  data_scale <- .check_and_list_data.scale(
    .call  = .call,
    data   = data,
    .envir = .envir,
    k      = k
  )
  data_random <- .check_and_list_data.random(
    .call                   = .call,
    data                    = data,
    .envir                  = .envir,
    k                       = k,
    random_effects_metadata = random_effects_metadata,
    random_group_covariance = random_group_covariance
  )

  data_mods  <- .check_and_list_data.coerce_character_predictors(data_mods)
  data_scale <- .check_and_list_data.scale_coerce_character_predictors(data_scale)
  if (!is.null(data_random[["data"]])) {
    data_random[["data"]] <- .check_and_list_data.coerce_character_predictors(data_random[["data"]])
  }

  ### Step 3: Apply subset (before NA handling)
  subset    <- .get_variable(.call, data, .envir, "subset", allow_NULL = TRUE)
  keep_rows <- rep(TRUE, k)

  if (!is.null(subset)) {
    # Validate and convert subset to logical
    subset <- .check_and_list_data.validate_subset(subset, k)
    keep_rows <- keep_rows & subset

    # Apply subsetting to all data frames
    data_outcome <- .check_and_list_data.subset(data_outcome, subset)
    data_mods    <- .check_and_list_data.subset(data_mods, subset)
    data_scale   <- .check_and_list_data.scale_subset(data_scale, subset)
    if (!is.null(data_random[["data"]])) {
      data_random[["data"]] <- .check_and_list_data.subset(data_random[["data"]], subset)
    }
  }

  ### Step 4: Handle NA dropping
  # Build list of data frames to check for NAs
  # Use only the columns specified by the outcome handler for NA checking
  data_list_for_na <- list(
    outcome_required = data_outcome[, na_check_cols, drop = FALSE]
  )
  if (!is.null(data_mods)) {
    data_list_for_na$mods <- data_mods
  }
  data_list_for_na <- c(
    data_list_for_na,
    .check_and_list_data.scale_na_frames(data_scale)
  )
  if (!is.null(data_random[["data"]])) {
    data_list_for_na$random <- data_random[["data"]]
  }

  # Get rows with NAs
  na_rows <- .check_and_list_data.check_na(data_list_for_na)

  # Drop NA rows if any
  n_dropped <- sum(na_rows)
  if (n_dropped > 0) {
    if (!allow_na_drop) {
      .check_and_list_data.stop_na(data_list_for_na)
    }

    warning(paste0(n_dropped, " observation(s) removed due to missing values."), call. = FALSE, immediate. = TRUE)

    keep_after_subset <- !na_rows
    current_rows <- which(keep_rows)
    keep_rows[current_rows] <- keep_after_subset

    data_outcome <- .check_and_list_data.subset(data_outcome, keep_after_subset)
    data_mods    <- .check_and_list_data.subset(data_mods, keep_after_subset)
    data_scale   <- .check_and_list_data.scale_subset(data_scale, keep_after_subset)
    if (!is.null(data_random[["data"]])) {
      data_random[["data"]] <- .check_and_list_data.subset(data_random[["data"]], keep_after_subset)
    }
  }

  ### Step 5: Final validation and processing
  k_final <- nrow(data_outcome)

  if (k_final == 0) {
    stop("No observations remaining after removing missing values.", call. = FALSE)
  }

  if (outcome_type == "norm" && effect_direction == "detect") {
    effect_direction <- if (stats::median(data_outcome[["yi"]], na.rm = TRUE) >= 0) "positive" else "negative"
  }

  # Validate predictor variables (after subsetting and NA dropping)
  # Skip validation when processing newdata (skip_validation = TRUE)
  .check_and_list_data.validate_predictors(data_mods, "mods", skip_validation)
  .check_and_list_data.scale_validate_predictors(data_scale, skip_validation)
  data_scale <- .check_and_list_data.validate_scale_random(data_scale, data_random)

  # Generate default study labels if not provided (after NA dropping)
  if (!slab_provided) {
    data_outcome$slab <- paste0("Study ", seq_len(k_final))
  }

  # Convert cluster to numeric factor for use as indices in fitting (after NA dropping)
  if (cluster_provided) {
    if (anyNA(data_outcome[["cluster"]])) {
      stop("The 'cluster' argument must not contain missing values.", call. = FALSE)
    }

    cluster_factor <- factor(data_outcome[["cluster"]], levels = unique(data_outcome[["cluster"]]))
    data_outcome[["cluster_label"]] <- as.character(cluster_factor)
    data_outcome[["cluster"]]       <- as.integer(cluster_factor)
  }

  ### Create output object
  known_V <- NULL
  if (!is.null(known_V_input)) {
    .check_and_list_data.mv_validate_hidden_inputs(
      hidden    = known_V_hidden,
      keep_rows = keep_rows
    )
    known_V <- .known_v_prepare(
      V                                   = known_V_input,
      keep_rows                           = keep_rows,
      known_v_parameterization            = known_v_parameterization,
      known_v_residual_fraction           = known_v_residual_fraction,
      known_v_residual_fraction_specified = known_v_residual_fraction_specified,
      known_v_is_scale                    = !is.null(data_scale)
    )
    data_outcome[["sei"]] <- sqrt(.known_v_diagonal(known_V))
  }

  data_list <- list(
    outcome  = data_outcome,
    mods     = data_mods,
    scale    = data_scale,
    location = .check_and_list_data.location(
      data_mods   = data_mods,
      data_random = data_random
    )
  )
  data_list <- .check_and_list_data.detach_formula_environments(data_list)

  class(data_list) <- "RoBMA_data"
  attr(data_list, "outcome_type")                       <- outcome_type
  attr(data_list, "measure")                            <- measure
  attr(data_list, "n_dropped")                          <- n_dropped
  attr(data_list, "k_final")                            <- k_final
  attr(data_list, "mods")                               <- !is.null(data_mods)
  attr(data_list, "scale")                              <- !is.null(data_scale)
  attr(data_list, "random")                             <- !is.null(data_random[["formula"]])
  attr(data_list, "weights")                            <- weights_provided
  attr(data_list, "known_V")                            <- !is.null(known_V)
  attr(data_list, "known_V_data")                       <- known_V
  attr(data_list, "slab")                               <- slab_provided
  attr(data_list, "cluster")                            <- cluster_provided
  attr(data_list, "standardize_continuous_predictors")  <- standardize_continuous_predictors
  attr(data_list, "set_contrast_factor_predictors")     <- set_contrast_factor_predictors
  attr(data_list, "effect_direction")                   <- effect_direction
  return(data_list)
}


.equal_within_double_roundoff <- function(x, y) {

  difference <- abs(x - y)
  scale      <- pmax(abs(x), abs(y))

  return(difference == 0 | difference <= .Machine$double.eps * scale)
}


.sampling_variance_matches_se <- function(vi, sei) {

  return(.equal_within_double_roundoff(sqrt(vi), sei))
}


# Fixed and scale formulas are later compiled by BayesTools::JAGS_formula().
# Validate its deliberately narrow replay grammar before model.frame() can
# materialize an input design that the fitting backend cannot reproduce.
.check_and_list_data.validate_fixed_formula <- function(formula, name) {

  if (!inherits(formula, "formula")) {
    stop("Internal error: 'formula' must be a formula.", call. = FALSE)
  }

  validate_expression <- function(expression) {

    if (is.symbol(expression)) {
      if (identical(as.character(expression), ".")) {
        stop(
          "Unsupported term '.' in the '", name,
          "' formula. List each data column explicitly.",
          call. = FALSE
        )
      }
      return(invisible(TRUE))
    }

    if (is.numeric(expression) && length(expression) == 1L &&
        is.finite(expression) && expression %in% c(0, 1)) {
      return(invisible(TRUE))
    }

    if (is.call(expression)) {
      call_name <- if (is.symbol(expression[[1L]])) {
        as.character(expression[[1L]])
      } else {
        ""
      }

      if (call_name %in% c("+", "-", "*", ":", "/", "^", "(")) {
        for (argument in as.list(expression)[-1L]) {
          validate_expression(argument)
        }
        return(invisible(TRUE))
      }

      expression_label <- paste(deparse(expression), collapse = " ")
      stop(
        "Unsupported call '", expression_label, "' in the '", name,
        "' formula. Precompute transformed predictors as columns in 'data' ",
        "and reference those columns by name.",
        call. = FALSE
      )
    }

    stop(
      "Unsupported expression '", paste(deparse(expression), collapse = " "),
      "' in the '", name, "' formula. Use literal data-column names.",
      call. = FALSE
    )
  }

  rhs_index <- if (length(formula) == 3L) 3L else 2L
  validate_expression(formula[[rhs_index]])
  invisible(TRUE)
}


# Internal function to extract and validate outcome variables for normal likelihood models
# Handles yi, vi/sei, ni, weights, cluster, slab, and yi as formula
#
# @param .call Matched call from the calling function
# @param data The data frame (can be NULL)
# @param .envir The enclosing environment
# @param effect_direction The effect direction ("positive", "negative", or "detect")
# @param skip_validation Whether to skip strict validation checks (for newdata).
#        When TRUE, allows sei/vi/ni = 0 (useful for prediction at zero SE).
#
# @return A list with:
#   - data_outcome: data.frame with yi, sei, ni, cluster, slab, weights
#   - k: number of observations
#   - mods_from_yi: moderators extracted from yi formula (or NULL)
#   - formula_yi: the original yi formula (or NULL)
#   - slab_provided: logical, whether slab was provided by user
#   - cluster_provided: logical, whether cluster was provided by user
#   - na_check_cols: character vector of column names to check for NAs
.check_and_list_data.outcome.norm <- function(.call, data, .envir, effect_direction, skip_validation = FALSE) {

  # Extract yi (may be a formula like yi ~ mod1 + mod2)
  yi <- .get_variable(.call, data, .envir, "yi", allow_NULL = FALSE)

  # Handle yi as formula (e.g., yi ~ mod1 + mod2)
  formula_yi   <- NULL
  mods_from_yi <- NULL

  if (inherits(yi, "formula")) {

    formula_yi <- yi
    .check_and_list_data.validate_fixed_formula(yi, "yi")

    # Check that mods is not also specified (would be ambiguous)
    mods_check <- .get_variable(.call, data, .envir, "mods", allow_NULL = TRUE)
    if (!is.null(mods_check)) {
      stop("Cannot specify 'mods' when 'yi' is a formula. Use either 'yi ~ mod1 + mod2' or 'yi = effect, mods = ~ mod1 + mod2', but not both.", call. = FALSE)
    }

    # Extract the model frame from the formula
    na_act <- getOption("na.action")
    options(na.action = "na.pass")
    on.exit(options(na.action = na_act), add = TRUE)

    full_mf <- stats::model.frame(yi, data = data)

    # Extract response (LHS)
    yi <- stats::model.response(full_mf)
    names(yi) <- NULL

    # Extract predictors (RHS)
    if (ncol(full_mf) > 1) {
      mods_from_yi <- full_mf[, -1, drop = FALSE]
    }
  }

  # Extract variance/standard error
  vi  <- .get_variable(.call, data, .envir, "vi",  allow_NULL = TRUE)
  sei <- .get_variable(.call, data, .envir, "sei", allow_NULL = TRUE)

  # Extract ni (sample sizes) - specific to normal likelihood models
  ni <- .get_variable(.call, data, .envir, "ni", allow_NULL = TRUE)

  ### Input validation

  # Validate yi
  BayesTools::check_real(yi, "yi", check_length = 0, allow_NULL = FALSE, allow_NA = TRUE)
  if (all(is.na(yi)))
    stop("The 'yi' argument must contain at least one non-NA value.", call. = FALSE)

  k <- length(yi)

  # Validate vi (sampling variance)
  # When skip_validation = TRUE (newdata), allow vi = 0 for prediction at zero variance
  if (!is.null(vi))
    BayesTools::check_real(vi, "vi", check_length = k, allow_NULL = TRUE, allow_NA = TRUE, lower = 0, allow_bound = skip_validation)

  # Validate sei (standard error)
  # When skip_validation = TRUE (newdata), allow sei = 0 for prediction at zero SE
  if (!is.null(sei))
    BayesTools::check_real(sei, "sei", check_length = k, allow_NULL = TRUE, allow_NA = TRUE, lower = 0, allow_bound = skip_validation)

  # Convert vi to sei or vice versa
  if (is.null(sei) && !is.null(vi)) {
    sei <- sqrt(vi)
  } else if (is.null(vi) && !is.null(sei)) {
    vi <- sei^2
  } else if (is.null(vi) && is.null(sei)) {
    stop("Either 'vi' (variance) or 'sei' (standard error) must be provided.", call. = FALSE)
  } else {
    # Both provided - check consistency
    if (any(!.sampling_variance_matches_se(vi, sei), na.rm = TRUE))
      stop("The provided 'vi' and 'sei' values are inconsistent.", call. = FALSE)
  }

  represented_vi  <- sei^2
  unrepresentable <- !is.na(sei) & sei > 0 &
    (!is.finite(represented_vi) | represented_vi == 0)
  if (any(unrepresentable)) {
    stop(
      "Positive 'sei' values must have positive finite squared sampling variances.",
      call. = FALSE
    )
  }

  # Validate ni (sample sizes)
  # When skip_validation = TRUE (newdata), allow ni = 0
  if (!is.null(ni))
    BayesTools::check_real(ni, "ni", check_length = k, allow_NULL = TRUE, allow_NA = TRUE, lower = 0, allow_bound = skip_validation)

  # Extract and validate common optional variables (weights, cluster, slab)
  optional <- .check_and_list_data.optional_vars(.call, data, .envir, k, "yi")

  ### Construct output data frame
  data_outcome <- data.frame(
    yi            = yi,
    sei           = sei,
    ni            = if (!is.null(ni))                  ni               else rep(NA_integer_, k),
    cluster       = if (!is.null(optional$cluster))    optional$cluster else rep(NA_character_, k),
    cluster_label = if (!is.null(optional$cluster))    as.character(optional$cluster) else rep(NA_character_, k),
    slab          = if (!is.null(optional$slab))       optional$slab    else rep(NA_character_, k),
    weights       = if (!is.null(optional$weights))    optional$weights else rep(NA, k),
    stringsAsFactors = FALSE
  )

  return(list(
    data_outcome       = data_outcome,
    k                  = k,
    mods_from_yi       = mods_from_yi,
    formula_yi         = formula_yi,
    weights_provided   = optional$weights_provided,
    slab_provided      = optional$slab_provided,
    cluster_provided   = optional$cluster_provided,
    na_check_cols      = c("yi", "sei"),  # Only check these columns for NAs
    effect_direction   = effect_direction,
    outcome_type       = "norm"
  ))
}


# Internal function to extract and validate outcome variables for normal models
# with a known working variance-covariance matrix.
.check_and_list_data.outcome.mv <- function(.call, data, .envir, effect_direction, skip_validation = FALSE) {

  yi <- .get_variable(.call, data, .envir, "yi", allow_NULL = FALSE)

  formula_yi   <- NULL
  mods_from_yi <- NULL

  if (inherits(yi, "formula")) {

    formula_yi <- yi
    .check_and_list_data.validate_fixed_formula(yi, "yi")

    mods_check <- .get_variable(.call, data, .envir, "mods", allow_NULL = TRUE)
    if (!is.null(mods_check)) {
      stop("Cannot specify 'mods' when 'yi' is a formula. Use either 'yi ~ mod1 + mod2' or 'yi = effect, mods = ~ mod1 + mod2', but not both.", call. = FALSE)
    }

    na_act <- getOption("na.action")
    options(na.action = "na.pass")
    on.exit(options(na.action = na_act), add = TRUE)

    full_mf <- stats::model.frame(yi, data = data)

    yi <- stats::model.response(full_mf)
    names(yi) <- NULL

    if (ncol(full_mf) > 1) {
      mods_from_yi <- full_mf[, -1, drop = FALSE]
    }
  }

  V   <- .get_variable(.call, data, .envir, "V", allow_NULL = TRUE)
  vi  <- .get_variable(.call, data, .envir, "vi", allow_NULL = TRUE)
  sei <- .get_variable(.call, data, .envir, "sei", allow_NULL = TRUE)
  ni  <- .get_variable(.call, data, .envir, "ni", allow_NULL = TRUE)

  BayesTools::check_real(yi, "yi", check_length = 0, allow_NULL = FALSE, allow_NA = TRUE)
  if (all(is.na(yi)))
    stop("The 'yi' argument must contain at least one non-NA value.", call. = FALSE)

  k             <- length(yi)
  known_V_input <- .check_and_list_data.mv_known_v_input(
    V   = V,
    vi  = vi,
    sei = sei,
    k   = k
  )
  if (.known_v_input_nrow(known_V_input[["V"]]) != k) {
    stop("The dimensions of 'V' must match the length of 'yi'.", call. = FALSE)
  }
  missing_for_na <- known_V_input[["missing_for_na"]]
  sei            <- numeric(k)
  sei[missing_for_na] <- NA_real_

  if (!is.null(ni))
    BayesTools::check_real(ni, "ni", check_length = k, allow_NULL = TRUE, allow_NA = TRUE, lower = 0, allow_bound = skip_validation)

  optional <- .check_and_list_data.optional_vars(.call, data, .envir, k, "yi")
  if (optional$weights_provided) {
    stop("'weights' are not supported in brma.mv().", call. = FALSE)
  }
  if (optional$cluster_provided) {
    stop(
      "'cluster' is not supported in brma.mv(); use the dedicated ",
      "'random' argument for multilevel structures.",
      call. = FALSE
    )
  }

  data_outcome <- data.frame(
    yi            = yi,
    sei           = sei,
    ni            = if (!is.null(ni))                  ni               else rep(NA_integer_, k),
    cluster       = rep(NA_character_, k),
    cluster_label = rep(NA_character_, k),
    slab          = if (!is.null(optional$slab))       optional$slab    else rep(NA_character_, k),
    weights       = rep(NA, k),
    stringsAsFactors = FALSE
  )

  return(list(
    data_outcome       = data_outcome,
    k                  = k,
    mods_from_yi       = mods_from_yi,
    formula_yi         = formula_yi,
    weights_provided   = FALSE,
    slab_provided      = optional$slab_provided,
    cluster_provided   = FALSE,
    na_check_cols      = c("yi", "sei"),
    effect_direction   = effect_direction,
    outcome_type       = "norm",
    known_V_input      = known_V_input[["V"]],
    known_V_hidden     = known_V_input[["hidden"]]
  ))
}


.check_and_list_data.mv_known_v_input <- function(V, vi, sei, k) {

  has_V   <- !is.null(V)
  has_vi  <- !is.null(vi)
  has_sei <- !is.null(sei)

  if (has_V && (has_vi || has_sei)) {
    stop("Use only one of 'V' and 'vi'/'sei' inputs in brma.mv().",
         call. = FALSE)
  }
  if (!has_V && !has_vi && !has_sei) {
    stop("For brma.mv(), provide 'V' or diagonal input 'vi'/'sei'.",
         call. = FALSE)
  }
  if (has_V) {
    return(list(
      V              = V,
      missing_for_na = is.na(.known_v_input_diagonal(V)),
      hidden         = NULL
    ))
  }

  if (has_vi) {
    .check_and_list_data.mv_hidden_input_structure(vi, "vi", k)
  }
  if (has_sei) {
    .check_and_list_data.mv_hidden_input_structure(sei, "sei", k)
  }

  missing_for_na <- rep(FALSE, k)
  if (has_vi) {
    missing_for_na <- missing_for_na | is.na(vi)
  }
  if (has_sei) {
    missing_for_na <- missing_for_na | is.na(sei)
  }

  return(list(
    V              = if (has_vi) vi else sei^2,
    missing_for_na = missing_for_na,
    hidden         = list(
      vi  = if (has_vi)  vi  else NULL,
      sei = if (has_sei) sei else NULL
    )
  ))
}


.check_and_list_data.mv_hidden_input_structure <- function(x, name, k) {

  if (!is.numeric(x) || !is.null(dim(x)) || length(x) != k) {
    stop(
      "brma.mv() input '", name,
      "' must be a numeric vector with the same length as 'yi'.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}


# Validate hidden diagonal inputs after row selection.
.check_and_list_data.mv_validate_hidden_inputs <- function(hidden, keep_rows) {

  if (is.null(hidden)) {
    return(invisible(TRUE))
  }
  if (!is.logical(keep_rows) || anyNA(keep_rows)) {
    stop("Internal error: invalid hidden known-V row selector.", call. = FALSE)
  }

  vi  <- hidden[["vi"]]
  sei <- hidden[["sei"]]
  if (!is.null(vi)) {
    vi <- vi[keep_rows]
    if (any(!is.finite(vi)) || any(vi <= 0)) {
      stop("brma.mv() input 'vi' must contain positive finite values.",
           call. = FALSE)
    }
  }
  if (!is.null(sei)) {
    sei <- sei[keep_rows]
    if (any(!is.finite(sei)) || any(sei <= 0)) {
      stop("brma.mv() input 'sei' must contain positive finite values.",
           call. = FALSE)
    }
  }
  if (!is.null(vi) && !is.null(sei)) {
    if (any(!.sampling_variance_matches_se(vi, sei))) {
      stop("brma.mv() inputs 'vi' and 'sei' must be consistent.",
           call. = FALSE)
    }
  }
  if (!is.null(sei)) {
    represented_vi <- sei^2
    if (any(!is.finite(represented_vi)) || any(represented_vi <= 0)) {
      stop(
        "brma.mv() input 'sei' must have positive finite squared ",
        "sampling variances.",
        call. = FALSE
      )
    }
  }

  return(invisible(TRUE))
}


# Internal helper to validate and canonicalize binomial arm counts.
.canonicalize_binomial_counts <- function(ai, bi = NULL, ci = NULL, di = NULL,
                                          n1i = NULL, n2i = NULL,
                                          skip_validation = FALSE) {

  if (is.null(ai)) {
    stop(
      "For GLMM models, provide either (ai, bi, ci, di) or ",
      "(ai, ci, n1i, n2i).",
      call. = FALSE
    )
  }

  k          <- length(ai)
  count_args <- list(ai = ai, bi = bi, ci = ci, di = di, n1i = n1i, n2i = n2i)
  for (arg_name in names(count_args)) {
    if (!is.null(count_args[[arg_name]])) {
      BayesTools::check_int(
        count_args[[arg_name]], arg_name,
        check_length = k, allow_NULL = FALSE, allow_NA = TRUE, lower = 0
      )
    }
  }

  has_cells  <- !is.null(ci) && !is.null(bi) && !is.null(di)
  has_totals <- !is.null(ci) && !is.null(n1i) && !is.null(n2i)

  if (!has_cells && !has_totals) {
    stop(
      "For GLMM models, provide either (ai, bi, ci, di) or ",
      "(ai, ci, n1i, n2i).",
      call. = FALSE
    )
  }

  if (!is.null(n1i) && any(ai > n1i, na.rm = TRUE)) {
    stop("Invalid data: some values of 'bi' (= n1i - ai) are negative.", call. = FALSE)
  }
  if (!is.null(n2i) && any(ci > n2i, na.rm = TRUE)) {
    stop("Invalid data: some values of 'di' (= n2i - ci) are negative.", call. = FALSE)
  }

  if (!is.null(bi) && !is.null(n1i) && any(n1i != ai + bi, na.rm = TRUE)) {
    stop("The provided 'n1i' values must equal 'ai + bi' when both are supplied.", call. = FALSE)
  }
  if (!is.null(di) && !is.null(n2i) && any(n2i != ci + di, na.rm = TRUE)) {
    stop("The provided 'n2i' values must equal 'ci + di' when both are supplied.", call. = FALSE)
  }

  if (!is.null(bi)) {
    n1i <- ai + bi
  }
  if (!is.null(di)) {
    n2i <- ci + di
  }

  if (!skip_validation) {
    if (any(n1i <= 0, na.rm = TRUE)) {
      stop("Invalid data: 'n1i' must contain positive arm totals.", call. = FALSE)
    }
    if (any(n2i <= 0, na.rm = TRUE)) {
      stop("Invalid data: 'n2i' must contain positive arm totals.", call. = FALSE)
    }
  }

  return(list(
    ai  = ai,
    ci  = ci,
    n1i = n1i,
    n2i = n2i,
    k   = k
  ))
}


# Internal function to extract and validate outcome variables for binomial GLMM models
# Handles ai, bi, ci, di, n1i, n2i, weights, cluster, slab (for measure = "OR")
#
# @param .call Matched call from the calling function
# @param data The data frame (can be NULL)
# @param .envir The enclosing environment
#
# @return A list with:
#   - data_outcome: data.frame with ai, ci, n1i, n2i, cluster, slab, weights
#   - k: number of observations
#   - mods_from_yi: NULL (GLMM does not support formula syntax for outcome)
#   - formula_yi: NULL
#   - slab_provided: logical, whether slab was provided by user
#   - cluster_provided: logical, whether cluster was provided by user
#   - na_check_cols: character vector of column names to check for NAs
# @param skip_validation Whether to skip strict validation checks (for newdata).
.check_and_list_data.outcome.bin <- function(.call, data, .envir, skip_validation = FALSE) {

  # Extract cell counts for 2x2 tables
  ai  <- .get_variable(.call, data, .envir, "ai",  allow_NULL = TRUE)
  bi  <- .get_variable(.call, data, .envir, "bi",  allow_NULL = TRUE)
  ci  <- .get_variable(.call, data, .envir, "ci",  allow_NULL = TRUE)
  di  <- .get_variable(.call, data, .envir, "di",  allow_NULL = TRUE)
  n1i <- .get_variable(.call, data, .envir, "n1i", allow_NULL = TRUE)
  n2i <- .get_variable(.call, data, .envir, "n2i", allow_NULL = TRUE)

  counts <- .canonicalize_binomial_counts(
    ai              = ai,
    bi              = bi,
    ci              = ci,
    di              = di,
    n1i             = n1i,
    n2i             = n2i,
    skip_validation = skip_validation
  )
  ai  <- counts[["ai"]]
  ci  <- counts[["ci"]]
  n1i <- counts[["n1i"]]
  n2i <- counts[["n2i"]]
  k   <- counts[["k"]]

  # Extract and validate common optional variables (weights, cluster, slab)
  optional <- .check_and_list_data.optional_vars(.call, data, .envir, k, "ai")

  ### Construct output data frame
  data_outcome <- data.frame(
    ai            = ai,
    ci            = ci,
    n1i           = n1i,
    n2i           = n2i,
    cluster       = if (!is.null(optional$cluster)) optional$cluster else rep(NA_character_, k),
    cluster_label = if (!is.null(optional$cluster)) as.character(optional$cluster) else rep(NA_character_, k),
    slab          = if (!is.null(optional$slab))      optional$slab      else rep(NA_character_, k),
    weights       = if (!is.null(optional$weights))   optional$weights   else rep(NA, k),
    stringsAsFactors = FALSE
  )

  return(list(
    data_outcome       = data_outcome,
    k                  = k,
    mods_from_yi       = NULL,  # GLMM does not support formula syntax
    formula_yi         = NULL,
    weights_provided   = optional$weights_provided,
    slab_provided      = optional$slab_provided,
    cluster_provided   = optional$cluster_provided,
    na_check_cols      = c("ai", "ci", "n1i", "n2i"),  # Check cell counts for NAs
    effect_direction   = "positive", # Non-consequential, for consistency with .norm
    outcome_type       = "bin"
  ))
}


# Internal function to extract and validate outcome variables for Poisson GLMM models
# Handles x1i, x2i, t1i, t2i, weights, cluster, slab (for measure = "IRR")
#
# @param .call Matched call from the calling function
# @param data The data frame (can be NULL)
# @param .envir The enclosing environment
# @param skip_validation Whether to skip strict validation checks (for newdata).
#        When TRUE, allows t1i/t2i = 0 (useful for prediction).
#
# @return A list with:
#   - data_outcome: data.frame with x1i, x2i, t1i, t2i, cluster, slab, weights
#   - k: number of observations
#   - mods_from_yi: NULL (GLMM does not support formula syntax for outcome)
#   - formula_yi: NULL
#   - slab_provided: logical, whether slab was provided by user
#   - cluster_provided: logical, whether cluster was provided by user
#   - na_check_cols: character vector of column names to check for NAs
.check_and_list_data.outcome.pois <- function(.call, data, .envir, skip_validation = FALSE) {

  # Extract event counts and person-time for Poisson models
  x1i <- .get_variable(.call, data, .envir, "x1i", allow_NULL = TRUE)
  x2i <- .get_variable(.call, data, .envir, "x2i", allow_NULL = TRUE)
  t1i <- .get_variable(.call, data, .envir, "t1i", allow_NULL = TRUE)
  t2i <- .get_variable(.call, data, .envir, "t2i", allow_NULL = TRUE)

  ### Validate that all required variables are provided
  if (is.null(x1i) || is.null(x2i) || is.null(t1i) || is.null(t2i)) {
    stop("For Poisson models (measure = 'IRR'), all of 'x1i', 'x2i', 't1i', and 't2i' must be provided.", call. = FALSE)
  }

  # Determine k from x1i
  k <- length(x1i)

  ### Validate inputs
  # x1i: event counts in treatment group
  BayesTools::check_int(x1i, "x1i", check_length = k, allow_NULL = FALSE, allow_NA = TRUE, lower = 0)

  # x2i: event counts in control group
  BayesTools::check_int(x2i, "x2i", check_length = k, allow_NULL = FALSE, allow_NA = TRUE, lower = 0)

  # t1i: person-time in treatment group
  # When skip_validation = TRUE (newdata), allow t1i = 0
  BayesTools::check_real(t1i, "t1i", check_length = k, allow_NULL = FALSE, allow_NA = TRUE, lower = 0, allow_bound = skip_validation)

  # t2i: person-time in control group
  # When skip_validation = TRUE (newdata), allow t2i = 0
  BayesTools::check_real(t2i, "t2i", check_length = k, allow_NULL = FALSE, allow_NA = TRUE, lower = 0, allow_bound = skip_validation)

  # Extract and validate common optional variables (weights, cluster, slab)
  optional <- .check_and_list_data.optional_vars(.call, data, .envir, k, "x1i")

  ### Construct output data frame
  data_outcome <- data.frame(
    x1i           = x1i,
    x2i           = x2i,
    t1i           = t1i,
    t2i           = t2i,
    cluster       = if (!is.null(optional$cluster)) optional$cluster else rep(NA_character_, k),
    cluster_label = if (!is.null(optional$cluster)) as.character(optional$cluster) else rep(NA_character_, k),
    slab          = if (!is.null(optional$slab))      optional$slab      else rep(NA_character_, k),
    weights       = if (!is.null(optional$weights))   optional$weights   else rep(NA, k),
    stringsAsFactors = FALSE
  )

  return(list(
    data_outcome       = data_outcome,
    k                  = k,
    mods_from_yi       = NULL,  # GLMM does not support formula syntax
    formula_yi         = NULL,
    weights_provided   = optional$weights_provided,
    slab_provided      = optional$slab_provided,
    cluster_provided   = optional$cluster_provided,
    na_check_cols      = c("x1i", "x2i", "t1i", "t2i"),  # Check Poisson data for NAs
    effect_direction   = "positive", # Non-consequential, for consistency with .norm
    outcome_type       = "pois"
  ))
}


# Internal function to extract and validate predictor variables (mods or scale)
# Handles formulas, matrices, and data.frames
#
# @param .call Matched call from the calling function
# @param data The data frame (can be NULL)
# @param .envir The enclosing environment
# @param name The name of the argument ("mods" or "scale")
# @param k The expected number of rows
#
# @return A data.frame with predictor variables (with "formula" attribute) or NULL
.check_and_list_data.predictors <- function(.call, data, .envir, name, k) {

  # Check if argument is in the call
  arg_index <- match(name, names(.call))

  if (is.na(arg_index))
    return(NULL)

  mods_expr <- .call[[arg_index]]

  if (is.null(mods_expr))
    return(NULL)

  # Evaluate the expression to get the formula, matrix, or data.frame
  if (!is.null(data) && is.data.frame(data)) {
    mods <- try(eval(mods_expr, data, .envir), silent = TRUE)
  } else {
    mods <- try(eval(mods_expr, .envir), silent = TRUE)
  }

  if (inherits(mods, "try-error"))
    stop(paste0("Cannot evaluate the '", name, "' argument: ",
                conditionMessage(attr(mods, "condition"))), call. = FALSE)

  .check_and_list_data.predictors_value(
    predictors = mods,
    data       = data,
    name       = name,
    k          = k
  )
}

.check_and_list_data.predictors_value <- function(predictors, data, name, k) {

  # Handle formula input
  if (inherits(predictors, "formula")) {

    original_formula <- predictors

    # Ensure formula has no LHS (response)
    if (length(predictors) == 3) {
      warning(paste0("The '", name, "' formula should not have a left-hand side. ",
                     "The LHS will be ignored."), call. = FALSE)
      predictors <- predictors[-2]
      original_formula <- predictors
    }

    .check_and_list_data.validate_fixed_formula(predictors, name)

    # Create model frame from formula
    if (!is.null(data) && is.data.frame(data)) {
      mf <- try(
        stats::model.frame(predictors, data = data, na.action = stats::na.pass),
        silent = TRUE
      )
    } else {
      mf <- try(
        stats::model.frame(predictors, na.action = stats::na.pass),
        silent = TRUE
      )
    }

    if (inherits(mf, "try-error"))
      stop(paste0("Cannot create model frame from '", name, "' formula: ",
                  conditionMessage(attr(mf, "condition"))), call. = FALSE)

    if (nrow(mf) != k)
      stop(paste0("The number of rows in '", name, "' (", nrow(mf),
                  ") must equal the number of effect sizes (", k, ")."), call. = FALSE)

    rownames(mf) <- NULL
    attr(mf, "formula") <- original_formula

    return(mf)

  } else if (is.matrix(predictors)) {

    mods_df <- as.data.frame(predictors)

    if (nrow(mods_df) != k)
      stop(paste0("The number of rows in '", name, "' (", nrow(mods_df),
                  ") must equal the number of effect sizes (", k, ")."), call. = FALSE)

    rownames(mods_df) <- NULL
    attr(mods_df, "formula") <- .create_formula_from_names(names(mods_df))

    return(mods_df)

  } else if (is.data.frame(predictors)) {

    mods_df <- predictors

    if (nrow(mods_df) != k)
      stop(paste0("The number of rows in '", name, "' (", nrow(mods_df),
                  ") must equal the number of effect sizes (", k, ")."), call. = FALSE)

    rownames(mods_df) <- NULL
    attr(mods_df, "formula") <- .create_formula_from_names(names(mods_df))

    return(mods_df)

  } else {
    stop(paste0("The '", name, "' argument must be a formula, matrix, or data.frame."),
         call. = FALSE)
  }
}

.check_and_list_data.scale <- function(.call, data, .envir, k) {

  arg_index <- match("scale", names(.call))
  if (is.na(arg_index)) {
    return(NULL)
  }

  scale_expr <- .call[[arg_index]]
  if (is.null(scale_expr)) {
    return(NULL)
  }

  scale <- if (!is.null(data) && is.data.frame(data)) {
    try(eval(scale_expr, data, .envir), silent = TRUE)
  } else {
    try(eval(scale_expr, .envir), silent = TRUE)
  }
  if (inherits(scale, "try-error")) {
    stop(
      "Cannot evaluate the 'scale' argument: ",
      conditionMessage(attr(scale, "condition")),
      call. = FALSE
    )
  }

  if (is.list(scale) && !is.data.frame(scale) && !inherits(scale, "formula")) {
    if (is.null(names(scale)) || any(!nzchar(names(scale)))) {
      stop(
        "Component-specific 'scale' lists must be named by random component.",
        call. = FALSE
      )
    }
    if (anyDuplicated(names(scale))) {
      stop(
        "Component-specific 'scale' list names must be unique.",
        call. = FALSE
      )
    }

    scale_components <- lapply(names(scale), function(component) {
      component_data <- .check_and_list_data.predictors_value(
        predictors = scale[[component]],
        data       = data,
        name       = paste0("scale$", component),
        k          = k
      )
      attr(component_data, "component") <- component
      attr(component_data, "source")    <- paste0("tau_", component)
      attr(component_data, "parameter") <- paste0("log_tau_", component)
      component_data
    })
    names(scale_components) <- names(scale)
    class(scale_components) <- c("RoBMA_scale_components", "list")

    return(scale_components)
  }

  .check_and_list_data.predictors_value(
    predictors = scale,
    data       = data,
    name       = "scale",
    k          = k
  )
}

.check_and_list_data.validate_scale_component_name <- function(component) {

  if (!grepl("^[A-Za-z][A-Za-z0-9_]*$", component)) {
    stop(
      "Scale component name '", component,
      "' is not a valid JAGS parameter-name fragment.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

.check_and_list_data.scale_component_aliases <- function(component, name,
                                                         scale_name) {

  aliases <- unique(c(component, name, scale_name))
  aliases[!is.na(aliases) & nzchar(aliases)]
}

.check_and_list_data.random <- function(.call, data, .envir, k,
                                        random_effects_metadata = NULL,
                                        random_group_covariance = NULL) {

  arg_index <- match("random", names(.call))
  if (is.na(arg_index)) {
    if (!is.null(random_group_covariance)) {
      stop("'R' requires a 'random' formula.", call. = FALSE)
    }
    if (!is.null(random_effects_metadata)) {
      if (!inherits(random_effects_metadata, "BayesTools_random_effects")) {
        stop("Internal error: random-effect metadata must be a BayesTools random-effect object.",
             call. = FALSE)
      }
      return(list(
        formula        = NULL,
        data           = NULL,
        terms          = random_effects_metadata[["terms"]],
        random_effects = random_effects_metadata
      ))
    }
    return(list(formula = NULL, data = NULL, terms = list()))
  }

  random_expr <- .call[[arg_index]]
  if (is.null(random_expr)) {
    if (!is.null(random_group_covariance)) {
      stop("'R' requires a non-NULL 'random' formula.", call. = FALSE)
    }
    return(list(formula = NULL, data = NULL, terms = list()))
  }

  random <- if (!is.null(data) && is.data.frame(data)) {
    try(eval(random_expr, data, .envir), silent = TRUE)
  } else {
    try(eval(random_expr, .envir), silent = TRUE)
  }
  if (inherits(random, "try-error")) {
    stop(
      "Cannot evaluate the 'random' argument: ",
      conditionMessage(attr(random, "condition")),
      call. = FALSE
    )
  }

  if (inherits(random, "BayesTools_random_effects")) {
    if (!is.null(random_group_covariance)) {
      stop(
        "'R' cannot be combined with a pre-parsed BayesTools random-effect ",
        "object; supply the random-effect formula to 'random' instead.",
        call. = FALSE
      )
    }
    random_effects <- random
  } else {
    random_effects <- tryCatch(
      BayesTools::random_effects_formula(
        random           = random,
        envir            = .envir,
        group_covariance = random_group_covariance
      ),
      error = function(e) {
        message <- conditionMessage(e)
        if (!is.null(random_group_covariance)) {
          message <- .brma_mv_group_covariance_error_message(message)
        }
        stop(message, call. = FALSE)
      }
    )
    .check_and_list_data.random_validate_formula_component(
      random_effects = random_effects,
      random         = random
    )
  }
  formula <- random_effects[["formula"]]
  terms   <- random_effects[["terms"]]
  .check_and_list_data.random_validate_terms(terms)
  .check_and_list_data.random_validate_group_covariance_terms(terms)

  variables <- all.vars(formula)
  random_data <- .check_and_list_data.random_variables(
    variables = variables,
    data      = data,
    .envir    = .envir,
    k         = k
  )
  attr(random_data, "formula") <- formula

  return(list(
    formula        = formula,
    data           = random_data,
    terms          = terms,
    random_effects = random_effects
  ))
}

.check_and_list_data.random_validate_group_covariance_terms <- function(terms) {

  for (term in terms) {
    if (!.random_effect_term_has_known_group_covariance(term)) {
      next
    }

    block <- term[["block_name"]]
    structure <- term[["structure"]]
    structure_label <- if (is.character(structure) && length(structure) == 1L &&
                           !is.na(structure) && nzchar(structure)) {
      structure
    } else {
      "unknown"
    }
    if (!is.character(structure) || length(structure) != 1L ||
        is.na(structure) || !tolower(structure) %in% c("id", "diag", "us")) {
      stop(
        "'R' currently supports random-intercept blocks only; ",
        "random-effect block '", block, "' uses structure '",
        structure_label, "'.",
        call. = FALSE
      )
    }
    if (!identical(term[["expr"]], 1)) {
      stop(
        "'R' currently supports random intercepts only; ",
        "random-effect block '", block, "' contains random slopes.",
        call. = FALSE
      )
    }
  }

  invisible(TRUE)
}

.check_and_list_data.random_validate_formula_component <- function(random_effects,
                                                                   random) {

  is_single_formula <- inherits(random, "formula") ||
    (is.list(random) && length(random) == 1L &&
       inherits(random[[1L]], "formula"))
  if (!is_single_formula) {
    return(invisible(TRUE))
  }

  terms <- random_effects[["terms"]]
  if (length(terms) == 0L) {
    return(invisible(TRUE))
  }
  component_visible <- vapply(terms, function(term) {
    isTRUE(term[["component_visible"]])
  }, logical(1))
  if (any(component_visible)) {
    return(invisible(TRUE))
  }
  if (length(terms) > 1L &&
      !.check_and_list_data.random_is_plain_nested(random)) {
    stop(
      "Bare 'random' formulas with multiple random-effect terms are ",
      "ambiguous. Use a named list such as ",
      "'random = list(component = ~ 1 | study)'.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

.check_and_list_data.random_is_plain_nested <- function(random) {

  rhs <- random[[if (length(random) == 3L) 3L else 2L]]
  while (is.call(rhs) && identical(rhs[[1L]], as.name("("))) {
    rhs <- rhs[[2L]]
  }

  is.call(rhs) &&
    identical(rhs[[1L]], as.name("|")) &&
    .check_and_list_data.random_group_has_nested_slash(rhs[[3L]])
}

.check_and_list_data.random_group_has_nested_slash <- function(expr) {

  if (!is.call(expr)) {
    return(FALSE)
  }
  if (identical(expr[[1L]], as.name("/"))) {
    return(TRUE)
  }

  any(vapply(as.list(expr[-1L]), .check_and_list_data.random_group_has_nested_slash,
             logical(1)))
}

.check_and_list_data.formula_rhs <- function(formula) {

  rhs_index <- if (length(formula) == 3L) 3L else 2L
  paste(
    deparse(formula[[rhs_index]], width.cutoff = 500L, backtick = TRUE),
    collapse = " "
  )
}

.check_and_list_data.random_validate_terms <- function(terms) {

  if (length(terms) == 0L) {
    stop("The 'random' formula must contain at least one random-effect term.",
         call. = FALSE)
  }

  for (term in terms) {
    is_plain <- !isTRUE(term[["explicit_special"]])
    if (is_plain && !identical(term[["expr"]], 1)) {
      stop(
        "Plain 'random' terms are supported only for random intercepts ",
        "such as '~ 1 | study'. Use an explicit covariance wrapper ",
        "such as 'diag()', 'us()', 'cs()', or 'ar1()', or the '||' ",
        "diagonal shorthand, for random slopes.",
        call. = FALSE
      )
    }
  }

  invisible(TRUE)
}

.check_and_list_data.random_variables <- function(variables, data, .envir, k) {

  if (length(variables) == 0L) {
    out <- data.frame(.RoBMA_random_intercept = rep(1, k))
    return(out[, FALSE, drop = FALSE])
  }

  out <- data.frame(row.names = seq_len(k))
  for (variable in variables) {
    value <- .check_and_list_data.random_variable(
      variable = variable,
      data     = data,
      .envir   = .envir
    )
    if (length(value) != k) {
      stop(
        "The random-effect variable '", variable,
        "' must have length ", k, " (same as 'yi').",
        call. = FALSE
      )
    }
    out[[variable]] <- value
  }
  rownames(out) <- NULL

  return(out)
}

.check_and_list_data.random_variable <- function(variable, data, .envir) {

  if (!is.null(data) && is.data.frame(data) && variable %in% names(data)) {
    return(data[[variable]])
  }
  if (exists(variable, envir = .envir, inherits = TRUE)) {
    return(get(variable, envir = .envir, inherits = TRUE))
  }

  stop(
    "Cannot find the random-effect variable ('", variable, "').",
    call. = FALSE
  )
}

.check_and_list_data.location <- function(data_mods, data_random) {

  if (is.null(data_mods) && is.null(data_random[["data"]])) {
    return(NULL)
  }

  if (!is.null(data_mods)) {
    location <- data_mods
  } else {
    location <- data.frame(
      row.names = seq_len(nrow(data_random[["data"]]))
    )
  }
  if (!is.null(data_random[["data"]])) {
    for (variable in names(data_random[["data"]])) {
      if (variable %in% names(location)) {
        if (!identical(location[[variable]], data_random[["data"]][[variable]])) {
          stop(
            "The variable '", variable, "' is used by both 'mods' and ",
            "'random' with different values. Use distinct names or make ",
            "both inputs refer to the same data column.",
            call. = FALSE
          )
        }
      } else {
        location[[variable]] <- data_random[["data"]][[variable]]
      }
    }
  }
  rownames(location) <- NULL

  attr(location, "formula") <- .check_and_list_data.location_formula(
    data_mods   = data_mods,
    data_random = data_random
  )
  attr(location, "random_effects") <- data_random[["random_effects"]]

  return(location)
}

.check_and_list_data.location_formula <- function(data_mods, data_random) {

  fixed_rhs <- if (!is.null(data_mods)) {
    .check_and_list_data.formula_rhs(attr(data_mods, "formula"))
  } else {
    "1"
  }

  if (is.null(data_random[["formula"]])) {
    return(stats::as.formula(paste("~", fixed_rhs)))
  }

  random_rhs <- .check_and_list_data.formula_rhs(data_random[["formula"]])

  formula <- stats::as.formula(
    paste("~", fixed_rhs, "+", random_rhs),
    env = environment(data_random[["formula"]])
  )
  attr(formula, "random_terms") <- data_random[["terms"]]
  if (!is.null(data_random[["random_effects"]][["components"]])) {
    attr(formula, "random_components") <- data_random[["random_effects"]][["components"]]
  }

  formula
}


# Internal function to coerce character predictors to factors
# Character moderators and scale predictors must be factors before downstream
# contrast handling builds design matrices.
#
# @param df A data.frame with predictor variables (can be NULL)
#
# @return The input data.frame with character columns converted to factors
.check_and_list_data.coerce_character_predictors <- function(df) {

  if (is.null(df))
    return(NULL)

  character_columns <- vapply(df, is.character, logical(1))

  if (!any(character_columns))
    return(df)

  for (col_name in names(df)[character_columns]) {
    df[[col_name]] <- factor(df[[col_name]])
  }

  return(df)
}

.check_and_list_data.scale_coerce_character_predictors <- function(data_scale) {

  if (.check_and_list_data.is_scale_components(data_scale)) {
    data_scale <- lapply(data_scale, .check_and_list_data.coerce_character_predictors)
    class(data_scale) <- c("RoBMA_scale_components", "list")
    return(data_scale)
  }

  .check_and_list_data.coerce_character_predictors(data_scale)
}


# Internal function to validate and convert subset argument to logical vector
#
# @param subset The subset argument (logical or numeric vector)
# @param k The expected length (number of observations)
#
# @return A logical vector of length k
.check_and_list_data.validate_subset <- function(subset, k) {

  if (is.logical(subset)) {

    if (length(subset) != k)
      stop(paste0("The 'subset' argument must have length ", k, " (same as the outcome) when logical."), call. = FALSE)

    # Check for NAs in subset
    if (any(is.na(subset)))
      stop("The 'subset' argument must not contain NA values.", call. = FALSE)

    return(subset)

  } else if (is.numeric(subset)) {

    # Check for NAs in subset
    if (any(is.na(subset)))
      stop("The 'subset' argument must not contain NA values.", call. = FALSE)

    if (any(!is.finite(subset)) || any(subset != floor(subset)))
      stop("The 'subset' argument must contain integer numeric indices.", call. = FALSE)

    if (any(subset < 1) || any(subset > k))
      stop(paste0("The 'subset' argument must contain values between 1 and ", k, "."), call. = FALSE)

    # Convert numeric indices to logical
    subset_logical <- rep(FALSE, k)
    subset_logical[subset] <- TRUE
    return(subset_logical)

  } else {
    stop("The 'subset' argument must be a logical or numeric vector.", call. = FALSE)
  }
}


# Internal function to apply subsetting to a data.frame
# Preserves formula attributes and drops unused factor levels
#
# @param df A data.frame (can be NULL)
# @param subset A logical vector indicating which rows to keep
#
# @return The subsetted data.frame (or NULL if input was NULL)
.check_and_list_data.subset <- function(df, subset) {

  if (is.null(df))
    return(NULL)

  # Preserve attributes
  saved_formula    <- attr(df, "formula")
  saved_formula_yi <- attr(df, "formula_yi")
  saved_component  <- attr(df, "component")
  saved_source     <- attr(df, "source")
  saved_parameter  <- attr(df, "parameter")

  # Apply subset
  df <- df[subset, , drop = FALSE]

  # Drop unused factor levels
  df <- droplevels(df)

  # Reset row names
  rownames(df) <- NULL

  # Restore attributes
  if (!is.null(saved_formula))
    attr(df, "formula") <- saved_formula
  if (!is.null(saved_formula_yi))
    attr(df, "formula_yi") <- saved_formula_yi
  if (!is.null(saved_component))
    attr(df, "component") <- saved_component
  if (!is.null(saved_source))
    attr(df, "source") <- saved_source
  if (!is.null(saved_parameter))
    attr(df, "parameter") <- saved_parameter

  return(df)
}

.check_and_list_data.scale_subset <- function(data_scale, subset) {

  if (.check_and_list_data.is_scale_components(data_scale)) {
    data_scale <- lapply(data_scale, .check_and_list_data.subset, subset = subset)
    class(data_scale) <- c("RoBMA_scale_components", "list")
    return(data_scale)
  }

  .check_and_list_data.subset(data_scale, subset)
}


# Internal function to validate predictor variables (mods or scale)
# Checks that:
#   - Continuous variables have sd > 0 (are not constant)
#   - Factor variables have more than one level
#
# These checks are appropriate for original data but should be skipped for
# newdata used in prediction, where single-value predictions are legitimate.
#
# @param df A data.frame with predictor variables (can be NULL)
# @param name The name of the predictor type ("mods" or "scale") for error messages
# @param skip_validation Logical; if TRUE, skip all validation checks.
#        Should be TRUE when processing newdata for prediction.
#
# @return NULL (invisibly); stops with error if validation fails
.check_and_list_data.validate_predictors <- function(df, name, skip_validation = FALSE) {

  if (is.null(df))
    return(invisible(NULL))

  # Skip validation for newdata (e.g., predicting at a single factor level or

  # single continuous value is valid for prediction)
  if (skip_validation)
    return(invisible(NULL))

  for (col_name in names(df)) {
    col <- df[[col_name]]

    if (is.factor(col) || is.character(col)) {
      # Factor/character: must have more than one unique level
      n_levels <- length(unique(col))
      if (n_levels < 2) {
        stop(paste0("The '", name, "' variable '", col_name,
                    "' has only one level. Factor predictors must have at least two levels."),
             call. = FALSE)
      }
    } else if (is.numeric(col)) {
      # Numeric: must have sd > 0 (not constant)
      col_sd <- stats::sd(col, na.rm = TRUE)
      if (is.na(col_sd) || col_sd == 0) {
        stop(paste0("The '", name, "' variable '", col_name,
                    "' has zero variance (all values are identical). ",
                    "Continuous predictors must have variation."),
             call. = FALSE)
      }
    }
    # For other types (logical, etc.), no specific check
  }

  return(invisible(NULL))
}

.check_and_list_data.scale_validate_predictors <- function(data_scale, skip_validation = FALSE) {

  if (.check_and_list_data.is_scale_components(data_scale)) {
    for (component in names(data_scale)) {
      .check_and_list_data.validate_predictors(
        df              = data_scale[[component]],
        name            = paste0("scale$", component),
        skip_validation = skip_validation
      )
    }
    return(invisible(NULL))
  }

  .check_and_list_data.validate_predictors(
    df              = data_scale,
    name            = "scale",
    skip_validation = skip_validation
  )
}

.check_and_list_data.scale_na_frames <- function(data_scale) {

  if (is.null(data_scale)) {
    return(list())
  }
  if (.check_and_list_data.is_scale_components(data_scale)) {
    out <- as.list(data_scale)
    names(out) <- paste0("scale$", names(out))
    return(out)
  }

  list(scale = data_scale)
}

.check_and_list_data.is_scale_components <- function(data_scale) {

  inherits(data_scale, "RoBMA_scale_components")
}

.check_and_list_data.validate_scale_random <- function(data_scale, data_random) {

  if (is.null(data_scale)) {
    return(NULL)
  }
  if (length(data_random[["terms"]]) == 0L) {
    if (.check_and_list_data.is_scale_components(data_scale)) {
      stop(
        "Component-specific 'scale' lists require a 'random' formula.",
        call. = FALSE
      )
    }
    return(data_scale)
  }

  components       <- .check_and_list_data.random_components(data_random[["terms"]])
  component_labels <- components[["label"]]

  if (.check_and_list_data.is_scale_components(data_scale)) {
    component_map <- stats::setNames(component_labels, component_labels)
    component_map[components[["name"]]] <- component_labels

    scale_names        <- names(data_scale)
    mapped_components  <- unname(component_map[scale_names])
    unknown_components <- names(data_scale)[is.na(mapped_components)]
    if (length(unknown_components) > 0L) {
      stop(
        "Component-specific 'scale' names must match top-level 'random' ",
        "components. Unknown: ",
        .check_and_list_data.collapse_or_none(unknown_components),
        ".",
        call. = FALSE
      )
    }
    if (anyDuplicated(mapped_components)) {
      stop(
        "Component-specific 'scale' names must uniquely match top-level ",
        "'random' components after BayesTools name normalization.",
        call. = FALSE
      )
    }

    scaled_components <- component_labels[component_labels %in% mapped_components]
    scale_names_by_component <- scale_names[match(scaled_components, mapped_components)]
    names(scale_names_by_component) <- scaled_components
    component_names <- stats::setNames(components[["name"]], component_labels)

    data_scale <- data_scale[match(scaled_components, mapped_components)]
    names(data_scale) <- scaled_components
    for (component in scaled_components) {
      .check_and_list_data.validate_scale_component_name(component)
      component_name <- component_names[[component]]
      scale_name     <- scale_names_by_component[[component]]
      attr(data_scale[[component]], "component")      <- component
      attr(data_scale[[component]], "component_name") <- component_name
      attr(data_scale[[component]], "scale_name")     <- scale_name
      attr(data_scale[[component]], "aliases")        <-
        .check_and_list_data.scale_component_aliases(
          component  = component,
          name       = component_name,
          scale_name = scale_name
        )
      attr(data_scale[[component]], "source")         <- paste0("tau_", component)
      attr(data_scale[[component]], "parameter")      <- paste0("log_tau_", component)
    }
    class(data_scale) <- c("RoBMA_scale_components", "list")

    return(data_scale)
  }

  if (length(component_labels) > 1L) {
    stop(
      "A single 'scale' formula is ambiguous when 'random' has multiple ",
      "top-level components. Use a named list such as ",
      "'scale = list(component = ~ x)'.",
      call. = FALSE
    )
  }

  attr(data_scale, "component")      <- component_labels[[1L]]
  attr(data_scale, "component_name") <- components[["name"]][[1L]]
  attr(data_scale, "aliases")        <- .check_and_list_data.scale_component_aliases(
    component  = component_labels[[1L]],
    name       = components[["name"]][[1L]],
    scale_name = NULL
  )

  data_scale
}

.check_and_list_data.random_components <- function(terms) {

  labels <- .check_and_list_data.random_component_labels(terms)
  names  <- vapply(seq_along(terms), function(i) {
    component <- terms[[i]][["component"]]
    if (is.null(component) || length(component) != 1L ||
        is.na(component) || !nzchar(component)) {
      return(labels[[i]])
    }
    component
  }, character(1))

  keep <- !duplicated(labels)
  data.frame(
    label = labels[keep],
    name  = names[keep],
    stringsAsFactors = FALSE
  )
}

.check_and_list_data.collapse_or_none <- function(x) {

  if (length(x) == 0L) {
    return("<none>")
  }

  paste(x, collapse = ", ")
}

.check_and_list_data.random_component_labels <- function(terms) {

  vapply(seq_along(terms), function(i) {
    label <- terms[[i]][["component_label"]]
    if (is.null(label) || length(label) != 1L ||
        is.na(label) || !nzchar(label)) {
      return(terms[[i]][["block_name"]])
    }
    label
  }, character(1))
}


# Internal function to check for NA values across a list of data.frames
# Returns a logical vector indicating which rows have at least one NA
#
# @param data_list A named list of data.frames to check
#
# @return A logical vector where TRUE indicates the row has at least one NA
.check_and_list_data.check_na <- function(data_list) {

  if (length(data_list) == 0)
    return(logical(0))

  # Combine all data.frames into one for vectorized NA check
  combined <- do.call(cbind, Filter(function(df) !is.null(df) && is.data.frame(df), data_list))

  # Vectorized NA check: TRUE if any column in that row has NA
  na_rows <- rowSums(is.na(combined)) > 0

  return(na_rows)
}

.check_and_list_data.stop_na <- function(data_list) {

  na_entries <- do.call(rbind, lapply(names(data_list), function(frame_name) {
    df <- data_list[[frame_name]]
    if (is.null(df) || !is.data.frame(df) || ncol(df) == 0L) {
      return(NULL)
    }

    indices <- which(is.na(df), arr.ind = TRUE)
    if (nrow(indices) == 0L) {
      return(NULL)
    }

    data.frame(
      row    = indices[, "row"],
      column = paste0(frame_name, "$", colnames(df)[indices[, "col"]]),
      stringsAsFactors = FALSE
    )
  }))

  rows    <- sort(unique(na_entries[["row"]]))
  columns <- unique(na_entries[["column"]])

  stop(
    "Prediction data must not contain missing values. Missing values found in row(s): ",
    paste(rows, collapse = ", "),
    "; column(s): ",
    paste(columns, collapse = ", "),
    ".",
    call. = FALSE
  )
}


# Internal helper to create a formula from column names
# Creates a one-sided formula like ~ var1 + var2 + var3
#
# @param col_names Character vector of column names
# @return A formula object
.create_formula_from_names <- function(col_names) {

  if (length(col_names) == 0) {
    return(stats::as.formula("~ 1", env = baseenv()))
  }

  # Backtick names that need protection (contain spaces, special chars, etc.)
  col_names_safe <- vapply(col_names, function(nm) {
    if (make.names(nm) != nm) {
      paste0("`", nm, "`")
    } else {
      nm
    }
  }, character(1), USE.NAMES = FALSE)

  # Create formula string and convert to formula
  formula_str <- paste("~", paste(col_names_safe, collapse = " + "))
  return(stats::as.formula(formula_str, env = baseenv()))
}


# Detach persisted formula metadata from caller evaluation frames. Predictor
# values have already been materialized in data frames at this boundary, so
# retaining caller environments only serializes unrelated session objects.
.check_and_list_data.detach_formula_environments <- function(x) {

  if (is.environment(x)) {
    return(x)
  }

  if (inherits(x, "formula")) {
    environment(x) <- baseenv()
  }
  if (inherits(x, "terms")) {
    attr(x, ".Environment") <- baseenv()
  }

  if (is.list(x)) {
    for (i in seq_along(x)) {
      x[i] <- list(.check_and_list_data.detach_formula_environments(x[[i]]))
    }
  }

  metadata <- attributes(x)
  if (!is.null(metadata)) {
    metadata_names <- intersect(
      names(metadata),
      c("formula", "formula_yi", "terms", "random_terms",
        "random_components", "random_effects")
    )
    for (name in metadata_names) {
      attr(x, name) <- .check_and_list_data.detach_formula_environments(
        attr(x, name, exact = TRUE)
      )
    }
  }

  return(x)
}


#' @title Print method for RoBMA_data objects
#'
#' @description Prints a concise summary of an RoBMA_data object.
#'
#' @param x an RoBMA_data object
#' @param n maximum number of rows to display for each component. Defaults to 6.
#' @param ... additional arguments (ignored)
#'
#' @return Invisibly returns the input object.
#'
#' @export
print.RoBMA_data <- function(x, n = 6, ...) {

  k_final   <- attr(x, "k_final")
  n_dropped <- attr(x, "n_dropped")

  # Outcome data header
  if (n_dropped > 0) {
    cat(sprintf("Outcome data (k = %d, dropped = %d)\n", k_final, n_dropped))
  } else {
    cat(sprintf("Outcome data (k = %d)\n", k_final))
  }
  print(utils::head(x$outcome, n = n))

  # Moderators (effect)
  if (!is.null(x$mods)) {
    cat("\nModerators (effect)\n")
    print(utils::head(x$mods, n = n))
  }

  # Moderators (scale)
  if (!is.null(x$scale)) {
    cat("\nModerators (scale)\n")
    print(utils::head(x$scale, n = n))
  }

  # Indicate if more rows exist

  if (k_final > n) {
    cat(sprintf("\n... with %d more rows\n", k_final - n))
  }

  invisible(x)
}


# Internal helper to normalize data.frame/list newdata inputs.
.prepare_newdata_as_data_frame <- function(newdata) {

  if (is.list(newdata) && !is.data.frame(newdata)) {
    newdata <- as.data.frame(newdata, stringsAsFactors = FALSE)
  }

  return(newdata)
}


# Internal helper to report missing newdata columns.
.prepare_newdata_stop_missing <- function(newdata, cols) {

  missing_cols <- setdiff(cols, names(newdata))
  if (length(missing_cols) == 0L) {
    return(invisible(TRUE))
  }

  stop(
    "The 'newdata' must contain columns: ",
    paste(missing_cols, collapse = ", "), ".",
    call. = FALSE
  )
}


# Parser placeholders keep the shared data-input path usable without making
# absent outcome variables visible to fitted moderator, scale, or random
# formulas. The namespace is reserved so user data cannot shadow a placeholder.
.prepare_newdata_placeholder_name <- function(name) {

  return(paste0(".RoBMA_newdata_parser_", name))
}


.prepare_newdata_reject_placeholder_columns <- function(newdata) {

  reserved <- startsWith(names(newdata), ".RoBMA_newdata_parser_")
  if (any(reserved)) {
    stop(
      "The 'newdata' contains a column reserved for internal prediction parsing: ",
      paste(names(newdata)[reserved], collapse = ", "), ".",
      call. = FALSE
    )
  }

  invisible(TRUE)
}


.prepare_newdata_add_placeholder <- function(newdata, name, value, n) {

  placeholder <- .prepare_newdata_placeholder_name(name)
  if (length(value) == 1L) {
    value <- rep(value, n)
  } else if (length(value) != n) {
    stop("Internal error: prediction placeholder length mismatch.", call. = FALSE)
  }
  newdata[[placeholder]] <- value
  return(newdata)
}


.prepare_newdata_parser_column <- function(newdata, name) {

  if (name %in% names(newdata)) {
    return(name)
  }

  placeholder <- .prepare_newdata_placeholder_name(name)
  if (placeholder %in% names(newdata)) {
    return(placeholder)
  }

  return(NULL)
}


# Internal helper to decide whether normal newdata needs sampling SE.
.prepare_newdata_needs_norm_sei <- function(object, type, bias_adjusted) {

  return(
    type == "response" ||
      (
        !bias_adjusted &&
          type %in% c("terms", "cluster", "estimate") &&
          (.is_PET(object) || .is_PEESE(object))
      )
  )
}


# Internal helper to add outcome placeholders for the shared data parser.
.prepare_newdata_outcome <- function(object, newdata, type, bias_adjusted) {

  outcome_type <- .outcome_type(object)
  n_new        <- nrow(newdata)

  if (outcome_type == "norm") {
    if (!"yi" %in% names(newdata)) {
      newdata <- .prepare_newdata_add_placeholder(
        newdata = newdata,
        name    = "yi",
        value   = 0,
        n       = n_new
      )
    }

    if (!("sei" %in% names(newdata) || "vi" %in% names(newdata))) {
      if (.prepare_newdata_needs_norm_sei(object, type, bias_adjusted)) {
        stop(
          "The 'newdata' must contain either 'sei' (standard error) or ",
          "'vi' (variance) column.",
          call. = FALSE
        )
      }
      newdata <- .prepare_newdata_add_placeholder(
        newdata = newdata,
        name    = "sei",
        value   = 0,
        n       = n_new
      )
    }

  } else if (outcome_type == "bin") {
    count_names      <- c("ai", "bi", "ci", "di", "n1i", "n2i")
    provided_counts  <- intersect(count_names, names(newdata))
    totals_available <- all(c("n1i", "n2i") %in% provided_counts)

    if (length(provided_counts) == 0L && type != "response") {
      for (name in c("ai", "ci", "n1i", "n2i")) {
        newdata <- .prepare_newdata_add_placeholder(
          newdata = newdata,
          name    = name,
          value   = 0L,
          n       = n_new
        )
      }
    } else {
      count_args <- lapply(count_names, function(name) {
        if (name %in% names(newdata)) {
          return(newdata[[name]])
        }
        if (totals_available && name %in% c("ai", "ci")) {
          return(rep(0L, n_new))
        }
        return(NULL)
      })
      names(count_args) <- count_names
      count_args[["skip_validation"]] <- TRUE
      counts <- do.call(.canonicalize_binomial_counts, count_args)

      for (name in c("ai", "ci", "n1i", "n2i")) {
        if (name %in% provided_counts) {
          newdata[[name]] <- counts[[name]]
        } else {
          newdata <- .prepare_newdata_add_placeholder(
            newdata = newdata,
            name    = name,
            value   = counts[[name]],
            n       = n_new
          )
        }
      }
    }

  } else if (outcome_type == "pois") {
    if (type == "response") {
      .prepare_newdata_stop_missing(newdata, c("t1i", "t2i"))
    }
    for (name in c("x1i", "x2i", "t1i", "t2i")) {
      if (!name %in% names(newdata)) {
        newdata <- .prepare_newdata_add_placeholder(
          newdata = newdata,
          name    = name,
          value   = 0,
          n       = n_new
        )
      }
    }
  }

  return(newdata)
}


# Internal helper to construct outcome arguments for `.check_and_list_data`.
.prepare_newdata_outcome_call_args <- function(outcome_type, newdata,
                                               prefer_vi = FALSE) {

  if (outcome_type == "norm") {
    call_args <- list(
      yi = as.name(.prepare_newdata_parser_column(newdata, "yi"))
    )

    sampling_names <- if (prefer_vi && "vi" %in% names(newdata)) {
      "vi"
    } else {
      c("vi", "sei")
    }
    for (name in sampling_names) {
      column <- .prepare_newdata_parser_column(newdata, name)
      if (!is.null(column)) {
        call_args[[name]] <- as.name(column)
      }
    }

  } else if (outcome_type == "bin") {
    call_args <- lapply(c("ai", "ci", "n1i", "n2i"), function(name) {
      as.name(.prepare_newdata_parser_column(newdata, name))
    })
    names(call_args) <- c("ai", "ci", "n1i", "n2i")

  } else if (outcome_type == "pois") {
    call_args <- lapply(c("x1i", "x2i", "t1i", "t2i"), function(name) {
      as.name(.prepare_newdata_parser_column(newdata, name))
    })
    names(call_args) <- c("x1i", "x2i", "t1i", "t2i")
  }

  return(call_args)
}


# Internal helper to validate formula variables before evaluation.
.prepare_newdata_validate_formula_vars <- function(newdata, formula, label) {

  missing_vars <- setdiff(all.vars(formula), names(newdata))
  if (length(missing_vars) == 0L) {
    return(invisible(TRUE))
  }

  stop(
    "The 'newdata' must contain all ", label, " variables. Missing: ",
    paste(missing_vars, collapse = ", "), ".",
    call. = FALSE
  )
}


.prepare_newdata_validate_formula_arg <- function(newdata, formula_arg, label) {

  if (is.null(formula_arg)) {
    return(invisible(TRUE))
  }
  if (inherits(formula_arg, "formula")) {
    .prepare_newdata_validate_formula_vars(newdata, formula_arg, label)
    return(invisible(TRUE))
  }
  if (inherits(formula_arg, "BayesTools_random_effects")) {
    .prepare_newdata_validate_formula_vars(newdata, formula_arg[["formula"]], label)
    return(invisible(TRUE))
  }
  if (is.list(formula_arg)) {
    for (i in seq_along(formula_arg)) {
      .prepare_newdata_validate_formula_vars(newdata, formula_arg[[i]], label)
    }
    return(invisible(TRUE))
  }

  stop("Internal error: unsupported formula argument for newdata.", call. = FALSE)
}


.prepare_newdata_formula_arg_variables <- function(formula_arg) {

  if (is.null(formula_arg)) {
    return(character())
  }
  if (inherits(formula_arg, "formula")) {
    return(all.vars(formula_arg))
  }
  if (inherits(formula_arg, "BayesTools_random_effects")) {
    return(all.vars(formula_arg[["formula"]]))
  }
  if (is.list(formula_arg)) {
    return(unique(unlist(
      lapply(formula_arg, .prepare_newdata_formula_arg_variables),
      use.names = FALSE
    )))
  }

  stop("Internal error: unsupported formula argument for newdata.", call. = FALSE)
}


.prepare_newdata_scale_formula_arg <- function(original_data) {

  if (!.is_data_scale(original_data)) {
    return(NULL)
  }

  if (!inherits(original_data[["scale"]], "RoBMA_scale_components")) {
    return(attr(original_data[["scale"]], "formula"))
  }

  scale_specs <- .data_scale_component_specs(original_data)
  out <- lapply(scale_specs, `[[`, "formula")
  names(out) <- vapply(scale_specs, function(scale_spec) {
    scale_name <- scale_spec[["scale_name"]]
    if (!is.null(scale_name) && length(scale_name) == 1L &&
        !is.na(scale_name) && nzchar(scale_name)) {
      return(scale_name)
    }
    scale_spec[["display_name"]]
  }, character(1))

  return(out)
}


.prepare_newdata_random_formula_arg <- function(original_data) {

  if (!.is_data_random(original_data)) {
    return(NULL)
  }

  random_effects <- attr(original_data[["location"]], "random_effects")
  if (is.null(random_effects) || is.null(random_effects[["formula"]])) {
    stop("Internal error: missing fitted random formula.", call. = FALSE)
  }

  return(random_effects)
}


.prepare_newdata_random_grouping_values <- function(newdata, random_effects,
                                                    protected_variables = character()) {

  terms <- random_effects[["terms"]]
  if (length(terms) == 0L) {
    return(newdata)
  }

  random_design_variables <- unique(unlist(lapply(terms, function(term) {
    all.vars(term[["expr"]])
  }), use.names = FALSE))
  protected_variables <- unique(c(protected_variables, random_design_variables))

  for (term in terms) {
    grouping_vars <- all.vars(term[["group_expr"]])
    missing_vars  <- setdiff(grouping_vars, names(newdata))
    if (length(missing_vars) == 0L) {
      next
    }
    if (.random_effect_term_has_known_group_covariance(term)) {
      stop(
        "Known-R new-effect prediction requires explicit fitted grouping ",
        "variable(s) for block '", term[["block_name"]], "': ",
        paste(missing_vars, collapse = ", "), ". Unseen known-R levels ",
        "require 'R_new', which is not supported.",
        call. = FALSE
      )
    }

    synthesis_vars <- setdiff(missing_vars, protected_variables)
    for (variable in synthesis_vars) {
      newdata[[variable]] <- paste0(
        ".RoBMA_new_", term[["block_name"]], "_", seq_len(nrow(newdata))
      )
    }
  }

  return(newdata)
}


# Internal helper function to prepare newdata for prediction
# Reuses `.check_and_list_data` by constructing appropriate call and environment
#
# The newdata data.frame must contain all variables used by the requested
# prediction. Moderator and scale variables are always required when referenced
# by the fitted formulas. Outcome columns are required only when the requested
# prediction uses the sampling distribution or PET/PEESE bias terms; otherwise
# internal dummy values are inserted to keep the shared parser path.
#
# @param object A fitted brma object
# @param newdata A data.frame with new data for prediction.
# @param type Prediction type: "terms", "effect", or "response"
# @param bias_adjusted Whether PET/PEESE terms should be omitted.
# @param include_scale Whether to replay the fitted scale formula.
# @param include_random Whether to replay the fitted random formula.
#
# @return A data list equivalent to `object[["data"]]` but for `newdata`
.prepare_newdata <- function(object, newdata, type, bias_adjusted = FALSE,
                             include_scale = type != "terms",
                             include_random = FALSE) {

  # extract settings from the original fitted object's data attributes
  original_data <- object[["data"]]
  set_contrast_factor_predictors    <- attr(original_data, "set_contrast_factor_predictors")
  standardize_continuous_predictors <- attr(original_data, "standardize_continuous_predictors")
  effect_direction                  <- .effect_direction(object)
  outcome_type                      <- .outcome_type(object)
  extra_env                         <- list()

  newdata <- .prepare_newdata_as_data_frame(newdata)
  .prepare_newdata_reject_placeholder_columns(newdata)
  newdata <- .prepare_newdata_outcome(
    object        = object,
    newdata       = newdata,
    type          = type,
    bias_adjusted = bias_adjusted
  )

  # build the synthetic call expression
  # start with the base call structure
  call_args <- list(quote(.prepare_newdata_call))

  # add data argument
  call_args[["data"]] <- quote(data)

  # add outcome arguments based on outcome_type
  call_args <- c(
    call_args,
    .prepare_newdata_outcome_call_args(
      outcome_type = outcome_type,
      newdata      = newdata,
      prefer_vi    = .is_data_known_v(original_data) && type == "response"
    )
  )

  # add mods formula if this is a regression model
  if (.is_mods(object)) {
    call_args[["mods"]] <- attr(original_data[["mods"]], "formula")
  }

  # add scale formula if present
  if (.is_scale(object) && include_scale) {
    scale_formula_arg <- .prepare_newdata_scale_formula_arg(original_data)
    call_args[["scale"]] <- quote(.RoBMA_scale)
    extra_env[[".RoBMA_scale"]] <- scale_formula_arg
  } else {
    scale_formula_arg <- NULL
  }

  # add random formula when the requested prediction needs random metadata
  if (include_random) {
    random_formula_arg <- .prepare_newdata_random_formula_arg(original_data)
    protected_variables <- character()
    if (.is_mods(object)) {
      protected_variables <- c(
        protected_variables,
        all.vars(attr(original_data[["mods"]], "formula"))
      )
    }
    if (.is_scale(object) && include_scale) {
      protected_variables <- c(
        protected_variables,
        .prepare_newdata_formula_arg_variables(scale_formula_arg)
      )
    }
    newdata <- .prepare_newdata_random_grouping_values(
      newdata             = newdata,
      random_effects      = random_formula_arg,
      protected_variables = protected_variables
    )
    call_args[["random"]] <- quote(.RoBMA_random)
    extra_env[[".RoBMA_random"]] <- random_formula_arg
  } else {
    random_formula_arg <- NULL
  }

  # add cluster structure for multilevel predictions
  if (.is_multilevel(object)) {
    if ("cluster" %in% names(newdata)) {
      call_args[["cluster"]] <- quote(cluster)
    } else {
      n_new <- nrow(newdata)
      newdata[[".RoBMA_cluster"]] <- seq_len(n_new)
      call_args[["cluster"]] <- quote(.RoBMA_cluster)
    }
  }

  # construct the call object
  .call <- as.call(call_args)

  # pre-validate that all formula variables exist in newdata
  # this prevents formulas from finding variables in parent environments
  if (.is_mods(object)) {
    mods_formula <- attr(original_data[["mods"]], "formula")
    .prepare_newdata_validate_formula_vars(newdata, mods_formula, "moderator")
  }

  if (.is_scale(object) && include_scale) {
    .prepare_newdata_validate_formula_arg(newdata, scale_formula_arg, "scale")
  }

  if (include_random) {
    .prepare_newdata_validate_formula_arg(newdata, random_formula_arg, "random-effect")
  }

  # create environment with newdata as "data"
  # use baseenv() as parent to prevent variable leakage from calling context
  .envir           <- new.env(parent = baseenv())
  .envir[["data"]] <- newdata
  for (name in names(extra_env)) {
    .envir[[name]] <- extra_env[[name]]
  }

  # Determine class and measure for .check_and_list_data.
  # brma.mv() known-V response newdata reaches this path after predict.brma()
  # has validated V_new and injected its diagonal as vi.
  measure    <- .measure(object)
  data_class <- switch(
    outcome_type,
    "norm" = "norm",
    "bin"  = "glmm",
    "pois" = "glmm"
  )

  # call .check_and_list_data with extracted settings
  # Use skip_validation = TRUE because:
  # 1. Newdata may have single-level factors (predicting for one group)
  # 2. Newdata may have zero-variance continuous predictors (predicting at one value)
  # 3. Factor levels are validated by BayesTools during formula evaluation
  # 4. Scaling uses original data's mean/sd stored in the fit object
  new_data <- .check_and_list_data(
    .call  = .call,
    .envir = .envir,
    class  = data_class,
    measure = measure,
    set_contrast_factor_predictors    = set_contrast_factor_predictors,
    standardize_continuous_predictors = standardize_continuous_predictors,
    effect_direction                  = effect_direction,
    skip_validation                   = TRUE,
    allow_na_drop                     = FALSE
  )

  return(new_data)
}
