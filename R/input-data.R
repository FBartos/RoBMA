### RoBMA 4.0.0

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
    if (is.atomic(out) && !inherits(out, "formula")) {
      out <- as.vector(out)
    }

    return(out)
  }

  # If direct evaluation failed, check if it's a string referring to a column
  if (is.character(mf_x) && length(mf_x) == 1 && !is.null(data) && is.data.frame(data) && mf_x %in% names(data)) {
    return(as.vector(data[[mf_x]]))
  }

  # If still not found, report error
  if (!allow_NULL) {
    stop(paste0("Cannot find the object/variable ('", deparse(mf_x),
                "') specified for the '", name, "' argument."), call. = FALSE)
  }

  return(NULL)
}


# Internal function to extract common optional variables shared by all outcome handlers
# Handles weights, study_ids, and slab extraction and validation
#
# @param .call Matched call from the calling function
# @param data The data frame (can be NULL)
# @param .envir The enclosing environment
# @param k The expected number of observations
# @param primary_var The name of the primary variable (for error messages, e.g., "yi" or "ai")
#
# @return A list with:
#   - weights: numeric vector or NULL
#   - study_ids: vector or NULL
#   - slab: character vector or NULL
#   - slab_provided: logical
#   - study_ids_provided: logical
.check_and_list_data.optional_vars <- function(.call, data, .envir, k, primary_var) {

  # Extract optional variables
  weights   <- .get_variable(.call, data, .envir, "weights",   allow_NULL = TRUE)
  study_ids <- .get_variable(.call, data, .envir, "study_ids", allow_NULL = TRUE)
  slab      <- .get_variable(.call, data, .envir, "slab",      allow_NULL = TRUE)

  # Track which optional fields were provided
  weights_provided   <- !is.null(weights)
  slab_provided      <- !is.null(slab)
  study_ids_provided <- !is.null(study_ids)

  # Validate weights
  if (!is.null(weights))
    BayesTools::check_real(weights, "weights", check_length = k, allow_NULL = TRUE, allow_NA = TRUE, lower = 0, allow_bound = FALSE)

  # Validate study_ids
  if (!is.null(study_ids) && length(study_ids) != k)
    stop(paste0("The 'study_ids' argument must have length ", k, " (same as '", primary_var, "')."), call. = FALSE)

  # Validate slab (study labels)
  if (!is.null(slab) && length(slab) != k)
    stop(paste0("The 'slab' argument must have length ", k, " (same as '", primary_var, "')."), call. = FALSE)

  return(list(
    weights            = weights,
    study_ids          = study_ids,
    slab               = slab,
    weights_provided   = weights_provided,
    slab_provided      = slab_provided,
    study_ids_provided = study_ids_provided
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
#
# @return A list with components:
#   - outcome: data.frame with outcome variables (columns depend on class)
#   - mods: moderator information (data.frame or NULL)
#   - scale: scale information (data.frame or NULL)
.check_and_list_data <- function(.call, .envir, class = "norm") {

  ### Extract the data argument first - other variables may reference columns within it
  data <- .get_variable(.call, NULL, .envir, "data", allow_NULL = TRUE)

  ### Step 1: Extract and validate outcome variables (dispatch based on class)
  outcome_result <- switch(
    class,
    "norm" = .check_and_list_data.outcome.norm(.call, data, .envir),
    "glmm" = .check_and_list_data.outcome.glmm(.call, data, .envir),
    stop(paste0("Unknown class '", class, "' for .check_and_list_data"), call. = FALSE)
  )

  data_outcome       <- outcome_result$data_outcome
  k                  <- outcome_result$k
  mods_from_yi       <- outcome_result$mods_from_yi
  formula_yi         <- outcome_result$formula_yi
  weights_provided   <- outcome_result$weights_provided
  slab_provided      <- outcome_result$slab_provided
  study_ids_provided <- outcome_result$study_ids_provided
  na_check_cols      <- outcome_result$na_check_cols

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

  data_scale <- .check_and_list_data.predictors(
    .call  = .call,
    data   = data,
    .envir = .envir,
    name   = "scale",
    k      = k
  )

  ### Step 3: Apply subset (before NA handling)
  subset <- .get_variable(.call, data, .envir, "subset", allow_NULL = TRUE)

  if (!is.null(subset)) {
    # Validate and convert subset to logical
    subset <- .check_and_list_data.validate_subset(subset, k)

    # Apply subsetting to all data frames
    data_outcome <- .check_and_list_data.subset(data_outcome, subset)
    data_mods    <- .check_and_list_data.subset(data_mods, subset)
    data_scale   <- .check_and_list_data.subset(data_scale, subset)
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
  if (!is.null(data_scale)) {
    data_list_for_na$scale <- data_scale
  }

  # Get rows with NAs
  na_rows <- .check_and_list_data.check_na(data_list_for_na)

  # Drop NA rows if any
  n_dropped <- sum(na_rows)
  if (n_dropped > 0) {

    warning(paste0(n_dropped, " observation(s) removed due to missing values."), call. = FALSE, immediate. = TRUE)

    keep_rows <- !na_rows
    data_outcome <- .check_and_list_data.subset(data_outcome, keep_rows)
    data_mods    <- .check_and_list_data.subset(data_mods, keep_rows)
    data_scale   <- .check_and_list_data.subset(data_scale, keep_rows)
  }

  ### Step 5: Final validation and processing
  k_final <- nrow(data_outcome)

  if (k_final == 0) {
    stop("No observations remaining after removing missing values.", call. = FALSE)
  }

  # Validate predictor variables (after subsetting and NA dropping)
  .check_and_list_data.validate_predictors(data_mods, "mods")
  .check_and_list_data.validate_predictors(data_scale, "scale")

  # Generate default study labels if not provided (after NA dropping)
  if (!slab_provided) {
    data_outcome$slab <- paste0("Study ", seq_len(k_final))
  }

  # Convert study_ids to numeric factor for use as indices in fitting (after NA dropping)
  if (study_ids_provided) {
    data_outcome[["study_ids"]] <- as.numeric(as.factor(data_outcome[["study_ids"]]))
  }

  data_list <- list(
    outcome = data_outcome,
    mods    = data_mods,
    scale   = data_scale
  )

  class(data_list) <- "RoBMA_data"
  attr(data_list, "n_dropped") <- n_dropped
  attr(data_list, "k_final")   <- k_final
  attr(data_list, "mods")      <- !is.null(data_mods)
  attr(data_list, "scale")     <- !is.null(data_scale)
  attr(data_list, "weights")   <- weights_provided
  attr(data_list, "slab")      <- slab_provided
  attr(data_list, "study_ids") <- study_ids_provided
  return(data_list)
}


# Internal function to extract and validate outcome variables for normal likelihood models
# Handles yi, vi/sei, ni, weights, study_ids, slab, and yi as formula
#
# @param .call Matched call from the calling function
# @param data The data frame (can be NULL)
# @param .envir The enclosing environment
#
# @return A list with:
#   - data_outcome: data.frame with yi, sei, ni, study_ids, slab, weights
#   - k: number of observations
#   - mods_from_yi: moderators extracted from yi formula (or NULL)
#   - formula_yi: the original yi formula (or NULL)
#   - slab_provided: logical, whether slab was provided by user
#   - study_ids_provided: logical, whether study_ids was provided by user
#   - na_check_cols: character vector of column names to check for NAs
.check_and_list_data.outcome.norm <- function(.call, data, .envir) {

  # Extract yi (may be a formula like yi ~ mod1 + mod2)
  yi <- .get_variable(.call, data, .envir, "yi", allow_NULL = FALSE)

  # Handle yi as formula (e.g., yi ~ mod1 + mod2)
  formula_yi   <- NULL
  mods_from_yi <- NULL

  if (inherits(yi, "formula")) {

    formula_yi <- yi

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
  if (!is.null(vi))
    BayesTools::check_real(vi, "vi", check_length = k, allow_NULL = TRUE, allow_NA = TRUE, lower = 0, allow_bound = FALSE)

  # Validate sei (standard error)
  if (!is.null(sei))
    BayesTools::check_real(sei, "sei", check_length = k, allow_NULL = TRUE, allow_NA = TRUE, lower = 0, allow_bound = FALSE)

  # Convert vi to sei or vice versa
  if (is.null(sei) && !is.null(vi)) {
    sei <- sqrt(vi)
  } else if (is.null(vi) && !is.null(sei)) {
    vi <- sei^2
  } else if (is.null(vi) && is.null(sei)) {
    stop("Either 'vi' (variance) or 'sei' (standard error) must be provided.", call. = FALSE)
  } else {
    # Both provided - check consistency
    if (any(abs(vi - sei^2) > 1e-10, na.rm = TRUE))
      stop("The provided 'vi' and 'sei' values are inconsistent.", call. = FALSE)
  }

  # Validate ni (sample sizes)
  if (!is.null(ni))
    BayesTools::check_real(ni, "ni", check_length = k, allow_NULL = TRUE, allow_NA = TRUE, lower = 0, allow_bound = FALSE)

  # Extract and validate common optional variables (weights, study_ids, slab)
  optional <- .check_and_list_data.optional_vars(.call, data, .envir, k, "yi")

  ### Construct output data frame
  data_outcome <- data.frame(
    yi        = yi,
    sei       = sei,
    ni        = if (!is.null(ni))                  ni                  else rep(NA_integer_, k),
    study_ids = if (!is.null(optional$study_ids))  optional$study_ids  else rep(NA_character_, k),
    slab      = if (!is.null(optional$slab))       optional$slab       else rep(NA_character_, k),
    weights   = if (!is.null(optional$weights))    optional$weights    else rep(NA, k),
    stringsAsFactors = FALSE
  )

  return(list(
    data_outcome       = data_outcome,
    k                  = k,
    mods_from_yi       = mods_from_yi,
    formula_yi         = formula_yi,
    weights_provided   = optional$weights_provided,
    slab_provided      = optional$slab_provided,
    study_ids_provided = optional$study_ids_provided,
    na_check_cols      = c("yi", "sei")  # Only check these columns for NAs
  ))
}


# Internal function to extract and validate outcome variables for GLMM models
# Handles ai, bi, ci, di, n1i, n2i, weights, study_ids, slab
#
# @param .call Matched call from the calling function
# @param data The data frame (can be NULL)
# @param .envir The enclosing environment
#
# @return A list with:
#   - data_outcome: data.frame with ai, bi, n1i, n2i, study_ids, slab, weights
#   - k: number of observations
#   - mods_from_yi: NULL (GLMM does not support formula syntax for outcome)
#   - formula_yi: NULL
#   - slab_provided: logical, whether slab was provided by user
#   - study_ids_provided: logical, whether study_ids was provided by user
#   - na_check_cols: character vector of column names to check for NAs
.check_and_list_data.outcome.glmm <- function(.call, data, .envir) {

  # Extract cell counts for 2x2 tables
  ai  <- .get_variable(.call, data, .envir, "ai",  allow_NULL = TRUE)
  bi  <- .get_variable(.call, data, .envir, "bi",  allow_NULL = TRUE)
  ci  <- .get_variable(.call, data, .envir, "ci",  allow_NULL = TRUE)
  di  <- .get_variable(.call, data, .envir, "di",  allow_NULL = TRUE)
  n1i <- .get_variable(.call, data, .envir, "n1i", allow_NULL = TRUE)
  n2i <- .get_variable(.call, data, .envir, "n2i", allow_NULL = TRUE)

  ### Determine input format and validate

  # Check which variables are provided
  has_abcd <- !is.null(ai) && !is.null(bi) && !is.null(ci) && !is.null(di)
  has_ac_n <- !is.null(ai) && !is.null(ci) && !is.null(n1i) && !is.null(n2i)

  if (!has_abcd && !has_ac_n) {
    stop("For GLMM models, provide either (ai, bi, ci, di) or (ai, ci, n1i, n2i).", call. = FALSE)
  }

  # Determine k from first available variable
  if (!is.null(ai)) {
    k <- length(ai)
  } else {
    stop("The 'ai' argument is required for GLMM models.", call. = FALSE)
  }

  # Validate and convert inputs
  if (has_abcd) {
    # Full 2x2 table specification
    BayesTools::check_int(ai, "ai", check_length = k, allow_NULL = FALSE, allow_NA = TRUE, lower = 0)
    BayesTools::check_int(bi, "bi", check_length = k, allow_NULL = FALSE, allow_NA = TRUE, lower = 0)
    BayesTools::check_int(ci, "ci", check_length = k, allow_NULL = FALSE, allow_NA = TRUE, lower = 0)
    BayesTools::check_int(di, "di", check_length = k, allow_NULL = FALSE, allow_NA = TRUE, lower = 0)

    # Compute marginal totals if not provided
    if (is.null(n1i)) n1i <- ai + bi
    if (is.null(n2i)) n2i <- ci + di

  } else {
    # ai, ci with marginal totals
    BayesTools::check_int(ai,  "ai",  check_length = k, allow_NULL = FALSE, allow_NA = TRUE, lower = 0)
    BayesTools::check_int(ci,  "ci",  check_length = k, allow_NULL = FALSE, allow_NA = TRUE, lower = 0)
    BayesTools::check_int(n1i, "n1i", check_length = k, allow_NULL = FALSE, allow_NA = TRUE, lower = 0)
    BayesTools::check_int(n2i, "n2i", check_length = k, allow_NULL = FALSE, allow_NA = TRUE, lower = 0)

    # Compute bi and di from marginals
    bi <- n1i - ai
    di <- n2i - ci

    # Check for negative values (invalid data)
    if (any(bi < 0, na.rm = TRUE))
      stop("Invalid data: some values of 'bi' (= n1i - ai) are negative.", call. = FALSE)
    if (any(di < 0, na.rm = TRUE))
      stop("Invalid data: some values of 'di' (= n2i - ci) are negative.", call. = FALSE)
  }

  # Extract and validate common optional variables (weights, study_ids, slab)
  optional <- .check_and_list_data.optional_vars(.call, data, .envir, k, "ai")

  ### Construct output data frame
  data_outcome <- data.frame(
    ai        = ai,
    bi        = bi,
    n1i       = n1i,
    n2i       = n2i,
    study_ids = if (!is.null(optional$study_ids)) optional$study_ids else rep(NA_character_, k),
    slab      = if (!is.null(optional$slab))      optional$slab      else rep(NA_character_, k),
    weights   = if (!is.null(optional$weights))   optional$weights   else rep(NA, k),
    stringsAsFactors = FALSE
  )

  return(list(
    data_outcome       = data_outcome,
    k                  = k,
    mods_from_yi       = NULL,  # GLMM does not support formula syntax
    formula_yi         = NULL,
    weights_provided   = optional$weights_provided,
    slab_provided      = optional$slab_provided,
    study_ids_provided = optional$study_ids_provided,
    na_check_cols      = c("ai", "bi", "n1i", "n2i")  # Check cell counts for NAs
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

  # Handle formula input
  if (inherits(mods, "formula")) {

    original_formula <- mods

    # Ensure formula has no LHS (response)
    if (length(mods) == 3) {
      warning(paste0("The '", name, "' formula should not have a left-hand side. ",
                     "The LHS will be ignored."), call. = FALSE)
      mods <- mods[-2]
      original_formula <- mods
    }

    # Create model frame from formula
    if (!is.null(data) && is.data.frame(data)) {
      mf <- try(
        stats::model.frame(mods, data = data, na.action = stats::na.pass),
        silent = TRUE
      )
    } else {
      mf <- try(
        stats::model.frame(mods, na.action = stats::na.pass),
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

  } else if (is.matrix(mods)) {

    mods_df <- as.data.frame(mods)

    if (nrow(mods_df) != k)
      stop(paste0("The number of rows in '", name, "' (", nrow(mods_df),
                  ") must equal the number of effect sizes (", k, ")."), call. = FALSE)

    rownames(mods_df) <- NULL
    attr(mods_df, "formula") <- .create_formula_from_names(names(mods_df))

    return(mods_df)

  } else if (is.data.frame(mods)) {

    mods_df <- mods

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

  return(df)
}


# Internal function to validate predictor variables (mods or scale)
# Checks that:
#   - Continuous variables have sd > 0 (are not constant)
#   - Factor variables have more than one level
#
# @param df A data.frame with predictor variables (can be NULL)
# @param name The name of the predictor type ("mods" or "scale") for error messages
#
# @return NULL (invisibly); stops with error if validation fails
.check_and_list_data.validate_predictors <- function(df, name) {

  if (is.null(df))
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


# Internal helper to create a formula from column names
# Creates a one-sided formula like ~ var1 + var2 + var3
#
# @param col_names Character vector of column names
# @return A formula object
.create_formula_from_names <- function(col_names) {

  if (length(col_names) == 0) {
    return(~ 1)
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
  return(stats::as.formula(formula_str))
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
