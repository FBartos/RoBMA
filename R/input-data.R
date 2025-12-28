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

# Internal function to check and list input data for RoBMA
# Handles flexible input where yi, vi, sei, ni, weights, study_ids, and slab can be passed as:
#   1. A variable name (unquoted)
#   2. A string name (quoted) referring to a column in `data`
#   3. A vector passed directly
#
# @param yi Vector of effect sizes (observed outcomes)
# @param vi Vector of sampling variances (optional if sei is provided)
# @param sei Vector of standard errors (optional if vi is provided)
# @param ni Vector of sample sizes (optional)
# @param weights Vector of weights (optional)
# @param mods Model matrix or formula for moderators (optional)
# @param scale Scale information (optional)
# @param study_ids Vector of study identifiers for multilevel modeling (optional)
# @param data Optional data frame containing the variables
# @param slab Vector of study labels (optional)
# @param subset Logical or numeric vector for subsetting (optional)
# @param .call Matched call from the calling function (required for NSE)
# @param .envir The environment where the calling function was invoked (required for NSE)
#
# @return A list with components:
#   - outcome: data.frame with yi, sei, ni, study_ids, slab, weights
#   - mods: moderator information
#   - scale: scale information
.check_and_list_data <- function(.call, .envir) {

  ### Extract variables using NSE (non-standard evaluation)
  # .call is the match.call() from the parent function (e.g., brma.uni)
  # .envir is the environment where the parent function was called

  # First, extract the data argument - it needs special handling
  # because other variables may reference columns within it
  data <- .get_variable(.call, NULL, .envir, "data", allow_NULL = TRUE)

  # Extract yi (may be a formula like yi ~ mod1 + mod2)
  yi <- .get_variable(.call, data, .envir, "yi", allow_NULL = FALSE)

  # Handle yi as formula (e.g., yi ~ mod1 + mod2)
  # This allows users to specify moderators directly via the yi argument
  formula_yi <- NULL
  mods_from_yi <- NULL

  if (inherits(yi, "formula")) {

    # Store the original formula for later reference
    formula_yi <- yi

    # Check that mods is not also specified (would be ambiguous)
    mods_check <- .get_variable(.call, data, .envir, "mods", allow_NULL = TRUE)
    if (!is.null(mods_check)) {
      stop("Cannot specify 'mods' when 'yi' is a formula. Use either 'yi ~ mod1 + mod2' or 'yi = effect, mods = ~ mod1 + mod2', but not both.", call. = FALSE)
    }

    # Extract the model frame from the formula (RHS becomes mods)
    # Use model.frame to preserve original variable types (not model.matrix which expands factors)
    na_act <- getOption("na.action")
    options(na.action = "na.pass")
    on.exit(options(na.action = na_act), add = TRUE)

    # Get the full model frame (includes response)
    full_mf <- stats::model.frame(yi, data = data)

    # Extract the response (LHS) - first column of model frame
    yi <- stats::model.response(full_mf)
    names(yi) <- NULL

    # Extract predictors (RHS) - remaining columns of model frame
    # Remove the response column to get just the moderators
    if (ncol(full_mf) > 1) {
      mods_from_yi <- full_mf[, -1, drop = FALSE]
    } else {
      # Formula was just ~ 1 (intercept only) or similar
      mods_from_yi <- NULL
    }
  }

  # Extract variance/standard error (one of vi or sei must be provided)
  vi  <- .get_variable(.call, data, .envir, "vi",  allow_NULL = TRUE)
  sei <- .get_variable(.call, data, .envir, "sei", allow_NULL = TRUE)

  # Extract optional variables
  ni        <- .get_variable(.call, data, .envir, "ni",        allow_NULL = TRUE)
  weights   <- .get_variable(.call, data, .envir, "weights",   allow_NULL = TRUE)
  study_ids <- .get_variable(.call, data, .envir, "study_ids", allow_NULL = TRUE)
  slab      <- .get_variable(.call, data, .envir, "slab",      allow_NULL = TRUE)
  subset    <- .get_variable(.call, data, .envir, "subset",    allow_NULL = TRUE)


  ### Input validation

  # Validate yi
  BayesTools::check_real(yi, "yi", check_length = 0, allow_NULL = FALSE, allow_NA = TRUE)
  if (all(is.na(yi)))
    stop("The 'yi' argument must contain at least one non-NA value.", call. = FALSE)

  # Determine number of studies from yi
  k <- length(yi)
  k_original <- k  # Store original k for moderator validation

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

  # Validate weights
  if (!is.null(weights))
    BayesTools::check_real(weights, "weights", check_length = k, allow_NULL = TRUE, allow_NA = TRUE, lower = 0, allow_bound = FALSE)

  # Validate study_ids
  if (!is.null(study_ids) && length(study_ids) != k)
    stop(paste0("The 'study_ids' argument must have length ", k, " (same as 'yi')."), call. = FALSE)

  # Validate slab (study labels)
  if (!is.null(slab)) {
    if (length(slab) != k) {
      stop(paste0("The 'slab' argument must have length ", k, " (same as 'yi')."), call. = FALSE)
    }
  } else {
    # Generate default study labels
    slab <- paste0("Study ", seq_len(k))
  }

  # Validate and apply subset
  if (!is.null(subset)) {
    if (is.logical(subset)) {
      if (length(subset) != k)
        stop(paste0("The 'subset' argument must have length ", k, " (same as 'yi') when logical."), call. = FALSE)
      # Subset is valid logical
    } else if (is.numeric(subset)) {
      if (any(subset < 1) || any(subset > k))
        stop(paste0("The 'subset' argument must contain values between 1 and ", k, "."), call. = FALSE)
      # Convert numeric indices to logical
      subset_logical <- rep(FALSE, k)
      subset_logical[subset] <- TRUE
      subset <- subset_logical
    } else {
      stop("The 'subset' argument must be a logical or numeric vector.", call. = FALSE)
    }

    # Apply subsetting
    yi        <- yi[subset]
    sei       <- sei[subset]
    ni        <- if (!is.null(ni))        ni[subset]        else NULL
    weights   <- if (!is.null(weights))   weights[subset]   else NULL
    study_ids <- if (!is.null(study_ids)) study_ids[subset] else NULL
    slab      <- slab[subset]
    k         <- sum(subset)
  }


  ### Construct output data frame
  data_outcome <- data.frame(
    yi        = yi,
    sei       = sei,
    ni        = if (!is.null(ni))        ni        else rep(NA_integer_, k),
    study_ids = if (!is.null(study_ids)) study_ids else rep(NA_character_, k),
    slab      = slab,
    weights   = if (!is.null(weights))   weights   else rep(NA, k),
    stringsAsFactors = FALSE
  )

  ### Clean outcome input
  # make sure that study ids can be used as indexes in fitting
  if (!is.null(study_ids)) {
    data_outcome[["study_ids"]] <- as.numeric(as.factor(data_outcome[["study_ids"]]))
  }


  ### Handle mods (moderator variables)
  # If yi was a formula, mods have already been extracted from it
  if (!is.null(mods_from_yi)) {

    # Convert model matrix to data.frame, applying subset if needed
    data_mods <- as.data.frame(mods_from_yi)

    # Apply subsetting
    if (!is.null(subset)) {
      data_mods <- data_mods[subset, , drop = FALSE]
      # Drop unused factor levels
      for (j in seq_len(ncol(data_mods))) {
        if (is.factor(data_mods[[j]])) {
          data_mods[[j]] <- droplevels(data_mods[[j]])
        }
      }
    }

    # Create formula attribute from the RHS of the original yi formula
    # (formula_yi[-2] removes the LHS, leaving just ~ rhs)
    attr(data_mods, "formula") <- formula_yi[-2]
    attr(data_mods, "formula_yi") <- formula_yi

  } else {
    # Extract mods from the 'mods' argument using standard extraction
    data_mods <- .extract_moderators(
      .call   = .call,
      data    = data,
      .envir  = .envir,
      name    = "mods",
      k       = k_original,
      subset  = subset
    )
  }

  ### Handle scale (scale moderator variables)
  data_scale <- .extract_moderators(
    .call   = .call,
    data    = data,
    .envir  = .envir,
    name    = "scale",
    k       = k_original,
    subset  = subset
  )


  return(list(
    outcome = data_outcome,
    mods    = data_mods,
    scale   = data_scale
  ))
}


# Internal helper function to extract moderator variables from a formula
# Handles formulas with variables from data frame or environment
# Preserves factor/numeric types and applies subsetting
#
# @param .call The matched call from the parent function
# @param data The data frame to search in (can be NULL)
# @param .envir The enclosing environment for evaluation
# @param name The name of the argument ("mods" or "scale")
# @param k The expected number of rows (before subsetting)
# @param subset Logical vector for subsetting (can be NULL)
# @return A data.frame with moderator variables or NULL
.extract_moderators <- function(.call, data, .envir, name, k, subset) {

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

    # Store the original formula before any modifications
    original_formula <- mods

    # Ensure formula has no LHS (response)
    if (length(mods) == 3) {
      warning(paste0("The '", name, "' formula should not have a left-hand side. ",
                     "The LHS will be ignored."), call. = FALSE)
      mods <- mods[-2]
      original_formula <- mods  # Update to the cleaned version (RHS only)
    }

    # Create model frame from formula
    # model.frame preserves factor levels and handles transformations
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

    # Check length matches yi
    if (nrow(mf) != k)
      stop(paste0("The number of rows in '", name, "' (", nrow(mf),
                  ") must equal the number of effect sizes (", k, ")."), call. = FALSE)

    # Apply subsetting if needed
    if (!is.null(subset)) {
      # Subset while preserving factor levels
      mf <- mf[subset, , drop = FALSE]
      # Drop unused factor levels after subsetting
      mf <- droplevels(mf)
    }

    # Reset row names
    rownames(mf) <- NULL

    # Preserve the original formula (with interactions, transformations, etc.)
    attr(mf, "formula") <- original_formula

    return(mf)

  } else if (is.matrix(mods)) {

    # Convert matrix to data.frame
    mods_df <- as.data.frame(mods)

    # Check length matches yi
    if (nrow(mods_df) != k)
      stop(paste0("The number of rows in '", name, "' (", nrow(mods_df),
                  ") must equal the number of effect sizes (", k, ")."), call. = FALSE)

    # Apply subsetting if needed
    if (!is.null(subset))
      mods_df <- mods_df[subset, , drop = FALSE]

    # Reset row names
    rownames(mods_df) <- NULL

    # Add formula attribute using column names
    attr(mods_df, "formula") <- .create_formula_from_names(names(mods_df))

    return(mods_df)

  } else if (is.data.frame(mods)) {

    # Use data.frame as-is
    mods_df <- mods

    # Check length matches yi
    if (nrow(mods_df) != k)
      stop(paste0("The number of rows in '", name, "' (", nrow(mods_df),
                  ") must equal the number of effect sizes (", k, ")."), call. = FALSE)

    # Apply subsetting if needed
    if (!is.null(subset)) {
      mods_df <- mods_df[subset, , drop = FALSE]
      mods_df <- droplevels(mods_df)
    }

    # Reset row names
    rownames(mods_df) <- NULL

    # Add formula attribute using column names
    attr(mods_df, "formula") <- .create_formula_from_names(names(mods_df))

    return(mods_df)

  } else {
    stop(paste0("The '", name, "' argument must be a formula, matrix, or data.frame."),
         call. = FALSE)
  }
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
