# ============================================================================ #
# bias-branch-design.R
# ============================================================================ #
#
# Design classes of bias-term ensembles.
#
# .get_model_matrix() returns one union design for a fitted object: the location
# columns plus a PET column (sei) and a PEESE column (sei^2) whenever any bias
# component contains that term. A RoBMA() ensemble mixes members with different
# bias terms, and no member contains both. The design-based diagnostics
# (hat values and everything built on them, VIF/GVIF, the typical sampling
# variance of the heterogeneity summary) therefore use, for each posterior
# draw, the design of the member that the draw belongs to. Members with the same
# design form one class:
#
#   1  location columns only (no bias component, or a weight function),
#   2  location columns and the PET column,
#   3  location columns and the PEESE column.
#
# ============================================================================ #


# ---------------------------------------------------------------------------- #
# .bias_design_member_classes
# ---------------------------------------------------------------------------- #
#
# Design class (1, 2, or 3) of every member of the bias prior, in the order of
# the member index that `bias_indicator` refers to.
#
# ---------------------------------------------------------------------------- #
.bias_design_member_classes <- function(priors) {

  members <- .selection_bias_priors(priors)

  return(vapply(members, function(prior) {
    if (BayesTools::is.prior.PET(prior)) {
      return(2L)
    }
    if (BayesTools::is.prior.PEESE(prior)) {
      return(3L)
    }
    return(1L)
  }, integer(1L)))
}


# ---------------------------------------------------------------------------- #
# .bias_class_design
# ---------------------------------------------------------------------------- #
#
# The union design restricted to one class: the location columns and the bias
# column of the class. The attributes follow what .get_model_matrix() returns
# for a model whose bias prior holds only that term.
#
# @param X union design from .get_model_matrix().
# @param is_PET,is_PEESE whether the union design carries the PET / PEESE column.
# @param design_class design class (1, 2, or 3).
#
# ---------------------------------------------------------------------------- #
.bias_class_design <- function(X, is_PET, is_PEESE, design_class) {

  bias_parameters <- c(if (is_PET) "PET", if (is_PEESE) "PEESE")
  n_location      <- ncol(X) - length(bias_parameters)
  bias_parameter  <- c(NA_character_, "PET", "PEESE")[design_class]
  bias_column     <- if (is.na(bias_parameter)) {
    integer()
  } else {
    n_location + match(bias_parameter, bias_parameters)
  }
  keep <- c(seq_len(n_location), bias_column)

  design <- X[, keep, drop = FALSE]

  assign <- attr(X, "assign")
  if (is.null(assign)) {
    assign <- seq_len(ncol(X)) - 1L
  }
  assign <- assign[seq_len(n_location)]
  if (length(bias_column) > 0L) {
    assign <- c(assign, max(c(0L, assign)) + 1L)
  }
  attr(design, "assign") <- assign

  term_labels <- attr(X, "term_labels")
  if (!is.null(term_labels)) {
    term_labels <- term_labels[seq_len(length(term_labels) - length(bias_parameters))]
    if (length(bias_column) > 0L) {
      term_labels <- c(term_labels, bias_parameter)
    }
    if (length(term_labels) > 0L) {
      attr(design, "term_labels") <- term_labels
    }
  }

  raw_colnames <- attr(X, "raw_colnames")
  if (!is.null(raw_colnames) && length(raw_colnames) == ncol(X)) {
    attr(design, "raw_colnames") <- raw_colnames[keep]
  }
  aliased <- attr(X, "aliased")
  if (!is.null(aliased) && length(aliased) == ncol(X)) {
    attr(design, "aliased") <- aliased[keep]
  }
  attr(design, "rank") <- qr(design)[["rank"]]

  return(design)
}


# ---------------------------------------------------------------------------- #
# .bias_design_classes
# ---------------------------------------------------------------------------- #
#
# Design class of every posterior draw of an ensemble whose members differ in
# their bias terms.
#
# @param object brma object.
# @param X union design from .get_model_matrix().
# @param posterior_samples optional posterior sample matrix; the rows of the
#   returned class vector are the rows of these samples (e.g. after known-V
#   thinning), or of the full posterior when omitted.
#
# @return NULL when one design applies to every draw (no bias term, a single
#   bias prior, or members of one class), which leaves the union design as it
#   is; otherwise a list with
#   - class: integer vector, the design class of each posterior row,
#   - designs: list of the three class designs (NULL for classes without a
#     member),
#   - n_draws: draws per class.
#
# ---------------------------------------------------------------------------- #
.bias_design_classes <- function(object, X = .get_model_matrix(object),
                                 posterior_samples = NULL) {

  member_class <- .bias_design_member_classes(object[["priors"]])
  if (length(unique(member_class)) < 2L) {
    return(NULL)
  }

  bias_indicator <- .extract_bias_indicator(
    object,
    posterior_samples = posterior_samples
  )
  if (any(is.na(bias_indicator)) ||
      any(bias_indicator < 1L | bias_indicator > length(member_class))) {
    stop("Invalid 'bias_indicator' values in posterior samples.", call. = FALSE)
  }
  draw_class <- member_class[bias_indicator]

  is_PET   <- .is_PET(object)
  is_PEESE <- .is_PEESE(object)
  designs  <- vector("list", 3L)
  for (design_class in unique(member_class)) {
    designs[[design_class]] <- .bias_class_design(
      X, is_PET, is_PEESE, design_class
    )
  }

  return(list(
    class   = draw_class,
    designs = designs,
    n_draws = tabulate(draw_class, nbins = 3L)
  ))
}
