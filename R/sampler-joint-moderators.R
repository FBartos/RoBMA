# Bind the fixed offset to the same coefficient coordinates and design
# columns that the compiler uses in the fitted location predictor.
.sampler_joint_bias_moderators <- function(prior_list, formula_design, data) {

  if (!identical(names(formula_design), "mu")) return(NULL)
  design <- formula_design[["mu"]]
  terms  <- design[["terms"]]
  labels <- attr(terms, "term.labels")
  if (length(labels) == 0L || !all(attr(terms, "order") == 1L) ||
      !all(vapply(as.list(attr(terms, "variables"))[-1L], is.name, logical(1))) ||
      length(attr(terms, "offset")) > 0L || isTRUE(design[["log_intercept"]]) ||
      length(design[["transformed_terms"]]) > 0L || length(design[["expression_specs"]]) > 0L ||
      length(design[["random_effects"]]) > 0L) return(NULL)

  map <- design[["name_map"]]
  fixed <- map[map[["kind"]] == "fixed" & map[["role"]] == "coefficient", , drop = FALSE]
  term_nodes <- fixed[["jags_name"]][match(c("intercept", labels), fixed[["term"]])]
  if (anyNA(term_nodes) || anyDuplicated(term_nodes) || !all(term_nodes %in% names(prior_list))) return(NULL)

  X      <- design[["model_matrix"]]
  assign <- design[["assign"]]
  n      <- data[["K"]]
  rows   <- data[["sel_joint_singleton_row"]]
  if (!is.matrix(X) || nrow(X) != n || length(rows) != n || !setequal(rows, seq_len(n)) ||
      length(assign) != ncol(X) || sum(assign == 0L) != 1L ||
      !setequal(assign, seq.int(0L, length(labels)))) return(NULL)

  references <- character(ncol(X))
  for (i in seq_along(labels)) {
    term <- labels[[i]]
    type <- design[["model_terms_type"]][[term]]
    if (!isTRUE(type %in% c("continuous", "factor"))) return(NULL)
    name <- term_nodes[[i + 1L]]
    if (!.sampler_joint_bias_moderator_prior(prior_list[[name]], type)) return(NULL)
    data_name <- design[["jags_data_names"]][[term]]
    if (length(data_name) != 1L) return(NULL)
    emitted <- data[[data_name]]
    columns <- which(assign == i)
    width   <- length(columns)
    if (type == "continuous") {
      if (width != 1L || length(emitted) != n || !is.null(dim(emitted))) return(NULL)
    } else if (!is.matrix(emitted) || nrow(emitted) != n || ncol(emitted) != width) {
      return(NULL)
    }
    if (!is.numeric(emitted) || !all(is.finite(emitted)) ||
        !identical(as.numeric(emitted), as.numeric(X[, columns, drop = FALSE]))) return(NULL)

    references[columns] <- if (width == 1L) name else paste0(name, "[", seq_len(width), "]")
  }

  moderators <- which(assign != 0L)
  return(list(intercept = term_nodes[[1L]], prior_names = term_nodes, coefficients = references[moderators],
              design = as.numeric(t(X[rows, moderators, drop = FALSE]))))
}


.sampler_joint_bias_moderator_prior <- function(prior, type) {

  if (is.null(attr(prior, "multiply_by", exact = TRUE)) &&
      BayesTools::is.prior.mixture(prior) && !BayesTools::is.prior.spike_and_slab(prior)) {
    return(all(vapply(prior, .sampler_joint_bias_moderator_prior, logical(1), type = type)))
  }
  if (type == "continuous") return(.sampler_joint_scalar_prior(prior))

  # The compiler takes the factor dimension from the bound design. The
  # underlying vector prior can still carry an unfilled K placeholder.
  parameters <- prior[["parameters"]]
  parameters[["K"]] <- NULL
  return(BayesTools::is.prior.factor(prior) && !BayesTools::is.prior.ordered(prior) &&
           !BayesTools::is.prior.R2(prior) && !BayesTools::is.prior.spike_and_slab(prior) &&
           !BayesTools::is.prior.regularized(prior) && !BayesTools::is.prior.regularized_term(prior) &&
           is.null(attr(prior, "multiply_by", exact = TRUE)) &&
           all(vapply(parameters, function(value) {
             is.numeric(value) && all(is.finite(value))
           }, logical(1))))
}
