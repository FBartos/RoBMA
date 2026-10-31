# Role and coverage guards for the fixed four-fixture ordered-source verifier.
# Sourcing this file does not load a product, read a fit or perform integration.
.ordered_source_fixture_role <- function(name, fit, context, specs, ns) {

  guard <- function(valid, detail) {
    if (!isTRUE(valid)) stop("Ordered fixture '", name, "' role: ", detail, ".", call. = FALSE)
  }
  guard(name %in% c("normal", "gamma", "selected", "shared"), "unknown fixture")
  guard(length(specs) > 0L, "a declared ordered owner is required")
  prior <- specs[[1L]][["total_prior"]]
  family <- if (name == "gamma") "gamma" else "normal"
  guard(identical(prior[["distribution"]], family), paste("declared", family, "total prior is required"))
  selected <- name %in% c("selected", "shared")
  declaration <- ns$.iwmde_selection_spec(context[["data"]], fit[["priors"]])
  guard(if (selected) !is.null(declaration) && identical(declaration[["mode"]], "step") else is.null(declaration),
    if (selected) "a declared step-selection law is required" else "an unselected declaration is required")
  if (name != "shared") return(invisible(NULL))

  records <- specs[[1L]][["allocations"]]
  shared <- lapply(records, function(record) {
    owners <- Filter(function(owner) {
      allocation <- owner[["allocations"]][[record[["key"]]]]
      total <- owner[["total_prior"]]
      live <- !(BayesTools::is.prior.point(total) && identical(total[["parameters"]][["location"]], 0))
      live && !is.null(allocation) && identical(allocation[["key"]], record[["key"]]) &&
        identical(allocation[["gamma_coordinates"]], record[["gamma_coordinates"]])
    }, specs)
    list(key = record[["key"]], gamma = record[["gamma_coordinates"]], owners = names(owners))
  })
  shared <- Filter(function(record) length(record[["gamma"]]) > 0L && length(record[["owners"]]) > 1L, shared)
  guard(length(shared) > 0L, "matching allocation key and Gamma coordinates must have multiple live declared owners")
  invisible(shared)
}

.ordered_source_block_role <- function(name, context, target, state, block, ns) {

  guard <- function(valid, detail) {
    if (!isTRUE(valid)) stop("Ordered fixture '", name, "', row ", state[["row_index"]], " role: ", detail, ".", call. = FALSE)
  }
  selected <- name %in% c("selected", "shared")
  guard(identical(isTRUE(state[["active_setup"]][["is_weightfunction"]]), selected), "the actual row-state selection route must match")
  guard(identical(isTRUE(block[["selected"]]), selected) &&
      identical(isTRUE(block[["state"]][["active_setup"]][["is_weightfunction"]]), selected),
    "the actual block selection route must match")
  if (selected) guard(!is.null(block[["selection_context"]]) &&
    all(block[["selection_context"]][["kernel_mode"]] == ns$SELKERNEL_STEP), "the actual selected block must use the step kernel")
  if (name != "shared") return(invisible(NULL))

  guard(identical(target[["ordered_point_plan"]][["supplier"]], "total") &&
    identical(block[["closed_total_numerator"]], FALSE), "a nonclosed total supplier is required")
  samples <- ns$.iwmde_ordered_block_samples(context, block, rep(0, 3L), c(.2, .5, .8))
  likelihood <- ns$.iwmde_log_lik_from_posterior_samples_sum_active_branch(
    context, samples, block[["state"]][["active_setup"]], "estimate")
  floor <- 64 * .Machine$double.eps * nrow(context[["data"]][["outcome"]]) * max(1, abs(likelihood))
  guard(length(likelihood) == 3L && all(is.finite(likelihood)) && diff(range(likelihood)) > floor,
    "the physical likelihood must retain a shared-allocation contribution at total zero")
  invisible(likelihood)
}

.ordered_source_verification_counts <- function(counts) {

  required <- c(source = 5L, grid = 4L, proof = 12L, physical_block = 12L, supplier = 2L)
  if (!identical(counts, required)) stop("Ordered verification coverage differs from the fixed target schedule: ",
    paste(names(counts), counts, sep = "=", collapse = ", "), ".", call. = FALSE)
  invisible(counts)
}


# Error estimates are relative to their own computed values P and R. Compare
# on the reference scale with r = P/R, keeping the existing accuracy ceiling.
.ordered_source_comparison <- function(log_producer, log_reference,
                                       producer_error, reference_error,
                                       floating_error, fixture, row, quantity) {

  fail <- function(detail) {
    stop("Ordered verification for fixture '", fixture, "', row ", row,
         ", ", quantity, ": ", detail, ".", call. = FALSE)
  }
  inputs <- list(log_producer, log_reference, producer_error, reference_error,
                 floating_error)
  if (!all(vapply(inputs, function(value) {
    is.numeric(value) && length(value) == 1L && is.finite(value)
  }, logical(1)))) {
    fail("comparison values and error estimates must be finite scalars")
  }
  if (any(c(producer_error, reference_error, floating_error) < 0)) {
    fail("error estimates must be nonnegative")
  }
  ratio <- exp(log_producer - log_reference)
  difference <- abs(expm1(log_producer - log_reference))
  allowance <- producer_error * ratio + reference_error +
    floating_error * max(1, ratio)
  if (!is.finite(ratio) || !is.finite(difference) || !is.finite(allowance)) {
    fail("discrepancy and comparison envelope must be finite")
  }
  if (allowance > .025) {
    fail("comparison envelope exceeds the 0.025 accuracy ceiling")
  }
  if (difference > allowance) {
    fail("discrepancy exceeds its comparison envelope")
  }
  return(list(difference = difference, allowance = allowance))
}


.ordered_source_reference_integral <- function(reference, fixture, row, quantity) {

  value <- reference[["value"]]
  error <- reference[["abs.error"]]
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) || value <= 0 ||
      !is.numeric(error) || length(error) != 1L || !is.finite(error) || error < 0) {
    stop("Ordered verification for fixture '", fixture, "', row ", row,
         ", ", quantity, ": reference integral must be finite and positive with a finite nonnegative error estimate.",
         call. = FALSE)
  }
  return(invisible(reference))
}
