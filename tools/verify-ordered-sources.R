# Focused verification against existing, immutable ordered fits. Run from the
# package directory with a manifest RDS: a named character vector of paths for
# normal, gamma, selected and shared. No fit or cache is created or refreshed.
# ROBMA_TEST_FILES_DIR remains the owner's external-cache control; paths in the
# manifest select the actual existing artifacts, independently of their location.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || !file.exists(args[[1L]])) {
  stop("Supply one existing ordered-fit manifest RDS.", call. = FALSE)
}
fixtures <- readRDS(args[[1L]])
required <- c("normal", "gamma", "selected", "shared")
if (!is.character(fixtures) || !identical(sort(names(fixtures)), sort(required)) ||
    anyNA(fixtures) || any(!file.exists(fixtures))) {
  stop("The manifest must name existing normal, gamma, selected and shared fit paths.", call. = FALSE)
}
options(warn = 1, width = 180)
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
pkgload::load_all(".", compile = FALSE, helpers = FALSE, quiet = TRUE)
ns <- asNamespace("RoBMA")
hashes <- vapply(fixtures, function(path) digest::digest(file = path, algo = "sha256"), character(1))
cat("External cache control: ", Sys.getenv("ROBMA_TEST_FILES_DIR"), "\n", sep = "")
print(lapply(c("BayesTools", "RoBMA"), function(package) {
  path <- getLoadedDLLs()[[package]][["path"]]
  list(package = package, version = as.character(utils::packageVersion(package)),
    namespace = getNamespaceInfo(asNamespace(package), "path"), dll = path,
    dll_sha256 = digest::digest(file = path, algo = "sha256"))
}))
records <- list()
for (name in required) {
  fit <- readRDS(fixtures[[name]])
  context <- ns$.iwmde_context(fit)
  # These are the original two-chain, 2000 samples/chain, twelve-study fixtures.
  # Artifact hashes below also preserve all seeds, adaptation and burn-in state.
  stopifnot(nrow(context[["posterior_samples"]]) == 4000L,
    length(unique(context[["chain_id"]])) == 2L, nrow(context[["data"]][["outcome"]]) == 12L)
  specs <- BayesTools::JAGS_ordered_parameter_spec(fit[["fit"]])
  coefficients <- specs[["mu_g"]][["coefficient_names"]]
  weights <- if (name == "gamma") list(first = stats::setNames(1, coefficients[[1L]]),
    total = stats::setNames(rep(1, 3L), coefficients)) else
    list(hi = stats::setNames(c(1, 1), coefficients[1:2]))
  for (target_name in names(weights)) {
    target <- ns$.iwmde_parameter_spec(context, paste(name, target_name),
      list(type = "linear", weights = weights[[target_name]]))
    state <- ns$.iwmde_ordered_row_state(context, target, 1021L)
    stopifnot(!is.null(state[["source_gaussian"]]))
    source <- ns$.iwmde_ordered_source_normalizer(context, state)
    stopifnot(!is.null(source))
    # The oracle integrates the complete original likelihood/prior kernel on
    # genuinely replayed candidates, independently of completion/CDF/GH formulas.
    physical_log <- function(values) {
      replaced <- ns$.iwmde_ordered_replace_source(context, state, values)
      out <- rep(-Inf, length(values))
      samples <- replaced[["samples"]][replaced[["valid"]], , drop = FALSE]
      if (nrow(samples)) out[replaced[["valid"]]] <-
        ns$.iwmde_log_lik_from_posterior_samples_sum_active_branch(context, samples, state[["active_setup"]], "estimate") +
        state[["prior_evaluator"]](samples) + state[["ordered_kernel"]](samples)
      out
    }
    support <- ns$.iwmde_prior_support(state[["ordered_chart"]][["prior"]])
    scale <- physical_log(state[["ordered_chart"]][["current"]])
    reference <- stats::integrate(function(value) exp(physical_log(value) - scale),
      support[[1L]], support[[2L]], subdivisions = 1000L, rel.tol = 1e-10, abs.tol = 0)
    reference_log <- log(reference[["value"]]) + scale
    difference <- abs(expm1(source[["log_normalizer"]] - reference_log))
    allowance <- source[["discretization_relative_error"]] + source[["exterior_relative_error"]] +
      reference[["abs.error"]] / reference[["value"]] +
      64 * .Machine$double.eps * 12 * source[["max_log_scale"]]
    stopifnot(difference <= allowance)
    record <- data.frame(fixture = name, target = target_name,
      source_family = state[["ordered_chart"]][["prior"]][["distribution"]],
      method = source[["method"]], relative_difference = difference, allowance = allowance)
    print(record)
    records[[paste(name, target_name)]] <- record
    if (target_name == "total") next
    states <- lapply(c(1L, 1021L, 3107L), function(row) ns$.iwmde_ordered_row_state(context, target, row))
    block <- ns$.iwmde_ordered_block_geometry(context, target, states[[2L]])
    candidates <- expand.grid(total = c(-.4, 0, .7), share = c(0, .5, 1))
    samples <- ns$.iwmde_ordered_block_samples(context, block, candidates[["total"]], candidates[["share"]])
    original <- ns$.iwmde_log_lik_from_posterior_samples_sum_active_branch(context, samples, state[["active_setup"]], "estimate")
    compiled <- ns$.iwmde_ordered_block_log_likelihood(context, block, candidates[["total"]], candidates[["share"]])
    floor <- 64 * .Machine$double.eps * 12 * pmax(1, abs(original), abs(compiled))
    stopifnot(all(abs(original - compiled) <= floor))
    produced <- ns$.iwmde_ordered_block_log_terms(context, target, states, target[["ordered_point_plan"]])
    for (proof in produced[["proofs"]]) {
      floating <- 64 * .Machine$double.eps * 12 * max(1, proof[["counts"]][["max_log_scale"]],
        abs(proof[["numerator"]][["log_value"]]), abs(proof[["denominator"]][["log_value"]]))
      e_N <- proof[["numerator"]][["relative_error"]] + expm1(floating)
      e_Z <- proof[["denominator"]][["relative_error"]] + expm1(floating)
      stopifnot(identical(proof[["numerator_relative_error"]], e_N),
        identical(proof[["denominator_relative_error"]], e_Z),
        identical(proof[["ratio_relative_error"]], (1 + e_N) / (1 - e_Z) - 1),
        proof[["ratio_relative_error"]] <= .025)
    }
    cat(name, "original/compiled maximum likelihood difference", max(abs(original - compiled)),
      "; three complete selected-row proofs passed.\n")
    if (!name %in% c("gamma", "shared")) next
    block[["counts"]][["max_log_scale"]] <- 1
    supplier <- target[["ordered_point_plan"]][["supplier"]]
    numerator <- ns$.iwmde_ordered_block_numerator(context, block, supplier)
    original_likelihood <- function(total, share) {
      samples <- ns$.iwmde_ordered_block_samples(context, block, total, share)
      ns$.iwmde_log_lik_from_posterior_samples_sum_active_branch(context, samples, state[["active_setup"]], "estimate")
    }
    scale <- numerator[["log_value"]]
    # Direct reciprocal-supplier integrals use the original likelihood, without
    # the inverse-Gamma or shifted-Beta numerator identities under verification.
    if (supplier == "share") {
      support <- ns$.iwmde_prior_support(block[["prior"]])
      reference <- stats::integrate(function(total) exp(original_likelihood(total, 0) +
        BayesTools::lpdf(block[["prior"]], total) + BayesTools::lpdf(block[["beta"]], 0) -
        log(total) - log(abs(block[["descriptor"]][["k"]])) - scale),
        support[[1L]], support[[2L]], rel.tol = 1e-10, abs.tol = 0, subdivisions = 1000L)
    } else {
      stopifnot(!block[["closed_total_numerator"]])
      reference <- stats::integrate(function(share) exp(original_likelihood(0, share) +
        BayesTools::lpdf(block[["prior"]], 0) + BayesTools::lpdf(block[["beta"]], share) -
        log(share) - log(abs(block[["descriptor"]][["k"]])) - scale),
        0, 1, rel.tol = 1e-10, abs.tol = 0, subdivisions = 1000L)
    }
    difference <- abs(expm1(numerator[["log_value"]] - log(reference[["value"]]) - scale))
    floating <- 64 * .Machine$double.eps * 12 * max(1,
      block[["counts"]][["max_log_scale"]], abs(numerator[["log_value"]]))
    allowance <- numerator[["relative_error"]] + expm1(floating) + reference[["abs.error"]] / reference[["value"]]
    stopifnot(difference <= allowance)
    print(data.frame(fixture = name, supplier = supplier, method = numerator[["method"]],
      relative_difference = difference, allowance = allowance))
  }
}
stopifnot(identical(hashes, vapply(fixtures, function(path) digest::digest(file = path, algo = "sha256"), character(1))))
cat("ORDERED_SOURCE_VERIFICATION_COMPLETE: five physical source normalizers, four original/compiled likelihood grids, twelve selected-row proofs and two physical supplier numerators passed; fit hashes unchanged.\n")
