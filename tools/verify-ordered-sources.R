# Focused verification against existing, immutable ordered fits. Run from the
# package directory with a manifest RDS: a named character vector of paths for
# normal, gamma, selected and shared. No fit or cache is created or refreshed.
# Manifest paths select owner-managed artifacts wherever they already reside.
# ROBMA_TEST_FILES_DIR is printed context only; it does not select these paths.

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
cat("ROBMA_TEST_FILES_DIR (context only; manifest paths select fits): ", Sys.getenv("ROBMA_TEST_FILES_DIR"), "\n", sep = "")
print(lapply(c("BayesTools", "RoBMA"), function(package) {
  path <- getLoadedDLLs()[[package]][["path"]]
  list(package = package, version = as.character(utils::packageVersion(package)),
    namespace = getNamespaceInfo(asNamespace(package), "path"), dll = path,
    dll_sha256 = digest::digest(file = path, algo = "sha256"))
}))
records <- list()
physical_blocks <- list()
physical_block_reference <- function(fit, context, target, state, proof) {
  # The public projection tensor and allocation metadata determine the physical
  # coordinates. No completed Gaussian geometry, production block likelihood,
  # GH rule or probability/CDF-panel formula enters this reference.
  specs <- BayesTools::JAGS_ordered_parameter_spec(fit[["fit"]])
  row <- state[["row"]]
  original <- matrix(row, 1L, dimnames = list(NULL, names(row)))
  projection <- BayesTools::JAGS_ordered_parameter_spec(fit[["fit"]],
    parameters = names(specs), weights = target[["weights"]], draws = original)
  terms <- Filter(function(term) any(term[["tensor"]][["coefficients"]] != 0), projection[["contractions"]])
  stopifnot(length(terms) == 1L, length(projection[["ordinary_weights"]]) == 0L)
  term <- terms[[1L]]
  tensor <- term[["tensor"]]
  stopifnot(length(tensor[["records"]]) == 1L)
  allocation <- tensor[["records"]][[1L]]
  coefficients <- tensor[["coefficients"]]
  J <- which(coefficients != 0)
  K <- which(coefficients == 0)
  slope <- unique(coefficients[J])
  stopifnot(length(slope) == 1L, is.finite(slope), slope != 0,
    length(J) > 0L, length(K) > 0L, allocation[["spec"]][["type"]] == "dirichlet")
  alpha <- allocation[["spec"]][["alpha"]]
  a <- sum(alpha[J])
  b <- sum(alpha[K])
  owner <- specs[[term[["parameter"]]]]
  prior <- owner[["total_prior"]]
  source <- owner[["total_names"]][[term[["slice"]]]]
  support <- ns$.iwmde_prior_support(prior)
  gamma <- allocation[["gamma_coordinates"]]
  H <- sum(row[gamma])
  r_J <- row[gamma[J]] / sum(row[gamma[J]])
  r_K <- row[gamma[K]] / sum(row[gamma[K]])
  stopifnot(is.finite(H), H > 0, all(is.finite(c(r_J, r_K))))
  replay <- function(total, share) {
    n <- max(length(total), length(share))
    total <- rep_len(total, n)
    share <- rep_len(share, n)
    samples <- matrix(rep(row, each = n), nrow = n, dimnames = list(NULL, names(row)))
    samples[, source] <- total
    # A shared allocation is changed once. Replay then updates every owner,
    # retaining the original within-group ratios and gamma sum H.
    samples[, gamma[J]] <- outer(H * share, r_J)
    samples[, gamma[K]] <- outer(H * (1 - share), r_K)
    ns$.iwmde_ordered_replay(context, samples, specs)
  }
  candidates <- replay(c(.4, .7, 1), c(.2, .5, .8))
  projected <- BayesTools::JAGS_ordered_parameter_spec(fit[["fit"]],
    parameters = names(specs), weights = target[["weights"]], draws = candidates)[["values"]]
  stopifnot(all(abs(projected - slope * c(.4, .7, 1) * c(.2, .5, .8)) <=
    64 * .Machine$double.eps * pmax(1, abs(projected))))
  likelihood <- function(total, share) {
    samples <- replay(total, share)
    value <- ns$.iwmde_log_lik_from_posterior_samples_sum_active_branch(
      context, samples, state[["active_setup"]], "estimate")
    stopifnot(length(value) == nrow(samples), all(is.finite(value)))
    value
  }
  # These positive scaled integrals retain finite log values even when their
  # unscaled physical values are outside double range. A zero reference never
  # passes as a physical zero. QUADPACK errors are explicit numerical estimates.
  scale_D <- proof[["denominator"]][["log_value"]]
  maximum_inner_error <- 0
  inner <- function(share) vapply(share, function(S) {
    value <- stats::integrate(function(T) exp(likelihood(T, S) +
      BayesTools::lpdf(prior, T) - scale_D), support[[1L]], support[[2L]],
      subdivisions = 1000L, rel.tol = 1e-9, abs.tol = 1e-11)
    stopifnot(is.finite(value[["value"]]), value[["value"]] > 0,
      is.finite(value[["abs.error"]]))
    maximum_inner_error <<- max(maximum_inner_error, value[["abs.error"]])
    value[["value"]]
  }, numeric(1))
  denominator <- stats::integrate(function(S) stats::dbeta(S, a, b) * inner(S),
    0, 1, subdivisions = 1000L, rel.tol = 1e-9, abs.tol = 1e-11)
  stopifnot(is.finite(denominator[["value"]]), denominator[["value"]] > 0)
  reference_D_error <- (denominator[["abs.error"]] + maximum_inner_error) / denominator[["value"]]
  log_D <- log(denominator[["value"]]) + scale_D
  scale_N <- proof[["numerator"]][["log_value"]]
  if (identical(prior[["distribution"]], "gamma")) {
    stopifnot(a == 1, support[[1L]] == 0, is.infinite(support[[2L]]),
      prior[["parameters"]][["shape"]] == 2)
    numerator <- stats::integrate(function(T) exp(likelihood(T, 0) + BayesTools::lpdf(prior, T) +
      stats::dbeta(0, a, b, log = TRUE) - log(T) - log(abs(slope)) - scale_N),
      support[[1L]], support[[2L]], subdivisions = 1000L, rel.tol = 1e-9, abs.tol = 1e-11)
  } else {
    stopifnot(identical(prior[["distribution"]], "normal"), a > 1)
    # Shared interactions still depend on S at T_main=0. Evaluate every S.
    numerator <- stats::integrate(function(S) exp(likelihood(0, S) + BayesTools::lpdf(prior, 0) +
      stats::dbeta(S, a, b, log = TRUE) - log(S) - log(abs(slope)) - scale_N),
      0, 1, subdivisions = 1000L, rel.tol = 1e-9, abs.tol = 1e-11)
  }
  stopifnot(is.finite(numerator[["value"]]), numerator[["value"]] > 0,
    is.finite(numerator[["abs.error"]]))
  reference_N_error <- numerator[["abs.error"]] / numerator[["value"]]
  log_N <- log(numerator[["value"]]) + scale_N
  difference_N <- abs(expm1(proof[["numerator"]][["log_value"]] - log_N))
  difference_D <- abs(expm1(proof[["denominator"]][["log_value"]] - log_D))
  difference_ratio <- abs(expm1(proof[["numerator"]][["log_value"]] -
    proof[["denominator"]][["log_value"]] - log_N + log_D))
  floating <- 64 * .Machine$double.eps * nrow(context[["data"]][["outcome"]]) *
    max(1, abs(c(log_N, log_D)), proof[["counts"]][["max_log_scale"]])
  allowance_N <- proof[["numerator_relative_error"]] + reference_N_error + expm1(floating)
  allowance_D <- proof[["denominator_relative_error"]] + reference_D_error + expm1(floating)
  reference_ratio_error <- (1 + reference_N_error) / (1 - reference_D_error) - 1
  allowance_ratio <- proof[["ratio_relative_error"]] + reference_ratio_error + expm1(floating)
  stopifnot(reference_D_error < 1, difference_N <= allowance_N,
    difference_D <= allowance_D, difference_ratio <= allowance_ratio)
  data.frame(row = state[["row_index"]], a = a, b = b, slope = slope,
    log_N = log_N, log_D = log_D, difference_N = difference_N, allowance_N = allowance_N,
    difference_D = difference_D, allowance_D = allowance_D,
    difference_ratio = difference_ratio, allowance_ratio = allowance_ratio,
    reference_N_error = reference_N_error, reference_D_error = reference_D_error)
}
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
    states <- lapply(c(1L, 1990L, 4000L), function(row) ns$.iwmde_ordered_row_state(context, target, row))
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
    for (i in seq_along(states)) {
      reference <- physical_block_reference(fit, context, target, states[[i]], produced[["proofs"]][[i]])
      reference <- cbind(fixture = name, reference)
      print(reference)
      physical_blocks[[paste(name, i)]] <- reference
    }
    cat(name, "original/compiled maximum likelihood difference", max(abs(original - compiled)),
      "; three N/D ratio arithmetic proofs and independent physical block references passed.\n")
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
cat("ORDERED_SOURCE_VERIFICATION_COMPLETE: five physical source normalizers, four original/compiled likelihood grids, twelve N/D ratio arithmetic proofs, twelve independent physical block N/D and ratio references and two physical supplier numerators passed; fit hashes unchanged.\n")
cat("LIMIT: selective original-source correctness evidence on rows 1/1990/4000; reference quadrature errors are estimates, not global integration or chain-calibration certification. No new native build or performance certification.\n")
