# TEMPORARY diagnostic for the Linux standard-profile failures; removed before merge.
suppressMessages(pkgload::load_all(".", quiet = TRUE))
options(width = 200)
print(La_library()); print(extSoftVersion()["BLAS"])

delta <- function(a, b) {
  a <- unlist(a); b <- unlist(b)
  c(max_abs = max(abs(a - b)), max_rel = max(abs(a - b) / pmax(abs(b), 1e-300)), max_b = max(abs(b)))
}

cat("\n== certificate fixture: actual - stored per field\n")
source("tests/testthat/helper-selection-certificate.R")
reference <- readRDS("tests/testthat/fixtures/selection-certificate-best-reference.rds")
factor_rules  <- .selection_joint_factor_quadrature_rules()
cluster_rules <- .selection_joint_cluster_quadrature_rules(SELNORM_CLUSTER_QUADRATURE_ORDERS)
fields <- c("log_density", "relative_mcse", "relative_change", "log_normalizer")
rows <- NULL
for (i in seq_along(reference[["blocks"]])) {
  b <- reference[["blocks"]][[i]]
  a <- .test_cert_native_best(b, factor_rules, cluster_rules)
  e <- reference[["results"]][[i]]
  d <- vapply(fields, function(f) as.numeric(a[[f]]) - as.numeric(e[[f]]), numeric(1))
  rows <- rbind(rows, data.frame(block = i, rank = b[["rank"]], rmcse = as.numeric(e[["relative_mcse"]]), t(d)))
}
print(rows, digits = 3)
cat("factor nodes digest:", rlang::hash(factor_rules), " cluster:", rlang::hash(cluster_rules), "\n")

cat("\n== block covariance matvec (block - dense)\n")
definitions <- function(file) {
  exprs <- parse(file, keep.source = FALSE)
  for (e in exprs) if (is.call(e) && identical(e[[1L]], as.name("<-"))) eval(e, globalenv())
}
definitions("tests/testthat/test-00-block-covariance.R")
definitions("tests/testthat/test-00-pr82-conditioned-streaming.R")
toy <- .block_covariance_toy()
vector <- { set.seed(1); stats::rnorm(toy[["K"]]) }
for (draw in seq_len(toy[["S"]])) {
  block <- .block_covariance_matvec(toy[["object"]], draw, vector)
  dense <- as.vector(toy[["dense"]][[draw]] %*% vector)
  bound <- 2 * toy[["K"]] * .Machine$double.eps * as.vector(abs(toy[["dense"]][[draw]]) %*% abs(vector))
  print(rbind(diff = block - dense, bound = bound))
}

cat("\n== pr82 dense sampling posterior: chunk 2 vs full\n")
setup <- .pr82_dense_sampling_setup()
posterior <- .selection_conditioned_sampling_posterior(setup)
local({
  withr::local_options(RoBMA.known_v_covariance_max_bytes = 4 * .known_v_covariance_peak_bytes(2L, setup[["K"]]))
  chunked <- .selection_conditioned_sampling_posterior(setup)
  for (n in c("e", "sources", "source_means")) { cat(n, ":\n"); print(delta(chunked[[n]], posterior[[n]])) }
})

cat("\n== pr82 global estimate targets: chunks 1..3 vs full\n")
setup <- .pr82_global_conditioned_setup()
expected <- .selection_conditioned_sampling_estimate_targets(setup, c("log_density", "mean"))
for (chunk_size in 1:3) local({
  withr::local_options(RoBMA.known_v_covariance_max_bytes = 4 * .known_v_covariance_peak_bytes(chunk_size, setup[["K"]]))
  actual <- .selection_conditioned_sampling_estimate_targets(setup, c("log_density", "mean"))
  for (n in names(expected)) { cat("chunk", chunk_size, n, ":\n"); print(delta(actual[[n]], expected[[n]])) }
  print(identical(attributes(actual), attributes(expected)))
})

cat("\n== residual numerics\n")
X <- matrix(c(1, 1, 1, 1, 6.072075e-7, -3.901161e5, 4.130463, 104.9131,
              -5.899230e-6, -1.276173e5, 4.659023e-6, 1270.221), nrow = 4L)
diagonal <- c(4.240359, 2.437391e-6, 1.313608, 6.024016)
rank_one <- c(-67.93884, 844.9597, .02247927, -2779.021)
W  <- .hat_precision_matrix(diagonal, rank_one, list(1:4))
WX <- W %*% X
B  <- .hat_solve_crossprod(crossprod(X, WX))
print(rbind(old = diag(W) - rowSums((WX %*% B) * WX),
            new = .hat_transformed_covariance_diag(W - WX %*% B %*% t(WX), diagonal, rank_one, list(1:4))), digits = 6)

cat("\n== zplot threshold quadrature vs R reduction\n")
definitions("tests/testthat/test-00-zplot-threshold-batch.R")
sei <- c(.3, .5, .3, .3, .8, .5); S <- 6L
mean <- matrix(seq(-.4, .7, length.out = S * length(sei)), S); mean[, 3L] <- mean[, 1L]
sd <- matrix(rep(sei + .2, each = S), S)
latent_sd <- matrix(seq(0, .3, length.out = length(mean)), S); latent_sd[, 3L] <- latent_sd[, 1L]
rule <- list(nodes = c(-3, -.7, 0, .2, 2.5), weights = c(0, .1, .2, .3, .4))
for (direction in c("positive", "negative")) for (normal in list(FALSE, TRUE, c(FALSE, TRUE))) {
  selection <- .zplot_threshold_batch_context(S, sei, direction, normal)
  selection$omega[2L, ] <- selection$omega[2L, ] * 1e280
  selection$omega[3L, ] <- selection$omega[3L, ] * 1e-280
  selection$omega[4L, ] <- c(1, 0, 0)
  selection$omega[5L, ] <- c(.1, .8, .2)
  args <- list(1.96, mean, sd, sei, selection, mean + .125, latent_sd, rule, FALSE)
  a <- do.call(.zplot_selnorm_threshold_quadrature, args)
  b <- do.call(.zplot_threshold_batch_reference, args)
  cat(direction, paste(normal, collapse = "/"), ":", vapply(names(b), function(n) {
    d <- abs(a[[n]] - b[[n]]) / pmax(abs(b[[n]]), 1e-300); d[!is.finite(d)] <- 0
    sprintf("%s %.2e", n, max(d))
  }, ""), "\n")
}

cat("\n== streaming with matprod = 'internal'\n")
local({
  withr::local_options(matprod = "internal")
  setup <- .pr82_global_conditioned_setup()
  expected <- .selection_conditioned_sampling_estimate_targets(setup, c("log_density", "mean"))
  for (chunk_size in 1:3) local({
    withr::local_options(RoBMA.known_v_covariance_max_bytes = 4 * .known_v_covariance_peak_bytes(chunk_size, setup[["K"]]))
    actual <- .selection_conditioned_sampling_estimate_targets(setup, c("log_density", "mean"))
    cat("chunk", chunk_size, "identical:", identical(actual, expected), "\n")
  })
})

cat("\n== qCMDE plot sample scales (bcg_meta-regression3b, alloc:year)\n")
try(local({
  Sys.setenv(NOT_CRAN = "true")
  testthat::source_test_helpers("tests/testthat", env = globalenv())
  record <- list()
  original <- .plot_brma_same_sample_scale
  testthat::local_mocked_bindings(.plot_brma_same_sample_scale = function(raw_samples, plotted_samples) {
    a <- as.numeric(raw_samples); b <- as.numeric(plotted_samples)
    record[[length(record) + 1L]] <<- c(n = length(a), max_abs = max(abs(a - b)), scale = max(abs(b)),
                                        identical = identical(a, b))
    original(raw_samples, plotted_samples)
  }, .package = "RoBMA")
  .local_mock_iwmde_estimate_success()
  testthat::local_mocked_bindings(plot_posterior = function(samples, parameter, ...) structure(list(), class = "mock_plot"),
                                  .package = "BayesTools")
  for (case in list(c("bcg_meta-regression3b", "alloc:year"), c("bcg_meta-regression4", "alloc:year_before1969"),
                    c("dat.lehmann2018_RoBMA_mods2", "Preregistered:Gender"))) {
    out <- try(plot(load_fit(case[[1]], validate = FALSE), parameter_mods = case[[2]], standardized_coefficients = TRUE,
                    plot_type = "ggplot", density_method = "qCMDE", density_control = list(n_points = 20, samples = 20)))
    cat(case, if (inherits(out, "try-error")) "ERROR" else "ok", "\n")
  }
  print(do.call(rbind, record))
}))
cat("\nCPU:\n"); try(cat(grep("model name|^flags", readLines("/proc/cpuinfo", n = 30), value = TRUE)[1:2], sep = "\n"))
for (v in c("CXX17", "CXX17FLAGS", "CXXFLAGS")) try(print(system2("R", c("CMD", "config", v), stdout = TRUE)))
