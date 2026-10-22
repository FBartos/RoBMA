context("Selection node memo")
source(testthat::test_path("common-functions.R"))
skip_on_cran()

# The factor and cluster selection nodes remember the density of their last
# deterministic evaluation (src/distributions/selnorm-node-memo.h). A node may
# answer a call only when every input equals the remembered one, so the key and
# the draws are what these tests pin.

.node_memo_repository_file <- function(...) {

  testthat::test_path("..", "..", ...)
}


.node_memo_compiler <- function() {

  r_bin  <- file.path(R.home("bin"), "R")
  config <- function(name) {
    output <- suppressWarnings(system2(r_bin, c("CMD", "config", name),
                                       stdout = TRUE, stderr = TRUE))
    trimws(output[length(output)])
  }
  list(command = strsplit(config("CXX17"), "[[:space:]]+")[[1L]],
       standard = config("CXX17STD"))
}


test_that("the node memo key covers every input and the map stays bounded", {

  header <- .node_memo_repository_file("src", "distributions", "selnorm-node-memo.h")
  skip_if_not(
    file.exists(header),
    "Repository build sources are not available in this installed-package context."
  )
  compiler <- .node_memo_compiler()
  skip_if(length(compiler$command) == 0L || !nzchar(Sys.which(compiler$command[[1L]])),
          "A C++ compiler is required.")

  work   <- normalizePath(withr::local_tempdir(), winslash = "/", mustWork = TRUE)
  probe  <- file.path(work, if (.Platform$OS.type == "windows") "probe.exe" else "probe")
  build  <- suppressWarnings(system2(
    compiler$command[[1L]],
    c(compiler$command[-1L], compiler$standard, "-O1",
      shQuote(paste0("-I", normalizePath(dirname(header), winslash = "/"))),
      shQuote(testthat::test_path("fixtures", "selnorm-node-memo-probe.cc")),
      "-o", shQuote(probe)),
    stdout = TRUE, stderr = TRUE
  ))
  expect_true(file.exists(probe), info = paste(build, collapse = "\n"))
  output <- system2(probe, stdout = TRUE)
  results <- strsplit(output, " ", fixed = TRUE)
  results <- stats::setNames(as.integer(vapply(results, `[[`, "", 2L)),
                             vapply(results, `[[`, "", 1L))

  # The observations and every parameter but the design of a four-row, rank-two
  # block hold 72 values; each one is mutated in place.
  expect_identical(results[["elements_checked"]], 72L)
  expected <- c("empty_misses", "stored_call_hits", "other_node_misses",
    "density_type_is_input", "type_restored_hits", "mutated_element_misses",
    "store_carries_mutation", "store_replaces_old_value", "restored_call_hits",
    "qmc_design_is_not_an_input", "length_change_misses", "new_layout_hits",
    "old_layout_misses", "limit_keys_all_hit", "update_at_limit_keeps_keys",
    "new_key_beyond_limit_clears")
  expect_setequal(setdiff(names(results), "elements_checked"), expected)
  expect_identical(unname(results[expected]), rep(1L, length(expected)))
})


# The memo switch is read once per process, so each setting fits in a fresh
# worker that loads the RoBMA under test: the checkout when the tests run under
# pkgload, the installed namespace otherwise. The worker's native library must
# be the parent's; a different build would test other code.
.node_memo_fit_in_worker <- function(memo) {

  withr::local_envvar(c(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep),
                        ROBMA_SELNORM_NODE_MEMO = memo))
  cluster <- parallel::makePSOCKcluster(1L, rscript_args = "--vanilla")
  on.exit(parallel::stopCluster(cluster), add = TRUE)
  namespace_path <- getNamespaceInfo(asNamespace("RoBMA"), "path")
  development    <- requireNamespace("pkgload", quietly = TRUE) &&
    pkgload::is_dev_package("RoBMA")
  parent_dll     <- normalizePath(getLoadedDLLs()[["RoBMA"]][["path"]],
                                  winslash = "/", mustWork = TRUE)

  fit_model <- function(paths, development, namespace_path) {

    .libPaths(paths)
    Sys.setenv(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep))
    if (development) {
      pkgload::load_all(namespace_path, compile = FALSE, helpers = FALSE, quiet = TRUE)
    } else {
      loadNamespace("RoBMA", lib.loc = dirname(namespace_path))
    }
    # Five studies of the assink2016 data: four rank-one clusters and one
    # two-factor block, so both selection nodes are evaluated.
    data("dat.assink2016", package = "metadat", envir = environment())
    dat <- dat.assink2016[dat.assink2016$study %in% c(2, 9, 12, 13, 15), ]
    V   <- metafor::vcalc(vi, cluster = study, type = deltype, obs = esid,
                          rho = c(.7, .5), data = dat)
    fit <- suppressWarnings(RoBMA::bselmodel.mv(
      yi = yi, V = V, random = ~ 1 | study / esid, data = dat, measure = "SMD",
      selection = RoBMA::selection_model(estimate_random_effects = "condition"),
      chains = 1L, adapt = 50L, burnin = 100L, sample = 200L, seed = 1L,
      parallel = FALSE, silent = TRUE,
      convergence_checks = RoBMA::set_convergence_checks(max_Rhat = NULL, min_ESS = NULL)
    ))
    list(
      dll     = normalizePath(getLoadedDLLs()[["RoBMA"]][["path"]],
                              winslash = "/", mustWork = TRUE),
      memo    = Sys.getenv("ROBMA_SELNORM_NODE_MEMO", unset = "<unset>"),
      methods = RoBMA:::.data_selection_execution_plan(fit[["data"]])[["block_methods"]],
      draws   = as.data.frame(RoBMA::as_draws_df(fit))
    )
  }
  # Under pkgload this helper descends from the RoBMA namespace; unserializing
  # such a closure would load the installed package in the worker before the
  # build under test, so the function is sent without its defining environment.
  environment(fit_model) <- globalenv()
  result <- parallel::clusterCall(cluster, fit_model, .libPaths(), development,
                                  namespace_path)[[1L]]
  same_build <- if (.Platform$OS.type == "windows") {
    identical(tolower(result[["dll"]]), tolower(parent_dll))
  } else identical(result[["dll"]], parent_dll)
  if (!same_build) {
    stop("The memo worker loaded a different RoBMA build than the tests: worker library '",
         result[["dll"]], "', tests '", parent_dll, "'.", call. = FALSE)
  }

  return(result)
}


test_that("a selection fit with factor and cluster blocks draws identically with the node memo on and off", {

  skip_if_not_installed("rjags")
  skip_if_not_installed("metadat")
  skip_if_not_installed("metafor")
  skip_if_not_installed("pkgload")

  with_memo    <- .node_memo_fit_in_worker("1")
  without_memo <- .node_memo_fit_in_worker("0")

  expect_identical(with_memo[["memo"]], "1")
  expect_identical(without_memo[["memo"]], "0")
  # Both node kinds are evaluated and memoized: the rank-one blocks are cluster
  # nodes, the two-factor block a factor node.
  expect_true(all(c("factor", "rank_one") %in% with_memo[["methods"]]))
  expect_identical(with_memo[["methods"]], without_memo[["methods"]])
  expect_identical(dim(with_memo[["draws"]]), c(200L, ncol(with_memo[["draws"]])))
  expect_identical(with_memo[["draws"]], without_memo[["draws"]])
})


# A step node may remember a density only when its kernel reports the result
# as deterministic: the ladder, the refinement tier or the sparse grid, which
# never read the QMC design that the memo key omits. The batch entries expose the
# flag as the attribute "deterministic" of their result (one value per row).
# A randomized result is never flagged, whatever its Monte Carlo standard error
# (which can be exactly zero), so it cannot be stored.
test_that("the kernels flag only ladder, tier and sparse-grid results as deterministic", {

  flagged <- function(result) attr(result, "deterministic")
  group   <- .test_cert_parallel_group

  # a smooth block: the ladder accepts a rung
  smooth <- .test_cert_block(list(group(4L, c(.6, .5), 1.1, .2)), c(1, .5))
  result <- .test_cert_native_factor(smooth)
  expect_false(is.na(.test_cert_probe_ladder(smooth)))
  expect_identical(flagged(result), TRUE)
  expect_true(.test_cert_deterministic(result))

  # the 80-row block that the shift certificate defers at every rung and the
  # refinement tier certifies
  fed80  <- .test_cert_block(list(group(80L, c(.496441, .542314), .246872, -2.125758)), c(1, .01))
  result <- .test_cert_native_factor(fed80)
  expect_true(is.na(.test_cert_probe_ladder(fed80)))
  expect_true(.test_cert_tier(fed80)[["accepted"]])
  expect_identical(flagged(result), TRUE)
  expect_true(.test_cert_deterministic(result))

  # rank eight with cyclically overlapping supports is not a forest and the node
  # budget affords one tensor rung: the sparse grid certifies it
  rank    <- 8L
  mean    <- seq(-.25, .35, length.out = rank)
  cyclic  <- matrix(0, rank, rank)
  for (column in seq_len(rank)) cyclic[c(column, column %% rank + 1L), column] <- .12
  sparse <- .test_cert_block(lapply(seq_len(rank), function(row) {
    group(1L, cyclic[row, ], .35, mean[row], sei = .25)
  }), c(.4, 1), c(stats::qnorm(.05, lower.tail = FALSE), -Inf))
  result <- .test_cert_native_factor(sparse)
  expect_true(is.na(.test_cert_probe_ladder(sparse)))
  expect_true(.test_cert_deterministic(result))
  expect_identical(flagged(result), TRUE)

  # the rank-one cluster ladder
  cluster <- .test_cert_block(list(group(4L, .8, 1.1, .3)), c(1, .5))
  result  <- .test_cert_native_cluster(cluster)
  expect_true(.test_cert_deterministic(result))
  expect_identical(flagged(result), TRUE)

  # a block forced to the randomized fallback by rules too short to certify it,
  # for the factor and the cluster kernel: flagged FALSE with a Monte Carlo error
  short  <- .selection_joint_cluster_quadrature_rules(c(3L, 5L, 7L))
  result <- .test_cert_native_factor(fed80, rules = short)
  expect_gt(result[["relative_mcse"]], 0)
  expect_identical(flagged(result), FALSE)
  sharp  <- .test_cert_block(list(group(100L, .5, .1, -.1)), c(1, .1))
  result <- .test_cert_native_cluster(sharp, rules = .selection_joint_cluster_quadrature_rules(c(15L, 31L)))
  expect_gt(result[["relative_mcse"]], 0)
  expect_identical(flagged(result), FALSE)

  # one flag per row
  expect_identical(length(flagged(.test_cert_native_factor(smooth))), 1L)
})
