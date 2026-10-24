context("Coarse surrogates of cluster and factor selection blocks")

# The corrected sampler (selection.sampler = "coarse_corrected") draws from the
# exact full conditional through a surrogate: for a cluster or factor selection
# block the exact numerator (Gaussian log density plus observed log weights)
# minus the log normalizer of rule j = min(selection.coarse_max_rules, rules) of
# the block's own Gauss-Hermite schedule at the actual parameters. These tests
# pin the surrogate value against independent base-R Gauss-Hermite rules, its
# fallbacks, its purity, and the validators the sampler splits into dynamic and
# control checks. They go through RoBMA_selnorm_surrogate_probe, which calls the
# distribution classes the sampler calls (helper-selection-coarse-blocks.R).

.coarse_rules_full <- function() .coarse_rules(c(7, 15, 31, 63, 127, 255, 511, 1023))
.coarse_rules_factor <- function() .coarse_rules(c(3, 5, 7, 9, 11, 15, 21, 31, 43, 63, 95))

# Equicorrelated rho = .8 cluster: residual SD sqrt(1 - rho), loading sqrt(rho).
.coarse_cluster_block <- function(...) {

  rho <- .8
  arguments <- list(kind = "cluster", mean = c(.1, -.2, .3, 0),
    residual_sd = rep(sqrt(1 - rho), 4), loading = rep(sqrt(rho), 4), sei = rep(1, 4),
    omega = c(1, .3), x = c(.2, -.5, .1, .4), obs_bin = c(1L, 2L, 1L, 2L),
    rules = .coarse_rules_full())
  do.call(.coarse_block, utils::modifyList(arguments, list(...)))
}

# Two latent axes, the first loading on every row and the second on rows 3-5.
.coarse_factor_block <- function(...) {

  arguments <- list(kind = "factor", mean = c(.1, -.2, .3, 0, .2), residual_sd = rep(.7, 5),
    loading = cbind(c(.5, .5, .4, .3, .6), c(0, 0, .4, .5, .3)), sei = rep(1, 5),
    omega = c(1, .3), x = c(.2, -.5, .1, .4, 0), obs_bin = c(1L, 2L, 1L, 2L, 1L),
    rules = .coarse_rules_factor())
  do.call(.coarse_block, utils::modifyList(arguments, list(...)))
}

# The log normalizer and numerator of a probe against the independent evaluation
# of the rule of order `order` (log scale, absolute).
.expect_coarse_value <- function(probe, block, order, tolerance = 1e-12, info = NULL) {

  reference_a <- .coarse_reference_log_normalizer(block, order)
  reference_n <- .coarse_reference_numerator(block)
  expect_false(probe[["exact_fallback"]], info = info)
  expect_lte(abs(probe[["log_normalizer"]] - reference_a), tolerance, label = paste(info, "log A_j"))
  expect_lte(abs(probe[["numerator"]] - reference_n), tolerance, label = paste(info, "numerator"))
  expect_lte(abs(probe[["surrogate"]] - (reference_n - reference_a)), tolerance,
             label = paste(info, "surrogate"))
  expect_lte(abs((probe[["numerator"]] - probe[["surrogate"]]) - reference_a), tolerance,
             label = paste(info, "numerator - surrogate"))
}


test_that("the cluster surrogate is the exact numerator minus the log normalizer of rule j", {

  block <- .coarse_cluster_block()
  for (max_rules in c(3L, 4L, 5L, 100L)) {
    probe <- .coarse_probe(block, max_rules)
    expect_identical(probe[["rule"]], min(max_rules, 8L))
    .expect_coarse_value(probe, block, .coarse_rule_order(block, max_rules),
                         info = paste("max_rules", max_rules))
    expect_identical(probe[["surrogate_error"]], "")
  }
  # a negative effect direction and three bins
  negative <- .coarse_cluster_block(sign = -1, z_lower = c(1.96, 0, -1e300),
    z_upper = c(1e300, 1.96, 0), omega = c(1, .6, .2), obs_bin = c(1L, 2L, 3L, 3L))
  .expect_coarse_value(.coarse_probe(negative, 3L), negative, 31L, info = "negative effect, three bins")
  # a crude rule 3 (three nodes): the surrogate is far from the exact density
  # and still equals its definition
  crude <- .coarse_cluster_block(rules = .coarse_rules(c(1, 2, 3, 7, 15, 31, 63, 127, 255, 511)))
  probe <- .coarse_probe(crude, 3L)
  .expect_coarse_value(probe, crude, 3L, info = "crude rule 3")
  expect_gt(abs(probe[["surrogate"]] - probe[["exact"]]), 1e-3)
  # a schedule shorter than the setting uses its last rule
  short <- .coarse_cluster_block(rules = .coarse_rules(c(7, 15)))
  probe <- .coarse_probe(short, 3L)
  expect_identical(probe[["rule"]], 2L)
  .expect_coarse_value(probe, short, 15L, info = "two-rule schedule")
})


test_that("the factor surrogate follows rule j of tensors, compacted columns and nested forests", {

  rank_two <- .coarse_factor_block()
  for (max_rules in c(3L, 4L, 6L)) {
    .expect_coarse_value(.coarse_probe(rank_two, max_rules), rank_two,
                         .coarse_rule_order(rank_two, max_rules),
                         info = paste("rank two, max_rules", max_rules))
  }

  # A declared third column without loadings is dropped before integration: the
  # surrogate is the rank-two one, bit for bit, and equals the two-axis rule.
  compacted <- .coarse_factor_block()
  compacted[["parameters"]][[3L]] <- c(compacted[["parameters"]][[3L]], rep(0, 5))
  compacted[["parameters"]][[15L]] <- as.double(BayesTools::selection_qmc_design(6L, 8L, 2L, 1L))
  compacted[["rank"]] <- 3L
  probe_compacted <- .coarse_probe(compacted, 3L)
  probe_two <- .coarse_probe(rank_two, 3L)
  expect_identical(probe_compacted[["surrogate"]], probe_two[["surrogate"]])
  expect_identical(probe_compacted[["log_normalizer"]], probe_two[["log_normalizer"]])
  reference_two <- .coarse_reference_log_normalizer(rank_two, 7L)
  expect_lte(abs(probe_compacted[["log_normalizer"]] - reference_two), 1e-12)

  # Rank three with a nested forest: one column on every row and two on disjoint
  # row subsets. The nested rule evaluates the same tensor sum.
  forest_loading <- cbind(c(.5, .4, .45, .35, .5, .3), c(.5, .4, 0, 0, 0, 0), c(0, 0, .4, .5, 0, 0))
  forest <- .coarse_block("factor", mean = c(.1, -.2, .3, 0, .2, -.1), residual_sd = rep(.6, 6),
    loading = forest_loading, sei = rep(1, 6), omega = c(1, .3),
    x = c(.2, -.5, .1, .4, 0, .3), obs_bin = c(1L, 2L, 1L, 2L, 1L, 1L),
    rules = .coarse_rules_factor())
  for (max_rules in c(3L, 4L)) {
    .expect_coarse_value(.coarse_probe(forest, max_rules), forest,
                         .coarse_rule_order(forest, max_rules),
                         info = paste("nested forest, max_rules", max_rules))
  }
})


test_that("blocks without an available rule fall back to the exact density", {

  same_as_exact <- function(probe, info) {

    expect_identical(probe[["surrogate_error"]], probe[["exact_error"]], info = info)
    expect_identical(probe[["surrogate"]], probe[["exact"]], info = info)
  }

  for (block in list(.coarse_cluster_block(), .coarse_factor_block())) {
    kind <- block[["kind"]]
    # the NORMAL kernel has no normalizer: the numerator is the density
    normal <- .coarse_probe(.coarse_block_set(block, 11L, 0), 3L)
    expect_false(normal[["exact_fallback"]], info = kind)
    expect_identical(normal[["log_normalizer"]], 0)
    same_as_exact(normal, paste(kind, "NORMAL kernel"))
    gaussian <- .coarse_reference_numerator(block) - sum(log(block[["parameters"]][[5L]][block[["parameters"]][[8L]]]))
    expect_lte(abs(normal[["surrogate"]] - gaussian), 1e-13)

    # a best-weight vector rule has no product rung
    for (vector_rule in c(1, 2)) {
      best <- .coarse_probe(.coarse_block_set(block, 20L, vector_rule), 3L)
      expect_true(best[["exact_fallback"]], info = kind)
      expect_true(is.nan(best[["log_normalizer"]]), info = kind)
      same_as_exact(best, paste(kind, "vector rule", vector_rule))
    }
  }

  # a rung the node budget does not afford: three axes, no forest (every column
  # loads on every row), rule 3 of order 70 (70^3 nodes exceed the budget)
  tensor_loading <- cbind(c(.5, .4, .45, .35), c(.3, -.3, .35, .25), c(.25, .3, -.3, .2))
  tensor <- .coarse_block("factor", mean = c(.1, -.2, .3, 0), residual_sd = rep(.6, 4),
    loading = tensor_loading, sei = rep(1, 4), omega = c(1, .3), x = c(.2, -.5, .1, .4),
    obs_bin = c(1L, 2L, 1L, 1L), rules = .coarse_rules(c(7, 31, 70, 95)),
    initial_points = 64, max_points = 512, scrambles = 8)
  probe <- .coarse_probe(tensor, 3L)
  expect_true(probe[["exact_fallback"]])
  same_as_exact(probe, "unaffordable rung")
  # the same block is on rung 3 once a rung that fits is the third
  affordable <- tensor
  affordable[["rules"]] <- .coarse_rules(c(7, 21, 31, 95))
  affordable[["parameters"]][c(12L, 13L, 14L)] <- list(
    as.double(affordable[["rules"]][["nodes"]]), as.double(affordable[["rules"]][["log_weights"]]),
    as.double(affordable[["rules"]][["orders"]]))
  .expect_coarse_value(.coarse_probe(affordable, 3L), affordable, 31L, info = "affordable rung")
})


test_that("dynamic-invalid inputs have a -Inf surrogate and fail the dynamic check only", {

  for (block in list(.coarse_cluster_block(), .coarse_factor_block())) {
    kind <- block[["kind"]]
    valid <- .coarse_probe(block)
    expect_true(valid[["valid"]] && valid[["dynamic_valid"]] && valid[["control_valid"]], info = kind)
    first <- function(index, value) {
      changed <- block[["parameters"]][[index]]
      changed[1L] <- value
      .coarse_block_set(block, index, changed)
    }
    invalid <- list(
      mean_infinite    = first(1L, Inf),
      mean_missing     = first(1L, NaN),
      sd_zero          = first(2L, 0),
      sd_negative      = first(2L, -1),
      sd_infinite      = first(2L, Inf),
      loading_infinite = first(3L, Inf),
      omega_negative   = .coarse_block_set(block, 5L, c(-1, .3)),
      omega_infinite   = .coarse_block_set(block, 5L, c(Inf, .3)),
      kernel_mode      = .coarse_block_set(block, 11L, 2)
    )
    for (case in names(invalid)) {
      probe <- .coarse_probe(invalid[[case]])
      info  <- paste(kind, case)
      expect_identical(probe[["surrogate"]], -Inf, info = info)
      expect_false(probe[["dynamic_valid"]], info = info)
      expect_true(probe[["control_valid"]], info = info)
      expect_false(probe[["valid"]], info = info)
    }
  }
})


test_that("surrogate values do not depend on the call history or on exact evaluations", {

  # The probe evaluates the exact density at each state too, so the exact node
  # memo is filled with these states between the surrogate calls. A crude rule 3
  # makes the surrogate differ from the exact density: a surrogate that read the
  # memo or kept state would return the exact value or change.
  crude <- list(cluster = .coarse_rules(c(1, 2, 3, 7, 15, 31, 63, 127, 255, 511)),
                factor  = .coarse_rules(c(1, 2, 3, 5, 7, 9, 11, 15, 21, 31, 43, 63, 95)))
  states <- list(
    cluster = list(.coarse_cluster_block(rules = crude$cluster),
                   .coarse_cluster_block(rules = crude$cluster, x = c(.9, -.1, .3, .2))),
    factor  = list(.coarse_factor_block(rules = crude$factor),
                   .coarse_factor_block(rules = crude$factor, x = c(.9, -.1, .3, .2, -.4))))
  for (kind in names(states)) {
    a <- states[[kind]][[1L]]
    b <- .coarse_block_set(states[[kind]][[2L]], 1L, states[[kind]][[2L]][["parameters"]][[1L]] + .1)
    sequence <- list(a, b, a, b, a, a)
    probes <- lapply(sequence, .coarse_probe)
    surrogate <- vapply(probes, `[[`, numeric(1L), "surrogate")
    exact     <- vapply(probes, `[[`, numeric(1L), "exact")
    expect_identical(surrogate[c(1L, 3L, 5L, 6L)], rep(surrogate[1L], 4L), info = kind)
    expect_identical(surrogate[c(2L, 4L)], rep(surrogate[2L], 2L), info = kind)
    expect_identical(exact[c(1L, 3L, 5L, 6L)], rep(exact[1L], 4L), info = kind)
    expect_gt(abs(surrogate[1L] - exact[1L]), 1e-3)
    expect_gt(abs(surrogate[2L] - exact[2L]), 1e-3)
    expect_true(all(is.finite(c(surrogate, exact))), info = kind)
  }
})


# The validators were split into a dynamic part (kernel mode, mean, residual SD,
# loadings, weights) and a control part (everything else) for the sampler's
# fixed-control certificate; checkParameterValue is both. Each validated
# parameter gets an invalid value that JAGS must refuse at compile time with
# the same message as before the split.
test_that("the split validators refuse every invalid dynamic and control value at compile time", {

  skip_if_not_installed("rjags")

  fixtures <- .native_selection_jags_fixtures()
  base     <- fixtures$base
  compile <- function(specification, data) {

    data   <- data[unique(c("y", specification$arguments))]
    syntax <- paste0("model { y[1:2] ~ ", specification$distribution, "(",
      paste(specification$arguments, collapse = ", "), ") }")
    connection <- textConnection(syntax)
    on.exit(close(connection), add = TRUE)
    rjags::jags.model(connection, data = data, n.chains = 1L, n.adapt = 0L, quiet = TRUE)
    invisible(TRUE)
  }
  at_first <- function(value) function(x) { x[1L] <- value; x }
  at_second <- function(value) function(x) { x[2L] <- value; x }
  replace  <- function(value) function(x) value
  # parameter, case, role, change
  table <- list(
    list("mu",       "infinite",            "dynamic", at_first(Inf)),
    list("sd",       "zero",                "dynamic", at_first(0)),
    list("sd",       "negative",            "dynamic", at_second(-1)),
    list("sd",       "infinite",            "dynamic", at_first(Inf)),
    list("loading",  "infinite",            "dynamic", at_second(Inf)),
    list("omega",    "negative",            "dynamic", at_second(-.1)),
    list("omega",    "infinite",            "dynamic", at_first(Inf)),
    list("mode",     "not a kernel mode",   "dynamic", at_first(2)),
    list("mode",     "fractional",          "dynamic", at_first(.5)),
    list("sei",      "zero",                "control", at_first(0)),
    list("sei",      "negative",            "control", at_second(-1)),
    list("sei",      "infinite",            "control", at_first(Inf)),
    list("lower",    "not below upper",     "control", replace(c(1, -1e300))),
    list("upper",    "not above lower",     "control", replace(c(0, 0))),
    list("lower",    "ascending partition", "control", replace(c(-1e300, 0))),
    list("bins",     "zero",                "control", at_first(0)),
    list("bins",     "beyond the weights",  "control", at_second(3)),
    list("bins",     "fractional",          "control", at_first(1.5)),
    list("sign",     "not plus or minus one", "control", at_first(2)),
    list("telescope","not a flag",          "control", at_first(2)),
    list("nodes",    "infinite",            "control", at_first(Inf)),
    list("log_weights", "infinite",         "control", at_first(Inf)),
    list("orders",   "lengths do not sum",  "control", replace(c(1, 3, 4))),
    list("orders",   "not increasing",      "control", replace(c(5, 3, 1))),
    list("orders",   "fractional",          "control", replace(c(1.5, 2.5, 5))),
    list("points",   "below four",          "control", at_first(3)),
    list("points",   "fractional",          "control", at_first(4.5)),
    list("points",   "above the maximum",   "control", at_first(8)),
    list("maximum",  "below four",          "control", at_first(3)),
    list("scrambles","below two",           "control", at_first(1)),
    list("tolerance","zero",                "control", at_first(0)),
    list("tolerance","infinite",            "control", at_first(Inf)),
    list("qmc",      "wrong length",        "control", function(x) x[-1L]),
    list("rule",     "beyond the rules",    "control", at_first(3)),
    list("rule",     "fractional",          "control", at_first(.5))
  )

  for (name in c("cluster", "factor")) {
    specification <- fixtures$cases[[name]]
    data <- utils::modifyList(base, specification$extra)
    expect_true(isTRUE(compile(specification, data)), info = paste(name, "baseline"))
    checked <- character()
    for (row in table) {
      field <- if (row[[1L]] == "loading" && name == "factor") "factor_loading" else row[[1L]]
      if (!field %in% specification$arguments) next
      changed <- data
      changed[[field]] <- row[[4L]](changed[[field]])
      error <- tryCatch(compile(specification, changed), error = identity)
      info  <- paste(name, row[[1L]], row[[2L]], row[[3L]])
      expect_s3_class(error, "error")
      expect_match(conditionMessage(error), "Invalid parent values", fixed = TRUE, info = info)
      checked <- c(checked, paste(row[[3L]], field))
    }
    # every validated dynamic parameter and every control parameter is covered
    expect_true(all(c("dynamic mu", "dynamic sd", "dynamic omega", "dynamic mode",
      "dynamic loading"[name == "cluster"], "dynamic factor_loading"[name == "factor"],
      "control sei", "control bins", "control sign", "control telescope", "control nodes",
      "control log_weights", "control orders", "control points", "control maximum",
      "control scrambles", "control tolerance", "control rule", "control lower",
      "control upper", "control qmc") %in% checked), info = name)
  }
})


# Fits of one cluster, one factor and one best-weight specification (Assink et
# al. subsets) under selection.sampler = "coarse_corrected". Studies 1, 2 and 4
# give two rank-one blocks and a singleton, studies 12, 4 and 5 one rank-two
# block and two singletons; the best-weight rule has no product rung.
.coarse_blocks_tiny_fit <- function(case) {

  data("dat.assink2016", package = "metadat", envir = environment())
  studies <- if (case == "factor") c(12, 4, 5) else c(1, 2, 4)
  dat     <- dat.assink2016[dat.assink2016$study %in% studies, ]
  dat[["study"]] <- as.integer(factor(dat[["study"]]))
  V <- if (case == "factor") {
    metafor::vcalc(vi, cluster = study, type = deltype, obs = esid, rho = c(.7, .5), data = dat)
  } else {
    metafor::vcalc(vi, cluster = study, obs = esid, rho = .6, data = dat)
  }
  selection <- if (case == "best") {
    RoBMA::selection_model(weight_rule = "best", group = "study")
  } else RoBMA::selection_model()
  fit <- suppressWarnings(RoBMA::bselmodel.mv(
    yi = yi, V = V, random = ~ 1 | study / esid, data = dat, measure = "SMD",
    selection = selection, chains = 1L, adapt = 50L, burnin = 50L, sample = 100L,
    seed = 1L, parallel = FALSE, silent = TRUE,
    convergence_checks = RoBMA::set_convergence_checks(max_Rhat = NULL, min_ESS = NULL)))
  list(fit = fit, samplers = fit[["fit"]][["samplers"]],
       methods = RoBMA:::.data_selection_execution_plan(fit[["data"]])[["block_methods"]])
}


test_that("cluster and factor blocks give their nodes the corrected sampler", {

  skip_on_cran()
  skip_if_not_installed("rjags")
  skip_if_not_installed("metadat")
  skip_if_not_installed("metafor")

  coarse   <- "RoBMA::CoarseCorrectedSlice"
  previous <- RoBMA::RoBMA.options()[c("selection.sampler")]
  on.exit(do.call(RoBMA::RoBMA.options, previous), add = TRUE)
  RoBMA::RoBMA.options(selection.sampler = "coarse_corrected")

  for (case in c("cluster", "factor")) {
    RoBMA::selection_sampler_info(clear = TRUE)
    tiny     <- .coarse_blocks_tiny_fit(case)
    samplers <- tiny[["samplers"]]
    counters <- RoBMA::selection_sampler_info()
    # a block of several rows makes its random-effect node and every global
    # node eligible; a node with only singleton children keeps the stock sampler
    expect_true(any(c("rank_one", "factor") %in% tiny[["methods"]]), info = case)
    is_z <- grepl("_Zx[", samplers[["Node"]], fixed = TRUE)
    expect_true(all(samplers[["Sampler"]][!is_z] == coarse), info = case)
    expect_identical(sum(samplers[["Sampler"]][is_z] == coarse),
                     sum(tiny[["methods"]] %in% c("rank_one", "factor")), info = case)
    expect_true(all(samplers[["Sampler"]][is_z & samplers[["Sampler"]] != coarse] == "base::RealSlicer"),
                info = case)
    n_coarse <- sum(samplers[["Sampler"]] == coarse)
    expect_gt(counters[["proposed"]], 0)
    expect_identical(counters[["proposed"]] %% n_coarse, 0)
    expect_identical(counters[["accepted"]] + counters[["correction_rejected"]],
                     counters[["proposed"]], info = case)
  }

  # blocks without a product rung (best weight rule) are not coarse-capable
  RoBMA::selection_sampler_info(clear = TRUE)
  best <- .coarse_blocks_tiny_fit("best")
  expect_true(any(c("rank_one", "factor") %in% best[["methods"]]))
  expect_true(all(best[["samplers"]][["Sampler"]] == "base::RealSlicer"))
  expect_identical(RoBMA::selection_sampler_info()[["proposed"]], 0)
})


# The node memo switch is read once per process, so each setting fits in a fresh
# worker that loads the RoBMA under test (as in test-00-selection-node-memo.R).
.coarse_blocks_memo_worker <- function(memo) {

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
    RoBMA::RoBMA.options(selection.sampler = "coarse_corrected")
    RoBMA::selection_sampler_info(clear = TRUE)
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
      dll      = normalizePath(getLoadedDLLs()[["RoBMA"]][["path"]],
                               winslash = "/", mustWork = TRUE),
      memo     = Sys.getenv("ROBMA_SELNORM_NODE_MEMO", unset = "<unset>"),
      methods  = RoBMA:::.data_selection_execution_plan(fit[["data"]])[["block_methods"]],
      samplers = fit[["fit"]][["samplers"]][["Sampler"]],
      counters = RoBMA::selection_sampler_info(),
      draws    = as.data.frame(RoBMA::as_draws_df(fit))
    )
  }
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


test_that("a corrected chain draws identically with the node memo on and off", {

  skip_on_cran()
  skip_if_not_installed("rjags")
  skip_if_not_installed("metadat")
  skip_if_not_installed("metafor")

  with_memo    <- .coarse_blocks_memo_worker("1")
  without_memo <- .coarse_blocks_memo_worker("0")

  expect_identical(with_memo[["memo"]], "1")
  expect_identical(without_memo[["memo"]], "0")
  expect_true(all(c("factor", "rank_one") %in% with_memo[["methods"]]))
  # the memo holds exact densities only, so the corrected chain, which asks for
  # them at every update, must not see any difference
  expect_true("RoBMA::CoarseCorrectedSlice" %in% with_memo[["samplers"]])
  expect_gt(with_memo[["counters"]][["proposed"]], 0)
  expect_identical(with_memo[["samplers"]], without_memo[["samplers"]])
  expect_identical(with_memo[["counters"]], without_memo[["counters"]])
  expect_identical(dim(with_memo[["draws"]]), c(200L, ncol(with_memo[["draws"]])))
  expect_identical(with_memo[["draws"]], without_memo[["draws"]])
})
