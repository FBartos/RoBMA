context("Cached fit catalog")

source(testthat::test_path("common-functions.R"))

test_that("GLMM test settings preserve full certification draws", {

  old_profile <- Sys.getenv("ROBMA_TEST_PROFILE", unset = NA_character_)
  on.exit({
    if (is.na(old_profile)) {
      Sys.unsetenv("ROBMA_TEST_PROFILE")
    } else {
      Sys.setenv(ROBMA_TEST_PROFILE = old_profile)
    }
  })

  Sys.setenv(ROBMA_TEST_PROFILE = "standard")
  standard <- test_glmm_fit_settings()
  Sys.setenv(ROBMA_TEST_PROFILE = "certification")
  certification <- test_glmm_fit_settings()

  expect_identical(standard, list(
    chains = 3L,
    sample = 1500L,
    burnin = 500L,
    adapt  = 500L
  ))
  expect_identical(certification, list(
    chains = 3L,
    sample = 5000L,
    burnin = 2000L,
    adapt  = 500L
  ))
})

test_that("test cache roots isolate profiles", {

  cache_root    <- tempfile("robma-profile-cache-")
  standard      <- .test_profile_cache_dir(cache_root, "standard")
  certification <- .test_profile_cache_dir(cache_root, "certification")

  expect_identical(basename(standard), "standard")
  expect_identical(basename(certification), "certification")
  expect_identical(dirname(standard), dirname(certification))

  named_root <- file.path(cache_root, "standard")
  expect_identical(
    .test_profile_cache_dir(named_root, "certification"),
    file.path(normalizePath(named_root, winslash = "/", mustWork = FALSE),
              "certification")
  )
})

.collect_save_fit_calls <- function(expr, source_file) {

  out <- data.frame(
    name        = character(),
    source_file = character(),
    stringsAsFactors = FALSE
  )

  if (is.call(expr) && identical(.source_hash_call_name(expr), "save_fit")) {
    fit_name <- expr[[2]]
    if (is.character(fit_name) && length(fit_name) == 1L) {
      out <- rbind(
        out,
        data.frame(
          name        = fit_name,
          source_file = source_file,
          stringsAsFactors = FALSE
        )
      )
    }
  }

  if (is.recursive(expr)) {
    children <- lapply(as.list(expr), .collect_save_fit_calls, source_file = source_file)
    children <- Filter(function(x) nrow(x) > 0L, children)
    if (length(children) > 0L) {
      out <- rbind(out, do.call(rbind, children))
    }
  }

  return(out)
}

.test_01_saved_fits <- function() {

  test_files <- list.files(
    testthat::test_path(),
    pattern    = "^test-01-.*\\.R$",
    full.names = FALSE
  )
  saved_fits <- lapply(test_files, function(source_file) {

    parsed <- parse(testthat::test_path(source_file))
    do.call(rbind, lapply(parsed, .collect_save_fit_calls, source_file = source_file))
  })
  saved_fits <- Filter(function(x) nrow(x) > 0L, saved_fits)

  if (length(saved_fits) == 0L) {
    return(data.frame(
      name        = character(),
      source_file = character(),
      stringsAsFactors = FALSE
    ))
  }

  return(do.call(rbind, saved_fits))
}

.fake_fit_for_catalog_entry <- function(name) {

  entry <- fit_catalog_entry(name)
  fit   <- list()

  if (isTRUE(entry[["has_loo"]])) {
    fit[["loo"]] <- list(estimate = "loo")
  }
  if (isTRUE(entry[["has_waic"]])) {
    fit[["waic"]] <- list(estimate = "waic")
  }
  if (isTRUE(entry[["has_marglik"]])) {
    fit[["marglik"]] <- list(estimate = "marglik")
  }

  class(fit) <- entry[["class"]]

  return(fit)
}

.fake_info_for_catalog_entry <- function(name) {

  entry <- fit_catalog_entry(name)

  if (isTRUE(entry[["has_metafor"]])) {
    return(list(metafor = list(reference = TRUE)))
  }

  return(list())
}

.valid_metadata_for_catalog_entry <- function(name) {

  entry    <- fit_catalog_entry(name)
  versions <- fit_cache_versions()

  return(list(
    schema_version                 = FIT_CACHE_SCHEMA_VERSION,
    name                           = name,
    saved_at                       = "synthetic",
    robma_version                  = versions[["robma_version"]],
    bayestools_version             = versions[["bayestools_version"]],
    r_version                      = versions[["r_version"]],
    jags_version                   = versions[["jags_version"]],
    fit_class                      = entry[["class"]],
    source_file                    = entry[["source_file"]],
    source_file_md5                = source_file_md5(entry[["source_file"]]),
    package_source_md5             = package_source_md5(),
    bayestools_backend_fingerprint = bayestools_backend_fingerprint(),
    has_loo                        = entry[["has_loo"]],
    has_waic                       = entry[["has_waic"]],
    has_marglik                    = entry[["has_marglik"]],
    has_metafor_info               = entry[["has_metafor"]]
  ))
}

.with_temp_fit_cache <- function(code) {

  old_test_files_dir    <- test_files_dir
  old_temp_fits_dir     <- temp_fits_dir
  old_temp_info_dir     <- temp_info_dir
  old_temp_metadata_dir <- temp_metadata_dir
  old_temp_temp_dir     <- temp_temp_dir
  old_env               <- Sys.getenv("ROBMA_TEST_FILES_DIR", unset = NA_character_)
  old_fit_names         <- ls(envir = .fit_object_cache)
  old_info_names        <- ls(envir = .info_object_cache)
  old_fit_cache         <- mget(old_fit_names, envir = .fit_object_cache, inherits = FALSE)
  old_info_cache        <- mget(old_info_names, envir = .info_object_cache, inherits = FALSE)

  cache_root <- tempfile("robma-fit-cache-")
  dir.create(cache_root, recursive = TRUE)

  test_files_dir <<- normalizePath(cache_root, winslash = "/", mustWork = TRUE)
  temp_fits_dir     <<- file.path(test_files_dir, "fits")
  temp_info_dir     <<- file.path(test_files_dir, "info")
  temp_metadata_dir <<- file.path(test_files_dir, "metadata")
  temp_temp_dir     <<- file.path(test_files_dir, "temp")
  dir.create(temp_fits_dir, recursive = TRUE)
  dir.create(temp_info_dir, recursive = TRUE)
  dir.create(temp_metadata_dir, recursive = TRUE)
  dir.create(temp_temp_dir, recursive = TRUE)

  .clear_fit_object_cache()
  .clear_info_object_cache()
  Sys.setenv(ROBMA_TEST_FILES_DIR = test_files_dir)

  on.exit({
    test_files_dir <<- old_test_files_dir
    temp_fits_dir <<- old_temp_fits_dir
    temp_info_dir <<- old_temp_info_dir
    temp_metadata_dir <<- old_temp_metadata_dir
    temp_temp_dir <<- old_temp_temp_dir
    .clear_fit_object_cache()
    .clear_info_object_cache()
    if (length(old_fit_cache) > 0L) {
      list2env(old_fit_cache, envir = .fit_object_cache)
    }
    if (length(old_info_cache) > 0L) {
      list2env(old_info_cache, envir = .info_object_cache)
    }
    if (is.na(old_env)) {
      Sys.unsetenv("ROBMA_TEST_FILES_DIR")
    } else {
      Sys.setenv(ROBMA_TEST_FILES_DIR = old_env)
    }
    unlink(cache_root, recursive = TRUE)
  }, add = TRUE)

  force(code)
}

# Structural carriers exercise cache provenance, never posterior calculations.
.with_cache_source_fixture <- function(code) {

  real_root <- .fit_cache_source_root()
  required  <- .fit_cache_source_files(package_root = real_root, relative = TRUE)
  skip_if_not(length(required) > 0L, "Complete package source tree is unavailable.")
  root <- tempfile("robma-producer-source-")
  paths <- unique(c(required, "DESCRIPTION", "tests/testthat/common-functions.R",
                    "tests/testthat/test-02-cache-catalog.R",
                    paste0("tests/testthat/", unique(fit_catalog()[["source_file"]]))))
  for (path in paths) {
    target <- file.path(root, path)
    dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
    stopifnot(file.copy(file.path(real_root, path), target))
  }
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  cache_names <- ls(envir = .package_source_md5_cache)
  old_cache   <- mget(cache_names, envir = .package_source_md5_cache)
  old_env     <- Sys.getenv("ROBMA_TEST_FILES_DIR", unset = NA_character_)
  on.exit({
    .clear_package_source_md5_cache()
    if (length(old_cache) > 0L) list2env(old_cache, .package_source_md5_cache)
    if (is.na(old_env)) Sys.unsetenv("ROBMA_TEST_FILES_DIR") else
      Sys.setenv(ROBMA_TEST_FILES_DIR = old_env)
  }, add = TRUE)
  .clear_package_source_md5_cache()
  Sys.setenv(ROBMA_TEST_FILES_DIR = file.path(root, "cache"))
  testthat::local_test_directory(file.path(root, "tests/testthat"),
                                 package = "RoBMA", .env = environment())
  e <- new.env(parent = environment(.with_cache_source_fixture))
  source(file.path(root, "tests/testthat/common-functions.R"), local = e)
  for (name in c(".fake_fit_for_catalog_entry", ".fake_info_for_catalog_entry")) {
    fun <- get(name)
    environment(fun) <- e
    assign(name, fun, e)
  }
  # The unchanged helper has no generation operation; initializing its actual
  # metadata here retains a genuine A-init/B-edit baseline for the regression.
  e$fit_cache_metadata("brma.mv_block_mvn_random",
                        e$.fake_fit_for_catalog_entry("brma.mv_block_mvn_random"))
  code(e, root)
}

test_that("producer generation refuses changed source before saving", {

  skip_on_cran()
  for (target in c("tests/testthat/test-01-brma.mv.R", "R/BMA.glmm.R")) {
    .with_cache_source_fixture(function(e, root) {
      name <- "brma.mv_block_mvn_random"
      fit  <- e$.fake_fit_for_catalog_entry(name)
      expect_true(e$save_fit(name, fit))
      paths <- unlist(e$fit_cache_paths(name)[c("fit", "metadata")])
      saved <- tools::md5sum(paths)
      old_metadata <- e$load_fit_metadata(name)
      path <- file.path(root, target)
      cat("\n.cache_generation_marker <- 'B'\n", file = path, append = TRUE)
      expect_error(e$save_fit(name, fit), "Reload.*rerun")
      expect_identical(tools::md5sum(paths), saved)
      expect_true(length(e$validate_cached_fit(name)) > 0L)
      source(file.path(root, "tests/testthat/common-functions.R"), local = e)
      expect_error(e$save_fit(name, fit), "Reload.*rerun")
      expect_identical(tools::md5sum(paths), saved)
      # Explicit new producer lifecycle is the only permitted recapture.
      if (exists(".fit_cache_generation_begin", e, inherits = FALSE)) {
        e$.fit_cache_generation_begin(force = TRUE)
      }
      expect_true(e$save_fit(name, fit))
      expect_equal(e$validate_cached_fit(name, deep = TRUE), character())
      changed_field <- if (startsWith(target, "R/")) "package_source_md5" else "source_file_md5"
      expect_false(identical(e$load_fit_metadata(name)[[changed_field]],
                             old_metadata[[changed_field]]))
    })
  }
})

test_that("producer generation captures all builders before later helper sources", {

  skip_on_cran()
  .with_cache_source_fixture(function(e, root) {
    first <- "brma.mv_block_mvn_random"
    later <- "bcg_meta-analysis"
    expect_true(e$save_fit(first, e$.fake_fit_for_catalog_entry(first)))
    cat("\n.cache_later_builder_marker <- 'B'\n",
        file = file.path(root, "tests/testthat/test-01-brma.norm.R"), append = TRUE)
    source(file.path(root, "tests/testthat/common-functions.R"), local = e)
    expect_error(e$save_fit(later, e$.fake_fit_for_catalog_entry(later),
                            e$.fake_info_for_catalog_entry(later)), "Reload.*rerun")
    expect_false(file.exists(e$fit_cache_paths(later)[["fit"]]))
    expect_equal(e$validate_cached_fit(first), character())
  })
})

test_that("helper initialization leaves product loading to the profile runner", {

  skip_on_cran()
  skip_if_not_installed("processx")
  .with_cache_source_fixture(function(e, root) {
    script <- file.path(root, "helper-startup.R")
    output <- file.path(root, "helper-startup.rds")
    writeLines(c(
      "main <- function() {",
      "  args <- commandArgs(trailingOnly = TRUE)",
      "  before <- 'RoBMA' %in% loadedNamespaces()",
      "  testthat::local_test_directory(file.path(args[[1L]], 'tests/testthat'), .env = environment())",
      "  source(file.path(args[[1L]], 'tests/testthat/common-functions.R'), local = new.env(parent = globalenv()))",
      "  saveRDS(list(before = before, after = 'RoBMA' %in% loadedNamespaces()), args[[2L]])",
      "}", "main()"
    ), script)
    result <- processx::run(
      file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"),
      c("--vanilla", script, root, output), error_on_status = FALSE
    )
    expect_equal(result[["status"]], 0L, info = result[["stderr"]])
    observed <- readRDS(output)
    expect_false(observed[["before"]])
    expect_false(observed[["after"]])
  })
})

.cache_dependency_parents <- function(e) {

  info_name <- "brma.mv_block_mvn_random"
  fit_name  <- "brma.mv_block_mvn_fixed_random_null"
  info <- list(data = data.frame(yi = 1), V = matrix(1), marker = "A")
  fit <- e$.fake_fit_for_catalog_entry(fit_name)
  fit[["marker"]] <- "A"
  e$save_fit(info_name, e$.fake_fit_for_catalog_entry(info_name), info)
  e$save_fit(fit_name, fit, e$.fake_info_for_catalog_entry(fit_name))
  return(c(info = info_name, fit = fit_name))
}

.consume_cache_dependency_parents <- function(e, parents) {
  list(info = e$load_info(parents[["info"]]), fit = e$load_fit(parents[["fit"]]))
}

.save_cache_dependency_children <- function(e) {
  children <- c("iwmde_known_v_tau_full", "iwmde_known_v_tau_null")
  for (name in children) e$save_fit(name, e$.fake_fit_for_catalog_entry(name))
  return(children)
}

test_that("derived cache entries declare and validate consumed artifacts", {

  skip_if_not_certification("These copied-source cache fixtures belong to the tooling case.")
  skip_on_cran()
  .with_cache_source_fixture(function(e, root) {
    parents <- .cache_dependency_parents(e)
    .consume_cache_dependency_parents(e, parents)
    children <- .save_cache_dependency_children(e)
    expected <- data.frame(name = unname(parents), artifact = names(parents))
    for (child in children) {
      expect_equal(e$fit_catalog_entry(child)[["dependencies"]][[1L]], expected)
      expect_equal(length(e$load_fit_metadata(child)[["dependencies"]]), 2L)
      expect_equal(e$validate_cached_fit(child, deep = TRUE), character())
    }
    # Rewriting identical bytes or a timestamp alone is the same consumed input.
    for (artifact in names(parents)) {
      path <- e$fit_cache_paths(parents[[artifact]])[[artifact]]
      writeBin(readBin(path, "raw", file.info(path)$size), path)
      metadata <- e$load_fit_metadata(parents[[artifact]])
      metadata[["saved_at"]] <- "different timestamp"
      saveRDS(metadata, e$fit_cache_paths(parents[[artifact]])[["metadata"]])
    }
    for (child in children) expect_equal(e$validate_cached_fit(child), character())
    child_paths <- unlist(lapply(children, function(name) e$fit_cache_paths(name)[c("fit", "metadata")]))
    child_bytes <- tools::md5sum(child_paths)
    stale_children <- function() {
      for (child in children) {
        expect_true(length(e$validate_cached_fit(child)) > 0L)
        expect_true(length(e$validate_cached_fit(child, deep = TRUE)) > 0L)
      }
      expect_identical(tools::md5sum(child_paths), child_bytes)
    }
    for (artifact in names(parents)) {
      path <- e$fit_cache_paths(parents[[artifact]])[[artifact]]
      bytes <- readBin(path, "raw", file.info(path)$size)
      value <- readRDS(path)
      value[["marker"]] <- "B"
      if (artifact == "info") value[["data"]][["yi"]] <- 2
      saveRDS(value, path)
      stale_children()
      writeBin(bytes, path)
      file.remove(path)
      stale_children()
      writeBin(bytes, path)
    }
    # The info dependency also requires its parent fit's normal provenance.
    path <- e$fit_cache_paths(parents[["info"]])[["fit"]]
    bytes <- readBin(path, "raw", file.info(path)$size)
    file.remove(path); stale_children(); writeBin(bytes, path)
    path <- file.path(root, "tests/testthat/test-01-brma.mv.R")
    bytes <- readBin(path, "raw", file.info(path)$size)
    cat("\n.cache_parent_builder_marker <- 'B'\n", file = path, append = TRUE)
    stale_children(); writeBin(bytes, path)
    for (child in children) {
      metadata <- e$load_fit_metadata(child)
      metadata[["dependencies"]] <- NULL
      saveRDS(metadata, e$fit_cache_paths(child)[["metadata"]])
      expect_true(any(grepl("dependenc", e$validate_cached_fit(child))))
    }
  })
})

test_that("derived cache saves require current generation consumption receipts", {

  skip_on_cran()
  .with_cache_source_fixture(function(e, root) {
    parents <- .cache_dependency_parents(e)
    children <- c("iwmde_known_v_tau_full", "iwmde_known_v_tau_null")
    for (child in children) {
      expect_error(e$save_fit(child, e$.fake_fit_for_catalog_entry(child)), "consum")
      expect_false(file.exists(e$fit_cache_paths(child)[["fit"]]))
    }
    .consume_cache_dependency_parents(e, parents)
    .save_cache_dependency_children(e)
    paths <- unlist(lapply(children, function(name) e$fit_cache_paths(name)[c("fit", "metadata")]))
    before <- tools::md5sum(paths)
    if (exists(".fit_cache_generation_begin", e, inherits = FALSE)) e$.fit_cache_generation_begin(force = TRUE)
    e$load_info(parents[["info"]], validate = FALSE)
    e$load_fit(parents[["fit"]], validate = FALSE)
    for (child in children) expect_error(e$save_fit(child, e$.fake_fit_for_catalog_entry(child)), "consum")
    expect_identical(tools::md5sum(paths), before)
    .consume_cache_dependency_parents(e, parents)
    source(file.path(root, "tests/testthat/common-functions.R"), local = e)
    expect_true(e$save_fit(children[[1L]], e$.fake_fit_for_catalog_entry(children[[1L]])))
    before <- tools::md5sum(paths)
    path <- e$fit_cache_paths(parents[["info"]])[["info"]]
    value <- readRDS(path); value[["data"]][["yi"]] <- 2
    saveRDS(value, path)
    parent_paths <- unlist(lapply(parents, function(name) e$fit_cache_paths(name)[c("fit", "info", "metadata")]))
    parent_paths <- parent_paths[file.exists(parent_paths)]
    parent_before <- tools::md5sum(parent_paths)
    for (child in children) expect_error(e$save_fit(child, e$.fake_fit_for_catalog_entry(child)), "changed")
    expect_identical(tools::md5sum(paths), before)
    expect_identical(tools::md5sum(parent_paths), parent_before)
  })
})

test_that("warmed dependency objects retain first consumption conflicts", {

  skip_if_not_certification("These copied-source cache fixtures belong to the tooling case.")
  skip_on_cran()
  for (artifact in c("fit", "info")) {
    .with_cache_source_fixture(function(e, root) {
      parents <- .cache_dependency_parents(e)
      consumed <- .consume_cache_dependency_parents(e, parents)
      expect_identical(consumed[[artifact]][["marker"]], "A")
      children <- .save_cache_dependency_children(e)
      paths <- unlist(lapply(children, function(name) e$fit_cache_paths(name)[c("fit", "metadata")]))
      before <- tools::md5sum(paths)
      path <- e$fit_cache_paths(parents[[artifact]])[[artifact]]
      bytes <- readBin(path, "raw", file.info(path)$size)
      value <- readRDS(path); value[["marker"]] <- "B"; saveRDS(value, path)
      reloaded <- if (artifact == "fit") e$load_fit(parents[[artifact]]) else e$load_info(parents[[artifact]])
      expect_identical(reloaded[["marker"]], "B")
      for (child in children) expect_error(e$save_fit(child, e$.fake_fit_for_catalog_entry(child)), "conflicting")
      # A later return to the first bytes cannot erase conflicting consumption.
      writeBin(bytes, path)
      .consume_cache_dependency_parents(e, parents)
      for (child in children) expect_error(e$save_fit(child, e$.fake_fit_for_catalog_entry(child)), "conflicting")
      expect_identical(tools::md5sum(paths), before)
    })
  }
})

test_that("declared cache dependency cycles fail validation", {

  skip_on_cran()
  .with_cache_source_fixture(function(e, root) {
    parents <- .cache_dependency_parents(e)
    .consume_cache_dependency_parents(e, parents)
    child <- .save_cache_dependency_children(e)[[1L]]
    catalog <- e$fit_catalog()
    if (is.null(catalog[["dependencies"]])) catalog[["dependencies"]] <- rep(list(NULL), nrow(catalog))
    catalog[["dependencies"]][[match(child, catalog[["name"]])]] <-
      data.frame(name = parents[["fit"]], artifact = "fit")
    catalog[["dependencies"]][[match(parents[["fit"]], catalog[["name"]])]] <-
      data.frame(name = child, artifact = "fit")
    expect_true(any(grepl("cycle", e$validate_cached_fit(child, catalog = catalog))))
  })
})

test_that("dependency reads detect artifact mutation during deserialization", {

  skip_on_cran()
  .with_cache_source_fixture(function(e, root) {
    parents <- .cache_dependency_parents(e)
    e$.clear_info_object_cache()
    read_info <- e$read_cached_info
    e$read_cached_info <- function(name, ...) {
      value <- read_info(name, ...)
      changed <- value; changed[["marker"]] <- "B"
      saveRDS(changed, e$fit_cache_paths(name)[["info"]])
      value
    }
    expect_error(e$load_info(parents[["info"]]), "changed.*read")
  })
})

test_that("present corrupt info cannot authorize a derived cache save", {

  skip_on_cran()
  .with_cache_source_fixture(function(e, root) {
    parents <- .cache_dependency_parents(e)
    path <- e$fit_cache_paths(parents[["info"]])[["info"]]
    bytes <- readBin(path, "raw", file.info(path)$size)
    expect_equal(e$load_info(parents[["info"]])[["data"]][["yi"]], 1)
    file.remove(path)
    expect_identical(e$load_info(parents[["info"]]), list())
    writeBin(bytes, path)
    .consume_cache_dependency_parents(e, parents)
    children <- .save_cache_dependency_children(e)
    e$.fit_cache_generation_begin(force = TRUE)
    e$load_fit(parents[["fit"]])
    writeBin(charToRaw("present but not an RDS artifact"), path)
    e$.clear_info_object_cache()
    paths <- unlist(lapply(c(parents, children), function(name) e$fit_cache_paths(name)[c("fit", "info", "metadata")]))
    paths <- paths[file.exists(paths)]
    before <- tools::md5sum(paths)
    expect_error(e$load_info(parents[["info"]]), "could not be read")
    for (child in children) expect_error(e$save_fit(child, e$.fake_fit_for_catalog_entry(child)), "consum")
    expect_identical(tools::md5sum(paths), before)
  })
})

test_that("cache inventory assertion rejects inactive and uncatalogued artifacts", {

  skip_if(is_certification_profile(),
          "The fixture switches to the standard profile, which changes the producer provenance of a certification process.")
  skip_on_cran()
  .with_cache_source_fixture(function(e, root) {
    withr::local_envvar(c(ROBMA_TEST_PROFILE = "standard", ROBMA_TEST_ACTIVE_FITS = ""))
    active <- "bcg_meta-analysis"
    inactive <- "bcg_glmm_3lvl_scale"
    stray <- "uncatalogued_cache_inventory_probe"
    e$save_fit(active, e$.fake_fit_for_catalog_entry(active), e$.fake_info_for_catalog_entry(active))
    expressions <- as.list(parse(file.path(root, "tests/testthat/test-02-cache-catalog.R")))
    assertion <- Filter(function(expr) is.call(expr) &&
      identical(expr[[1L]], as.name("test_that")) &&
      identical(expr[[2L]], "cached fits are catalogued and valid"), expressions)[[1L]]
    run_assertion <- function() {
      reporter <- testthat::ListReporter$new()
      testthat::with_reporter(reporter, {
        reporter$start_file("test-02-cache-catalog.R")
        testthat::context("Cached fit catalog")
        eval(assertion, e)
        reporter$end_file()
      })
      as.data.frame(reporter$get_results())
    }
    control <- run_assertion()
    expect_equal(sum(control[["failed"]]), 0L)
    expect_false(any(control[["error"]] | control[["skipped"]]))
    e$save_fit(inactive, e$.fake_fit_for_catalog_entry(inactive), e$.fake_info_for_catalog_entry(inactive))
    expect_true(sum(run_assertion()[["failed"]]) > 0L)
    file.remove(e$fit_cache_paths(inactive)[["fit"]])
    saveRDS(list(structural_probe = TRUE), e$fit_cache_paths(stray)[["fit"]])
    expect_true(sum(run_assertion()[["failed"]]) > 0L)
    expect_equal(e$list_fits(validate = FALSE), active)
  })
})

test_that("cache profile membership includes actual positive and negative members", {

  withr::local_envvar(c(ROBMA_TEST_PROFILE = "standard", ROBMA_TEST_ACTIVE_FITS = ""))
  expect_true("bcg_meta-analysis" %in% active_fit_catalog()[["name"]])
  expect_false("bcg_glmm_3lvl_scale" %in% active_fit_catalog()[["name"]])
  Sys.setenv(ROBMA_TEST_PROFILE = "certification")
  expect_true("bcg_glmm_3lvl_scale" %in% active_fit_catalog()[["name"]])
  expect_true("bcg_meta-analysis" %in% active_fit_catalog()[["name"]])
})


test_that("cache source hash tracks only cache-affecting fitting sources", {

  source_files <- .fit_cache_source_files(relative = TRUE)
  skip_if_not(
    length(source_files) > 0L,
    "Complete package R/src source tree is not available under source-build checks."
  )

  expect_equal(anyDuplicated(source_files), 0L)
  expect_true(all(c(
    # fit-time and stored-result code that is easy to overlook: the selection
    # runtime cache read by fit.R, the block covariance and known-V helpers of
    # the stored LOO and marginal likelihood, the allocation-inclusion priors of
    # the preflight, and the native sampler and tables of the JAGS module
    "R/block-covariance.R",
    "R/brma-mv-targets.R",
    "R/distributions.R",
    "R/known-v-gls.R",
    "R/outcome-helpers.R",
    "R/pdf.R",
    "R/predict-brma-mv.R",
    "R/random-inclusion.R",
    "R/random-marginal-vcov.R",
    "R/rng.R",
    "R/selection-cache.R",
    "R/selection-normalizer-grid.R",
    "src/glmm-binomial-loglik.h",
    "src/samplers/CoarseCorrectedSlice.cc",
    "src/samplers/CoarseCorrectedSlice.h",
    "src/selnorm/selnorm-cdf-table-data.h",
    "src/selnorm/selnorm-mixture-compression.h",
    "src/selnorm/selnorm-mv.cc",
    "R/brma-mv-known-r.R",
    "R/covariance-factorization.R",
    "R/fit.R",
    "R/formula-design.R",
    "R/glmm-aghq.R",
    "R/input-data.R",
    "R/input-object.R",
    "R/input-priors.R",
    "R/input-priors-assignment.R",
    "R/input-priors-check-list.R",
    "R/input-priors-documentation.R",
    "R/input-priors-formula.R",
    "R/input-priors-heterogeneity-allocation.R",
    "R/jags-formula-args.R",
    "R/known-v-preflight.R",
    "R/known-v-representation.R",
    "R/random-effects-compile.R",
    "R/selection-mapping.R",
    "R/marglik.R",
    "R/loo.R",
    "src/Makevars.in",
    "src/Makevars.ucrt",
    "src/Makevars.win",
    "src/RoBMA.cc",
    "src/distributions/DWN.cc",
    "src/glmm-aghq.cc",
    "src/glmm-aghq.h",
    "src/r-native-api.h",
    "src/r-native-boundary.h",
    "src/r-glmm.cc",
    "src/r-selnorm.cc",
    "src/r-selnorm-common.cc.inc",
    "src/r-selnorm-funnel-zcurve.cc.inc",
    "src/r-selnorm-kernel.cc.inc",
    "src/r-selnorm-loglik.cc.inc",
    "src/selnorm/selnorm-api.cc.inc",
    "src/selnorm/selnorm-boundary.cc.inc",
    "src/selnorm/selnorm-parallel.h",
    "src/selnorm/selnorm-phack.cc.inc",
    "src/selnorm/selnorm-probability.cc.inc",
    "src/selnorm/selnorm-step.cc.inc"
  ) %in% source_files))
  expect_false(any(c(
    "DESCRIPTION",
    "NAMESPACE",
    "R/as_draws.R",
    "R/effect-transformations.R",
    "R/plot.R",
    "R/predict.R",
    "R/summary.R",
    "src/r-regplot.cc",
    # post-fit methods and plot or documentation carriers
    "R/forest.R",
    "R/hypothesis.R",
    "R/hypothesis-plan.R",
    "R/iwmde-estimate.R",
    "R/marginal_means.R",
    "R/random-effect-formula-tags.R",
    "R/random-priors.R",
    "R/regplot-data.R",
    "R/zplot-api.R",
    "R/zplot-compute.R",
    "src/plot-root.h",
    "src/r-native-probe.cc.inc"
  ) %in% source_files))
})

test_that("in-memory fit caches are synchronized with save and clean", {

  skip_on_cran()

  .with_temp_fit_cache({
    name <- "bcg_meta-analysis"

    old_fit <- .fake_fit_for_catalog_entry(name)
    old_fit[["cache_marker"]] <- "old"
    old_info <- .fake_info_for_catalog_entry(name)
    old_info[["cache_marker"]] <- "old"

    new_fit <- .fake_fit_for_catalog_entry(name)
    new_fit[["cache_marker"]] <- "new"
    new_info <- .fake_info_for_catalog_entry(name)
    new_info[["cache_marker"]] <- "new"

    save_fit(name, old_fit, old_info)
    expect_identical(load_fit(name)[["cache_marker"]], "old")
    expect_identical(load_info(name)[["cache_marker"]], "old")

    save_fit(name, new_fit, new_info)
    expect_identical(load_fit(name)[["cache_marker"]], "new")
    expect_identical(load_info(name)[["cache_marker"]], "new")

    clean_cached_fits(name)
    expect_false(exists(name, envir = .fit_object_cache, inherits = FALSE))
    expect_false(exists(name, envir = .info_object_cache, inherits = FALSE))

    fit_names <- c("bcg_meta-analysis", "bcg_meta-regression")
    for (fit_name in fit_names) {
      save_fit(
        fit_name,
        .fake_fit_for_catalog_entry(fit_name),
        .fake_info_for_catalog_entry(fit_name)
      )
    }
    clean_cached_fits(fit_names)
    expect_false(any(file.exists(vapply(
      fit_names,
      function(fit_name) fit_cache_paths(fit_name)[["fit"]],
      character(1)
    ))))
  })
})

test_that("cache source hash requires a complete source root", {

  package_root <- tempfile("robma-incomplete-source-")
  dir.create(file.path(package_root, "R"), recursive = TRUE)
  writeLines("Package: RoBMA", file.path(package_root, "DESCRIPTION"), useBytes = TRUE)
  writeLines("fit <- function() NULL", file.path(package_root, "R", "fit.R"), useBytes = TRUE)
  on.exit(unlink(package_root, recursive = TRUE), add = TRUE)

  expect_equal(
    .fit_cache_source_files(package_root = package_root, relative = TRUE),
    character()
  )
})

test_that("package source hash refreshes after source edits", {

  real_root    <- .fit_cache_source_root()
  source_files <- .fit_cache_source_files(package_root = real_root, relative = TRUE)
  skip_if_not(
    length(source_files) > 0L,
    "Complete package R/src source tree is not available under source-build checks."
  )

  package_root <- tempfile("robma-complete-source-")
  dir.create(package_root, recursive = TRUE)
  on.exit(unlink(package_root, recursive = TRUE), add = TRUE)

  for (source_file in source_files) {
    from <- file.path(real_root, source_file)
    to   <- file.path(package_root, source_file)
    dir.create(dirname(to), recursive = TRUE, showWarnings = FALSE)
    file.copy(from, to, overwrite = TRUE)
  }

  old_cache_names <- ls(envir = .package_source_md5_cache)
  old_cache       <- mget(old_cache_names, envir = .package_source_md5_cache,
                          inherits = FALSE)
  on.exit({
    .clear_package_source_md5_cache()
    if (length(old_cache) > 0L) {
      list2env(old_cache, envir = .package_source_md5_cache)
    }
  }, add = TRUE)

  .clear_package_source_md5_cache()
  assign(
    "source_root",
    normalizePath(package_root, winslash = "/", mustWork = TRUE),
    envir = .package_source_md5_cache
  )

  hash_before <- package_source_md5()
  fit_file    <- file.path(package_root, "R", "fit.R")

  writeLines(
    c(readLines(fit_file, warn = FALSE), "# cache hash refresh probe"),
    fit_file,
    useBytes = TRUE
  )
  expect_identical(package_source_md5(), hash_before)

  writeLines(
    c(readLines(fit_file, warn = FALSE), ".robma_cache_hash_probe <- 1"),
    fit_file,
    useBytes = TRUE
  )
  expect_false(identical(package_source_md5(), hash_before))

  for (header in c("src/r-native-api.h", "src/r-native-boundary.h",
                   "src/selnorm/selnorm-parallel.h",
                   "src/glmm-binomial-loglik.h",
                   "src/samplers/CoarseCorrectedSlice.h",
                   "src/selnorm/selnorm-cdf-table-data.h",
                   "src/selnorm/selnorm-mixture-compression.h")) {
    before_header_edit <- package_source_md5()
    header_path <- file.path(package_root, header)
    writeLines(
      c(readLines(header_path, warn = FALSE), "#define ROBMA_NATIVE_CACHE_TEST 1"),
      header_path,
      useBytes = TRUE
    )
    expect_false(identical(package_source_md5(), before_header_edit), info = header)
  }

  # R sources read while a model is constructed or its stored results are
  # computed change the hash, a native source of the JAGS module as well
  for (source_file in c("R/selection-cache.R", "R/block-covariance.R",
                        "R/random-inclusion.R", "R/known-v-gls.R",
                        "R/outcome-helpers.R", "R/rng.R",
                        "src/samplers/CoarseCorrectedSlice.cc")) {
    before_edit <- package_source_md5()
    source_path <- file.path(package_root, source_file)
    probe <- if (startsWith(source_file, "R/")) {
      ".robma_cache_hash_probe_added <- 1"
    } else {
      "int robma_cache_hash_probe_added = 1;"
    }
    writeLines(
      c(readLines(source_path, warn = FALSE), probe), source_path, useBytes = TRUE
    )
    expect_false(identical(package_source_md5(), before_edit), info = source_file)
  }

  # A post-fit-only source is outside the hash, wherever it is edited
  for (post_fit in c("R/forest.R", "R/zplot-api.R", "R/zplot-compute.R",
                   "src/r-regplot.cc")) {
    before_edit <- package_source_md5()
    post_fit_path <- file.path(package_root, post_fit)
    dir.create(dirname(post_fit_path), recursive = TRUE, showWarnings = FALSE)
    file.copy(file.path(real_root, post_fit), post_fit_path, overwrite = TRUE)
    probe <- if (startsWith(post_fit, "R/")) {
      ".robma_post_fit_probe <- 1"
    } else {
      "int robma_post_fit_probe = 1;"
    }
    writeLines(
      c(readLines(post_fit_path, warn = FALSE), probe), post_fit_path, useBytes = TRUE
    )
    expect_identical(package_source_md5(), before_edit, info = post_fit)
  }
})

test_that("GLMM fit settings enter the hash of the files that read them", {

  old_profile <- Sys.getenv("ROBMA_TEST_PROFILE", unset = NA_character_)
  on.exit({
    if (is.na(old_profile)) {
      Sys.unsetenv("ROBMA_TEST_PROFILE")
    } else {
      Sys.setenv(ROBMA_TEST_PROFILE = old_profile)
    }
  })
  hash_files <- function(profile) {

    Sys.setenv(ROBMA_TEST_PROFILE = profile)
    c(
      glmm  = source_file_md5("test-01-brma.glmm.R"),
      other = source_file_md5("test-01-brma.norm.R")
    )
  }

  standard      <- hash_files("standard")
  certification <- hash_files("certification")

  # test-01-brma.glmm.R only calls test_glmm_fit_settings(); the draws it
  # returns differ between the profiles
  expect_false(identical(standard[["glmm"]], certification[["glmm"]]))
  expect_identical(standard[["other"]], certification[["other"]])

  # the value is hashed for a file that references the helper, not for others
  Sys.setenv(ROBMA_TEST_PROFILE = "standard")
  expect_identical(.source_hash_fit_settings(FALSE), character())
  expect_identical(
    .source_hash_fit_settings(TRUE),
    "test_glmm_fit_settings() = list(chains = 3L, sample = 1500L, burnin = 500L, adapt = 500L)"
  )
})

test_that("R source hashes ignore comments and formatting", {

  file_a <- tempfile("robma-source-a-", fileext = ".R")
  file_b <- tempfile("robma-source-b-", fileext = ".R")
  file_c <- tempfile("robma-source-c-", fileext = ".R")
  on.exit(unlink(c(file_a, file_b, file_c)), add = TRUE)

  writeLines(c(
    "# comment before code",
    "f <- function(x) {",
    "  x + 1",
    "}"
  ), file_a, useBytes = TRUE)
  writeLines(c(
    "f<-function(x){x+1} # trailing comment"
  ), file_b, useBytes = TRUE)
  writeLines(c(
    "f <- function(x) {",
    "  x + 2",
    "}"
  ), file_c, useBytes = TRUE)

  expect_identical(
    .fit_cache_source_file_md5(file_a),
    .fit_cache_source_file_md5(file_b)
  )
  expect_false(identical(
    .fit_cache_source_file_md5(file_a),
    .fit_cache_source_file_md5(file_c)
  ))
})


test_that("fit catalog is internally consistent", {

  catalog <- fit_catalog()

  expect_equal(anyDuplicated(catalog[["name"]]), 0L,
               info = "fit catalog names must be unique")
  expect_true(all(file.exists(testthat::test_path(catalog[["source_file"]]))),
              info = "all catalog source files must exist")
  no_loo_fits <- c(
    "brma.mv_block_mvn_random_mods_scale",
    "nielweise2008_glmm_effect_null",
    "dat.lehmann2018-3PSM_effect_null",
    "hypothesis_plan_treatment",
    "hypothesis_plan_meandif",
    "hypothesis_plan_orthonormal",
    "hypothesis_plan_independent",
    "hypothesis_plan_ordered_2",
    "hypothesis_plan_ordered_3",
    "hypothesis_plan_ordered_4",
    "hypothesis_plan_scaled_interaction",
    "hypothesis_plan_location_scale",
    "hypothesis_plan_model_averaged",
    "hypothesis_plan_robma_mixture",
    "hypothesis_plan_two_scale",
    "hypothesis_plan_meandif_mt",
    "marginal_means_expression",
    "marginal_means_zero_ordinate",
    "lane_carrier_treatment",
    "lane_carrier_meandif",
    "lane_carrier_nonzero",
    "lane_carrier_averaged",
    "lane_carrier_ordinary",
    "lane_carrier_gated",
    "lane_contribution_constant2",
    "lane_contribution_constant1",
    "lane_contribution_tau",
    "lane_contribution_point",
    "lane_contribution_discovery",
    "lane_exp_affine",
    "lane_factor_averaged_reference",
    "lane_factor_ordered_labels",
    "lane_factor_numeric_labels",
    "lane_scalar_standardized",
    "lane_scalar_raw",
    "lane_scalar_convolution",
    "lane_scalar_averaged",
    "lane_scalar_averaged_default",
    "lane_known_v_factor_scale",
    "lane_random_correlation",
    "lane_random_correlation_direct",
    "lane_random_correlation_gated",
    "lane_nested_total",
    "lane_selection_reference",
    "lane_glmm_point_baserate",
    "lane_glmm_boundary_beta",
    "lane_glmm_truncated_beta",
    "lane_glmm_truncated_lograte",
    "bias_branch_PP_intercept",
    "bias_branch_bPET_mods",
    "bias_branch_known_v_ensemble",
    "plot_prior_random_allocated",
    "plot_prior_random_scale"
  )
  expect_equal(catalog[["name"]][!catalog[["has_loo"]]], no_loo_fits,
               info = "selected brma.mv cache fixtures omit pre-computed LOO")
  # the weighted cluster LOO fixture stores both criteria: its test compares
  # the LOO and WAIC deletion scores
  expect_identical(catalog[["name"]][catalog[["has_waic"]]], "weighted_gaussian_cluster_loo",
                   info = "WAIC is pre-computed only for the weighted cluster LOO fixture")
  expect_true(all(!catalog[catalog[["class"]] %in% c("BMA.norm", "BMA.glmm", "BMA.mv", "RoBMA", "RoBMA.mv"), "has_marglik"]),
              info = "product-space model averaging fits must not cache marginal likelihoods")
  v14_marglik_fits <- c(
    "brma.mv_v14_konstantopoulos2011_cs",
    "brma.mv_v14_assink2016_nested",
    "brma.mv_v14_ishak2007_har",
    "brma.mv_v14_begg1989_study_treatment",
    "iwmde_known_v_tau_full",
    "iwmde_known_v_tau_null",
    "fixed_nonzero_brma_mv"
  )
  expect_equal(catalog[["name"]][catalog[["class"]] == "brma.mv" & catalog[["has_marglik"]]],
               v14_marglik_fits,
               info = "selected brma.mv cache fixtures pre-compute marginal likelihoods")
  no_marglik_single_fits <- c(
    "hypothesis_plan_treatment",
    "hypothesis_plan_meandif",
    "hypothesis_plan_orthonormal",
    "hypothesis_plan_independent",
    "hypothesis_plan_ordered_2",
    "hypothesis_plan_ordered_3",
    "hypothesis_plan_ordered_4",
    "hypothesis_plan_scaled_interaction",
    "hypothesis_plan_location_scale",
    "hypothesis_plan_meandif_mt",
    "marginal_means_expression",
    "marginal_means_zero_ordinate",
    "lane_carrier_treatment",
    "lane_carrier_meandif",
    "lane_carrier_nonzero",
    "lane_carrier_ordinary",
    "lane_carrier_gated",
    "lane_contribution_constant2",
    "lane_contribution_constant1",
    "lane_contribution_tau",
    "lane_contribution_point",
    "lane_contribution_discovery",
    "lane_exp_affine",
    "lane_factor_ordered_labels",
    "lane_factor_numeric_labels",
    "lane_scalar_standardized",
    "lane_scalar_raw",
    "lane_scalar_convolution",
    "lane_known_v_factor_scale",
    "lane_random_correlation",
    "lane_random_correlation_direct",
    "lane_selection_reference",
    "lane_glmm_point_baserate",
    "lane_glmm_boundary_beta",
    "lane_glmm_truncated_beta",
    "lane_glmm_truncated_lograte",
    "bias_branch_bPET_mods",
    "weighted_gaussian_cluster_loo"
  )
  expect_true(all(!catalog[catalog[["name"]] %in% no_marglik_single_fits, "has_marglik"]))
  expect_true(all(catalog[!catalog[["class"]] %in% c("BMA.norm", "BMA.glmm", "BMA.mv", "RoBMA", "RoBMA.mv", "brma.mv") &
                           !catalog[["name"]] %in% no_marglik_single_fits, "has_marglik"]),
              info = "implemented single-model cached fits include marginal likelihoods")
  expect_true(all(catalog[catalog[["class"]] %in% c("brma.glmm", "BMA.glmm"), "family"] == "glmm"),
              info = "GLMM classes must be tagged as glmm family")
  expect_true(all(catalog[!catalog[["class"]] %in% c("brma.glmm", "BMA.glmm"), "family"] == "norm"),
              info = "normal-outcome classes must be tagged as norm family")
  expect_equal(length(catalog[["features"]]), nrow(catalog),
               info = "each catalog row must have feature tags")
  expect_true(all(vapply(catalog[["features"]], length, integer(1)) > 0L),
              info = "feature tags cannot be empty")
  expect_setequal(unique(catalog[["profile"]]), TEST_PROFILES)
  expect_setequal(
    catalog[["name"]][catalog[["profile"]] == "standard"],
    c(
      "hypothesis_plan_treatment",
      "hypothesis_plan_meandif",
      "hypothesis_plan_orthonormal",
      "hypothesis_plan_independent",
      "hypothesis_plan_ordered_2",
      "hypothesis_plan_ordered_3",
      "hypothesis_plan_ordered_4",
      "hypothesis_plan_scaled_interaction",
      "hypothesis_plan_location_scale",
      "hypothesis_plan_model_averaged",
      "hypothesis_plan_robma_mixture",
      "hypothesis_plan_two_scale",
      "hypothesis_plan_meandif_mt",
      "marginal_means_expression",
      "marginal_means_zero_ordinate",
      "bcg_meta-analysis",
      "bcg_meta-regression",
      "bcg_meta-regression2",
      "bcg_meta-regression2b",
      "bcg_meta-regression3",
      "bcg_meta-regression3b",
      "bcg_meta-regression4",
      "bcg_meta-regression4b",
      "bangertdrowns2004_location-scale",
      "konstantopoulos2011_3lvl",
      "konstantopoulos2011_3lvl2",
      "bcg_glmm",
      "bcg_glmm_reg",
      "nielweise2008_glmm",
      "dat.lehmann2018-PET",
      "dat.lehmann2018-PEESE",
      "dat.lehmann2018-3PSM",
      "dat.lehmann2018_BMA.norm",
      "dat.lehmann2018_BMA.norm_mods",
      "bcg_BMA.glmm",
      "dat.lehmann2018_RoBMA",
      "dat.lehmann2018_RoBMA_mods",
      "dat.lehmann2018_RoBMA_mods2",
      "dat.lehmann2018_RoBMA_3lvl_mods_scale",
      "BMA.mv_random_components",
      "RoBMA.mv_marg_product_space",
      "bselmodel.mv_marg_random",
      "bPET.mv_random",
      "bPEESE.mv_random",
      "brma.mv_latent",
      "brma.mv_whitened",
      "brma.mv_block_mvn",
      "brma.mv_block_mvn_fixed_random_null",
      "brma.mv_block_mvn_random",
      "vif_parity_brma",
      "vif_parity_brma_mv",
      "fixed_null_brma",
      "fixed_nonzero_brma",
      "fixed_nonzero_brma_mv",
      "fixed_null_brma_glmm",
      "fixed_null_brma_glmm_pois",
      "fixed_null_bPET",
      "fixed_null_bPEESE",
      "fixed_null_bselmodel",
      "lane_carrier_treatment",
      "lane_carrier_meandif",
      "lane_carrier_nonzero",
      "lane_carrier_averaged",
      "lane_carrier_ordinary",
      "lane_carrier_gated",
      "lane_contribution_constant2",
      "lane_contribution_constant1",
      "lane_contribution_tau",
      "lane_contribution_point",
      "lane_contribution_discovery",
      "lane_exp_affine",
      "lane_factor_averaged_reference",
      "lane_factor_ordered_labels",
      "lane_factor_numeric_labels",
      "lane_scalar_standardized",
      "lane_scalar_raw",
      "lane_scalar_convolution",
      "lane_scalar_averaged",
      "lane_scalar_averaged_default",
      "lane_known_v_factor_scale",
      "lane_random_correlation",
      "lane_random_correlation_direct",
      "lane_random_correlation_gated",
      "lane_nested_total",
      "lane_selection_reference",
      "lane_glmm_point_baserate",
      "lane_glmm_boundary_beta",
      "lane_glmm_truncated_beta",
      "lane_glmm_truncated_lograte",
      "bias_branch_PP_mods",
      "bias_branch_PP_intercept",
      "bias_branch_bPET_mods",
      "bias_branch_known_v_ensemble",
      "plot_prior_random_allocated",
      "plot_prior_random_scale",
      "weighted_gaussian_cluster_loo"
    )
  )
})

test_that("fit catalog mirrors test-01 save_fit calls", {

  catalog    <- fit_catalog()
  saved_fits <- .test_01_saved_fits()

  expect_equal(anyDuplicated(saved_fits[["name"]]), 0L,
               info = "each cached fit name is saved once")
  expect_equal(setdiff(saved_fits[["name"]], catalog[["name"]]), character(),
               info = "all save_fit() names must be catalogued")
  expect_equal(setdiff(catalog[["name"]], saved_fits[["name"]]), character(),
               info = "all catalogued fits must be saved by a test-01 file")

  saved_source <- saved_fits[["source_file"]][match(catalog[["name"]], saved_fits[["name"]])]
  expect_equal(saved_source, catalog[["source_file"]],
               info = "catalog source_file must match the test-01 file calling save_fit()")
})

test_that("BayesTools cache fingerprint fails closed", {

  fingerprint <- BayesTools::fit_backend_fingerprint()
  expect_identical(bayestools_backend_fingerprint(), fingerprint)

  invalid_values <- list(NULL, NA_character_, "", "not-a-fingerprint", 1)
  for (value in invalid_values) {
    expect_error(
      bayestools_backend_fingerprint(value = value),
      "unavailable or invalid",
      fixed = TRUE
    )
    expect_identical(
      bayestools_backend_fingerprint(required = FALSE, value = value),
      NA_character_
    )
  }
})


test_that("cache validation rejects corrupted synthetic metadata", {

  name     <- "bcg_meta-analysis"
  fit      <- .fake_fit_for_catalog_entry(name)
  info     <- .fake_info_for_catalog_entry(name)
  metadata <- .valid_metadata_for_catalog_entry(name)

  expect_equal(
    validate_cached_fit(
      name        = name,
      fit         = fit,
      info        = info,
      metadata    = metadata,
      check_files = FALSE,
      deep        = TRUE
    ),
    character(),
    info = "valid synthetic metadata passes without cache files"
  )
  expect_identical(names(metadata), FIT_CACHE_METADATA_FIELDS)
  expect_identical(
    metadata[["robma_version"]],
    as.character(utils::packageVersion("RoBMA"))
  )
  expect_identical(
    metadata[["bayestools_version"]],
    as.character(utils::packageVersion("BayesTools"))
  )
  expect_identical(
    metadata[["r_version"]],
    paste(R.version$major, R.version$minor, sep = ".")
  )
  expect_true(
    is.na(metadata[["jags_version"]]) || nzchar(metadata[["jags_version"]])
  )

  invalid_metadata <- validate_cached_fit(
    name        = name,
    fit         = fit,
    info        = info,
    metadata    = "not-a-list",
    check_files = FALSE,
    deep        = TRUE
  )
  expect_true("metadata must be a list" %in% invalid_metadata)

  stale_schema <- metadata
  stale_schema[["schema_version"]] <- FIT_CACHE_SCHEMA_VERSION + 1L
  expect_true(
    "cache schema changed" %in% validate_cached_fit(
      name        = name,
      fit         = fit,
      info        = info,
      metadata    = stale_schema,
      check_files = FALSE,
      deep        = TRUE
    )
  )

  missing_schema <- metadata
  missing_schema[["schema_version"]] <- NULL
  missing_schema_problems <- validate_cached_fit(
    name        = name,
    fit         = fit,
    info        = info,
    metadata    = missing_schema,
    check_files = FALSE,
    deep        = TRUE
  )
  expect_true("cache metadata fields changed" %in% missing_schema_problems)
  expect_true("cache schema changed" %in% missing_schema_problems)

  version_messages <- c(
    robma_version      = "RoBMA version changed",
    bayestools_version = "BayesTools version changed",
    r_version          = "R version changed",
    jags_version       = "JAGS version changed"
  )
  for (field in names(version_messages)) {
    stale_version <- metadata
    stale_version[[field]] <- "stale-version"
    expect_true(
      unname(version_messages[[field]]) %in% validate_cached_fit(
        name         = name,
        fit          = fit,
        info         = info,
        metadata     = stale_version,
        check_source = FALSE,
        check_files  = FALSE,
        deep         = TRUE
      ),
      info = field
    )
  }

  unexpected_field <- metadata
  unexpected_field[["obsolete"]] <- TRUE
  expect_true(
    "cache metadata fields changed" %in% validate_cached_fit(
      name        = name,
      fit         = fit,
      info        = info,
      metadata    = unexpected_field,
      check_files = FALSE,
      deep        = TRUE
    )
  )

  stale_name <- metadata
  stale_name[["name"]] <- "wrong-name"
  expect_true(
    "metadata name mismatch" %in% validate_cached_fit(
      name        = name,
      fit         = fit,
      info        = info,
      metadata    = stale_name,
      check_files = FALSE,
      deep        = TRUE
    )
  )

  stale_source_file <- metadata
  stale_source_file[["source_file"]] <- "test-01-wrong.R"
  expect_true(
    "metadata source file mismatch" %in% validate_cached_fit(
      name        = name,
      fit         = fit,
      info        = info,
      metadata    = stale_source_file,
      check_files = FALSE,
      deep        = TRUE
    )
  )

  unavailable_backend <- validate_cached_fit(
    name                   = name,
    fit                    = fit,
    info                   = info,
    metadata               = metadata,
    check_files            = FALSE,
    deep                   = TRUE,
    bayestools_fingerprint = NA_character_
  )
  expect_true(
    "BayesTools backend fingerprint unavailable" %in% unavailable_backend
  )

  unavailable_source_validator <- validate_cached_fit
  environment(unavailable_source_validator) <- list2env(
    list(source_file_md5 = function(...) NA_character_),
    parent = environment(validate_cached_fit)
  )
  unavailable_source <- unavailable_source_validator(
    name = name, metadata = metadata, check_files = FALSE
  )
  expect_true("source file fingerprint unavailable" %in% unavailable_source)
  unavailable_package <- validate_cached_fit(
    name = name, metadata = metadata, check_files = FALSE,
    package_md5 = NA_character_
  )
  expect_true("cache source fingerprint unavailable" %in% unavailable_package)

  malformed_backend <- validate_cached_fit(
    name                   = name,
    fit                    = fit,
    info                   = info,
    metadata               = metadata,
    check_files            = FALSE,
    deep                   = TRUE,
    bayestools_fingerprint = "malformed"
  )
  expect_true(
    "BayesTools backend fingerprint unavailable" %in% malformed_backend
  )

  stale_source <- metadata
  stale_source[["source_file_md5"]] <- "not-current"
  expect_true(
    "source file hash changed" %in% validate_cached_fit(
      name        = name,
      fit         = fit,
      info        = info,
      metadata    = stale_source,
      check_files = FALSE,
      deep        = TRUE
    ),
    info = "test source hash is enforced"
  )

  stale_package <- metadata
  stale_package[["package_source_md5"]] <- "not-current"
  if (!is.na(package_source_md5())) {
    expect_true(
      "cache source hash changed" %in% validate_cached_fit(
        name        = name,
        fit         = fit,
        info        = info,
        metadata    = stale_package,
        check_files = FALSE,
        deep        = TRUE
      ),
      info = "cache source hash is enforced"
    )
  }

  stale_bayestools <- metadata
  stale_bayestools[["bayestools_backend_fingerprint"]] <- "not-current"
  if (!is.na(BayesTools::fit_backend_fingerprint())) {
    expect_true(
      "BayesTools backend changed" %in% validate_cached_fit(
        name        = name,
        fit         = fit,
        info        = info,
        metadata    = stale_bayestools,
        check_files = FALSE,
        deep        = TRUE
      ),
      info = "BayesTools backend fingerprint is enforced"
    )

    missing_bayestools <- metadata
    missing_bayestools[["bayestools_backend_fingerprint"]] <- NULL
    expect_true(
      "BayesTools backend changed" %in% validate_cached_fit(
        name        = name,
        fit         = fit,
        info        = info,
        metadata    = missing_bayestools,
        check_files = FALSE,
        deep        = TRUE
      ),
      info = "missing BayesTools backend fingerprint is rejected"
    )
  }

  missing_marglik <- metadata
  missing_marglik[["has_marglik"]] <- FALSE
  expect_true(
    "metadata reports missing marginal likelihood" %in% validate_cached_fit(
      name        = name,
      fit         = fit,
      info        = info,
      metadata    = missing_marglik,
      check_files = FALSE,
      deep        = TRUE
    ),
    info = "required marginal likelihood metadata is enforced"
  )

  missing_metafor <- metadata
  missing_metafor[["has_metafor_info"]] <- FALSE
  expect_true(
    "metadata reports missing metafor reference" %in% validate_cached_fit(
      name        = name,
      fit         = fit,
      info        = info,
      metadata    = missing_metafor,
      check_files = FALSE,
      deep        = TRUE
    ),
    info = "required metafor metadata is enforced"
  )

  product_name     <- "dat.lehmann2018_RoBMA"
  product_fit      <- .fake_fit_for_catalog_entry(product_name)
  product_info     <- .fake_info_for_catalog_entry(product_name)
  product_metadata <- .valid_metadata_for_catalog_entry(product_name)
  product_metadata[["has_marglik"]] <- TRUE
  product_fit[["marglik"]] <- list(estimate = "unexpected")

  problems <- validate_cached_fit(
    name        = product_name,
    fit         = product_fit,
    info        = product_info,
    metadata    = product_metadata,
    check_files = FALSE,
    deep        = TRUE
  )
  expect_true(any(grepl("unexpected marginal likelihood", problems)),
              info = "product-space marginal likelihoods are forbidden")
})

test_that("cached fits are catalogued and valid", {

  skip_if(is_certification_profile(),
          "A certification case activates only its own fits; the standard profile checks the cache inventory.")
  skip_on_cran()

  cached_names <- list_fits(validate = FALSE, active_only = FALSE)
  if (length(cached_names) == 0) {
    skip("No cached fits available for catalog validation.")
  }

  catalog_names <- fit_catalog()[["name"]]
  uncatalogued  <- setdiff(cached_names, catalog_names)
  inactive      <- setdiff(cached_names, active_fit_catalog()[["name"]])

  expect_equal(uncatalogued, character(),
               info = "all cached fits must be listed in fit_catalog()")
  expect_equal(inactive, character(),
               info = "the active cache must not contain fits from another profile")

  # One package-source fingerprint serves the whole batch, as in list_fits().
  package_md5 <- package_source_md5()
  problems <- lapply(cached_names, validate_cached_fit, deep = TRUE,
                     package_md5 = package_md5)
  names(problems) <- cached_names
  problems <- problems[lengths(problems) > 0]
  problem_text <- if (length(problems) == 0) {
    ""
  } else {
    paste(
      paste0(names(problems), ": ", vapply(problems, paste, "", collapse = "; ")),
      collapse = "\n"
    )
  }

  expect_true(
    length(problems) == 0,
    info = paste(
      "cached fits must match catalog requirements; rerun test-01 after cache-affecting fitting changes",
      problem_text,
      sep = "\n"
    )
  )
})

.review_binary_text_md5 <- function(text) {

  path <- tempfile("robma-binary-reference-")
  on.exit(unlink(path), add = TRUE)
  writeBin(charToRaw(enc2utf8(text)), path)
  unname(tools::md5sum(path))
}


test_that("cache text hashes use binary UTF-8 LF and preserve empty semantics", {

  root <- tempfile("robma-source-text-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  testthat::local_test_directory(root, .env = environment())
  for (text in c("a\nb\n", "\u00e9\n\u010d\n", "", "\n")) {
    paths <- file.path(root, c("lf.cc", "crlf.cc"))
    writeBin(charToRaw(enc2utf8(text)), paths[[1L]])
    writeBin(charToRaw(enc2utf8(gsub("\n", "\r\n", text, fixed = TRUE))), paths[[2L]])
    original <- tools::md5sum(paths)
    expected <- .review_binary_text_md5(text)
    expect_identical(unname(vapply(paths, .fit_cache_source_file_md5, character(1L))),
      rep(expected, 2L))
    expect_identical(tools::md5sum(paths), original)
  }
  values <- c("1.0000000000000000", "1.0000000000000002")
  paths <- file.path(root, c("first.R", "second.R"))
  for (j in 1:2) writeBin(charToRaw(paste0("value <- ", values[[j]], "\n")), paths[[j]])
  parsed <- vapply(paths, function(path) eval(parse(path)[[1L]][[3L]]), numeric(1L))
  expect_false(identical(parsed[[1L]], parsed[[2L]]))
  expect_false(identical(.fit_cache_source_file_md5(paths[[1L]]),
    .fit_cache_source_file_md5(paths[[2L]])))
  expect_false(identical(source_file_md5("first.R"), source_file_md5("second.R")))
  for (text in c("", "# comment only\n")) {
    path <- file.path(root, "empty.R")
    writeBin(charToRaw(text), path)
    expect_identical(.fit_cache_source_file_md5(path), .review_binary_text_md5(""))
    expect_identical(source_file_md5("empty.R"), .review_binary_text_md5(""))
  }
})


test_that("package source aggregates hash the declared LF bytes without unavailable values", {

  skip_on_cran()
  .with_cache_source_fixture(function(e, root) {
    files <- e$.fit_cache_source_files(package_root = root)
    original <- lapply(files, readLines, warn = FALSE, encoding = "UTF-8")
    observed <- character()
    for (ending in c("\n", "\r\n")) {
      for (j in seq_along(files)) {
        text <- if (length(original[[j]])) paste0(paste(original[[j]], collapse = ending), ending) else ""
        writeBin(charToRaw(enc2utf8(text)), files[[j]])
      }
      e$.clear_package_source_md5_cache()
      hashes <- vapply(files, e$.fit_cache_source_file_md5, character(1L))
      relative <- substring(files, nchar(root) + 2L)
      expected <- .review_binary_text_md5(paste0(paste(paste0(relative, ":", hashes), collapse = "\n"), "\n"))
      observed <- c(observed, e$package_source_md5())
      expect_identical(tail(observed, 1L), expected)
    }
    expect_identical(observed[[1L]], observed[[2L]])
    e$.clear_package_source_md5_cache()
    assign(".fit_cache_source_file_md5", function(path) NA_character_, e)
    expect_identical(e$package_source_md5(), NA_character_)
  })
})


test_that("text references retain scalar CRLF and missing-candidate policy", {

  root <- tempfile("robma-reference-text-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  baseline <- file.path(root, "reference.txt")
  writeBin(charToRaw("a\nb\n"), baseline)
  before <- readBin(baseline, "raw", n = file.size(baseline))
  test_reference_text("a\r\nb", "reference.txt", print_dir = root)
  expect_identical(readBin(baseline, "raw", n = file.size(baseline)), before)
  writeBin(charToRaw("a\r\nb\r\n"), baseline)
  test_reference_text("a\nb", "reference.txt", print_dir = root)
  for (invalid in list(c("a", "ignored"), character(), NA_character_, 1)) {
    expect_error(test_reference_text(invalid, "reference.txt", print_dir = root),
      "'text' must be one non-missing character string.", fixed = TRUE)
  }
  missing <- tryCatch(test_reference_text("new", "missing.txt", print_dir = root), skip = identity)
  expect_s3_class(missing, "skip")
  expect_false(file.exists(file.path(root, "missing.txt")))
  candidate <- .reference_candidate_path(file.path(root, "missing.txt"))
  expect_true(file.exists(candidate))
  expect_match(conditionMessage(missing), paste0("Reference file missing.txt not found. ",
    "Reference candidate cached at ", candidate, ". Review it with 'review_test_snapshots()'."),
    fixed = TRUE)
  if (file.exists(candidate)) expect_identical(readLines(candidate), "new")
  candidate <- .reference_candidate_path(baseline)
  writeLines("old candidate", candidate)
  before <- readBin(baseline, "raw", n = file.size(baseline))
  test_reference_text("a\nb", "reference.txt", print_dir = root)
  expect_false(file.exists(candidate))
  expect_identical(readBin(baseline, "raw", n = file.size(baseline)), before)
})


test_that("owned cache cleanup checks unsuccessful filesystem operations", {

  .with_temp_fit_cache({
    name <- "bcg_meta-analysis"
    paths <- unlist(fit_cache_paths(name), use.names = FALSE)
    for (path in paths) {
      dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
      writeLines("owned", path)
    }
    neighbor <- file.path(test_files_dir, "unrelated.txt")
    writeLines("retained", neighbor)
    clean <- clean_cached_fits
    environment(clean) <- list2env(list(file.remove = function(paths) rep(FALSE, length(paths))),
      parent = environment(clean))
    expect_error(clean(name), "Could not remove cached artifacts", fixed = TRUE)
    expect_true(all(file.exists(paths)))
    clean_cached_fits(name)
    expect_false(any(file.exists(paths)))
    expect_true(file.exists(neighbor))
    clean <- clean_cached_fits
    environment(clean) <- list2env(list(unlink = function(...) 1L), parent = environment(clean))
    expect_error(clean(), "Could not remove all owned cache directories", fixed = TRUE)
    clean_cached_fits()
    expect_true(all(dir.exists(c(temp_fits_dir, temp_info_dir, temp_metadata_dir, temp_temp_dir))))
    expect_identical(readLines(neighbor), "retained")
  })
})
