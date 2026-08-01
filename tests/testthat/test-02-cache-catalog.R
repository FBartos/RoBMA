context("Cached fit catalog")

source(testthat::test_path("common-functions.R"))

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

  entry <- fit_catalog_entry(name)

  return(list(
    schema_version                 = FIT_CACHE_SCHEMA_VERSION,
    name                           = name,
    saved_at                       = "synthetic",
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


test_that("cache source hash tracks only cache-affecting fitting sources", {

  source_files <- .fit_cache_source_files(relative = TRUE)
  skip_if_not(
    length(source_files) > 0L,
    "Complete package R/src source tree is not available under source-build checks."
  )

  expect_equal(anyDuplicated(source_files), 0L)
  expect_true(all(c(
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
    "src/r-glmm.cc",
    "src/r-selnorm.cc",
    "src/r-selnorm-common.cc.inc",
    "src/r-selnorm-funnel-zcurve.cc.inc",
    "src/r-selnorm-kernel.cc.inc",
    "src/r-selnorm-loglik.cc.inc",
    "src/selnorm/selnorm-api.cc.inc",
    "src/selnorm/selnorm-boundary.cc.inc",
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
    "src/r-regplot.cc"
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
    "dat.lehmann2018-3PSM_effect_null"
  )
  expect_equal(catalog[["name"]][!catalog[["has_loo"]]], no_loo_fits,
               info = "selected brma.mv cache fixtures omit pre-computed LOO")
  expect_true(all(!catalog[["has_waic"]]),
              info = "WAIC is not pre-computed for cached test fits")
  expect_true(all(!catalog[catalog[["class"]] %in% c("BMA.norm", "BMA.glmm", "RoBMA"), "has_marglik"]),
              info = "product-space model averaging fits must not cache marginal likelihoods")
  v14_marglik_fits <- c(
    "brma.mv_v14_konstantopoulos2011_cs",
    "brma.mv_v14_assink2016_nested",
    "brma.mv_v14_ishak2007_har",
    "brma.mv_v14_begg1989_study_treatment",
    "iwmde_known_v_tau_full",
    "iwmde_known_v_tau_null"
  )
  expect_equal(catalog[["name"]][catalog[["class"]] == "brma.mv" & catalog[["has_marglik"]]],
               v14_marglik_fits,
               info = "selected brma.mv cache fixtures pre-compute marginal likelihoods")
  expect_true(all(catalog[!catalog[["class"]] %in% c("BMA.norm", "BMA.glmm", "RoBMA", "brma.mv"), "has_marglik"]),
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
      "brma.mv_latent",
      "brma.mv_whitened",
      "brma.mv_block_mvn",
      "brma.mv_block_mvn_fixed_random_null",
      "brma.mv_block_mvn_random",
      "vif_parity_brma",
      "vif_parity_brma_mv"
    )
  )
  expect_true(all(active_fit_catalog()[["profile"]] %in%
                    if (is_certification_profile()) TEST_PROFILES else "standard"))
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

  skip_on_cran()

  cached_names <- list_fits(validate = FALSE)
  if (length(cached_names) == 0) {
    skip("No cached fits available for catalog validation.")
  }

  catalog_names <- fit_catalog()[["name"]]
  uncatalogued  <- setdiff(cached_names, catalog_names)

  expect_equal(uncatalogued, character(),
               info = "all cached fits must be listed in fit_catalog()")

  problems <- lapply(cached_names, validate_cached_fit, deep = TRUE)
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
