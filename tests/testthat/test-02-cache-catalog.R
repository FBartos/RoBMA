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
    version            = FIT_CACHE_VERSION,
    name               = name,
    saved_at           = "synthetic",
    fit_class          = entry[["class"]],
    source_file        = entry[["source_file"]],
    source_file_md5    = source_file_md5(entry[["source_file"]]),
    package_source_md5 = package_source_md5(),
    has_loo            = entry[["has_loo"]],
    has_waic           = entry[["has_waic"]],
    has_marglik        = entry[["has_marglik"]],
    has_metafor_info   = entry[["has_metafor"]]
  ))
}


test_that("cache source hash tracks only cache-affecting fitting sources", {

  source_files <- .fit_cache_source_files(relative = TRUE)

  expect_equal(anyDuplicated(source_files), 0L)
  expect_true(all(c(
    "R/fit.R",
    "R/input-data.R",
    "R/input-object.R",
    "R/input-priors.R",
    "R/selection-mapping.R",
    "R/marglik.R",
    "R/loo.R",
    "src/RoBMA.cc",
    "src/distributions/DWN.cc",
    "src/r-glmm.cc",
    "src/r-selnorm.cc"
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
  expect_true(all(catalog[["has_loo"]]),
              info = "all cached test fits include pre-computed LOO")
  expect_true(all(!catalog[["has_waic"]]),
              info = "WAIC is not pre-computed for cached test fits")
  expect_true(all(!catalog[catalog[["class"]] %in% c("BMA.norm", "BMA.glmm", "RoBMA"), "has_marglik"]),
              info = "product-space model averaging fits must not cache marginal likelihoods")
  expect_true(all(catalog[!catalog[["class"]] %in% c("BMA.norm", "BMA.glmm", "RoBMA"), "has_marglik"]),
              info = "single-model cached fits include marginal likelihoods")
  expect_true(all(catalog[catalog[["class"]] %in% c("brma.glmm", "BMA.glmm"), "family"] == "glmm"),
              info = "GLMM classes must be tagged as glmm family")
  expect_true(all(catalog[!catalog[["class"]] %in% c("brma.glmm", "BMA.glmm"), "family"] == "norm"),
              info = "normal-outcome classes must be tagged as norm family")
  expect_equal(length(catalog[["features"]]), nrow(catalog),
               info = "each catalog row must have feature tags")
  expect_true(all(vapply(catalog[["features"]], length, integer(1)) > 0L),
              info = "feature tags cannot be empty")
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
