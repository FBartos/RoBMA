vignette_cache <- function(name, objects, packages = character(), root = "../models") {

  stopifnot(
    is.character(name), length(name) == 1L,
    is.character(objects), length(objects) > 0L,
    is.character(packages),
    is.character(root), length(root) == 1L
  )

  object_names <- names(objects)
  file_names   <- unname(objects)

  if (is.null(object_names)) {
    object_names <- rep("", length(objects))
  }

  unnamed_objects <- !nzchar(object_names)
  object_names[unnamed_objects] <- sub(
    "\\.[Rr][Dd][Ss]$",
    "",
    basename(file_names[unnamed_objects])
  )

  missing_extension <- !grepl("\\.[Rr][Dd][Ss]$", file_names)
  file_names[missing_extension] <- paste0(file_names[missing_extension], ".RDS")

  paths <- file.path(root, name, file_names)
  names(paths) <- object_names

  structure(
    paths,
    class    = "vignette_cache",
    name     = name,
    packages = packages,
    root     = root
  )
}

vignette_cache_checking <- function() {

  ("CheckExEnv" %in% search()) ||
    any(c("_R_CHECK_TIMINGS_", "_R_CHECK_LICENSE_") %in% names(Sys.getenv()))
}

vignette_cache_missing <- function(cache) {

  names(cache)[!file.exists(unname(cache))]
}

vignette_cache_missing_packages <- function(cache) {

  packages <- attr(cache, "packages")
  if (is.null(packages)) {
    packages <- character()
  }

  has_package <- vapply(
    packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )

  packages[!has_package]
}

vignette_cache_eval <- function(cache) {

  !vignette_cache_checking() &&
    !length(vignette_cache_missing_packages(cache)) &&
    !length(vignette_cache_missing(cache))
}

vignette_cache_load <- function(cache, envir = parent.frame()) {

  missing <- vignette_cache_missing(cache)
  if (length(missing)) {
    stop(
      "Missing cached vignette fits: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  for (i in seq_along(cache)) {
    assign(
      names(cache)[i],
      readRDS(file = unname(cache[i])),
      envir = envir
    )
  }

  invisible(cache)
}

vignette_cache_save <- function(cache, envir = parent.frame(), compress = "xz") {

  for (cache_dir in unique(dirname(unname(cache)))) {
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  }

  has_object <- vapply(
    names(cache),
    exists,
    logical(1),
    envir    = envir,
    inherits = FALSE
  )

  if (any(!has_object)) {
    stop(
      "Missing objects to cache: ",
      paste(names(cache)[!has_object], collapse = ", "),
      call. = FALSE
    )
  }

  for (i in seq_along(cache)) {
    saveRDS(
      get(names(cache)[i], envir = envir, inherits = FALSE),
      file     = unname(cache[i]),
      compress = compress
    )
  }

  invisible(cache)
}
