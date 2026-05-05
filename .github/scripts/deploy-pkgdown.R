excluded_pkgdown_pages <- c(
  "AGENTS.html"
)

pkgdown_destination_path <- function(pkg, dest_dir = "docs") {

  if (inherits(pkg, "pkgdown") && !is.null(pkg$dst_path)) {
    return(pkg$dst_path)
  } else {
    return(pkgdown:::as_pkgdown(pkg, override = list(destination = dest_dir))$dst_path)
  }
}

remove_excluded_pkgdown_pages <- function(pkg, dest_dir = "docs") {

  dst_path <- pkgdown_destination_path(pkg, dest_dir = dest_dir)

  paths <- file.path(dst_path, excluded_pkgdown_pages)
  paths <- paths[file.exists(paths)]

  if (length(paths) > 0) {
    unlink(paths, force = TRUE)
    message("Removed excluded pkgdown pages: ", paste(basename(paths), collapse = ", "))
  }

  invisible(paths)
}

remove_excluded_sitemap_urls <- function(pkg, dest_dir = "docs") {

  sitemap_path <- file.path(pkgdown_destination_path(pkg, dest_dir = dest_dir), "sitemap.xml")

  if (!file.exists(sitemap_path)) {
    return(invisible(FALSE))
  }

  sitemap_lines <- readLines(sitemap_path, warn = FALSE)
  page_pattern  <- paste(gsub(".", "\\.", excluded_pkgdown_pages, fixed = TRUE), collapse = "|")
  keep_lines    <- !grepl(page_pattern, sitemap_lines)

  if (!all(keep_lines)) {
    writeLines(sitemap_lines[keep_lines], sitemap_path, useBytes = TRUE)
    message("Removed excluded pkgdown pages from sitemap.")
  }

  invisible(TRUE)
}

patch_pkgdown_deploy_cleanup <- function() {

  pkgdown_namespace <- asNamespace("pkgdown")
  original_build    <- get("build_site_github_pages", envir = pkgdown_namespace)

  patched_build <- function(pkg = ".", ..., dest_dir = "docs", clean = TRUE, install = FALSE, new_process = FALSE) {

    original_build(
      pkg         = pkg,
      ...,
      dest_dir    = dest_dir,
      clean       = clean,
      install     = install,
      new_process = new_process
    )

    remove_excluded_pkgdown_pages(pkg, dest_dir = dest_dir)
    remove_excluded_sitemap_urls(pkg, dest_dir = dest_dir)
    invisible()
  }

  unlockBinding("build_site_github_pages", pkgdown_namespace)
  assign("build_site_github_pages", patched_build, envir = pkgdown_namespace)
  lockBinding("build_site_github_pages", pkgdown_namespace)

  invisible()
}

patch_pkgdown_deploy_cleanup()
pkgdown::deploy_to_branch(new_process = FALSE)
