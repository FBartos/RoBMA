context("JAGS 4.x build requirements")

# RoBMA's JAGS module implements the JAGS 4 module interface only. These tests
# run the build scripts of a source checkout against fake JAGS installations:
# headers (Console.h, module/Module.h, and a version.h declaring JAGS_MAJOR)
# and empty library directories. Nothing is linked against them.

.jags_build_repository_file <- function(...) {

  testthat::test_path("..", "..", ...)
}


.jags_build_write_lf <- function(lines, path) {

  writeBin(charToRaw(paste0(paste(lines, collapse = "\n"), "\n")), path)
}


.jags_build_fake_tree <- function(root, header_major,
                                  include_dir = file.path("include", "JAGS"),
                                  lib_dir = "lib") {

  headers <- file.path(root, include_dir)
  dir.create(file.path(headers, "module"), recursive = TRUE, showWarnings = FALSE)
  file.create(file.path(headers, "Console.h"))
  file.create(file.path(headers, "module", "Module.h"))
  .jags_build_write_lf(
    c("#ifndef JAGS_VERSION_H_", "#define JAGS_VERSION_H_", "",
      paste("#define JAGS_MAJOR", header_major), "", "#endif"),
    file.path(headers, "version.h")
  )
  dir.create(file.path(root, lib_dir), recursive = TRUE, showWarnings = FALSE)

  normalizePath(root, winslash = "/", mustWork = TRUE)
}


.jags_build_env_unset <- function() {

  c(JAGS_ROOT = NA, JAGS_PREFIX = NA, JAGS_INCLUDEDIR = NA, JAGS_LIBDIR = NA,
    JAGS_VERSION = NA, JAGS_MAJOR_VERSION = NA, PKG_CONFIG_PATH = NA)
}


.jags_build_run_tool <- function(command, args, directory, env) {

  withr::local_envvar(env)
  withr::local_dir(directory)
  output <- suppressWarnings(system2(command, args, stdout = TRUE, stderr = TRUE))
  status <- attr(output, "status", exact = TRUE)

  list(
    status = if (is.null(status)) 0L else as.integer(status),
    output = paste(output, collapse = "\n")
  )
}


test_that("configure accepts only JAGS 4.x from pkg-config and the headers", {

  skip_on_cran()
  configure_file <- .jags_build_repository_file("configure")
  makevars_file  <- .jags_build_repository_file("src", "Makevars.in")
  skip_if_not(
    file.exists(configure_file) && file.exists(makevars_file),
    "Repository build sources are not available in this installed-package context."
  )
  sh <- unname(Sys.which("sh"))
  skip_if(!nzchar(sh), "A POSIX shell is required to run configure.")

  work <- normalizePath(withr::local_tempdir(), winslash = "/", mustWork = TRUE)
  skip_if(grepl("[[:space:]]", work), "Fake JAGS paths are space-separated shell words.")
  package_dir <- file.path(work, "pkg")
  dir.create(file.path(package_dir, "src"), recursive = TRUE)
  file.copy(configure_file, file.path(package_dir, "configure"))
  file.copy(makevars_file, file.path(package_dir, "src", "Makevars.in"))

  jags4 <- .jags_build_fake_tree(file.path(work, "jags4"), 4L)
  jags5 <- .jags_build_fake_tree(file.path(work, "jags5"), 5L)

  # pkg-config stubs placed first on PATH: one without JAGS metadata (so a
  # system jags.pc cannot leak into the legacy configuration) and ones that
  # report a JAGS version with the given headers.
  stub_none <- file.path(work, "stub-none")
  dir.create(stub_none)
  .jags_build_write_lf(c("#!/bin/sh", "exit 1"), file.path(stub_none, "pkg-config"))
  Sys.chmod(file.path(stub_none, "pkg-config"), mode = "0755")
  stub_jags <- function(version, atleast, tree) {

    stub <- file.path(work, paste0("stub-", version, "-", basename(tree)))
    dir.create(stub, showWarnings = FALSE)
    .jags_build_write_lf(c(
      "#!/bin/sh",
      "case \"$*\" in",
      paste0("  *--atleast-version*) exit ", if (atleast) 0L else 1L, " ;;"),
      "  *--exists*) exit 0 ;;",
      paste0("  *--modversion*) echo \"", version, "\" ;;"),
      paste0("  *--cflags*) echo \"-I", tree, "/include/JAGS\" ;;"),
      paste0("  *--libs*) echo \"-L", tree, "/lib -ljags\" ;;"),
      paste0("  *--variable*) echo \"", tree, "/lib\" ;;"),
      "esac"
    ), file.path(stub, "pkg-config"))
    Sys.chmod(file.path(stub, "pkg-config"), mode = "0755")
    stub
  }

  run_configure <- function(args = character(), stub = stub_none) {

    unlink(file.path(package_dir, "src", "Makevars"))
    path <- paste(stub, Sys.getenv("PATH"), sep = .Platform$path.sep)
    result <- .jags_build_run_tool(
      sh, c("./configure", args), package_dir,
      c(.jags_build_env_unset(), LIBnn = "lib", PATH = path)
    )
    result$makevars <- file.exists(file.path(package_dir, "src", "Makevars"))
    result
  }
  requirement <- "RoBMA requires JAGS 4.x (>= 4.3.1, < 5.0.0); "
  expect_rejected <- function(result, reported) {

    expect_false(identical(result$status, 0L), info = result$output)
    expect_false(result$makevars, info = result$output)
    expect_match(result$output, paste0(requirement, reported), fixed = TRUE)
  }
  # The fake libraries cannot be linked, so an accepted JAGS 4 installation
  # passes every version check and stops only at the link sanity check.
  expect_version_accepted <- function(result) {

    expect_false(grepl(requirement, result$output, fixed = TRUE), info = result$output)
    expect_match(result$output, "cannot link to JAGS library", fixed = TRUE)
  }

  expect_rejected(
    run_configure(stub = stub_jags("5.0.0", TRUE, jags4)),
    "pkg-config reported 5.0.0"
  )
  expect_rejected(
    run_configure(stub = stub_jags("4.2.0", FALSE, jags4)),
    "pkg-config reported 4.2.0"
  )
  expect_rejected(
    run_configure(stub = stub_jags("4.3.2", TRUE, jags5)),
    "the JAGS headers in"
  )
  expect_version_accepted(run_configure(stub = stub_jags("4.3.2", TRUE, jags4)))

  # Without pkg-config metadata the legacy configuration reads the headers.
  expect_rejected(
    run_configure(args = paste0("--with-jags-prefix=", jags5)),
    paste0("the JAGS headers in '", jags5, "/include/JAGS' declare JAGS_MAJOR 5")
  )
  expect_version_accepted(run_configure(args = paste0("--with-jags-prefix=", jags4)))
})


test_that("Windows builds select only JAGS 4.x installations", {

  skip_on_cran()
  makevars_files <- c(
    ucrt = .jags_build_repository_file("src", "Makevars.ucrt"),
    win  = .jags_build_repository_file("src", "Makevars.win")
  )
  skip_if_not(
    all(file.exists(makevars_files)),
    "Repository build sources are not available in this installed-package context."
  )
  make <- unname(Sys.which("make"))
  skip_if(!nzchar(make), "GNU make is required to evaluate the Windows Makevars.")

  work <- normalizePath(withr::local_tempdir(), winslash = "/", mustWork = TRUE)
  skip_if(grepl("[[:space:]]", work), "Candidate JAGS roots are space-separated make words.")
  install_root <- function(name, header_major) {

    .jags_build_fake_tree(
      file.path(work, "JAGS", name), header_major,
      include_dir = "include", lib_dir = file.path("x64", "bin")
    )
  }
  roots <- c(
    "4.2.0"  = install_root("JAGS-4.2.0", 4L),
    "4.3.1"  = install_root("JAGS-4.3.1", 4L),
    "4.3.2"  = install_root("JAGS-4.3.2", 4L),
    "5.0.0"  = install_root("JAGS-5.0.0", 5L),
    custom4  = install_root("custom4", 4L),
    custom5  = install_root("custom5", 5L)
  )

  for (variant in names(makevars_files)) {
    file.copy(makevars_files[[variant]], file.path(work, paste0("Makevars.", variant)))
    .jags_build_write_lf(c(
      "print-selection:",
      "\t@echo \"selected $(notdir $(JAGS_ROOT)) version $(JAGS_VERSION) link -ljags-$(JAGS_MAJOR) assumed $(JAGS_MAJOR_ASSUMED)\"",
      paste0("include Makevars.", variant)
    ), file.path(work, paste0("probe-", variant, ".mk")))

    run_make <- function(candidates, vars = character()) {

      .jags_build_run_tool(
        make,
        c("-s", "-f", paste0("probe-", variant, ".mk"),
          shQuote(paste0("JAGS_INSTALL_ROOTS=", paste(roots[candidates], collapse = " "))),
          vapply(vars, shQuote, character(1)), "print-selection"),
        work,
        .jags_build_env_unset()
      )
    }
    expect_selected <- function(result, selection) {

      expect_identical(result$status, 0L, info = paste(variant, result$output))
      expect_match(result$output, selection, fixed = TRUE, info = variant)
    }
    expect_stopped <- function(result, message) {

      expect_false(identical(result$status, 0L), info = paste(variant, result$output))
      expect_match(result$output, message, fixed = TRUE, info = variant)
    }
    requirement <- "RoBMA requires JAGS 4.x (>= 4.3.1, < 5.0.0)"

    expect_selected(
      run_make(c("4.2.0", "4.3.1", "5.0.0")),
      "selected JAGS-4.3.1 version 4.3.1 link -ljags-4 assumed 4"
    )
    expect_selected(
      run_make(c("4.3.1", "4.3.2", "5.0.0")),
      "selected JAGS-4.3.2 version 4.3.2 link -ljags-4"
    )
    expect_stopped(
      run_make("5.0.0"),
      paste0(requirement, ", but only other JAGS versions were found: JAGS-5.0.0")
    )
    expect_stopped(
      run_make("4.2.0"),
      paste0(requirement, "; JAGS_VERSION/JAGS_ROOT reported 4.2.0")
    )
    expect_stopped(
      run_make("4.3.1", paste0("JAGS_ROOT=", roots[["5.0.0"]])),
      paste0(requirement, "; JAGS_VERSION/JAGS_ROOT reported 5.0.0")
    )
    expect_stopped(
      run_make("4.3.1", paste0("JAGS_ROOT=", roots[["custom5"]])),
      "declare JAGS_MAJOR 5"
    )
    expect_selected(
      run_make("4.3.1", paste0("JAGS_ROOT=", roots[["custom4"]])),
      "selected custom4 version  link -ljags-4"
    )
    expect_stopped(
      run_make("4.3.1", "JAGS_VERSION=5.0.0"),
      paste0(requirement, "; JAGS_VERSION/JAGS_ROOT reported 5.0.0")
    )
    expect_stopped(
      run_make("4.3.1", "JAGS_MAJOR_VERSION=5"),
      paste0(requirement, "; JAGS_MAJOR_VERSION is 5")
    )
    expect_selected(
      run_make("4.3.1", "JAGS_MAJOR_VERSION=4"),
      "selected JAGS-4.3.1 version 4.3.1 link -ljags-4"
    )
  }
})


test_that("the JAGS version header guard rejects JAGS 5 at compile time", {

  skip_on_cran()
  header_file <- .jags_build_repository_file("src", "jagsversions.h")
  skip_if_not(
    file.exists(header_file),
    "Repository build sources are not available in this installed-package context."
  )
  r_bin <- file.path(R.home("bin"), "R")
  cxx <- suppressWarnings(system2(r_bin, c("CMD", "config", "CXX"), stdout = TRUE, stderr = TRUE))
  cxx <- strsplit(trimws(cxx[length(cxx)]), "[[:space:]]+")[[1L]]
  skip_if(length(cxx) == 0L || !nzchar(Sys.which(cxx[[1L]])), "A C++ compiler is required.")

  work <- normalizePath(withr::local_tempdir(), winslash = "/", mustWork = TRUE)
  jags4 <- .jags_build_fake_tree(file.path(work, "jags4"), 4L)
  jags5 <- .jags_build_fake_tree(file.path(work, "jags5"), 5L)
  source_file <- file.path(work, "guard.cc")
  .jags_build_write_lf(
    c(paste0("#include \"", normalizePath(header_file, winslash = "/"), "\""), "int main() { return 0; }"),
    source_file
  )
  compile <- function(tree, forced = 0L, assumed = 0L) {

    .jags_build_run_tool(
      cxx[[1L]],
      c(cxx[-1L], "-fsyntax-only", paste0("-I", tree, "/include/JAGS"),
        paste0("-DJAGS_MAJOR_FORCED=", forced), paste0("-DJAGS_MAJOR_ASSUMED=", assumed),
        source_file),
      work,
      character()
    )
  }
  guard <- "RoBMA requires JAGS 4.x (>= 4.3.1, < 5.0.0); JAGS 5 is not supported."

  accepted <- compile(jags4, assumed = 4L)
  expect_identical(accepted$status, 0L, info = accepted$output)
  rejected <- compile(jags5)
  expect_false(identical(rejected$status, 0L), info = rejected$output)
  expect_match(rejected$output, guard, fixed = TRUE)
  forced <- compile(jags4, forced = 5L)
  expect_false(identical(forced$status, 0L), info = forced$output)
  expect_match(forced$output, guard, fixed = TRUE)
})


test_that("macOS CI builds against a cached JAGS 4 installer, never Homebrew's jags", {

  workflow_dir <- .jags_build_repository_file(".github", "workflows")
  skip_if_not(
    dir.exists(workflow_dir),
    "Repository workflow sources are not available in this installed-package context."
  )

  for (path in list.files(workflow_dir, pattern = "\\.ya?ml$", full.names = TRUE)) {
    workflow <- readLines(path, warn = FALSE, encoding = "UTF-8")
    brew_installs <- workflow[grepl("brew install", workflow, fixed = TRUE)]
    expect_false(any(grepl("\\bjags\\b", brew_installs)), info = basename(path))
    if (!any(grepl("name: Install JAGS (macOS)", workflow, fixed = TRUE))) {
      next
    }

    text <- paste(workflow, collapse = "\n")
    expect_match(text, "path: ~/jags-installer", fixed = TRUE, info = basename(path))
    expect_match(text, "~/jags-installer/JAGS-4\\.3\\.[0-9]+\\.pkg", info = basename(path))
    installer <- regexpr("sudo installer -pkg", text, fixed = TRUE)
    expect_true(installer > 0L, info = basename(path))
    # Homebrew links into /usr/local on Intel runners, so a Homebrew jags is
    # removed before the package installs JAGS 4 there, never afterwards.
    removal <- regexpr("brew uninstall --ignore-dependencies jags", text, fixed = TRUE)
    expect_true(removal > 0L && removal < installer, info = basename(path))
  }
})


test_that("DESCRIPTION declares the supported JAGS range", {

  description_file <- .jags_build_repository_file("DESCRIPTION")
  skip_if_not(
    file.exists(description_file),
    "Repository DESCRIPTION is not available in this installed-package context."
  )

  requirements <- unname(read.dcf(description_file, fields = "SystemRequirements")[1, 1])
  expect_match(requirements, "JAGS (>= 4.3.1, < 5.0.0)", fixed = TRUE)
})
