context("JAGS module build and runtime requirements")

test_that("the compiled RoBMA module matches the linked JAGS runtime", {

  flags <- rjags::pkg.compile.flags("ALL")
  major <- as.integer(flags$JAGS_VERSION[, 1L])
  compiled <- .C("getjagsversions", forced = integer(1), assumed = integer(1),
    detected = integer(1), used = integer(1), PACKAGE = "RoBMA")

  expect_true(major >= 5L)
  expect_identical(compiled$detected, major)
  expect_identical(compiled$used, major)
  expect_true("RoBMA" %in% rjags::list.modules())
})

test_that("both package builds get their JAGS flags from rjags", {

  root <- testthat::test_path("..", "..")
  skip_if_not(file.exists(file.path(root, "src", "Makevars")), "Package build sources are unavailable.")
  for (name in c("Makevars", "Makevars.win")) {
    source <- readLines(file.path(root, "src", name), warn = FALSE)
    for (flag in c("CPPFLAGS", "CXXFLAGS", "LIBS")) {
      expect_true(any(grepl(paste0("rjags::pkg.compile.flags('", flag, "')"), source, fixed = TRUE)),
        info = paste(name, flag))
    }
  }
})

test_that("DESCRIPTION declares the required JAGS and rjags versions", {

  description <- utils::packageDescription("RoBMA")
  expect_match(description$SystemRequirements, "JAGS (>= 5.0.0)", fixed = TRUE)
  expect_match(description$Imports, "rjags (>= 5-5)", fixed = TRUE)
})
