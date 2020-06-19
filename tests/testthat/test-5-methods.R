context("(5) print and summary functions")
skip_on_cran()
# Make sure that the print window is streched as much as possible, errors might emerge because the prints messages get folded.

# test objects - assuming that the fit function worked properly
saved_fits    <- readRDS(file = "saved_fits.RDS")
saved_methods <- readRDS(file = "saved_methods.RDS")

test_that("Print function works", {
  for(i in 1:length(saved_fits)){
    expect_equal(capture.output(saved_fits[[i]]), saved_methods$print[[i]])
  }
})

test_that("Summary function works", {
  for(i in 1:length(saved_fits)){
    expect_equal(capture.output(summary(saved_fits[[i]], conditional = TRUE, include_theta = TRUE)), saved_methods$summary[[i]])
  }

  expect_equal(capture.output(summary(saved_fits[[1]], conditional = TRUE,  include_theta = TRUE, probs = c(.15), logBF = TRUE, digits_estimates = 5)), saved_methods$summary_add[[1]])
  expect_equal(capture.output(summary(saved_fits[[1]], conditional = FALSE, include_theta = FALSE, probs = c(.66, .33, .50), logBF = TRUE, BF01 = TRUE)), saved_methods$summary_add[[2]])
})

test_that("Diagnostics summary works", {
  for(i in 1:length(saved_fits)){
    expect_equal(capture.output(summary(saved_fits[[i]], type = "models", diagnostics = TRUE)), saved_methods$diagnostics[[i]])
  }
})

test_that("Individual summary works", {
  for(i in 1:length(saved_fits)){
    expect_equal(capture.output(summary(saved_fits[[i]], type = "individual")), saved_methods$individual[[i]])
  }
})


#### creating / updating the test settings ####
if(FALSE){
  saved_fits       <- readRDS(file = "tests/testthat/saved_fits.RDS")
  print_fits       <- list()
  summary_fits     <- list()
  diagnostics_fits <- list()
  individual_fits  <- list()
  for(i in 1:length(saved_fits)){
    print_fits[[i]]       <- capture.output(saved_fits[[i]])
    summary_fits[[i]]     <- capture.output(summary(saved_fits[[i]], conditional = TRUE, include_theta = TRUE))
    diagnostics_fits[[i]] <- capture.output(summary(saved_fits[[i]], type = "models", diagnostics = TRUE))
    individual_fits[[i]]  <- capture.output(summary(saved_fits[[i]], type = "individual"))
  }
  summary_add <- list(
    capture.output(summary(saved_fits[[1]], conditional = TRUE,  include_theta = TRUE, probs = c(.15), logBF = TRUE, digits_estimates = 5)),
    capture.output(summary(saved_fits[[1]], conditional = FALSE, include_theta = FALSE, probs = c(.66, .33, .50), logBF = TRUE, BF01 = TRUE))
  )
  saved_methods <- list(
    print       = print_fits,
    summary     = summary_fits,
    summary_add = summary_add,
    diagnostics = diagnostics_fits,
    individual  = individual_fits
  )
  saveRDS(saved_methods, file = "tests/testthat/saved_methods.RDS")
}

