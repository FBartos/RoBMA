context("Native C++ unwinding across R conditions")

.native_unwind_probe <- function(action, callback = NULL, threads = 1L) {

  .Call("RoBMA_native_unwind_probe", action, callback, threads, PACKAGE = "RoBMA")
}

test_that("R and C++ failures release all native probe workspaces", {

  for (mode in c("cpp", "bad_alloc", "allocation", "r_error")) {
    before <- .native_unwind_probe("counts")
    expect_error(.native_unwind_probe(mode))
    after <- .native_unwind_probe("counts")
    expect_identical(after[[1L]], 0L)
    expect_identical(after[[2L]] - before[[2L]], 2L)
  }
  before <- .native_unwind_probe("counts")
  expect_identical(.native_unwind_probe("normal"), 42L)
  after <- .native_unwind_probe("counts")
  expect_identical(after[[1L]], 0L)
  expect_identical(after[[2L]] - before[[2L]], 2L)
})

test_that("malformed R objects and protection-stack failures also unwind", {

  before <- .native_unwind_probe("counts")
  expect_error(.native_unwind_probe("malformed_data", new.env()), "environment")
  expect_error(.native_unwind_probe("protect_stack"), "protection stack overflow")
  after <- .native_unwind_probe("counts")
  expect_identical(after[[1L]], 0L)
  expect_identical(after[[2L]] - before[[2L]], 4L)
})

test_that("native boundaries preserve R errors and interrupt conditions", {

  for (kind in c("error", "interrupt")) {
    condition <- structure(list(message = paste("probe", kind), call = NULL,
                                marker = 731L),
                           class = c("native_probe_condition", kind, "condition"))
    before <- .native_unwind_probe("counts")
    caught <- tryCatch(.native_unwind_probe("callback", function() stop(condition)),
                       native_probe_condition = identity)
    expect_identical(caught, condition)
    after <- .native_unwind_probe("counts")
    expect_identical(after[[1L]], 0L)
    expect_identical(after[[2L]] - before[[2L]], 2L)
  }
  # Repeat the actual R allocation-error route; native state must not
  # accumulate across R tryCatch recovery or a later successful .Call.
  for (iteration in seq_len(30L)) {
    tryCatch(.native_unwind_probe("allocation"), error = function(e) NULL)
  }
  gc()
  expect_identical(.native_unwind_probe("counts")[[1L]], 0L)
  expect_identical(.native_unwind_probe("normal"), 42L)
})

test_that("nested native calls restore the outer continuation boundary", {

  condition <- structure(list(message = "nested condition", call = NULL),
                         class = c("nested_probe_error", "error", "condition"))
  before <- .native_unwind_probe("counts")
  caught <- tryCatch(.native_unwind_probe("callback", function() {
    .native_unwind_probe("normal")
    stop(condition)
  }), nested_probe_error = identity)
  expect_identical(caught, condition)
  after <- .native_unwind_probe("counts")
  expect_identical(after[[1L]], 0L)
  expect_identical(after[[2L]] - before[[2L]], 4L)
})

test_that("native interrupt checks unwind without leaving live workspaces", {

  on.exit(setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE), add = TRUE)
  before <- .native_unwind_probe("counts")
  expect_error(.native_unwind_probe("time_limit", function() {
    setTimeLimit(elapsed = .02, transient = TRUE)
  }), "elapsed time limit")
  setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
  after <- .native_unwind_probe("counts")
  expect_identical(after[[1L]], 0L)
  expect_identical(after[[2L]] - before[[2L]], 2L)
})

test_that("worker failures return to the main boundary without R API calls", {

  for (threads in 1:2) {
    expect_error(.native_unwind_probe("worker", threads = threads),
                 "native probe worker failure")
    expect_identical(.native_unwind_probe("counts")[[1L]], 0L)
  }
  if (isTRUE(.native_unwind_probe("openmp"))) {
    expect_error(.native_unwind_probe("worker_api", threads = 2L),
                 "outside the native main-thread boundary")
    expect_identical(.native_unwind_probe("counts")[[1L]], 0L)
  }
  expect_identical(.native_unwind_probe("altrep", as.double(1:16), 2L),
                   as.double(2:17))
})

test_that("interrupted native simulations save consumed R random draws", {

  withr::local_seed(987L)
  expected <- runif(2L)[[2L]]
  set.seed(987L)
  expect_error(.native_unwind_probe("rng_callback", function() stop("rng probe")),
               "rng probe")
  expect_identical(runif(1L), expected)
  expect_identical(.native_unwind_probe("counts")[[1L]], 0L)
})

test_that("nested native simulations share their parent's active RNG stream", {

  withr::local_seed(741L)
  expected <- runif(4L)
  set.seed(741L)
  .native_unwind_probe("rng_callback", function() {
    .native_unwind_probe("rng_success")
    .native_unwind_probe("rng_success")
  })
  expect_identical(runif(1L), expected[[4L]])

  set.seed(741L)
  expect_error(.native_unwind_probe("rng_callback", function() {
    .native_unwind_probe("rng_callback", function() stop("nested RNG failure"))
  }), "nested RNG failure")
  expect_identical(runif(1L), expected[[3L]])
  expect_identical(.native_unwind_probe("counts")[[1L]], 0L)
})

test_that("paired RNG calls borrow and release only their own native scope", {

  spec <- .test_step_spec(0, 1)
  draw <- function() .selnorm_kernel_rng_matrix(matrix(0), matrix(1), 1,
    matrix(1, 1L, spec$n_bins), spec, kernel_mode = 0L)[[1L]]
  withr::local_seed(971L)
  expected <- runif(4L)
  set.seed(971L)
  nested <- NULL
  .native_unwind_probe("rng_callback", function() nested <<- draw())
  expect_equal(nested, qnorm(expected[[3L]]), tolerance = 1e-14)
  expect_identical(runif(1L), expected[[4L]])

  # A child is the owner when the parent has not acquired RNG state.
  set.seed(971L)
  .native_unwind_probe("callback", function() nested <<- draw())
  expect_equal(nested, qnorm(expected[[2L]]), tolerance = 1e-14)
  expect_identical(runif(1L), expected[[3L]])
})

test_that("R warnings promoted to errors unwind native state", {

  withr::local_options(warn = 2)
  before <- .native_unwind_probe("counts")
  expect_error(.native_unwind_probe("warning"), "native probe warning")
  after <- .native_unwind_probe("counts")
  expect_identical(after[[1L]], 0L)
  expect_identical(after[[2L]] - before[[2L]], 2L)
})

test_that("successful results remain rooted across R API calls and RNG saving", {

  gctorture(TRUE)
  on.exit(gctorture(FALSE), add = TRUE)
  normal <- .native_unwind_probe("normal")
  rng <- .native_unwind_probe("rng_success")
  counts <- .native_unwind_probe("counts")
  gctorture(FALSE)
  expect_identical(normal, 42L)
  expect_identical(rng, 42L)
  expect_identical(counts[[1L]], 0L)
})
