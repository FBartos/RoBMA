context("Session cache of hypothesis plans")

# The cache of an object is found by the hash of the object's whole content. The
# objects here are plain lists: the cache does not read what it is a cache of.


test_that("an object has one cache, found by its content", {

  .hypothesis_plan_cache_clear()
  withr::defer(.hypothesis_plan_cache_clear())
  object  <- list(fit = list(draws = c(1, 2, 3)), data = "a", priors = "p")
  copy    <- list(fit = list(draws = c(1, 2, 3)), data = "a", priors = "p")
  other   <- list(fit = list(draws = c(1, 2, 4)), data = "a", priors = "p")
  changed <- list(fit = list(draws = c(1, 2, 3)), data = "a", priors = "q")

  cache <- .hypothesis_plan_cache(object)
  expect_true(is.environment(cache))
  expect_identical(.hypothesis_plan_cache(object), cache)
  expect_identical(.hypothesis_plan_cache(copy), cache)
  # Other draws, or other priors, are other objects.
  expect_false(identical(.hypothesis_plan_cache(other), cache))
  expect_false(identical(.hypothesis_plan_cache(changed), cache))
  # Without an object the cache is that of the call.
  expect_false(identical(.hypothesis_plan_cache(), .hypothesis_plan_cache()))
  expect_false(isTRUE(.hypothesis_plan_cache()[[".bounded"]]))
  expect_true(isTRUE(cache[[".bounded"]]))

  # An object that cannot be hashed has the cache of its call only.
  testthat::local_mocked_bindings(
    hash = function(x, ...) stop("cannot hash"),
    .package = "rlang"
  )
  expect_false(identical(.hypothesis_plan_cache(object),
                         .hypothesis_plan_cache(object)))
})


test_that("a plan is computed once per object, statement and arguments", {

  .hypothesis_plan_cache_clear()
  withr::defer(.hypothesis_plan_cache_clear())
  builds <- 0L
  build  <- function() {
    builds <<- builds + 1L
    list(value = builds)
  }
  cache <- .hypothesis_plan_cache(list(a = 1))
  expect_identical(.hypothesis_plan_cached(cache, "one", build), list(value = 1L))
  expect_identical(.hypothesis_plan_cached(cache, "one", build), list(value = 1L))
  expect_identical(.hypothesis_plan_cached(cache, "two", build), list(value = 2L))
  expect_identical(builds, 2L)
  # The cache of another object has its own entries, and a failing computation
  # keeps nothing.
  other <- .hypothesis_plan_cache(list(a = 2))
  expect_identical(.hypothesis_plan_cached(other, "one", build), list(value = 3L))
  expect_error(.hypothesis_plan_cached(cache, "failing", function() {
    stop("no plan")
  }), "no plan")
  expect_false(exists("failing", envir = cache, inherits = FALSE))
  # Without a cache the plan is computed each time.
  expect_identical(.hypothesis_plan_cached(NULL, "one", build), list(value = 4L))
  expect_identical(.hypothesis_plan_cached(NULL, "one", build), list(value = 5L))
})


test_that("the cache of an object keeps a bounded number of entries and no entry above its size", {

  .hypothesis_plan_cache_clear()
  withr::defer(.hypothesis_plan_cache_clear())
  cache <- .hypothesis_plan_cache(list(a = 1))
  limit <- .hypothesis_plan_cache_limit()
  for (i in seq_len(limit + 5L)) {
    .hypothesis_plan_cached(cache, paste0("entry-", i), function() list(i))
  }
  kept <- setdiff(ls(cache, all.names = TRUE), c(".bounded", ".order", ".sizes"))
  expect_length(kept, limit)
  # The oldest entries are the ones dropped.
  expect_false(exists("entry-1", envir = cache, inherits = FALSE))
  expect_true(exists(paste0("entry-", limit + 5L), envir = cache,
                     inherits = FALSE))
  expect_identical(cache[[".order"]], paste0("entry-", 6:(limit + 5L)))

  # Entries above the size bound are computed again, not kept.
  testthat::local_mocked_bindings(
    .hypothesis_plan_cache_max_bytes = function() 4000,
    .package = "RoBMA"
  )
  cache  <- .hypothesis_plan_cache(list(a = 2))
  builds <- 0L
  large  <- function() {
    builds <<- builds + 1L
    numeric(1000L)
  }
  .hypothesis_plan_cached(cache, "large", large)
  .hypothesis_plan_cached(cache, "large", large)
  expect_identical(builds, 2L)
  expect_false(exists("large", envir = cache, inherits = FALSE))
  # The entries together stay under the size bound, the oldest out.
  small <- function() numeric(200L)
  for (i in 1:6) {
    .hypothesis_plan_cached(cache, paste0("small-", i), small)
  }
  expect_lte(sum(cache[[".sizes"]]), 4000)
  expect_true(exists("small-6", envir = cache, inherits = FALSE))
  expect_false(exists("small-1", envir = cache, inherits = FALSE))
})


test_that("the package keeps the caches of a few objects, the oldest out", {

  .hypothesis_plan_cache_clear()
  withr::defer(.hypothesis_plan_cache_clear())
  limit   <- .hypothesis_plan_registry_limit()
  objects <- lapply(seq_len(limit + 2L), function(i) list(draws = i))
  caches  <- lapply(objects, .hypothesis_plan_cache)
  registry <- setdiff(ls(.hypothesis_plan_registry, all.names = TRUE), ".order")
  expect_length(registry, limit)
  # The last objects' caches are kept and found again; the first is new again.
  expect_identical(.hypothesis_plan_cache(objects[[limit + 2L]]),
                   caches[[limit + 2L]])
  expect_false(identical(.hypothesis_plan_cache(objects[[1L]]), caches[[1L]]))

  .hypothesis_plan_cache_clear()
  expect_length(ls(.hypothesis_plan_registry, all.names = TRUE), 0L)
})
