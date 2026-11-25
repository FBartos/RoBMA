context("Persisted outcome fingerprint compatibility")

test_that("native outcome hash preserves the exact persisted byte recurrence", {

  reference <- function(bytes) {

    hash1 <- 5381
    hash2 <- 0
    for (byte in as.integer(bytes)) {
      hash1 <- (hash1 * 33 + byte) %% 2147483647
      hash2 <- (hash2 * 65599 + byte) %% 2147483629
    }
    c(hash1, hash2)
  }
  for (bytes in list(raw(), as.raw(0:255), charToRaw("v1:outcomes"),
                    rep(as.raw(255:0), 17L))) {
    actual <- .Call("RoBMA_outcome_hash_raw", bytes, PACKAGE = "RoBMA")
    expect_identical(actual, reference(bytes))
  }
  expect_error(.Call("RoBMA_outcome_hash_raw", 0, PACKAGE = "RoBMA"),
               "raw vector")
})
