# The native probe compares against the original scalar double expression and
# observes owner identity without retaining a process-wide test registry.
.test_double_rule_weights <- function(nodes, log_weights, shift = 0, evict = 0L) {
  .Call("RoBMA_selnorm_double_rule_weights_check", as.double(nodes),
        as.double(log_weights), as.double(shift), as.integer(evict), PACKAGE = "RoBMA")
}

test_that("double rule weights preserve arithmetic and exact owned keys", {
  nodes <- seq(-3, 3, length.out = 31)
  logs <- -nodes^2 / 2 - 2
  for (shift in c(0, -0, 0.25, -0.125)) {
    result <- .test_double_rule_weights(nodes, logs, shift)
    expect_identical(result$mismatches, 0L)
    expect_true(result$reused)
    expect_true(all(result$mutation_distinct))
    expect_true(result$retained)
    expect_equal(result$values, exp(logs - shift * nodes - 0.5 * shift^2),
                 tolerance = 1e-15)
  }
  # Same prefix, different count: the longer entry cannot answer the short rule.
  short <- .test_double_rule_weights(nodes[1:15], logs[1:15], 0.25)
  expect_length(short$values, 15)
  expect_identical(short$mismatches, 0L)
  expect_equal(short$values, exp(logs[1:15] - 0.25 * nodes[1:15] - 0.5 * 0.25^2),
               tolerance = 1e-15)
})

test_that("double rule owners survive entry and payload eviction", {
  for (count in c(31L, 12000L)) {
    nodes <- seq(-4, 4, length.out = count)
    logs <- -nodes^2 / 2 - 3
    # Four large entries exceed 1 MiB, while fewer than 32 keys are used.
    evictions <- if (count == 31L) 40L else 4L
    result <- .test_double_rule_weights(nodes, logs, 0.25, evict = evictions)
    expect_identical(result$mismatches, 0L)
    expect_true(result$reused)
    expect_true(result$retained)
    expect_true(result$evicted)
  }
})

test_that("oversized double rules compute without retention", {
  nodes <- seq(-4, 4, length.out = 70000)
  result <- .test_double_rule_weights(nodes, -nodes^2 / 2 - 3, 0.25)
  expect_identical(result$mismatches, 0L)
  expect_false(result$reused)
  expect_true(result$retained)
  expect_length(result$values, 70000)
})

test_that("extreme and nonfinite inputs preserve fresh double results", {
  nodes <- c(-2, 0, 2)
  for (shift in c(1e100, -1e100, Inf, -Inf, NaN)) {
    result <- .test_double_rule_weights(nodes, c(-1000, -1, 0), shift)
    expect_identical(result$mismatches, 0L)
    expect_true(result$retained)
    if (!is.finite(shift)) expect_false(result$reused)
  }
  result <- .test_double_rule_weights(c(-Inf, 0, Inf, NaN),
                                c(-Inf, Inf, NaN, -1), 0.25)
  expect_identical(result$mismatches, 0L)
  expect_true(result$retained)
})

test_that("private double rule probe validates its inputs", {
  probe <- function(nodes, logs, shift, evict) {
    .Call("RoBMA_selnorm_double_rule_weights_check", nodes, logs, shift, evict,
          PACKAGE = "RoBMA")
  }
  expect_error(probe(1, c(1, 2), 0, 0L), "equal nonempty")
  expect_error(probe(numeric(), numeric(), 0, 0L), "equal nonempty")
  expect_error(probe(1, 1, c(0, 1), 0L), "one numeric")
  expect_error(probe(1, 1, 0, NA_integer_), "integer between")
  expect_error(probe(1, 1, 0, -1L), "integer between")
  expect_error(probe(1, 1, 0, 65L), "integer between")
  expect_error(probe(1, 1, 0, 1.5), "integer between")
})
