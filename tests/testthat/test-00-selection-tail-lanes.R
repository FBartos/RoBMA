.test_tail_lanes <- function(locations, offset = 0, inverse_sd = -1,
                             last = 0.1, step = 0.9) {
  .Call("RoBMA_selnorm_tail_lanes_check", as.double(locations), as.double(offset),
        as.double(inverse_sd), as.double(last), as.double(step), PACKAGE = "RoBMA")
}

.test_tail_lane_agreement <- function(result, count) {
  expect_identical(is.null(result$avx2), !result$avx2_available)
  expect_identical(is.null(result$avx512), !result$avx512_available)
  for (lane in c("avx2", "avx512", "automatic")) {
    if (!is.null(result[[lane]])) expect_true(result[[lane]]$guards)
  }
  if (result$avx2_available && result$avx512_available) {
    expect_identical(result$avx512, result$avx2)
  }
  expected <- if (result$avx512_available && count >= 8L) {
    result$avx512
  } else if (result$avx2_available) {
    result$avx2
  } else {
    NULL
  }
  if (!is.null(expected)) expect_identical(result$automatic, expected)
}

test_that("direct tail lanes preserve finite values and scalar accuracy", {
  locations <- seq(-37, 37, length.out = 20001)
  result <- .test_tail_lanes(locations)
  .test_tail_lane_agreement(result, length(locations))
  reference <- pnorm(locations, lower.tail = FALSE)
  expect_lt(max(abs(result$automatic$tail - reference) / reference), 1e-14)
  expect_lt(max(abs(result$scalar - reference) / reference), 1e-14)
  expect_lt(max(abs(result$automatic$tail - result$scalar) / reference), 1e-14)
  expect_true(result$automatic$valid)
})

test_that("tail lane blocks and remainders preserve boundary values", {
  cuts <- c(0.67448975, sqrt(32), 37.5193)
  # Exact cuts and close values on both sides, rotated across lane positions.
  scores <- c(0, -0, cuts, -cuts, cuts * (1 - 1e-15), cuts * (1 + 1e-15),
              -cuts * (1 - 1e-15), -cuts * (1 + 1e-15))
  for (count in c(0:17, 31L, 63L, 127L)) {
    for (phase in 0:7) {
      locations <- rep(scores, length.out = count + phase)[seq_len(count) + phase]
      result <- .test_tail_lanes(locations)
      .test_tail_lane_agreement(result, count)
      expect_length(result$automatic$tail, count)
      expect_length(result$automatic$normalizers, count)
    }
  }
})

test_that("mixed nonfinite lanes preserve classifications and validation", {
  scores <- c(NaN, Inf, -Inf, 0, -0, 1e300, -1e300, 0.2, 8, -8)
  for (phase in 0:7) {
    locations <- rep(scores, length.out = 17 + phase)[seq_len(17) + phase]
    for (offset in c(0, Inf, NaN)) {
      result <- .test_tail_lanes(locations, offset)
      .test_tail_lane_agreement(result, length(locations))
      expect_false(result$automatic$valid)
      if (result$avx2_available && result$avx512_available) {
        expect_identical(is.nan(result$avx2$tail), is.nan(result$avx512$tail))
        expect_identical(is.finite(result$avx2$normalizers),
                         is.finite(result$avx512$normalizers))
      }
    }
  }
})

test_that("tail lanes agree on cancellation and subnormal rejection", {
  settings <- list(c(1, -1), c(0, 1), c(.Machine$double.xmin / 2, 0))
  for (count in c(1L, 7L, 8L, 9L, 15L, 17L, 31L)) {
    for (setting in settings) {
      locations <- rep(c(-40, 40), length.out = count)
      # A single positive score must also reject the zero-base case.
      if (setting[1] == 0) locations <- abs(locations)
      result <- .test_tail_lanes(locations, last = setting[1], step = setting[2])
      .test_tail_lane_agreement(result, count)
      expect_false(result$automatic$valid)
    }
  }
})

test_that("private tail lane probe validates shape before evaluation", {
  probe <- function(x, offset = 0, scale = -1, last = 0.1, step = 0.9) {
    .Call("RoBMA_selnorm_tail_lanes_check", x, offset, scale, last, step,
          PACKAGE = "RoBMA")
  }
  expect_error(probe(1L), "numeric vector")
  expect_error(probe(1, offset = c(0, 1)), "one numeric")
  expect_error(probe(1, scale = numeric()), "one numeric")
  expect_error(probe(1, last = TRUE), "one numeric")
})
