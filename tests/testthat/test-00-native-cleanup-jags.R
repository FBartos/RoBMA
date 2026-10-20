test_that("JAGS selected distributions refuse unsupported unobserved sampling", {

  skip_if_not_installed("rjags")
  .load_RoBMA_module()
  fixtures <- .native_selection_jags_fixtures()
  compile <- function(specification, data, observed) {
    used <- unique(c(if (observed) "y", specification$arguments))
    syntax <- paste0("model { y[1:2] ~ ", specification$distribution, "(",
      paste(specification$arguments, collapse = ", "), ") }")
    connection <- textConnection(syntax)
    on.exit(close(connection), add = TRUE)
    model <- rjags::jags.model(connection, data = data[used],
      inits = list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = 73071L),
      n.chains = 1L, n.adapt = 0L, quiet = TRUE)
    if (!observed) rjags::coda.samples(model, variable.names = "y",
      n.iter = 1L, progress.bar = "none")
    model
  }
  for (specification in fixtures$cases) {
    data <- utils::modifyList(fixtures$base, specification$extra)
    for (mode in 0:1) {
      data$mode <- mode
      observed <- compile(specification, data, TRUE)
      expect_s3_class(observed, "jags")
      result <- tryCatch(compile(specification, data, FALSE), error = identity)
      expect_s3_class(result, "error")
      message <- paste0("'", specification$distribution,
        "' supports observed nodes only; sampling unobserved selected nodes is unsupported.")
      expect_match(conditionMessage(result), message, fixed = TRUE)
      expect_false(grepl("Invalid parent values", conditionMessage(result), fixed = TRUE))
    }
  }
})

test_that("known-V JAGS sampling retains seeded PD draws and refuses non-PD covariance", {

  skip_if_not_installed("rjags")
  .load_RoBMA_module()
  compile <- function(covariance) {
    connection <- textConnection("model { y[1:2] ~ dknown_v_mnorm(mu, tau2, V_lower) }")
    on.exit(close(connection), add = TRUE)
    rjags::jags.model(connection, data = list(mu = c(.03, -.01), tau2 = c(.01, .04),
      V_lower = covariance[lower.tri(covariance, diag = TRUE)]),
      inits = list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = 73079L),
      n.chains = 1L, n.adapt = 0L, quiet = TRUE)
  }
  sample <- function() as.matrix(rjags::coda.samples(compile(matrix(c(.04, .015, .015, .09), 2L)),
    variable.names = "y", n.iter = 5L, progress.bar = "none"))
  draws <- sample()
  expect_true(all(is.finite(draws)))
  expect_identical(sample(), draws)
  result <- tryCatch(rjags::coda.samples(compile(matrix(c(1, 2, 2, 1), 2L)),
    variable.names = "y", n.iter = 1L, progress.bar = "none"), error = identity)
  expect_s3_class(result, "error")
  expect_match(conditionMessage(result), "Known-V sampling covariance is not positive definite.", fixed = TRUE)
  expect_false(grepl("Invalid parent values", conditionMessage(result), fixed = TRUE))
})

test_that("factor and zplot ingress mirror shared mode partition and order checks", {

  rule <- .selection_joint_factor_quadrature_rules()
  factor_args <- list(c(.1, -.1), matrix(c(0, 0), 1L), matrix(c(1, 1), 1L),
    matrix(c(.2, .2, .1, -.1), 1L), c(1, 1), matrix(c(1, .2), 1L),
    c(0, -Inf), c(Inf, 0), c(1L, 2L), 1L, TRUE, 0L,
    rule$nodes, rule$log_weights, as.double(rule$orders),
    rep(.5, 64L), 4L, 8L, 2L, .005, FALSE, 0L)
  factor <- function(args) do.call(.Call, c(list("RoBMA_selnorm_factor_step_loglik_batch"),
    args, list(PACKAGE = "RoBMA")))
  reference <- factor(factor_args)
  expect_true(all(is.finite(reference$log_density)))
  for (mode in c(NA_integer_, -1L, 2L, 3L, 999L)) {
    changed <- factor_args; changed[[12L]] <- mode
    expect_error(factor(changed), "Multivariate selection requires a normal or step kernel.", fixed = TRUE)
  }
  for (order in c(Inf, NaN, 2.5, .Machine$integer.max + 1, length(rule$nodes) + 1)) {
    changed <- factor_args; changed[[15L]][1L] <- order
    expect_error(factor(changed), "'quadrature_orders' must contain increasing positive integers matching the quadrature rules.", fixed = TRUE)
  }
  changed <- factor_args; changed[[7L]] <- c(-Inf, 0); changed[[8L]] <- c(0, Inf)
  expect_error(factor(changed), "Selection bounds do not define a descending step partition.", fixed = TRUE)
  expect_identical(factor(factor_args), reference)
  zplot_args <- list(c(.1, -.1), matrix(c(0, 0), 1L), matrix(c(1, .2, 1), 1L),
    c(1, 1), matrix(c(1, .2), 1L), c(0, -Inf), c(Inf, 0), c(1L, 2L),
    1L, TRUE, 0L, rep(.5, 32L), 4L, 2L, .005, c(-.1, .1), FALSE, NULL, 0L)
  zplot <- function(args) do.call(.Call, c(list("RoBMA_selnorm_mnorm_zplot_batch"),
    args, list(PACKAGE = "RoBMA")))
  reference <- zplot(zplot_args)
  expect_true(all(is.finite(reference$density)))
  changed <- zplot_args; changed[[11L]] <- 3L
  expect_error(zplot(changed), "Multivariate selection requires a normal or step kernel.", fixed = TRUE)
  changed <- zplot_args; changed[[6L]] <- c(-Inf, 0); changed[[7L]] <- c(0, Inf)
  expect_error(zplot(changed), "Selection bounds do not define a descending step partition.", fixed = TRUE)
  expect_identical(zplot(zplot_args), reference)
})

test_that("response RNG interrupts between rows and releases its consumed stream", {

  skip_on_cran()
  withr::local_seed(73081L)
  rows <- 50000L
  means <- matrix(0, rows, 2L)
  covariance <- array(rep(diag(2), each = rows), c(rows, 2L, 2L))
  weights <- matrix(1, rows, 1L)
  generate <- function() .Call("RoBMA_selnorm_mnorm_step_rng_batch", means,
    covariance, c(1, 1), weights, -Inf, Inf, 1L, 0L, list(1:2), 1L, 0L,
    PACKAGE = "RoBMA")
  initial <- .Random.seed
  complete <- generate()
  expect_identical(complete$failure_code, 0L)
  complete_seed <- .Random.seed
  assign(".Random.seed", initial, envir = .GlobalEnv)
  # The generator repeats inside the callback until the elapsed limit fires, so
  # the test does not depend on one call outlasting the limit (a fast machine
  # completes 50,000 rows in under 20 ms). The bound only stops a broken
  # interrupt check from running forever; reaching it fails the assertions.
  calls <- 0L
  interrupt <- function() {
    on.exit(setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE), add = TRUE)
    .Call("RoBMA_native_unwind_probe", "callback", function() {
      setTimeLimit(elapsed = .02, transient = TRUE)
      while (calls < 2000L) {
        generate()
        calls <<- calls + 1L
      }
    }, 1L, PACKAGE = "RoBMA")
  }
  result <- tryCatch(interrupt(), error = identity)
  expect_s3_class(result, "error")
  if (inherits(result, "error")) expect_match(conditionMessage(result), "elapsed time limit", fixed = TRUE)
  expect_false(identical(.Random.seed, initial))
  expect_false(identical(.Random.seed, complete_seed))
  expect_identical(.Call("RoBMA_native_unwind_probe", "state", NULL, 1L,
    PACKAGE = "RoBMA"), rep(FALSE, 3L))
  assign(".Random.seed", initial, envir = .GlobalEnv)
  expect_identical(generate(), complete)
  expect_identical(.Random.seed, complete_seed)
})
