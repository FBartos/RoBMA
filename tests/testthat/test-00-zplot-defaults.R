test_that("zplot density defaults follow the model class and preserve explicit budgets", {

  grid <- sort(unique(c(seq(-6, 6, .05), stats::qnorm(.975))))
  calls <- list()
  grids <- list()
  intervals <- list()
  summaries <- numeric()
  record <- function(max_samples, route){

    calls[[length(calls) + 1L]] <<- list(max_samples = max_samples, route = route)
    matrix(.1, 2L, length(grid))
  }
  testthat::local_mocked_bindings(
    .zplot_bins = function(priors, from, to, by, length.out, type){
      grids[[length(grids) + 1L]] <<- list(from = from, to = to, by = by,
        length.out = length.out, type = type)
      grid
    },
    .zplot_fun.brma = function(object, max_samples, ...){
      record(max_samples, "single")
    },
    .zplot_density_pair = function(object, max_samples, ...){
      density <- record(max_samples, "paired")
      list(fitted = density, extrapolated = density)
    },
    .zplot_density_data = function(z_sequence, z_density, probs){
      intervals[[length(intervals) + 1L]] <<- probs
      data.frame(x = z_sequence, y = .1, y_lCI = .05, y_uCI = .15)
    },
    hist.zplot_brma = function(...) invisible(NULL),
    as_zplot = function(object, max_samples, ...){
      summaries <<- c(summaries, max_samples)
      object
    },
    .package = "RoBMA"
  )
  testthat::local_mocked_bindings(
    lines = function(...) invisible(NULL),
    polygon = function(...) invisible(NULL),
    .package = "graphics"
  )

  model_classes <- list(c("brma.norm", "brma"),
    c("brma.mv", "brma.norm", "brma"),
    c("bselmodel.mv", "bselmodel", "brma.mv", "brma.norm", "brma"))
  default_budgets <- c(10000, 1000, 1000)
  for(model_index in seq_along(model_classes)){
    x <- structure(list(priors = list(), zplot = list(data = list(
      conditioning_depth = "marginal", integration_control = set_selection_likelihood_control()))),
      class = c("zplot_brma", model_classes[[model_index]]))
    expected <- default_budgets[[model_index]]
    calls <- list()
    plot.zplot_brma(x)
    expect_identical(calls[[1L]]$route, "single")
    expect_equal(calls[[1L]]$max_samples, expected)
    plot.zplot_brma(x, plot_extrapolation = TRUE)
    expect_identical(calls[[2L]]$route, "paired")
    expect_equal(calls[[2L]]$max_samples, expected)
    expect_s3_class(lines.zplot_brma(x, as_data = TRUE), "data.frame")
    expect_equal(calls[[3L]]$max_samples, expected)

    zplot(x)
    expect_equal(tail(calls, 1L)[[1L]]$max_samples, expected)
    zplot.brma(x)
    expect_equal(tail(calls, 1L)[[1L]]$max_samples, expected)
    expect_equal(tail(summaries, 1L), 10000)

    for(budget in c(137, Inf)){
      plot.zplot_brma(x, max_samples = budget)
      expect_equal(tail(calls, 1L)[[1L]]$max_samples, budget)
      lines.zplot_brma(x, max_samples = budget, as_data = TRUE)
      expect_equal(tail(calls, 1L)[[1L]]$max_samples, budget)
    }
  }
  expect_length(grid, 242L)
  expect_true(all(vapply(grids, identical, logical(1L),
    list(from = -6, to = 6, by = .05, length.out = NULL, type = "dens"))))
  expect_true(all(vapply(intervals, identical, logical(1L), c(.025, .975))))
})


test_that("the EDR summary budget remains independent of density plot defaults", {

  expect_identical(formals(as_zplot.brma)$max_samples, 10000)
  expect_identical(formals(zplot.brma)$summary_max_samples, 10000)
})


test_that("zplot displays only the fitted curve unless extrapolation is requested", {

  expect_identical(formals(plot.zplot_brma)$plot_fit, TRUE)
  expect_identical(formals(plot.zplot_brma)$plot_extrapolation, FALSE)
  expect_identical(formals(lines.zplot_brma)$extrapolate, FALSE)
})


test_that("the credible band is the one two separate quantile passes returned", {

  # `probs` reaches this helper already checked as two probabilities; one pass
  # per probability is what the band used to be built from, including the names
  # it gave the result, and `length.out` lets a caller ask for a grid of one
  # point or of none at all.
  probs <- c(.025, .975)
  reference <- function(z_sequence, z_density) {

    data.frame(
      x     = z_sequence,
      y     = colMeans(z_density),
      y_lCI = apply(z_density, 2, stats::quantile, probs = probs[[1L]]),
      y_uCI = apply(z_density, 2, stats::quantile, probs = probs[[2L]])
    )
  }

  set.seed(20260918)
  grids <- list(
    "many points"          = matrix(stats::rnorm(40L * 13L), nrow = 40L),
    "one grid point"       = matrix(stats::rnorm(40L), ncol = 1L),
    "no grid point"        = matrix(numeric(0), nrow = 40L, ncol = 0L),
    "one draw"             = matrix(stats::rnorm(7L), nrow = 1L),
    "two draws"            = matrix(stats::rnorm(14L), nrow = 2L),
    "one draw, one point"  = matrix(1.5, 1L, 1L),
    "two draws, one point" = matrix(c(1.5, 2.5), 2L, 1L),
    "tied draws"           = matrix(rep(c(1, 1, 1, 2, 2), each = 6L), nrow = 6L),
    "infinite entries"     = matrix(c(-Inf, 0, 1, Inf, 2, 3), nrow = 3L),
    "named grid"           = matrix(stats::rnorm(12L), nrow = 4L,
                                    dimnames = list(NULL, c("a", "b", "c")))
  )
  for (label in names(grids)) {
    z_density  <- grids[[label]]
    z_sequence <- seq_len(ncol(z_density))
    expect_identical(
      .zplot_density_data(z_sequence, z_density, probs),
      reference(z_sequence, z_density),
      info = label
    )
  }
})
