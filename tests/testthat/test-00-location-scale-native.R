context("Heterogeneity scale move of the location block sampler")
source(testthat::test_path("common-functions.R"))
skip_on_cran()

# The location block sampler can own the heterogeneity sources of the
# random-effect term whose latents it samples (src/samplers/LocationScale.cc):
# their ordinary samplers run embedded, then one scale-group move of
# (a, w, Z) that leaves every child density unchanged. The hand models and the
# independent references are in helper-location-scale.R.


test_that("the scale owner claims the declared sources and embeds their ordinary samplers", {

  skip_if_not_installed("rjags")
  setup <- .scale_allocation_setup(shapes = c(1.5, 0.8))
  model <- .scale_allocation_model(setup)
  inits <- list(I = 1, a = 0.4, eta = c(0.7, 1.1))
  owned <- .scale_compile(model, inits)
  expect_identical(owned$decisions$status, "owned")
  expect_identical(owned$decisions$layout, "allocation")
  expect_identical(owned$decisions$check, "passed")
  members <- c("a", "eta[1]", "eta[2]", "g", "slab", paste0("Z[", 1:4, "]"))
  expect_true(all(owned$samplers[members] == "RoBMA::LocationBlock"))
  expect_false(owned$samplers[["I"]] == "RoBMA::LocationBlock")

  # declined, every node keeps its sampler; the embedded samplers are those
  declined <- .scale_compile(model, inits, decline = 1L)
  expect_identical(declined$decisions$status, "declined")
  expect_identical(declined$decisions$reason, "test switch")
  expect_true(all(declined$samplers[c("g", "slab", paste0("Z[", 1:4, "]"))] == "RoBMA::LocationBlock"))
  expect_false(any(declined$samplers[c("a", "eta[1]", "eta[2]")] == "RoBMA::LocationBlock"))
  expect_identical(sort(strsplit(owned$decisions$embedded, ",", fixed = TRUE)[[1L]]),
                   sort(unname(declined$samplers[c("a", "eta[1]", "eta[2]")])))

  # the clustered layout with gated sources
  clustered <- .scale_compile(.scale_clustered_model(setup),
                              list(tau_indicator = 2, tau_component_2 = 0.3, rho_component_1 = 0.4))
  expect_identical(clustered$decisions$status, "owned")
  expect_identical(clustered$decisions$layout, "clustered")
  expect_true(all(clustered$samplers[c("tau_component_2", "rho_component_1", paste0("gamma[", 1:4, "]"))] ==
                    "RoBMA::LocationBlock"))
  expect_false(any(clustered$samplers[c("tau_indicator", "rho_indicator")] == "RoBMA::LocationBlock"))
})


test_that("frozen states: the scale move keeps the children, the gamma total and its group law, and its density is Liu and Sabatti's", {

  skip_if_not_installed("rjags")
  setup <- .scale_allocation_setup()
  z <- c(0.8, -1.3, 0.4, 1.7)
  cases <- list(
    list(model = .scale_allocation_model(.scale_allocation_setup(shapes = c(1.5, 0.8))), s = 2L,
         inits = list(I = 1), x = c(z, 0.45, 0.7, 1.1), active = TRUE),
    list(model = .scale_allocation_model(.scale_allocation_setup(shapes = c(0.6, 2.5, 1.2))), s = 2L,
         inits = list(I = 1), x = c(z, 0.45, 0.3, 1.9, 0.8), active = TRUE),
    list(model = .scale_allocation_model(.scale_allocation_setup(shapes = c(0.6, 2.5, 1.2))), s = 2L,
         inits = list(I = 0), x = c(z, 0.45, 0.3, 1.9, 0.8), active = FALSE),
    list(model = .scale_clustered_model(setup), s = 1L,
         inits = list(tau_indicator = 2), x = c(z, 0.35, 0.3), active = TRUE),
    list(model = .scale_clustered_model(setup), s = 1L,
         inits = list(tau_indicator = 1), x = c(z, 0.35, 0.3), active = FALSE)
  )
  close <- function(a, b, tolerance) max(abs(a - b) / (1 + abs(b)))
  for (case in cases) {
    model <- case$model
    info <- paste(model$layout, "L =", model$L, "active =", case$active)
    compiled <- .scale_compile(model, case$inits)
    expect_identical(compiled$decisions$status, "owned", info = info)
    owner <- compiled$owner
    .scale_hook("set", owner, 0L, case$x)
    start <- .scale_hook("state", owner, 0L)
    expect_identical(start$active, case$active, info = info)
    expect_equal(start$prior[[2L]], start$prior[[1L]], tolerance = 1e-12, info = info)
    for (u in c(-0.9, 0.35)) {
      expect_true(.scale_hook("transform", owner, 0L, u), info = info)
      moved <- .scale_hook("state", owner, 0L)
      # every child density unchanged, sigma_s scaled by c, the R map
      expect_lt(close(moved$children, start$children), 1e-10)
      expect_lt(abs(moved$sigma - exp(u) * start$sigma), 1e-12 * (1 + start$sigma))
      expect_lt(close(moved$values, .scale_transform(case$x, u, model$layout, model$L, model$K, case$s)), 1e-13)
      if (model$layout == "allocation") {
        total <- sum(case$x[-seq_len(model$K + 1L)])
        expect_lt(abs(sum(moved$values[-seq_len(model$K + 1L)]) - total), 1e-14 * total)
      }
      # T_c T_r = T_(c r)
      expect_true(.scale_hook("transform", owner, 0L, 0.5), info = info)
      composed <- .scale_hook("state", owner, 0L)
      expect_lt(close(composed$values, .scale_transform(case$x, u + 0.5, model$layout, model$L, model$K, case$s)),
                1e-13)
      .scale_hook("set", owner, 0L, case$x)
    }
    # log f(u) against pi(T_u x) |det J_u(x)| from the node densities and a
    # numerical Jacobian of the source map
    u <- c(-1.2, -0.4, 0.3, 0.9)
    hook <- .scale_hook("density", owner, 0L, c(0, u))
    direct <- .scale_direct_density(case$x, u, model, case$s)
    expect_lt(max(abs((hook[-1L] - hook[1L]) - direct) / pmax(1, abs(direct))), 1e-8)
  }

  # each mutant breaks the identity it is meant to break
  model <- cases[[1L]]$model
  x <- cases[[1L]]$x
  u <- c(-0.8, 0.6)
  direct <- .scale_direct_density(x, u, model, 2L)
  for (mutant in 1:4) {
    compiled <- .scale_compile(model, list(I = 1), mutant = mutant)
    owner <- compiled$owner
    .scale_hook("set", owner, 0L, x)
    start <- .scale_hook("state", owner, 0L)
    if (mutant == 3L) {
      # a transformation that changes the children: the move's check undoes it
      # with the classed error naming the child
      expect_error(.scale_hook("transform", owner, 0L, u[[1L]]), "RoBMA_location_scale_invariance", fixed = TRUE)
      expect_identical(.scale_hook("state", owner, 0L)$values, start$values)
    } else {
      hook <- .scale_hook("density", owner, 0L, c(0, u))
      expect_gt(max(abs((hook[-1L] - hook[1L]) - direct)), 1e-3, label = paste("mutant", mutant))
    }
  }
})


test_that("the production slice leaves a broad bounded target invariant and a fixed cap per side does not", {

  # one update of the move's slice (width 0.5, a budget of 40 expansions) from
  # exact draws of a uniform target on (0, 100), which no interval reaches
  set.seed(5113L)
  n <- 200000L
  edges <- function(points) {
    p <- c(mean(points < 2), mean(points > 98), mean(points > 49 & points < 51))
    (p - 0.02) / sqrt(0.02 * 0.98 / n)
  }
  valid <- edges(.scale_hook("slice", 0L, 0L, c(0, 100, stats::runif(n, 0, 100))))
  expect_lt(max(abs(valid)), 4)
  fixed <- edges(.scale_hook("slice", 1L, 0L, c(0, 100, stats::runif(n, 0, 100))))
  expect_lt(min(fixed[1:2]), -6)
})


test_that("a truncated total SD and extreme shares keep the move and the embedded updates valid", {

  skip_if_not_installed("rjags")
  truncated <- .scale_allocation_model(.scale_allocation_setup(a_upper = 0.6))
  compiled <- .scale_compile(truncated, list(I = 1))
  owner <- compiled$owner
  x <- c(0.8, -1.3, 0.4, 1.7, 0.59, 0.4, 1.6)
  .scale_hook("set", owner, 0L, x)
  # a' = a sqrt(1 - w_s + c^2 w_s) beyond the bound: no density, no move
  expect_identical(.scale_hook("density", owner, 0L, 0.5), -Inf)
  expect_false(.scale_hook("transform", owner, 0L, 0.5))
  expect_identical(.scale_hook("state", owner, 0L)$values, x)
  rjags::adapt(compiled$jags, 0L, end.adaptation = TRUE)
  draws <- rjags::coda.samples(compiled$jags, "a", n.iter = 300L, progress.bar = "none")
  expect_lt(max(as.matrix(draws)), 0.6)

  # w_s at 1e-8 and 1 - 1e-8: outside the move's set ([2^-20, 1 - 2^-20]),
  # where the move is the identity; the embedded updates run
  model <- .scale_allocation_model(.scale_allocation_setup())
  for (eta in list(c(1, 1e-8), c(1e-8, 1))) {
    compiled <- .scale_compile(model, list(I = 1))
    x <- c(0.8, -1.3, 0.4, 1.7, 0.45, eta)
    .scale_hook("set", compiled$owner, 0L, x)
    expect_identical(.scale_hook("density", compiled$owner, 0L, c(0, 0.3)), c(-Inf, -Inf))
    expect_false(.scale_hook("transform", compiled$owner, 0L, 0.3))
    expect_identical(.scale_hook("state", compiled$owner, 0L)$values, x)
    rjags::adapt(compiled$jags, 0L, end.adaptation = TRUE)
    expect_silent(stats::update(compiled$jags, n.iter = 50L, progress.bar = "none"))
  }
  clustered <- .scale_clustered_model(.scale_allocation_setup())
  for (rho in c(1e-8, 1 - 1e-8)) {
    compiled <- .scale_compile(clustered, list(tau_indicator = 2))
    x <- c(0.8, -1.3, 0.4, 1.7, 0.35, rho)
    .scale_hook("set", compiled$owner, 0L, x)
    expect_identical(.scale_hook("density", compiled$owner, 0L, 0), -Inf)
    expect_false(.scale_hook("transform", compiled$owner, 0L, -0.3))
    expect_identical(.scale_hook("state", compiled$owner, 0L)$values, x)
    rjags::adapt(compiled$jags, 0L, end.adaptation = TRUE)
    expect_silent(stats::update(compiled$jags, n.iter = 50L, progress.bar = "none"))
  }
})


test_that("the move acts on weights within [2^-20, 1 - 2^-20] with a stable chart and is the identity beyond", {

  skip_if_not_installed("rjags")
  # The reviewers' frozen states. (1) Gamma sources (1e-16, 1): the sampled
  # weight rounds to one and 1 + w_s expm1(2u) cancelled to zero; the state is
  # outside the move's set. Inside it, at a weight one step from the bound,
  # d = (eta_other + c^2 eta_s) / T is formed from the sources and log f
  # agrees with the direct density in both directions.
  model <- .scale_allocation_model(.scale_allocation_setup())
  compiled <- .scale_compile(model, list(I = 1, a = 0.4, eta = c(0.7, 1.1), g = 2, slab = 0.3))
  owner <- compiled$owner
  outside <- c(rep(0, 4L), 0.45, 1e-16, 1)
  .scale_hook("set", owner, 0L, outside)
  # no density at the state itself and no move from it; the orbit point at
  # u = -19.5 lies in the set and has a finite density (d = 1.1e-16, which the
  # cancelling form gave as zero)
  density <- .scale_hook("density", owner, 0L, c(0, -19.5))
  expect_identical(density[1L], -Inf)
  expect_true(is.finite(density[2L]))
  expect_false(.scale_hook("transform", owner, 0L, -19.5))
  expect_identical(.scale_hook("state", owner, 0L)$values, outside)
  edge <- c(0.8, -1.3, 0.4, 1.7, 0.45, 2^-19, 1 - 2^-19)
  .scale_hook("set", owner, 0L, edge)
  u <- -3
  hook <- .scale_hook("density", owner, 0L, c(0, u))
  expect_lt(abs((hook[2L] - hook[1L]) - .scale_direct_density(edge, u, model, 2L)), 1e-8)
  expect_true(.scale_hook("transform", owner, 0L, u))
  moved <- .scale_hook("state", owner, 0L)$values
  expect_lt(max(abs(moved - .scale_transform(edge, u, "allocation", 2L, 4L, 2L)) / (1 + abs(moved))), 1e-13)
  back <- .scale_hook("density", owner, 0L, c(0, -u))
  expect_lt(abs((back[2L] - back[1L]) + (hook[2L] - hook[1L])), 1e-8)
  expect_true(.scale_hook("transform", owner, 0L, -u))
  expect_lt(max(abs(.scale_hook("state", owner, 0L)$values - edge) / (1 + abs(edge))), 1e-13)
  # a transformed point beyond the bound has no density on the orbit
  expect_identical(.scale_hook("density", owner, 0L, 3), -Inf)

  # (2) The clustered share one ULP below one with SEs near the data domain's
  # lower limit: the rebuilt within-cluster variance changed the likelihood.
  # The state is outside the move's set; one step inside it, the move keeps
  # every child to rounding.
  setup <- .scale_allocation_setup(a_sd = 300)
  z <- c(0.8, -1.3, 0.4, 1.7)
  setup$y <- 100 * z[setup$cluster] + 1e-5 * setup$y
  setup$variance[] <- 1.1e-12
  clustered <- .scale_compile(.scale_clustered_model(setup, rho_shapes = c(5e15, 1)),
                              list(tau_indicator = 2, tau_component_2 = 100, rho_component_1 = 0.4, g = 1L, slab = 0))
  owner <- clustered$owner
  x <- c(z, 100, 1 - .Machine$double.eps / 2)
  .scale_hook("set", owner, 0L, x)
  expect_identical(.scale_hook("density", owner, 0L, log(0.8)), -Inf)
  expect_false(.scale_hook("transform", owner, 0L, log(0.8)))
  expect_identical(.scale_hook("state", owner, 0L)$values, x)
  inside <- c(z, 100, 1 - 2^-19)
  .scale_hook("set", owner, 0L, inside)
  before <- .scale_hook("state", owner, 0L)
  expect_true(.scale_hook("transform", owner, 0L, log(1.2)))
  after <- .scale_hook("state", owner, 0L)
  expect_lt(abs(sum(after$children) - sum(before$children)), 1e-6)
  expect_identical(.scale_hook("density", owner, 0L, log(2)), -Inf)
})


test_that("ownership does not depend on initial values, and a move that changes the likelihood is an error", {

  skip_if_not_installed("rjags")
  # (1) A truncated total SD: from a start near its bound the factory's probe
  # (now at interior values) decides as from a start far from it.
  truncated <- .scale_allocation_model(.scale_allocation_setup(a_upper = 0.6))
  for (a in c(0.59, 0.05)) {
    compiled <- .scale_compile(truncated, list(I = 1, a = a))
    expect_identical(compiled$decisions$status, "owned", info = paste("a =", a))
    expect_identical(compiled$decisions$check, "passed", info = paste("a =", a))
  }

  # (2) A hand-tagged graph whose integrated variance also depends on the
  # sampled SD through an undeclared indicator h: the same decision from
  # either initial value of h. With h = 1 the move would change the
  # likelihood: the classed error names the child and the state is restored.
  model <- .scale_allocation_model(.scale_allocation_setup())
  model$syntax <- sub("  integrated <- ([^\n]+)", "  h ~ dbern(0.5)\n  integrated <- \\1 + h * 0.5 * pow(sigma, 2)",
                      model$syntax)
  starts <- list(I = 1, a = 0.4, eta = c(0.7, 1.1), g = 2, slab = 0.3)
  zero <- .scale_compile(model, c(starts, list(h = 0)))
  one  <- .scale_compile(model, c(starts, list(h = 1)))
  expect_identical(zero$decisions$status, "owned")
  expect_identical(one$decisions[c("status", "check")], zero$decisions[c("status", "check")])
  owner <- one$owner
  x <- c(0.8, -1.3, 0.4, 1.7, 0.45, 0.7, 1.1)
  .scale_hook("set", owner, 0L, x)
  before <- .scale_hook("state", owner, 0L)
  expect_error(.scale_hook("transform", owner, 0L, log(1.3)), "RoBMA_location_scale_invariance", fixed = TRUE)
  after <- .scale_hook("state", owner, 0L)
  expect_identical(after$values, x)
  expect_identical(after$children, before$children)
  # in the chain the error arises at the first move with h = 1 and the
  # inclusion on, and the move's entry is restored (the graph refused nothing)
  rjags::adapt(one$jags, 0L, end.adaptation = TRUE)
  error <- tryCatch({
    stats::update(one$jags, n.iter = 200L, progress.bar = "none")
    NULL
  }, error = conditionMessage)
  expect_match(error, "RoBMA_location_scale_invariance", fixed = TRUE)
  expect_match(error, "y[", fixed = TRUE)
})


test_that("a dependence on the sampled SD that only the open gate shows is an error at the first move", {

  skip_if_not_installed("rjags")
  # The factory's probe checks the charts at values the declaration fixes; the
  # likelihood's invariance is checked at every move, so the extra dependence
  # raises the classed error once the gate is open, from either initial value
  # of the gate, with the move's entry restored.
  model <- .scale_allocation_model(.scale_allocation_setup(), extra = TRUE)
  for (included in 0:1) {
    compiled <- .scale_compile(model, list(I = included, a = 0.4))
    expect_identical(compiled$decisions$status, "owned")
    expect_identical(compiled$decisions$check, "passed")
    .scale_hook("set", compiled$owner, 0L, c(0.8, -1.3, 0.4, 1.7, 0.45, 0.7, 1.1))
    before <- .scale_hook("state", compiled$owner, 0L)
    if (before$active) {
      expect_error(.scale_hook("transform", compiled$owner, 0L, 0.2), "RoBMA_location_scale_invariance", fixed = TRUE)
      expect_identical(.scale_hook("state", compiled$owner, 0L)$values, before$values)
    }
    rjags::adapt(compiled$jags, 0L, end.adaptation = TRUE)
    error <- tryCatch({
      stats::update(compiled$jags, n.iter = 300L, progress.bar = "none")
      NULL
    }, error = conditionMessage)
    expect_match(error, "RoBMA_location_scale_invariance", fixed = TRUE)
  }
})


test_that("after adaptation the move's invariance is also checked in log-likelihood units", {

  skip_if_not_installed("rjags")
  # A relative tolerance on parent values has no unit-free meaning. Two
  # hand-tagged graphs whose means also depend on the sampled SD, in small
  # units: (1) the allocation fixture in units of 1e-5 (row SEs about 1.5e-6,
  # which RoBMA admits) with + 0.5 sigma, where a move with u = 0.2 changes
  # the rows' log density by -0.63; (2) row SEs of 1e-6 around 1e-5 with
  # + 0.1 sigma, where c = 1.8 changes it by +1.40. Every mean moves by less
  # than the parent check's 1e-6, so only the check in log-likelihood units
  # sees them: the classed error at the first active move after adaptation,
  # with the move's entry restored.
  first <- function() {
    k <- 1e-5
    setup <- .scale_allocation_setup()
    setup$y <- setup$y * k
    setup$variance <- setup$variance * k^2
    setup$a_sd <- setup$a_sd * k
    setup$prior_precision <- 1 / k^2
    model <- .scale_allocation_model(setup)
    model$syntax <- sub("mu[i] <- b + sigma * Z[cluster[i]]", "mu[i] <- b + sigma * Z[cluster[i]] + 0.5 * sigma",
                        model$syntax, fixed = TRUE)
    list(compiled = .scale_compile(model, list(I = 1L, a = 0.4 * k, g = 1L)),
         x = c(0.8, -1.3, 0.4, 1.7, 0.45 * k, 0.7, 1.1), u = 0.2)
  }
  second <- function() {
    setup <- .scale_allocation_setup()
    setup$variance[] <- 1e-12
    setup$y[] <- 1e-5
    model <- .scale_allocation_model(setup)
    model$syntax <- sub("sigma * Z[cluster[i]]", "sigma * Z[cluster[i]] + 0.1 * sigma", model$syntax, fixed = TRUE)
    list(compiled = .scale_compile(model, list(I = 1, a = 1e-5, eta = c(1, 1), g = 2, slab = 0, Z = rep(0, 4))),
         x = c(rep(0, 4), 1e-5, 1, 1), u = log(1.8))
  }
  for (case in list(first = first, second = second)) {
    graph <- case()
    owner <- graph$compiled$owner
    expect_identical(graph$compiled$decisions$status, "owned")
    rjags::adapt(graph$compiled$jags, 0L, end.adaptation = TRUE)
    .scale_hook("set", owner, 0L, graph$x)
    before <- .scale_hook("state", owner, 0L)
    expect_true(before$active)
    error <- tryCatch({
      .scale_hook("transform", owner, 0L, graph$u)
      ""
    }, error = conditionMessage)
    expect_match(error, paste0("Location scale move [RoBMA_location_scale_invariance]: the log density of this ",
                               "likelihood node changed"), fixed = TRUE)
    after <- .scale_hook("state", owner, 0L)
    expect_identical(after$values, graph$x)
    expect_identical(after$children, before$children)
  }

  # An ordinary clustered graph inside the move's chart (tau = 0.001, the share
  # 1 - 2^-20, SEs 1e-6, outcomes 0.08 to 0.22): the rows are tens of
  # thousands of their SDs from their means, and the move's check refuses the
  # state with the precision guard's error before writing
  setup <- .scale_allocation_setup(a_sd = .5)
  setup$variance[] <- 1e-12
  setup$y <- seq(.08, .22, length.out = length(setup$y))
  clustered <- .scale_compile(.scale_clustered_model(setup), list(g = 1, slab = 0, gamma = c(.8, -1.3, .4, 1.7),
    tau_indicator = 2, tau_component_2 = 1e-3, rho_component_1 = 1 - 2^-20))
  rjags::adapt(clustered$jags, 0L, end.adaptation = TRUE)
  values <- c(.8, -1.3, .4, 1.7, 1e-3, 1 - 2^-20)
  .scale_hook("set", clustered$owner, 0L, values)
  before <- .scale_hook("state", clustered$owner, 0L)
  error <- tryCatch({
    .scale_hook("transform", clustered$owner, 0L, log(.8))
    ""
  }, error = conditionMessage)
  expect_match(error, "Location scale move [RoBMA_likelihood_not_representable]", fixed = TRUE)
  expect_identical(.scale_hook("state", clustered$owner, 0L)$values, before$values)

  # The same ordinary clustered graph with its outcomes times 0.05 is inside
  # the chart at M = 1.94e8: moves that stretch it by 0.65 to 0.9 change the
  # rows' summed log density by rounding (up to 0.006, below the tolerance
  # 1e-9 M = 0.19). Each is evaluated, counted and kept, without an error.
  setup$y <- setup$y * .05
  clustered <- .scale_compile(.scale_clustered_model(setup), list(g = 1, slab = 0, gamma = c(.8, -1.3, .4, 1.7),
    tau_indicator = 2, tau_component_2 = 1e-3, rho_component_1 = 1 - 2^-20))
  rjags::adapt(clustered$jags, 0L, end.adaptation = TRUE)
  for (factor in c(.65, .7, .75, .8, .85, .9)) {
    .scale_hook("set", clustered$owner, 0L, values)
    expect_true(.scale_hook("transform", clustered$owner, 0L, log(factor)))
  }
  counters <- RoBMA:::.location_scale_report_read()$counters
  expect_identical(counters[["likelihood_checks"]], 6)
  expect_gt(counters[["likelihood_changed"]], 0)
  expect_gt(counters[["likelihood_change_max"]], 1e-4)
  expect_lt(counters[["likelihood_change_max"]], 0.19)

  # Every move after adaptation is evaluated, not the first and every 64th:
  # in the second hand-tagged graph a move that stretches by 1 + 1e-10
  # changes the rows by 1.8e-10 and is kept, the next (c = 1.8) changes them
  # by +1.40 and fails with the entry of that move restored.
  graph <- second()
  rjags::adapt(graph$compiled$jags, 0L, end.adaptation = TRUE)
  .scale_hook("set", graph$compiled$owner, 0L, graph$x)
  expect_true(.scale_hook("transform", graph$compiled$owner, 0L, log(1 + 1e-10)))
  kept <- .scale_hook("state", graph$compiled$owner, 0L)
  expect_false(identical(kept$values, graph$x))
  error <- tryCatch({
    .scale_hook("transform", graph$compiled$owner, 0L, log(1.8))
    ""
  }, error = conditionMessage)
  expect_match(error, "[RoBMA_location_scale_invariance]: the log density of this likelihood node changed", fixed = TRUE)
  expect_identical(.scale_hook("state", graph$compiled$owner, 0L)$values, kept$values)

  # the valid fixtures pass the check at every one of 20,000 moves after
  # adaptation; their changes are rounding
  allocation <- .scale_allocation_model(.scale_allocation_setup(include = 1))
  clustered  <- .scale_clustered_model(.scale_allocation_setup())
  for (case in list(list(model = allocation, inits = list(I = 1), monitor = "a", draws = 20500L),
                    list(model = clustered, inits = list(tau_indicator = 2), monitor = "tau_component_2",
                         draws = 30000L))) {
    run <- .scale_run(.scale_compile(case$model, case$inits), case$monitor, draws = case$draws, burnin = 0L)
    counters <- run$report$counters
    # every move that was not exhausted wrote its point and was checked
    expect_identical(counters[["likelihood_checks"]], counters[["moves"]] - counters[["exhausted"]],
                     info = case$monitor)
    expect_gte(counters[["likelihood_checks"]], 20000, label = case$monitor)
    expect_lt(counters[["likelihood_change_max"]], 1e-10, label = case$monitor)
  }
})


test_that("the scale move's likelihood correction samples a graph whose change stays within the tolerance", {

  skip_if_not_installed("rjags")
  # A hand-tagged allocation graph whose likelihood the orbit changes by less
  # than the move's tolerance: a child v of sigma alone, mean 1e-7 sigma^2 /
  # (0.09 + sigma^2) (below the parent check's 1e-6) at precision 1.7e14, so
  # its log density spans (-0.85, 0]; and a child far from its constant mean
  # (log density -9e8, unchanged by the move) that sets M, and so the
  # tolerance 1e-9 M, near 0.9. Every move's change stays below it, and the
  # move corrects it: the posterior mean of sigma agrees with the stock
  # samplers' (the scale declined) within Monte Carlo error (|z| < 4, batch
  # means), and without the correction (mutant 5) it does not.
  setup <- .scale_allocation_setup(include = 1)
  model <- .scale_allocation_model(setup)
  model$syntax <- sub("  for (k in 1:K) { Z[k] ~ dnorm(0, 1) }",
                      paste("  for (k in 1:K) { Z[k] ~ dnorm(0, 1) }", "  far <- 0 * a + far_mean",
                            "  y_far ~ dnorm(far, 1)", "  v_mean <- 1e-7 * pow(sigma, 2) / (0.09 + pow(sigma, 2))",
                            "  v ~ dnorm(v_mean, v_precision)", sep = "\n"),
                      model$syntax, fixed = TRUE)
  model$data <- c(model$data, list(far_mean = sqrt(2 * 9e8), y_far = 0, v = 0, v_precision = 1.7e14))
  run <- function(mutant = 0L, decline = 0L) {
    compiled <- .scale_compile(model, list(I = 1), mutant = mutant, decline = decline, seed = 9157L)
    out <- .scale_run(compiled, "sigma", draws = 100000L, burnin = 500L, adapt = 500L)
    sigma <- out$draws[, "sigma"]
    means <- colMeans(matrix(sigma, length(sigma) / 50L, 50L))
    list(status = compiled$decisions$status, mean = mean(sigma), mcse = stats::sd(means) / sqrt(50),
         counters = out$report$counters)
  }
  z <- function(a, b) (a$mean - b$mean) / sqrt(a$mcse^2 + b$mcse^2)
  stock     <- run(decline = 1L)
  corrected <- run()
  expect_identical(stock$status, "declined")
  expect_identical(corrected$status, "owned")
  expect_lt(abs(z(corrected, stock)), 4)
  # the moves changed the likelihood, within the tolerance, and the
  # correction rejected some of them
  expect_gt(corrected$counters[["likelihood_changed"]], 0.99 * corrected$counters[["likelihood_checks"]])
  expect_gt(corrected$counters[["likelihood_rejected"]], 0)
  expect_gt(corrected$counters[["likelihood_change_max"]], 0.1)
  expect_lt(corrected$counters[["likelihood_change_max"]], 0.9)
  uncorrected <- run(mutant = 5L)
  expect_identical(uncorrected$counters[["likelihood_rejected"]], 0)
  expect_gt(abs(z(uncorrected, stock)), 4)
})


test_that("declined and undeclared scales draw exactly as before", {

  skip_if_not_installed("rjags")
  setup <- .scale_allocation_setup()
  monitor <- c("a", "eta", "I", "g", "slab", "Z")
  draw <- function(model, decline) {
    .scale_run(.scale_compile(model, list(I = 1), decline = decline), monitor, draws = 300L, burnin = 100L,
               adapt = 100L)$draws
  }
  undeclared <- draw(.scale_allocation_model(setup, tag = FALSE), 0L)
  declined   <- draw(.scale_allocation_model(setup), 1L)
  owned      <- draw(.scale_allocation_model(setup), 0L)
  expect_identical(declined, undeclared)
  expect_false(identical(owned, undeclared))
})


test_that("the heterogeneity block keeps its corrected record under the scale owner of a RoBMA.mv fit", {

  skip_if_not_installed("rjags")
  options <- RoBMA.options()[c("sampler.joint")]
  on.exit(do.call(RoBMA.options, options), add = TRUE)
  RoBMA.options(sampler.joint = "auto")
  dat <- data.frame(
    study = rep(1:8, each = 3L),
    esid  = seq_len(24L),
    yi    = 0.25 + rep(c(0.21, -0.18, 0.32, -0.27, 0.05, 0.24, -0.15, 0.08), each = 3L) +
      rep(c(0.06, -0.04, 0.02), 8L),
    vi    = rep(c(0.010, 0.014, 0.012), 8L)
  )
  fit <- function(decline) {
    previous <- .scale_config(decline = decline)
    on.exit(.scale_config(decline = previous[["decline"]]), add = TRUE)
    suppressWarnings(RoBMA.mv(yi = yi, V = vi, random = ~ 1 | study / esid, data = dat, measure = "SMD",
                              adapt = 60L, burnin = 50L, sample = 100L, chains = 1L, seed = 3L,
                              autofit = FALSE, silent = TRUE))
  }
  owned <- fit(0L)
  declined <- fit(1L)
  scale <- attr(owned[["fit"]], "sampler_declarations")[["robma"]][["joint"]][["scale"]]
  expect_identical(scale$decisions$status, "owned")
  expect_identical(scale$decisions$check, "passed")
  expect_gt(scale$counters[["moves"]], 0)
  table <- declined[["fit"]][["samplers"]]
  sources <- c("mu__xRE_ALLOCx_heterogeneity__allocation_sd",
               paste0("prior_par_eta_mu__xRE_ALLOCx_component_1_split__weight[", 1:2, "]"))
  expect_identical(unique(table[["Sampler"]][table[["Node"]] %in% sources]), scale$decisions$embedded)
  owned_table <- owned[["fit"]][["samplers"]]
  expect_true(all(owned_table[["Sampler"]][owned_table[["Node"]] %in% sources] == "RoBMA::LocationBlock"))
  expect_identical(attr(declined[["fit"]], "sampler_declarations")[["robma"]][["joint"]][["scale"]]$decisions$reason,
                   "test switch")
  # the same corrected samplers, with their record
  record <- function(fit) attr(fit[["fit"]], "selection_sampler", exact = TRUE)
  expect_identical(record(owned)[["samplers"]], record(declined)[["samplers"]])
  expect_identical(record(owned)[["fallback_samplers"]], record(declined)[["fallback_samplers"]])
})


test_that("a fit's scale record starts empty although an earlier report was not read", {

  skip_if_not_installed("rjags")
  # an unrelated live scale owner leaves decisions and moves in the process store
  unrelated <- .scale_compile(.scale_clustered_model(.scale_allocation_setup()),
                              list(tau_indicator = 2, tau_component_2 = 0.3, rho_component_1 = 0.4))
  rjags::adapt(unrelated$jags, 0L, end.adaptation = TRUE)
  stats::update(unrelated$jags, n.iter = 20L, progress.bar = "none")
  pending <- .Call("RoBMA_location_scale_report", FALSE, PACKAGE = "RoBMA")
  expect_gt(pending$counters[["moves"]], 0)
  dat <- data.frame(study = rep(1:8, each = 3L), esid = seq_len(24L),
                    yi = 0.25 + rep(c(0.21, -0.18, 0.32, -0.27, 0.05, 0.24, -0.15, 0.08), each = 3L) +
                      rep(c(0.06, -0.04, 0.02), 8L),
                    vi = rep(c(0.010, 0.014, 0.012), 8L))
  fit <- function() RoBMA.mv(yi = yi, V = vi, random = ~ 1 | study / esid, data = dat, measure = "SMD",
                             adapt = 60L, burnin = 50L, sample = 100L, chains = 2L, seed = 3L,
                             autofit = FALSE, silent = TRUE)
  record <- function(fit) attr(fit$fit, "sampler_declarations")[["robma"]][["joint"]][["scale"]]
  first  <- fit()
  second <- fit()
  expect_false("clustered" %in% record(first)$decisions$layout)
  expect_identical(record(first)$decisions, record(second)$decisions)
  expect_identical(record(first)$counters[c("owners", "moves", "identity", "outside")],
                   record(second)$counters[c("owners", "moves", "identity", "outside")])
  expect_identical(first$fit$mcmc, second$fit$mcmc)
})
