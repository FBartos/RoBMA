skip_on_cran()

test_that("clustered opt-in owns the gate, full Gaussian slab and all study effects", {

  skip_if_not_installed("rjags")
  options <- RoBMA.options()[c("sampler.joint", "selection.sampler")]
  on.exit(do.call(RoBMA.options, options), add = TRUE)
  RoBMA.options(selection.sampler = "slice")

  dat <- data.frame(
    y = c(.15, .20, -.15, .31, .17, -.21, .08, .26, -.19),
    v = seq(.04, .18, length.out = 9),
    study = rep(letters[1:3], c(2L, 3L, 4L)),
    group = factor(rep(letters[1:3], 3), levels = c("c", "a", "b"))
  )
  fit_cluster <- function(direction = "positive", ...) {
    RoBMA(yi = y, vi = v, cluster = study, mods = ~ group, data = dat,
          measure = "SMD", effect_direction = direction, adapt = 60L, burnin = 50L,
          sample = 100L, chains = 2L, seed = 37L, autofit = FALSE, silent = TRUE, ...)
  }
  owned <- function(fit) {
    a <- fit[["fit"]][["samplers"]]
    a[["Node"]][a[["Sampler"]] == "RoBMA::JointCluster"]
  }
  group <- c("mu_group_indicator", "mu_intercept_component_2",
             "mu_group_component_2", paste0("gamma[", 1:3, "]"))
  scales <- c("tau_component_2", "rho_component_1")

  RoBMA.options(sampler.joint = "clustered")
  joint <- fit_cluster()
  expect_setequal(owned(joint), group)

  RoBMA.options(sampler.joint = "clustered.scale")
  scaled <- fit_cluster()
  expect_setequal(owned(scaled), c(group, scales))
  request <- attr(scaled[["fit"]], "sampler_declarations", exact = TRUE)[["robma"]][["joint"]]
  expect_identical(request[["status"]], "active")
  expect_setequal(request[["nodes"]], c(group, scales))
  expect_identical(colnames(scaled[["fit"]][["mcmc"]][[1L]]),
                   colnames(joint[["fit"]][["mcmc"]][[1L]]))
  expect_true(all(is.finite(as.matrix(scaled[["fit"]][["mcmc"]]))))

  ordinary <- fit_cluster("negative")
  expect_length(owned(ordinary), 0L)
  expect_false(grepl("robma_cluster_setup", ordinary[["fit"]][["model"]], fixed = TRUE))

  one_bias <- BayesTools::prior_weightfunction("one-sided", .05, BayesTools::wf_fixed(c(1, .4)))
  singleton <- fit_cluster(prior_bias = one_bias, prior_bias_null = NULL)
  expect_true(singleton[["fit"]][["has_posterior"]])
  expect_length(owned(singleton), 0L)
  request <- attr(singleton[["fit"]], "sampler_declarations", exact = TRUE)[["robma"]][["joint"]]
  expect_identical(request[["status"]], "graph_declined")

  # A malformed native declaration with two markers naming the same study
  # slot must fall through as an entire block. No stochastic source is lost.
  original <- scaled[["fit"]]
  broken <- sub("robma_cluster_member\\(gamma\\[1\\], [^,]+, [^)]+\\)",
                "robma_cluster_member(gamma[1], 4, 2)", original[["model"]])
  expect_false(identical(broken, original[["model"]]))
  con <- textConnection(broken)
  on.exit(close(con), add = TRUE)
  live <- rjags::jags.model(con, data = original[["data"]], inits = original[["end_state"]],
                          n.chains = 2L, n.adapt = 0L, quiet = TRUE)
  ordinary_samplers <- rjags::list.samplers(live)
  expect_false("RoBMA::JointCluster" %in% names(ordinary_samplers))
  expect_true(all(c(group, scales) %in% unlist(ordinary_samplers, use.names = FALSE)))
})

.cluster_provider_model <- function(scale = FALSE, declare_scales = scale, observed = list(),
                                    second_study = "u[2]", location = "a", tags = TRUE,
                                    scale_aliases = "", heterogeneity_gate = "2",
                                    heterogeneity_label = 2L, scale_members = c("s", "r"),
                                    cadence = 5L, initial = list(), pet_setup = "0") {

  scales <- if (scale) "s ~ dnorm(0, 1) T(0,)
                       r ~ dbeta(2, 2)" else ""
  if (scale && !nzchar(scale_aliases)) scale_aliases <- "tau_value <- s
                                                      rho_value <- r"
  between <- if (scale) "tau_value * sqrt(rho_value)" else ".3"
  precision <- if (scale) "1 / (1 + pow(tau_value, 2) * (1 - rho_value))" else "1"
  physical <- if (scale) "tau_value, rho_value" else ".5, .36"
  declaration <- if (tags) paste0("
    m1 <- robma_cluster_member(g, 1, 1)
    m2 <- robma_cluster_member(a, 2, 1)
    m3 <- robma_cluster_member(b, 3, 1)
    m4 <- robma_cluster_member(u[1], 4, 1)
    m5 <- robma_cluster_member(", second_study, ", 4, 2)\n",
    if (declare_scales) paste0("m6 <- robma_cluster_member(", scale_members[1L], ", 5, 1)
                               m7 <- robma_cluster_member(", scale_members[2L], ", 6, 1)\n"),
    "setup <- robma_cluster_setup(g, X, group, y, ones, ones, ones, -1, 1, 1, 0, 0, 1, 1, 2,
                                  ", physical, ", ", heterogeneity_gate, ", ", heterogeneity_label,
                                  ", ", pet_setup, ", 0, 2, ", cadence, ", ", as.integer(scale), ")\n"
  ) else ""
  model <- paste0("model {
    g ~ dcat(c(1, 1))
    a ~ dnorm(0, 1)
    b[1:2] ~ dmnorm(zero, P)
    ", scales, "
    ", scale_aliases, "
    for (j in 1:2) { u[j] ~ dnorm(0, 1) }
    for (i in 1:4) {
      y[i] ~ dnorm(", location, " + equals(g, 2) * (b[1]*X[i] + b[2]*X[i+4]) + ",
                     between, " * u[group[i]], ", precision, ")
    }
    ", declaration, "
  }")
  data <- c(list(zero = c(0, 0), P = diag(2), group = c(1, 1, 2, 2),
                 y = c(.3, .25, -.13, .08), X = c(0, 1, 0, 1, 1, 0, 1, 0)), observed)
  if (tags) data[["ones"]] <- rep(1, 4)
  inits <- utils::modifyList(list(a = 1e-6, g = 2L, .RNG.name = "base::Mersenne-Twister",
                                .RNG.seed = 7301L), initial)
  if (!is.null(inits[[".RNG.state"]])) inits[[".RNG.seed"]] <- NULL
  con <- textConnection(model)
  on.exit(close(con))
  rjags::jags.model(con, data = data, inits = inits, n.chains = 1L, n.threads = 1L,
                    n.adapt = 0L, quiet = TRUE)
}

test_that("a clustered provider must supply distinct free roots for its requested scale mode", {

  skip_if_not_installed("rjags")
  roots <- c("g", "a", "b", "u[1]", "u[2]")
  own <- function(model) rjags::list.samplers(model)[["RoBMA::JointCluster"]]
  decline <- function(model, free) {
    samplers <- rjags::list.samplers(model)
    expect_null(samplers[["RoBMA::JointCluster"]])
    expect_setequal(unlist(samplers, use.names = FALSE), free)
  }
  expect_setequal(own(.cluster_provider_model()), roots)
  expect_setequal(own(.cluster_provider_model(scale = TRUE)), c(roots, "s", "r"))
  decline(.cluster_provider_model(second_study = "u[1]"), roots)
  decline(.cluster_provider_model(scale = TRUE, declare_scales = FALSE), c(roots, "s", "r"))
  decline(.cluster_provider_model(scale = TRUE, observed = list(s = .5)), c(roots, "r"))
  decline(.cluster_provider_model(scale = TRUE, observed = list(r = .36)), c(roots, "s"))
  decline(.cluster_provider_model(scale = TRUE, observed = list(s = .5, r = .36)), roots)
})

test_that("a failed clustered update restores every owned coordinate", {

  skip_if_not_installed("rjags")
  state <- function(model, internal = FALSE) {
    model$sync()
    model$state(internal = internal)[[1L]]
  }
  step <- function(model) {
    tryCatch({
      stats::update(model, n.iter = 1L, progress.bar = "none")
      NULL
    }, error = conditionMessage)
  }
  cases <- list(
    list(location = "sqrt(a)"),
    list(location = "a + sqrt(u[2])", initial = list(u = c(.1, .05))),
    list(scale = TRUE, location = "a + sqrt(r - .35)",
         initial = list(s = .5, r = .36, u = c(.1, .05)))
  )
  for (case in cases) {
    live <- do.call(.cluster_provider_model, case)
    owned <- c("g", "a", "b", "u[1]", "u[2]", if (isTRUE(case$scale)) c("s", "r"))
    expect_named(rjags::list.samplers(live), "RoBMA::JointCluster")
    expect_setequal(rjags::list.samplers(live)[[1L]], owned)

    # Freeze the scalar widths at their initial values. A replay can then
    # start from the current RNG state without reconstructing adaptation.
    rjags::adapt(live, n.iter = 0L, end.adaptation = TRUE)
    error <- NULL
    for (iteration in 1:5) {
      entry <- state(live)
      error <- step(live)
      if (!is.null(error)) break
    }
    expect_match(error, "Invalid parent values")
    expect_identical(state(live), entry)
    # Both models take the same refresh-only path before the grouped move
    # is due, even if the fresh model starts at a different cadence tick.
    if (iteration < 5L) {
      case$initial <- state(live, internal = TRUE)
      replay <- do.call(.cluster_provider_model, case)
      rjags::adapt(replay, n.iter = 0L, end.adaptation = TRUE)
      expect_identical(step(live), step(replay))
      expect_identical(state(live, internal = TRUE), state(replay, internal = TRUE))
    }
  }
})

test_that("a failed grouped proposal keeps the clustered move due on the next sweep", {

  skip_if_not_installed("rjags")
  model <- function(cadence, initial = list(u = c(.1, 1), .RNG.seed = 210L)) {
    live <- .cluster_provider_model(location = "a + sqrt(u[2] + .85)",
                                    cadence = cadence, initial = initial)
    rjags::adapt(live, n.iter = 0L, end.adaptation = TRUE)
    live
  }
  state <- function(model, internal = TRUE) {
    model$sync()
    model$state(internal = internal)[[1L]]
  }
  advance <- function(model, n = 1L) {
    stats::update(model, n.iter = n, progress.bar = "none")
  }

  live <- model(cadence = 3L)
  refresh_only <- model(cadence = 100L)
  advance(live, 2L)
  advance(refresh_only, 2L)
  expect_identical(state(live), state(refresh_only))
  entry <- state(live, internal = FALSE)
  expect_no_error(advance(refresh_only))
  expect_error(advance(live), "Invalid parent values")
  expect_identical(state(live, internal = FALSE), entry)

  # Rebuild from the consumed RNG stream with and without a grouped move due.
  initial <- state(live)
  due <- model(cadence = 1L, initial = initial)
  not_due <- model(cadence = 100L, initial = initial)
  advance(due)
  advance(not_due)
  expect_false(identical(state(due), state(not_due)))
  advance(live)
  expect_identical(state(live), state(due))
})

test_that("clustered scales resolve the declared physical aliases exactly", {

  skip_if_not_installed("rjags")
  roots <- c("g", "a", "b", "u[1]", "u[2]", "s", "r")
  aliases <- "h ~ dcat(c(1, 1))
              rho_gate ~ dcat(c(1))
              tau_value <- 0 * (h == 1) + s * (h == 2)
              rho_value <- r * (rho_gate == 1)"
  for (active in 1:2) {
    syntax <- if (active == 1L) {
      sub("0 * (h == 1) + s * (h == 2)", "s * (h == 1) + 0 * (h == 2)",
          aliases, fixed = TRUE)
    } else aliases
    live <- .cluster_provider_model(scale = TRUE, scale_aliases = syntax,
      heterogeneity_gate = "h", heterogeneity_label = active, initial = list(h = 3L - active))
    expect_setequal(rjags::list.samplers(live)[["RoBMA::JointCluster"]], roots)
  }

  common <- list(scale = TRUE, heterogeneity_gate = "h", scale_aliases = aliases)
  declined <- list(
    utils::modifyList(common, list(heterogeneity_label = 1L)),
    utils::modifyList(common, list(heterogeneity_gate = "other_h",
      scale_aliases = paste(aliases, "other_h ~ dcat(c(1, 1))", sep = "\n"),
      location = "a + (other_h - other_h)")),
    utils::modifyList(common, list(scale_aliases = sub("0 *", ".1 *", aliases, fixed = TRUE))),
    utils::modifyList(common, list(scale_aliases = sub("rho_gate == 1", "h == 1", aliases, fixed = TRUE))),
    utils::modifyList(common, list(scale_members = c("aux_s", "r"),
      scale_aliases = paste(aliases, "aux_s ~ dnorm(0, 1) T(0,)", sep = "\n"),
      location = "a + (aux_s - aux_s)", initial = list(aux_s = .5))),
    utils::modifyList(common, list(scale_members = c("s", "aux_r"),
      scale_aliases = paste(aliases, "aux_r ~ dbeta(2, 2)", sep = "\n"),
      location = "a + (aux_r - aux_r)", initial = list(aux_r = .36))),
    utils::modifyList(common, list(scale_members = c("r", "s")))
  )
  for (case in declined) {
    tagged <- do.call(.cluster_provider_model, case)
    stock <- do.call(.cluster_provider_model, c(case, list(tags = FALSE)))
    expect_identical(rjags::list.samplers(tagged), rjags::list.samplers(stock))
    draw <- function(model) {
      rjags::adapt(model, n.iter = 100L, end.adaptation = TRUE)
      rjags::coda.samples(model, c("g", "a", "b", "u", "s", "r", "h"), n.iter = 40L, thin = 2L)
    }
    expect_identical(draw(tagged), draw(stock))
    expect_identical(tagged$state(internal = TRUE), stock$state(internal = TRUE))
  }
})

test_that("declined clustered tags preserve stock samplers, draws and RNG state", {

  skip_if_not_installed("rjags")
  for (config in list(list(scale = FALSE, second_study = "u[1]"),
                      list(scale = TRUE, declare_scales = FALSE),
                      list(scale = TRUE, observed = list(s = .5)),
                      list(scale = TRUE, observed = list(s = .5, r = .36)))) {
    tagged <- do.call(.cluster_provider_model, config)
    stock <- do.call(.cluster_provider_model, c(config, list(tags = FALSE)))
    expect_identical(rjags::list.samplers(tagged), rjags::list.samplers(stock))
    expect_true("bugs::ConjugateNormal" %in% names(rjags::list.samplers(stock)))
    draw <- function(model) {
      rjags::adapt(model, n.iter = 100L, end.adaptation = TRUE)
      rjags::coda.samples(model, c("g", "a", "b", "u"), n.iter = 40L, thin = 2L)
    }
    expect_identical(draw(tagged), draw(stock))
    expect_identical(tagged$state(internal = TRUE), stock$state(internal = TRUE))
  }
})

test_that("clustered proposal arguments must be current in the sampling graph", {

  skip_if_not_installed("rjags")
  make <- function(location, pet_setup = "pet_alias", tags = TRUE) {
    .cluster_provider_model(location = location, pet_setup = pet_setup, tags = tags,
      scale_aliases = "p ~ dnorm(0, 1)
                       pet_alias <- p * 1",
      initial = list(p = .1))
  }
  draw <- function(model) {
    rjags::adapt(model, n.iter = 100L, end.adaptation = TRUE)
    rjags::coda.samples(model, c("g", "a", "b", "u", "p"), n.iter = 40L, thin = 2L)
  }

  # The stochastic slope is informative, but its private alias is not.
  tagged <- make("a + p")
  stock <- make("a + p", tags = FALSE)
  expect_identical(rjags::list.samplers(tagged), rjags::list.samplers(stock))
  expect_identical(draw(tagged), draw(stock))
  expect_identical(tagged$state(internal = TRUE), stock$state(internal = TRUE))

  # The same alias is safe once the likelihood itself uses it.
  live <- make("a + pet_alias")
  direct <- make("a + pet_alias", pet_setup = "p")
  owners <- rjags::list.samplers(live)
  expect_setequal(owners[["RoBMA::JointCluster"]], c("g", "a", "b", "u[1]", "u[2]"))
  expect_identical(owners, rjags::list.samplers(direct))
  expect_identical(draw(live), draw(direct))
  expect_identical(live$state(internal = TRUE), direct$state(internal = TRUE))
})
