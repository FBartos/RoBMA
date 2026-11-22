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
                                    second_study = "u[2]", nonlinear = FALSE, tags = TRUE) {

  location <- if (nonlinear) "sqrt(a)" else "a"
  scales <- if (scale) "s ~ dnorm(0, 1) T(0,)
                       r ~ dbeta(2, 2)" else ""
  between <- if (scale) "s * sqrt(r)" else ".3"
  precision <- if (scale) "1 / (1 + pow(s, 2) * (1 - r))" else "1"
  physical <- if (scale) "s, r" else ".5, .36"
  declaration <- if (tags) paste0("
    m1 <- robma_cluster_member(g, 1, 1)
    m2 <- robma_cluster_member(a, 2, 1)
    m3 <- robma_cluster_member(b, 3, 1)
    m4 <- robma_cluster_member(u[1], 4, 1)
    m5 <- robma_cluster_member(", second_study, ", 4, 2)\n",
    if (declare_scales) "m6 <- robma_cluster_member(s, 5, 1)
                         m7 <- robma_cluster_member(r, 6, 1)\n",
    "setup <- robma_cluster_setup(g, X, group, y, ones, ones, ones, -1, 1, 1, 0, 0, 1, 1, 2,
                                  ", physical, ", 2, 2, 0, 0, 2, 5, ", as.integer(scale), ")\n"
  ) else ""
  model <- paste0("model {
    g ~ dcat(c(1, 1))
    a ~ dnorm(0, 1)
    b[1:2] ~ dmnorm(zero, P)
    ", scales, "
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
  inits <- list(a = 1e-6, g = 2L, .RNG.name = "base::Mersenne-Twister", .RNG.seed = 7301L)
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

test_that("a failed clustered scalar refresh restores its entry value", {

  skip_if_not_installed("rjags")
  live <- .cluster_provider_model(nonlinear = TRUE)
  expect_setequal(rjags::list.samplers(live)[["RoBMA::JointCluster"]],
                   c("g", "a", "b", "u[1]", "u[2]"))
  stats::update(live, n.iter = 1L, progress.bar = "none")
  live$sync()
  entry <- live$state()[[1L]][["a"]]
  expect_true(entry > 0)
  expect_error(stats::update(live, n.iter = 1L, progress.bar = "none"), "Invalid parent values")
  live$sync()
  expect_identical(live$state()[[1L]][["a"]], entry)
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
