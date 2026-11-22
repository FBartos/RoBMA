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
  fit_cluster <- function(direction = "positive") {
    RoBMA(yi = y, vi = v, cluster = study, mods = ~ group, data = dat,
          measure = "SMD", effect_direction = direction, adapt = 60L, burnin = 50L,
          sample = 100L, chains = 2L, seed = 37L, autofit = FALSE, silent = TRUE)
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
  expect_identical(colnames(scaled[["fit"]][["mcmc"]][[1L]]),
                   colnames(joint[["fit"]][["mcmc"]][[1L]]))
  expect_true(all(is.finite(as.matrix(scaled[["fit"]][["mcmc"]]))))

  ordinary <- fit_cluster("negative")
  expect_length(owned(ordinary), 0L)
  expect_false(grepl("drobma_cluster_setup", ordinary[["fit"]][["model"]], fixed = TRUE))

  # A malformed native declaration with two markers naming the same study
  # slot must fall through as an entire block. No stochastic source is lost.
  original <- scaled[["fit"]]
  broken <- sub("drobma_cluster_member\\(gamma\\[1\\], [^,]+, [^)]+\\)",
                "drobma_cluster_member(gamma[1], 4, 2)", original[["model"]])
  expect_false(identical(broken, original[["model"]]))
  con <- textConnection(broken)
  on.exit(close(con), add = TRUE)
  live <- rjags::jags.model(con, data = original[["data"]], inits = original[["end_state"]],
                          n.chains = 2L, n.adapt = 0L, quiet = TRUE)
  ordinary_samplers <- rjags::list.samplers(live)
  expect_false("RoBMA::JointCluster" %in% names(ordinary_samplers))
  expect_true(all(c(group, scales) %in% unlist(ordinary_samplers, use.names = FALSE)))
})
