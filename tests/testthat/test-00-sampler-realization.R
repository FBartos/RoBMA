test_that("a replayed sampler report resolves singleton aliases on each worker", {

  owner <- "RoBMA::GaussianJoint"
  first <- data.frame(Node = c("gate", "vector", "matrix", "exact[1]"),
                      Sampler = c(owner, "bugs::DMNorm", "bugs::DMNorm", "base::Finite"))
  second <- first
  second[["Node"]][[1L]] <- "gate[1]"
  normalized <- c("gate[1]", "gate", "gate[2]", "vector[1]", "matrix[1]", "exact[1]")
  requested <- list(joint = "gaussian", specialized = TRUE, nested_weights = FALSE,
                    probabilities = c(" gate [ 1 ] ", normalized[-1L]))
  state <- list(gate = 2, vector = c(1, 2), matrix = matrix(1, 1L, 1L))
  fit <- structure(list(samplers = first, end_state = list(state, state)),
                   sampler_declarations = list(joint = list(), probabilities = data.frame(node = normalized)))
  fitted <- .sampler_request_record(fit, requested)

  second[["Sampler"]][[1L]] <- "base::Finite"
  fitted[["samplers_by_worker"]] <- list(list(chains = 1L, samplers = first),
                                        list(chains = 2L, samplers = second))
  record <- attr(.sampler_request_record(fitted), "sampler_declarations")[["robma"]]
  expect_identical(record[["request"]], requested)
  expect_identical(record[["joint"]][["status"]], "mixed")
  expect_identical(record[["probabilities"]][["node"]], normalized)
  expect_identical(record[["probabilities"]][["requested_node"]], requested[["probabilities"]])
  expect_identical(record[["probabilities"]][["sampler"]],
                   c("mixed", "mixed", NA, NA, NA, "base::Finite"))
  expect_true(all(record[["probabilities"]][["status"]] == "unavailable"))
})


test_that("serialized updates refresh ownership with the original sampler request", {

  skip_on_cran()
  previous <- RoBMA.options()
  on.exit(do.call(RoBMA.options, previous), add = TRUE)
  owner <- "RoBMA::IndependentBiasJoint"
  .load_RoBMA_module()
  factories <- rjags::list.factories("sampler")
  previous_factory <- factories[["status"]][match(owner, factories[["factory"]])]
  on.exit(rjags::set.factory(owner, "sampler", previous_factory), add = TRUE)
  rjags::set.factory(owner, "sampler", TRUE)
  file <- tempfile(fileext = ".rds")
  on.exit(unlink(file), add = TRUE)
  data("Bem2011", package = "RoBMA", envir = environment())
  fit_small <- function() {
    RoBMA(yi = d, sei = se, data = Bem2011, measure = "SMD", chains = 2L,
          adapt = 500L, burnin = 50L, sample = 100L, seed = 708L,
          convergence_checks = set_convergence_checks(max_Rhat = NULL, min_ESS = NULL))
  }
  read_saved <- function(fit) {
    saveRDS(fit, file)
    readRDS(file)
  }
  declaration <- function(fit) attr(fit[["fit"]], "sampler_declarations", exact = TRUE)
  aliases <- c("bias_indicator", " mu_indicator [ 1 ] ", "tau_indicator[1]")
  RoBMA.options(sampler.joint = "auto", sampler.nested.weights = TRUE, selection.sampler = "slice",
                sampler.probabilities = aliases)
  aliased <- read_saved(fit_small())
  record <- declaration(aliased)[["robma"]][["probabilities"]]
  expect_identical(record[["node"]], c("bias_indicator", "mu_indicator[1]", "tau_indicator[1]"))
  expect_identical(record[["requested_node"]], aliases)
  expect_identical(record[["sampler"]], c(owner, owner, "BayesTools::Finite"))
  expect_identical(record[["status"]], c("unavailable", "unavailable", "active"))
  expect_null(aliased[["fit"]][["sampler_probabilities"]][["mu_indicator[1]"]])

  RoBMA.options(sampler.probabilities = "tau_indicator[1]")
  original <- read_saved(fit_small())
  original_plan <- declaration(original)
  expect_identical(original_plan[["robma"]][["joint"]][["status"]], "active")
  RoBMA.options(sampler.joint = "off", sampler.nested.weights = FALSE, selection.sampler = "jags.default",
                sampler.probabilities = character())
  rjags::set.factory(owner, "sampler", FALSE)

  expect_warning(unchanged <- update(aliased, sample_extend = 20L, adapt = 500L),
                 "incompatible probability trace")
  expect_identical(unchanged[["fit"]][["mcmc"]], aliased[["fit"]][["mcmc"]])
  expect_identical(declaration(unchanged)[["robma"]], declaration(aliased)[["robma"]])

  disabled <- update(original, sample_extend = 20L, adapt = 500L)
  actual <- disabled[["fit"]][["samplers"]]
  expect_identical(actual[["Sampler"]][match("bias_indicator", actual[["Node"]])], "base::Finite")
  expect_identical(declaration(disabled)[["robma"]][["joint"]][["status"]], "native_declined")
  expect_length(declaration(disabled)[["robma"]][["joint"]][["nodes"]], 0L)
  expect_true(declaration(disabled)[["robma"]][["nested_weights"]][["realized"]])
  expect_identical(declaration(disabled)[["robma"]][["request"]], original_plan[["robma"]][["request"]])
  expect_identical(declaration(disabled)[["joint"]], original_plan[["joint"]])
  expect_equal(nrow(disabled[["fit"]][["mcmc"]][[1L]]), 120L)
  expect_identical(as.matrix(disabled[["fit"]][["mcmc"]][[1L]])[seq_len(100L), , drop = FALSE],
                   as.matrix(original[["fit"]][["mcmc"]][[1L]]))

  rjags::set.factory(owner, "sampler", TRUE)
  restored <- update(read_saved(disabled), sample_extend = 20L, adapt = 500L)
  expect_identical(declaration(restored)[["robma"]][["joint"]], original_plan[["robma"]][["joint"]])
  expect_identical(declaration(restored)[["robma"]][["probabilities"]], original_plan[["robma"]][["probabilities"]])
  expect_equal(nrow(restored[["fit"]][["sampler_probabilities"]][["tau_indicator[1]"]][[1L]]), 140L)

  ordinary <- read_saved(fit_small())
  expect_null(declaration(ordinary)[["robma"]])
  RoBMA.options(sampler.joint = "auto", sampler.probabilities = "tau_indicator")
  quiet <- update(ordinary, sample_extend = 20L, adapt = 500L)
  expect_null(declaration(quiet)[["robma"]])
})
