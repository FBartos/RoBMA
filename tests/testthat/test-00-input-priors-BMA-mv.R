skip_on_cran()

test_that("component scale overrides retain real global SD and inclusion consumers", {

  dat <- data.frame(yi = c(0.1, 0.2, 0.3, 0.4), x = c(-1, 0, 1, 2),
    study = factor(c("a", "a", "b", "b")), obs = factor(seq_len(4)))
  normal <- BayesTools::prior("normal", list(0, 1))
  global <- BayesTools::prior("gamma", list(shape = 2, rate = 2), prior_weights = 3)
  zero <- BayesTools::prior("point", list(0), prior_weights = 1)
  args <- list(yi = quote(yi), V = diag(0.04, 4),
    random = list(study_effect = ~ 1 | study, observation_effect = ~ 1 | obs),
    scale = list(study = ~ x, obs = ~ 1),
    prior_scale = list(study = list(intercept = BayesTools::prior("point", list(0.3)), x = normal),
      obs = list(intercept = BayesTools::prior("point", list(0.4)))),
    data = dat, measure = "GEN", prior_effect = normal, only_priors = TRUE, silent = TRUE)
  check_binding <- function(object, probability = NULL) {
    expect_s3_class(object, "brma.mv")
    design <- .fitted_formula_design(object, "mu", required = TRUE)
    for (term in design[["random_effects"]]) {
      expect_s3_class(term[["sd_binding"]][["source"]], "random_sd_source")
      expect_match(term[["sd_binding"]][["source"]][["name"]], "^tau_")
    }
    expect_identical(object$priors$scale$log_tau_study_effect$intercept$parameters$location, 0.3)
    expect_identical(object$priors$scale$log_tau_observation_effect$intercept$parameters$location, 0.4)
    if (!is.null(probability)) {
      allocations <- object$priors$random$allocation
      expect_length(allocations, 2L)
      expect_true(all(vapply(allocations, function(a) is.null(a$sd), logical(1))))
      expect_equal(unlist(lapply(allocations, function(a) lapply(a$inclusion, mean))),
        c(study_effect = probability, observation_effect = probability))
    }
  }
  expect_silent(check_binding(do.call(brma.mv, args)))
  for (name in c("omitted", "weighted", "null_NULL", "alt_FALSE")) {
    selected <- args
    expected <- switch(name, omitted = 0.5, weighted = 0.75, null_NULL = 1, alt_FALSE = 0)
    if (name %in% c("weighted", "null_NULL")) selected$prior_heterogeneity <- global
    if (name == "weighted") selected$prior_heterogeneity_null <- zero
    if (name == "null_NULL") selected["prior_heterogeneity_null"] <- list(NULL)
    if (name == "alt_FALSE") selected$prior_heterogeneity <- FALSE
    expect_silent(check_binding(do.call(BMA.mv, selected), expected))
  }
  selected <- args
  selected$prior_heterogeneity <- global
  expect_error(do.call(brma.mv, selected), "prior_heterogeneity", class = "RoBMA_prior_intercept_conflict", fixed = TRUE)
  selected$scale <- list(study = ~ x)
  selected$prior_scale <- selected$prior_scale["study"]
  expect_silent(partial <- do.call(brma.mv, selected))
  design <- .fitted_formula_design(partial, "mu", required = TRUE)
  sd_names <- unique(design$random_effects[[2L]]$sd_parameter_names)
  binding_args <- .create_jags_formula_args(partial$data, partial$priors)
  compiled <- BayesTools::JAGS_formula(
    formula = binding_args$formula_list$mu, parameter = "mu", data = binding_args$formula_data_list$mu,
    prior_list = binding_args$formula_prior_list$mu, formula_scale = binding_args$formula_scale_list$mu,
    prior_random = binding_args$formula_random_prior_list$mu,
    random_effects_compile = binding_args$formula_random_effects_compile_list$mu)
  expect_length(sd_names, 1L)
  expect_identical(compiled$prior_list[[sd_names]]$distribution, global$distribution)
  expect_identical(compiled$prior_list[[sd_names]]$parameters, global$parameters)
  expect_identical(partial$priors$scale$log_tau_study_effect$intercept$parameters$location, 0.3)
})

.bma_mv_input_data <- function() {

  data.frame(
    yi    = c(0.10, 0.25, 0.15, 0.30),
    study = factor(c("s1", "s1", "s2", "s2")),
    esid  = factor(seq_len(4L)),
    obs   = factor(seq_len(4L)),
    group = factor(c("a", "b", "a", "b")),
    x     = c(0, 1, 0, 1),
    z     = c(-1, -1, 1, 1)
  )
}

.bma_mv_input_V <- function() {

  diag(c(0.04, 0.06, 0.05, 0.08))
}

.bma_mv_gate_priors <- function(prior_random) {

  allocations <- prior_random[["allocation"]]
  unlist(lapply(allocations, function(allocation) {
    allocation[["inclusion"]]
  }), recursive = FALSE)
}


test_that("BMA.mv stages multivariate product-space objects", {

  data <- .bma_mv_input_data()
  V    <- .bma_mv_input_V()
  object <- BMA.mv(
    yi = yi, V = V, random = ~ 1 | obs, data = data,
    measure = "GEN", prior_unit_information_sd = 1,
    only_priors = TRUE
  )

  expect_identical(
    class(object),
    c("only_priors.brma", "BMA.mv", "RoBMA", "brma.mv", "brma.norm", "brma")
  )
  expect_null(object[["priors"]][["outcome"]][["bias"]])
  expect_s3_class(object[["priors"]][["random"]], "prior_random")

  allocation <- object[["priors"]][["random"]][["allocation"]][[1L]]
  expect_length(allocation[["terms"]], 1L)
  expect_null(allocation[["weights"]])
  expect_equal(mean(allocation[["inclusion"]][[1L]]), 0.5)
})


test_that("BMA.mv uses independent default gates without renormalizing allocations", {

  data <- .bma_mv_input_data()
  object <- BMA.mv(
    yi = yi, V = .bma_mv_input_V(),
    random = list(study = ~ 1 | study, observation = ~ 1 | obs),
    data = data, measure = "GEN", prior_unit_information_sd = 1,
    only_priors = TRUE
  )
  allocation <- object[["priors"]][["random"]][["allocation"]][[1L]]
  gates      <- .bma_mv_gate_priors(object[["priors"]][["random"]])

  expect_equal(unname(allocation[["terms"]]), c("study", "observation"))
  expect_s3_class(allocation[["weights"]], "prior.simplex")
  expect_equal(allocation[["weights"]][["parameters"]][["alpha"]], c(1, 1))
  expect_equal(unname(vapply(gates, mean, numeric(1))), c(0.5, 0.5))

  formula_design <- .compile_prior_random_formula_design(
    data           = object[["data"]],
    prior_location = object[["priors"]][["location"]],
    prior_random   = object[["priors"]][["random"]]
  )
  allocation_compiled <- formula_design[["random_allocations"]][[1L]]
  indicators <- vapply(
    allocation_compiled[["inclusion"]],
    `[[`,
    character(1),
    "indicator_name"
  )
  expect_length(unique(indicators), 2L)
  expect_setequal(names(indicators), c("study", "observation"))
})


test_that("BMA.mv decomposes heterogeneity mixture weights into gate and slab", {

  data <- .bma_mv_input_data()
  alternatives <- list(
    BayesTools::prior(
      "normal", parameters = list(mean = 0, sd = 0.25),
      truncation = list(0, Inf), prior_weights = 3
    ),
    BayesTools::prior(
      "normal", parameters = list(mean = 0, sd = 0.75),
      truncation = list(0, Inf), prior_weights = 2
    )
  )
  null <- BayesTools::prior(
    "spike", parameters = list(location = 0), prior_weights = 5
  )
  object <- BMA.mv(
    yi = yi, V = .bma_mv_input_V(), random = ~ 1 | obs,
    prior_heterogeneity = alternatives,
    prior_heterogeneity_null = null,
    data = data, measure = "GEN", prior_unit_information_sd = 1,
    only_priors = TRUE
  )
  allocation <- object[["priors"]][["random"]][["allocation"]][[1L]]
  slab       <- allocation[["sd"]]
  gate       <- allocation[["inclusion"]][[1L]]

  expect_s3_class(slab, "prior.mixture")
  expect_equal(attr(slab, "prior_weights"), c(3, 2))
  expect_equal(mean(gate), 0.5)
})


test_that("BMA.mv supports fixed inclusion and exclusion gates", {

  data <- .bma_mv_input_data()
  included <- BMA.mv(
    yi = yi, V = .bma_mv_input_V(), random = ~ 1 | obs,
    prior_heterogeneity_null = FALSE,
    data = data, measure = "GEN", prior_unit_information_sd = 1,
    only_priors = TRUE
  )
  excluded <- BMA.mv(
    yi = yi, V = .bma_mv_input_V(), random = ~ 1 | obs,
    prior_heterogeneity = FALSE,
    data = data, measure = "GEN", prior_unit_information_sd = 1,
    only_priors = TRUE
  )

  expect_equal(mean(.bma_mv_gate_priors(included[["priors"]][["random"]])[[1L]]), 1)
  expect_equal(mean(.bma_mv_gate_priors(excluded[["priors"]][["random"]])[[1L]]), 0)

  included_null <- BMA.mv(
    yi = yi, V = .bma_mv_input_V(), random = ~ 1 | obs,
    prior_heterogeneity_null = NULL,
    data = data, measure = "GEN", prior_unit_information_sd = 1,
    only_priors = TRUE
  )
  excluded_null <- BMA.mv(
    yi = yi, V = .bma_mv_input_V(), random = ~ 1 | obs,
    prior_heterogeneity = NULL,
    data = data, measure = "GEN", prior_unit_information_sd = 1,
    only_priors = TRUE
  )

  expect_equal(mean(.bma_mv_gate_priors(included_null[["priors"]][["random"]])[[1L]]), 1)
  expect_equal(mean(.bma_mv_gate_priors(excluded_null[["priors"]][["random"]])[[1L]]), 0)
})


test_that("BMA.mv gates nested structures before their within-component split", {

  data <- .bma_mv_input_data()
  object <- BMA.mv(
    yi = yi, V = .bma_mv_input_V(), random = ~ 1 | study / esid,
    data = data, measure = "GEN", prior_unit_information_sd = 1,
    only_priors = TRUE
  )
  allocations <- object[["priors"]][["random"]][["allocation"]]

  expect_length(allocations, 2L)
  expect_null(allocations[[1L]][["weights"]])
  expect_length(allocations[[1L]][["inclusion"]], 1L)
  expect_s3_class(allocations[[2L]][["weights"]], "prior.simplex")
  expect_length(allocations[[2L]][["terms"]], 2L)
  expect_equal(
    allocations[[2L]][["parent"]][["allocation"]],
    allocations[[1L]][["name"]]
  )

  inclusion_map <- .random_component_inclusion_map(object)
  gate <- unique(unlist(inclusion_map, use.names = FALSE))
  expect_length(gate, 1L)
  expect_identical(inclusion_map[["study"]], gate)
  expect_identical(inclusion_map[["esid_study"]], gate)

  conditioning <- .random_component_conditioning_parameters(
    object     = object,
    components = c("study", "esid_study", "total")
  )
  expect_identical(conditioning[["study"]], gate)
  expect_identical(conditioning[["esid_study"]], gate)
  expect_identical(conditioning[["total"]], gate)
  expect_error(
    .random_component_conditioning_parameters(
      object     = object,
      components = "missing"
    ),
    "No inclusion gate is available for random-effect component 'missing'.",
    fixed = TRUE
  )
})


test_that("BMA.mv gates a random-coefficient block before its SD split", {

  data <- .bma_mv_input_data()
  object <- BMA.mv(
    yi = yi, V = .bma_mv_input_V(),
    random = ~ us(0 + group | study),
    prior_heterogeneity = BayesTools::prior_random(
      study = BayesTools::random_block(
        contrasts = c(group = "independent")
      )
    ),
    data = data, measure = "GEN", prior_unit_information_sd = 1,
    only_priors = TRUE
  )
  allocations <- object[["priors"]][["random"]][["allocation"]]

  expect_length(allocations, 2L)
  expect_length(allocations[[1L]][["inclusion"]], 1L)
  expect_null(allocations[[1L]][["weights"]])
  expect_identical(allocations[[2L]][["target"]], "sd_component")
  expect_identical(allocations[[2L]][["scale"]], "mean_variance")
  expect_equal(
    allocations[[2L]][["weights"]][["parameters"]][["alpha"]],
    c(1, 1)
  )
})


test_that("BMA.mv hides a sole implicit random-component name", {

  data <- .bma_mv_input_data()
  object <- BMA.mv(
    yi = yi, V = .bma_mv_input_V(),
    random = ~ 1 | study / obs,
    data = data, measure = "GEN", prior_unit_information_sd = 1,
    only_priors = TRUE
  )
  allocations <- object[["priors"]][["random"]][["allocation"]]

  expect_identical(allocations[[1L]][["component_names"]], "")
  expect_identical(allocations[[2L]][["display_name"]], "")
  expect_false(any(grepl(
    "component 1",
    unlist(lapply(allocations, `[[`, "display_name")),
    fixed = TRUE
  )))
})


test_that("BMA.mv gates component-specific scale regressions independently", {

  data <- .bma_mv_input_data()
  object <- BMA.mv(
    yi = yi, V = .bma_mv_input_V(),
    random = list(study = ~ 1 | study, observation = ~ 1 | obs),
    scale = list(study = ~ x, observation = ~ 1),
    data = data, measure = "GEN", prior_unit_information_sd = 1,
    only_priors = TRUE
  )
  allocations <- object[["priors"]][["random"]][["allocation"]]

  expect_identical(
    .data_scale_formula_sources(object[["data"]]),
    c(study = "tau_study", observation = "tau_observation")
  )
  expect_length(allocations, 2L)
  expect_true(all(vapply(
    allocations,
    function(allocation) length(allocation[["inclusion"]]) == 1L,
    logical(1)
  )))
  expect_true(all(vapply(
    allocations,
    function(allocation) is.null(allocation[["weights"]]),
    logical(1)
  )))
  expect_named(
    .random_component_inclusion_map(object),
    c("study", "observation")
  )
})


test_that("BMA.mv accepts block scale aliases for separately named components", {

  data <- .bma_mv_input_data()
  object <- BMA.mv(
    yi = yi, V = .bma_mv_input_V(),
    random = list(
      study_effect = ~ 1 | study,
      observation_effect = ~ 1 | obs
    ),
    scale = list(study = ~ x, obs = ~ 1),
    data = data, measure = "GEN", prior_unit_information_sd = 1,
    only_priors = TRUE
  )
  allocations <- object[["priors"]][["random"]][["allocation"]]

  expect_identical(
    .data_scale_formula_sources(object[["data"]]),
    c(
      study_effect       = "tau_study_effect",
      observation_effect = "tau_observation_effect"
    )
  )
  expect_identical(
    vapply(allocations, function(allocation) {
      allocation[["sd_source"]][["name"]]
    }, character(1)),
    c("tau_study_effect", "tau_observation_effect")
  )
})


test_that("BMA.mv preserves structural prior_random overrides", {

  data <- .bma_mv_input_data()
  override <- BayesTools::prior_random(
    study = BayesTools::random_block(parameterization = "noncentered")
  )
  object <- BMA.mv(
    yi = yi, V = .bma_mv_input_V(), random = ~ us(1 + x | study),
    prior_heterogeneity = override,
    data = data, measure = "GEN", prior_unit_information_sd = 1,
    only_priors = TRUE
  )

  expect_identical(
    object[["priors"]][["random"]][["blocks"]][["study"]][["parameterization"]],
    "noncentered"
  )
  expect_length(.bma_mv_gate_priors(object[["priors"]][["random"]]), 1L)

  expect_error(
    BMA.mv(
      yi = yi, V = .bma_mv_input_V(), random = ~ 1 | obs,
      prior_heterogeneity = BayesTools::prior_random(
        sd = BayesTools::prior(
          "normal", list(mean = 0, sd = 1), truncation = list(0, Inf)
        )
      ),
      data = data, measure = "GEN", prior_unit_information_sd = 1,
      only_priors = TRUE
    ),
    "cannot add component gates to a user-supplied random-effect scale architecture",
    fixed = TRUE
  )
})


test_that("BMA.mv singular known V requires structurally active regularization", {

  data <- data.frame(
    yi  = c(0.10, 0.25),
    obs = seq_len(2L)
  )
  singular_V <- matrix(0.04, nrow = 2L, ncol = 2L)

  suppressWarnings(expect_error(
    BMA.mv(
      yi = yi, V = singular_V, random = ~ 1 | obs,
      known_v_parameterization = "block_mvn",
      data = data, measure = "GEN", prior_unit_information_sd = 1,
      only_priors = TRUE
    ),
    "singular",
    ignore.case = TRUE
  ))
  suppressWarnings(expect_no_error(BMA.mv(
    yi = yi, V = singular_V, random = ~ 1 | obs,
    known_v_parameterization = "block_mvn",
    prior_heterogeneity_null = FALSE,
    data = data, measure = "GEN", prior_unit_information_sd = 1,
    only_priors = TRUE
  )))
})


test_that("BMA.mv rejects publication-bias arguments", {

  data <- .bma_mv_input_data()
  expect_error(
    BMA.mv(
      yi = yi, V = .bma_mv_input_V(), random = ~ 1 | obs,
      prior_bias = BayesTools::prior_none(),
      data = data, measure = "GEN", prior_unit_information_sd = 1,
      only_priors = TRUE
    ),
    "Selection/publication-bias priors are not supported in BMA.mv().",
    fixed = TRUE
  )
})

.bma_mv_collision_object <- function(case) {

  dat <- data.frame(
    yi = c(.08, .13, .18, .20, .01, .05),
    study = factor(rep(c("s1", "s2", "s3"), each = 2)),
    esid = factor(seq_len(6)), lab = factor(rep(c("l1", "l2"), 3)),
    x = c(-1, 1, -1, 1, -1, 1)
  )
  random <- list(study = ~ random(1 | study, name = "study") +
    random(1 | esid, name = "esid"))
  if (case == "combined") names(random) <- "combined"
  if (case %in% c("scaled", "suffix")) {
    names(random) <- "a"
    random$a_split <- ~ random(1 | lab, name = "lab")
  }
  if (case == "suffix") {
    random$a_split_1 <- ~ random(1 | lab, name = "suffix_one")
    random$a_split_2 <- ~ random(1 | lab, name = "suffix_two")
  }
  args <- list(yi = quote(yi), V = diag(c(.04, .05, .04, .05, .04, .05)),
    random = random, data = dat, measure = "GEN", prior_unit_information_sd = 1,
    only_priors = TRUE, silent = TRUE)
  if (case %in% c("scaled", "suffix")) {
    args$scale <- stats::setNames(rep(list(~ x), length(random)), names(random))
  }
  do.call(if (case == "brma") brma.mv else BMA.mv, args)
}

.expect_bma_mv_allocation_graph <- function(object) {

  setup <- review_compile_object(object)
  allocations <- setup$compiled$mu$formula_design$random_allocations
  labels <- vapply(allocations, `[[`, character(1), "label")
  expect_length(unique(labels), length(labels))
  expect_identical(names(allocations), unname(labels))
  for (allocation in allocations) {
    parent <- allocation$parent
    if (is.null(parent)) next
    expect_true(parent$allocation %in% labels)
    expect_true(parent$component %in% names(allocations[[parent$allocation]]$terms))
  }
  terms <- setup$compiled$mu$formula_design$random_effects
  expect_true(all(vapply(terms, function(term) !is.null(term$sd_binding), logical(1))))
  allocations
}

test_that("a consumed component keeps its public key when it names an actual block", {

  object <- .bma_mv_collision_object("study")
  allocations <- .expect_bma_mv_allocation_graph(object)
  root <- allocations$heterogeneity
  expect_identical(names(root$terms), "study")
  expect_identical(root$component_names, "study")
  expect_false(unname(root$terms) %in% c("study", "esid"))
  expect_identical(allocations$study_split$parent$allocation, root$label)
  expect_identical(allocations$study_split$parent$component, "study")
  expect_identical(allocations$study_split$display_name, "study")
  expect_identical(unname(allocations$study_split$terms), c("study", "esid"))
  expect_equal(mean(root$inclusion$study$prior), .5)
})

test_that("valid allocation names, parent keys and priors survive noncolliding controls", {

  ordinary <- .expect_bma_mv_allocation_graph(.bma_mv_collision_object("brma"))
  expect_identical(names(ordinary), "heterogeneity")
  expect_identical(unname(ordinary$heterogeneity$terms), c("study", "esid"))
  expect_identical(ordinary$heterogeneity$display_name, "")
  combined <- .expect_bma_mv_allocation_graph(.bma_mv_collision_object("combined"))
  expect_identical(names(combined), c("heterogeneity", "combined_split"))
  expect_identical(combined$heterogeneity$terms, c(combined = "combined"))
  expect_identical(combined$combined_split$parent$allocation, "heterogeneity")
  expect_identical(combined$combined_split$parent$component, "combined")
  expect_identical(combined$combined_split$display_name, "combined")
  expect_equal(mean(combined$heterogeneity$inclusion$combined$prior), .5)
  prior <- .bma_mv_collision_object("combined")$priors$random$allocation[[1L]]$sd
  expect_identical(prior$distribution, "normal")
  expect_equal(prior$parameters$sd, .25)
  expect_equal(.bma_mv_collision_object("combined")$priors$random$allocation[[2L]]$weights$parameters$alpha, c(1, 1))
})

for (.collision_case in c("scaled", "suffix")) {
  local({
    case <- .collision_case
    test_that(paste("component scale graph avoids declared allocation collisions:", case), {

      object <- .bma_mv_collision_object(case)
      allocations <- .expect_bma_mv_allocation_graph(object)
      split <- Filter(function(a) identical(unname(a$terms), c("study", "esid")), allocations)
      expect_length(split, 1L)
      expect_identical(split[[1L]]$parent$allocation, "a")
      expect_identical(split[[1L]]$parent$component, "a")
      expect_identical(split[[1L]]$display_name, "a")
      # A preferred allocation ID remains valid when it merely names a public component.
      expect_identical(allocations$a$label, "a")
      labels <- vapply(allocations, `[[`, character(1), "label")
      expect_identical(split[[1L]]$label, "a_split")
      for (owner in if (case == "suffix") c("a_split", "a_split_1", "a_split_2") else "a_split") {
        root <- Filter(function(a) identical(a$display_name, owner) && is.null(a$parent), allocations)
        expect_length(root, 1L)
        if (owner == "a_split") {
          expect_false(root[[1L]]$label %in% (if (case == "suffix") c("a_split", "a_split_1", "a_split_2") else "a_split"))
        }
        expect_identical(names(root[[1L]]$terms), owner)
        expect_equal(mean(root[[1L]]$inclusion[[owner]]$prior), .5)
      }
      expect_identical(names(.data_scale_formula_sources(object$data)),
        if (case == "suffix") c("a", "a_split", "a_split_1", "a_split_2") else c("a", "a_split"))
    })
  })
}

test_that("an unused heterogeneity root does not reserve a component allocation ID", {

  data <- .bma_mv_input_data()
  object <- brma.mv(yi = yi, V = .bma_mv_input_V(), data = data,
    random = list(other = ~ random(1 | study, name = "study"),
      heterogeneity = ~ random(1 | group, name = "group") + random(1 | esid, name = "esid")),
    scale = list(other = ~ x, heterogeneity = ~ x),
    measure = "GEN", prior_unit_information_sd = 1, only_priors = TRUE, silent = TRUE)
  allocations <- .expect_bma_mv_allocation_graph(object)
  expect_identical(names(allocations), "heterogeneity")
  expect_identical(allocations$heterogeneity$label, "heterogeneity")
})



test_that("BMA.mv mean-difference factor slopes are exchangeable units with one SD component", {

  dat <- data.frame(
    yi    = c(-0.2, 0.1, 0.3, 0.5, 0.0, 0.2, 0.4, -0.1, 0.25),
    group = factor(rep(c("a", "b", "c"), 3L)),
    pair  = factor(rep(c("u", "v", "u"), 3L)),
    study = factor(rep(c("s1", "s2", "s3"), each = 3L))
  )
  stage <- function(mods, random){
    BMA.mv(
      yi = yi, V = diag(0.04, nrow(dat)), mods = mods, random = random,
      data = dat, measure = "GEN", prior_unit_information_sd = 1,
      only_priors = TRUE
    )
  }
  random_term <- function(object){
    .fitted_formula_design(object, "mu", required = TRUE)[["random_effects"]][[1L]]
  }
  component_weights <- function(object){
    allocations <- object[["priors"]][["random"]][["allocation"]]
    components  <- Filter(function(allocation) {
      identical(allocation[["target"]], "sd_component")
    }, allocations)
    expect_length(components, 1L)
    components[[1L]][["weights"]][["parameters"]][["alpha"]]
  }

  # three levels: two contrast columns that share one SD, so the block has
  # two SD components (the intercept and the factor)
  for (random in list(~ diag(1 + group | study), ~ us(1 + group | study))) {
    object <- stage(~ group, random)
    term   <- random_term(object)
    expect_identical(term[["n_columns"]], 3L, info = deparse(random))
    expect_identical(term[["sd_leaves"]][["n_leaves"]], 2L, info = deparse(random))
    expect_identical(
      unique(term[["sd_parameter_names"]]),
      # the allocated SD component nodes, one per leaf
      paste0(term[["parameter_stem"]], c("_intercept", "_group[1]")),
      info = deparse(random)
    )
    expect_equal(component_weights(object), c(1, 1), info = deparse(random))
  }
  expect_identical(
    random_term(stage(~ group, ~ us(1 + group | study)))[["correlation"]][["unit_columns"]],
    2:3
  )

  # a block that is one unit has one SD component and no correlation
  object <- stage(~ group, ~ us(0 + group | study))
  term   <- random_term(object)
  expect_identical(term[["sd_leaves"]][["n_leaves"]], 1L)
  expect_null(term[["correlation"]])
  expect_false(any(vapply(object[["priors"]][["random"]][["allocation"]], function(allocation) {
    identical(allocation[["target"]], "sd_component")
  }, logical(1))))

  # a two-level factor is one column: an ordinary column of the LKJ block
  term <- random_term(stage(~ pair, ~ us(1 + pair | study)))
  expect_identical(term[["n_columns"]], 2L)
  expect_identical(term[["correlation"]][["unit_columns"]], integer())

  # the SD component of a column after the unit is found at its own column
  dat$x <- seq(-1, 1, length.out = nrow(dat))
  term <- random_term(stage(~ group + x, ~ diag(1 + group + x | study)))
  expect_identical(term[["sd_leaves"]][["leaf_index_by_column"]], c(1L, 2L, 2L, 3L))
  expect_identical(vapply(c("intercept", "group", "x"), function(component) {
    .brma_random_parameter_component_column(term, list(random_component = component))
  }, integer(1)), c(intercept = 1L, group = 2L, x = 4L))
})
