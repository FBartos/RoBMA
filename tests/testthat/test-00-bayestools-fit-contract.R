context("BayesTools fitted metadata contract")


test_that("fitted formula identities come from the versioned name map", {

  fit <- structure(list(), class = "BayesTools_fit")
  object <- list(fit = fit)
  checked <- NULL
  name_map <- data.frame(
    encoded_name      = c("BT1_intercept", "BT1_first", "BT1_second"),
    jags_name         = c("mu_intercept", "opaque_backend_1", "opaque_backend_2"),
    kind              = "fixed",
    formula_parameter = "mu",
    term              = c("intercept", "x:a_b", "x__xXx__a_b"),
    role              = "coefficient",
    stringsAsFactors  = FALSE
  )
  testthat::local_mocked_bindings(
    JAGS_validate_fit_contract = function(fit, requires) {
      checked <<- requires
      invisible(TRUE)
    },
    JAGS_formula_name_map = function(fit, parameter) name_map,
    JAGS_parameter_names = function(...) {
      stop("backend names must not be reconstructed")
    },
    .package = "BayesTools"
  )

  observed <- list()
  add <- function(parameter, component, term, aliases,
                  source, formula_parameter) {
    observed[[length(observed) + 1L]] <<- list(
      parameter         = parameter,
      component         = component,
      term              = term,
      aliases           = aliases,
      source            = source,
      formula_parameter = formula_parameter
    )
  }
  expect_silent(.brma_parameter_catalog_terms(
    object            = object,
    model_parameter   = "mu",
    component         = "mods",
    formula_parameter = "mu",
    source            = "mods",
    add               = add
  ))

  expect_setequal(
    checked,
    c(
      "name_encoding",
      "formula_name_map",
      "formula_design",
      "parameter_map"
    )
  )
  expect_identical(
    vapply(observed, `[[`, character(1), "parameter"),
    c("opaque_backend_1", "opaque_backend_2")
  )
  expect_identical(
    vapply(observed, `[[`, character(1), "term"),
    c("x:a_b", "x__xXx__a_b")
  )
})


test_that("fitted metadata consumers fail closed on stale contracts", {

  # A BayesTools fit without the fit contract of this version (such as a
  # fit of BayesTools 0.3.0): BayesTools' refit error passes through with
  # its classes.
  object <- list(fit = structure(list(), class = "BayesTools_fit"))
  error  <- tryCatch(.fitted_formula_name_map(object, "mu"), error = identity)

  expect_true(inherits(error, "BayesTools_refit_required"))
  expect_false(inherits(error, "RoBMA_refit_required"))
  expect_match(conditionMessage(error), "Refit the model", fixed = TRUE)
})


test_that("RoBMA refit errors inherit BayesTools_refit_required", {

  refit <- c("RoBMA_refit_required", "BayesTools_refit_required",
             "error", "condition")

  # The helper builds the message as stop() does and raises no call.
  error <- tryCatch(
    .stop_refit_required("Metadata of '", c("a", "b"), "' are missing: ", 2L, "."),
    error = identity
  )
  expect_identical(class(error), refit)
  expect_identical(conditionMessage(error), "Metadata of 'ab' are missing: 2.")
  expect_null(conditionCall(error))
  # A handler of BayesTools' refit class catches RoBMA's refit errors.
  expect_identical(
    tryCatch(
      .stop_refit_required("Refit the model."),
      BayesTools_refit_required = function(error) "caught"
    ),
    "caught"
  )

  # RoBMA's refit stops, with their complete messages.
  catch <- function(expr) tryCatch(expr, error = identity)
  fit   <- structure(list(), class = "BayesTools_fit")
  attr(fit, "formula_design") <- list(log_tau = list(parameter = "log_tau"))
  testthat::local_mocked_bindings(
    JAGS_validate_fit_contract = function(...) invisible(TRUE),
    JAGS_formula_name_map      = function(...) NULL,
    .package = "BayesTools"
  )
  errors <- list(
    gate = catch(.brma_validate_fit_contract(
      list(fit = NULL),
      requires = "parameter_map"
    )),
    name_map = catch(.fitted_formula_name_map(list(fit = fit), "mu")),
    design = catch(.fitted_formula_design(list(fit = fit), "mu")),
    term_labels = catch(.formula_design_term_labels(list(
      parameter         = "mu",
      model_terms       = c("x", "g"),
      model_term_labels = "x"
    ))),
    grouped_scale = catch(.brma_parameter_catalog_formula_quantity(
      catalog            = NULL,
      map_row            = NULL,
      coordinates        = data.frame(fitted_scale = c("a", "b"),
                                      display_scale = "a"),
      semantic_component = "mods"
    )),
    random_parts = catch(.brma_random_parameter_io_parts(list(NULL))),
    random_support = catch(.brma_random_parameter_catalog_support(list(
      entry = list(parameter = "tau", selection = list(quantities = NULL))
    ))),
    allocated_sd = catch(.marginalized_random_effect_allocated_sd_nodes(list(
      block_name = "study",
      sd_binding = list(
        true_allocation = TRUE,
        allocations     = list(list(
          source               = list(shape = "scalar"),
          target               = "sd_component",
          leaf_names           = character(),
          leaf_index_by_column = integer()
        ))
      )
    ))),
    summary_rows = catch(.summary_scale_row_labels(
      estimates      = matrix(1, dimnames = list("x", "Mean")),
      object         = NULL,
      formula_prefix = NULL
    ))
  )
  messages <- c(
    gate = paste0(
      "Current BayesTools fitted metadata are unavailable. Refit the model ",
      "with the current RoBMA/BayesTools build."
    ),
    name_map = paste0(
      "Fitted formula name-map metadata for parameter 'mu' are missing. ",
      "Refit the model with the current RoBMA/BayesTools build."
    ),
    design = paste0(
      "Fitted formula design metadata for parameter 'mu' is missing. ",
      "Refit the model with the current RoBMA/BayesTools build."
    ),
    term_labels = paste0(
      "Formula design metadata of 'mu' do not label its model terms. Refit ",
      "the model with the current RoBMA/BayesTools build."
    ),
    grouped_scale = paste0(
      "Grouped formula coefficient coordinates have inconsistent scale ",
      "metadata. Refit the model with the current RoBMA/BayesTools build."
    ),
    random_parts = paste0(
      "Random-effect catalog quantities have no random-effect label parts. ",
      "Refit the model with the current RoBMA/BayesTools build."
    ),
    random_support = paste0(
      "Random-effect quantity 'tau' has no catalog support metadata. Refit ",
      "the model with the current BayesTools version."
    ),
    allocated_sd = paste0(
      "Allocated random-effect SD nodes of block 'study' are unavailable. ",
      "Refit the model with the current RoBMA/BayesTools build."
    ),
    summary_rows = paste0(
      "Scale summary rows have no catalog label parts. Refit the model with ",
      "the current RoBMA/BayesTools build."
    )
  )
  for (site in names(messages)) {
    expect_identical(class(errors[[site]]), refit, info = site)
    expect_identical(conditionMessage(errors[[site]]), messages[[site]],
                     info = site)
  }
})


# String constants of the package sources that contain the word "refit" in
# any case ("refitting" and "refitted" are other words): the file, line, and
# whether a call of .stop_refit_required() encloses the string.
.refit_message_sites <- function(files) {

  sites <- lapply(files, function(file) {
    data <- utils::getParseData(parse(file, keep.source = TRUE, encoding = "UTF-8"))
    strings <- data[
      data[["token"]] == "STR_CONST" &
        grepl("\\brefit\\b", data[["text"]], ignore.case = TRUE, perl = TRUE),
      ,
      drop = FALSE
    ]
    if (nrow(strings) == 0L) {
      return(NULL)
    }
    helper <- data[["parent"]][
      data[["token"]] == "SYMBOL_FUNCTION_CALL" &
        data[["text"]] == ".stop_refit_required"
    ]
    helper_calls <- data[["parent"]][match(helper, data[["id"]])]
    enclosed <- vapply(strings[["id"]], function(id) {
      while (length(id) == 1L && id > 0L) {
        if (id %in% helper_calls) {
          return(TRUE)
        }
        id <- data[["parent"]][data[["id"]] == id]
      }
      FALSE
    }, logical(1))
    data.frame(
      file     = basename(file),
      line     = strings[["line1"]],
      enclosed = enclosed,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, sites)
}


test_that("every message that asks for a refit is raised through .stop_refit_required()", {

  # The scanner finds refit messages of every call form and their enclosure.
  planted <- tempfile(fileext = ".R")
  on.exit(unlink(planted), add = TRUE)
  writeLines(c(
    "f <- function() {",
    "  stop(\"Refit the model.\", call. = FALSE)",
    "  .stop_refit_required(\"Refit the model.\")",
    "  .stop_refit_required(paste0(\"a\", \" or refit the model.\"))",
    "  .hypothesis_refusal(paste0(\"Please \", \"REFIT.\"), \"target\")",
    "  message(\"refitting and refitted are other words\")",
    "}"
  ), planted)
  planted_sites <- .refit_message_sites(planted)
  expect_identical(planted_sites[["line"]], c(2L, 3L, 4L, 5L))
  expect_identical(planted_sites[["enclosed"]], c(FALSE, TRUE, TRUE, FALSE))

  package_r_dir <- testthat::test_path("..", "..", "R")
  skip_if_not(dir.exists(package_r_dir), "Package sources are not available.")
  sites <- .refit_message_sites(
    list.files(package_r_dir, pattern = "[.][Rr]$", full.names = TRUE)
  )
  expect_gt(sum(sites[["enclosed"]]), 0L)
  labels <- paste0(sites[["file"]], ":", sites[["line"]])
  expect_identical(labels[!sites[["enclosed"]]], character())
})
