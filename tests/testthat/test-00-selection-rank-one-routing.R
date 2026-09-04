context("Exact selection rank-one routing")
skip_on_cran()

test_that("nested intercept random effects use the exact rank-one kernel", {

  dat <- data.frame(
    yi    = c(0.10, 0.20, 0.05, 0.15),
    vi    = rep(0.01, 4L),
    study = factor(c("a", "a", "b", "b")),
    esid  = factor(seq_len(4L))
  )
  object <- bselmodel.mv(
    yi                        = yi,
    vi                        = vi,
    random                    = ~ 1 | study / esid,
    data                      = dat,
    measure                   = "SMD",
    prior_unit_information_sd = 1,
    selection_likelihood      = "exact",
    only_priors               = TRUE,
    silent                    = TRUE
  )
  setup       <- object[["selection_likelihood"]]
  exact_setup <- .data_exact_selection_setup(object[["data"]])
  syntax      <- .create_model_syntax(object[["data"]], object[["priors"]])

  expect_identical(setup[["exactness"]], "E1")
  expect_identical(
    exact_setup[["random_covariance"]][["representation"]],
    "diagonal_factor"
  )
  expect_match(syntax, "dselnorm_cluster_step", fixed = TRUE)
  expect_false(grepl("dselnorm_mnorm_step", syntax, fixed = TRUE))
  expect_false(grepl("sel_exact_qmc_", syntax, fixed = TRUE))
})


test_that("two-coefficient random covariance uses the same rank-one route", {

  dat <- data.frame(
    yi      = c(-0.2, 0.1, 0.3, 0.5),
    vi      = rep(0.04, 4L),
    outcome = factor(c("sensitivity", "specificity",
                       "sensitivity", "specificity")),
    study   = factor(c("s1", "s1", "s2", "s2"))
  )
  object <- bselmodel.mv(
    yi                        = yi,
    vi                        = vi,
    random                    = ~ us(0 + outcome | study),
    prior_heterogeneity       = BayesTools::prior_random(
      study = BayesTools::random_block(
        contrasts = c(outcome = "independent")
      )
    ),
    data                      = dat,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    selection_likelihood      = "exact",
    only_priors               = TRUE,
    silent                    = TRUE
  )
  setup       <- object[["selection_likelihood"]]
  exact_setup <- .data_exact_selection_setup(object[["data"]])
  syntax      <- .create_model_syntax(object[["data"]], object[["priors"]])

  expect_identical(setup[["exactness"]], "E1")
  expect_identical(
    exact_setup[["random_covariance"]][["representation"]],
    "diagonal_factor"
  )
  expect_match(syntax, "_factor_basis", fixed = TRUE)
  expect_match(syntax, "dselnorm_cluster_step", fixed = TRUE)
  expect_false(grepl("dselnorm_mnorm_step", syntax, fixed = TRUE))
})


test_that("certified factors route by structural rank and dense V stays dense", {

  dat <- data.frame(
    yi    = c(0.10, 0.20, 0.05, 0.15, 0.12, 0.22),
    study = factor(rep(c("a", "b"), each = 3L)),
    esid  = factor(seq_len(6L)),
    time  = rep(1:3, 2L)
  )
  V <- diag(rep(0.01, 6L))
  V[1L, 2L] <- V[2L, 1L] <- 0.002
  V[2L, 3L] <- V[3L, 2L] <- 0.001
  V[4L, 5L] <- V[5L, 4L] <- 0.002
  V[5L, 6L] <- V[6L, 5L] <- 0.001

  correlated <- bselmodel.mv(
    yi                        = yi,
    V                         = V,
    random                    = ~ 1 | study / esid,
    data                      = dat,
    measure                   = "SMD",
    prior_unit_information_sd = 1,
    selection_likelihood      = "exact",
    only_priors               = TRUE,
    silent                    = TRUE
  )
  correlated_setup <- correlated[["selection_likelihood"]]
  correlated_exact_setup <- .data_exact_selection_setup(correlated[["data"]])
  correlated_syntax <- .create_model_syntax(
    correlated[["data"]],
    correlated[["priors"]]
  )
  expect_identical(correlated_setup[["exactness"]], "E2")
  expect_identical(
    correlated_exact_setup[["random_covariance"]][["representation"]],
    "dense"
  )
  expect_match(correlated_syntax, "dselnorm_mnorm_step", fixed = TRUE)
  expect_false(grepl("dselnorm_cluster_step", correlated_syntax, fixed = TRUE))

  higher_rank <- bselmodel.mv(
    yi                        = yi,
    vi                        = rep(0.01, 6L),
    random                    = ~ har(time | study),
    data                      = dat,
    measure                   = "SMD",
    prior_unit_information_sd = 1,
    selection_likelihood      = "exact",
    only_priors               = TRUE,
    silent                    = TRUE
  )
  higher_rank_setup <- higher_rank[["selection_likelihood"]]
  higher_rank_exact_setup <- .data_exact_selection_setup(higher_rank[["data"]])
  higher_rank_syntax <- .create_model_syntax(
    higher_rank[["data"]],
    higher_rank[["priors"]]
  )
  expect_identical(higher_rank_setup[["exactness"]], "EF")
  expect_identical(
    higher_rank_exact_setup[["random_covariance"]][["representation"]],
    "diagonal_factor"
  )
  expect_identical(
    higher_rank_exact_setup[["factor_ranks"]],
    c(2L, 2L)
  )
  expect_match(higher_rank_syntax, "dselnorm_factor_step", fixed = TRUE)
  expect_false(grepl("dselnorm_mnorm_step", higher_rank_syntax, fixed = TRUE))
  expect_false(grepl("dselnorm_cluster_step", higher_rank_syntax, fixed = TRUE))
})


test_that("crossed random dependencies combine through the generic factor plan", {

  dat <- data.frame(
    yi       = c(-.2, .1, .3, .5),
    vi       = rep(.04, 4L),
    study    = factor(c("s1", "s1", "s2", "s2")),
    observer = factor(c("a", "b", "a", "b"))
  )
  object <- bselmodel.mv(
    yi                        = yi,
    vi                        = vi,
    random                    = list(
      study    = ~ 1 | study,
      observer = ~ 1 | observer
    ),
    data                      = dat,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    selection_likelihood      = "exact",
    only_priors               = TRUE,
    silent                    = TRUE
  )
  exact_setup <- .data_exact_selection_setup(object[["data"]])
  syntax <- .create_model_syntax(object[["data"]], object[["priors"]])

  expect_identical(exact_setup[["row_blocks"]], list(1:4))
  expect_identical(exact_setup[["exactness"]], "EF")
  expect_identical(
    exact_setup[["factor_ranks"]],
    4L
  )
  expect_match(syntax, "dselnorm_factor_step", fixed = TRUE)
  expect_false(grepl("dselnorm_mnorm_step", syntax, fixed = TRUE))
})


test_that("diagonal block-list V retains generic factor routing", {

  dat <- data.frame(
    yi    = c(-.2, .1, .3, .5),
    study = factor(c("s1", "s1", "s2", "s2"))
  )
  object <- bselmodel.mv(
    yi                        = yi,
    V                         = list(diag(c(.02, .03)), diag(c(.04, .05))),
    random                    = ~ 1 | study,
    data                      = dat,
    measure                   = "SMD",
    prior_unit_information_sd = 1,
    selection_likelihood      = "exact",
    only_priors               = TRUE,
    silent                    = TRUE
  )
  exact_setup <- .data_exact_selection_setup(object[["data"]])

  expect_identical(exact_setup[["exactness"]], "E1")
  expect_identical(
    exact_setup[["block_methods"]],
    rep("rank_one", 2L)
  )
  expect_identical(
    exact_setup[["random_covariance"]][["representation"]],
    "diagonal_factor"
  )
})
