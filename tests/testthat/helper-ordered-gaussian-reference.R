# Independent declared-model and rounded-input references are documented beside
# the exact binary fixture. These helpers only check captured input identities
# and the explicit physical mean/variance perturbation; they do not evaluate an
# alternative likelihood or change the original output-scale allowance.
.ordered_gaussian_reference <- function(case) {

  readRDS(testthat::test_path("fixtures", "ordered-gaussian-reference.rds"))[[case]]
}


.ordered_gaussian_reference_inputs <- function(reference, block) {

  input <- reference[["inputs"]]
  expect_identical(block[["setup"]][["yi"]], input[["yi"]])
  for (field in c("mu0", "mu_a", "mu_b", "mu_c", "variance")) {
    expect_identical(as.numeric(block[[field]]), as.numeric(input[[field]]))
  }
  expect_identical(as.numeric(block[["weights"]]), as.numeric(input[["powers"]]))
}


.ordered_gaussian_rounding_bound <- function(reference, weights, variance, selection = NULL) {

  delta <- abs(reference[["delta_mean"]])
  residual <- reference[["residual_declared"]]
  v <- matrix(variance, nrow(delta), ncol(delta), byrow = TRUE)
  dv <- reference[["delta_variance"]]
  v2 <- v + dv
  # |L(m+delta,v)-L(m,v)| plus the variance change at the perturbed mean.
  gaussian <- abs(residual) * delta / v + delta^2 / (2 * v) +
    .5 * (abs(log1p(dv / v)) + (residual - reference[["delta_mean"]])^2 * abs(dv) / (v * v2))
  selected <- 0
  if (!is.null(selection)) {
    omega <- as.numeric(selection[["omega"]])
    variation <- sum(abs(diff(omega)))
    minimum <- min(omega)
    lower_v <- pmin(v, v2)
    # Z >= min(omega); telescoping step-CDF sensitivities bound both changes.
    selected <- variation * delta / (sqrt(2 * pi * lower_v) * minimum) +
      variation * abs(dv) / (2 * lower_v * sqrt(2 * pi * exp(1)) * minimum)
  }
  rowSums(sweep(gaussian + selected, 2L, weights, "*"))
}
