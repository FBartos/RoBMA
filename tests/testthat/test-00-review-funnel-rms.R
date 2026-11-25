test_that("plug-in RMS spread is stable at zero and extreme finite scales", {

  local_mocked_bindings(
    .funnel_joint_model_groups = function(x, posterior_samples) rep(1L, nrow(posterior_samples)),
    .funnel_setup_from_samples = function(x, posterior_samples, tau_samples, ...) tau_samples,
    .package = "RoBMA"
  )
  for (scale in c(0, 1e-200, 1, 1e200)) {
    common <- list(posterior_samples = matrix(0, 2L, 1L), tau = scale * c(1, 3))
    actual <- .funnel_sampling_setup(list(), TRUE, FALSE, Inf, "plugin", common)
    expect_equal(unname(actual / if (scale == 0) 1 else scale),
                  if (scale == 0) 0 else sqrt(5), tolerance = 1e-14)
  }
})
