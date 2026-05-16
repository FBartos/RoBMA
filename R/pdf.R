# ============================================================================ #
# Log-Likelihood File Map
# ============================================================================ #
#
# The old monolithic pdf.R implementation is split by concern:
# - pdf-utils.R: shared numerical, weighting, native-availability, and quadrature helpers
# - pdf-outcome-normal.R: normal and selected-normal estimate likelihood kernels
# - pdf-outcome-glmm.R: binomial and Poisson GLMM estimate likelihood kernels
# - log-lik.R: estimate-unit dispatch, setup, and posterior-sample evaluation
# - log-lik-cluster.R: cluster-unit dispatch and setup
# - log-lik-cluster-normal.R: normal cluster-unit likelihood implementations
# - log-lik-cluster-glmm.R: GLMM cluster-unit likelihood implementations
#
# This file is intentionally kept as an index for maintainers and old mental maps.
# ============================================================================ #
