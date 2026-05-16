# ============================================================================ #
# IWMDE File Map
# ============================================================================ #

# The old monolithic iwmde.R implementation is split by concern:
# - iwmde-api.R: exported diagnostics and top-level orchestration
# - iwmde-parameters.R: parameter specs, active components, and focal priors
# - iwmde-grid.R: support transforms, normalization grids, and display grids
# - iwmde-density.R: row states, density aggregation, Chen weights, and MCSE helpers
# - iwmde-q-grid.R: q-grid batch dispatch and GLMM conditional batches
# - iwmde-predictor.R: predictor/location fast paths and candidate construction
# - iwmde-replacement.R: scalar, linear, and indexed replacement maps
# - iwmde-likelihood.R: likelihood and prior ordinates for replacement rows
# - iwmde-active.R: active prior branches, row parameters, selection, and formula inputs
# - iwmde-plot.R: KDE/histogram data and base plotting
#
# This file is intentionally kept as an index for maintainers and old mental maps.
# ============================================================================ #

