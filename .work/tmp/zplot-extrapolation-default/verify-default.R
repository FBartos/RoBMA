library(RoBMA)

set.seed(1)

# Small selection model so the fitted and extrapolated curves genuinely differ.
data("dat.lehmann2018", package = "metadat")
fit <- bselmodel(
  yi      = dat.lehmann2018$yi,
  vi      = dat.lehmann2018$vi,
  measure = "SMD",
  sample  = 2000, burnin = 1000, adapt = 500, chains = 2,
  seed    = 1, silent = TRUE
)

.geom_line_count <- function(plot) {

  sum(vapply(plot$layers, function(layer) inherits(layer$geom, "GeomLine"), logical(1)))
}

.geom_line_data <- function(plot) {

  lapply(
    Filter(function(layer) inherits(layer$geom, "GeomLine"), plot$layers),
    function(layer) layer$mapping
  )
}

zp <- as_zplot(fit, max_samples = 200)

default_plot       <- zplot(zp, plot_type = "ggplot", max_samples = 200)
explicit_off_plot  <- zplot(zp, plot_type = "ggplot", max_samples = 200, plot_extrapolation = FALSE)
explicit_on_plot   <- zplot(zp, plot_type = "ggplot", max_samples = 200, plot_extrapolation = TRUE)

cat("GeomLine layers, default             :", .geom_line_count(default_plot), "\n")
cat("GeomLine layers, plot_extrapolation=F:", .geom_line_count(explicit_off_plot), "\n")
cat("GeomLine layers, plot_extrapolation=T:", .geom_line_count(explicit_on_plot), "\n")

stopifnot(
  identical(.geom_line_count(default_plot), .geom_line_count(explicit_off_plot)),
  .geom_line_count(default_plot) == 1L,
  .geom_line_count(explicit_on_plot) == 2L
)

# The fitted curve the default draws must equal the fitted curve of the old
# (extrapolating) display, and must differ from the extrapolated curve.
fitted_density <- lines(zp, as_data = TRUE, max_samples = 200, extrapolate = FALSE, length.out = 41)
extrap_density <- lines(zp, as_data = TRUE, max_samples = 200, extrapolate = TRUE,  length.out = 41)

cat("\nmax |fitted - extrapolated| density:",
    max(abs(fitted_density$y - extrap_density$y)), "\n")
stopifnot(max(abs(fitted_density$y - extrap_density$y)) > 1e-6)

# Base rendering must run for every combination.
file <- tempfile(fileext = ".png")
grDevices::png(file)
zplot(zp, max_samples = 200)
zplot(zp, max_samples = 200, plot_extrapolation = TRUE)
zplot(zp, max_samples = 200, plot_fit = FALSE, plot_extrapolation = TRUE)
zplot(zp, max_samples = 200, plot_fit = FALSE, plot_extrapolation = FALSE)
grDevices::dev.off()
unlink(file)

cat("\nformals(plot.zplot_brma)$plot_extrapolation = ",
    format(formals(RoBMA:::plot.zplot_brma)$plot_extrapolation), "\n")
cat("formals(lines.zplot_brma)$extrapolate       = ",
    format(formals(RoBMA:::lines.zplot_brma)$extrapolate), "\n")

cat("\nAll zplot extrapolation-default checks passed.\n")
