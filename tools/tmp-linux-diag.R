# TEMPORARY diagnostic for the Linux standard-profile failures; removed before merge.
suppressMessages(pkgload::load_all(".", quiet = TRUE))
options(width = 200)
definitions <- function(file) {
  for (e in parse(file, keep.source = FALSE)) if (is.call(e) && identical(e[[1L]], as.name("<-"))) eval(e, globalenv())
}
zplot_check <- function(label) {
  definitions("tests/testthat/test-00-zplot-threshold-batch.R")
  sei <- c(.3, .5, .3, .3, .8, .5); S <- 6L
  mean <- matrix(seq(-.4, .7, length.out = S * length(sei)), S); mean[, 3L] <- mean[, 1L]
  sd <- matrix(rep(sei + .2, each = S), S)
  latent_sd <- matrix(seq(0, .3, length.out = length(mean)), S); latent_sd[, 3L] <- latent_sd[, 1L]
  rule <- list(nodes = c(-3, -.7, 0, .2, 2.5), weights = c(0, .1, .2, .3, .4))
  kind <- function(x) ifelse(is.nan(x), "NaN", ifelse(is.na(x), "NA", ifelse(is.finite(x), "", "Inf")))
  for (direction in c("positive", "negative")) for (normal in list(FALSE, TRUE, c(FALSE, TRUE))) for (fo in c(FALSE, TRUE)) {
    selection <- .zplot_threshold_batch_context(S, sei, direction, normal)
    selection$omega[2L, ] <- selection$omega[2L, ] * 1e280
    selection$omega[3L, ] <- selection$omega[3L, ] * 1e-280
    selection$omega[4L, ] <- c(1, 0, 0)
    selection$omega[5L, ] <- c(.1, .8, .2)
    args <- list(1.96, mean, sd, sei, selection, mean + .125, latent_sd, rule, fo)
    a <- do.call(.zplot_selnorm_threshold_quadrature, args); b <- do.call(.zplot_threshold_batch_reference, args)
    if (!identical(a, b)) {
      cat(label, direction, paste(normal, collapse = "/"), fo, "NOT identical\n")
      for (n in names(b)) if (!identical(a[[n]], b[[n]])) {
        print(rbind(native = a[[n]], reference = b[[n]]), digits = 17)
        cat("kinds native:", kind(a[[n]]), " reference:", kind(b[[n]]), "\n")
        cat("bits differ:", which(!mapply(identical, a[[n]], b[[n]])), "\n")
      }
    }
  }
  cat(label, "zplot check done\n")
}
zplot_check("fresh")
Sys.setenv(NOT_CRAN = "true")
r <- as.data.frame(testthat::test_file("tests/testthat/test-00-zplot-threshold-batch.R", reporter = "summary",
                                       package = "RoBMA", load_package = "none"))
cat("zplot test file alone: failed", sum(r$failed), "passed", sum(r$passed), "\n")
cat("MXCSR-sensitive check: 2^-1074 * 1 =", 2^-1074 * 1, " .Machine$double.xmin/2 =", .Machine$double.xmin / 2, "\n")
