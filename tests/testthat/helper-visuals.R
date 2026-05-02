expect_brma_plot_snapshot <- function(name, plot) {

  vdiffr::expect_doppelganger(name, plot)
}

expect_side_by_side_plot <- function(name, metafor, brma) {

  vdiffr::expect_doppelganger(name, function() {
    old_par <- graphics::par(mfrow = c(1, 2))
    on.exit(graphics::par(old_par), add = TRUE)
    metafor()
    brma()
  })
}

