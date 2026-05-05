# R Package Plotting Function Architecture

**Domain**: R package development, visualization
**Confidence**: High
**Source**: RoBMA package patterns (funnel.R, radial.R, regplot.R)

## Pattern

When creating plotting functions for R packages that need both base R and ggplot2 support, use a **3-layer architecture**:

### Layer 1: Main Dispatcher Function

```r
#' @export
myplot <- function(x, ...) UseMethod("myplot")

#' @export
myplot.myclass <- function(x, plot_type = "base", as_data = FALSE, ...) {

  # Input validation

BayesTools::check_char(plot_type, "plot_type", allow_values = c("base", "ggplot"))
  BayesTools::check_bool(as_data, "as_data")

  # Set up graphical parameters with defaults
  dots <- .set_dots_myplot(...)

  # Generate plot data (pure computation)
  plot_data <- .myplot_data(x = x, dots = dots)

  # Return data for programmatic access
  if (isTRUE(as_data)) {
    return(plot_data)
  }

  # Dispatch to rendering function
  if (plot_type == "ggplot") {
    return(.myplot_plot_ggplot(plot_data, dots))
  } else {
    .myplot_plot_base(plot_data, dots)
    return(invisible(NULL))
  }
}
```

### Layer 2: Data Generation Function

Pure computation, no graphics code. Returns structured list:

```r
.myplot_data <- function(x, dots) {

  # Extract data from object
  yi  <- .outcome_data_yi(x)
  sei <- .outcome_data_sei(x)

  # Compute all necessary coordinates
  # ...computation logic...

  # Return structured list
  return(list(
    points     = data.frame(x = ..., y = ...),
    line       = data.frame(x = ..., y = ...),
    polygon    = data.frame(x = ..., y = ...),
    xlim       = xlim,
    ylim       = ylim,
    xlab       = xlab,
    ylab       = ylab
  ))
}
```

### Layer 3: Rendering Functions

Separate functions for each backend:

```r
# Base R graphics
.myplot_plot_base <- function(data, dots) {

  # Set up empty plot
  graphics::plot(NA, NA,
    xlim = data$xlim, ylim = data$ylim,
    xlab = data$xlab, ylab = data$ylab,
    type = "n"
  )

  # Draw elements in order (back to front)
  graphics::polygon(data$polygon$x, data$polygon$y, col = dots$shade)
  graphics::lines(data$line$x, data$line$y, col = dots$lcol)
  graphics::points(data$points$x, data$points$y, pch = dots$pch)

  return(invisible(NULL))
}

# ggplot2 graphics
.myplot_plot_ggplot <- function(data, dots) {

  out <- ggplot2::ggplot() +
    ggplot2::geom_polygon(
      mapping = ggplot2::aes(x = data$polygon$x, y = data$polygon$y),
      fill = dots$shade
    ) +
    ggplot2::geom_line(
      mapping = ggplot2::aes(x = data$line$x, y = data$line$y),
      colour = dots$lcol
    ) +
    ggplot2::geom_point(
      mapping = ggplot2::aes(x = data$points$x, y = data$points$y),
      shape = dots$pch
    ) +
    ggplot2::scale_x_continuous(limits = data$xlim, name = data$xlab) +
    ggplot2::scale_y_continuous(limits = data$ylim, name = data$ylab)

  return(out)
}
```

### Defaults Helper Function

```r
.set_dots_myplot <- function(...) {

  dots <- list(...)

  # Point styling
  if (is.null(dots[["pch"]]))   dots[["pch"]]   <- 21
  if (is.null(dots[["col"]]))   dots[["col"]]   <- "black"
  if (is.null(dots[["bg"]]))    dots[["bg"]]    <- "gray"
  if (is.null(dots[["cex"]]))   dots[["cex"]]   <- 1
  if (is.null(dots[["size"]]))  dots[["size"]]  <- 2  # ggplot equivalent

  # Line styling
  if (is.null(dots[["lcol"]]))  dots[["lcol"]]  <- "black"
  if (is.null(dots[["lwd"]]))   dots[["lwd"]]   <- 2

  # Region styling
  if (is.null(dots[["shade"]])) dots[["shade"]] <- "gray90"

  # Title
  if (is.null(dots[["main"]]))  dots[["main"]]  <- NULL

  return(dots)
}
```

## Benefits

1. **Separation of concerns** - computation vs rendering cleanly separated
2. **Testability** - `as_data = TRUE` enables programmatic testing of computed values
3. **Maintainability** - changes to one backend don't affect the other
4. **Consistency** - both backends produce identical plots from same data
5. **Flexibility** - users can choose their preferred graphics system

## When to Use

- Creating visualization functions for R packages
- Need to support both base R and ggplot2 users
- Complex plots with multiple elements (polygons, lines, points)
- Plots that may need programmatic data access

## Examples in RoBMA

- `funnel.brma()` - funnel plots
- `radial.brma()` - radial/Galbraith plots
- `regplot.brma()` - regression/bubble plots
