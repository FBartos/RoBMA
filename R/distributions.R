# Native interface helpers ----------------------------------------------------

# Compute sqrt(x^2 + y^2) without intermediate overflow or underflow.
.root_sum_squares <- function(x, y) {

  scale    <- pmax(abs(x), abs(y))
  positive <- !is.na(scale) & is.finite(scale) & scale > 0
  out      <- scale
  x_scaled <- x / scale
  y_scaled <- y / scale
  out[positive] <- scale[positive] * sqrt(
    x_scaled[positive]^2 + y_scaled[positive]^2
  )
  zero <- !is.na(scale) & scale == 0
  out[zero] <- 0
  return(out)
}

.native_numeric_matrix <- function(x) {

  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }
  if (!is.double(x)) {
    storage.mode(x) <- "double"
  }
  return(x)
}

.native_numeric_vector <- function(x) {

  if (!is.double(x)) {
    x <- as.double(x)
  }
  return(x)
}

.native_integer_vector <- function(x) {

  if (!is.integer(x)) {
    x <- as.integer(x)
  }
  return(x)
}
