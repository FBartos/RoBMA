# Native interface helpers ----------------------------------------------------

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
