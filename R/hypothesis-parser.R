.hypothesis_brma_point_refs <- function(hypothesis, parameter) {

  rows <- list()
  row_i <- 1L

  for (text in hypothesis) {
    text <- .hypothesis_brma_normalize_level_references(text)
    sides <- unlist(strsplit(text, "\\s+[Vv][Ss]\\s+", perl = TRUE),
                    use.names = FALSE)

    for (side in sides) {
      if (grepl("[&|]", side, perl = TRUE)) {
        next
      }

      relation <- .hypothesis_brma_find_relation(side)
      if (is.null(relation) ||
          !relation[["operator"]] %in% c("=", "==", "!=")) {
        next
      }

      lhs <- trimws(substr(side, 1L, relation[["start"]] - 1L))
      rhs <- trimws(substr(side, relation[["end"]] + 1L, nchar(side)))
      symbol <- .hypothesis_brma_direct_symbol(lhs)
      if (is.null(symbol)) {
        stop(
          "qCMDE/IWMDE precomputed ordinates require direct parameter or ",
          "level point hypotheses; unsupported point expression: '",
          lhs, "'.",
          call. = FALSE
        )
      }

      value <- .hypothesis_brma_numeric_rhs(rhs)
      ref <- .hypothesis_brma_level_ref(symbol)
      if (is.null(ref) && identical(symbol, parameter)) {
        rows[[row_i]] <- data.frame(
          symbol       = symbol,
          level        = NA_character_,
          value        = value,
          stringsAsFactors = FALSE
        )
        row_i <- row_i + 1L
      } else if (!is.null(ref) && identical(ref[["parameter"]], parameter)) {
        rows[[row_i]] <- data.frame(
          symbol       = symbol,
          level        = ref[["level"]],
          value        = value,
          stringsAsFactors = FALSE
        )
        row_i <- row_i + 1L
      }
    }
  }

  if (length(rows) == 0L) {
    return(data.frame(
      symbol       = character(),
      level        = character(),
      value        = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  out <- unique(do.call(rbind, rows))
  rownames(out) <- NULL
  return(out)
}


.hypothesis_brma_find_relation <- function(text) {

  chars       <- strsplit(text, "", fixed = TRUE)[[1L]]
  in_backtick <- FALSE
  depth       <- 0L
  found       <- list()
  two_char    <- c("<=", ">=", "==", "!=")
  one_char    <- c("<", ">", "=")

  i <- 1L
  while (i <= length(chars)) {
    ch <- chars[[i]]
    if (ch == "`") {
      in_backtick <- !in_backtick
      i <- i + 1L
      next
    }
    if (!in_backtick) {
      if (ch == "(") {
        depth <- depth + 1L
      } else if (ch == ")") {
        depth <- max(0L, depth - 1L)
      } else if (depth == 0L) {
        op <- NULL
        if (i < length(chars)) {
          candidate <- paste0(chars[[i]], chars[[i + 1L]])
          if (candidate %in% two_char) {
            op <- candidate
          }
        }
        if (is.null(op) && ch %in% one_char) {
          op <- ch
        }
        if (!is.null(op)) {
          found[[length(found) + 1L]] <- list(
            operator = op,
            start    = i,
            end      = i + nchar(op) - 1L
          )
          i <- i + nchar(op)
          next
        }
      }
    }
    i <- i + 1L
  }

  if (length(found) == 0L) {
    return(NULL)
  }
  if (length(found) > 1L) {
    return(NULL)
  }

  return(found[[1L]])
}


.hypothesis_brma_direct_symbol <- function(text) {

  parsed <- try(parse(text = text, keep.source = FALSE), silent = TRUE)
  if (inherits(parsed, "try-error") || length(parsed) != 1L) {
    return(NULL)
  }
  expr <- parsed[[1L]]
  if (!is.name(expr)) {
    return(NULL)
  }

  return(as.character(expr))
}


.hypothesis_brma_numeric_rhs <- function(text) {

  parsed <- try(parse(text = text, keep.source = FALSE), silent = TRUE)
  if (inherits(parsed, "try-error") || length(parsed) != 1L) {
    stop("Right side of a point hypothesis must be numeric.", call. = FALSE)
  }
  symbols <- setdiff(
    all.names(parsed[[1L]], functions = TRUE, unique = TRUE),
    c("(", "+", "-", "*", "/", "^", "Inf", "NaN", "NA", "TRUE", "FALSE")
  )
  if (length(symbols) > 0L) {
    stop("Right side of a point hypothesis must be numeric.", call. = FALSE)
  }

  value <- eval(parsed[[1L]], envir = baseenv())
  BayesTools::check_real(value, "hypothesis value", check_length = 1, allow_NA = FALSE)
  if (!is.finite(value)) {
    stop("Hypothesis value must be finite.", call. = FALSE)
  }

  return(value)
}
