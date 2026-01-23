---
description: R code style guidelines for the RoBMA package
---

# R Code Style Guidelines

When writing or modifying R code in this package, follow these conventions:

## Assignment Arrow Alignment

When there are multiple related assignments in a block, align the assignment arrows (`<-`) for readability:

### Correct
```r
is_multilevel     <- .is_multilevel(x)
is_weightfunction <- .is_weightfunction(x)
is_PET            <- .is_PET(x)
is_PEESE          <- .is_PEESE(x)
effect_direction  <- .effect_direction(x)
```

### Incorrect
```r
is_multilevel <- .is_multilevel(x)
is_weightfunction <- .is_weightfunction(x)
is_PET <- .is_PET(x)
is_PEESE <- .is_PEESE(x)
effect_direction <- .effect_direction(x)
```

## Function Arguments Alignment

When a function call spans multiple lines, align arguments for clarity:

```r
ci_bounds <- .get_funnel_quantiles(
  x                 = x,
  se_sequence       = se_sequence,
  tau               = tau,
  sampling_bias     = sampling_bias,
  is_weightfunction = is_weightfunction,
  is_PET            = is_PET,
  is_PEESE          = is_PEESE,
  effect_direction  = effect_direction
)
```

## Return Statement Alignment

Align elements in return lists:

```r
return(list(
  points       = df_points,
  funnel       = df_funnel,
  funnel_edge1 = df_funnel_edge1,
  funnel_edge2 = df_funnel_edge2,
  background   = df_background,
  x_range      = x_range,
  y_range      = se_range,
  xlab         = xlab,
  ylab         = ylab
))
```

## General Spacing

- Use spaces around assignment operators: `x <- 1` not `x<-1`
- Use spaces after commas: `c(1, 2, 3)` not `c(1,2,3)`
- Use 2 spaces for indentation
- Add an empty line after the opening brace `{` of a function:
  ```r
  my_function <- function(x) {
    
    # code starts after an empty line
    return(x)
  }
  ```
