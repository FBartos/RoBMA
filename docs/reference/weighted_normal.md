# Weighted normal distribution

Density, distribution function, quantile function and random generation
for the weighted normal distribution with `mean`, standard deviation
`sd`, steps `steps` (or critical values) `crit_x`), and weights `omega`.

## Usage

``` r
dwnorm(
  x,
  mean,
  sd,
  steps = if (!is.null(crit_x)) NULL,
  omega,
  crit_x = if (!is.null(steps)) NULL,
  type = "two.sided",
  log = FALSE
)

pwnorm(
  q,
  mean,
  sd,
  steps = if (!is.null(crit_x)) NULL,
  omega,
  crit_x = if (!is.null(steps)) NULL,
  type = "two.sided",
  lower.tail = TRUE,
  log.p = FALSE
)

qwnorm(
  p,
  mean,
  sd,
  steps = if (!is.null(crit_x)) NULL,
  omega,
  crit_x = if (!is.null(steps)) NULL,
  type = "two.sided",
  lower.tail = TRUE,
  log.p = FALSE
)

rwnorm(
  n,
  mean,
  sd,
  steps = if (!is.null(crit_x)) NULL,
  omega,
  crit_x = if (!is.null(steps)) NULL,
  type = "two.sided"
)
```

## Arguments

- x, q:

  vector of quantiles.

- mean:

  mean

- sd:

  standard deviation.

- steps:

  vector of steps for the weight function.

- omega:

  vector of weights defining the probability of observing a t-statistics
  between each of the two steps.

- crit_x:

  vector of critical values defining steps (if `steps` are not
  supplied).

- type:

  type of weight function (defaults to `"two.sided"`).

- log, log.p:

  logical; if `TRUE`, probabilities `p` are given as `log(p)`.

- lower.tail:

  logical; if `TRUE` (default), probabilities are \\P\[X \le x\]\\,
  otherwise, \\P\[X \ge x\]\\.

- p:

  vector of probabilities.

- n:

  number of observations. If length(n) \> 1, the length is taken to be
  the number required.

## Value

`dwnorm` gives the density, `dwnorm` gives the distribution function,
`qwnorm` gives the quantile function, and `rwnorm` generates random
deviates.

## Details

The `mean`, `sd`, `steps`, `omega` can be supplied as a vectors (`mean`,
`sd`) or matrices (`steps`, `omega`) with length / number of rows equal
to `x`/`q`/ `p`. Otherwise, they are recycled to the length of the
result.

## See also

[Normal](https://rdrr.io/r/stats/Normal.html)
