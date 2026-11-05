# Interpret brma Results

Creates a concise textual interpretation of fitted RoBMA `brma` objects.

## Usage

``` r
interpret(object, ...)

# Default S3 method
interpret(object, ...)

# S3 method for class 'brma'
interpret(
  object,
  output_measure = NULL,
  transform = NULL,
  conditional = FALSE,
  scope = "core",
  probs = c(0.025, 0.975),
  central = NULL,
  priors = FALSE,
  digits = 3,
  ...
)

# S3 method for class 'interpret.brma'
print(x, ...)
```

## Arguments

- object:

  a fitted model object.

- ...:

  additional arguments passed to methods. The `brma` method reserves
  `...`; unused arguments error, except deprecated `output_scale`, which
  is accepted with a warning.

- output_measure:

  optional effect-size measure used for the pooled effect only.
  Supported conversion targets include `"SMD"`, `"COR"`, `"ZCOR"`, and
  `"OR"` when the input scale allows conversion.

- transform:

  optional display transformation. `"EXP"` exponentiates log-scale
  `"OR"`, `"RR"`, `"HR"`, and `"IRR"` pooled effects; aliases for no
  transform are accepted.

- conditional:

  whether to summarize conditional estimates for RoBMA product-space
  objects. Defaults to `FALSE`.

- scope:

  character vector specifying sections to include. Use `"core"` for the
  default concise interpretation, `"all"` for all sections, or any of
  `"components"`, `"estimates"`, `"moderators"`, `"scale"`, and
  `"bias"`.

- probs:

  two quantiles used for credible intervals. Defaults to
  `c(.025, .975)`.

- central:

  whether estimates are described by posterior mean or median. Defaults
  to `NULL`, which uses the posterior mean, except for
  `transform = "EXP"` pooled effects where the posterior mode is used.

- priors:

  whether to print prior distributions after the interpretation.
  Defaults to `FALSE`.

- digits:

  number of digits after the decimal point.

- x:

  an `interpret.brma` object.

## Value

A character vector with class `"interpret.brma"`. The normalized
BayesTools interpretation records are stored in the `"records"`
attribute; the `brma` method also stores `"scope"`, `"conditional"`, and
optionally `"priors"`.
