# Interprets results of a RoBMA model.

`interpret` creates a brief textual summary of a fitted RoBMA object.

## Usage

``` r
interpret(object, output_scale = NULL)
```

## Arguments

- object:

  a fitted RoBMA object

- output_scale:

  transform the meta-analytic estimates to a different scale. Defaults
  to `NULL` which returns the same scale as the model was estimated on.

## Value

`interpret` returns a character.
