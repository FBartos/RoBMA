# 39 study rows on household chaos and child executive functions from a meta-analysis by Andrews et al. (2021)

The data set contains raw correlation coefficients and sampling
information using metafor-compatible column names (`ri`, `vi`, `sei`,
`ni`, and `slab`), together with study-level moderators from the
original data set. Three source rows have missing effect-size inputs and
are omitted automatically by model-fitting functions that require
complete outcomes. The original data set assessed the effect of
household chaos on child executive functions (Andrews et al. 2021) and
was used as an example in Bartoš et al. (2025) .

## Usage

``` r
Andrews2021
```

## Format

A data.frame with 15 columns and 39 observations:

- `study_id`:

  Source row identifier.

- `slab`:

  Study label.

- `ri`:

  Raw correlation coefficient.

- `vi`:

  Sampling variance of the raw correlation coefficient.

- `sei`:

  Sampling standard error of the raw correlation coefficient.

- `ni`:

  Sample size.

- `year`:

  Publication year.

- `percent_female`:

  Percentage of female children in the sample.

- `age`:

  Mean age of the children in years.

- `assessment_interval_months`:

  Time between household-chaos and executive-function assessments in
  months.

- `dissertation`:

  Whether the study was a dissertation.

- `percent_minority`:

  Percentage of minority participants in the sample.

- `percent_low_parental_education`:

  Percentage of parents with high school, GED, or lower education.

- `hc_dimension`:

  Household-chaos dimension using the source coding.

- `measure`:

  Executive-function assessment type, direct or informant.

## References

Andrews K, Atkinson L, Harris M, Gonzalez A (2021). “Examining the
effects of household chaos on child executive functions: A
meta-analysis.” *Psychological Bulletin*, **147**(1), 16–32.
[doi:10.1037/bul0000311](https://doi.org/10.1037/bul0000311) .  
  
Bartoš F, Maier M, Stanley TD, Wagenmakers E (2025). “Robust Bayesian
meta-regression: Model-averaged moderation analysis in the presence of
publication bias.” *Psychological Methods*.
[doi:10.1037/met0000737](https://doi.org/10.1037/met0000737) .
