# 106 analyst-team effect-size estimates of the relation between religiosity and well-being from Hoogeveen et al. (2023)

The data set contains standardized effect-size estimates, standard
errors, 95% confidence interval bounds, sample sizes, effect-size types,
sequential analyst-team identifiers, and team-level characteristics
describing psychological dependent-variable coding and self-reported
methodological knowledge. The Many-Analysts Religion Project examined
whether religious people self-report greater well-being (Hoogeveen et
al. 2023) .

## Usage

``` r
Hoogeveen2023
```

## Format

A data.frame with 9 columns and 106 observations:

- `yi`:

  Standardized effect-size estimate.

- `sei`:

  Standard error of the effect-size estimate.

- `lci`:

  Lower bound of the reported 95% confidence interval.

- `uci`:

  Upper bound of the reported 95% confidence interval.

- `ni`:

  Analysis sample size.

- `type`:

  Reported effect-size type: standardized regression coefficient beta or
  Cohen's d.

- `team`:

  Sequential analyst-team identifier in this example.

- `team_psych`:

  Whether the effect estimate used the psychological dependent-variable
  coding.

- `team_knowledge`:

  Team's self-reported methodological-knowledge rating, from 1 (no
  knowledge) to 5 (expert).

## Source

<https://osf.io/gkxqy/>

## References

Hoogeveen S, Sarafoglou A, Aczel B, Aditya Y, Alayan AJ, others (2023).
“A many-analysts approach to the relation between religiosity and
well-being.” *Religion, Brain & Behavior*, **13**(3), 237–283.
[doi:10.1080/2153599X.2022.2070255](https://doi.org/10.1080/2153599X.2022.2070255)
.
