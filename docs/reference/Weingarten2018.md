# 582 effect sizes examining the ease-of-retrieval effect from a meta-analysis by Weingarten and Hutchinson (2018)

The data set contains correlation coefficients between the manipulation
and outcome variable (r_xy), sample sizes (N), and various study
characteristics including publication status, country (USA vs. other),
number of few and many examples requested, whether the trial targeted
episodic memory, paradigm type (standard vs. other), dataset type
(proximal vs. distal), and mediation variables (r_xm, r_my). The
meta-analysis examined whether subjective ease mediates the
ease-of-retrieval effect, where participants list either few or many
examples and then make judgments (Weingarten and Hutchinson 2018) .

## Usage

``` r
Weingarten2018
```

## Format

A data.frame with 12 columns and 582 observations:

- `r_xy`:

  Correlation between manipulation and outcome.

- `N`:

  Sample size.

- `paper_id`:

  Paper identifier.

- `published`:

  Publication-status indicator.

- `USA`:

  Whether the study was conducted in the USA.

- `number_of_few`:

  Number of examples in the few-examples condition.

- `number_of_many`:

  Number of examples in the many-examples condition.

- `episodic_memory`:

  Whether the trial targeted episodic memory.

- `standard_paradigm`:

  Whether the standard paradigm was used.

- `proximal_dataset`:

  Whether the data set was proximal.

- `r_xm`:

  Correlation between manipulation and mediator.

- `r_my`:

  Correlation between mediator and outcome.

## References

Weingarten E, Hutchinson JW (2018). “Does ease mediate the
ease-of-retrieval effect? A meta-analysis.” *Psychological Bulletin*,
**144**(3), 227–283.
[doi:10.1037/bul0000122](https://doi.org/10.1037/bul0000122) .
