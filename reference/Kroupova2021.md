# 881 estimates from 69 studies of a relationship between employment and educational outcomes collected by Kroupova et al. (2024)

The data set contains partial correlation coefficients, standard errors,
study labels, sample sizes, type of the educational outcome, intensity
of the employment, gender of the student population, study location,
study design, whether the study controlled for endogeneity, and whether
the study controlled for motivation. The original data set including
additional variables and the publication are available from the project
page. (Note that some standard errors and employment intensities are
missing.)

## Usage

``` r
Kroupova2021
```

## Format

A data.frame with 11 columns and 881 observations:

- `r`:

  Partial correlation coefficient.

- `se`:

  Standard error of r.

- `study`:

  Study label.

- `sample_size`:

  Sample size.

- `education_outcome`:

  Type of educational outcome.

- `employment_intensity`:

  Employment intensity.

- `students_gender`:

  Gender composition of the student sample.

- `location`:

  Study location.

- `design`:

  Study design.

- `endogenity_control`:

  Whether endogeneity was controlled.

- `motivation_control`:

  Whether motivation was controlled.

## Source

<https://meta-analysis.cz/students/>

## References

Kroupova K, Havranek T, Irsova Z (2024). “Student employment and
education: A meta-analysis.” *Economics of Education Review*, **100**,
102539.
[doi:10.1016/j.econedurev.2024.102539](https://doi.org/10.1016/j.econedurev.2024.102539)
.
