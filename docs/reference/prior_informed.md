# Informed Prior

Create empirical informed prior distributions.

## Usage

``` r
prior_informed(name, parameter = NULL, type = "smd")
```

## Arguments

- name:

  character. Informed prior name. Supported examples include
  `"van Erp"`, `"Oosterwijk"`, `"Cochrane"`, and medicine subfields from
  [`BayesTools::prior_informed_medicine_names`](https://fbartos.github.io/BayesTools/reference/prior_informed_medicine_names.html).

- parameter:

  character. Parameter subset. Required for medicine priors and not
  needed for psychology priors.

- type:

  character. Effect-size type. Medicine priors use values such as
  `"smd"`, `"logOR"`, `"logRR"`, `"RD"`, and `"logHR"`.

## Value

An object of class `prior`.

## Details

Further details can be found in van Erp et al. (2017) , Gronau et al.
(2017) , Bartoš et al. (2021) , and Bartoš et al. (2023) .

## References

Bartoš F, Gronau QF, Timmers B, Otte WM, Ly A, Wagenmakers E (2021).
“Bayesian model-averaged meta-analysis in medicine.” *Statistics in
Medicine*, **40**(30), 6743–6761.
[doi:10.1002/sim.9170](https://doi.org/10.1002/sim.9170) .\
\
Bartoš F, Otte WM, Gronau QF, Timmers B, Ly A, Wagenmakers E (2023).
“Empirical prior distributions for Bayesian meta-analyses of binary and
time-to-event outcomes.”
[doi:10.48550/arXiv.2306.11468](https://doi.org/10.48550/arXiv.2306.11468)
. Preprint available at https://doi.org/10.48550/arXiv.2306.11468.\
\
Gronau QF, Van Erp S, Heck DW, Cesario J, Jonas KJ, Wagenmakers E
(2017). “A Bayesian model-averaged meta-analysis of the power pose
effect with informed and default priors: The case of felt power.”
*Comprehensive Results in Social Psychology*, **2**(1), 123–138.
[doi:10.1080/23743603.2017.1326760](https://doi.org/10.1080/23743603.2017.1326760)
.\
\
van Erp S, Verhagen J, Grasman RP, Wagenmakers E (2017). “Estimates of
between-study heterogeneity for 705 meta-analyses reported in
Psychological Bulletin from 1990–2013.” *Journal of Open Psychology
Data*, **5**(1), 1–5.
[doi:10.5334/jopd.33](https://doi.org/10.5334/jopd.33) .
