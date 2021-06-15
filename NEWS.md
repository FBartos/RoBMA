## version 1.3.0
### Features
- adding PET and PEESE style models for adjusting for publication bias
- adding functions for effect sizes and standard error transformations
- the prior distributions can be specified on a different scale ('prior_scale') than is used for fitting the effect sizes ('transform'). Defaults to Cohen's d and Fisher's z respectively. (The estimates are transformed back to a scale with majority of effect size estimates.)
- additional prior distribution: lognormal
- transforming the results between different effect sizes with the print.RoBMA, summary.RoBMA, and plot.RoBMA functions

### Changes to the RoBMA function
- the default ensemble now uses 36 models - the publication bias adjustment part was extended to 6 weight functions, PET, and PEESE style models. It corresponds to RoBMA 6wPP model introduce in out new [preprint](ADD LINK)

### Changes to the data input
- studies with different effect size measures can be supplied simultaneously - RoBMA will internally transform them into a common scale ('transform'). The common scale defaults to Fisher's z since it stabilizes variances. Not stabilizing variances leads to overestimating the evidence for publication bias in case that that PET and PEESE style models are used.
- one-sample t-tests are no longer supported as input - the corresponding effect sizes and standard errors must be manually computed and passed as 'd' and 'se'.
- see metafor::escalc for more details about calculating effect sizes and standard errors


## version 1.2.1
### Fixes
- check_setup function not working at all

## version 1.2.0
### Changes
- the studies's true effects are now marginalized out of the random effects models and are no longer estimated (see Appendix A of our [prerint](https://psyarxiv.com/u4cns/) for more details). As a results, arguments referring to the true effects are now disabled.
- all models are now being estimated using the likelihood of effect sizes (instead of test-statistics as usually defined). We reproduced the simulation study that we used to evaluate the method performance and it achieved identical results (up to MCMC error, before marginalizing out the true effects). A big advantage of using the normal likelihood for effect sizes is a considerable speed up of the whole estimation process.
- as a results of these two changes, the results of the models will differ to those of pre 1.2.0 version

### Fixes
- autofit being turn on if any control argument was specified


## version 1.1.2
### Fixes
- vdiffr not being used conditionally in unit tests

## version 1.1.1
### Fixes
- inability to fit a model without specifying a seed
- inability to produce individual model plots due to incompatibility with the newer versions of ggplot2  

## version 1.1.0
### Features
- parallel within and between model fitting using the parallel package with 'parallel = TRUE' argument


## version 1.0.5
### Fixes:
- models being fitted automatically until reaching R-hat lower than 1.05 without setting max_rhat and autofit control parameters
- bug preventing to draw a bivariate plot of mu and tau
- range for parameter estimates from individual models no containing 0 (or 1 in case of OR measured effect sizes)
- inability to fit a model with only null mu distributions if correlation or OR measured effect sizes were specified
- ordering of the estimated and observed effects when both of them are requested simultaneously
- formatting of this file (NEWS.md)

### Improvements:
- priors plot: parameter specification, default plotting range, clearer x-axis labels in cases when the parameter is defined on transformed scale
- parameters plots: probability scale always ends at the same spot as is the last tick on the density scale
- adding warnings if any of the specified models has Rhat higher than 1.05 or the specified value
- grouping the same warnings messages together

## version 1.0.4
### Fixes:
- inability to run models without the silent = TRUE control

## version 1.0.3
### Features:
- x-axis rescaling for the weight function plot (by setting 'rescale_x = TRUE' in the 'plot.RoBMA' function)
- setting expected direction of the effect in for RoBMA function

### Fixes:
- marginal likelihood calculation for models with spike prior distribution on mean parameter which location was not set to 0
- some additional error messages 

### CRAM requested changes:
- changing information messages from 'cat' to 'message' from plot related functions
- saving and returning the 'par' settings to the user defined one in the base plot functions

## version 1.0.2
### Fixes:
- the summary and plot function now shows quantile based confidence intervals for individual models instead of the HPD provided before (this affects only 'summary'/'plot' with 'type = "individual"', all other confidence intervals were quantile based before)

## version 1.0.1
### Fixes:
- summary function returning median instead of mean

## version 1.0.0 (vs the osf version)
### Fixes:
- incorrectly weighted theta estimates
- models with non-zero point prior distribution incorrectly plotted using when "models" option in case that the mu parameter was transformed

### Additional features:
- analyzing OR
- distributions implemented using boost library (helps with convergence issues)
- ability to mute the non-suppressible "precision not achieved" warning messages by using "silent" = TRUE inside of the control argument
- vignettes

### Notable changes:
- the way how the seed is set before model fitting (the simulation study will not be reproducible with the new version of the package)
