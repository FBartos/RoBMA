# version 1.0.2
Fixes:
- the summary and plot function now shows quantile based confidence intervals for individual models instead of the HPD provided before (this affects only 'summary'/'plot' with 'type = "individual"', all other confidence intervals were quantile based before)

# version 1.0.1
Fixes:
- summary function returning median instead of mean

# version 1.0.0 (vs the osf version)
Fixes:
- incorrectly weighted theta estimates
- models with non-zero point prior distribution incorrectly plotted using when "models" option in case that the mu parameter was transformed

Additional features:
- analyzing OR
- distributions implemented using boost library (helps with convergence issues)
- ability to mute the non-suppressible "precision not achieved" warning messages by using "silent" = TRUE inside of the control argument
- vignettes

Notable changes:
- the way how the seed is set before model fitting (the simulation study will not be reproducible with the new version of the package)
