# version 1.0.4
Fixes:
- inability to run models without the silent = TRUE control

# version 1.0.3
Features:
- x-axis rescaling for the weight function plot (by setting 'rescale_x = TRUE' in the 'plot.RoBMA' function)
- setting expected direction of the effect in for RoBMA function

Fixes:
- marginal likelihood calculation for models with spike prior distribution on mean parameter which location was not set to 0
- some additional error messages 

CRAM requested changes:
- changing information messages from 'cat' to 'message' from plot related functions
- saving and returning the 'par' settings to the user defined one in the base plot functions

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
