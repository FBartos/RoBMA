% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check-input-and-settings.R
\name{RoBMA_control}
\alias{RoBMA_control}
\alias{set_autofit_control,}
\alias{set_convergence_checks}
\alias{set_autofit_control}
\title{Control MCMC fitting process}
\usage{
set_autofit_control(
  max_Rhat = 1.05,
  min_ESS = 500,
  max_error = NULL,
  max_SD_error = NULL,
  max_time = list(time = 60, unit = "mins"),
  sample_extend = 1000,
  restarts = 10,
  max_extend = 10
)

set_convergence_checks(
  max_Rhat = 1.05,
  min_ESS = 500,
  max_error = NULL,
  max_SD_error = NULL,
  remove_failed = FALSE,
  balance_probability = TRUE
)
}
\arguments{
\item{max_Rhat}{maximum value of the R-hat diagnostic.
Defaults to \code{1.05}.}

\item{min_ESS}{minimum estimated sample size.
Defaults to \code{500}.}

\item{max_error}{maximum value of the MCMC error.
Defaults to \code{NULL}. Be aware that PEESE publication bias
adjustment can have estimates on different scale than the rest of
the output, resulting in relatively large max MCMC error.}

\item{max_SD_error}{maximum value of the proportion of MCMC error
of the estimated SD of the parameter.
Defaults to \code{NULL}.}

\item{max_time}{list with the time and unit specifying the maximum
autofitting process per model. Passed to \link[base]{difftime} function
(possible units are \code{"secs", "mins", "hours", "days", "weeks", "years"}).
Defaults to \code{list(time = 60, unit = "mins")}.}

\item{sample_extend}{number of samples to extend the fitting process if
the criteria are not satisfied.
Defaults to \code{1000}.}

\item{restarts}{number of times new initial values should be generated in case a
model fails to initialize. Defaults to \code{10}.}

\item{max_extend}{number of times after which the automatic fitting function is stopped.}

\item{remove_failed}{whether models not satisfying the convergence checks should
be removed from the inference. Defaults to \code{FALSE} - only a warning is raised.}

\item{balance_probability}{whether prior model probability should be balanced
across the combinations of models with the same H0/H1 for effect / heterogeneity / bias
in the case of non-convergence. Defaults to \code{TRUE}.}
}
\value{
\code{set_autofit_control} returns a list of autofit control settings
and \code{set_convergence_checks} returns a list of convergence checks settings.
}
\description{
Controls settings for the autofit
process of the MCMC JAGS sampler (specifies termination
criteria), and values for the convergence checks.
}
\seealso{
\link{RoBMA}, \link{update.RoBMA}
}
