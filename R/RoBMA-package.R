##' RoBMA: Robust Bayesian Meta-Analysis
##'
##' RoBMA provides Bayesian meta-analysis, meta-regression, multilevel
##' meta-analysis, model averaging, and publication-bias adjustment. The main
##' user-facing fitters are \code{\link{RoBMA}}, \code{\link{BMA}},
##' \code{\link{brma}}, \code{\link{brma.glmm}}, \code{\link{bselmodel}},
##' \code{\link{bPET}}, and \code{\link{bPEESE}}.
##'
##' @name RoBMA-package
##' @author Frantisek Bartos \email{f.bartos96@@gmail.com}
##' @keywords package
##' @aliases RoBMA-package RoBMA_package RoBMA.package
##' @docType package
##' @section User guide:
##' See \insertCite{bartos2021no;textual}{RoBMA},
##' \insertCite{maier2020robust;textual}{RoBMA},
##' \insertCite{bartos2020adjusting;textual}{RoBMA}, and
##' \insertCite{bartos2025robust;textual}{RoBMA} for the RoBMA methodology.
##' Use \code{vignette(package = "RoBMA")} to list installed vignettes.
##'
##' @references \insertAllCited{}
##' @importFrom BayesTools is.prior is.prior.none is.prior.point is.prior.simple is.prior.factor is.prior.PET is.prior.PEESE is.prior.weightfunction
##' @importFrom BayesTools is.prior.independent is.prior.spike_and_slab is.prior.mixture
##' @importFrom Rdpack reprompt
##' @importFrom rlang .data
"_PACKAGE"
