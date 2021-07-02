##' RoBMA: Robust Bayesian meta-analysis
##'
##' RoBMA: Bayesian model-averaged meta-analysis with adjustments for publication
##' bias and ability to specify informed prior distributions and draw inference with
##' inclusion Bayes factors.
##'
##' \code{\link{RoBMA}} fits parametric models ....
##'
##'
##' @name RoBMA-package
##' @aliases RoBMA-package
##' @docType package
##' @section User guide: The \bold{RoBMA user guide} vignette explains the
##' methods in detail, and gives several worked examples.  A further vignette
##' \bold{RoBMA-examples} gives a few more complicated examples, and users
##' are encouraged to submit their own.
##' @author František Bartoš \email{f.bartos96@@gmail.com}
##' @references
##' @keywords package
##' @importFrom BayesTools prior prior_none prior_PET prior_PEESE prior_weightfunction
##' @importFrom BayesTools is.prior is.prior.none is.prior.point is.prior.simple is.prior.PET is.prior.PEESE is.prior.weightfunction
##' @export prior
##' @export prior_none
##' @export prior_PET
##' @export prior_PEESE
##' @export prior_weightfunction
"_PACKAGE"

#//##' @useDynLib RoBMA, .registration = TRUE
#//.onUnload <- function(libpath) {
#//  library.dynam.unload("RoBMA", libpath)
#//}

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL

