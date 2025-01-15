##' RoBMA: Robust Bayesian meta-analysis
##'
##' RoBMA: Bayesian model-averaged meta-analysis with adjustments for publication
##' bias and ability to specify informed prior distributions and draw inference with
##' inclusion Bayes factors.
##'
##'
##' @name RoBMA-package
##' @author František Bartoš \email{f.bartos96@@gmail.com}
##' @keywords package
##' @aliases RoBMA-package RoBMA_package RoBMA.package
##' @docType package
##' @section
##' User guide: See \insertCite{bartos2021no;textual}{RoBMA},
##' \insertCite{maier2020robust;textual}{RoBMA}, and
##' \insertCite{bartos2020adjusting;textual}{RoBMA} for details regarding the RoBMA
##' methodology.
##'
##' More details regarding customization of the model ensembles are provided in the
##' \href{../doc/ReproducingBMA.html}{\bold{Reproducing BMA}},
##' \href{../doc/MedicineBMA.html}{\bold{BMA in Medicine}}, and
##' \href{../doc/CustomEnsembles.html}{\bold{Fitting Custom Meta-Analytic Ensembles}}
##' vignettes. Please, use the "Issues" section in the GitHub repository to ask any
##' further questions.
##'
##' @references \insertAllCited{}
##' @importFrom BayesTools is.prior is.prior.none is.prior.point is.prior.simple is.prior.factor is.prior.PET is.prior.PEESE is.prior.weightfunction
##' @importFrom BayesTools is.prior.independent is.prior.spike_and_slab is.prior.mixture
##' @importFrom Rdpack reprompt
##' @importFrom rlang .data
"_PACKAGE"

