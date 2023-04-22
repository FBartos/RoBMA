#' @title Combines different effect sizes into a common metric
#'
#' @description \code{combine_data} combines different effect sizes
#' into a common measure specified in \code{transformation}. Either
#' a data.frame \code{data} with columns named corresponding to the
#' arguments or vectors with individual values can be passed.
#'
#' @param d a vector of effect sizes measured as Cohen's d
#' @param r a vector of effect sizes measured as correlations
#' @param logOR a vector of effect sizes measured as log odds ratios
#' @param z a vector of effect sizes measured as Fisher's z
#' @param t a vector of t/z-statistics
#' @param y a vector of unspecified effect sizes (note that effect size
#' transformations are unavailable with this type of input)
#' @param se a vector of standard errors of the effect sizes
#' @param v a vector of variances of the effect sizes
#' @param n a vector of overall sample sizes
#' @param lCI a vector of lower bounds of confidence intervals
#' @param uCI a vector of upper bounds of confidence intervals
#' @param study_names an optional argument with the names of the studies
#' @param study_ids an optional argument specifying dependency between the
#' studies (for using a multilevel model). Defaults to \code{NULL} for
#' studies being independent.
#' @param data a data frame with column names corresponding to the
#' variable names used to supply data individually
#' @param transformation transformation to be applied to the supplied
#' effect sizes before fitting the individual models. Defaults to
#' \code{"fishers_z"}. We highly recommend using \code{"fishers_z"}
#' transformation since it is the only variance stabilizing measure
#' and does not bias PET and PEESE style models. The other options are
#' \code{"cohens_d"}, correlation coefficient \code{"r"} and \code{"logOR"}.
#' Supplying \code{"none"} will treat the effect sizes as unstandardized and
#' refrain from any transformations.
#' @param weight specifies likelihood weights of the individual estimates.
#' Notes that this is an untested experimental feature.
#' @param return_all whether data frame containing all filled values should be
#' returned. Defaults to \code{FALSE}
#'
#' @details The aim of the function is to combine different, already calculated,
#' effect size measures. In order to obtain effect size measures from raw values,
#' e.g, mean differences, standard deviations, and sample sizes, use
#' \link[metafor]{escalc} function.
#'
#' The function checks the input values and in transforming the input into a common
#' effect size measure in the following fashion:
#' 1) obtains missing standard errors by squaring variances
#' 2) obtains missing standard errors from confidence intervals (after transformation to
#' Fisher's z scale for \code{d} and \code{r}).
#' 3) obtains missing sample sizes (or standard errors for logOR) from t-statistics
#' and effect sizes
#' 4) obtains missing standard errors from sample sizes and effect sizes
#' 5) obtains missing sample sizes from standard errors and effect sizes
#' 6) obtains missing t-statistics from sample sizes and effect sizes
#' (or standard errors and effect sizes for logOR)
#' 7) changes the effect sizes direction to be positive
#' 8) transforms effect sizes into the common effect size
#' 9) transforms standard errors into the common metric
#'
#' If the \code{transforms} is \code{NULL} or an unstandardized effect size \code{y} is
#' supplied, steps 4-9 are skipped.
#'
#'
#' @return \code{combine_data} returns a data.frame.
#'
#' @seealso [RoBMA()], [check_setup()], [effect_sizes()], [standard_errors()], and [sample_sizes()]
#' @export
combine_data  <- function(d = NULL, r = NULL, z = NULL, logOR = NULL, t = NULL, y = NULL, se = NULL, v = NULL, n = NULL, lCI = NULL, uCI = NULL, study_names = NULL, study_ids = NULL, weight = NULL, data = NULL, transformation = "fishers_z", return_all = FALSE){

  # settings & input  check
  BayesTools::check_char(transformation, "transformation")
  BayesTools::check_bool(return_all, "return_all")
  BayesTools::check_real(d[!is.na(d)],           "d",      allow_NULL = TRUE, check_length = FALSE)
  BayesTools::check_real(r[!is.na(r)],           "r",      allow_NULL = TRUE, check_length = FALSE, lower = -1, upper = 1, allow_bound = FALSE)
  BayesTools::check_real(z[!is.na(z)],           "z",      allow_NULL = TRUE, check_length = FALSE)
  BayesTools::check_real(logOR[!is.na(logOR)],  "logOR",   allow_NULL = TRUE, check_length = FALSE, allow_bound = FALSE)
  BayesTools::check_real(t[!is.na(t)],           "t",      allow_NULL = TRUE, check_length = FALSE)
  BayesTools::check_real(y[!is.na(y)],           "y",      allow_NULL = TRUE, check_length = FALSE)
  BayesTools::check_real(se[!is.na(se)],         "se",     allow_NULL = TRUE, check_length = FALSE, lower = 0, allow_bound = FALSE)
  BayesTools::check_real(v[!is.na(v)],           "v",      allow_NULL = TRUE, check_length = FALSE, lower = 0, allow_bound = FALSE)
  BayesTools::check_int( n[!is.na(n)],           "n",      allow_NULL = TRUE, check_length = FALSE, lower = 0, allow_bound = FALSE)
  BayesTools::check_real(lCI[!is.na(lCI)],       "lCI",    allow_NULL = TRUE, check_length = FALSE)
  BayesTools::check_real(uCI[!is.na(uCI)],       "uCI",    allow_NULL = TRUE, check_length = FALSE)
  BayesTools::check_real(weight[!is.na(weight)], "weight", allow_NULL = TRUE, check_length = FALSE, lower = 0, allow_bound = FALSE)
  BayesTools::check_char(study_names[!is.na(study_names)], "study_names", allow_NULL = TRUE, check_length = FALSE)


  transformation <- .transformation_var(transformation)


  # forward information about the original measure when re-transformating the data
  if(inherits(data, "data.RoBMA")){
    colnames(data)[colnames(data) == "y"] <- attr(data, "effect_measure")
    original_measure                      <- attr(data, "original_measure")
  }else{
    original_measure <- NULL
  }

  input_variables <- c("d", "r", "z", "logOR", "y", "se", "v", "n", "lCI", "uCI", "t", "study_names", "study_ids", "weight")

  if(!is.null(data)){
    if(!is.data.frame(data))
      stop("Data must be passed as a data.frame.")
    if(any(!colnames(data) %in% input_variables))
      stop(paste0("The following variables do not correspond to any effect size/variability measure: ", paste(colnames(data)[!colnames(data) %in% input_variables], collapse = ", ")))
    data <- data[,colnames(data) %in% input_variables]
  }else{
    data <- data.frame(do.call(cbind, list(d = d, r = r, z = z, logOR = logOR, t = t, y = y, se = se, v = v, n = n, lCI = lCI, uCI = uCI, study_names = study_names, study_ids = study_ids, weight = weight)))
  }

  if(is.null(original_measure)){
    original_measure <- rep(NA, nrow(data))
  }

  ### add the remaining columns
  for(var in input_variables){
    if(!any(colnames(data) == var)){
      data <- cbind(data, NA)
      colnames(data)[ncol(data)] <- var
    }
  }

  ### into numeric
  for(var in c("d", "r", "z", "logOR", "y", "se", "v", "n", "lCI", "uCI", "t", "weight")){
    data[,var] <- as.numeric(as.character(data[,var]))
  }

  ### create holder of the output
  output <- data.frame(
    y  = rep(NA, nrow(data)),
    se = rep(NA, nrow(data)),
    study_names = rep(NA, nrow(data)),
    study_ids   = rep(NA, nrow(data)),
    weight      = rep(NA, nrow(data))
  )

  ### check for sufficient input
  if(all(is.na(data[, c("d", "r", "z", "logOR", "y")])))
    stop("The data do not contain any effect size measure.")
  if(all(is.na(data[, c("se", "v", "n", "lCI", "uCI", "t")])))
    stop("The data do not contain any variability measure.")

  # logical check for confidence intervals
  if(nrow(stats::na.omit(data[,c("lCI", "uCI")]) > 0)){
    if(any(stats::na.omit(data[,"lCI"]) > stats::na.omit(data[,"uCI"])))
      stop("'lCI' must be lower than 'uCI'.")
    for(var in c("d", "r", "z", "logOR", "y")){
      if(any(!is.na(data[,var]))){
        if((any(data[!is.na(data[,var]) & !is.na(data[,"lCI"]),var] < data[!is.na(data[,var]) & !is.na(data[,"lCI"]),"lCI"]) | any(data[!is.na(data[,var]) & !is.na(data[,"uCI"]),var] > data[!is.na(data[,var]) & !is.na(data[,"uCI"]),"uCI"])))
          stop("All effect sizes must be within the CI intervals.")
      }
    }
  }

  if(any(rowSums(!is.na(data[,c("d", "r", "z", "logOR", "y")])) > 1))
    stop("Only one effect size measure per study can be supplied.")

  if(any(rowSums(!is.na(data[,c("se", "v", "n", "t", "lCI")])) > 1))
    stop("Only one variability measure per study can be supplied.")

  if(any(rowSums(!is.na(data[,c("d", "r", "z", "logOR", "y")])) != 1))
    stop("At least one effect size measure per study must be supplied.")

  if(any(!rowSums(!is.na(data[,c("lCI", "uCI")])) %in% c(0,2)))
    stop("Either none or both CI bounds must be supplied.")

  if(any(rowSums(!is.na(data[,c("t", "se", "v", "n", "lCI")])) < 1))
    stop("At least one variability measure per study must be supplied.")


  ### store the original effect size measure
  original_measure[!is.na(data[,"d"])]     <- "d"
  original_measure[!is.na(data[,"r"])]     <- "r"
  original_measure[!is.na(data[,"z"])]     <- "z"
  original_measure[!is.na(data[,"logOR"])] <- "logOR"
  original_measure[!is.na(data[,"y"])]     <- "none"

  if(anyNA(original_measure))
    stop("At least one effect size measure per study must be supplied.")
  if(!(sum(original_measure == "none") == 0 | sum(original_measure == "none") == nrow(data)))
    stop("Standardized and general effect sizes cannot be combined.")

  # transform variance to standard errors
  data[is.na(data[,"se"]),"se"] <- sqrt(data[is.na(data[,"se"]),"v"])

  # add study names if missing
  if(all(is.na(data[,"study_names"]))){
    data[,"study_names"] <- paste0("Study ", 1:nrow(data))
  }

  # remove indicators from independent studies
  data[,"study_ids"][!data[,"study_ids"] %in% data[,"study_ids"][duplicated(data[,"study_ids"])]] <- NA
  # assign factor levels
  data[,"study_ids"] <- as.integer(as.factor(data[,"study_ids"]))

  # add weights if missing
  if(all(is.na(data[,"weight"]))){
    data[,"weight"] <- NA
  }

  ### deal with general 'unstandardized' input
  if(!anyNA(data[,"y"])){

    if(transformation != "y")
      stop("Effect sizes cannot be transformed in a presence of unstandardized measure ('y').")

    data[is.na(data[,"se"]),"se"] <- (data[is.na(data[,"se"]),"uCI"] - data[is.na(data[,"se"]),"lCI"]) / (2*stats::qnorm(.975))

    ### t-statistics from effect sizes and standard errors
    data[is.na(data[,"t"]),"t"] <- data[is.na(data[,"t"]),"y"] / data[is.na(data[,"t"]),"se"]

    ### add missing standard errors
    data[is.na(data[,"se"]),"se"] <- data[is.na(data[,"se"]),"y"] / data[is.na(data[,"se"]),"t"]

    if(return_all){
      data[,"lCI"] <- data[,"y"] + stats::qnorm(0.025) * data[,"se"]
      data[,"uCI"] <- data[,"y"] + stats::qnorm(0.975) * data[,"se"]
      return(data)
    }else{
      output$y  <- data[,"y"]
      output$se <- data[,"se"]
      output$study_names <- data[,"study_names"]
      output$study_ids   <- data[,"study_ids"]
      output$weight      <- data[,"weight"]
      attr(output, "effect_measure")   <- transformation
      attr(output, "original_measure") <- original_measure
      attr(output, "all_independent")  <- all(is.na(data[,"study_ids"]))
      attr(output, "weighted")         <- !all(is.na(data[,"weight"]))
      class(output) <- c(class(output), "data.RoBMA")

      return(output)
    }
  }

  ### calculate standard errors from CI using Fisher's z (for Cohen's d and r) and on original scale for log(OR) and Fisher's z
  data[is.na(data[,"se"]) & !.row_NA(data[,c("d", "lCI", "uCI")]), "se"] <- se_z2se_d(
    ( d2z(data[is.na(data[,"se"]) & !.row_NA(data[,c("d", "lCI", "uCI")]),"uCI"]) -
        d2z(data[is.na(data[,"se"]) & !.row_NA(data[,c("d", "lCI", "uCI")]),"lCI"]))/(2*stats::qnorm(.975)),
    d2z(data[is.na(data[,"se"]) & !.row_NA(data[,c("d", "lCI", "uCI")]), "d"])
  )
  data[is.na(data[,"se"]) & !.row_NA(data[,c("r", "lCI", "uCI")]), "se"] <- se_z2se_r(
    ( r2z(data[is.na(data[,"se"]) & !.row_NA(data[,c("r", "lCI", "uCI")]),"uCI"]) -
        r2z(data[is.na(data[,"se"]) & !.row_NA(data[,c("r", "lCI", "uCI")]),"lCI"]))/(2*stats::qnorm(.975)),
    r2z(data[is.na(data[,"se"]) & !.row_NA(data[,c("r", "lCI", "uCI")]), "r"])
  )
  data[is.na(data[,"se"]) & !.row_NA(data[,c("logOR", "lCI", "uCI")]), "se"] <-
    (data[is.na(data[,"se"]) & !.row_NA(data[,c("logOR", "lCI", "uCI")]),"uCI"] -
       data[is.na(data[,"se"]) & !.row_NA(data[,c("logOR", "lCI", "uCI")]),"lCI"])/(2*stats::qnorm(.975))
  data[is.na(data[,"se"]) & !.row_NA(data[,c("z", "lCI", "uCI")]), "se"] <-
    (data[is.na(data[,"se"]) & !.row_NA(data[,c("z", "lCI", "uCI")]),"uCI"] -
       data[is.na(data[,"se"]) & !.row_NA(data[,c("z", "lCI", "uCI")]),"lCI"])/(2*stats::qnorm(.975))


  ### calculate sample sizes
  # based on effect sizes and standard errors
  data[is.na(data[,"n"]) & !.row_NA(data[,c("d", "se")]),"n"] <- n_d(
    data[is.na(data[,"n"]) & !.row_NA(data[,c("d", "se")]),"d"], data[is.na(data[,"n"]) & !.row_NA(data[,c("d", "se")]),"se"])
  data[is.na(data[,"n"]) & !.row_NA(data[,c("z", "se")]),"n"] <- n_z(
    data[is.na(data[,"n"]) & !.row_NA(data[,c("z", "se")]),"se"])
  data[is.na(data[,"n"]) & !.row_NA(data[,c("r", "se")]),"n"] <- n_r(
    data[is.na(data[,"n"]) & !.row_NA(data[,c("r", "se")]),"r"], data[is.na(data[,"n"]) & !.row_NA(data[,c("r", "se")]),"se"])

  # based on effect sizes and t-statistics
  data[is.na(data[,"n"]) & !.row_NA(data[,c("d", "t")]),"n"] <- .n_dt(
    data[is.na(data[,"n"]) & !.row_NA(data[,c("d", "t")]),"d"], data[is.na(data[,"n"]) & !.row_NA(data[,c("d", "t")]),"t"])
  data[is.na(data[,"n"]) & !.row_NA(data[,c("z", "t")]),"n"] <- .n_zt(
    data[is.na(data[,"n"]) & !.row_NA(data[,c("z", "t")]),"z"], data[is.na(data[,"n"]) & !.row_NA(data[,c("z", "t")]),"t"])
  data[is.na(data[,"n"]) & !.row_NA(data[,c("r", "t")]),"n"] <- .n_rt(
    data[is.na(data[,"n"]) & !.row_NA(data[,c("r", "t")]),"r"], data[is.na(data[,"n"]) & !.row_NA(data[,c("r", "t")]),"t"])

  # check whether all of them are positive
  if(any(data[!is.na(data[,"n"]), "n"] < 0)){
    stop("One of the effect sizes and standard errors implies a negative sample size. Please, check the input.")
  }

  ### calculate standard errors on effect sizes and sample sizes (pushes some remaining info from t-stats)
  data[is.na(data[,"se"]) & !.row_NA(data[,c("d", "n")]),"se"] <- se_d(
    data[is.na(data[,"se"]) & !.row_NA(data[,c("d", "n")]),"d"], data[is.na(data[,"se"]) & !.row_NA(data[,c("d", "n")]),"n"])
  data[is.na(data[,"se"]) & !.row_NA(data[,c("z", "n")]),"se"] <- se_z(
    data[is.na(data[,"se"]) & !.row_NA(data[,c("z", "n")]),"n"])
  data[is.na(data[,"se"]) & !.row_NA(data[,c("r", "n")]),"se"] <- se_r(
    data[is.na(data[,"se"]) & !.row_NA(data[,c("r", "n")]),"r"], data[is.na(data[,"se"]) & !.row_NA(data[,c("r", "n")]),"n"])
  # or add standard errors directly on t-statistics for logOR
  data[is.na(data[,"se"]) & !.row_NA(data[,c("logOR", "t")]),"se"] <-
    data[is.na(data[,"se"]) & !.row_NA(data[,c("logOR", "t")]),"logOR"] / data[is.na(data[,"se"]) & !.row_NA(data[,c("logOR", "t")]),"t"]


  ### compute test statistics based on sample sizes and standard errors / sample sizes
  data[is.na(data[,"t"]) & !.row_NA(data[,c("d", "n")]),"t"] <- .t_dn(
    data[is.na(data[,"t"]) & !.row_NA(data[,c("d", "n")]),"d"], data[is.na(data[,"t"]) & !.row_NA(data[,c("d", "n")]),"n"])
  data[is.na(data[,"t"]) & !.row_NA(data[,c("z", "n")]),"t"] <- .t_zn(
    data[is.na(data[,"t"]) & !.row_NA(data[,c("z", "n")]),"z"], data[is.na(data[,"t"]) & !.row_NA(data[,c("z", "n")]),"n"])
  data[is.na(data[,"t"]) & !.row_NA(data[,c("r", "n")]),"t"] <- .t_rn(
    data[is.na(data[,"t"]) & !.row_NA(data[,c("r", "n")]),"r"], data[is.na(data[,"t"]) & !.row_NA(data[,c("r", "n")]),"n"])
  data[is.na(data[,"t"]) & !.row_NA(data[,c("logOR", "se")]),"t"] <-
    data[is.na(data[,"t"]) & !.row_NA(data[,c("logOR", "se")]),"logOR"] / data[is.na(data[,"t"]) & !.row_NA(data[,c("logOR", "se")]),"se"]

  # transform effect sizes and standard errors to the required metric
  if(transformation == "d"){
    data[is.na(data[,"d"]) & !is.na(data[,"r"]),    "d"] <-     r2d(data[is.na(data[,"d"]) & !is.na(data[,"r"]),        "r"])
    data[is.na(data[,"d"]) & !is.na(data[,"z"]),    "d"] <-     z2d(data[is.na(data[,"d"]) & !is.na(data[,"z"]),        "z"])
    data[is.na(data[,"d"]) & !is.na(data[,"logOR"]),"d"] <- logOR2d(data[is.na(data[,"d"]) & !is.na(data[,"logOR"]),"logOR"])

    data[!is.na(data[,"se"]) & !is.na(data[,"r"]),    "se"] <-     se_r2se_d(
      data[!is.na(data[,"se"]) & !is.na(data[,"r"]),        "se"],
      data[!is.na(data[,"se"]) & !is.na(data[,"r"]),        "r"])
    data[!is.na(data[,"se"]) & !is.na(data[,"z"]),    "se"] <-     se_z2se_d(
      data[!is.na(data[,"se"]) & !is.na(data[,"z"]),        "se"],
      data[!is.na(data[,"se"]) & !is.na(data[,"z"]),        "z"])
    data[!is.na(data[,"se"]) & !is.na(data[,"logOR"]),"se"] <- se_logOR2se_d(
      data[!is.na(data[,"se"]) & !is.na(data[,"logOR"]),    "se"])

    # add CI for export
    data[,"lCI"] <- z2d(d2z(data[,"d"]) + stats::qnorm(0.025) * se_d2se_z(data[,"se"], data[,"d"]))
    data[,"uCI"] <- z2d(d2z(data[,"d"]) + stats::qnorm(0.975) * se_d2se_z(data[,"se"], data[,"d"]))

  }else if(transformation == "z"){
    data[is.na(data[,"z"]) & !is.na(data[,"d"]),    "z"] <-     d2z(data[is.na(data[,"z"]) & !is.na(data[,"d"]),        "d"])
    data[is.na(data[,"z"]) & !is.na(data[,"r"]),    "z"] <-     r2z(data[is.na(data[,"z"]) & !is.na(data[,"r"]),        "r"])
    data[is.na(data[,"z"]) & !is.na(data[,"logOR"]),"z"] <- logOR2z(data[is.na(data[,"z"]) & !is.na(data[,"logOR"]),"logOR"])

    data[!is.na(data[,"se"]) & !is.na(data[,"d"]),    "se"] <-     se_d2se_z(
      data[!is.na(data[,"se"]) & !is.na(data[,"d"]),        "se"],
      data[!is.na(data[,"se"]) & !is.na(data[,"d"]),        "d"])
    data[!is.na(data[,"se"]) & !is.na(data[,"r"]),    "se"] <-     se_r2se_z(
      data[!is.na(data[,"se"]) & !is.na(data[,"r"]),        "se"],
      data[!is.na(data[,"se"]) & !is.na(data[,"r"]),        "r"])
    data[!is.na(data[,"se"]) & !is.na(data[,"logOR"]),"se"] <- se_logOR2se_z(
      data[!is.na(data[,"se"]) & !is.na(data[,"logOR"]),    "se"],
      data[!is.na(data[,"se"]) & !is.na(data[,"logOR"]), "logOR"])

    # add CI for export
    data[,"lCI"] <- data[,"z"] + stats::qnorm(0.025) * data[,"se"]
    data[,"uCI"] <- data[,"z"] + stats::qnorm(0.975) * data[,"se"]

  }else if(transformation == "r"){
    data[is.na(data[,"r"]) & !is.na(data[,"d"]),    "r"] <-     d2r(data[is.na(data[,"r"]) & !is.na(data[,"d"]),        "d"])
    data[is.na(data[,"r"]) & !is.na(data[,"z"]),    "r"] <-     z2r(data[is.na(data[,"r"]) & !is.na(data[,"z"]),        "z"])
    data[is.na(data[,"r"]) & !is.na(data[,"logOR"]),"r"] <- logOR2r(data[is.na(data[,"r"]) & !is.na(data[,"logOR"]),"logOR"])

    data[!is.na(data[,"se"]) & !is.na(data[,"d"]),    "se"] <-     se_d2se_r(
      data[!is.na(data[,"se"]) & !is.na(data[,"d"]),        "se"],
      data[!is.na(data[,"se"]) & !is.na(data[,"d"]),        "d"])
    data[!is.na(data[,"se"]) & !is.na(data[,"z"]),    "se"] <-     se_z2se_r(
      data[!is.na(data[,"se"]) & !is.na(data[,"z"]),        "se"],
      data[!is.na(data[,"se"]) & !is.na(data[,"z"]),        "z"])
    data[!is.na(data[,"se"]) & !is.na(data[,"logOR"]),"se"] <- se_logOR2se_r(
      data[!is.na(data[,"se"]) & !is.na(data[,"logOR"]),   "se"],
      data[!is.na(data[,"se"]) & !is.na(data[,"logOR"]),"logOR"])

    # add CI for export
    data[,"lCI"] <- z2r(r2z(data[,"r"]) + stats::qnorm(0.025) * se_r2se_z(data[,"se"], data[,"r"]))
    data[,"uCI"] <- z2r(r2z(data[,"r"]) + stats::qnorm(0.975) * se_r2se_z(data[,"se"], data[,"r"]))

  }else if(transformation == "logOR"){
    data[is.na(data[,"logOR"]) & !is.na(data[,"d"]), "logOR"] <- d2logOR(data[is.na(data[,"logOR"]) & !is.na(data[,"d"]), "d"])
    data[is.na(data[,"logOR"]) & !is.na(data[,"r"]), "logOR"] <- r2logOR(data[is.na(data[,"logOR"]) & !is.na(data[,"r"]), "r"])
    data[is.na(data[,"logOR"]) & !is.na(data[,"z"]), "logOR"] <- z2logOR(data[is.na(data[,"logOR"]) & !is.na(data[,"z"]), "z"])

    data[!is.na(data[,"se"]) & !is.na(data[,"d"]), "se"] <- se_d2se_logOR(
      data[!is.na(data[,"se"]) & !is.na(data[,"d"]),        "se"])
    data[!is.na(data[,"se"]) & !is.na(data[,"r"]), "se"] <- se_r2se_logOR(
      data[!is.na(data[,"se"]) & !is.na(data[,"r"]),        "se"],
      data[!is.na(data[,"se"]) & !is.na(data[,"r"]),         "r"])
    data[!is.na(data[,"se"]) & !is.na(data[,"z"]), "se"] <- se_z2se_logOR(
      data[!is.na(data[,"se"]) & !is.na(data[,"z"]),        "se"],
      data[!is.na(data[,"se"]) & !is.na(data[,"z"]),        "z"])

    # add CI for export
    data[,"lCI"] <- data[,"logOR"] + stats::qnorm(0.025) * data[,"se"]
    data[,"uCI"] <- data[,"logOR"] + stats::qnorm(0.975) * data[,"se"]
  }

  if(return_all){
    return(data)
  }else{
    output$y           <- data[,transformation]
    output$se          <- data[,"se"]
    output$study_names <- data[,"study_names"]
    output$study_ids   <- data[,"study_ids"]
    output$weight      <- data[,"weight"]
    attr(output, "effect_measure")   <- transformation
    attr(output, "original_measure") <- original_measure
    attr(output, "all_independent")  <- all(is.na(data[,"study_ids"]))
    attr(output, "weighted")         <- !all(is.na(data[,"weight"]))
    class(output) <- c(class(output), "data.RoBMA")

    if(anyNA(data[,"se"]) | anyNA(data[,"se"])){
      stop("One of the effect sizes and/or standard errors could not be transformed. Please, check the input.")
    }

    return(output)
  }
}



.transform_posterior       <- function(object, current_scale, output_scale){

  for(type in c("posteriors", "posteriors_conditional")){

    if(!is.null(object$RoBMA[[type]][["mu"]])){
      object$RoBMA[[type]][["mu"]] <- .transform_mu(
        object$RoBMA[[type]][["mu"]],
        current_scale,
        output_scale
      )
    }

    if(!is.null(object$RoBMA[[type]][["tau"]])){
      object$RoBMA[[type]][["tau"]] <- .scale(
        object$RoBMA[[type]][["tau"]],
        current_scale,
        output_scale
      )
    }

    if(!is.null(object$RoBMA[[type]][["PEESE"]])){
      # the transformation for PEESE is inverse
      object$RoBMA[[type]][["PEESE"]] <- .scale(
        object$RoBMA[[type]][["PEESE"]],
        output_scale,
        current_scale
      )
    }
  }

  for(type in c("posteriors_predictors", "posteriors_predictors_conditional")){

    for(i in seq_along(object$RoBMA[[type]])){

      if(inherits(object$RoBMA[[type]][[i]], "mixed_posteriors.factor")){
        for(j in 1:ncol(object$RoBMA[[type]][[i]])){
          object$RoBMA[[type]][[i]][,j] <- .transform_mu(
            object$RoBMA[[type]][[i]][,j],
            current_scale,
            output_scale
          )
        }
      }else if(inherits(object$RoBMA[[type]][[i]], "mixed_posteriors.simple")){
        object$RoBMA[[type]][[i]] <- .transform_mu(
          object$RoBMA[[type]][[i]],
          current_scale,
          output_scale
        )
      }

    }

  }


  object$add_info[["output_scale"]] <- output_scale

  return(object)
}
.transform_model_posterior <- function(model,  current_scale, output_scale){

  if(length(model[["fit"]]) == 0){
    return(model)
  }

  for(chain in seq_along(model[["fit"]][["mcmc"]])){
    for(par in colnames(model[["fit"]][["mcmc"]][[chain]])){

      if(par == "mu"){
        model[["fit"]][["mcmc"]][[chain]][,par] <- .transform_mu(
          model[["fit"]][["mcmc"]][[chain]][,par],
          current_scale,
          output_scale
        )
      }else if(par == "tau"){
        model[["fit"]][["mcmc"]][[chain]][[par]] <- .scale(
          model[["fit"]][["mcmc"]][[chain]][,par],
          current_scale,
          output_scale
        )
      }else if(par == "PEESE"){
        # the transformation for PEESE is inverse
        model[["fit"]][["mcmc"]][[chain]][[par]] <- .scale(
          model[["fit"]][["mcmc"]][[chain]][,par],
          output_scale,
          current_scale
        )
      }
    }

  }

  model[["output_scale"]] <- output_scale

  return(model)
}
#### effect sizes transformation ####
# based on Borenstein, M., Hedges, L. V., Higgins, J. P., & Rothstein, H. R. (2011). Introduction to meta-analysis. John Wiley & Sons.

#' @title Effect size transformations
#'
#' @description Functions for transforming between different
#' effect size measures.
#'
#' @param d Cohen's d.
#' @param r correlation coefficient.
#' @param z Fisher's z.
#' @param logOR log(odds ratios).
#' @param OR offs ratios.
#'
#' @details All transformations are based on
#' \insertCite{borenstein2011introduction}{RoBMA}. In case that
#' a direct transformation is not available, the transformations
#' are chained to provide the effect size of interest.
#'
#' @export d2r
#' @export d2z
#' @export d2logOR
#' @export d2OR
#' @export r2d
#' @export r2z
#' @export r2logOR
#' @export r2OR
#' @export z2d
#' @export z2r
#' @export z2logOR
#' @export z2OR
#' @export logOR2d
#' @export logOR2r
#' @export logOR2z
#' @export logOR2OR
#' @export OR2d
#' @export OR2r
#' @export OR2z
#' @export OR2logOR
#'
#' @name effect_sizes
#'
#' @references
#' \insertAllCited{}
#' @seealso [standard_errors()], [sample_sizes()]
NULL

#' @title Standard errors transformations
#'
#' @description Functions for transforming between
#' standard errors of different effect size measures.
#'
#' @param se_d standard error of Cohen's d
#' @param se_r standard error of correlation coefficient
#' @param se_z standard error of Fisher's z
#' @param se_logOR standard error of log(odds ratios)
#' @param d Cohen's d
#' @param r correlation coefficient
#' @param z Fisher's z
#' @param logOR log(odds ratios)
#'
#' @details Transformations for Cohen's d, Fisher's z, and log(OR) are
#' based on \insertCite{borenstein2011introduction}{RoBMA}. Calculations
#' for correlation coefficient were modified to make the standard error
#' corresponding to the computed on Fisher's z scale under the same sample
#' size (in order to make all other transformations consistent). In case that
#' a direct transformation is not available, the transformations
#' are chained to provide the effect size of interest.
#'
#' It is important to keep in mind that the transformations are only
#' approximations to the true values. From our experience,
#' \code{se_d2se_z} works well for values of se(Cohen's d) < 0.5. Do
#' not forget that the effect sizes are standardized and variance of
#' Cohen's d = 1. Therefore, a standard error of study cannot be larger
#' unless the participants provided negative information (of course, the
#' variance is dependent on the effect size as well, and, can therefore be
#' larger).
#'
#' When setting prior distributions, do NOT attempt to transform a standard
#' normal distribution on Cohen's d (mean = 0, sd = 1) to a normal
#' distribution on Fisher's z with mean 0 and sd = \code{se_d2se_z(0, 1)}.
#' The approximation does NOT work well in this range of values. Instead,
#' approximate the sd of distribution on Fisher's z using samples in this way:
#' \code{sd(d2z(rnorm(10000, 0, 1)))} or, specify the distribution on Cohen's d
#' directly.
#'
#'
#' @export se_d2se_r
#' @export se_d2se_z
#' @export se_d2se_logOR
#' @export se_r2se_d
#' @export se_r2se_z
#' @export se_r2se_logOR
#' @export se_z2se_d
#' @export se_z2se_r
#' @export se_z2se_logOR
#' @export se_logOR2se_d
#' @export se_logOR2se_r
#' @export se_logOR2se_z
#'
#' @name standard_errors
#'
#' @references
#' \insertAllCited{}
#' @seealso [effect_sizes()], [sample_sizes()]
NULL

#' @title Sample sizes to standard errors calculations
#'
#' @description Functions for transforming between standard
#' errors and sample sizes (assuming equal sample sizes per group).
#'
#' @param d Cohen's d
#' @param r correlation coefficient
#' @param se standard error of the corresponding effect size
#' @param n sample size of the corresponding effect size
#'
#' @details Calculations for Cohen's d, Fisher's z, and log(OR) are
#' based on \insertCite{borenstein2011introduction}{RoBMA}. Calculations
#' for correlation coefficient were modified to make the standard error
#' corresponding to the computed on Fisher's z scale under the same sample
#' size (in order to make all other transformations consistent). In case that
#' a direct transformation is not available, the transformations
#' are chained to provide the effect size of interest.
#'
#' Note that sample size and standard error calculation for log(OR)
#' is not available. The standard error is highly dependent on the
#' odds within the groups and sample sizes for individual events are
#' required. Theoretically, the sample size could be obtained by
#' transforming the effect size and standard error to a different measure
#' and obtaining the sample size using corresponding function, however,
#' it leads to a very poor approximation and it is not recommended.
#'
#' @export se_d
#' @export se_r
#' @export se_z
#' @export n_d
#' @export n_r
#' @export n_z
#'
#' @name sample_sizes
#'
#' @references
#' \insertAllCited{}
#' @seealso [effect_sizes()], [standard_errors()]
NULL

#' @rdname effect_sizes
d2r     <- function(d) .d2r$fun(d)
#' @rdname effect_sizes
d2z     <- function(d) .d2z$fun(d)
#' @rdname effect_sizes
d2logOR <- function(d) .d2logOR$fun(d)
#' @rdname effect_sizes
d2OR    <- function(d) .d2OR$fun(d)

#' @rdname effect_sizes
r2d     <- function(r) .r2d$fun(r)
#' @rdname effect_sizes
r2z     <- function(r) .r2z$fun(r)
#' @rdname effect_sizes
r2logOR <- function(r) .r2logOR$fun(r)
#' @rdname effect_sizes
r2OR    <- function(r) .r2OR$fun(r)

#' @rdname effect_sizes
z2r     <- function(z) .z2r$fun(z)
#' @rdname effect_sizes
z2d     <- function(z) .z2d$fun(z)
#' @rdname effect_sizes
z2logOR <- function(z) .z2logOR$fun(z)
#' @rdname effect_sizes
z2OR    <- function(z) .z2OR$fun(z)

#' @rdname effect_sizes
logOR2r     <- function(logOR) .logOR2r$fun(logOR)
#' @rdname effect_sizes
logOR2z     <- function(logOR) .logOR2z$fun(logOR)
#' @rdname effect_sizes
logOR2d     <- function(logOR) .logOR2d$fun(logOR)
#' @rdname effect_sizes
logOR2OR    <- function(logOR) .logOR2OR$fun(logOR)

#' @rdname effect_sizes
OR2r     <- function(OR) .OR2r$fun(OR)
#' @rdname effect_sizes
OR2z     <- function(OR) .OR2z$fun(OR)
#' @rdname effect_sizes
OR2logOR <- function(OR) .OR2logOR$fun(OR)
#' @rdname effect_sizes
OR2d     <- function(OR) .OR2d$fun(OR)


# sample size / standard errors calculations
#' @rdname sample_sizes
se_d     <- function(d, n) sqrt(4/n + d^2/(2*n))
#' @rdname sample_sizes
n_d      <- function(d, se) (d^2 + 8) / (2 * se^2)
#' @rdname sample_sizes
se_r     <- function(r, n){
  # sqrt((1-r^2)^2/(n-1)) : according to the Borenstein, however, it is not consistent with the remaining transformations
  d    <- r2d(r)
  se_d <- se_d(d, n)
  se_d2se_r(se_d, d)
}
#' @rdname sample_sizes
n_r      <- function(r, se){
  # (r^4 - 2*r^2 + se^2 + 1)/se^2 : according to the Borenstein, however, it is not consistent with the remaining
  d    <- r2d(r)
  se_d <- se_r2se_d(se, r)
  n_d(d, se_d)
}
#' @rdname sample_sizes
se_z     <- function(n) sqrt(1/(n-3))
#' @rdname sample_sizes
n_z      <- function(se) (3*se^2 + 1)/se^2

# sample size / standard errors calculations based on chaining the transformations among the remaining effect sizes (introduces more error due to multiple approximations and missing detailed information - terrible, don't use)
# se_logOR <- function(logOR, n){
#   d    <- logOR2d(logOR)
#   se_d <- se_d(d, n)
#   se_d2se_logOR(se_d)
# }
# n_logOR  <- function(logOR, se){
#   (3*logOR^2 + 8*pi^2)/(6*se^2) # backsolved by WolframAlpha
# }

# transformation between standard errors of effect sizes
#' @rdname standard_errors
se_d2se_logOR <- function(se_d, logOR) sqrt(se_d^2 * pi^2 / 3)
#' @rdname standard_errors
se_d2se_r     <- function(se_d, d) sqrt((4^2*se_d^2)/(d^2+4)^3)
#' @rdname standard_errors
se_r2se_d     <- function(se_r, r) sqrt((4*se_r^2)/(1-r^2)^3)
#' @rdname standard_errors
se_logOR2se_d <- function(se_logOR, logOR) sqrt(se_logOR^2 * 3 / pi^2)

# compound transformations
#' @rdname standard_errors
se_d2se_z     <- function(se_d, d){
  n <- n_d(d, se_d)
  se_z(n)
}
#' @rdname standard_errors
se_r2se_z     <- function(se_r, r){
  n <- n_r(r, se_r)
  se_z(n)
}
#' @rdname standard_errors
se_r2se_logOR <- function(se_r, r){
  se_d <- se_r2se_d(se_r, r)
  se_d2se_logOR(se_d)
}
#' @rdname standard_errors
se_logOR2se_r <- function(se_logOR, logOR){
  se_d <- se_logOR2se_d(se_logOR)
  d    <- logOR2d(logOR)
  se_d2se_r(se_d, d)
}
#' @rdname standard_errors
se_logOR2se_z <- function(se_logOR, logOR){
  se_d <- se_logOR2se_d(se_logOR)
  d    <- logOR2d(logOR)
  se_d2se_z(se_d, d)
}
#' @rdname standard_errors
se_z2se_d     <- function(se_z, z){
  n <- n_z(se_z)
  d <- z2d(z)
  se_d(d, n)
}
#' @rdname standard_errors
se_z2se_r     <- function(se_z, z){
  r <- z2r(z)
  n <- n_z(se_z)
  se_r(r, n)
}
#' @rdname standard_errors
se_z2se_logOR <- function(se_z, z){
  se_d <- se_z2se_d(se_z, z)
  se_d2se_logOR(se_d)
}


# sample size based on effect sizes and test statistics
.n_rt <- function(r, t){
  (r^2 * -(t^2) + 2 * r^2 + t^2 ) / r^2
}
.n_dt <- function(d, t){
  4*t^2/d^2
}
.n_zt <- function(z, t){
  n_z(z/t)
}
# test statistics based on effect sizes and sample sizes
.t_rn <- function(r, n){
  r * sqrt((n - 2)/(1 - r^2))
}
.t_dn <- function(d, n){
  d * sqrt(n)/2
}
.t_zn <- function(z, n){
  z/se_z(n)
}
# effect sizes based on test statistics and sample sizes
.r_tn <- function(t, n){
  t / sqrt(n + t^2 - 2)
}
.d_tn <- function(t, n){
  2 * t / sqrt(n)
}
.z_tn <- function(t, n){
  t * se_z(n)
}
# cutoffs based on effect sizes, standard errors, and t-statistics
# 'c_z' and 'c_logOR' do not require the effect sizes - keeping them for consistency when creating calls to the functions automatically
.c_r     <- function(r, se, t){
  .r_tn(t, n_r(r, se))
}
.c_d     <- function(d, se, t){
  .d_tn(t, n_d(d, se))
}
.c_z     <- function(z, se, t){
  t * se
}
.c_logOR <- function(logOR, se, t){
  t * se
}

.get_cutoffs <- function(y, se, prior, original_measure, effect_measure){

  crit_y <- matrix(ncol = length(prior$parameters$steps), nrow = length(y))

  if(effect_measure == "y"){

    crit_t <- matrix(ncol = 0, nrow = length(y))

    for(step in prior$parameters$steps){
      if(prior$distribution == "one.sided"){
        crit_t <- cbind(crit_t, stats::qnorm(step,   lower.tail = FALSE))
      }else if(prior$distribution == "two.sided"){
        crit_t <- cbind(crit_t, stats::qnorm(step/2, lower.tail = FALSE))
      }
    }

    crit_y <- crit_t * matrix(se, ncol = ncol(crit_t), nrow = nrow(crit_t))

  }else{

    for(i in 1:length(y)){

      backtransformed_y    <- do.call(
        .get_transformation(effect_measure, original_measure[i]),
        args = list(y[i]))

      backtransformed_y_se <- do.call(
        .get_transformation_se(effect_measure, original_measure[i]),
        args = list(se[i], y[i]))

      backtransformed_df   <- switch(
        original_measure[i],
        "d"     = n_d(backtransformed_y, backtransformed_y_se) - 2,
        "r"     = n_r(backtransformed_y, backtransformed_y_se) - 2,
        "z"     = n_z(backtransformed_y_se) - 2,
        "logOR" = Inf
      )

      if(grepl("one.sided", prior[["distribution"]])){
        backtransformed_t    <- stats::qt(prior$parameters[["steps"]],   backtransformed_df, lower.tail = FALSE)
      }else if(grepl("two.sided", prior[["distribution"]])){
        backtransformed_t    <- stats::qt(prior$parameters[["steps"]]/2, backtransformed_df, lower.tail = FALSE)
      }

      backtransformed_c    <- do.call(
        eval(parse(text = paste0(".c_", original_measure[i]))),
        args = list(backtransformed_y, backtransformed_y_se, backtransformed_t))

      # transform back the backtransformed measures
      crit_y[i,] <- do.call(
        .get_transformation(original_measure[i], effect_measure),
        args = list(backtransformed_c))
    }
  }

  return(crit_y)
}


.transformation_var    <- function(name){
  if(!name %in% c("fishers_z", "cohens_d", "r", "logOR", "none"))
    stop("Unknown effect size / transformation. The available options are 'fishers_z', 'cohens_d', 'r', 'logOR', and 'none'.")
  switch(
    name,
    "fishers_z" = "z",
    "cohens_d"  = "d",
    "r"         = "r",
    "logOR"     = "logOR",
    "none"      = "y"
  )
}
.transformation_invar  <- function(name){
  if(!name %in% c("z", "d", "r", "logOR", "y"))
    stop("Unknown effect size / transformation shortcut. The available options are 'z', 'd', 'r', 'logOR', and 'y'.")
  switch(
    name,
    "z"     = "fishers_z",
    "d"     = "cohens_d",
    "r"     = "r",
    "logOR" = "logOR",
    "y"     = "none"
  )
}
.transformation_names  <- function(var){
  if(!var %in% c("z", "d", "r", "logOR", "y"))
    stop("Unknown variable type. The available options are 'z', 'd', 'r', 'logOR', and 'y'.")
  switch(
    var,
    "z"       = "Fisher's z",
    "d"       = "Cohen's d",
    "r"       = "correlation",
    "logOR"   = "log(OR)",
    "y"       = "none"
  )
}

.get_transformation    <- function(from, to){
  if(any(c(from, to) == "y"))
    stop("Prior / effect size transformations are not available for unstandardized effect sizes.")
  if(from == to){
    return(function(x)x)
  }else{
    return(eval(parse(text = paste0(from, "2", to))))
  }
}
.get_transformation_se <- function(from, to){
  if(any(c(from, to) == "y"))
    stop("Prior / effect size transformations are not available for unstandardized effect sizes.")
  if(from == to){
    return(function(se_x, x)se_x)
  }else{
    return(eval(parse(text = paste0("se_", from, "2", "se_", to))))
  }
}

# approximate linear transformations used for prior distribution - as in metaBMA package
scale_d2r     <- function(d) .scale_d2r$fun(d)
scale_d2z     <- function(d) .scale_d2z$fun(d)
scale_d2logOR <- function(d) .scale_d2logOR$fun(d)

scale_r2d     <- function(r) .scale_r2d$fun(r)
scale_r2z     <- function(r) .scale_r2z$fun(r)
scale_r2logOR <- function(r) .scale_r2logOR$fun(r)

scale_z2r     <- function(z) .scale_z2r$fun(z)
scale_z2d     <- function(z) .scale_z2d$fun(z)
scale_z2logOR <- function(z) .scale_z2logOR$fun(z)

scale_logOR2r     <- function(logOR) .scale_logOR2r$fun(logOR)
scale_logOR2z     <- function(logOR) .scale_logOR2z$fun(logOR)
scale_logOR2d     <- function(logOR) .scale_logOR2d$fun(logOR)


.get_scale   <- function(from, to){
  if(any(c(from, to) %in% c("y")))
    stop("Prior rescaling is not available for unstandardized effect sizes.")
  if(from == to){
    return(function(x)x)
  }else{
    return(eval(parse(text = paste0("scale_", from, "2", to))))
  }
}
.scale       <- function(x, from, to){
  do.call(
    .get_scale(from, to),
    args = list(x))
}
.get_scale_b <- function(from, to){
  return(switch(
    paste0(from, "2", to),
    "d2logOR"  = pi/sqrt(3),
    "logOR2d"  = 1/(pi/sqrt(3)),
    "d2z"      = 1/2,
    "z2d"      = 2,
    "d2r"      = 1/2,
    "r2d"      = 2,
    "z2logOR"  = 2 * (pi/sqrt(3)),
    "logOR2z"  = (1/(pi/sqrt(3))) * (1/2),
    "z2r"      = (1/2) * 2,
    "r2z"      = 2 * (1/2),
    "logOR2r"  = (1/(pi/sqrt(3))) * (1/2),
    "r2logOR"  = 2 * (pi/sqrt(3))
  ))
}

.scale_note <- function(prior_scale, output_scale){
  return(sprintf(
    "The estimates are summarized on the %1$s scale (priors were specified on the %2$s scale).",
    .transformation_names(output_scale),
    .transformation_names(prior_scale)))
}



.transform_mu    <- function(mu, from, to){
  do.call(
    .get_transformation(from, to),
    args = list(mu))
}
.transform_tau   <- function(tau, mu, from, to){
  if(all(c(from, to) %in% c("d", "logOR"))){
    return(do.call(
      .get_transformation_se(from, to),
      args = list(tau)))
  }else{
    return(do.call(
      .get_transformation_se(from, to),
      args = list(tau, if(!is.null(mu)) mu else 0)))
  }
}
# .transform_PEESE <- function(PEESE, mu, from, to){
#   do.call(
#     .get_transformation_PEESE(from, to),
#     args = list(PEESE, if(!is.null(mu)) mu else 0))
# }




#### helper functions ####
# helper functions
.row_NA <- function(data){
  if(!is.data.frame(data))
    stop("data must be a data.frame")
  apply(data, 1, anyNA)
}
