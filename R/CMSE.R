#' Compute Experimental CMSE Elements
#'
#' Computes c.experimental and a.experimental values for the CMSE computation.
#'
#' @param dataset the actual data as a data frame or mxDataStatic object
#' @param intervention character string containing the name of the intervention column
#' @param posttest character string containing the name of the posttest outcome
#' @param outcomes character vector including names of all additional outcomes
#' @param covariates character vector including the names of all covariates to be included
#'
#' @return a list containing elements:
#'         - `$a` contains the experimentally-determined effect of the intervention on posttest
#'         - `$c` contains a named list of all outcomes; each outcome contains the experimentally-determined effect of the intervention on that outcome
#' 
#' @importFrom stats as.formula coef lm na.exclude
#' @export
computeExperimental <- function(dataset, intervention, posttest, outcomes, covariates=NULL) {
  posttestFormula <- as.formula(paste(posttest, "~", paste(c(intervention, covariates), collapse=' + ')))
  # Standardized:
  # posttestFormula <- as.formula(paste("scale(", posttest, ") ~ ", paste("scale(", c(intervention, covariates), ")", collapse=' + ')))
  regressPosttest <- lm(posttestFormula, data=dataset, na.action =na.exclude)
  a.experimental <- coef(regressPosttest)[intervention]
  # Standardized:
  #a.experimental <- coef(regressPosttest)[paste0("scale(", intervention, ")")]
  
  c.experimental <- as.list(rep(NA, length(outcomes)))
  names(c.experimental) <- outcomes
  # browser()
  for(anOutcome in outcomes) {
    outcomeFormula <- as.formula(paste(anOutcome, "~", paste(c(intervention, covariates), collapse=' + ')))
    regressOutcome <- lm(outcomeFormula, data=dataset)
    c.experimental[[anOutcome]] <- unname(coef(regressOutcome)[intervention])
  }
  return(list(a=unname(a.experimental), c=unlist(c.experimental)))
}

#' Compute Model-implied CMSE Elements
#'
#' Computes b.nonexperimental and a.nonexperimental values for the CMSE computation.
#'
#' @param ctl the control group only of the experimental data set
#' @param models a list of mxModel objects (like those returned from makeARPanelModels) for fitting
#' @param posttest character string containing the name of the posttest true score (see details), or a vector of the same length as models containing the true scores for each
#' @param outcomes character vector including names of all additional outcomes, or a list the same length as models with outcomes for each model
#'
#' @return a list containing the model-implied effect of the posttest on each outcome
#'
#' @import MICr
#' 
#' @importFrom methods is
#' 
#' @export
computeNonExperimental <- function(ctl, models, posttest, outcomes) {
  if(is.data.frame(ctl)) {
    mxctl <- mxData(ctl, type="raw")
  } else if(methods::is(ctl, "MxDataStatic")) {
    mxctl <- ctl
    ctl <- mxctl$observed
    nrows <- mxctl$numObs
  } else if(is.matrix(ctl) && !is.na(nrows)) {
    warning("Convering covariance matrix into mxData object.")
    mxctl <- mxData(ctl, type="cov", numObs = nrows)
  } else {
    stop("ctl and ixn must be data frames, mxData objects, or covariance matrices with nrows")
  }
  if(!is.list(models)){
    models <- list(models)
  }
  # browser()
  ctlModels <- lapply(models, function(x) {OpenMx::mxTryHard(mxModel(x, mxctl), extraTries = 50)})
  b.table <- lapply(ctlModels, "MICr::MICTable", from = posttest, to=outcomes, splitByType = FALSE)
  b.nonexperimental <- rep(list(list(b=NA)), length(ctlModels))
  for(modelNo in seq_along(b.table)) {
    aModel <- data.frame(b.table[[modelNo]]) 
    names(b.nonexperimental)[modelNo] <- names(aModel[3])
    b.nonexperimental[[modelNo]]$b <- aModel[,3]
    names(b.nonexperimental[[modelNo]]$b) <- aModel[,2]
  }
  return(b.nonexperimental)
}

#' CMSE: Compute the Causal Mean Squared Error
#'
#' `Causal Mean Squared Error (CMSE)` is a fit index that quantifies misfit 
#' between the model-implied causal effects of an intervention and the 
#' empirically determined effects.
#' 
#' This is quantified as the difference between the (regression-determined)
#' empirical effects of the intervention on distal outcomes, and the model's 
#' predictions about those distal effects given the empirical proximal effect.
#' See details for more.
#'
#' @param dataset the actual data as a data frame or mxDataStatic object
#' @param intervention character string containing the name of the intervention column
#' @param models a list of mxModel objects (like those returned from makeARPanelModels) for fitting
#' @param posttest character string containing the name of the posttest outcome
#' @param outcomes character vector including names of all additional outcomes
#' @param covariates character vector including the names of all covariates to be included
#' @param ... Does not accept arguments; only there so that later arguments must be named
#' @param nrows number of data rows in the data set; Required if dataset is a covariance matrix without means
#' @param latentPosttest the latent true posttest score; Required if the posttest outcome has a measurement model (e.g. for the RI-AR model)
#'
#' @return a list containing elements:
#'         - `$CMSE` a named data frame of CMSE scores for each model and each outcome, as well as model-mean CMSEs
#'         - `$experimental` contains the experimentally-determined effect of the intervention on posttest and outcomes
#'         - `$nonexperimental` contains the model-implied effect of the intervention on posttest and each outcome
#'
#' @export
CMSE <- function(dataset, intervention, models, posttest, outcomes, covariates=NULL, ..., nrows=NA, latentPosttest=NA) {
  # browser()
  
  # Default to posttest-is-latent-posttest
  if(is.na(latentPosttest)) {
    latentPosttest=posttest
  }
  
  # Handle single model case
  if(!is.list(models)) {
    models <- list(models)
  }
  
  # Handle no-names case
  if(is.null(names(models))) {
    names(models) <- lapply(models, function(x) {getElement(x, "name")})
  }
  
  experimental <- computeExperimental(dataset, intervention, posttest, outcomes, covariates)
  # Estimated on Control group only; hence dataset$intervention==FALSE
  nonExperimental <- computeNonExperimental(dataset[dataset$intervention==FALSE,], models, latentPosttest, outcomes)
  
  # browser()
  outcomes <- as.list(rep(NA, length(models)))
  names(outcomes) <- names(models)
  CMSE <- outcomes
  for(aModel in models) {
    mName <- aModel$name
    nonExperimental[[mName]]$c <- nonExperimental[[mName]]$b * experimental$a
    cdiff <- nonExperimental[[mName]]$c - experimental$c
    CMSE[[mName]] <- c(mean=mean(cdiff^2), cdiff^2)
  }
  
  if(length(CMSE)==1) {
    CMSE <- CMSE[[1]]
    nonExperimental <- nonExperimental[[1]]
  }
  
  return(list(CMSE=CMSE, experimental=experimental, nonExperimental=nonExperimental))
}