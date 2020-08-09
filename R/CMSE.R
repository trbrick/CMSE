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
#' 
#' @export
computeExperimental <- function(dataset, intervention, posttest, outcomes, covariates=NULL) {
  # posttestFormula <- as.formula(paste(posttest, "~", paste(c(intervention, covariates), collapse=' + ')))

  # Semi-standardized:
  posttestFormula <- as.formula(paste("scale(", posttest, ") ~ ", paste(c(intervention, covariates), collapse=' + ')))
  regressPosttest <- lm(posttestFormula, data=dataset, na.action =na.exclude)
  a.experimental <- coef(regressPosttest)[intervention]

  c.experimental <- as.list(rep(NA, length(outcomes)))
  names(c.experimental) <- outcomes

  for(anOutcome in outcomes) {
    outcomeFormula <- as.formula(paste("scale(", anOutcome, ") ~", paste(c(intervention, covariates), collapse=' + ')))
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
#' @param intervention character vector containing the name of the intervention.  If present, the intervention will be isolated to avoid causing misfit.
#' @param ... additional arguments to \link{mxTryHard()}
#' @param returnModels (default FALSE) whether to return the models themselves in addition to the CMSE elements.
#' 
#' @return a list containing the model-implied effect of the (standardized) posttest on each (standardized) outcome
#' 
#' @details This function returns the (standardized) model-implied causal effect of a change to the 
#' posttest change of one (standardized) unit.  CMSE computation requires this to be performed only
#' on the control group.  If the intervention exists in the model, it is isolated by removing all
#' outgoing paths.  This prevents it from influencing the paths.  Its variance is bounded above zero
#' to ensure that the model-implied covariance matrix does not become non-positive-definite.
#' 
#'#' @import MICr
#' 
#' @importFrom methods is
#' @importFrom OpenMx mxData mxModel mxTryHard mxPath
#' 
#' @export
#' 
computeNonExperimental <- function(ctl, models, posttest, outcomes, intervention=NA, 
                                   ..., returnModels=FALSE) {
  if(is.data.frame(ctl)) {
    if(intervention %in% names(ctl)) {
      columns <- intersect(c(posttest, outcomes), names(ctl))
      ctl[,columns] <- lapply(ctl[,columns], 
                              function(x) { return(x/sd(x[ctl$intervention==0], na.rm = TRUE))}
                            )
    }
    mxctl <- mxData(ctl, type="raw")
  } else if(methods::is(ctl, "MxDataStatic") && ctl$type == "raw") {
    mxctl <- ctl
    ctl <- mxctl$observed
    nrows <- mxctl$numObs
    if(intervention %in% names(ctl)) {
      columns <- intersect(c(posttest, outcomes), names(ctl))
      ctl[,columns] <- lapply(ctl[,columns], 
                              function(x) { return(x/sd(x[ctl$intervention==0], na.rm = TRUE))}
      )
    }
    mxctl <- mxData(ctl, type="raw", numObs=nrows)
  } else if(methods::is(ctl, "MxDataStatic")) {
    mxctl <- ctl
    ctl <- mxctl$observed
    nrows <- mxctl$numObs
    warning(paste0("Cannot confirm appropriate standardization.\n",
    "CMSE requires posttests and outcomes be scaled so that",
    "the standard deviation of the control group is 1.0."))
  } else if(is.matrix(ctl) && !is.na(nrows)) {
    warning(paste0("Convering covariance matrix into mxData object. ",
            "Cannot confirm appropriate standardization.\n",
            "CMSE requires posttests and outcomes be scaled so that ",
            "the standard deviation of the control group is 1.0."))
    mxctl <- mxData(ctl, type="cov", numObs = nrows)
  } else {
    stop("ctl and ixn must be data frames, mxData objects, or covariance matrices with nrows")
  }
  if(!is.list(models)){
    models <- list(models)
  }
  # browser()
  ctlModels <- lapply(models, function(x) {
          ctlOnly <- mxModel(x, mxctl)
          if(!is.na(intervention) &&
             intervention %in% ctlOnly$manifestVars &&
             methods::is(ctlOnly, "MxRAMModel")) {
             if(any(x$A$values[intervention,] !=0) ||
                any(x$A$free[intervention,] !=FALSE) ||
                any(x$S$values[intervention,colnames(x$S$values) != intervention] !=0) ||
                any(x$S$free[intervention,colnames(x$S$values) != intervention] !=FALSE) ) {
                    warning(paste("CMSE is not designed for models",
                       "with an endogenous intervention;",
                       "It may not function as expected."))}
          # If there's a manifest intervention, isolate it.
              allVars <- c(x$manifestVars, x$latentVars)
              ixnVar <- x$A$values[intervention, intervention]
              # No loadings from intervention for control group.
              ctlOnly <- mxModel(ctlOnly,
                          mxPath(from=intervention, to=allVars),
                          mxPath(from=intervention, 
                                 to=setdiff(allVars, intervention), 
                                 arrows=2),
                            remove=TRUE)
              ctlOnly <- mxModel(ctlOnly, 
                                 mxPath(from=intervention, arrows=2,
                                        values=max(ixnVar, .01), free=FALSE))
          }
          suppressMessages(
          OpenMx::mxTryHard( ctlOnly,
            extraTries = 50, 
            silent = TRUE)
          )
        }
      )
  
  b.table <- lapply(ctlModels, MICr::MICTable, from = posttest, to=outcomes, 
                      splitByType = FALSE, print=FALSE, standardize=TRUE)
  b.nonexperimental <- rep(list(list(b=NA)), length(ctlModels))
  for(modelNo in seq_along(b.table)) {
    aModel <- data.frame(b.table[[modelNo]]) 
    names(b.nonexperimental)[modelNo] <- names(aModel[3])
    b.nonexperimental[[modelNo]]$b <- aModel[,3]
    names(b.nonexperimental[[modelNo]]$b) <- aModel[,2]
  }
  if(returnModels) return(list(b.nonexperimental, ctlModels))
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
#' 
CMSE <- function(dataset, intervention, models, posttest, outcomes, covariates=NULL, ..., nrows=NA, latentPosttest=NA) {
  # browser()
  if(methods::is(dataset, "MxDataStatic")) {
    dataset <- dataset$observed
  } 
  if(!is.data.frame(dataset)){
    stop("dataset must be a data frame, or raw mxData object")
  }
  
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
  nonExperimental <- computeNonExperimental(dataset[dataset$intervention==FALSE,], 
                                            models, latentPosttest, outcomes,
                                            intervention)
  
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