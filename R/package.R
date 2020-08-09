
#' CMSE: Causal Mean Squared Error index
#'
#' CMSE: Causal Mean Squared Error index
#'
#' @docType package
#' @name CMSE
#' 
#' @description CMSE computes the Causal Mean Squared Error (CMSE).
#' CMSE indexes the difference between the expected effect of the intervention and
#' the real impact of the intervention.  Specifically, the CMSE computes the true
#' effect of the intervention on each outcome using a simple linear regression
#' (controlling appropriately for covariates) and compares that true effect
#' to the model-implied effect on that outcome.  By focusing on the specific
#' causal predictions of the model, CMSE captures a different slice of information
#' than other fit statistics.
#' 
#' See Wan, et al., for more.
NULL
