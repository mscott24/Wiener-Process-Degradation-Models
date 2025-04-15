#' WPDM: Wiener Process Degradation Model Estimation
#'
#' Fits a Wiener process degradation model to longitudinal data using various estimation procedures
#'
#' @param df data.frame containing the input data with subject IDs, time, and outcome variables
#' Default: no default - required field
#' 
#' @param idvar string indicating the subject column name in 'df'
#' Default: 'Patient'
#' 
#' @param timevar string indicating the time column name in 'df'
#' Default: 't'
#' 
#' @param Yvar string indicating the outcome column name in 'df'
#' Default: 'Y'
#' 
#' @param estimator string specifying the estimation method:
#'   - 'optim': maximum likelihood using optimization
#'   - 'mle': profile likelihood based on MLE
#'   - 'u': profile likelihood based on empirically unbiased 
#' Default: 'optim'
#' 
#' @param modtype string indicating the model type:
#'   - 'm1': c(mu, sigma_sq)
#'   - 'm2': c(mu, sigma_sq, sigma_sq_mu) 
#'   - 'm3': c(mu, sigma_sq, sigma_sq_mu, sigma_sq_ep)
#' Default: 'm3'
#' 
#' @param params0 numeric vector of initial parameter values for estimation, 
#' depending on the optimizer:
#'   - 'optim' and 'm1': {mu, sigma_sq}
#'   - 'optim' and 'm2': {mu, sigma_sq, sigma_sq_mu}
#'   - 'optim' and 'm3': {mu, sigma_sq, sigma_sq_mu, sigma_sq_ep}
#'   - 'mle' and 'm1': none required
#'   - 'mle' and 'm2': {sigma_sq, sigma_sq_mu}
#'   - 'mle' and 'm3': {psi}
#'   - 'u' and 'm2': none required
#'   - 'u' and 'm3': {psi}
#' Default: mu=0, sigma_sq=0.1, sigma_sq_mu=0.1, sigma_sq_ep=0.1, psi=1
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{model}{model and estimator}
#'   \item{stats}{estimates, standard errors (SE), and Wald 95% confidence intervals (CIs).}
#'   \item{ll}{log-likelihood}
#' }
#'
#' @details This function preprocesses the data into a difference format and applies one of three
#' estimation strategies to fit a Wiener process degradation model. Internal utility files
#' are sourced dynamically, and confidence intervals are computed from the Hessian.
#'
#' @import MASS rlang tidyr dplyr pbapply parallel
#' @examples
#' \dontrun{
#' result <- WPDM(df = mydata, idvar = "Patient", timevar = "t", Yvar = "Y", 
#'                estimator = "optim", modtype = "m3")
#' print(result$stats)
#' }

WPDM <- function(df, idvar, timevar, Yvar, estimator, modtype, params0, ...) {
  
  #default values
  if (missing(idvar)) {idvar <- 'Patient'}
  if (missing(timevar)) {timevar <- 't'}
  if (missing(Yvar)) {Yvar <- 'Y'}
  if (missing(estimator)) {estimator <- 'optim'} 
  if (missing(modtype)) {modtype <- 'm3'} 
  
  require(MASS)
  require(rlang)
  require(tidyr)
  require(dplyr)
  require(pbapply)
  require(parallel)
  
  files <- c("profile_u.R", "profile_mle.R", "optim_ll.R", "utility.R", "simulations.R", "derivatives.R")
  invisible(lapply(files, function(f) source(f, echo = FALSE)))
  
  df <- FormatDiffDF(df = df, 
                     idvar = idvar, 
                     timevar = timevar, 
                     Yvar = Yvar)
  
  if (estimator=='optim') {
    if (missing(params0)) {
      results <- run_optim_ll(df = df, modtype = modtype)
    }
    else {
      results <- run_optim_ll(df = df, modtype = modtype, params0 = params0)
    }
  }
  
  if (estimator=='mle') {
    if (missing(params0)) {
      results <- run_profile_mle(df = df, modtype = modtype)
    }
    else {
      results <- run_profile_mle(df = df, modtype = modtype, params0 = params0)
    }
    results$psi <- NULL
  }
  
  if (estimator=='u') {
    if (missing(params0)) {
      results <- run_profile_u(df = df, modtype = modtype)
    }
    else {
      results <- run_profile_u(df = df, modtype = modtype, params0 = params0)
    }
    results$psi <- NULL
  }
  
  out <- list()
  out$model <- paste("Wiener process degradation model using", modtype, "with", estimator)
  out$stats <- calcCIs(df=df, results=results)
  out$ll <- ll(df=df, params=as.numeric(results))
  return(out)
}

