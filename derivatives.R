#VARIANCE AND DERIVATIVE FUNCTIONS#
calcHessian <- function(df, params, ...) {
  
  n_parms <- length(params)
  if (length(params) < 4) {params[(length(params)+1):4] <- 0}
  patient_df <- split(df, df$Patient, drop = TRUE)
  
  compute_hess_i <- function(data, params) {
    S <- diag(data$tau, nrow = length(data$tau), ncol = length(data$tau))
    n <- length(data$tau)
    A <- toeplitz(c(2, -1, rep(0, n - 2)))
    Omega <- params[2] * S +
      params[3] * data$tau %*% t(data$tau) +
      params[4] * A
    
    if (det(Omega) < 1e-4) {  
      epsilon <- 1e-5
      Omega <- Omega + diag(epsilon, nrow(Omega))
    }
    
    Omega_inv <- solve(Omega)
    ders <- list(S, data$tau %*% t(data$tau), A)
    
    hess_i <- matrix(nrow=n_parms, ncol=n_parms)
    hess_i[1,1] <- t(data$tau) %*% Omega_inv %*% data$tau 
    hess_i[1,2:n_parms] <- hess_i[2:n_parms,1] <- 0
    for (i in 2:n_parms) { 
      for (j in 2:n_parms) {
        hess_i[i,j] <- hess_i[j,i] <- 0.5 * sum(diag(Omega_inv %*% ders[[i-1]] %*% Omega_inv %*% ders[[j-1]]))
      }
    }
    return(-hess_i)
  }
  
  hess_elements <- lapply(patient_df, function(data) compute_hess_i(data, params))
  hess <-  Reduce(`+`, hess_elements)
  
  return(hess)
}
calcGrad <- function(df, params, ...) {
  
  n_parms <- length(param)
  if (length(params) < 4) {params[(length(params)+1):4] <- 0}
  
  patient_df <- split(df, df$Patient, drop = TRUE)
  
  compute_grad_i <- function(data, params) {
    S <- diag(data$tau, nrow = length(data$tau), ncol = length(data$tau))
    n <- length(data$tau)
    A <- toeplitz(c(2, -1, rep(0, n - 2)))
    Omega <- params[2] * S +
      params[3] * data$tau %*% t(data$tau) +
      params[4] * A
    
    if (det(Omega) < 1e-4) {  
      epsilon <- 1e-5
      Omega <- Omega + diag(epsilon, nrow(Omega))
    }
    
    Omega_inv <- solve(Omega)
    ders <- list(S,data$tau %*% t(data$tau),A)
    
    grad_i <- matrix(nrow=n_parms, ncol=1)
    grad_i[1] <- t(data$tau) %*% Omega_inv %*% (data$V - params[1] *data$tau)
    for (p in 2:(n_parms)) {
      grad_i[p] <-  (-1) * 0.5 * (sum(diag(Omega_inv %*% ders[[p-1]])) - 
                                    t(data$V - params[1] *data$tau) %*% 
                                    Omega_inv %*% ders[[p-1]] %*% Omega_inv %*%
                                    (data$V - params[1] *data$tau))
    }
    return(grad_i)
  }
  
  grad_elements <- lapply(patient_df, function(data) compute_grad_i(data, param))
  grad <- Reduce(`+`, grad_elements)
  
  return(grad)
}
calcCIs <- function(df, results) {
  
  par <- as.numeric(results)
  fim <- -calcHessian(df=df, params = par)
  varcov <- solve(fim)
  se <- sqrt(diag(varcov))
  
  z <- qnorm(0.975)
  ci_lb <- numeric(length(par))
  ci_ub <- numeric(length(par))
  
  for (i in seq_along(par)) {
    if (i == 1) {
      ci_lb[i] <- par[i] - z * se[i]
      ci_ub[i] <- par[i] + z * se[i]
    } else {
      log_est <- log(par[i])
      se_log <- se[i] / par[i]
      log_ci_lb <- log_est - z * se_log
      log_ci_ub <- log_est + z * se_log
      ci_lb[i] <- exp(log_ci_lb)
      ci_ub[i] <- exp(log_ci_ub)
    }
  }
  
  out <- data.frame(Estimate = par, SE = se, CI_lower = ci_lb, CI_upper = ci_ub)
  
  return(out)
}

