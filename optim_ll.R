#estimation of degradation models using joint loglikelihood
run_optim_ll <- function(df, modtype, params0, ...) { 
  
  if (!all(c("V", "tau", "Patient") %in% names(df))) {
    stop('Missing columns (V, tau, or Patient) in df')
  }
  
  if (modtype=='m3') {
    
    if (missing(params0)) {
      params0 <- c(0, 0.1, 0.1, 0.1)
    }
    
    result <- optim(fn = ll,
                    par = params0,
                    df = df, 
                    method = "L-BFGS-B",
                    lower = c(-Inf, 1e-4, 1e-4, 1e-4), 
                    upper = c(Inf, Inf, Inf, Inf),
                    control = list(fnscale = -1))
    stats <- list()
    stats$mu <- result$par[1]
    stats$sigma_sq <- result$par[2]
    stats$sigma_sq_mu <- result$par[3]
    stats$sigma_sq_ep <- result$par[4]
  } else if (modtype=='m2') {
    
    if (missing(params0)) {
      params0 <- c(0, 0.1, 0.1)
    }
    
    result <- optim(fn = ll,
                    par = params0,
                    df = df, 
                    method = "L-BFGS-B",
                    lower = c(-Inf, 1e-4, 1e-4), upper = c(Inf, Inf, Inf),
                    control = list(fnscale = -1))
    stats <- list()
    stats$mu <- result$par[1]
    stats$sigma_sq <- result$par[2]
    stats$sigma_sq_mu <- result$par[3]
  } else if (modtype=='m1') {
    
    if (missing(params0)) {
      params0 <- c(0, 0.1)
    }
    
    result <- optim(fn = ll,
                    par = params0,
                    df = df, 
                    method = "L-BFGS-B",
                    lower = c(-Inf, 1e-4), upper = c(Inf, Inf),
                    control = list(fnscale = -1))
    stats <- list()
    stats$mu <- result$par[1]
    stats$sigma_sq <- result$par[2]
  } else {
    stop('model type not specified.')
  }
  
  return(stats)
}
ll <- function(df, params, ...) {
  
  if (length(params) < 4) {params[(length(params)+1):4] <- 0}
  
  patient_df <- split(df, df$Patient, drop = TRUE)
  func <- lapply(patient_df, function(data) {
    S <- diag(data$tau)
    n <- length(data$tau)
    A <- toeplitz(c(2, -1, rep(0, n-2)))
    Sigma <- params[2] * S + 
      params[3] * data$tau %*% t(data$tau) + 
      params[4] * A
    
    ll <- as.numeric((-1) * (n/2) * log(2*pi) - (1/2) * log(det(Sigma)) - 
                       (1/2) * t(data$V - params[1] * data$tau) %*% solve(Sigma) %*% (data$V - params[1] * data$tau)) 
  })
  
  ll <- sum(sapply(func, function(x) x))
  
  return(ll)
}
