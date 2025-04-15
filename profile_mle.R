#interface function that runs user specified models
run_profile_mle <- function(df, modtype, params0, ...) { 

  if (!all(c("V", "tau", "Patient") %in% names(df))) {
    stop('Missing columns (V, tau, or Patient) in df')
  }
  
  if (modtype=='m1') {
    stats <- mle_m1(df=df)
  } else if (modtype=='m2') {
    if (missing(params0)) {params0 <- c(0.1, 0.1)}
    stats <- mle_m2(df=df, params0=params0)
  } else if (modtype=='m3') {
    if (missing(params0)) {params0 <- 1}
    stats <- mle_m3(df = df, psi0=params0)
  } else {
    stop('model type not specified.')
  }
  
  return(stats)
}

#M1
mle_m1 <- function(df, ...) {
  
  patient_df <- split(df, df$Patient, drop = TRUE)
  
  mu_func <- lapply(patient_df, function(data) {
    list(num_mu = sum(data$V), denom_mu = sum(data$tau))
  })
  
  num_mu <- sapply(mu_func, function(x) x$num_mu)
  denom_mu <- sapply(mu_func, function(x) x$denom_mu)
  mu_hat <- sum(num_mu) / sum(denom_mu)
  
  sigma_sq_func <- lapply(patient_df, function(data) {
    S_i <- diag(data$tau)
    num_sigma_sq <- t(data$V - mu_hat * data$tau) %*% solve(S_i) %*%  (data$V - mu_hat * data$tau)
    denom_sigma_sq <- nrow(data)
    list(num_sigma_sq = num_sigma_sq, denom_sigma_sq = denom_sigma_sq)
  })
  
  num_sigma_sq <- sapply(sigma_sq_func, function(x) x$num_sigma_sq)
  demon_sigma_sq <- sapply(sigma_sq_func, function(x) x$denom_sigma_sq)
  sigma_sq_hat <- sum(num_sigma_sq) / sum(demon_sigma_sq)
  
  return(list(mu = mu_hat,
              sigma_sq = sigma_sq_hat))
}

#M2
mu_mle_m2 <- function(df, params,  ...) {
  patient_df <- split(df, df$Patient, drop = TRUE)
  func <- lapply(patient_df, function(data) {
    Tmax <- sum(data$tau)
    w <- Tmax / (params[1] + params[2]*Tmax)
    mu <- sum(data$V) / Tmax
    list(mu = mu, w = w)
  })
  
  mus <- sapply(func, function(x) x$mu)
  ws <- sapply(func, function(x) x$w)
  mu_hat <- sum(mus  * ws/sum(ws))
  return(mu_hat)
}
profileL_mle_m2 <- function(df, params, ...) {
  #params[1] <- wiener process
  #params[2] <- random slope
  #this function will run a profile likelihood over the two variance comps
  
  mu <- mu_mle_m2(df, params)
  
  patient_df <- split(df, df$Patient, drop = TRUE)
  func <- lapply(patient_df, function(data) {
    data <- tidyr::drop_na(data, V, tau)
    S <- diag(data$tau)
    S_inv <- S
    diag(S_inv) <- 1/diag(S)
    Sigma <- params[1] * S + params[2] * data$tau %*% t(data$tau)
    n <- nrow(data)
    Tmax <- sum(data$tau)
    Sigma_det <- params[1]^(n-1) * (params[1]  + params[2] * Tmax) * prod(diag(S)) 
    Sigma_inv <- 1/params[1]  * (S_inv - (params[2])/(params[1]  + params[2] * Tmax) 
                                 * S_inv %*% data$tau %*% t(data$tau) %*% S_inv)
    
    ll <- (-n / 2)*log(2*pi) - 
      (1/2)*log(Sigma_det) - 
      (1/2) * (t(data$V - mu * data$tau) %*% Sigma_inv %*% (data$V - mu * data$tau))
  })
  
  ll <- sum(sapply(func, function(x) x))
}
mle_m2 <- function(df, params0, ...) {
  
  if (missing(params0)) {
    params0 <- c(0.1, 0.1)
  }
  
  df <- tidyr::drop_na(df, V, tau)
  result <- optim(fn = profileL_mle_m2,
                  par = params0,
                  df = df, 
                  method = "L-BFGS-B",
                  lower = c(1e-4, 1e-4), upper = c(Inf, Inf),
                  control = list(fnscale = -1))
  sigma_sq_hat <- result$par[1]
  sigma_sq_mu_hat <- result$par[2]
  
  stats_m2 <- list()
  stats_m2$mu <- mu_mle_m2(df, c(result$par[1], result$par[2]))
  stats_m2$sigma_sq <- result$par[1]
  stats_m2$sigma_sq_mu <- result$par[2]
  return(stats_m2)
}

#M3
mu_mle_m3 <- function(df, psi, ...){
  patient_df <- split(df, df$Patient, drop = TRUE)
  n <- length(patient_df)
  mu_f <- lapply(patient_df, function(data) {
    Omega <- diag(data$tau)
    FF <- toeplitz(c(2, -1, rep(0, nrow(data)-2)))
    D <- Omega + psi*FF
    lamda_i <- (t(data$V) %*% solve(D) %*% data$tau) / 
      (t(data$tau) %*% solve(D) %*% data$tau)
  })
  return(sum(sapply(mu_f, function(x) x))/n)
}
sigma_sq_mle_m3 <- function(df, psi, ...){
  patient_df <- split(df, df$Patient, drop = TRUE)
  n <- length(patient_df)
  sigma_sq_f <- lapply(patient_df, function(data) {
    Omega <- diag(data$tau)
    FF <- toeplitz(c(2, -1, rep(0, nrow(data)-2)))
    D <- Omega + psi*FF
    lamda_i <- as.numeric((t(data$V) %*% solve(D) %*% data$tau) / 
                            (t(data$tau) %*% solve(D) %*% data$tau))
    sigma_sq_i <- (t(data$V - lamda_i) %*% solve(D) %*% (data$V - lamda_i)) / (nrow(data)-1)
  })
  return(sum(sapply(sigma_sq_f, function(x) x))/n)
}
sigma_sq_mu_mle_m3 <- function(df, psi, ...){
  patient_df <- split(df, df$Patient, drop = TRUE)
  n <- length(patient_df)
  sigma_sq_mu2_f <- lapply(patient_df, function(data) {
    Omega <- diag(data$tau)
    FF <- toeplitz(c(2, -1, rep(0, nrow(data)-2)))
    D <- Omega + psi*FF
    lamda_i <- as.numeric((t(data$V) %*% solve(D) %*% data$tau) / 
                            (t(data$tau) %*% solve(D) %*% data$tau))
    sigma_sq_mu_i <- as.numeric((t(data$V - lamda_i*data$tau) %*% solve(D) %*% (data$V - lamda_i*data$tau)) / 
                                  ((nrow(data)-1) * (t(data$tau) %*% solve(D) %*% (data$tau))))
    return(list(sigma_sq_mu_i=sigma_sq_mu_i, lamda_i=lamda_i))
  })
  sigma_sq_mu_2 <- sum(sapply(sigma_sq_mu2_f, function(x) x$sigma_sq_mu_i))/n
  lamda <- sapply(sigma_sq_mu2_f, function(x) x$lamda_i)
  mu <- sum(sapply(sigma_sq_mu2_f, function(x) x$lamda_i))/n
  sigma_sq_mu_1 <- sum((mu-lamda)^2)/(n) #n-1 for U estimator 
  sigma_sq_mu <- sigma_sq_mu_1 - sigma_sq_mu_2
  return(sigma_sq_mu)
}
psi_ll_mle_m3 <- function(df, psi, ...) {
  mu <- mu_mle_m3(df=df, psi=psi)
  sigma_sq <- sigma_sq_mle_m3(df=df, psi=psi)
  sigma_sq_mu <-  sigma_sq_mu_mle_m3(df=df, psi=psi)
  sigma_sq_ep <- psi * sigma_sq
  
  patient_df <- split(df, df$Patient, drop = TRUE)
  n <- length(patient_df)
  ll_f <- lapply(patient_df, function(data) {
    Omega <- diag(data$tau)
    FF <- toeplitz(c(2, -1, rep(0, nrow(data)-2)))
    D <- Omega + psi*FF
    Sigma <- sigma_sq_mu * data$tau %*% t(data$tau) + sigma_sq * D
    ll <- (-1)* log(2*pi)/2 * length(data$tau) - 
      (1/2) * log(det(Sigma)) - 
      (1/2) * t(data$V - mu * data$tau) %*% solve(Sigma) %*% (data$V - mu * data$tau)
    return(ll[1])
  })
  ll <- sum(sapply(ll_f, function(x) x))
}
mle_m3 <- function(df, psi0, ...) {
  
  stats_m3 <- list()
  tryCatch({
    result <- optim(par = psi0, 
                    fn = psi_ll_mle_m3, df = df, 
                    method = "L-BFGS-B", 
                    lower = 1e-3, upper = Inf,
                    control = list(fnscale = -1))
    
    stats_m3$psi <- result$par
    stats_m3$mu <- mu_mle_m3(df = df, psi = stats_m3$psi)
    stats_m3$sigma_sq <- sigma_sq_mle_m3(df = df, psi = stats_m3$psi)
    stats_m3$sigma_sq_mu <- sigma_sq_mu_mle_m3(df = df, psi = stats_m3$psi)
    stats_m3$sigma_sq_ep <- stats_m3$psi * stats_m3$sigma_sq
    
  }, error = function(e) {
    stats_m3$psi <- NULL
    stats_m3$mu <- NULL
    stats_m3$sigma_sq <- NULL
    stats_m3$sigma_sq_mu <- NULL
    stats_m3$sigma_sq_ep <- NULL
    
    message("optim failed: ", e$message)
  })
  return(stats_m3)
}