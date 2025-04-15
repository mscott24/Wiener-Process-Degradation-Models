#interface function that runs user specified models
run_profile_u <- function(df, modtype, params0, ...) { 
  
  if (!all(c("V", "tau", "Patient") %in% names(df))) {
    stop('Missing columns (V, tau, or Patient) in df')
  }
  
  if (modtype=='m2') {
    stats <- u_m2(df=df)
  } else if (modtype=='m3') {
    if (missing(params0)) {params0 <- 1}
    stats <- u_m3(df = df, psi0=params0)
  } else {
    stop('model type not specified.')
  }
  
  return(stats)
}

#M2
u_m2 <- function(df, ...) {
  patient_df <- split(df, df$Patient, drop = TRUE)
  n <- length(patient_df)
  lamda_sigma_f <- lapply(patient_df, function(data) {
    lamda_i <- sum(data$V) / sum(data$tau)
    sigma_sq_i <- sum((data$V - lamda_i * data$tau)^2 / (((nrow(data)-1)) * data$tau))
    return(list(lamda_i = lamda_i, sigma_sq_i = sigma_sq_i))
  })
  
  sigma_sq_mu_2_f <- lapply(patient_df, function(data) {
    lamda_i <- sum(data$V) / sum(data$tau)
    sigma_sq_i <- sum((data$V - lamda_i * data$tau)^2 / (((nrow(data)-1)*sum(data$tau)) * data$tau))
  })
  
  sigma_sq <- sum(sapply(lamda_sigma_f, function(x) x$sigma_sq_i))/n
  mu <- sum(sapply(lamda_sigma_f, function(x) x$lamda_i))/n
  lamda <- sapply(lamda_sigma_f, function(x) x$lamda_i)
  sigma_sq_mu_1 <- sum((mu-lamda)^2)/(n-1) 
  sigma_sq_mu_2 <- sum(sapply(sigma_sq_mu_2_f, function(x) x))/n
  sigma_sq_mu <- sigma_sq_mu_1 - sigma_sq_mu_2
  
  stats_m2 <- list()
  stats_m2$mu <- mu
  stats_m2$sigma_sq <- sigma_sq
  stats_m2$sigma_sq_mu <- sigma_sq_mu
  return(stats_m2)
}

#M3
mu_u_m3 <- function(df, psi, ...){
  patient_df <- split(df, df$Patient, drop = TRUE)
  n <- length(patient_df)
  mu_f <- lapply(patient_df, function(data) {
    Omega <- diag(data$tau)
    FF <- toeplitz(c(2, -1, rep(0, nrow(data)-2)))
    #FF[1,1] <- 1
    D <- Omega + psi*FF
    lamda_i <- (t(data$V) %*% solve(D) %*% data$tau) / 
      (t(data$tau) %*% solve(D) %*% data$tau)
  })
  return(sum(sapply(mu_f, function(x) x))/n)
}
sigma_sq_u_m3 <- function(df, psi, ...){
  patient_df <- split(df, df$Patient, drop = TRUE)
  n <- length(patient_df)
  sigma_sq_f <- lapply(patient_df, function(data) {
    Omega <- diag(data$tau)
    FF <- toeplitz(c(2, -1, rep(0, nrow(data)-2)))
    #FF[1,1] <- 1
    D <- Omega + psi*FF
    lamda_i <- as.numeric((t(data$V) %*% solve(D) %*% data$tau) / 
                            (t(data$tau) %*% solve(D) %*% data$tau))
    sigma_sq_i <- (t(data$V - lamda_i) %*% solve(D) %*% (data$V - lamda_i)) / (nrow(data)-1)
  })
  return(sum(sapply(sigma_sq_f, function(x) x))/n)
}
sigma_sq_mu_u_m3 <- function(df, psi, ...){
  patient_df <- split(df, df$Patient, drop = TRUE)
  n <- length(patient_df)
  sigma_sq_mu2_f <- lapply(patient_df, function(data) {
    Omega <- diag(data$tau)
    FF <- toeplitz(c(2, -1, rep(0, nrow(data)-2)))
    #FF[1,1] <- 1
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
  sigma_sq_mu_1 <- sum((mu-lamda)^2)/(n-1) 
  sigma_sq_mu <- sigma_sq_mu_1 - sigma_sq_mu_2
  return(sigma_sq_mu)
}
psi_ll_u_m3 <- function(df, psi, ...) {
  mu <- mu_u_m3(df=df, psi=psi)
  sigma_sq <- sigma_sq_u_m3(df=df, psi=psi)
  sigma_sq_mu <-  sigma_sq_mu_u_m3(df=df, psi=psi)
  sigma_sq_ep <- psi * sigma_sq
  
  patient_df <- split(df, df$Patient, drop = TRUE)
  n <- length(patient_df)
  ll_f <- lapply(patient_df, function(data) {
    Omega <- diag(data$tau)
    FF <- toeplitz(c(2, -1, rep(0, nrow(data)-2)))
    #FF[1,1] <- 1
    D <- Omega + psi*FF
    Sigma <- sigma_sq_mu * data$tau %*% t(data$tau) + sigma_sq * D
    ll <- (-1)* log(2*pi)/2 * length(data$tau) - 
      (1/2) * log(det(Sigma)) - 
      (1/2) * t(data$V - mu * data$tau) %*% solve(Sigma) %*% (data$V - mu * data$tau)
    return(ll[1])
  })
  ll <- sum(sapply(ll_f, function(x) x))
}
u_m3 <- function(df, psi0, ...) {
  
  if (missing(psi0)) {
    psi0=1 
  } 
  
  stats_m3 <- list()
  df<-tidyr::drop_na(df, V, t)
  tryCatch({
    iter <- optim(par = psi0, 
                  fn = psi_ll_u_m3, df = df, 
                  method = "L-BFGS-B", 
                  lower = 1e-3, upper = Inf,
                  control = list(fnscale = -1))
    
    stats_m3$psi <- iter$par
    stats_m3$mu <- mu_u_m3(df = df, psi = stats_m3$psi)
    stats_m3$sigma_sq <- sigma_sq_u_m3(df = df, psi = stats_m3$psi)
    stats_m3$sigma_sq_mu <- sigma_sq_mu_u_m3(df = df, psi = stats_m3$psi)
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