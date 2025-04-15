#data generation tools
DefaultSetup <- function(...) {
  setup <- list()
  setup$n_sims <- 1000  #number of simulations     
  setup$n_subjs <- 100 #number of subjects
  setup$alpha <- 0.05 #coverage probability
  setup$mu <- -1 #slope mean
  setup$sigma_sq <- 0.25 #brownian motion
  setup$sigma_sq_mu <- 0.25 #slope variance
  setup$sigma_sq_ep <- 0.25 #measurement error
  setup$fu_lb <- 3 #lower bound for number of follow-up visits in unequal case
  setup$fu_ub <- 10 #upper bound for number of follow-up visits  in unequal case
  setup$n <- 10 #equal FU case
  setup$home_dir <- '/project/jointage/matt/DM/paper'
  setup$batch_size <- min(100, ceiling(setup$n_sims / (2 * parallel::detectCores())))
  setup$u_lb <- 0.5
  setup$u_ub <- 1.5
  setup$p <- 0
  setup$tau <- 0.05
  return(setup)
}
GenTimes <- function(bounds, lambda, p, tau, ...) {
  if (missing(p)) {p=0}
  if (missing(tau)) {tau=0.05}
  
  if (missing(lambda) & length(bounds)==1) {
    #equal FU case
    return(0:(bounds[1])) 
  }
  
  else {
    #different number of FU visits and duration
    n_samp <- sample((bounds[1]+1):(bounds[2]+1), size = 1)
    t <- cumsum(runif(n_samp, min = 0.5, max = 1.5))
    if (p>0) {
      ind_close <- which(runif(n_samp) < p)
      for (i in ind_close) {
        t[i + 1] <- t[i] + runif(1, 0, tau)
      }
    }
    t <- t - t[1] 
    return(t)
  }
}
GenData <- function(setup, ff, opt, ...) {
  if (!dir.exists(ff)) {dir.create(ff)}
  if (missing(opt)) {opt <- 'all'}
  
  #equal FU visits and duration
  ff_eq <- paste(ff, '/eq.RData', sep='')
  if (!file.exists(ff_eq) & (opt %in% c('all', 'eq'))) {
    print('Generating eq')
    set.seed(2024)
    time_eq <- matrix(rep(t(GenTimes(bounds=setup$n)), setup$n_subjs), 
                      nrow=(setup$n+1), ncol=setup$n_subjs)
    eq <- list()
    eq$setup <- setup
    sims <- MakeSim(setup=setup, time=time_eq)
    eq$sim_m1 <- sims$sim_m1
    eq$sim_m2 <- sims$sim_m2
    eq$sim_m3 <- sims$sim_m3
    eq$sim_m1$dfs <- lapply(eq$sim_m1$dfs, as.data.table)
    eq$sim_m2$dfs <- lapply(eq$sim_m2$dfs, as.data.table)
    eq$sim_m3$dfs <- lapply(eq$sim_m3$dfs, as.data.table)
    save(eq, file=ff_eq)
    rm(eq)
  }
  
  #unequal FU time and number of visits
  ff_ue <- paste(ff, '/ue.RData', sep='')
  if (!file.exists(ff_ue) & (opt %in% c('all', 'ue'))) {
    print('Generating ue')
    set.seed(2024)
    time_ue <- PadList(replicate(setup$n_subjs, 
                                 GenTimes(bounds=c(setup$fu_lb, setup$fu_ub), lambda=setup$lambda),
                                 simplify = FALSE))
    ue <- list()
    ue$setup <- setup
    sims <- MakeSim(setup=setup, time=time_ue)
    ue$sim_m1 <- sims$sim_m1
    ue$sim_m2 <- sims$sim_m2
    ue$sim_m3 <- sims$sim_m3
    ue$sim_m1$dfs <- lapply(ue$sim_m1$dfs, as.data.table)
    ue$sim_m2$dfs <- lapply(ue$sim_m2$dfs, as.data.table)
    ue$sim_m3$dfs <- lapply(ue$sim_m3$dfs, as.data.table)
    save(ue, file=ff_ue)
    rm(ue)
  }
}
RunSim <- function(params, t, ...) {
  n <- length(t) - 1 
  q <- outer(t, t, FUN = pmin)
  iden <- diag(n+1) 
  Sigma <- params[2] * q +
    params[3] * (t %*% t(t))+ 
    params[4] * iden
  mu <- params[1] * t
  return(mvrnorm(mu=mu, Sigma=Sigma))
}
MakeSim <- function(setup, time, ...) {
  
  sims <- list()
  sims$sim_m3 <- list()
  sims$sim_m2 <- list()
  sims$sim_m1 <- list()
  sims$sim_m3$dfs <- list()
  sims$sim_m2$dfs <- list()
  sims$sim_m1$dfs <- list()
  
  params_m3 <- c(setup$mu,
                 setup$sigma_sq,
                 setup$sigma_sq_mu,
                 setup$sigma_sq_ep)
  params_m2 <- c(setup$mu,
                 setup$sigma_sq,
                 setup$sigma_sq_mu,
                 0)
  params_m1 <- c(setup$mu,
                 setup$sigma_sq,
                 0,
                 0)
  
  max_fu <- nrow(time)
  for (k in 1:setup$n_sims) {
    Y_k_m3 <- matrix(nrow=max_fu, ncol = setup$n_subjs)
    Y_k_m2 <- matrix(nrow=max_fu, ncol = setup$n_subjs)
    Y_k_m1 <- matrix(nrow=max_fu, ncol = setup$n_subjs)
    for (i in 1:setup$n_subjs) {
      time_i <- as.numeric(na.omit(time[,i]))
      Y_k_m3[1:length(time_i),i] <- RunSim(params = params_m3, t=time_i)
      Y_k_m2[1:length(time_i),i] <- RunSim(params = params_m2, t=time_i)
      Y_k_m1[1:length(time_i),i] <- RunSim(params = params_m1, t=time_i)
    }
    sims$sim_m3$dfs[[k]] <- AddDiff(MakeLong(Y_k_m3, time), Y='Y', Patient='Patient', t='t')
    sims$sim_m2$dfs[[k]] <- AddDiff(MakeLong(Y_k_m2, time), Y='Y', Patient='Patient', t='t')
    sims$sim_m1$dfs[[k]] <- AddDiff(MakeLong(Y_k_m1, time), Y='Y', Patient='Patient', t='t')
  }
  
  return(sims)
}

#joint likelihood estimators using LGBF-S from optim
mc_optim_ll  <- function(simdat, setup, modtype, ...) {
  num_cores <- detectCores() - 1
  cl <- makeCluster(num_cores)
  clusterExport(cl, c("ll", "optim_ll", 
                      "simdat", "modtype"), envir = environment())
  batch_size <- setup$batch_size
  ests <- NULL
  print(paste('simulations', modtype, ' optim running'))
  for (i in seq(1, setup$n_sims, by = batch_size)) {
    print(i)
    indices <- i:min(i + batch_size - 1, setup$n_sims)
    results <- pblapply(indices, function(idx) {
      optim_ll(simdat$dfs[[idx]], modtype)
    }, cl = cl)
    ests <- rbind(ests, do.call(rbind, results))
    gc()
  }
  stopCluster(cl)
  gc()
  
  stats <- list()
  stats$mu <- unlist(ests[,1])
  stats$sigma_sq <- unlist(ests[,2])
  
  if (modtype=='m2') {
    stats$sigma_sq_mu <- unlist(ests[,3])
  }
  
  if (modtype=='m3') {
    stats$sigma_sq_mu <- unlist(ests[,3])
    stats$sigma_sq_ep <- unlist(ests[,4])
  }
  return(stats)
}

#profile likelihood estimators - MLEs
mc_mle_m1 <- function(sim_m1, ...) {
  print('m1 mle simulations running')
  ests <- pblapply(sim_m1$df, mle_m1)
  stats_m1 <- list()
  stats_m1$mu <- sapply(ests, function(x) x$mu)
  stats_m1$sigma_sq <- sapply(ests, function(x) x$sigma_sq)
  return(stats_m1)
}
mc_mle_m2 <- function(sim_m2, setup, ...) {
  
  num_cores <- detectCores() - 1
  cl <- makeCluster(num_cores)
  clusterExport(cl, c("mle_m2", "mu_mle_m2", 
                      "sim_m2", "profileL_mle_m2"), envir = environment())
  batch_size <- setup$batch_size
  ests <- NULL
  print('m2 mle simulations running')
  for (i in seq(1, setup$n_sims, by = batch_size)) {
    print(i)
    indices <- i:min(i + batch_size - 1, setup$n_sims)
    results <- pblapply(indices, function(idx) {
      mle_m2(sim_m2$dfs[[idx]])
    }, cl = cl)
    ests <- rbind(ests, do.call(rbind, results))
    gc()
  }
  stopCluster(cl)
  gc()
  stats_m2 <- list()
  stats_m2$mu <- unlist(ests[,1])
  stats_m2$sigma_sq <- unlist(ests[,2])
  stats_m2$sigma_sq_mu <- unlist(ests[,3])
  return(stats_m2)
}
mc_mle_m3 <- function(sim_m3, setup, ...) {
  
  print('m3 mle simulations running')
  num_cores <- detectCores() - 1
  cl <- makeCluster(num_cores)
  clusterExport(cl, c("mle_m3", "psi_ll_mle_m3", 
                      "sigma_sq_mu_mle_m3", "sigma_sq_mle_m3", "mu_mle_m3",
                      "sim_m3"), envir = environment())
  batch_size <- setup$batch_size
  ests <- NULL
  for (i in seq(1, setup$n_sims, by = batch_size)) {
    print(i)
    indices <- i:min(i + batch_size - 1, setup$n_sims)
    results <- pblapply(indices, function(idx) {
      mle_m3(df=sim_m3$df[[idx]]) #using default psi value
    }, cl = cl)
    ests <- rbind(ests, do.call(rbind, results))
    gc()
  }
  
  stopCluster(cl)
  gc()
  
  stats_m3 <- list()
  stats_m3$psi <- unlist(ests[,1])
  stats_m3$mu <- unlist(ests[,2])
  stats_m3$sigma_sq <-  unlist(ests[,3])
  stats_m3$sigma_sq_mu <-  unlist(ests[,4])
  stats_m3$sigma_sq_ep <-  unlist(ests[,5])
  return(stats_m3)
}

#profile likelihood estimators - MLEs
mc_u_m2 <- function(sim_m2, ...) {
  print('m2 u simulations running')
  tmp <- lapply(sim_m2$df, u_m2)
  stats_m2 <- list()
  stats_m2$mu <- sapply(tmp, function(x) x$mu)
  stats_m2$sigma_sq <- sapply(tmp, function(x) x$sigma_sq)
  stats_m2$sigma_sq_mu <- sapply(tmp, function(x) x$sigma_sq_mu)
  return(stats_m2)
}
mc_u_m3 <- function(sim_m3, setup, ...) {
  print('m3 u simulations running')
  num_cores <- detectCores() - 1
  cl <- makeCluster(num_cores)
  clusterExport(cl, c("u_m3", "psi_ll_u_m3", 
                      "sigma_sq_mu_u_m3", "sigma_sq_u_m3", "mu_u_m3",
                      "sim_m3"), envir = environment())
  batch_size <- setup$batch_size
  ests <- NULL
  for (i in seq(1, setup$n_sims, by = batch_size)) {
    print(i)
    indices <- i:min(i + batch_size - 1, setup$n_sims)
    results <- pblapply(indices, function(idx) {
      u_m3(sim_m3$df[[idx]])
    }, cl = cl)
    ests <- rbind(ests, do.call(rbind, results))
    gc()
  }
  
  stopCluster(cl)
  gc()
  
  stats_m3 <- list()
  stats_m3$phi <- unlist(ests[,1])
  stats_m3$mu <- unlist(ests[,2])
  stats_m3$sigma_sq <-  unlist(ests[,3])
  stats_m3$sigma_sq_mu <-  unlist(ests[,4])
  stats_m3$sigma_sq_ep <-  unlist(ests[,5])
  return(stats_m3)
}