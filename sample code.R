# Here are a few examples of implementing the WPDM.R function to run a Wiener 
# process degradation model.

source('functions/WPDM.R')
load('data/example_data.RData')

# model 3 with L-BFGS optimizer (default)
WPDM(df=example_data$m3, idvar='Patient', timevar='t', Yvar='Y')

# model 2 with profile likelihood using MLE (params0=(sigma_sq, sigma_sq_mu))
WPDM(df=example_data$m2, idvar='Patient', timevar='t', Yvar='Y', 
     estimator='mle', modtype='m2', params0=c(0.1, 0.1))

# model 3 with profile likelihood using empirically unbiased estimator (params0=psi0)
WPDM(df=example_data$m3, idvar='Patient', timevar='t', Yvar='Y', 
     estimator='u', modtype='m3', params0=0.1)
