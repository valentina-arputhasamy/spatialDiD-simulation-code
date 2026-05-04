## CV SVC DiD Fit onto CV Classical ##
## Model Selection ## 

rm(list=ls())

## Load in Libraries ##
## ----------------------------------------------------------------------------- 

library(loo)

## Set Working Directory ## 
## ----------------------------------------------------------------------------- 

setwd("/Users/tinaarputhasamy/Desktop/spatialDiD_simulations/dgp-cls-homoskedastic/test")

## Read in Posterior Samples ##
## ----------------------------------------------------------------------------- 

samps <- read.table("../outfiles/jagsOutput/jagsOut_svcHomoskedastic2.txt",
                    header = T,
                    check.names = F)

## Read in Data ##
## ----------------------------------------------------------------------------- 

data <- read.table("../data/simData_clsHomoskedastic.txt", 
                   header = T)

## Calculate WAIC ##
## ----------------------------------------------------------------------------- 

extract <- as.matrix(samps)

beta.df <- as.data.frame(extract[, 
                                 colnames(extract)[grep("beta", 
                                                        colnames(extract))]])
colnames(beta.df) <- c("beta0", "beta1", "beta2", "beta3")

omega.df <- extract[, colnames(extract)[grep("omega", colnames(extract))]]

## separate each of the omega effects into 4 dataframes ##

omega0.indices <- grep("^omega0\\[", colnames(omega.df))
omega0.df <- omega.df[, omega0.indices]

omega1.indices <- grep("^omega1\\[", colnames(omega.df))
omega1.df <- omega.df[, omega1.indices]

omega2.indices <- grep("^omega2\\[", colnames(omega.df))
omega2.df <- omega.df[, omega2.indices]

omega3.indices <- grep("^omega3\\[", colnames(omega.df))
omega3.df <- omega.df[, omega3.indices]

## calculate treatxtime ##

treatXtime <- data$time * data$treatment

## calculate fixed mu component ##

fixed.mu <- sapply(seq_along(data$treatment), function(x) 
  beta.df$beta0 + beta.df$beta1 * as.numeric(data$treatment[x]) + 
    beta.df$beta2 * as.numeric(data$time[x]) + beta.df$beta3 * 
    as.numeric(treatXtime[x]))

## calculate spatial/random mu components ##

# create function for calculating spatial components

compute.random.mu <- function(i, region.vec, omega0, omega1, omega2, omega3, 
                              treat, time, treatxtime) {
  
  # Vectorized multiplication for each row `i`
  
    omega0[i, region.vec] + 
    omega1[i, region.vec]*treat + 
    omega2[i, region.vec]*time +
    omega3[i, region.vec]*treatxtime
}

random.mu <- t(sapply(seq_len(nrow(omega.df)), function(i){
  
  compute.random.mu(i, data$region, omega0.df, omega1.df, omega2.df, 
                    omega3.df, data$treatment, data$time, treatXtime)
  
}))

# combine the components #

total.mu <- fixed.mu + random.mu

# Extract tausqs

tausq.df <- extract[, "tausq"]

# Create log likelihood calculation function

calc.log.likelihood <- function(y, mu, tausq){
  
  dnorm(y, mu, sqrt(tausq), log = T)
}

# Calculate Log likelihoods

log.likelihood.mat <- t(sapply(seq_len(nrow(beta.df)), function(i){
  
  calc.log.likelihood(data$y, total.mu[i, ], tausq.df[i])
}))


waic <- loo(log.likelihood.mat)

# Computed from 30000 by 2000 log-likelihood matrix.
# 
# Estimate   SE
# elpd_loo  -3541.5 31.2
# p_loo        48.0  2.1
# looic      7083.1 62.4
# ------
#   MCSE of elpd_loo is 0.0.
# MCSE and ESS estimates assume independent draws (r_eff=1).
# 
# All Pareto k estimates are good (k < 0.7).
# See help('pareto-k-diagnostic') for details.

## Using # n.samp = 50000, n.chains = 3, n.adapt = 2500, burnin = 5000 ## 
## -----------------------------------------------------------------------------

# Computed from 150000 by 2000 log-likelihood matrix.
# 
# Estimate   SE
# elpd_loo  -3538.8 31.2
# p_loo        34.3  1.5
# looic      7077.5 62.5
# ------
#   MCSE of elpd_loo is 0.0.
# MCSE and ESS estimates assume independent draws (r_eff=1).
# 
# All Pareto k estimates are good (k < 0.7).
# See help('pareto-k-diagnostic') for details.

## Calculate DIC ##
## -----------------------------------------------------------------------------

#dic <- dic.samples(m1, n.iter = NITER, type = "pD")
