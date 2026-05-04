## CV SVC DiD Fit onto CV SVC ##
## Model Selection ## 

rm(list=ls())

## Load in Libraries ##
## ----------------------------------------------------------------------------- 

library(loo)

### Set Working Directory ### 
## ----------------------------------------------------------------------------- 

setwd("/Users/tinaarputhasamy/Desktop/spatialDiD_simulations/dgp-svc-homoskedastic/test")

## Read in Posterior Samples ##
## ----------------------------------------------------------------------------- 

samps <- read.table("../outfiles/jagsOutput/jagsOut_svcHeteroskedastic2.txt",
                    header = T,
                    check.names = F)

## Read in Data ##
## ----------------------------------------------------------------------------- 

data <- read.table("../data/simData_svcHomoskedastic", 
                   header = T)

## Calculate WAIC ##
## ----------------------------------------------------------------------------- 

# Get log likelihood matrix #

extract <- as.matrix(samps)

beta.df <- as.data.frame(extract[, 
                                 colnames(extract)[grep("beta", 
                                                        colnames(extract))]])
colnames(beta.df) <- c("beta0", "beta1", "beta2", "beta3")

omega.df <- extract[, colnames(extract)[grep("omega", colnames(extract))]]

## Separate each of the omega effects into 4 dataframes ##

omega0.indices <- grep("^omega0\\[", colnames(omega.df))
omega0.df <- omega.df[, omega0.indices]

omega1.indices <- grep("^omega1\\[", colnames(omega.df))
omega1.df <- omega.df[, omega1.indices]

omega2.indices <- grep("^omega2\\[", colnames(omega.df))
omega2.df <- omega.df[, omega2.indices]

omega3.indices <- grep("^omega3\\[", colnames(omega.df))
omega3.df <- omega.df[, omega3.indices]

## Calculate treatxtime ##

treatXtime <- data$time * data$treatment

## Calculate fixed mu component ##

fixed.mu <- sapply(seq_along(data$treatment), function(x) 
  beta.df$beta0 + beta.df$beta1 * as.numeric(data$treatment[x]) + 
    beta.df$beta2 * as.numeric(data$time[x]) + beta.df$beta3 * 
    as.numeric(treatXtime[x]))

## Calculate spatial/random mu components ##

# Create function for calculating spatial components

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

# Combine the components #

total.mu <- fixed.mu + random.mu

# Extract tausqs

tausq.df <- extract[, grep("\\btausq.\\b", colnames(extract))]
tausq.pre <- tausq.df[,"tausq[1]"]
tausq.post <- tausq.df[,"tausq[2]"]

# Split observations into pre- and post-exposure groups
indexPre <- which(data$time == 0)    # Pre-exposure
indexPost <- which(data$time == 1)   # Post-exposure

# Create Empty Log-Likelihood Matrix

log.likelihood.mat <- matrix(nrow = nrow(samps), ncol = length(data$y))

# Create log likelihood calculation function

calc.log.likelihood <- function(y, mu, tausq){
  
  dnorm(y, mu, sqrt(tausq), log = T)
}

# Fill in pre-index values into log-likelihood matrix

log.likelihood.mat[,indexPre] <- t(sapply(seq_len(nrow(beta.df)), function(i){
  
  calc.log.likelihood(data$y[indexPre], total.mu[i, indexPre],
                      tausq.pre[i])
}, simplify = "matrix"))

# Fill in post-index values into log-likelihood matrix

log.likelihood.mat[,indexPost] <- t(sapply(seq_len(nrow(beta.df)), function(i){
  
  calc.log.likelihood(data$y[indexPost], total.mu[i, indexPost],
                      tausq.post[i])
}, simplify = "matrix"))

# Use log-likelihood to calculate waic

waic <- loo(log.likelihood.mat, is_method="psis")

# Computed from 30000 by 2000 log-likelihood matrix.
# 
# Estimate   SE
# elpd_loo  -3609.6 30.9
# p_loo       132.5  6.7
# looic      7219.2 61.8
# ------
#   MCSE of elpd_loo is NA.
# MCSE and ESS estimates assume independent draws (r_eff=1).
# 
# Pareto k diagnostic values:
#   Count Pct.    Min. ESS
# (-Inf, 0.7]   (good)     1987  99.4%   667     
# (0.7, 1]   (bad)        12   0.6%   <NA>    
#   (1, Inf)   (very bad)    1   0.0%   <NA>    
#   See help('pareto-k-diagnostic') for details.

## -----------------------------------------------------------------------------

waic(log.likelihood.mat)

# Computed from 210000 by 2000 log-likelihood matrix.
# 
# Estimate   SE
# elpd_waic  -3600.6 30.6
# p_waic       122.4  5.6
# waic        7201.2 61.3

# using waic function, calculated w/, 
# n.samp = 70000, n.chains = 3, n.adapt = 3500, burnin = 7000


## Calculate DIC ##
## -----------------------------------------------------------------------------

#dic <- dic.samples(m1, n.iter = NITER, type = "pD")
