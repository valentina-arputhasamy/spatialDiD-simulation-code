## CV US DiD Fit onto CV SVC ##
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

samps <- read.table("../outfiles/jagsOutput/jagsOut_usHomoskedastic2.txt",
                    header = T,
                    check.names = F)

## Read in Data ##
## ----------------------------------------------------------------------------- 

data <- read.table("../data/simData_svcHomoskedastic", 
                   header = T)

## Calculate WAIC ##
## ----------------------------------------------------------------------------- 

# get log likelihood matrix #

extract <- as.matrix(samps)

beta.df <- as.data.frame(extract[, 
                                 colnames(extract)[grep("beta", 
                                                        colnames(extract))]])
colnames(beta.df) <- c("beta0", "beta1", "beta2", "beta3")

omega.df <- extract[, colnames(extract)[grep("omega", colnames(extract))]]

treatXtime <- as.numeric(data$time) * as.numeric(data$treatment)

fixed.mu <- sapply(seq_along(data$treatment), function(x) 
  beta.df$beta0 + beta.df$beta1 * as.numeric(data$treatment[x]) + 
    beta.df$beta2 * as.numeric(data$time[x]) + beta.df$beta3 * 
    as.numeric(treatXtime[x]))


# function to compute random component of mu

compute.random.mu <- function(i, region, cov) {
  # Vectorized multiplication for each row `i`
  omega.df[i, region] * cov
}

# apply the function across all rows using sapply
random.mu <- t(sapply(seq_len(nrow(beta.df)), 
                      function(i) compute.random.mu(i, data$region,
                                                    treatXtime)))

total.mu <- fixed.mu + random.mu

# extract tausqs

tausq.df <- extract[, "tausq"]

# create log likelihood calculation function

calc.log.likelihood <- function(y, mu, tausq){
  
  dnorm(y, mu, sqrt(tausq), log = T)
}

# create Empty Log-Likelihood Matrix

log.likelihood.mat <- matrix(nrow = nrow(samps), ncol = length(data$y))

# calculate Log likelihoods

log.likelihood.mat <- t(sapply(seq_len(nrow(beta.df)), function(i){
  
  calc.log.likelihood(data$y, total.mu[i, ], tausq.df[i])
}))

# use log-likelihood to calculate waic

waic <- loo(log.likelihood.mat)

# Computed from 30000 by 2000 log-likelihood matrix.
# 
# Estimate   SE
# elpd_loo  -4579.1 34.7
# p_loo        25.7  2.6
# looic      9158.2 69.4
# ------
#   MCSE of elpd_loo is NA.
# MCSE and ESS estimates assume independent draws (r_eff=1).
# 
# Pareto k diagnostic values:
#   Count Pct.    Min. ESS
# (-Inf, 0.7]   (good)     1994  99.7%   815     
# (0.7, 1]   (bad)         6   0.3%   <NA>    
#   (1, Inf)   (very bad)    0   0.0%   <NA>    
#   See help('pareto-k-diagnostic') for details.

## -----------------------------------------------------------------------------

## Using 1.08, n.samp = 40000, n.chains = 3, n.adapt = 2000, burnin = 4000

# Computed from 120000 by 2000 log-likelihood matrix.
# 
# Estimate   SE
# elpd_loo  -4579.5 34.7
# p_loo        26.0  2.7
# looic      9159.1 69.4
# ------
#   MCSE of elpd_loo is NA.
# MCSE and ESS estimates assume independent draws (r_eff=1).
# 
# Pareto k diagnostic values:
#   Count Pct.    Min. ESS
# (-Inf, 0.7]   (good)     1994  99.7%   1756    
# (0.7, 1]   (bad)         6   0.3%   <NA>    
#   (1, Inf)   (very bad)    0   0.0%   <NA>    
#   See help('pareto-k-diagnostic') for details.

## Calculate DIC ##
## -----------------------------------------------------------------------------

# dic <- dic.samples(m1, n.iter = NITER, type = "pD")
