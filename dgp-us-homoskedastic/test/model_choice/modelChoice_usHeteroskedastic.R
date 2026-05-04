## TV US DiD Fit onto CV US DiD ##
## Model Selection ## 

rm(list=ls())

## Load in Libraries ##
## ----------------------------------------------------------------------------- 

library(loo)

### Set Working Directory ### 
## ----------------------------------------------------------------------------- 

setwd("/Users/tinaarputhasamy/Desktop/spatialDiD_simulations/dgp-us-homoskedastic/test")

## Read in Posterior Samples ##
## ----------------------------------------------------------------------------- 

samps <- read.table("../outfiles/jagsOutput/jagsOut_usHeteroskedastic2.txt",
                    header = T,
                    check.names = F)

## Read in Data ##
## ----------------------------------------------------------------------------- 

data <- read.table("../data/simData_usHomoskedastic.txt", 
                   header = T)

## Calculate WAIC ##
## ----------------------------------------------------------------------------- 

# Get log-likelihood matrix #

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

# Function to compute random component of mu
compute.random.mu <- function(i, omega, region.vec, cov) {
  # Vectorized multiplication for each row `i`
  omega[i, region.vec] * cov
}

# Apply the function across all rows using sapply
random.mu <- t(sapply(seq_len(nrow(omega.df)), 
                      function(i) compute.random.mu(i, omega.df, 
                                                    data$region,
                                                    treatXtime)))

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

waic <- loo(log.likelihood.mat, is_method = "psis")

# Computed from 30000 by 2000 log-likelihood matrix.
# 
# Estimate   SE
# elpd_loo  -3556.7 32.0
# p_loo        40.0  4.0
# looic      7113.4 64.1
# ------
#   MCSE of elpd_loo is NA.
# MCSE and ESS estimates assume independent draws (r_eff=1).
# 
# Pareto k diagnostic values:
#   Count Pct.    Min. ESS
# (-Inf, 0.7]   (good)     1999  100.0%  368     
# (0.7, 1]   (bad)         1    0.0%  <NA>    
#   (1, Inf)   (very bad)    0    0.0%  <NA>    
#   See help('pareto-k-diagnostic') for details.

## -----------------------------------------------------------------------------

# using n.samp = 20000, n.chains = 3 n.adapt = 1000, burnin = 2000

# Computed from 60000 by 2000 log-likelihood matrix.
# 
# Estimate   SE
# elpd_loo  -3556.9 32.0
# p_loo        40.1  4.0
# looic      7113.9 64.0
# ------
#   MCSE of elpd_loo is NA.
# MCSE and ESS estimates assume independent draws (r_eff=1).
# 
# Pareto k diagnostic values:
#   Count Pct.    Min. ESS
# (-Inf, 0.7]   (good)     1999  100.0%  919     
# (0.7, 1]   (bad)         1    0.0%  <NA>    
#   (1, Inf)   (very bad)    0    0.0%  <NA>    
#   See help('pareto-k-diagnostic') for details.

## Calculate DIC ##
## -----------------------------------------------------------------------------

# dic <- dic.samples(m1, n.iter = NITER, type = "pD")