## CV Classical DiD Fit onto CV Classical ##
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

samps <- read.table("../outfiles/jagsOutput/jagsOut_clsHeteroskedastic.txt",
                    header = T,
                    check.names = F)

## Read in Data ##
## ----------------------------------------------------------------------------- 

data <- read.table("../data/simData_clsHomoskedastic.txt", 
                   header = T)

## Calculate WAIC ##
## ----------------------------------------------------------------------------- 

# get log likelihood matrix #

extract <- as.matrix(samps)

beta.df <- as.data.frame(extract[, 
                                 colnames(extract)[grep("beta", 
                                                        colnames(extract))]])
colnames(beta.df) <- c("beta0", "beta1", "beta2", "beta3")

treatXtime <- as.numeric(data$time) * as.numeric(data$treatment)

mu <- sapply(seq_along(data$treatment), function(x) 
  beta.df$beta0 + beta.df$beta1 * as.numeric(data$treatment[x]) + 
    beta.df$beta2 * as.numeric(data$time[x]) + beta.df$beta3 * 
    as.numeric(treatXtime[x]))

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

# fill in pre-index values into log-likelihood matrix

log.likelihood.mat[,indexPre] <- t(sapply(seq_len(nrow(beta.df)), function(i){
  
  calc.log.likelihood(data$y[indexPre], mu[i, indexPre],
                      tausq.pre[i])
}, simplify = "matrix"))

# fill in post-index values into log-likelihood matrix

log.likelihood.mat[,indexPost] <- t(sapply(seq_len(nrow(beta.df)), function(i){
  
  calc.log.likelihood(data$y[indexPost], mu[i, indexPost],
                      tausq.post[i])
}, simplify = "matrix"))

# use log-likelihood to calculate waic

waic <- loo(log.likelihood.mat, is_method = "psis")

# Computed from 30000 by 2000 log-likelihood matrix.
# 
# Estimate   SE
# elpd_loo  -3535.9 31.2
# p_loo         5.9  0.2
# looic      7071.9 62.4
# ------
#   MCSE of elpd_loo is 0.0.
# MCSE and ESS estimates assume independent draws (r_eff=1).
# 
# All Pareto k estimates are good (k < 0.7).
# See help('pareto-k-diagnostic') for details.

## Calculate DIC ##
## -----------------------------------------------------------------------------

#dic <- dic.samples(m1, n.iter = NITER, type = "pD")
