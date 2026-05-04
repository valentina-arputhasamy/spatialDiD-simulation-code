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

samps <- read.table("../outfiles/jagsOutput/jagsOut_clsHomoskedastic.txt",
                    header = T,
                    check.names = F)

## Read in Data ##
## ----------------------------------------------------------------------------- 

data <- read.table("../data/simData_clsHomoskedastic.txt", 
                   header = T)

## Calculate WAIC ##
## ----------------------------------------------------------------------------- 

## get log-likelihood matrix ##

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

## extract tausqs ##

tausq.df <- extract[, "tausq"]

## create log likelihood calculation function ##

calc.log.likelihood <- function(y, mu, tausq){
  
  dnorm(y, mu, sqrt(tausq), log = T)
}

## create empty log-likelihood matrix ##

log.likelihood.mat <- matrix(nrow = nrow(samps), ncol = length(data$y))

## calculate log-likelihoods ##

log.likelihood.mat <- t(sapply(seq_len(nrow(beta.df)), function(i){
  
  calc.log.likelihood(data$y, mu[i, ], tausq.df[i])
}))

## use log-likelihood to calculate waic ##

waic <- loo(log.likelihood.mat)

# Computed from 30000 by 2000 log-likelihood matrix.
# 
# Estimate   SE
# elpd_loo  -3536.2 31.2
# p_loo         4.9  0.2
# looic      7072.4 62.3
# ------
#   MCSE of elpd_loo is 0.0.
# MCSE and ESS estimates assume independent draws (r_eff=1).
# 
# All Pareto k estimates are good (k < 0.7).
# See help('pareto-k-diagnostic') for details.

## Calculate DIC ##
## -----------------------------------------------------------------------------

#dic <- dic.samples(m1, n.iter = NITER, type = "pD")
