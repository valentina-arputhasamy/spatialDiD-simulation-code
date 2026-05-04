## CV Classical DiD onto CV SVC ##
## Model Fit ## 


rm(list=ls())

## Load in Libraries ##
## ----------------------------------------------------------------------------- 

library(rjags)
library(coda)
library(loo)
library(tidyverse)
library(sf)
library(maps)
library(tigris)
library(spdep)

### Set Working Directory ### 
## ----------------------------------------------------------------------------- 

setwd("/Users/tinaarputhasamy/Desktop/spatialDiD_simulations/dgp-svc-homoskedastic/test")

## Read in Data Prep Files ## 
## ----------------------------------------------------------------------------- 

source("../src/dataPrep_svcHomoskedastic.R")

## Set Initial Values ##
## ----------------------------------------------------------------------------- 

## set initial values for chain 1 ## 

beta.init.1 = rnorm(q, 0, 0.01)
tausq.inv.init.1 = rep(1, times = 1)

initial.values.chain.1 = list(beta = beta.init.1, 
                              tausq.inv=tausq.inv.init.1) 

## set initial values for chain 2 ## 

beta.init.2 = rnorm(q, 0, 0.01)
tausq.inv.init.2 <- rep(1, times = 1)

initial.values.chain.2 = list(beta = beta.init.2, 
                              tausq.inv=tausq.inv.init.2)

## set initial values for chain 3 ## 

beta.init.3 = rnorm(q, 0, 0.01)
tausq.inv.init.3 <- rep(1, times = 1)

initial.values.chain.3 = list(beta = beta.init.2, 
                              tausq.inv=tausq.inv.init.2)

## combine initial values into a single list ##

initial.values.all = list(initial.values.chain.1, 
                          initial.values.chain.2,
                          initial.values.chain.3)

## Set MCMC Settings ##
## ----------------------------------------------------------------------------- 

n.samp = 10000
n.chains = 3
n.adapt = 500
burnin = 1000

## Assign Monitored Parameters ##
## ----------------------------------------------------------------------------- 

model.parameters = c("beta", "tausq.inv", "tausq")

## Load in Model File ##
## ----------------------------------------------------------------------------- 

model.file <- "../src/modelFile_clsHomoskedastic.txt"

## Compile Model ##
## ----------------------------------------------------------------------------- 

m1 <- jags.model(model.file, data = lt.Data.cls, inits = initial.values.all,
                 n.chains = n.chains, n.adapt = n.adapt)

## Posterior Sampling ##
## ----------------------------------------------------------------------------- 

#set.seed(327)
set.seed(712)
update(m1, n.iter = burnin)
m1.out <- coda.samples(model = m1, variable.names = model.parameters,
                       n.iter = n.samp)

## Posterior Summary Statistics ##
## -----------------------------------------------------------------------------

sm <- summary(m1.out)
chains <- as.mcmc.list(m1.out)
# Gelman–Rubin PSRFs
print(gelman.diag(chains))
# Multivariate psrf
# 1, n.samp = 10000, n.chains = 3, n.adapt = 500, burnin = 1000

samps <- do.call(rbind, m1.out)

# calculate mcse (another measure of convergence) #

# manual calc 
mcmc_obj <- as.mcmc(as.matrix(samps)) 
# turn your matrix into an mcmc object
ess_vec <- effectiveSize(mcmc_obj) 
# get posterior SD per parameter
sd_vec  <- apply(as.matrix(samps), 2, sd)
# compute MCSE per parameter
mcse_vec <- sd_vec / sqrt(ess_vec)
# summarize
mcse_vec
# median mcse
median_mcse_manual <- median(mcse_vec)
# 0.001302808

# # using posterior package 
# mcse_all <- mcse_mean(draws)      # named numeric vector, length 6
# mcse_all <- apply(draws, 2, function(x) mcse_mean(x))
# median_mcse <- median(mcse_all)
# # 0.001316831

# get effective sample size using coda package
draws <- as_draws_matrix(as.matrix(samps))
ess <- apply(draws, 2, function(x) effectiveSize(x))
median(ess)
# 3471.141

## Save Full Posterior Samples ##
## ----------------------------------------------------------------------------- 

write.table(samps, 
            file = "../outfiles/jagsOutput/jagsOut_clsHomoskedastic.txt", 
            row.names=FALSE, col.names=TRUE)
