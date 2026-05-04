## TV Classical DiD onto CV Classical ##
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

setwd("/Users/tinaarputhasamy/Desktop/spatialDiD_simulations/dgp-cls-homoskedastic/test/")

## Read in Data Prep Files ## 
## ----------------------------------------------------------------------------- 

source("../src/dataPrep_clsHomoskedastic.R")

## Set Initial Values ##
## ----------------------------------------------------------------------------- 

## set initial values for chain 1 ## 

beta.init.1 = rnorm(q, 0, 0.01)
tausq.inv.init.1 = rep(1, times = 2)

initial.values.chain.1 = list(beta = beta.init.1, 
                              tausq.inv=tausq.inv.init.1) 

## set initial values for chain 2 ## 

beta.init.2 = rnorm(q, 0, 0.01)
tausq.inv.init.2 <- rep(1, times = 2)

initial.values.chain.2 = list(beta = beta.init.2, 
                              tausq.inv=tausq.inv.init.2)

## set initial values for chain 3 ## 

beta.init.3 = rnorm(q, 0, 0.01)
tausq.inv.init.3 <- rep(1, times = 2)

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

model.file <- "../src/modelFile_clsHeteroskedastic.txt"

## Compile Model ##
## ----------------------------------------------------------------------------- 

m1 <- jags.model(model.file, data = lt.Data.cls, inits = initial.values.all,
                 n.chains = n.chains, n.adapt = n.adapt)

## Posterior Sampling ##
## ----------------------------------------------------------------------------- 

set.seed(327)
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
# 
# 1

## Save Full Posterior Samples ##
## ----------------------------------------------------------------------------- 

samps <- do.call(rbind, m1.out)
write.table(samps, 
            file = "../outfiles/jagsOutput/jagsOut_clsHeteroskedastic.txt", 
            row.names=FALSE, col.names=TRUE)