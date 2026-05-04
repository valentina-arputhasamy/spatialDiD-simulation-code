## TV SVC DiD Fit onto CV SVC ##
## Model Fit ## 
## COMPLETELY CORRECTED CODE ## 

rm(list=ls())

## Load in Libraries ##
## ----------------------------------------------------------------------------- 

library(rjags)
library(coda)
library(loo)
library(tidyverse)

## Set Working Directory ## 
## ----------------------------------------------------------------------------- 

setwd("/Users/tinaarputhasamy/Desktop/spatialDiD_simulations/dgp-svc-homoskedastic/test")

## Read in Data and JAGS List ## 
## ----------------------------------------------------------------------------- 

source("../src/dataPrep_svcHomoskedastic.R")

## Set Initial Values ##
## ----------------------------------------------------------------------------- 

## set initial values for chain 1 ## 

beta.init.1 = rnorm(q, 0, 0.01)
tausq.inv.init.1 = rep(1, times = 2)
sigmasq.inv.init.1 = rep(1, q)

initial.values.chain.1 = list(beta = beta.init.1, 
                              tausq.inv=tausq.inv.init.1, 
                              sigmasq.inv=sigmasq.inv.init.1) 

## set initial values for chain 2 ##

beta.init.2 = rnorm(q, 0, 0.01)
tausq.inv.init.2 <- rep(1, times = 2)
sigmasq.inv.init.2 = rep(1, q)

initial.values.chain.2 = list(beta = beta.init.2,
                              tausq.inv=tausq.inv.init.2,
                              sigmasq.inv=sigmasq.inv.init.2)

## set initial values for chain 3 ##

beta.init.3 = rnorm(q, 0, 0.01)
tausq.inv.init.3 <- rep(1, times = 2)
sigmasq.inv.init.3 = rep(1, q)

initial.values.chain.3 = list(beta = beta.init.3,
                              tausq.inv=tausq.inv.init.3,
                              sigmasq.inv=sigmasq.inv.init.3)

## combine initial values into a single list ##

initial.values.all = list(initial.values.chain.1,
                          initial.values.chain.2,
                          initial.values.chain.3)

## Set MCMC Settings ##
## ----------------------------------------------------------------------------- 

n.samp = 70000
n.chains = 3
n.adapt = 3500
burnin = 7000

## Assign Monitored Parameters ##
## ----------------------------------------------------------------------------- 

model.parameters = c("beta", "tausq.inv", "tausq", "sigmasq.inv", "sigsq", 
                     "omega0", "omega1", "omega2", "omega3")

## Load in Model File ##
## ----------------------------------------------------------------------------- 

model.file <- "../src/modelFile_svcHeteroskedastic.txt"

## Compile Model ##
## ----------------------------------------------------------------------------- 

m1 <- jags.model(model.file, data = lt.Data.svc, inits = initial.values.all,
                 n.chains = n.chains, n.adapt = n.adapt)

## Posterior Sampling ##
## ----------------------------------------------------------------------------- 

set.seed(602)

update(m1, n.iter = burnin)
m1.out <- coda.samples(model = m1, variable.names = model.parameters,
                       n.iter = n.samp)

## Posterior Summary Statistics ##
## ----------------------------------------------------------------------------- 

sm <- summary(m1.out)
chains <- as.mcmc.list(m1.out)
# Gelman–Rubin PSRFs
print(gelman.diag(chains))
# 1.06, n.samp = 70000, n.chains = 3, n.adapt = 3500, burnin = 7000

## Save Full Posterior Samples ##
## ----------------------------------------------------------------------------- 

samps <- do.call(rbind, m1.out)
write.table(samps, 
            file = "../outfiles/jagsOutput/jagsOut_svcHeteroskedastic2.txt", 
            row.names=FALSE, col.names=TRUE)
