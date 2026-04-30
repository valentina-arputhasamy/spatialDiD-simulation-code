## SVC (Time-Varying) DiD Data Preparation for JAGS ##

rm(list=ls())

## Load Libraries ##
## ----------------------------------------------------------------------------- 

library(tidyverse)

### Set Working Directory ### 
## ----------------------------------------------------------------------------- 

setwd("/Users/tinaarputhasamy/Desktop/spatialDiD_simulations/dgp-svc-homoskedastic/test")

## Read in Data ##
## ----------------------------------------------------------------------------- 

data <- read.table("../data/simData_svcHomoskedastic", header = T)
# data

X <- data %>% dplyr::select(-c(region, y))
# design matrix

W <- read.table("../data/adjMatrix_svcHomoskedastic") %>%
  as.matrix()
# adjacency matrix

## Save Some Regression Related Parameters ##
## ----------------------------------------------------------------------------- 

c <- length(unique(data$region))
# number of counties
q <- ncol(X) 
# number of regression parameters
N <- nrow(data)

## Save Some CAR Related Parameters ##
## -----------------------------------------------------------------------------

D <- diag(c(W %*% rep(1, nrow(W)))) 
# number of neighbors for each county
phi <- 0.99
car.precision.unscaled <- D - phi*W

car.mu <- rep(0, c)
# mean vector of car

## Get Lists Ready for Jags ##
## ----------------------------------------------------------------------------- 

lt.Data.svc <- list(y = data$y, trt = X$treatment, time = X$time, car.mu = car.mu, 
                car.precision.unscaled = car.precision.unscaled, q=q, N = N,
                c = c, region = data$region)
# svc tv or cv models

lt.Data.cls <- list(y = data$y, trt = X$treatment, time = X$time, q = q, N = N)
# classical did model

lt.Data.uspatial <- list(y = data$y, trt = X$treatment, time = X$time, 
                         car.mu = car.mu, 
                         car.precision.unscaled = car.precision.unscaled,
                         q = q, N = N, c = c, region = data$region)
# uniform spatial model