### Generate Homoskedastic SVC DiD Data ###

### Set Seed and Load Libraries ###
## ----------------------------------------------------------------------------- 

set.seed(1234)

library(maps)
library(mapdata)
#library(maptools)
library(spdep)
library(tidyverse)
library(readxl)
library(MASS)
library(Matrix)

### Set Working Directory ### 
## ----------------------------------------------------------------------------- 

setwd("~/Desktop/Desktop - Tina’s MacBook Pro/spatialDiD_simulations/dgp-svc-homoskedastic/test")

### Load California County Data ###
## ----------------------------------------------------------------------------- 

ca_counties <- maps::map("county", "california", fill=TRUE, plot=FALSE)

### Get Adjacency Matrix ### 
## ----------------------------------------------------------------------------- 

county.ID <- sapply(strsplit(ca_counties$names, ":"), function(x) x[1])
ca.poly = map2SpatialPolygons(ca_counties, IDs=county.ID)
ca.nb = poly2nb(ca.poly)
ca.adj.mat = nb2mat(ca.nb, style="B")
# adjacency matrix (W)
c <- nrow(ca.adj.mat) 
# number of counties
county.names <- sapply(strsplit(county.ID, ","), function(x) x[2])
county.names <- as.factor(str_to_title(county.names))
# county names

### Get Population Count by County ###
## ----------------------------------------------------------------------------- 

cpop_data <- read_excel("../data/ca_county_pop.xlsx")

cpop_data_2018 <- cpop_data %>% 
  filter(Year == 2018) %>%
  dplyr::group_by(County) %>%
  dplyr::summarize(pop_total = sum(Population)) %>% 
  mutate(pop_prop = pop_total/sum(pop_total)) 

### Assign Subject Regions Based on County Proportions ###
## ----------------------------------------------------------------------------- 

n.subjects <- 1000
# number of subjects
n.obs <- 2 
# number of times each subject is observed 
N <- n.subjects * n.obs
# total number of observations
probs <- cpop_data_2018$pop_prop
# proportion of people in each county
initial.assignments <- 1:c
# ensures each county is assigned at least once
set.seed(1234)
remaining.assignments <- sample(1:c, size = n.subjects - c, prob = probs, 
                                replace = TRUE)
# fills in remaining assignments based on probability vector
region <- c(initial.assignments, remaining.assignments) 
# combine assignment lists
set.seed(1234)
region <- sample(region)
# shuffle assignment list
region <- rep(region, each = n.obs)
# repeating each region twice for data set
county.region.assignments <- data.frame(unique(cpop_data_2018$County), 
                                        sort(unique(region))) 
colnames(county.region.assignments) <- c("county", "region")
# creating a county-region dataset

## NOTE: assumes each county is represented in sample

### Generate Components of Spatial Random Effects Vector ###
## ----------------------------------------------------------------------------- 

W <- ca.adj.mat
# adjacency matrix
D <- diag(c(W %*% rep(1, nrow(W)))) 
# number of neighbors for each county
phi <- 0.99
car.precision.raw <- D - phi*W

### Create Spatial Random Effects Vectors ###
## ----------------------------------------------------------------------------- 

## to vary intercept (beta0.gen) by region ##
set.seed(1234)
sigmasq0.gen <- 2
car.sigma0 <- (solve(car.precision.raw))*sigmasq0.gen
omega0.gen <- mvrnorm(1, rep(0, c), car.sigma0)
# vector of spatial random effects for each county for intercept


## to vary coefficient of treatment (beta1.gen) by region ##
set.seed(1234)
sigmasq1.gen <- 3
car.sigma1 <- (solve(car.precision.raw))*sigmasq1.gen
omega1.gen <- mvrnorm(1, rep(0, c), car.sigma1)

## to vary coefficient of time (beta2.gen) by region ##
set.seed(1234)
sigmasq2.gen <- 5
car.sigma2 <- (solve(car.precision.raw))*sigmasq2.gen
omega2.gen <- mvrnorm(1, rep(0, c), car.sigma2)

## to vary coefficient of treatXtime (beta3.gen) by region ##
set.seed(1234)
sigmasq3.gen <- 4
car.sigma3 <- (solve(car.precision.raw))*sigmasq3.gen
omega3.gen <- mvrnorm(1, rep(0, c), car.sigma3)

### Assign Values to Fixed Parameters  ### 
## ----------------------------------------------------------------------------- 

beta0.gen <- 2
beta1.gen <- 4
beta2.gen <- 7
beta3.gen <- 12
beta.gen.vec <- c(beta0.gen, beta1.gen, beta2.gen, beta3.gen)

tausq.gen <- 2
# variance of epsilon

### Create Design Matrix ###
## ----------------------------------------------------------------------------- 

prob <- 0.5 
# probability of treatment
set.seed(1234)
treatment <- rep(rbinom(n.subjects, 1, p = prob), each = n.obs)
# randomize treatment assignment 
time <- rep(c(1,0), n.subjects)
treatXtime <- treatment * time
# interaction term creation
X <- data.frame(intercept = rep(1, n.obs * n.subjects), treatment = treatment,
                time = time, treatXtime = treatXtime)

# design matrix
p <- ncol(X) - 1
# number of predictors in model
q <- ncol(X)
# number of regression coefficients estimated 

## Generate Outcomes ##
## ----------------------------------------------------------------------------- 

mu.spat0 <- sapply(region, function(x) omega0.gen[x])
mu.spat1 <- sapply(region, function(x) omega1.gen[x]) * treatment
mu.spat2 <- sapply(region, function(x) omega2.gen[x]) * time
mu.spat3 <- sapply(region, function(x) omega3.gen[x]) * treatXtime

mu.spat.total <- mu.spat0 + mu.spat1 + mu.spat2 + mu.spat3

mu.vec <- (as.matrix(X) %*% beta.gen.vec) + mu.spat.total

set.seed(123)
y <- sapply(mu.vec, function(x) rnorm(1, x, sqrt(tausq.gen)))

## Save Data and Adjacency Matrix ##
## -----------------------------------------------------------------------------

df <- data.frame(cbind(y, X, region))
write.table(df, file="../data/simData_svcHomoskedastic", row.names=FALSE)
write.table(W, file = "../data/adjMatrix_svcHomoskedastic", row.names = TRUE)
write.table(county.region.assignments, file = "../data/simData_countyRegions",
            col.names = TRUE, row.names = F)
