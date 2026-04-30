### Generate Homoskedastic Classical DiD Data ###

### Set Seed and Load Libraries ###
## ----------------------------------------------------------------------------- 

set.seed(1234)
library(maps)
library(mapdata)
library(maptools)
library(spdep)
library(tidyverse)
library(maptools)
library(readxl)
library(MASS)
library(Matrix)

### Set Working Directory ### 
## ----------------------------------------------------------------------------- 

setwd("/Users/tinaarputhasamy/Desktop/spatialDiD_simulations/dgp-cls-homoskedastic/test")

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
county.names
# county names

### Get Population Count by County ###
## ----------------------------------------------------------------------------- 

cpop_data <- read_excel("../data/ca_county_pop.xlsx")

cpop_data_2018 <- cpop_data %>% 
  filter(Year == 2018) %>%
  dplyr::group_by(County) %>%
  dplyr::summarize(pop_total = sum(Population)) %>% 
  mutate(pop_prop = pop_total/sum(pop_total)) 

### Get Population Count by County ###
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
remaining.assignments <- sample(1:c, size = n.subjects - c, prob = probs, 
                                replace = TRUE)
# fills in remaining assignments based on probability vector
region <- c(initial.assignments, remaining.assignments) 
# combine assignment lists
region <- sample(region)
# shuffle assignment list
region <- rep(region, each = n.obs)
# repeating each region twice for data set

## NOTE: assumes each county is represented in sample

### Assign Values to Fixed Parameters  ### 
## ----------------------------------------------------------------------------- 

beta0.gen <- 2
beta1.gen <- 4
beta2.gen <- 7
beta3.gen <- 12
beta.gen.vec <- c(beta0.gen, beta1.gen, beta2.gen, beta3.gen)

tausq.gen <- 2 
# variance of epislon

### Create Design Matrix ###
## ----------------------------------------------------------------------------- 

prob <- 0.5 
# probability of being in treatment group
treatment <- rep(rbinom(n.subjects, 1, prob), each = n.obs) 
# treatment assignment
time <- rep(c(1,0), times = n.subjects)
# time period assignment
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

mu.vec <- as.matrix(X) %*% beta.gen.vec

set.seed(123)

y <- sapply(mu.vec, function(x) rnorm(1, x, sqrt(tausq.gen)))

### Generate Components of Spatial Random Effects Vector For US/SVC Testing ###
## ----------------------------------------------------------------------------- 

W <- ca.adj.mat
# adjacency matrix
D <- diag(c(W %*% rep(1, nrow(W)))) 
# number of neighbors for each county
phi <- 0.99
car.precision.unscaled <- D - phi*W

## Save Files ## 
## ----------------------------------------------------------------------------- 

df <- data.frame(cbind(y, X, region))
write.table(df, file = "../data/simData_clsHomoskedastic.txt", 
            row.names = FALSE)
write.table(W, file = "../data/adjMatrix_clsHomoskedastic", row.names = TRUE)
