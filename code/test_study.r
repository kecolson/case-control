################################################################################################
# NAME: test_study.r
# AUTHORS: Ellie Matthay, Catherine Li, Chris Rowe
# DATE STARTED: 2/20/18
# PURPOSE: This script is to allow for testing of the study() function, which applies
#            various survey sampling and case control procedures to 
#            test the performance of different approaches to selecting controls from external 
#            data sources with complex sampling methods.
# UPDATES: [date]: XX
################################################################################################

# Clear workspace
rm(list=ls())

# Load Packages. Be sure your local computer has these installed before running. No need to re-install each time you run. 
library(Epi) # for case-cohort and density sampling designs
library(survival) # for clogit analysis
library(dplyr) # for data management
library(parallel) # for setting random seeds

# Set Working Directory
#setwd("~/Documents/PhD/Ahern GSR/Case Control Simulation") # Chris's directory
#setwd("C:/Users/kecolson/Google Drive/simulation/case-control-other") # Ellie's directory
setwd("C:/Users/Catherine/Documents/GitHub/case-control_master") # Catherine's directory

# Source the study function
source('code/study.r')

# Bring in data and true parameters
pop <- read.csv("data/population.data.csv", stringsAsFactors = F)

# Sample down the dataset for testing purposes
#data <- pop[sample(1:nrow(pop), 4000, replace=F),]

# Full dataset for testing purposes
data <- pop

# Generate random seed to feed into functino
RNGkind("L'Ecuyer-CMRG")
set.seed(1)
s <- .Random.seed
seeds <- list(s)

# Select names of exposure, outcome, and time variables to study
outcome <- "Y.02.A.5"
exposure <- "A.5"
timevar <- "time"

####
## Test function
####

# Warning: density sampling a long time to run.

study1 <- study(iteration=1, cctype = "cumulative", samp = "srs", ratio = 1, data=data,
                exposure=exposure, outcome=outcome, timevar=timevar, seeds=seeds)
round(c(study1$est, study1$lower, study1$upper, study1$truth),3) 


study2 <- study(iteration=1, cctype = "density", samp = "srs", ratio = 1, data=data,
                exposure=exposure, outcome=outcome, timevar=timevar, seeds=seeds) 
round(c(study2$est, study2$lower, study2$upper, study2$truth),3) 

# END