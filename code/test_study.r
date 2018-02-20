################################################################################################
# NAME: test_study.r
# AUTHORS: Ellie Matthay, Catherine Li, Chris Rowe
# DATE STARTED: 12/15/2017    
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
library(dplyr)
library(geepack) # for modified poisson

# Set Working Directory
#setwd("~/Documents/PhD/Ahern GSR/Case Control Simulation") # Chris's directory
#setwd("C:/Users/kecolson/Google Drive/simulation/case-control-other") # Ellie's directory
setwd("C:/Users/Catherine/Desktop/Case Control GSR") # Catherine's directory

# Bring in data and true parameters
pop <- read.csv("data/population.data.csv", stringsAsFactors = F)

# Sample down the dataset for testing purposes
#data <- pop[sample(1:nrow(pop), 4000, replace=F),]

# Full dataset for testing purposes
data <- pop

####
## Test function
####

# Warning: density sampling a long time to run.

study1 <- study(seed=1, cctype = "cumulative", samp = "srs", ratio = 1, data=data)
round(c(study1$est, study1$lower, study1$upper),3) 

study2 <- study(seed=1, cctype = "density", samp = "srs", ratio = 1, data=data) # warning has tied failure times. Look into this/check with Patrick.
round(c(study3$est, study3$lower, study3$upper),3) 

# END