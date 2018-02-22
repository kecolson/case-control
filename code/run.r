################################################################################################
# NAME: sim.r
# AUTHORS: Ellie Matthay, Catherine Li, Chris Rowe
# DATE STARTED: 2/21/2018
# PURPOSE: Script to run the various case-control simulations
# UPDATES: [date]: XX
################################################################################################

# Clear workspace
rm(list=ls())

######
# Load Packages and functions
######

library("parallel") # For the random seed generation (see below)
library("Epi") # for case-cohort and density sampling designs
library("survival") # for clogit analysis
library("dplyr") # for data management operations

# Source the study, sim, and perfomance functions
source('code/study.r')
source('code/sim.r')
source('code/performance.r')


######
## Run simulations of interest to us
######

# Simple cumulative case-control, simple random sampling, ratio of 1
sim1 <- sim(nsims=100, cluster=T, cctype="cumulative", samp="srs", ratio=1, 
            exposure="A.5", outcome="Y.02.A.5", timevar="time")
sim1$est.lower.upper # view the point estimates and CIs
summary(sim1$results[[1]]$mod) # summarize the model object from the first simulation
summary(sim1$results[[3]]$sample) # summarize the sample data object from the third simulation
sim1$results[[1]]$truth # view the true OR

# Simple density sampled case-control, simple random sampling, ratio of 1
# sim2 <- sim(nsims=3, cluster=F, cctype="density", samp="srs", ratio=1, 
#             exposure="A.5", outcome="Y.02.A.5", timevar="time")

######
## Aggregate/summarize performance metrics
######

performance(sim1)
#performance(sim2)

######
## Save everything and close
######

save(sim1, file = "results/cumulative_SRS_ratio1_test.rdata")
#save(sim2, file = "results/density_SRS_ratio1_test.rdata")

# Close MPI. This is only for the cluster. try() tells r to ignore it and continue if it fails (which it will, if running on a local computer)
try(mpi.quit(), silent=T)

# END