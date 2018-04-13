################################################################################################
# NAME: sim.r
# AUTHORS: Ellie Matthay, Catherine Li, Chris Rowe
# DATE STARTED: 2/21/2018
# PURPOSE: Script to run the various case-control simulations
# UPDATES: [date]: XX
################################################################################################

# Clear workspace
rm(list=ls())

options(warn = 2, show.error.messages=T) # To help with debuggin

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
## Save and aggregate/summarize performance metrics after each simulation, in case the next one fails. 
## 
######

# Options:
# cctype:
# "cumulative"  for cumulative case-control
# "density" for density-sampled
# samp:
# "srs" for simple random sample
# "sps" for simple probability sample with know probability of selection for each individual
# "clustered1" for single stage clustered design
# "clustered2" for two-stage clustered design 
# "stratified" for single stage stratified design

sim1 <- sim(nsims=1000, cluster=T, cctype="cumulative", samp="srs",        ratio=4, exposure="A.5", outcome="Y.02.A.5", timevar="time")
save(sim1, file = "results/cumulative_SRS_ratio4_1000sims.rdata")
performance(sim1)

sim2 <- sim(nsims=1000, cluster=T, cctype="cumulative", samp="sps",        ratio=4, exposure="A.5", outcome="Y.02.A.5", timevar="time")
save(sim2, file = "results/cumulative_SPS_ratio4_1000sims.rdata")
performance(sim2)

sim3 <- sim(nsims=1000, cluster=T, cctype="cumulative", samp="clustered1", ratio=4, exposure="A.5", outcome="Y.02.A.5", timevar="time")
save(sim3, file = "results/cumulative_clustered1_ratio4_1000sims.rdata")
performance(sim3)

sim4 <- sim(nsims=1000, cluster=T, cctype="cumulative", samp="clustered2", ratio=4, exposure="A.5", outcome="Y.02.A.5", timevar="time")
save(sim4, file = "results/cumulative_clustered2_ratio4_1000sims.rdata")
performance(sim4)

sim5 <- sim(nsims=1000, cluster=T, cctype="cumulative", samp="stratified", ratio=4, exposure="A.5", outcome="Y.02.A.5", timevar="time")
save(sim5, file = "results/cumulative_stratified_ratio4_1000sims.rdata")
performance(sim5)


######
## Close
######

# Close MPI. This is only for the cluster. try() tells r to ignore it and continue if it fails (which it will, if running on a local computer)
try(mpi.quit(), silent=T)

# END