################################################################################################
# NAME: sim.r
# AUTHORS: Ellie Matthay, Catherine Li, Chris Rowe
# DATE STARTED: 2/21/2018
# PURPOSE: Script to run the various case-control simulations
# UPDATES: 4/17/2018: CR updated code to write individual performance measures to csv; updated
#                     exposure/outcome inputs from "A.5" to "A.50" and time input from "time"
#                     to "time.Y.02.A.50"; updated output files to include exposure and outcome
#                     frequencies in names; added comments for new sampling schemes
#          4/23/2018: CL added code to source in ccwc.weights function
#          4/24/2018: Updated arguments to include svysize and method
#          5/1/2018:  ECM minor edits to make code compatible with cluster
#          11/7/2018: CR added svycase study.r argument throughout
#          5/22/2019: CR added stratified.bias survey design 
#          1/13/2020: CR added matching functionality
#          5/19/2020: Added population argument ("standard", "strong_c", "single_c")
################################################################################################

# Clear workspace
rm(list=ls())

#options(warn = 2, show.error.messages=T) # To help with debugging

######
# Load Packages and functions
######

library("parallel") # For the random seed generation (see below)
library("Epi") # for case-cohort and density sampling designs
library("survival") # for clogit analysis
library("dplyr") # for data management operations

# Source the study, ccwc.weights, sim, and perfomance functions
source('code/study.r')
source('code/ccwc.weights.R')
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
  # "density weights" for density-sampled with weights in control selection
# samp:
  # "srs" for simple random sample
  # "sps" for simple probability sample with know probability of selection for each individual
  # "exp.ps" for probability sample with sampling probability corresponds to exposure mechanism
  # "exp.ps2" for probability sample where exposed have sampling probability of 0.75 and unexposed have sampling probability of 0.25
  # "clustered1" for single stage clustered design
  # "clustered2" for two-stage clustered design 
  # "stratified.bias" for stratified design with bias
  # "stratified" for single stage stratified design
  # "age.stratified" for age stratified design
  # "race.stratified" for race stratified design
# svysize:
  # "benchmark" sufficient respondents to allow for unique controls
  # "small" 1000 respondents
  # "medium" 5000 respondents
  # "large" 15000 respondents
# svycase:
  # FALSE: cases are excluded from survey from which controls are obtained
  # TRUE: cases are included in survey from which controls are obtained
# method (NOTE: Only relevant for svysize="small", "medium", or "large"):
  # "expand" expand the sampled controls by their survey weights
  # "sample" sample from controls with replacement and with sampling probability proportional to weights
# match_var (NOTE: Only relevant for cctype = "density" and method = "expand")
  # input should be either a single variable (match_var =  "race"), or multiple (match_var = "list(race, educ, age, sex)")
  # always a character class object, if no matching is desired, set match_var = NULL
# population: 
  # "standard" original population generated for case-control simulations, multiple weak confounders
  # "strong_c" modified population generated for matching scenarios, contains multiple weak confounders and one very strong confounder (sex)
  # "single_c" modified population generated for matching scenarios, contains single strong confounder (sex)



sim1 <- sim(nsims=1, cluster=F, cctype="density", samp="srs", svysize="benchmark", svycase=FALSE, method="expand", ratio=4, 
            match_var = "sex", exposure="A.50", outcome="Y.02.A.50", timevar="time.Y.02.A.50", population = "single_c")
save(sim1, file = "results/cumulative_SRS_bm_ratio4_A50_Y02_1000sims.rdata")
perf1 <- performance(sim1)
write.csv(perf1, "results/perf_cumulative_SRS_bm_ratio4_A50_Y02_1000sims.csv", row.names=F)

sim2 <- sim(nsims=1000, cluster=T, cctype="cumulative", samp="sps", svysize="benchmark", svycase=FALSE, method="sample", ratio=4, 
            match_var = NULL,  exposure="A.50", outcome="Y.02.A.50", timevar="time.Y.02.A.50", population = "standard")
save(sim2, file = "results/cumulative_SPS_bm_ratio4_A50_Y02_1000sims.rdata")
perf2 <- performance(sim2)
write.csv(perf2, "results/perf_cumulative_SPS_bm_ratio4_A50_Y02_1000sims.csv", row.names=F)

sim3 <- sim(nsims=1000, cluster=T, cctype="cumulative", samp="clustered1", svysize="benchmark", svycase=FALSE, method="sample", ratio=4,
            match_var = NULL,  exposure="A.50", outcome="Y.02.A.50", timevar="time.Y.02.A.50", population = "standard")
save(sim3, file = "results/cumulative_clustered1_bm_ratio4_A50_Y02_1000sims.rdata")
perf3 <- performance(sim3)
write.csv(perf3, "results/perf_cumulative_clustered1_bm_ratio4_A50_Y02_1000sims.csv", row.names=F)

sim4 <- sim(nsims=1000, cluster=T, cctype="cumulative", samp="clustered2", svysize="benchmark", svycase=FALSE, method="sample", ratio=4, 
            match_var = NULL,  exposure="A.50", outcome="Y.02.A.50", timevar="time.Y.02.A.50", population = "standard")
save(sim4, file = "results/cumulative_clustered2_bm_ratio4_A50_Y02_1000sims.rdata")
perf4 <- performance(sim4)
write.csv(perf4, "results/perf_cumulative_clustered2_ratio4_A50_Y02_1000sims.csv", row.names=F)

sim5 <- sim(nsims=1000, cluster=T, cctype="cumulative", samp="stratified", svysize="small", svycase=FALSE, method="sample", ratio=4, 
            match_var = NULL,  exposure="A.50", outcome="Y.02.A.50", timevar="time.Y.02.A.50", population = "standard")
save(sim5, file = "results/cumulative_stratified_bm_ratio4_A50_Y02_1000sims.rdata")
perf5 <- performance(sim5)
write.csv(perf5, "results/perf_cumulative_stratified_bm_ratio4_A50_Y02_1000sims.csv", row.names=F)

######
## Close
######

# Close MPI. This is only for the cluster. try() tells r to ignore it and continue if it fails (which it will, if running on a local computer)
try(mpi.quit(), silent=T)

# END