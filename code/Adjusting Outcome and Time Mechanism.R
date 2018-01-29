################################################################################################
# NAME: Case Control Simulation: Adjusting Outcome and Time Mechanism
# AUTHORS: Ellie Matthay, Catherine Li, Chris Rowe*
# DATE STARTED: 01/22/2018    
# PURPOSE: To play around with exposure and outcome mechanism, and confounding in "true" population
# UPDATES: 01/22/2018: CL adapted code from main Case Control Simulation Program to play around 
#                      with confounding and outcome
################################################################################################

####### PROGRAM START

# Clear workspace
rm(list=ls())

# Load Packages. Be sure your local computer has these installed before running. No need to re-install each time you run. 
library(readstata13)
library(locfit) # for expit function
library(Epi) # for case-cohort and density sampling designs
library(survival) # for clogit analysis
library(dplyr)

# Set Working Directory
#setwd("~/Documents/PhD/Ahern GSR/Case Control Simulation") # Chris's directory
#setwd("C:/Users/kecolson/Google Drive/simulation/case-control") # Ellie's directory
setwd("C:/Users/Catherine/Desktop/Case Control GSR") # Catherine's directory

# Read in population data
pop <- read.csv("data/population.data.csv", stringsAsFactors = F)

N <- 1000000

# Generate Time Variable (X)
set.seed(1)

pop$X <- rnorm(N, 365.25*5, 365.25*1) # Think about other distributions as well

# Generate Exposure
baseline_A <- log(0.5)
set.seed(2)
pop$A <- rbinom(N, size = 1, prob = expit(baseline_A + 
                                            (log(20)/5)*pop$black + 
                                            (log(0.6))*pop$asian + 
                                            (log(0.8))*pop$hispanic +
                                            (log(1.2))*pop$otherrace + 
                                            (log(2))*pop$age_18_24 +
                                            (log(3))*pop$age_25_34 +
                                            (log(2.5))*pop$age_35_44 +
                                            (log(2))*pop$age_45_54 +
                                            (log(1.5))*pop$age_55_64 +
                                            (log(25))*pop$age_over64 +
                                            (log(20)/5)*pop$male +
                                            (log(20)/5)*pop$educ_ged +
                                            (log(0.9))*pop$educ_hs + 
                                            (log(0.8))*pop$educ_somecollege +
                                            (log(0.7))*pop$educ_associates +
                                            (log(0.6))*pop$educ_bachelors +
                                            (log(0.5))*pop$educ_advdegree))

# Check Exposure Frequency
summary(pop$A)

# Generate Outcome as Function of Exposure and Covariates - Socio-Demographic Patterns Loosely Based on Incidence of Lung Cancer in U.S. (https://www.cdc.gov/cancer/lung/statistics/race.htm)
trueOR <- log(2)
baseline_Y <- log(0.005)
set.seed(3)
pop$Y <- rbinom(N, size = 1, prob = expit(baseline_Y + 
                                            trueOR*pop$A +
                                            (log(20)/5)*pop$black + 
                                            (log(0.8)/1000)*pop$asian + 
                                            (log(0.6)/1000)*pop$hispanic +
                                            (log(1)/1000)*pop$otherrace + 
                                            (log(1.2)/1000)*pop$age_18_24 +
                                            (log(1.5)/1000)*pop$age_25_34 +
                                            (log(2)/1000)*pop$age_35_44 +
                                            (log(2.5)/1000)*pop$age_45_54 +
                                            (log(3)/1000)*pop$age_55_64 +
                                            (log(25)/5)*pop$age_over64 +
                                            (log(20)/5)*pop$male +
                                            (log(20)/5)*pop$educ_ged +
                                            (log(.8)/1000)*pop$educ_hs + 
                                            (log(.7)/1000)*pop$educ_somecollege +
                                            (log(.6)/1000)*pop$educ_associates +
                                            (log(.5)/1000)*pop$educ_bachelors +
                                            (log(.4)/1000)*pop$educ_advdegree +
                                            (log(1.2)/2000)*pop$X))

# Set all non-cases to time = end of study 
pop$X[pop$Y==0] <- 365.25*10

# Check Outcome Frequency
summary(pop$Y)

# Logistic Model to Check Fidelity of Exposure/Outcome Relationship in Total Population
mod <- glm(Y ~ A + black + asian + hispanic + otherrace + age_18_24 + age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_over64 + male + educ_ged + educ_hs + educ_somecollege + educ_associates + educ_bachelors + educ_advdegree, data=pop, family='binomial')
summary(mod)
exp(coef(mod))

# Check for presence of confounding in total population
crude_mod <- glm(Y ~ A, data = pop, family = "binomial")
summary(crude_mod)
exp(coef(crude_mod))

# True Odds Ratio in total population
trueOR_mod <- glm(Y ~ A + black + asian + hispanic + otherrace + age_18_24 + age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_over64 + male + educ_ged + educ_hs + educ_somecollege + educ_associates + educ_bachelors + educ_advdegree, data=pop, family='binomial')
trueOR <- as.numeric(exp(coef(trueOR_mod)["A"])); trueOR

# True Cumulative Incidence Ratio in total population
# log binomial model
trueCIR_lbmod <- glm(Y ~ A + black + asian + hispanic + otherrace + age_18_24 + age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_over64 + male + educ_ged + educ_hs + educ_somecollege + educ_associates + educ_bachelors + educ_advdegree, data=pop, family=binomial(log))
trueCIR_lb <- as.numeric(exp(coef(trueCIR_lbmod)["A"])); trueCIR_lb
# poisson model
trueCIR_pmod <- glm(Y ~ A + black + asian + hispanic + otherrace + age_18_24 + age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_over64 + male + educ_ged + educ_hs + educ_somecollege + educ_associates + educ_bachelors + educ_advdegree, data=pop, family='poisson')
trueCIR_p <- as.numeric(exp(coef(trueCIR_pmod)["A"])); trueCIR_p

# True Incidence Density Raio in total population
trueIDR_mod <- glm(Y ~ A + black + asian + hispanic + otherrace + age_18_24 + age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_over64 + male + educ_ged + educ_hs + educ_somecollege + educ_associates + educ_bachelors + educ_advdegree, offset = log(X), data=pop, family='poisson')
trueIDR <- as.numeric(exp(coef(trueIDR_mod)["A"])); trueIDR

## To save population data
#write.csv(pop, "data/population.data.csv", row.names=F)
  
# END
