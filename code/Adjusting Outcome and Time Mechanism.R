################################################################################################
# NAME: Case Control Simulation: Adjusting Outcome and Time Mechanism
# AUTHORS: Ellie Matthay, Catherine Li, Chris Rowe*
# DATE STARTED: 01/22/2018    
# PURPOSE: To play around with exposure and outcome mechanism, and confounding in "true" population
# UPDATES: 01/22/2018: CL adapted code from main Case Control Simulation Program to play around 
#                      with confounding and outcome
#          02/05/2018: CL added function to generate different outcome frequencies
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
library(MASS)
library(AER) # for dispersion test https://www.rdocumentation.org/packages/AER/versions/1.2-5/topics/dispersiontest

# Set Working Directory
#setwd("~/Documents/PhD/Ahern GSR/Case Control Simulation") # Chris's directory
#setwd("C:/Users/kecolson/Google Drive/simulation/case-control") # Ellie's directory
setwd("C:/Users/Catherine/Desktop/Case Control GSR") # Catherine's directory

# Read in population data
pop <- read.csv("data/population.data.csv", stringsAsFactors = F)

# Turn dummy variables into factors
pop$race <- factor(apply(pop[,8:12], 1, function(x) which(x == 1)), labels = colnames(pop[,8:12]))
pop$age <- factor(apply(pop[,14:19], 1, function(x) which(x == 1)), labels = colnames(pop[,14:19]))
pop$education <- factor(apply(pop[,20:26], 1, function(x) which(x == 1)), labels = colnames(pop[,20:26]))

N <- 1000000

# Generate Time Variable (X)
set.seed(1)

pop$time <- runif(N, 0, 365.25*10)
#pop$time <- rnorm(N, 365.25*5, 365.25*1)# Think about other distributions as well

# Generate Exposure
baseline_A <- log(0.4)
set.seed(2)
pop$A <- rbinom(N, size = 1, prob = expit(baseline_A + 
                                            (log(1))*pop$black + 
                                            (log(1.1))*pop$asian + 
                                            (log(0.8))*pop$hispanic +
                                            (log(0.9))*pop$otherrace + 
                                            (log(1.1))*pop$age_25_34 +
                                            (log(1.2))*pop$age_35_44 +
                                            (log(1.3))*pop$age_45_54 +
                                            (log(1.4))*pop$age_55_64 +
                                            (log(3))*pop$age_over64 +
                                            (log(3))*pop$male +
                                            (log(3))*pop$educ_ged +
                                            (log(1.2))*pop$educ_hs + 
                                            (log(1.1))*pop$educ_somecollege +
                                            (log(1))*pop$educ_associates +
                                            (log(0.9))*pop$educ_bachelors +
                                            (log(0.8))*pop$educ_advdegree))

# Check Exposure Frequency
summary(pop$A)

# Generate Outcome Function
generate_outcome <- function(OR, baselineY) {
  OR <- log(OR)
  baseline_Y <- log(baselineY) # set to 0.06 for 20% (2.0 vs 2.8), 0.025 for 10% (2.0 vs. 2.9), 0.005 for 2% (2.0 vs. 3.09)
  set.seed(3)
  pop$Y <- rbinom(N, size = 1, prob = expit(baseline_Y + 
                                              OR*pop$A +
                                              (log(1.1))*pop$black + 
                                              (log(0.9))*pop$asian + 
                                              (log(1.1))*pop$hispanic +
                                              (log(1))*pop$otherrace + 
                                              (log(1.1))*pop$age_25_34 +
                                              (log(1.2))*pop$age_35_44 +
                                              (log(1.3))*pop$age_45_54 +
                                              (log(1.4))*pop$age_55_64 +
                                              (log(3))*pop$age_over64 +
                                              (log(3))*pop$male +
                                              (log(3))*pop$educ_ged +
                                              (log(1.1))*pop$educ_hs + 
                                              (log(1))*pop$educ_somecollege +
                                              (log(.9))*pop$educ_associates +
                                              (log(.8))*pop$educ_bachelors +
                                              (log(.7))*pop$educ_advdegree)) 
                                              #+ (log(1.00000000000001))*pop$time)) No longer including time in outcome mechanims
  
}

pop$Y.20 <- generate_outcome(OR = 2, baselineY = 0.06) # for ~20% outcome frequency
pop$Y.10 <- generate_outcome(OR = 2, baselineY = 0.025) # for ~10% outcome frequency
pop$Y.2 <- generate_outcome(OR = 2, baselineY = 0.005) # for ~2% outcome frequency
pop$Y.05 <- generate_outcome(OR = 2, baselineY = 0.0011) # for ~0.05% outcome frequency

# Set all non-cases to time = end of study 
pop$time.20 <- ifelse(pop$Y.20==0, 365.25*10, pop$time)
pop$time.10 <- ifelse(pop$Y.10==0, 365.25*10, pop$time)
pop$time.2 <- ifelse(pop$Y.2==0, 365.25*10, pop$time)
pop$time.05 <- ifelse(pop$Y.05==0, 365.25*10, pop$time)

# Check Outcome Frequency
summary(pop$Y.20)
summary(pop$Y.10)
summary(pop$Y.2)
summary(pop$Y.05)

# Set X and Y to run following models
pop$Y <- pop$Y.2
pop$time <- pop$time.2

# Logistic Model to Check Fidelity of Exposure/Outcome Relationship in Total Population
mod <- glm(Y ~ A + black + asian + hispanic + otherrace + age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_over64 + male + educ_ged + educ_hs + educ_somecollege + educ_associates + educ_bachelors + educ_advdegree, data=pop, family='binomial')
summary(mod)
exp(coef(mod))

# Check for presence of confounding in total population
crude_mod <- glm(Y ~ A, data = pop, family = "binomial")
summary(crude_mod)
exp(coef(crude_mod))

# True Odds Ratio in total population
trueOR_mod <- glm(Y ~ A + black + asian + hispanic + otherrace + age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_over64 + male + educ_ged + educ_hs + educ_somecollege + educ_associates + educ_bachelors + educ_advdegree, data=pop, family='binomial')
trueOR <- as.numeric(exp(coef(trueOR_mod)["A"])); trueOR

# True Cumulative Incidence Ratio in total population
# log binomial model --> more likely to converge when outcome frequency is low?
trueCIR_lbmod <- glm(Y ~ A + black + asian + hispanic + otherrace + age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_over64 + male + educ_ged + educ_hs + educ_somecollege + educ_associates + educ_bachelors + educ_advdegree, data=pop, family=binomial(log))
trueCIR_lb <- as.numeric(exp(coef(trueCIR_lbmod)["A"])); trueCIR_lb
# poisson model
trueCIR_pmod <- glm(Y ~ A + black + asian + hispanic + otherrace + age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_over64 + male + educ_ged + educ_hs + educ_somecollege + educ_associates + educ_bachelors + educ_advdegree, data=pop, family='poisson')
trueCIR_p <- as.numeric(exp(coef(trueCIR_pmod)["A"])); trueCIR_p

# True Incidence Density Ratio in total population
# Poisson Model with X offset
trueIDR_mod <- glm(Y ~ A + black + asian + hispanic + otherrace + age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_over64 + male + educ_ged + educ_hs + educ_somecollege + educ_associates + educ_bachelors + educ_advdegree, offset = log(time), data=pop, family='poisson')
trueIDR <- as.numeric(exp(coef(trueIDR_mod)["A"])); trueIDR
# Quasi-poisson model with X offset
trueIDR_mod.qp <- glm(Y ~ A + black + asian + hispanic + otherrace + age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_over64 + male + educ_ged + educ_hs + educ_somecollege + educ_associates + educ_bachelors + educ_advdegree, offset = log(time), data=pop, family='quasipoisson')
trueIDR.qp <- as.numeric(exp(coef(trueIDR_mod.qp)["A"])); trueIDR.qp
# Negative Binomial Model
trueIDR_mod.nb <- glm.nb(Y ~ A + black + asian + hispanic + otherrace + age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_over64 + male + educ_ged + educ_hs + educ_somecollege + educ_associates + educ_bachelors + educ_advdegree + offset(log(time)), data=pop, control=glm.control(maxit=50))
trueIDR.nb <- as.numeric(exp(coef(trueIDR_mod.nb)["A"]))

# Looking at dispersion
dispersiontest(trueCIR_pmod)
dispersiontest(trueIDR_mod)
pchisq(2 * (logLik(trueIDR_mod) - logLik(trueIDR_mod.nb)), df = 1, lower.tail = FALSE)

 ## To save population data
#write.csv(pop, "data/population.data.csv", row.names=F)

# Frequency Tables
table <- xtabs(~ male + race, data = pop) 
#table <- xtabs(~ male + education, data = pop)
#table <- xtabs(~ male + age, data = pop)
#table <- xtabs(~ male + A, data = pop)
#table <- xtabs(~ race + education, data = pop)
#table <- xtabs(~ race + age, data = pop)
#table <- xtabs(~ race + A, data = pop)
#table <- xtabs(~ education + age, data = pop)
#table <- xtabs(~ education + A, data = pop)
#table <- xtabs(~ age + A, data = pop)

ftable(table)
prop.table(table)
summary(table)
# END