################################################################################################
# NAME: Case Control Simulation: Adjusting Outcome and Time Mechanism
# AUTHORS: Ellie Matthay, Catherine Li, Chris Rowe*
# DATE STARTED: 01/22/2018    
# PURPOSE: To play around with exposure and outcome mechanism, and confounding in "true" population
# UPDATES: 01/22/2018: CL adapted code from main Case Control Simulation Program to play around 
#                      with confounding and outcome
#          02/05/2018: CL added function to generate different outcome frequencies
#          03/13/2018: CL updated code to reflect new mechanism for generating outcome
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

# Generate Exposure
exposure_gen <- function(BaselineA){
  set.seed(2)
  
  exposure <- rbinom(N, size = 1, prob = expit(BaselineA + 
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
  print(summary(exposure))
  return(exposure)
}


pop$A.50 <- exposure_gen(log(0.4)) # 48%
pop$A.20 <- exposure_gen(log(0.094)) # 20%
pop$A.10 <- exposure_gen(log(0.041)) # 10.2%
pop$A.05 <- exposure_gen(log(0.019)) # 5.2%

# Generate Outcome as Function of Exposure and Covariates - Socio-Demographic Patterns Loosely Based on Incidence of Lung Cancer in U.S. (https://www.cdc.gov/cancer/lung/statistics/race.htm)
# Using exponential distribution

outcome_gen <- function(trueRate = log(2),
                        BaselineHazard, 
                        expfreq, # vector of exposure frequency names (format "A.50")
                        outcomefreq # format "Y.02"
){
  a <- length(expfreq) # to get the number of exposure frequencies
  output <- NULL
  for(A in 1:a){
    set.seed(3)
    lambda <- exp(BaselineHazard[A] + 
                    trueRate*pop[,expfreq[A]] +
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
                    (log(0.9))*pop$educ_associates +
                    (log(0.8))*pop$educ_bachelors +
                    (log(0.7))*pop$educ_advdegree)
    set.seed(4)
    time <- rexp(N, lambda)
    
    # Set those with event before 10 years as a case
    outcome <- rep(NA, N)
    outcome[time <= 365.25*10] <- 1
    outcome[time > 365.25*10] <- 0
    
    # Set event time for those did not get the event to end of follow up
    time[outcome == 0] <- 365.25*10
    
    # Check Outcome Frequency
    print(summary(outcome))
    
    # Return Outcome and Time
    newoutput <- cbind(outcome, time)
    colnames(newoutput) <- c(paste(outcomefreq, expfreq[A], sep = "."), #column name for outcome 
                             paste("time", outcomefreq, expfreq[A], sep = ".")) # column name for time
    output <- cbind(output, newoutput)
  }
  
  return(output)
}
outcome.20 <- outcome_gen(trueRate = log(2), 
                          BaselineHazard = c(log(0.0000145), log(0.0000185), log(0.00002), log(0.000021)),
                          expfreq = c("A.50", "A.20", "A.10", "A.05"), 
                          outcomefreq = "Y.20")

outcome.10 <- outcome_gen(trueRate = log(2), 
                          BaselineHazard = c(log(0.0000065), log(0.0000081), log(0.000009), log(0.0000095)),
                          expfreq = c("A.50", "A.20", "A.10", "A.05"), 
                          outcomefreq = "Y.10")

outcome.02 <- outcome_gen(trueRate = log(2), 
                          BaselineHazard = c(log(0.0000012), log(0.0000015), log(0.00000166), log(0.00000176)),
                          expfreq = c("A.50", "A.20", "A.10", "A.05"), 
                          outcomefreq = "Y.02")

outcome.005 <- outcome_gen(trueRate = log(2), 
                          BaselineHazard = c(log(0.000000295), log(0.00000036), log(0.0000004), log(0.000000425)),
                          expfreq = c("A.50", "A.20", "A.10", "A.05"), 
                          outcomefreq = "Y.005")

# add outcomes and times to pop data
pop <- data.frame(cbind(pop, outcome.02))


# make nested for loop to check for fidelity of exposure outcome relationship in total population
# and check for confounding
exposurefreq <- list("A.50", "A.20", "A.10", "A.05")
outcomefreq <- list("Y.02")
#outcomefreq <- list("Y.20", "Y.10", "Y.02", "Y.005")
TrueMeasures <- matrix(NA, ncol = length(exposurefreq), nrow = length(outcomefreq))
for(i in exposurefreq) {
  for(j in outcomefreq){
    print(paste(j, i, sep = "."))
    print(paste("time", j, i, sep = "."))
    
    Y <- paste(j, i, sep = ".")
    A <- i
    time <- paste("time", j, i, sep = ".")
    
    # Check Fidelity of Exposure/Outcome Relationship in Total Population
    # Logistic Model
    log_mod <- glm(pop[,Y] ~ pop[,A] + black + asian + hispanic + otherrace + age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_over64 + male + 
                     educ_ged + educ_hs + educ_somecollege + educ_associates + educ_bachelors + educ_advdegree, data=pop, family='binomial')
    
    summary(log_mod)
    
    # Poisson Model
    pos_mod <- glm(pop[,Y] ~ pop[,A] + black + asian + hispanic + otherrace + age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_over64 + male + 
                     educ_ged + educ_hs + educ_somecollege + educ_associates + educ_bachelors + educ_advdegree, offset = log(pop[,time]), data=pop, family='poisson')
    summary(pos_mod)
    
    # Check for presence of confounding in total population
    # Logistic Model
    logcrude_mod <- glm(pop[,Y] ~ pop[,A], data = pop, family = "binomial")
    summary(logcrude_mod)
    exp(coef(logcrude_mod))
    # Poisson Model
    poscrude_mod <- glm(pop[,Y] ~ pop[,A], offset = log(pop[,time]), data=pop, family='poisson')
    summary(poscrude_mod)
    exp(coef(poscrude_mod))
  }
}
# Check Fidelity of Exposure/Outcome Relationship in Total Population
# Logistic Model
log_mod <- glm(Y.02.A.5 ~ A.5 + black + asian + hispanic + otherrace + age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_over64 + male + 
                 educ_ged + educ_hs + educ_somecollege + educ_associates + educ_bachelors + educ_advdegree, data=pop, family='binomial')

summary(log_mod)
print(pop$trueOR.Y.02.A.5 <- exp(coef(log_mod)[2]))
# Poisson Model
pos_mod <- glm(Y.02.A.5 ~ A.5 + black + asian + hispanic + otherrace + age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_over64 + male + 
                 educ_ged + educ_hs + educ_somecollege + educ_associates + educ_bachelors + educ_advdegree, offset = log(time), data=pop, family='poisson')
summary(pos_mod)
print(pop$trueIDR.Y.02.A.5 <- exp(coef(pos_mod)[2]))

# Check for presence of confounding in total population
# Logistic Model
logcrude_mod <- glm(Y.02.A.5 ~ A.5, data = pop, family = "binomial")
summary(logcrude_mod)
exp(coef(logcrude_mod)) # crude OR = 2.9728 vs. true OR = 1.9495
# Poisson Model
poscrude_mod <- glm(Y.02.A.5 ~ A.5, offset = log(time), data=pop, family='poisson')
summary(poscrude_mod)
exp(coef(poscrude_mod)) # crude IDR = 2.9422 vs. true IDR = 1.9342

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