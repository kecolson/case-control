################################################################################################
# NAME: pop_gen_strongconfounder.r
# AUTHORS: Ellie Matthay, Catherine Li, Chris Rowe
# DATE STARTED: 2/20/2018    
# PURPOSE: Produce an artifical population based on the 2010 California 
#          American Community Survey data for use in various case control simulations. 
#          Includes changes to original pop_gen.r to incorporate stronger confounding by gender
# UPDATES: 03/13/2018: CL added functions to add different exposure and outcome frequencies
#                      CL added function to calculate true measure of association and confounding
#                      for all exposure and outcome frequency combinations
################################################################################################

# ASSUMPTIONS
  # Perfect follow-up
  # No competing risks

# NOTES
  # Analyzing ACS data: https://stats.idre.ucla.edu/other/mult-pkg/faq/sample-setups-for-commonly-used-survey-data-sets/

# Clear workspace
rm(list=ls())

# Load Packages.
library(readstata13)
library(locfit) # for expit function

# Set Working Directory
setwd("~/Documents/PhD/github/case-control/") # Chris's directory
#setwd("C:/Users/kecolson/Google Drive/simulation/case-control-other") # Ellie's directory
#setwd("/Users/cxli/Documents/case-control/") # Catherine's directory

# Import Raw ACS Data - California 2010-2013
raw <- read.dta13("data/usa_00005_10_13.dta")

# Subset Data - Relevant Variables and Single Year (2010)
data <- raw[raw$year==2010 & as.numeric(raw$age) >= 18, # Year 2010 and age 18 or older only
             c('cluster','strata','serial','perwt',   # Survey design variables: HH cluster, HH strata, HH serial #, person weight
               'county','city','puma','cpuma0010',    # Geographic identifiers
               'sex','age','educd','race','hispan')]  # Individual-level covariates

# Expand Raw Dataset Using Individual Person Weights
fulldata <- data[rep(row.names(data),data$perwt),]

# Set Desired Size for Total Population (N) and Sample Total Population from Full Expanded Dataset
N <- 1000000
set.seed(1)
pop <- fulldata[sample(1:nrow(fulldata), size=N, replace=F),]

# Clean Covariates

  # Race/Ethnicity (Five Categories: White; Black/African-American, Non-Hispanic; Asian/Pacific Islander, Non-Hispanic; Other Race or Two or More Races, Non-Hispanic; Hispanic)
  # (NOTE: American Indian/Alaskan Native is included in "Other" Category)
  levels(pop$race) <- c('White, Non-Hispanic', 
                        'Black/African-American, Non-Hispanic', 
                        'Other or Two or More Races, Non-Hispanic', 
                        'Asian/Pacific Islander, Non-Hispanic', 
                        'Asian/Pacific Islander, Non-Hispanic', 
                        'Asian/Pacific Islander, Non-Hispanic', 
                        'Other or Two or More Races, Non-Hispanic', 
                        'Other or Two or More Races, Non-Hispanic', 
                        'Other or Two or More Races, Non-Hispanic', 
                        'Hispanic')
  
  pop$race[pop$hispan=='mexican' | pop$hispan=='puerto rican' | pop$hispan=='cuban' | pop$hispan=='other'] <- 'Hispanic'
  
  pop$white <- ifelse(pop$race=='White, Non-Hispanic',1,0)
  pop$black <- ifelse(pop$race=='Black/African-American, Non-Hispanic',1,0)
  pop$asian <- ifelse(pop$race=='Asian/Pacific Islander, Non-Hispanic',1,0)
  pop$hispanic <- ifelse(pop$race=='Hispanic',1,0)
  pop$otherrace <- ifelse(pop$race=='Other or Two or More Races, Non-Hispanic',1,0)

  # Sex (Two Categories: Male; Female)
  pop$male <- ifelse(pop$sex=='male',1,0)

  # Age (Six Categories: 18-24; 25-34; 35-44; 45-54; 55-64; 65 and older)
  pop$agenum <- as.numeric(pop$age)
  #pop$age_u18   <- ifelse(pop$agenum<18,1,0)
  pop$age_18_24 <- ifelse(pop$agenum>=18 & pop$agenum<=24,1,0)
  pop$age_25_34 <- ifelse(pop$agenum>=25 & pop$agenum<=34,1,0)
  pop$age_35_44 <- ifelse(pop$agenum>=35 & pop$agenum<=44,1,0)
  pop$age_45_54 <- ifelse(pop$agenum>=45 & pop$agenum<=54,1,0)
  pop$age_55_64 <- ifelse(pop$agenum>=55 & pop$agenum<=64,1,0)
  pop$age_over64 <- ifelse(pop$agenum>=65,1,0)

  # Education (Seven Categores: Less than high school; GED; High School Diploma; Some College; Associate's Degree; Bachelor's Degree; Advanced Degree)
  pop$educ_lesshs <- ifelse(pop$educd %in% c('n/a or no schooling','n/a','no schooling completed','nursery school to grade 4',
                                             'nursery school, preschool','kindergarten','grade 1, 2, 3, or 4','grade 1','grade 2','grade 3',
                                             'grade 4','grade 5, 6, 7, or 8','grade 5 or 6','grade 5','grade 6','grade 7 or 8','grade 7','grade 8',
                                             'grade 9','grade 10','grade 11','grade 12','12th grade, no diploma'),1,0)
  pop$educ_ged <- ifelse(pop$educd=='ged or alternative credential',1,0)
  pop$educ_hs <- ifelse(pop$educd=='regular high school diploma',1,0)
  pop$educ_somecollege <- ifelse(pop$educd=='some college, but less than 1 year' | pop$educd=='1 or more years of college credit, no degree',1,0)
  pop$educ_associates <- ifelse(pop$educd=="associate's degree, type not specified",1,0)
  pop$educ_bachelors <- ifelse(pop$educd=="bachelor's degree",1,0)
  pop$educ_advdegree <- ifelse(pop$educd=="master's degree" | pop$educd=="professional degree beyond a bachelor's degree" | pop$educd=="doctoral degree",1,0)
  
  
# Keep only Final Variables
pop <- pop[,c('cluster','strata','serial', # Survey design variables: HH cluster, HH strata, HH serial #, person weight
              'county','city','puma','cpuma0010', # Geographic identifiers
              'white','black','asian','hispanic','otherrace','male', # Individual-level covariates
              'age_18_24','age_25_34','age_35_44','age_45_54','age_55_64','age_over64',
              'educ_lesshs','educ_ged','educ_hs','educ_somecollege','educ_associates','educ_bachelors','educ_advdegree')]

# Generate Exposure - Socio-Demographic Patterns Loosely Based on Cigarette Smoking in U.S. (https://www.cdc.gov/tobacco/campaign/tips/resources/data/cigarette-smoking-in-united-states.html)
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
                                          (log(8))*pop$male +
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
  
  
pop$A.50 <- exposure_gen(log(0.27)) # 50.3%
pop$A.20 <- exposure_gen(log(0.047)) # 20.02%
pop$A.10 <- exposure_gen(log(0.019)) # 10.09%
pop$A.05 <- exposure_gen(log(0.00845)) # 5.01%

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
                    (log(8))*pop$male +
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
                          BaselineHazard = c(log(0.0000064), log(0.0000084), log(0.0000093), log(0.0000098)),
                          expfreq = c("A.50", "A.20", "A.10", "A.05"), 
                          outcomefreq = "Y.20")

outcome.10 <- outcome_gen(trueRate = log(2), 
                          BaselineHazard = c(log(0.0000028), log(0.0000036), log(0.0000040), log(0.0000043)),
                          expfreq = c("A.50", "A.20", "A.10", "A.05"), 
                          outcomefreq = "Y.10")

outcome.02 <- outcome_gen(trueRate = log(2), 
                          BaselineHazard = c(log(0.00000050), log(0.00000063), log(0.00000071), log(0.00000077)),
                          expfreq = c("A.50", "A.20", "A.10", "A.05"), 
                          outcomefreq = "Y.02")

outcome.005 <- outcome_gen(trueRate = log(2), 
                           BaselineHazard = c(log(0.000000119), log(0.000000152), log(0.000000175), log(0.000000188)),
                           expfreq = c("A.50", "A.20", "A.10", "A.05"), 
                           outcomefreq = "Y.005")

# Add the outcome frequencies and time to the population data
pop <- data.frame(cbind(pop, outcome.20, outcome.10, outcome.02, outcome.005))


# make nested for loop to check for fidelity of exposure outcome relationship in total population
# and check for confounding

exposurefreq <- list("A.50", "A.20", "A.10", "A.05")
outcomefreq <- list("Y.20", "Y.10", "Y.02", "Y.005")

TrueMeasures_OR <- matrix(NA, ncol = length(exposurefreq), nrow = length(outcomefreq))
TrueMeasures_IDR <- matrix(NA, ncol = length(exposurefreq), nrow = length(outcomefreq))
rownames(TrueMeasures_OR) <- rownames(TrueMeasures_IDR) <- c("A.50", "A.20", "A.10", "A.05")
colnames(TrueMeasures_OR) <- colnames(TrueMeasures_IDR) <- c("Y.20", "Y.10", "Y.02", "Y.005")

Crude_OR <- matrix(NA, ncol = length(exposurefreq), nrow = length(outcomefreq))
Crude_IDR <- matrix(NA, ncol = length(exposurefreq), nrow = length(outcomefreq))
rownames(Crude_OR) <- rownames(Crude_IDR) <- c("A.50", "A.20", "A.10", "A.05")
colnames(Crude_OR) <- colnames(Crude_IDR) <- c("Y.20", "Y.10", "Y.02", "Y.005")

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
    
    #print(summary(log_mod))
    TrueMeasures_OR[i, j] <- exp(coef(log_mod)[2])
    print(pop[,paste("trueOR", j, i, sep = ".")] <- exp(coef(log_mod)[2]))
    
    # Poisson Model
    pos_mod <- glm(pop[,Y] ~ pop[,A] + black + asian + hispanic + otherrace + age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_over64 + male + 
                     educ_ged + educ_hs + educ_somecollege + educ_associates + educ_bachelors + educ_advdegree, offset = log(pop[,time]), data=pop, family='poisson')
    #print(summary(pos_mod))
    TrueMeasures_IDR[i, j] <- exp(coef(pos_mod)[2])
    print(pop[,paste("trueIDR", j, i, sep = ".")] <- exp(coef(pos_mod)[2]))
    
    # Check for presence of confounding in total population
    # Logistic Model
    logcrude_mod <- glm(pop[,Y] ~ pop[,A], data = pop, family = "binomial")
    #summary(logcrude_mod)
    Crude_OR[i,j] <- exp(coef(logcrude_mod)[2])
    
    # Poisson Model
    poscrude_mod <- glm(pop[,Y] ~ pop[,A], offset = log(pop[,time]), data=pop, family='poisson')
    #summary(poscrude_mod)
    Crude_IDR[i,j] <- exp(coef(poscrude_mod)[2])
  }
}


### KEEP POPULATION FIXED - TRUTH

## Save population data
write.csv(pop, "data/population.data.strongconfounder.csv", row.names=F)

# END
