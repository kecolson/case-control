################################################################################################
# NAME: Case Control Simulation Program
# AUTHORS: Ellie Matthay, Catherine Li, Chris Rowe
# DATE STARTED: 12/15/2017
# PURPOSE: Produce an artifical population based on the 2010 California 
#          American Community Survey data and conduct various case control simulations 
#          to mimic selection of controls from external data sources with complex sampling methods.
# UPDATES: 12/15/2017: EM created basic code frame during meeting with CR and CL.
#          12/19/2017: CR cleaned up & drafted code for final population generation (incl.
#                      covariates, exposure, and outcome).
#          12/28/2017: EM made further updates 
################################################################################################

# COMMENTS AND PENDING QUESTIONS
  # Which type(s) of multistage surveys to mimic?
  # Any other covariates we want to pull?
  # Create exposure so it mimics a real relationship in real data, using the covariates we've used here? Same for outcome. 
  # Make sure true RR is on the right scale
  # Set seeds throughout
  # Think through true RR instead of OR
  # Add in: time of case occurence, control occurence  -- time of event / censoring and event indicator -- so that we can do density sampling
  # Want to build this out so it can do all three types of case control studies (+/- matching) 
  # And also multiple types of sampling - simple random, cluster sampling, stratified, etc. 
  # Assumptions: perfect follow up, no competing risks
  # Finalize covariates we are interested in (i.e., ones that we sample based on such as cluster and race and ones that we will control for but not sample on such as education)
  # Do we want to only include individuals 18 and older? How do we deal with educaiton covariate for youngest age group?
  # Do we ultimately want ot make it very easy for users to alter true odds ratio and outcome prevalence, or should we just generate an assortment of datasets that have a range of combinations of true odds ratio and prevalences. That is, will the code to create the population be published or will the final population be provided as a dataset.

# NEXT STEPS FROM 12/15/17 MEETING
  # Chris to try to hash out the covariates and exposure and outcome mechanisms and then we will check in
  # Catherine working on lit review
  # Ellie to work on setting us up fo

# PROGRAM START

# Install and Load Packages
install.packages("readstata13")
library(readstata13)
install.packages("locfit")
library(locfit)

# Set Working Directory
setwd("~/Documents/PhD/Fall 2017/GSR/Case Control Simulation/Raw Data")

# Import Raw ACD Data - California 2010-2013
raw <- read.dta13("usa_00005_10_13.dta")

# Subset Data - Relevant Variables and Single Year (2010)
data <- raw[raw$year==2010, 
             c('cluster','county','metro','city','citypop','puma','strata','sex',
               'age','educd','race','hispan','perwt')]

# Generate Synthetic Household ID (CR: Waiting until we confirm which geographic areas we will utilize for survey sampling simulation)
#data$hh <- 

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
  
  pop$black <- ifelse(pop$race=='Black/African-American, Non-Hispanic',1,0)
  pop$asian <- ifelse(pop$race=='Asian/Pacific Islander, Non-Hispanic',1,0)
  pop$hispanic <- ifelse(pop$race=='Hispanic',1,0)
  pop$otherrace <- ifelse(pop$race=='Other or Two or More Races, Non-Hispanic',1,0)

  # Sex (Two Categories: Male; Female)
  pop$male <- ifelse(pop$sex=='male',1,0)

  # Age (Seven Categories: Under 18; 18-24; 25-34; 35-44; 45-54; 55-64; 65 and older)
  pop$agenum <- as.numeric(pop$age)
  pop$age_18_24 <- ifelse(pop$agenum>=18 & pop$agenum<=24,1,0)
  pop$age_25_34 <- ifelse(pop$agenum>=25 & pop$agenum<=34,1,0)
  pop$age_35_44 <- ifelse(pop$agenum>=35 & pop$agenum<=44,1,0)
  pop$age_45_54 <- ifelse(pop$agenum>=45 & pop$agenum<=54,1,0)
  pop$age_55_64 <- ifelse(pop$agenum>=55 & pop$agenum<=64,1,0)
  pop$age_over64 <- ifelse(pop$agenum>=65,1,0)

  # Education (Seven Categores: Less than high school; GED; High School Diploma; Some College; Associate's Degree; Bachelor's Degree; Advanced Degree)
  pop$educ_ged <- ifelse(pop$educd=='ged or alternative credential',1,0)
  pop$educ_hs <- ifelse(pop$educd=='regular high school diploma',1,0)
  pop$educ_somecollege <- ifelse(pop$educd=='some college, but less than 1 year' | pop$educd=='1 or more years of college credit, no degree',1,0)
  pop$educ_associates <- ifelse(pop$educd=="associate's degree, type not specified",1,0)
  pop$educ_bachelors <- ifelse(pop$educd=="bachelor's degree",1,0)
  pop$educ_advdegree <- ifelse(pop$educd=="master's degree" | pop$educd=="professional degree beyond a bachelor's degree" | pop$educd=="doctoral degree",1,0)
 
# Keep only Final Variables
pop <- pop[,c('cluster','county','metro','city','citypop','puma','strata','black', 'asian', 'hispanic', 'otherrace', 'male', 'age_18_24', 'age_25_34', 'age_35_44', 'age_45_54', 'age_55_64', 'age_over64', 'educ_ged', 'educ_hs', 'educ_somecollege', 'educ_associates', 'educ_bachelors', 'educ_advdegree')]

# Generate Exposure - Socio-Demographic Patterns Loosely Based on Cigarette Smoking in U.S. (https://www.cdc.gov/tobacco/campaign/tips/resources/data/cigarette-smoking-in-united-states.html)
baseline_A <- log(0.5)
set.seed(2)
pop$A <- rbinom(N, size = 1, prob = expit(baseline_A + 
                                          (log(1))*pop$black + 
                                          (log(0.6))*pop$asian + 
                                          (log(0.8))*pop$hispanic +
                                          (log(1.2))*pop$otherrace + 
                                          (log(2))*pop$age_18_24 +
                                          (log(3))*pop$age_25_34 +
                                          (log(2.5))*pop$age_35_44 +
                                          (log(2))*pop$age_45_54 +
                                          (log(1.5))*pop$age_55_64 +
                                          (log(1.2))*pop$age_over64 +
                                          (log(1.5))*pop$male +
                                          (log(2))*pop$educ_ged +
                                          (log(0.9))*pop$educ_hs + 
                                          (log(0.8))*pop$educ_somecollege +
                                          (log(0.7))*pop$educ_associates +
                                          (log(0.6))*pop$educ_bachelors +
                                          (log(0.5))*pop$educ_advdegree))

  # Check Exposure Frequency
  summary(pop$A)

# Generate Outcome as Function of Exposure and Covariates - Socio-Demographic Patterns Loosely Based on Incidence of Lung Cancer in U.S. (https://www.cdc.gov/cancer/lung/statistics/race.htm)
trueOR <- log(2)
baseline_Y <- log(0.03)
set.seed(3)
pop$Y <- rbinom(N, size = 1, prob = expit(baseline_Y + 
                                          trueOR*pop$A +
                                          (log(1.2))*pop$black + 
                                          (log(0.8))*pop$asian + 
                                          (log(0.6))*pop$hispanic +
                                          (log(1))*pop$otherrace + 
                                          (log(1.2))*pop$age_18_24 +
                                          (log(1.5))*pop$age_25_34 +
                                          (log(2))*pop$age_35_44 +
                                          (log(2.5))*pop$age_45_54 +
                                          (log(3))*pop$age_55_64 +
                                          (log(3.5))*pop$age_over64 +
                                          (log(1.5))*pop$male +
                                          (log(.9))*pop$educ_ged +
                                          (log(.8))*pop$educ_hs + 
                                          (log(.7))*pop$educ_somecollege +
                                          (log(.6))*pop$educ_associates +
                                          (log(.5))*pop$educ_bachelors +
                                          (log(.4))*pop$educ_advdegree))

  # Check Outcome Frequency
  summary(pop$Y)

# Logistic Model to Check Fidelity of Exposure/Outcome Relationship in Total Population
mod <- glm(Y ~ A + black + asian + hispanic + otherrace + age_18_24 + age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_over64 + male + educ_ged + educ_hs + educ_somecollege + educ_associates + educ_bachelors + educ_advdegree, data=pop, family='binomial')
summary(mod)
exp(coef(mod))

# Generate Case Occurrence/Censoring Time Variable (10 Years of Follow-Up, Measured in Days)
pop$time[pop$Y==1] <- runif(N, 0, 365.25*10)
pop$time[pop$Y==0] <- 365.25*10

### KEEP POPULATION FIXED - TRUTH
### CR Stopped Here 12/17/2017 ###















###### SAMPLING SCHEME 1
# Take all cases
# SRS controls

ratioControlsCases <- 4

cases <- data2[data2$Y==1,]
controls <- data2[data2$Y==0,]

## Assess distribution of estimates
nsims <- 5

results <- data.frame(sim=1:nsims, est = NA, SE = NA)

for (i in 1:nsims) {
  print(paste0("iteration ",i))
  set.seed(i)
  
  # Phase 1: source of controls -- SRS, complex, etc. 
  control.samp <- controls[sample(1:nrow(controls), size = nrow(cases)*ratioControlsCases, replace=F),] # simple random
  
  # Phase 2: implement case-control - cumulative, density sampled, etc. 
  sample <- rbind(cases, control.samp) # cumulative - taking all controls
  
  sample <- ccwc() # density sample
  
  
  # Run model
  mod <- glm(Y ~ age + sex + A, data=sample, family=binomial)

  # Save point estimates and SE/CIs
  results$est[results$sim==i] <- summary(mod)$coefficients[4,1]
  results$SE [results$sim==i] <- summary(mod)$coefficients[4,2]
}


# Summarize distribution of point estimates and CIs
summary(results)
exp(results$est)




# END