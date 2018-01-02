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
#          12/28/2017: EM cleaned up population generation; added function to apply the three basic
#                      case-control designs with simple random sampling; wrapped the function in 
#                      a test simulation loop. 
#          01/02/2018: CL added bias, variance and MSE to test distribution of estimates
################################################################################################

# PENDING QUESTIONS TO CHECK WITH JEN AND PATRICK
  # Which type(s) of multistage surveys to mimic? Start with ACS?
  # Current covariates we have are XX. Covariates include variables we will want to sample based on (e.g. cluster), and ones we will want to control for in the analysis (e.g. education). Exposure and outcome mechanisms have been honed based on these, but are there others that would be critical to include? Some others we could consider: household income, food stamps, household composition, marital status, foreign born, citizenship status, years in US, language spoken, health insurance, school attendance, employment status, occupation, household or personal income, poverty status, residential transience, disabilities/difficulties of living, veteran status.
  # Specifications: We created a time of case occurrence drawn from a uniform distribution; are assuming 10 years of perfect follow up for controls; population of 1 million, outcome frequency of 6%, true OR of 2, etc. Sound good?
  # Patrick - For case-cohort, cox model is the right one to use? How to deal with tied failure times in density sampled?

# NEXT STEPS 
  # Continue to develop and test sampling and analysis part of code for the three different types of case control designs (cumulative, case-cohort, density-sampled), with controls drawn randomly from the population. Confirm that each of these designs recovers the parameter/quantity we expect.

# QUESTIONS FOR US TO DISCUSS AT OUR NEXT MEETING
  # In the case-cohort design, the controls are not sampled separately from the cases. Rather, you sample from a baseline cohort irrespective of later case/control status. Think more about whether the case-cohort design is relevant to what we're doing and if so, how it would mesh with the concept of drawing controls from complex surveys. 
  # Do we want to only include individuals 18 and older? How do we deal with education covariate for youngest age group?
  # Do we ultimately want to make it very easy for users to alter true odds ratio and outcome prevalence, or should we just generate an assortment of datasets that have a range of combinations of true odds ratio and prevalences. That is, will the code to create the population be published or will the final population be provided as a dataset.

# FEATURES WE WILL/MAY WANT TO INCORPORATE LATER
  # Different types of sampling of controls: cluster sampling, stratified, ACS, NHANES, Ad Health, NHIS, BRFSS
  # Matching 

# ASSUMPTIONS
  # Perfect follow-up
  # No competing risks

# NOTES
  # Analyzing ACS data: https://stats.idre.ucla.edu/other/mult-pkg/faq/sample-setups-for-commonly-used-survey-data-sets/



####### PROGRAM START

# Clear workspace
rm(list=ls())

# Load Packages. Be sure your local computer has these installed before running. No need to re-install each time you run. 
library(readstata13)
library(locfit) # for expit function
library(Epi) # for case-cohort and density sampling designs
library(survival) # for clogit analysis

# Set Working Directory
# setwd("~/Documents/PhD/Fall 2017/GSR/Case Control Simulation/Raw Data") # Chris's directory
#setwd("C:/Users/kecolson/Google Drive/simulation/case-control") # Ellie's directory
setwd("C:/Users/Catherine/Desktop/Case Control GSR") # Catherine's directory

########################
# Create and Save Population Data
########################

if (1==2) { # Use to skip over this section after we've run it for the first time. 1 never = 2, so this section will be skipped. 

# Import Raw ACS Data - California 2010-2013
raw <- read.dta13("data/usa_00005_10_13.dta")

# Subset Data - Relevant Variables and Single Year (2010)
data <- raw[raw$year==2010, 
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

  # Age (Seven Categories: Under 18; 18-24; 25-34; 35-44; 45-54; 55-64; 65 and older)
  pop$agenum <- as.numeric(pop$age)
  pop$age_u18   <- ifelse(pop$agenum<18,1,0)
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
              'age_u18','age_18_24','age_25_34','age_35_44','age_45_54','age_55_64','age_over64',
              'educ_ged','educ_hs','educ_somecollege','educ_associates','educ_bachelors','educ_advdegree')]

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
set.seed(4)
pop$time[pop$Y==1] <- runif(sum(pop$Y), 0, 365.25*10)
pop$time[pop$Y==0] <- 365.25*10

### KEEP POPULATION FIXED - TRUTH

## Save population data
write.csv(pop, "data/population.data.csv", row.names=F)

rm(list=setdiff(ls(), "pop")) # clear everything except population data

}


########################
# Test Designs
########################

# Bring in data
pop <- read.csv("data/population.data.csv", stringsAsFactors = F)

# Sample down the dataset for testing purposes
#data <- pop[sample(1:nrow(pop), 4000, replace=F),]

####
# Create a function that will apply the different study designs and analyses
####

study <- function(seed, # random seed to make sampling replicable
                  cctype = "cumulative", # options will be "cumulative", "cch" for case-cohort, and "density" for density-sampled
                  samp   = "srs", # sampling of controls. Options will be "srs" for simple random sample, and others TBD
                  ratio = 1, # ratio of controsl to cases
                  data, # argument to provide the population data. Population data will take the format of the pop data we've created. 
                  subcohort.size # size of subcohort for case-cohort design, if relevant
                   # any other options here TBD
) {
  
  ####### PHASE 1: source of cases and controls -- SRS, complex survey, etc. 
  set.seed(seed)
  
  allcases    <- data[data$Y==1,]
  allcontrols <- data[data$Y==0,]
  
  if (samp=="srs") { # simple random sample of controls
    
    control.samp <- allcontrols[sample(1:nrow(allcontrols), size = nrow(allcases)*ratio, replace=F),] 
    
  } else if (samp=="clustered") {
    
  } else if (samp=="ACS") {
    
  } else if (samp=="NHANES") {
    
  }
  
  ####### PHASE 2: implement case-control - cumulative, density sampled, or case-cohort and analyse the data appropriately
  set.seed(seed+1)
  
  # CUMULATIVE CASE-CONTROL
  if (cctype=="cumulative") {
    
    # Apply design
    sample <- rbind(allcases, control.samp) 
    
    # Run model
    mod <- glm(Y ~ A + black + asian + hispanic + otherrace + age_18_24 + age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_over64 + male + educ_ged + educ_hs + educ_somecollege + educ_associates + educ_bachelors + educ_advdegree, 
               data=sample, family='binomial')
    
    # Pull the main point estimate and CI
    est <- exp(coef(mod)[2])
    lower <- exp(coef(mod)[2] - 1.96*summary(mod)$coefficients[2,2])
    upper <- exp(coef(mod)[2] + 1.96*summary(mod)$coefficients[2,2])

  # CASE-COHORT
  } else if (cctype=="cch") {
    
    # Apply design
    # Give people an ID
    data$id <- 1:nrow(data)
      
    # Identify a subcohort from baseline (randomly drawn from base population). 
    subcohort <- data[sample(1:nrow(data), subcohort.size, replace=F),]
    subcohort$subcohort <- 1
    
    # Take as the analysis data, this subcohort, plus and additional cases not in the subcohort
    extracases <- data[!data$id %in% subcohort$id & data$Y==1,]
    extracases$subcohort <- 0
    sample <- rbind(subcohort,extracases)
    
    # Run model - cox regression for case-cohort design. I'm not 100% sure this is correct
    mod <- cch(Surv(time, Y) ~ A + black + asian + hispanic + otherrace + age_18_24 + age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_over64 + male + educ_ged + educ_hs + educ_somecollege + educ_associates + educ_bachelors + educ_advdegree, 
               data = sample,
               subcoh =~subcohort, id=~id, 
               cohort.size=nrow(data))
    
    # Pull the main point estimate and CI
    est <- exp(coef(mod)[1])
    lower <- exp(coef(mod)[1] - 1.96*summary(mod)$coefficients[1,2])
    upper <- exp(coef(mod)[1] + 1.96*summary(mod)$coefficients[1,2])
    
  # DENSITY-SAMPLED CASE-CONTROL
  } else if (cctype=="density") {
    
    # Apply design
    sample <- ccwc(entry=0, exit=time, fail=Y, origin=0, controls=ratio, 
                   #match=list(), # use this argument for variables we want to match on
      include=list(A,black,asian,hispanic,otherrace,age_18_24,age_25_34,age_35_44,
                   age_45_54,age_55_64,age_over64,male,educ_ged,educ_hs,educ_somecollege,
                   educ_associates,educ_bachelors,educ_advdegree), data=data, silent=FALSE)
    
    # Run model
    mod <- clogit(Fail ~ A + black + asian + hispanic + otherrace + age_18_24 + age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_over64 + male + educ_ged + educ_hs + educ_somecollege + educ_associates + educ_bachelors + educ_advdegree + strata(Set), data = sample)

    # Pull the main point estimate and CI
    est <- exp(coef(mod)[1])
    lower <- exp(coef(mod)[1] - 1.96*summary(mod)$coefficients[1,3])
    upper <- exp(coef(mod)[1] + 1.96*summary(mod)$coefficients[1,3])
  }
  
  # Return the sampled data, and model object
  return(list(sample=sample, mod=mod, est=est, lower=lower, upper=upper))
}


####
## Test function
####

# Warning: the second and third take a long time to run.

study1 <- study(seed=1, cctype = "cumulative", samp = "srs", ratio = 1, data=pop)
round(c(study1$est, study1$lower, study1$upper),3) 

study2 <- study(seed=1, cctype = "cch", samp = "srs", ratio = 1, data=pop, subcohort.size=10000) 
round(c(study2$est, study2$lower, study2$upper),3) 

study3 <- study(seed=1, cctype = "density", samp = "srs", ratio = 1, data=pop) # warning has tied failure times. Look into this/check with Patrick.
round(c(study3$est, study3$lower, study3$upper),3) 

# ECM 12/28/2017: these all appear to be producing estimates right around 2, the true OR. This is good.


####
## Test assessing distribution of estimates
####

# We will want to parallelize this for 1000-2000 sims and more complex designs. 

nsims <- 50 # I run out of memory when I try to do this on more than 70 sims.

results <- data.frame(sim=1:nsims, est = NA, lower = NA, upper=NA)
samples <- models <- NULL

for (i in 1:nsims) {
  print(paste0("iteration ",i))
  
  run <- study(seed=i, cctype = "cumulative", samp   = "srs", ratio = 1, data = pop)
  
  samples[[i]] <- run$sample
  models [[i]] <- run$mod
  
  results$est  [results$sim==i] <- run$est
  results$lower[results$sim==i] <- run$lower
  results$upper[results$sim==i] <- run$upper
  
  rm(run)
}


# Summarize distribution of point estimates and CIs
summary(results)

hist(results$est)

# Calculate 95% CI coverage - % of calculated CIs that include the true OR
results$CIcover <- as.numeric(results$lower<=2 & results$upper>=2)
round(mean(results$CIcover, na.rm=T)*100,1)

# Calculate Bias - average distance of point estimates away from true OR in repeated simulations
bias <- mean(results$est - 2) #should assign trueOR, is currently in line 164 which gets skipped
bias

# Calculate Variance - variance of point estimate of repeated simulations
variance <- var(results$est)
variance

# Calculate MSE - mean squared error of point estimate of repeated simulations
MSE <- mean((results$est - 2)^2)
MSE


# END