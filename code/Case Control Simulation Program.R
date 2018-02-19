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
#          01/16/2018: CR added simple probability and clustered sampling designs
#          01/17/2018: CL added crude OR calcuation for total population, true OR, CIR, and IDR values
#                      for the total population
#          01/31/2018: CL added code to look at dispersion
#          02/05/2018: CL changed total population so that only those aged 18 or older subset from ACS,
#                      changed code throughout to reflect this change (18-24 now reference group)
#          02/06/2018: CR added second cluster design, calculated sampling weights for each design, 
#                      and incorporated weights into cumulative case control and incidence density 
#                      sampling analyses.
#          02/14/2018: CL changed outcome mechansim to come from exponential distribution and revised 
#                      code for true OR and IDR accordingly
#          02/19/2018: EM cleaned up commenting throughout; added lines to save true OR and IDR; 
#                      corrected srs weights for controls to be 1/pr(samp) instead of 1; removed
#                      case-cohort
################################################################################################

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
library(dplyr)
library(MASS) # for negative binomial models
library(AER) # for dispersion test https://www.rdocumentation.org/packages/AER/versions/1.2-5/topics/dispersiontest
library(pscl) 
library(geepack) # for modified poisson

# Set Working Directory
#setwd("~/Documents/PhD/Ahern GSR/Case Control Simulation") # Chris's directory
#setwd("C:/Users/kecolson/Google Drive/simulation/case-control-other") # Ellie's directory
setwd("C:/Users/Catherine/Desktop/Case Control GSR") # Catherine's directory

########################
# Create and Save Population Data
########################

if (1==2) { # Use to skip over this section after we've run it for the first time. 1 never = 2, so this section will be skipped. 

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

  # Age (Seven Categories: Under 18; 18-24; 25-34; 35-44; 45-54; 55-64; 65 and older)
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
  summary(pop$A) # 48%

# Generate Outcome as Function of Exposure and Covariates - Socio-Demographic Patterns Loosely Based on Incidence of Lung Cancer in U.S. (https://www.cdc.gov/cancer/lung/statistics/race.htm)
# Using exponential distribution
trueRate <- log(2)
#baseline_hazard <- log(0.00005) # for 50% outcome freq
baseline_hazard <- log(0.0000012) #for 2% outcome freq
set.seed(3)
lambda <- exp(baseline_hazard + 
              trueRate*pop$A +
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
pop$time <- rexp(N, lambda)

# Set those with event before 10 years as a case
pop$Y[pop$time <= 365.25*10] <- 1
pop$Y[pop$time > 365.25*10] <- 0

# Set event time for those did not get the event to end of follow up
pop$time[pop$Y==0] <- 365.25*10

  # Check Outcome Frequency
  summary(pop$Y) # 2%

# Check Fidelity of Exposure/Outcome Relationship in Total Population
# Logistic Model
log_mod <- glm(Y ~ A + black + asian + hispanic + otherrace + age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_over64 + male + educ_ged + educ_hs + educ_somecollege + educ_associates + educ_bachelors + educ_advdegree, data=pop, family='binomial')
summary(log_mod)
exp(coef(log_mod))
# Poisson Model
pos_mod <- glm(Y ~ A + black + asian + hispanic + otherrace + age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_over64 + male + educ_ged + educ_hs + educ_somecollege + educ_associates + educ_bachelors + educ_advdegree, offset = log(time), data=pop, family='poisson')
summary(pos_mod)
exp(coef(pos_mod))

# Check for presence of confounding in total population
# Logistic Model
logcrude_mod <- glm(Y ~ A, data = pop, family = "binomial")
summary(logcrude_mod)
exp(coef(logcrude_mod)) # crude OR = 2.9728 vs. true OR = 1.9495
# Poisson Model
poscrude_mod <- glm(Y ~ A, offset = log(time), data=pop, family='poisson')
summary(poscrude_mod)
exp(coef(poscrude_mod)) # crude IDR = 2.9422 vs. true IDR = 1.9342

### KEEP POPULATION FIXED - TRUTH

## Save population data
write.csv(pop, "data/population.data.csv", row.names=F)
write.csv(exp(coef(log_mod)["A"]), "data/trueOR.csv", row.names=F)
write.csv(exp(coef(pos_mod)["A"]), "data/trueIDR.csv", row.names=F)

rm(list=setdiff(ls(), "pop")) # clear everything except population data

}


########################
# Test Designs
########################

# Bring in data and true parameters
pop <- read.csv("data/population.data.csv", stringsAsFactors = F)

trueOR  <- read.csv("data/trueOR.csv",  stringsAsFactors = F)
trueIDR <- read.csv("data/trueIDR.csv", stringsAsFactors = F)

# Sample down the dataset for testing purposes
#data <- pop[sample(1:nrow(pop), 4000, replace=F),]

# Full dataset for testing purposes
data <- pop

####
# Create a function that will apply the different study designs and analyses
####

study <- function(seed, # random seed to make sampling replicable
                  cctype = "cumulative", # options will be "cumulative" and "density" for density-sampled
                  samp   = "srs", # sampling of controls. Options will be:
                    # "srs" for simple random sample
                    # "sps" for simple probability sample with know probability of selection for each individual
                    # "clustered1" for single stage clustered design
                    # "clustered2" for two-stage clustered design 
                  ratio = 1, # ratio of controls to cases
                  data # argument to provide the population data. Population data will take the format of the pop data we've created. 
                   # any other options here TBD
) {
  
  ####### PHASE 1: source of cases and controls -- SRS, complex survey, etc. 
  set.seed(seed)
  
  # Select cases and controls
  allcases    <- data[data$Y==1,]
  allcontrols <- data[data$Y==0,]

  # Create Case Weights
  allcases$sampweight <- 1
  
  if (samp=="srs") { # simple random sample of controls
    
    control.samp <- allcontrols[sample(1:nrow(allcontrols), size = nrow(allcases)*ratio, replace=F),] 
    control.samp$sampweight <- nrow(allcases)*ratio / nrow(allcontrols)
    
  } else if (samp=="sps") {  # simple probability sample of controls with known probability of selection for each individual
    
    allcontrols$sampprob <- runif(nrow(allcontrols), 0, 1) # generate probability of being selected
    control.samp <- allcontrols[sample(1:nrow(allcontrols), size = nrow(allcases)*ratio, prob = allcontrols$sampprob, replace=F),] # Sample control units
    control.samp$sampweight <- 1/control.samp$sampprob # Calculate Weights
    control.samp <- subset(control.samp, select = -sampprob) # Remove unneeded column  
    allcontrols <- subset(allcontrols, select = -sampprob) # Remove unneeded column 
 
  } else if (samp=="clustered1") { # single state cluster design in which clusters are sampled and all individuals within selected clusters are selected.          
    cluster <- aggregate(data.frame(popsize = allcontrols$cluster), list(cluster = allcontrols$cluster), length) # Calculate cluster (i.e. cluster) population size to determine cluster sampling probability (proportional to cluster population size)
    cluster$cls.sampprob <- cluster$popsize/nrow(allcontrols) # Calculate cluster sampling probability
    cluster.samp <- cluster[sample(1:nrow(cluster), size = round((nrow(allcases)*ratio/mean(table(data$cluster)))/2,0), prob = cluster$cls.sampprob, replace=F),] # Sample clusters using cluster sampling probability; note difficulty in arriving at desired sample size
    control.samp <- allcontrols[allcontrols$cluster %in% cluster.samp[,"cluster"],] # Sample all controls from each of the randomly sampled clusters
    control.samp <- merge(control.samp, cluster.samp, by="cluster") # Merge cluster characteristics with sampled controls
    control.samp$sampweight <- 1/(control.samp$cls.sampprob) # Calculate sampling weight
    control.samp <- subset(control.samp, select = -c(popsize, cls.sampprob)) # Remove unneeded column
    control.samp <- control.samp[sample(1:nrow(control.samp)), ] # Order randomly
    rm(cluster,cluster.samp) # Remove unneeded objects      
    
  } else if (samp=="clustered2") { # two stage cluster design in which cluster are sampled and individuals are sampled from within selected clusters.

    puma <- aggregate(data.frame(popsize = allcontrols$puma), list(puma = allcontrols$puma), length) # Calculate cluster (i.e. PUMA) population size to determine cluster sampling probability (proportional to cluster population size)
    puma$cls.sampprob <- puma$popsize/nrow(allcontrols) # Calculate cluster sampling probability
    puma.samp <- puma[sample(1:nrow(puma), size = 150, prob = puma$cls.sampprob, replace=F),] # Sample 150 clusters using cluster sampling probability
    control.samp <- allcontrols[allcontrols$puma %in% puma.samp[,"puma"],] %>% group_by(puma) %>% sample_n(150)# Randomly sample 150 controls from each of the 150 selected clusters
    control.samp <- merge(control.samp, puma.samp, by="puma") # Merge cluster characteristics with sampled controls
    control.samp$sampprob <- 150/control.samp$popsize # Calculate individual within-cluster sampling probability (i.e. 150 divided by cluster population size)
    control.samp$sampweight <- 1/(control.samp$cls.sampprob*control.samp$sampprob) # Calculate Sampling Weight
    control.samp <- subset(control.samp, select = -c(popsize,cls.sampprob, sampprob)) # Remove unneeded columns
    control.samp <- control.samp[sample(1:nrow(control.samp)), ] # Order randomly
    rm(puma,puma.samp) # Remove unneeded objects

  } else if (samp=="ACS") {
    
  } else if (samp=="NHANES") {
    
  }
  
  ####### PHASE 2: implement case-control - cumulative or density sampled, and analyse the data appropriately
  set.seed(seed+1)
  
  # CUMULATIVE CASE-CONTROL
  if (cctype=="cumulative") {
    
    # Apply design
    sample <- rbind(allcases, control.samp) 
    
    # Run model
    mod <- glm(Y ~ A + black + asian + hispanic + otherrace + age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_over64 + male + 
                 educ_ged + educ_hs + educ_somecollege + educ_associates + educ_bachelors + educ_advdegree, 
               data=sample, family='binomial', weights = sampweight)
    
    # Pull the main point estimate and CI
    est <- exp(coef(mod)[2])
    lower <- exp(coef(mod)[2] - 1.96*summary(mod)$coefficients[2,2])
    upper <- exp(coef(mod)[2] + 1.96*summary(mod)$coefficients[2,2])

    
  # DENSITY-SAMPLED CASE-CONTROL
  } else if (cctype=="density") {
    
    # Apply design
    presample <- rbind(allcases, control.samp) 
    sample <- ccwc(entry=0, exit=time, fail=Y, origin=0, controls=ratio, 
                   #match=list(), # use this argument for variables we want to match on
      include=list(A,black,asian,hispanic,otherrace,age_25_34,age_35_44,
                   age_45_54,age_55_64,age_over64,male,educ_ged,educ_hs,educ_somecollege,
                   educ_associates,educ_bachelors,educ_advdegree,sampweight), data=presample, silent=FALSE)
    
    # Run model
    mod <- clogit(Fail ~ A + black + asian + hispanic + otherrace + age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_over64 + male + 
                    educ_ged + educ_hs + educ_somecollege + educ_associates + educ_bachelors + educ_advdegree + strata(Set), 
                  data = sample, weights = sampweight, method = "efron")

    # Pull the main point estimate and CI
    est <- exp(coef(mod)[1])
    lower <- exp(coef(mod)[1] - 1.96*summary(mod)$coefficients[1,3])
    upper <- exp(coef(mod)[1] + 1.96*summary(mod)$coefficients[1,3])
  }
  
  # Return the sampled data, model object, point estimate, and CI
  return(list(sample=sample, mod=mod, est=est, lower=lower, upper=upper))
}


####
## Test function
####

# Warning: density sampling a long time to run.

study1 <- study(seed=1, cctype = "cumulative", samp = "srs", ratio = 1, data=pop)
round(c(study1$est, study1$lower, study1$upper),3) 

study2 <- study(seed=1, cctype = "density", samp = "srs", ratio = 1, data=pop) # warning has tied failure times. Look into this/check with Patrick.
round(c(study3$est, study3$lower, study3$upper),3) 


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
results$CIcover <- as.numeric(results$lower<=trueOR & results$upper>=trueOR)
round(mean(results$CIcover, na.rm=T)*100,1)

# Calculate Bias - average distance of point estimates away from true OR in repeated simulations
bias <- mean(results$est - trueOR) 
bias

# Calculate Variance - variance of point estimate of repeated simulations
variance <- var(results$est)
variance

# Calculate MSE - mean squared error of point estimate of repeated simulations
MSE <- mean((results$est - trueOR)^2)
MSE


# END