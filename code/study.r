################################################################################################
# NAME: study.r
# AUTHORS: Ellie Matthay, Catherine Li, Chris Rowe
# DATE STARTED: 2/20/18 
# PURPOSE: Define a function to apply various survey sampling and case control procedures to 
#            test the performance of different approaches to selecting controls from external 
#            data sources with complex sampling methods.
# UPDATES: 2/26/2018: CR updated 2-stage clustered design to flexibly handle any outcome frequency
#          2/27/2018  CR updated 2-stage clustered design to incorporate control:case ratio
#          2/27/2018 CR initialized main results with NA values and wrapped ccwc/clogit in try()
#          3/08/2018 CR fixed SRS sampling weights
#          3/14/2018 CR Added code to expand sample by person weights prior to risk set sampling
#          3/28/2018: CR normalized/scaled weights for population expansion prior to risks set
#                     sampling
#          4/12/2018: EM added extra print statements to help with debugging
#          4/16/2018: CR adding new sampling schemes related to exposure/covariates (exposure-
#                     based probability sample, age-stratified, race-stratified)
#          4/23/2018: CL added ccwc.weights function and updated function to have density sampling
#                     designs that use weight expansion or weights in control selection
################################################################################################

####
# Create a function that will apply the different study designs and analyses
####

study <- function(iteration, # iteration number for indexing runs and seeds
                  cctype = "cumulative", # case control type. Options will be:
                    # "cumulative"  for cumulative case-control
                    # "density expansion" for density-sampled with weight expansion
                    # "density weights" for density-sampled with weights in control selection
                  samp   = "srs", # sampling of controls. Options will be:
                    # "srs" for simple random sample
                    # "sps" for simple probability sample with known probability of selection for each individual
                    # "exp.ps" for probabiliy sample where sampling probability corresponds to exposure mechanism
                    # "clustered1" for single stage clustered design
                    # "clustered2" for two-stage clustered design 
                    # "stratified" for single stage stratified design
                    # "age.stratified" for age stratified design
                    # "race.stratified" for race stratified design
                  ratio = 1, # ratio of controls to cases
                  data, # argument to provide the population data. Population data will take the format of the pop data we've created. 
                  exposure, # name of exposure variable
                  outcome,  # name of outcome variable
                  timevar,   # name of the time variable corresponding to the given outcome frequency
                  seeds # rigorously produced random seeds
                   # any other options here TBD
) {
  
  print(paste0("Running iteration number ",iteration))
  
  print("Loading packages...")
  
  # Load packages on each compute node
  library("Epi") # for case-cohort and density sampling designs
  library("survival") # for clogit analysis
  library("dplyr") # for data management
  library("parallel") # for setting seeds
  
  print("Packages loaded. Setting seed...")
  
  # Set seed for entire process
  RNGkind("L'Ecuyer-CMRG")
  .Random.seed <- seeds[[iteration]]
  
  print("Seed set. Putting variables into standardized names...")
  
  ####### PHASE 1: source of cases and controls -- SRS, complex survey, etc. 
  
  # Put variables into standardized names
  data$Y <- data[[outcome]]
  data$A <- data[[exposure]]
  data$time <- data[[timevar]]
  
  print("Variables put in standardized names. Pulling all cases and controls and creating case weights...")
  
  print("Data summary:")
  print(summary(data))
  
  # Select cases and controls
  allcases    <- data[data$Y==1,]
  allcontrols <- data[data$Y==0,]

  # Create Case Weights
  allcases$sampweight <- 1
  
  print("Summary of allcases and allcontrols:")
  print(summary(allcases))
  print(summary(allcontrols))
  
  print("Cases and controls pulled. Case weights created. Entering control selection...")
  
  if (samp=="srs") { # simple random sample of controls
    
    control.samp <- allcontrols[sample(1:nrow(allcontrols), size = nrow(allcases)*ratio, replace=F),] 
    control.samp$sampweight <- 1/(nrow(allcases)*ratio / nrow(allcontrols))
    
  } else if (samp=="sps") {  # simple probability sample of controls with known probability of selection for each individual
    
    print("Simple probability sample. Generating probability of being selected for each control...")
    
    allcontrols$sampprob <- runif(nrow(allcontrols), 0, 1) # generate probability of being selected
    
    print("summary of allcontrols:")
    print(summary(allcontrols))
    
    print("Probability of selection generated. Sampling based on this probability...") 
    
    control.samp <- allcontrols[sample(1:nrow(allcontrols), size = nrow(allcases)*ratio, prob = allcontrols$sampprob, replace=F),] # Sample control units
    
    print("summary of control.samp:")
    print(summary(control.samp))
    
    print("Controls sampled. Creating sampling weights....")
    
    control.samp$sampweight <- 1/control.samp$sampprob # Calculate Weights
    
    print("Sampling weights created. Removing unneeded columns...") 
    
    control.samp <- subset(control.samp, select = -sampprob) # Remove unneeded column  
    allcontrols <- subset(allcontrols, select = -sampprob) # Remove unneeded column 
    
    print("Unneeded columns removed. Proceeding to analysis...")
 
  } else if (samp=="exp.ps") {  # probability sample where sampling probability corresponds exposure mechanism
 
    print("Exposure probability sample. Generating probability of being selected for each control...")
    
    allcontrols$sampprob <- plogis((log(1))*allcontrols$black + 
                                    (log(1.1))*allcontrols$asian + 
                                    (log(0.8))*allcontrols$hispanic +
                                    (log(0.9))*allcontrols$otherrace + 
                                    (log(1.1))*allcontrols$age_25_34 +
                                    (log(1.2))*allcontrols$age_35_44 +
                                    (log(1.3))*allcontrols$age_45_54 +
                                    (log(1.4))*allcontrols$age_55_64 +
                                    (log(3))*allcontrols$age_over64 +
                                    (log(3))*allcontrols$male +
                                    (log(3))*allcontrols$educ_ged +
                                    (log(1.2))*allcontrols$educ_hs + 
                                    (log(1.1))*allcontrols$educ_somecollege +
                                    (log(1))*allcontrols$educ_associates +
                                    (log(0.9))*allcontrols$educ_bachelors +
                                    (log(0.8))*allcontrols$educ_advdegree) # generate probability of being selected
    
    print("summary of allcontrols:")
    print(summary(allcontrols))
    
    print("Probability of selection generated. Sampling based on this probability...") 
    
    control.samp <- allcontrols[sample(1:nrow(allcontrols), size = nrow(allcases)*ratio, prob = allcontrols$sampprob, replace=F),] # Sample control units
    
    print("summary of control.samp:")
    print(summary(control.samp))
    
    print("Controls sampled. Creating sampling weights....")
    
    control.samp$sampweight <- 1/control.samp$sampprob # Calculate Weights
    
    print("Sampling weights created. Removing unneeded columns...") 
    
    control.samp <- subset(control.samp, select = -sampprob) # Remove unneeded column  
    allcontrols <- subset(allcontrols, select = -sampprob) # Remove unneeded column 
    
    print("Unneeded columns removed. Proceeding to analysis...")       
    
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
    puma.samp <- puma[sample(1:nrow(puma), size = ceiling(sqrt(nrow(allcases)*ratio)), prob = puma$cls.sampprob, replace=F),] # Sample clusters using cluster sampling probability
    control.samp <- allcontrols[allcontrols$puma %in% puma.samp[,"puma"],] %>% group_by(puma) %>% sample_n(ceiling(sqrt(nrow(allcases)*ratio)))# Randomly sample controls from each of the selected clusters
    control.samp <- merge(control.samp, puma.samp, by="puma") # Merge cluster characteristics with sampled controls
    control.samp$sampprob <- ceiling(sqrt(nrow(allcases)*ratio))/control.samp$popsize # Calculate individual within-cluster sampling probability (i.e. 150 divided by cluster population size)
    control.samp$sampweight <- 1/(control.samp$cls.sampprob*control.samp$sampprob) # Calculate Sampling Weight
    control.samp <- subset(control.samp, select = -c(popsize,cls.sampprob, sampprob)) # Remove unneeded columns
    control.samp <- control.samp[sample(1:nrow(control.samp)), ] # Order randomly
    rm(puma,puma.samp) # Remove unneeded objects

  } else if (samp=="stratified") {    
    
    allcontrols$strata2 <- as.numeric(cut(allcontrols$county, unique(quantile(allcontrols$county,seq(0,1,.1))), include.lowest=TRUE)) # Split counties into 8 strata
    stratainfo <- data.frame(table(allcontrols$strata2)) # Create dataframe for strata info for calculating sampling weights later
    stratainfo$size <- round((stratainfo$Freq/nrow(allcontrols))*(nrow(allcases)*ratio)) # Calculate sample size for each strata that is proportional to strata size
    colnames(stratainfo) <- c("strata2", "stratasize", "stratasampsize") # Rename strata data colunms for merging with sampled controls
    control.samp <- allcontrols[0,] # Create empty data.frame for samples
    for(i in 1:length(unique(allcontrols$strata2))) { # Sample controls proportional to strata size
      controls.strata <- allcontrols[allcontrols$strata2==i,]
      control.samp.strata <- controls.strata[sample(1:nrow(controls.strata), size = stratainfo$stratasampsize[i], replace=F),]   
      control.samp <- rbind(control.samp, control.samp.strata)
    }
    control.samp <- merge(control.samp, stratainfo, by="strata2") # Merge in strata info to calculate weights
    control.samp$sampweight <- 1/(control.samp$stratasampsize/control.samp$stratasize) # Calculate Sampling Weight 
    control.samp <- subset(control.samp, select = -c(strata2, stratasize, stratasampsize)) # Remove unneeded columns
    control.samp <- control.samp[sample(1:nrow(control.samp)), ] # Order randomly
    rm(control.samp.strata, controls.strata, stratainfo, i) # Remove unneeded objects
  
  } else if (samp=="age.stratified") {    
    
    print("Age stratified sample. Generating probability of being selected for each control...")
    
    allcontrols$sampprob[allcontrols$age_18_24==1] <- 1/(sum(allcontrols$age_18_24)/nrow(allcontrols)) # Calculate sampling probabilities proportional to inverse of age group frequency (i.e. rarer age groups are sampled more)
    allcontrols$sampprob[allcontrols$age_25_34==1] <- 1/(sum(allcontrols$age_25_34)/nrow(allcontrols))
    allcontrols$sampprob[allcontrols$age_35_44==1] <- 1/(sum(allcontrols$age_35_44)/nrow(allcontrols))
    allcontrols$sampprob[allcontrols$age_45_54==1] <- 1/(sum(allcontrols$age_45_54)/nrow(allcontrols))
    allcontrols$sampprob[allcontrols$age_55_64==1] <- 1/(sum(allcontrols$age_55_64)/nrow(allcontrols))
    allcontrols$sampprob[allcontrols$age_over64==1] <- 1/(sum(allcontrols$age_over64)/nrow(allcontrols))
    allcontrols$sampprob <- allcontrols$sampprob/sum(allcontrols$sampprob) # Scale sampling probabilities to fractions summing to 1
 
    print("summary of allcontrols:")
    print(summary(allcontrols))
    
    print("Probability of selection generated. Sampling based on this probability...") 
    
    control.samp <- allcontrols[sample(1:nrow(allcontrols), size = nrow(allcases)*ratio, prob = allcontrols$sampprob, replace=F),] # Sample control units
    
    print("summary of control.samp:")
    print(summary(control.samp))
    
    print("Controls sampled. Creating sampling weights....")
    
    control.samp$sampweight <- 1/control.samp$sampprob # Calculate Weights
    
    print("Sampling weights created. Removing unneeded columns...") 
    
    control.samp <- subset(control.samp, select = -sampprob) # Remove unneeded column  
    allcontrols <- subset(allcontrols, select = -sampprob) # Remove unneeded column 
    
    print("Unneeded columns removed. Proceeding to analysis...")       
    
  } else if (samp=="race.stratified") {    
    
    print("Race stratified sample. Generating probability of being selected for each control...")
    
    allcontrols$sampprob[allcontrols$white==1] <- 1/(sum(allcontrols$white)/nrow(allcontrols)) # Calculate sampling probabilities proportional to inverse of racial group frequency (i.e. rarer racial groups are sampled more)
    allcontrols$sampprob[allcontrols$black==1] <- 1/(sum(allcontrols$black)/nrow(allcontrols))
    allcontrols$sampprob[allcontrols$asian==1] <- 1/(sum(allcontrols$asian)/nrow(allcontrols))
    allcontrols$sampprob[allcontrols$hispanic==1] <- 1/(sum(allcontrols$hispanic)/nrow(allcontrols))
    allcontrols$sampprob[allcontrols$otherrace==1] <- 1/(sum(allcontrols$otherrace)/nrow(allcontrols))
    allcontrols$sampprob <- allcontrols$sampprob/sum(allcontrols$sampprob) # Scale sampling probabilities to fractions summing to 1
    
    print("summary of allcontrols:")
    print(summary(allcontrols))
    
    print("Probability of selection generated. Sampling based on this probability...") 
    
    control.samp <- allcontrols[sample(1:nrow(allcontrols), size = nrow(allcases)*ratio, prob = allcontrols$sampprob, replace=F),] # Sample control units
    
    print("summary of control.samp:")
    print(summary(control.samp))
    
    print("Controls sampled. Creating sampling weights....")
    
    control.samp$sampweight <- 1/control.samp$sampprob # Calculate Weights
    
    print("Sampling weights created. Removing unneeded columns...") 
    
    control.samp <- subset(control.samp, select = -sampprob) # Remove unneeded column  
    allcontrols <- subset(allcontrols, select = -sampprob) # Remove unneeded column 
    
    print("Unneeded columns removed. Proceeding to analysis...")       
    
  }
  
  
  
  ####### PHASE 2: implement case-control - cumulative or density sampled, and analyse the data appropriately
  
  print("Starting Phase 2. Initialize main results with NA values...")
  
  # Initialize main results with NA values
  est <- as.numeric(NA)
  lower <- as.numeric(NA)
  upper <- as.numeric(NA)
  truth <- as.numeric(NA)
  
  # CUMULATIVE CASE-CONTROL
  if (cctype=="cumulative") {
    
    print("Values initialized. rbinding cases and controls into sample data.frame...")
    
    # Apply design
    sample <- rbind(allcases, control.samp) 
    
    print("Combined data frame created. Running model...")

    print("Summary of data to be input to model:")
    print(summary(sample))
    
    # Run model
    mod <- glm(Y ~ A + black + asian + hispanic + otherrace + age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_over64 + male +
                 educ_ged + educ_hs + educ_somecollege + educ_associates + educ_bachelors + educ_advdegree,
               data=sample, family='binomial', weights = sampweight)
    
    print("Summary and coef of model:")
    print(summary(mod))
    
    print("Model run successfully. Pulling point estimate and CI...")
    
    # Pull the main point estimate and CI
    est <- exp(coef(mod)[2])
    lower <- exp(coef(mod)[2] - 1.96*summary(mod)$coefficients[2,2])
    upper <- exp(coef(mod)[2] + 1.96*summary(mod)$coefficients[2,2])
    
    print("Pulled point estimate and CI. Pulling true parameter...")
    
    # pull the relevant true parameter
    truth <- data[[paste0("trueOR.",outcome)]][1]
    
    print("True parameter pulled. Returning results...")
    
  # DENSITY-SAMPLED CASE-CONTROL - WEIGHT EXPANSION
  } else if (cctype=="density expansion") {
    
    # Apply design
    
    # Normalize & Scale Control Sampling Weights to Number of Controls in Population
    control.samp$sampweight = (control.samp$sampweight/sum(control.samp$sampweight))*nrow(allcontrols)
    
    # Combine Cases and Sampled Controls
    presample <- rbind(allcases, control.samp) 
    
    # Expand Sample Using Individual Person Weights
    presample.expnd <- presample[rep(row.names(presample),round(presample$sampweight,0)),]
    
    # Risk Set Sampling
    try(suppressWarnings(sample <- ccwc.weights(entry=0, exit=time, fail=Y, origin=0, controls=ratio, 
                   #match=list(), # use this argument for variables we want to match on
      include=list(A,black,asian,hispanic,otherrace,age_25_34,age_35_44,
                   age_45_54,age_55_64,age_over64,male,educ_ged,educ_hs,educ_somecollege,
                   educ_associates,educ_bachelors,educ_advdegree), data=presample.expnd, silent=FALSE)))
    
    # Run model
    try(mod <- clogit(Fail ~ A + black + asian + hispanic + otherrace + age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_over64 +
                    male + educ_ged + educ_hs + educ_somecollege + educ_associates + educ_bachelors + educ_advdegree + 
                    strata(Set),
                  data = sample, method = "efron"))

    # Pull the main point estimate and CI
    est <- exp(coef(mod)[1])
    lower <- exp(coef(mod)[1] - 1.96*summary(mod)$coefficients[1,3])
    upper <- exp(coef(mod)[1] + 1.96*summary(mod)$coefficients[1,3])
    
    # pull the relevant true parameter
    truth <- data[[paste0("trueIDR.",outcome)]][1]
  
  # DENSITY-SAMPLED CASE-CONTROL - WEIGHTS IN CONTROL SELECTION 
  } else if (cctype == "density weights") {
    
    # Apply Design
    
    #Combine Cases and Selected Controls
    presample <- rbind(allcases, control.samp) 
    
    # Risk Set Sampling with Weights
    try(suppressWarnings(sample <- ccwc.weights(entry=0, exit=time, fail=Y, origin=0, controls=ratio, weights = sampweight,
                                                #match=list(), # use this argument for variables we want to match on
                                                include=list(A,black,asian,hispanic,otherrace,age_25_34,age_35_44,
                                                             age_45_54,age_55_64,age_over64,male,educ_ged,educ_hs,educ_somecollege,
                                                             educ_associates,educ_bachelors,educ_advdegree, sampweight), data=presample, silent=FALSE)))
    
    # Run model - DOUBLE CHECK TO SEE IF WEIGHTS NEED TO BE INCLUDED IN REGRESSION
    try(mod <- clogit(Fail ~ A + black + asian + hispanic + otherrace + age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_over64 +
                        male + educ_ged + educ_hs + educ_somecollege + educ_associates + educ_bachelors + educ_advdegree + 
                        strata(Set),
                      data = sample, method = "efron"))
    
    # Pull the main point estimate and CI
    est <- exp(coef(mod)[1])
    lower <- exp(coef(mod)[1] - 1.96*summary(mod)$coefficients[1,3])
    upper <- exp(coef(mod)[1] + 1.96*summary(mod)$coefficients[1,3])
    
    # pull the relevant true parameter
    truth <- data[[paste0("trueIDR.",outcome)]][1]
  }
  
  # Return the sampled data, model object, point estimate, and CI
  return(list(est=est, lower=lower, upper=upper, truth=truth))
}

# END