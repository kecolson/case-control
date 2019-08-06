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
#          4/24/2018 CR updated sampling and cumulative code to allow for finite control survey
#                    sample sizes (i.e., argument "svysize") and method of obtaining sufficient
#                    number of controls when control survey sample size is less than number of
#                    controls needed.
#          4/23/2018: CL added ccwc.weights function and updated function to have density sampling
#                     designs that use weight expansion or weights in control selection
#          4/25/2018: CL updated density sampling designs to include 4 different methods (expand, 
#                     sample, model, unweighted)
#		       4/26/2018: CL updated cumulative sampling designs to include model and unweighted methods
#          5/1/2018: CR updated code so control weights are always rescaled to sum to number of controls
#                    in population; updated code to remove data, allcases, and allcontrols dataframes
#                    when not needed.
#          6/3/2018: CR updated cumulative code analysis methods (e.g., expand, sample, etc.)
#          11/7/2018: CR Modified density model and unweighted methods to risk-set sample only 
#                     controls, not cases and controls; added survey sampling design where
#                     exposed controls sampled with probability 0.75, unexposed controls sampled
#                     with 0.25; added function argument allowing for cases to be included in
#                     survey; capture SE as part of results
#          5/24/2019: CR added stratified.bias survey design
#          08//06/2019: Updated exp.ps2 design to explicilty sample 75% exposed and 25% unexposed
################################################################################################

####
# Create a function that will apply the different study designs and analyses
####

study <- function(iteration, # iteration number for indexing runs and seeds
                  cctype = "cumulative", # case control type. Options will be:
                    # "cumulative"  for cumulative case-control
                    # "density" for density-sampled (nested case-control design)
                  samp   = "srs", # sampling of controls. Options will be:
                    # "srs" for simple random sample
                    # "sps" for simple probability sample with known probability of selection for each individual
                    # "exp.ps" for probabiliy sample where sampling probability corresponds to exposure mechanism
                    # "exp.ps2" for probability sample where exposed have sampling probability of 0.75 and unexposed have sampling probability of 0.25
                    # "clustered1" for single stage clustered design
                    # "clustered2" for two-stage clustered design 
                    # "stratified.bias" for stratified design with bias
                    # "stratified" for single stage stratified design
                    # "age.stratified" for age stratified design
                    # "race.stratified" for race stratified design
                  svysize = "benchmark", # Size of survey from which controls are sampled
                    # "benchmark" sufficient respondents to allow for unique controls
                    # "small" 1000 respondents
                    # "medium" 5000 respondents
                    # "large" 15000 respondents
                    # "xl" 100,000 respondents
                    # "xxl" 400,000 respondents
                  svycase = FALSE, # Cases are excluded from survey used to obtain controls
                    # TRUE, cases are included in survey used to obtain controls
                  method = "expand", # Method of obtaining controls for when survey size is smaller than number of controls needed
                    # "expand" expand the sampled controls by their survey weights
                    # "sample" sample from controls with replacement and with sampling probability proportional to weights
                    # "model"  weights are incoporated into the model
                    # "unweighted"  weights from survey are not considered at all
                  ratio = 1, # ratio of controls to cases
                  data, # argument to provide the population data. Population data will take the format of the pop data we've created. 
                  exposure, # name of exposure variable
                  outcome,  # name of outcome variable
                  timevar,   # name of the time variable corresponding to the given outcome frequency
                  seeds # rigorously produced random seeds
                   # any other options here TBD
) {
  
  print(paste0("Running iteration number ",iteration))
  
  print("Loading packages and ccwc.weights() function ...")
  
  # Load packages on each compute node
  library("Epi") # for case-cohort and density sampling designs
  library("survival") # for clogit analysis
  library("dplyr") # for data management
  library("parallel") # for setting seeds
  
  source('code/ccwc.weights.R')
  
  print("Packages and function loaded. Setting seed...")
  
  # Set seed for entire process
  RNGkind("L'Ecuyer-CMRG")
  .Random.seed <- seeds[[iteration]]
  
  print("Seed set. Putting variables into standardized names...")
  
  ####### PHASE 1: source of cases and controls -- SRS, complex survey, etc. 
  
  # Put variables into standardized names
  data$Y <- data[[outcome]]
  data$A <- data[[exposure]]
  data$time <- data[[timevar]]
  
  data <- data[,!names(data) %in% c(outcome, exposure, timevar)]

  print("Variables put in standardized names. Pulling all cases and controls and creating case weights...")
  
  print("Data summary:")
  print(summary(data))
  
  # Select cases and controls
  allcases    <- data[data$Y==1,]
  allcontrols <- data[data$Y==0,]
  
  # Extract Truth
  if (cctype=="cumulative") { 
    truth <- data[[paste0("trueOR.",outcome)]][1]
  } else if (cctype=="density") {
    truth <- data[[paste0("trueIDR.",outcome)]][1]
  }
  
  # Extract number of cases and controls
  Ncases <- nrow(allcases)
  Ncontrols <- nrow(allcontrols)

  # Create Case Weights
  allcases$sampweight <- 1
  
  print("Summary of allcases and allcontrols:")
  print(summary(allcases))
  print(summary(allcontrols))
  
  print("Cases and controls pulled. Case weights created. Entering control selection...")
  
  if (svysize=="benchmark") { # Obtain size of survey from svysize argument
    svysizenum <- Ncases*ratio
  } else if (svysize=="small") {
    svysizenum <- 1000   
  } else if (svysize=="medium") {  
    svysizenum <- 5000   
  } else if (svysize=="large") {  
    svysizenum <- 15000   
  } else if (svysize=="xl") {
    svysizenum <- 100000
  } else if (svysize=="xxl") {
    svysizenum <- 400000
  }
  
  if (svycase==TRUE) { # Determine whether cases should be included or excluded from survey from svycase argument
    svypop <- data # Cases included in survey
  } else if (svycase==FALSE) { 
    svypop <- allcontrols # Cases excluded from survey   
  }
    
    # Generate Survey Population Size
    Nsvypop <- nrow(svypop)
    
    # Remove Full Dataset - No Longer Needed
    rm(data)
    
    # Remove allcontrols dataframe - no longer needed
    rm(allcontrols)  
    
  if (samp=="srs") { # simple random sample of controls
    
    control.samp <- svypop[sample(1:Nsvypop, size = svysizenum, replace=F),] 
    control.samp$sampweight <- 1/(nrow(control.samp)/Nsvypop)
    
  } else if (samp=="sps") {  # simple probability sample of controls with known probability of selection for each individual
    
    print("Simple probability sample. Generating probability of being selected for each control...")
    
    svypop$sampprob <- runif(Nsvypop, 0, 1) # generate probability of being selected
    
    print("summary of svypop:")
    print(summary(svypop))
    
    print("Probability of selection generated. Sampling based on this probability...") 
    
    control.samp <- svypop[sample(1:Nsvypop, size = svysizenum, prob = svypop$sampprob, replace=F),] # Sample control units
    
    print("summary of control.samp:")
    print(summary(control.samp))
    
    print("Controls sampled. Creating sampling weights....")
    
    control.samp$sampweight <- 1/control.samp$sampprob # Calculate Weights
    
    print("Sampling weights created. Removing unneeded columns...") 
    
    control.samp <- subset(control.samp, select = -sampprob) # Remove unneeded column  
    svypop <- subset(svypop, select = -sampprob) # Remove unneeded column 
    
    print("Unneeded columns removed. Proceeding to analysis...")
 
  } else if (samp=="exp.ps") {  # probability sample where sampling probability corresponds exposure mechanism
 
    print("Exposure probability sample. Generating probability of being selected for each control...")
    
    svypop$sampprob <- plogis((log(1))*svypop$black + 
                                    (log(1.1))*svypop$asian + 
                                    (log(0.8))*svypop$hispanic +
                                    (log(0.9))*svypop$otherrace + 
                                    (log(1.1))*svypop$age_25_34 +
                                    (log(1.2))*svypop$age_35_44 +
                                    (log(1.3))*svypop$age_45_54 +
                                    (log(1.4))*svypop$age_55_64 +
                                    (log(3))*svypop$age_over64 +
                                    (log(3))*svypop$male +
                                    (log(3))*svypop$educ_ged +
                                    (log(1.2))*svypop$educ_hs + 
                                    (log(1.1))*svypop$educ_somecollege +
                                    (log(1))*svypop$educ_associates +
                                    (log(0.9))*svypop$educ_bachelors +
                                    (log(0.8))*svypop$educ_advdegree) # generate probability of being selected
    
    print("summary of svypop:")
    print(summary(svypop))
    
    print("Probability of selection generated. Sampling based on this probability...") 
    
    control.samp <- svypop[sample(1:Nsvypop, size = svysizenum, prob = svypop$sampprob, replace=F),] # Sample control units
    
    print("summary of control.samp:")
    print(summary(control.samp))
    
    print("Controls sampled. Creating sampling weights....")
    
    control.samp$sampweight <- 1/control.samp$sampprob # Calculate Weights
    
    print("Sampling weights created. Removing unneeded columns...") 
    
    control.samp <- subset(control.samp, select = -sampprob) # Remove unneeded column  
    svypop <- subset(svypop, select = -sampprob) # Remove unneeded column 
    
    print("Unneeded columns removed. Proceeding to analysis...")       

  } else if (samp=="exp.ps2") {  # for probability sample where exposed have sampling probability of 0.75 and unexposed have sampling probability of 0.25
    
    print("Exposure probability sample. Generating number of exposed and unexposed to be sample...")
    n.exp <- round(0.75*svysizenum)
    n.unexp <- svysizenum - n.exp

    print("Numbers to be sampled generated. Sampling based on these numbers...") 
    
    id.exp <- sample(which(svypop$A==1),size=n.exp,replace=F) # Sample exposed rows
    id.unexp <- sample(which(svypop$A==0),size=n.unexp,replace=F) # Sample unexposed rows
    control.samp <- svypop[c(id.exp,id.unexp),] # Construct sample
    
    print("summary of control.samp:")
    print(summary(control.samp))
    
    print("Controls sampled. Generating sampling weights....")
    #Note: weights are calculated as pr(A=a in source population)/Pr(A=a in sample)
    control.samp$sampweight[control.samp$A==1] <- prop.table(table(svypop$A))[2]/prop.table(table(control.samp$A))[2] # Calculate Weights for exposed
    control.samp$sampweight[control.samp$A==0] <- prop.table(table(svypop$A))[1]/prop.table(table(control.samp$A))[1] # Calculate Weights for unexposed
    
    print("Sampling weights created. Proceeding to analysis...") 
    
    rm(n.exp,n.unexp,id.exp,id.unexp) # Remove unneeded objects      
    
  } else if (samp=="clustered1") { # single state cluster design in which clusters are sampled and all individuals within selected clusters are selected.          
    cluster <- aggregate(data.frame(popsize = svypop$cluster), list(cluster = svypop$cluster), length) # Calculate cluster (i.e. cluster) population size to determine cluster sampling probability (proportional to cluster population size)
    cluster$cls.sampprob <- cluster$popsize/Nsvypop # Calculate cluster sampling probability
    cluster.samp <- cluster[sample(1:nrow(cluster), size = round((svysizenum/mean(table(svypop$cluster)))/1.84,0), prob = cluster$cls.sampprob, replace=F),] # Sample clusters using cluster sampling probability; note difficulty in arriving at desired sample size
    control.samp <- svypop[svypop$cluster %in% cluster.samp[,"cluster"],] # Sample all controls from each of the randomly sampled clusters
    control.samp <- merge(control.samp, cluster.samp, by="cluster") # Merge cluster characteristics with sampled controls
    control.samp$sampweight <- 1/(control.samp$cls.sampprob) # Calculate sampling weight
    control.samp <- subset(control.samp, select = -c(popsize, cls.sampprob)) # Remove unneeded column
    control.samp <- control.samp[sample(1:nrow(control.samp)), ] # Order randomly
    rm(cluster,cluster.samp) # Remove unneeded objects      
    
  } else if (samp=="clustered2") { # two stage cluster design in which cluster are sampled and individuals are sampled from within selected clusters.

    puma <- aggregate(data.frame(popsize = svypop$puma), list(puma = svypop$puma), length) # Calculate cluster (i.e. PUMA) population size to determine cluster sampling probability (proportional to cluster population size)
    puma$cls.sampprob <- puma$popsize/Nsvypop # Calculate cluster sampling probability
    puma.samp <- puma[sample(1:nrow(puma), size = ceiling(sqrt(svysizenum)/3), prob = puma$cls.sampprob, replace=F),] # Sample clusters using cluster sampling probability
    control.samp <- svypop[svypop$puma %in% puma.samp[,"puma"],] %>% group_by(puma) %>% sample_n(ceiling(sqrt(svysizenum)*3))# Randomly sample controls from each of the selected clusters
    control.samp <- merge(control.samp, puma.samp, by="puma") # Merge cluster characteristics with sampled controls
    control.samp$sampprob <- ceiling(sqrt(svysizenum)*3)/control.samp$popsize # Calculate individual within-cluster sampling probability (i.e. 150 divided by cluster population size)    
    control.samp$sampweight <- 1/(control.samp$cls.sampprob*control.samp$sampprob) # Calculate Sampling Weight
    control.samp <- subset(control.samp, select = -c(popsize,cls.sampprob, sampprob)) # Remove unneeded columns
    control.samp <- control.samp[sample(1:nrow(control.samp)), ] # Order randomly
    rm(puma,puma.samp) # Remove unneeded objects

  } else if (samp=="stratified.bias") { # two stage stratified design in which individuals are sampled from within strata inversely proportional to size of strata
    
    print("Biased stratified sample. Sample equal number of individuals from each strata")
    
    # Identify number of individuals to be sampled per strata (balanced across strata)
    samp.per.strata <- round(svysizenum/length(unique(svypop[[paste0("strata.bias.",exposure)]])))
    svypop$id <- 1:nrow(svypop) # Create ID variable for all individuals survey pool
    
    # Sample equal number of individuals from each strata; note this step just samples IDs
    control.samp.id <- c()
    for(i in unique(svypop[[paste0("strata.bias.",exposure)]])){
      sample <- sample(svypop$id[svypop[[paste0("strata.bias.",exposure)]]==i],size=samp.per.strata,replace=F)
      control.samp.id <- c(control.samp.id,sample)
    }
    
    print("IDs sampled. Calculate sampling probability for each strata.")
    
    # Calculate the sampling probability for individuals, by strata, equal to number of individuals sampled from each strata divided by strata size
    strata <- data.frame(table(svypop[[paste0("strata.bias.",exposure)]]))
    strata$sampprob <- samp.per.strata/strata$Freq
    strata$Freq <- NULL
    names(strata) <- c(paste0("strata.bias.",exposure),"sampprob")

    
    
    print("Sampling probabilities calculated for later merging. Construct sample from sampled IDs") 
    
    # Construct sample from previously sampled IDs
    control.samp <- svypop[svypop$id %in% control.samp.id,] # Sample control units
    
    print("Controls sampled. Merging in sampling probabilities and Creating sampling weights....")
    
    # Merge in strata-specific sampling probabilities
    control.samp <- merge(control.samp,strata,by=paste0("strata.bias.",exposure))
    control.samp$sampweight <- 1/control.samp$sampprob # Calculate Weights
    
    print("Sampling weights created. Removing unneeded columns...") 
    
    control.samp <- subset(control.samp, select = -c(sampprob,id)) # Remove unneeded column  
    svypop <- subset(svypop, select = -id) # Remove unneeded column 
    
    print("Unneeded columns removed. Proceeding to analysis...")      

  } else if (samp=="stratified") {    
    
    svypop$strata2 <- as.numeric(cut(svypop$county, unique(quantile(svypop$county,seq(0,1,.1))), include.lowest=TRUE)) # Split counties into 8 strata
    stratainfo <- data.frame(table(svypop$strata2)) # Create dataframe for strata info for calculating sampling weights later
    stratainfo$size <- round((stratainfo$Freq/Nsvypop)*(svysizenum)) # Calculate sample size for each strata that is proportional to strata size
    colnames(stratainfo) <- c("strata2", "stratasize", "stratasampsize") # Rename strata data colunms for merging with sampled controls
    control.samp <- svypop[0,] # Create empty data.frame for samples
    for(i in 1:length(unique(svypop$strata2))) { # Sample controls proportional to strata size
      controls.strata <- svypop[svypop$strata2==i,]
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
    
    svypop$sampprob[svypop$age_18_24==1] <- 1/(sum(svypop$age_18_24)/Nsvypop) # Calculate sampling probabilities proportional to inverse of age group frequency (i.e. rarer age groups are sampled more)
    svypop$sampprob[svypop$age_25_34==1] <- 1/(sum(svypop$age_25_34)/Nsvypop)
    svypop$sampprob[svypop$age_35_44==1] <- 1/(sum(svypop$age_35_44)/Nsvypop)
    svypop$sampprob[svypop$age_45_54==1] <- 1/(sum(svypop$age_45_54)/Nsvypop)
    svypop$sampprob[svypop$age_55_64==1] <- 1/(sum(svypop$age_55_64)/Nsvypop)
    svypop$sampprob[svypop$age_over64==1] <- 1/(sum(svypop$age_over64)/Nsvypop)
    svypop$sampprob <- svypop$sampprob/sum(svypop$sampprob) # Scale sampling probabilities to fractions summing to 1
 
    print("summary of svypop:")
    print(summary(svypop))
    
    print("Probability of selection generated. Sampling based on this probability...") 
    
    control.samp <- svypop[sample(1:Nsvypop, size = svysizenum, prob = svypop$sampprob, replace=F),] # Sample control units
    
    print("summary of control.samp:")
    print(summary(control.samp))
    
    print("Controls sampled. Creating sampling weights....")
    
    control.samp$sampweight <- 1/control.samp$sampprob # Calculate Weights
    
    print("Sampling weights created. Removing unneeded columns...") 
    
    control.samp <- subset(control.samp, select = -sampprob) # Remove unneeded column  
    svypop <- subset(svypop, select = -sampprob) # Remove unneeded column 
    
    print("Unneeded columns removed. Proceeding to analysis...")       
    
  } else if (samp=="race.stratified") {    
    
    print("Race stratified sample. Generating probability of being selected for each control...")
    
    svypop$sampprob[svypop$white==1] <- 1/(sum(svypop$white)/Nsvypop) # Calculate sampling probabilities proportional to inverse of racial group frequency (i.e. rarer racial groups are sampled more)
    svypop$sampprob[svypop$black==1] <- 1/(sum(svypop$black)/Nsvypop)
    svypop$sampprob[svypop$asian==1] <- 1/(sum(svypop$asian)/Nsvypop)
    svypop$sampprob[svypop$hispanic==1] <- 1/(sum(svypop$hispanic)/Nsvypop)
    svypop$sampprob[svypop$otherrace==1] <- 1/(sum(svypop$otherrace)/Nsvypop)
    svypop$sampprob <- svypop$sampprob/sum(svypop$sampprob) # Scale sampling probabilities to fractions summing to 1
    
    print("summary of svypop:")
    print(summary(svypop))
    
    print("Probability of selection generated. Sampling based on this probability...") 
    
    control.samp <- svypop[sample(1:Nsvypop, size = svysizenum, prob = svypop$sampprob, replace=F),] # Sample control units
    
    print("summary of control.samp:")
    print(summary(control.samp))
    
    print("Controls sampled. Creating sampling weights....")
    
    control.samp$sampweight <- 1/control.samp$sampprob # Calculate Weights
    
    print("Sampling weights created. Removing unneeded columns...") 
    
    control.samp <- subset(control.samp, select = -sampprob) # Remove unneeded column  
    svypop <- subset(svypop, select = -sampprob) # Remove unneeded column 
    
    print("Unneeded columns removed. Proceeding to analysis...")       
    
  }
  
  # Normalize & Scale Control Sampling Weights to Number of people in underlying survey Population
  control.samp$sampweight = (control.samp$sampweight/sum(control.samp$sampweight))*(Nsvypop)  
  
  # Change control sample Y to zero and time to 3652.5 (this is only relevant when svycase=TRUE and cases are included in the survey
  control.samp$time <- svypop$time[svypop$Y==0][1]  
  control.samp$Y <- 0

  
  ####### PHASE 2: implement case-control - cumulative or density sampled, and analyse the data appropriately
  
  print("Starting Phase 2. Initialize main results with NA values...")
  
  # Initialize main results with NA values
  est <- as.numeric(NA)
  lower <- as.numeric(NA)
  upper <- as.numeric(NA)
  se <- as.numeric(NA)
  
  
  # CUMULATIVE CASE-CONTROL
  if (cctype=="cumulative") {
    
    print("Values initialized. Processing Controls...")
    
    if (method=="expand") { # Obtain sufficient number of controls 

      print("Expanding Control Sample to obtain sufficient controls for study design...")
      
      # Expand Control Sample using sampling weights
      control.samp.expnd <- control.samp[rep(row.names(control.samp),round(control.samp$sampweight,0)),]
      
      # Sample Controls
      new.control.samp <- control.samp.expnd[sample(1:nrow(control.samp.expnd), size = Ncases*ratio, replace=F),]
      
      # Clear Sampling Weights
      new.control.samp$sampweight <- 1
      
      # Combine Cases/Controls
      sample <- rbind(allcases, new.control.samp) 
      
    } else if (method=="sample") {
      
      print("Sampling from control sample with replacement to obtain sufficient controls for study design...")
        
      # Sample Controls with replacement and sampling probability proportional to weights
      new.control.samp <- control.samp[sample(1:nrow(control.samp), size = Ncases*ratio, prob=control.samp$sampweight, replace=T),]
      
      # Clear Sampling Weights
      new.control.samp$sampweight <- 1

      # Combine Cases/Controls
      sample <- rbind(allcases, new.control.samp) 
      
    } else if (method=="model") {

    	print("No control expansion/resampling, weights will be accounted for in analysis...")

    	# Combine Cases/Controls
    	sample <- rbind(allcases, control.samp)

    } else if (method=="unweighted") {

    	print("Survey weights removed for unweighted analysis...")

    	# Make all control weights 1
    	control.samp$sampweight <- 1

    	# Combine Cases/Controls
      sample <- rbind(allcases, control.samp)

    }
  
    # Remove allcases dataframe - no longer needed
    rm(allcases)
    
    
    # Apply design

    
    print("Combined case/control data frame created. Running model...")

    print("Summary of data to be input to model:")
    print(summary(sample))
    
    # Run model
    mod <- glm(Y ~ A + black + asian + hispanic + otherrace + age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_over64 + male +
                 educ_ged + educ_hs + educ_somecollege + educ_associates + educ_bachelors + educ_advdegree,
               data=sample, family='binomial', weights = sampweight)
    
    print("Summary and coef of model:")
    print(summary(mod))
    
    print("Model run successfully. Pulling point estimate and CI...")
    
    # Pull the main point estimate and CI and SE
    est <- exp(coef(mod)[2])
    lower <- exp(coef(mod)[2] - 1.96*summary(mod)$coefficients[2,2])
    upper <- exp(coef(mod)[2] + 1.96*summary(mod)$coefficients[2,2])
    se <- summary(mod)$coefficients[2,2]
    
    print("Pulled point estimate and CI. Returning results...")
    
    
  # DENSITY-SAMPLED CASE-CONTROL - WEIGHT EXPANSION
  } else if (cctype=="density") {
    
    if (method=="expand"){
      
      print("Expanding Control Sample to obtain sufficient controls for study design...")
      
      # Combine Cases and Sampled Controls
      presample <- rbind(allcases, control.samp) 
      
      # Remove allcases dataframe - no longer needed
      rm(allcases)
      
      # Expand Sample Using Individual Person Weights
      presample.expnd <- presample[rep(row.names(presample),round(presample$sampweight,0)),]
      
      print("Expansion complete. Running ccwc...")
      
      # Risk Set Sampling
      try(sample <- ccwc.weights(entry=0, exit=time, fail=Y, origin=0, controls=ratio, weights=NULL,
                                                  #match=list(), # use this argument for variables we want to match on
                                                  include=list(A,black,asian,hispanic,otherrace,age_25_34,age_35_44,
                                                               age_45_54,age_55_64,age_over64,male,educ_ged,educ_hs,educ_somecollege,
                                                               educ_associates,educ_bachelors,educ_advdegree), data=presample.expnd, silent=FALSE))
      
      print("ccwc complete. Running unweighted model...")
      
      # Run model
      try(mod <- clogit(Fail ~ A + black + asian + hispanic + otherrace + age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_over64 +
                          male + educ_ged + educ_hs + educ_somecollege + educ_associates + educ_bachelors + educ_advdegree + 
                          strata(Set),
                        data = sample, method = "efron"))
      
      print("model completed. Storing results...")
      
      # Pull the main point estimate and CI and SE
      est <- exp(coef(mod)[1])
      lower <- exp(coef(mod)[1] - 1.96*summary(mod)$coefficients[1,3])
      upper <- exp(coef(mod)[1] + 1.96*summary(mod)$coefficients[1,3])
      se <- summary(mod)$coefficients[1,3]
      
    } else if (method=="sample") {
      
      print("running ccwc with resampling of controls proportional to survey weights")
      
      #Combine Cases and Selected Controls
      presample <- rbind(allcases, control.samp) 
      
      # Remove allcases dataframe - no longer needed
      rm(allcases)
      
      # Risk Set Sampling with Weights
      try(sample <- ccwc.weights(entry=0, exit=time, fail=Y, origin=0, controls=ratio, weights = sampweight,
                                                  #match=list(), # use this argument for variables we want to match on
                                                  include=list(A,black,asian,hispanic,otherrace,age_25_34,age_35_44,
                                                               age_45_54,age_55_64,age_over64,male,educ_ged,educ_hs,educ_somecollege,
                                                               educ_associates,educ_bachelors,educ_advdegree, sampweight), data=presample, silent=FALSE))
      
      print("ccwc complete. Running unweighted model...")
      
      # Run model
      try(mod <- clogit(Fail ~ A + black + asian + hispanic + otherrace + age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_over64 +
                          male + educ_ged + educ_hs + educ_somecollege + educ_associates + educ_bachelors + educ_advdegree + 
                          strata(Set),
                        data = sample, method = "efron"))
      
      print("model completed. Storing results...")
      
      # Pull the main point estimate and CI and SE
      est <- exp(coef(mod)[1])
      lower <- exp(coef(mod)[1] - 1.96*summary(mod)$coefficients[1,3])
      upper <- exp(coef(mod)[1] + 1.96*summary(mod)$coefficients[1,3])
      se <- summary(mod)$coefficients[1,3]
      
    } else if (method=="model") {
      
      print("risk-set sampling with method=='model' (i.e. without regard for sampling weights)")
      
      # Create risk set strata for cases
      allcases$Set <- 1:Ncases
      
      # Risk set sample controls (which is equivalent to randomly selecting controls with replacement)
      rss.controls <- control.samp[sample(1:nrow(control.samp),size=Ncases*ratio,replace=T),]
      
      # Create risk set strata for controls
      rss.controls$Set <- rep(1:Ncases,ratio)
      
      # Combine cases and risk-set sample controls to create sample
      sample <- rbind(allcases, rss.controls)
      
      # Create "Fail" Variable
      sample$Fail <- sample$Y

      # Remove allcases dataframe - no longer needed
      rm(allcases, rss.controls)
      
      print("risk set sampling complete. Running weighted model...")
      
      # Run model
      try(mod <- clogit(Fail ~ A + black + asian + hispanic + otherrace + age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_over64 +
                          male + educ_ged + educ_hs + educ_somecollege + educ_associates + educ_bachelors + educ_advdegree + 
                          strata(Set), weights = sampweight,
                        data = sample, method = "efron"))
      
      print("model completed. Storing results...")
      
      # Pull the main point estimate and CI and SE
      est <- exp(coef(mod)[1])
      lower <- exp(coef(mod)[1] - 1.96*summary(mod)$coefficients[1,3])
      upper <- exp(coef(mod)[1] + 1.96*summary(mod)$coefficients[1,3])
      se <- summary(mod)$coefficients[1,3]
      
    } else if (method=="unweighted") {
      
      print("risk-set sampling with method=='unweighted' (i.e. without regard for sampling weights)")
      
      # Create risk set strata for cases
      allcases$Set <- 1:Ncases
      
      # Risk set sample controls (which is equivalent to randomly selecting controls with replacement)
      rss.controls <- control.samp[sample(1:nrow(control.samp),size=Ncases*ratio,replace=T),]
      
      # Create risk set strata for controls
      rss.controls$Set <- rep(1:Ncases,ratio)
      
      # Combine cases and risk-set sample controls to create sample
      sample <- rbind(allcases, rss.controls)
      
      # Create "Fail" Variable
      sample$Fail <- sample$Y
      
      # Remove allcases dataframe - no longer needed
      rm(allcases, rss.controls)
      
      print("risk set sampling complete. Running unweighted model...")
      
      # Run model
      try(mod <- clogit(Fail ~ A + black + asian + hispanic + otherrace + age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_over64 +
                          male + educ_ged + educ_hs + educ_somecollege + educ_associates + educ_bachelors + educ_advdegree + 
                          strata(Set),
                        data = sample, method = "efron"))
      
      print("model completed. Storing results...")
      
      # Pull the main point estimate and CI and SE
      est <- exp(coef(mod)[1])
      lower <- exp(coef(mod)[1] - 1.96*summary(mod)$coefficients[1,3])
      upper <- exp(coef(mod)[1] + 1.96*summary(mod)$coefficients[1,3])
      se <- summary(mod)$coefficients[1,3]
      
    }
    
  }
  
  # Return the sampled data, model object, point estimate, and CI
  return(list(est=est, lower=lower, upper=upper, se=se, truth=truth))
}

# END