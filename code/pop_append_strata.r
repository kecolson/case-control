################################################################################################
# NAME: pop_append.r
# AUTHORS: Ellie Matthay, Catherine Li, Chris Rowe
# DATE STARTED: 5/22/2019 
# PURPOSE: Assign each record in populataion.data.csv to a strata, where smaller strata have
#          higher proportion of exposed individuals (in order to generate bias when weights are 
#         not incorporated in analysis)
# UPDATES: 
################################################################################################

# Clear workspace
rm(list=ls())

# Load Packages.
library(lpSolve)

# Set Working Directory
#setwd("~/Documents/GitHub/case-control") # Chris's Directory
#setwd("C:/Users/kecolson/Google Drive/simulation/case-control-other") # Ellie's directory
setwd("/Users/cxli/Documents/case-control") # Catherine's directory


# Import Population Data from local drive (most updated is on GitHub)
data <- read.csv("data/population.data.csv", stringsAsFactors = F)


# Generage biased stratas

# Set Seed
set.seed(1234)
popsize <- nrow(data)

# Generate stratas Randomly of Population Sizes Ranging from 20,000 to 80,000 for ~20 stratas total (i.e., strata population sizes are uniformly selected from 20,000 to 80,000 for population of 1M)
strata.pops <- c()
sum <- 0
while (sum<popsize){
  single.strata <- sample(20000:80000,1)
  strata.pops <- c(strata.pops, single.strata)
  sum <- sum(strata.pops)
  if (sum>popsize){
    diff <- sum - popsize
    strata.pops[length(strata.pops)] <- strata.pops[length(strata.pops)] - diff
  }
}

# Remove unnecessary data
rm(single.strata,sum,diff,popsize)

# Organize Data
strata.pops <- sort(strata.pops,decreasing=F) # Sort in descending order
strata.data <- data.frame(1:length(strata.pops),strata.pops) # Assign strata ID 
names(strata.data) <- c("strata","popsize")
rm(strata.pops)

# Calculate appropriate exposed and unexposed sampling probabilities for each strata using constrained linear optimization
# Overview of constraints:
  # Weighted sum of strata sizes (where weights are sampling probabilities of exposed in each strata) equals total number of exposed in full population
  # Largest strata sampling probability must be between total exposure frequency and 0.99.
  # Smallest strata sampling probability must be between 0.01 and total exposure frequency.
  # The ratio of sampling probabilities for a strata and the subsequent strata (in descending order of size) must be at least 1.01 and less than 1.10


  # Set up Vector of Coefficients as strata Population Sizes
  f.obj  <- strata.data[,2]
  
  # Set up constraints
  
    # Total sum and bound constraints
    f.con <- matrix(c(strata.data[,2], # Weighted sum equals number of exposed
                      c(1,rep(0,length(strata.data[,2])-1)), # lower bound for sampling probability of exposed within largest strata 
                      c(1,rep(0,length(strata.data[,2])-1)), # upper bound for sampling probability of exposed within largest strata 
                      c(rep(0,length(strata.data[,2])-1),1), # upper bound for sampling probability of exposed within smallest strata 
                      c(rep(0,length(strata.data[,2])-1),1)), nrow=5,byrow=TRUE) # lower bound for sampling probability of exposed within smallest strata 
  
    # Initialize constraints that larger stratas have larger sampling probability of exposed (between 1% and 10% larger)
    order.constraints1 <- matrix(rep(0,length(strata.data[,2])*(length(strata.data[,2])-1)),nrow=nrow(strata.data)-1,ncol=nrow(strata.data))
    for(i in 1:nrow(order.constraints1)){
      order.constraints1[i,i] <- 1
      order.constraints1[i,i+1] <- -1.05 # Ratio of exposed sampling probabilities in one strata to its smaller neighbor is at least 1.01
    }
  
    order.constraints2 <- matrix(rep(0,length(strata.data[,2])*(length(strata.data[,2])-1)),nrow=nrow(strata.data)-1,ncol=nrow(strata.data))
    for(i in 1:nrow(order.constraints2)){
      order.constraints2[i,i] <- 1
      order.constraints2[i,i+1] <- -1.06 # Ratio of exposed sampling probabilities in one strata to its smaller neighbor is less than 1.1
    }
    
    # Combine two sets of constraints into single matrix
    f.con <- rbind(f.con,order.constraints1,order.constraints2)
  
  # Set up constraint directions
  f.dir <- c("=",">","<","<",">",rep(">",length(strata.data[,2])-1),rep("<",length(strata.data[,2])-1))
  
  # Initialize constraint values and calculate/save sampling probabilities for each exposure frequency
  
    # Initialize empty matrix for all exposed and unexposed sampling probabilities
    sp.exp <- matrix(rep(NA,length(strata.data[,2])*4),nrow=length(strata.data[,2]),ncol=4)
    
    # Loop through exposure frequencies and calculate sampling probabilities for exposed and unexposed for each strata size
    i = 1
    for(exp in c("A.50", "A.20", "A.10", "A.05")){
    
      # Set Up constraint values
      f.rhs <- c(sum(data[,exp]),sum(data[,exp])/nrow(data),0.99,sum(data[,exp])/nrow(data),0.01,rep(0,length(strata.data[,2])-1),rep(0,length(strata.data[,2])-1))
      
      # Run and save sampling probabilities
      sp.exp[,i] <- lp("max", f.obj, f.con, f.dir, f.rhs)$solution
      
      i = i + 1
    }

    # Remove data
    rm(f.con,order.constraints1,order.constraints2,exp,f.dir,f.obj,f.rhs,i)

# Assess solution
plot(sp.exp[,1], type="l")
plot(sp.exp[,2], type="l")
plot(sp.exp[,3], type="l")
plot(sp.exp[,4], type="l")

# Calculate exposed and unexposed sampling sizes for each strata and exposure frequency
clpop.exp <- sp.exp*strata.data[,2]
  clpop.exp[1:5,] <- ceiling(clpop.exp[1:5,])
  clpop.exp[6:nrow(sp.exp),] <- round(clpop.exp[6:nrow(sp.exp),])
clpop.unexp <- strata.data[,2] - clpop.exp
rm(sp.exp)


# For each exposure frequency, sample Individuals into stratas using calculated sampling probabiltiies

  # A.50 Exposure Frequency
  set.seed(3422)
  data$id <- as.numeric(rownames(data))
  pool.exp <- data$id[data$A.50==1]
  pool.unexp <- data$id[data$A.50==0]
  allobs.exp <- c()
  allobs.unexp <- c()
  strata.bias.A.50.exp <- data.frame(pool.exp,rep(NA,length(pool.exp)))
    names(strata.bias.A.50.exp) <- c("id","strata.bias.A.50")
  strata.bias.A.50.unexp <- data.frame(pool.unexp,rep(NA,length(pool.unexp)))
    names(strata.bias.A.50.unexp) <- c("id","strata.bias.A.50")
  
  for(i in strata.data[,1]){
    
    # Updating Sampling Pool to exclude individuals already sampled
    pool.exp <- (pool.exp)[!((pool.exp) %in% allobs.exp)]
    pool.unexp <- (pool.unexp)[!((pool.unexp) %in% allobs.unexp)]
    
    # Sample
    if((clpop.exp[i,1]>length(pool.exp)) | (clpop.exp[i,1]<length(pool.exp) & i == nrow(strata.data))){
      clpop.exp[i,1] <- length(pool.exp)
    }
    if((clpop.unexp[i,1]>length(pool.unexp)) | (clpop.unexp[i,1]<length(pool.unexp) & i == nrow(strata.data))){
      clpop.unexp[i,1] <- length(pool.unexp)
    }
   
    if((clpop.exp[i,1])==1){ 
      strata.obs.exp <- pool.exp
    } else {
      strata.obs.exp <- sample(pool.exp,size=clpop.exp[i,1])
    }
    
    if((clpop.unexp[i,1])==1){ 
      strata.obs.unexp <- pool.unexp
    } else {
      strata.obs.unexp <- sample(pool.unexp,size=clpop.unexp[i,1])
    }
    
    # Assign strata
    strata.bias.A.50.exp[strata.bias.A.50.exp$id %in% strata.obs.exp,2] <- i
    strata.bias.A.50.unexp[strata.bias.A.50.unexp$id %in% strata.obs.unexp,2] <- i
    
    # Finish
    allobs.exp <- c(allobs.exp,strata.obs.exp)
    allobs.unexp <- c(allobs.unexp,strata.obs.unexp)
    print(paste0("strata ",i," Complete"))
    print(paste0("    Exposed Pool: ",length(pool.exp)))
    print(paste0("    Unexposed Pool: ",length(pool.unexp)))
  }  
   
  # Append Exposed and Unexposed strata assignments
  strata.bias.A.50 <- rbind(strata.bias.A.50.exp,strata.bias.A.50.unexp)    
     
  rm(allobs.exp,allobs.unexp,pool.exp,pool.unexp,strata.obs.exp,strata.obs.unexp,strata.bias.A.50.exp,strata.bias.A.50.unexp,i)
  
  # A.20 Exposure Frequency
  set.seed(3422)
  data$id <- as.numeric(rownames(data))
  pool.exp <- data$id[data$A.20==1]
  pool.unexp <- data$id[data$A.20==0]
  allobs.exp <- c()
  allobs.unexp <- c()
  strata.bias.A.20.exp <- data.frame(pool.exp,rep(NA,length(pool.exp)))
  names(strata.bias.A.20.exp) <- c("id","strata.bias.A.20")
  strata.bias.A.20.unexp <- data.frame(pool.unexp,rep(NA,length(pool.unexp)))
  names(strata.bias.A.20.unexp) <- c("id","strata.bias.A.20")
  
  for(i in strata.data[,1]){
    
    # Updating Sampling Pool to exclude individuals already sampled
    pool.exp <- (pool.exp)[!((pool.exp) %in% allobs.exp)]
    pool.unexp <- (pool.unexp)[!((pool.unexp) %in% allobs.unexp)]
    
    # Sample
    if((clpop.exp[i,2]>length(pool.exp)) | (clpop.exp[i,2]<length(pool.exp) & i == nrow(strata.data))){
      clpop.exp[i,2] <- length(pool.exp)
    }
    
    if((clpop.unexp[i,2]>length(pool.unexp)) | (clpop.unexp[i,2]<length(pool.unexp) & i == nrow(strata.data))){
      clpop.unexp[i,2] <- length(pool.unexp)
    }
    
    if((clpop.exp[i,2])==1){ 
      strata.obs.exp <- pool.exp
    } else {
      strata.obs.exp <- sample(pool.exp,size=clpop.exp[i,2])
    }
    
    if((clpop.unexp[i,2])==1){ 
      strata.obs.unexp <- pool.unexp
    } else {
      strata.obs.unexp <- sample(pool.unexp,size=clpop.unexp[i,2])
    }
    
    # Assign strata
    strata.bias.A.20.exp[strata.bias.A.20.exp$id %in% strata.obs.exp,2] <- i
    strata.bias.A.20.unexp[strata.bias.A.20.unexp$id %in% strata.obs.unexp,2] <- i
    
    # Finish
    allobs.exp <- c(allobs.exp,strata.obs.exp)
    allobs.unexp <- c(allobs.unexp,strata.obs.unexp)
    print(paste0("strata ",i," Complete"))
    print(paste0("    Exposed Pool: ",length(pool.exp)))
    print(paste0("    Unexposed Pool: ",length(pool.unexp)))
  }  
  
  # Append Exposed and Unexposed strata assignments
  strata.bias.A.20 <- rbind(strata.bias.A.20.exp,strata.bias.A.20.unexp)    
  
  # Remove data
  rm(allobs.exp,allobs.unexp,pool.exp,pool.unexp,strata.obs.exp,strata.obs.unexp,strata.bias.A.20.exp,strata.bias.A.20.unexp,i)

  # A.10 Exposure Frequency
  set.seed(3422)
  data$id <- as.numeric(rownames(data))
  pool.exp <- data$id[data$A.10==1]
  pool.unexp <- data$id[data$A.10==0]
  allobs.exp <- c()
  allobs.unexp <- c()
  strata.bias.A.10.exp <- data.frame(pool.exp,rep(NA,length(pool.exp)))
  names(strata.bias.A.10.exp) <- c("id","strata.bias.A.10")
  strata.bias.A.10.unexp <- data.frame(pool.unexp,rep(NA,length(pool.unexp)))
  names(strata.bias.A.10.unexp) <- c("id","strata.bias.A.10")
  
  for(i in strata.data[,1]){
    
    # Updating Sampling Pool to exclude individuals already sampled
    pool.exp <- (pool.exp)[!((pool.exp) %in% allobs.exp)]
    pool.unexp <- (pool.unexp)[!((pool.unexp) %in% allobs.unexp)]
    
    # Sample
    if((clpop.exp[i,3]>length(pool.exp)) | (clpop.exp[i,3]<length(pool.exp) & i == nrow(strata.data))){
      clpop.exp[i,3] <- length(pool.exp)
    }
    if((clpop.unexp[i,3]>length(pool.unexp)) | (clpop.unexp[i,3]<length(pool.unexp) & i == nrow(strata.data))){
      clpop.unexp[i,3] <- length(pool.unexp)
    }

    if((clpop.exp[i,3])==1){ 
      strata.obs.exp <- pool.exp
    } else {
      strata.obs.exp <- sample(pool.exp,size=clpop.exp[i,3])
    }
    
    if((clpop.unexp[i,3])==1){ 
      strata.obs.unexp <- pool.unexp
    } else {
      strata.obs.unexp <- sample(pool.unexp,size=clpop.unexp[i,3])
    }
    
    # Assign strata
    strata.bias.A.10.exp[strata.bias.A.10.exp$id %in% strata.obs.exp,2] <- i
    strata.bias.A.10.unexp[strata.bias.A.10.unexp$id %in% strata.obs.unexp,2] <- i
    
    # Finish
    allobs.exp <- c(allobs.exp,strata.obs.exp)
    allobs.unexp <- c(allobs.unexp,strata.obs.unexp)
    print(paste0("strata ",i," Complete"))
    print(paste0("    Exposed Pool: ",length(pool.exp)))
    print(paste0("    Unexposed Pool: ",length(pool.unexp)))
  }  
  
  # Append Exposed and Unexposed strata assignments
  strata.bias.A.10 <- rbind(strata.bias.A.10.exp,strata.bias.A.10.unexp)    
  
  # Remove data
  rm(allobs.exp,allobs.unexp,pool.exp,pool.unexp,strata.obs.exp,strata.obs.unexp,strata.bias.A.10.exp,strata.bias.A.10.unexp,i)
  
  # A.05 Exposure Frequency
  set.seed(3422)
  data$id <- as.numeric(rownames(data))
  pool.exp <- data$id[data$A.05==1]
  pool.unexp <- data$id[data$A.05==0]
  allobs.exp <- c()
  allobs.unexp <- c()
  strata.bias.A.05.exp <- data.frame(pool.exp,rep(NA,length(pool.exp)))
    names(strata.bias.A.05.exp) <- c("id","strata.bias.A.05")
  strata.bias.A.05.unexp <- data.frame(pool.unexp,rep(NA,length(pool.unexp)))
    names(strata.bias.A.05.unexp) <- c("id","strata.bias.A.05")
  
  for(i in strata.data[,1]){
    
    # Updating Sampling Pool to exclude individuals already sampled
    pool.exp <- (pool.exp)[!((pool.exp) %in% allobs.exp)]
    pool.unexp <- (pool.unexp)[!((pool.unexp) %in% allobs.unexp)]
    
    # Sample
    if((clpop.exp[i,4]>length(pool.exp)) | (clpop.exp[i,4]<length(pool.exp) & i == nrow(strata.data))){
      clpop.exp[i,4] <- length(pool.exp)
    }
   
    if((clpop.unexp[i,4]>length(pool.unexp)) | (clpop.unexp[i,4]<length(pool.unexp) & i == nrow(strata.data))){
      clpop.unexp[i,4] <- length(pool.unexp)
    }
    
    if((clpop.exp[i,4])==1){ 
      strata.obs.exp <- pool.exp
    } else {
        strata.obs.exp <- sample(pool.exp,size=clpop.exp[i,4])
    }
    
    if((clpop.unexp[i,4])==1){ 
      strata.obs.unexp <- pool.unexp
    } else {
        strata.obs.unexp <- sample(pool.unexp,size=clpop.unexp[i,4])
    }
    
    # Assign strata
    strata.bias.A.05.exp[strata.bias.A.05.exp$id %in% strata.obs.exp,2] <- i
    strata.bias.A.05.unexp[strata.bias.A.05.unexp$id %in% strata.obs.unexp,2] <- i
    
    # Finish
    allobs.exp <- c(allobs.exp,strata.obs.exp)
    allobs.unexp <- c(allobs.unexp,strata.obs.unexp)
    print(paste0("strata ",i," Complete"))
    print(paste0("    Exposed Pool: ",length(pool.exp)))
    print(paste0("    Unexposed Pool: ",length(pool.unexp)))
  }  
  
  # Append Exposed and Unexposed strata assignments
  strata.bias.A.05 <- rbind(strata.bias.A.05.exp,strata.bias.A.05.unexp)    
  
  # Remove data
  rm(allobs.exp,allobs.unexp,pool.exp,pool.unexp,strata.obs.exp,strata.obs.unexp,strata.bias.A.05.exp,strata.bias.A.05.unexp,i)
  
  
# Merge in strata Assignments
data <- merge(data,strata.bias.A.50,by="id")
data <- merge(data,strata.bias.A.20,by="id")
data <- merge(data,strata.bias.A.10,by="id")
data <- merge(data,strata.bias.A.05,by="id")
data$id <- NULL

# Remove data
rm(strata.bias.A.50,strata.bias.A.20,strata.bias.A.10,strata.bias.A.05,clpop.exp,clpop.unexp,strata.data)

# Check Data

  # Confirm no missing strata assignments
  sum(is.na(data$strata.bias.A.50))
  sum(is.na(data$strata.bias.A.20))
  sum(is.na(data$strata.bias.A.10))
  sum(is.na(data$strata.bias.A.05))

  # Confirm increasing strata size
  plot(prop.table(table(data$strata.bias.A.50)), type="l") # Proportion
  plot(table(data$strata.bias.A.50), type="l") # Size
  
  # Confirm descending probability of sampling for exposued as strata size increases
  plot(prop.table(table(data$strata.bias.A.50, data$A.50),margin=1)[,2], type="l")
  plot(prop.table(table(data$strata.bias.A.20, data$A.20),margin=1)[,2], type="l")
  plot(prop.table(table(data$strata.bias.A.10, data$A.10),margin=1)[,2], type="l")
  plot(prop.table(table(data$strata.bias.A.05, data$A.05),margin=1)[,2], type="l")

## Save population data
write.csv(data, "data/population.data.append.csv", row.names=F)
