################################################################################################
# NAME: sim.r
# AUTHORS: Ellie Matthay, Catherine Li, Chris Rowe
# DATE STARTED: 2/21/2018
# PURPOSE: Define a function that will parallelize and run the different case-control simulations
# UPDATES: [date]: XX
################################################################################################


sim <- function(nsims, # Number of simulations to run. Probably 3 for testing, 500-1000 for initial, 2000 for final results
                cluster, # Set to TRUE to run on grizzlybear; set to FALSE to run locally
                cctype, samp, ratio, exposure, outcome, timevar) {
  
  # For testing: nsims <- 3; cluster <- F; cctype <- "cumulative"; samp <- "srs"; ratio <- 1; exposure <- "A.5"; outcome <- "Y.02.A.5"; timevar <- "time"
  
  library("parallel") # For setting random seeds
  library("Epi") # for density sampling design
  library("survival") # for clogit analysis
  library("dplyr") # for data management
  
  index <- 1:nsims
  
  # Set Working Directory
  if (cluster==F) {
    #setwd("~/Documents/PhD/Ahern GSR/Case Control Simulation") # Chris's directory
    #setwd("C:/Users/kecolson/Google Drive/simulation/case-control-other") # Ellie's directory
    setwd("C:/Users/Catherine/Documents/GitHub/case-control_master") # Catherine's directory
  }

  # Bring in data and true parameters
  data <- read.csv("data/population.data.csv", stringsAsFactors = F)
  
  # Sample down the dataset for testing purposes
  #data <- data[sample(1:nrow(data), 4000, replace=F),]

  # Create #nsims# random seeds that will be really independent and not likely to loop 
  # back through the same random number generator. This is good practice for doing 
  # simulation work.
  RNGkind("L'Ecuyer-CMRG")
  set.seed(1)
  s <- .Random.seed
  seeds <- list(s)  
  if (nsims>=2) {
    for (i in 2:nsims) {
      s <- nextRNGStream(s)
      seeds[[i]] <- s
    }
  }
  
  if (cluster==T) {
    # Set up the worker nodes
    library("Rmpi")
    mpi.spawn.Rslaves() # # Spawn as many worker nodes as possible
    
    # In case R exits unexpectedly, have it automatically clean up
    # resources taken up by Rmpi (slaves, memory, etc...)
    .Last <- function(){
      if (is.loaded("mpi_initialize")){
        if (mpi.comm.size(1) > 0){
          print("Please use mpi.close.Rslaves() to close slaves.")
          mpi.close.Rslaves()
        }
        print("Please use mpi.quit() to quit R")
        .Call("mpi_finalize")
      }
    }
    
    # Tell all worker nodes to return a message identifying themselves
    mpi.remote.exec(paste("I am",mpi.comm.rank(),"of",mpi.comm.size()))
    
    # Send all necessary information and functions to the worker nodes
    mpi.bcast.Robj2slave(index)
    mpi.bcast.Robj2slave(study)
    mpi.bcast.Robj2slave(cctype)
    mpi.bcast.Robj2slave(samp)
    mpi.bcast.Robj2slave(ratio)
    mpi.bcast.Robj2slave(data)
    mpi.bcast.Robj2slave(exposure)
    mpi.bcast.Robj2slave(outcome)
    mpi.bcast.Robj2slave(timevar)
    mpi.bcast.Robj2slave(seeds)
    
    
    # Send the work to the worker nodes and run the simulations
    results <- mpi.parLapply(index, study, cctype=cctype, samp=samp, ratio=ratio, data=data, 
                             exposure=exposure, outcome=outcome, timevar=timevar, seeds=seeds) 
    
    # Turn off worker nodes
    mpi.close.Rslaves()
    
  } else { # if not on cluster, run these simulations locally
    results <- lapply(index, study, cctype=cctype, samp=samp, ratio=ratio, data=data, 
                      exposure=exposure, outcome=outcome, timevar=timevar, seeds=seeds) 
  }
  
  # What comes out of study(): list(sample=sample, # the sampled data
  #                                 mod=mod, # the model object
  #                                 est=est, lower=lower, upper=upper, # the point estimate and CI
  #                                 truth=truth) # the true OR or IDR as relevant
    
  # Collapse est, lower, and upper results; leave samples and model objects in list form in "results" object
  est.lower.upper <- do.call(rbind, lapply(index, function(x) results[[x]][c('est','lower','upper')] ))
  
  return(list(results=results, est.lower.upper=est.lower.upper))
}



# END