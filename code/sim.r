################################################################################################
# NAME: sim.r
# AUTHORS: Ellie Matthay, Catherine Li, Chris Rowe
# DATE STARTED: 2/21/2018
# PURPOSE: Define a function that will parallelize and run the different case-control simulations
# UPDATES: 4/24/2018: CR add svysize and method study.r arguments throughout.
#          11/7/2018: Added svycase argument to funciton, allowing choice of including/excluding
#                     cases from survey; capture SE as part of results
################################################################################################


sim <- function(nsims, # Number of simulations to run. Probably 3 for testing, 500-1000 for initial, 2000 for final results
                cluster, # Set to TRUE to run on grizzlybear; set to FALSE to run locally
                cctype, samp, svysize, svycase, method, ratio, exposure, outcome, timevar) {
  
  # For testing: nsims <- 3; cluster <- F; cctype <- "cumulative"; samp <- "exp.ps"; svysize <- "small"; svycase <- FALSE; method <- "expand"; ratio <- 1; exposure <- "A.50"; outcome <- "Y.02.A.50"; timevar <- "time.Y.02.A.50"
  
  library("parallel") # For setting random seeds
  library("Epi") # for density sampling design
  library("survival") # for clogit analysis
  library("dplyr") # for data management
  
  index <- 1:nsims
  
  # Set Working Directory
  if (cluster==F) {
    setwd("~/Documents/GitHub/case-control")# Chris's directory
    #setwd("C:/Users/kecolson/Google Drive/simulation/case-control-other") # Ellie's directory
    #setwd("C:/Users/Catherine/Documents/GitHub/case-control_master") # Catherine's directory
  }

  # Bring in data and true parameters
  data <- read.csv("data/population.data.csv", stringsAsFactors = F)
  
  # Save memory by keeping only variables of interest
  data <- data[,c('cluster','strata','serial','county','city','puma','cpuma0010',
		  'white','black','asian','hispanic','otherrace','male',
		  'age_18_24','age_25_34','age_35_44','age_45_54','age_55_64','age_over64',
		  'educ_lesshs','educ_ged','educ_hs','educ_somecollege','educ_associates',
		  'educ_bachelors','educ_advdegree',
		  exposure, outcome, timevar, paste0("trueOR.",outcome),paste0('trueIDR.',outcome) )]
  
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
    mpi.bcast.Robj2slave(svysize)
    mpi.bcast.Robj2slave(svycase)
    mpi.bcast.Robj2slave(method)
    mpi.bcast.Robj2slave(ratio)
    mpi.bcast.Robj2slave(data)
    mpi.bcast.Robj2slave(exposure)
    mpi.bcast.Robj2slave(outcome)
    mpi.bcast.Robj2slave(timevar)
    mpi.bcast.Robj2slave(seeds)
    
    print("sim.r: running sims...")
    
    # Send the work to the worker nodes and run the simulations
    results <- mpi.parLapply(index, study, cctype=cctype, samp=samp, svysize=svysize, svycase=svycase, method=method, ratio=ratio, data=data, 
                             exposure=exposure, outcome=outcome, timevar=timevar, seeds=seeds)
    
    print("sim.r: sims completed. closing workers...")
    
    try(mpi.close.Rslaves())

    print("workers closed. Collapsing results....")
    
  } else { # if not on cluster, run these simulations locally
    results <- lapply(index, study, cctype=cctype, samp=samp, svysize=svysize, svycase=svycase, method=method, ratio=ratio, data=data, 
                      exposure=exposure, outcome=outcome, timevar=timevar, seeds=seeds) 
  }
  
  # What comes out of study(): list(est=est, lower=lower, upper=upper, se=se, # the point estimate and CI and SE
  #                                 truth=truth) # the true OR or IDR as relevant
  
  # Print some things, to help debug
  print("length of results object to be collapsed:")
  print(length(results))

  print("summary of results object to be collapsed:")
  print(summary(results))

  print("head of first object to be collapsed:")
  print(head(results[[1]]))
    
  # Collapse est, lower, upper, and se results; leave samples and model objects in list form in "results" object
  est.lower.upper <- do.call(rbind, lapply(index, function(x) results[[x]][c('est','lower','upper','se','truth')] ))
  
  print("sim.r: Results collapsed. Returning results...")
  
  return(list(est.lower.upper=est.lower.upper))
}



# END