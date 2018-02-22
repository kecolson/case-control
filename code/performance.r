################################################################################################
# NAME: performance.r
# AUTHORS: Ellie Matthay, Catherine Li, Chris Rowe
# DATE STARTED: 2/21/2018
# PURPOSE: Script to define a function to calculate performance metrics for the simulations
# UPDATES: [date]: XX
################################################################################################

performance <- function(sim) {
  print(paste0("RESULTS FROM ", nrow(sim$est.lower.upper)," SIMULATIONS:"))
  
  # Format the point estimates and CIs
  ests <- as.data.frame(sim$est.lower.upper)
  for (vv in names(ests)) ests[[vv]] <- as.numeric(ests[[vv]])
  
  # Pull the truth
  truth <- sim$results[[1]]$truth
  
  # Summarize distribution of point estimates and CIs
  hist(ests$est, main="Distribution of Point Estimates", xlab="Point Estimate", breaks=30)
  
  # Calculate 95% CI coverage - % of calculated CIs that include the true OR
  CIcover <- as.numeric(ests$lower<=truth & ests$upper>=truth)
  print(paste0(round(mean(CIcover, na.rm=T)*100,1),"% of confidence intervals include the true measure of association."))
  
  # Calculate Bias - average distance of point estimates away from true OR in repeated simulations
  bias <- round(mean(ests$est - truth),4) 
  print(paste0("Mean bias = ",bias))
  
  # Calculate Variance - variance of point estimate of repeated simulations
  variance <- round(var(ests$est),4)
  print(paste0("Variance = ",variance))
  
  # Calculate MSE - mean squared error of point estimate of repeated simulations
  MSE <- round(mean((ests$est - truth)^2),4)
  print(paste0("MSE = ",MSE))
}

# END