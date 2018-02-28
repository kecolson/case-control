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
  truth <- ests$truth[1]
  
  # Summarize distribution of point estimates and CIs
  hist(ests$est, main="Distribution of Point Estimates", xlab="Point Estimate", breaks=30)
  
  # Calculate 95% CI coverage - % of calculated CIs that include the true OR
  CIcover <- as.numeric(ests$lower<=truth & ests$upper>=truth)
  print(paste0(round(mean(CIcover, na.rm=T)*100,5),"% of confidence intervals include the true measure of association."))
  
  # Calculate mean estimate
  avg <- round(mean(ests$est))
  print(paste0("Mean estimate = "),avg)
  
  # Calculate Bias - average distance of point estimates away from true OR in repeated simulations
  bias <- round(mean(ests$est - truth),10) 
  print(paste0("Mean bias = ",bias))
  
  # Calculate Relative Bias - average distance of point estimates away from true OR in repeated simulations
  relbias <- round(mean((ests$est - truth)/truth*100),10)
  print(paste0("Mean relative bias (%) = ",relbias))
  
  # Calculate Variance - variance of point estimate of repeated simulations
  variance <- round(var(ests$est),10)
  print(paste0("Variance = ",variance))
  
  # Calculate MSE - mean squared error of point estimate of repeated simulations
  MSE <- round(mean((ests$est - truth)^2),10)
  print(paste0("MSE = ",MSE))
}

# END