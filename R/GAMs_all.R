####################################################################################
#     Modelling fish habitat selection using functional regression 🐟
#     Created by :    Jeremie Boudreault (Jeremie.Boudreault@ete.inrs.ca)
#     Current file :  Models_GAM_all.R
#     Aim :           Model fish abundance with GAM models
####################################################################################

#### Part 0 : Libraries, functions and variable definition

# Import library
library('dplyr') # Data.frame manipulation
library('ggplot2') # Data visualisation
library('mgcv') # For fitting Gams

# Load goodness-of-fit functions
source('R/functions/Goodness_of_fit.R')

# Import data per site (SMR)
SMR.per.site <- readRDS('data/SMR_2017_clean_per_site.Rda')

# Calculate the number of observations
nObs <- nrow(SMR.per.site)

# Ensure we obtain the same resuls
set.seed('2912')

#### Part 1 : Models fitting with different scenarios

# Create a data.frame of different scenarios (0 = GLM, 3:6 = degrees of freedom for GAM)
Scenarios <- expand.grid(k = c(0, c(3:6)), family = c('NBinomial'), Y = names(SMR.per.site[, c('nFry_M', 'nParr1_M')]))

# Add other informations we want to save
Scenarios$Dev <- NA
Scenarios$DevExpl <- NA
Scenarios$Rsq <- NA
Scenarios$Rsqadj <- NA
Scenarios$RMSELoo <- NA
Scenarios$x0 <- NA
Scenarios$x1 <- NA
Scenarios$x2 <- NA
Scenarios$x3 <- NA
Scenarios$x4 <- NA

# Number of scenarios
nScenarios <- nrow(Scenarios)

# We run each scenarios
for (Scenario.i in 1:nScenarios) {

  # If we want to focus on only one scenario
  # Scenario.i <- 3
  
  # Take the information of the current scenario
  Scenario <- Scenarios[Scenario.i, ]
  
  # Degrees of freedom of the smooth effect
  k <- Scenario$k
  
  # Family of the response variable
  if (Scenario$family == 'Normal')  {
    Family.gam <- gaussian(link = "identity")
  } else if (Scenario$family == 'Poisson') {
    Family.gam <- poisson(link = "log")
  } else if (Scenario$family == 'NBinomial') {
    Family.gam <- nb(link = 'log')
  }
  
  # Response variable (many will be tested)
  Y <- as.character(Scenario$Y)
  YObs <- as.vector(unlist(SMR.per.site[, Y]))
  
  # Create the formula for the GAM
  if (k != 0) {
    Formula.gam <- as.formula(paste0(Y, ' ~ s(Depth, k=', k, ') + s(Velocity, k=', k, ') + s(D50, k=', k, ') + s(Temp, k=', k, ')'))
  } else {
    Formula.gam <- as.formula(paste0(Y, ' ~ Depth + Velocity + D50 + Temp'))
  }
  
  # Generate a vector to save the leave-one-out prediction
  YHatLoo <- rep(NA, nObs)
  
  # Loop on every observations for the leave-one-out procedure
  for (Obs.i in 1:nObs) {
    
    # Fit the GAM model
    model.gam.loo <- gam(formula=Formula.gam, 
                         data = SMR.per.site[-Obs.i, ], 
                         family=Family.gam,
                         method = 'GCV.Cp')
    
    # Prediction of the model
    YHatLoo[Obs.i]  <- predict(model.gam.loo, newdata=SMR.per.site[Obs.i, ], type='response')
    
    # End of loop on each prediction
  }
  
  # Save these predictions
  saveRDS(YHatLoo, paste0('out/models/all/tmp/GAM/YHatLoo_', Y, '_k_', k, '_family_',Scenario$family, ".Rds"))
  
  # Fit the full GAM model
  model.gam <- gam(formula=Formula.gam, 
                   data = SMR.per.site, 
                   family=Family.gam, 
                   method = 'GCV.Cp')
  
  # Save this model
  saveRDS(model.gam, paste0('out/models/all/tmp/GAM/Model_', Y, '_k_', k, '_family_',Scenario$family, ".Rds"))
  
  # Summary of the model
  (a <- summary(model.gam))
  
  # Extract the p-value of the coefficient
  Scenarios[Scenario.i, paste0('x', 0:4)] <-  c(a$p.pv, a$s.pv)
  
  # Prediction of the model
  YHat <- predict(model.gam, type='response')
  
  # Calculate the deviance of the FRM
  Dev.res <- CalcDevianceNegBin(model.gam, YObs)
  
  # Save some Goodness-of-fit criteria
  Scenarios[Scenario.i, 'Dev'] <-  Dev.res$ModelDeviance
  Scenarios[Scenario.i, 'DevExpl'] <-  Dev.res$ExplDeviance
  Scenarios[Scenario.i, 'Rsq'] <- CalcRsquare(YObs, YHat)
  Scenarios[Scenario.i, 'Rsqadj'] <- CalcRsquareAdj(YObs, YHat, p = 4)
  Scenarios[Scenario.i, 'RMSELoo'] <- CalcMSE(YObs, YHatLoo, root = T)
  
  # Plotting the coefficient
  # plot(model.gam)
  
}

# Save results to .csv file
write.csv2(Scenarios, 'out/models/all/GAMs_results.csv')
