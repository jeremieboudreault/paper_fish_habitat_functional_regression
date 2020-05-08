####################################################################################
#     Modelling fish habitat selection using functional regression 🐟
#     Created by :    Jeremie Boudreault (Jeremie.Boudreault@ete.inrs.ca)
#     Current file :  Models_GAMs_Predictions.R
#     Aim :           Calculate the leave-one-out predictions with the best GAMs
####################################################################################

#### Part 0 : Libraries, functions and variable definition

# Import library
library('tidyr') # To pivot a data.frame
library('dplyr') # Data.frame manipulation
library('ggplot2') # Data visualisation
library('mgcv') # For fitting Gams

# Import data per site of the calibration river (SMR)
SMR.per.site <- readRDS('data/SMR_2017_clean_per_site.Rda')

# Calculate the number of observations
nObs <- nrow(SMR.per.site)

# Import GAMs models
GAM.fry <- readRDS('out/models/GAM_nFry_M.rda')
GAM.parr1 <- readRDS('out/models/GAM_nParr1_M.rda')

#### Part 1 : Models fitting with different scenarios

# For each fish
for (Y in c('nFry_M', 'nParr1_M'))
{
  
  # Extract the formula and the family of the GAM
  if (Y == 'nFry_M')  {
    Formula.gam <- GAM.fry$formula
    Family.gam <- nb(link = 'log')
  } else if (Y == 'nParr1_M') {
    Formula.gam <- GAM.parr1$formula
    Family.gam <- nb(link = 'log')
  }

  # Generate a vector to save the leave-one-out prediction
  Y.loo <- rep(NA, nObs)

  # Loop on every observations
  for (Obs.i in 1:nObs) {
    
    # Fit the GAM model
    model.gam.loo <- gam(formula=Formula.gam, data = SMR.per.site[-Obs.i, ], family=Family.gam)
    
    # Prediction of the model
    Y.loo[Obs.i]  <- predict(model.gam.loo, newdata=SMR.per.site[Obs.i, ], type='response')
    
    # End of loop on each prediction
  }
  
  # Calculate an error measure on the Leave-one-out predictions
  saveRDS(Y.loo, paste0('out/predictions/YHatLoo_GAM_', Y, '.rda'))
  
}

