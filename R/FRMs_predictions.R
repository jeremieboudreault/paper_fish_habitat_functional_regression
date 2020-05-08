####################################################################################
#     Modelling fish habitat selection using functional regression 🐟
#     Created by :    Jeremie Boudreault (Jeremie.Boudreault@ete.inrs.ca)
#     Current file :  Predictions_FRM.R
#     Aim :           Calculate the leave-one-out predictions with the FRMs
####################################################################################

#### Part 0 : Libraries, functions and variable definition

# Import library
library('dplyr') # Data.frame manipulation
library('ggplot2') # Data visualisation
library('FDboost') # For fitting Frms

# Import data per site of the calibration river (SMR)
SMR.per.site <- readRDS('data/SMR_2017_clean_per_site.Rda')

# Import the functional observations
SMR.funobs <- readRDS('data/SMR_2017_clean_fun_obs.Rda')

# Count the number of observations
nObs <- nrow(SMR.per.site)

# Import the models
FRM.fry <- readRDS('out/models/FRM_nFry_M.rda')
FRM.parr1 <- readRDS('out/models/FRM_nParr1_M.rda')

# Function to eval any string character
evalStr <- function(...)
  eval(parse(text=paste0(..., collapse = "", sep="")))

# Function to rbind and create the vector if it doesn't exist
rbind2 <- function(matname, vect, index=2, name="") {
  if (index == 1) {
    suppressWarnings(remove(list = matname, envir = .GlobalEnv))
    assign(matname, vect, envir = .GlobalEnv)
  } else {
    assign(matname, rbind(evalStr(matname), vect), envir= .GlobalEnv)
  }
  evalStr(matname)
}

#### Part 1 : Leave-one-out predictions for Fry

# The Y variable
Y <- 'nFry_M'
  
# Parameter of the model
Family.frm <- FRM.fry$family
mstop <- FRM.fry$control$mstop
nu <-  FRM.fry$control$nu
Formula.frm <- as.formula(FRM.fry$formulaFDboost)

# Saving the Y.loo
Y.loo <- rep(NA, nObs)
  
# Loop on every observations
for (Obs.i in 1:nObs) {
  
  # If we want to do only one observations
  # Obs.i <- 1
  
  # Create the list of the Y variable
  Subset <- list(Y=as.vector(unlist(SMR.per.site[, Y]))[-Obs.i],
                 Velocity = t(as.matrix(SMR.funobs[which(SMR.funobs$Var=='Velocity'), 4:29]))[-Obs.i, ],
                 Velocity.x = SMR.funobs[which(SMR.funobs$Var=='Velocity'), 2],
                 D50 = t(as.matrix(SMR.funobs[which(SMR.funobs$Var=='D50'), 4:29]))[-Obs.i, ],
                 D50.x = SMR.funobs[which(SMR.funobs$Var=='D50'), 2],
                 Temp = t(as.matrix(SMR.funobs[which(SMR.funobs$Var=='Temp'), 4:29]))[-Obs.i, ],
                 Temp.x = SMR.funobs[which(SMR.funobs$Var=='Temp'), 2])
  
  # Adjust the functional regression model
  model.frm.loo <- FDboost(Formula.frm ,  timeformula=NULL, family=Family.frm, data=Subset,  control=boost_control(mstop = mstop, nu = nu))
  
  # Create a new subset for the prediction
  Subset_new <- list(Velocity = rbind(t(as.matrix(SMR.funobs[which(SMR.funobs$Var=='Velocity'), 4:29]))[Obs.i, ]),
                     Velocity.x = SMR.funobs[which(SMR.funobs$Var=='Velocity'), 2],
                     D50 = rbind(t(as.matrix(SMR.funobs[which(SMR.funobs$Var=='D50'), 4:29]))[Obs.i, ]),
                     D50.x = SMR.funobs[which(SMR.funobs$Var=='D50'), 2],
                     Temp = rbind(t(as.matrix(SMR.funobs[which(SMR.funobs$Var=='Temp'), 4:29]))[Obs.i, ]),
                     Temp.x = SMR.funobs[which(SMR.funobs$Var=='Temp'), 2])
  
  # Predictions versus observation
  Y.loo[Obs.i] <- predict(model.frm.loo, newdata=Subset_new , type='response')
  
  # Extract the coefficient
  modelcoef <- coef(model.frm.loo)
  
  # On enregistre les coefficients
  for (beta.i in 1:3) {
    rbind2(paste0("Beta", beta.i), t(modelcoef$smterms[[beta.i]]$value), Obs.i)
    assign(paste0("Beta", beta.i, ".x"),  modelcoef$smterms[[beta.i]]$x)
  }
  
  # End of loop on each prediction
}

# We are saving the results
saveRDS(Y.loo, paste0('out/predictions/YHatLoo_FRM_', Y, '.rda'))

# We rearrange the coefficient to save them
Coef <- data.frame(var='Velocity', x=Beta1.x, t(Beta1))
Coef <- rbind(Coef, data.frame(var='D50', x=Beta2.x, t(Beta2)))
Coef <- rbind(Coef, data.frame(var='Temp', x=Beta3.x, t(Beta3)))

# Save the leave-one-out coefficients
saveRDS(Coef, paste0('out/coefficients/Coef_LOO_FRM_', Y, '.rda'))

#### Part 2 : Leave-one-out prediction for Parr 1+

# The Y variable
Y <- 'nParr1_M'

# Parameter of the model
Family.frm <- FRM.parr1$family
mstop <- FRM.parr1$control$mstop
nu <-  FRM.parr1$control$nu
Formula.frm <- as.formula(FRM.parr1$formulaFDboost)

# Saving the Y.loo
Y.loo <- rep(NA, nObs)

# Loop on every observations
for (Obs.i in 1:nObs) {
  
  # If we want to do only one observations
  # Obs.i <- 1
  
  # Create the list of the Y variable
  Subset <- list(Y=as.vector(unlist(SMR.per.site[, Y]))[-Obs.i],
                 Depth = t(as.matrix(SMR.funobs[which(SMR.funobs$Var=='Depth'), 4:29]))[-Obs.i, ],
                 Depth.x = SMR.funobs[which(SMR.funobs$Var=='Depth'), 2],
                 D50 = t(as.matrix(SMR.funobs[which(SMR.funobs$Var=='D50'), 4:29]))[-Obs.i, ],
                 D50.x = SMR.funobs[which(SMR.funobs$Var=='D50'), 2],
                 Temp = t(as.matrix(SMR.funobs[which(SMR.funobs$Var=='Temp'), 4:29]))[-Obs.i, ],
                 Temp.x = SMR.funobs[which(SMR.funobs$Var=='Temp'), 2])
  
  # Adjust the functional regression model
  model.frm.loo <- FDboost(Formula.frm ,  timeformula=NULL, family=Family.frm, data=Subset,  control=boost_control(mstop = mstop, nu = nu))
  
  # Create a new subset for the prediction
  Subset_new <- list(Depth = rbind(t(as.matrix(SMR.funobs[which(SMR.funobs$Var=='Depth'), 4:29]))[Obs.i, ]),
                     Depth.x = SMR.funobs[which(SMR.funobs$Var=='Depth'), 2],
                     D50 = rbind(t(as.matrix(SMR.funobs[which(SMR.funobs$Var=='D50'), 4:29]))[Obs.i, ]),
                     D50.x = SMR.funobs[which(SMR.funobs$Var=='D50'), 2],
                     Temp = rbind(t(as.matrix(SMR.funobs[which(SMR.funobs$Var=='Temp'), 4:29]))[Obs.i, ]),
                     Temp.x = SMR.funobs[which(SMR.funobs$Var=='Temp'), 2])
  
  # Predictions versus observation
  Y.loo[Obs.i] <- predict(model.frm.loo, newdata=Subset_new , type='response')
  
  # Extract the coefficient
  modelcoef <- coef(model.frm.loo)
  
  # On enregistre les coefficients
  for (beta.i in 1:3) {
    rbind2(paste0("Beta", beta.i), t(modelcoef$smterms[[beta.i]]$value), Obs.i)
    assign(paste0("Beta", beta.i, ".x"),  modelcoef$smterms[[beta.i]]$x)
  }
  
  # End of loop on each prediction
}

# We are saving the results
saveRDS(Y.loo, paste0('out/predictions/YHatLoo_FRM_', Y, '.rda'))

# We rearrange the coefficient to save them
Coef <- data.frame(var='Depth', x=Beta1.x, t(Beta1))
Coef <- rbind(Coef, data.frame(var='D50', x=Beta2.x, t(Beta2)))
Coef <- rbind(Coef, data.frame(var='Temp', x=Beta3.x, t(Beta3)))

# Save the leave-one-out coefficients
saveRDS(Coef, paste0('out/coefficients/Coef_LOO_FRM_', Y, '.rda'))
