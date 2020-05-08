####################################################################################
#     Modelling fish habitat selection using functional regression 🐟
#     Created by :   Jeremie Boudreault (Jeremie.Boudreault@ete.inrs.ca)
#     Current file : FRMs_all.R
#     Aim :          Model fish abundance with FRM models
####################################################################################

#### Part 0 : Libraries, functions and variable definition

# Import library
library('dplyr') # Data.frame manipulation
library('ggplot2') # Data visualisation
library('FDboost') # For fitting Frms

# Load goodness-of-fit functions
source('R/functions/Goodness_of_fit.R')

# Load the leave-one-out fitting of the FDboost
source('R/functions/FDboost_loo.R')

# Import data per site for the SMR
SMR.per.site <- readRDS('data/SMR_2017_clean_per_site.Rda')

# Import the functional observations
SMR.funobs <- readRDS('data/SMR_2017_clean_fun_obs.Rda')

# Count the number of observations
nObs <- nrow(SMR.per.site)

# Ensure we obtain the same resuls
set.seed('2912')

# Set a variable to stop the iteration
mstop_max <- 1000

#### Part 1 : Models fitting with different scenarios

# Create a data.frame of different scenarios
Scenarios <- expand.grid(df = 4, family = c('nbinomial'), Y = names(SMR.per.site[, c('nFry_M', 'nParr1_M')]))

# Add other informations we want to save
Scenarios$mstop <- NA
Scenarios$nu <- NA
Scenarios$Dev <- NA
Scenarios$DevExpl <- NA
Scenarios$Rsq <- NA
Scenarios$Rsqadj <- NA
Scenarios$RMSELoo <- NA
Scenarios$x1 <- NA
Scenarios$x2 <- NA
Scenarios$x3 <- NA
Scenarios$x4 <- NA

# Number of scenarios
nScenarios <- nrow(Scenarios)

# We run each scenarios
for (Scenario.i in 1:nScenarios) {
  
  # Specific scenario to run
  # Scenario.i <- 1
  
  # Print progression
  print(paste0('Scenario #', Scenario.i, ' (out of n=', nScenarios, ')'))

  # Take the information of the current scenario
  Scenario <- Scenarios[Scenario.i, ]
  
  # The degree of freedom of the functional observation
  df <- Scenario$df
  
  # Family of the response variable
  if (Scenario$family == 'gaussian')  {
    Family.frm <- Gaussian() 
    nu <- 0.1
  } else if (Scenario$family == 'poisson') {
    Family.frm <- Poisson() 
    nu <- 0.01
  } else if (Scenario$family == 'nbinomial') {
    Family.frm <- NBinomial() 
    nu <- 0.01
  }
  
  # The Y variable
  Y <- as.character(Scenario$Y)
  YObs <- as.vector(unlist(SMR.per.site[, Y]))
  
  # The explanatory variables with want to use in our model
  Vars <- c(names(SMR.per.site)[c(12, 13, 14, 16)])
  
  # Create the list of the Y variable for the model with complete observations
  Subset <- list(Y=as.vector(unlist(SMR.per.site[, Y])),
                 Depth = t(as.matrix(SMR.funobs[which(SMR.funobs$Var=='Depth'), 4:29])),
                 Depth.x = SMR.funobs[which(SMR.funobs$Var=='Depth'), 2],
                 Velocity = t(as.matrix(SMR.funobs[which(SMR.funobs$Var=='Velocity'), 4:29])),
                 Velocity.x = SMR.funobs[which(SMR.funobs$Var=='Velocity'), 2],
                 D50 = t(as.matrix(SMR.funobs[which(SMR.funobs$Var=='D50'), 4:29])),
                 D50.x = SMR.funobs[which(SMR.funobs$Var=='D50'), 2],
                 Temp = t(as.matrix(SMR.funobs[which(SMR.funobs$Var=='Temp'), 4:29])),
                 Temp.x = SMR.funobs[which(SMR.funobs$Var=='Temp'), 2])
  
  # The formula for the model
  Formula.frm <- as.formula(paste0("Y~", paste0("bsignal(", Vars, ", ", Vars, ".x", ", df=", df, ", knots=6, degree=3, difference=1)", collapse="+")))
  
  # Generate a matrix of mstop x nobs observations to store the leave-one-out prediction
  LooResultTmp <- FDboost_loo(Formula.frm, Family.frm, data.list = Subset, mstop_max = mstop_max, nu=nu)
    
  # We will save the matrix of leave-one-out estimates for future purpose
  saveRDS(LooResultTmp, paste0('out/models/all/tmp/FRM/YHatLoo_', Y, '_df_', df, '_family_',Scenario$family, ".Rds"))
  
  # We will calculate the MSE on the leave-one-out predictions
  MSE <- colSums((t(LooResultTmp)-YObs)^2)/nObs
  
  # Defining the optimal mstop_best from the MSE
  mstop_best <- which.min(MSE)
  
  # We will save a MVP of the visualisation of the MSE and the number of iteration
  pdf(paste0('out/models/all/tmp/FRM/Boosting_it_', Y, '_df_', df, '_family_', Scenario$family, ".pdf"), width=5, height=4)
  print(ggplot(data=data.frame(mstop=1:mstop_max, MSE=MSE), aes(x=mstop, y=MSE)) + geom_line() + 
     geom_vline(aes(xintercept=mstop_best), lty=3) + ggtitle('Boosting iterations versus MSE'))
  dev.off()

  # Adjust the functional regression model with all covariates
  model.frm <- FDboost(Formula.frm ,  timeformula=NULL, family=Family.frm, data=Subset,  control=boost_control(mstop = mstop_best, nu = nu))
  
  # Save this model 
  saveRDS(model.frm, paste0('out/models/all/tmp/FRM/Model_', Y, '_df_', df, '_family_',Scenario$family, ".Rds"))
  
  # Save the information
  Scenarios[Scenario.i, 'mstop'] <- mstop_best
  Scenarios[Scenario.i, 'nu'] <- nu
  
  # Extract the summary information
  (a <- summary(model.frm))

  # Extract the probability of selection for each baselearners
  Selprob <- a$selprob
  
  # Extract the var names associated with each selected var
  Selprob.vars <- sapply(names(Selprob), function(w) substr(w, 13, regexpr ( ',' , w)-1))
    
  # Save the information about the probability of selection
  Scenarios[Scenario.i, 'x1'] <- ifelse(length(Selprob[which(Selprob.vars=='Depth')]) == 0, 0, Selprob[which(Selprob.vars=='Depth')])
  Scenarios[Scenario.i, 'x2'] <- ifelse(length(Selprob[which(Selprob.vars=='Velocity')]) == 0, 0, Selprob[which(Selprob.vars=='Velocity')])
  Scenarios[Scenario.i, 'x3'] <- ifelse(length(Selprob[which(Selprob.vars=='D50')]) == 0, 0, Selprob[which(Selprob.vars=='D50')])
  Scenarios[Scenario.i, 'x4'] <- ifelse(length(Selprob[which(Selprob.vars=='Temp')]) == 0, 0, Selprob[which(Selprob.vars=='Temp')])
  
  # Calculate the number of covariates in the model
  p <- sum(Scenarios[Scenario.i, c('x1', 'x2', 'x3', 'x4')] > 0.05)
  
  # Predictions versus observation
  YHat <- predict(model.frm, type='response')
  
  # Calculate the deviance of the FRM
  Dev.res <- CalcDevianceNegBin(model.frm, YObs)
  
  # Save some Goodness-of-fit criteria
  Scenarios[Scenario.i, 'Dev'] <-  Dev.res$ModelDeviance
  Scenarios[Scenario.i, 'DevExpl'] <-  Dev.res$ExplDeviance
  Scenarios[Scenario.i, 'Rsq'] <- CalcRsquare(YObs, YHat)
  Scenarios[Scenario.i, 'Rsqadj'] <- CalcRsquareAdj(YObs, YHat, p = p)
  Scenarios[Scenario.i, 'RMSELoo'] <- CalcMSE(YObs, LooResultTmp[mstop_best, ], root = T)
  
  # Plot the coefficients results
  pdf(paste0('out/models/all/tmp/FRM/Coef_', Y, '_df_', df, '_family_',Scenario$family, ".pdf"))
  plot(model.frm, ask = FALSE)
  dev.off()
  
# End of the loop for each scenario
}

# Save results to .csv file
write.csv2(Scenarios, 'out/models/all/FRMs_results.csv')



