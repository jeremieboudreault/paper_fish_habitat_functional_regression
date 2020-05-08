####################################################################################
#     Modelling fish habitat selection using functional regression 🐟
#     Created by :    Jeremie Boudreault (Jeremie.Boudreault@ete.inrs.ca)
#     Current file :  Models_GAMs_best.R
#     Aim :           Find the best GAMS to model fish abundance
####################################################################################

#### Part 0 : Libraries, functions and variable definition

# Import library
library('dplyr') # Data.frame manipulation
library('ggplot2') # Data visualisation
library('mgcv') # For fitting Gams

# Import data per site of the SMR
SMR.per.site <- readRDS('data/SMR_2017_clean_per_site.Rda')

# Import results for Models_GAMs_all
GAMs_results <- read.csv2('out/models/all/GAMS_results.csv')[, -1]

# Extract best models / parameter for Fry / Parr
GAMs_best <- GAMs_results %>% 
  group_by(Y) %>%
  top_n(-1, RMSELoo)

#### Part 1 : Generate the best model for Fry

# The response variable 
Y <- 'nFry_M'

# The original value of the model
GAMs_best[which(GAMs_best$Y == Y), ]
k <- GAMs_best[which(GAMs_best$Y == Y), 'k']

# We will paste the original p-value to the article
(pval_original <- round(GAMs_best[which(GAMs_best$Y == Y), paste0('x', 0:4)], 3))

# Full model
Full_model_GAM <-  gam(   
                   formula= as.formula(paste0(Y, ' ~ s(Depth, k=', k, ') + s(Velocity, k=', k, ') + s(D50, k=', k, ') + s(Temp, k=', k, ')')),
                   data = SMR.per.site, 
                   family=nb(link = 'log'),
                   method = 'GCV.cp'
                   ) 

# Summary of the full model
(a.full <- summary(Full_model_GAM))

# We remove velocity from the model
Reduced_model_GAM_1 <-  gam(   
  formula= as.formula(paste0(Y, ' ~ s(Depth, k=', k, ') + s(D50, k=', k, ') + s(Temp, k=', k, ')')),
  data = SMR.per.site, 
  family=nb(link = 'log'),
  method = 'GCV.cp'
) 

# Summary of the reduced model 1
(a.r1 <- summary(Reduced_model_GAM_1))

# We remove D50 from the model
Reduced_model_GAM_2 <-  gam(   
  formula= as.formula(paste0(Y, '~ s(Depth, k=', k, ') + s(Temp, k=', k, ')')),
  data = SMR.per.site, 
  family=nb(link = 'log'),
  method = 'GCV.cp'
) 

# Summary of the reduced model 2
(a.r2 <- summary(Reduced_model_GAM_2))

# We remove Depth from the model
Reduced_model_GAM_3 <-  gam(   
  formula= as.formula(paste0(Y, ' ~ s(Temp, k=', k, ')')),
  data = SMR.per.site, 
  family=nb(link = 'log'),
  method = 'GCV.cp'
) 

# Summary of the reduced model 3 (All covariates significant)
(a.r3 <- summary(Reduced_model_GAM_3))

# Spot check on the AIC (Selected model have the lowest AIC)
Full_model_GAM$aic
Reduced_model_GAM_1$aic 
Reduced_model_GAM_2$aic  
Reduced_model_GAM_3$aic # The best model according to AIC

# Save the model
saveRDS(Reduced_model_GAM_3, paste0('out/models/GAM_', Y, '.rda'))

#### Part 3 : Generate the best model for Parr (1+)

# The response variable 
Y <- 'nParr1_M'

# The original value of the model
GAMs_best[which(GAMs_best$Y == Y), ]
k <- GAMs_best[which(GAMs_best$Y == Y), 'k']

# We will paste the original p-value to the article
(pval_original <- round(GAMs_best[which(GAMs_best$Y == Y), paste0('x', 0:4)], 3))

# Full model
Full_model_GAM <-  gam(   
  formula= as.formula(paste0(Y, ' ~ s(Depth, k=', k, ') + s(Velocity, k=', k, ') + s(D50, k=', k, ') + s(Temp, k=', k, ')')),
  data = SMR.per.site, 
  family=nb(link = 'log'),
  method = 'GCV.cp'
) 

# Summary of the full model
(a.full <- summary(Full_model_GAM))

# We remove velocity from the model
Reduced_model_GAM_1 <-  gam(   
  formula= as.formula(paste0(Y, ' ~ s(Depth, k=', k, ') + s(D50, k=', k, ') + s(Temp, k=', k, ')')),
  data = SMR.per.site, 
  family=nb(link = 'log'),
  method = 'GCV.cp'
) 

# Summary of the reduced model 1
(a.r1 <- summary(Reduced_model_GAM_1))

# If we were to remove Depth from the model
Reduced_model_GAM_2 <-  gam(   
  formula= as.formula(paste0(Y, ' ~ s(D50, k=', k, ') + s(Temp, k=', k, ')')),
  data = SMR.per.site, 
  family=nb(link = 'log'),
  method = 'GCV.cp'
) 

# Summary of the reduced model 2 (All covariates significative)
(a.r2 <- summary(Reduced_model_GAM_2))

# Just an overview of the AIC (very close to each other)
Full_model_GAM$aic
Reduced_model_GAM_1$aic # Lowest AIC is model 1
Reduced_model_GAM_2$aic 

# Save the best (model)
saveRDS(Reduced_model_GAM_1, paste0('out/models/GAM_', Y, '.rda')) # Based on AIC

