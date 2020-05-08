####################################################################################
#     Modelling fish habitat selection using functional regression 🐟
#     Created by :   Jeremie Boudreault (Jeremie.Boudreault@ete.inrs.ca)
#     Current file : Models_performance.R
#     Aim :          Calculate models performance for GAM and FRM
####################################################################################

#### Part 0 : Libraries, functions and variable definition

# Import library
library('tidyr') # 
library('dplyr') # Data.frame manipulation
library('ggplot2') # Data visualisation
library('mgcv') #GAMs
library('FDboost') #FRMs

# Import function to calculate deviance
source('R/functions/Goodness_of_fit.R')

# Import per site data
SMR.per.site <- readRDS('data/SMR_2017_clean_per_site.Rda')

# Import GAMs models
GAM.fry <- readRDS('out/models/GAM_nFry_M.rda')
GAM.parr1 <- readRDS('out/models/GAM_nParr1_M.rda')

# Import FRMs models
FRM.fry <- readRDS('out/models/FRM_nFry_M.rda')
FRM.parr1 <- readRDS('out/models/FRM_nParr1_M.rda')

# Import Leave-one-out predictions
GAM.fry.loo <- readRDS('out/predictions/YHatLoo_GAM_nFry_M.rda')
GAM.parr1.loo <- readRDS('out/predictions/YHatLoo_GAM_nParr1_M.rda')
FRM.fry.loo <- readRDS('out/predictions/YHatLoo_FRM_nFry_M.rda')
FRM.parr1.loo <- readRDS('out/predictions/YHatLoo_FRM_nParr1_M.rda')

### Part 1 : Basic exploration of the results

# Create a data.frame with the response for GAM / FRM
ResultsDF <- data.frame(
  Model = rep(c('FRM', 'GAM'), each=26*2),
  Yvar = rep(c('Fry', 'Parr1', 'Fry', 'Parr1'), each=26),
  YObs = c(SMR.per.site$nFry_M, SMR.per.site$nParr1_M, SMR.per.site$nFry_M, SMR.per.site$nParr1_M),
  YHat = c(predict(FRM.fry, type='response'), predict(FRM.parr1, type='response'), predict(GAM.fry, type = 'response'), predict(GAM.parr1, type = 'response')),
  YHatLoo = c(FRM.fry.loo, FRM.parr1.loo, GAM.fry.loo, GAM.parr1.loo)
)

# Melt to have the YHat and YHatLoo on the same graph
ResultsDFMelt <- pivot_longer(ResultsDF, cols = c(YHat, YHatLoo))

# Plot the results against observed values
ggplot(data=ResultsDFMelt, aes(x=YObs, y=value)) +
  geom_abline(aes(intercept=0, slope=1)) + 
  geom_point(aes(color=Model)) +
  xlab('Observed abundance (number of fish)') +
  ylab('Predicted abundance (number of fish)') + 
  facet_wrap(name~Yvar, nrow=2, ncol=2, scales='free')

#### Part 2 : Extract the selection probability for the FRM

# Fry functional model
round(summary(FRM.fry)$selprob, 3)

# Parr 1+ functional model
round(summary(FRM.parr1)$selprob, 3)

# Fry generalised additive model
GAM.fry
round(c(summary(GAM.fry)$p.pv, summary(GAM.fry)$s.pv), 4)

# Parr 1+ generalised additive model
GAM.parr1
round(c(summary(GAM.parr1)$p.pv, summary(GAM.parr1)$s.pv), 4)

#### Part 3 : Calculate the goodness-of-fit for the models

# Calculate the Dsq for each model
D2DF <- data.frame(Model=c(rep('FRM', 2), rep('GAM', 2)),
                        Yvar=rep(c('Fry', 'Parr1'), length.out=2 * 2),
                        D2=c(CalcDevianceNegBin(FRM.fry, SMR.per.site$nFry_M)$ExplDeviance,
                                  CalcDevianceNegBin(FRM.parr1, SMR.per.site$nParr1_M)$ExplDeviance,
                                  CalcDevianceNegBin(GAM.fry, SMR.per.site$nFry_M)$ExplDeviance,
                                  CalcDevianceNegBin(GAM.parr1, SMR.per.site$nParr1_M)$ExplDeviance))

# Calculate all other goodness-of-fit criteria
GofDF <-  ResultsDF %>% 
          group_by(Model, Yvar) %>%
          summarize(RMSE = CalcMSE(YObs, YHat, T), RMSELoo = CalcMSE(YObs, YHatLoo, T)) %>%
          left_join(D2DF, by=c("Model", "Yvar")) %>%
          select(Yvar, Model, D2, RMSE, RMSELoo) %>%
          mutate(D2 = round(D2, 3), RMSE = round(RMSE, 1), RMSELoo = round(RMSELoo, 1)) %>%
          arrange(desc(Model)) %>%
          pivot_wider(id_cols = Yvar, names_from='Model', values_from = c('D2', 'RMSE', 'RMSELoo'))

View(GofDF)
