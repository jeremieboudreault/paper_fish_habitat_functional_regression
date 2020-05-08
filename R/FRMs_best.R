####################################################################################
#     Modelling fish habitat selection using functional regression 🐟
#     Created by :    Jeremie Boudreault (Jeremie.Boudreault@ete.inrs.ca)
#     Current file :  FRMs_best.R
#     Aim :           Find the best FRM models for fish abundance
####################################################################################

#### Part 0 : Libraries, functions and variable definition

# Import library
library('dplyr') # Data.frame manipulation
library('ggplot2') # Data visualisation
library('FDboost') # For fitting Frms

# Load the leave-one-out fitting of the FDboost
source('R/functions/FDboost_loo.R')

# Import data per site of the SMR
SMR.per.site <- readRDS('data/SMR_2017_clean_per_site.Rda')

# Import the functional observations
SMR.funobs <- readRDS('data/SMR_2017_clean_fun_obs.Rda')

# Import results for Models_FRMs_all
FRMs_best <- read.csv2('out/models/all/FRMs_results.csv')[, -1]

# Ensure we obtain the same resuls
set.seed('2912')

#### Part 1 : Generate the best model for Fry

# Set the response variable for Fry
Y <- 'nFry_M'

# The degree of freedom of the functional observation
df <- FRMs_best$df[FRMs_best$Y == Y]

# The family of the best model
family <- FRMs_best$family[FRMs_best$Y == Y]

# Probability of selection
round(FRMs_best[FRMs_best$Y == Y, c('x1', 'x2', 'x3', 'x4')], 3)

# Read the Fry model
FRM.fry <- readRDS(paste0('out/models/all/tmp/FRM/Model_nFry_M_df_', df, '_family_', family, '.Rds'))

# We will just refit the model without Depth as it had a probability of selection of 0.000
Vars <- c('Velocity', 'D50', 'Temp')
Subset <- list(Y=as.vector(unlist(SMR.per.site[, Y])),
               Velocity = t(as.matrix(SMR.funobs[which(SMR.funobs$Var=='Velocity'), 4:29])),
               Velocity.x = SMR.funobs[which(SMR.funobs$Var=='Velocity'), 2],
               D50 = t(as.matrix(SMR.funobs[which(SMR.funobs$Var=='D50'), 4:29])),
               D50.x = SMR.funobs[which(SMR.funobs$Var=='D50'), 2],
               Temp = t(as.matrix(SMR.funobs[which(SMR.funobs$Var=='Temp'), 4:29])),
               Temp.x = SMR.funobs[which(SMR.funobs$Var=='Temp'), 2])
formula.frm <- as.formula(paste0("Y~", paste0("bsignal(", Vars, ", ", Vars, ".x", ", df=", df, ", knots=6, degree=3, difference=1)", collapse="+")))
LooResultTmp  <- FDboost_loo(formula = formula.frm,
                             data.list = Subset,
                             family = NBinomial(),
                             mstop_max = 1000,
                             nu = 0.01)

# We will calculate the MSE on the leave-one-out predictions
YObs <- unlist(SMR.per.site[, Y])
MSE <- colSums((t(LooResultTmp)-YObs)^2)/nObs

# Defining the optimal mstop_best from the MSE
mstop_best <- which.min(MSE)

# We will save a MVP of the visualisation of the MSE and the number of iteration
pdf(paste0('out/models/all/tmp/FRM/Boosting_it_', Y, '_df_', df, '_family_', family, "_minus_Depth.pdf"), width=5, height=4)
print(ggplot(data=data.frame(mstop=1:mstop_max, MSE=MSE), aes(x=mstop, y=MSE)) + geom_line() + 
        geom_vline(aes(xintercept=mstop_best), lty=3) + ggtitle('Boosting iterations versus MSE') +
        annotate(geom='text', label=mstop_best, x=mstop_best, y=Inf, size=3, vjust=1.2))
dev.off()

# Adjust the final FRM with all observations
model.frm <- FDboost(formula.frm,  timeformula=NULL, family= NBinomial(), data=Subset,  control=boost_control(mstop = mstop_best, nu = 0.01))

# 
summary(model.frm)$selprob

# Predictions versus observation
YHat <- predict(model.frm, type='response')

# Calculate the deviance of the FRM
Dev.res <- CalcDevianceNegBin(model.frm, YObs)

# Plot the coefficients results
pdf(paste0('out/models/all/tmp/FRM/Coef_', Y, '_df_', df, '_family_', family, "_minus_Depth.pdf"))
plot(model.frm, ask = FALSE)
dev.off()

# Save some Goodness-of-fit criteria
Dev.res$ModelDeviance
Dev.res$ExplDeviance
CalcRsquare(YObs, YHat)
CalcRsquareAdj(YObs, YHat, p = p)
CalcMSE(YObs, LooResultTmp[mstop_best, ], root = T)

# Save the best model model
saveRDS(model.frm, paste0('out/models/FRM_', Y, '.rda'))

#### Part 2 : Generate the best model for 1 + parr

# Set the response variable for Fry
Y <- 'nParr1_M'

# The degree of freedom of the functional observation
df <- FRMs_best$df[FRMs_best$Y == Y]

# The family of the best model
family <- FRMs_best$family[FRMs_best$Y == Y]

# Probability of selection
round(FRMs_best[FRMs_best$Y == Y, c('x1', 'x2', 'x3', 'x4')], 3)

# Read the 1+ parr model model
FRM.parr1 <- readRDS(paste0('out/models/all/tmp/FRM/Model_nParr1_M_df_', df, '_family_', family, '.Rds'))

# We will test a model without velocity
Vars <- c('Depth', 'D50', 'Temp')
Subset <- list(Y=as.vector(unlist(SMR.per.site[, Y])),
               Depth = t(as.matrix(SMR.funobs[which(SMR.funobs$Var=='Depth'), 4:29])),
               Depth.x = SMR.funobs[which(SMR.funobs$Var=='Depth'), 2],
               D50 = t(as.matrix(SMR.funobs[which(SMR.funobs$Var=='D50'), 4:29])),
               D50.x = SMR.funobs[which(SMR.funobs$Var=='D50'), 2],
               Temp = t(as.matrix(SMR.funobs[which(SMR.funobs$Var=='Temp'), 4:29])),
               Temp.x = SMR.funobs[which(SMR.funobs$Var=='Temp'), 2])
formula.frm <- as.formula(paste0("Y~", paste0("bsignal(", Vars, ", ", Vars, ".x", ", df=", df, ", knots=6, degree=3, difference=1)", collapse="+")))
LooResultTmp  <- FDboost_loo(formula = formula.frm,
                             data.list = Subset,
                             family = NBinomial(),
                             mstop_max = 1000,
                             nu = 0.01)

# We will calculate the MSE on the leave-one-out predictions
YObs <- unlist(SMR.per.site[, Y])
MSE <- colSums((t(LooResultTmp)-YObs)^2)/nObs

# Defining the optimal mstop_best from the MSE
mstop_best <- which.min(MSE)

# We will save a MVP of the visualisation of the MSE and the number of iteration
pdf(paste0('out/models/all/tmp/FRM/Boosting_it_', Y, '_df_', df, '_family_', family, "_minus_Velocity.pdf"), width=5, height=4)
print(ggplot(data=data.frame(mstop=1:mstop_max, MSE=MSE), aes(x=mstop, y=MSE)) + geom_line() + 
        geom_vline(aes(xintercept=mstop_best), lty=3) + ggtitle('Boosting iterations versus MSE') +
        annotate(geom='text', label=mstop_best, x=mstop_best, y=Inf, size=3, vjust=1.2))
dev.off()

# Adjust the functional regression model with all observations
model.frm <- FDboost(formula.frm,  timeformula=NULL, family= NBinomial(), data=Subset,  control=boost_control(mstop = mstop_best, nu = 0.01))

# Selection probability of the new model
summary(model.frm)$selprob

# Plot the coefficients results
pdf(paste0('out/models/all/tmp/FRM/Coef_', Y, '_df_', df, '_family_', family, "_minus_Velocity.pdf"))
plot(model.frm, ask = FALSE)
dev.off()

# Save the model
saveRDS(model.frm, paste0('out/models/FRM_', Y, '.rda'))


