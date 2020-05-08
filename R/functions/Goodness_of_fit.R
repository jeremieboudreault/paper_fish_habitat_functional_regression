# Requires packages
library('mgcv')
library('FDboost')

# Function to calculate the log-likelihood of the negative binomial
CalcLLNegBin <- function(yi, mui, theta)
{
  # If only one mu parameter is given (classical case)
  if (length(mui) == 1){
    sum(log(dnbinom(yi, mu=mui, size=theta)))
    # If multiple mu parameters are giver (given regression model such where yihat=mui)
  } else {
    sum(sapply(1:length(yi), function(w) log(dnbinom(yi[w], mu=mui[w], size=theta))))
  }
}

# Function to calculated the deviance of a NB model
CalcDevianceNegBin <- function(FittedModel, y)
{
  
  # Calculate null model (the log-likelihood of the null model will be the same for both GAM and FRM. GAM is used here)
  NullModel <- gam(y ~ 1, data=data.frame(y=y), family=nb(link = 'log'))
  
  # Spot check that the deviance is equal to 0 for this model
  if (summary(NullModel)$dev.expl > 0.001)
    stop('Null model has a deviance greater than 0')
  
  # Calculate the LogLikelihood for this null model
  NullModelLL <- logLik(NullModel)
  
  # Calculate the LL of a null model estimated by MOM
  mu <- mean(y) # Estimate of mu with MOM
  theta <- mu^2/(var(y)-mu) # Estimate of theta with MOM
  NullModelLL_2 <- CalcLLNegBin(y, mu, theta) # Calculate the log-lik of this model
  
  # Calculate the difference in % between the two calculated NullModelDeviance
  NullModelLL_diff <- as.vector(abs((NullModelLL-NullModelLL_2)/NullModelLL))
  
  # If the difference is greater than 1%
  if (NullModelLL_diff > 0.01)
    warning(paste0('Log-likelihood of the null models differ from ', round(NullModelLL_diff*100, 1), '%'))
  
  # Calculate saturated model (the log-likelihood of the full model will be the same for both GAM and FRM. GAM is used here)
  SaturatedModel <- gam(y ~ 0 + log(x), data=data.frame(y=y, x=y), family=nb(link = 'log'))
  
  # Spot check that the deviance is 100% for this model
  if (summary(SaturatedModel)$dev.expl < 0.999)
    stop('Saturated model has a deviance lower than 99.9%')
  
  # Calculate the LogLikelihood for this saturated model
  SaturatedModelLL <- logLik(SaturatedModel)
  
  # Calculate the loglikelihood of a saturated model in another way
  SaturatedModelLL_2 <- CalcLLNegBin(y, y, theta=999999) # mui = yi (all predictions = observation) + theta parameter tend to Inf as mui = yi
  
  # Calculate the difference between the two saturated model loglikelihood
  SaturatedModelLL_diff <-   as.vector(abs((SaturatedModelLL-SaturatedModelLL_2)/SaturatedModelLL))
  
  # If the difference is greater than 1%
  if (SaturatedModelLL_diff > 0.01)
    warning(paste0('Log-likelihood of the saturated models differ from ', round(SaturatedModelLL_diff*100, 1), '%'))
  
  # Calculate the log likelihood of the fitted model
  FittedModelLL <- logLik(FittedModel)
  
  # Calculate the null deviance
  NullDeviance <- as.vector(2 * (SaturatedModelLL - NullModelLL))
  
  # Calculate the deviance of the model
  ModelDeviance <- as.vector(2 * (FittedModelLL - NullModelLL))
  
  # Calculate the model deviance
  ResidualDeviance <- as.vector(2 * (SaturatedModelLL- FittedModelLL))
  
  # Calculate the explained deviance by the model
  ExplDeviance <- as.vector(1 - ResidualDeviance / NullDeviance)
  
  # Return the final result
  list(NullDeviance = NullDeviance, ResidualDeviance = ResidualDeviance, ModelDeviance = ModelDeviance, ExplDeviance = ExplDeviance)
}

# Function to calculate the Rsquare
CalcRsquare <- function(yobs, yhat){
  if (length(yobs) != length(yhat))
    stop('yobs and yhat do not have the same length')
  1 - sum((yobs-yhat)^2)/sum((yobs-mean(yobs))^2)
}

# Function to calculate the adjusted R2
CalcRsquareAdj <- function(yobs, yhat, p) {
  R2 <- CalcRsquare(yobs, yhat)
  n <- length(yobs)
  1 - (1-R2) * (n-1)/(n-p-1)
}

# Function to calculate the MSE / RMSE
CalcMSE <- function(yobs, yhat, root=FALSE){
  if (length(yobs) != length(yhat))
    stop('yobs and yhat do not have the same length')
  if (root) {
    exponent <- 1/2
    name <- 'RMSE'
  } else {
    exponent <- 1
    name <- 'MSE'
  }
  n <- length(yobs)
  r <- (1/n * sum((yhat-yobs)^2))^exponent
  names(r) <- name
  r
}
