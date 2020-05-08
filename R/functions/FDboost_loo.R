# Parallel packages required
library('parallel') # For parallel computing of the function
library('doParallel') # For parallel computing of the function
library('foreach') # For parallel computing of the function
library('FDboost') # For FDboost fitting

# Function to be performed a single leave-one-out estimate of the model
FDboost_loo_i <- function(i, # The observation to remove
                       formula, # The formula of the model
                       family, # The family of the model
                       data.list,  # The data for the model
                       mstop_max, # The boosting parameters
                       nu) {
  
  # The data.list should be formatting as follows :
  # * The first element data.list[[1]] is contains the (n x 1) Y observations variable
  # * The ith explanatory variable is stored as follows :
  #    * data.list[[i*2]] is a matrix of the (n x t) function observations (each row j is the function observation associated jth response variable y)
  #    * data.list[[i*2 + 1 ]] is a vector of (t x 1) giving the value at which the functional observations are made
  
  # Hence we can deduct the number of variables from 
  nvars <- (length(data.list)-1) / 2
  
  # And the variable names as liste
  vars.names <- names(data.list)[(1:nvars)*2]
    
  # Create a new subset of (n-1) observations
  subset_loo <- vector(mode='list', length=nvars*2 + 1)
  
  # Store the Y variable
  subset_loo[[1]] <- data.list[[1]][-i]
  
  # Store the functional observations Xi[-i, ]
  for (var.i in 1:nvars)
    subset_loo[[var.i * 2]] <- data.list[[var.i * 2]][-i, ]
  
  # Store the values were the functional observations are observed
  for (var.i in 1:nvars)
    subset_loo[[var.i * 2 + 1]] <- data.list[[var.i * 2 + 1]]
  
  # Add the names
  names(subset_loo) <- names(data.list)
  
  # Create a new subset for prediction of (1) observations
  subset_new <- vector(mode='list', length=nvars*2)
  
  # Store the functional observation Xi[i, ]
  for (var.i in 1:nvars)
    subset_new[[var.i * 2 - 1]] <- matrix(data.list[[var.i * 2]][i, ], nrow=1)
  
  # Store the values where the functional observations are observed
  for (var.i in 1:nvars)
    subset_new[[var.i * 2]] <- data.list[[var.i * 2 + 1]]
  
  # Store the names
  names(subset_new) <- names(data.list)[-1]
  
  # Adjust the functional regression model
  model.frm.loo <- FDboost(formula ,  timeformula=NULL, family=family, data=subset_loo,  control=boost_control(mstop = mstop_max, nu = nu))
  
  # Create a temp vector to store the result
  LooResultTmp <- rep(NA, mstop_max)
  
  # Start looping on each iteration of the model
  for (mstop in mstop_max:1)
  {
    
    # Update the model
    model.frm.loo.mstop <- model.frm.loo[mstop]
    
    # Generate the prediction
    LooResultTmp[mstop] <- predict(model.frm.loo.mstop, newdata=subset_new , type='response')
    
    # End of loop on each iteration of the model
  }
  
  return(LooResultTmp)
 
  # End of the function 
}

# Wrapper function for the FDboost package to perform leave-one-out estimate
FDboost_loo <- function(formula, # The formula of the model
                        family, # The family of the model
                        data.list,  # The data for the model
                        mstop_max=100, # The boosting parameters
                        nu=0.01) {

    # Calculate the number of cores
    no_cores <- detectCores() - 1
  
    # Create a cluster with all cores except one
    cl <- makeCluster(no_cores - 1)
          
    # Create a cluster with these cores
    registerDoParallel(cl) # register the cluster
    
    # Count the number of observation
    nObs <- length(data.list[[1]])
    
    # Use the foreach and %dopar% to calculate parallel GEV distributions
    res <- foreach(
      i          =  1:nObs,
      .combine   = "cbind",
      .packages  = c("FDboost"),
      .export = 'FDboost_loo_i') %dopar% {
      # Run the estimate of the FDboost on (n-ith) observation
      FDboost_loo_i(            i, # The observation to remove
                                formula = formula, # The formula of the model
                                family = family, # The family of the model
                                data.list = data.list,  # The data for the model
                                mstop_max=mstop_max, # The boosting parameters
                                nu=nu)
      }
              
# Shut down the cluster
stopCluster(cl)

# Return the result
res

# End of function
}
  
