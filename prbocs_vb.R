# Function for PRBOCS-VB-TS-optim

prbocs_vb_optim <- function(data, evalBudget, n_iter, n_vars, xTrain, xTrain_in, theta_current, order){
  vb_data <- data[,-1]
  
  # Initialize a data frame to store iteration results
  prbocs_vb_result <- matrix(0, evalBudget, n_vars)
  
  # browser()
  theta_current <-  rep(0.5, ncol(xTrain))
  
  # Find duplicate columns
  duplicate_cols <- which(duplicated(as.list(vb_data)))
  
  # Save the removed columns in a separate data frame
  removed_columns <- vb_data[, duplicate_cols, drop = FALSE]
  
  # Keep only unique columns
  data_reduced <- vb_data[, !duplicated(as.list(vb_data))]
  
  # Prepare data
  X <- as.matrix(data_reduced)  # Assuming the first column is the response variable
  Y <- data[, 1]              # Response variable
  
  # Fit the variational Bayesian model
  vb_model <- svb.fit(
    X = X,
    Y = Y,
    family = "linear",     # For linear regression
    slab = "laplace",      # Default slab prior
    intercept = TRUE       # Include intercept in the model
  )
  
  for (t in 1:n_iter) {

    stat_model <- function(theta) {
      thompson_sam_svb(theta, vb_model = vb_model, duplicate_cols, vb_data, order)
    }
    
    min_acq <- optim(theta_current, stat_model, method='L-BFGS-B', lower=1e-8, upper=0.99999)
    
    expected_val <- min_acq$par
    cat("expected_val", expected_val, "\n")
    x_new <- rbinom(length(expected_val), 1, expected_val)
    cat("New evaluation point", x_new, "\n")
    
    # browser()
    #Append new point to existing x_vals
    x_vals_updated <- rbind(xTrain, x_new)
    # Evaluate model objective at the new evaluation point
    x_new <- matrix(x_new, nrow = 1, ncol = n_vars)
    y_new <- model(x_new, seed)
    
    x_new_in_comb <- order_effects(x_new, order)
    x_new_in <- x_new_in_comb$xTrain_in
    
    data_new <- data.frame(y = y_new, x_new_in)
    data <- rbind(data, data_new)
    
    theta_current <- expected_val
    
    prbocs_vb_result[t,] <- expected_val
    
    vb_data <- data[,-1]
    
    # Find duplicate columns
    duplicate_cols <- which(duplicated(as.list(vb_data)))
    
    # Save the removed columns in a separate data frame
    removed_columns <- vb_data[, duplicate_cols, drop = FALSE]
    
    # Keep only unique columns
    data_reduced <- vb_data[, !duplicated(as.list(vb_data))]
    
    # Prepare data
    X <- as.matrix(data_reduced)  # Assuming the first column is the response variable
    Y <- data[, 1]              # Response variable
    
    # Fit the variational Bayesian model
    vb_model <- svb.fit(
      X = X,
      Y = Y,
      family = "linear",     # For linear regression
      slab = "laplace",      # Default slab prior
      intercept = TRUE       # Include intercept in the model
    )
  }
  result <- list(solution = tail(vb_data, n =1), data = data, 
                 model = vb_model) 
  return(result)
}
