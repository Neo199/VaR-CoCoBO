# Test problem for a constrained combinatorial Bayesian optimisation
# PRBOCS-VB method is used to solve the problem here.
# Approaching the constraints by PRing them before sampling and adding gumbel trick
# GA (genetic algorithm) for comparison and true answer

# This code is for the 0/1 KNAPSACK test problem.  
# MAXIMISATION PROBELM
# BENCHMARK INSTANCE
# https://people.sc.fsu.edu/~jburkardt/datasets/knapsack_01/knapsack_01.html
# p08 - solution
# 110111000110100100000111
# 
# 10 instances with n_ini = 5, 10, 20

# Author: Niyati Seth
# Date:  October 2025

# ---------------------------------------------------------
# LOAD Functions and Libraries
# ---------------------------------------------------------

source("~/Projects:Codes/VaR-CoCoBO/sample_models.R")
source("~/Projects:Codes/VaR-CoCoBO/thompson_svb.R")
library(GA)
library(psych)
library(lpSolve)
library(ompr)
library(ompr.roi)
library(ggplot2)
library(dplyr)
library(ROI)
library(GPfit)
library(rstanarm)
library(rstan)
library(bayesplot)
library(sparsevb)
library(selectiveInference)


# ---------------------------------------------------------
# Knapsack Problem
# ---------------------------------------------------------

knapsack <- function(x) {
  # x: binary vector (1 = facility installed)
  # runlength, seed: for reproducibility if stochastic
  n <- length(x)
  stopifnot(n == n_vars)
  
  fn <- sum(values * x)
  
  # Only include items that are taken (x == 1)
  constraint <- rep(0, n)
  for (k in 1:n) {
    constraint[k] <- weights[k] * x[k]
  }    # For penalty
  
  total_constraint <- sum(weights * x)  # same measure, but treated as constraint
  limit <- W
  
  return(list(fn = fn, constraint = constraint, limit = limit, total_constraint = total_constraint))
}

knapsack_prob <- function(x_mat) {
  num_inputs <- nrow(x_mat)
  cost <- numeric(num_inputs)
  
  for (i in seq_len(num_inputs)) {
    res <- knapsack(x_mat[i, ])
    cost[i] <- res$fn
  }
  return(cost)
}


# ---------------------------------------------------------
# ORDER EFFECTS/THETA INTERACTION - to account interaction
# ---------------------------------------------------------
# # Function to generate interaction terms up to a given order
order_effects <- function(xTrain, order) {
  # Find dimensions of the matrix
  n_samp <- nrow(xTrain)
  n_vars <- ncol(xTrain)
  
  # Initialize result matrix
  xTrain_in <- xTrain
  
  # Generate interaction terms for each order up to ord_t
  for (ord_i in 2:order) {
    # Generate all combinations of indices (without diagonals)
    # cat("ord_i", ord_i, "\n")
    # cat("n_vars", n_vars, "\n")
    offdProd <- combn(n_vars, ord_i)
    
    # Generate products of input variables
    x_comb <- array(dim = c(n_samp, ncol(offdProd)))
    for (j in 1:ncol(offdProd)) {
      x_comb[, j] <- apply(xTrain[, offdProd[, j], drop = FALSE], 1, prod)
    }
    
    # Append interaction terms to the result matrix
    xTrain_in <- cbind(xTrain_in, x_comb)
  }
  
  return(list(xTrain_in = xTrain_in, combos = offdProd))
}

theta_interaction <- function(theta, order, n_vars){
  
  # Initialize result matrix
  theta <- matrix(theta, nrow = 1, ncol = n_vars)
  theta_in <- theta
  
  for (ord_i in 2:order) {
    # Generate all combinations of indices (without diagonals)
    # cat("ord_i", ord_i, "\n")
    # cat("n_vars", n_vars, "\n")
    offdProd <- combn(n_vars, ord_i)
    # Generate products of input variables
    theta_comb <- array(dim = c(1, ncol(offdProd)))
    for (j in 1:ncol(offdProd)) {
      theta_comb[, j] <- apply(theta[, offdProd[, j], drop = FALSE], 1, prod)
    }
    
    # Append interaction terms to the result matrix
    theta_in <- cbind(theta_in, theta_comb)
  }
  
  return(theta_in)
}

# ---------------------------------------------------------
# ACQUISITION FUNCTION
# ---------------------------------------------------------
# Function to compute thompson sampling
# this function is designed to be used with sparse
# variational bayes implementation

thompson_sam_svb <- function(theta_current, vb_model, duplicate_cols, vb_data, order){
  
  # browser()
  coeffs <- list()
  
  #Create a full mu vector of correct length
  full_mu <- numeric(ncol(vb_data))  # Initialize with zeros
  
  # Identify non-duplicate (i.e., retained) columns
  kept_cols <- setdiff(seq_along(full_mu), duplicate_cols)
  
  #Fill in mu values from reduced model
  full_mu[kept_cols] <- vb_model$mu
  
  # Copy values for removed (duplicate) columns
  for (col in duplicate_cols) {
    dup_column <- vb_data[, col]
    
    # Find the original column that matches the binary values of the duplicate column
    original_col <- which(apply(vb_data, 2, function(x) all(x == dup_column)) & !(seq_along(vb_data) %in% duplicate_cols))
    
    if (length(original_col) == 1) {
      full_mu[col] <- full_mu[original_col]
    } else {
      warning(paste("Could not uniquely match duplicate column:", names(vb_data)[col]))
    }
  }
  # Add a column of 1s to 'theta_current_in' to account for the intercept
  theta_current_in <- theta_interaction(theta_current, order, n_vars)
  theta_current_in <- c(1, theta_current_in)
  
  coeffs <- full_mu 
  coeffs <- c(vb_model$intercept, coeffs)
  # cat("Estimated coeffs", coeffs, "\n")
  y_pred <-  sum(theta_current_in * coeffs)
  
  return(y_pred = y_pred)
}

# ---------------------------------------------------------
# SET INPUTS
# ---------------------------------------------------------
# -------------------------
# Problem / knapsack setup
# -------------------------
n_vars <- 24                       # number of items / binary variables

# item values (v) and weights (w) - change to your instance
weights <- c(382745,
             799601,
             909247,
             729069,
             467902,
             44328,
             34610,
             698150,
             823460,
             903959,
             853665,
             551830,
             610856,
             670702,
             488960,
             951111,
             323046,
             446298,
             931161,
             31385,
             496951,
             264724,
             224916,
             169684)

values <- c( 825594,
             1677009,
             1676628,
             1523970,
             943972,
             97426,
             69666,
             1296457,
             1679693,
             1902996,
             1844992,
             1049289,
             1252836,
             1319836,
             953277,
             2067538,
             675367,
             853655,
             1826027,
             65731,
             901489,
             577243,
             466257,
             369261)
W <- 6404180                      # knapsack capacity

# ---------------------------------------------------------
# SINGLE EXPERIMENT RUN
# ---------------------------------------------------------
run_single_experiment <- function(n_init, instance_id, n_vars) {
  
  cat("Running n_init =", n_init, ", instance =", instance_id, "\n")
  
  evalBudget <- 200
  lambda <- 10
  minSpend <- 30
  order <- 2
  seed <- instance_id + 100
  set.seed(seed)
  # -----------------------------------
  # Generate initial x values (random)
  # -----------------------------------
  x_vals <- sample_models(n_init, n_vars)
  
  # -----------------------------------
  # Compute initial y
  # -----------------------------------
  y_vals <- knapsack_prob(x_vals)
  
  model <- function(x_vals) {
    knapsack_prob(x_vals)
  }
  
  n_init_actual <- nrow(x_vals)
  n_iter <- evalBudget - n_init_actual
  
  xTrain <- x_vals
  yTrain <- y_vals
  
  xTrain_in_comb <- order_effects(xTrain, order)
  xTrain_in <- xTrain_in_comb$xTrain_in
  inter_combos <- xTrain_in_comb$combos
  
  cDims <- dim(xTrain_in)
  nSamps <- cDims[1]
  nCoeffs <- cDims[2]
  # cat("nSamps, nCoeffs", c(nSamps, nCoeffs), "\n")
  
  hs_ss_sd <- sd(xTrain)
  n <- nrow(xTrain_in)
  p <- ncol(xTrain_in)
  
  #Setup dataframe for stan_glam training 
  X <- xTrain_in
  y <- y_vals
  
  #Create a dataframe
  data <- data.frame(y=y , X)
  data_history <- data.frame(y = y, x_vals)
  colnames(data_history) <- c("y", paste0("X", 1:n_vars))
  
  
  # Initialize a data frame to store iteration results
  optim_result <- matrix(0, evalBudget, n_vars)
  
  # browser()
  theta_ini <- rep(0.5, ncol(xTrain))
  
  costs <- rep(1, n_vars)
  theta_current <- rep(0.5, n_vars)
  x_history <- list()
  y_history <- numeric(evalBudget)
  feasible <- rep(NA, evalBudget)
  total_weight <- rep(NA, evalBudget)
  
  vb_data <- data[,-1]
  
  # Initialize a data frame to store iteration results
  prbocs_vb_result <- matrix(0, evalBudget, n_vars)
  
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
  
  
  for (t in 1:evalBudget-n_init) {
    print(paste("prbocsvb_iteration_",t))
    
    # browser()
    # Sample new x using Gumbel rounding with constraints
    # x_new <- gumbel_acquisition(theta_current, vb_model, lambda, order, n_vars)
    
    stat_model <- function(theta) {
      thompson_sam_svb(theta, vb_model = vb_model, duplicate_cols, vb_data, order)
    }
    
    # Here you can still optimize theta if desired
    min_acq <- optim(theta_current, stat_model, method='L-BFGS-B', lower=1e-8, upper=0.999,
                     control = list(fnscale = -1))
    expected_val <- min_acq$par
    cat("expected_val", expected_val, "\n")
    theta_current <- expected_val
    
    # Gumbel noise
    g <- -log(-log(runif(n_vars)))
    
    x_theta_current <- rbinom(n_vars, 1, theta_current)
    cat("X before constraint", x_theta_current, "\n")
    
    score <- numeric(n_vars)
    knapsack_res <- knapsack(x_theta_current)
    constr <- knapsack_res$constraint
    limit <- knapsack_res$limit
    
    cumsum_weight <- 0
    
    for (j in 1:n_vars) {
      if (x_theta_current[j] == 1) {
        projected_weight <- cumsum_weight + weights[j]
        penalty <- max(0, projected_weight - W)
        score[j] <- theta_current[j] + g[j] - lambda * penalty
        cumsum_weight <- projected_weight
      } else {
        score[j] <- theta_current[j] + g[j]
      }
    }
    cat("Score:", score, "\n")
    # browser()
    # Rank indices
    idx <- order(score, decreasing = TRUE)
    
    # Initialize binary vector
    x_new <- rep(0, n_vars)
    
    # Initialize constraint tracker
    constraint_satisfied <- FALSE
    
    # Start flipping bits in order of rank
    for (k in seq_along(idx)) {
      x_new[idx[k]] <- 1
      # browser()
      # Evaluate using provided function
      eval_res <- knapsack(x_new)
      
      # Check satisfaction
      if (sum(eval_res$constraint) >= eval_res$limit) {
        # Violates constraint — undo last change and stop
        x_new[idx[k]] <- 0
        constraint_satisfied <- TRUE
      }
    }
    
    
    
    # Evaluate knapsack for x_new
    res <- knapsack(x_new)
    y_new <- res$fn
    x_new <- matrix(x_new, nrow = 1)
    
    # Save history
    data_hnew <- data.frame(y = y_new, x_new)
    colnames(data_hnew) <- colnames(data_history)
    data_history <- rbind(data_history, data_hnew)
    
    x_history[[t]] <- x_theta_current
    y_history[t] <- y_new
    feasible[t] <- ifelse(res$total_constraint <= res$limit, 1, 0)
    total_weight[t] <- res$total_constraint
    
    
    # Update dataset (optional: append new x_new, y_new)
    # data <- rbind(data, c(y_new, theta_interaction(x_new, order, n_vars)))
    
    cat(sprintf("Iteration %d: x = %s, y = %.3f\n", t, paste(x_new, collapse=""), y_new))
    
    x_new <- matrix(x_new, nrow = 1, ncol = n_vars)
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
  
  # Save whatever you want
  result_list <- list(
    n_init = n_init,
    instance = instance_id,
    data = data_history
  )
  
  return(result_list)
}

# ---------------------------------------------------------
# MAIN WRAPPER FOR ALL RUNS
# ---------------------------------------------------------
run_all_experiments <- function() {
  
  n_inis <- c(5, 10, 20)
  n_instances <- 10
  
  base_dir <- "~/Projects:Codes/P3Compute/knapsack_gumbel_results"
  dir.create(base_dir, showWarnings = FALSE)
  
  for (n_init in n_inis) {
    for (inst in 1:n_instances) {
      
      # Folder: results/n_ini_5/instance_1/
      dir_path <- file.path(base_dir, paste0("n_ini_", n_init))
      dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
      # -----------------------
      # TIMING START
      # -----------------------
      t_start <- Sys.time()
      # Run your experiment
      res <- run_single_experiment(
        n_init = n_init,
        instance_id = inst,
        n_vars = n_vars
      )
      # -----------------------
      # TIMING END
      # -----------------------
      t_end <- Sys.time()
      runtime_sec <- as.numeric(difftime(t_end, t_start, units = "secs"))
      
      # Add timing inside result object
      res$runtime_sec <- runtime_sec
      # Save result to .RData
      save(res, file = file.path(dir_path, paste0("instance_", inst, ".RData")))
      
      cat("Saved:", dir_path, "\n")
    }
  }
}

# Run everything
run_all_experiments()


