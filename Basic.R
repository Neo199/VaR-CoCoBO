# Test problem for a constrained combinatorial Bayesian optimisation
# PRBOCS-VB method is used to solve the problem here.
# Approaching the constraints by PRing them before sampling and adding gumbel trick
# GA (genetic algorithm) for comparison and true answer

# This code is for the CONTAMINATION test problem in BOCS

# Author: Niyati Seth
# Date:  October 2025
# ---------------------------------------------------------
# LOAD Functions and Libraries
# ---------------------------------------------------------

source("~/Projects:Codes/P3Compute/sample_models.R")
source("~/Projects:Codes/P3Compute/thompson_svb.R")
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

set.seed(1)

# ---------------------------------------------------------
# Contamination Problem
# ---------------------------------------------------------

COST <- function(x, runlength = 1, seed = 10) {
  # x: binary vector (1 = facility installed)
  # runlength, seed: for reproducibility if stochastic
  n <- length(x)
  
  if (!is.numeric(x) || any(is.na(x)) ||
      any(x > 1) || any(x < 0) ||
      !is.numeric(runlength) || runlength <= 0 || runlength != round(runlength) ||
      !is.numeric(seed) || seed <= 0 || seed != round(seed)) {
    stop("Invalid inputs: x numeric in [0,1], runlength positive integer, seed positive integer.")
  }
  
  set.seed(seed)
  cost <- rep(15, n)  # cost for each facility
  constraint <- rep(0, n)
 
  # Only include facilities that are installed (x == 1)
  fn <- sum(cost * x)             # total cost of installed facilities
  for (k in 1:n) {
    constraint[k] <- cost[k] * x[k]
  }
  total_constraint <- sum(cost * x)     # same measure, but treated as constraint
  
  limit <- 100                    # budget limit
  
  return(list(fn = fn, constraint = constraint, limit = limit))
}



cost_prob <- function(x_mat, gamma = 10) {
  # x_mat: matrix where each row is one binary vector of length n
  # returns: vector of penalised objective values (cost - gamma * sum(constraint_violation_positive?))
  num_inputs <- nrow(x_mat)
  cost <- numeric(num_inputs)
  
  for (i in seq_len(num_inputs)) {
    res <- COST(x_mat[i, ], n_samples, seed = 10)
    cost[i] <- res$fn
  }
  return(cost)
}

cost_prob_lang <- function(x_mat, n_samples, gamma = 10) {
  # x_mat: matrix where each row is one binary vector of length n
  # returns: vector of penalised objective values (cost - gamma * sum(constraint_violation_positive?))
  num_inputs <- nrow(x_mat)
  out <- numeric(num_inputs)
  
  # Iterate over each input sample
  for (i in 1:num_inputs) {
    # Run contamination study
    cost_res <- COST(x_mat[i, ])
    cost <- cost_res$fn
    constraint <- sum(cost_res$constraint)
    
    # Compute total output
    out[i] <- cost - sum(gamma * (constraint - cost_res$limit))
  }
  return(out)
}

langmodel <- function(x_vals){
  cost_prob_lang(x_vals) 
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
n_vars <-10
evalBudget <-100
n_init <- 5
lambda <- 10
minSpend <- 30
order <- 2
seed <- 1


# ---------------------------------------------------------
# INITIAL SAMPLES FOR STATISTICAL MODELS
# ---------------------------------------------------------
x_vals <- sample_models(n_init, n_vars)
y_vals <- cost_prob_lang(x_vals)

x_vals
y_vals

# ---------------------------------------------------------
# DEFINE TRUE MODEL
# ---------------------------------------------------------
model <- function(x_vals){
  cost_prob(x_vals) 
}

# ---------------------------------------------------------
# RUN PROCS
# ---------------------------------------------------------

# Find number of iterations based on total budget
n_init <- nrow(x_vals)
n_iter <- evalBudget - n_init

# Train initial statistical model
# setup data for training
xTrain <- x_vals
yTrain <- y_vals

# Generate interaction terms for the input data
xTrain_in_comb <- order_effects(xTrain, order)
xTrain_in <- xTrain_in_comb$xTrain_in
inter_combos <- xTrain_in_comb$combos
# cat("xTrain with interaction terms","\n")
# print(xTrain_in)

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

# Initialize a data frame to store iteration results
optim_result <- matrix(0, evalBudget, n_vars)

# browser()
theta_ini <- rep(0.5, ncol(xTrain))

costs <- rep(1, n_vars)
n_vars <- length(costs)
theta_current <- rep(0.5, n_vars)
x_history <- list()
y_history <- numeric(evalBudget)

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


for (t in 1:evalBudget) {
  print(paste("prbocsvb_iteration_",t))
  
  # browser()
  # Sample new x using Gumbel rounding with constraints
  # x_new <- gumbel_acquisition(theta_current, vb_model, lambda, order, n_vars)
  
  stat_model <- function(theta) {
    thompson_sam_svb(theta, vb_model = vb_model, duplicate_cols, vb_data, order)
  }
  
  # Here you can still optimize theta if desired
  min_acq <- optim(theta_current, stat_model, method='L-BFGS-B', lower=1e-8, upper=0.999)
  expected_val <- min_acq$par
  cat("expected_val", expected_val, "\n")
  theta_current <- expected_val
  
  # Gumbel noise
  g <- -log(-log(runif(n_vars)))
  
  x_theta_current <- rbinom(n_vars, 1, theta_current)
  cat("X before constraint", x_theta_current, "\n")
  
  score <- numeric(n_vars)
  cost_res <- COST(x_theta_current)
  constr <- cost_res$constraint
  limit <- cost_res$limit
  
  for (j in seq_len(n_vars)) {
    # Penalize only violated constraints
    violations <- pmin(0, constr[j])
    penalty <- -violations  # positive if constraint violated
    score[j] <- theta_current[j] + g[j] + lambda * penalty
  }
  cat("Score:", score, "\n")
  
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
    eval_res <- COST(x_new)
    
    # Check satisfaction
    if (sum(eval_res$constraint) < eval_res$limit && sum(x_new)>=1) {
      # Violates constraint — undo last change and stop
      constraint_satisfied <- TRUE
      break
    }
  }
  
  # Final evaluation (guaranteed feasible)
  y_new <- COST(x_new)$fn
  
  # Store history and update data
  x_history[[t]] <- x_new
  y_history[t] <- y_new
  
  
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
  Y <- data_reduced[, 1]              # Response variable
  
  # Fit the variational Bayesian model
  vb_model <- svb.fit(
    X = X,
    Y = Y,
    family = "linear",     # For linear regression
    slab = "laplace",      # Default slab prior
    intercept = TRUE       # Include intercept in the model
  )
  
}

ga_model <- function(x) {
  x_mat <- matrix(x, nrow = 1)  # convert vector to 1-row matrix
  langmodel(x_mat)
}

GA_run <- ga(type = "binary", fitness = ga_model, nBits = n_vars,
             popSize = 100, maxiter = 1000, run = 100, monitor = FALSE)
ga_result <- list(solution = GA_run@solution, fitness_value = GA_run@fitnessValue)

ga_result$solution


