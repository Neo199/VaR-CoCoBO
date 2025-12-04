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

Contamination <- function(x, runlength, seed = 10) {
  # x: binary or {0,1} vector of length n (prevention decisions)
  # runlength: number of independent generations (positive integer)
  # seed: integer seed for reproducibility
  n <- length(x)
  # browser()
  # Input validation
  if (!is.numeric(x) || any(is.na(x)) ||
      any(x > 1) || any(x < 0) ||
      !is.numeric(runlength) || runlength <= 0 || runlength != round(runlength) ||
      !is.numeric(seed) || seed <= 0 || seed != round(seed)) {
    stop("Invalid inputs: x numeric in [0,1], runlength positive integer, seed positive integer.")
  }
  
  # Parameters
  nGen <- runlength         # number of independent generations
  u <- x                    # prevention binary decision variable
  X <- matrix(0, n, nGen)   # fraction contaminated at each stage for each generation
  epsilon <- rep(0.05, n)   # error probability
  p <- rep(0.1, n)          # proportion limit
  cost <- rep(1, n)         # cost for prevention at stage i
  
  # Beta parameters for initial contamination, contamination rate, restoration rate
  initialAlpha <- 1
  initialBeta <- 30
  contamAlpha <- 1
  contamBeta <- 17 / 3
  restoreAlpha <- 1
  restoreBeta <- 3 / 7
  
  
  # RNG
  set.seed(as.integer(seed))
  
  # Generate initial contamination per generation (vector length nGen)
  initialX <- rbeta(nGen, initialAlpha, initialBeta)
  
  # Generate contamination and restoration rates: matrices n x nGen
  Lambda <- matrix(rbeta(n * nGen, contamAlpha, contamBeta), nrow = n, ncol = nGen)
  Gamma  <- matrix(rbeta(n * nGen, restoreAlpha, restoreBeta), nrow = n, ncol = nGen)
  
  # Fill X: rows = 1..n variables (stages), columns = 1..nGen generations
  X[1, ] <- Lambda[1, ] * (1 - u[1]) * (1 - initialX) + (1 - Gamma[1, ] * u[1]) * initialX
  if (n >= 2) {
    for (i in 2:n) {
      X[i, ] <- Lambda[i, ] * (1 - u[i]) * (1 - X[i - 1, ]) + (1 - Gamma[i, ] * u[i]) * X[i - 1, ]
    }
  }
  
  # Limit and cost of contamination control
  limit <- 1 - epsilon
  fn <- sum(cost * u)
  
  # Constraint checking
  con <- matrix(0, nGen, n)
  for (j in 1:nGen) {
    con[j, ] <- X[, j] <= p
  }
  le <- sum(rowSums(con) == n)
  constraint <- rep(0, n)
  for (k in 1:n) {
    constraint[k] <- (sum(con[, k]) / runlength)
  }
  ConstraintCov <- cov(con)
  
  return(list(fn = fn, constraint = constraint, ConstraintCov = ConstraintCov, limit = limit))
}


contamination_prob <- function(x_mat, n_samples, gamma = 10) {
  # x_mat: matrix where each row is one binary vector of length n
  # returns: vector of penalised objective values (cost - gamma * sum(constraint_violation_positive?))
  num_inputs <- nrow(x_mat)
  cost <- numeric(num_inputs)
  
  for (i in seq_len(num_inputs)) {
    res <- Contamination(x_mat[i, ], n_samples, seed = 10)
    cost[i] <- res$fn
    }
  return(cost)
  }
    
# contamination_prob <- function(x_mat, n_samples, gamma = 10) {
#   # x_mat: matrix where each row is one binary vector of length n
#   # returns: vector of penalised objective values (cost - gamma * sum(constraint_violation_positive?))
#   num_inputs <- nrow(x_mat)
#   out <- numeric(num_inputs)
#   
#   for (i in seq_len(num_inputs)) {
#     res <- Contamination(x_mat[i, ], n_samples, seed = 10)
#     cost <- res$fn
#     constr <- res$constraint   # length n, positive = fraction_satisfied - limit
#     # If you want to penalize *violations* only (i.e. when fraction_satisfied < limit), use:
#     violations <- pmin(0, constr)  # negative values indicate violations
#     penalty <- -sum(violations)    # positive penalty = total shortfall fraction
#     out[i] <- cost + gamma * penalty
#     # If instead you want cost - gamma * constraint (as earlier), uncomment below:
#     # out[i] <- cost - sum(gamma * constr)
#   }
#   return(out)
# }

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
y_vals <- contamination_prob_lang(x_vals, 100)

x_vals
y_vals

# ---------------------------------------------------------
# DEFINE TRUE MODEL
# ---------------------------------------------------------
model <- function(x_vals){
  contamination_prob(x_vals, 100) 
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
  
  g <- -log(-log(runif(n_vars)))
  x_theta_current <- rbinom(n_vars, 1, theta_current)
  
  score <- numeric(n_vars)
  contam_res <- Contamination(x_theta_current, 100)
  constr <- contam_res$constraint
  limit <- contam_res$limit
  
  for (j in seq_len(n_vars)) {
    # Penalize only violated constraints
    violations <- pmin(0, constr[j])
    penalty <- -violations  # positive if constraint violated
    score[j] <- theta_current[j] + g[j] + lambda * penalty
  }
  cat("Score:", score, "\n")
  
  # Sort indices descending by score
  idx <- order(score, decreasing = TRUE)
  
  
  # ---------------------------------------------------------
  # Construct x_new incrementally until contamination constraint satisfied
  # ---------------------------------------------------------
  
  x_new <- rep(0, n_vars)
  
  for (k in seq_along(idx)) {
    # Turn on the next best bit
    x_new[idx[k]] <- 1
    
    # Evaluate contamination constraints for the current x
    contam_res <- Contamination(x_new, 100)
    constr <- contam_res$constraint
    limit <- contam_res$limit
    
    # Suppose Contamination()$constraint gives vector of constraint values
    # and all constraints are satisfied if they are >= 0
    if (all(constr > limit)) {
      message("All constraints satisfied at step ", k)
      break
    }
  }
  
  # Evaluate objective for this final feasible x
  y_new <- contam_res$fn
  
  # Store results
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

contamination_prob_lang <- function(x_mat, n_samples, gamma = 10) {
  # x_mat: matrix where each row is one binary vector of length n
  # returns: vector of penalised objective values (cost - gamma * sum(constraint_violation_positive?))
  num_inputs <- nrow(x_mat)
  out <- numeric(num_inputs)

  # Iterate over each input sample
  for (i in 1:num_inputs) {
    # Run contamination study
    contamination_result <- Contamination(x[i, ], n_samples, seed)
    cost <- contamination_result$fn
    constraint <- contamination_result$constraint
    
    # Compute total output
    out[i] <- cost - sum(gamma * constraint)
  }
  return(out)
}

langmodel <- function(x_vals){
  contamination_prob_lang(x_vals, 100) 
}

ga_model <- function(x) {
  x_mat <- matrix(x, nrow = 1)  # convert vector to 1-row matrix
  langmodel(x_mat)
}

GA_run <- ga(type = "binary", fitness = ga_model, nBits = n_vars,
             popSize = 100, maxiter = 1000, run = 100, monitor = FALSE)
ga_result <- list(solution = GA_run@solution, fitness_value = GA_run@fitnessValue)

ga_result$solution

# ---------------------------------------------------------
# PLOTS
# ---------------------------------------------------------
# 
# res <- prbocs_vb_result$data
# 
# # Plot objective function values versus iterations for normal optimization
# plot(1:nrow(res), res$y , type = "l", xlab = "Iterations", 
#      ylab = "Objective Function Value", main = "Objective Function vs Iterations", 
#      col = "red", xlim = c(1, nrow(res)))
