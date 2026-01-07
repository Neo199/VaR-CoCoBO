# Test problem for a constrained combinatorial Bayesian optimisation
# Var-CoCoBO method is used to solve the problem here.
# Approaching the constraints by PRing them before sampling and adding gumbel trick
# GA (genetic algorithm) for comparison and true answer
# MINIMISATION PROBELM
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
  set.seed(seed)
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
  #constraint[k] = proportion of generations where stage k satisfies X[i] ≤ p, we want constraint > limit
  #limit = 1 − epsilon
  
  # Violation vector (how much safety falls short)
  violation_vec <- pmax(0, limit - constraint) # positive if constraint < limit
  # aggregate measure to flag infeasibility (max violation across stages)
  total_constraint_violation <- max(violation_vec)
  
  return(list(fn = fn, constraint = constraint, ConstraintCov = ConstraintCov, limit = limit,
              total_constraint_violation = total_constraint_violation))
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
n_vars <-25  
evalBudget <-100
n_init <- 10
lambda <- 10
order <- 2
seed <- 1


# ---------------------------------------------------------
# INITIAL SAMPLES FOR STATISTICAL MODELS
# ---------------------------------------------------------
x_vals <- sample_models(n_init, n_vars)
res <- Contamination(x_vals, 100)

x_vals
num_inputs <- nrow(x_vals)
y_vals <- numeric(num_inputs)

for (i in seq_len(num_inputs)) {
  res <- Contamination(x_vals[i, ], 100, seed = 10)
  y_vals[i] <- res$fn
  feasible[i] <- all(res$constraint > res$limit)
}

y_vals
feasible


# ---------------------------------------------------------
# DEFINE TRUE MODEL
# ---------------------------------------------------------
model <- function(x_vals){
  Contamination(x_vals, 100)
}

# ---------------------------------------------------------
# RUN VaR-CoCoBO
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

#Setup dataframe for stan_glm training 
X <- xTrain_in
y <- y_vals

#Create a dataframe
data <- data.frame(y=y , X)
data_history <- data.frame(y = y, x_vals, feasibility = data_feasible)

# Initialize a data frame to store iteration results
optim_result <- matrix(0, evalBudget, n_vars)

# browser()
theta_ini <- rep(0.5, ncol(xTrain))

costs <- rep(1, n_vars)
n_vars <- length(costs)
theta_current <- rep(0.5, n_vars)
x_history <- list()
y_history <- numeric(evalBudget)
data_feasible <- rep(NA, evalBudget)
data_feasible[1:n_init] <- as.logical(feasible)
total_violation <- rep(NA, evalBudget)

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

# ---------- Main optimization loop ----------
for (t in 1:evalBudget) {
  print(paste("prbocsvb_iteration_",t))

  stat_model <- function(theta) {
    thompson_sam_svb(theta, vb_model = vb_model, duplicate_cols, vb_data, order)
  }
  
  # Here you can still optimize theta if desired
  min_acq <- optim(theta_current, stat_model, method='L-BFGS-B', lower=1e-8, upper=0.999)
  expected_val <- min_acq$par
  cat("expected_val", expected_val, "\n")
  theta_current <- expected_val
  
  # Gumbel + sample
  g <- -log(-log(runif(n_vars)))
  x_theta_current <- rbinom(n_vars, 1, theta_current)
  cat("X before constraint", x_theta_current, "\n")
  
  # compute per-stage violations for sampled x (violation = pmax(0, limit - constraint))
  score <- numeric(n_vars)
  contam_res <- Contamination(x_theta_current, 100)
  constr <- contam_res$constraint     # proportion SAFE
  limit <- contam_res$limit           # required safety probability (1-epsilon)
  violation_sample <- pmax(0, limit - constr)  # positive when unsafe

  for (j in seq_len(n_vars)) {
    penalty <- lambda * sum(violation_sample)      # penalty per-coordinate: how much stage falls short of required safety
    # subtract penalty: higher penalty reduces score
    score[j] <- theta_current[j] + g[j] - penalty
  }
  cat("Score:", score, "\n")
  
  # Sort indices descending by score
  idx <- order(score, decreasing = TRUE)
  
  
  # greedily build a feasible x_new: add bit if resulting design remains feasible
  
  x_new <- rep(0, n_vars)
  
  
  for (k in seq_along(idx)) {
    # Check if this stage is currently unsafe
    contam_try <- Contamination(x_new, 100)
    constr_try <- contam_try$constraint
    limit_try <- contam_try$limit
    # browser()
    if (constr_try[idx[k]] < limit_try[idx[k]]) {
      # Stage k is unsafe → apply prevention
      x_new[idx[k]] <- 1
    }
    # else stage is already safe → leave prevention off
  }
  browser()

  # Evaluate contamination for x_new
  res <- Contamination(x_new,100)
  y_new <- res$fn
  x_new <- matrix(x_new, nrow = 1)
  data_feasible[t + n_init] <- as.logical(all(res$constraint > res$limit))
  
  # browser()
  # Save history
  data_hnew <- data.frame(y = y_new, x_new, feasibility = data_feasible[t + n_init])
  colnames(data_hnew) <- colnames(data_history)
  data_history <- rbind(data_history, data_hnew)
  
  x_history[[t]] <- x_new
  y_history[t] <- y_new
  
  cat(sprintf("Iteration %d: x = %s, y = %.3f\n", t, paste(x_new, collapse=""), y_new))
  
  # append new observation to dataset (include interactions)
  x_new <- matrix(x_new, nrow = 1, ncol = n_vars)
  x_new_in_comb <- order_effects(x_new, order)
  x_new_in <- x_new_in_comb$xTrain_in
  
  data_new <- data.frame(y = y_new, x_new_in)
  data <- rbind(data, data_new)
  
  theta_current <- expected_val
  
  prbocs_vb_result[t,] <- expected_val
  
  # prepare for next SVB fit: rebuild design and detect duplicates
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

test_feasible <- function(x) {
  res <- Contamination(x, 100)
  cat("\nx =", x, "\n")
  cat("constraint =", round(res$constraint, 3), "\n")
  cat("limit      =", round(res$limit, 3), "\n")
  cat("feasible?  =", all(res$constraint >= res$limit), "\n")
  cat("shortfall  =", round(res$limit - res$constraint, 3), "\n")
}

for (i in 1:length(x_history)) test_feasible(x_history[[i]])

# # ---------------------------------------------------------
# # Penalised wrapper used by GA 
# # ---------------------------------------------------------
# contamination_prob_lang <- function(x_mat, n_samples, gamma = 10, seed = 10) {
#   if (is.vector(x_mat)) x_mat <- matrix(x_mat, nrow = 1)
#   num_inputs <- nrow(x_mat)
#   out <- numeric(num_inputs)
#   
#   for (i in seq_len(num_inputs)) {
#     contamination_result <- Contamination(x_mat[i, ], n_samples, seed = seed)
#     cost <- contamination_result$fn
#     constraint <- contamination_result$constraint    # proportion SAFE
#     limit_vec <- contamination_result$limit
#    
#     # MINIMISATION objective: add penalty proportional to shortfall
#     out[i] <- cost - sum(gamma * constraint)
#   }
#   return(out)
# }
# 
langmodel <- function(x_vals) {
  num_inputs <- nrow(x_vals)
  out <- numeric(num_inputs)

  res <- Contamination(x_vals, 100)
  gamma <- 10
  cost <- res$fn
  # penalty must be based on SHORTFALL, not raw constraint
  shortfall <- pmax(0, res$limit - res$constraint)

  out[1] <- cost + gamma * sum(shortfall)
  return(out)
}


# GA expects a fitness to MAXIMISE, convert minimisation -> negative
ga_model <- function(x) {
  x_mat <- matrix(x, nrow = 1)
  -langmodel(x_mat)
}

GA_run <- ga(type = "binary", fitness = ga_model, nBits = n_vars,
             popSize = 100, maxiter = 1000, run = 100, monitor = FALSE)
ga_result <- list(solution = GA_run@solution, fitness_value = GA_run@fitnessValue)

print(ga_result)


# ---------------------------------------------------------
# PLOTS
# ---------------------------------------------------------
# ---------------------------------------------------------
# Optimisation Trace Plot
# ---------------------------------------------------------

# df_trace <- data.frame(
#   iter = 1:length(y_history),
#   y = y_history,
#   feasible = sapply(1:length(y_history), function(i) {
#     x_i <- optim_result[i, ]
#     cont_res <- Contamination(x_i,100)
#     cont_res$total_constraint <= cont_res$limit
#   })
# )
# 
# df_trace$best_feas <- NA
# best_so_far <- -Inf
# for (i in 1:nrow(df_trace)) {
#   if (df_trace$feasible[i]) {
#     best_so_far <- max(best_so_far, df_trace$y[i])
#   }
#   df_trace$best_feas[i] <- best_so_far
# }
# 
# library(ggplot2)
# 
# ggplot(df_trace, aes(x = iter)) +
#   geom_line(aes(y = best_feas), linewidth = 1, colour = "blue") +
#   geom_point(aes(y = y, colour = feasible)) +
#   scale_colour_manual(values = c("red", "darkgreen")) +
#   labs(
#     title = "PRBOCS-VB Optimisation Trace",
#     y = "Objective Value",
#     colour = "Feasible?"
#   ) +
#   theme_minimal()
