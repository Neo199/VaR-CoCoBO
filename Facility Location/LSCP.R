# Application for a constrained combinatorial Bayesian optimisation
# PRBOCS-VB method is used to solve the problem here.
# Approaching the constraints by PRing them before sampling and adding gumbel trick
# GA (genetic algorithm) for comparison and true answer

# This code is for the FACILITY LOCATION test problem
# LSCP
# MINIMISATION

# Author: Niyati Seth
# Date:  December 2025
# ---------------------------------------------------------
# LOAD Functions and Libraries
# ---------------------------------------------------------

source("~/Projects:Codes/VaR-CoCoBO/sample_models.R")
source("~/Projects:Codes/VaR-CoCoBO/thompson_svb.R")

library(tidyr)
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
library(maxcovr)
library(dplyr)

# ---------------------------------------------------------
# LOAD Data
# ---------------------------------------------------------

demand_data <- read.csv("PhD-Data/SF_data/SF_demand_205_centroid_uniform_weight.csv")
facility_loc <- read.csv("PhD-Data/SF_data/SF_store_site_16_longlat.csv")
distance_matrix <- read.csv("PhD-Data/SF_data/SF_network_distance_candidateStore_16_censusTract_205_new.csv")

demand_data
facility_loc
distance_matrix

# ---------------------------------------------------------
# LSCP DATA PREPARATION
# ---------------------------------------------------------

service_radius <- 5000

distance_matrix$covered <- as.integer(distance_matrix$distance <= service_radius)

A_df <- distance_matrix %>%
  dplyr::select(DestinationName, name, covered) %>%
  tidyr::pivot_wider(
    names_from  = name,
    values_from = covered,
    values_fill = list(covered = 0)
  )

A <- as.matrix(A_df[, -1])
n_demand <- nrow(A)
n_vars   <- ncol(A)

n_vars

# ---------------------------------------------------------
# LSCP PROBLEM
# ---------------------------------------------------------

lscp <- function(x) {
  # x: binary vector of length n_vars (facility selection)
  stopifnot(length(x) == n_vars)
  
  # Objective 
  fn <- sum(x)
  
  # Coverage per demand point
  coverage <- A %*% x    # n_demand x 1 vector
  
  # Constraint violations (same style as knapsack)
  # coverage >= 1  <=>  1 - coverage <= 0
  constraint <- 1 - coverage
  
  return(list(
    fn = fn,
    constraint = as.vector(constraint)
  ))
}

lscp_prob <- function(x_mat) {
  num_inputs <- nrow(x_mat)
  cost <- numeric(num_inputs)
  
  for (i in seq_len(num_inputs)) {
    res <- lscp(x_mat[i, ])
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
# Other Initialisations
# -------------------------
evalBudget <-250
n_init <- 5
lambda <- 1
minSpend <- 30
order <- 2
seed <- 1


# ---------------------------------------------------------
# INITIAL SAMPLES FOR STATISTICAL MODELS
# ---------------------------------------------------------
x_vals <- sample_models(n_init, n_vars)
y_vals <- lscp_prob(x_vals)

x_vals
y_vals

# ---------------------------------------------------------
# DEFINE TRUE MODEL
# ---------------------------------------------------------
model <- function(x_vals){
  lscp_prob(x_vals)
}

# ---------------------------------------------------------
# SAVE DIRECTORY SETUP
# ---------------------------------------------------------
# Create results directory if it doesn't exist
results_dir <- "results_LSCP"
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}

# Create unique identifier for this run
run_id <- format(Sys.time(), "%Y%m%d_%H%M%S")

# Start timing the optimization
start_time <- Sys.time()
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
theta_current <- rep(0.5, n_vars)
x_history <- list()
y_history <- numeric(evalBudget)
data_history <- data.frame(y = y, x_vals)

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
  
  x_theta_current <- rbinom(n_vars, 1, theta_current)
  cat("X before constraint", x_theta_current, "\n")
  
  # -------------------------
  # Gumbel-based ranking (LSCP)
  # -------------------------
  
  # Gumbel noise
  g <- -log(-log(runif(n_vars)))
  
  # SCORING STEP: compute score for each facility
  score <- numeric(n_vars)
  
  marginal_coverage <- colSums(A)
  
  # Standardize to mean=0, sd=1
  normalized_coverage <- (marginal_coverage - mean(marginal_coverage)) / sd(marginal_coverage)
  
  for (j in seq_len(n_vars)) {
    reward <- lambda * normalized_coverage[j]
    score[j] <- theta_current[j] + g[j] + reward
  }
  
  cat("Score:", score, "\n")
  
  # RANKING STEP: rank facilities by score (higher is better)
  idx <- order(score, decreasing = TRUE)
  
  # GREEDY CONSTRUCTION: build solution following the ranking
  x_new <- rep(0, n_vars)
  
  for (k in seq_along(idx)) {
    facility_id <- idx[k]
    
    # Check current coverage with x_new so far
    coverage_try <- as.vector(A %*% x_new)
    uncovered_try <- which(coverage_try < 1)
    
    if (length(uncovered_try) > 0) {
      if (sum(A[uncovered_try, facility_id]) > 0) {
        x_new[facility_id] <- 1
      }
    }
    
    # Early stopping
    coverage_final <- as.vector(A %*% x_new)
    if (all(coverage_final >= 1)) {
      break
    }
  }
  
  # Evaluate the solution
  res <- lscp(x_new)
  y_new <- res$fn
  
  data_hnew <- data.frame(y = y_new, t(x_new))
  colnames(data_hnew) <- colnames(data_history)
  data_history <- rbind(data_history, data_hnew)
  
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

# End timing the optimization
end_time <- Sys.time()
optimization_time <- as.numeric(difftime(end_time, start_time, units = "secs"))

cat(sprintf("\nOptimization completed in %.2f seconds (%.2f minutes)\n", 
            optimization_time, optimization_time/60))

# ---------------------------------------------------------
# SAVE RESULTS AFTER OPTIMIZATION
# ---------------------------------------------------------

# ---------------------------------------------------------
# SAVE RESULTS AFTER OPTIMIZATION
# ---------------------------------------------------------

# Save all results as RData
save(
  x_history,
  y_history,
  prbocs_vb_result,
  data_history,
  theta_current,
  optimization_time,
  file = file.path(results_dir, paste0("LSCP_BO_results_", run_id, ".RData"))
)

# Convert results to dataframes for CSV export
y_history_df <- data.frame(
  iteration = 1:evalBudget,
  objective = y_history,
  best_so_far = cummin(y_history)
)

theta_history_df <- as.data.frame(prbocs_vb_result)
colnames(theta_history_df) <- paste0("theta_", 1:n_vars)
theta_history_df$iteration <- 1:evalBudget

# Save as CSV
write.csv(y_history_df,
          file = file.path(results_dir, paste0("objective_history_", run_id, ".csv")),
          row.names = FALSE)

write.csv(theta_history_df,
          file = file.path(results_dir, paste0("theta_history_", run_id, ".csv")),
          row.names = FALSE)

write.csv(data_history,
          file = file.path(results_dir, paste0("solution_history_", run_id, ".csv")),
          row.names = FALSE)

# ---------------------------------------------------------
# SOLVE WITH EXACT SOLVER AND SAVE
# ---------------------------------------------------------

library(ompr)
library(ompr.roi)
library(ROI.plugin.glpk)
library(ROI)

# Time the exact solver
lp_start_time <- Sys.time()

# Pure Location Set Covering Problem (LSCP)
model <- MIPModel() %>%
  # Decision variables: which facilities to open
  add_variable(x[j], j = 1:n_vars, type = "binary") %>%

  # Constraints: each demand point i must be covered by at least one facility
  add_constraint(sum_expr(A[i,j] * x[j], j = 1:n_vars) >= 1, i = 1:nrow(A)) %>%

  # Objective: minimize number of facilities
  set_objective(sum_expr(x[j], j = 1:n_vars), sense = "min")

result <- solve_model(model, with_ROI(solver = "glpk"))
lp_end_time <- Sys.time()
lp_time <- as.numeric(difftime(lp_end_time, lp_start_time, units = "secs"))
lp_objective <- result$objective_value

# Save optimal solution
lp_solution <- data.frame(
  facility_id = 1:n_vars,
  selected = result$solution,
  lp_objective = lp_objective
)

write.csv(lp_solution,
          file = file.path(results_dir, paste0("lp_solution_", run_id, ".csv")),
          row.names = FALSE)

# ---------------------------------------------------------
# SAVE COMPARISON SUMMARY
# ---------------------------------------------------------

# Best BO solution
best_bo_iteration <- which.min(y_history)
best_bo_objective <- min(y_history)
best_bo_solution <- x_history[[best_bo_iteration]]

comparison_summary <- data.frame(
  method = c("GLPK", "Best_BO", "Final_BO"),
  objective = c(lp_objective, best_bo_objective, y_history[evalBudget]),
  iteration = c(NA, best_bo_iteration, evalBudget),
  gap_from_optimal = c(0,
                       best_bo_objective - lp_objective,
                       y_history[evalBudget] - lp_objective),
  time_seconds = c(lp_time, optimization_time, optimization_time)
)

write.csv(comparison_summary,
          file = file.path(results_dir, paste0("comparison_summary_", run_id, ".csv")),
          row.names = FALSE)

print(comparison_summary)

# ---------------------------------------------------------
# SAVE METADATA
# ---------------------------------------------------------

metadata <- list(
  run_id = run_id,
  run_date = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  evalBudget = evalBudget,
  n_init = n_init,
  lambda = lambda,
  order = order,
  seed = seed,
  n_vars = n_vars,
  n_demand = n_demand,
  service_radius = service_radius,
  lp_objective = lp_objective,
  lp_time_seconds = lp_time,
  best_bo_objective = best_bo_objective,
  convergence_iteration = best_bo_iteration,
  optimization_time_seconds = optimization_time,
  optimization_time_minutes = optimization_time / 60,
  time_per_iteration = optimization_time / evalBudget,
  speedup_vs_lp = lp_time / optimization_time
)

save(metadata,
        file = file.path(results_dir, paste0("metadata_", run_id, ".RData")))

library(GA)

# ---------------------------------------------------------
# GLOBAL COUNTER FOR FITNESS EVALUATIONS
# ---------------------------------------------------------
objective_counter <<- 0

# ---------------------------------------------------------
# Fitness function for GA (binary LSCP)
# ---------------------------------------------------------
lscp_fitness <- function(x) {
  
  # Count every evaluation
  objective_counter <<- objective_counter + 1
  
  x <- round(x)   # enforce binary (0/1)
  
  res <- lscp(x)  # your LSCP function from earlier
  
  obj  <- res$fn
  cons <- res$constraint
  
  # Penalty: number of uncovered demand points × large constant
  penalty <- sum(cons > 0) * 1000
  
  # GA maximizes → return negative penalized objective
  return(-(obj + penalty))
}

# ---------------------------------------------------------
# RUN GA WITH DEFAULTS + EARLY STOP (no improvement for 40 generations)
# ---------------------------------------------------------

ga_start_time <- Sys.time()

ga_LSCP <- ga(
  type    = "binary",
  fitness = lscp_fitness,
  nBits   = n_vars,
  run     = 100       # stop if no improvement for 40 generations
)

ga_end_time <- Sys.time()
ga_time_seconds <- as.numeric(difftime(ga_end_time, ga_start_time, units = "secs"))


# ---------------------------------------------------------
# RESULTS
# ---------------------------------------------------------

ga_best_solution <- ga_LSCP@solution[1, ]
ga_best_value    <- -ga_LSCP@fitnessValue

cat("\n====== GA LSCP RESULTS ======\n")
cat("Best number of facilities:", ga_best_value, "\n")
cat("Selected facilities (0/1):\n")
print(ga_best_solution)

# Optional: feasibility check
coverage_final <- as.vector(A %*% ga_best_solution)
cat("\nDemand points covered:", sum(coverage_final >= 1), "of", nrow(A), "\n")
cat("Feasible:", all(coverage_final >= 1), "\n")

# ---------------------------------------------------------
# FITNESS EVALUATION COUNT
# ---------------------------------------------------------

cat("\n===== Fitness Evaluations =====\n")
cat("Counted via fitness():   ", objective_counter, "\n")


# ---------------------------------------------------------
# TIME TAKEN
# ---------------------------------------------------------

cat("\n===== GA Runtime =====\n")
cat("Time taken (seconds): ", ga_time_seconds, "\n")
cat("Time taken (minutes): ", ga_time_seconds / 60, "\n")
