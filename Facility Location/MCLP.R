# Application for a constrained combinatorial Bayesian optimisation
# PRBOCS-VB method is used to solve the problem here.
# Approaching the constraints by PRing them before sampling and adding gumbel trick
# GA (genetic algorithm) for comparison and true answer

# This code is for the FACILITY LOCATION test problem
# MCLP
# MAXIMISATION

# Author: Niyati Seth
# Date:  January 2026
# ---------------------------------------------------------
# LOAD Functions and Libraries
# ---------------------------------------------------------

source("~/Projects:Codes/P3Compute/sample_models.R")
source("~/Projects:Codes/P3Compute/thompson_svb.R")

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

# ---------------------------------------------------------
# MCLP DATA PREPARATION
# ---------------------------------------------------------

service_radius <- 5000
max_facilities <- 4  # Maximum number of facilities to open

# Mark which demand-facility pairs are within service radius
distance_matrix$covered <- as.integer(distance_matrix$distance <= service_radius)

# Create coverage matrix A and population vector
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

# Extract population for each demand point (census tract)
# Match demand points to population
demand_lookup <- demand_data %>%
  dplyr::select(NAME, POP2000) %>%
  rename(DestinationName = NAME, population = POP2000)

population_df <- A_df %>%
  dplyr::select(DestinationName) %>%
  left_join(demand_lookup, by = "DestinationName")

population <- population_df$population

cat("Number of facilities:", n_vars, "\n")
cat("Number of demand points:", n_demand, "\n")
cat("Total population:", sum(population), "\n")
cat("Max facilities allowed:", max_facilities, "\n")

# ---------------------------------------------------------
# MCLP PROBLEM
# ---------------------------------------------------------

mclp <- function(x) {
  # x: binary vector of length n_vars (facility selection)
  stopifnot(length(x) == n_vars)
  
  # Check facility constraint
  n_facilities <- sum(x)
  
  # Coverage: which demand points are covered (1 if covered by at least one facility)
  coverage_vector <- pmin(A %*% x, 1)  # Binary: 1 if covered, 0 otherwise
  
  # Objective: total population covered (MAXIMIZE)
  # For minimization in optimization, return negative
  total_covered_pop <- sum(coverage_vector * population)
  
  # Constraint violation: number of facilities must be <= max_facilities
  # constraint: n_facilities - max_facilities <= 0
  constraint <- n_facilities - max_facilities
  
  return(list(
    fn = total_covered_pop,  
    coverage = total_covered_pop,  # Actual coverage for reporting
    constraint = constraint,
    n_facilities = n_facilities
  ))
}

mclp_prob <- function(x_mat) {
  num_inputs <- nrow(x_mat)
  cost <- numeric(num_inputs)
  
  for (i in seq_len(num_inputs)) {
    res <- mclp(x_mat[i, ])
    cost[i] <- res$fn
  }
  return(cost)
}

# ---------------------------------------------------------
# ORDER EFFECTS/THETA INTERACTION
# ---------------------------------------------------------

order_effects <- function(xTrain, order) {
  n_samp <- nrow(xTrain)
  n_vars <- ncol(xTrain)
  
  xTrain_in <- xTrain
  
  for (ord_i in 2:order) {
    offdProd <- combn(n_vars, ord_i)
    x_comb <- array(dim = c(n_samp, ncol(offdProd)))
    for (j in 1:ncol(offdProd)) {
      x_comb[, j] <- apply(xTrain[, offdProd[, j], drop = FALSE], 1, prod)
    }
    xTrain_in <- cbind(xTrain_in, x_comb)
  }
  
  return(list(xTrain_in = xTrain_in, combos = offdProd))
}

theta_interaction <- function(theta, order, n_vars){
  theta <- matrix(theta, nrow = 1, ncol = n_vars)
  theta_in <- theta
  
  for (ord_i in 2:order) {
    offdProd <- combn(n_vars, ord_i)
    theta_comb <- array(dim = c(1, ncol(offdProd)))
    for (j in 1:ncol(offdProd)) {
      theta_comb[, j] <- apply(theta[, offdProd[, j], drop = FALSE], 1, prod)
    }
    theta_in <- cbind(theta_in, theta_comb)
  }
  
  return(theta_in)
}

# ---------------------------------------------------------
# ACQUISITION FUNCTION
# ---------------------------------------------------------

thompson_sam_svb <- function(theta_current, vb_model, duplicate_cols, vb_data, order){
  coeffs <- list()
  
  full_mu <- numeric(ncol(vb_data))
  kept_cols <- setdiff(seq_along(full_mu), duplicate_cols)
  full_mu[kept_cols] <- vb_model$mu
  
  for (col in duplicate_cols) {
    dup_column <- vb_data[, col]
    original_col <- which(apply(vb_data, 2, function(x) all(x == dup_column)) & 
                            !(seq_along(vb_data) %in% duplicate_cols))
    if (length(original_col) == 1) {
      full_mu[col] <- full_mu[original_col]
    } else {
      warning(paste("Could not uniquely match duplicate column:", names(vb_data)[col]))
    }
  }
  
  theta_current_in <- theta_interaction(theta_current, order, n_vars)
  theta_current_in <- c(1, theta_current_in)
  
  coeffs <- full_mu 
  coeffs <- c(vb_model$intercept, coeffs)
  y_pred <- sum(theta_current_in * coeffs)
  
  return(y_pred = y_pred)
}

# ---------------------------------------------------------
# SET INPUTS
# ---------------------------------------------------------

evalBudget <- 250
n_init <- 5
lambda <- 1
order <- 2
seed <- 1

# ---------------------------------------------------------
# INITIAL SAMPLES
# ---------------------------------------------------------

set.seed(seed)

# Generate random initial samples (no constraint enforcement)
x_vals <- sample_models(n_init, n_vars)

y_vals <- mclp_prob(x_vals)

cat("Initial solutions:\n")
print(x_vals)
cat("Number of facilities in each initial solution:", apply(x_vals, 1, sum), "\n")
cat("\nInitial objectives (negative population covered):\n")
print(y_vals)


# ---------------------------------------------------------
# SAVE DIRECTORY SETUP
# ---------------------------------------------------------

results_dir <- "results_MCLP"
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}

run_id <- format(Sys.time(), "%Y%m%d_%H%M%S")
start_time <- Sys.time()

# ---------------------------------------------------------
# RUN OPTIMIZATION
# ---------------------------------------------------------

n_iter <- evalBudget - n_init

xTrain <- x_vals
yTrain <- y_vals

xTrain_in_comb <- order_effects(xTrain, order)
xTrain_in <- xTrain_in_comb$xTrain_in

X <- xTrain_in
y <- y_vals
data <- data.frame(y = y, X)

theta_current <- rep(0.5, n_vars)
x_history <- list()
y_history <- numeric(evalBudget)
coverage_history <- numeric(evalBudget)
data_history <- data.frame(y = y, x_vals)

# Store initial results
for (i in 1:n_init) {
  res <- mclp(x_vals[i,])
  y_history[i] <- res$fn
  coverage_history[i] <- res$coverage
}

prbocs_vb_result <- matrix(0, evalBudget, n_vars)

vb_data <- data[, -1]
duplicate_cols <- which(duplicated(as.list(vb_data)))
data_reduced <- vb_data[, !duplicated(as.list(vb_data))]

X_vb <- as.matrix(data_reduced)
Y_vb <- data[, 1]

vb_model <- svb.fit(
  X = X_vb,
  Y = Y_vb,
  family = "linear",
  slab = "laplace",
  intercept = TRUE
)

for (t in 1:evalBudget) {
  print(paste("VaRCBO Iteration", t))
  
  stat_model <- function(theta) {
    -thompson_sam_svb(theta, vb_model = vb_model, duplicate_cols, vb_data, order)
  }
  
  max_acq <- optim(theta_current, stat_model, method = 'L-BFGS-B', 
                   lower = 1e-8, upper = 0.999)
  expected_val <- max_acq$par
  theta_current <- expected_val
  
  x_theta_current <- rbinom(n_vars, 1, theta_current)
  cat("X before constraint", x_theta_current, "\n")
  
  # ---------------------------------------------------------
  # Gumbel-based selection for MCLP
  # ---------------------------------------------------------
  
  # Gumbel noise
  g <- -log(-log(runif(n_vars)))
  
  # SCORING STEP: MCLP-specific scoring
  score <- numeric(n_vars)
  
  # Each facility's value: population it can cover
  marginal_coverage <- colSums(A * population)
  
  # Standardize to match gumbel/theta scale (0 mean, unit variance)
  normalized_coverage <- (marginal_coverage - mean(marginal_coverage)) / 
    (sd(marginal_coverage))
  
  for (j in seq_len(n_vars)) {
    # Score: exploration (theta + gumbel) + exploitation (coverage potential)
    score[j] <- theta_current[j] + g[j] + lambda * normalized_coverage[j]
  }
  
  cat("Score:", score, "\n")
  
  # RANKING STEP
  idx <- order(score, decreasing = TRUE)
  
  # CONSTRUCTION: select top max_facilities
  x_new <- rep(0, n_vars)
  x_new[idx[1:max_facilities]] <- 1
  
  # Evaluate
  res <- mclp(x_new)
  y_new <- res$fn
  coverage_new <- res$coverage
  
  cat(sprintf("Iteration %d: x = %s, y = %.3f\n", t, paste(x_new, collapse=""), y_new))
  cat(sprintf("Iteration %d: Facilities = %d, Coverage = %.0f people (%.1f%%)\n", 
              t, sum(x_new), coverage_new, 100 * coverage_new / sum(population)))
  
  data_hnew <- data.frame(y = y_new, t(x_new))
  colnames(data_hnew) <- colnames(data_history)
  data_history <- rbind(data_history, data_hnew)
  
  x_history[[t]] <- x_new
  y_history[t] <- y_new
  coverage_history[t] <- coverage_new
  
  # Update model
  x_new_matrix <- matrix(x_new, nrow = 1, ncol = n_vars)
  x_new_in_comb <- order_effects(x_new_matrix, order)
  x_new_in <- x_new_in_comb$xTrain_in
  
  data_new <- data.frame(y = y_new, x_new_in)
  data <- rbind(data, data_new)
  
  prbocs_vb_result[t, ] <- expected_val
  
  vb_data <- data[, -1]
  duplicate_cols <- which(duplicated(as.list(vb_data)))
  data_reduced <- vb_data[, !duplicated(as.list(vb_data))]
  
  X_vb <- as.matrix(data_reduced)
  Y_vb <- data[, 1]
  
  vb_model <- svb.fit(
    X = X_vb,
    Y = Y_vb,
    family = "linear",
    slab = "laplace",
    intercept = TRUE
  )
}

end_time <- Sys.time()
optimization_time <- as.numeric(difftime(end_time, start_time, units = "secs"))

cat(sprintf("\nOptimization completed in %.2f seconds (%.2f minutes)\n", 
            optimization_time, optimization_time/60))

# ---------------------------------------------------------
# SAVE RESULTS
# ---------------------------------------------------------

save(
  x_history,
  y_history,
  coverage_history,
  prbocs_vb_result,
  data_history,
  theta_current,
  optimization_time,
  file = file.path(results_dir, paste0("MCLP_BO_results_", run_id, ".RData"))
)

results_df <- data.frame(
  iteration = 1:evalBudget,
  objective = y_history,
  coverage = coverage_history,
  coverage_pct = 100 * coverage_history / sum(population),
  best_coverage = cummax(coverage_history)
)

write.csv(results_df, 
          file = file.path(results_dir, paste0("objective_history_", run_id, ".csv")),
          row.names = FALSE)

theta_history_df <- as.data.frame(prbocs_vb_result)
colnames(theta_history_df) <- paste0("theta_", 1:n_vars)
theta_history_df$iteration <- 1:evalBudget

write.csv(theta_history_df, 
          file = file.path(results_dir, paste0("theta_history_", run_id, ".csv")),
          row.names = FALSE)

write.csv(data_history, 
          file = file.path(results_dir, paste0("solution_history_", run_id, ".csv")),
          row.names = FALSE)

# ---------------------------------------------------------
# SOLVE WITH EXACT SOLVER
# ---------------------------------------------------------

library(ROI.plugin.glpk)

lp_start_time <- Sys.time()

# Maximum Coverage Location Problem (MCLP)
model <- MIPModel() %>%
  add_variable(x[j], j = 1:n_vars, type = "binary") %>%
  add_variable(y[i], i = 1:n_demand, type = "binary") %>%
  
  # Coverage constraint: y[i] = 1 only if at least one covering facility is open
  add_constraint(y[i] <= sum_expr(A[i,j] * x[j], j = 1:n_vars), i = 1:n_demand) %>%
  
  # Facility limit constraint
  add_constraint(sum_expr(x[j], j = 1:n_vars) <= max_facilities) %>%
  
  # Objective: maximize covered population
  set_objective(sum_expr(population[i] * y[i], i = 1:n_demand), sense = "max")

result <- solve_model(model, with_ROI(solver = "glpk"))
lp_end_time <- Sys.time()
lp_time <- as.numeric(difftime(lp_end_time, lp_start_time, units = "secs"))
lp_objective <- result$objective_value

lp_solution <- data.frame(
  facility_id = 1:n_vars,
  selected = as.integer(get_solution(result, x[j])$value)
)

lp_summary <- data.frame(
  optimal_coverage = lp_objective,
  optimal_pct = 100 * lp_objective / sum(population),
  n_facilities = sum(lp_solution$selected)
)

write.csv(lp_solution, 
          file = file.path(results_dir, paste0("lp_solution_", run_id, ".csv")),
          row.names = FALSE)

cat(sprintf("\nOptimal solution covers %.0f people (%.1f%%)\n", 
            lp_objective, 100 * lp_objective / sum(population)))

# ---------------------------------------------------------
# COMPARISON SUMMARY
# ---------------------------------------------------------

best_bo_iteration <- which.max(coverage_history)
best_bo_coverage <- max(coverage_history)
best_bo_solution <- x_history[[best_bo_iteration]]

comparison_summary <- data.frame(
  method = c("GLPK_Optimal", "Best_BO", "Final_BO"),
  coverage = c(lp_objective, best_bo_coverage, coverage_history[evalBudget]),
  coverage_pct = c(100 * lp_objective / sum(population),
                   100 * best_bo_coverage / sum(population),
                   100 * coverage_history[evalBudget] / sum(population)),
  iteration = c(NA, best_bo_iteration, evalBudget),
  gap_from_optimal = c(0, 
                       lp_objective - best_bo_coverage,
                       lp_objective - coverage_history[evalBudget]),
  gap_pct = c(0,
              100 * (lp_objective - best_bo_coverage) / lp_objective,
              100 * (lp_objective - coverage_history[evalBudget]) / lp_objective),
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
  problem_type = "MCLP",
  evalBudget = evalBudget,
  n_init = n_init,
  lambda = lambda,
  order = order,
  seed = seed,
  n_vars = n_vars,
  n_demand = n_demand,
  max_facilities = max_facilities,
  service_radius = service_radius,
  total_population = sum(population),
  lp_objective = lp_objective,
  lp_coverage_pct = 100 * lp_objective / sum(population),
  lp_time_seconds = lp_time,
  best_bo_coverage = best_bo_coverage,
  best_bo_coverage_pct = 100 * best_bo_coverage / sum(population),
  convergence_iteration = best_bo_iteration,
  optimization_time_seconds = optimization_time,
  optimization_time_minutes = optimization_time / 60,
  time_per_iteration = optimization_time / evalBudget
)

save(metadata, file = file.path(results_dir, paste0("metadata_", run_id, ".RData")))

cat("\n=== MCLP Optimization Complete ===\n")
cat(sprintf("Total population: %d\n", sum(population)))
cat(sprintf("Optimal coverage: %.0f (%.1f%%)\n", 
            lp_objective, metadata$lp_coverage_pct))
cat(sprintf("Best BO coverage: %.0f (%.1f%%)\n", 
            best_bo_coverage, metadata$best_bo_coverage_pct))
cat(sprintf("Gap from optimal: %.1f%%\n", 
            comparison_summary$gap_pct[2]))

library(GA)

# ---------------------------------------------------------
# SINGLE GA COMPARATOR FOR MCLP
# ---------------------------------------------------------

run_ga_mclp_single <- function(evalBudget = 250, seed = 123) {
  
  cat("\n==============================================\n")
  cat("Running GA Comparator\n")
  cat("==============================================\n")
  
  set.seed(seed)
  
  start_time <- Sys.time()
  
  # Fitness function for GA (maximize coverage)
  fitness_mclp <- function(x) {
    # Check constraint: exactly max_facilities
    if (sum(x) != max_facilities) {
      # Penalize heavily for constraint violation
      return(-1e10)
    }
    
    # Evaluate coverage
    res <- mclp(x)
    return(res$fn)  # Return coverage (positive value to maximize)
  }
  
  # Run GA
  ga_result <- ga(
    type = "binary",
    fitness = fitness_mclp,
    nBits = n_vars,
    maxiter = evalBudget,  # Match BO budget
    popSize = 50,
    pcrossover = 0.8,
    pmutation = 0.1,
    elitism = base::max(1, round(50 * 0.05)),
    seed = seed,
    monitor = TRUE
  )
  
  end_time <- Sys.time()
  runtime <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  # Extract best solution
  best_solution <- as.vector(ga_result@solution[1, ])
  best_fitness <- ga_result@fitnessValue
  # Get convergence history
  fitness_history <- ga_result@summary[, "max"]  # Best fitness per generation
  
  cat(sprintf("\nGA completed in %.2f seconds (%.2f minutes)\n", 
              runtime, runtime/60))
  cat(sprintf("Best coverage: %.0f people (%.2f%%)\n", 
              best_fitness, 100 * best_fitness / sum(population)))
  cat(sprintf("Best solution: %s\n", paste(best_solution, collapse="")))
  
  # Store results
  ga_results <- list(
    seed = seed,
    best_solution = best_solution,
    best_coverage = best_fitness,
    best_coverage_pct = 100 * best_fitness / sum(population),
    convergence_history = fitness_history,
    runtime_sec = runtime,
    runtime_min = runtime / 60,
    ga_object = ga_result
  )
  
  # Save results
  ga_dir <- "results_MCLP/GA_results"
  dir.create(ga_dir, showWarnings = FALSE, recursive = TRUE)
  
  save(ga_results, 
       file = file.path(ga_dir, "ga_single_run.RData"))
  
  # Save summary
  ga_summary <- data.frame(
    method = "GA",
    best_coverage = best_fitness,
    best_coverage_pct = 100 * best_fitness / sum(population),
    runtime_sec = runtime,
    runtime_min = runtime / 60
  )
  
  write.csv(ga_summary, 
            file = file.path(ga_dir, "ga_summary.csv"),
            row.names = FALSE)
  
  return(ga_results)
}

# ---------------------------------------------------------
# RUN GA
# ---------------------------------------------------------

ga_result <- run_ga_mclp_single(evalBudget = 250, seed = 123)

