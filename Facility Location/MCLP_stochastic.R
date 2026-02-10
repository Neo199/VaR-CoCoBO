# Application for a constrained combinatorial Bayesian optimisation
# PRBOCS-VB method for STOCHASTIC MCLP
# Population values are drawn from Gaussian distributions (single sample per evaluation)

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
# STOCHASTIC MCLP DATA PREPARATION
# ---------------------------------------------------------

service_radius <- 5000
max_facilities <- 4
population_cv <- 0.15  # Coefficient of variation for population (15% noise)

# Mark which demand-facility pairs are within service radius
distance_matrix$covered <- as.integer(distance_matrix$distance <= service_radius)

# Create coverage matrix A
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

# Extract MEAN population for each demand point
demand_lookup <- demand_data %>%
  dplyr::select(NAME, POP2000) %>%
  rename(DestinationName = NAME, population_mean = POP2000)

population_df <- A_df %>%
  dplyr::select(DestinationName) %>%
  left_join(demand_lookup, by = "DestinationName")

population_mean <- population_df$population_mean
population_sd <- population_mean * population_cv  # SD = CV * mean

cat("Number of facilities:", n_vars, "\n")
cat("Number of demand points:", n_demand, "\n")
cat("Total mean population:", sum(population_mean), "\n")
cat("Max facilities allowed:", max_facilities, "\n")
cat("Population CV (noise level):", population_cv, "\n")

# ---------------------------------------------------------
# STOCHASTIC MCLP PROBLEM
# ---------------------------------------------------------

# Single noisy evaluation - sample population once
mclp_stochastic <- function(x) {
  stopifnot(length(x) == n_vars)
  
  n_facilities <- sum(x)
  constraint <- n_facilities - max_facilities
  
  # Sample population once (adds noise)
  population_sample <- pmax(0, rnorm(n_demand, mean = population_mean, sd = population_sd))
  
  # Evaluate coverage with this noisy population
  coverage_vector <- pmin(A %*% x, 1)
  total_coverage <- sum(coverage_vector * population_sample)
  
  return(list(
    fn = total_coverage,
    coverage = total_coverage,
    constraint = constraint,
    n_facilities = n_facilities
  ))
}

mclp_prob <- function(x_mat) {
  num_inputs <- nrow(x_mat)
  cost <- numeric(num_inputs)
  
  for (i in seq_len(num_inputs)) {
    res <- mclp_stochastic(x_mat[i, ])
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
  full_mu <- numeric(ncol(vb_data))
  kept_cols <- setdiff(seq_along(full_mu), duplicate_cols)
  full_mu[kept_cols] <- vb_model$mu
  
  for (col in duplicate_cols) {
    dup_column <- vb_data[, col]
    original_col <- which(apply(vb_data, 2, function(x) all(x == dup_column)) & 
                            !(seq_along(vb_data) %in% duplicate_cols))
    if (length(original_col) == 1) {
      full_mu[col] <- full_mu[original_col]
    }
  }
  
  theta_current_in <- theta_interaction(theta_current, order, n_vars)
  theta_current_in <- c(1, theta_current_in)
  
  coeffs <- c(vb_model$intercept, full_mu)
  y_pred <- sum(theta_current_in * coeffs)
  
  return(y_pred = y_pred)
}

# ---------------------------------------------------------
# SINGLE EXPERIMENT RUN (BO)
# ---------------------------------------------------------

run_single_bo_experiment <- function(instance_id, evalBudget, n_init, lambda, order) {
  
  cat("\n==============================================\n")
  cat("Running BO Instance", instance_id, "\n")
  cat("==============================================\n")
  
  seed <- instance_id + 100
  set.seed(seed)
  
  # Initial samples
  x_vals <- sample_models(n_init, n_vars)
  y_vals <- mclp_prob(x_vals)  # Each evaluation gets noisy population
  
  cat("Initial objectives (noisy coverage):\n")
  print(y_vals)
  
  start_time <- Sys.time()
  
  # Setup
  xTrain <- x_vals
  yTrain <- y_vals
  
  xTrain_in_comb <- order_effects(xTrain, order)
  xTrain_in <- xTrain_in_comb$xTrain_in
  
  data <- data.frame(y = y_vals, xTrain_in)
  
  theta_current <- rep(0.5, n_vars)
  x_history <- list()
  y_history <- numeric(evalBudget)
  data_history <- data.frame(y = y_vals, x_vals)
  
  # Store initial results
  for (i in 1:n_init) {
    y_history[i] <- y_vals[i]
  }
  
  prbocs_vb_result <- matrix(0, evalBudget, n_vars)
  
  vb_data <- data[, -1]
  duplicate_cols <- which(duplicated(as.list(vb_data)))
  data_reduced <- vb_data[, !duplicated(as.list(vb_data))]
  
  vb_model <- svb.fit(
    X = as.matrix(data_reduced),
    Y = data[, 1],
    family = "linear",
    slab = "laplace",
    intercept = TRUE
  )
  
  # Optimization loop
  for (t in 1:evalBudget) {
    print(paste("BO Iteration", t))
    
    stat_model <- function(theta) {
      -thompson_sam_svb(theta, vb_model = vb_model, duplicate_cols, vb_data, order)
    }
    
    max_acq <- optim(theta_current, stat_model, method = 'L-BFGS-B', 
                     lower = 1e-8, upper = 0.999)
    expected_val <- max_acq$par
    theta_current <- expected_val
    
    # Gumbel-based selection
    g <- -log(-log(runif(n_vars)))
    
    # Sample from theta
    x_theta_current <- rbinom(n_vars, 1, theta_current)
    
    # Adjust to max_facilities
    if (sum(x_theta_current) > max_facilities) {
      selected <- which(x_theta_current == 1)
      keep <- sample(selected, max_facilities)
      x_theta_current <- rep(0, n_vars)
      x_theta_current[keep] <- 1
    } else if (sum(x_theta_current) < max_facilities) {
      unselected <- which(x_theta_current == 0)
      add <- sample(unselected, max_facilities - sum(x_theta_current))
      x_theta_current[add] <- 1
    }
    
    # Scoring based on mean population (use mean for scoring, not noisy samples)
    score <- numeric(n_vars)
    coverage_sample <- pmin(as.vector(A %*% x_theta_current), 1)
    uncovered_pop <- population_mean * (1 - coverage_sample)
    
    for (j in seq_len(n_vars)) {
      covers_uncovered <- sum(A[, j] * uncovered_pop)
      normalized_coverage <- covers_uncovered / sum(population_mean)
      score[j] <- theta_current[j] + g[j] + lambda * normalized_coverage
    }
    
    # Ranking and construction
    idx <- order(score, decreasing = TRUE)
    x_new <- rep(0, n_vars)
    x_new[idx[1:max_facilities]] <- 1
    
    # Evaluate with noisy population (single sample)
    res <- mclp_stochastic(x_new)
    y_new <- res$fn
    
    x_history[[t]] <- x_new
    y_history[t] <- y_new
    
    data_hnew <- data.frame(y = y_new, t(x_new))
    colnames(data_hnew) <- colnames(data_history)
    data_history <- rbind(data_history, data_hnew)
    
    cat(sprintf("Iteration %d: Coverage = %.0f\n", t, y_new))
    
    # Update model
    x_new_matrix <- matrix(x_new, nrow = 1, ncol = n_vars)
    x_new_in <- order_effects(x_new_matrix, order)$xTrain_in
    
    data <- rbind(data, data.frame(y = y_new, x_new_in))
    prbocs_vb_result[t, ] <- expected_val
    
    vb_data <- data[, -1]
    duplicate_cols <- which(duplicated(as.list(vb_data)))
    data_reduced <- vb_data[, !duplicated(as.list(vb_data))]
    
    vb_model <- svb.fit(
      X = as.matrix(data_reduced),
      Y = data[, 1],
      family = "linear",
      slab = "laplace",
      intercept = TRUE
    )
  }
  
  end_time <- Sys.time()
  runtime <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  return(list(
    instance_id = instance_id,
    seed = seed,
    x_history = x_history,
    y_history = y_history,
    data_history = data_history,
    best_coverage = max(y_history),
    best_iteration = which.max(y_history),
    runtime_sec = runtime
  ))
}

# ---------------------------------------------------------
# SINGLE EXPERIMENT RUN (GA)
# ---------------------------------------------------------

run_single_ga_experiment <- function(instance_id, evalBudget) {
  
  cat("\n==============================================\n")
  cat("Running GA Instance", instance_id, "\n")
  cat("==============================================\n")
  
  seed <- instance_id + 100
  set.seed(seed)
  
  start_time <- Sys.time()
  
  # Fitness function (noisy evaluation)
  fitness_mclp <- function(x) {
    if (sum(x) != max_facilities) {
      return(-1e10)
    }
    res <- mclp_stochastic(x)  # Noisy evaluation
    return(res$fn)
  }
  
  # Run GA
  ga_result <- ga(
    type = "binary",
    fitness = fitness_mclp,
    nBits = n_vars,
    maxiter = evalBudget,
    popSize = 50,
    pcrossover = 0.8,
    pmutation = 0.1,
    elitism = max(1, round(50 * 0.05)),
    seed = seed,
    monitor = TRUE
  )
  
  end_time <- Sys.time()
  runtime <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  best_solution <- as.vector(ga_result@solution[1, ])
  best_fitness <- ga_result@fitnessValue
  fitness_history <- ga_result@summary[, "max"]
  
  cat(sprintf("GA completed: Best coverage = %.0f\n", best_fitness))
  
  return(list(
    instance_id = instance_id,
    seed = seed,
    best_solution = best_solution,
    best_coverage = best_fitness,
    convergence_history = fitness_history,
    runtime_sec = runtime
  ))
}

# ---------------------------------------------------------
# RUN ALL EXPERIMENTS
# ---------------------------------------------------------

run_all_stochastic_experiments <- function() {
  
  n_instances <- 10
  evalBudget <- 250
  n_init <- 5
  lambda <- 1
  order <- 2
  
  base_dir <- "results_MCLP_Stochastic"
  dir.create(base_dir, showWarnings = FALSE, recursive = TRUE)
  
  bo_results <- list()
  ga_results <- list()
  
  # Run BO experiments
  cat("\n########## RUNNING BO EXPERIMENTS ##########\n")
  for (inst in 1:n_instances) {
    result <- run_single_bo_experiment(
      instance_id = inst,
      evalBudget = evalBudget,
      n_init = n_init,
      lambda = lambda,
      order = order
    )
    
    bo_results[[inst]] <- result
    
    # Save individual instance - use temporary variable
    save(result, 
         file = file.path(base_dir, paste0("bo_instance_", inst, ".RData")))
    
    cat("Saved BO instance", inst, "\n")
  }
  
  # Run GA experiments
  cat("\n########## RUNNING GA EXPERIMENTS ##########\n")
  for (inst in 1:n_instances) {
    result <- run_single_ga_experiment(
      instance_id = inst,
      evalBudget = evalBudget
    )
    
    ga_results[[inst]] <- result
    
    # Save individual instance - use temporary variable
    save(result, 
         file = file.path(base_dir, paste0("ga_instance_", inst, ".RData")))
    
    cat("Saved GA instance", inst, "\n")
  }
  
  # Create summaries
  bo_summary <- data.frame(
    instance = 1:n_instances,
    best_coverage = sapply(bo_results, function(x) x$best_coverage),
    best_iteration = sapply(bo_results, function(x) x$best_iteration),
    runtime_sec = sapply(bo_results, function(x) x$runtime_sec)
  )
  
  ga_summary <- data.frame(
    instance = 1:n_instances,
    best_coverage = sapply(ga_results, function(x) x$best_coverage),
    runtime_sec = sapply(ga_results, function(x) x$runtime_sec)
  )
  
  bo_summary$runtime_min <- bo_summary$runtime_sec / 60
  ga_summary$runtime_min <- ga_summary$runtime_sec / 60
  
  write.csv(bo_summary, file.path(base_dir, "bo_summary.csv"), row.names = FALSE)
  write.csv(ga_summary, file.path(base_dir, "ga_summary.csv"), row.names = FALSE)
  
  # Comparison
  comparison <- data.frame(
    instance = 1:n_instances,
    BO_coverage = bo_summary$best_coverage,
    GA_coverage = ga_summary$best_coverage,
    Diff = bo_summary$best_coverage - ga_summary$best_coverage,
    BO_runtime_min = bo_summary$runtime_min,
    GA_runtime_min = ga_summary$runtime_min
  )
  
  write.csv(comparison, file.path(base_dir, "comparison.csv"), row.names = FALSE)
  
  cat("\n=== SUMMARY ===\n")
  cat("BO Mean Coverage:", mean(bo_summary$best_coverage), "± SD:", sd(bo_summary$best_coverage), "\n")
  cat("GA Mean Coverage:", mean(ga_summary$best_coverage), "± SD:", sd(ga_summary$best_coverage), "\n")
  cat("Mean Difference:", mean(comparison$Diff), "\n")
  
  print(comparison)
  
  # Save all results together
  save(bo_results, ga_results, 
       file = file.path(base_dir, "all_results.RData"))
  
  return(list(bo_results = bo_results, ga_results = ga_results))
}

# RUN EVERYTHING
results <- run_all_stochastic_experiments()