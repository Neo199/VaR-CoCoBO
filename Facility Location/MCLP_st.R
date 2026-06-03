# Application for a constrained combinatorial Bayesian optimisation - STRESS TEST
# PRBOCS-VB method compared with GLPK (optimal) and GA
# Multiple perturbations of initial population

# Author: Niyati Seth
# Date:  January 2026
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
library(ROI.plugin.glpk)
library(bayesplot)
library(sparsevb)
library(selectiveInference)
library(maxcovr)

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

mclp <- function(x, pop = population) {
  # pop parameter allows using perturbed population data
  stopifnot(length(x) == n_vars)
  
  n_facilities <- sum(x)
  coverage_vector <- pmin(A %*% x, 1)
  total_covered_pop <- sum(coverage_vector * pop)
  constraint <- n_facilities - max_facilities
  
  return(list(
    fn = total_covered_pop,  
    coverage = total_covered_pop,
    constraint = constraint,
    n_facilities = n_facilities
  ))
}

mclp_prob <- function(x_mat, pop = population) {
  num_inputs <- nrow(x_mat)
  cost <- numeric(num_inputs)
  
  for (i in seq_len(num_inputs)) {
    res <- mclp(x_mat[i, ], pop = pop)
    cost[i] <- res$fn
  }
  return(cost)
}

# ---------------------------------------------------------
# GLPK OPTIMAL SOLVER
# ---------------------------------------------------------

solve_mclp_optimal <- function(pop = population) {
  
  cat("\n==============================================\n")
  cat("SOLVING MCLP OPTIMALLY WITH GLPK\n")
  cat("==============================================\n")
  
  start_time <- Sys.time()
  
  model <- MIPModel() %>%
    # Decision variables: x[j] = 1 if facility j is opened
    add_variable(x[j], j = 1:n_vars, type = "binary") %>%
    # Auxiliary variables: y[i] = 1 if demand point i is covered
    add_variable(y[i], i = 1:n_demand, type = "binary") %>%
    # Maximize covered population
    set_objective(sum_expr(pop[i] * y[i], i = 1:n_demand), "max") %>%
    # Constraint: at most max_facilities can be opened
    add_constraint(sum_expr(x[j], j = 1:n_vars) <= max_facilities) %>%
    # Linking constraint: y[i] can only be 1 if at least one covering facility is open
    add_constraint(y[i] <= sum_expr(A[i, j] * x[j], j = 1:n_vars), i = 1:n_demand)
  
  result <- solve_model(model, with_ROI(solver = "glpk", verbose = TRUE))
  
  end_time <- Sys.time()
  solve_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  # Extract solution
  x_optimal <- get_solution(result, x[j])$value
  optimal_coverage <- objective_value(result)
  
  cat("\nOptimal solution found!\n")
  cat("Facilities opened:", which(x_optimal == 1), "\n")
  cat("Optimal coverage:", optimal_coverage, "\n")
  cat("Solve time:", solve_time, "seconds\n")
  
  return(list(
    x_optimal = x_optimal,
    optimal_coverage = optimal_coverage,
    solve_time = solve_time
  ))
}

# ---------------------------------------------------------
# GENETIC ALGORITHM
# ---------------------------------------------------------

run_ga <- function(seed, pop = population, maxiter = 250, popSize = 50) {
  
  cat("\n==============================================\n")
  cat("RUNNING GENETIC ALGORITHM\n")
  cat("==============================================\n")
  
  set.seed(seed)
  
  # Fitness function for GA (must return single value to maximize)
  fitness_ga <- function(x) {
    # Penalize if constraint violated
    if (sum(x) != max_facilities) {
      return(-1e9)  # Large penalty
    }
    return(mclp(x, pop = pop)$coverage)
  }
  
  start_time <- Sys.time()
  
  ga_result <- ga(
    type = "binary",
    fitness = fitness_ga,
    nBits = n_vars,
    maxiter = maxiter,
    popSize = popSize,
    pcrossover = 0.8,
    pmutation = 0.1,
    elitism = max(1, round(popSize * 0.05)),
    seed = seed,
    monitor = FALSE
  )
  
  end_time <- Sys.time()
  ga_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  best_solution <- as.vector(ga_result@solution[1, ])
  best_fitness <- ga_result@fitnessValue
  
  cat("\nGA completed!\n")
  cat("Best solution:", best_solution, "\n")
  cat("Best fitness:", best_fitness, "\n")
  cat("GA time:", ga_time, "seconds\n")
  
  return(list(
    x_best = best_solution,
    best_coverage = best_fitness,
    ga_time = ga_time,
    convergence_history = ga_result@summary[, "max"]
  ))
}

# ---------------------------------------------------------
# PERTURBATION FUNCTION
# ---------------------------------------------------------

perturb_demand_data <- function(base_population, perturbation_rate = 0.1, seed) {
  # Perturb the demand (population) at each census tract
  # This simulates uncertainty in population estimates
  
  set.seed(seed)
  
  perturbed_population <- base_population
  
  # Add multiplicative noise to each demand point
  # Each population value is multiplied by (1 + noise)
  # where noise ~ Uniform[-perturbation_rate, +perturbation_rate]
  
  for (i in 1:length(base_population)) {
    noise <- runif(1, min = -perturbation_rate, max = perturbation_rate)
    perturbed_population[i] <- base_population[i] * (1 + noise)
    
    # Ensure non-negative population
    perturbed_population[i] <- max(0, perturbed_population[i])
  }
  
  # Round to integers (population counts)
  perturbed_population <- round(perturbed_population)
  
  cat("\nOriginal total population:", sum(base_population), "\n")
  cat("Perturbed total population:", sum(perturbed_population), "\n")
  cat("Population change:", 
      sprintf("%.2f%%", 100 * (sum(perturbed_population) - sum(base_population)) / sum(base_population)), "\n")
  
  return(perturbed_population)
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
# SINGLE BO EXPERIMENT RUN
# ---------------------------------------------------------

run_single_bo_experiment <- function(instance_id, init_population, pop, evalBudget, lambda, order) {
  
  cat("\n--- Running BO Instance", instance_id, "---\n")
  
  # Set seed for reproducibility
  seed <- instance_id + 1000
  set.seed(seed)
  
  # Use provided initial population
  x_vals <- init_population
  n_init <- nrow(x_vals)
  y_vals <- mclp_prob(x_vals, pop = pop)
  
  cat("Initial objectives (population covered):", y_vals, "\n")
  
  # ---------------------------------------------------------
  # OPTIMIZATION SETUP
  # ---------------------------------------------------------
  
  start_time <- Sys.time()
  
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
    res <- mclp(x_vals[i,], pop = pop)
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
  
  # ---------------------------------------------------------
  # OPTIMIZATION LOOP
  # ---------------------------------------------------------
  
  for (t in 1:n_iter) {
    
    stat_model <- function(theta) {
      -thompson_sam_svb(theta, vb_model = vb_model, duplicate_cols, vb_data, order)
    }
    
    max_acq <- optim(theta_current, stat_model, method = 'L-BFGS-B', 
                     lower = 1e-8, upper = 0.999)
    expected_val <- max_acq$par
    theta_current <- expected_val
    
    # Gumbel noise
    g <- -log(-log(runif(n_vars)))
    
    # Sample from theta to understand uncovered areas
    x_theta_current <- rbinom(n_vars, 1, theta_current)

    # Adjust to exactly max_facilities
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

    # SCORING STEP: prioritize facilities that cover UNCOVERED population
    score <- numeric(n_vars)

    # What's uncovered in the sampled solution?
    coverage_sample <- pmin(as.vector(A %*% x_theta_current), 1)
    uncovered_pop <- pop * (1 - coverage_sample)

    for (j in seq_len(n_vars)) {
      # How much uncovered population does facility j cover?
      covers_uncovered <- sum(A[, j] * uncovered_pop)

      # Normalize: as fraction of total population (0 to 1 scale)
      normalized_coverage <- covers_uncovered / sum(pop)

      # Score: theta + gumbel + reward for covering uncovered areas
      score[j] <- theta_current[j] + g[j] + lambda * normalized_coverage
    }

    # RANKING STEP
    idx <- order(score, decreasing = TRUE)

    # CONSTRUCTION: select top max_facilities
    x_new <- rep(0, n_vars)
    x_new[idx[1:max_facilities]] <- 1
    
    # Evaluate
    res <- mclp(x_new, pop = pop)
    y_new <- res$fn
    coverage_new <- res$coverage
    
    iter_idx <- n_init + t
    data_hnew <- data.frame(y = y_new, t(x_new))
    colnames(data_hnew) <- colnames(data_history)
    data_history <- rbind(data_history, data_hnew)
    
    x_history[[iter_idx]] <- x_new
    y_history[iter_idx] <- y_new
    coverage_history[iter_idx] <- coverage_new
    
    # Update model
    x_new_matrix <- matrix(x_new, nrow = 1, ncol = n_vars)
    x_new_in_comb <- order_effects(x_new_matrix, order)
    x_new_in <- x_new_in_comb$xTrain_in
    
    data_new <- data.frame(y = y_new, x_new_in)
    data <- rbind(data, data_new)
    
    prbocs_vb_result[iter_idx, ] <- expected_val
    
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
  
  # ---------------------------------------------------------
  # EXTRACT BEST RESULTS (EXCLUDING INITIAL SAMPLES)
  # ---------------------------------------------------------
  
  # Coverage from BO iterations only (exclude initial samples)
  bo_coverage_only <- coverage_history[(n_init + 1):evalBudget]
  
  # Overall best (including initial samples for comparison)
  best_coverage_overall <- max(coverage_history)
  best_iteration_overall <- which.max(coverage_history)
  
  # Best from BO iterations only
  best_coverage_bo <- max(bo_coverage_only)
  best_iteration_bo <- which.max(bo_coverage_only) + n_init  # Adjust index
  
  cat(sprintf("  Best overall: %.0f (iteration %d)\n", best_coverage_overall, best_iteration_overall))
  cat(sprintf("  Best from BO: %.0f (iteration %d)\n", best_coverage_bo, best_iteration_bo))
  
  # ---------------------------------------------------------
  # RETURN RESULTS
  # ---------------------------------------------------------
  
  result_list <- list(
    instance_id = instance_id,
    seed = seed,
    coverage_history = coverage_history,
    coverage_history_bo_only = bo_coverage_only,
    optimization_time = optimization_time,
    best_coverage = best_coverage_bo,  # Use BO-only best
    best_coverage_overall = best_coverage_overall,  # For reference
    best_iteration = best_iteration_bo,
    best_iteration_overall = best_iteration_overall,
    n_init = n_init
  )
  
  return(result_list)
}

# ---------------------------------------------------------
# STRESS TEST: SINGLE PERTURBATION TRIAL
# ---------------------------------------------------------

run_single_perturbation_trial <- function(trial_id, base_population, base_demand, perturbation_rate = 0.1) {
  
  cat("\n##################################################\n")
  cat("PERTURBATION TRIAL", trial_id, "\n")
  cat("##################################################\n")
  
  # Generate perturbed DEMAND data
  perturbed_seed <- trial_id * 100
  perturbed_demand <- perturb_demand_data(base_demand, perturbation_rate, perturbed_seed)
  
  # Use the same initial solutions (not perturbed)
  init_pop <- base_population
  
  cat("\nInitial solutions:\n")
  print(init_pop)
  cat("Initial objectives with perturbed demand:", mclp_prob(init_pop, pop = perturbed_demand), "\n")
  
  # Parameters
  evalBudget <- 250
  lambda <- 1
  order <- 2
  n_bo_instances <- 10
  
  # ---------------------------------------------------------
  # 1. RUN GLPK (once per perturbation)
  # ---------------------------------------------------------
  
  glpk_result <- solve_mclp_optimal(pop = perturbed_demand)
  
  # ---------------------------------------------------------
  # 2. RUN GA (once per perturbation)
  # ---------------------------------------------------------
  
  ga_seed <- trial_id * 1000
  ga_result <- run_ga(seed = ga_seed, pop = perturbed_demand, maxiter = evalBudget, popSize = 50)
  
  # ---------------------------------------------------------
  # 3. RUN BO (10 instances per perturbation)
  # ---------------------------------------------------------
  
  bo_results <- list()
  
  for (bo_inst in 1:n_bo_instances) {
    bo_instance_id <- trial_id * 100 + bo_inst
    
    bo_result <- run_single_bo_experiment(
      instance_id = bo_instance_id,
      init_population = init_pop,
      pop = perturbed_demand,
      evalBudget = evalBudget,
      lambda = lambda,
      order = order
    )
    
    bo_results[[bo_inst]] <- bo_result
    
    cat(sprintf("BO instance %d/%d: Best coverage = %.0f\n", 
                bo_inst, n_bo_instances, bo_result$best_coverage))
  }
  
  # ---------------------------------------------------------
  # COMPILE TRIAL RESULTS
  # ---------------------------------------------------------
  
  trial_summary <- data.frame(
    trial_id = trial_id,
    glpk_coverage = glpk_result$optimal_coverage,
    glpk_time = glpk_result$solve_time,
    ga_coverage = ga_result$best_coverage,
    ga_time = ga_result$ga_time,
    bo_mean_coverage = mean(sapply(bo_results, function(x) x$best_coverage)),
    bo_sd_coverage = sd(sapply(bo_results, function(x) x$best_coverage)),
    bo_min_coverage = min(sapply(bo_results, function(x) x$best_coverage)),
    bo_max_coverage = max(sapply(bo_results, function(x) x$best_coverage)),
    bo_mean_time = mean(sapply(bo_results, function(x) x$optimization_time))
  )
  
  # Calculate optimality gaps
  trial_summary$ga_gap <- 100 * (glpk_result$optimal_coverage - ga_result$best_coverage) / glpk_result$optimal_coverage
  trial_summary$bo_gap_mean <- 100 * (glpk_result$optimal_coverage - trial_summary$bo_mean_coverage) / glpk_result$optimal_coverage
  
  return(list(
    trial_summary = trial_summary,
    glpk_result = glpk_result,
    ga_result = ga_result,
    bo_results = bo_results,
    perturbed_demand = perturbed_demand
  ))
}

# ---------------------------------------------------------
# MAIN STRESS TEST WRAPPER
# ---------------------------------------------------------

run_stress_test <- function() {
  
  cat("\n##################################################\n")
  cat("STARTING STRESS TEST\n")
  cat("##################################################\n")
  
  # Parameters
  n_trials <- 5
  n_init <- 5
  perturbation_rate <- 0.1  # 10% perturbation in demand
  
  # Generate base initial population (same for all trials)
  set.seed(42)
  base_population <- sample_models(n_init, n_vars)
  
  # Store original demand (will be perturbed for each trial)
  base_demand <- population
  
  cat("\nBase initial population (same across all trials):\n")
  print(base_population)
  cat("Base objectives (with original demand):", mclp_prob(base_population, pop = base_demand), "\n")
  cat("Original total demand:", sum(base_demand), "\n")
  
  # Create results directory
  base_dir <- "results_MCLP/stress_test"
  dir.create(base_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Store all trial results
  all_trial_results <- list()
  all_trial_summaries <- list()
  
  # ---------------------------------------------------------
  # RUN ALL TRIALS (each with different perturbed demand)
  # ---------------------------------------------------------
  
  for (trial in 1:n_trials) {
    
    trial_result <- run_single_perturbation_trial(
      trial_id = trial,
      base_population = base_population,
      base_demand = base_demand,
      perturbation_rate = perturbation_rate
    )
    
    all_trial_results[[trial]] <- trial_result
    all_trial_summaries[[trial]] <- trial_result$trial_summary
    
    # Save individual trial
    save(
      trial_result,
      file = file.path(base_dir, paste0("trial_", trial, ".RData"))
    )
    
    cat("\nTrial", trial, "completed!\n")
    print(trial_result$trial_summary)
  }
  
  # ---------------------------------------------------------
  # AGGREGATE RESULTS ACROSS ALL TRIALS
  # ---------------------------------------------------------
  
  summary_df <- do.call(rbind, all_trial_summaries)
  
  cat("\n==============================================\n")
  cat("STRESS TEST SUMMARY ACROSS ALL TRIALS\n")
  cat("==============================================\n")
  print(summary_df)
  
  # Overall statistics
  overall_stats <- data.frame(
    metric = c("GLPK Coverage", "GA Coverage", "BO Mean Coverage",
               "GA Gap (%)", "BO Gap (%)"),
    mean = c(
      mean(summary_df$glpk_coverage),
      mean(summary_df$ga_coverage),
      mean(summary_df$bo_mean_coverage),
      mean(summary_df$ga_gap),
      mean(summary_df$bo_gap_mean)
    ),
    sd = c(
      sd(summary_df$glpk_coverage),
      sd(summary_df$ga_coverage),
      sd(summary_df$bo_mean_coverage),
      sd(summary_df$ga_gap),
      sd(summary_df$bo_gap_mean)
    )
  )
  
  cat("\nOverall Statistics:\n")
  print(overall_stats)
  
  # Save results
  write.csv(summary_df, 
            file = file.path(base_dir, "stress_test_summary.csv"),
            row.names = FALSE)
  
  write.csv(overall_stats,
            file = file.path(base_dir, "overall_statistics.csv"),
            row.names = FALSE)
  
  save(all_trial_results, summary_df, overall_stats, base_population, base_demand,
       file = file.path(base_dir, "stress_test_complete.RData"))
  
  # ---------------------------------------------------------
  # CREATE VISUALIZATION
  # ---------------------------------------------------------
  
  # Convergence plot comparing methods
  create_convergence_plot(all_trial_results, base_dir)
  
  # Boxplot comparison
  create_boxplot_comparison(summary_df, base_dir)
  
  return(list(
    all_trial_results = all_trial_results,
    summary_df = summary_df,
    overall_stats = overall_stats
  ))
}

# ---------------------------------------------------------
# VISUALIZATION FUNCTIONS
# ---------------------------------------------------------

create_convergence_plot <- function(all_trial_results, base_dir) {
  
  # Extract convergence data for first trial as example
  trial_1 <- all_trial_results[[1]]
  
  # Prepare data for plotting
  n_iter <- 250
  
  # BO convergence (take mean across 10 instances)
  bo_convergence <- matrix(0, nrow = 10, ncol = n_iter)
  for (i in 1:10) {
    bo_convergence[i, ] <- trial_1$bo_results[[i]]$coverage_history
  }
  bo_mean <- colMeans(bo_convergence)
  bo_max <- apply(bo_convergence, 2, max)
  bo_min <- apply(bo_convergence, 2, min)
  
  # GA convergence
  ga_convergence <- trial_1$ga_result$convergence_history
  
  # GLPK optimal (constant line)
  glpk_optimal <- rep(trial_1$glpk_result$optimal_coverage, n_iter)
  
  # Create plot
  pdf(file.path(base_dir, "convergence_comparison_trial1.pdf"), width = 10, height = 6)
  
  plot(1:n_iter, glpk_optimal, type = "l", col = "red", lwd = 2, lty = 2,
       xlab = "Iteration", ylab = "Population Covered",
       main = "Convergence Comparison: GLPK vs GA vs BO (Trial 1)",
       ylim = c(min(c(bo_min, ga_convergence)), max(glpk_optimal) * 1.05))
  
  lines(1:length(ga_convergence), ga_convergence, col = "blue", lwd = 2)
  lines(1:n_iter, bo_mean, col = "darkgreen", lwd = 2)
  
  # Add BO confidence band
  polygon(c(1:n_iter, rev(1:n_iter)), 
          c(bo_max, rev(bo_min)),
          col = rgb(0, 0.5, 0, 0.2), border = NA)
  
  legend("bottomright", 
         legend = c("GLPK (Optimal)", "GA", "BO (Mean)", "BO (Range)"),
         col = c("red", "blue", "darkgreen", rgb(0, 0.5, 0, 0.2)),
         lty = c(2, 1, 1, 1), lwd = c(2, 2, 2, 10),
         cex = 0.8)
  
  dev.off()
  
  cat("Convergence plot saved!\n")
}

create_boxplot_comparison <- function(summary_df, base_dir) {
  
  # Prepare data for boxplot
  plot_data <- data.frame(
    Coverage = c(summary_df$glpk_coverage, 
                 summary_df$ga_coverage,
                 summary_df$bo_mean_coverage),
    Method = rep(c("GLPK", "GA", "BO"), each = nrow(summary_df))
  )
  
  pdf(file.path(base_dir, "method_comparison_boxplot.pdf"), width = 8, height = 6)
  
  p <- ggplot(plot_data, aes(x = Method, y = Coverage, fill = Method)) +
    geom_boxplot() +
    geom_jitter(width = 0.1, alpha = 0.5) +
    labs(title = "Coverage Comparison Across Methods",
         subtitle = paste0(nrow(summary_df), " perturbation trials"),
         y = "Population Covered",
         x = "Method") +
    theme_minimal() +
    theme(legend.position = "none")
  
  print(p)
  
  dev.off()
  
  cat("Boxplot saved!\n")
}

# ---------------------------------------------------------
# RUN THE STRESS TEST
# ---------------------------------------------------------

stress_test_results <- run_stress_test()

cat("\n##################################################\n")
cat("STRESS TEST COMPLETED!\n")
cat("##################################################\n")


############################################################
# RE-RUN GA ONLY FOR EACH TRIAL — RECORD OBJ EVAL COUNTS
############################################################

library(GA)

ga_eval_results <- data.frame(
  trial_id = integer(),
  ga_evals = integer(),
  ga_best_coverage = numeric(),
  ga_time_sec = numeric()
)

n_trials          <- 5
perturbation_rate <- 0.1
n_instances       <- 10

for (trial_id in 1:n_trials) {
  
  cat("\n==============================================\n")
  cat("GA TRIAL", trial_id, "\n")
  cat("==============================================\n")
  
  # Same perturbed demand as stress test
  perturbed_demand <- perturb_demand_data(population, perturbation_rate, seed = trial_id * 100)
  
  ga_coverages <- numeric(n_instances)
  ga_evals     <- numeric(n_instances)
  
  for (ga_inst in 1:n_instances) {
    
    ga_seed <- trial_id * 1000 + ga_inst
    eval_count <- 0L
    
    fitness_ga <- function(x) {
      eval_count <<- eval_count + 1L
      if (sum(x) != max_facilities) return(-1e9)
      mclp(x, pop = perturbed_demand)$coverage
    }
    
    set.seed(ga_seed)
    
    ga_out <- ga(
      type       = "binary",
      fitness    = fitness_ga,
      nBits      = n_vars,
      monitor    = FALSE
    )
    
    ga_coverages[ga_inst] <- ga_out@fitnessValue
    ga_evals[ga_inst]     <- eval_count
    
    cat(sprintf("  Instance %2d | Coverage: %.0f | Evals: %d\n",
                ga_inst, ga_out@fitnessValue, eval_count))
  }
  
  cat(sprintf("Trial %d SUMMARY | Mean: %.0f (%.2f%%) | Std: %.2f%% | Mean evals: %.1f\n",
              trial_id,
              mean(ga_coverages),
              100 * mean(ga_coverages) / sum(population),
              100 * sd(ga_coverages) / sum(population),
              mean(ga_evals)))
}

