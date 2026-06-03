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
  stopifnot(length(x) == n_vars)
  
  n_facilities <- sum(x)
  coverage_vector <- pmin(A %*% x, 1)
  total_covered_pop <- sum(coverage_vector * population)
  constraint <- n_facilities - max_facilities
  
  return(list(
    fn = total_covered_pop,  
    coverage = total_covered_pop,
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
# SINGLE EXPERIMENT RUN
# ---------------------------------------------------------

run_single_experiment <- function(instance_id, evalBudget, n_init, lambda, order) {
  
  cat("\n==============================================\n")
  cat("Running Instance", instance_id, "\n")
  cat("==============================================\n")
  
  # Set seed for reproducibility
  seed <- instance_id + 100
  set.seed(seed)
  
  # ---------------------------------------------------------
  # INITIAL SAMPLES
  # ---------------------------------------------------------
  
  x_vals <- sample_models(n_init, n_vars)
  y_vals <- mclp_prob(x_vals)
  
  cat("Initial solutions:\n")
  print(x_vals)
  cat("Number of facilities in each initial solution:", apply(x_vals, 1, sum), "\n")
  cat("\nInitial objectives (population covered):\n")
  print(y_vals)
  
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
  
  # ---------------------------------------------------------
  # OPTIMIZATION LOOP
  # ---------------------------------------------------------
  
  for (t in 1:evalBudget) {
    print(paste("VaRCBO Iteration", t))
    
    stat_model <- function(theta) {
      -thompson_sam_svb(theta, vb_model = vb_model, duplicate_cols, vb_data, order)
    }
    
    max_acq <- optim(theta_current, stat_model, method = 'L-BFGS-B', 
                     lower = 1e-8, upper = 0.999)
    expected_val <- max_acq$par
    theta_current <- expected_val
    
    # x_theta_current <- rbinom(n_vars, 1, theta_current)
    # cat("X before constraint", x_theta_current, "\n")
    # 
    # Gumbel noise
    g <- -log(-log(runif(n_vars)))
    # 
    # # CONSTRAINT PENALTY φᵢ
    # # For cardinality constraint: sum(x) = max_facilities
    # # Penalize facilities that would violate this constraint
    # 
    # # Expected number of facilities if we sample from theta
    # expected_facilities <- sum(theta_current)
    # 
    # # Compute penalty for each facility
    # phi <- numeric(n_vars)
    # 
    # if (expected_facilities > max_facilities) {
    #   # We're over budget - penalize all facilities proportionally
    #   # Facilities with higher theta contribute more to violation
    #   for (j in seq_len(n_vars)) {
    #     phi[j] <- theta_current[j] * (expected_facilities - max_facilities)
    #   }
    # } else {
    #   # We're under budget - no penalty (or small penalty to encourage reaching max)
    #   for (j in seq_len(n_vars)) {
    #     phi[j] <- 0
    #     # Optional: small penalty for being under-budget
    #     # phi[j] <- -theta_current[j] * (max_facilities - expected_facilities)
    #   }
    # }
    # 
    # # LINE 5: Compute scores with PENALTY (subtract phi)
    # score <- theta_current + g - lambda * phi
    # 
    # # LINE 6: Ranking
    # idx <- order(score, decreasing = TRUE)
    # 
    # # LINES 7-14: Construction with constraint checking
    # x_new <- rep(0, n_vars)
    # 
    # for (h in 1:n_vars) {
    #   facility_idx <- idx[h]
    #   x_new[facility_idx] <- 1
    #   
    #   if (sum(x_new) <= max_facilities) {
    #     # Constraint satisfied, keep it
    #   } else {
    #     # Constraint violated, reject it
    #     x_new[facility_idx] <- 0
    #   }
    #   
    #   if (sum(x_new) == max_facilities) {
    #     break
    #   }
    # }
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

    cat("X sampled from theta:", x_theta_current, "\n")

    # SCORING STEP: prioritize facilities that cover UNCOVERED population
    score <- numeric(n_vars)

    # What's uncovered in the sampled solution?
    coverage_sample <- pmin(as.vector(A %*% x_theta_current), 1)
    uncovered_pop <- population * (1 - coverage_sample)

    for (j in seq_len(n_vars)) {
      # How much uncovered population does facility j cover?
      covers_uncovered <- sum(A[, j] * uncovered_pop)

      # Normalize: as fraction of total population (0 to 1 scale)
      normalized_coverage <- covers_uncovered / sum(population)

      # Score: theta + gumbel + reward for covering uncovered areas
      score[j] <- theta_current[j] + g[j] + lambda * normalized_coverage
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
    
    cat(sprintf("Iteration %d: x = %s, coverage = %.0f\n", t, paste(x_new, collapse=""), coverage_new))
    
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
  
  cat(sprintf("\nInstance %d completed in %.2f seconds (%.2f minutes)\n", 
              instance_id, optimization_time, optimization_time/60))
  
  # ---------------------------------------------------------
  # RETURN RESULTS
  # ---------------------------------------------------------
  
  result_list <- list(
    instance_id = instance_id,
    seed = seed,
    x_history = x_history,
    y_history = y_history,
    coverage_history = coverage_history,
    prbocs_vb_result = prbocs_vb_result,
    data_history = data_history,
    theta_current = theta_current,
    optimization_time = optimization_time,
    best_coverage = max(coverage_history),
    best_iteration = which.max(coverage_history)
  )
  
  return(result_list)
}

# ---------------------------------------------------------
# MAIN WRAPPER FOR ALL INSTANCES
# ---------------------------------------------------------

run_all_experiments <- function() {
  
  # Parameters
  n_instances <- 20
  evalBudget <- 250
  n_init <- 5
  lambda <- 1
  order <- 2
  
  # Create results directory
  base_dir <- "results_MCLP/multiple_instances"
  dir.create(base_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Store all results
  all_results <- list()
  
  for (inst in 1:n_instances) {
    
    cat("\n##################################################\n")
    cat("STARTING INSTANCE", inst, "of", n_instances, "\n")
    cat("##################################################\n")
    
    # Run single experiment
    result <- run_single_experiment(
      instance_id = inst,
      evalBudget = evalBudget,
      n_init = n_init,
      lambda = lambda,
      order = order
    )
    
    # Save individual instance
    save(
      result,
      file = file.path(base_dir, paste0("instance_", inst, ".RData"))
    )
    
    # Store in list
    all_results[[inst]] <- result
    
    cat("Saved instance", inst, "\n")
    cat("Best coverage:", result$best_coverage, "at iteration", result$best_iteration, "\n")
  }
  
  # ---------------------------------------------------------
  # AGGREGATE SUMMARY ACROSS INSTANCES
  # ---------------------------------------------------------
  
  summary_df <- data.frame(
    instance = 1:n_instances,
    best_coverage = sapply(all_results, function(x) x$best_coverage),
    best_iteration = sapply(all_results, function(x) x$best_iteration),
    runtime_sec = sapply(all_results, function(x) x$optimization_time)
  )
  
  summary_df$best_coverage_pct <- 100 * summary_df$best_coverage / sum(population)
  
  cat("\n==============================================\n")
  cat("SUMMARY ACROSS ALL INSTANCES\n")
  cat("==============================================\n")
  print(summary_df)
  
  cat("\nMean best coverage:", mean(summary_df$best_coverage), 
      "(", mean(summary_df$best_coverage_pct), "%)\n")
  cat("Std best coverage:", sd(summary_df$best_coverage), "\n")
  cat("Mean runtime:", mean(summary_df$runtime_sec), "seconds\n")
  
  write.csv(summary_df, 
            file = file.path(base_dir, "summary_all_instances.csv"),
            row.names = FALSE)
  
  # Save all results
  save(all_results, 
       file = file.path(base_dir, "all_instances.RData"))
  
  return(all_results)
}

# ---------------------------------------------------------
# RUN ALL EXPERIMENTS
# ---------------------------------------------------------

all_results <- run_all_experiments() 


# ---------------------------------------------------------
# PLOTS & TABLES
# ---------------------------------------------------------

library(xtable)

# Read data
summary <- read.csv(file = "P3Compute/Facility Location/results_MCLP/multiple_instances/summary_all_instances.csv")

# Convert runtime to minutes
summary$runtime_min <- summary$runtime_sec / 60

# Reorder columns for better presentation
summary_table <- summary[, c("instance", "best_coverage", "best_coverage_pct", 
                             "best_iteration", "runtime_min")]

# Rename columns for better display
colnames(summary_table) <- c("Instance", "Best Coverage", "Coverage \\%", 
                             "Best Iteration", "Runtime (min)")

# Create xtable
latex_table <- xtable(summary_table, 
                      caption = "MCLP Optimization Results Across 10 Instances",
                      label = "tab:mclp_results",
                      digits = c(0, 0, 0, 2, 0, 2))  # Control decimal places

# Print LaTeX code
print(latex_table, 
      include.rownames = FALSE,
      caption.placement = "top",
      booktabs = TRUE,
      sanitize.text.function = function(x){x})  # Preserve LaTeX commands like \%


############################################################
# GA FOR MCLP — INDEPENDENT VERIFICATION + EVAL COUNTER
############################################################

library(GA)

cat("\n==============================================\n")
cat("Running GA for MCLP (Independent Check)\n")
cat("==============================================\n")

# ---------------------------------------------------------
# GLOBAL EVALUATION COUNTER
# ---------------------------------------------------------
ga_eval_count <<- 0

# ---------------------------------------------------------
# GA FITNESS FUNCTION
# ---------------------------------------------------------
ga_mclp_fitness <- function(x) {
  
  # Count evaluations
  ga_eval_count <<- ga_eval_count + 1
  
  # Enforce binary chromosome
  x <- round(x)
  
  # Evaluate the true MCLP model
  res <- mclp(x)
  
  # Penalize violating max_facilities
  penalty <- 0
  if (res$constraint > 0) {
    penalty <- -1e6 * res$constraint   # strong penalty for > max facilities
  }
  
  # MAXIMISE covered population
  return(res$coverage + penalty)
}

# ---------------------------------------------------------
# RUN GA (default settings except nBits)
# ---------------------------------------------------------
ga_start <- Sys.time()

ga_mclp <- ga(
  type    = "binary",
  fitness = ga_mclp_fitness,
  nBits   = n_vars
  # run     = 100         # early stopping if no improvement for 40 generations
)

ga_end <- Sys.time()
ga_time_sec <- as.numeric(difftime(ga_end, ga_start, units = "secs"))

# ---------------------------------------------------------
# RESULTS
# ---------------------------------------------------------
ga_sol <- round(ga_mclp@solution[1, ])
ga_obj <- mclp(ga_sol)$coverage

cat("\n========== GA RESULTS ==========\n")
cat("Best GA coverage: ", ga_obj, "\n")
cat("Facilities selected (0/1):\n")
print(ga_sol)
cat("Total facilities: ", sum(ga_sol), "\n")

cat("\n========== GA EVALUATION COUNT ==========\n")
cat("Number of objective evaluations: ", ga_eval_count, "\n")
cat("GA runtime (seconds): ", ga_time_sec, "\n")
cat("GA runtime (minutes): ", ga_time_sec / 60, "\n")


# ---------------------------------------------------------
# GA FOR MCLP — 20 INSTANCES, CONSOLE OUTPUT ONLY
# ---------------------------------------------------------

library(GA)

# Parameters (matching BO setup)
n_instances  <- 20
evalBudget   <- 250

all_coverage <- numeric(n_instances)
all_evals    <- numeric(n_instances)

for (inst in seq_len(n_instances)) {
  
  seed <- inst + 100
  set.seed(seed)
  
  eval_count <- 0L
  
  fitness_fn <- function(x) {
    eval_count <<- eval_count + 1L
    x   <- round(x)
    res <- mclp(x)
    penalty <- if (res$constraint > 0) -1e6 * res$constraint else 0
    res$coverage + penalty
  }
  
  ga_res <- ga(
    type    = "binary",
    fitness = fitness_fn,
    nBits   = n_vars,
    monitor = FALSE
  )
  
  best_sol      <- round(ga_res@solution[1, ])
  best_coverage <- mclp(best_sol)$coverage
  best_pct      <- 100 * best_coverage / sum(population)
  
  all_coverage[inst] <- best_coverage
  all_evals[inst]    <- eval_count
  
  cat(sprintf("Instance %2d | Coverage: %7.0f (%5.2f%%) | Evals: %d | Facilities: %s\n",
              inst, best_coverage, best_pct, eval_count, paste(best_sol, collapse = "")))
}

cat("\n--- SUMMARY ---\n")
cat(sprintf("Mean coverage : %.0f (%.2f%%)\n", mean(all_coverage), 100 * mean(all_coverage) / sum(population)))
cat(sprintf("Std  coverage : %.0f\n", sd(all_coverage)))
cat(sprintf("Min  coverage : %.0f\n", min(all_coverage)))
cat(sprintf("Max  coverage : %.0f\n", max(all_coverage)))
cat(sprintf("Mean evals    : %.1f\n", mean(all_evals)))

