# Test problem for a constrained combinatorial Bayesian optimisation
# PRBOCS-VB method is used to solve the problem here.
# Approaching the constraints by PRing them before sampling and adding gumbel trick
# GA (genetic algorithm) for comparison and true answer

# MINIMISATION PROBELM
# This code is for the CONTAMINATION test problem in BOCS
# 
# 10 instances with lambda = 0.1, 1, 10,100

# Author: Niyati Seth
# Date:  Decmber 2025

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
# Contamination Problem
# ---------------------------------------------------------

Contamination <- function(x, runlength, seed) {
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
  epsilon <- rep(0.35, n)   # error probability
  p <- rep(0.2, n)          # proportion limit
  cost <- rep(1, n)         # cost for prevention at stage i
  
  # Beta parameters for initial contamination, contamination rate, restoration rate
  initialAlpha <- 1
  initialBeta <- 30
  contamAlpha <- 1
  contamBeta <- 17 / 3
  restoreAlpha <- 5
  restoreBeta <- 3 / 7
  
  
  # RNG
  set.seed(as.integer(seed))
  # 
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
evalBudget <-250
n_init <- 10
order <- 2


# ---------------------------------------------------------
# SINGLE EXPERIMENT RUN
# ---------------------------------------------------------
run_single_experiment <- function(instance_id, n_vars, lambda) {
  
  cat("Running instance =", instance_id, "\n")
  
  seed <- instance_id + 100
  set.seed(seed)
  
  # ---------------------------------------------------------
  # INITIAL SAMPLES (shared across all three conditions)
  # ---------------------------------------------------------
  x_vals    <- sample_models(n_init, n_vars)
  num_inputs <- nrow(x_vals)
  y_vals    <- numeric(num_inputs)
  feasible  <- logical(num_inputs)
  
  for (i in seq_len(num_inputs)) {
    res        <- Contamination(x_vals[i, ], 100, seed = seed)
    y_vals[i]  <- res$fn
    feasible[i] <- all(res$constraint > res$limit)
  }
  
  n_init_actual <- nrow(x_vals)
  n_iter        <- evalBudget - n_init_actual
  
  # Helper: build initial SVB data structures from x_vals / y_vals
  build_svb <- function(xv, yv) {
    xin       <- order_effects(xv, order)$xTrain_in
    df        <- data.frame(y = yv, xin)
    vb_d      <- df[, -1]
    dup       <- which(duplicated(as.list(vb_d)))
    red       <- vb_d[, !duplicated(as.list(vb_d))]
    mdl       <- svb.fit(X = as.matrix(red), Y = yv,
                         family = "linear", slab = "laplace", intercept = TRUE)
    list(data = df, vb_data = vb_d, dup = dup, reduced = red, model = mdl)
  }
  
  # Helper: refit SVB after appending one new row
  refit_svb <- function(df_old, y_new, x_new_vec) {
    x_new_in  <- order_effects(matrix(x_new_vec, nrow = 1), order)$xTrain_in
    df_new    <- rbind(df_old, data.frame(y = y_new, x_new_in))
    vb_d      <- df_new[, -1]
    dup       <- which(duplicated(as.list(vb_d)))
    red       <- vb_d[, !duplicated(as.list(vb_d))]
    mdl       <- svb.fit(X = as.matrix(red), Y = df_new[, 1],
                         family = "linear", slab = "laplace", intercept = TRUE)
    list(data = df_new, vb_data = vb_d, dup = dup, reduced = red, model = mdl)
  }
  
  # Helper: extract posterior mean as probability vector (length n_vars)
  extract_mean_prob <- function(vb_model, vb_data, dup) {
    full_mu <- numeric(ncol(vb_data))
    kept    <- setdiff(seq_along(full_mu), dup)
    full_mu[kept] <- vb_model$mu
    for (col in dup) {
      dup_vals <- vb_data[, col]
      orig     <- which(apply(vb_data, 2, function(x) all(x == dup_vals)) &
                          !(seq_along(vb_data) %in% dup))
      if (length(orig) == 1) full_mu[col] <- full_mu[orig]
    }
    # first n_vars entries are main-effect coefficients; sigmoid → probabilities
    1 / (1 + exp(-full_mu[seq_len(n_vars)]))
  }
  
  # Helper: greedy constraint repair given a score vector
  greedy_repair <- function(score, seed_val) {
    idx   <- order(score, decreasing = TRUE)
    x_new <- rep(0, n_vars)
    for (k in seq_along(idx)) {
      ct <- Contamination(x_new, 100, seed = seed_val)
      if (ct$constraint[idx[k]] < ct$limit[idx[k]])
        x_new[idx[k]] <- 1
    }
    x_new
  }
  
  # ---------------------------------------------------------
  # INITIALISE ALL THREE CONDITIONS FROM THE SAME STARTING DATA
  # ---------------------------------------------------------
  s_full <- build_svb(x_vals, y_vals)
  s_ab1  <- build_svb(x_vals, y_vals)   # Ablation 1: Thompson only, no Gumbel
  s_ab2  <- build_svb(x_vals, y_vals)   # Ablation 2: Gumbel + E[theta], no Thompson
  
  # Shared data_history (feasibility fixed-length from the start)
  data_feasible_full <- rep(NA, evalBudget)
  data_feasible_ab1  <- rep(NA, evalBudget)
  data_feasible_ab2  <- rep(NA, evalBudget)
  data_feasible_full[1:n_init_actual] <- as.logical(feasible)
  data_feasible_ab1 [1:n_init_actual] <- as.logical(feasible)
  data_feasible_ab2 [1:n_init_actual] <- as.logical(feasible)
  
  data_history_full <- data.frame(y = y_vals, x_vals,
                                  feasibility = data_feasible_full[1:n_init_actual])
  data_history_ab1  <- data_history_full
  data_history_ab2  <- data_history_full
  colnames(data_history_full) <- c("y", paste0("X", 1:n_vars), "Feasibility")
  colnames(data_history_ab1)  <- colnames(data_history_full)
  colnames(data_history_ab2)  <- colnames(data_history_full)
  
  theta_full <- rep(0.5, n_vars)
  theta_ab1  <- rep(0.5, n_vars)
  
  y_history_full <- numeric(n_iter)
  y_history_ab1  <- numeric(n_iter)
  y_history_ab2  <- numeric(n_iter)
  
  # ---------------------------------------------------------
  # MAIN LOOP — all three conditions run each iteration
  # ---------------------------------------------------------
  for (t in 1:n_iter) {
    cat(sprintf("Instance %d | lambda %.1f | iter %d\n", instance_id, lambda, t))
    
    # ── FULL METHOD (Thompson + Gumbel) ──────────────────────
    stat_model_full <- function(theta)
      thompson_sam_svb(theta, s_full$model, s_full$dup, s_full$vb_data, order)
    
    opt_full   <- optim(theta_full, stat_model_full, method = "L-BFGS-B",
                        lower = 1e-8, upper = 0.999,
                        control = list(fnscale = -1))
    theta_full <- opt_full$par
    
    g_full         <- -log(-log(runif(n_vars)))
    x_tmp_full     <- rbinom(n_vars, 1, theta_full)
    cr_full        <- Contamination(x_tmp_full, 100, seed)
    viol_full      <- pmax(0, cr_full$limit - cr_full$constraint)
    score_full     <- theta_full + g_full - lambda * sum(viol_full)
    x_new_full     <- greedy_repair(score_full, seed)
    res_full       <- Contamination(x_new_full, 100, seed)
    y_new_full     <- res_full$fn
    feas_full      <- all(res_full$constraint > res_full$limit)
    
    data_feasible_full[t + n_init_actual] <- as.logical(feas_full)
    y_history_full[t] <- y_new_full
    dh_full <- data.frame(y = y_new_full,
                          matrix(x_new_full, nrow = 1),
                          feasibility = feas_full)
    colnames(dh_full) <- colnames(data_history_full)
    data_history_full <- rbind(data_history_full, dh_full)
    s_full <- refit_svb(s_full$data, y_new_full, x_new_full)
    
    # ── ABLATION 1: Thompson only, no Gumbel ─────────────────
    # g_j = 0; only randomness is the posterior sample theta_ab1
    stat_model_ab1 <- function(theta)
      thompson_sam_svb(theta, s_ab1$model, s_ab1$dup, s_ab1$vb_data, order)
    
    opt_ab1  <- optim(theta_ab1, stat_model_ab1, method = "L-BFGS-B",
                      lower = 1e-8, upper = 0.999,
                      control = list(fnscale = -1))
    theta_ab1 <- opt_ab1$par
    
    # no Gumbel draw — score uses theta directly
    x_tmp_ab1  <- rbinom(n_vars, 1, theta_ab1)
    cr_ab1     <- Contamination(x_tmp_ab1, 100, seed)
    viol_ab1   <- pmax(0, cr_ab1$limit - cr_ab1$constraint)
    score_ab1  <- theta_ab1 - lambda * sum(viol_ab1)   # g_j = 0
    x_new_ab1  <- greedy_repair(score_ab1, seed)
    res_ab1    <- Contamination(x_new_ab1, 100, seed)
    y_new_ab1  <- res_ab1$fn
    feas_ab1   <- all(res_ab1$constraint > res_ab1$limit)
    
    data_feasible_ab1[t + n_init_actual] <- as.logical(feas_ab1)
    y_history_ab1[t] <- y_new_ab1
    dh_ab1 <- data.frame(y = y_new_ab1,
                         matrix(x_new_ab1, nrow = 1),
                         feasibility = feas_ab1)
    colnames(dh_ab1) <- colnames(data_history_ab1)
    data_history_ab1 <- rbind(data_history_ab1, dh_ab1)
    s_ab1 <- refit_svb(s_ab1$data, y_new_ab1, x_new_ab1)
    
    # ── ABLATION 2: Gumbel + E[theta_j], no Thompson ─────────
    # theta is fixed at posterior mean mu_j; all randomness from Gumbel
    mu_prob_ab2 <- extract_mean_prob(s_ab2$model, s_ab2$vb_data, s_ab2$dup)
    
    g_ab2      <- -log(-log(runif(n_vars)))
    x_tmp_ab2  <- rbinom(n_vars, 1, mu_prob_ab2)
    cr_ab2     <- Contamination(x_tmp_ab2, 100, seed)
    viol_ab2   <- pmax(0, cr_ab2$limit - cr_ab2$constraint)
    score_ab2  <- mu_prob_ab2 + g_ab2 - lambda * sum(viol_ab2)
    x_new_ab2  <- greedy_repair(score_ab2, seed)
    res_ab2    <- Contamination(x_new_ab2, 100, seed)
    y_new_ab2  <- res_ab2$fn
    feas_ab2   <- all(res_ab2$constraint > res_ab2$limit)
    
    data_feasible_ab2[t + n_init_actual] <- as.logical(feas_ab2)
    y_history_ab2[t] <- y_new_ab2
    dh_ab2 <- data.frame(y = y_new_ab2,
                         matrix(x_new_ab2, nrow = 1),
                         feasibility = feas_ab2)
    colnames(dh_ab2) <- colnames(data_history_ab2)
    data_history_ab2 <- rbind(data_history_ab2, dh_ab2)
    s_ab2 <- refit_svb(s_ab2$data, y_new_ab2, x_new_ab2)
  }
  
  # ---------------------------------------------------------
  # RETURN all three conditions
  # ---------------------------------------------------------
  list(
    instance  = instance_id,
    lambda    = lambda,
    n_init    = n_init_actual,
    # Full method
    full = list(
      data       = data_history_full,
      y_history  = y_history_full,
      feasible   = data_feasible_full
    ),
    # Ablation 1: Thompson only
    ablation_no_gumbel = list(
      data       = data_history_ab1,
      y_history  = y_history_ab1,
      feasible   = data_feasible_ab1
    ),
    # Ablation 2: Gumbel + posterior mean
    ablation_gumbel_mean = list(
      data       = data_history_ab2,
      y_history  = y_history_ab2,
      feasible   = data_feasible_ab2
    )
  )
}

# ---------------------------------------------------------
# MAIN WRAPPER FOR ALL RUNS
# ---------------------------------------------------------
run_all_experiments <- function() {
  
  n_instances <- 10
  lambda_grid <- c(0.1, 1, 10, 100)
  
  base_dir <- "~/Projects:Codes/VaR-CoCoBO/Contamination/results_corrections"
  dir.create(base_dir, showWarnings = FALSE, recursive = TRUE)
  
  for (lam in lambda_grid) {
    
    cat("Running experiments for lambda =", lam, "\n")
    
    lambda_dir <- file.path(base_dir, paste0("lambda_", lam))
    dir.create(lambda_dir, showWarnings = FALSE)
    
    for (inst in 1:n_instances) {
      
      cat("  Instance =", inst, "\n")
      
      t_start <- Sys.time()
      
      res <- run_single_experiment(
        instance_id = inst,
        n_vars = n_vars,
        lambda = lam
      )
      
      t_end <- Sys.time()
      res$runtime_sec <- as.numeric(difftime(t_end, t_start, units = "secs"))
      
      save(
        res,
        file = file.path(lambda_dir, paste0("instance_", inst, ".RData"))
      )
      
      cat("  Saved instance", inst, "for lambda", lam, "\n")
    }
  }
}


# Run everything
run_all_experiments()


