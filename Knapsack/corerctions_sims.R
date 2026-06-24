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
# KNAPSACK PROBLEM SETUP
# ---------------------------------------------------------

n_vars <- 24

weights <- c(382745, 799601, 909247, 729069, 467902,  44328,
             34610, 698150, 823460, 903959, 853665, 551830,
             610856, 670702, 488960, 951111, 323046, 446298,
             931161,  31385, 496951, 264724, 224916, 169684)

values  <- c( 825594, 1677009, 1676628, 1523970,  943972,   97426,
              69666, 1296457, 1679693, 1902996, 1844992, 1049289,
              1252836, 1319836,  953277, 2067538,  675367,  853655,
              1826027,   65731,  901489,  577243,  466257,  369261)

W       <- 6404180
optimal <- 13549094

# ---------------------------------------------------------
# KNAPSACK EVALUATION
# ---------------------------------------------------------

knapsack <- function(x) {
  fn               <- sum(values * x)
  total_constraint <- sum(weights * x)
  limit            <- W
  feasible         <- total_constraint <= limit
  list(fn               = fn,
       total_constraint = total_constraint,
       limit            = limit,
       feasible         = feasible)
}

# ---------------------------------------------------------
# ORDER EFFECTS / THETA INTERACTION
# ---------------------------------------------------------

order_effects <- function(xTrain, order) {
  n_samp   <- nrow(xTrain)
  n_vars   <- ncol(xTrain)
  xTrain_in <- xTrain
  for (ord_i in 2:order) {
    offdProd <- combn(n_vars, ord_i)
    x_comb   <- array(dim = c(n_samp, ncol(offdProd)))
    for (j in 1:ncol(offdProd))
      x_comb[, j] <- apply(xTrain[, offdProd[, j], drop = FALSE], 1, prod)
    xTrain_in <- cbind(xTrain_in, x_comb)
  }
  list(xTrain_in = xTrain_in, combos = offdProd)
}

theta_interaction <- function(theta, order, n_vars) {
  theta    <- matrix(theta, nrow = 1, ncol = n_vars)
  theta_in <- theta
  for (ord_i in 2:order) {
    offdProd  <- combn(n_vars, ord_i)
    theta_comb <- array(dim = c(1, ncol(offdProd)))
    for (j in 1:ncol(offdProd))
      theta_comb[, j] <- apply(theta[, offdProd[, j], drop = FALSE], 1, prod)
    theta_in <- cbind(theta_in, theta_comb)
  }
  theta_in
}

# ---------------------------------------------------------
# ACQUISITION FUNCTION
# ---------------------------------------------------------

thompson_sam_svb <- function(theta_current, vb_model, duplicate_cols, vb_data, order) {
  full_mu   <- numeric(ncol(vb_data))
  kept_cols <- setdiff(seq_along(full_mu), duplicate_cols)
  full_mu[kept_cols] <- vb_model$mu
  for (col in duplicate_cols) {
    dup_column   <- vb_data[, col]
    original_col <- which(
      apply(vb_data, 2, function(x) all(x == dup_column)) &
        !(seq_along(vb_data) %in% duplicate_cols)
    )
    if (length(original_col) == 1) full_mu[col] <- full_mu[original_col]
  }
  theta_current_in <- theta_interaction(theta_current, order, n_vars)
  theta_current_in <- c(1, theta_current_in)
  coeffs           <- c(vb_model$intercept, full_mu)
  sum(theta_current_in * coeffs)
}

# ---------------------------------------------------------
# SHARED HELPERS
# ---------------------------------------------------------

build_svb <- function(xv, yv, order) {
  xin  <- order_effects(xv, order)$xTrain_in
  df   <- data.frame(y = yv, xin)
  vb_d <- df[, -1]
  dup  <- which(duplicated(as.list(vb_d)))
  red  <- vb_d[, !duplicated(as.list(vb_d))]
  mdl  <- svb.fit(X = as.matrix(red), Y = yv,
                  family = "linear", slab = "laplace", intercept = TRUE)
  list(data = df, vb_data = vb_d, dup = dup, reduced = red, model = mdl)
}

refit_svb <- function(df_old, y_new, x_new_vec, order) {
  x_new_in <- order_effects(matrix(x_new_vec, nrow = 1), order)$xTrain_in
  df_new   <- rbind(df_old, data.frame(y = y_new, x_new_in))
  vb_d     <- df_new[, -1]
  dup      <- which(duplicated(as.list(vb_d)))
  red      <- vb_d[, !duplicated(as.list(vb_d))]
  mdl      <- svb.fit(X = as.matrix(red), Y = df_new[, 1],
                      family = "linear", slab = "laplace", intercept = TRUE)
  list(data = df_new, vb_data = vb_d, dup = dup, reduced = red, model = mdl)
}

extract_mean_prob <- function(vb_model, vb_data, dup, n_vars) {
  full_mu          <- numeric(ncol(vb_data))
  kept             <- setdiff(seq_along(full_mu), dup)
  full_mu[kept]    <- vb_model$mu
  for (col in dup) {
    dup_vals <- vb_data[, col]
    orig     <- which(
      apply(vb_data, 2, function(x) all(x == dup_vals)) &
        !(seq_along(vb_data) %in% dup)
    )
    if (length(orig) == 1) full_mu[col] <- full_mu[orig]
  }
  1 / (1 + exp(-full_mu[seq_len(n_vars)]))
}

# Greedy knapsack selection given a score vector
# Ranks items by score and adds them greedily while capacity allows
greedy_knapsack <- function(score) {
  idx   <- order(score, decreasing = TRUE)
  x_new <- rep(0, n_vars)
  cum_w <- 0
  for (k in seq_along(idx)) {
    j <- idx[k]
    if (cum_w + weights[j] <= W) {
      x_new[j] <- 1
      cum_w    <- cum_w + weights[j]
    }
  }
  x_new
}

# phi_j for knapsack: cumulative weight overshoot penalty
# computed on a candidate x before greedy selection
knapsack_phi <- function(x_candidate) {
  phi   <- numeric(n_vars)
  cum_w <- 0
  for (j in seq_len(n_vars)) {
    if (x_candidate[j] == 1) {
      cum_w  <- cum_w + weights[j]
      phi[j] <- max(0, cum_w - W)
    }
  }
  phi
}

# ---------------------------------------------------------
# SINGLE EXPERIMENT RUN
# ---------------------------------------------------------

run_single_experiment <- function(n_init, instance_id, n_vars) {
  
  cat("Running n_init =", n_init, ", instance =", instance_id, "\n")
  
  evalBudget <- 200
  lambda     <- 10
  order      <- 2
  seed       <- instance_id + 100
  set.seed(seed)
  
  # ---------------------------------------------------------
  # INITIAL SAMPLES — shared across all conditions
  # ---------------------------------------------------------
  x_vals    <- sample_models(n_init, n_vars)
  y_vals    <- sapply(seq_len(nrow(x_vals)),
                      function(i) knapsack(x_vals[i, ])$fn)
  feasible  <- sapply(seq_len(nrow(x_vals)),
                      function(i) knapsack(x_vals[i, ])$feasible)
  
  n_init_actual <- nrow(x_vals)
  n_iter        <- evalBudget - n_init_actual
  
  # Build initial SVB — shared starting point for all surrogate conditions
  s_full <- build_svb(x_vals, y_vals, order)
  s_ab1  <- build_svb(x_vals, y_vals, order)   # Thompson only
  s_ab2  <- build_svb(x_vals, y_vals, order)   # Gumbel + E[theta]
  
  # Initialise data histories
  make_history <- function() {
    df <- data.frame(y = y_vals, x_vals,
                     Feasibility = feasible)
    colnames(df) <- c("y", paste0("X", 1:n_vars), "Feasibility")
    df
  }
  
  data_history_full   <- make_history()
  data_history_ab1    <- make_history()
  data_history_ab2    <- make_history()
  data_history_greedy <- make_history()
  
  y_history_full   <- numeric(n_iter)
  y_history_ab1    <- numeric(n_iter)
  y_history_ab2    <- numeric(n_iter)
  y_history_greedy <- numeric(n_iter)
  
  feas_full   <- rep(NA, n_iter)
  feas_ab1    <- rep(NA, n_iter)
  feas_ab2    <- rep(NA, n_iter)
  feas_greedy <- rep(NA, n_iter)
  
  theta_full <- rep(0.5, n_vars)
  theta_ab1  <- rep(0.5, n_vars)
  
  # Pure greedy score: fixed value-to-weight ratio, computed once
  # This never changes — greedy has no surrogate and no noise
  vw_ratio    <- values / weights
  x_greedy    <- greedy_knapsack(vw_ratio)
  res_greedy  <- knapsack(x_greedy)
  
  # ---------------------------------------------------------
  # MAIN LOOP
  # ---------------------------------------------------------
  
  for (t in 1:n_iter) {
    cat(sprintf("Instance %d | n_init %d | iter %d\n",
                instance_id, n_init, t))
    
    # ── FULL METHOD: Thompson + Gumbel ───────────────────────
    
    stat_model_full <- function(theta)
      thompson_sam_svb(theta, s_full$model, s_full$dup, s_full$vb_data, order)
    
    opt_full   <- optim(theta_full, stat_model_full, method = "L-BFGS-B",
                        lower = 1e-8, upper = 0.999,
                        control = list(fnscale = -1))
    theta_full <- opt_full$par
    
    g_full        <- -log(-log(runif(n_vars)))         # Gumbel draw — new each iteration
    x_tmp_full    <- rbinom(n_vars, 1, theta_full)
    phi_full      <- knapsack_phi(x_tmp_full)
    score_full    <- theta_full + g_full + lambda * phi_full
    x_new_full    <- greedy_knapsack(score_full)
    res_full      <- knapsack(x_new_full)
    y_new_full    <- res_full$fn
    feas_full[t]  <- res_full$feasible
    
    y_history_full[t] <- y_new_full
    dh <- data.frame(y = y_new_full, matrix(x_new_full, nrow = 1),
                     Feasibility = feas_full[t])
    colnames(dh) <- colnames(data_history_full)
    data_history_full <- rbind(data_history_full, dh)
    s_full <- refit_svb(s_full$data, y_new_full, x_new_full, order)
    
    # ── ABLATION 1: Thompson only, no Gumbel ─────────────────
    # g_j = 0; score = theta_hat + lambda * phi only
    # demonstrates Thompson sampling alone is insufficient
    
    stat_model_ab1 <- function(theta)
      thompson_sam_svb(theta, s_ab1$model, s_ab1$dup, s_ab1$vb_data, order)
    
    opt_ab1   <- optim(theta_ab1, stat_model_ab1, method = "L-BFGS-B",
                       lower = 1e-8, upper = 0.999,
                       control = list(fnscale = -1))
    theta_ab1 <- opt_ab1$par
    
    x_tmp_ab1   <- rbinom(n_vars, 1, theta_ab1)
    phi_ab1     <- knapsack_phi(x_tmp_ab1)
    score_ab1   <- theta_ab1 + lambda * phi_ab1        # g_j = 0
    x_new_ab1   <- greedy_knapsack(score_ab1)
    res_ab1     <- knapsack(x_new_ab1)
    y_new_ab1   <- res_ab1$fn
    feas_ab1[t] <- res_ab1$feasible
    
    y_history_ab1[t] <- y_new_ab1
    dh <- data.frame(y = y_new_ab1, matrix(x_new_ab1, nrow = 1),
                     Feasibility = feas_ab1[t])
    colnames(dh) <- colnames(data_history_ab1)
    data_history_ab1 <- rbind(data_history_ab1, dh)
    s_ab1 <- refit_svb(s_ab1$data, y_new_ab1, x_new_ab1, order)
    
    # ── ABLATION 2: Gumbel + E[theta], no Thompson ───────────
    # theta fixed at posterior mean; all randomness from Gumbel only
    # demonstrates Gumbel contribution independently of posterior uncertainty
    
    mu_prob_ab2 <- extract_mean_prob(s_ab2$model, s_ab2$vb_data,
                                     s_ab2$dup, n_vars)
    
    g_ab2       <- -log(-log(runif(n_vars)))            # Gumbel draw — new each iteration
    x_tmp_ab2   <- rbinom(n_vars, 1, mu_prob_ab2)
    phi_ab2     <- knapsack_phi(x_tmp_ab2)
    score_ab2   <- mu_prob_ab2 + g_ab2 + lambda * phi_ab2
    x_new_ab2   <- greedy_knapsack(score_ab2)
    res_ab2     <- knapsack(x_new_ab2)
    y_new_ab2   <- res_ab2$fn
    feas_ab2[t] <- res_ab2$feasible
    
    y_history_ab2[t] <- y_new_ab2
    dh <- data.frame(y = y_new_ab2, matrix(x_new_ab2, nrow = 1),
                     Feasibility = feas_ab2[t])
    colnames(dh) <- colnames(data_history_ab2)
    data_history_ab2 <- rbind(data_history_ab2, dh)
    s_ab2 <- refit_svb(s_ab2$data, y_new_ab2, x_new_ab2, order)
    
    # ── PURE GREEDY: value-to-weight ratio, no surrogate ─────
    # fixed ranking — same solution every iteration
    # this is the degenerate special case with no perturbations
    # included as practical external baseline per reviewer correction 1
    
    y_history_greedy[t] <- res_greedy$fn
    feas_greedy[t]      <- res_greedy$feasible
    dh <- data.frame(y = res_greedy$fn,
                     matrix(x_greedy, nrow = 1),
                     Feasibility = res_greedy$feasible)
    colnames(dh) <- colnames(data_history_greedy)
    data_history_greedy <- rbind(data_history_greedy, dh)
  }
  
  # ---------------------------------------------------------
  # RETURN all four conditions
  # ---------------------------------------------------------
  list(
    instance = instance_id,
    n_init   = n_init_actual,
    # Full method
    full = list(
      data      = data_history_full,
      y_history = y_history_full,
      feasible  = feas_full
    ),
    # Ablation 1: Thompson only
    ablation_no_gumbel = list(
      data      = data_history_ab1,
      y_history = y_history_ab1,
      feasible  = feas_ab1
    ),
    # Ablation 2: Gumbel + posterior mean
    ablation_gumbel_mean = list(
      data      = data_history_ab2,
      y_history = y_history_ab2,
      feasible  = feas_ab2
    ),
    # Pure greedy: value-to-weight ratio baseline
    greedy = list(
      data      = data_history_greedy,
      y_history = y_history_greedy,
      feasible  = feas_greedy
    )
  )
}

# ---------------------------------------------------------
# MAIN WRAPPER
# ---------------------------------------------------------

run_all_experiments <- function() {
  
  n_inis      <- c(5, 10, 20)
  n_instances <- 10
  
  base_dir <- "~/Projects:Codes/VaR-CoCoBO/knapsack/results_corrections"
  dir.create(base_dir, showWarnings = FALSE, recursive = TRUE)
  
  for (n_init in n_inis) {
    
    cat("Running experiments for n_init =", n_init, "\n")
    
    init_dir <- file.path(base_dir, paste0("n_init_", n_init))
    dir.create(init_dir, showWarnings = FALSE)
    
    for (inst in 1:n_instances) {
      
      cat("  Instance =", inst, "\n")
      
      t_start <- Sys.time()
      
      res <- run_single_experiment(
        n_init      = n_init,
        instance_id = inst,
        n_vars      = n_vars
      )
      
      t_end          <- Sys.time()
      res$runtime_sec <- as.numeric(difftime(t_end, t_start, units = "secs"))
      
      save(res,
           file = file.path(init_dir,
                            paste0("instance_", inst, ".RData")))
      
      cat("  Saved instance", inst, "for n_init", n_init, "\n")
    }
  }
}

run_all_experiments()