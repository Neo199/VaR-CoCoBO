
library(ggplot2)
library(dplyr)

num_instances <- 10
vd <- load("/Users/niyati/Projects:Codes/P3Compute/knapsack_gumbel_results/n_ini_5/instance_1.RData")

# Create a list to hold all instances
dt <- vector("list", num_instances)

data <- list()
time <- list()

# Loop over instance files
for (i in 1:num_instances) {
  file_path <- paste0("/Users/niyati/Projects:Codes/P3Compute/knapsack_gumbel_results/n_ini_5/instance_", i, ".RData")
  load(file_path)
  dt[[i]] <- res
}

for (i in 1:num_instances) {
  data[[i]] <- dt[[i]]$data
  time[[i]] <- dt[[i]]$runtime_sec
}

get_mean_df <- function(data) {
  num_rows <- nrow(data[[1]])
  y_means <- numeric(num_rows)
  y_lowers <- numeric(num_rows)
  y_uppers <- numeric(num_rows)
  y_max <- numeric(num_rows)
  
  y_vals <- sapply(data, function(df) df$y)
  for (i in 1:num_instances) {
    for (j in 1:num_rows) {
      y_means[j] <- rowMeans(y_vals)
      y_max[j] <- cummax(y_vals[,i])
      stderr <- sd(y_vals) / sqrt(length(y_vals))
      error_margin <- qt(0.975, df = length(y_vals) - 1) * stderr
      y_lowers[j] <- y_max[j] - error_margin
      y_uppers[j] <- y_max[j] + error_margin
    }
  }
  
  
  data.frame(
    row = 1:num_rows,
    mean_y = y_means,
    max_y = y_max,
    ci_lower = y_lowers,
    ci_upper = y_uppers
  )
}

df_trace <- data.frame(
  iter = 1:length(y_history),
  y = y_history,
  feasible = sapply(1:length(y_history), function(i) {
    x_i <- optim_result[i, ]
    knapsack_res <- knapsack(x_i)
    knapsack_res$total_constraint <= knapsack_res$limit
  })
)

df_trace$best_feas <- NA
best_so_far <- -Inf
for (i in 1:nrow(df_trace)) {
  if (df_trace$feasible[i]) {
    best_so_far <- max(best_so_far, df_trace$y[i])
  }
  df_trace$best_feas[i] <- best_so_far
}


ggplot(df_trace, aes(x = iter)) +
  geom_line(aes(y = best_feas), linewidth = 1, colour = "blue") +
  geom_point(aes(y = y, colour = feasible)) +
  scale_colour_manual(values = c("red", "darkgreen")) +
  labs(
    title = "PRBOCS-VB Optimisation Trace",
    y = "Objective Value",
    colour = "Feasible?"
  ) +
  theme_minimal()

# ---------------------------------------------------------
# Knapsack Problem
# ---------------------------------------------------------

knapsack <- function(x, seed = 10) {
  # x: binary vector (1 = facility installed)
  # runlength, seed: for reproducibility if stochastic
  n <- length(x)
  stopifnot(n == n_vars)
  
  fn <- sum(values * x)
  
  # Only include items that are taken (x == 1)
  constraint <- rep(0, n)
  for (k in 1:n) {
    constraint[k] <- weights[k] * x[k]
  }    # For penalty
  
  total_constraint <- sum(weights * x)  # same measure, but treated as constraint
  limit <- W
  
  return(list(fn = fn, constraint = constraint, limit = limit, total_constraint = total_constraint))
}


# Load all instances (each has $data and $runtime_sec)
load_all_instances <- function(n_ini, num_instances) {
  dt <- vector("list", num_instances)
  
  for (i in 1:num_instances) {
    path <- sprintf(
      "/Users/niyati/Projects:Codes/P3Compute/knapsack_gumbel_results/n_ini_%s/instance_%s.RData",
      n_ini, i
    )
    load(path)   # loads `res`
    
    # ---- FIX: compute feasibility column immediately ----
    df <- res$data
    df$feasible <- apply(df[, -1], 1, function(x) {
      out <- knapsack(x)
      out$total_constraint <= out$limit
    })
    res$data <- df
    # -----------------------------------------------------
    
    dt[[i]] <- res
  }
  dt
}

# Compute cumulative best feasible value for *one* instance
compute_best_feasible <- function(df) {
  best <- -Inf
  best_vec <- numeric(nrow(df))
  
  for (i in seq_len(nrow(df))) {
    if (df$feasible[i]) {
      best <- max(best, df$y[i])
    }
    best_vec[i] <- best
  }
  best_vec
}

summarise_instances <- function(dt) {
  num_instances <- length(dt)
  num_iter <- nrow(dt[[1]]$data)
  
  BF <- matrix(NA, nrow = num_iter, ncol = num_instances)
  
  for (i in 1:num_instances) {
    df <- dt[[i]]$data
    BF[, i] <- compute_best_feasible(df)
  }
  
  mean_bf <- rowMeans(BF)
  sd_bf   <- apply(BF, 1, sd)
  stderr  <- sd_bf / sqrt(num_instances)
  err     <- qt(0.975, df = num_instances - 1) * stderr
  
  summary_df <- data.frame(
    iter = 1:num_iter,
    mean = mean_bf,
    lower = mean_bf - err,
    upper = mean_bf + err
  )
  
  # ---- final best objective summary ----
  final_vals <- BF[num_iter, ]
  
  final_mean <- mean(final_vals)
  final_sd   <- sd(final_vals)
  final_se   <- final_sd / sqrt(num_instances)
  final_err  <- qt(0.975, df = num_instances - 1) * final_se
  
  runtime_vec <- sapply(dt, function(x) x$runtime_sec)
  
  summary_table <- data.frame(
    n_ini  = nrow(dt[[1]]$data[dt[[1]]$data$init == TRUE,]),
    best_obj_mean = final_mean,
    best_obj_lower = final_mean - final_err,
    best_obj_upper = final_mean + final_err,
    runtime_mean = mean(runtime_vec),
    runtime_sd   = sd(runtime_vec)
  )
  
  list(summary_df = summary_df, summary_table = summary_table)
}


all_results <- list()
all_summary_tables <- list()
n_inis <- c(5, 10, 20)

for (n in n_inis) {
  dt <- load_all_instances(n, num_instances = 10)
  sm <- summarise_instances(dt)
  
  all_results[[as.character(n)]] <- list(
    summary_df = sm$summary_df,
    summary_table = sm$summary_table,
    runtimes = sapply(dt, function(x) x$runtime_sec)
  )
  
  all_summary_tables[[as.character(n)]] <- sm$summary_table
}

final_table <- do.call(rbind, all_summary_tables)
final_table

#Add n_ini tag to summary df
tag_results <- function(all_results) {
  bind_rows(lapply(names(all_results), function(n) {
    df <- all_results[[n]]$summary_df
    df$n_ini <- as.integer(n)
    df
  }))
}
summary_long <- tag_results(all_results)
initial_counts <- c("5" = 5, "10" = 10, "20" = 20)



#---------PLOTS-----------------------------------------------
true_x <- as.integer(strsplit("110111000110100100000111", "")[[1]])
true_opt <- knapsack(true_x)$objective


plot_iter_vs_y <- function(dt_summary, df_first_instance, true_opt) {
  
  ggplot() +
    
    # Mean best-feasible trace
    geom_line(data = dt_summary,
              aes(x = iter, y = mean),
              linewidth = 1, colour = "blue") +
    
    geom_ribbon(data = dt_summary,
                aes(x = iter, ymin = lower, ymax = upper),
                alpha = 0.2, fill = "blue") +
    
    # All points from the first instance
    geom_point(data = df_first_instance,
               aes(x = iter, y = y, color = feasible),
               size = 2, alpha = 0.9) +
    
    # The true optimum
    geom_hline(yintercept = true_opt,
               linetype = "dashed",
               linewidth = 1,
               colour = "darkgreen") +
    
    scale_color_manual(values = c("FALSE" = "red", "TRUE" = "black"),
                       labels = c("Infeasible", "Feasible")) +
    
    theme_bw(base_size = 13) +
    labs(
      title = "Optimisation Trace (Iteration vs y)",
      x = "Iteration",
      y = "Objective Value y",
      color = "Feasible?",
      caption = sprintf("True optimum y = %.3f", true_opt)
    )
}
# x string -> vector
true_x <- as.integer(strsplit("110111000110100100000111", "")[[1]])
true_opt <- knapsack(true_x)$fn

# Load instance 1 for n_ini = 10 (example)
df1 <- load_all_instances(10, 10)[[1]]$data

# Add iteration index (needed for plotting)
df1$iter <- seq_len(nrow(df1))

# Plot
plot_iter_vs_y(all_results[["10"]]$summary_df, df1, true_opt)

#Combined comparison plot (all curves together)
plot_combined <- ggplot(summary_long,
                        aes(x = iter, y = mean, color = factor(n_ini))) +
  
  geom_line(linewidth = 1) +
  
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = factor(n_ini)),
              alpha = 0.15, colour = NA) +
  
  # Vertical lines for init sample region
  geom_vline(data = data.frame(n_ini = c(5,10,20),
                               x = c(5,10,20)),
             aes(xintercept = x, color = factor(n_ini)),
             linetype = "dotted", linewidth = 0.9) +
  
  # True optimum
  geom_hline(yintercept = true_opt,
             linetype = "dashed", linewidth = 1.1, color = "black") +
  
  scale_color_brewer(palette = "Dark2", name = "n_ini") +
  scale_fill_brewer(palette = "Dark2", name = "n_ini") +
  
  theme_bw(base_size = 14) +
  labs(
    title = "Comparison of Best-Feasible Traces for Different n_ini",
    x = "Iteration",
    y = "Best Feasible y",
    subtitle = sprintf("True optimum y = %.3f (dashed)", true_opt)
  )
ggsave(
  filename = "/Users/niyati/Projects:Codes/P3Compute/knapsack_gumbel_results/comparison.pdf",
  plot = plot_combined,
  device = cairo_pdf,
  width = 11, height = 6, units = "in",
  dpi = 1200, bg = "white"
)

#Multi-panel figure for n_ini
plot_faceted <- ggplot(summary_long,
                       aes(x = iter, y = mean)) +
  
  geom_line(color = "blue", linewidth = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              fill = "blue", alpha = 0.2) +
  
  # Vertical line separating initial samples
  geom_vline(aes(xintercept = n_ini),
             linetype = "dotted", linewidth = 1, colour = "black") +
  
  # True optimum
  geom_hline(yintercept = true_opt,
             linetype = "dashed", linewidth = 1, color = "darkgreen") +
  
  facet_wrap(~ n_ini, scales = "free_y") +
  
  theme_bw(base_size = 14) +
  labs(
    title = "Multi-Panel Comparison of Best-Feasible Traces",
    x = "Iteration",
    y = "Best Feasible y"
  )

#Multi-panel with initial points displayed
df_first_instances <- bind_rows(lapply(c(5,10,20), function(n) {
  df <- load_all_instances(n, 10)[[1]]$data
  df$iter <- seq_len(nrow(df))
  df$init <- df$iter <= n
  df$n_ini <- n
  df
}))

# create data frame for vertical lines
vline_df <- data.frame(
  n_ini = c(5, 10, 20),
  xintercept = c(5, 10, 20)
)

plot_faceted_with_initials <- ggplot() +
  
  # Best feasible traces
  geom_line(data = summary_long,
            aes(x = iter, y = mean),
            color = "blue", linewidth = 1) +
  
  geom_ribbon(data = summary_long,
              aes(x = iter, ymin = lower, ymax = upper),
              fill = "blue", alpha = 0.2) +
  
  # All points from first instance
  geom_point(data = df_first_instances,
             aes(x = iter, y = y, color = feasible, shape = init),
             size = 2) +
  
  scale_color_manual(values = c("FALSE" = "red", "TRUE" = "black")) +
  scale_shape_manual(values = c(`TRUE` = 17, `FALSE` = 16)) +
  
  # --- FIXED VERTICAL LINES ---
  geom_vline(data = vline_df,
             aes(xintercept = xintercept),
             linetype = "dotted", linewidth = 1, colour = "black") +
  
  # true optimum
  geom_hline(yintercept = true_opt,
             linetype = "dashed", linewidth = 1, color = "darkgreen") +
  
  facet_wrap(~ n_ini, scales = "free_y") +
  
  theme_bw(base_size = 14) +
  labs(
    title = "Optimisation Traces for Different Initial Samples",
    x = "Iteration",
    y = "Objective Value y",
    color = "Feasible",
    shape = "Initial?"
  )

plot_faceted_with_initials
ggsave(
  filename = "/Users/niyati/Projects:Codes/P3Compute/knapsack_gumbel_results/faceted_with_initials.pdf",
  plot = plot_faceted_with_initials,
  device = cairo_pdf,
  width = 11, height = 6, units = "in",
  dpi = 1200, bg = "white"
)

#---------------------TIME-----------------------------------

make_latex_table <- function(all_results) {
  tab <- bind_rows(lapply(names(all_results), function(n) {
    t <- all_results[[n]]$summary_table
    t$n_ini <- n
    t
  }))
  
  tab$runtime_min <- tab$runtime_mean / 60
  tab$runtime_sd_min <- tab$runtime_sd / 60
  
  tab2 <- tab %>%
    select( n_ini,
      best_obj_mean,
      best_obj_lower,
      best_obj_upper,
      runtime_min,
      runtime_sd_min
    )
  
  print(xtable::xtable(tab2, digits = 3), include.rownames = FALSE)
}

make_latex_table(all_results)

