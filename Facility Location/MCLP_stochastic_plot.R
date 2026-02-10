library(ggplot2)
library(dplyr)
library(tidyr)

# ---------------------------------------------------------
# OKABE-ITO COLOR PALETTE (colorblind-friendly)
# ---------------------------------------------------------

okabe_ito <- c(
  orange = "#E69F00",
  sky_blue = "#56B4E9",
  bluish_green = "#009E73",
  yellow = "#F0E442",
  blue = "#0072B2",
  vermillion = "#D55E00",
  reddish_purple = "#CC79A7",
  black = "#000000"
)

# For BO vs GA comparison - extract as unnamed values
bo_color <- unname(okabe_ito["blue"])
ga_color <- unname(okabe_ito["vermillion"])

# ---------------------------------------------------------
# LOAD RESULTS
# ---------------------------------------------------------

base_dir <- "results_MCLP_Stochastic"
load(file.path(base_dir, "all_results.RData"))

bo_summary <- read.csv(file.path(base_dir, "bo_summary.csv"))
ga_summary <- read.csv(file.path(base_dir, "ga_summary.csv"))

# ---------------------------------------------------------
# 1. CONVERGENCE PLOT: BO vs GA with uncertainty bands
# ---------------------------------------------------------

# Extract convergence histories for all BO instances
n_instances <- length(bo_results)
evalBudget <- length(bo_results[[1]]$y_history)

# For BO: use cummax of y_history
bo_convergence <- matrix(0, nrow = evalBudget, ncol = n_instances)
for (i in 1:n_instances) {
  bo_convergence[, i] <- cummax(bo_results[[i]]$y_history)
}

# For GA: need to filter out penalty values (-1e10) and track best VALID solution
ga_maxiter <- length(ga_results[[1]]$convergence_history)

ga_convergence <- matrix(0, nrow = ga_maxiter, ncol = n_instances)
for (i in 1:n_instances) {
  raw_history <- ga_results[[i]]$convergence_history
  
  # Filter out penalty values (anything less than 0 or suspiciously large negative)
  # Replace penalties with previous best valid value
  best_so_far <- 0
  cleaned_history <- numeric(length(raw_history))
  
  for (j in 1:length(raw_history)) {
    if (raw_history[j] > 0) {  # Valid solution
      best_so_far <- max(best_so_far, raw_history[j])
      cleaned_history[j] <- best_so_far
    } else {  # Penalty value
      cleaned_history[j] <- best_so_far
    }
  }
  
  ga_convergence[, i] <- cleaned_history
}

# Make sure both have same length
if (nrow(bo_convergence) != nrow(ga_convergence)) {
  cat("WARNING: BO and GA have different lengths!\n")
  min_length <- min(nrow(bo_convergence), nrow(ga_convergence))
  bo_convergence <- bo_convergence[1:min_length, ]
  ga_convergence <- ga_convergence[1:min_length, ]
  evalBudget <- min_length
}

# Calculate mean and standard deviation across instances
bo_mean <- rowMeans(bo_convergence)
bo_sd <- apply(bo_convergence, 1, sd)
bo_se <- bo_sd / sqrt(n_instances)

ga_mean <- rowMeans(ga_convergence)
ga_sd <- apply(ga_convergence, 1, sd)
ga_se <- ga_sd / sqrt(n_instances)

# Debug
cat("BO final mean:", tail(bo_mean, 1), "\n")
cat("GA final mean:", tail(ga_mean, 1), "\n")
cat("BO range:", range(bo_mean), "\n")
cat("GA range:", range(ga_mean), "\n")

# Create data frame for plotting
convergence_data <- data.frame(
  iteration = rep(1:evalBudget, 2),
  method = rep(c("BO", "GA"), each = evalBudget),
  mean_coverage = c(bo_mean, ga_mean),
  lower_sd = c(bo_mean - bo_sd, ga_mean - ga_sd),
  upper_sd = c(bo_mean + bo_sd, ga_mean + ga_sd),
  lower_se = c(bo_mean - 1.96 * bo_se, ga_mean - 1.96 * ga_se),
  upper_se = c(bo_mean + 1.96 * bo_se, ga_mean + 1.96 * ga_se)
)

# Plot with shaded uncertainty bands
p1 <- ggplot(convergence_data, aes(x = iteration, y = mean_coverage, color = method, fill = method)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = lower_se, ymax = upper_se), alpha = 0.2, color = NA) +
  labs(
    title = "Convergence: BO vs GA on Stochastic MCLP",
    subtitle = "Shaded area shows 95% confidence interval across 10 runs",
    x = "Iteration",
    y = "Best Coverage (Population)",
    color = "Method",
    fill = "Method"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(size = 12)
  ) +
  scale_color_manual(values = c("BO" = bo_color, "GA" = ga_color)) +
  scale_fill_manual(values = c("BO" = bo_color, "GA" = ga_color))

ggsave(file.path(base_dir, "stcohastic_convergence_comparison.pdf"), p1, width = 10, height = 6)
print(p1)

# ---------------------------------------------------------
# 2. INDIVIDUAL RUNS PLOT: Show all 10 trajectories
# ---------------------------------------------------------

all_bo_runs <- data.frame()
for (i in 1:n_instances) {
  run_data <- data.frame(
    iteration = 1:evalBudget,
    coverage = cummax(bo_results[[i]]$y_history)[1:evalBudget],
    instance = as.factor(i),
    method = "BO"
  )
  all_bo_runs <- rbind(all_bo_runs, run_data)
}

all_ga_runs <- data.frame()
for (i in 1:n_instances) {
  # Clean GA data same way
  raw_history <- ga_results[[i]]$convergence_history[1:evalBudget]
  best_so_far <- 0
  cleaned_history <- numeric(length(raw_history))
  
  for (j in 1:length(raw_history)) {
    if (raw_history[j] > 0) {
      best_so_far <- max(best_so_far, raw_history[j])
      cleaned_history[j] <- best_so_far
    } else {
      cleaned_history[j] <- best_so_far
    }
  }
  
  run_data <- data.frame(
    iteration = 1:evalBudget,
    coverage = cleaned_history,
    instance = as.factor(i),
    method = "GA"
  )
  all_ga_runs <- rbind(all_ga_runs, run_data)
}

all_runs <- rbind(all_bo_runs, all_ga_runs)

p2 <- ggplot(all_runs, aes(x = iteration, y = coverage, group = instance, color = method)) +
  geom_line(alpha = 0.4, linewidth = 0.5) +
  facet_wrap(~method, ncol = 2) +
  labs(
    title = "Individual Run Trajectories (10 instances each)",
    subtitle = "Shows variability across different initializations and stochastic evaluations",
    x = "Iteration",
    y = "Best Coverage (Population)"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 14),
    strip.text = element_text(face = "bold", size = 12)
  ) +
  scale_color_manual(values = c("BO" = bo_color, "GA" = ga_color))

ggsave(file.path(base_dir, "individual_runs.pdf"), p2, width = 12, height = 5)
print(p2)

# ---------------------------------------------------------
# 3. BOXPLOT: Final performance comparison
# ---------------------------------------------------------

final_comparison <- data.frame(
  coverage = c(bo_summary$best_coverage, ga_summary$best_coverage),
  method = rep(c("BO", "GA"), each = n_instances)
)

p3 <- ggplot(final_comparison, aes(x = method, y = coverage, fill = method)) +
  geom_boxplot(alpha = 0.7, width = 0.5) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 2) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 4, fill = "white") +
  labs(
    title = "Final Best Coverage: BO vs GA",
    subtitle = "Diamond shows mean, box shows quartiles, points show individual runs",
    x = "Method",
    y = "Best Coverage (Population)"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 14)
  ) +
  scale_fill_manual(values = c("BO" = bo_color, "GA" = ga_color))

ggsave(file.path(base_dir, "final_performance_boxplot.pdf"), p3, width = 7, height = 6)
print(p3)

# ---------------------------------------------------------
# 4. VARIANCE ACROSS ITERATIONS
# ---------------------------------------------------------

variance_data <- data.frame(
  iteration = 1:evalBudget,
  bo_variance = apply(bo_convergence, 1, var),
  ga_variance = apply(ga_convergence, 1, var)
) %>%
  pivot_longer(cols = c(bo_variance, ga_variance), 
               names_to = "method", 
               values_to = "variance") %>%
  mutate(method = ifelse(method == "bo_variance", "BO", "GA"))

p4 <- ggplot(variance_data, aes(x = iteration, y = variance, color = method)) +
  geom_line(linewidth = 1) +
  labs(
    title = "Variance in Performance Across Instances",
    subtitle = "Lower variance indicates more consistent performance despite stochasticity",
    x = "Iteration",
    y = "Variance in Coverage",
    color = "Method"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 14)
  ) +
  scale_color_manual(values = c("BO" = bo_color, "GA" = ga_color))

ggsave(file.path(base_dir, "variance_plot.pdf"), p4, width = 10, height = 6)
print(p4)

# ---------------------------------------------------------
# 5. NOISY OBJECTIVE PLOT
# ---------------------------------------------------------

example_instance <- 1
raw_evaluations <- bo_results[[example_instance]]$y_history
best_so_far <- cummax(raw_evaluations)

noise_example <- data.frame(
  iteration = 1:evalBudget,
  raw = raw_evaluations,
  best = best_so_far
) %>%
  pivot_longer(cols = c(raw, best), names_to = "type", values_to = "coverage")

p5 <- ggplot(noise_example, aes(x = iteration, y = coverage, color = type)) +
  geom_line(data = subset(noise_example, type == "best"), linewidth = 1.2) +
  geom_point(data = subset(noise_example, type == "raw"), alpha = 0.3, size = 1) +
  labs(
    title = "Effect of Stochastic Objective (Example BO Run)",
    subtitle = "Points show noisy evaluations, line shows best-so-far",
    x = "Iteration",
    y = "Coverage (Population)",
    color = ""
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 14)
  ) +
  scale_color_manual(
    values = c("raw" = unname(okabe_ito["black"]), "best" = bo_color),
    labels = c("Best-so-far", "Raw evaluations")
  )

ggsave(file.path(base_dir, "stochastic_objective_example.pdf"), p5, width = 10, height = 6)
print(p5)

# ---------------------------------------------------------
# 6. PERFORMANCE TABLE
# ---------------------------------------------------------

performance_table <- data.frame(
  Method = c("BO", "GA"),
  Mean_Coverage = c(mean(bo_summary$best_coverage), mean(ga_summary$best_coverage)),
  SD_Coverage = c(sd(bo_summary$best_coverage), sd(ga_summary$best_coverage)),
  Min_Coverage = c(min(bo_summary$best_coverage), min(ga_summary$best_coverage)),
  Max_Coverage = c(max(bo_summary$best_coverage), max(ga_summary$best_coverage)),
  Mean_Runtime_Min = c(mean(bo_summary$runtime_min), mean(ga_summary$runtime_min))
)

print(performance_table)
write.csv(performance_table, 
          file.path(base_dir, "performance_summary.csv"), 
          row.names = FALSE)

# ---------------------------------------------------------
# 8. COMBINED SUMMARY PLOT
# ---------------------------------------------------------

library(gridExtra)

pdf(file.path(base_dir, "combined_summary.pdf"), width = 14, height = 10)
grid.arrange(p1, p3, p4, p5, ncol = 2)
dev.off()

cat("\n=== Plots saved to:", base_dir, "===\n")

# ---------------------------------------------------------
# DETERMINISTIC EVALUATION FUNCTION
# ---------------------------------------------------------

# Evaluate solution on TRUE MEAN population (no noise)
mclp_deterministic <- function(x) {
  stopifnot(length(x) == n_vars)
  
  coverage_vector <- pmin(A %*% x, 1)
  # Use MEAN population (deterministic)
  total_coverage <- sum(coverage_vector * population_mean)
  
  return(total_coverage)
}

# ---------------------------------------------------------
# 1. EVALUATE ALL SOLUTIONS DETERMINISTICALLY
# ---------------------------------------------------------

cat("Evaluating all solutions on deterministic (mean) population...\n")

# For BO: get best solution from each run
bo_deterministic_scores <- numeric(n_instances)
for (i in 1:n_instances) {
  best_idx <- which.max(bo_results[[i]]$y_history)
  best_solution <- bo_results[[i]]$x_history[[best_idx]]
  bo_deterministic_scores[i] <- mclp_deterministic(best_solution)
}

# For GA: get best solution from each run
ga_deterministic_scores <- numeric(n_instances)
for (i in 1:n_instances) {
  ga_deterministic_scores[i] <- mclp_deterministic(ga_results[[i]]$best_solution)
}

cat("BO mean deterministic score:", mean(bo_deterministic_scores), "\n")
cat("GA mean deterministic score:", mean(ga_deterministic_scores), "\n")

# ---------------------------------------------------------
# 2. COMPARE STOCHASTIC VS DETERMINISTIC PERFORMANCE
# ---------------------------------------------------------

comparison_stoch_vs_det <- data.frame(
  instance = rep(1:n_instances, 2),
  method = rep(c("BO", "GA"), each = n_instances),
  stochastic_score = c(bo_summary$best_coverage, ga_summary$best_coverage),
  deterministic_score = c(bo_deterministic_scores, ga_deterministic_scores)
)

comparison_stoch_vs_det$score_difference <- comparison_stoch_vs_det$stochastic_score - 
  comparison_stoch_vs_det$deterministic_score

# Plot: Stochastic vs Deterministic
p_stoch_det <- ggplot(comparison_stoch_vs_det, 
                      aes(x = stochastic_score, y = deterministic_score, 
                          color = method, shape = method)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  labs(
    title = "Stochastic vs Deterministic Performance",
    subtitle = "Points above line indicate stochastic noise helped; below line indicates noise hurt",
    x = "Best Score (Stochastic Evaluation)",
    y = "Score (Deterministic Evaluation on Mean Population)",
    color = "Method",
    shape = "Method"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 14)
  ) +
  scale_color_manual(values = c("BO" = bo_color, "GA" = ga_color))

ggsave(file.path(base_dir, "stochastic_vs_deterministic.pdf"), p_stoch_det, width = 8, height = 7)
print(p_stoch_det)

# ---------------------------------------------------------
# 3. NOISE IMPACT ANALYSIS
# ---------------------------------------------------------

# How much does noise affect each method?
noise_impact <- data.frame(
  method = c("BO", "GA"),
  mean_stochastic = c(mean(bo_summary$best_coverage), mean(ga_summary$best_coverage)),
  mean_deterministic = c(mean(bo_deterministic_scores), mean(ga_deterministic_scores)),
  sd_stochastic = c(sd(bo_summary$best_coverage), sd(ga_summary$best_coverage)),
  sd_deterministic = c(sd(bo_deterministic_scores), sd(ga_deterministic_scores))
)

noise_impact$noise_effect <- noise_impact$mean_stochastic - noise_impact$mean_deterministic
noise_impact$cv_stochastic <- noise_impact$sd_stochastic / noise_impact$mean_stochastic
noise_impact$cv_deterministic <- noise_impact$sd_deterministic / noise_impact$mean_deterministic

print("=== Noise Impact Analysis ===")
print(noise_impact)

write.csv(noise_impact, 
          file.path(base_dir, "noise_impact_analysis.csv"), 
          row.names = FALSE)

# ---------------------------------------------------------
# 4. ROBUSTNESS TO NOISE: Coefficient of Variation
# ---------------------------------------------------------

robustness_data <- data.frame(
  method = rep(c("BO", "GA"), 2),
  evaluation_type = rep(c("Stochastic", "Deterministic"), each = 2),
  cv = c(noise_impact$cv_stochastic[1], noise_impact$cv_stochastic[2],
         noise_impact$cv_deterministic[1], noise_impact$cv_deterministic[2])
)

p_robustness <- ggplot(robustness_data, aes(x = method, y = cv, fill = evaluation_type)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.7, width = 0.6) +
  labs(
    title = "Robustness to Noise: Coefficient of Variation",
    subtitle = "Lower CV indicates more consistent performance",
    x = "Method",
    y = "Coefficient of Variation (SD / Mean)",
    fill = "Evaluation Type"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 14)
  ) +
  scale_fill_manual(values = c("Stochastic" = "gray60", "Deterministic" = "gray30"))

ggsave(file.path(base_dir, "robustness_analysis.pdf"), p_robustness, width = 8, height = 6)
print(p_robustness)

# ---------------------------------------------------------
# 5. PAIRED COMPARISON: Deterministic Performance
# ---------------------------------------------------------

# Boxplot of deterministic scores
det_comparison <- data.frame(
  coverage = c(bo_deterministic_scores, ga_deterministic_scores),
  method = rep(c("BO", "GA"), each = n_instances)
)

p_det_boxplot <- ggplot(det_comparison, aes(x = method, y = coverage, fill = method)) +
  geom_boxplot(alpha = 0.7, width = 0.5) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 2) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 4, fill = "white") +
  labs(
    title = "True Optimization Performance (Deterministic Evaluation)",
    subtitle = "All solutions evaluated on mean population (no noise)",
    x = "Method",
    y = "Coverage on Mean Population"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 14)
  ) +
  scale_fill_manual(values = c("BO" = bo_color, "GA" = ga_color))

ggsave(file.path(base_dir, "deterministic_performance.pdf"), p_det_boxplot, width = 7, height = 6)
print(p_det_boxplot)
