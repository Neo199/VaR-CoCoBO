library(dplyr)
library(stringr)
library(ggplot2)

# ---------------------------------------------------------
# LOAD ALL DATA FROM NEW STRUCTURE
# ---------------------------------------------------------

base_dir <- "/Users/niyati/Projects:Codes/VaR-CoCoBO/Contamination/results_corrections"

lambda_dirs <- list.dirs(base_dir, recursive = FALSE)

# Will collect each condition separately
all_full      <- NULL
all_ab1       <- NULL
all_ab2       <- NULL
runtime_data  <- NULL

condition_labels <- c(
  "full"                  = "Full method (Thompson + Gumbel)",
  "ablation_no_gumbel"    = "Ablation 1: Thompson only",
  "ablation_gumbel_mean"  = "Ablation 2: Gumbel + E[\u03b8]"
)

for (ldir in lambda_dirs) {
  
  cat("Processing lambda dir:", ldir, "\n")
  
  lambda_val <- as.numeric(str_extract(basename(ldir), "[0-9.]+"))
  
  files <- list.files(ldir, pattern = "instance_.*\\.RData", full.names = TRUE)
  cat("  Found", length(files), "files\n")
  
  for (f in files) {
    
    cat("  Loading:", f, "\n")
    load(f)  # loads object `res`
    
    # runtime
    runtime_data <- bind_rows(
      runtime_data,
      data.frame(
        lambda      = lambda_val,
        instance    = res$instance,
        runtime_sec = res$runtime_sec
      )
    )
    
    # helper to extract and label one condition
    extract_condition <- function(cond_list, cond_name) {
      df <- cond_list$data
      if (!is.data.frame(df)) {
        warning("Not a data.frame for condition ", cond_name, " in ", f)
        return(NULL)
      }
      df$iter      <- seq_len(nrow(df))
      df$lambda    <- lambda_val
      df$instance  <- res$instance
      df$n_init    <- res$n_init
      df$condition <- cond_name
      df
    }
    
    all_full <- bind_rows(all_full,
                          extract_condition(res$full, "full"))
    all_ab1  <- bind_rows(all_ab1,
                          extract_condition(res$ablation_no_gumbel, "ablation_no_gumbel"))
    all_ab2  <- bind_rows(all_ab2,
                          extract_condition(res$ablation_gumbel_mean, "ablation_gumbel_mean"))
  }
}

# replace this:
# all_data <- bind_rows(all_full, all_ab1, all_ab2)

# with this:
all_data <- bind_rows(all_full, all_ab1)
all_data <- bind_rows(all_data, all_ab2)
# ---------------------------------------------------------
# COMPUTE BEST FEASIBLE SO FAR — do this BEFORE factor conversion
# ---------------------------------------------------------

all_data <- all_data %>%
  arrange(condition, lambda, instance, iter) %>%
  group_by(condition, lambda, instance) %>%
  mutate(
    best_so_far = {
      bf           <- rep(NA_real_, n())
      current_best <- Inf
      for (i in seq_len(n())) {
        if (!is.na(Feasibility[i]) && Feasibility[i]) {
          current_best <- min(current_best, y[i])
        }
        if (is.finite(current_best)) bf[i] <- current_best
      }
      bf
    }
  ) %>%
  ungroup()

# ---------------------------------------------------------
# NOW convert to factor for plotting
# ---------------------------------------------------------

all_data <- all_data %>%
  mutate(condition = factor(
    condition,
    levels = c(
      "full",
      "ablation_no_gumbel",
      "ablation_gumbel_mean"
    ),
    labels = c(
      "Thompson + Gumbel",
      "Thompson only",
      "Gumbel + E[\u03b8]"
    )
  ))

# confirm
print(unique(all_data$condition))
print(all_data %>% count(condition))

# ---------------------------------------------------------
# COLOUR PALETTE (Okabe-Ito)
# ---------------------------------------------------------

okabe_ito <- c(
  "#000000",
  "#E69F00",
  "#56B4E9",
  "#009E73",
  "#F0E442",
  "#0072B2",
  "#D55E00",
  "#CC79A7"
)

lambda_labeller <- labeller(
  lambda = function(x) paste0("\u03BB = ", x)
)

# ---------------------------------------------------------
# PLOT 1: LAMBDA COMPARISON — full method only
# ---------------------------------------------------------

full_data <- all_data %>%
  filter(condition == "Thompson + Gumbel")

summary_full <- full_data %>%
  group_by(lambda, iter) %>%
  summarise(
    mean_best = mean(best_so_far, na.rm = TRUE),
    sd_best   = sd(best_so_far,   na.rm = TRUE),
    n         = sum(!is.na(best_so_far)),
    se_best   = sd_best / sqrt(n),
    lower     = mean_best - 1.96 * se_best,
    upper     = mean_best + 1.96 * se_best,
    .groups   = "drop"
  )

plot_lambda <- ggplot(
  summary_full,
  aes(
    x     = iter,
    y     = mean_best,
    color = factor(lambda),
    fill  = factor(lambda)
  )
) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper),
    alpha = 0.25,
    color = NA
  ) +
  geom_line(linewidth = 1) +
  geom_vline(
    xintercept  = 10,
    linetype    = "dotted",
    linewidth   = 1,
    inherit.aes = FALSE
  ) +
  facet_wrap(
    ~ lambda,
    scales   = "free_y",
    labeller = lambda_labeller
  ) +
  scale_color_manual(values = okabe_ito, name = expression(lambda)) +
  scale_fill_manual( values = okabe_ito, name = expression(lambda)) +
  labs(
    x = "Iteration",
    y = "Best feasible objective so far"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text      = element_text(size = 13)
  )

ggsave(
  filename = file.path(base_dir, "contan_faceted_trace.pdf"),
  plot     = plot_lambda,
  device   = cairo_pdf,
  width = 11, height = 6, units = "in",
  dpi = 1200, bg = "white"
)

# ---------------------------------------------------------
# PLOT 2: ABLATION COMPARISON — all three conditions
# ---------------------------------------------------------

# rep_lambda <- 1
# 
# summary_ablation <- all_data %>%
#   filter(lambda == rep_lambda) %>%
#   group_by(condition, iter) %>%
#   summarise(
#     mean_best = mean(best_so_far, na.rm = TRUE),
#     sd_best   = sd(best_so_far,   na.rm = TRUE),
#     n         = sum(!is.na(best_so_far)),
#     se_best   = sd_best / sqrt(n),
#     lower     = mean_best - 1.96 * se_best,
#     upper     = mean_best + 1.96 * se_best,
#     .groups   = "drop"
#   )
# 
# plot_ablation <- ggplot(
#   summary_ablation,
#   aes(
#     x     = iter,
#     y     = mean_best,
#     color = condition,
#     fill  = condition
#   )
# ) +
#   geom_ribbon(
#     aes(ymin = lower, ymax = upper),
#     alpha = 0.25,
#     color = NA
#   ) +
#   geom_line(linewidth = 1) +
#   geom_vline(
#     xintercept  = 10,
#     linetype    = "dotted",
#     linewidth   = 1,
#     inherit.aes = FALSE
#   ) +
#   facet_wrap(
#     ~ condition,
#     scales = "free_y"
#   ) +
#   scale_color_manual(
#     values = c(
#       "Full method (Thompson + Gumbel)" = "#009E73",
#       "Ablation 1: Thompson only"       = "#D55E00",
#       "Ablation 2: Gumbel + E[\u03b8]" = "#0072B2"
#     )
#   ) +
#   scale_fill_manual(
#     values = c(
#       "Full method (Thompson + Gumbel)" = "#009E73",
#       "Ablation 1: Thompson only"       = "#D55E00",
#       "Ablation 2: Gumbel + E[\u03b8]" = "#0072B2"
#     )
#   ) +
#   labs(
#     x = "Iteration",
#     y = "Best feasible objective so far"
#   ) +
#   theme_bw(base_size = 14) +
#   theme(
#     legend.position = "none",
#     strip.text      = element_text(size = 13)
#   )
# 
# ggsave(
#   filename = file.path(base_dir, "contan_ablation_trace.pdf"),
#   plot     = plot_ablation,
#   device   = cairo_pdf,
#   width    = 13, height = 5, units = "in",
#   dpi      = 1200, bg = "white"
# )

# ---------------------------------------------------------
# PLOT 2: ABLATION COMPARISON — all conditions, faceted by lambda
# Fixed y-axis for direct comparison
# ---------------------------------------------------------

summary_ablation <- all_data %>%
  group_by(condition, lambda, iter) %>%
  summarise(
    mean_best = mean(best_so_far, na.rm = TRUE),
    sd_best   = sd(best_so_far,   na.rm = TRUE),
    n         = sum(!is.na(best_so_far)),
    se_best   = ifelse(n > 1, sd_best / sqrt(n), NA_real_),
    lower     = ifelse(n > 1, mean_best - 1.96 * se_best, NA_real_),
    upper     = ifelse(n > 1, mean_best + 1.96 * se_best, NA_real_),
    .groups   = "drop"
  ) %>%
  mutate(
    mean_best = ifelse(is.nan(mean_best), NA_real_, mean_best),
    lower     = ifelse(is.nan(lower),     NA_real_, lower),
    upper     = ifelse(is.nan(upper),     NA_real_, upper)
  )

# compute shared y limits across all conditions and lambdas
y_min <- min(summary_ablation$lower,    na.rm = TRUE)
y_max <- max(summary_ablation$upper,    na.rm = TRUE)

plot_ablation <- ggplot(
  summary_ablation,
  aes(
    x     = iter,
    y     = mean_best,
    color = condition,
    fill  = condition,
    linetype = condition
  )
) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper),
    alpha = 0.15,
    color = NA
  ) +
  geom_line(linewidth = 1) +
  geom_vline(
    xintercept  = 10,
    linetype    = "dotted",
    color       = "grey40",
    linewidth   = 0.8,
    inherit.aes = FALSE
  ) +
  facet_wrap(
    ~ lambda,
    labeller = lambda_labeller
  ) +
  coord_cartesian(ylim = c(y_min, y_max)) +
  scale_color_manual(
    name   = "Method",
    values = c(
      "Thompson + Gumbel" = "#009E73",
      "Thompson only"     = "#D55E00",
      "Gumbel + E[\u03b8]" = "#0072B2"
    )
  ) +
  scale_fill_manual(
    name   = "Method",
    values = c(
      "Thompson + Gumbel" = "#009E73",
      "Thompson only"     = "#D55E00",
      "Gumbel + E[\u03b8]" = "#0072B2"
    )
  ) +
  scale_linetype_manual(
    name   = "Method",
    values = c(
      "Thompson + Gumbel" = "solid",
      "Thompson only"     = "dashed",
      "Gumbel + E[\u03b8]" = "dotdash"
    )
  ) +
  labs(
    x = "Iteration",
    y = "Best feasible objective so far"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "bottom",
    strip.text      = element_text(size = 13),
    legend.title    = element_text(size = 12),
    legend.text     = element_text(size = 11)
  )

ggsave(
  filename = file.path(base_dir, "contan_ablation_trace.pdf"),
  plot     = plot_ablation,
  device   = cairo_pdf,
  width    = 13, height = 8, units = "in",
  dpi      = 1200, bg = "white"
)
# ---------------------------------------------------------
# PLOT 3a: FEASIBILITY TRACE — one instance, fixed y-axis
# ---------------------------------------------------------

rep_inst   <- 4
rep_lambda <- 1

inst_data <- all_data %>%
  filter(lambda == rep_lambda, instance == rep_inst) %>%
  mutate(
    init     = iter <= n_init,
    feasible = !is.na(Feasibility) & Feasibility
  )

# fixed y limits from this instance across all conditions
y_min_inst <- min(inst_data$y, na.rm = TRUE)
y_max_inst <- max(inst_data$y, na.rm = TRUE)

plot_feasibility <- ggplot() +
  geom_line(
    data = inst_data,
    aes(x = iter, y = best_so_far),
    color     = "darkgreen",
    linewidth = 1
  ) +
  geom_point(
    data = inst_data,
    aes(
      x     = iter,
      y     = y,
      color = feasible,
      shape = init
    ),
    size  = 2,
    alpha = ifelse(!is.na(inst_data$feasible) & inst_data$feasible, 1, 0.4)
  ) +
  geom_vline(
    data = inst_data %>% distinct(condition, n_init),
    aes(xintercept = n_init),
    linetype  = "dotted",
    linewidth = 1
  ) +
  facet_wrap(~ condition) +
  coord_cartesian(ylim = c(y_min_inst, y_max_inst)) +
  scale_color_manual(
    values = c("FALSE" = "red", "TRUE" = "black"),
    labels = c("FALSE" = "No",  "TRUE" = "Yes"),
    name   = "Feasible"
  ) +
  scale_shape_manual(
    values = c(`TRUE` = 17, `FALSE` = 16),
    labels = c(`TRUE` = "Initial samples", `FALSE` = "BO samples"),
    name   = ""
  ) +
  labs(
    x = "Iteration",
    y = "Objective value"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "right",
    strip.text      = element_text(size = 13)
  )

ggsave(
  filename = file.path(base_dir, "contan_ablation_feasibility_4.pdf"),
  plot     = plot_feasibility,
  device   = cairo_pdf,
  width    = 13, height = 5, units = "in",
  dpi      = 1200, bg = "white"
)

# ---------------------------------------------------------
# PLOT 3b: AVERAGE FEASIBILITY RATE — all instances
# Rolling mean across iterations, faceted by lambda
# Shows proportion of feasible solutions per method
# ---------------------------------------------------------

library(zoo)

feas_summary <- all_data %>%
  mutate(feasible_num = as.numeric(!is.na(Feasibility) & Feasibility)) %>%
  group_by(condition, lambda, iter) %>%
  summarise(
    feas_rate = mean(feasible_num, na.rm = TRUE),
    .groups   = "drop"
  ) %>%
  group_by(condition, lambda) %>%
  mutate(
    feas_rate_smooth = rollmean(feas_rate, k = 10, fill = NA, align = "right")
  ) %>%
  ungroup()

plot_feas_rate <- ggplot(
  feas_summary,
  aes(
    x        = iter,
    y        = feas_rate_smooth,
    color    = condition,
    linetype = condition
  )
) +
  geom_line(linewidth = 1, na.rm = TRUE) +
  geom_vline(
    xintercept  = 10,
    linetype    = "dotted",
    color       = "grey40",
    linewidth   = 0.8,
    inherit.aes = FALSE
  ) +
  facet_wrap(
    ~ lambda,
    labeller = lambda_labeller
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_color_manual(
    name   = "Method",
    values = c(
      "Thompson + Gumbel"  = "#009E73",
      "Thompson only"      = "#D55E00",
      "Gumbel + E[\u03b8]" = "#0072B2"
    )
  ) +
  scale_linetype_manual(
    name   = "Method",
    values = c(
      "Thompson + Gumbel"  = "solid",
      "Thompson only"      = "dashed",
      "Gumbel + E[\u03b8]" = "dotdash"
    )
  ) +
  labs(
    x = "Iteration",
    y = "Proportion feasible (rolling mean, k = 10)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "bottom",
    strip.text      = element_text(size = 13),
    legend.title    = element_text(size = 12),
    legend.text     = element_text(size = 11)
  )

ggsave(
  filename = file.path(base_dir, "contan_feasibility_rate.pdf"),
  plot     = plot_feas_rate,
  device   = cairo_pdf,
  width    = 13, height = 8, units = "in",
  dpi      = 1200, bg = "white"
)
# ---------------------------------------------------------
# RUNTIME TABLE — full method only, by lambda
# ---------------------------------------------------------

library(xtable)

latex_tab <- runtime_data %>%
  group_by(lambda) %>%
  summarise(
    runtime_min    = mean(runtime_sec) / 60,
    runtime_sd_min = sd(runtime_sec)   / 60,
    .groups        = "drop"
  ) %>%
  mutate(
    runtime = sprintf("%.2f $\\pm$ %.2f", runtime_min, runtime_sd_min)
  ) %>%
  select(lambda, runtime)

print(
  xtable(latex_tab),
  include.rownames       = FALSE,
  sanitize.text.function = identity
)

