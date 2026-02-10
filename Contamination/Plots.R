library(dplyr)
library(stringr)
library(ggplot2)

all_data <- NULL   # will collect everything here

lambda_dirs <- list.dirs(
  "/Users/niyati/Projects:Codes/VaR-CoCoBO/Contamination/results",
  recursive = FALSE
)

for (ldir in lambda_dirs) {
  
  cat("Processing lambda dir:", ldir, "\n")
  
  # extract lambda value from directory name
  lambda_val <- str_extract(basename(ldir), "[0-9.]+")
  lambda_val <- as.numeric(lambda_val)
  
  # list instance files
  files <- list.files(
    ldir,
    pattern = "instance_.*\\.RData",
    full.names = TRUE
  )
  
  cat("  Found", length(files), "files\n")
  
  for (f in files) {
    
    cat("Loading:", f, "\n")
    
    load(f)   # loads object `res`
    
    ## ---- DEBUG POINT 1 ----
    ## Check what was loaded
    # ls()
    # str(res)
    
    df <- res$data
    
    ## ---- DEBUG POINT 2 ----
    if (!is.data.frame(df)) {
      warning("res$data is not a data.frame in ", f)
      next
    }
    
    df$iter <- seq_len(nrow(df))
    
    df <- df %>%
      mutate(
        lambda   = lambda_val,
        instance = res$instance,
        n_init   = res$n_init
      )
    
    ## ---- DEBUG POINT 3 ----
    # str(df)
    # head(df)
    
    # bind to master data
    if (is.null(all_data)) {
      all_data <- df
    } else {
      all_data <- bind_rows(all_data, df)
    }
  }
}

all_data <- all_data %>%
  arrange(lambda, instance, iter) %>%
  group_by(lambda, instance) %>%
  mutate(
    best_so_far = {
      bf <- rep(NA_real_, n())
      current_best <- Inf
      
      for (i in seq_len(n())) {
        if (Feasibility[i]) {
          current_best <- min(current_best, y[i])
        }
        if (is.finite(current_best)) {
          bf[i] <- current_best
        }
      }
      bf
    }
  ) %>%
  ungroup()


summary_trace <- all_data %>%
  group_by(lambda, iter) %>%
  summarise(
    mean_best = mean(best_so_far, na.rm = TRUE),
    sd_best   = sd(best_so_far, na.rm = TRUE),
    n         = sum(!is.na(best_so_far)),
    se_best   = sd_best / sqrt(n),
    lower     = mean_best - 1.96 * se_best,
    upper     = mean_best + 1.96 * se_best,
    .groups = "drop"
  )

lambda_labeller <- labeller(
  lambda = function(x) paste0("\u03BB = ", x)  # Greek lambda
)

okabe_ito <- c(
  "#000000", # black
  "#E69F00", # orange
  "#56B4E9", # sky blue
  "#009E73", # bluish green
  "#F0E442", # yellow
  "#0072B2", # blue
  "#D55E00", # vermillion
  "#CC79A7"  # reddish purple
)


trace_plot <- ggplot(
  summary_trace,
  aes(
    x = iter,
    y = mean_best,
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
    xintercept = 10,  #n_init
    linetype = "dotted",
    linewidth = 1,
    inherit.aes = FALSE
  ) +
  
  facet_wrap(
    ~ lambda,
    scales = "free_y",
    labeller = lambda_labeller
  ) +
  
  scale_color_manual(
    values = okabe_ito,
    name = expression(lambda)
  ) +
  
  scale_fill_manual(
    values = okabe_ito,
    name = expression(lambda)
  ) +
  
  labs(
    x = "Iteration",
    y = "Best feasible objective so far"
  ) +
  
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",   # redundant with facets
    strip.text = element_text(size = 13)
  )

ggsave(
  filename = "/Users/niyati/Projects:Codes/VaR-CoCoBO/Contamination/contan_faceted_trace.pdf",
  plot = trace_plot,
  device = cairo_pdf,
  width = 11, height = 6, units = "in",
  dpi = 1200, bg = "white"
)

#-----------------FEASIBILITY--------------------------

inst_id <- 4

inst_data <- all_data %>%
  filter(instance == inst_id) %>%
  mutate(
    init = iter <= n_init,
    feasible = Feasibility  # assuming this already exists
  )

inst_data <- inst_data %>%
  mutate(trace = "Best feasible optimum")

feasible <- ggplot() +
  
  # --- Best feasible trace ---
  geom_line(
    data = inst_data,
    aes(
      x = iter,
      y = best_so_far,
      linetype = trace
    ),
    color = "darkgreen",
    linewidth = 1
  ) +
  
  scale_linetype_manual(
    values = c("Best feasible optimum" = "solid"),
    name = NULL
  ) +
  
  # --- All evaluated points ---
  geom_point(
    data = inst_data,
    aes(
      x = iter,
      y = y,
      color = feasible,
      shape = init
    ),
    size = 2,
    alpha = ifelse(inst_data$feasible, 1, 0.4)
  ) +
  
  # --- Initial design separator ---
  geom_vline(
    data = inst_data %>% distinct(n_init),
    aes(xintercept = n_init),
    linetype = "dotted",
    linewidth = 1
  ) +
  
  # --- Facet by lambda ---
  facet_wrap(~ lambda, scales = "free_y",
             labeller = lambda_labeller) +
  
  scale_color_manual(
    values = c("FALSE" = "red", "TRUE" = "black"),
    labels = c("FALSE" = "No", "TRUE" = "Yes"),
    name = "Feasible"
  ) +
  scale_shape_manual(
    values = c(`TRUE` = 17, `FALSE` = 16),
    labels = c(`TRUE` = "Initial Samples", `FALSE` = "BO Samples"),
    name = ""
  ) +
  
  labs(
    x = "Iteration",
    y = "Objective value",
    title = paste("Feasibility and optimisation trace for one instance")
  ) +
  
  theme_bw(base_size = 14) +
  theme(
    legend.position = "right",
    strip.text = element_text(size = 13)
  )
ggsave(
  filename = "/Users/niyati/Projects:Codes/VaR-CoCoBO/Contamination/contan_faceted_feasibility.pdf",
  plot = feasible,
  device = cairo_pdf,
  width = 11, height = 6, units = "in",
  dpi = 1200, bg = "white"
)

#-------------Time--------------------------

runtime_data <- NULL

for (ldir in lambda_dirs) {
  lambda_val <- as.numeric(str_extract(basename(ldir), "[0-9.]+"))
  files <- list.files(ldir, pattern = "instance_.*\\.RData", full.names = TRUE)
  
  for (f in files) {
    load(f)  # loads res
    runtime_data <- bind_rows(
      runtime_data,
      data.frame(
        lambda = lambda_val,
        instance = res$instance,
        runtime_sec = res$runtime_sec
      )
    )
  }
}

latex_tab <- runtime_data %>%
  group_by(lambda) %>%
  summarise(
    runtime_min    = mean(runtime_sec) / 60,
    runtime_sd_min = sd(runtime_sec) / 60,
    .groups = "drop"
  )

latex_tab_fmt <- latex_tab %>%
  mutate(
    runtime = sprintf(
      "%.2f $\\pm$ %.2f",
      runtime_min,
      runtime_sd_min
    )
  )

library(xtable)
print(
  xtable(latex_tab_fmt),
  include.rownames = FALSE,
  sanitize.text.function = identity
)


