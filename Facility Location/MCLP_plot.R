# PLOTTING FOR MCLP
# Load libraries
library(sf)
library(ggplot2)
library(osmdata)
library(dplyr)

# ---------------------------------------------------------
# LOAD Data
# ---------------------------------------------------------

demand_data <- read.csv("PhD-Data/SF_data/SF_demand_205_centroid_uniform_weight.csv")
facility_loc <- read.csv("PhD-Data/SF_data/SF_store_site_16_longlat.csv")
distance_matrix <- read.csv("PhD-Data/SF_data/SF_network_distance_candidateStore_16_censusTract_205_new.csv")

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

load(file = "P3Compute/Facility Location/results_MCLP/multiple_instances/instance_1.RData")

results_df <- data.frame(
  iteration = 1:250,
  objective = result$y_history,
  coverage = result$coverage_history,
  coverage_pct = 100 * (result$coverage_history) / sum(population),
  best_coverage = cummax(result$coverage_history)
)

# Create spatial objects
facility_sf <- st_as_sf(
  facility_loc,
  coords = c("long", "lat"),
  crs = 4326
)

demand_sf <- st_as_sf(
  demand_data,
  coords = c("long", "lat"),
  crs = 4326
)


best_bo_iteration <- which.max(result$coverage_history)
best_bo_coverage <- max(result$coverage_history)
best_bo_solution <- result$x_history[[best_bo_iteration]]

# Assign facility selection status
facility_sf$selected <- best_bo_solution

# # Coverage assignment
coverage <- A %*% best_bo_solution
demand_sf$covered <- coverage > 0
# 
# # Calculate coverage statistics
# total_demand <- sum(demand_data$POP2000)  # assuming you have weights
# covered_demand <- sum(demand_data$POP2000[demand_sf$covered])
coverage_pct <- max(results_df$coverage_pct)

# Calculate bounding box
bbox_all <- st_bbox(st_union(
  st_geometry(facility_sf),
  st_geometry(demand_sf)
))

pad_x <- (bbox_all["xmax"] - bbox_all["xmin"]) * 0.05
pad_y <- (bbox_all["ymax"] - bbox_all["ymin"]) * 0.05

bbox_expanded <- bbox_all
bbox_expanded[c("xmin","xmax")] <- bbox_expanded[c("xmin","xmax")] + c(-pad_x, pad_x)
bbox_expanded[c("ymin","ymax")] <- bbox_expanded[c("ymin","ymax")] + c(-pad_y, pad_y)

# Fetch OpenStreetMap data
osm_map <- opq(bbox = bbox_expanded) %>%
  add_osm_feature(key = "highway") %>%
  osmdata_sf()

streets <- osm_map$osm_lines

# Okabe-Ito color palette
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

# Prepare data with proper factor levels for legend control
demand_sf <- demand_sf %>%
  mutate(coverage_status = factor(
    ifelse(covered, "Covered", "Uncovered"),
    levels = c("Covered", "Uncovered")
  ))

facility_sf <- facility_sf %>%
  mutate(selection_status = factor(
    ifelse(selected == 1, "Selected", "Not selected"),
    levels = c("Selected", "Not selected")
  ))

# Create the plot
plot_mclp <- ggplot() +
  
  # Street network
  geom_sf(
    data = streets,
    colour = "grey",
    linewidth = 0.3,
    alpha = 0.3
  ) +
  
  # Demand nodes
  geom_sf(
    data = demand_sf,
    aes(colour = coverage_status, shape = coverage_status),
    size = 2.5,
    alpha = 0.8
  ) +
  
  # Facilities
  geom_sf(
    data = facility_sf,
    aes(fill = selection_status),
    shape = 23,
    size = 4,
    colour = "grey",
    stroke = 0.8
  ) +
  
  # Demand colour scale
  scale_colour_manual(
    name = "Demand Status",
    values = c(
      "Covered"   = okabe_ito[4],
      "Uncovered" = okabe_ito[7]
    ),
    guide = guide_legend(order = 1)
  ) +
  
  # Demand shape scale
  scale_shape_manual(
    name = "Demand Status",
    values = c(
      "Covered"   = 16,
      "Uncovered" = 1
    ),
    guide = "none"
  ) +
  
  # Facility fill scale
  scale_fill_manual(
    name = "Facility Status",
    values = c(
      "Selected"     = okabe_ito[6],
      "Not selected" = okabe_ito[8]
    ),
    guide = guide_legend(
      order = 2, 
      override.aes = list(shape = 23, size = 4, colour = okabe_ito[1], stroke = 0.8)
    )
  ) +
  
  coord_sf(
    xlim = c(bbox_expanded["xmin"], bbox_expanded["xmax"]),
    ylim = c(bbox_expanded["ymin"], bbox_expanded["ymax"]),
    expand = FALSE
  ) +
  
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.spacing.y = unit(0.5, "cm"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey90", linewidth = 0.2),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.caption = element_text(hjust = 0.5, size = 11)
  ) +
  
  labs(
    title = "MCLP Solution: Maximum Coverage with Limited Facilities",
    caption = paste0(
      "Facilities selected: ", sum(best_bo_solution), 
      " | Demand covered: ", sum(demand_sf$covered), "/", nrow(demand_sf),
      " (", round(coverage_pct,2), "% of total demand)"
    )
  )

# Save the plot
ggsave(
  filename = "/Users/niyati/Projects:Codes/P3Compute/Facility Location/MCLP.pdf",
  plot = plot_mclp,
  device = cairo_pdf,
  width = 11, height = 10, units = "in",
  dpi = 1200, bg = "white"
)

# Optional: Print summary statistics
cat("\n=== MCLP Solution Summary ===\n")
cat("Facilities selected:", sum(best_bo_solution), "\n")
cat("Demand nodes covered:", sum(demand_sf$covered), "/", nrow(demand_sf), "\n")
cat("Coverage percentage:", coverage_pct, "%\n")
cat("Uncovered demand nodes:", sum(!demand_sf$covered), "\n")

# -----------------------------
# GLPK solution 
# -----------------------------

best_lp_solution <- lp_solution$selected


# -----------------------------
# Coverage from LP solution
# -----------------------------
coverage_lp <- A %*% best_lp_solution
covered_lp  <- as.vector(coverage_lp > 0)
coverage_pct_lp <- round(lp_summary$optimal_pct, 2)

facility_sf_lp <- st_as_sf(
  facility_loc,
  coords = c("long", "lat"),
  crs = 4326
)

demand_sf_lp <- st_as_sf(
  demand_data,
  coords = c("long", "lat"),
  crs = 4326
)

facility_sf_lp$selected <- best_lp_solution
demand_sf_lp$covered <- covered_lp

demand_sf_lp <- demand_sf_lp %>%
  mutate(coverage_status = factor(
    ifelse(covered, "Covered", "Uncovered"),
    levels = c("Covered", "Uncovered")
  ))

facility_sf_lp <- facility_sf_lp %>%
  mutate(selection_status = factor(
    ifelse(selected == 1, "Selected", "Not selected"),
    levels = c("Selected", "Not selected")
  ))

bbox_all_lp <- st_bbox(st_union(
  st_geometry(facility_sf_lp),
  st_geometry(demand_sf_lp)
))

pad_x <- (bbox_all_lp["xmax"] - bbox_all_lp["xmin"]) * 0.05
pad_y <- (bbox_all_lp["ymax"] - bbox_all_lp["ymin"]) * 0.05

bbox_expanded_lp <- bbox_all_lp
bbox_expanded_lp[c("xmin","xmax")] <- bbox_expanded_lp[c("xmin","xmax")] + c(-pad_x, pad_x)
bbox_expanded_lp[c("ymin","ymax")] <- bbox_expanded_lp[c("ymin","ymax")] + c(-pad_y, pad_y)

osm_map_lp <- opq(bbox = bbox_expanded_lp) %>%
  add_osm_feature(key = "highway") %>%
  osmdata_sf()

streets_lp <- osm_map_lp$osm_lines

plot_mclp_lp <- ggplot() +
  
  geom_sf(data = streets_lp,
          colour = "grey",
          linewidth = 0.3,
          alpha = 0.3) +
  
  geom_sf(
    data = demand_sf_lp,
    aes(colour = coverage_status, shape = coverage_status),
    size = 2.5,
    alpha = 0.8
  ) +
  
  geom_sf(
    data = facility_sf_lp,
    aes(fill = selection_status),
    shape = 23,
    size = 4,
    colour = "grey",
    stroke = 0.8
  ) +
  
  scale_colour_manual(
    name = "Demand Status",
    values = c(
      "Covered"   = okabe_ito[4],
      "Uncovered" = okabe_ito[7]
    ),
    guide = guide_legend(order = 1)
  ) +
  
  scale_shape_manual(
    name = "Demand Status",
    values = c("Covered" = 16, "Uncovered" = 1),
    guide = "none"
  ) +
  
  scale_fill_manual(
    name = "Facility Status",
    values = c(
      "Selected"     = okabe_ito[6],
      "Not selected" = okabe_ito[8]
    ),
    guide = guide_legend(
      order = 2,
      override.aes = list(
        shape = 23, size = 4,
        colour = okabe_ito[1], stroke = 0.8
      )
    )
  ) +
  
  coord_sf(
    xlim = c(bbox_expanded_lp["xmin"], bbox_expanded_lp["xmax"]),
    ylim = c(bbox_expanded_lp["ymin"], bbox_expanded_lp["ymax"]),
    expand = FALSE
  ) +
  
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.spacing.y = unit(0.5, "cm"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey90", linewidth = 0.2),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.caption = element_text(hjust = 0.5, size = 11)
  ) +
  
  labs(
    title = "MCLP Solution (GLPK)",
    caption = paste0(
      "Facilities selected: ", sum(best_lp_solution),
      " | Demand covered: ", sum(demand_sf_lp$covered), "/", nrow(demand_sf_lp),
      " (", coverage_pct_lp, "% of total demand)"
    )
  )

# # -----------------------------
# # Distance Coverage Plot
# # -----------------------------
# 
# service_dist_plot <- ggplot(
#   data = inst_data,
#   aes(x = iter, y = service_dist)
# ) +
#   geom_line(linewidth = 1) +
#   geom_point(size = 1.8) +
#   labs(
#     x = "Iteration",
#     y = "Service distance"
#   ) +
#   theme_bw()
# 
# service_dist_plot

