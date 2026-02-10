# PLOTTING FOR LSCP
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

demand_data
facility_loc
distance_matrix

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

# Assign facility selection status
facility_sf$selected <- best_bo_solution

# Coverage assignment
coverage <- A %*% best_bo_solution
demand_sf$covered <- coverage > 0

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

#  Fetch OpenStreetMap data
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
plot_lscp <- ggplot() +
  
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
    plot.caption = element_text(hjust = 0.5)
  ) +
  
  labs(
    title = "LSCP Solution: Facilities and Demand Coverage",
    caption = paste0(
      "Selected facilities: ", sum(best_bo_solution), 
      " | Covered demand: ", sum(demand_sf$covered), 
      "/", nrow(demand_sf)
    )
  )


ggsave(
  filename = "/Users/niyati/Projects:Codes/P3Compute/Facility Location/LSCP.pdf",
  plot = plot_lscp,
  device = cairo_pdf,
  width = 11, height = 10, units = "in",
  dpi = 1200, bg = "white"
)
