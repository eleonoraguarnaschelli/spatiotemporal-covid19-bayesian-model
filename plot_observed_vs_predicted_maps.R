# -------------------------------------------------------------
# Observed vs predicted spatial maps for COVID-19 incidence
#
# This script:
#   1. Loads Italian regional boundaries
#   2. Imports region-level observed and predicted summaries
#      for each epidemic phase
#   3. Merges model outputs with spatial geometries
#   4. Builds observed, predicted, and difference maps
#   5. Combines all phase-specific maps into a single figure
# -------------------------------------------------------------


# -------------------------------------------------------------
# 0. WORKING ENVIRONMENT
# Set the working directory containing the prediction files and outputs
setwd("/Users/sguar/Desktop/DESKTOP ONE DRIVE/ELEONORA/STATISTICA COMPUTAZIONALE/MODELLO_LOG")

# -------------------------------------------------------------
# 1. LIBRARIES
# Load required packages for spatial data handling, manipulation,
# plotting, and multi-panel figure composition
library(sf)
library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)
library(purrr)
library(patchwork)

# -------------------------------------------------------------
# 2. ITALIAN REGIONAL GEOMETRIES
# Import the GeoJSON file containing the boundaries of Italian regions
italy_regions <- st_read("https://raw.githubusercontent.com/openpolis/geojson-italy/master/geojson/limits_IT_regions.geojson")

italy_regions <- italy_regions %>%
  mutate(
    regione = recode(reg_name,
                     "Valle d'Aosta/Vallée d'Aoste" = "Valle d'Aosta",
                     "Friuli-Venezia Giulia" = "Friuli Venezia Giulia",
                     "Trentino-Alto Adige/Südtirol" = "Trentino-Alto Adige")
  ) %>%
  filter(!regione %in% c("Sicilia", "Sardegna"))

# -------------------------------------------------------------
# 3. HELPER FUNCTION TO LOAD AND AGGREGATE PHASE-SPECIFIC RESULTS
# Read a CSV file containing observed and predicted regional means,
# harmonize region names, and aggregate values by region
load_phase_data <- function(file_path, phase_label) {
  df <- read_csv(file_path) %>%
    mutate(denominazione_regione = recode(denominazione_regione,
                                          "P.A. Trento" = "Trento",
                                          "P.A. Bolzano" = "Bolzano")) %>%
    mutate(denominazione_regione = ifelse(denominazione_regione %in% c("Trento", "Bolzano"),
                                          "Trentino-Alto Adige", denominazione_regione)) %>%
    group_by(denominazione_regione) %>%
    summarise(y_obs_mean = mean(y_obs_mean_regione, na.rm = TRUE),
              y_pred_mean = mean(y_pred_mean_regione, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(phase = phase_label)
  
  return(df)
}

# -------------------------------------------------------------
# 4. LOAD PHASE-SPECIFIC FILES AND COMBINE THEM INTO A SINGLE TABLE
file_names <- paste0("previsioni_F", 1:4, ".csv")
phase_labels <- paste0("Fase ", 1:4)
all_data <- map2_dfr(file_names, phase_labels, load_phase_data)

# -------------------------------------------------------------
# 5. MERGE MODEL OUTPUTS WITH REGIONAL GEOMETRIES
map_data <- italy_regions %>%
  left_join(all_data, by = c("regione" = "denominazione_regione"))

# -------------------------------------------------------------
# 6. COMPUTE COMMON COLOR SCALE LIMITS
# Use shared limits across phases to ensure visual comparability
map_data <- map_data %>%
  mutate(difference = y_obs_mean - y_pred_mean)

max_value <- max(c(map_data$y_obs_mean, map_data$y_pred_mean), na.rm = TRUE)
max_diff <- max(abs(map_data$difference), na.rm = TRUE)

# -------------------------------------------------------------
# 7. GENERATE OBSERVED, PREDICTED, AND DIFFERENCE MAPS FOR EACH PHASE
all_maps <- map(phase_labels, function(phase_label) {
  phase_data <- map_data %>% filter(phase == phase_label)
  
  map_obs <- ggplot(phase_data) +
    geom_sf(aes(fill = y_obs_mean), color = "white") +
    scale_fill_viridis_c(option = "C", limits = c(0, max_value), name = "Osservato", na.value = "grey90") +
    theme_void() +
    ggtitle(paste(phase_label, "- Osservato"))
  
  map_pred <- ggplot(phase_data) +
    geom_sf(aes(fill = y_pred_mean), color = "white") +
    scale_fill_viridis_c(option = "C", limits = c(0, max_value), name = "Predetto", na.value = "grey90") +
    theme_void() +
    ggtitle(paste(phase_label, "- Predetto"))
  
  map_diff <- ggplot(phase_data) +
    geom_sf(aes(fill = difference), color = "white") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                         limits = c(-max_diff, max_diff),
                         midpoint = 0, name = "O - P") +
    theme_void() +
    ggtitle(paste(phase_label, "- Differenza (firmata)"))
  
  list(map_obs, map_pred, map_diff)
}) %>% flatten()

# -------------------------------------------------------------
# 8. COMBINE ALL MAPS INTO A SINGLE MULTI-PANEL FIGURE
# Layout: 4 rows (phases) x 3 columns (observed, predicted, difference)
final_plot <- wrap_plots(all_maps, ncol = 3)

# -------------------------------------------------------------
# 9. SAVE THE FINAL FIGURE
ggsave(
  filename = "mappe_confronto_tutte_le_fasi.png",
  plot = final_plot,
  width = 15,
  height = 20
)