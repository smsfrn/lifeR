# Using eBird Status & Trends products to map cumulative potential lifers. Note that running this for a big region (like the whole US) or for many species requires *lots* of working memory.

# Load packages
library(tidyverse)
library(here)
library(ebirdst)
library(rebird)
library(terra)
library(sf)
library(rnaturalearth)
library(rnaturalearthhires)
library(tidyterra)
library(magick)

# Set parameters
region <- "US" # "US" for whole continental US or "US-NY" for an individual state (replace "NY" with any 2-letter state code).
user <- "Sam Safran" # enter how you want to be identified in the map caption.
your_ebird_dat <- here("Data", "sam", "ebird_1705018505439", "MyEBirdData.csv") # path to where your personal eBird data are stored
needs_list_to_use <- "global" # set to "global" if you want to map true lifers (species you haven't observed anywhere); set to "regional" if you'd like to map needs for the specified region.
resolution <- "9km" # "3km", "9km", or "27km"

# Make directories for user & region
user_file <- tolower(str_replace(user, " ", ""))
mainDir <- here("Results")
dir.create(file.path(mainDir, user_file, region, needs_list_to_use, resolution), recursive = TRUE, showWarnings = TRUE)
outputDir <- here("Results", user_file, region, needs_list_to_use, resolution)
subdirectories <- c("Weekly_maps", "Animated_map")
lapply(file.path(mainDir, user_file, region, needs_list_to_use, resolution, subdirectories), function(x) if (!dir.exists(x)) dir.create(x))

# Save region info in df for later use
region_info <- data.frame(region = region) %>%
  separate(region, into = c("country", "state"), sep = "-", remove = FALSE)

# Get full eBird taxonomy
sp_all <- rebird::ebirdtaxonomy("species") %>%
  rename(Common.Name = comName)

# User global life list
sp_user_all <- read.csv(your_ebird_dat) %>%
  select(Common.Name) %>%
  unique() %>%
  left_join(sp_all) %>%
  filter(category == "species")

# User regional list
sp_user_region <- read.csv(your_ebird_dat) %>%
  separate(State.Province, into = c("country", "state"), sep = "-", remove = FALSE) %>%
  filter(country == "US")

if (region != "US") {
  sp_user_region <- filter(sp_user_region, State.Province == region)
}

sp_user_region <- sp_user_region %>%
  select(Common.Name) %>%
  unique() %>%
  left_join(sp_all) %>%
  filter(category == "species")

# All species observed in region (by anyone)
sp_region <- ebirdregionspecies(region) %>%
  left_join(sp_all) %>%
  drop_na(Common.Name)

# Species user needs in region
if (needs_list_to_use == "global") {
  sp_needed <- setdiff(sp_region$Common.Name, sp_user_all$Common.Name) %>%
    as.data.frame() %>%
    rename(Common.Name = ".") %>%
    left_join(sp_all)
}

if (needs_list_to_use == "regional") {
  sp_needed <- setdiff(sp_region$Common.Name, sp_user_region$Common.Name) %>%
    as.data.frame() %>%
    rename(Common.Name = ".") %>%
    left_join(sp_all)
}

# All species with ebst data
sp_ebst <- ebirdst_runs %>%
  rename(Common.Name = common_name)

# All needed species with ebst data
sp_ebst_for_run <- inner_join(sp_ebst, sp_needed) %>%
  filter(!species_code %in% c("laugul", "yebsap-example")) # not sure why Laughing Gull is tossing an error

# Download data for needed sp and append paths to a species df
species_list_paths <- sapply(sp_ebst_for_run$species_code, ebirdst_download_status, download_abundance = FALSE, download_occurrence = TRUE, pattern = paste0("occurrence_median_", resolution), USE.NAMES = FALSE)
species_df <- data.frame(sp_ebst_for_run, species_list_paths)

# Load occurrence rasters for all species in species_list
occ_combined <- sapply(sp_ebst_for_run$species_code, load_raster, product = "occurrence", period = "weekly", metric = "median", resolution = resolution)

# Vector data for region
study_area <- ne_states(iso_a2 = "US", returnclass = "sf")
if (!region %in% c("US-HI", "US-AK")) {
  study_area <- filter(study_area, !postal %in% c("HI", "AK"))
} # if region is US only mapping conttinental US
if (region != "US") {
  study_area <- filter(study_area, postal == region_info$state)
}
study_area <- st_transform(study_area, st_crs(occ_combined[[1]]))

# Crop the occurrence rasters using the vector extent
occ_crop_combined <- sapply(occ_combined, crop, y = study_area, USE.NAMES = FALSE, overwrite = TRUE)
occ_crop_combined <- sapply(occ_crop_combined, trim)

view_sp <- function(x){sp_max <- (raster::raster(max(occ_crop_combined[[x]])))
                     mapview::mapview(sp_max)}
# view_sp("chclon")

# some of the included species never have occurrence porbabilities >0 in the region. save resources by filtering them out and not processing their layers.
sp_maxvals <- lapply(occ_crop_combined, minmax) %>% sapply(max)
sp_ebst_for_run <- bind_cols(sp_ebst_for_run, as.data.frame(sp_maxvals))
sp_ebst_for_run_in_region <- filter(sp_ebst_for_run, sp_maxvals > 0) 
occ_crop_combined_ <- occ_crop_combined[names(occ_crop_combined) %in% sp_ebst_for_run_in_region$species_code]

# Define occurrence threhsold for when a species is "possible"
possible_occurrence_threshold <- 0.01 # minimum occurrence probability for a species to be considered "possible" at a given time/location.

# Sum number of "possible" species in each cell based on occurrence probability and the defined threshold. Each week stored as a single-layer SpatRaster in a list.
possible_lifers <- list()
for (i in 1:52) {
  week_slice <- lapply(occ_crop_combined, subset, subset = i)
  week_slice <- rast(week_slice)
  week_slice <- ifel(week_slice > possible_occurrence_threshold, 1, 0)
  week_slice <- sum(week_slice, na.rm = TRUE)
  possible_lifers[[i]] <- week_slice
}

# Reproject, mask, and trim weekly rasters for plotting
possible_lifers <- sapply(possible_lifers, project, y = "epsg:5070", method = "near")
possible_lifers <- sapply(possible_lifers, mask, mask = project(vect(study_area), y = "epsg:5070"))
possible_lifers <- sapply(possible_lifers, trim)

# Get maximum lifer count (across all cells and weeks). Needed for fill scale.
max_val_possible <- lapply(possible_lifers, minmax) %>% sapply(max) %>% max()

# Generate weekly maps
bg_color <- "azure"
week_plots_possible <- list()
for (i in 1:length(possible_lifers)) {
  date <- occ_crop_combined[[1]]@cpp[["names"]][i]
  if(needs_list_to_use == "global"){legend_lab <- "Your potential lifers"}
  if(needs_list_to_use == "regional"){legend_lab <- "Your regional needs"} 
  week_plot <- ggplot() +
    geom_spatraster(data = possible_lifers[[i]]) +
    geom_sf(data = study_area, fill = NA, color = alpha("white", .3)) +
    scale_fill_viridis_c(
      limits = c(0, max_val_possible), na.value = "transparent", option = "turbo",
      guide = guide_colorbar(title.position = "top", title.hjust = 0.2),
      # get the units (species) after the last value in the legend
      labels = function(x) {
        lab <- " species"
        chars <- nchar(paste0(tail(x, n = 1), lab))
        x_last <- as.character(paste0(tail(x, n = 1), lab))
        x_last_pad <- str_pad(x_last, nchar(x_last) * 2 + nchar(tail(x, n = 1)) + 1, side = "left")
        c(x[1:length(x) - 1], x_last_pad)
      }
    ) +
    labs(
      title = "Lifer finder: mapping the birds you've yet to meet",
      tag = paste0(format(ymd(date), format = "%b-%d")),
      # subtitle = "Mapping the birds you've yet to meet",
      fill = legend_lab,
      caption = paste0("Lifers mapped for: ", user, ". \nLifer analysis and map by Sam Safran.\nA candidate lifer is considered `possible` if the species has a >", round(possible_occurrence_threshold * 100, 0), "% modeled occurrence probability at the location and date.\n
Data from 2022 eBird Status & Trends products (https://ebird.org/science/status-and-trends): Fink, D., T. Auer, A. Johnston, M. Strimas-Mackey, S. Ligocki, O. Robinson,\n W. Hochachka, L. Jaromczyk, C. Crowley, K. Dunham, A. Stillman, I. Davies, A. Rodewald, V. Ruiz-Gutierrez, C. Wood. 2023. eBird Status and Trends, Data Version:\n2022; Released: 2023. Cornell Lab of Ornithology, Ithaca, New York. https://doi.org/10.2173/ebirdst.2022. This material uses data from the eBird Status and Trends\n Project at the Cornell Lab of Ornithology, eBird.org. Any opinions, findings, and conclusions or recommendations expressed in this material are those of the author(s)\nand do not necessarily reflect the views of the Cornell Lab of Ornithology.")
    ) +
    ggthemes::theme_fivethirtyeight() +
    theme(
      rect = element_rect(linetype = 1, colour = NA),
      plot.title = element_text(
        size = 10, hjust = 0, face = "plain", color = "azure3",
        margin = margin(0, 0, 15, 0)
      ),
      plot.tag = element_text(size = 10, face = "bold"),
      axis.title.y = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      plot.caption = element_text(size = 5, hjust = 0),
      plot.background = element_rect(fill = bg_color),
      panel.background = element_rect(fill = bg_color),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      legend.title = element_text(size = 10, face = "bold"),
      legend.title.align = 1,
      legend.direction = "horizontal",
      legend.key.width = unit(.34, "inch"),
      legend.key.height = unit(.06, "inch"),
      legend.position = c(.773, 1.03),
      plot.tag.position = c(0.05, .94),
      legend.background = element_rect(colour = NA, fill = NA, linetype = "solid"),
      legend.text.align = 0.5
    )
  week_plot
  ggsave(filename = here(outputDir, "Weekly_maps", paste0(region, "_", date, ".jpg")), plot = week_plot, bg = "white", height = 5.35, width = 6.6)
  week_plots_possible[[i]] <- week_plot
}

# Generate animated gif
imgs <- list.files(here(outputDir, "Weekly_maps"), full.names = T)
img_joined <- image_join(lapply(imgs, image_read))
img_animated <- image_animate(img_joined, fps = 5)
image_write(image = img_animated, path = here(outputDir, "Animated_map", paste0(region, "_Animated_map_annual_hires.gif")))
hires <- image_read(here(outputDir, "Animated_map", paste0(region, "_Animated_map_annual_hires.gif")))
lores <- image_scale(hires, geometry_size_percent(width = 38, height = NULL))
image_write(image = lores, path = here(outputDir, "Animated_map", paste0(region, "_Animated_map_annual_lores.gif")))
