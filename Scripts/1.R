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
your_ebird_dat <- here("Data", "ebird_1704596325086", "MyEBirdData.csv") # path to where your personal eBird data are stored
needs_list_to_use <- "global" # set to global if you want to map true lifers (species you haven't observed anywhere); set to "regional" if you'd like to map needs for the specified region.

# Make directories for region
mainDir <- here("Results")
dir.create(file.path(mainDir, region), showWarnings = TRUE)
outputDir <- here("Results", region)
subdirectories <- c("Expected_lifers_weekly_maps", "Expected_lifers_animation", "Possible_lifers_weekly_maps", "Possible_lifers_animation")
lapply(file.path(mainDir,region,subdirectories), function(x) if(!dir.exists(x)) dir.create(x))

# Save region info in df for late ruse
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
  filter(category == 'species')

# User regional list
sp_user_region <- read.csv(your_ebird_dat) %>%
  separate(State.Province, into = c("country", "state"), sep = "-", remove = FALSE) %>%
  filter(country == "US")
if(region != "US"){sp_user_region <- filter(sp_user_region, State.Province == region)}
sp_user_region <- sp_user_region %>%
  select(Common.Name) %>%
  unique() %>%
  left_join(sp_all) %>%
  filter(category == 'species')

# All species observed in region (by anyone)
sp_region <- ebirdregionspecies(region) %>%
  left_join(sp_all) %>%
  drop_na(Common.Name)

# Species user needs in region
if(needs_list_to_use == "global"){
sp_needed <- setdiff(sp_region$Common.Name, sp_user_all$Common.Name) %>%
  as.data.frame() %>%
  rename(Common.Name = ".") %>%
  left_join(sp_all)}

if(needs_list_to_use == "regional"){
  sp_needed <- setdiff(sp_region$Common.Name, sp_user_region$Common.Name) %>%
    as.data.frame() %>%
    rename(Common.Name = ".") %>%
    left_join(sp_all)}
  
# All species with ebst data
sp_ebst <- ebirdst_runs %>%
  rename(Common.Name = common_name)

# All needed species with ebst data
sp_ebst_for_run <- inner_join(sp_ebst, sp_needed) %>%
  filter(species_code != "laugul") # not sure why Laughing Gull is tossing an error

# Download data for needed sp and append paths to a species df
species_list_paths <- sapply(sp_ebst_for_run$species_code, ebirdst_download_status, download_abundance = FALSE, download_occurrence = TRUE, pattern = "occurrence_median_27km", USE.NAMES = FALSE)
species_df <- data.frame(sp_ebst_for_run, species_list_paths)

# Load occurrence rasters for all species in species_list
occ_combined <- sapply(sp_ebst_for_run$species_code,load_raster, product = "occurrence", 
                       period = "weekly", metric = "median", 
                       resolution = "27km")

crs_rasters <- st_crs(occ_combined[[1]])

# Vector data for region
study_area <- ne_states(iso_a2 = "US", returnclass = "sf")
if(!region %in% c("US-HI", "US-AK")){study_area <- filter(study_area, !postal %in% c("HI", "AK"))} # if region is US only mapping conttinental US
if(region != "US"){study_area <- filter(study_area, postal == region_info$state)} 
study_area <- st_transform(study_area, crs_rasters)
mapview::mapview(study_area)  

# Crop the occurrence rasters using the vector extent
occ_crop_combined <- sapply(occ_combined, crop, y = study_area, USE.NAMES = FALSE, overwrite = TRUE)
occ_crop_combined <- sapply(occ_crop_combined, trim)

# Sum occurrence probability of each species in each cell for each week to get "expected" lifer count. Each week stored as a single-layer SpatRaster in a list.
expected_lifers <- list()
for(i in 1:52){
  week_slice <- lapply(occ_crop_combined, subset, subset = i)
  week_slice <- rast(week_slice)
  week_slice <- sum(week_slice, na.rm = TRUE)
  expected_lifers[[i]] <- week_slice
}

# Reproject, mask, and trim weekly rasters for plotting
expected_lifers <- sapply(expected_lifers, project, y = "epsg:5070", method = "near")
expected_lifers <- sapply(expected_lifers, mask, mask = project(vect(study_area), y = "epsg:5070"))
expected_lifers <- sapply(expected_lifers, trim)
# plot(expected_lifers[[26]])

# Get maximum lifer count (across all cells and weeks. Needed for fill scale.
r_c <- rast(expected_lifers) 
max_val <- max(minmax(r_c))

# Generate weekly maps
bg_color <- "azure"
week_plots <- list()
for(i in 1:length(expected_lifers)){
  date <- occ_crop_combined[[1]]@cpp[["names"]][i]
  week_plot <-  ggplot() +
    geom_spatraster(data = expected_lifers[[i]]) +
    geom_sf(data = study_area, fill = NA, color = alpha("white", .4)) +   
    scale_fill_viridis_c(limits = c(0, max_val), na.value = "transparent", option = "turbo") +
    labs(title = "Lifer finder",
         subtitle = paste0(format(ymd(date), format = "%b-%d")),
         caption = paste0("Lifers mapped for: ", user,"\nData from 2022 eBird Status & Trends products. Analysis and map by Sam Safran.\nThe number of expected lifers is obtained by summing the occurrence probabilities of all lifer species at each locaiton and date.")) +
    ggthemes::theme_fivethirtyeight() +
    theme(
      rect = element_rect(linetype = 1, colour = NA),
      plot.title = element_text(size = 10, hjust = 0.05, face = "bold", color = "azure3"),
      plot.subtitle = element_text(size = 12, hjust = 0.05, face = "bold"),
      axis.title.y = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      plot.caption = element_text(size = 6, hjust = 0),
      plot.background = element_rect(fill = bg_color),
      panel.background = element_rect(fill = bg_color),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      legend.title = element_text(size = 12, face = "bold"),
      legend.direction="horizontal",
      legend.key.width = unit(.4, "inch"),
      legend.position = c(.67,1.045),
      legend.background = element_rect(colour = NA, fill = NA, linetype='solid')) +
    labs(fill = "Expected lifers")
  week_plot
  ggsave(filename = here("Results", region, "Expected_lifers_weekly_maps", paste0(region,"_",date,".jpg")), plot = week_plot, bg="white", height = 7, width = 7)
  week_plots[[i]] <-  week_plot 
}

# Generate animated gif
imgs <- list.files(here("Results", region, "Expected_lifers_weekly_maps"), full.names=T)
img_joined <- image_join(lapply(imgs, image_read))
img_animated <- image_animate(img_joined, fps = 5)
image_write(image = img_animated,
            path = here("Results", region, "Expected_lifers_animation", paste0(region,"_Expected_lifers_annual.gif")))

##############################
# Repeat for "possible" lifers

possible_occurrence_threshold <- 0.05 # minimum occurrence probabilty for a species to be considered "possible" at a given time/location.

possible_lifers <- list()
for(i in 1:52){
  week_slice <- lapply(occ_crop_combined, subset, subset = i)
  week_slice <- rast(week_slice)
  week_slice <- ifel(week_slice > possible_occurrence_threshold, 1, 0)
  week_slice <- sum(week_slice, na.rm = TRUE)
  possible_lifers[[i]] <- week_slice
}

possible_lifers <- sapply(possible_lifers, project, y = "epsg:5070", method = "near")
possible_lifers <- sapply(possible_lifers, mask, mask = project(vect(study_area), y = "epsg:5070"))
possible_lifers <- sapply(possible_lifers, trim)
# plot(possible_lifers[[26]])

r_c_possible <- rast(possible_lifers) 
max_val_possible <- max(minmax(r_c_possible))

bg_color = "azure"
week_plots_possible <- list()
for(i in 1:length(possible_lifers)){
  date <- occ_crop_combined[[1]]@cpp[["names"]][i]
  week_plot <-  ggplot() +
    geom_spatraster(data = possible_lifers[[i]]) +
    geom_sf(data = study_area, fill = NA, color = alpha("white", .4)) +   
    scale_fill_viridis_c(limits = c(0, max_val_possible), na.value = "transparent", option = "turbo") +
    labs(title = "Lifer finder",
         subtitle = paste0(format(ymd(date), format = "%b-%d")),
         caption = paste0("Lifers mapped for: ", user, "\nData from 2022 eBird Status & Trends products. Analysis and map by Sam Safran.\nA canididate lifer is considered `possible` if the species has a >", round(possible_occurrence_threshold*100,0),"% modeled occurrence probability at the location and date.")) +
    ggthemes::theme_fivethirtyeight() +
    theme(
      rect = element_rect(linetype = 1, colour = NA),
      plot.title = element_text(size = 10, hjust = 0.05, face = "bold", color = "azure3"),
      plot.subtitle = element_text(size = 12, hjust = 0.05, face = "bold"),
      axis.title.y = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      plot.caption = element_text(size = 6, hjust = 0),
      plot.background = element_rect(fill = bg_color),
      panel.background = element_rect(fill = bg_color),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      legend.title = element_text(size = 12, face = "bold"),
      legend.direction="horizontal",
      legend.key.width = unit(.4, "inch"),
      legend.position = c(.67,1.045),
      legend.background = element_rect(colour = NA, fill = NA, linetype='solid')) +
  labs(fill = "Possible lifers")
  week_plot
  ggsave(filename = here("Results", region, "Possible_lifers_weekly_maps", paste0(region,"_",date,".jpg")), plot = week_plot, bg="white", height = 7, width = 7)
  week_plots_possible[[i]] <-  week_plot 
}

imgs <- list.files(here("Results", region, "Possible_lifers_weekly_maps"), full.names=T)
img_joined <- image_join(lapply(imgs, image_read))
img_animated <- image_animate(img_joined, fps = 5)
image_write(image = img_animated,
            path = here("Results", region, "Possible_lifers_animation", paste0(region,"_Possible_lifers_annual.gif")))
