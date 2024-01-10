library(tidyverse)
library(here)
library(ebirdst)
library(rebird)
library(terra)
library(sf)
library(rnaturalearth)
library(tidyterra)

# region of interest
region <- "US-AZ" #US or US-MN (for example, any 2-letter state code will work)

region_info <- data.frame(region = region) %>%
  separate(region, into = c("country", "state"), sep = "-", remove = FALSE)

# all_sp
sp_all <- rebird::ebirdtaxonomy("species") %>%
  rename(Common.Name = comName)

# load user ebird dat
sp_user_all <- read.csv(here("Data", "ebird_1704596325086", "MyEBirdData.csv")) %>%
  select(Common.Name) %>%
  unique() %>%
  left_join(sp_all) %>%
  filter(category == 'species')

# sp_user_region <- read.csv(here("Data", "ebird_1704596325086", "MyEBirdData.csv")) %>%
#   separate(State.Province, into = c("country", "state"), sep = "-", remove = FALSE) %>%
#   filter(country == "US") %>%
#   if(region != "US"){filter(State.Province == region)} %>% 
#   select(Common.Name) %>%
#   unique() %>%
#   left_join(sp_all) %>%
#   filter(category == 'species')

# all sp in region
sp_region <- ebirdregionspecies(region) %>%
  left_join(sp_all) %>%
  drop_na(Common.Name)

# needed sp in region
sp_needed <- setdiff(sp_region$Common.Name, sp_user_all$Common.Name) %>%
  as.data.frame() %>%
  rename(Common.Name = ".") %>%
  left_join(sp_all)
  
# set_ebirdst_access_key(overwrite = TRUE, "16f8fmeqjqn2")

#Global variables
analysis_projection <- "+proj=moll +lon_0=-90 +x_0=0 +y_0=0 +ellps=WGS84 +units=km"
raster_projection <- "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m"

#focal_species_test <- focal_species_test[1:10,]

# all sp with ebst
sp_ebst <- ebirdst_runs %>%
  rename(Common.Name = common_name)

# all needed sp with ebst
sp_ebst_for_run <- inner_join(sp_ebst, sp_needed)

# download data for needed sp
species_list_paths <- sapply(sp_ebst_for_run$species_code, ebirdst_download_status, download_abundance = FALSE, download_occurrence = TRUE, pattern = "occurrence_median_27km", USE.NAMES = FALSE)

#species_df <- data.frame(sp_ebst_for_run,species_list_paths)

#load occurrence rasters for all species in species_list
occ_combined <- sapply(sp_ebst_for_run$species_code,load_raster, product = "occurrence", 
                       period = "weekly", metric = "median", 
                       resolution = "27km")

crs_rasters <- st_crs(occ_combined[[1]])

crs_target <- st_crs("ESPG:5070") 

study_area <- ne_states(iso_a2 = "US", returnclass = "sf")
study_area <- filter(study_area, !postal %in% c("HI", "AK"))
if(region != "US"){study_area <- filter(study_area, postal == region_info$state)} 
study_area <- st_transform(study_area, crs_rasters)
mapview::mapview(study_area)  

# crop and mask the occurence rasters using the vector extent
occ_crop_combined <- sapply(occ_combined, crop, y = study_area, USE.NAMES = FALSE, overwrite = TRUE)

occ_crop_combined <- sapply(occ_crop_combined, trim)

occ_crop_combined_project <- sapply(occ_crop_combined, project, y = "epsg:5070", method = "near")
                               

occ_crop_max <- lapply(occ_crop_combined_project, app, max)

# plot occurrence
occ_bins <- seq(0, 1, by = 0.1)
par(mar = c(0, 0, 0, 2), cex = 0.9)
plot(occ_crop_max[[1]], 
     breaks = occ_bins, 
     col = abundance_palette(length(occ_bins) - 1, season = "weekly"),
     axes = FALSE, box = FALSE, 
     #maxpixels = ncell(occ_crop),
     legend.width = 2, legend.shrink = 0.97)

week_list <- list()
for(i in 1:52){
test <- lapply(occ_crop_combined_project, subset, subset = i)
test_c <- rast(test)
test_s <- sum(test_c)
week_list[[i]] <- test_s
}


r_c <- rast(week_list) 
animate(r_c)
max_val <- max(minmax(r_c))


ggplot() +
  geom_spatraster(data = subset(r_c, 1)) +
  scale_fill_viridis_c(limits = c(0, max_val))


week_plots <- list()
for(i in 1:length(week_list)){
week_plot <-  ggplot() +
  geom_spatraster(data = week_list[[i]]) +
  geom_sf(data = study_area, fill = NA, color = "white") +   
  scale_fill_viridis_c(limits = c(0, max_val)) +
  labs(title = paste0("Expected lifers: ", region),
       subtitle = occ_crop_combined_project[["aplfal"]]@ptr[["names"]][i])
week_plot
ggsave(filename = here("Results", paste0(region,"_week",i,".jpg")), plot = week_plot)
week_plots[[i]] <-  week_plot 
}

library(magick)
imgs <- list.files(here("Results"), full.names=T)
img_list <- lapply(imgs, image_read)
img_joined <- image_join(img_list)

## animate at 2 frames per second
img_animated <- image_animate(img_joined, fps = 10)

image_write(image = img_animated,
            path = here("Results", "animation.gif"))

# library(av)
# list_of_frames <- list.files(here("Results"), full.names=T)
# av::av_encode_video(list_of_frames, framerate = 16,
#                     output = here("Results", "animation.mp4"))
