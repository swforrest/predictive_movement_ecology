#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#| warning=FALSE

library(tidyverse)
packages <- c("amt", "sf", "terra", "RColorBrewer", "leaflet")
walk(packages, require, character.only = T)

#
#
#
#
#

buffalo_data <- read_csv("data/buffalo.csv") 

# remove individuals that have poor data quality or less than about 3 months of data. 
# The "2014.GPS_COMPACT copy.csv" string is a duplicate of ID 2024, so we exclude it
buffalo_data <- buffalo_data %>% filter(!node %in% c("2014.GPS_COMPACT copy.csv", 
                                           2029, 2043, 2265, 2284, 2346))

buffalo_data <- buffalo_data %>%  
  group_by(node) %>% 
  arrange(DateTime, .by_group = T) %>% 
  distinct(DateTime, .keep_all = T) %>% 
  arrange(node) %>% 
  mutate(ID = node)

buffalo_clean <- buffalo_data[, c(12, 2, 4, 3)]
colnames(buffalo_clean) <- c("id", "time", "lon", "lat")
attr(buffalo_clean$time, "tzone") <- "Australia/Queensland"
head(buffalo_clean)
tz(buffalo_clean$time)

buffalo_ids <- unique(buffalo_clean$id)

#
#
#
#
#
#
#
#| label: make_track
#| code-summary: "Create a trajectory object"

buffalo_all <- buffalo_clean %>% mk_track(id = id,
                                           lon,
                                           lat, 
                                           time, 
                                           all_cols = T,
                                           crs = 4326) %>% 
  transform_coords(crs_to = 3112, crs_from = 4326) # Transformation to GDA94 / 
# Geoscience Australia Lambert (https://epsg.io/3112)

#
#
#
#
#

buffalo_all %>%
  ggplot(aes(x = x_, y = y_, colour = id)) +
  geom_point(alpha = 0.5, size = 0.1) + 
  coord_fixed() +
  scale_x_continuous("Easting (m)") +
  scale_y_continuous("Northing (m)") +
  scale_colour_viridis_d() +
  theme_classic() +
  theme(legend.position = "right") 

# ggsave("outputs/data_prep/buffalo_djelk_map.png",
#        width = 150, height = 150, units = "mm",  dpi = 600)

#
#
#
#
#
#| message: false
#| warning: false

# Use the original longitude/latitude data for mapping
leaflet(data = buffalo_clean) %>%
  addProviderTiles(providers$Esri.WorldImagery) %>%
  addCircleMarkers(
    lng = ~lon, lat = ~lat,
    color = "red", radius = 2, opacity = 0.5,
    popup = ~paste("ID:", id, "<br>", "Time:", time)
  ) %>%
  addLayersControl(
    baseGroups = c("Satellite"),
    options = layersControlOptions(collapsed = TRUE)
  )

#
#
#
#
#

which_buffalo <- "2022" # select a single buffalo ID

buffalo_id <- buffalo_all %>% filter(id == which_buffalo)

buffalo_id %>%
  ggplot(aes(x = x_, y = y_, colour = t_)) +
  geom_point(alpha = 0.5, size = 0.1) + 
  coord_fixed() +
  scale_x_continuous("Easting (m)") +
  scale_y_continuous("Northing (m)") +
  theme_classic()

#
#
#
#
#
#
#

buffalo_id_steps <- buffalo_id %>% 
  steps()

head(buffalo_id_steps)

#
#
#
#
#
#
#

ndvi <- rast("mapping/ndvi_aug_2018.tif")
plot(ndvi, main = "NDVI August 2018")
points(buffalo_id$x_, buffalo_id$y_, col = "red", pch = 16, cex = 0.5)
ndvi

#
#
#

canopy <- rast("mapping/canopy_cover.tif")
herby <- rast("mapping/veg_herby.tif")
slope <- rast("mapping/slope_raster.tif")

spatial_covs <- c(ndvi, canopy, herby, slope)
names(spatial_covs) <- c("ndvi", "canopy", "herby", "slope")
plot(spatial_covs)

#
#
#
#
#

buffalo_id_steps <- buffalo_id_steps %>% 
  extract_covariates(spatial_covs, where = "end")

head(buffalo_id_steps)

#
#
#
#
#

buffalo_id_steps <- buffalo_id_steps %>% 
  mutate(t1_ = lubridate::with_tz(buffalo_id_steps$t1_, tzone = "Australia/Darwin"),
         t2_ = lubridate::with_tz(buffalo_id_steps$t2_, tzone = "Australia/Darwin"))

buffalo_id_steps <- buffalo_id_steps %>%
  mutate(
    # id_num = as.numeric(factor(id)), 
        #  step_id = step_id_, 
         x1 = x1_, x2 = x2_, 
         y1 = y1_, y2 = y2_, 
         t1 = t1_, 
         t1_rounded = round_date(t1_, "hour"), 
         hour_t1 = hour(t1_rounded),
         t2 = t2_, 
         t2_rounded = round_date(t2_, "hour"), 
         hour_t2 = hour(t2_rounded),
         hour = ifelse(hour_t2 == 0, 24, hour_t2),
         yday = yday(t1_),
         year = year(t1_), 
         month = month(t1_),
         sl = sl_, 
         log_sl = log(sl_), 
         ta = ta_, 
         cos_ta = cos(ta_),
         canopy_01 = canopy/100)

head(buffalo_id_steps) 

#
#
#
#
#
#
#
#
