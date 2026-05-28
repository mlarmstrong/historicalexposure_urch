# set wd
setwd("~/Desktop/urbanurchins")
# packages
library(sf)
library(dplyr)
library(ggrepel)
library(ggspatial)
library(leaflet)
library(ggplot2)
library(maps) 
library(tools)
library(rnaturalearth)
library(forcats)
library(lubridate)
library(tidyr)
library(terra)
library(tidyterra)
library(maptiles)

##################

# create maps
#############
# prelim visualization of sites of interest in google maps
m <- leaflet()
m <- addTiles(m)
m <- addCircleMarkers(m, long=long, lat=lat, radius =2, opacity = 1, 
                      label = sites$name, labelOptions = labelOptions(noHide = T))
m 

# set object
world <- ne_countries(scale = "medium", returnclass = "sf")
sf_use_s2(FALSE)
world_crop <- st_crop(world, c(xmin =-117.2 , xmax = 110, ymin = -60, ymax = 60))
class(world)
class(world_crop)

# add state borders for the sake of clarity
states <- st_as_sf(map("state", plot = FALSE, fill = TRUE))
head(states)
states <- cbind(states, st_coordinates(st_centroid(states)))
states$ID <- toTitleCase(states$ID)
head(states)

#urban and nonurban larvae sites####
larv<-read.csv("ch2_sites.csv", header=TRUE, sep=",")

#old map
  ggplot(data = world) +
  theme_bw() + 
  geom_sf(data = world_crop, fill = 'antiquewhite1') +
  geom_sf(data = states, fill = 'antiquewhite1') +
  geom_point(data = larv, aes(x = Longitude, y = Latitude, color = dev), size = 5) +
 # geom_text_repel(data = larv, aes(x = Longitude, y = Latitude, label = Site),
#                  size = 5, nudge_x = c(-0.5, 0.5), fontface = "bold") +
  coord_sf(xlim = c(-119, -117.5), ylim = c(33.4, 34), expand = FALSE) +
  scale_color_manual(values = c("urban" = "#5C5649", "nonurban" = "#99D1EC")) +
  theme(
    plot.title = element_text(size = 24),
    panel.grid.major = element_line(color = gray(0.5), linetype = "dashed", size = 0.5),
    panel.background = element_rect(fill = "aliceblue"),
    legend.position = 'none'
  )
  
  
##fancier map
  # Define your area as an sf object
bbox_sf <- st_bbox(c(xmin = -119, ymin = 33.4, xmax = -118.5, ymax = 34),
                     crs = st_crs(4326)) |> st_as_sfc()
  
# Download tiles (no API key needed)
tiles <- get_tiles(bbox_sf, provider = "OpenStreetMap.HOT",zoom = 10,crop = TRUE)
#"Esri.NatGeoWorldMap" also good

ch2map<-
ggplot() +
  geom_spatraster_rgb(data = tiles) +   # basemap tiles
  geom_sf(data = states, fill = NA, color = "white", linewidth = 0.4) +
  geom_point(data = larv, 
             aes(x = Longitude, y = Latitude, color = dev), 
             size = 5) +
    annotation_scale(             
    location = "bl",                     
    width_hint = 0.3,                       
    text_cex = 0.8,
    bar_cols = c("black", "white")
  ) +
  annotation_north_arrow(                 
    location = "bl",                        
    pad_y = unit(0.5, "cm"),                 
    which_north = "true",
    style = north_arrow_fancy_orienteering() # options below
  ) +
  coord_sf(xlim = c(-118.6, -117.5), ylim = c(33.4, 34), expand = FALSE) +
  scale_color_manual(values = c("urban" = "#5C5649", "nonurban" = "#99D1EC")) +
  theme_bw() +
  theme(
    plot.title       = element_text(size = 24),
    panel.grid.major = element_line(color = gray(0.5), linetype = "dashed", size = 0.5),
    legend.position  = "none"
  )

ggsave("histexp_map.png", ch2map, width=15, height=10, units = "cm")
