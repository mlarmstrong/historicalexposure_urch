# set wd
setwd("~/Desktop/urbanurchins")
# packages
libraxry(sf)
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
##################
treat.colors<-(c("C"="#A6CEE3" , "NP"="#006D2C"))

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

ch2map<-
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

ggsave("histexp_map.png", ch2map, width=20, height=20, units = "cm")
