#Map
library(rgdal)
library(sf)
library(ggplot2)
library(ggspatial)
library(ggrepel)


setwd("C:/Users/istas/OneDrive/Documents/Dryas Research/Dryas 2.0")

states = st_read('alaska bounds/cb_2018_us_state_500k.shp')
alaska = states[states$NAME == "Alaska",]
crs = st_crs(alaska)

sites1 = read.csv('sites.csv')
sites = st_as_sf(sites1, coords = c("X", "Y"))
st_crs(sites) <- crs

labels = read.csv('labels.csv')

dev.new(width = 6, height = 6, unit = 'in')
ggplot() + 
  geom_sf(data = alaska) + 
  geom_sf(data = sites, shape =19) +
  geom_text_repel(data = sites1, aes(x = X, y = Y, label = Name),
                  nudge_x = c(2, 1, -2,1,-3, -.5),
                  nudge_y = c(.5, -1, -.25, 1, 1,-1)) +
  coord_sf(xlim = c(-170,-140), ylim = c(55,73), expand = F) +
  ggtitle("Collection Locations") + 
  theme_bw() +
  annotation_scale(location = 'br') +
  annotation_north_arrow(location = 'bl', which_north = "true")+
  xlab("Longitude") + ylab("Latitude")



ggplot()+
  geom_sf(data = sites)


geom_text(data = labels, aes(X,Y, label = Name)) +