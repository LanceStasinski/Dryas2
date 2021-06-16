library(sf)
library(raster)
library(dplyr)
library(spData)
library(spDataLarge)
library(tmap)
library(ggplot2)
library(rgdal)
library(elevatr)
library(sp)
library(grid)

setwd("C:/Users/istas/OneDrive/Documents/Dryas Research/Dryas 2.0")

################################################################################
#North America shapefile
################################################################################
#shapefile from https://www.naturalearthdata.com/downloads/10m-cultural-vectors/

countries = st_read('alaska bounds/ne_10m_admin_0_countries.shp')
na_countries = c('Canada', 'Mexico', 'United States of America')
na = countries[countries$NAME %in% na_countries,]
crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
na2 = st_transform(na, crs)
bounds_na = st_bbox(c(xmin = -170, xmax = -55, ymin = 15, ymax = 78), crs = crs)

#create site rectangle
x1 = -150
x2 = -144
y1 = 63.5
y2 = 65.75

rect = Polygon(cbind(c(x1,x1,x2,x2,x1),c(y1,y2,y2,y1,y1)))
rect1 = Polygons(list(rect), ID = 'A')
sp.rect = SpatialPolygons(list(rect1))
proj4string(sp.rect) = crs
df = matrix(data = c(0))
rownames(df) = "A"
rect_spp = SpatialPolygonsDataFrame(sp.rect, data = as.data.frame(df))

#map
na_map = tm_shape(na2[1], bbox = bounds_na) + tm_fill() + tm_borders()
na_map_box = na_map + tm_shape(rect_spp) + tm_borders(col = 'red', lwd = 2)
################################################################################
#Alaska shapefile
################################################################################

states = st_read('alaska bounds/cb_2018_us_state_500k.shp')
alaska = states[states$NAME == "Alaska",]
crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
st_crs(alaska) <- crs
bounds = st_bbox(c(xmin = -170, xmax = -140, ymax = 71, ymin = 55), crs = crs)

#create site rectangle
x1 = -150
x2 = -144
y1 = 63.5
y2 = 65.75

rect = Polygon(cbind(c(x1,x1,x2,x2,x1),c(y1,y2,y2,y1,y1)))
rect1 = Polygons(list(rect), ID = 'A')
sp.rect = SpatialPolygons(list(rect1))
proj4string(sp.rect) = crs
df = matrix(data = c(0))
rownames(df) = "A"
rect_spp = SpatialPolygonsDataFrame(sp.rect, data = as.data.frame(df))

#map
alaska_map = tm_shape(alaska, bbox = bounds) + tm_fill()
alaska_map_box = alaska_map + tm_shape(rect_spp) + tm_borders()

################################################################################
#Points
################################################################################

sites1 = read.csv('sites.csv')
sites = st_as_sf(sites1, coords = c("X", "Y"))
st_crs(sites) <- crs

#map
points = tm_shape(sites) + tm_dots(size = .3, shape = 18) +
  tm_text('Name', auto.placement = .7, size = 1.15)

################################################################################
#DEM
################################################################################
examp_df <- data.frame(x = seq(from = -144.0, to = -150.0, length.out = 1000),
                       y = seq(from = 63.5, to = 65.75, length.out = 1000))

# get raster
elevation_df <- get_elev_raster(examp_df, prj = crs, z = 8)
slope = terrain(elevation_df, opt = 'slope')
aspect = terrain(elevation_df, opt = 'aspect')
hill = hillShade(slope, aspect)
#map

dem = 
  tm_shape(elevation_df) + 
  tm_raster(n=6, palette = 'YlGnBu', title = "Elevation (m)") +
  tm_shape(hill) +
  tm_raster(palette = gray(0:100 / 100), n = 100, legend.show = FALSE, alpha = 0.6) +
  tm_legend(outside = F, position = c(.8023, .265), bg.color = 'white', frame = T) + 
  tm_compass(type = '4star', position = c("right","top")) +
  tm_scale_bar(breaks = c(0, 50, 100, 200), text.size = 1, position = c('left', 'bottom')) +
  tm_grid(alpha = .1) +
  tm_xlab("Longitude") +
  tm_ylab("Latitude")

################################################################################
#Final Map
################################################################################
map_final = dem + points

map_final
print(na_map_box, vp = viewport(0.76, 0.195, width = 0.21, height = 0.3))
