# Dryas topo map
# https://cran.r-project.org/web/packages/elevatr/vignettes/introduction_to_elevatr.html
# https://www.neonscience.org/dc-plot-raster-data-r

library(elevatr)
library(sp)

setwd("C:/Users/istas/OneDrive/Documents/Dryas Research/Dryas 2.0")


# Create a data.frame
examp_df <- data.frame(x = seq(from = -144.5, to = -149.5, length.out = 500), y = seq(from = 63.6, to = 65.6, length.out = 500))
examp_df <- data.frame(x = seq(from = -130.0, to = -180.0, length.out = 1000), y = seq(from = 51, to = 72, length.out = 1000)) # AKmap
# set data
prj_dd <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
# get raster
elevation_df <- get_elev_raster(examp_df, prj = prj_dd, z = 8)
# plot
plot(elevation_df, xlim=c(-180,-130))
#plot with breaks
plot(elevation_df, 
breaks = 10, 
col = terrain.colors(10),
main="10")

# 9 breaks
plot(elevation_df, 
     breaks = c(0,400,800,1200,1600,2000,3000,4000), 
     col = c("#00A600FF", "#35B800FF", "#E6E600FF", "#E8C32EFF", "#EAB64EFF", "#ECB17EFF", "#F0C9C0FF", "#F2F2F2FF"))

# 9 breaks, 760 m
plot(elevation_df, 
     breaks = c(0,400,760,1200,1600,2000,3000,4000), 
     col = c("#00A600FF", "#35B800FF", "#E6E600FF", "#E8C32EFF", "#EAB64EFF", "#ECB17EFF", "#F0C9C0FF", "#F2F2F2FF"))

# 10 breaks
plot(elevation_df, 
     breaks = c(0,300,600,900,1200,1500,1800,2100,3000,4000), 
     col = c("#00A600FF", "#35B800FF", "#85CF00FF", "#E6E600FF", "#E8C32EFF", "#EAB64EFF", "#EEB99FFF", "#F0C9C0FF", "#F2F2F2FF"))

# 11 breaks
plot(elevation_df, 
     breaks = c(0,200,400,600,800,1000,1500,2000,2500,3000,4000), 
     col = c("#00A600FF", "#19AF00FF", "#35B800FF", "#E6E600FF", "#EBB25EFF", "#EBB16EFF", "#EDB48EFF", "#EFC0AFFF", "#F1D5D0FF", "#F1E2E1FF", "#F2F2F2FF"))

terrain.colors(30)
[1] "#00A600FF" "#0CAA00FF" "#19AF00FF" "#26B300FF" "#35B800FF" "#43BD00FF" "#53C100FF" "#63C600FF" "#74CA00FF"
[10] "#85CF00FF" "#97D300FF" "#AAD800FF" "#BDDC00FF" "#D1E100FF" "#E6E600FF" "#E6D80FFF" "#E7CC1FFF" "#E8C32EFF"
[19] "#E9BB3EFF" "#EAB64EFF" "#EBB25EFF" "#EBB16EFF" "#ECB17EFF" "#EDB48EFF" "#EEB99FFF" "#EFC0AFFF" "#F0C9C0FF"
[28] "#F1D5D0FF" "#F1E2E1FF" "#F2F2F2FF"

#generate site points
sites <- as.data.frame(read.csv("sites.csv"))
site <- sites[,-1]
# Create an example SpatialPoints
examp_sp <- SpatialPoints(site, proj4string = CRS(prj_dd))

plot(examp_sp, add = TRUE, pch = 13, col = "black", 
     lwd = 2, cex = 2)

####################### AK map
library(rnaturalearth)
library(ggplot2)
world <- ne_countries(scale = "medium", returnclass = "sf")

m <- ggplot(data = world) +
  geom_sf(fill="white") +
  coord_sf(xlim = c(-170, -140), ylim = c(51, 72), expand = FALSE) + 
  theme(panel.background = element_rect(fill = "aliceblue"))
m
