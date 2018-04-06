###############################################
#### Exploratory analysis for the Scotland data
###############################################

setwd("C:/Users/Nina/Documents/ULB Assistant Stat/GEOMED 2017 - Workshop 1 - CARBayes Course/GEOMED 2017 - Workshop 1 - CARBayes Course")
#####################
#### Read in the data
#####################
#### Data
dat <- read.csv(file="Scotland spatial data.csv")
head(dat)


#### Shapefiles
library(shapefiles)
dbf <- read.dbf(dbf.name = "ScotlandIZ.dbf")
shp <- read.shp(shp.name = "ScotlandIZ.shp")


#############################################
#### Create a SpatialPolygonsDataFrame object
#############################################
#### Create the object
#### The rownames of the data set must be contained within the first column of the dbf$dbf file
library(sp)
library(CARBayes)
rownames(dat) <- dat$IZ
dat$IZ <- NULL
head(dbf$dbf)
sp.dat <- combine.data.shapefile(data=dat, shp=shp, dbf=dbf)


#### Plot the object 
head(sp.dat@data)
plot(sp.dat)


#### Add the SMR to the object
sp.dat@data$smr <- sp.dat@data$Y / sp.dat@data$E
head(sp.dat@data)



##############################################################
#### Create the neighbourhood matrix W based on border sharing
##############################################################
library(spdep)
W.nb <- poly2nb(sp.dat, row.names = rownames(sp.dat@data))
W <- nb2mat(W.nb, style = "B")
class(W)
dim(W)



################################################
#### Compute the Moran's I statistic for the SMR
################################################
W.list <- nb2listw(W.nb, style = "B")
moran.mc(x = sp.dat@data$smr, listw = W.list, nsim = 10000)



##########################
#### Mapping using ggplot2
##########################
library(ggplot2)
library(rgeos)
library(maptools)

#### Turn the SpatialPolygonsDataFrame object into a data.frame object that contains the spaital information
sp.dat@data$id <- rownames(sp.dat@data)
temp1 <- fortify(sp.dat, region = "id")
sp.dat2 <- merge(temp1, sp.dat@data, by = "id")


#### Divide the coordinates by 1000 to get them in KM and not metres.
sp.dat2$long <- sp.dat2$long / 1000
sp.dat2$lat <- sp.dat2$lat / 1000

#### Create a basic map
ggplot(data = sp.dat2, aes(x=long, y=lat, goup=group, fill = smr)) + 
    geom_polygon()


#### Create a nicer map
ggplot(data = sp.dat2, aes(x=long, y=lat, goup=group, fill = smr)) + 
    geom_polygon() + 
    coord_equal() + 
    xlab("Easting (km)") + 
    ylab("Northing (km)") + 
    labs(title = "SMR for respiratory disease in 2011", fill = "SMR") +  
    theme(title = element_text(size=14))


#### Change the colour scale
library(RColorBrewer)
ggplot(data = sp.dat2, aes(x=long, y=lat, goup=group, fill = smr)) + 
    geom_polygon() + 
    coord_equal() + 
    xlab("Easting (km)") + 
    ylab("Northing (km)") + 
    labs(title = "SMR for respiratory disease in 2011", fill = "SMR") +  
    theme(title = element_text(size=14)) + 
    scale_fill_gradientn(colors=brewer.pal(n=9, name="Reds"))




################################
#### Add a Google map underneath
################################
#### Create a new spatial data object with a long-lat coordinate system
library(rgdal)
sp.dat3 <- sp.dat
proj4string(sp.dat3) <- CRS("+init=epsg:27700")                   
sp.dat3 <- spTransform(sp.dat3, CRS("+init=epsg:4326"))  


#### Turn the new spatial data object into a dataframe
sp.dat3@data$id <- rownames(sp.dat3@data)
temp1 <- fortify(sp.dat3, region = "id")
sp.dat4 <- merge(temp1, sp.dat3@data, by = "id")


#### Define the centre and extent of the data set and get a map
library(ggmap)
extent <- bbox(sp.dat3)
centre <- apply(extent, 1,mean)
myMap <- get_map(location=centre, maptype="roadmap", zoom=10)


#### Plot the map
ggmap(myMap) + 
    geom_polygon(data=sp.dat4, aes(x=long, y=lat, group=group, fill=smr), alpha=0.8) +
    xlab("Longitude") + 
    ylab("Latitude") + 
    labs(title = "SMR for respiratory disease in 2011", fill = "SMR") +  
    theme(title = element_text(size=14)) + 
    scale_fill_gradientn(colors=brewer.pal(n=9, name="Reds"))



