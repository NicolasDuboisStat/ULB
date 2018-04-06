############################################
#### Spatial modelling for the Scotland data
############################################


#####################
#### Read in the data
#####################
#### Data
setwd("C:/Users/Nina/Documents/ULB Assistant Stat/GEOMED 2017 - Workshop 1 - CARBayes Course/GEOMED 2017 - Workshop 1 - CARBayes Course")
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
library(sp)
library(CARBayes)
rownames(dat) <- dat$IZ
dat$IZ <- NULL
sp.dat <- combine.data.shapefile(data=dat, shp=shp, dbf=dbf)



##############################################################
#### Create the neighbourhood matrix W based on border sharing
##############################################################
library(spdep)
W.nb <- poly2nb(sp.dat, row.names = rownames(sp.dat@data))
W <- nb2mat(W.nb, style = "B")


#############################################
#### Fit the S.CARleroux model using CARBayes
#############################################
#### Specify the formula object
formula <- Y ~ offset(log(E)) + jsa + ethnic + no2


#### Fit the model
model <- S.CARleroux(formula=formula, family="poisson", data=sp.dat@data, W=W, 
        burnin=20000, n.sample=100000, verbose=TRUE)
print(model)


#### See what elements the model object contains
summary(model)


#############################################
#### Checking convergence of the MCMC samples
#############################################
summary(model$samples)
plot(model$samples$rho)
plot(model$samples$beta)
plot(model$samples$tau2)
     

#################################
#### Residual spatial correlation
#################################
W.list <- nb2listw(W.nb, style = "B")
moran.mc(x = residuals(model, type="pearson"), listw = W.list, nsim = 10000)


######################
#### Covariate effects
######################
exp(model$summary.results[2:4 , 1:3])


################################
#### Compute the estimated risks
################################
sp.dat@data$risk <- model$fitted.values / sp.dat@data$E


###################################################
#### Compute the posterior exceedence probabilities
###################################################
m <- nrow(model$samples$fitted)
risk <- model$samples$fitted / matrix(rep(sp.dat@data$E,m), 
                                nrow=m, byrow=T)
sp.dat@data$PEP <- as.numeric(summarise.samples(risk, exceedences=1)
                              $exceedences)


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



#### Plot the estimated risks
library(RColorBrewer)
ggplot(data = sp.dat2, aes(x=long, y=lat, goup=group, fill = risk)) + 
    geom_polygon() + 
    coord_equal() + 
    xlab("Easting (km)") + 
    ylab("Northing (km)") + 
    labs(title = "Estimated risks for respiratory disease in 2011", fill = "Risks") +  
    theme(title = element_text(size=14)) + 
    scale_fill_gradientn(colors=brewer.pal(n=9, name="Reds"))


#### Plot the PEP
ggplot(data = sp.dat2, aes(x=long, y=lat, goup=group, fill = PEP)) + 
  geom_polygon() + 
  coord_equal() + 
  xlab("Easting (km)") + 
  ylab("Northing (km)") + 
  labs(title = "Posterior probabilities the risks are greater than 1", fill = "PEP") +  
  theme(title = element_text(size=14)) + 
  scale_fill_gradientn(colors=brewer.pal(n=9, name="Reds"))

