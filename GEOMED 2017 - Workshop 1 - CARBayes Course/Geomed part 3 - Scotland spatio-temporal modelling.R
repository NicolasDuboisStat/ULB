####################################################
#### Spatio-temporal modelling for the Scotland data
####################################################


#####################
#### Read in the data
#####################
#### Data
dat <- read.csv(file="Scotland space time data.csv")
head(dat)


#### Shapefiles
library(shapefiles)
shp <- read.shp(shp.name = "ScotlandIZ.shp")
dbf <- read.dbf(dbf.name = "ScotlandIZ.dbf")


#### Add the SMR
dat$smr <- dat$Y / dat$E


########################################
#### Plot the temporal trend in the data
########################################
library(ggplot2)
ggplot(dat, aes(x=jitter(year), y=smr)) +
    geom_point(colour="red") + 
    xlab("Year") + 
    scale_x_continuous(breaks=c(2002, 2005, 2009, 2012 )) + 
    ylab("SMR") + 
    ggtitle("SMR for coronary heart disease") +
    theme(text=element_text(face="bold", size=14))


#############################################
#### Create a SpatialPolygonsDataFrame object
#############################################
library(sp)
library(CARBayes)
dat.2002 <- dat[dat$year==2002, ]
rownames(dat.2002) <- dat.2002$IZ
sp.dat <- combine.data.shapefile(data=dat.2002, shp=shp, dbf=dbf)


##############################################################
#### Create the neighbourhood matrix W based on border sharing
##############################################################
library(spdep)
W.nb <- poly2nb(sp.dat, row.names = rownames(sp.dat@data))
W <- nb2mat(W.nb, style = "B")
W.list <- nb2listw(W.nb, style = "B")


########################################
#### Compute Moran's I for a sample year
########################################
moran.mc(dat$smr[dat$year==2002], listw = W.list, nsim = 10000)



###############################################
#### Fit the ST.CARanova model using CARBayesST
###############################################
library(CARBayesST)
model1 <- ST.CARanova(formula=Y~offset(log(E)), family="poisson", data=dat, 
                      W=W, burnin=10000, n.sample=50000, verbose=TRUE)
print(model1)


#############################################
#### Checking convergence of the MCMC samples
#############################################
plot(model1$samples$rho)
plot(model1$samples$beta)
plot(model1$samples$tau2)
plot(model1$sample$delta[ ,1:2])


#################################################
#### Compute the estimated temporal trend in risk
#################################################
theta.trend <- as.data.frame(array(NA, c(11,4)))
colnames(theta.trend) <- c("year", "median", "LCI", "UCI")
theta.trend$year <- 2002:2012
    for(i in 1:11)
    {
    temp <- quantile(model1$samples$beta[ ,1] + model1$samples$delta[ ,i], 
                     c(0.5, 0.025, 0.975))
    theta.trend[i,2:4] <- exp(temp)
    }


###############################################
#### Add the temporal trend to the raw SMR plot
###############################################
ggplot(dat, aes(x=jitter(year), y=smr)) +
    geom_point(colour="red") + 
    xlab("Year") + 
    scale_x_continuous(breaks=c(2002, 2005, 2009, 2012 )) + 
    ylab("SMR") + 
    ggtitle("Estimated risk for coronary heart disease") +
    theme(text=element_text(face="bold", size=14)) + 
    geom_line(mapping=aes(x=year, y=median), data=theta.trend) +
    geom_line(mapping=aes(x=year, y=LCI), colour="blue", data=theta.trend) + 
    geom_line(mapping=aes(x=year, y=UCI), colour="blue", data=theta.trend)




###############################################
#### Fit the ST.CARanova model using CARBayesST
###############################################
model2 <- ST.CARar(formula=Y~offset(log(E)), family="poisson", data=dat, 
                   W=W, burnin=10000, n.sample=50000, verbose=TRUE)
print(model2)


#############################################
#### Checking convergence of the MCMC samples
#############################################
plot(model2$samples$rho)
plot(model2$samples$beta)
plot(model2$samples$tau2)
plot(model2$sample$phi[ ,1:2])


###################################################################
#### Compute the estimated risks and add to the spatial data object
###################################################################
risk.median.all <- model2$fitted.values / dat$E
risk.2002 <- risk.median.all[dat$year==2002]
risk.2005 <- risk.median.all[dat$year==2005]
risk.2009 <- risk.median.all[dat$year==2009]
risk.2012 <- risk.median.all[dat$year==2012]

sp.dat@data$risk.2002 <- risk.2002 
sp.dat@data$risk.2005 <- risk.2005
sp.dat@data$risk.2009 <- risk.2009
sp.dat@data$risk.2012 <- risk.2012


#########################
#### Prepare for plotting
#########################
#### Load the libraries required
library(rgeos)
library(maptools)

#### Turn into a data.frame
sp.dat@data$id <- rownames(sp.dat@data)
temp1 <- fortify(sp.dat, region = "id")
sp.dat2 <- merge(temp1, sp.dat@data, by = "id")


#### Change the scale to kilometres
sp.dat2$long <- sp.dat2$long / 1000
sp.dat2$lat <- sp.dat2$lat / 1000


##################
#### Plot the maps
##################
#### Create each map separately
library(RColorBrewer)

#### Map 2002
map2002 <-  ggplot(data = sp.dat2, aes(x=long, y=lat, group=group, 
                                       fill = risk.2002)) + 
    geom_polygon() + 
    coord_equal() + 
    xlab("Easting (km)") + 
    ylab("Northing (km)") + 
    labs(title = "(A) - 2002", fill = "Proportion") +  
    theme(title = element_text(face="bold", size=14)) + 
    scale_fill_gradientn(colors=brewer.pal(n=9, name="Reds"), 
                         limits=c(0.4,2.7))

#### Map 2005
map2005 <-  ggplot(data = sp.dat2, aes(x=long, y=lat, group=group, 
                                       fill = risk.2005)) + 
    geom_polygon() + 
    coord_equal() + 
    xlab("Easting (km)") + 
    ylab("Northing (km)") + 
    labs(title = "(B) - 2005", fill = "Proportion") +  
    theme(title = element_text(face="bold", size=14)) + 
    scale_fill_gradientn(colors=brewer.pal(n=9, name="Reds"), 
                         limits=c(0.4,2.7))

#### Map 2009
map2009 <-  ggplot(data = sp.dat2, aes(x=long, y=lat, group=group, 
                                       fill = risk.2009)) + 
    geom_polygon() + 
    coord_equal() + 
    xlab("Easting (km)") + 
    ylab("Northing (km)") + 
    labs(title = "(C) - 2009", fill = "Proportion") +  
    theme(title = element_text(face="bold", size=14)) + 
    scale_fill_gradientn(colors=brewer.pal(n=9, name="Reds"), 
                         limits=c(0.4,2.7))

#### Map 2012
map2012 <-  ggplot(data = sp.dat2, aes(x=long, y=lat, group=group, 
                                       fill = risk.2012)) + 
    geom_polygon() + 
    coord_equal() + 
    xlab("Easting (km)") + 
    ylab("Northing (km)") + 
    labs(title = "(D) - 2012", fill = "Proportion") +  
    theme(title = element_text(face="bold", size=14)) + 
    scale_fill_gradientn(colors=brewer.pal(n=9, name="Reds"), 
                         limits=c(0.4,2.7))


#### Plot them together
library(gridExtra)
grid.arrange(map2002, map2005, map2009, map2012, ncol=2, nrow=2)

