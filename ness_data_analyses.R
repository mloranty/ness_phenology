######################
# create ndvi and evi 
# stacks from planet
# 4band surf refl
#
# ML 03/13/18
#####################
require(raster)
require(spatial)
require(rgdal)
require(lubridate)
rm(list=ls())

setwd('L:/data_repo/gis_data/planet_cherskii/PlanetScope_4band_with_SR/')

#########################################################################
############# define functions ###############
# also accounts for scaling factor of 0.0001
#ndvi
ndvi <- function(x){nir <- x[[4]]*0.0001
                    r <- x[[3]]*0.0001
                    vi <- (nir-r)/(nir+r)
                    return(vi)}

#evi using coefficients from MODIS/Landsat algorithms
evi <- function(x){nir <- x[[4]]*0.0001
                   r <- x[[3]]*0.0001
                   b <- x[[1]]*0.0001
                   vi <- 2.5*((nir-r)/(nir+6*r+7.5*b+1))
                   return(vi)}

# ndwi after McFeeters 1996 for water delineation
ndwi <- function(x){nir <- x[[4]]*0.0001
                    g <- x[[2]]*0.0001
                    vi <- (g-nir)/(g+nir)
                    return(vi)}

#########################################################################
############# calculate indices ###############

# list all surface reflectance files
sr <- list.files(pattern='SR.tif',
                 full.names=T,
                 recursive=T)

# create vi layers for each SR file
for(i in 1:length(sr))
{
  x <- stack(sr[i])
  n <- ndvi(x)
  writeRaster(n,
              filename=gsub('.tif','_NDVI.tif',sr[i]),
              overwrite=T)
  e <- evi(x)
  writeRaster(e,
              filename=gsub('.tif','_EVI.tif',sr[i]),
              overwrite=T)
  w <- ndwi(x)
  writeRaster(w,
              filename=gsub('.tif','_NDWI.tif',sr[i]),
              overwrite=T)
}

#########################################################################
############# crop files to study area immediately surrounding ness ###############

# list all of the NDVI files
n <- list.files(pattern=glob2rx('*NDVI.tif'),
                full.names=T,
                recursive=T)

# list all of the EVI files
e <- list.files(pattern=glob2rx('*EVI.tif'),
                full.names=T,
                recursive=T)

# list all of the NDWI files
w <- list.files(pattern=glob2rx('*NDWI.tif'),
                full.names=T,
                recursive=T)

# read roi shapefile #
roi <- readOGR('C:/Users/mloranty/Documents/GitHub/ness_phenology/data',
             layer='study_area')
#roi <- extent(c(595893,602349,7624983,7630554))

# crop all VI files to just the roi
# reclassify ndvi and evi to include only positive values
for(i in 1:length(n))
#for(i in 25:26) #process new files - need to code better
  {
    # see if the extents overlap -
    t <- tryCatch(!is.null(crop(raster(n[i]),roi)), error=function(e) return(FALSE))

    if(t==F) {
    print('no overlap')}
    else{
      #ndvi
      nc <- crop(raster(n[i]),roi,
                 filename=gsub('_NDVI.tif','_NDVI_ness.tif',n[i]),
                 overwrite=T)
      nc <- extend(nc,roi,
                   value=NA,
                   filename=gsub('_NDVI.tif','_NDVI_ness.tif',n[i]),
                   overwrite=T)
      reclassify(nc,c(-Inf,0,NA),
                filename=gsub('_NDVI.tif','_NDVI_ness.tif',n[i]),
                overwrite=T)
      #evi
      ec <- crop(raster(e[i]),roi,
              filename=gsub('_EVI.tif','_EVI_ness.tif',e[i]),
              overwrite=T)
      ec <- extend(ec,roi,
                   value=NA,
                   filename=gsub('_EVI.tif','_EVI_ness.tif',e[i]),
                   overwrite=T)
      reclassify(ec,c(-Inf,0,NA),
                 filename=gsub('_EVI.tif','_EVI_ness.tif',e[i]),
                 overwrite=T)

      #ndwi
      wc <- crop(raster(w[i]),roi,
                 filename=gsub('_NDWI.tif','_NDWI_ness.tif',w[i]),
                 overwrite=T)
      wc <- extend(wc,roi,
                   value=NA,
                   filename=gsub('_NDWI.tif','_NDWI_ness.tif',w[i]),
                   overwrite=T)
        }
  }

#########################################################################
############# create water mask from NDWI and NDVI ###############
w <- list.files(pattern=glob2rx('*NDWI_ness.tif'),
                full.names=T,
                recursive=T)

n <- list.files(pattern=glob2rx('*NDVI_ness.tif'),
                full.names=T,
                recursive=T)

#NDWI
# stack June-Sept images (too much snow in May)
wf <- stack(w[5:37])
wf.18 <- stack(w[42:59])

wf.mean <- calc(wf,mean,na.rm=T,
                filename='2017_3B_AnalyticMS_SR_NDWI_ness_MEAN.tif')
wf.mean <- calc(wf.18,mean,na.rm=T,
                filename='2018_3B_AnalyticMS_SR_NDWI_ness_MEAN.tif')

wf.max <- calc(wf,max,na.rm=T,progress=T,
                filename='2017_3B_AnalyticMS_SR_NDWI_ness_MAX.tif')
wf.max <- calc(wf.18,max,na.rm=T,
                filename='2018_3B_AnalyticMS_SR_NDWI_ness_MAX.tif')

wf.min <- calc(wf,min,na.rm=T,
               filename='2017_3B_AnalyticMS_SR_NDWI_ness_MIN.tif')
wf.min <- calc(wf.18,min,na.rm=T,
               filename='2018_3B_AnalyticMS_SR_NDWI_ness_MIN.tif')
#NDVI
nf <- stack(n)
nf.18 <- stack(n[42:59])

nf.mean <- calc(nf,mean,na.rm=T,
                filename='2017_3B_AnalyticMS_SR_NDVI_ness_MEAN.tif')
nf.mean <- calc(nf.18,mean,na.rm=T,
                filename='2018_3B_AnalyticMS_SR_NDVI_ness_MEAN.tif')

nf.max <- calc(nf,max,na.rm=T,
               filename='2017_3B_AnalyticMS_SR_NDVI_ness_MAX.tif')
nf.max <- calc(nf.18,max,na.rm=T,
               filename='2018_3B_AnalyticMS_SR_NDVI_ness_MAX.tif')

# note water mask threshold values based on visual interpretation
# used seasonal mean because if differences in coverage and hydrology
reclassify(wf.mean,c(-0.44,Inf,NA),
                   filename='2017_3B_AnalyticMS_SR_ness_water_mask.tif')
reclassify(wf.mean,c(-0.44,Inf,NA),
                   filename='2018_3B_AnalyticMS_SR_ness_water_mask.tif')
#########################################################################
############# apply water mask to all files ###############
m <- stack('2017_3B_AnalyticMS_SR_ness_water_mask.tif',
           '2018_3B_AnalyticMS_SR_ness_water_mask.tif')

# list all of the NDVI files
n <- list.files(pattern=glob2rx('*NDVI_ness.tif'),
                full.names=T,
                recursive=T)

# list all of the EVI files
e <- list.files(pattern=glob2rx('*EVI_ness.tif'),
                full.names=T,
                recursive=T)

# create index of year vectors
yr <- as.numeric(substr(e,3,6))-2016

for(i in 1:length(n))
#for(i in 19:20)
{
         mask(raster(n[i]),m[[yr[i]]],overwrite=T,filename=gsub('_NDVI_ness.tif','_NDVI_ness_h2o_mask.tif',n[i]))

         mask(raster(e[i]),m[[yr[i]]],overwrite=T,filename=gsub('_EVI_ness.tif','_EVI_ness_h2o_mask.tif',e[i]))

}

## NOW MOSAIC ALL FILES FROM THE SAME DAY
# list all of the NDVI files
n <- list.files(pattern=glob2rx('*_NDVI_ness_h2o_mask.tif'),
                full.names=T,
                recursive=T)

# list all of the EVI files
e <- list.files(pattern=glob2rx('*_EVI_ness_h2o_mask.tif'),
                full.names=T,
                recursive=T)

# strip timestamp from filenames
ts <- strsplit(n,"/")
ts <- sapply(ts,'[[',3)
ts <- substr(ts,1,15)
ts <- strptime(ts,"%Y%m%d_%H%M%S",tz="GMT")

# now correct for local time
ts <- with_tz(ts,tzone="Asia/Magadan")

d <- date(ts)
d.ag <- unique(d)

#mosaic by date
ns <- stack(n)
es <- stack(e)

for (i in 1:length(d.ag))
{
  r <- which(d == d.ag[i])
  
  if(length(r)<2) {
    print('only one file')}
  
  else{
  s <- as.list(ns[[r]])
  s$fun = mean
  s$NA.RM=TRUE
  s$filename = paste('ness_vi_mosaic/',d.ag[i],'_mosaic_3B_AnalyticMS_SR_NDVI_ness_h2o_mask.tif',sep='')
  s$overwrite=T
  do.call(mosaic,s)
  
  s <- as.list(es[[r]])
  s$fun = mean
  s$NA.RM=TRUE
  s$filename = paste('ness_vi_mosaic/',d.ag[i],'_mosaic_3B_AnalyticMS_SR_EVI_ness_h2o_mask.tif',sep='')
  s$overwrite=T
  do.call(mosaic,s)
  }
}

#########################################################################
#########################################################################
# START HERE FOR ANALYSES - EVERYTHING ABOVE IS PLANET PRE-PROCESSING
#########################################################################
#########################################################################

#########################################################################
########### CREATE CANOPY COVER MAP FROM FIELD OBSERVATIONS #############
#########################################################################


# read panchromatic data
pan <- raster("L:/data_repo/gis_data/siberia_high_resolution/chersky_area/WV02_20150308004730_103001003F9F4900_15MAR08004730-P1BS-500379646120_01_P002_u16ns3413.tif")

# read stand data
den <- read.csv("L:/projects/ness_phenology/stand_data.csv", header=T)

d <- SpatialPointsDataFrame(den[,5:4],den,
                            proj4string = CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'))

f <- readOGR("L:/projects/ness_phenology/gis",
             layer="forest")
t <- extract(pan,d,buffer=15)
t2 <- extract(pan,d,buffer=15,fun=sd)

dn <- unlist(lapply(t,FUN="mean"))
dn.sd <- unlist(lapply(t,FUN="sd"))

d.norm <- ((dn-450)/(1000-450))+1
cc <- d$cc.pct
lb <- d$agb
m1 <- nls(cc ~ a*exp(-b*d.norm),
            start=list(a=15000,b=4))
m2 <- nls(lb ~ a*exp(-b*d.norm),
          start=list(a=15000,b=4))

#########################################################################
############# create canopy cover map from field obs ###############

# just upland forest area
# f.pl <- rasterize(f,m,
#                   mask=T,
#                   filename='ness_forest_roi_planet.tif',
#                   overwrite=T)

f.pl <- raster('ness_forest_roi_planet.tif')

# reproject panchromatic to planet resolution
# pan.pl <- projectRaster(pan,m,
#                         method='bilinear',
#                         filename = 'WV02_20150308004730_planet_ness.tif',
#                         overwrite=T)

# mask to only water & forest
# pan.pl <- mask(pan.pl,f.pl,
#                filename = 'WV02_20150308004730_planet_ness_h2o_mask.tif',
#                overwrite=T)

pan.pl <- raster('WV02_20150308004730_planet_ness_h2o_mask.tif')

#cc.lm <- exp(coefficients(m2)[1]+log(pan.pl)*coefficients(m2)[2])
cc.exp <- coefficients(m1)[1]*exp(-coefficients(m1)[2]*(((pan.pl-450)/(1000-450))+1))
agb.pl <- coefficients(m2)[1]*exp(-coefficients(m2)[2]*(((pan.pl-450)/(1000-450))+1))

writeRaster(agb.pl,
            filename='agb_planet_exp_model.tif',
            overwrite=T)

cc.exp <- reclassify(cc.exp,c(100,Inf,100),
                    filename='canopy_cover_planet_exp_model.tif',
                    overwrite=T)
					
cc.bin <- reclassify(cc.exp,cbind(seq(0,90,10),seq(10,100,10),seq(5,95,10)),
                    filename='canopy_cover_planet_bin.tif',
                    overwrite=T)
#########################################################################
############# create a figure of model fits ###############
pdf(file='L:/projects/ness_phenology/figures/fig4_canopy_cover_exp_model.pdf',5,5)
par(cex=1,cex.axis=1.25,cex.lab=1.25)
plot(dn,cc,pch=16,
     ylab = "",yaxt="n",
     xlab = "DN",
     ylim=c(0,100),
     xlim=c(400,1000))
lines(400:1000,lwd=2,lty='dashed',
      coefficients(m1)[1]*exp(-coefficients(m1)[2]*((((400:1000)-450)/(1000-450))+1)))

axis(2,labels=T,tick=T,las=2)	
mtext("Canopy Cover (%)",side=2,line=3,cex=1.25)

segments(dn,d$cc.pct-d$SE,dn,d$cc.pct+d$SE,lwd=0.5)
segments(dn-dn.sd,d$cc.pct,dn+dn.sd,d$cc.pct,lwd=0.5)
dev.off()

pdf(file='L:/projects/ness_phenology/figures/larch_biomass_exp_model.pdf',5,5)
par(cex=1,cex.axis=1.25,cex.lab=1.25)
plot(dn,lb,pch=16,
     ylab = expression(paste("Larch Biomass (g ",C^-1," ",m^-2,")", sep="")),
     xlab = "DN",
     ylim=c(0,3000),
     xlim=c(400,1000))
lines(450:1000,lwd=2,lty='dashed',
      coefficients(m2)[1]*exp(-coefficients(m2)[2]*((((450:1000)-450)/(1000-450))+1)))
dev.off()

##################################################################################
################ canopy cover model fit figure #################
po <- lm(predict(m1)~cc)
cr <- lm(predict(m1)-cc~cc)
inst <- 0.05
pdf(file='L:/projects/ness_phenology/figures/fig6_canopy_cover_model_fit.pdf',5,8)
par(cex=1,cex.axis=1.25,cex.lab=1.25,mfrow=c(2,1),mar=c(0,0,0,0),oma=c(5,5,5,2))
plot(cc,predict(m1),pch=16,
     xlim=c(0,100),
     ylim=c(0,100),
     xaxt="n",
     yaxt="n",ylab="")
axis(2,labels=T,tick=T,las=2)	
axis(1,labels=F,tick=T)
mtext("Predicted Canopy Cover (%)",side=2,line=3,cex=1.25)
abline(lm(predict(m1)~cc),lwd=2,lty='dashed')
abline(0,1,lwd=0.5)
text(20,70,paste("r2 =",round(summary(po)$adj.r.squared,2)))
legend("bottomright","A",bg="white",box.col="white",cex=1.25,inset=inst)

## residuals ##
plot(cc,predict(m1)-cc,pch=16,
     xlim=c(0,100),
     ylim=c(-42,32),
     xaxt="n",
     yaxt="n")
axis(2,labels=T,tick=T,las=2)	
axis(1,labels=T,tick=T)
mtext("Residual Canopy Cover (%)",side=2,line=3,cex=1.25)
mtext("Observed Canopy Cover (%)",side=1,line=3,cex=1.25)
abline(lm(predict(m1)-cc~cc),lwd=2,lty='dashed')

text(70,20,paste("r2 =",round(summary(cr)$adj.r.squared,2)))
#points(cc,exp(predict(m2))-cc,col='gray50',pch=17)
abline(h=0,lwd=0.5)
legend("bottomright","B",bg="white",box.col="white",cex=1.25,inset=inst)
dev.off()

###
pdf(file='L:/projects/ness_phenology/figures/larch_biomass_obs_vs_pred.pdf',5,5)
par(cex=1,cex.axis=1.25,cex.lab=1.25)
plot(na.omit(lb),predict(m2),pch=16,
     xlab = expression(paste("Observed Biomass (g ",C^-1," ",m^-2,")", sep="")),
     ylab = expression(paste("Predicted Biomass (g ",C^-1," ",m^-2,")", sep="")))

abline(lm(predict(m2)~na.omit(lb)),lwd=2,lty='dashed')

text(450,2000,paste("r2 =",round(cor(predict(m2),na.omit(lb)),2)))
dev.off()
###########################################################################


# #########################################################################
# #########################################################################
# START HERE ONCE CANOPY COVER MODEL IS FINALIZED
# #########################################################################
# #########################################################################

# read canopy cover map
cc.exp <- raster('canopy_cover_planet_exp_model.tif')
cc.bin <- raster('canopy_cover_planet_bin.tif')
agb.pl <- raster('agb_planet_exp_model.tif')
# read stand data
den <- read.csv("L:/projects/ness_phenology/stand_data.csv", header=T)

d <- SpatialPointsDataFrame(den[,5:4],den,
                            proj4string = CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'))

# read vi files
e <- c(list.files(path='ness_vi_mosaic/',
                pattern=glob2rx('*EVI_ness_h2o_mask.tif'),
                full.names=T,
                recursive=T),
       list.files(pattern=glob2rx('20170831*EVI_ness_h2o_mask.tif'),
                  full.names=T,
                  recursive=T),
       list.files(pattern=glob2rx('20170615*EVI_ness_h2o_mask.tif'),
                  full.names=T,
                  recursive=T))

n <- c(list.files(path='ness_vi_mosaic/',
                pattern=glob2rx('*NDVI_ness_h2o_mask.tif'),
                full.names=T,
                recursive=T),
       list.files(pattern=glob2rx('20170831*NDVI_ness_h2o_mask.tif'),
                  full.names=T,
                  recursive=T),
       list.files(pattern=glob2rx('20170615*NDVI_ness_h2o_mask.tif'),
                  full.names=T,
                  recursive=T))

# create raster stacks
nf <- stack(n)
ef <- stack(e)

# extract data for sample stands using 15m buffer from center of stand
nv <- extract(nf[[1]],d,buffer=15,na.rm=T)
ev <- extract(ef[[1]],d,buffer=15,na.rm=T)

# make a new matrix of means
nvs <- unlist(lapply(nv,FUN="mean"))
nvsd <- unlist(lapply(nv,FUN="sd"))
evs <- unlist(lapply(ev,FUN="mean"))
evsd  <- unlist(lapply(ev,FUN="sd"))

# extract for each layer/timstep
for(i in 2:nlayers(nf))
{
  nv <- extract(nf[[i]],d,buffer=15,na.rm=T) 
  ev <- extract(ef[[i]],d,buffer=15,na.rm=T)
  
  nvs <- cbind(nvs,
               unlist(lapply(nv,FUN="mean")))
  nvsd <- cbind(nvsd,
               unlist(lapply(nv,FUN="sd")))
  evs <- cbind(evs,
               unlist(lapply(ev,FUN="mean")))
  evsd <- cbind(evsd,
               unlist(lapply(ev,FUN="sd")))
  
}

# strip timestamp from filenames
ts <- strsplit(n,"/")
ts <- sapply(ts,'[[',2)
ts <- substr(ts,1,10)
ts <- c(strptime(ts[1:19],"%Y-%m-%d",tz="GMT"),
        strptime(substr(ts[20:21],1,8),"%Y%m%d",tz="GMT"))

# correct for local time
ts <- with_tz(ts,tzone="Asia/Magadan")

# calculate fractional julian day
fjday <- julian(ts,origin = as.POSIXct('2017-01-01', tz = "Asia/Magadan"))

# set colnames with timestamp
colnames(nvs) <- as.character(ts)
colnames(evs) <- as.character(ts)

#####################################################
## subset to exclude May & sept (e.g. only June-Aug)
## this also excludes data with haze or other atmospheric effects as determined via visual inspection of imagery
s <- c(21,5,6,7,20)
ts1 <- ts[s]
fjday1 <- fjday[s]
nvs1 <- nvs[,s]
nvsd1 <- nvsd[,s]
evs1 <- evs[,s]
evsd1 <- evsd[,s]
nf1 <- nf[[s]]
ef1 <- ef[[s]]

###################
# analyze some data

############################## Figure 3 canopy cover vs. biomass ############################
z <- lm(d$agb~d$cc.pct)
pdf(file='L:/projects/ness_phenology/figures/fig3_larch_biomass_vs_canopy_cov.pdf',5,5)
par(cex=1,cex.axis=1.25,cex.lab=1.25,mar=c(5,5,4,2))
plot(d$cc.pct,d$agb,pch=16,
     xlab = "Canopy Cover (%)",
     yaxt="n",ylab="")

axis(2,labels=T,tick=T,las=2)	

mtext(expression(paste("Larch Biomass C (g C"," ",m^-2,")", sep="")),side=2,line=3.5,cex=1.25)

segments(d$cc.pct,d$agb-d$agb.se,d$cc.pct,d$agb+d$agb.se,lwd=0.5)
segments(d$cc.pct-d$SE,d$agb,d$cc.pct+d$SE,d$agb,lwd=0.5)

abline(lm(d$agb~d$cc.pct),lwd=2,lty='dashed')

text(20,2500,paste("r2 =",round(summary(z)$adj.r.squared,2)))
dev.off()

####################################################################################
# zonal stats for entire watershed
ws <- zonal(nf1,cc.bin,fun='mean',na.rm=T)
ws.sd <- zonal(nf1,cc.bin,fun='sd',na.rm=T)
colnames(ws) <- c("canopy",fjday1)

wse <- zonal(ef1,cc.bin,fun='mean',na.rm=T)
wse.sd <- zonal(ef1,cc.bin,fun='sd',na.rm=T)
colnames(wse) <- c("canopy",fjday1)


#######################################
# examine linear relationships between 
# canopy cover and NDVI/EVI
# write stats to table and plot
# significant relationships in Figure 6
#######################################

# set up tables to hold vars of interest
cc.n.obs <- matrix(0,6,5)

row.names(cc.n.obs) <- c(colnames(nvs1),"mean")
colnames(cc.n.obs) <- c('p','r2','f','int','slope')

cc.e.obs <- cc.n.obs
cc.n.mod <- cc.n.obs
cc.e.mod <- cc.n.obs

for(i in 1:5)
{
  l <- lm(nvs1[,i]~d$cc.pct)
  x <- summary(l)
  cc.n.obs[i,] <- c(x$coefficients[2,4],
                    x$adj.r.squared,
                    x$fstatistic[1],
                    x$coefficients[,1])
  
  l <- lm(evs1[,i]~d$cc.pct)
  x <- summary(l)
  cc.e.obs[i,] <- c(x$coefficients[2,4],
                    x$adj.r.squared,
                    x$fstatistic[1],
                    x$coefficients[,1])
  
  l <- lm(ws[,i+1]~ws[,1])
  x <- summary(l)
  cc.n.mod[i,] <- c(x$coefficients[2,4],
                    x$adj.r.squared,
                    x$fstatistic[1],
                    x$coefficients[,1])
  
  l <- lm(wse[,i+1]~wse[,1])
  x <- summary(l)
  cc.e.mod[i,] <- c(x$coefficients[2,4],
                    x$adj.r.squared,
                    x$fstatistic[1],
                    x$coefficients[,1])
}

# now for seasonal mean values

l <- lm(rowMeans(nvs1)~d$cc.pct)
x <- summary(l)
cc.n.obs[6,] <- c(x$coefficients[2,4],
                  x$adj.r.squared,
                  x$fstatistic[1],
                  x$coefficients[,1])

l <- lm(rowMeans(evs1)~d$cc.pct)
x <- summary(l)
cc.e.obs[6,] <- c(x$coefficients[2,4],
                  x$adj.r.squared,
                  x$fstatistic[1],
                  x$coefficients[,1])

l <- lm(rowMeans(ws[,2:6])~ws[,1])
x <- summary(l)
cc.n.mod[6,] <- c(x$coefficients[2,4],
                  x$adj.r.squared,
                  x$fstatistic[1],
                  x$coefficients[,1])

l <- lm(rowMeans(wse[,2:6])~wse[,1])
x <- summary(l)
cc.e.mod[6,] <- c(x$coefficients[2,4],
                  x$adj.r.squared,
                  x$fstatistic[1],
                  x$coefficients[,1])
rm(x)

write.csv(cc.n.obs,row.names = T,
          file='L:/projects/ness_phenology/canopy_vs_ndvi_obs.csv' )
write.csv(cc.n.mod,row.names = T,
          file='L:/projects/ness_phenology/canopy_vs_ndvi_mod.csv' )
write.csv(cc.e.obs,row.names = T,
          file='L:/projects/ness_phenology/canopy_vs_evi_obs.csv' )
write.csv(cc.e.mod,row.names = T,
          file='L:/projects/ness_phenology/canopy_vs_evi_mod.csv' )

###########################################################################
#### PLOT SEASONAL VI TRAJECTORIES ########################################
###########################################################################
######## Figure 8 - six panels with observed and modeled canopy cover #####
######## scatterplots of observed Canopy Cover and PLanet NDVI
###########################################################################

# plot vars
cl <- c('black','blue','red','orange','purple')
inst <- 0.05
wc <- 0.5

pdf(file="L:/projects/ness_phenology/figures/fig8.pdf",7,10)

par(cex=1,cex.axis=1.25,cex.lab=1.25,mar=c(0,0,0,0),mfcol=c(3,2),oma=c(5,5,5,2))

# panel A obs canopy vs. NDVI
plot(0,0,col=0,
     ylim=c(0,0.8), xlim=c(0,100),
	 xaxt = "n",
	 yaxt= "n")
axis(2,labels=T,tick=T,las=2)	
axis(1,labels=F,tick=T)
mtext("NDVI",side=2,line=3)
for(i in 1:5)
{
  points(d$cc.pct,nvs1[,i],pch=16,col=cl[i])
  segments(d$cc.pct,nvs1[,i]-nvsd[,i],d$cc.pct,nvs1[,i]+nvsd[,i],col=cl[i],lwd=wc)
  segments(d$cc.pct-d$SE,nvs1[,i],d$cc.pct+d$SE,nvs1[,i],col=cl[i],lwd=wc)
}

#### need to figure out how to plot significant regression lines
# also plot obs vs predicted canopy cover
# and planet vs landsat NDVI/EVI, and perhaps field obs of canopy cover vs. landsat ndvi
r <- which(cc.n.obs[1:5,1]<0.05)

for(i in 1:length(r))
{
     abline(a=cc.n.obs[r[i],4],b=cc.n.obs[r[i],5],col=cl[r[i]])
}

grid()
legend("bottomright","A",bg="white",box.col="white",cex=1.25,inset=inst)
legend("bottomleft", bg="white",box.col="white",ncol=3,
       c("15 June", "27 June", "26 July", "16 Aug", "31 Aug"),
       fill =cl,inset=inst,cex=1.25)



# panel B obs canopy vs. EVI

plot(0,0,col=0,
     ylim=c(0,0.43), xlim=c(0,100),
     xaxt = "n",
     yaxt="n",
     ylab = "EVI")

axis(2,labels=T,tick=T,las=2)	
axis(1,labels=F,tick=T)
mtext("EVI",side=2,line=3)

for(i in 1:5)
{
  points(d$cc.pct,evs1[,i],pch=16,col=cl[i])
  segments(d$cc.pct,evs1[,i]-evsd[,i],d$cc.pct,evs1[,i]+evsd[,i],col=cl[i],lwd=wc)
  segments(d$cc.pct-d$SE,evs1[,i],d$cc.pct+d$SE,evs1[,i],col=cl[i],lwd=wc)
}

r <- which(cc.n.obs[1:5,1]<0.05)

for(i in 1:length(r))
{
  abline(a=cc.e.obs[r[i],4],b=cc.n.obs[r[i],5],col=cl[r[i]])
}

legend("bottomright","B",bg="white",box.col="white",cex=1.25,inset=inst)

grid()

# panel C obs canopy vs mean growing season NDVI and EVI
plot(d$cc.pct,rowMeans(nvs1),pch=16,
     ylim=c(0,0.8),xlim=c(0,100),
     xlab='Canopy Cover (%)',
     yaxt="n",
     ylab='Mean Seasonal VI')
points(d$cc.pct,rowMeans(evs1),pch=17)
axis(2,labels=T,tick=T,las=2)	
mtext("Season Mean VI",side=2,line=3)
mtext("Observed Canopy Cover (%)",side=1,line=3)
## only if these results are significant ##
abline(a=cc.n.obs[6,4],b=cc.n.obs[6,5],col="black")
abline(a=cc.e.obs[6,4],b=cc.e.obs[6,5],col="black")

# not sure if these are necessary - sd for VI is conveyed in panels above
# segments(d$cc.pct-d$SE,rowMeans(nvs1[,1:4]),d$cc.pct+d$SE,rowMeans(nvs1[,1:4]))
# segments(d$cc.pct-d$SE,rowMeans(evs1[,1:4]),d$cc.pct+d$SE,rowMeans(evs1[,1:4]))
# 
# segments(d$cc.pct,rowMeans(nvs1[,1:4])-apply(nvs1[,1:4],1,sd),d$cc.pct,rowMeans(nvs1[,1:4])+apply(nvs1[,1:4],1,sd))
# segments(d$cc.pct,rowMeans(evs1[,1:4])-apply(evs1[,1:4],1,sd),d$cc.pct,rowMeans(evs1[,1:4])+apply(evs1[,1:4],1,sd))

legend("bottomleft",bty="n",
       c("NDVI","EVI"),
       pch=c(16,17),inset=inst,cex=1.25)
legend("bottomright","C",bg="white",box.col="white",cex=1.25,inset=inst)

grid()

## modeled data 
# panel D modeled canopy  NDVI 

plot(0,0,col=0,
     ylim=c(0,0.8), xlim=c(0,100),
     xaxt = "n",
     yaxt = "n")
axis(2,labels=F,tick=T)	
axis(1,labels=F,tick=T)

for(i in 1:5)
{
  polygon(c(ws[,1],rev(ws[,1])),c(ws[,i+1]-ws.sd[,i+1],rev(ws[,i+1]+ws.sd[,i+1])),
          col=rgb(t(col2rgb(cl[i],alpha=F)),alpha = 64,maxColorValue = 255),density=NULL,border=NA) 
  
  lines(ws[,1],ws[,i+1],type="b",col=cl[i],pch=16)
}

r <- which(cc.n.mod[1:5,1]<0.05)

for(i in 1:length(r))
{
  abline(a=cc.n.mod[r[i],4],b=cc.n.mod[r[i],5],col=cl[r[i]])
}

legend("bottomright","D",bg="white",box.col="white",cex=1.25,inset=inst)

grid()

####################################################
# panel E modeled canopy  EVI 

plot(0,0,col=0,
     ylim=c(0,0.43), xlim=c(0,100),
     xaxt="n",
     yaxt="n")
axis(2,labels=F,tick=T)	
axis(1,labels=F,tick=T)

for(i in 1:5)
{
  polygon(c(wse[,1],rev(wse[,1])),c(wse[,i+1]-wse.sd[,i+1],rev(wse[,i+1]+wse.sd[,i+1])),
          col=rgb(t(col2rgb(cl[i],alpha=F)),alpha = 64,maxColorValue = 255),density=NULL,border=NA) 
  
  lines(wse[,1],wse[,i+1],type="b",col=cl[i],pch=16)
}

r <- which(cc.e.mod[1:5,1]<0.05)

for(i in 1:length(r))
{
  abline(a=cc.e.mod[r[i],4],b=cc.e.mod[r[i],5],col=cl[r[i]])
}

legend("bottomright","E",bg="white",box.col="white",cex=1.25,inset=inst)

grid()

# panel F model mean season VI
plot(ws[,1],rowMeans(ws[,2:6]),pch=16,type="b",
     ylim=c(0,0.8),xlim=c(0,100),
     xlab='Canopy Cover (%)',
     yaxt="n",
     ylab='Mean Seasonal VI')
lines(wse[,1],rowMeans(wse[,2:6]),type="b",pch=17)
axis(2,labels=F,tick=T,las=2)	
mtext("Modeled Canopy Cover (%)",side=1,line=3)
# not sure if these are necessary - sd for VI is conveyed in panels above


legend("bottomleft",bty="n",
       c("NDVI","EVI"),
       pch=c(16,17),inset=inst,cex=1.25)
legend("bottomright","F",bg="white",box.col="white",cex=1.25,inset=inst)

grid()

dev.off()
###################################################################



############# compare Landsat NDVI with Berner map ###############
# Berner et al 2012 Larch Biomass Map
agb <- raster("L:/data_repo/gis_data/Berner_2011_Kolyma_fire_biomass/kolyma_landsat5_larch_AGB_gm2_2007.tif")

# mostly cloud free landsat map from 28 July 2017
ls <- stack(list.files("L:/data_repo/gis_data/landsat/LC081050122017072801T1-SC20180524140512/",pattern="sr_band",full.names = T))

lsm <- raster("L:/data_repo/gis_data/landsat/LC081050122017072801T1-SC20180524140512/LC08_L1TP_105012_20170728_20170810_01_T1_pixel_qa.tif")
lsm <- raster("L:/data_repo/gis_data/landsat/LC08_L1TP_105012_20170728_20170810_01_T1/LC08_L1TP_105012_20170728_20170810_01_T1_BQA.TIF")

## mask using clear pixel values from landsat BQA band: https://landsat.usgs.gov/collectionqualityband
ls <- mask(ls,lsm,inverse=T,maskvalue=2720,updatevalue=NA)
# calculate ndvi
ls.n <- (ls[[5]]-ls[[4]])/(ls[[5]]+ls[[4]])
ls.e <- evi(ls[[2:5]])
ls.e <- reclassify(ls.e,matrix(c(-Inf,0.1,NA,1,Inf,NA),nrow = 2,byrow = T),
                   filename="L:/data_repo/gis_data/landsat/LC081050122017072801T1-SC20180524140512/LC081050122017072801T1_EVI.tif",
                   overwrite = T)
ls.e <- projectRaster(ls.e,agb,
                      filename = "L:/data_repo/gis_data/landsat/LC081050122017072801T1-SC20180524140512/LC081050122017072801T1_EVI_aea.tif",
                      overwrite = T)

ls.w <- ndwi(ls[[2:5]])
ls.n <- reclassify(ls.n,matrix(c(-Inf,0.5,NA,1,Inf,NA),nrow = 2,byrow = T),
                   filename="L:/data_repo/gis_data/landsat/LC081050122017072801T1-SC20180524140512/LC081050122017072801T1_NDVI.tif",
                   overwrite = T)

ls.n <- projectRaster(ls.n,agb,
                      filename = "L:/data_repo/gis_data/landsat/LC081050122017072801T1-SC20180524140512/LC081050122017072801T1_NDVI_aea.tif",
                      overwrite = T)

## load processed landsat files ##
ls.n <- raster("L:/data_repo/gis_data/landsat/LC081050122017072801T1-SC20180524140512/LC081050122017072801T1_NDVI_aea.tif")
ls.e <- raster("L:/data_repo/gis_data/landsat/LC081050122017072801T1-SC20180524140512/LC081050122017072801T1_EVI_aea.tif")

#r30 <- lm(getValues(ls.n)~getValues(agb))

pl.ls.n <- projectRaster(nf1[[3]],ls.n,filename = "planet_07262018_ndvi_30m.tif",overwrite=T)
pl.ls.e <- projectRaster(ef1[[3]],ls.n,filename = "planet_07262018_evi_30m.tif",overwrite=T)


## extract landsat data for study sites ##
ls.f <- extract(ls.n,d,buffer=15,na.rm=T,small=T)
ls.fe <- extract(ls.e,d,buffer=15,na.rm=T,small=T)
###########################################################################
###########################################################################
######## Figure 9 - Landsat VI vs Canopy Cover                        #####
###########################################################################
### plot Landsat VI vs planet VA

pdf(file='L:/projects/ness_phenology/figures/fig9_landsat_vs_planet_VI.pdf',5,5)
par(cex=1,cex.axis=1.25,cex.lab=1.25
    #mfrow=c(2,1),mar=c(0,0,0,0),oma=c(5,5,5,2)
    )
plot(getValues(ls.n),getValues(pl.ls),,pch=16,
     xlim=c(0.4,0.9),
     ylim=c(0.4,0.9),
     xlab="Landsat NDVI",
     yaxt="n",ylab="Planet NDVI")
axis(2,labels=T,tick=T,las=2)	
#axis(1,labels=F,tick=T)
abline(0,1,lty="dashed")

dev.off()
### plot Landsat VI vs canopy cover
pdf(file="L:/projects/ness_phenology/figures/figure9.pdf",6,6)
plot(d$cc.pct,ls.f,pch=16,
     ylim=c(0,0.9),xlim=c(0,100),
     xlab='Observed Canopy Cover (%)',
     yaxt="n",
     ylab='Landsat VI')
points(d$cc.pct,ls.fe,pch=17)
segments(d$cc.pct-d$SE,unlist(ls.f),d$cc.pct+d$SE,unlist(ls.f))
segments(d$cc.pct-d$SE,unlist(ls.fe),d$cc.pct+d$SE,unlist(ls.fe))
axis(2,labels=T,tick=T,las=2)	

abline(summary(lm(unlist(ls.fe)~d$cc.pct)),lty="dashed",lwd=2)
summary(lm(unlist(ls.f)~d$cc.pct))

legend("bottomleft",bty="n",
       c("NDVI","EVI"),
       pch=c(16,17),inset=inst,cex=1.25)

dev.off()

###########################################################################

plot(ls.f,nvs1[,3])



plot(d$agb,nvs1[,6],pch=16,
     ylim=c(0.55,0.85))
points(d$agb,ls.f,pch=16,col='red',
       xlim=c(0.5,0.8))

l.r <- lm(ls.f~d$cc.pct)
p.r <- lm(nvs1[,6]~d$cc.pct)
l.a <- lm(ls.f~d$agb)
p.a <- lm(nvs1[,6]~d$agb)

pdf(file='L:/projects/ness_phenology/figures/landsat_vs_planet.pdf',5,5)

par(cex=1,cex.axis=1.25,cex.lab=1.25,mar=c(5,5,4,2))
plot(ls.f,nvs1[,3],pch=16,
     xlim=c(0.6,0.85),
     ylim=c(0.6,0.85),
     xlab = "Landsat8 NDVI",
     ylab = "PlanetScope NDVI")
	 abline(0,1,lwd=2,lty='dashed')
	 dev.off()
	 
########################################################################
