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

# list all of the EVI files
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

wf.mean <- calc(wf,mean,na.rm=T,
                filename='2017_3B_AnalyticMS_SR_NDWI_ness_MEAN.tif')

wf.max <- calc(wf,max,na.rm=T,progress=T,
                filename='2017_3B_AnalyticMS_SR_NDWI_ness_MAX.tif')

wf.min <- calc(wf,min,na.rm=T,
               filename='2017_3B_AnalyticMS_SR_NDWI_ness_MIN.tif')

#NDVI
nf <- stack(n)

nf.mean <- calc(nf,mean,na.rm=T,
                filename='2017_3B_AnalyticMS_SR_NDVI_ness_MEAN.tif')

nf.max <- calc(nf,max,na.rm=T,
               filename='2017_3B_AnalyticMS_SR_NDVI_ness_MAX.tif')

# note water mask values based on visual interpretation
# used seasonal mean because if differences in coverage and hydrology
mask <- reclassify(wf.mean,c(-0.44,Inf,NA),
                   filename='2017_3B_AnalyticMS_SR_ness_water_mask.tif')

#########################################################################
############# apply water mask to all files ###############
mask <- raster('2017_3B_AnalyticMS_SR_ness_water_mask.tif')

# list all of the NDVI files
n <- list.files(pattern=glob2rx('*NDVI_ness.tif'),
                full.names=T,
                recursive=T)

# list all of the EVI files
e <- list.files(pattern=glob2rx('*EVI_ness.tif'),
                full.names=T,
                recursive=T)

for(i in 1:length(n))
{
  mask(raster(n[i]),mask,
       filename=gsub('_NDVI_ness.tif','_NDVI_ness_h2o_mask.tif',n[i]))
  mask(raster(e[i]),mask,
       filename=gsub('_EVI_ness.tif','_EVI_ness_h2o_mask.tif',e[i]))
}

#########################################################################
############# create canopy cover map from field obs ###############

# read panchromatic data
pan <- raster("L:/data_repo/gis_data/siberia_high_resolution/chersky_area/WV02_20150308004730_103001003F9F4900_15MAR08004730-P1BS-500379646120_01_P002_u16ns3413.tif")

# read stand data
den <- read.csv("L:/projects/ness_phenology/stand_data.csv", header=T)

d <- SpatialPointsDataFrame(den[,5:4],den,
                            proj4string = CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'))

f <- readOGR("L:/projects/ness_phenology/gis",
             layer="forest")
t <- extract(pan,d,buffer=15)

dn <- unlist(lapply(t,FUN="mean"))


d.norm <- ((dn-400)/(1400-400))+1
cc <- d$cc.pct
m1 <- nls(cc ~ a*exp(-b*d.norm),
            start=list(a=15000,b=4))
m2 <- lm(log(cc+0.1)~log(dn))
plot(log(cc+0.1)~log(dn))
#########################################################################
############# create canopy cover map from field obs ###############

# just upland forest area
f.pl <- rasterize(f,mask,
                  mask=T,
                  filename='ness_forest_roi_planet.tif',
                  overwrite=T)

# reproject panchromatic to planet resolution
pan.pl <- projectRaster(pan,mask,
                        method='bilinear',
                        filename = 'WV02_20150308004730_planet_ness.tif',
                        overwrite=T)

# mask to only water & forest
pan.pl <- mask(pan.pl,f.pl,
               filename = 'WV02_20150308004730_planet_ness_h2o_mask.tif',
               overwrite=T)

cc.lm <- exp(coefficients(m2)[1]+log(pan.pl)*coefficients(m2)[2])
cc.exp <- coefficients(m1)[1]*exp(-coefficients(m1)[2]*(((pan.pl-400)/(1400-400))+1))

cc.lm <- reclassify(cc.lm,c(100,Inf,100),
                    filename='canopy_cover_planet_linear_model.tif')
cc.exp <- reclassify(cc.exp,c(100,Inf,100),
                    filename='canopy_cover_planet_exp_model.tif')
############# extract subset of data for phenology modeling ###############

#evi files
e <- list.files(pattern=glob2rx('*EVI_ness.tif'),
                full.names=T,
                recursive=T)

# create raster stacks
nf <- stack(n)
ef <- stack(e)


# strip timestamp from filenames
ts <- strsplit(n,"/")
ts <- sapply(ts,'[[',4)
ts <- substr(ts,1,15)
ts <- strptime(ts,"%Y%m%d_%H%M%S",tz="GMT")

# now correct for local time
ts <- with_tz(ts,tzone="Asia/Magadan")

# calculate fractional julian day
fjday <- julian(ts,origin = as.POSIXct('2017-01-01', tz = "Asia/Magadan"))

# read in sample points for preliminary analyses
s <- readOGR('L:/projects/ness_phenology/gis',
             layer='cherskiy_veg_pheno_sample')

# something wonky with original layer - but OK for now
# looks like there was a columns for elevation that we don't want/need
r <- s@coords[,1:2]
v <- SpatialPoints(r,proj4string = CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'))


# extract the ndvi and combine date/site info
nv <- extract(nf,v)
ev <- extract(ef,v)

colnames(nv) <- fjday

# vector of vegetation types for preliminary analyses
veg <- c(rep('fs',2),rep('s',14),rep('f',15),rep('fp',10))

# create data frame for modeling
# as vector extracts by column, and here each column is an image (jday), and each row is a site
mod.dat <- as.data.frame(as.vector(nv))
mod.dat$evi <- as.vector(ev)
mod.dat$doy <- rep(fjday,each=nrow(nv))
mod.dat$site <-  rep(s@data$Name,41)
mod.dat$veg <- rep(veg,41)
names(mod.dat) <- c("ndvi","evi","doy","site","veg")

mod.dat <- na.omit(mod.dat)
# write to file
write.csv(mod.dat,file="C:/Users/mloranty/Documents/GitHub/ness_phenology/planet_veg_class_sample.csv",row.names = F)

write.csv(mod.dat,file="L:/projects/ness_phenology/planet_veg_class_sample.csv",row.names = F)

############# create data frame for full phenology model ###############


