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

# test on mac
# setwd('/Users/mloranty/Desktop/20170611/')

# function to calculate various RS indices
# note this also accounts for scaling factor of 0.0001
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

# list all surface reflectance files
sr <- list.files(pattern='SR',
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

################################################################

# now reproject each file in order to be able to stack them for TS analyses

# list all of the NDVI files
n <- list.files(pattern=glob2rx('*NDVI.tif'),
                full.names=T,
                recursive=T)

# list all of the EVI files
e <- list.files(pattern=glob2rx('*EVI.tif'),
                full.names=T,
                recursive=T)

# the files have the same projection and origin, but different extents
# need to match extents to stack
x <- extent(raster(n[1]))

for(i in 2:length(n))
{
  y <- extent(raster(n[i]))
  x <- merge(x,y)
}
rm(y)

# now use the exent to extend files
for(i in 1:length(n))
{
  extend(raster(n[i]),x,
        value=NA,
        filename=gsub('_NDVI.tif','_NDVI_RS.tif',n[i]),
        overwrite=T)
  
  extend(raster(e[i]),x,
         value=NA,
         filename=gsub('_EVI.tif','_EVI_RS.tif',e[i]),
         overwrite=T)
}



# now stack the files and extract data to test the model

n <- list.files(pattern=glob2rx('*NDVI_RS.tif'),
                full.names=T,
                recursive=T)

# list all of the EVI files
e <- list.files(pattern=glob2rx('*EVI_RS.tif'),
                full.names=T,
                recursive=T)

# extract ndvi files
nf <- stack(n)
ef <- stack(e)

################################
# strip timestamp from filenames
ts <- strsplit(n,"/")
ts <- sapply(ts,'[[',4)
ts <- substr(ts,1,15)
ts <- strptime(ts,"%Y%m%d_%H%M%S",tz="GMT")

# now correct for local time
ts <- with_tz(ts,tzone="Asia/Magadan")

# calculate fractional julian day
fjday <- julian(ts,origin = as.POSIXct('2017-01-01', tz = "Asia/Magadan"))

#######################
# read in sample points
s <- readOGR('L:/data_repo/gis_data',
             layer='cherskiy_veg_pheno_sample')
# something wonky with original layer - but OK for now
# looks like there was a columns for elevation that we don't want/need
r <- s@coords[,1:2]
v <- SpatialPoints(r,proj4string = CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'))

##############################################
# extract the ndvi and combine date/site info
nv <- extract(nf,v)
ev <- extract(ef,v)

colnames(nv) <- fjday

veg <- c(rep('fs',2),rep('s',14),rep('f',15),rep('fp',10))

# as vector extracts by column, and here each column is an image (jday), and each row is a site
mod.dat <- as.data.frame(as.vector(nv))
mod.dat$evi <- as.vector(ev)
mod.dat$doy <- rep(fjday,each=41)
mod.dat$site <-  rep(s@data$Name,41)
mod.dat$veg <- rep(veg,41)
names(mod.dat) <- c("ndvi","evi","doy","site","veg")
write.csv(mod.dat,file="C:/Users/mloranty/Documents/GitHub/ness_phenology/planet_veg_class_sample.csv",
          row.names = F)




