######################
# collate and analyze 
# phenology data 
# from cherskii, russia
#
# ML 03/20/18
#####################
require(raster)
require(spatial)
require(rgdal)
require(plyr)

rm(list=ls())

##########################################################
# read and collate decagon NDVI data from field sites
##########################################################

setwd('L:/data_repo/field_data/viperSensor/decagon/ndvi/')

# read raw reflectances and alpha values from upward facing reference sensors
# above-canopy sensors at the low-density site need to be corrected using
# ref sensor from high-density site to acount for differences nir/r
# all other sites have painir uppward facing reference sensors

# read ndvi data, get rid of gaps in data using na.omit
d.ndvi <- na.omit(read.csv('ndvi.SRSnr.csv', header=T))

# read in raw radiance/irradiance and alpha values
# get rid of gaps in data using na.omit
d.alpha <- na.omit(read.csv('alpha.SRSni.csv', header=T))
r.630 <- na.omit(read.csv('630nm.SRSnr.csv', header=T))
r.800 <- na.omit(read.csv('800nm.SRSnr.csv', header=T))
i.630 <- na.omit(read.csv('630nm.SRSni.csv',header=T))
i.800 <- na.omit(read.csv('800nm.SRSni.csv', header=T))

# calculate red and nir reflectance
#################################
# join incident and reflected red
red <- join(r.630,i.630,by=c('doy','year','hour','site',
                             'sensorUnit','sensorZ',
                             'loggerName'),type='left')

# join incoming data from high-density with reflected data from low-density
red.cor <- join(r.630[which(r.630$loggerName=='LDF2RS'),],
                i.630[which(i.630$loggerName=='DavCnpy'),],
                by = c('doy','year','hour',
                       'sensorUnit','sensorZ'))

# get rid of extra columns
red.cor <- red.cor[,-c(16,17)]

# correct column names
colnames(red.cor) <- colnames(red)

# replace in joined data frame and remove the extra one
red <- rbind(red[which(red$loggerName!='LDF2RS'),],red.cor)
rm(red.cor)

# calculate red reflectance
red$refl <- red$X630nm.SRSnr/red$X630nm.SRSni

#################################
# join incident and reflected nir
nir <- join(r.800,i.800,by=c('doy','year','hour','site',
                             'sensorUnit','sensorZ',
                             'loggerName'),type='left')

# join incoming data from high-density with reflected data from low-density
nir.cor <- join(r.800[which(r.800$loggerName=='LDF2RS'),],
                i.800[which(i.800$loggerName=='DavCnpy'),],
                by = c('doy','year','hour',
                       'sensorUnit','sensorZ'))

# get rid of extra columns
nir.cor <- nir.cor[,-c(16,17)]

# correct column names
colnames(nir.cor) <- colnames(nir)

# replace in joined data frame and remove the extra one
nir <- rbind(nir[which(nir$loggerName!='LDF2RS'),],nir.cor)
rm(nir.cor)

# calculate nir reflectance
nir$refl <- nir$X800nm.SRSnr/nir$X800nm.SRSni

#################################
# correct ndvi for low-density with alpha values from high-density
# see for http://library.metergroup.com/Manuals/14597_SRS_Web.pdf for info on method

# select low-density overstory reflectances
r.630 <- r.630[which(r.630$loggerName=='LDF2RS'),]
r.800 <- r.800[which(r.800$loggerName=='LDF2RS'),]
# select high-density overstory alpha
d.alpha <- d.alpha[which(d.alpha$loggerName=='DavCnpy'),]

# join all three
ref <- join(r.630,r.800,by=c('doy','year','hour','site',
                             'sensorUnit','sensorLoc','sensorZ',
                             'loggerName','sensorUID'))
ref <- join(ref,d.alpha,by=c('doy','year','hour'))

# calculate corrected decagon ndvi = (alpha*nir-r)/(alpha*nir+r)
ref$ndvi.SRSnr <- (ref$alpha.SRSni*ref$X800nm.SRSnr-ref$X630nm.SRSnr)/
                  (ref$alpha.SRSni*ref$X800nm.SRSnr+ref$X630nm.SRSnr)

# pull out relevant variables
ref <- ref[,c(1:3,25,9,12)]

#join with original ndvi vars from that logger, with ndvi omitted
ref <- join(ref,d.ndvi[which(d.ndvi$loggerName=='LDF2RS'),-4])

# add back to orignial ndvi data frame
d.ndvi <- rbind(d.ndvi[which(d.ndvi$loggerName!='LDF2RS'),],ref)

# write corrected ndvi file to GitHub repo
write.csv(d.ndvi,
          file='C:/Users/mloranty/Documents/GitHub/ness_phenology/field_ndvi.csv',
          row.names = F)
rm(ref)
##########################################################################################
#screen for precip
##########################################################################################
# read daily precip from Cherskii airport
prcp <- read.csv('L:/data_repo/field_data/viperSensor/airport/airport.csv',
                 header=T)[,-c(3,5)]

nir <- join(nir,prcp)
nir[which(nir$Pr.mm>=1),c(4,13,18)] <- NA
red <- join(red,prcp)
red[which(red$Pr.mm>=1),c(4,13,18)] <- NA
d.ndvi <- join(d.ndvi,prcp)
d.ndvi$ndvi.SRSnr[which(d.ndvi$Pr.mm>=1)] <- NA
##########################################################################################
##########################################################################################
# data now ready for analyses
# get daily values - using several approaches
##########################################################################################
##########################################################################################
# subset for growing season (May 15th to Sept 15th)
d.ndvi <- d.ndvi[which(d.ndvi$doy>134 & d.ndvi$doy<259),]
#subset to include only data from 10am to 6pm
# night time data are no good, and there are artifacts around dusk
d.ndvi <- d.ndvi[which(d.ndvi$hour>=10 & d.ndvi$hour<=18),]
# screen by quantiles
q <- quantile(d.ndvi$ndvi.SRSnr,probs=c(0.025,0.975),na.rm=T)
d.ndvi$ndvi.SRSnr[which(d.ndvi$ndvi.SRSnr>q[2] | d.ndvi$ndvi.SRSnr<q[1])] <- NA
# daily maximum value
ndvi.max <- aggregate(d.ndvi$ndvi.SRSnr,
                      list(d.ndvi$doy,d.ndvi$year,d.ndvi$site,
                           d.ndvi$sensorLoc,d.ndvi$sensorUID),
                      max)

# value at local solar noon
ndvi.noon <- d.ndvi[which(d.ndvi$hour==14),]

# mean value from 12-2pm daily
ndvi.mid <- d.ndvi[which(d.ndvi$hour >= 12 & d.ndvi$hour <= 14),]
ndvi.avg <- aggregate(ndvi.mid$ndvi.SRSnr,
                      list(ndvi.mid$doy,ndvi.mid$year,ndvi.mid$site,
                           ndvi.mid$sensorLoc, ndvi.mid$sensorUID),
                      mean)
rm(ndvi.mid)
# set column names
colnames(ndvi.avg) <- c('doy','year','stand','sensorLoc','sensorUID','ndvi')
colnames(ndvi.max) <- c('doy','year','stand','sensorLoc','sensorUID','ndvi')

# choose a dataset for analyses
ndvi <- ndvi.avg

# create subsets with unique sensor IDs
# 8 & 9 = overstory high density
# 18 & 19 = understory high density
# 12 & 13 = overstory low density
# 22 & 23 = understory low density

oh1.16 <- which(ndvi$year==2016 & 
                ndvi$sensorUID==8)
oh2.16 <- which(ndvi$year==2016 & 
                ndvi$sensorUID==9)
ol1.16 <- which(ndvi$year==2016 & 
                ndvi$sensorUID==12)
ol2.16 <- which(ndvi$year==2016 & 
                ndvi$sensorUID==13)
uh1.16 <- which(ndvi$year==2016 & 
                  ndvi$sensorUID==18)
uh2.16 <- which(ndvi$year==2016 & 
                  ndvi$sensorUID==19)
ul1.16 <- which(ndvi$year==2016 & 
                  ndvi$sensorUID==22)
ul2.16 <- which(ndvi$year==2016 & 
                  ndvi$sensorUID==23)

oh1.17 <- which(ndvi$year==2017 & 
                  ndvi$sensorUID==8)
oh2.17 <- which(ndvi$year==2017 & 
                  ndvi$sensorUID==9)
ol1.17 <- which(ndvi$year==2017 & 
                  ndvi$sensorUID==12)
ol2.17 <- which(ndvi$year==2017 & 
                  ndvi$sensorUID==13)
uh1.17 <- which(ndvi$year==2017 & 
                  ndvi$sensorUID==18)
uh2.17 <- which(ndvi$year==2017 & 
                  ndvi$sensorUID==19)
ul1.17 <- which(ndvi$year==2017 & 
                  ndvi$sensorUID==22)
ul2.17 <- which(ndvi$year==2017 & 
                  ndvi$sensorUID==23)

# function to calculate running mean
rnmn <- function(x){for(i in 6:length(x))
                    {
                      mean(x[i-5:i])
                    }
  
}
#declare plot vars
xl <- c(130,265)
yl <- c(0,0.8)
pw <- 0.75
lw <- 3

#############
# make a four panel plot
tiff(file='C:/Users/mloranty/Documents/GitHub/ness_phenology/field_ndvi.tiff',
     width=10,height=10,units="in",res=300,compression="lzw",bg="white")
par(mfcol=c(2,2))

#plot  2016 overstory data
plot(ndvi$doy[oh1.16],ndvi$ndvi[oh1.16],
     pch=1,main='2016',
     col='black',
     lwd=pw,
     xlim = xl,
     ylim = yl,
     xlab = 'DOY',
     ylab='NDVI')
lines(ndvi$doy[oh1.16],filter(ndvi$ndvi[oh1.16],rep(1,5),method="convolution",sides=2)/5,
      lwd=lw,
      col='black',
      lty='dotted')
# draw grid lines
grid()
points(ndvi$doy[oh2.16],ndvi$ndvi[oh2.16],
       pch=16,
       col='black',
       lwd=pw)
lines(ndvi$doy[oh2.16],filter(ndvi$ndvi[oh2.16],rep(1,5),method="convolution",sides=2)/5,
      col='black',
      lwd=lw)
#
points(ndvi$doy[ol1.16],ndvi$ndvi[ol1.16],
       pch=16,
       col='blue',
       lwd=pw)
lines(ndvi$doy[ol1.16],filter(ndvi$ndvi[ol1.16],rep(1,5),method="convolution",sides=2)/5,
      col='blue',
      lwd=lw)

points(ndvi$doy[ol2.16],ndvi$ndvi[ol2.16],
       pch=1,
       col='blue',
       lwd=pw)
lines(ndvi$doy[ol2.16],filter(ndvi$ndvi[ol2.16],rep(1,5),method="convolution",sides=2)/5,
      col='blue',
      lwd=lw,
      lty='dotted')
legend('topleft',c('High Density','Low Density'),fill=c('black','blue'),bty='n')
legend('top','Overstory',bty='n')
###########################################################################
#par(mfcol=c(2,2))
plot(ndvi$doy[uh1.16],ndvi$ndvi[uh1.16],
     pch=1,
     col='black',
     lwd=pw,
     xlim = xl,
     ylim = yl,
     xlab = 'DOY',
     ylab='NDVI')
lines(ndvi$doy[uh1.16],filter(ndvi$ndvi[uh1.16],rep(1,5),method="convolution",sides=2)/5,
      lwd=lw,
      col='black',
      lty='dotted')
# draw grid lines
grid()
points(ndvi$doy[uh2.16],ndvi$ndvi[uh2.16],
       pch=16,
       col='black',
       lwd=pw)
lines(ndvi$doy[uh2.16],filter(ndvi$ndvi[uh2.16],rep(1,5),method="convolution",sides=2)/5,
      col='black',
      lwd=lw)
points(ndvi$doy[ul1.16],ndvi$ndvi[ul1.16],
       pch=16,
       col='blue',
       lwd=pw)
lines(ndvi$doy[ul1.16],filter(ndvi$ndvi[ul1.16],rep(1,5),method="convolution",sides=2)/5,
      col='blue',
      lwd=lw)

points(ndvi$doy[ul2.16],ndvi$ndvi[ul2.16],
       pch=1,
       col='blue',
       lwd=pw)
lines(ndvi$doy[ul2.16],filter(ndvi$ndvi[ul2.16],rep(1,5),method="convolution",sides=2)/5,
      col='blue',
      lwd=lw,
      lty='dotted')
legend('top','Understory',bty='n')
##################### 2017 ##############################
#plot  2017 overstory data
plot(ndvi$doy[oh1.17],ndvi$ndvi[oh1.17],
     pch=1,main='2017',
     col='black',
     lwd=pw,
     xlim = xl,
     ylim = yl,
     xlab = 'DOY',
     ylab='NDVI')
lines(ndvi$doy[oh1.17],filter(ndvi$ndvi[oh1.17],rep(1,5),method="convolution",sides=2)/5,
      lwd=lw,
      col='black',
      lty='dotted')
# draw grid lines
grid()
points(ndvi$doy[oh2.17],ndvi$ndvi[oh2.17],
       pch=16,
       col='black',
       lwd=pw)
lines(ndvi$doy[oh2.17],filter(ndvi$ndvi[oh2.17],rep(1,5),method="convolution",sides=2)/5,
      col='black',
      lwd=lw)
#
points(ndvi$doy[ol1.17],ndvi$ndvi[ol1.17],
       pch=16,
       col='blue',
       lwd=pw)
lines(ndvi$doy[ol1.17],filter(ndvi$ndvi[ol1.17],rep(1,5),method="convolution",sides=2)/5,
      col='blue',
      lwd=lw)

points(ndvi$doy[ol2.17],ndvi$ndvi[ol2.17],
       pch=1,
       col='blue',
       lwd=pw)
lines(ndvi$doy[ol2.17],filter(ndvi$ndvi[ol2.17],rep(1,5),method="convolution",sides=2)/5,
      col='blue',
      lwd=lw,
      lty='dotted')
legend('top','Overstory',bty='n')
## plot understory 2017##
plot(ndvi$doy[uh1.17],ndvi$ndvi[uh1.17],
     pch=1,
     col='black',
     lwd=pw,
     xlim = xl,
     ylim = yl,
     xlab = 'DOY',
     ylab='NDVI')
lines(ndvi$doy[uh1.17],filter(ndvi$ndvi[uh1.17],rep(1,5),method="convolution",sides=2)/5,
      lwd=lw,
      col='black',
      lty='dotted')
# draw grid lines
grid()
points(ndvi$doy[uh2.17],ndvi$ndvi[uh2.17],
       pch=16,
       col='black',
       lwd=pw)
lines(ndvi$doy[uh2.17],filter(ndvi$ndvi[uh2.17],rep(1,5),method="convolution",sides=2)/5,
      col='black',
      lwd=lw)
points(ndvi$doy[ul1.17],ndvi$ndvi[ul1.17],
       pch=16,
       col='blue',
       lwd=pw)
lines(ndvi$doy[ul1.17],filter(ndvi$ndvi[ul1.17],rep(1,5),method="convolution",sides=2)/5,
      col='blue',
      lwd=lw)

points(ndvi$doy[ul2.17],ndvi$ndvi[ul2.17],
       pch=1,
       col='blue',
       lwd=pw)
lines(ndvi$doy[ul2.17],filter(ndvi$ndvi[ul2.17],rep(1,5),method="convolution",sides=2)/5,
      col='blue',
      lwd=lw,
      lty='dotted')
legend('top','Understory',bty='n')
dev.off()


# write csv files for output
write.csv(ndvi.avg[which(ndvi.avg$sensorLoc=='overstory'),],
          file='C:/Users/mloranty/Documents/GitHub/ness_phenology/field_ndvi_avg.csv',
          row.names = F)

write.csv(ndvi.noon[which(ndvi.noon$sensorLoc=='overstory'),],
          file='C:/Users/mloranty/Documents/GitHub/ness_phenology/field_ndvi_noon.csv',
          row.names = F)

write.csv(ndvi.max[which(ndvi.max$sensorLoc=='overstory'),],
          file='C:/Users/mloranty/Documents/GitHub/ness_phenology/field_ndvi_max.csv',
          row.names = F)



##########################################################
# read and extract planet NDVI data
##########################################################

setwd('L:/data_repo/gis_data/planet_cherskii/PlanetScope_4band_with_SR/')

