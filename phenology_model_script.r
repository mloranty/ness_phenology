##############################################
##script for fitting curves and analysis    ##
##############################################
library(plyr)
library(rjags)
library(coda)
library(mcmcplots)

###############################################
###############read in data####################
###############################################
#cubesat data
datP <- read.csv("z:\\projects\\ness_phenology\\planet_veg_class_sample.csv")
modDI <- "z:\\projects\\ness_phenology\\sample_model\\run1"


###############################################
###############set up data for run ############
###############################################
###now look at data plot to make sure this curve makes sense
plot(datP$doy,datP$ndvi)
min(datP$ndvi)
abline(v=193)
#day 193 does not look reliable, clouds? exclude for now

datP1 <- datP[floor(datP$doy)!=193,]
plot(datP1$doy,datP1$ndvi)


#organize a table of IDs, just start with vegetation

vegDF <- unique(data.frame(veg=datP$veg))

vegDF$vegID <- seq(1,dim(vegDF)[1])

#join back into the dataframe
datP2 <- join(datP1,vegDF, by="veg", type="left")


#normalize day timeframe
doyMin <- 135
doyMax <- 273

datP2$doyN <- (datP2$doy-doyMin)/(doyMax-doyMin)

plot(datP2$doyN, datP2$ndvi)


###############################################
###############set up model run    ############
###############################################

datalist <- list(Nobs=dim(datP2)[1], ndvi=datP2$ndvi, vegID=datP2$vegID,doyN=datP2$doyN,
				Nveg=dim(vegDF)[1])

parms <- c("ndvi.rep", "sig.ndvi","base","Amp","slopeG","slopeB","G1","B1", "half.G","half.B")

samp.modI<-jags.model(file="c:\\Users\\hkropp\\Documents\\GitHub\\ness_phenology\\phenology_model_code.r",
						data=datalist,
						n.adapt=5000,
						n.chains=3)

samp.sample <- coda.samples(samp.modI,variable.names=parms,
                       n.iter=60000, thin=30)	
			   
mcmcplot(samp.sample, parms=c("sig.ndvi","base","Amp","slopeG","slopeB","G1","B1"),
			dir=paste0(modDI,"\\history"))	
					   
mod.out <- summary(samp.sample)

write.table(mod.out$statistics,paste0(modDI,"\\comp_mod_stats.csv"),
			sep=",",row.names=TRUE)
write.table(mod.out$quantiles,paste0(modDI,"\\comp_mod_quant.csv"),
			sep=",",row.names=TRUE)



###############################################
##############look at curves###################
###############################################

#negative
phenC1 <- function(amp,base,c1,slope1,doy){

	base+(amp/((1+exp(slope1*(doy-c1)))))

}
#positive
phenC2 <- function(amp,base,c1,slope1,doy){

	base+(amp/((1+exp(slope1*(c1-doy)))))

}

#both
phenC3 <- function(amp,base,c1,slope1,c2,slope2,doy){

	base+(amp/((1+exp(slope1*(c1-doy)))*(1+exp(slope2*(doy-c2)))))

}


xval <- seq(0,1, by=.1)
yval1 <- phenC1(.5,.3,.4,10,xval)

yval3 <- phenC3(.5,.3,0.4,2,.6,10,xval)

plot(xval,yval1, type="l", lwd=2, col="cornflowerblue")
points(xval,yval3, type="l",lwd=2, col="tomato3")
abline(v=.4)
abline(v=.6)

###################################################