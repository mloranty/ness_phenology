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
modDI <- "z:\\projects\\ness_phenology\\sample_model\\run3"
#set to 1 if runing model and not just plotting
modRun <-1

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

#remove forest shrub mix because too few observations to look at

datP1 <- datP1[datP1$veg!="fs",]


#organize a table of IDs, just start with vegetation

vegDF <- unique(data.frame(veg=datP1$veg))

vegDF$vegID <- seq(1,dim(vegDF)[1])


#join back into the dataframe
datP2 <- join(datP1,vegDF, by="veg", type="left")

#pixel vege id
pixDF <- unique(data.frame(site=datP2$site,vegID=datP2$vegID))

pixDF$pixID <- seq(1,dim(pixDF)[1])

datP3 <- join(datP2, pixDF, by=c("site","vegID"),type="left")

#normalize day timeframe
doyMin <- 135
doyMax <- 273

datP3$doyN <- (datP3$doy-doyMin)/(doyMax-doyMin)

plot(datP3$doyN, datP3$ndvi)



###############################################
###############set up model run    ############
###############################################
if(modRun==1){
datalist <- list(Nobs=dim(datP3)[1], ndvi=datP3$ndvi, vegID=datP3$vegID,doyN=datP3$doyN,
					pixID=datP3$pixID,NpixID=dim(pixDF)[1], vegIDP=pixDF$vegID,
				Nveg=dim(vegDF)[1])

parms <- c("ndvi.rep", "sig.ndvi","base","Amp","slopeG","slopeB","G1","B1", "mu.Amp","mu.slopeG","mu.slopeB","mu.G1","mu.B1","half.G","half.B",
			"sig.base","sig.Amp","sig.slopeG","sig.slopeB","sig.G1","sig.B1")

samp.modI<-jags.model(file="c:\\Users\\hkropp\\Documents\\GitHub\\ness_phenology\\phenology_model_code.r",
						data=datalist,
						n.adapt=50000,
						n.chains=3)

samp.sample <- coda.samples(samp.modI,variable.names=parms,
                       n.iter=200000, thin=200)	
			   
mcmcplot(samp.sample, parms=c("sig.ndvi","base","Amp","slopeG","slopeB","G1","B1",
							"mu.Amp","mu.slopeG","mu.slopeB","mu.G1","mu.B1","half.G","half.B",
							"sig.base","sig.Amp","sig.slopeG","sig.slopeB","sig.G1","sig.B1"),
			dir=paste0(modDI,"\\history"))	
					   
mod.out <- summary(samp.sample)

write.table(mod.out$statistics,paste0(modDI,"\\comp_mod_stats.csv"),
			sep=",",row.names=TRUE)
write.table(mod.out$quantiles,paste0(modDI,"\\comp_mod_quant.csv"),
			sep=",",row.names=TRUE)
}

###############################################
##############look at results##################
###############################################			
#read in model results
modS <- read.csv(paste0(modDI,"\\comp_mod_stats.csv"))		
modQ <- read.csv(paste0(modDI,"\\comp_mod_quant.csv"))				

datC <- cbind(modS,modQ)
			
dexps<-"\\[*[[:digit:]]*\\]"

datC$parms <- gsub(dexps,"", rownames(datC))	

#check model fit
mfit <- lm(datC$Mean[datC$parms=="ndvi.rep"]~datP2$ndvi)
jpeg(paste0(modDI,"\\modelfit.jpeg"), width=700, height=700, units="px",quality=100)
	plot(datP2$ndvi, datC$Mean[datC$parms=="ndvi.rep"], xlim=c(0,1), ylim=c(0,1), xlab="observed ndvi", ylab="predicted ndvi",pch=19)
	text(.3,.9, paste("pred=",round(summary(mfit)$coefficients[1,1],3), "+",round(summary(mfit)$coefficients[2,1],2), "*obs"))
	text(.3,.8, paste("r2=",round(summary(mfit)$r.squared,3)))
	abline(0,1,lwd=2,col="red")
	abline(mfit, lwd=2,lty=3)
dev.off()


#check if vegetation is different
jpeg(paste0(modDI,"\\halfmax.jpeg"), width=700, height=700, units="px",quality=100)
	par(mfrow=c(1,2))
	plot(seq(1,3), datC$Mean[datC$parms=="half.G"], xaxt="n", xlab="vegetation type", ylab="halfway max green", pch=19, ylim=c(0,.01))
	axis(1,seq(1,3), vegDF$veg)
	arrows(seq(1,3),datC$X2.5[datC$parms=="half.G"],seq(1,3),datC$X97.5[datC$parms=="half.G"], code=0)
	plot(seq(1,3), datC$Mean[datC$parms=="half.B"], xaxt="n", xlab="vegetation type", ylab="halfway max brown", pch=19, ylim=c(0.1,.15))
	axis(1,seq(1,3), vegDF$veg)
	arrows(seq(1,3),datC$X2.5[datC$parms=="half.B"],seq(1,3),datC$X97.5[datC$parms=="half.B"], code=0)
dev.off()

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

yval3 <- phenC3(.5,.3,0.4,10,.6,10,xval)

plot(xval,yval1, type="l", lwd=2, col="cornflowerblue")
points(xval,yval3, type="l",lwd=2, col="tomato3")
abline(v=.4)
abline(v=.6)

###################################################