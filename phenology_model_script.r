##############################################
##script for fitting curves and analysis    ##
##############################################

###############read in data###################
#cubesat data
datP <- read.csv("c:\\Users\\hkropp\\Documents\\GitHub\\ness_phenology\\planet_veg_class_sample.csv")


##############look at curves###################
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


###now look at data plot to make sure this curve makes sense
plot(datP$doy,datP$ndvi)

#steps

#remove negative values and potentially cloudy day
#normalize days

