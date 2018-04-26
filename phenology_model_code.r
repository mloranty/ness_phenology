##############################################
##model code for fitting curves and analysis##
##############################################

model{
	for(i in 1:Nobs){
		ndvi[i]~dnorm(mu.ndvi[i],tau.ndvi)
		mu.ndvi[i] <- base[pixelID[i]]+ (Amp[pixelID[i]]/(Greening[i]*Browning[i]))
		Greening[i] <- 1+exp(slopeG[pixelID[i]]*(G1[pixelID[i]]-doy[i]))
		Browning[i] <- 1+exp(slopeB[pixelID[i]]*(doy[i]-B1[pixelID[i]])
	}

	for(i in 1:Npixels){
		base[i] ~ dunif(0,.6)
		Amp[i] ~ dunif(0,1)
		slopeG[i] ~ dunif(0,300)
		slopeB[i] ~ dunif(0,300)
		G1[i] ~ dunif(0,.5)
		B1[i] ~ dunif(0.5,1)
	}

}