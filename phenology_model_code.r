##############################################
##model code for fitting curves and analysis##
##############################################

model{
	#likelihood
	for(i in 1:Nobs){
		ndvi[i] ~ dnorm(mu.ndvi[i],tau.ndvi[vegID[i]])
		mu.ndvi[i] <- base[vegID[i]]+ (Amp[vegID[i]]/(Greening[i]*Browning[i]))
		Greening[i] <- 1+exp(slopeG[vegID[i]]*(G1[vegID[i]]-doyN[i]))
		Browning[i] <- 1+exp(slopeB[vegID[i]]*(doyN[i]-B1[vegID[i]]))
		ndvi.rep[i] ~ dnorm(mu.ndvi[i],tau.ndvi[vegID[i]]) 
	}
	#priors
	for(i in 1:Nveg){
		base[i] ~ dunif(0,.6)
		Amp[i] ~ dunif(0,1)
		slopeG[i] ~ dunif(0,300)
		slopeB[i] ~ dunif(0,300)
		G1[i] ~ dunif(0,.5)
		B1[i] ~ dunif(0.5,1)
		tau.ndvi[i] <- pow(sig.ndvi[i],-2)
		sig.ndvi[i] ~ dunif(0,1)
		half.G[i] <- G1[i]/slopeG[i]
		half.B[i] <- B1[i]/slopeB[i]
	}

}