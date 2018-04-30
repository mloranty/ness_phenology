##############################################
##model code for fitting curves and analysis##
##############################################

model{
	#likelihood
	for(i in 1:Nobs){
		ndvi[i] ~ dnorm(mu.ndvi[i],tau.ndvi[vegID[i]])
		mu.ndvi[i] <- base[pixID[i]]+ (Amp[pixID[i]]/(Greening[i]*Browning[i]))
		Greening[i] <- 1+exp(slopeG[pixID[i]]*(G1[pixID[i]]-doyN[i]))
		Browning[i] <- 1+exp(slopeB[pixID[i]]*(doyN[i]-B1[pixID[i]]))
		ndvi.rep[i] ~ dnorm(mu.ndvi[i],tau.ndvi[vegID[i]]) 
	}
	#priors
	for(i in 1:NpixID){
		base[i] ~ dnorm(mu.base[vegIDP[i]],tau.base)T(0,1)
		Amp[i] ~ dnorm(mu.Amp[vegIDP[i]],tau.Amp)T(0,1)
		slopeG[i] ~ dnorm(mu.slopeG[vegIDP[i]],tau.slopeG)T(0,)
		slopeB[i] ~ dnorm(mu.slopeB[vegIDP[i]],tau.slopeB)T(0,)
		G1[i] ~ dnorm(mu.G1[vegIDP[i]],tau.G1)T(0,.5)
		B1[i] ~ dnorm(mu.B1[vegIDP[i]],tau.B1)T(.5,1)

	}
	#hyper priors
	for(i in 1:Nveg){
		mu.base[i] ~ dunif(0,.6)
		mu.Amp[i] ~ dunif(0,1)
		mu.slopeG[i] ~ dunif(0,300)
		mu.slopeB[i] ~ dunif(0,300)
		mu.G1[i] ~ dunif(0,.5)
		mu.B1[i] ~ dunif(0.5,1)
		tau.ndvi[i] <- pow(sig.ndvi[i],-2)
		sig.ndvi[i] ~ dunif(0,1)
		half.G[i] <- mu.G1[i]/mu.slopeG[i]
		half.B[i] <- mu.B1[i]/mu.slopeB[i]
	
		
	}
	tau.base <- pow(sig.base,-2)
	sig.base ~ dunif(0,1)
	
	tau.Amp <- pow(sig.Amp,-2)
	sig.Amp ~ dunif(0,1)	
	
	tau.slopeG <- pow(sig.slopeG,-2)
	sig.slopeG ~ dunif(0,1)
	
	tau.slopeB <- pow(sig.slopeB,-2)
	sig.slopeB ~ dunif(0,1)
	
	tau.G1 <- pow(sig.G1,-2)
	sig.G1 ~ dunif(0,1)
	
	tau.B1 <- pow(sig.B1,-2)
	sig.B1 ~ dunif(0,1)
	
}

