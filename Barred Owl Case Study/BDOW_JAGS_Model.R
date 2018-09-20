model {
  #Priors
  lambda2 ~ dunif(-5,5) # Log transformation on prior appears to help initialize the model
  log(lambda) <- lambda2 
  a0 ~ dnorm(0,0.5)  # intercept on gains
  a1 ~ dnorm(0,0.5) # slope on gains for neighborhood effect
  b0 ~ dnorm(0,0.5) # effect of crepuscular on detection
  b1 ~ dnorm(0,0.5) # effect of day on detection
  b2 ~ dnorm(0,0.5) # effect of night on detection
  c0 ~ dnorm(0,0.5) # intercept on survival
  c1 ~ dnorm(0,0.5) # slope on survival for habitat covariate
  
  for (i in 1:nSites){
    N[i,1] ~ dpois(lambda) #intialize counts at each site
    for (t in 2:nYears){
      logit(omega[i,t-1]) <- c0 + c1*habcov[i,t-1]  #estimate survial
      S[i,t-1] ~ dbin(omega[i,t-1], N[i,t-1])       #estimate number that survive
      G[i,t-1] ~ dpois(gamma[t-1])                  #estimate gains
      N[i,t] <- S[i,t-1] + G[i,t-1]               #calculate counts
    }
  }
#Detection model for occupancy data
for (k in 1:nSamples){
  occ.p[k] <- 1-pow((1-p[k]),N[site[k],year[k]]) # probability of not detecting any individuals
  y[k] ~ dbern(occ.p[k])                             
  logit(p[k])<-b0+b1*day[k]+b2*night[k] # estimate detection of a single indivudal
}

# covariate for gamma
for (t in 2:nYears){
   log(gamma[t-1]) <- a0 + a1*N.mean[t-1] 
   N.mean[t-1] <- mean(N[,t-1])          # per site mean from previous year
}
}