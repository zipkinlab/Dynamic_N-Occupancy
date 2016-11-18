model {
  #Priors
  lambda2 ~ dunif(-5,5) # Log transformation on prior (helps initialize the model)
  a0 ~ dnorm(0,0.5)  # intercept on gains
  a1 ~ dnorm(0,0.5) # slope on gains for “neighborhood effect”
  b0 ~ dnorm(0,0.5) # effect of crepuscular time period on detection
  b1 ~ dnorm(0,0.5) # effect of day on detection
  b2 ~ dnorm(0,0.5) # effect of night on detection
  c0 ~ dnorm(0,0.5) # intercept on survival
  c1 ~ dnorm(0,0.5) # slope on survival for habitat covariate
  
  for (i in 1:nSites){				# number of sites
    log(lambda[i]) <- lambda2 			# transform log prior
    N[i,1] ~ dpois(lambda[i]) 			# initialize counts at each site
    for (t in 2:nYears){
      logit(omega[i,t]) <- c0 + c1*habcov[i,t-1]  # covariate on survival 
      S[i,t-1] ~ dbin(omega[i,t], N[i,t-1])       # estimate number of individuals                        
      # that survive
      G[i,t-1] ~ dpois(gamma[t])                  # estimate number of individuals 
      # gained
      N[i,t] <- S[i,t-1] + G[i,t-1]               # calculate local abundance
    }
  }
  
  #Detection model 
  for (t in 1:nYears) { 				# number of years
    for (i in 1:nSites){ 			      # number of sites 
      for (j in 1:rep.count[i,t]){ 	# number of survey reps at site i in year t
        occ.p[i,j,t] <- 1-pow((1- theta[i,j,t]),N[i,t])
        # probability of detecting at least one individual
        logit(theta[i,j,t])<-    
          b0 + b1*day[i,j,t] + b2*night[i,j,t] 
        # estimate covariate effects on detection 
        y[site[i,t],j,t] ~ dbern(occ.p[i,j,t])                              
        # detection model
      }}}
  
  # covariate for gamma
  for (t in 2:nYears){
    log(gamma[t]) <- a0 + a1*N.mean[t-1]  # covariates on gains 
    N.mean[t-1] <- mean(N[,t-1]) – 0.75   # per site mean from previous year, 0.75 was
    # subtracted to center the covariate
    
    