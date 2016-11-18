############################################################
# Data generation
############################################################

# "True" values
lam <- runif(1,0.5,3)
omega <- runif(1,0,1)
gamma <- runif(1,0,3)
p <- runif(1,0,1)
nYears <- 10
nReps <- 3
nSites<- 75

#Simulate true abundances, N, for each location
N <- matrix(NA, nSites, nYears)
S <- G <- matrix(NA, nSites, nYears-1)

#First year of sampling follows a Poisson distribution
N[,1] <- rpois(nSites, lam)

#subsequent years follow the birth-death-immigration process
for(t in 2:nYears) {
  S[,t-1] <- rbinom(nSites, N[,t-1], omega)
  G[,t-1] <- rpois(nSites, gamma)
  N[,t] <- S[,t-1] + G[,t-1] 
}

# Generate data vector y for the counts
y <- array(NA, c(nSites, nYears, nReps))
for(t in 1:nYears) {
  for(j in 1:nReps) {
    y[,t,j] <- rbinom(nSites, N[,t], p)
  }
}

#Now assume that there are a vector of sites, x, have only have occupancy data
#And change their data to detection/nondetection data
for (i in 1:nSites) {
  for (j in 1:nReps) {
    a = which(y[i,,j]>0)
    y[i,a,j] = 1
  }
}

############################################################
# JAGS model
###########################################################

#Jags model
sink("dyn_N_Occ.txt")
cat("
    model {
    #Priors
    lambda2 ~ dunif(-10,10)       
    log(lambda) <- lambda2
    gamma2 ~  dunif(-10,10)
    log(gamma) <- gamma2
    omega ~ dunif(0, 1)
    p ~ dunif(0, 1)
    
    #Likelihood - Biological process model
    for(i in 1:nSites) {
    #First year of sampling - process and observation components
    N[i,1] ~ dpois(lambda)
    
    #All other years of sampling - process and observation components   
    for(t in 2:nYears) {
    S[i,t-1] ~ dbin(omega, N[i,t-1])
    G[i,t-1] ~ dpois(gamma)
    N[i,t] <- S[i,t-1] + G[i,t-1] 
    }
    }
    
    #Detection model for occupancy data
    for (i in 1:nSites) {
    for (t in 1:nYears) {
    occ.p[i,t] <- 1-pow( (1-p),N[i,t] )
    for (j in 1:nReps) {  
    y[i,t,j] ~ dbern(occ.p[i,t])
    }}}
    }
    ",fill = TRUE)
sink()


############################################################
# Code to run the JAGS model
############################################################

# Format data 
jags.data <- list(nSites=nSites, 
                  nYears=nYears, 
                  y=y,
                  nReps=nReps)

# Parameters monitored
params <- c("lambda", "gamma", "omega", "p", "N")

#Initial values (may need to be adjusted in certain scenarios)
Ni <- y[,,1]+20
Si <- S
Si[] <- 2
Gi<-matrix(10,nrow=nSites,ncol=(nYears-1))
Ni[,-1] <- NA
inits <- function() list(N=Ni,
                         S=Si,
                         G=Gi)

model=normalizePath("dyn_N_Occ.txt")

#Load the correct library
library(jagsUI)
# Compile the model
jags.out<-jags(jags.data, inits, params, model.file=model, store.data=TRUE,
               n.chains=3, n.iter=25000, n.burnin=5000, n.thin=50, parallel=TRUE)


#######################################################
# Estimate colonization and extinction dynamics
#######################################################

N.post<-jags.out$sims.list$N
nIter<-dim(N.post)[1]
col<-array(NA,dim=dim(N.post))
ext<-array(NA,dim=dim(N.post))

for (j in 1:nIter){
  for (t in 2:nYears){
    for (i in 1:nSites){
      #colonization
      if(N.post[j,i,t-1]==0){
        if(N.post[j,i,t]>0){
          col[j,i,t]<-1       # colonization event
        }else{
          col[j,i,t]<-0       # did not colonize when possible (site was not previously occupied)
        }}
      #extinction
      if(N.post[j,i,t-1]>0){
        if(N.post[j,i,t]==0){
          ext[j,i,t]<-1       # extinction event
        }else{
          ext[j,i,t]<-0       # no extinction event when possible (site was previously occupied)
        }}
    }}}

#calculate probability of colonization for each year then take credible interval across all iterations
col.summary<-apply(apply(col, c(1,3), mean, na.rm=TRUE), 2,quantile,na.rm=TRUE,probs=c(0.025,0.25,0.5,0.75,0.975))

#calculate P of extinction for each year then take credible interval across #all iterations
ext.summary<-apply(apply(ext, c(1,3), mean, na.rm=TRUE), 2,quantile,na.rm=TRUE,probs=c(0.025,0.25,0.5,0.75,0.975)) 
