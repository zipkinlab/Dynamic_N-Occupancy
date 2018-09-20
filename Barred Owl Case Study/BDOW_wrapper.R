rm(list=ls())

# Set working directory where files are located
#setwd(dirname(file.choose())) # select BDOW.csv or Habitat.csv

# Load data
habcov<-read.csv("Habitat.csv") # Habitat covariate
nYears<-ncol(habcov) # Years in data
nSites<-nrow(habcov) # Sites in data
nReps<-8             # Max number of reps

BDOW.raw<-read.csv("BDOW.csv") # Barred oWl data
day.raw<-read.csv("day.csv")       # Variable indicating day sampling events
night.raw<-read.csv("night.csv")   # Variable indicating night sampling events
#Format data into 3D array
reshape.seq<-seq(0,ncol(BDOW.raw),8)
BDOW<-night<-day<-array(NA, dim=c(nSites,nReps,nYears))
for (t in 1:(length(reshape.seq)-1)){
  BDOW[,,t]<-as.matrix(BDOW.raw[,(reshape.seq[t]+1):reshape.seq[(t+1)]])
  day[,,t]<-as.matrix(day.raw[,(reshape.seq[t]+1):reshape.seq[(t+1)]])
  night[,,t]<-as.matrix(night.raw[,(reshape.seq[t]+1):reshape.seq[(t+1)]])
}

# Convert arrays to vectors with index variables
year<-site<-array(NA, dim=c(nSites, nReps, nYears))
for (t in 1:nYears){
  year[,,t] <- t
}
for (i in 1:nSites){
  site[i,,] <- i
}
BDOW<-c(BDOW)
site<-c(site)[!is.na(BDOW)]
year<-c(year)[!is.na(BDOW)]
day<-c(day)[!is.na(BDOW)]
night<-c(night)[!is.na(BDOW)]
BDOW<-c(BDOW)[!is.na(BDOW)]


#Parameters to save
params <- c("a0", "a1", "b0", "b1", "b2", "c0", "c1","S","G", "lambda","N")

# Compile JAGS data
jags.data<-list(y=BDOW,
                nYears=nYears,
                nSites=nSites,
                nSamples=length(BDOW),
                year=year,
                site=site,
                day=day,
                night=night,
                habcov=habcov)

# Generate initial values
inits<-function() list(N=matrix(c(rep(10,nSites),rep(NA,nSites*(nYears-1))),nrow=nSites,ncol=nYears),
                       S=matrix(5,nrow=nSites,ncol=(nYears-1)),
                       G=matrix(2,nrow=nSites,ncol=(nYears-1))
                       )
library(jagsUI)
  ni=10000
  na=5000
  nb=5000
  nc=3
  nt=20

posterior<-((ni-nb)*nc)/nt
model<-"BDOW_JAGS_Model.R"
jags.out<-jags(jags.data, inits, params, model.file=model, store.data = TRUE, n.adapt=na,
                   n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nt, parallel = TRUE)


