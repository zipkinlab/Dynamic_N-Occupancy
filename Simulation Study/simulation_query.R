# This script allows users to query our simulation results for the "Dynamic N-Occupancy"
# model containing more than 30,000 simulations across a wide parameter space. 

# load simulation study
simulations<-read.table('https://github.com/zipkinlab/Dynamic_N-Occupancy/raw/master/Simulation%20Study/dynNocc_simulation.txt',
           fill=TRUE, row.names=NULL)

# filter simulations (values below are the maximum range)
lambda.range<-c(0.5, 3) # simualted range = 0.5 to 3
gamma.range<-c(0, 3)    # simualted range = 0 to 3
omega.range<-c(0, 1)    # simualted range = 0 to 1
p.range<-c(0, 1)        # simualted range = 0 to 1
occ.sites<-c(25, 50,  75, 100, 150)  # levels = c(25, 50, 75, 100, 150)

# remove unwanted data
sims <- simulations[which(simulations$true_lam>lambda.range[1]&simulations$true_lam<lambda.range[2]&
                            simulations$true_gamma>gamma.range[1]&simulations$true_gamma<gamma.range[2]&
                            simulations$true_omega>omega.range[1]&simulations$true_omega<omega.range[2]&
                            simulations$true_p>p.range[1]&simulations$true_p<p.range[2]&
                            simulations$nOcc%in%occ.sites),]

# Produce plot with 95% high density range
library(hdrcde)
param1 <- "true_omega"       #parameter for x and y axis (see 'colnames(simulations)' for options)
param2 <- "m.omega"          #parameter for y axis

sim.plot<-hdr.2d(sims[,param1], sims[,param2], prob=0.95)
plot.hdr2d(sim.plot, show.points = TRUE, xlab=paste(param1), ylab=paste(param2))
cor(sims[,param1],sims[,param2])
