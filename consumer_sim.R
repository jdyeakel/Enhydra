#Simulated consumer code
library(ellipse)
library(RColorBrewer)


#Import prey data
prey <- read.csv("Prey.csv",header=TRUE)
nprey <- length(prey[,1])
ellip_prey <- list(nprey); CI <- 0.95
#build ellipses for prey
for (i in 1:nprey) {
  #Sigma <- matrix(c(prey$CSD[i],0,0,prey$NSD[i]),2,2)
  mu <- c(prey$CM[i],prey$NM[i])
  ellip_prey[[i]] <- ellipse(x=0,scale = c(prey$CSD[i],prey$NSD[i]),centre=mu,level=CI)
}

#Number of consumers
N = 1

#Body size of consumers (kg)
bmass <- rep(20,N)

#Time-steps
t_term <- 10000

#Matrix for saving consumer C values
c_m <- matrix(0,N,t_term)

#Matrix for saving consumer N values
n_m <- matrix(0,N,t_term)

#Matrix cor saving consumer prey bouts
d_m <- matrix(0,N,t_term)

#Initial consumer isotope values
#For now, the initial value will just be the means of the prey
c_init <- rep(mean(prey$CM),N)
n_init <- rep(mean(prey$NM),N)

c_m[,1] <- c_init
n_m[,1] <- n_init

#0 indicates generalist; 1 indicates specialist
#theta <- rep(1,N)
theta <- rep(1,N)

#Which prey item does each consumer specialize on?
s_prey <- sample(nprey,N,replace=TRUE)
#cumulative distribution
#plot(sapply(seq(1,12),function(x){length(which(s_prey<x))/length(s_prey)}))

#Loop over time
for (t in 1:(t_term-1)) {
  
  #Loop over consumers
  for (i in 1:N) {
    #Define the carbon, nitrogen isotope values of body at time t
    cb <- c_m[i,t]
    nb <- n_m[i,t]
    
    #Define body mass of consumer
    mb <- bmass[i]
    
    #Determine next prey item
    
    #With probability equal to e_gen[i], they will specialize
    #With probability equal to 1-e_gen[i], they will draw from prey randomly
    
    #Draw random value
    rdraw <- runif(1)
    if (rdraw < theta[i]) {
      #If specialist, select it's preferred prey
      next_prey <- s_prey[i]
    } else {
      #If generalist, randomly select from all prey
      next_prey <- sample(nprey,1,replace=TRUE)
    }
    
    #Prey biomass
    #set to one if each prey is to be equally weighted
    #(assume 1 kg of each thing is eaten rather than at individual level)
    mp <- 1 #prey$Biomass[next_prey]
    
    #Randomly draw prey isotope values from known mean and sd
    cp_mean <- prey$CM[next_prey]
    cp_sd <- prey$CSD[next_prey]
    cp <- rnorm(1,cp_mean,cp_sd)
    
    np_mean <- prey$NM[next_prey]
    np_sd <- prey$NSD[next_prey]
    np <- rnorm(1,np_mean,np_sd)
    
    #weights for body size
    f <- mb/(mb + mp)
    
    cb_next <- f*cb + (1-f)*cp
    nb_next <- f*nb + (1-f)*np
    
    c_m[i,t+1] <- cb_next
    n_m[i,t+1] <- nb_next
    
  } #end i
  
} #end t


colors <- rep(brewer.pal(8,"Set1"),round(N/9)+1)

#Plot prey ellipses
plot(ellip_prey[[1]],type="l",xlim=c(-22,-10),ylim=c(6,18),col="gray",xlab="d13C",ylab="d15N")
for (i in 2:nprey) {
  lines(ellip_prey[[i]],col="gray")
}
for (i in 1:N) {
ind <- i
points(c_m[ind,],n_m[ind,],pch=16,cex=0.25,col=colors[i])
lines(c_m[ind,],n_m[ind,],pch=16,cex=0.25,col=colors[i])
}

ind <- 1
plot(c_m[ind,],pch=16,cex=0.25) #; lines(c_m[ind,])

par(mfrow=c(2,1))

#Analytical approximation for pure specialist, SD=0
analyticE <- sapply(seq(1,t_term),function(x){f^x*(c_init - cp_mean) + cp_mean})
plot(c_m[ind,2:10000],pch=16,cex=0.5,xlab="time",ylab="d13C",col="gray")
lines(analyticE)

binsize = 1000
analyticSD <- sapply(seq(1,t_term),function(x){sqrt(0.5*cp_sd^2*(f-1)*(exp(2*(f-1)*x)-1))})
bins <- seq(binsize+1,t_term,by=binsize)
sd_bin <- numeric(length(bins)-1)
for (i in 1:(length(bins)-1)) {
  bin0 <- bins[i]
  bin1 <- bins[i+1]
  sd_bin[i] <- sd(c_m[ind,bin0:bin1])
}
#sapply(seq(binsize+1,t_term,by=binsize),function(x){sd(c_m[ind,x-binsize:x])})
#I don't know why ^ doesn't work!
plot(bins[-1],sd_bin,
     pch=16,cex=0.5,xlim=c(0,tail(bins,n=1)),ylim=c(0,0.5),xlab="time",ylab="d13C",col="gray")
lines(analyticSD)



