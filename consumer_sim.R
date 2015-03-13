#Simulated consumer code
library(ellipse)


#Import prey data
prey <- read.csv("Prey.csv",header=TRUE)
nprey <- length(prey[,1])
ellip_prey <- list(nprey); CI <- 0.95
#build ellipses for prey
for (i in 1:nprey) {
  Sigma <- matrix(c(prey$CSD[i],0,0,prey$NSD[i]),2,2)
  mu <- c(prey$CM[i],prey$NM[i])
  ellip_prey[[i]] <- ellipse(Sigma,centre=mu,level=CI)
}

#Number of consumers
N = 10

#Body size of consumers (kg)
bmass <- rep(20,N)

#Time-steps
t_term <- 1000

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
e_gen <- rep(1,N)
#Which prey item does each consumer specialize on?
s_prey <- sample(nprey,N,replace=TRUE)

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
    if (e_gen[i] == 0) {
      #If generalist, randomly select from all prey
      next_prey <- sample(nprey,1,replace=TRUE)
    } else {
      #If specialist, select it's preferred prey
      next_prey <- s_prey[i]
    }
    
    #Randomly draw prey isotope values from known mean and sd
    cp_mean <- prey$CM[next_prey]
    cp_sd <- prey$CSD[next_prey]
    cp <- rnorm(1,cp_mean,cp_sd)
    
    np_mean <- prey$NM[next_prey]
    np_sd <- prey$NSD[next_prey]
    np <- rnorm(1,np_mean,np_sd)
    
    #weights for body size
    f <- mb/(mb + 1)
    
    cb_next <- f*cb + (1-f)*cp
    nb_next <- f*nb + (1-f)*np
    
    c_m[i,t+1] <- cb_next
    n_m[i,t+1] <- nb_next
    
  } #end i
  
} #end t

#Plot prey ellipses
plot(ellip_prey[[1]],type="l",xlim=c(-22,-10),ylim=c(5,20),col="gray")
for (i in 2:nprey) {
  lines(ellip_prey[[i]],col="gray")
}
for (i in 1:N) {
ind <- i
points(c_m[ind,],n_m[ind,],pch=16,cex=0.25)
}
