#Simulated consumer code

#Import prey data
prey <- read.csv("Prey.csv",header=TRUE)

#Number of consumers
N = 5

#Time-steps
t_term <- 100

#Matrix for saving consumer C values
c_m <- matrix(0,N,t_term)

#Matrix for saving consumer N values
n_m <- matrix(0,N,t_term)

#Matrix cor saving consumer prey bouts
d_m <- matrix(0,N,t_term)

#Initial consumer isotope values
