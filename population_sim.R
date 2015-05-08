rm(list=c(ls()))

#Simulated consumer code
library(ellipse)
library(gtools)
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
N = 3

#Body size of consumers (kg)
bmass <- rep(20,N)

#Time-steps
t_term <- 5000

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
#theta <- rep(1,N)

#Dirichlet distribution for p_vector
#p_vector is the proportional contribution of prey to the diet of consumer
#p_vector determines a per-diem vector of proportional contributions
#given the Dirichlet from which it is drawn
Dir_param <- matrix(1,N,nprey)

Dir_param[1,7] <- 10000
Dir_param[2,6] <- 10000
Dir_param[3,6] <- 10


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
    #mb <- bmass[i]
    
    #Determine next prey item
    
    #With probability equal to e_gen[i], they will specialize
    #With probability equal to 1-e_gen[i], they will draw from prey randomly
    
    #RANDOM DRAW VERSION
    #Draw random value
    #     rdraw <- runif(1)
    #     if (rdraw < theta[i]) {
    #       #If specialist, select it's preferred prey
    #       next_prey <- s_prey[i]
    #     } else {
    #       #If generalist, randomly select from all prey
    #       next_prey <- sample(nprey,1,replace=TRUE)
    #     }
    #     #Randomly draw prey isotope values from known mean and sd
    #     cp_mean <- prey$CM[next_prey]
    #     cp_sd <- prey$CSD[next_prey]
    #     cp <- rnorm(1,cp_mean,cp_sd)
    #     np_mean <- prey$NM[next_prey]
    #     np_sd <- prey$NSD[next_prey]
    #     np <- rnorm(1,np_mean,np_sd)
    
    
    #DIRICHLET VERSION
    #draw proportional contribution vector from random dirichlet
    p_vec <- numeric(nprey)
    #if (N == 1) {
    #  p_vec[1:nprey-1] <- rdirichlet(1,Dir_param)
    #} else {
    p_vec[1:nprey] <- rdirichlet(1,Dir_param[i,])
    #}
    #p_vec[nprey] <- 1 - sum(p_vec)
    #Draw prey values from each prey
    cp_vec <- sapply(seq(1,nprey),function(x){rnorm(1,prey$CM[x],prey$CSD[x])})
    np_vec <- sapply(seq(1,nprey),function(x){rnorm(1,prey$NM[x],prey$NSD[x])})
    #Calculate weighted average
    cp <- p_vec %*% cp_vec
    np <- p_vec %*% np_vec
    
    #Prey biomass
    #set to one if each prey is to be equally weighted
    #(assume 1 kg of each thing is eaten rather than at individual level)
    #mp <- 1 #prey$Biomass[next_prey]
    
    #Define incorporation rate
    incorp_rate <- 0.05
    
    #weights for body size
    #f <- mb/(mb + mp)
    f <- 1 - incorp_rate
    
    cb_next <- f*cb + (1-f)*cp
    nb_next <- f*nb + (1-f)*np
    
    c_m[i,t+1] <- cb_next
    n_m[i,t+1] <- nb_next
    
  } #end i
  
} #end t



#Calculate the analytical expectation and variance trajectory for EACH INDIVIDUAL
E_lim_c <- numeric(N)
E_lim_n <- numeric(N)
Var_lim_c <- numeric(N)
Var_lim_n <- numeric(N)

for (i in 1:N) {
  
  #ind <- 1
  #plot(c_m[ind,],pch=16,cex=0.25) #; lines(c_m[ind,])
  
  #Dirichlet-Normal Approximation for expectation and variance
  a0 <- sum(Dir_param[i,1:nprey])
  
  #Expectation of the sum(p_i*X_i) quantity
  EDir_c <- (Dir_param[i,1:nprey]/a0) %*% prey$CM
  EDir_n <- (Dir_param[i,1:nprey]/a0) %*% prey$NM
  
  #Expectation of p_i
  Ep <- Dir_param[i,1:nprey]/a0
  
  #Variance of p_i
  VarDir_p <- sapply(Dir_param[i,],function(x){(x*(a0-x)) / (a0^2*(a0+1))})
  
  #Covariance of p_i,p_j
  CovDir_p <- matrix(0,nprey,nprey)
  for (k in 1:nprey) {
    for (j in 1:nprey) {
      CovDir_p[k,j] <- -(Dir_param[i,k]*Dir_param[i,j]) / (a0^2*(a0 + 1))
    }
  } 
  diag(CovDir_p) <- 0
  
  #Variance of the sum(p_i*X_i) quantity for CARBON
  VarDir_c_vec <- numeric(nprey)
  for (j in 1:nprey) {
    VarDir_c_vec[j] <- (prey$CSD[j]^2)*VarDir_p[j] + Ep[j]^2*(prey$CSD[j])^2  + VarDir_p[j]*prey$CM[j]^2 + prey$CM[j]*sum(CovDir_p[j,-j] * prey$CM[-j])       
  }
  VarDir_c <- sum(VarDir_c_vec)
  
  #Variance of the sum(p_i*X_i) quantity for NITROGEN
  VarDir_n_vec <- numeric(nprey)
  for (j in 1:nprey) {
    VarDir_n_vec[j] <- (prey$NSD[j]^2)*VarDir_p[j] + Ep[j]^2 * prey$NSD[j]^2  + VarDir_p[j]*prey$NM[j]^2 + prey$NM[j]*sum(CovDir_p[j,-j] * prey$NM[-j])       
  }
  VarDir_n <- sum(VarDir_n_vec)
  
  #Individual expectations as a function of time.
  analyticEDir_c <- sapply(seq(1,t_term),function(x){f^x*(c_init - EDir_c) + EDir_c})
  analyticEDir_n <- sapply(seq(1,t_term),function(x){f^x*(n_init - EDir_n) + EDir_n})
  #Individual SD as a function of time
  analyticDirSD_c <- sapply(seq(1,t_term),function(x){sqrt(0.5*VarDir_c*(f-1)*(exp(2*(f-1)*x)-1))})
  analyticDirSD_n <- sapply(seq(1,t_term),function(x){sqrt(0.5*VarDir_n*(f-1)*(exp(2*(f-1)*x)-1))})
  
  #Limit
  E_lim_c[i] <- analyticEDir_c[t_term]
  E_lim_n[i] <- analyticEDir_n[t_term]
  Var_lim_c[i] <- analyticDirSD_c[t_term]^2
  Var_lim_n[i] <- analyticDirSD_n[t_term]^2
  
}


#Combine individual measurements to get population statistics (at the limit)
Exp_lim_pop <- c(rep((1/N),N) %*% E_lim_c,rep((1/N),N) %*% E_lim_n)

Var_lim_pop_vec_c <- numeric(N)
Var_lim_pop_vec_n <- numeric(N)
for (i in 1:N) {
  Var_lim_pop_vec_c[i] <- N*Var_lim_c[i] + (N-1)*E_lim_c[i]^2 - sum(E_lim_c[i]*E_lim_c[-i])
  Var_lim_pop_vec_n[i] <- N*Var_lim_n[i] + (N-1)*E_lim_n[i]^2 - sum(E_lim_n[i]*E_lim_n[-i])
}
Var_lim_pop_c <- (1/(N^2))*sum(Var_lim_pop_vec_c)
Var_lim_pop_n <- (1/(N^2))*sum(Var_lim_pop_vec_n)

#NEED COVARIANCE BTW C AND N FOR THE POPULATIONS

cov_pop <- ((rep((1/N),N) %*% (E_lim_c*E_lim_n)) - (Exp_lim_pop[1]*Exp_lim_pop[2]))
pop_Sigma <- matrix(c(Var_lim_pop_c,cov_pop,cov_pop,Var_lim_pop_n),2,2)

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
  lines(ellipse(x=0,scale=c(sqrt(Var_lim_c[i]),sqrt(Var_lim_n[i])),centre=c(E_lim_c[i],E_lim_n[i]),level=CI))
}
points(E_lim_c,E_lim_n,col="black",pch=16,cex=0.7)
#Population Ellipse
#pop_ellipse <- ellipse(x=0,scale=c(sqrt(Var_lim_pop_c),sqrt(Var_lim_pop_n)),centre=Exp_lim_pop)
pop_ellipse <- ellipse(pop_Sigma,centre=Exp_lim_pop)
points(Exp_lim_pop[1],Exp_lim_pop[2],pch=16,col="black")
lines(pop_ellipse,col="black",lty=3,lwd=3)



#CARBON
#Analytical approximation for pure specialist, SD=0

par(mfrow=c(2,1))
# analyticE <- sapply(seq(1,t_term),function(x){f^x*(c_init - cp_mean) + cp_mean})
# plot(c_m[ind,2:10000],pch=16,cex=0.5,xlab="time",ylab="d13C",col="gray")
# lines(analyticE)

#Plotting observed and expected values for the expectation CARBON
ind <- 2
analyticEDir <- sapply(seq(1,t_term),function(x){f^x*(c_init[init] - EDir_c) + EDir_c})
plot(c_m[ind,2:t_term],pch=16,cex=0.5,xlab="time",ylab="d13C",col="darkgray")
lines(analyticEDir[ind,])

#Plotting observed and expected values for the variance CARBON
binsize = 100
#analyticSD <- sapply(seq(1,t_term),function(x){sqrt(0.5*cp_sd^2*(f-1)*(exp(2*(f-1)*x)-1))})
analyticDirSD_c <- sapply(seq(1,t_term),function(x){sqrt(0.5*VarDir_c*(f-1)*(exp(2*(f-1)*x)-1))})
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
     pch=16,cex=0.5,xlim=c(0,tail(bins,n=1)),ylim=c(0,0.5),xlab="time",ylab="d13C SD",col="black")
lines(analyticDirSD_c)
#lines(analyticSD)


#NITROGEN
#Analytical approximation for pure specialist, SD=0
# par(mfrow=c(2,1))
# analyticE <- sapply(seq(1,t_term),function(x){f^x*(n_init - np_mean) + np_mean})
# plot(n_m[ind,2:t_term],pch=16,cex=0.5,xlab="time",ylab="d15N",col="gray")
# lines(analyticE)  

#Plotting observed and expected values for the variance NITROGEN
analyticEDir <- sapply(seq(1,t_term),function(x){f^x*(n_init - EDir_n) + EDir_n})
plot(n_m[ind,2:t_term],pch=16,cex=0.5,xlab="time",ylab="d15N",col="gray")
lines(analyticEDir)

#Plotting observed and expected values for the variance NITROGEN
binsize = 5000
#analyticSD <- sapply(seq(1,t_term),function(x){sqrt(0.5*np_sd^2*(f-1)*(exp(2*(f-1)*x)-1))})
analyticDirSD_n <- sapply(seq(1,t_term),function(x){sqrt(0.5*VarDir_n*(f-1)*(exp(2*(f-1)*x)-1))})
bins <- seq(binsize+1,t_term,by=binsize)
sd_bin <- numeric(length(bins)-1)
for (i in 1:(length(bins)-1)) {
  bin0 <- bins[i]
  bin1 <- bins[i+1]
  sd_bin[i] <- sd(n_m[ind,bin0:bin1])
}
#sapply(seq(binsize+1,t_term,by=binsize),function(x){sd(c_m[ind,x-binsize:x])})
#I don't know why ^ doesn't work!
plot(bins[-1],sd_bin,
     pch=16,cex=0.5,xlim=c(0,tail(bins,n=1)),ylim=c(0,0.5),xlab="time",ylab="d15N SD",col="black")
lines(analyticDirSD_n)


#Exploring the population-level parameter range
#specialization measure

epsilon_v <- seq(0,1,length.out=10)
avec <- matrix(0,10,nprey)
generalist <- rep(1/nprey,nprey)
specialist <- numeric(nprey); specialist[1] <- 1

#Vary from generalist to a specialist
for (i in 1:10) {
  avec[i,] <- epsilon_v[i]*generalist + (1-epsilon_v[i])*specialist
}


sim_matrix <- matrix(1,N,N)
diff_matrix <- matrix(0,N,N); diag(diff_matrix) <- 1
#similarity measure
psi_v <- seq(0,1,length.out=10)

sim_list <- list()
for (i in 1:length(psi)) {
  sim_list[[i]] <- psi_v[i]*sim_matrix + (1-psi_v[i])*diff_matrix
}

a_matrix <- matrix(0,N,nprey)

#For first individual, given value of epsilon, randomly sample alpha vector
#for second individual, given value of epsilon, and similarity between ind2 and ind1, randomly sample alpha vector
#for third individual, given value of epsilon, and similarity between ind3 and [ind1, ind2], randomly sample alpha vector
#Do this through N individuals
#Calculate Pop-Level Variance for each [epsilon, phi]

#Do this L times, and determine E[Pop-Level Variance] over L

#Plot

for (i in 1:length(epsilon)) {
  for (j in 1:length(psi)) {
    eps <- epsilon_v[i]
    psi <- psi_v[j]
    
    
    
    
    
  }
}

#Or just sample the shit out of the parameter space
L <- 10000
Epsi <- numeric(L)
Eepsilon <- numeric(L)
gen_vec <- rep(1/nprey,nprey)
spec_vec <- numeric(nprey); spec_vec[1] <- 1

for (l in 1:L) {
  #Build population p matrix
  #p_m <- matrix(runif(N*nprey),N,nprey)
  p_m <- matrix(rgamma(N*nprey,runif(1,1,10),runif(1,1,10)),N,nprey)
  p_m <- p_m/apply(p_m,1,sum) #ensures rows sum to 1
  
  #Build similarity matrix
  S_m <- matrix(0,N,N)
  for (i in 1:N) {
    for (j in 1:N) {
      v1 <- p_m[i,]
      v2 <- p_m[j,]
      mag1 <- sqrt(v1 %*% v1)
      mag2 <- sqrt(v2 %*% v2)
      S_m[i,j] <- (v1/mag1) %*% (v2/mag2)
    }
  }
  
  #Calculate <epsilon> and psi
  psi <- diag(ginv(sim_matrix - diff_matrix) %*% (S_m - diff_matrix))
  Epsi[l] <- 1 - mean(psi)
  
  M <- as.matrix(cbind(gen_vec,spec_vec))
  Minv <- ginv(M)
  Eepsilon[l] <- min(mean(apply(p_m,1,function(x){Minv %*% x})[1,]),1)
}

plot(Epsi,Eepsilon,xlim=c(0,1),ylim=c(0,1))







