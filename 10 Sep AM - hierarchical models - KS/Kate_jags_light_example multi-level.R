########################################################################################
model{
  ### likelihood
  # Data model
    for(s in 1:3){
      for(i in 1:n){
        y[i,s] ~ dnorm(mu[i,s], tau.o) #Note JAGS uses precision, not variance
      }  
    }
  # process model
    for(s in 1:3){
      for(i in 1:n){
        mu[i,s] ~ dnorm(mu2[i,s], tau.p) #Note JAGS uses precision, not variance
    
        mu2[i,s] <- a[s] * (x[i,s]-c) / ((a[s]/b)+(x[i,s]-c))
        }
    }
  
  # priors
  for(s in 1:3){
    a[s] ~ dnorm(mu.a, tau.a)
  }
  mu.a ~ dunif(0,100)
  c ~ dunif(-10,10)
  b ~ dgamma(0.01,0.01)
  
  sigma.o ~ dnorm(5, 1/(0.5*0.5)) ## assume prior knowledge of observation error (5 with SD of 0.5)
  sigma.p ~ dunif(0, 50)
  sigma.a ~ dunif(0, 50)
  
  # derived quantities: convert precisions to standard deviations
  tau.o <- 1/(sigma.o * sigma.o)
  tau.p <- 1/(sigma.p * sigma.p)
  tau.a <- 1/(sigma.a * sigma.a)
  
} # end of model
#########################################################################################

