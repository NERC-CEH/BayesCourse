# JAGS model

# y is the observed growth rate
# x is the measurement of light, L

# a is alpha, the max. growth rate at infinite light
# c is light level at which growth is zero
# b is rate at which curve tails off (gamma)

########################################################################################
model{
  ### likelihood
  # Data model
  for (i in 1:n)
  {
    y[i] ~ dnorm(mu[i], tau.o) #Note JAGS uses precision, not variance
  }  
  
  # process model
  for (i in 1:n)
  {
    mu[i] ~ dnorm(mu2[i], tau.p) #Note JAGS uses precision, not variance
    
    mu2[i] <- a * (x[i]-c) / ((a/b)+(x[i]-c))
    
  }
  
  # priors
  a ~ dgamma(0.01,0.01)
  c ~ dunif(-10,10)
  b ~ dgamma(0.01,0.01)
  
  sigma.o ~ dnorm(5, 1/(0.5*0.5)) ## assume prior knowledge of observation error (5 with SD of 0.5)
  sigma.p ~ dunif(0, 50)
  
  # derived quantities: convert precisions to standard deviations
  tau.o <- 1/(sigma.o * sigma.o)
  tau.p <- 1/(sigma.p * sigma.p)
  
} # end of model
#########################################################################################

