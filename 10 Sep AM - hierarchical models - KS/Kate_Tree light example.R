
setwd("N://STATS/CEH Bayes Course/")

library(rjags)
library(jagsUI)

rm(list=ls(all=TRUE)) # clear the R slate!

#############################################################################
### load in data
tree.data <- read.csv("Hemlock-light-data-simple.csv")
str(tree.data)
head(tree.data)
#############################################################################


#############################################################################
### do some plotting to look at data

plot(tree.data$Light, tree.data$Observed.growth.rate)

#############################################################################


#############################################################################
# Plot priors to visualise, using 'hist' or 'curve'
# EG
par(mfrow=c(1,2))
curve(dgamma(x, shape=0.001, rate=0.001), xlab="p",
     ylab="Probability Density",xlim=c(0,0.5), lwd=2, col="green")
curve(dgamma(x, shape=0.01, rate=0.01), xlab="p",
     ylab="Probability Density",xlim=c(0,0.5), lwd=2, col="red")
#############################################################################


#############################################################################

## set up for running model in JAGS

# Initialize the MCMC chains. This needs to be a list.
inits <- function() list(a=runif(1,35,40), b=runif(1,1,3), c=runif(1,-5,5), 
                         sigma.o=1, sigma.p=1) #small tau makes variance big
inits

# Specify data; must be a list. 
data=list(
	n=nrow(tree.data),
	x=as.numeric(tree.data$Light),
	y=as.numeric(tree.data$Observed.growth.rate )
	)
data

# Parameters monitored (i.e., for which estimates are saved)
params <- c("a","b","c","sigma.o","sigma.p")

# MCMC settings
ni <- 75000   ;   nt <- 1   ;   nb <- 50000   ;  nc <- 3

### fit the model in JAGS
model_1 <- jags(data=data, inits, params, 
             model.file = "Kate_jags_light_example.R", n.chains = nc, 
             n.thin = nt, n.iter = ni, n.burnin = nb)

summary(model_1)
print(model_1, 6) # the 6 sets the number of significant digits to print out
plot(model_1) # check convergence


#############################################################################
## extract the MCMC estimates for each model parameter
names(model_1)
str(model_1$sims.list)
mymodel.dat <- as.data.frame(model_1$sims.list)
head(mymodel.dat)
## plot observed data against predictions
par(mfrow=c(1,1))
plot(tree.data$Light, tree.data$Observed.growth.rate)
#############################################################################


#############################################################################
## plot observed against predicted using all MCMC output

## Simulate the range of the environmental variable (Light)
x.sim <- seq(min(tree.data$Light), max(tree.data$Light), by = 0.1)
#############################################################################
## Calculate prediction over all mcmc samples for each parameter
int.sim <- matrix(rep(NA, nrow(mymodel.dat)*length(x.sim)), nrow = nrow(mymodel.dat))
for(i in 1:length(x.sim)){
  int.sim[, i] <- mymodel.dat$a * (x.sim[i]-mymodel.dat$c) / 
    ((mymodel.dat$a/mymodel.dat$b)+(x.sim[i]-mymodel.dat$c))
  }
bayes.c.eff.mean <- apply(int.sim, 2, mean)
bayes.c.eff.lower <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.025)))
bayes.c.eff.upper <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.975)))
## add to plot
lines(x.sim, bayes.c.eff.mean, col="red")
lines(x.sim, bayes.c.eff.lower, lty="dashed", col="red")
lines(x.sim, bayes.c.eff.upper, lty="dashed", col="red")
############################################################################


