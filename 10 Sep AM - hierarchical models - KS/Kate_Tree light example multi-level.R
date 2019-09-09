
setwd("N://STATS/CEH Bayes Course/")

library(rjags)
library(jagsUI)

rm(list=ls(all=TRUE)) # clear the R slate!

#############################################################################
### load in data
tree.data <- read.csv("Hemlock-light-data-hierarchical_2.csv")
str(tree.data)
#############################################################################


#############################################################################
### do some plotting to look at data and fit simple non-linear model

colour <- c("red", "blue", "green")
col.list <- rep(0, length(tree.data$Site))
col.list[tree.data$Site=="1"] <- 1 
col.list[tree.data$Site=="2"] <- 2
col.list[tree.data$Site=="3"] <- 3
plot(tree.data$Light, tree.data$Observed.growth.rate, col=c(colour[col.list]))

#############################################################################


#############################################################################
## set up for running model in JAGS

# Initialize the MCMC chains. This needs to be a list.
inits <- function() list(a=rep(30,3), b=runif(1,1,3), c=runif(1,-5,5), 
                         mu.a=30, sigma.o=1, sigma.p=1, sigma.a=1) #small tau makes variance big
inits

## create dataframe for y data (observed growth rate)
head(tree.data)
y.dat <- array(NA, c(77,3))
y.dat[1:77,1] <- tree.data$Observed.growth.rate[which(tree.data$Site==1)]
y.dat[1:77,2] <- tree.data$Observed.growth.rate[which(tree.data$Site==2)]
y.dat[1:77,3] <- tree.data$Observed.growth.rate[which(tree.data$Site==3)]
head(y.dat)

## create dataframe for x data (Light)
head(tree.data)
x.dat <- array(NA, c(77,3))
x.dat[1:77,1] <- tree.data$Light[which(tree.data$Site==1)]
x.dat[1:77,2] <- tree.data$Light[which(tree.data$Site==2)]
x.dat[1:77,3] <- tree.data$Light[which(tree.data$Site==3)]
head(x.dat)

# Specify data; must be a list. 
data=list(
  n=77,
  x=x.dat,
  y=y.dat
)
data

# Parameters monitored (i.e., for which estimates are saved)
params <- c("a","b","c","sigma.o","sigma.p", "sigma.a", "mu.a")

# MCMC settings
ni <- 75000   ;   nt <- 1   ;   nb <- 50000   ;  nc <- 3

### fit the model in JAGS
model_2 <- jags(data=data, inits, params, 
                model.file = "Kate_jags_light_example multi-level.R", n.chains = nc, 
                n.thin = nt, n.iter = ni, n.burnin = nb)

summary(model_2)
print(model_2, 6) # the 6 sets the number of significant digits to print out
plot(model_2) # check convergence


#############################################################################
## extract the MCMC estimates for each model parameter
names(model_2)
str(model_2$sims.list)
mymodel.dat <- as.data.frame(model_2$sims.list)
head(mymodel.dat)
## plot observed data against predictions
par(mfrow=c(1,1))
plot(tree.data$Light, tree.data$Observed.growth.rate, col=c(colour[col.list]))
#############################################################################


#############################################################################
## plot observed against predicted using all MCMC output

## Simulate the range of the environmental variable (Light)
x.sim <- seq(min(tree.data$Light), max(tree.data$Light), by = 0.1)
#############################################################################
## Calculate prediction over all mcmc samples for each parameter
## site 1
int.sim <- matrix(rep(NA, nrow(mymodel.dat)*length(x.sim)), nrow = nrow(mymodel.dat))
for(i in 1:length(x.sim)){
  int.sim[, i] <- mymodel.dat$a.1 * (x.sim[i]-mymodel.dat$c) / 
    ((mymodel.dat$a.1/mymodel.dat$b)+(x.sim[i]-mymodel.dat$c))
}
bayes.c.eff.mean <- apply(int.sim, 2, mean)
bayes.c.eff.lower <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.025)))
bayes.c.eff.upper <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.975)))
## add to plot
lines(x.sim, bayes.c.eff.mean, col="red")
lines(x.sim, bayes.c.eff.lower, lty="dashed", col="red")
lines(x.sim, bayes.c.eff.upper, lty="dashed", col="red")
############################################################################
## site 2
for(i in 1:length(x.sim)){
  int.sim[, i] <- mymodel.dat$a.2 * (x.sim[i]-mymodel.dat$c) / 
    ((mymodel.dat$a.2/mymodel.dat$b)+(x.sim[i]-mymodel.dat$c))
}
bayes.c.eff.mean <- apply(int.sim, 2, mean)
bayes.c.eff.lower <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.025)))
bayes.c.eff.upper <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.975)))
## add to plot
lines(x.sim, bayes.c.eff.mean, col="blue")
lines(x.sim, bayes.c.eff.lower, lty="dashed", col="blue")
lines(x.sim, bayes.c.eff.upper, lty="dashed", col="blue")
############################################################################
## site 3
for(i in 1:length(x.sim)){
  int.sim[, i] <- mymodel.dat$a.3 * (x.sim[i]-mymodel.dat$c) / 
    ((mymodel.dat$a.3/mymodel.dat$b)+(x.sim[i]-mymodel.dat$c))
}
bayes.c.eff.mean <- apply(int.sim, 2, mean)
bayes.c.eff.lower <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.025)))
bayes.c.eff.upper <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.975)))
## add to plot
lines(x.sim, bayes.c.eff.mean, col="green")
lines(x.sim, bayes.c.eff.lower, lty="dashed", col="green")
lines(x.sim, bayes.c.eff.upper, lty="dashed", col="green")
############################################################################