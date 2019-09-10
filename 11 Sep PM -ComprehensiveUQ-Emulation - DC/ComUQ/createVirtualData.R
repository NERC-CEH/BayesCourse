library(BayesianTools)

# load reference parameter definition (upper, lower prior)
refPars <- VSEMgetDefaults()

# this adds one additional parameter for the likelihood standard deviation (see below)
refPars[12,] <- c(0.1, 0.001, 0.5) 
rownames(refPars)[12] <- "error-sd"

# create some virtual data 
set.seed(123)

ndays = 2048
PAR   <- VSEMcreatePAR(1:ndays)

referenceData <- VSEM(refPars$best[1:11], PAR) # model predictions with reference parameters  
referenceData[,1] = 1000 * referenceData[,1] # unit change for NEE from kg to g

nvar <- 12
obs                   <- referenceData + rnorm(length(referenceData),
                                  sd = (abs(referenceData) + 1E-7) * refPars$best[nvar])

## Set the minimum sd for NEE to 0.5 gC /m^2 /day.
## This is to avoid very small uncertainty for NEE values close to zero
obs[,1] <- referenceData[,1] + rnorm(length(referenceData[,1]),
                                     sd = pmax((abs(referenceData[,1]) + 1E-7) * refPars$best[nvar],0.5))

obsSel <- c(1,202,390,550,750,920)*2.0

save(file="calibrationData.RData", PAR, ndays, nvar, obs, obsSel, referenceData)

## Plot the calibration data
oldpar <- par(mfrow = c(2,2))
plot(obs[,1], main = colnames(referenceData)[1],pch=3)
plot(obs[obsSel,2], main = colnames(referenceData)[2],pch=3)
plot(obs[,3], main = colnames(referenceData)[3],pch=3)
