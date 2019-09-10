## ----rendering, eval=FALSE, echo=FALSE-----------------------------------
## library(rmarkdown)
## system.time(render("combining_Data_Practical.Rmd", output_file = "combining_Data_Practical.html"))
## system.time(render("combining_Data_Practical.Rmd", output_file = "combining_Data_Practical.pdf"))
## system.time(render("combining_Data_Practical.Rmd", output_file = "combining_Data_Practical.docx"))
## knitr::purl("combining_Data_Practical.Rmd")

## ----start_up, eval=TRUE, echo=FALSE, include=FALSE----------------------
rm(list=ls(all=TRUE))
#install.packages("BayesianTools")
library(knitr)
load(file = "landUseChange.RData", verbose = TRUE)

## ---- eval=TRUE, echo=FALSE----------------------------------------------
knitr::kable(m_B_true, caption = "Area changing from land-use in row $i$ to land-use in column $j$ in km$^2$.")

## ---- eval=FALSE, echo=TRUE----------------------------------------------
##   dA_net <- colSums(m_B) - rowSums(m_B)

## ---- eval=FALSE, echo=TRUE----------------------------------------------
## A_t1 <- rowSums(m_B)
## A_t2 <- colSums(m_B)

## ---- eval=FALSE, echo=TRUE----------------------------------------------
## dA_gain_obs <- colSums(m_B) - diag(m_B)
## dA_loss_obs <- rowSums(m_B) - diag(m_B)

## ---- eval=TRUE, echo=TRUE-----------------------------------------------
# predict the net change in area for each land use from a matrix
getAreaNetChange_fromBeta <- function(v_B){
  m_B <- matrix(v_B, n_u, n_u)
  dA_net <- colSums(m_B) - rowSums(m_B)
  return(dA_net)  
}

# calculate the gross changes in area for each land use in a Beta transition matrix
getAreaGrossChange_fromBeta <- function(v_B){
  m_B <- matrix(v_B, n_u, n_u)
  # subtract the unchanged area (diagonal) to get gross gains and losses
  A_gain <- colSums(m_B) - diag(m_B)
  A_loss <- rowSums(m_B) - diag(m_B)
  return(rbind(A_gain, A_loss))  
}

# calculate the area at t1 and t2 for each land use in a Beta transition matrix
getAreas_fromBeta <- function(v_B){
  m_B <- matrix(v_B, n_u, n_u)
  # rows add to the time1 area; cols add to the time2 area
  A_t1 <- rowSums(m_B)
  A_t2 <- colSums(m_B)
  return(rbind(A_t1, A_t2))  
}

## ----load_data, eval=TRUE, echo=TRUE-------------------------------------
library(BayesianTools)
set.seed(448)
load(file = "landUseChange.RData", verbose = TRUE)

## ---- eval=TRUE, echo=TRUE-----------------------------------------------
v_B_prior <- as.vector(m_B_prior)
prior <- createUniformPrior(
  best = v_B_prior, 
  lower = rep(0, n_u^2), 
  upper = rep(50, n_u^2))
  
# Initial values for chains
m_starter <- matrix(rep(m_B_prior, 3), nrow = 3, byrow = TRUE)
m_starter[1,] <- m_starter[1,] * runif(n_u^2, 0.5, 1.5)
m_starter[3,] <- m_starter[3,] * runif(n_u^2, 0.5, 1.5)

## ---- eval=TRUE, echo=TRUE-----------------------------------------------
# calculates the log-likelihood of B matrix given observations 
# of net area change 
getLogLikelihood_net <- function(v_B){
  # net change
  dA_net_pred <- getAreaNetChange_fromBeta(v_B)
  # use constant CV for sigma in obs for simple example
  loglik_net <- dnorm(dA_net_obs, mean = dA_net_pred, 
                sd = max(1, 0.1*abs(dA_net_obs)), log = T)
  return(sum(loglik_net, na.rm = TRUE))  
}

## ---- eval=FALSE, echo=TRUE----------------------------------------------
## loglik_net <- dnorm(dA_net_obs, mean = dA_net_pred,
##                 sd = max(1, 0.1*abs(dA_net_obs)), log = T)

## ---- eval=TRUE, echo=FALSE----------------------------------------------
x <- seq(-2, 2, by = 0.1)
plot(x, dnorm(x, mean = 0, sd = 0.1, log = T), type = "l",
  xlab = "Difference between observation and prediction",
  ylab = "Log-likelihood")

## ---- eval=TRUE, echo=TRUE-----------------------------------------------
setUp <- createBayesianSetup(getLogLikelihood_net, prior = prior)

n_iter <- 9000 # set the number of MCMC iterations
thin <- round(max(1, n_iter/1000))
n_chains <- 3  # use three chains
settings <- list(iterations = n_iter,  thin = thin,  nrChains = n_chains, startValue = m_starter, message = FALSE)
system.time(out <- runMCMC(bayesianSetup = setUp, sampler = "DEzs", settings = settings))

## ---- eval=TRUE, echo=TRUE-----------------------------------------------
plot(out, which = 1:2)
#marginalPlot(out, whichParameters = 1:2)
# plot(out, which = 5:8)
# plot(out, which = 9:12)
# plot(out, which = 13:16)
#gelmanDiagnostics(out, plot = TRUE)

## ---- eval=TRUE, echo=TRUE-----------------------------------------------
v_B_map <- MAP(out)$parametersMAP
getAreaNetChange_fromBeta(v_B_map)
dA_net_obs
# save MAP prediction values to array
m_B_map <- matrix(v_B_map, n_u, n_u)  
round(m_B_map, 2)
round(m_B_true, 2)
plot(m_B_true, m_B_map)
abline(0, 1)

## ---- eval=TRUE, echo=TRUE-----------------------------------------------
# calculates the log-likelihood of B matrix given observations 
# of areas at t1 and t2
getLogLikelihood_Areas_t1t2 <- function(v_B){
  # get total area for each land use at time 1 and 2
  A_pred <- getAreas_fromBeta(v_B)
  loglik_A_t1  <-  dnorm(A_obs[1,], mean = A_pred[1,], sd = max(1, 0.1*abs(A_obs[1,])), log = TRUE)
  loglik_A_t2  <-  dnorm(A_obs[2,], mean = A_pred[2,], sd = max(1, 0.1*abs(A_obs[2,])), log = TRUE)
  return(sum(loglik_A_t1, loglik_A_t2, na.rm = TRUE))  
}

# calculates the log-likelihood of B matrix given observations 
# of net area change and areas as t1 and t2
getLogLikelihood_net_t1t2 <- function(v_B){
  loglik_net <- getLogLikelihood_net(v_B)
  loglik_t1t2 <- getLogLikelihood_Areas_t1t2(v_B)
  return(sum(loglik_net, loglik_t1t2, na.rm = TRUE))  
}

## ---- eval=TRUE, echo=TRUE-----------------------------------------------
setUp <- createBayesianSetup(getLogLikelihood_net_t1t2, prior = prior)
system.time(out <- runMCMC(bayesianSetup = setUp, sampler = "DEzs", settings = settings))
plot(out, which = 1:2)
v_B_map <- MAP(out)$parametersMAP
getAreaNetChange_fromBeta(v_B_map)
dA_net_obs
m_B_map <- matrix(v_B_map, n_u, n_u)  
round(m_B_map, 2)
round(m_B_true, 2)
plot(m_B_true, m_B_map)
abline(0, 1)

## ---- eval=TRUE, echo=TRUE-----------------------------------------------
# calculates the log-likelihood of B matrix given observations 
# of gross area change 
getLogLikelihood_gross <- function(v_B){
  # gross changes
  dA_gross_pred <- getAreaGrossChange_fromBeta(v_B)
  loglik_gain  <-  dnorm(dA_gain_obs, mean = dA_gross_pred[1,], 
                     sd = max(1, 0.1*abs(dA_gain_obs)), log = TRUE)
  loglik_loss  <-  dnorm(dA_loss_obs, mean = dA_gross_pred[2,], 
                     sd = max(1, 0.1*abs(dA_loss_obs)), log = TRUE)
  return(sum(loglik_gain, loglik_loss, na.rm = TRUE))  
}

# calculates the log-likelihood of B matrix given observations 
# of net area change and areas as t1 and t2 and gross change
getLogLikelihood_net_t1t2_gross <- function(v_B){
  loglik_net <- getLogLikelihood_net(v_B)
  loglik_t1t2 <- getLogLikelihood_Areas_t1t2(v_B)
  loglik_gross <- getLogLikelihood_gross(v_B)
  return(sum(loglik_net, loglik_t1t2, loglik_gross, na.rm = TRUE))  
}

## ---- eval=TRUE, echo=TRUE-----------------------------------------------
setUp <- createBayesianSetup(getLogLikelihood_net_t1t2_gross, prior = prior)
system.time(out <- runMCMC(bayesianSetup = setUp, sampler = "DEzs", settings = settings))
plot(out, which = 1:2)
plot(out, which = 3:4)
#plot(out, which = 5:8)
v_B_map <- MAP(out)$parametersMAP
getAreaNetChange_fromBeta(v_B_map)
dA_net_obs
m_B_map <- matrix(v_B_map, n_u, n_u)  
round(m_B_map, 2)
round(m_B_true, 2)
plot(m_B_true, m_B_map)
abline(0, 1)

