## ----rendering, eval=FALSE, echo=FALSE-----------------------------------
## library(rmarkdown)
## system.time(render("uncertaintyProp_Practical1.Rmd", output_file = "uncertaintyProp_Practical1.html"))
## system.time(render("uncertaintyProp_Practical1.Rmd", output_file = "uncertaintyProp_Practical1.pdf"))
## knitr::purl("uncertaintyProp_Practical1.Rmd")

## ----load_data, eval=TRUE, echo=TRUE, warnings=FALSE, messages=FALSE-----
rm(list=ls(all=TRUE))
library(ggplot2)
library(openair)
library(rstanarm)
load(file = "streamflow.RData", verbose = TRUE)

## ---- eval=TRUE, echo=TRUE-----------------------------------------------
# df is the common generic data frame name
summary(df)
timePlot(df, pollutant = c("P_stream")) 

## ---- eval=TRUE, echo=TRUE-----------------------------------------------
summary(df_calib)
plot(Q ~ P_stream, data = df_calib)
m_lm <- lm(Q ~ P_stream, data = df_calib)
abline(m_lm, col = "blue")

## ---- eval=TRUE, echo=TRUE-----------------------------------------------
plot(sqrt(Q) ~ P_stream, data = df_calib)

## Fit a linear model
m_lm <- lm(sqrt(Q) ~ P_stream, data = df_calib)
abline(m_lm, col = "blue")

## ---- eval=TRUE, echo=TRUE-----------------------------------------------
df$Qpred_lm <- predict(m_lm, newdata = df)^2

## ---- eval=TRUE, echo=TRUE-----------------------------------------------
timePlot(df, pollutant = c("P_stream", "Qpred_lm"), scales = "free")
summary(df)

## ---- eval=FALSE, echo=TRUE----------------------------------------------
## m <- stan_glm(sqrt(Q) ~ P_stream, data = df_calib)

## ---- eval=TRUE, echo=FALSE, include=FALSE, warnings=FALSE, messages=FALSE----
m <- stan_glm(sqrt(Q) ~ P_stream, data = df_calib)

## ---- eval=TRUE, echo=TRUE-----------------------------------------------
prior_summary(m)

## ---- eval=TRUE, echo=TRUE-----------------------------------------------
plot(m, "trace")  

## ---- eval=TRUE, echo=TRUE-----------------------------------------------
plot(m, "hist")

## ---- eval=FALSE, echo=TRUE----------------------------------------------
## launch_shinystan(m)

## ---- eval=TRUE, echo=TRUE-----------------------------------------------
p <- posterior_vs_prior(m, group_by_parameter = TRUE, 
  color_by = "vs", facet_args = list(scales = "free"))
p <- p + 
 geom_hline(yintercept = 0, size = 0.3, linetype = 3) + 
 coord_flip() + 
 ggtitle("Comparing the prior and posterior")
p <- p + geom_pointrange(size = 0.2)
p

posterior_interval(m, 0.90)

## ---- eval=TRUE, echo=TRUE-----------------------------------------------
# make a new data frame without missing values
df_pred <- na.omit(df)
# and predict with the parameter distributions estimated by stan_glm
Q_pred <- posterior_predict(m, newdata = df_pred)^2

## ---- eval=TRUE, echo=TRUE-----------------------------------------------
df_pred$Q_pred1 <- Q_pred[1,] 
df_pred$Q_pred2 <- Q_pred[2,] 
df_pred$Q_pred3 <- Q_pred[3,] 
timePlot(df_pred, pollutant = c("Q_pred1", "Q_pred2", "Q_pred3"), ylab = "Q_flow", xlab = "Date, 2017")
summary(df_pred)

## ---- eval=TRUE, echo=TRUE-----------------------------------------------
df_pred$Q_mean   <- apply(Q_pred, 2, mean)
df_pred$quant_05 <- apply(Q_pred, 2, quantile, c(0.05))
df_pred$quant_50 <- apply(Q_pred, 2, quantile, c(0.50))
df_pred$quant_95 <- apply(Q_pred, 2, quantile, c(0.95))
p <- ggplot(df_pred, aes(date, y = Q_mean, ymin = quant_05, ymax = quant_95))
p <- p + geom_ribbon(fill = "blue")
p <- p + geom_line(colour = "red")
p

## ---- eval=TRUE, echo=TRUE-----------------------------------------------
# conventional mean and uncertainty
mean(df$Qpred_lm, na.rm = TRUE)
# SE in mean = sd/sqrt(n)
se <- sd(df$Qpred_lm, na.rm = TRUE) / sqrt(sum(!is.na(df$Qpred_lm)))
se; se^2 # ^2 to give in variance units

## ---- eval=TRUE, echo=TRUE-----------------------------------------------
# get the mean of each realisation
Q_pred_means <- apply(Q_pred, 1, mean)
hist(Q_pred_means)

## ---- eval=TRUE, echo=TRUE-----------------------------------------------
sd(Q_pred_means)
# a coefficient of variation
sd(Q_pred_means)/ mean(Q_pred_means) *100

## ---- eval=TRUE, echo=TRUE-----------------------------------------------
# get cumulative sum for each realisation
Q_pred_sum <- apply(Q_pred, 1, cumsum) / 1000 # l to m3

df_pred$Q_sum        <- apply(Q_pred_sum, 1, mean)
df_pred$sum_quant_05 <- apply(Q_pred_sum, 1, quantile, c(0.05))
df_pred$sum_quant_50 <- apply(Q_pred_sum, 1, quantile, c(0.50))
df_pred$sum_quant_95 <- apply(Q_pred_sum, 1, quantile, c(0.95))
p <- ggplot(df_pred, aes(date, y = Q_sum, ymin = sum_quant_05, ymax = sum_quant_95))
p <- p + geom_ribbon(fill = "light blue")
p <- p + geom_line(colour = "red")
p

## ---- eval=TRUE, echo=TRUE-----------------------------------------------
hist(Q_pred_sum[dim(Q_pred_sum)[1],], xlab = "Q_sum, m3 y-1", main = 
 "Posterior distribution of cumulative annual stream flow, 2017")
quantile(Q_pred_sum[dim(Q_pred_sum)[1],], c(0.025, 0.5, 0.975))
sd(Q_pred_sum[dim(Q_pred_sum)[1],]) /
mean(Q_pred_sum[dim(Q_pred_sum)[1],]) *100

