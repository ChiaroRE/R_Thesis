#Descriptive Time series analysis of the S&P500 Index from approximately 2013 - 2022###################################

# IMPORT DATA #################################################################################################
setwd("C:/Users/edoar/OneDrive/Desktop/Thesis/R_Thesis/Data")
SP500_data = read.csv("S&P500, 24-12-2012 to 21-12-2022_CLEAN.csv", header = TRUE)

## PRICES, TIME AND LENGHT ####################################################################################
Pt_SP <- SP500_data[,6] 
Time_SP <- as.Date(as.character(SP500_data[,1]), "%m/%d/%Y")
N_SP <- length(Pt_SP)
pt_SP <- log(Pt_SP)            # log prices for computation of log returns

### RETURNS ###################################################################################################
Rt_SP <- diff(Pt_SP, lag = 1)/Pt_SP[-1]
n_SP <- length(Rt_SP)

plot(Time_SP[-1],Rt_SP,
     type = "l",
     main = "S&P500 Returns: 25/12/2012 - 21/12/2022",
     xlab = "Time",
     ylab = "Returns")

summary(Rt_SP)

#### LOG RETURNS ##############################################################################################
rt_SP <- diff(pt_SP, lag = 1)

plot(Time_SP[-1],rt_SP,
     type = "l",
     main = "S&P500 Log Returns: 25/12/2012 - 21/12/2022",
     xlab = "Time",
     ylab = "Returns")

summary(rt_SP)
cor(Rt_SP,rt_SP)                # correlation between Returns and Log Returns, expected good approximation for +-5%
plot(Rt_SP,rt_SP)

##### MOMENTS #################################################################################################
library(moments)
sd_SP <- sd(rt_SP)
skew_SP <- skewness(rt_SP)         
kurt_SP <- kurtosis(Rt_SP)

# both skewness and kurtosis are signs of non-normality since skewness != 0 and kurtosis > 3, both significantly
# still less than what we can observe in the EUROSTOXX (especially excess kurtosis)

##### TEST FOR MEAN = 0 #######################################################################################
t.test(rt_SP)
# the small p-value indicates the refusal of the null hypothesis of true mean = 0

###### TESTING FOR NORMALITY ##################################################################################
library(tseries)
jarque.bera.test(rt_SP)
# the extremely small p-value signifies that we can refuse the null hypothesis of normality

qqnorm(rt_SP,
       datax = TRUE,
       main = "Normal qq plot for LogReturns")
qqline(rt_SP,
       datax = TRUE)
# as we could expect, while the sample and theoretical quantiles are approximately equal in the central range
# (around 0.00) the tails deviate from the theoretical quantiles. This indicates fatter tails than a 
# normal distributions as is often the case with financial data. We now proceed to test the hypothesis of a 
# t-distribution at different levels of degrees of freedom

####### TESTING FOR T-DISTRIBUTION ############################################################################
grid <- (1:n_SP)/(n_SP+1)
par(mfrow=c(1,4))

qqplot(rt_SP, qt(grid,df=5), 
       main="t-plot, df = 5 ",
       xlab="Data",
       ylab="t-quantiles")

lmfit <- lm(qt(c(.25,.75),df=5) ~ quantile(rt_SP,c(.25,.75)) )
abline(lmfit)

qqplot(rt_SP, qt(grid,df=4), 
       main="t-plot, df = 4 ",
       xlab="Data",
       ylab="t-quantiles")

lmfit <- lm(qt(c(.25,.75),df=4) ~ quantile(rt_SP,c(.25,.75)) )
abline(lmfit)

qqplot(rt_SP, qt(grid,df=3),
       main="t-plot, df = 3",
       xlab="Data",
       ylab="t-quantiles")
lmfit <- lm(qt(c(.25,.75),df=3) ~ quantile(rt_SP,c(.25,.75)) )
abline(lmfit)

qqplot(rt_SP, qt(grid,df=2),
       main="t-plot, df = 2",
       xlab="Data",
       ylab="t-quantiles")
lmfit <- lm(qt(c(.25,.75),df=2) ~ quantile(rt_SP,c(.25,.75)) )
abline(lmfit)

# Between the three alternatives, it seems that a t-distribution with 4 degrees of freedom is the best approximation
# for lower quantiles while a t-distribution with 3 degrees of freedom is a better suit for higher quantiles
# It becomes then difficult to establish the best approximation. Both will be tested during the KDE.

####### KERNEL DENSITY ESTIMATION ##############################################################################
# We now try to perform a Kernel density estimation and comparing it with random data generated according to a
# t-distribution with 4 degrees of freedom and a t-distribution with 3 degrees of freedom

par(mfrow=c(1,3))
KDE_SP <- density(rt_SP)

plot(KDE_SP,
     main = "KDE",
     col = "red")

sim1_SP <- rt(n_SP, df = 4) * sqrt(var(rt_SP) * 1/2) + mean(rt_SP) #simulate t dist with rt parameters rt*sqrt(var * df - 2/df ) + mean
KDEsim1_SP <- density(sim1_SP)
sim2_SP <- rt(n_SP, df = 3) * sqrt(var(rt_SP) * 1/3) + mean(rt_SP) #simulate t dist with rt parameters rt*sqrt(var * df - 2/df ) + mean
KDEsim2_SP <- density(sim2_SP)

plot(KDEsim1_SP,
     main = "KDE Simulation",
     col = "red")
plot(KDEsim2_SP,
     main = "KDE Simulation",
     col = "red")

# We try to generate random data based on a t-distribution (shifted by the log returns parameters).
# It is still difficult to determine which df are better suited for this distribution as both seem to replicate
# different aspects of the original distriution

######## DISTRIBUTION ESTIMATION ##############################################
#We now use another approach to estimate the parameters of a t-distribution to our
#dataset. Then we compare the data plot with the estimation
library(MASS)
library(fitdistrplus) #fitting t-dist
library(metRology) #addition of t-scaled distrbiuton
fit_t_SP <- fitdist(rt_SP,"t.scaled", start = list(df = 3, mean = mean(rt_SP), sd = sd(rt_SP)))
summary(fit_t_SP)
df_t_SP = fit_t_SP$estimate[1]

#Let's try to test this estimate and compare it with the previous chosen 3 dfs
par(mfrow=c(1,2))

qqplot(rt_SP, qt(grid,df=3), 
       main="t-plot, df = 3 ",
       xlab="Data",
       ylab="t-quantiles")

lmfit <- lm(qt(c(.25,.75),df=3) ~ quantile(rt_SP,c(.25,.75)) )
abline(lmfit)

qqplot(rt_SP, qt(grid,df= df_t_SP), 
       main="t-plot, df estimated ",
       xlab="Data",
       ylab="t-quantiles")

lmfit <- lm(qt(c(.25,.75),df= df_t_SP) ~ quantile(rt_SP,c(.25,.75)) )
abline(lmfit)

#visually, it seems that the estimated dfs overperform the 3 dfs. However, we should still note that the theoretical quantiles
#not seem to capture the lower tail of the distribution
#we proceed with a KDE to close the test
par(mfrow=c(1,3))

plot(KDE_SP,
     main = "KDE",
     col = "red")

sim1_est_SP <- rt(n_SP, df = 3) * sqrt(var(rt_SP) * 1/3) + mean(rt_SP) #simulate t dist with rt parameters rt*sqrt(var * df - 2/df ) + mean
KDEsim1_est_SP <- density(sim1_est_SP)
sim2_est_SP <- rt(n_SP, df = df_t) * sqrt(var(rt_SP) * (df_t - 2)/df_t) + mean(rt_SP) #simulate t dist with rt parameters rt*sqrt(var * df - 2/df ) + mean
KDEsim2_est_SP <- density(sim2_est_SP)

plot(KDEsim1_est_SP,
     main = "KDE 3 df",
     col = "red")
plot(KDEsim2_est_SP,
     main = "KDE df estimated",
     col = "blue")


######### TEST FOR SERIAL CORRELATION #########################################################################
par(mfrow=c(1,2))
acf(rt_SP)
pacf(rt_SP)       
Box.test(rt_SP, lag = 5, type = c("Ljung-Box")) 
Box.test(rt_SP, lag = 7, type = c("Ljung-Box")) 
Box.test(rt_SP, lag = 9, type = c("Ljung-Box")) 
Box.test(rt_SP, lag = 21, type = c("Ljung-Box")) 
Box.test(rt_SP, lag = 23, type = c("Ljung-Box")) 
Box.test(rt_SP, lag = 30, type = c("Ljung-Box")) 

# Using the Ljung-Box test, we test the null hypotheis of serial independence. As we can see from the 
# correlogram, there are different lags at which the serial correlation could be considered statistically
# significant. 
# We can as a matter of fact refuse the null hypothesis of serial independence
# 1% level: 
# lag 5 - lag 7 - lag 9 - lag 21 - lag 23 - lag 30 
# This could be a starting point for eventual modelling

rm(list = ls())
