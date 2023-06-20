# IMPORT DATA #################################################################################################
setwd("C:/Users/edoar/OneDrive/Desktop/Thesis/R_Thesis/Data/")
EUROSTOX_data = read.csv("EUROSTOXX50, 24-12-2012 to 21-12-2022_CLEAN.csv", header = TRUE)

## PRICES, TIME AND LENGHT ####################################################################################
Pt_EUR <- EUROSTOX_data[,6] 
Time_EUR <- as.Date(as.character(EUROSTOX_data[,1]), "%d/%m/%Y")
N_EUR <- length(Pt_EUR)
pt_EUR <- log(Pt_EUR)            # log prices for computation of log returns

### RETURNS ###################################################################################################
Rt_EUR <- diff(Pt_EUR, lag = 1)/Pt_EUR[-1]
n_EUR <- length(Rt_EUR)

plot(Time_EUR[-1],Rt_EUR,
     type = "l",
     main = "EUROSTOXX50 Returns: 25/12/2012 - 21/12/2022",
     xlab = "Time",
     ylab = "Returns")

summary(Rt_EUR)

#### LOG RETURNS ##############################################################################################
rt_EUR <- diff(pt_EUR, lag = 1)

plot(Time_EUR[-1],rt_EUR,
     type = "l",
     main = "EUROSTOXX50 Log Returns: 25/12/2012 - 21/12/2022",
     xlab = "Time",
     ylab = "Returns")

summary(rt_EUR)
cor(Rt_EUR,rt_EUR)                # correlation between Returns and Log Returns, expected good approximation for +-5%
plot(Rt_EUR,rt_EUR)

##### MOMENTS #################################################################################################
library(moments)
sd_EUR <- sd(rt_EUR)
skew_EUR <- skewness(rt_EUR)         
kurt_EUR <- kurtosis(rt_EUR)

# both skewness and kurtosis are signs of non-normality since skewness != 0 and kurtosis > 3, both significantly

##### TEST FOR MEAN = 0 #######################################################################################
t.test(rt_EUR)
# the high p-value indicates the impossibility to refuse the null hypothesis of true mean = 0

###### TESTING FOR NORMALITY ##################################################################################
library(tseries)
jarque.bera.test(rt_EUR)
# the extremely small p-value signifies that we can refuse the null hypothesis of normality

qqnorm(rt_EUR,
       datax = TRUE,
       main = "Normal qq plot for LogReturns")
qqline(rt_EUR,
       datax = TRUE)
# as we could expect, while the sample and theoretical quantiles are approximately equal in the central range
# (around 0.00) the tails deviate from the theoretical quantiles. This indicates fatter tails than a 
# normal distributions as is often the case with financial data. We now proceed to test the hypothesis of a 
# t-distribution at different levels of degrees of freedom

####### TESTING FOR T-DISTRIBUTION ############################################################################
grid <- (1:n)/(n+1)
par(mfrow=c(1,4))

qqplot(rt_EUR, qt(grid,df=5), 
       main="t-plot, df = 5 ",
       xlab="Data",
       ylab="t-quantiles")

lmfit <- lm(qt(c(.25,.75),df=5) ~ quantile(rt_EUR,c(.25,.75)) )
abline(lmfit)

qqplot(rt_EUR, qt(grid,df=4), 
       main="t-plot, df = 4 ",
       xlab="Data",
       ylab="t-quantiles")

lmfit <- lm(qt(c(.25,.75),df=4) ~ quantile(rt_EUR,c(.25,.75)) )
abline(lmfit)

qqplot(rt_EUR, qt(grid,df=3),
       main="t-plot, df = 3",
       xlab="Data",
       ylab="t-quantiles")
lmfit <- lm(qt(c(.25,.75),df=3) ~ quantile(rt_EUR,c(.25,.75)) )
abline(lmfit)

qqplot(rt_EUR, qt(grid,df=2),
       main="t-plot, df = 2",
       xlab="Data",
       ylab="t-quantiles")
lmfit <- lm(qt(c(.25,.75),df=2) ~ quantile(rt_EUR,c(.25,.75)) )
abline(lmfit)

# Between the four alternatives, it seems that a t-distribution with 4 degrees of freedom is the best approximation
# However we should notice that the higher quantiles seems to deviate significantly, still the smallest deviation
# compared to other degrees of freedom and the normal distribution

####### KERNEL DENSITY ESTIMATION ##############################################################################
# We now try to perform a Kernel density estimation and comparing it with random data generated according to a
# t-distribution with 4 degrees of freedom

par(mfrow=c(1,2))
KDE_EUR <- density(rt_EUR)

plot(KDE_EUR,
     main = "KDE",
     col = "red")

sim1_EUR <- rt(n_EUR, df = 4) * sqrt(var(rt_EUR) * 2/4) + mean(rt_EUR) #simulate t dist with rt parameters rt*sqrt(var * df - 2/df ) + mean
KDEsim_EUR <- density(sim1_EUR)

plot(KDEsim_EUR,
     main = "KDE Simulation",
     col = "red")

# We try to generate random data based on a t-distribution (shifted by the log returns parameters).
# It seems that the simulation gives a good approximation, it should be noticed however that the kurtosis
# of the simulated data seems to be different from our original density estimation

######## DISTRIBUTION ESTIMATION ##############################################
#We now use another approach to estimate the parameters of a t-distribution to our
#dataset. Then we compare the data plot with the estimation
library(MASS)
library(fitdistrplus) #fitting t-dist
library(metRology) #addition of t-scaled distrbiuton
fit_t_EUR <- fitdist(rt_EUR,"t.scaled", start = list(df = 4, mean = mean(rt_EUR), sd = sd(rt_EUR)))
summary(fit_t_EUR)
df_t_EUR = fit_t_EUR$estimate[1]

#Let's try to test this estimate and compare it with the previous chosen 3 dfs
par(mfrow=c(1,2))

qqplot(rt_EUR, qt(grid,df=4), 
       main="t-plot, df = 4 ",
       xlab="Data",
       ylab="t-quantiles")

lmfit <- lm(qt(c(.25,.75),df=4) ~ quantile(rt_EUR,c(.25,.75)) )
abline(lmfit)

qqplot(rt_EUR, qt(grid,df= df_t_EUR), 
       main="t-plot, df estimated ",
       xlab="Data",
       ylab="t-quantiles")

lmfit <- lm(qt(c(.25,.75),df= df_t_EUR) ~ quantile(rt_EUR,c(.25,.75)) )
abline(lmfit)

#visually, it seems that the estimated dfs overperform the 3 dfs. However, we should still note that the theoretical quantiles
#seem to better capture the lower tail of the distribution. These are usually our quantiles of interest
#we proceed with a KDE to close the test
par(mfrow=c(1,3))

plot(KDE_EUR,
     main = "KDE",
     col = "red")

sim1_est_EUR <- rt(n_EUR, df = 4) * sqrt(var(rt_EUR) * 2/4) + mean(rt_EUR) #simulate t dist with rt parameters rt*sqrt(var * df - 2/df ) + mean
KDEsim1_est_EUR<- density(sim1_est_EUR)
sim2_est_EUR <- rt(n_EUR, df = df_t_EUR) * sqrt(var(rt_EUR) * (df_t_EUR - 2)/df_t_EUR) + mean(rt_EUR) #simulate t dist with rt parameters rt*sqrt(var * df - 2/df ) + mean
KDEsim2_est_EUR <- density(sim2_est_EUR)

plot(KDEsim1_est_EUR,
     main = "KDE 4 df",
     col = "red")
plot(KDEsim2_est_EUR,
     main = "KDE df estimated",
     col = "blue")



######## TEST FOR SERIAL CORRELATION #########################################################################
par(mfrow=c(1,2))
acf(rt_EUR)
pacf(rt_EURt)       
Box.test(rt_EUR, lag = 7, type = c("Ljung-Box")) 
Box.test(rt_EUR, lag = 13, type = c("Ljung-Box")) 
Box.test(rt_EUR, lag = 15, type = c("Ljung-Box")) 
Box.test(rt_EUR, lag = 23, type = c("Ljung-Box")) 

# Using the Ljung-Box test, we test the null hypotheis of serial independence. As we can see from the 
# correlogram, there are different lags at which the serial correlation could be considered statistically
# significant. 
# We can as a matter of fact refuse the null hypothesis of serial independence
# 5% level:
# lag 7 - lag 13
# 1% level: 
# lag 15 - lag 23
# This could be a starting point for eventual modelling

rm(list = ls())

