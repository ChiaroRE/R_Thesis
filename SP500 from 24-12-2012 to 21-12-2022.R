#Descriptive Time series analysis of the S&P500 Index from approximately 2013 - 2022###################################

# IMPORT DATA #################################################################################################
setwd("C:/Users/edoar/Documents/R_Thesis/Data")
SP500_data = read.csv("S&P500, 24-12-2012 to 21-12-2022_CLEAN.csv", header = TRUE)

## PRICES, TIME AND LENGHT ####################################################################################
Pt <- SP500_data[,6] 
Time <- as.Date(as.character(SP500_data[,1]), "%m/%d/%Y")
N <- length(Pt)
pt <- log(Pt)            # log prices for computation of log returns

### RETURNS ###################################################################################################
Rt <- diff(Pt, lag = 1)/Pt[-1]
n <- length(Rt)

plot(Time[-1],Rt,
     type = "l",
     main = "S&P500 Returns: 25/12/2012 - 21/12/2022",
     xlab = "Time",
     ylab = "Returns")

summary(Rt)

#### LOG RETURNS ##############################################################################################
rt <- diff(pt, lag = 1)

plot(Time[-1],rt,
     type = "l",
     main = "S&P500 Log Returns: 25/12/2012 - 21/12/2022",
     xlab = "Time",
     ylab = "Returns")

summary(rt)
cor(Rt,rt)                # correlation between Returns and Log Returns, expected good approximation for +-5%
plot(Rt,rt)

##### MOMENTS #################################################################################################
library(moments)
sd <- sd(rt)
skew <- skewness(rt)         
kurt <- kurtosis(rt)

# both skewness and kurtosis are signs of non-normality since skewness != 0 and kurtosis > 3, both significantly
# still less than what we can observe in the EUROSTOXX (especially excess kurtosis)

##### TEST FOR MEAN = 0 #######################################################################################
t.test(rt)
# the small p-value indicates the refusal of the null hypothesis of true mean = 0

###### TESTING FOR NORMALITY ##################################################################################
library(tseries)
jarque.bera.test(rt)
# the extremely small p-value signifies that we can refuse the null hypothesis of normality

qqnorm(rt,
       datax = TRUE,
       main = "Normal qq plot for LogReturns")
qqline(rt,
       datax = TRUE)
# as we could expect, while the sample and theoretical quantiles are approximately equal in the central range
# (around 0.00) the tails deviate from the theoretical quantiles. This indicates fatter tails than a 
# normal distributions as is often the case with financial data. We now proceed to test the hypothesis of a 
# t-distribution at different levels of degrees of freedom

####### TESTING FOR T-DISTRIBUTION ############################################################################
grid <- (1:n)/(n+1)
par(mfrow=c(1,4))

qqplot(rt, qt(grid,df=5), 
       main="t-plot, df = 5 ",
       xlab="Data",
       ylab="t-quantiles")

lmfit <- lm(qt(c(.25,.75),df=5) ~ quantile(rt,c(.25,.75)) )
abline(lmfit)

qqplot(rt, qt(grid,df=4), 
       main="t-plot, df = 4 ",
       xlab="Data",
       ylab="t-quantiles")

lmfit <- lm(qt(c(.25,.75),df=4) ~ quantile(rt,c(.25,.75)) )
abline(lmfit)

qqplot(rt, qt(grid,df=3),
       main="t-plot, df = 3",
       xlab="Data",
       ylab="t-quantiles")
lmfit <- lm(qt(c(.25,.75),df=3) ~ quantile(rt,c(.25,.75)) )
abline(lmfit)

qqplot(rt, qt(grid,df=2),
       main="t-plot, df = 2",
       xlab="Data",
       ylab="t-quantiles")
lmfit <- lm(qt(c(.25,.75),df=2) ~ quantile(rt,c(.25,.75)) )
abline(lmfit)

# Between the three alternatives, it seems that a t-distribution with 4 degrees of freedom is the best approximation
# for lower quantiles while a t-distribution with 3 degrees of freedom is a better suit for higher quantiles
# It becomes then difficult to establish the best approximation. Both will be tested during the KDE.

####### KERNEL DENSITY ESTIMATION ##############################################################################
# We now try to perform a Kernel density estimation and comparing it with random data generated according to a
# t-distribution with 4 degrees of freedom and a t-distribution with 3 degrees of freedom

par(mfrow=c(1,3))
KDE <- density(rt)

plot(KDE,
     main = "KDE",
     col = "red")

sim1 <- rt(n, df = 4) * sqrt(var(rt) * 2/4) + mean(rt) #simulate t dist with rt parameters rt*sqrt(var * df - 2/df ) + mean
KDEsim1 <- density(sim1)
sim2 <- rt(n, df = 3) * sqrt(var(rt) * 1/3) + mean(rt) #simulate t dist with rt parameters rt*sqrt(var * df - 2/df ) + mean
KDEsim2 <- density(sim2)

plot(KDEsim1,
     main = "KDE Simulation",
     col = "red")
plot(KDEsim2,
     main = "KDE Simulation",
     col = "red")

# We try to generate random data based on a t-distribution (shifted by the log returns parameters).
# It is still difficult to determine which df are better suited for this distribution as both seem to replicate
# different aspects of the original distriution 

######## TEST FOR SERIAL CORRELATION #########################################################################
par(mfrow=c(1,2))
acf(rt)
pacf(rt)       
Box.test(rt, lag = 5, type = c("Ljung-Box")) 
Box.test(rt, lag = 7, type = c("Ljung-Box")) 
Box.test(rt, lag = 9, type = c("Ljung-Box")) 
Box.test(rt, lag = 21, type = c("Ljung-Box")) 
Box.test(rt, lag = 23, type = c("Ljung-Box")) 
Box.test(rt, lag = 30, type = c("Ljung-Box")) 

# Using the Ljung-Box test, we test the null hypotheis of serial independence. As we can see from the 
# correlogram, there are different lags at which the serial correlation could be considered statistically
# significant. 
# We can as a matter of fact refuse the null hypothesis of serial independence
# 1% level: 
# lag 5 - lag 7 - lag 9 - lag 21 - lag 23 - lag 30 
# This could be a starting point for eventual modelling

rm(list = ls())
