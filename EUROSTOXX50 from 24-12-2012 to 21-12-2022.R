#Descriptive Time series analysis of the EURO STOXX50 Index from approximately 2013 - 2022###################################

# IMPORT DATA #################################################################################################
setwd("C:/Users/edoar/Documents/R_Thesis/Data")
EUROSTOX_data = read.csv("EUROSTOXX50, 24-12-2012 to 21-12-2022_CLEAN.csv", header = TRUE)

## PRICES, TIME AND LENGHT ####################################################################################
Pt <- EUROSTOX_data[,6] 
Time <- as.Date(as.character(EUROSTOX_data[,1]), "%d/%m/%Y")
N <- length(Pt)
pt <- log(Pt)            # log prices for computation of log returns

### RETURNS ###################################################################################################
Rt <- diff(Pt, lag = 1)/Pt[-1]
n <- length(Rt)

plot(Time[-1],Rt,
     type = "l",
     main = "EUROSTOXX50 Returns: 25/12/2012 - 21/12/2022",
     xlab = "Time",
     ylab = "Returns")

summary(Rt)

#### LOG RETURNS ##############################################################################################
rt <- diff(pt, lag = 1)

plot(Time[-1],rt,
     type = "l",
     main = "EUROSTOXX50 Log Returns: 25/12/2012 - 21/12/2022",
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

##### TEST FOR MEAN = 0 #######################################################################################
t.test(rt)
# the high p-value indicates the impossibility to refuse the null hypothesis of true mean = 0

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
par(mfrow=c(1,3))

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

# Between the three alternatives, it seems that a t-distribution with 4 degrees of freedom is the best approximation
# However we should notice that the higher quantiles seems to deviate significantly, still the smallest deviation
# compared to other degrees of freedom and the normal distribution

####### KERNEL DENSITY ESTIMATION ##############################################################################
# We now try to perform a Kernel density estimation and comparing it with random data generated according to a
# t-distribution with 4 degrees of freedom

par(mfrow=c(1,2))
KDE <- density(rt)

plot(KDE,
     main = "KDE",
     col = "red")

sim1 <- rt(n, df = 4) * sqrt(var(rt) * 2/4) + mean(rt) #simulate t dist with rt parameters rt*sqrt(var * df - 2/df ) + mean
KDEsim <- density(sim1)

plot(KDEsim,
     main = "KDE Simulation",
     col = "red")

# We try to generate random data based on a t-distribution (shifted by the log returns parameters).
# It seems that the simulation gives a good approximation, it should be noticed however that the kurtosis
# of the simulated data seems to be different from our original density estimation

######## TEST FOR SERIAL CORRELATION #########################################################################
par(mfrow=c(1,2))
acf(rt)
pacf(rt)       
Box.test(rt, lag = 7, type = c("Ljung-Box")) 
Box.test(rt, lag = 13, type = c("Ljung-Box")) 
Box.test(rt, lag = 15, type = c("Ljung-Box")) 
Box.test(rt, lag = 23, type = c("Ljung-Box")) 

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

