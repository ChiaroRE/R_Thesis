#Descriptive Time series analysis of the S&P500 Index from approximately 2007 - 2009###################################

# IMPORT DATA #################################################################################################
setwd("C:/Users/edoar/OneDrive/Desktop/Thesis/R_Thesis/Data")
SP500_data = read.csv("S&P500, 01-01-2007 to 31-12-2009_CLEAN.csv", header = TRUE)

## PRICES, TIME AND LENGHT ####################################################################################
Pt_SP <- rev(SP500_data[,5])
Pt_SP
Time_SP <- rev(as.Date(as.character(SP500_data[,1]), "%d/%m/%Y"))
N_SP <- length(Pt_SP)
pt_SP <- log(Pt_SP)            # log prices for computation of log returns

### RETURNS ###################################################################################################
Rt_SP <- diff(Pt_SP, lag = 1)/Pt_SP[-1]
n_SP <- length(Rt_SP)

par(mfrow = c(1,1))
plot(Time_SP[-1],Rt_SP,
     type = "l",
     main = "S&P500 Returns: 25/12/2012 - 21/12/2022",
     xlab = "Time",
     ylab = "Returns")

summary(Rt_SP)

#outliers##
outliers_SP <- sum(Rt_SP > 0.15 | Rt_SP < -0.15)

z_scores <- (Rt_SP - mean(Rt_SP)) / sd(Rt_SP)

# Step 3: Set threshold (e.g., 2 or 3)
threshold <- 5

# Step 4: Identify outliers
outliers <- Rt_SP[abs(z_scores) > threshold]

# Print the outliers
print(outliers)


#### LOG RETURNS ##############################################################################################
rt_SP <- diff(pt_SP, lag = 1)

plot(Time_SP[-1],rt_SP,
     type = "l",
     main = "S&P500 Log-Returns: 25/12/2012 - 21/12/2022",
     xlab = "Time",
     ylab = "Returns",
     col = "blue")

summary(rt_SP)
cor(Rt_SP,rt_SP)                # correlation between Returns and Log Returns, expected good approximation for +-5%
plot(Rt_SP,rt_SP)

library(plotrix)
std.error(rt_SP)


range_SP <- max(rt_SP) - min(rt_SP)
range_SP

IQR_SP <- IQR(rt_SP)
##### STATIONARITY ###########################################################################################
library(timeSeries)
library(aTSA)
ADF_SP <- adf.test(rt_SP)
ADF_SP$p.value
ADF_SP


##### MOMENTS #################################################################################################
library(moments)
sd_SP <- sd(rt_SP)
se_sd_SP <- sd_SP/sqrt(2*n_SP - 2)
skew_SP <- skewness(rt_SP)
V_skew = 6 * n_SP * (n_SP - 1) / ((n_SP - 2) * (n_SP + 1) * (n_SP + 3))
sqrt(V_skew)
kurt_SP <- kurtosis(rt_SP)
V_kur = 4 * (n_SP^2 - 1) * V_skew / ((n_SP - 3) * (n_SP + 5))
sqrt(V_kur)


# both skewness and kurtosis are signs of non-normality since skewness != 0 and kurtosis > 3, both significantly
# still less than what we can observe in the SPOSTOXX (especially excess kurtosis)

##### TEST FOR MEAN = 0 #######################################################################################
t.test(rt_SP)
# the big p-value indicates the faliure to refuse the null hypothesis of true mean = 0

###### TESTING FOR NORMALITY ##################################################################################
library(tseries)
jarque.bera.test(rt_SP)
# the extremely small p-value signifies that we can refuse the null hypothesis of normality

qqnorm(rt_SP,
       datax = TRUE,
       main = "Normal qq plot for SP LogReturns")
qqline(rt_SP,
       datax = TRUE)
# as we could expect, while the sample and theoretical quantiles are approximately equal in the central range
# (around 0.00) the tails deviate from the theoretical quantiles. This indicates fatter tails than a 
# normal distribution as is often the case with financial data. We now proceed to test the hypothesis of a 
# t-distribution at different levels of degrees of freedom

####### TESTING FOR T-DISTRIBUTION ############################################################################
grid <- (1:n_SP)/(n_SP+1)
par(mfrow=c(1,3))

qqplot(rt_SP, qt(grid,df=4), 
       main="t-plot, df = 4 ",
       xlab="Data",
       ylab="t-quantiles")

lmfit <- lm(qt(c(.25,.75),df=4) ~ quantile(rt_SP,c(.25,.75)) )
abline(lmfit)

qqplot(rt_SP, qt(grid,df=3), 
       main="t-plot, df = 3 ",
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


# Between the three alternatives, it seems that a t-distribution with 2 degrees of freedom is the best approximation
#however, it doesn't appropriately capture the tails and especially the high quantiles

####### KERNEL DENSITY ESTIMATION ##############################################################################
# We now try to perform a Kernel density estimation and comparing it with random data generated according to a
# t-distribution with 2 degrees of freedom

par(mfrow=c(1,2))
KDE_SP <- density(rt_SP)

plot(KDE_SP,
     main = "KDE",
     col = "red")

sim2_SP <- rt(n_SP, df = 3) * sqrt(var(rt_SP) * (3 - 1)/3) + mean(rt_SP) #simulate t dist with rt parameters rt*sqrt(var * df - 2/df ) + mean
KDEsim2_SP <- density(sim2_SP)

plot(KDEsim2_SP,
     main = "KDE Simulation",
     col = "red")

# We try to generate random data based on a t-distribution (shifted by the log returns parameters).
#Some differences appears but the distribution seems to capture the bandwith

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

qqplot(rt_SP, qt(grid,df=2), 
       main="t-plot, df = 2 ",
       xlab="Data",
       ylab="t-quantiles")

lmfit <- lm(qt(c(.25,.75),df=2) ~ quantile(rt_SP,c(.25,.75)) )
abline(lmfit)

qqplot(rt_SP, qt(grid,df= df_t_SP), 
       main="t-plot, df estimated ",
       xlab="Data",
       ylab="t-quantiles")

lmfit <- lm(qt(c(.25,.75),df= df_t_SP) ~ quantile(rt_SP,c(.25,.75)) )
abline(lmfit)

#visually, it seems that the estimated dfs overperform the 2 dfs. However, we should still note that the theoretical quantiles
#not seem to capture the lower tail and expecially the upper tail of the distribution
#we proceed with a KDE to close the test
par(mfrow=c(1,3))

plot(KDE_SP,
     main = "KDE",
     col = "red")

sim1_est_SP <- rt(n_SP, df = 3) * sqrt(var(rt_SP) * 1/3) + mean(rt_SP) #simulate t dist with rt parameters rt*sqrt(var * df - 2/df ) + mean
KDEsim1_est_SP <- density(sim1_est_SP)
sim2_est_SP <- rt(n_SP, df = df_t_SP) * sqrt(var(rt_SP) * (df_t_SP - 2)/df_t_SP) + mean(rt_SP) #simulate t dist with rt parameters rt*sqrt(var * df - 2/df ) + mean
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
Box.test(rt_SP, lag = 10, type = c("Ljung-Box")) 
Box.test(rt_SP, lag = 20, type = c("Ljung-Box")) 
Box.test(rt_SP, lag = 30, type = c("Ljung-Box"))

cor(rt_SP,rt_SP)


# Using the Ljung-Box test, we test the null hypotheis of serial independence. As we can see from the 
# correlogram, there are different lags at which the serial correlation could be considered statistically
# significant. 
# We can as a matter of fact refuse the null hypothesis of serial independence
# 1% level: 
# lag 5 - lag  - lag 10 - lag 20 - lag 30 
# This could be a starting point for eventual modelling

##ARCH effect############################################################
sq_SP <- rt_SP^2
par(mfrow=c(1,2))
acf(sq_SP)
pacf(rt_SP)
pacf(sq_SP)

##ARMA model########################################################################
library(broom)
niter_order = 5
AIC = matrix(0, nrow = niter_order, ncol = niter_order)
bic = matrix(0, nrow = niter_order, ncol = niter_order)
LB = matrix(0, nrow = niter_order, ncol = niter_order)
AIC_nm = matrix(0, nrow = niter_order, ncol = niter_order)
bic_nm = matrix(0, nrow = niter_order, ncol = niter_order)
LB_nm = matrix(0, nrow = niter_order, ncol = niter_order)
for(i in c(0:4))
{
  for(j in c(0:4)){
    fitAR <- arima(rt_SP, order = c(i,0,j))
    AIC[i + 1,j + 1] = fitAR$aic
    bic[i + 1,j + 1] = BIC(fitAR)
    y = Box.test(fitAR$residuals, lag = 10,fitdf = (i + j), type = c("Ljung-Box"))
    LB[i + 1,j + 1] = y$p.value 
    fitAR_nm <- arima(rt_SP, order = c(i,0,j), include.mean = FALSE)
    AIC_nm[i + 1,j + 1] = fitAR_nm$aic
    bic_nm[i + 1,j + 1] = BIC(fitAR_nm)
    y_nm = Box.test(fitAR_nm$residuals, lag = 10, fitdf = (i + j), type = c("Ljung-Box"))
    LB_nm[i + 1,j + 1] = y_nm$p.value
  }
}

#finding minimum variance
vector_sigma <- as.vector(sigma)
sorted_sigma <- sort(vector_sigma)
min_sigma <- head(sorted_sigma, 5)
min_sigma_order <- list()
for(i in c(1:5))
{
  min_sigma_order[[i]] <- which(sigma == min_sigma[i], arr.ind = TRUE)
}
min_sigma_order

vector_sigma_nm <- as.vector(sigma_nm)
sorted_sigma_nm <- sort(vector_sigma_nm)
min_sigma_nm <- head(sorted_sigma_nm, 5)
min_sigma_order_nm <- list()
for(i in c(1:5))
{
  min_sigma_order_nm[[i]] <- which(sigma_nm == min_sigma_nm[i], arr.ind = TRUE)
}
min_sigma_order_nm


#finding minimum AIC
vector_AIC <- as.vector(AIC)
sorted_AIC <- sort(vector_AIC)
min_AIC <- head(sorted_AIC, 10)
min_AIC_order <- list()
for(i in c(1:10))
{
  min_AIC_order[[i]] <- which(AIC == min_AIC[i], arr.ind = TRUE)
}
min_AIC_order

vector_AIC_nm <- as.vector(AIC_nm)
sorted_AIC_nm <- sort(vector_AIC_nm)
min_AIC_nm <- head(sorted_AIC_nm, 10)
min_AIC_order_nm <- list()
for(i in c(1:10))
{
  min_AIC_order_nm[[i]] <- which(AIC_nm == min_AIC_nm[i], arr.ind = TRUE)
}
min_AIC_order_nm



##finding minimum BIC
vector_bic <- as.vector(bic)
sorted_bic <- sort(vector_bic)
min_bic <- head(sorted_bic, 5)
min_bic_order <- list()
for(i in c(1:5))
{
  min_bic_order[[i]] <- which(bic == min_bic[i], arr.ind = TRUE)
}
min_bic_order

vector_bic_nm <- as.vector(bic_nm)
sorted_bic_nm <- sort(vector_bic_nm)
min_bic_nm <- head(sorted_bic_nm, 5)
min_bic_order_nm <- list()
for(i in c(1:5))
{
  min_bic_order_nm[[i]] <- which(bic_nm == min_bic_nm[i], arr.ind = TRUE)
}
min_bic_order_nm


##finding maximum LB
vector_LB <- as.vector(LB)
sorted_LB <- sort(vector_LB)
max_LB <- tail(sorted_LB, 5)
max_LB_order <- list()
for(i in c(1:5))
{
  max_LB_order[[i]] <- which(LB == max_LB[i], arr.ind = TRUE)
}
max_LB_order

vector_LB_nm <- as.vector(LB_nm)
sorted_LB_nm <- sort(vector_LB_nm)
max_LB_nm <- tail(sorted_LB_nm, 5)
max_LB_order_nm <- list()
for(i in c(1:5))
{
  max_LB_order_nm[[i]] <- which(LB_nm == max_LB_nm[i], arr.ind = TRUE)
}
max_LB_order_nm

##ARMA-GARCH######################################
niter_order = 5
AIC = matrix(0, nrow = niter_order, ncol = niter_order)
BIC = matrix(0, nrow = niter_order, ncol = niter_order)
SIC = matrix(0, nrow = niter_order, ncol = niter_order)
HQIC = matrix(0, nrow = niter_order, ncol = niter_order)
LB = matrix(0, nrow = niter_order, ncol = niter_order)
for(i in 1:niter_order)
{
  for(j in 1:niter_order){
    uspec <- ugarchspec(variance.model = list(model = "sGARCH",
                                              garchOrder = c(i,j)), 
                        mean.model = list(armaOrder = c(3,3), include.mean = FALSE),
                        distribution.model = "std")
    fit <- ugarchfit(uspec, data = rt_SP)
    crit <- infocriteria(fit)
    AIC[i,j] = crit[1]
    BIC[i,j] = crit[2]
    SIC[i,j] = crit[3]
    HQIC[i,j] = crit[4]
    y = Box.test(fit@fit$z, lag = 10, type = c("Ljung-Box"))
    LB[i,j] = y$p.value 
  }
}



#finding minimum AIC
vector_AIC <- as.vector(AIC)
sorted_AIC <- sort(vector_AIC)
min_AIC <- head(sorted_AIC, 5)
min_AIC_order <- list()
for(i in c(1:5))
{
  min_AIC_order[[i]] <- which(AIC == min_AIC[i], arr.ind = TRUE)
}
min_AIC_order


##finding minimum BIC
vector_BIC <- as.vector(BIC)
sorted_BIC <- sort(vector_BIC)
min_BIC <- head(sorted_BIC, 5)
min_BIC_order <- list()
for(i in c(1:5))
{
  min_BIC_order[[i]] <- which(BIC == min_BIC[i], arr.ind = TRUE)
}
min_BIC_order

##finding minimum SIC
vector_SIC <- as.vector(SIC)
sorted_SIC <- sort(vector_SIC)
min_SIC <- head(sorted_SIC, 5)
min_SIC_order <- list()
for(i in c(1:5))
{
  min_SIC_order[[i]] <- which(SIC == min_SIC[i], arr.ind = TRUE)
}
min_SIC_order

##finding minimum HQIC
vector_HQIC <- as.vector(HQIC)
sorted_HQIC <- sort(vector_HQIC)
min_HQIC <- head(sorted_HQIC, 5)
min_HQIC_order <- list()
for(i in c(1:5))
{
  min_HQIC_order[[i]] <- which(HQIC == min_HQIC[i], arr.ind = TRUE)
}
min_HQIC_order

cor(rt_EUR, rt_SP)

##finding maximum LB
vector_LB <- as.vector(LB)
sorted_LB <- sort(vector_LB)
max_LB <- tail(sorted_LB, 5)
max_LB_order <- list()
for(i in c(1:5))
{
  max_LB_order[[i]] <- which(LB == max_LB[i], arr.ind = TRUE)
}
max_LB_order

library(rugarch)
SP_uspec <- ugarchspec(variance.model = list(model = "sGARCH",
                                             garchOrder = c(1,1)),
                       mean.model = list(armaOrder = c(3,3), include.mean = FALSE),
                       distribution.model = "std")
SP_fit <- ugarchfit(SP_uspec, data = rt_SP)
SP_fit@fit$matcoef
par(mfrow = c(1,1))
acf(SP_fit@fit$z)
Box.test(SP_fit@fit$z, lag = 30, type = c("Ljung-Box"))
plot(SP_fit)

infocriteria(SP_fit)

z_SP <- SP_fit@fit$z

rm(list = ls())
