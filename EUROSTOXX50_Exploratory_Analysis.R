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

#outliers###
outliers_EUR <- sum(Rt_EUR > 0.15 | Rt_EUR < -0.15)
plot(rt_EUR)

z_scores <- (Rt_EUR - mean(Rt_EUR)) / sd(Rt_EUR)

# Step 3: Set threshold (e.g., 2 or 3)
threshold <- 5

# Step 4: Identify outliers
outliers <- Rt_EUR[abs(z_scores) > threshold]

# Print the outliers
print(outliers)

#### LOG RETURNS ##############################################################################################
rt_EUR <- diff(pt_EUR, lag = 1)

plot(Time_EUR[-1],rt_EUR,
     type = "l",
     main = "EUROSTOXX50 Log-Returns: 25/12/2012 - 21/12/2022",
     xlab = "Time",
     ylab = "Returns",
     col = "blue")

summary(rt_EUR)
cor(Rt_EUR,rt_EUR)                # correlation between Returns and Log Returns, expected good approximation for +-5%
plot(Rt_EUR,rt_EUR)

library(plotrix)
std.error(rt_EUR)


range_EUR <- max(rt_EUR) - min(rt_EUR)
range_EUR

IQR_EUR <- IQR(rt_EUR)
IQR_EUR

##### MOMENTS #################################################################################################
library(moments)
sd_EUR <- sd(rt_EUR)
se_sd_EUR <- sd_EUR/sqrt(2*n_EUR- 2)
skew_EUR <- skewness(rt_EUR)         
kurt_EUR <- kurtosis(rt_EUR)
V_skew = 6 * n_EUR * (n_EUR - 1) / ((n_EUR - 2) * (n_EUR + 1) * (n_EUR + 3))
sqrt(V_skew)
kurt_EUR <- kurtosis(rt_EUR)
V_kur = 4 * (n_EUR^2 - 1) * V_skew / ((n_EUR - 3) * (n_EUR + 5))
sqrt(V_kur)

# both skewness and kurtosis are signs of non-normality since skewness != 0 and kurtosis > 3, both significantly

##### TEST FOR MEAN = 0 #######################################################################################
t.test(rt_EUR)
# the high p-value indicates the impossibility to refuse the null hypothesis of true mean = 0

###### TESTING FOR NORMALITY ##################################################################################
library(tseries)
JB_EUR <- jarque.bera.test(rt_EUR)
# the extremely small p-value signifies that we can refuse the null hypothesis of normality

qqnorm(rt_EUR,
       datax = TRUE,
       main = "Normal qq plot for EUR LogReturns")
qqline(rt_EUR,
       datax = TRUE)
# as we could expect, while the sample and theoretical quantiles are approximately equal in the central range
# (around 0.00) the tails deviate from the theoretical quantiles. This indicates fatter tails than a 
# normal distributions as is often the case with financial data. We now proceed to test the hypothesis of a 
# t-distribution at different levels of degrees of freedom

####### TESTING FOR T-DISTRIBUTION ############################################################################
grid <- (1:n)/(n+1)
par(mfrow=c(1,3))

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
pacf(rt_EUR)       
Box.test(rt_EUR, lag = 5, type = c("Ljung-Box")) 
Box.test(rt_EUR, lag = 10, type = c("Ljung-Box")) 
Box.test(rt_EUR, lag = 20, type = c("Ljung-Box")) 
b <- Box.test(rt_EUR, lag = 30, type = c("Ljung-Box"))  

# Using the Ljung-Box test, we test the null hypotheis of serial independence. As we can see from the 
# correlogram, there are different lags at which the serial correlation could be considered statistically
# significant. 
# We can as a matter of fact refuse the null hypothesis of serial independence
# 5% level:
# lag 7 - lag 13
# 1% level: 
# lag 15 - lag 23
# This could be a starting point for eventual modelling

##ARCH effect#########################################################################################################
sq_EUR <- rt_EUR^2
par(mfrow = c(1,2))
acf(sq_EUR)
pacf(rt_EUR)
pacf(sq_EUR)

##ARMA model################################################################################################################
library(broom)
niter_order = 5
sigma = matrix(0, nrow = niter_order, ncol = niter_order)
AIC = matrix(0, nrow = niter_order, ncol = niter_order)
bic = matrix(0, nrow = niter_order, ncol = niter_order)
LB = matrix(0, nrow = niter_order, ncol = niter_order)
for(i in c(0:4))
{
  for(j in c(0:4)){
    fitAR <- arima(rt_EUR, order = c(i,0,j))
    sigma[i + 1,j + 1] = fitAR$sigma2
    AIC[i + 1,j + 1] = fitAR$aic
    bic[i + 1,j + 1] = BIC(fitAR)
    y = Box.test(fitAR$residuals, lag = 10, type = c("Ljung-Box"))
    LB[i + 1,j + 1] = y$p.value 
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

#ARMA minimum AIC
EUR_ARMA_33 <- arima(rt_EUR, order = c(3,0,3), include.mean = FALSE)
acf(EUR_ARMA_33$residuals)
EUR_ARMA_33
coeftest(EUR_ARMA_33) ##perfect
Box.test(EUR_ARMA_33$residuals,  lag = 30, type = c("Ljung-Box")) ##definitely more independent than autoarima (0.12 vs 0.006)
BIC(EUR_ARMA_33) ##BIC higher than autoarima
EUR_ARMA_33$aic  ##AIC lower than autoarima


EUR_ARMA_32 <- arima(rt_EUR, order = c(3,0,2))
acf(EUR_ARMA_32$residuals)
EUR_ARMA_32
coeftest(EUR_ARMA_32) ##residuals strongly not correlated but AR(3) strongly insignificant (0.91), intercept too
LB[4,3]

EUR_ARMA_32_nm <- arima(rt_EUR, order = c(3,0,2), include.mean = FALSE)
acf(EUR_ARMA_32_nm$residuals)
EUR_ARMA_32_nm
coeftest(EUR_ARMA_32_nm) ##residuals strongly not correlated but AR(3) strongly insignificant (0.91), intercept too
Box.test(EUR_ARMA_32_nm$residuals,  lag = 30, type = c("Ljung-Box"))

EUR_ARMA_23 <- arima(rt_EUR, order = c(2,0,3))
acf(EUR_ARMA_23$residuals)
EUR_ARMA_23
coeftest(EUR_ARMA_23) ##residuals strongly not correlated but MA(3) strongly insignificant (0.87), intercept too
LB[3,4]

EUR_ARMA_00 <- arima(rt_EUR, order = c(0,0,0))
acf(EUR_ARMA_00$residuals)
EUR_ARMA_00
coeftest(EUR_ARMA_00) ##residuals mildly correlated (5%), intercept strongly insignificant
LB[1,1]

EUR_ARMA_10 <- arima(rt_EUR, order = c(1,0,0))
acf(EUR_ARMA_10$residuals)
EUR_ARMA_10
coeftest(EUR_ARMA_10) ##residuals mildly correlated, intercept and AR(1) strongly insignificant
LB[2,1]

EUR_ARMA_01 <- arima(rt_EUR, order = c(0,0,1))
acf(EUR_ARMA_01$residuals)
EUR_ARMA_01
coeftest(EUR_ARMA_01) ##residuals mildly correlated, intercept and MA(1) strongly insignificant
LB[1,2]

EUR_ARMA_42 <- arima(rt_EUR, order = c(4,0,2))
acf(EUR_ARMA_42$residuals)
EUR_ARMA_42
coeftest(EUR_ARMA_42) ##residuals mildly correlated (10%), intercept, AR(3) and AR(4) strongly insignificant
LB[5,3]

EUR_ARMA_43 <- arima(rt_EUR, order = c(4,0,3))
acf(EUR_ARMA_43$residuals)
EUR_ARMA_43
coeftest(EUR_ARMA_43) ##residuals not correlated, intercept, only 2/8 parameters are significant
LB[5,4]

EUR_ARMA_34 <- arima(rt_EUR, order = c(3,0,4))
acf(EUR_ARMA_34$residuals)
EUR_ARMA_34
coeftest(EUR_ARMA_34) ##residuals not correlated, intercept, only 2/8 parameters are significant
LB[4,5]

EUR_ARMA_02 <- arima(rt_EUR, order = c(0,0,2))
acf(EUR_ARMA_02$residuals)
EUR_ARMA_02
coeftest(EUR_ARMA_02) ##residuals mildly correlated (10%), all parameters are insignificant
LB[1,3]

EUR_ARMA_20 <- arima(rt_EUR, order = c(2,0,0))
acf(EUR_ARMA_20$residuals)
EUR_ARMA_20
coeftest(EUR_ARMA_20) ##residuals mildly correlated (10%), all parameters are insignificant
LB[3,1]

EUR_ARMA_00_nm <- arima(rt_EUR, order = c(0,0,0), include.mean = FALSE)
acf(EUR_ARMA_00_nm$residuals)
EUR_ARMA_00_nm
coeftest(EUR_ARMA_00_nm) ##residuals mildly correlated (10%), all parameters are insignificant
Box.test(EUR_ARMA_00_nm$residuals,  lag = 30, type = c("Ljung-Box"))
BIC(EUR_ARMA_00_nm)

auto.arima(rt_EUR)

acf(EUR_ARMA_33$residuals^2)

#GARCH model
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
                        mean.model = EUR_ARMA_33)
    fit <- ugarchfit(uspec, data = EUR_ARMA_33$residuals)
    crit <- infocriteria(fit)
    AIC[i,j] = crit[1]
    BIC[i,j] = crit[2]
    SIC[i,j] = crit[3]
    HQIC[i,j] = crit[4]
    y = Box.test(fit@fit$residuals, lag = 10, type = c("Ljung-Box"))
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


EUR_uspec <- ugarchspec(variance.model = list(model = "sGARCH",
                                          garchOrder = c(1,1)), 
                    mean.model = EUR_ARMA_33)
EUR_fit <- ugarchfit(EUR_uspec, data = EUR_ARMA_33$residuals)
EUR_fit@fit$matcoef
acf(EUR_fit@fit$residuals)
Box.test(EUR_fit@fit$residuals, lag = 10, type = c("Ljung-Box"))
plot(EUR_fit)
infocriteria(EUR_fit)

#it appears that the best fit is a (3,3) or (1,1). Both are present in the top 5 across all criterion.
#(1,1) since all parameters are significant, while critertions give different conclusions, this one is the most probable
#across all regarding independence and (3,3) has two parameters insignificant (one at 10% and one even higher)

rm(list = ls())

