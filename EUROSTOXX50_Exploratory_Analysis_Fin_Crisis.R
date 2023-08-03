# IMPORT DATA #################################################################################################
setwd("C:/Users/edoar/OneDrive/Desktop/Thesis/R_Thesis/Data/")
EUROSTOX_data = read.csv("EUROSTOXX50, 01-01-2007  to 31-12-2009_CLEAN.csv", header = TRUE)

## PRICES, TIME AND LENGHT ####################################################################################
Pt_EUR <- rev(EUROSTOX_data[,5]) 
Time_EUR <- rev(as.Date(as.character(EUROSTOX_data[,1]), "%m/%d/%Y"))
N_EUR <- length(Pt_EUR)
pt_EUR <- log(Pt_EUR)            # log prices for computation of log returns

### RETURNS ###################################################################################################
Rt_EUR <- diff(Pt_EUR, lag = 1)/Pt_EUR[-1]
n_EUR <- length(Rt_EUR)
rt_EUR <- diff(pt_EUR, lag = 1)

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
     main = "EUROSTOXX50 Log-Returns: 03/01/2020 - 30/12/2022",
     xlab = "Time",
     ylab = "Returns",
     col = "blue")

summary <- as.matrix(summary(rt_EUR))
cor(Rt_EUR,rt_EUR)                # correlation between Returns and Log Returns, expected good approximation for +-5%
plot(Rt_EUR,rt_EUR)
print(summary)
library(plotrix)
std.error(rt_EUR)


range_EUR <- max(rt_EUR) - min(rt_EUR)
range_EUR

IQR_EUR <- IQR(rt_EUR)
IQR_EUR
##### STATIONARITY ###########################################################################################
library(timeSeries)
ADF_EUR <- adf.test(rt_EUR)
ADF_EUR$p.value
ADF_EUR
##### MOMENTS #################################################################################################
library(moments)
mean_EUR <- mean(rt_EUR)
std.error(rt_EUR)
sd_EUR <- sd(rt_EUR)
sd_EUR
se_sd_EUR <- sd_EUR/sqrt(2*n_EUR- 2)
se_sd_EUR
skew_EUR <- skewness(rt_EUR) 
skew_EUR
kurt_EUR <- kurtosis(rt_EUR)
kurt_EUR
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
JB_EUR
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

grid <- (1:n_EUR)/(n_EUR+1)
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

#It seems that the best distribution is the df 3

####### KERNEL DENSITY ESTIMATION ##############################################################################
# We now try to perform a Kernel density estimation and comparing it with random data generated according to a
# t-distribution with 3 degrees of freedom

par(mfrow=c(1,2))
KDE_EUR <- density(rt_EUR) 

plot(KDE_EUR,
     main = "KDE",
     col = "red")

sim1_EUR <- rt(n_EUR, df = 3) * sqrt(var(rt_EUR) * 1/3) + mean(rt_EUR) #simulate t dist with rt parameters rt*sqrt(var * df - 2/df ) + mean
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
fit_t_EUR <- fitdist(rt_EUR,"t.scaled", start = list(df = 3, mean = mean(rt_EUR), sd = sd(rt_EUR)))
summary(fit_t_EUR)
df_t_EUR = fit_t_EUR$estimate[1]

#Let's try to test this estimate and compare it with the previous chosen 3 dfs
par(mfrow=c(1,2))

qqplot(rt_EUR, qt(grid,df=3), 
       main="t-plot, df = 3 ",
       xlab="Data",
       ylab="t-quantiles")

lmfit <- lm(qt(c(.25,.75),df=3) ~ quantile(rt_EUR,c(.25,.75)) )
abline(lmfit)

qqplot(rt_EUR, qt(grid,df= df_t_EUR), 
       main="t-plot, df estimated ",
       xlab="Data",
       ylab="t-quantiles")

lmfit <- lm(qt(c(.25,.75),df= df_t_EUR) ~ quantile(rt_EUR,c(.25,.75)) )
abline(lmfit)

#visually, it seems that the estimated dfs underperforms the 3 dfs. However, we should still note that the theoretical quantiles
#seem to better capture the tails of the distribution. These are usually our quantiles of interest
#we proceed with a KDE to close the test


sim1_est_EUR <- rt(n_EUR, df = 3) * sqrt(var(rt_EUR) * 1/3) + mean(rt_EUR) #simulate t dist with rt parameters rt*sqrt(var * df - 2/df ) + mean
KDEsim1_est_EUR<- density(sim1_est_EUR)
sim2_est_EUR <- rt(n_EUR, df = df_t_EUR) * sqrt(var(rt_EUR) * (df_t_EUR - 2)/df_t_EUR) + mean(rt_EUR) #simulate t dist with rt parameters rt*sqrt(var * df - 2/df ) + mean
KDEsim2_est_EUR <- density(sim2_est_EUR)

par(mfrow=c(1,3))

plot(KDE_EUR,
     main = "KDE",
     col = "black")

plot(KDEsim1_est_EUR,
     main = "KDE 3 df",
     col = "red")
plot(KDEsim2_est_EUR,
     main = "KDE df estimated",
     col = "blue")



######## TEST FOR SERIAL CORRELATION #########################################################################
par(mfrow=c(1,2))

acf(rt_EUR,
    main = "ACF LogReturns EUR")
pacf(rt_EUR,
     main = "PACF LogReturns EUR")       
Box.test(rt_EUR, lag = 5, type = c("Ljung-Box")) 
Box.test(rt_EUR, lag = 10, type = c("Ljung-Box")) 
Box.test(rt_EUR, lag = 20, type = c("Ljung-Box")) 
Box.test(rt_EUR, lag = 30, type = c("Ljung-Box"))  

# Using the Ljung-Box test, we test the null hypotheis of serial independence. As we can see from the 
# correlogram, there are different lags at which the serial correlation could be considered statistically
# significant. 
# We can as a matter of fact refuse the null hypothesis of serial independence
# 1% level: 
# lag 5,10
# 5% level:
# lag 20,30
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
AIC = matrix(0, nrow = niter_order, ncol = niter_order)
bic = matrix(0, nrow = niter_order, ncol = niter_order)
LB = matrix(0, nrow = niter_order, ncol = niter_order)
AIC_nm = matrix(0, nrow = niter_order, ncol = niter_order)
bic_nm = matrix(0, nrow = niter_order, ncol = niter_order)
LB_nm = matrix(0, nrow = niter_order, ncol = niter_order)
for(i in c(0:4))
{
  for(j in c(0:4)){
    fitAR <- arima(rt_EUR, order = c(i,0,j))
    AIC[i + 1,j + 1] = fitAR$aic
    bic[i + 1,j + 1] = BIC(fitAR)
    y = Box.test(fitAR$residuals, lag = 10, fitdf = (i + j), type = c("Ljung-Box"))
    LB[i + 1,j + 1] = y$p.value
    fitAR_nm <- arima(rt_EUR, order = c(i,0,j), include.mean = FALSE)
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

par(mfrow = c(1,1))


library(forecast)
auto.arima(rt_EUR)
#best (2,3) zero mean

#GARCH model
library(rugarch)
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
                        mean.model = list(armaOrder = c(2,3), include.mean = FALSE),
                        distribution.model = "std")
    fit <- ugarchfit(uspec, data = rt_EUR)
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
                        mean.model =list(armaOrder = c(2,3), include.mean = FALSE),
                        distribution.model = "std")
EUR_fit <- ugarchfit(EUR_uspec, data = rt_EUR)
EUR_fit@fit$matcoef
par(mfrow = c(1,1))
acf(EUR_fit@fit$z)
Box.test(EUR_fit@fit$z, lag = 30, type = c("Ljung-Box"))
plot(EUR_fit)
infocriteria(EUR_fit)

#ARMA(2,3)zm - GARCH(1,1) seems to be the best fit. second best ARMA(3,3) - GARCH(2,1)

EUR_fit@fit$matcoef  

z_EUR <- EUR_fit@fit$z
