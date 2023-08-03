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

SP500_data = read.csv("S&P500, 24-12-2012 to 21-12-2022_CLEAN.csv", header = TRUE)

## PRICES, TIME AND LENGHT ####################################################################################
Pt_SP <- SP500_data[,6] 
Time_SP <- as.Date(as.character(SP500_data[,1]), "%m/%d/%Y")
N_SP <- length(Pt_SP)
pt_SP <- log(Pt_SP)            # log prices for computation of log returns

### RETURNS ###################################################################################################
Rt_SP <- diff(Pt_SP, lag = 1)/Pt_SP[-1]
n_SP <- length(Rt_SP)

#Summary statistics on daily return##########################################


#mean
SP_mean = mean(Rt_SP)
SP_se_mean = sd(Rt_SP)/sqrt(n_SP)
EUR_mean = mean(Rt_EUR)
EUR_se_mean = sd(Rt_EUR)/sqrt(n_EUR)
#negative mean for both approximately equal to 0 as expected

#standard deviation
SP_sd = sd(Rt_SP)
SP_se_sd = SP_sd / sqrt(2 * (n_SP - 1))
EUR_sd = sd(Rt_EUR)
EUR_se_sd = EUR_sd / sqrt(2 * (n_EUR - 1))
#sligthly higher average volatility in European markets

#skewness
library(moments)
SP_skew = skewness(Rt_SP)
SP_se_skew = sqrt(6 * n_SP * (n_SP - 1) / ((n_SP - 2) * (n_SP + 1) * (n_SP + 3)))
EUR_skew = skewness(Rt_EUR)
EUR_se_skew = sqrt(6 * n_EUR * (n_EUR - 1) / ((n_EUR - 2) * (n_EUR + 1) * (n_EUR + 3)))
#as expected, both positively skewed

#kurtosis
SP_kurt = kurtosis(Rt_SP)
SP_se_kurt = sqrt((24 * n_SP * (n_SP - 2) * (n_SP - 3)) / ((n_SP + 1) * (n_SP + 3) * (n_SP + 5) * (n_SP + 7))) 
EUR_kurt = kurtosis(Rt_EUR)
EUR_se_kurt = sqrt((24 * n_EUR * (n_EUR - 2) * (n_EUR - 3)) / ((n_EUR + 1) * (n_EUR + 3) * (n_EUR + 5) * (n_EUR + 7))) 
#again as expected, fatter tails than normal

#Normality
SP_JB = jarque.test(Rt_SP)
SP_JB
EUR_JB = jarque.test(Rt_EUR)
EUR_JB
#reject the null hypothesis of normality, as confirmed previously by skewness and kurtosis

#Autocorrelation
##Ljung_Box S&P500
par(mfrow=c(1,2))
acf(Rt_SP)
pacf(Rt_SP) 
LB_SP_1 = Box.test(Rt_SP, lag = 1, type = c("Ljung-Box")) 
LB_SP_5 = Box.test(Rt_SP, lag = 5, type = c("Ljung-Box")) 
LB_SP_7 = Box.test(Rt_SP, lag = 7, type = c("Ljung-Box")) 
LB_SP_8 = Box.test(Rt_SP, lag = 8, type = c("Ljung-Box"))
LB_SP_22 = Box.test(Rt_SP, lag = 22, type = c("Ljung-Box")) 
LB_SP_24 = Box.test(Rt_SP, lag = 24, type = c("Ljung-Box")) 
LB_SP_30 = Box.test(Rt_SP, lag = 30, type = c("Ljung-Box"))

###Ljung-Box statistics
LB_SP_1_st = LB_SP_1$statistic
LB_SP_5_st = LB_SP_5$statistic
LB_SP_7_st = LB_SP_7$statistic
LB_SP_8_st = LB_SP_8$statistic
LB_SP_22_st = LB_SP_22$statistic
LB_SP_24_st = LB_SP_24$statistic
LB_SP_30_st = LB_SP_30$statistic

###Ljung-Box p-value
LB_SP_1_p = LB_SP_1$p.value
LB_SP_5_p = LB_SP_5$p.value
LB_SP_7_p = LB_SP_7$p.value
LB_SP_8_p = LB_SP_8$p.value
LB_SP_22_p = LB_SP_22$p.value
LB_SP_24_p = LB_SP_24$p.value
LB_SP_30_p = LB_SP_30$p.value

##Ljung_Box EUROSTOXX50
par(mfrow=c(1,2))
acf(Rt_EUR)
pacf(Rt_EUR) 
LB_EUR_7 = Box.test(Rt_EUR, lag = 7, type = c("Ljung-Box")) 
LB_EUR_13 = Box.test(Rt_EUR, lag = 13, type = c("Ljung-Box")) 
LB_EUR_15 = Box.test(Rt_EUR, lag = 15, type = c("Ljung-Box")) 
LB_EUR_23 = Box.test(Rt_EUR, lag = 23, type = c("Ljung-Box"))
LB_EUR_28 = Box.test(Rt_EUR, lag = 28, type = c("Ljung-Box"))
LB_EUR_32 = Box.test(Rt_EUR, lag = 32, type = c("Ljung-Box")) 

###Ljung-Box statistics
LB_EUR_7_st = LB_EUR_7$statistic
LB_EUR_13_st = LB_EUR_13$statistic
LB_EUR_15_st = LB_EUR_15$statistic
LB_EUR_23_st = LB_EUR_23$statistic
LB_EUR_28_st = LB_EUR_28$statistic
LB_EUR_32_st = LB_EUR_32$statistic

###Ljung-Box p-value
LB_EUR_7_st = LB_EUR_7$p.value
LB_EUR_13_st = LB_EUR_13$p.value
LB_EUR_15_st = LB_EUR_15$p.value
LB_EUR_23_st = LB_EUR_23$p.value
LB_EUR_28_st = LB_EUR_28$p.value
LB_EUR_32_st = LB_EUR_32$p.value

#Correlation matrix
Rt_matrix <- cbind(Rt_SP, Rt_EUR)
cor_matrix <- cor(Rt_matrix)

#Transformation by CDF
par(mfrow=c(1,1))
library(spgs)
library(ggplot2)
norm_SP <- pnorm(Rt_SP, mean = SP_mean, sd = SP_sd)
ks.test(norm_SP, "punif", min(norm_SP), max(norm_SP))
mean(norm_SP)
sd(norm_SP)

seq(0.05,1,length(norm_SP))
plot(1:n_SP/(n_SP+1),sort(norm_SP))

plot(ecdf(norm_SP))
curve(punif(x, min(norm_SP), max(norm_SP)), add=TRUE, col="red")

hist(norm_SP, breaks = 20, freq = FALSE)
lines(c(0, 1), c(1.25, 1.25), col = "red", lwd = 2)
lines(c(0, 1), c(1,1), col = "blue", lwd = 2)
lines(c(0, 1), c(0.75, 0.75), col = "red", lwd = 2)

norm_SP <- pnorm(norm_SP, mean = mean(norm_SP), sd = sd(norm_SP))
ks.test(norm_SP, "punif", min(norm_SP), max(norm_SP))
chisq.unif.test(norm_SP)

seq(0.05,1,length(norm_SP))
plot(1:n_SP/(n_SP+1),sort(norm_SP))

plot(ecdf(norm_SP))
curve(punif(x, min(norm_SP), max(norm_SP)), add=TRUE, col="red")

hist(norm_SP, breaks = 20, freq = FALSE)
lines(c(0, 1), c(1.25, 1.25), col = "red", lwd = 2)
lines(c(0, 1), c(1,1), col = "blue", lwd = 2)
lines(c(0, 1), c(0.75, 0.75), col = "red", lwd = 2)

mean(norm_SP)
sd(norm_SP)

standardized_t_SP <- (Rt_SP - SP_mean)/SP_sd
t_SP <- pt(standardized_t_SP, df = 3, lower.tail = TRUE)
hist(t_SP, breaks = 20)
t_SP <- pt(t_SP, df = 3, lower.tail = TRUE)
hist(t_SP, breaks = 20)


norm_EUR <- pnorm(Rt_EUR, mean = EUR_mean, sd = EUR_sd)
hist(norm_EUR, breaks = 20)
chisq.unif.test(norm_EUR)

seq(0.05,1,length(norm_EUR))
plot(1:n_EUR/(n_EUR+1),sort(norm_EUR))

plot(ecdf(norm_EUR))
curve(punif(x, min(norm_EUR), max(norm_EUR)), add=TRUE, col="red")

norm_EUR <- pnorm(norm_EUR, mean = mean(norm_EUR), sd = sd(norm_EUR))
hist(norm_EUR, breaks = 20)
chisq.unif.test(norm_EUR)

seq(0.05,1,length(norm_SP))
plot(1:n_SP/(n_SP+1),sort(norm_SP))

plot(ecdf(norm_SP))
curve(punif(x, min(norm_SP), max(norm_SP)), add=TRUE, col="red")

standardized_t_EUR <- (Rt_EUR - EUR_mean)/EUR_sd
t_EUR <- pt(standardized_t_EUR, df = 4, lower.tail = TRUE)
hist(t_EUR, breaks = 20)





rm(list = ls())
