#correlation plot
par(mfrow = c(1,2))
cor(rt_EUR, rt_SP)
plot(rt_SP, rt_EUR)
cor(z_EUR, z_SP)
plot(z_SP, z_EUR)

#no ties are present in the residuals
duplicates <- duplicated(z_EUR)
print(z_EUR[duplicates])
duplicates <- duplicated(z_SP)
print(z_SP[duplicates])


#denisty estimation################################################################################################################
#normal
library(spgs)
par(mfrow = c(1,1))
norm_z_SP <- pnorm(z_SP, mean = mean(z_SP), sd = sd(z_SP))
hist(norm_z_SP,
     breaks = 20,
     col = "blue")
chisq.unif.test(norm_z_SP)

norm_z_EUR <- pnorm(z_EUR, mean = mean(z_EUR), sd = sd(z_EUR))
hist(norm_z_EUR,
     breaks = 20,
     col = "blue")
chisq.unif.test(norm_z_EUR)


#t-student with 3 and 4 df
t3_z_SP <- pt(z_SP, df = 3, lower.tail = TRUE)
hist(t3_z_SP, 
     breaks = 20,
     col = "blue")

chisq.unif.test(t3_z_SP)

t4_z_SP <- pt(z_SP, df = 4, lower.tail = TRUE)
hist(t4_z_SP, 
     breaks = 20,
     col = "blue")

chisq.unif.test(t3_z_SP)

t4_z_EUR <- pt(z_EUR, df = 4, lower.tail = TRUE)
hist(t4_z_EUR, 
     breaks = 20,
     col = "blue")

chisq.unif.test(t4_z_EUR)

#standardized t
par(mfrow = c(1,1))
library(fGarch)
stdFit(z_EUR)
EUR_fit@fit$matcoef
std_z_EUR <- pstd(z_EUR, mean = mean(z_EUR), sd = sd(z_EUR) , nu = 4.024836)  
hist(std_z_EUR, 
     main = "F(z) EUR Estimated Standardized t",
     breaks = 10,
     col = "blue")

chisq.unif.test(std_z_EUR) #appears to be the correct distribution

stdFit(z_SP)
SP_fit@fit$matcoef
std_z_SP <- pstd(z_SP, mean = mean(z_SP), sd = sd(z_SP) , nu = 9.387547)  
hist(std_z_SP,
     main = "F(z) SP Estimated Standardized t",
     breaks = 10,
     col = "blue")

chisq.unif.test(std_z_SP)


#EMP#
##since none of the parametric distributions seems to work, we will try to use the empirical distribution function
EMP_SP <- ecdf(z_SP)
u_SP <- EMP_SP(z_SP)
hist(u_SP, 
     breaks = 20,
     col = "blue")
chisq.unif.test(u_SP)


EMP_EUR <- ecdf(z_EUR)
u_EUR <- EMP_EUR(z_EUR)
hist(u_EUR, 
     breaks = 20,
     col = "blue")
chisq.unif.test(u_EUR)

#Copula-GARCH model################################################################################################################
library(npcp)
library(copula)
library(lattice)
Z <- data.frame(SP = z_SP, EUR = z_EUR)
pobs(Z)
U_EMP <- as.matrix(data.frame(SP = u_SP, EUR = u_EUR))
U <- pobs(Z)


par(mfrow = c(1,3))
plot(rt_SP, rt_EUR)
plot(z_SP,z_EUR)
plot(U)


par(mfrow = c(1,1))
plot(U)

set.seed(2013)
cpCop <- (cpCopula(Z, b = 1, method = "nonseq"))
cpCop 
plot(U)
?cpCopula  #time invariant copula at the p-value: 0.1943

#Methods of moments based on tau Normal Copula ###################################################################
tau.n <- cor(Z[,1], Z[,2], method = "kendall") ## The sample version of Kendall’s tau

fitNorm_t <- (fitCopula(normalCopula(dim = 2, dispstr = "un"), data = U, method = "itau")) ## The corresponding copula parameter estimate
summary(fitNorm_t)

library(sfsmisc)
## Extract parameter estimates
len_nt <- length(coef(fitNorm_t))
p_nt <- coef(fitNorm_t)

## Define the H_0 copula
cop.n <- normalCopula(param = p_nt, dim = ncol(U))
## Build the array of pairwise H_0-transformed data columns
cu.u.n <- pairwiseCcop(U, cop.n)
## Compute pairwise matrix (d by d) of p-values and corresponding colors
set.seed(1389) # for reproducibility
pw.indep.n <- pairwiseIndepTest(cu.u.n, N = 742, verbose = FALSE)
p.val.n <- pviTest(pw.indep.n) # matrix of p-values
## Pairwise Rosenblatt plot (test for the fitted normal copula)
title <- list("P. R. transformed pseudo-observations rho Normal")
cols <- gray(seq(0.1, 1, length.out = 6))
pairsRosenblatt(cu.u.n, pvalueMat = p.val.n, method = "QQchisq", pch=".",
                xaxt = "n", yaxt = "n",
                colList = pairsColList(p.val.n, bucketCols = NULL,
                                       BWcutoff = 170, cols = c(1,1)),
                main.centered = TRUE, main = title)


#Methods of moments based on rho Normal Copula##################################
fitNorm_r <- (fitCopula(normalCopula(dim = 2, dispstr = "un"), data = U, method = "irho")) ## The corresponding copula parameter estimate
summary(fitNorm_r)

library(sfsmisc)
## Extract parameter estimates
len_nr <- length(coef(fitNorm_r))
p_nr <- coef(fitNorm_r)

## Define the H_0 copula
cop.n <- normalCopula(param = p_nr, dim = ncol(U))
## Build the array of pairwise H_0-transformed data columns
cu.u.n <- pairwiseCcop(U, cop.n)
## Compute pairwise matrix (d by d) of p-values and corresponding colors
set.seed(1389) # for reproducibility
pw.indep.n <- pairwiseIndepTest(cu.u.n, N = 742, verbose = FALSE)
p.val.n <- pviTest(pw.indep.n) # matrix of p-values
## Pairwise Rosenblatt plot (test for the fitted normal copula)
title <- list("P. R. transformed pseudo-observations rho Normal")
cols <- gray(seq(0.1, 1, length.out = 6))
pairsRosenblatt(cu.u.n, pvalueMat = p.val.n, method = "QQchisq", pch=".",
                xaxt = "n", yaxt = "n",
                colList = pairsColList(p.val.n, bucketCols = NULL,
                                       BWcutoff = 170, cols = c(1,1)),
                main.centered = TRUE, main = title)


#T-copula MOME_MPL#################################################################################
fitT_t <- (fitCopula(tCopula(dim = 2, dispstr = "un"), data = U, method = "itau.mpl")) ## The corresponding copula parameter estimate
summary(fitT_t)
fitT_t@copula@parameters

## Extract parameter estimates
len <- length(coef(fitT_t))
p_tt <- coef(fitT_t)[1]
df_tt <- coef(fitT_t)[2]

## Define the H_0 copula
cop.t <- tCopula(param = p_tt, df = df_tt, dim = ncol(U))
## Build the array of pairwise H_0-transformed data columns
cu.u.t <- pairwiseCcop(U, cop.t)
## Compute pairwise matrix (d by d) of p-values and corresponding colors
set.seed(1389) # for reproducibility
pw.indep.t <- pairwiseIndepTest(cu.u.t, N = 742, verbose = FALSE)
p.val.t <- pviTest(pw.indep.t) # matrix of p-values
## Pairwise Rosenblatt plot (test for the fitted normal copula)
title <- list("P.R. transformed pseudo-observations tau.mpl t")
cols <- gray(seq(0.1, 1, length.out = 6))
pairsRosenblatt(cu.u.t, pvalueMat = p.val.t, method = "QQchisq", pch=".",
                xaxt = "n", yaxt = "n",
                colList = pairsColList(p.val.t, bucketCols = NULL,
                                       BWcutoff = 0),
                main.centered = TRUE, main = title)



#Frank Copula MOM tau#########################################################################################
fitFr_t <- (fitCopula(frankCopula(dim = 2), data = U, method = "itau")) ## The corresponding copula parameter estimate
summary(fitFr_t)
fitFr_t@copula@parameters


## Extract parameter estimates
len_frt <- length(coef(fitFr_t))
p_frt <- fitFr_t@copula@parameters[1]



#Frank Copula MOM rho#########################################################################################
fitFr_r <- (fitCopula(frankCopula(dim = 2), data = U, method = "irho")) ## The corresponding copula parameter estimate
summary(fitFr_r)
fitFr_r@copula@parameters

## Extract parameter estimates
len_frr<- length(coef(fitFr_r))
p_frr <- fitFr_r@copula@parameters[1]



#Gumebl-copula MOME tau #################################################################################
fitG_t <- (fitCopula(gumbelCopula(dim = 2), data = U, method = "itau")) ## The corresponding copula parameter estimate
summary(fitG_t)
fitG_t@copula@parameters

## Extract parameter estimates
len_Gt <- length(coef(fitG_t))
p_Gt <- coef(fitG_t)[1]



#Gumbel-copula MOME rho #################################################################################
fitG_r <- (fitCopula(gumbelCopula(dim = 2), data = U, method = "irho")) ## The corresponding copula parameter estimate
summary(fitG_r)
fitG_r@copula@parameters

## Extract parameter estimates
len_Cr <- length(coef(fitC_r))
p_Cr <- coef(fitC_r)[1]



#Maximum pseudo likelihood estimation and Rosenblatt GOF###########################################################
#Normal copula MPL#############################################################################
fitNorm_mpl <- (fitCopula(normalCopula(dim = 2, dispstr = "un"), data = U, method = "mpl")) ## The corresponding copula parameter estimate
summary(fitNorm_mpl)



library(sfsmisc)
## Extract parameter estimates
len_nm <- length(coef(fitNorm_mpl))
p_nm <- coef(fitNorm_mpl)

## Define the H_0 copula
cop.n <- normalCopula(param = p_nm, dim = ncol(U))
## Build the array of pairwise H_0-transformed data columns
cu.u.n <- pairwiseCcop(U, cop.n)
## Compute pairwise matrix (d by d) of p-values and corresponding colors
set.seed(1389) # for reproducibility
pw.indep.n <- pairwiseIndepTest(cu.u.n, N = 742, verbose = FALSE)
p.val.n <- pviTest(pw.indep.n) # matrix of p-values
## Pairwise Rosenblatt plot (test for the fitted normal copula)
title <- list("P.R. transformed pseudo-observations MPL Normal")
cols <- gray(seq(0.1, 1, length.out = 6))
pairsRosenblatt(cu.u.n, pvalueMat = p.val.n, method = "QQchisq", pch=".",
                xaxt = "n", yaxt = "n",
                colList = pairsColList(p.val.n, bucketCols = NULL,
                                       BWcutoff = 0),
                main.centered = TRUE, main = title)

?gof

#T-Copula MPL##################################################################################################
fitT_m <- (fitCopula(tCopula(dim = 2, dispstr = "un"), data = U, method = "mpl")) ## The corresponding copula parameter estimate
summary(fitT_m)
fitT_m@copula@parameters


## Extract parameter estimates
len <- length(coef(fitT_m))
p_tm <- fitT_m@copula@parameters[1]
df_tm <- fitT_m@copula@parameters[2]

## Define the H_0 copula
cop.t <- tCopula(param = p_tm, df = df_tm, dim = ncol(U))
## Build the array of pairwise H_0-transformed data columns
cu.u.t <- pairwiseCcop(U, cop.t)
## Compute pairwise matrix (d by d) of p-values and corresponding colors
set.seed(1389) # for reproducibility
pw.indep.t <- pairwiseIndepTest(cu.u.t, N = 742, verbose = FALSE)
p.val.t <- pviTest(pw.indep.t) # matrix of p-values
## Pairwise Rosenblatt plot (test for the fitted normal copula)
title <- list("P.R. transformed pseudo-observations MPL t")
cols <- gray(seq(0.1, 1, length.out = 6))
pairsRosenblatt(cu.u.t, pvalueMat = p.val.t, method = "QQchisq", pch=".",
                xaxt = "n", yaxt = "n",
                colList = pairsColList(p.val.t, bucketCols = NULL,
                                       BWcutoff = 0),
                main.centered = TRUE, main = title)



gofCopula(tCopula(param = p_tm, df = df_tm), U, method = "SnC")
gofCopula(normalCopula(param = p_nm), U, method = "SnC")

?gofCopula

#Frank Copula MPL################################################################################################################
fitFr_m <- (fitCopula(frankCopula(dim = 2), data = U, method = "mpl")) ## The corresponding copula parameter estimate
summary(fitFr_m)
fitFr_m@copula@parameters

## Extract parameter estimates
len_frm<- length(coef(fitFr_m))
p_frm <- fitFr_m@copula@parameters[1]



#Gumbel-copula MPL#################################################################################
fitG_m <- (fitCopula(gumbelCopula(dim = 2), data = U, method = "mpl")) ## The corresponding copula parameter estimate
summary(fitG_m)
fitG_m@copula@parameters

## Extract parameter estimates
len_Gm <- length(coef(fitG_m))
p_Gm <- coef(fitG_m)[1]

gofCopula(gumbelCopula(param = p_Gm), U)

#CVC###################################################################################################################
CVC_gum <- xvCopula(gumbelCopula(), U)
CVC_fr <- xvCopula(frankCopula(), U)
CVC_joe <- xvCopula(joeCopula(), U)
CVC_cla <- xvCopula(claytonCopula(), U)
CVC_nor <- xvCopula(normalCopula(), U)
CVC_t <- xvCopula(tCopula(), U)
CVC_pla <- xvCopula(plackettCopula(), U)

CVC_t
CVC_nor
CVC_fr
CVC_gum

library(Copula.surv)
e <- matrix(1, ncol = 755, nrow = 1)
U_1 <- U[,1]
U_1
U_2 <- U[,2]
param_SG <- U1.Gumbel(U_1, U_2,e,e, lower = 0.01, upper = 10)
param_SG
testg <- Test.Gumbel(U_1, U_2, e, e) 

h <- simu.Gumbel(755, 1)
plot(h)

n=755 

Test.Gumbel(x.obs,y.obs,dx,dy)

##best fit seems mpl t-dist########################################################
#visual test with mpl t
library(lattice)
library(UsingR)
# needs to be an integer here (only) because of pCopula()
tc <- gumbelCopula(param = p_Gm)
set.seed(271)
X <- rCopula(755, copula = tc) # sample from the t copula
wireframe2(tc, FUN = dCopula, delta = 0.025) # density
contourplot2(tc, FUN = pCopula) # copula
contourplot2(tc, FUN = dCopula, n.grid = 42, cuts = 27) # density
par(mfrow = c(1,2))
plot(X, xlab = quote(X[1]), ylab = quote(X[2])) # scatter plot
plot(U)
plot(h)



##visual test with mpl normal###################################################
nc <- normalCopula(param = 0.5563)
set.seed(271)
X <- rCopula(742, copula = nc) # sample from the normal copula
wireframe2(nc, FUN = dCopula, delta = 0.025, main = "denisty plot 2008 Financial crisis") # density
contourplot2(nc, FUN = pCopula) # copula
contourplot2(nc, FUN = dCopula, n.grid = 42, cuts = 33, lwd = 1/2) # density
par(mfrow = c(1,2))
plot(X, xlab = quote(U[1]), ylab = quote(U[2])) # scatter plot
plot(U)

rm(list = ls())
