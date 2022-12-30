#EXAMPLES OF COPULAS TO BE USED IN THE THESIS#############################
#Independence copula######################################################
d <- 2
ic <- indepCopula(dim = d)
set.seed(2008)
u <-runif(d)
(Pi <- pCopula(u, copula = ic)) #simple multiplication
wireframe2(ic, FUN = pCopula, col.4 = adjustcolor("black", alpha.f = 0.25)) #surface plot
contourplot2(ic, FUN = pCopula) #contour plot

##Gaussian copula######################################################
nc <- normalCopula(iTau(normalCopula(), tau = 0.5))
set.seed(271)
U <- rCopula(1000, copula = nc) # sample from the normal copula
wireframe2(nc, FUN = dCopula, delta = 0.025) # density
contourplot2(nc, FUN = pCopula) # copula
contourplot2(nc, FUN = dCopula, n.grid = 42, cuts = 33, lwd = 1/2) # density
plot(U, xlab = quote(U[1]), ylab = quote(U[2])) # scatter plot

###T-Copula############################################################
nu <- 4 # needs to be an integer here (only) because of pCopula()
tc <- tCopula(iTau(tCopula(df = nu), tau = 0.5), df = nu)
set.seed(271)
X <- rCopula(1000, copula = tc) # sample from the t copula
wireframe2(tc, FUN = dCopula, delta = 0.025) # density
contourplot2(tc, FUN = pCopula) # copula
contourplot2(tc, FUN = dCopula, n.grid = 42, cuts = 27) # density
plot(X, xlab = quote(U[1]), ylab = quote(X[2])) # scatter plot

#generator function for the Frank########################################
library(copula)
x = seq(0.000001,1,length = 500)
frank = iPsi(copula = archmCopula(family = "frank", param = 1),x)
plot(x,frank,
     type = "l",
     lwd = 3,
     ylab = expression(phi(x)))
abline(h = 0)
abline(v = 0)

##Frank Copula###########################################################