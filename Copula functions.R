#EXAMPLES OF COPULAS TO BE USED IN THE THESIS#############################
library(copula)
library(lattice)
library(UsingR)
#Independence copula######################################################
d <- 2
ic <- indepCopula(dim = d)
set.seed(2008)
u <-runif(d)
(Pi <- pCopula(u, copula = ic)) #simple multiplication
wireframe2(ic, FUN = pCopula, col.4 = adjustcolor("black", alpha.f = 0.25)) #surface plot
contourplot2(ic, FUN = pCopula) #contour plot

#Co-monotonic and Counter-monotonic copula#####################################################
U <- runif(200)
u <- seq(0, 1, length.out = 40)
u12 <- expand.grid("u[1]" = u, "u[2]" = u)
W <- pmax(u12[,1] + u12[,2] - 1, 0)
M <- pmin(u12[,1], u12[,2])
val.W <- cbind(u12, "Cminus(u[1],u[2])" = W)
val.M <- cbind(u12, "Cplus(u[1],u[2])" = M)
par(mfrow = c(2, 2), mar = c(4, 4, 2, 2), oma = c(0, 0, 2, 0))
plot(cbind(U, 1-U), xlab = quote(U[1]), ylab = quote(U[2]))
plot(cbind(U, U), xlab = quote(U[1]), ylab = quote(U[2]))
wireframe2(val.W, col = "blue", size = 0.5)
wireframe2(val.M, col = "red", size = 0.5)
contourplot2(val.W, xlim = 0:1, ylim = 0:1)
contourplot2(val.M, xlim = 0:1, ylim = 0:1)

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
set.seed(5640)
theta = c(-100, -50, -10, -1, 0, 5, 20, 50, 500)
par(mfrow = c(3,3), cex.axis = 1.2, cex.lab = 1.2, cex.main = 1.2)
par(mar=c(4,4,2,1))
for(i in 1:9){
  U = rCopula(n = 500, 
              copula = archmCopula(family = "frank", param = theta[i]))
  plot(U, xlab = expression(u[1]), ylab = expression(u[2]), main = eval(substitute(expression(paste(theta,"=",j)),
                                                                                   list(j = as.character(theta[i])))))
}
##Clayton Copula########################################################
set.seed(5640)
theta = c(-0.98, -0.7, -0.3, -0.1, 0.1, 1, 5, 15, 100)
par(mfrow = c(3,3), cex.axis = 1.2, cex.lab = 1.2, cex.main = 1.2)
par(mar=c(4,4,2,1))
for(i in 1:9){
  U = rCopula(n = 500, 
              copula = archmCopula(family = "clayton", param = theta[i]))
  plot(U, xlab = expression(u[1]), ylab = expression(u[2]), main = eval(substitute(expression(paste(theta,"=",j)),
                                                                                   list(j = as.character(theta[i])))))
}
##Gumbel Copula###########################################################
set.seed(5640)
theta = c(1.1, 1.5, 2, 4, 8, 50)
par(mfrow = c(2,3), cex.axis = 1.2, cex.lab = 1.2, cex.main = 1.2)
for(i in 1:6){
  U = rCopula(n = 500, 
              copula = archmCopula(family = "gumbel", param = theta[i]))
  plot(U, xlab = expression(u[1]), ylab = expression(u[2]), main = eval(substitute(expression(paste(theta,"=",j)),
                                                                                   list(j = as.character(theta[i])))))
}
##Joe Copula##########################################################
set.seed(5640)
theta = c(1.1, 1.5, 2, 4, 8, 50)
par(mfrow = c(2,3), cex.axis = 1.2, cex.lab = 1.2, cex.main = 1.2)
for(i in 1:6){
  U = rCopula(n = 500, 
              copula = archmCopula(family = "joe", param = theta[i]))
  plot(U, xlab = expression(u[1]), ylab = expression(u[2]), main = eval(substitute(expression(paste(theta,"=",j)),
                                                                                   list(j = as.character(theta[i])))))
}
###LogNormal correlation max and min####################################
CorBoundLN <- function(s,bound = c("max","min"))
{
  if(!is.matrix(s)) s <- rbind(s)
  bound <- match.arg(bound)
  if(bound == "min") s[,2] <- -s[,2]
  (exp((s[,1] + s[,2])^2/2)-exp((s[,1]^2 + s[,2]^2)/2))/
    sqrt(expm1(s[,1]^2) * exp(s[,1]^2) * expm1(s[,2]^2) * exp(s[,2]^2))
}

s <- seq(0.01, 5, length.out = 20)
s12 <- expand.grid("sigma[1]" = s,"sigma[2]" = s)
wireframe2(cbind(s12, '(Cor_min)(sigma[1],sigma[2])'= CorBoundLN(s12, bound = "min")))
wireframe2(cbind(s12, '(Cor_max)(sigma[1],sigma[2])'= CorBoundLN(s12)))

rm(list =ls())
