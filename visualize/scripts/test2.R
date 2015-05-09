library(mvtnorm) # References rmvnorm()
library(ellipse) # References ellipse()
set.seed(17)
 
# Set the covariance matrix
sigma2 <- matrix(c(9, 0, 0, 1), ncol=2)
 
# Set the means
mu <- c(0,0)
 
# Get the correlation matrix
#P <- cov2cor(sigma2)
P <- sigma2
 
# Generate the data
p <- rmvnorm(n=1000, mean=mu, sigma=sqrt(sigma2))
 
# Plot the data
plot(p,xlim=c(-6,6),ylim=c(-6,6))
 
# Plot the ellipse
#lines( ellipse( P, centre = c(0,0)) , col='red')

evals <- eigen(P)$values
evecs <- eigen(P)$vectors
# Angles of a circle
a <- seq(0, 2*pi, len=100)
 
# Get critical value
c2 <- qchisq(0.39, 2)
c <- sqrt(c2)

# Get the distances
xT <- c * sqrt(evals[1]) * cos(a)
yT <- c * sqrt(evals[2]) * sin(a)
 
M <- cbind(xT, yT)

# Covert the coordinates
transM <- evecs %*% t(M)
transM <- t(transM)

lines(transM + mu)
