# Make the magnetism data from the raw data.
# The raw data is from
# Creer, Irving and Nairn (1959),
# Geophysical J. of the Royal Astron. Soc., 2, 306--323.

Ds <- c(
176, 185, 196, 179, 186,
192, 183, 203, 193, 171,
191, 190, 187, 203, 196,
205, 182, 190, 183, 190,
183, 179, 190, 177, 195,
182, 182, 190, 183, 194,
186, 192, 192, 180)

Is <- c(
-10,  -6,  -5,  -4,  -5,
  3,  -1, -11, -28, -14,
-14, -20, -15, -26, -19,
  6,   9, -13,   1,  15,
  7,  11,  13, -10,  11,
 11,   7, -15, -20,  -7,
 -6,  -9,   0,   0)

# Eq.(1.3) in Kent (1982), J.Roy. Statist. Soc. B, 44, 71--80.
thetas <- (pi/180) * (90 + Is)
phis <- (pi/180) * (360 - Ds)
x1 <- cos(thetas)
x2 <- sin(thetas) * cos(phis)
x3 <- sin(thetas) * sin(phis)

#x <- cbind(x1,x2,x3)
x <- cbind(x2,x3,x1)
n <- nrow(x)

xmeansq <- (t(x)%*%x) / n
xmean <- apply(x,2,mean)

SA <- round(xmeansq,3)
SB <- round(xmean,3)

#> SA
#       x1     x2     x3
#x1  0.045 -0.075  0.014
#x2 -0.075  0.921 -0.122
#x3  0.014 -0.122  0.034
#> SB
#    x1     x2     x3 
# 0.082 -0.959  0.131 

# wrong data!
# SBkent <- c(0.083,-0.959,-0.131)

