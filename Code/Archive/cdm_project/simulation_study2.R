# simulation study
library(CircStats)
rm(list = ls())

# 1 ####
# Choose a set of parameter values for two different locations. 
# These parameters would be the regression coefficients, not the pde constants.

b_01 <- 4.885; b_11 <- 0.022; b_12 <- -0.194; b_13 <- 0.191   
b_02 <- 7.780; b_21 <- -0.009; b_22 <- 5.089; b_23 <- -0.475

# 2 ####
# Decide on a simple g function for each location: the von mises density and choose some parameter values by hand.

mu1 <- pi/4; mu2 <- (5*pi)/4  #  location parameter μ is in [0, 2 π)
kappa1  <- 5; kappa2  <- 2  # concentration parameter κ > 0
# dvm(theta, mu, kappa)  # this is the form of the density function

# 3 ####
# Fix a relationship between two source locations i.e choose two points in 2d space (good separation)

x1 <- 0; y1 <- 0; x2 <- 10; y2 <- 10  # Location 1 is (0.0), location 2 is (1000,1000)
sources <- data.frame( "x" = c(x1,x2), "y" = c(y1,y2))
#plot(sources)

# 4 ####
# For a few sample sizes (maybe n = 30, 50, 80)
## a
##  Fix sample size and proportion of sample that will come from each model – some imbalance would be good, like a proportion around 0.2 or 0.3.

N <- 100  # 50, 100, 200
prop1 <- 0.7 * N  # proportion 1
prop2 <- N - prop1  # proportion 2
sim <- 1000
# store the source 1 coefficients in vectors so that I can get mean and sd
intercept.vec.1 <- rep(0, sim)
time.vec.1 <- rep(0, sim)
sinb1.vec.1 <- rep(0, sim)
sinb2.vec.1 <- rep(0, sim)
# store the source 2 coefficients in vectors so that I can get mean and sd
intercept.vec.2 <- rep(0, sim)
time.vec.2 <- rep(0, sim)
sinb1.vec.2 <- rep(0, sim)
sinb2.vec.2 <- rep(0, sim)
# ####
for (count in 1:sim){
## b time grid ####
## Generate two time grids, one for each sample size (random set of times from a discrete uniform distribution) 

time_grid1 <- runif(prop1, min = 50, max = 150)  # daily time grids
time_grid2 <- runif(prop2, min = 50, max = 150)

## c angle grid ####
## Generate two angle grids, one for each sample size (sample from uniform(0, 2pi))

angle_grid1 <- runif(prop1, min = 0, max = 2*pi)
angle_grid2 <- runif(prop2, min = 0, max = 2*pi)

## d y ####
## Generate data according to the two models by inputting the corresponding time grid, angle grid, and coefficients, 
## and adding random normal errors for each model. This will produce a set of responses log(1 + r/g(phi)). 
## Careful about the choice of error variance – shouldn’t be too big or too small.

y_generated1 <- rep(0, length(prop1))  # wrt source 1
for (i in 1:prop1) {
y_generated1[i] <- b_01 + (b_11 * time_grid1[i]) + (b_12 * sin(angle_grid1[i])) + 
                          (b_13 * (sin(angle_grid1[i] + pi/4))) + rnorm(1, 0, 0.5)  # μ = 0, var =  1 or 0.5
}

y_generated2 <- rep(0, length(prop2))  # wrt source 2
for (i in 1:prop2) {
  y_generated2[i] <- b_02 + (b_21 * time_grid2[i]) + (b_22 * sin(angle_grid2[i])) + 
                            (b_23 * (sin(angle_grid2[i] + pi/4))) + rnorm(1, 0, 0.5)  # μ = 0, var =  1 or 0.5
}

## e r ####
## Calculate r by back-transforming, and convert the phi, r pairs to Cartesian coordinates. 
## Use the actual g(phi) functions, no need to estimate them

# rho for source ####
# proportion for source 1
#
r_1 <- rep(0, prop1)  # rho distance of each point in proportion 1 from source 1
g_phi1 <- dvm(angle_grid1, mu1, kappa1)  # calculate g(phi) by applying the density function to the angle grid

for (i in 1:prop1){
r_1[i] <- ( g_phi1[i]*( (exp(y_generated1[i])) - 1))  
}
r_1 <- r_1

# rho for source 2 ####
# proportion for source 2
#
r_2 <- rep(0, prop2)  # rho distance of each point in proportion 2 from source 2
g_phi2 <- dvm(angle_grid2, mu2, kappa2)
for (i in 1:prop2){
  r_2[i] <- ( g_phi2[i] * ( (exp(y_generated2[i])) - 1))  # some are negative. I want positive distances
}
r_2 <- r_2 #+ 10

# x and y catersian for source 1 ####
x_1 <- rep(0, prop1)  
y_1 <- rep(0, prop1)  

for (i in 1:prop1){
  x_1[i] <- ( r_1[i] * cos(angle_grid1[i]) )   # x vector of each point in proportion 1
  y_1[i] <- ( r_1[i] * sin(angle_grid1[i]) )  # y vector of each point in proportion 1
}

# x and y catersian for source 2 ####
x_2 <- rep(0, prop2)  
y_2 <- rep(0, prop2)  
 
for (i in 1:prop2){
  x_2[i] <- ( r_2[i] * cos(angle_grid2[i]) )   # x vector of each point in proportion 2
  y_2[i] <- ( r_2[i] * sin(angle_grid2[i]) )  # y vector of each point in proportion 2
  }


# convert (r_1, theta) -> (x,y) -> (r2, theta2). 
# theta here is angle_grid1, r2 and theta2 are to be calculated
#
x_1_2 <- rep(0, prop1)  
y_1_2 <- rep(0, prop1)  
r_1_2 <- rep(0, prop1)  
theta_1_2 <- rep(0, prop1)  
for (i in 1:prop2){
  x_1_2[i] <- ( r_1[i] * cos(angle_grid1[i]) ) + 10  # 1000  # x vector of each point in proportion 1
  y_1_2[i] <- ( r_1[i] * sin(angle_grid1[i]) ) + 10  # 1000  # y vector of each point in proportion 1
  r_1_2[i] <- sqrt( (x_1_2[i])^2 + (y_1_2[i])^2 )  # rho distance of each point in proportion 1 from source 2
  theta_1_2[i] <- atan2(x_1_2[i], y_1_2[i])  # arctan (x/y) # angle of each point in proportion 1 from source 2
}
g_phi1_2 <- dvm(theta_1_2, mu1, kappa1)  # calculate g(phi) of each point in proportion 1 from source 2

# convert (r_2, theta) -> (x,y) -> (r1, theta1). 
#
x_2_1 <- rep(0, prop2)
y_2_1 <- rep(0, prop2)
r_2_1 <- rep(0, prop2)
theta_2_1 <- rep(0, prop2)
for (i in 1:prop2){
  x_2_1[i] <- ( r_2[i] * cos(angle_grid2[i]) ) - 10  # 1000
  y_2_1[i] <- ( r_2[i] * sin(angle_grid2[i]) ) - 10  # 1000
  r_2_1[i] <- sqrt( (x_2_1[i])^2 + (y_2_1[i])^2 )
  theta_2_1[i] <- atan2(x_2_1[i], y_2_1[i])
}
g_phi2_1 <- dvm(theta_2_1, mu2, kappa2)

## f ####
## Pool the simulated Cartesian coordinates and times and then convert to two sets of polar coordinates, one for each source
source1 <- data.frame("rho" = r_1, 
                       "theta" = angle_grid1,
                       "g.theta" = g_phi1,
                       "sinb1" = sin(angle_grid1), 
                       "sinb2" = sin(angle_grid1 + pi/4), 
                       "response" = y_generated1,
                       "time" = time_grid1,
                       "rho.src2" = r_1_2,
                       "theta.src2" = theta_1_2,
                       "g.theta.src2" = g_phi1_2) 

rho.src1 <- c(r_1, r_2_1); rho.src2 <- c(r_1_2, r_2)
g.phi.1 <- c(g_phi1, g_phi2_1); g.phi.2 <- c(g_phi2, g_phi1_2)

initdata <- data.frame(z.reg.1 = rho.src1/g.phi.1,  # first source regressor to initialize weights)
                       z.reg.2 = rho.src2/g.phi.2,  # second source regressor to initialize weights
                       z.fit = rho.src1 < rho.src2)  # is first source closer than second source?

source2 <- data.frame("rho" = r_2, 
                       "theta" = angle_grid2,
                       "g.theta" = g_phi2,
                       "sinb1" = sin(angle_grid2), 
                       "sinb2" = sin(angle_grid2 + pi/4), 
                       "response" = y_generated2,
                       "time" = time_grid2,
                       "rho.src1" = r_2_1,
                       "theta.src1" = theta_2_1,
                       "g.theta.src1" = g_phi2_1)

## g
## Apply the fitting procedure, and output the estimated coefficients.
source('fit_twosource_fn.R')
fit_out <- fit_ts(source1, 
                  source2, 
                  initdata, 
                  isotropy = F)
fit_src1 <- fit_out$fit1; confidence.1 <- fit_out$conf1
fit_src2 <- fit_out$fit2; confidence.2 <- fit_out$conf2
pr_z <- fit_out$prz
#print(fit_src1$coefficients)
#print(fit_src2$coefficients)

# store the source 1 coefficients in vectors so that I can get mean and sd
intercept.vec.1[count] <- as.numeric(fit_src1$coefficients[1])
time.vec.1[count] <- as.numeric(fit_src1$coefficients[2])
sinb1.vec.1[count] <- as.numeric(fit_src1$coefficients[3])
sinb2.vec.1[count] <- as.numeric(fit_src1$coefficients[4])

# store the source 2 coefficients in vectors so that I can get mean and sd
intercept.vec.2[count] <- as.numeric(fit_src2$coefficients[1])
time.vec.2[count] <- as.numeric(fit_src2$coefficients[2])
sinb1.vec.2[count] <- as.numeric(fit_src2$coefficients[3])
sinb2.vec.2[count] <- as.numeric(fit_src2$coefficients[4])
}

# calculate the means source 1
print(mean(intercept.vec.1, na.rm = T))
print(mean(time.vec.1, na.rm = T))
print(mean(sinb1.vec.1, na.rm = T))
print(mean(sinb2.vec.1, na.rm = T))

# source 2
print(mean(intercept.vec.2, na.rm = T))
print(mean(time.vec.2, na.rm = T))
print(mean(sinb1.vec.2, na.rm = T))
print(mean(sinb2.vec.2, na.rm = T))

# calculate the sd source 1
print(sd(intercept.vec.1), na.rm = T)
print(sd(time.vec.1, na.rm = T))
print(sd(sinb1.vec.1, na.rm = T))
print(sd(sinb2.vec.1, na.rm = T))

# source 2
print(sd(intercept.vec.2, na.rm = T))
print(sd(time.vec.2, na.rm = T))
print(sd(sinb1.vec.2, na.rm = T))
print(sd(sinb2.vec.2, na.rm = T))

# plots
#dev.off()
plot(x_1, y_1, col = 2, ylim = c(-50,100), xlim = c(-50,100), xlab = "x", ylab = "y", main = "Kappa = 5,5, VAR = 0.5")
points(sources$x, sources$y, col = 1, pch = 8)
points(x_2 + 10, y_2 + 10, col = 3)


  