# simulation study
library(CircStats)
rm(list = ls())

# 1 ####
# Choose a set of parameter values for two different locations. 
# These parameters would be the regression coefficients, not the pde constants.

b_01 <- 4.89; b_11 <- 0.02; b_12 <- -0.19; b_13 <- 0.3  # rounded off true values 
b_02 <- 18.6; b_21 <- 0.03; b_22 <- 6.89; b_23 <- 18.3

# 2 ####
# random wind
wind.data.source1 <- runif(24, min = 0, max = 2*pi)
wind.data.source2 <- runif(24, min = 0, max = 2*pi)

# 3 ####
# Fix a relationship between two source locations i.e choose two points in 2d space (good separation)

x1 <- 0; y1 <- 0; x2 <- 10; y2 <- 10  # Location 1 is (0.0), location 2 is (10,10)
sources <- data.frame( "x" = c(x1, x2), "y" = c(y1, y2))

# 4 ####
# For a few sample sizes (maybe n = 30, 50, 80)
## a
##  Fix sample size and proportion of sample that will come from each model – some imbalance would be good, like a proportion around 0.2 or 0.3.

N <- 100  # 50, 100, 200
prop1 <- 0.7 * N  # proportion 1
prop2 <- N - prop1  # proportion 2
sim <- 3
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
  
  ## e) r ####
  ## Calculate r by back-transforming, and convert the phi, r pairs to Cartesian coordinates. 
  ## Use the actual g(phi) functions, no need to estimate them
  

  # estimate g_phi for source 1 ####
  g_phi1 <- rep(0, prop1)
  require(circular)
  for (i in 1:prop1){
  g_phi1[i] <- density.circular(x = wind.data.source1, 
                             z = angle_grid1[i], 
                             bw = 50)$y  
  }
  
  # rho for source ####
  # proportion for source 1
  r_1 <- rep(0, prop1)  # rho distance of each point in proportion 1 from source 1
  for (i in 1:prop1){
    r_1[i] <- ( g_phi1[i]*( (exp(y_generated1[i])) - 1))  # I want positive distances
  }
  
  # estimate g_phi for source 2 ####
  g_phi2 <- rep(0, prop2)
  require(circular)
  for (i in 1:prop2){
    g_phi2[i] <- density.circular(x = wind.data.source2, 
                                  z = angle_grid2[i], 
                                  bw = 50)$y  
  }
  
  # rho for source 2 ####
  # proportion for source 2
  #
  r_2 <- rep(0, prop2)  # rho distance of each point in proportion 2 from source 2
  for (i in 1:prop2){
    r_2[i] <- ( g_phi2[i] * ( (exp(y_generated2[i])) - 1))  # some are negative. I want positive distances
  }
  
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
  g_phi1_2 <- rep(0, prop1)
  for (i in 1:prop2){
        x_1_2[i] <- ( r_1[i] * cos(angle_grid1[i]) ) + 10  # 1000  # x vector of each point in proportion 1
        y_1_2[i] <- ( r_1[i] * sin(angle_grid1[i]) ) + 10  # 1000  # y vector of each point in proportion 1
        r_1_2[i] <- sqrt( (x_1_2[i])^2 + (y_1_2[i])^2 )  # rho distance of each point in proportion 1 from source 2
        theta_1_2[i] <- atan2(x_1_2[i], y_1_2[i])  # arctan (x/y) # angle of each point in proportion 1 from source 2
        g_phi1_2[i] <- density.circular(x = wind.data.source1, 
                                       z = theta_1_2[i], 
                                       bw = 50)$y    # calculate g(phi) of each point in proportion 1 from source 2
  }
  # convert (r_2, theta) -> (x,y) -> (r1, theta1). 
  #
  x_2_1 <- rep(0, prop2)
  y_2_1 <- rep(0, prop2)
  r_2_1 <- rep(0, prop2)
  theta_2_1 <- rep(0, prop2)
  g_phi2_1 <- rep(0, prop2)
  for (i in 1:prop2){
        x_2_1[i] <- ( r_2[i] * cos(angle_grid2[i]) ) - 10  # 1000
        y_2_1[i] <- ( r_2[i] * sin(angle_grid2[i]) ) - 10  # 1000
        r_2_1[i] <- sqrt( (x_2_1[i])^2 + (y_2_1[i])^2 )
        theta_2_1[i] <- atan2(x_2_1[i], y_2_1[i])
        g_phi2_1 <- density.circular(x = wind.data.source2, 
                                     z = theta_2_1[i], 
                                     bw = 50)$y    # calculate g(phi) of each point in proportion 1 from source 2
  }

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
print(mean(intercept.vec.1))
print(mean(time.vec.1))
print(mean(sinb1.vec.1))
print(mean(sinb2.vec.1))

# source 2
print(mean(intercept.vec.2))
print(mean(time.vec.2))
print(mean(sinb1.vec.2))
print(mean(sinb2.vec.2))

# calculate the sd source 1
print(sd(intercept.vec.1))
print(sd(time.vec.1))
print(sd(sinb1.vec.1))
print(sd(sinb2.vec.1))

# source 2
print(sd(intercept.vec.2))
print(sd(time.vec.2))
print(sd(sinb1.vec.2))
print(sd(sinb2.vec.2))

# plots
#dev.off()
# plot(x_1, y_1, col = time_grid1, ylim = c(-100,100), xlim = c(-100,100), xlab = "x", ylab = "y", main = "Kappa = 2,2, VAR = 0.5")
# points(sources$x, sources$y, col = 1, pch = 8)
# points(x_2 + 10, y_2 + 10, col = time_grid2)

plot.df.source1 <- data.frame(x1 = x_1, y1 = y_1, Shape = 19, Day = time_grid1)
plot.df.source2 <- data.frame(x2 = x_2, y2 = y_2, Shape = 17, Day = time_grid2)
mylegend <- c()

require(ggplot2)
p <- ggplot() + 
  geom_point(aes(x = x1, y = y1, col = Day, shape = Shape), data = plot.df.source1) + 
  geom_point(aes(x = x2 + 10, y = y2 + 10, col = Day, shape = Shape), data = plot.df.source2) + ylim(c(-300,300)) + xlim(c(-300,300)) +
  scale_colour_gradientn(colours = rev(rainbow(3))) + scale_shape_identity() +
  geom_point(aes(x = sources$x, y = sources$y), shape = 15, data = sources) +
  theme_bw() +
  ggtitle("Estimated g(phi)") +
  theme(plot.title = element_text(hjust = 0.5))
print(p)
