# simulation study
library(CircStats)
library(ggplot2)
rm(list = ls())

#set.seed(10)

# 1 ####
# Choose a set of parameter values for two different locations. 
# These parameters would be the regression coefficients, not the pde constants.

# b_01 <- 4.89; b_11 <- 0.02; b_12 <- -0.19; b_13 <- 0.3  # rounded off true values 
b_01 <- 4.85; b_11 <- 0.03; b_12 <- 0.5; b_13 <- 0   
# b_02 <- 5; b_21 <- 0.03; b_22 <- 0.689; b_23 <- 1.83
b_02 <- 3.75; b_21 <- 0.02; b_22 <- 0; b_23 <- -1

# 2 ####
# Decide on a simple g function for each location: the von mises density and choose some parameter values by hand.

mu1 <- pi/4; mu2 <- (5*pi)/4  #  location parameter μ is in [0, 2 π)
# kappa1  <- 5; kappa2  <- 5  # concentration parameter κ > 0
kappa1  <- 2; kappa2  <- 2; 

# 3 ####
# Fix a relationship between two source locations i.e choose two points in 2d space (good separation)

source_dis <- 2000
x1 <- 0; y1 <- 0; x2 <- source_dis; y2 <- source_dis  # Location 1 is (0.0), location 2 is (1000,1000)
sources <- data.frame( "x" = c(x1,x2), "y" = c(y1,y2))

# 4 ####
# For a few sample sizes (maybe n = 30, 50, 80)
N <- 50  # 50, 100, 200
prop1 <- 0.7 * N  # proportion 1
prop2 <- N - prop1  # proportion 2
sim <- 100
source_id <- c(rep(1, prop1), rep(2, prop2))
#
# var <- 1
var <- 0.5

# 4a) store stuff #####
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
# store probabilities
probability <- rep(0, sim)
probability.source1 <- rep(0, sim)
probability.source2 <- rep(0, sim)

# 5 run 1000 simulations ####
for (count in 1:sim){
  ## b time grid ####
  ## Generate two time grids, one for each sample size (random set of times from a discrete uniform distribution) 
  
  # count <- 1
  
  time_grid1 <- runif(prop1, min = 50, max = 150)  # daily time grids
  time_grid2 <- runif(prop2, min = 50, max = 150)
  # time_grid2 <- runif(prop2, min = 100, max = 150) # late time 2 grid
  
  # missing time ####
  missing <- 0.2 # 0.2 0.3 0.4 0.5
  
  # missing.prop1 <- sample(time_grid1, missing*prop1)
  # for (l in missing.prop1) {
  # time_grid1[time_grid1 == l] <- 0
  # }
  # missing.prop2 <- sample(time_grid2, missing*prop2)
  # for (t in missing.prop2) {
  #   time_grid2[time_grid2 == t] <- 0
  # }
  # c angle grid
  # Generate two angle grids, one for each sample size (sample from uniform(0, 2pi))
  
  angle_grid1 <- runif(prop1, min = 0, max = 2*pi)
  angle_grid2 <- runif(prop2, min = 0, max = 2*pi)
  
  # missing angle ####
  # missing.prop1.angle <- sample(angle_grid1, missing*prop1)
  # for (p in missing.prop1.angle) {
  #   angle_grid1[angle_grid1 == p] <- 0
  # }
  # missing.prop2.angle <- sample(angle_grid2, missing*prop2)
  # for (q in missing.prop2.angle) {
  #   angle_grid2[angle_grid2 == q] <- 0
  # }
  # 
  
  ## d generate y ####
  ## Generate data according to the two models by inputting the corresponding time grid, angle grid, and coefficients
  
  y_generated1 <- rep(0, length(prop1))  # wrt source 1
  for (i in 1:prop1) {
    y_generated1[i] <- b_01 + (b_11 * time_grid1[i]) + (b_12 * sin(angle_grid1[i])) + 
      (b_13 * (sin(angle_grid1[i] + pi/4))) + rnorm(1, 0, var)  # μ = 0, var =  1 or 0.5
  }
  
  y_generated2 <- rep(0, length(prop2))  # wrt source 2
  for (i in 1:prop2) {
    y_generated2[i] <- b_02 + (b_21 * time_grid2[i]) + (b_22 * (sin(angle_grid2[i]))) + 
      (b_23 * (sin(angle_grid2[i] + pi/4))) + rnorm(1, 0, var)  
    # μ = 0, var =  1 or 0.5
  }
  
  ## e calculate rho ####
  ## Calculate r by back-transforming, and convert the phi, r pairs to Cartesian coordinates. 
  ## Use the actual g(phi) functions, no need to estimate them
  
  ### 1. rho for source ####
  
  r_1 <- rep(0, prop1)  # rho distance of each point in proportion 1 from source 1
  g_phi1 <- dvm(angle_grid1, mu1, kappa1)  # calculate g(phi) by applying the density function to the angle grid
  
  for (i in 1:prop1){
    r_1[i] <- ( g_phi1[i]*( (exp(y_generated1[i])) - 1))  # I want positive distances
  }
  
  ### 2. rho for source 2 ####
  
  r_2 <- rep(0, prop2)  # rho distance of each point in proportion 2 from source 2
  g_phi2 <- dvm(angle_grid2, mu2, kappa2)
  for (i in 1:prop2){
    r_2[i] <- ( g_phi2[i] * ( (exp(y_generated2[i])) - 1))  # some are negative. I want positive distances
  }
  #r_2 <- r_2 + source_dis
  
  #### 3. x and y catersian for source 1 ####
  x_1 <- rep(0, prop1)  
  y_1 <- rep(0, prop1)  
  
  for (i in 1:prop1){
    x_1[i] <- ( r_1[i] * cos(angle_grid1[i]) )   # x vector of each point in proportion 1
    y_1[i] <- ( r_1[i] * sin(angle_grid1[i]) )  # y vector of each point in proportion 1
  }
  
  #### 4. x and y catersian for source 2 ####
  x_2 <- rep(0, prop2)  
  y_2 <- rep(0, prop2)  
  
  for (i in 1:prop2){
    x_2[i] <- ( r_2[i] * cos(angle_grid2[i]) )   # x vector of each point in proportion 2
    y_2[i] <- ( r_2[i] * sin(angle_grid2[i]) )  # y vector of each point in proportion 2
  }
  
  #### 5 the x, y rho, theta and of points in prop 1 from source 2 ####
  
  x_1_2 <- rep(0, prop1)  # the x coordinate of points in prop 1 from source 2
  y_1_2 <- rep(0, prop1)  # the y coordinate of points in prop 1 from source 2
  r_1_2 <- rep(0, prop1)  # the rho of points in prop 1 from source 2
  theta_1_2 <- rep(0, prop1)  # the angle of points in prop 1 from source 2
  
  for (i in 1:prop1){
    x_1_2[i] <- ( (r_1[i]) * cos(angle_grid1[i]) ) + source_dis  # x vector of each point in proportion 1
    y_1_2[i] <- ( (r_1[i]) * sin(angle_grid1[i]) ) + source_dis  # y vector of each point in proportion 1
    r_1_2[i] <- sqrt( (x_1_2[i])^2 + (y_1_2[i])^2 )  # rho distance of each point in proportion 1 from source 2
    theta_1_2[i] <- atan2(x_1_2[i], y_1_2[i])  # arctan (x/y) # angle of each point in proportion 1 from source 2
  }
  g_phi1_2 <- dvm(theta_1_2, mu1, kappa1)  # calculate g(phi) of each point in proportion 1 from source 2
  
  #### 6 the x, y rho, theta and of points in prop 2 from source 1 ####
  
  x_2_1 <- rep(0, prop2)
  y_2_1 <- rep(0, prop2)
  r_2_1 <- rep(0, prop2)
  theta_2_1 <- rep(0, prop2)
  
  for (i in 1:prop2){
    x_2_1[i] <- ( (r_2[i] ) * cos(angle_grid2[i]) ) - source_dis 
    y_2_1[i] <- ( (r_2[i] ) * sin(angle_grid2[i]) ) - source_dis  
    r_2_1[i] <- sqrt( (x_2_1[i])^2 + (y_2_1[i])^2 )
    theta_2_1[i] <- atan2(x_2_1[i], y_2_1[i])  
  }
  g_phi2_1 <- dvm(theta_2_1, mu2, kappa2)
  
  ## f prepare for estimation ####
  # source 1 - first entry in the dataframe
  
  # source labels
  labels_src <- c(rep(1, prop1), rep(2, prop2))
  
  # coordinates relative to source 1
  rho.src1 <- c(r_1, r_2_1); 
  theta1 <- c(angle_grid1, theta_2_1); 
  g.phi1 <- c(g_phi1, g_phi2_1)
  
  sinb1a <- as.numeric(sin(theta1)); 
  sinb1b <- as.numeric(sin(theta1 + pi/4))
  response1 <- c(y_generated1, log(1 + r_2_1/g_phi2_1));
  time1 <- c(time_grid1, time_grid2)
  
  # coordinates relative to source 2
  rho.src2 <- c(r_1_2, r_2); 
  theta.src2 <- c(theta_1_2, angle_grid2); 
  g.phi1_2 <- c(g_phi1_2, g_phi2)
  
  source1 <- data.frame("rho" = rho.src1, 
                        "theta" = theta1,
                        "g.theta" = g.phi1,
                        "sinb1" = sinb1a, 
                        "sinb2" = sinb1b, 
                        "response" = response1,
                        "time" = time1) 
  
  # source 2 - second entry in the dataframe
  
  rho.src2 <- c(r_1_2, r_2); 
  theta2 <- c(theta_1_2, angle_grid2); 
  g.phi2 <- c(g_phi1_2, g_phi2)
  
  sinb2a <- as.numeric(sin(theta2)); sinb2b <- as.numeric(sin(theta2 + pi/4))
  response2 <- c(log(1 + r_1_2/g_phi1_2), y_generated2)
  time2 <- c(time_grid1, time_grid2)
  theta.src1 <- c(angle_grid1, theta_2_1); 
  g.phi2_1 <- c(g_phi1, g_phi2_1) 
  
  source2 <- data.frame("rho" = rho.src2, 
                        "theta" = theta2,
                        "g.theta" = g.phi2,
                        "sinb1" = sinb2a, 
                        "sinb2" = sinb2b, 
                        "response" = response2,
                        "time" = time2)
  # initial data
  initdata <- data.frame(z.reg.1 = rho.src1/g.phi1,  # first source regressor to initialize weights)
                         z.reg.2 = rho.src2/g.phi2,  # second source regressor to initialize weights
                         z.fit = rho.src1 < rho.src2)  # is first source closer than second source?
  ## g estimation ####
  ## Apply the fitting procedure, and output the estimated coefficients.
  source('fit_twosource_fn.R')
  fit_out <- fit_ts(source1, 
                    source2, 
                    initdata, 
                    isotropy = F)
  fit_src1 <- fit_out$fit1; confidence.1 <- fit_out$conf1
  fit_src2 <- fit_out$fit2; confidence.2 <- fit_out$conf2
  pr_z <- fit_out$prz
  
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
  
  # MODIFICATION
  estimated.labels <- as.numeric(pr_z < 0.5) + 1
  estimated.labels1 <- estimated.labels[1:prop1] 
  estimated.labels2 <- estimated.labels[(prop1+1):N]
  probability.source1[count] <- pct.correct.source1 <- mean(estimated.labels1 == 1)
  probability.source2[count] <- pct.correct.source2 <- mean(estimated.labels2 == 2)
  
  # maybe just store this and report an average?
  probability[count] <- (mean((pr_z < 0.5) + 1 == source_id)) # looks good
}

# calculate the means source 1 ####
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

# probabilities
plot(probability)

# plot ####
plot(probability.source1, xlim = c(0, sim), ylim = c(0, 1.2), col = "red", ylab = "Proportion correctly classified", xlab = "Simulation", pch = 16, cex = .4)
points(probability.source2, col = "blue", pch = 16, cex = .4)
legend("topleft", legend = c("Source 1", "Source 2"),
       col = c("red", "blue"), pch = 16, cex = .8)

plot.df.source1 <- data.frame(x1 = x_1, y1 = y_1, Shape = 19, Day = time_grid1)
plot.df.source2 <- data.frame(x2 = x_2, y2 = y_2, Shape = 17, Day = time_grid2)
mylegend <- c()

require(ggplot2)
p <- ggplot() + 
  geom_point(aes(x = x1, y = y1, col = Day, shape = Shape), data = plot.df.source1) + 
  geom_point(aes(x = x2 + source_dis, y = y2 + source_dis, col = Day, shape = Shape),
             data = plot.df.source2) + ylim(c(min(y_1, y_2), max(y_1, y_2))) + xlim(c(min(x_1, x_2), max(x_1, x_2))) +
  scale_colour_gradientn(colours = rev(rainbow(3))) + 
  scale_shape_identity(guide = "legend", breaks = c(17,19), 
                       labels = c("(2000,2000)", "(0,0)"), name = "Source" ) + 
  geom_point(aes(x = sources$x, y = sources$y), shape = 15, data = sources) +
  theme_bw() + xlab("x") + ylab("y") +
  ggtitle(paste('kappa1 = ', kappa1, ',', 'kappa2 = ', kappa2, ',', 'var = ', var, sep = " ")) +
  theme(plot.title = element_text(hjust = 0.5))
print(p)
