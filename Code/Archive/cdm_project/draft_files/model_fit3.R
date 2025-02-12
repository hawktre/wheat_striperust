fit_pred_twosource <- function(yr, training_prop = 0.8, nb = 1) {
  
  # load("~/Documents/Trevor/postprocessed_cdm_data.rdata")
  filename <- paste("cdm_multisource_test_data", yr, ".RData", sep="")
  load(filename)
  
  n <- length(DATA$report_date)
  # nb <- 1
  DATA$month <- as.numeric(unlist(strsplit(as.character(DATA$report_date),
                                           "-")))[seq(2, 3*n-1, 3)]
  DATA$day <- as.numeric(unlist(strsplit(as.character(DATA$report_date),
                                         "-")))[seq(3, 3*n, 3)]
  dm <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
  DATA$time <- rep(0, n)
  for(i in 1:n) {
    DATA$time[i] <- sum(dm[1:(DATA$month[i]-1)]) + DATA$day[i]
  }
  DATA$response1 <- log(1 + DATA$rho1/DATA$g.theta.1)
  DATA$response2 <- log(1 + DATA$rho2/DATA$g.theta.2)
  # nb <- 1
  sinb1 <- matrix(0, nrow = n, ncol = 2*nb)
  cosb1 <- matrix(0, nrow = n, ncol = 2*nb)
  sinb2 <- matrix(0, nrow = n, ncol = 2*nb)
  cosb2 <- matrix(0, nrow = n, ncol = 2*nb)
  #forml <- "response ~ time "
  for(k in seq(1, 2*nb-1, 2)){
    sinb1[, k] <- sin(k*as.numeric(DATA$theta1))
    sinb1[, (k+1)] <- sin(k*as.numeric(DATA$theta1) + pi/4)
    sinb2[, k] <- sin(k*as.numeric(DATA$theta2))
    sinb2[, (k+1)] <- sin(k*as.numeric(DATA$theta2) + pi/4)
    #  sinb[, (k+2)] <- sin(k*as.numeric(DATA$theta) + 3*pi/4)
    cosb1[, k] <- cos(k*as.numeric(DATA$theta1))
    cosb1[, (k+1)] <- cos(k*as.numeric(DATA$theta1) + pi/4)
    cosb2[, k] <- cos(k*as.numeric(DATA$theta2))
    cosb2[, (k+1)] <- cos(k*as.numeric(DATA$theta2) + pi/4)
    #  cosb[, (k+2)] <- cos(k*as.numeric(DATA$theta) + pi/2)
    #  forml <- paste(forml, "+ sinb.", k, "+ cosb.", k, sep="")
  }
  DATA$sinb1 <- sinb1
  DATA$cosb1 <- cosb1
  DATA$sinb2 <- sinb2
  DATA$cosb2 <- cosb2
  
  # yr = 2008
#  data1 <- DATA[DATA$year == yr, ]
  data1 <- DATA
  # training_prop <- 0.8
  n1 <- length(data1$report_date)
  n2 <- floor(training_prop*n1)
  data2 <- data1[1:n2, ]
  z_reg1 <- DATA$rho1/DATA$g.theta.1
  z_reg2 <- DATA$rho2/DATA$g.theta.2
  z_fit <- rep(1, n2)
  for(i in 1:n2) {
    if(DATA$rho1[i] > DATA$rho2[i]) {
      z_fit[i] <- 0
    }
  }
  data2$z_fit <- z_fit
  data2$z_reg1 <- z_reg1[1:n2]
  data2$z_reg2 <- z_reg2[1:n2]
  z_glm <- glm(z_fit ~ z_reg1 + z_reg2, data = data2, family = binomial())
  pr_z <- fitted(z_glm)
  
  err <- 1
  while(err > 0.001) {
  lmod1 <- lm(response1 ~ time + sinb1, data = data1[1:n2, ],
              weights = pr_z)
  mean1 <- residuals(lmod1)/summary(lmod1)$sigma 
  lmod2 <- lm(response2 ~ time + sinb2, data = data1[1:n2, ],
              weights = (1 - pr_z))
  mean2 <- residuals(lmod2)/summary(lmod2)$sigma 
  pr_z1 <- (dnorm(mean1)*pr_z)/(dnorm(mean1)*pr_z +
                                   dnorm(mean2)*(1 - pr_z))
  err <- mean(abs(pr_z - pr_z1))
  pr_z <- pr_z1
  }
  summary(lmod1)
  summary(lmod2)
  
  datapr <- data1[(n2+1):n1, ]
  z_fit <- rep(1, (n1 - n2))
  for(i in (n2+1):n1) {
    if(DATA$rho1[i] > DATA$rho2[i]) {
      z_fit[(i-n2)] <- 0
    }
  }
  datapr$z_fit <- z_fit
  datapr$z_reg1 <- z_reg1[(n2+1):n1]
  datapr$z_reg2 <- z_reg2[(n2+1):n1]
  z_glm_pr <- glm(z_fit ~ z_reg1 + z_reg2, data = datapr, family = binomial())
  pr_z2 <- fitted(z_glm_pr)
  pr_value1 <- predict(lmod1, datapr, weights = pr_z2)
  pr_value2 <- predict(lmod2, datapr, weights = (1 - pr_z2))
  true_value1 <- datapr$response1
  true_value2 <- datapr$response2
  pr_value1*round(pr_z2)
  true_value1*round(pr_z2)
  pr_value2*round(1 - pr_z2)
  true_value2*round(1 - pr_z2)
  
# radial_dist <- data1$rho1[(n2+1):n1]
# latlong <- cbind(data1$lat[(n2+1):n1], data1$long[(n2+1):n1])
  
  #plot(datapr$response, predict(lmod1, datapr))
  #abline(a=0,b=1)
  return(list(lmod1 = lmod1, lmod2 = lmod2, pr_value1 = pr_value1,  
              pr_value2 = pr_value2, true_value1 = true_value1, 
              true_value2 = true_value2))
}


