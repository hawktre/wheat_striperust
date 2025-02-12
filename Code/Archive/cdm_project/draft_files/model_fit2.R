fit_pred <- function(yr, training_prop = 0.64, nb = 1) {

# load("~/Documents/Trevor/postprocessed_cdm_data.rdata")

n <- length(DATA$year)
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
DATA$response <- log(1 + DATA$rho/DATA$g.theta)
sinb <- matrix(0, nrow = n, ncol = 2*nb)
cosb <- matrix(0, nrow = n, ncol = 2*nb)
#forml <- "response ~ time "
for(k in seq(1, 2*nb-1, 2)){
  sinb[, k] <- sin(k*as.numeric(DATA$theta))
  sinb[, (k+1)] <- sin(k*as.numeric(DATA$theta) + pi/4)
  #  sinb[, (k+2)] <- sin(k*as.numeric(DATA$theta) + 3*pi/4)
  cosb[, k] <- cos(k*as.numeric(DATA$theta))
  cosb[, (k+1)] <- cos(k*as.numeric(DATA$theta) + pi/4)
  #  cosb[, (k+2)] <- cos(k*as.numeric(DATA$theta) + pi/2)
  #  forml <- paste(forml, "+ sinb.", k, "+ cosb.", k, sep="")
}
DATA$sinb <- sinb
DATA$cosb <- cosb

# yr = 2009
data1 <- DATA[DATA$year == yr, ]
n1 <- length(data1$year)
n2 <- floor(training_prop*n1)
lmod <- lm(response ~ time + sinb, data = data1[1:n2, ])
# summary(lmod1)
datapr <- data1[(n2+1):n1, ]
pr_value <- predict(lmod, datapr)
true_value <- datapr$response
radial_dist <- data1$rho[(n2+1):n1]
latlong <- cbind(data1$lat[(n2+1):n1], data1$long[(n2+1):n1])

#plot(datapr$response, predict(lmod1, datapr))
#abline(a=0,b=1)
return(list(lmod = lmod, pr_value = pr_value, true_value = true_value,
            radial_dist = radial_dist, latlong = latlong))
}

pred_res <- function(lmod, date, theta, nb = 1) {
  
  n <- length(as.numeric(theta))
  month <- as.numeric(unlist(strsplit(as.character(date),
                                           "-")))[seq(2, 3*n-1, 3)]
  day <- as.numeric(unlist(strsplit(as.character(date),
                                         "-")))[seq(3, 3*n, 3)]
  dm <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
  time <- rep(0, n)
  for(i in 1:n) {
    time[i] <- sum(dm[1:(month[i]-1)]) + day[i]
  }
  sinb <- matrix(0, nrow = n, ncol = 2*nb)
  cosb <- matrix(0, nrow = n, ncol = 2*nb)
  for(k in seq(1, 2*nb-1, 2)){
    sinb[, k] <- sin(k*as.numeric(theta))
    sinb[, (k+1)] <- sin(k*as.numeric(theta) + pi/4)
    cosb[, k] <- cos(k*as.numeric(theta))
    cosb[, (k+1)] <- cos(k*as.numeric(theta) + pi/4)
  }
#  data1 <- data.frame(time = time, sinb = sinb, cosb = cosb)
#  pr_value <- predict(lmod, newdata = data1)
  pr_value <- lmod$coefficients[1] + lmod$coefficients[2]*time +
    lmod$coefficients[3]*sinb[,1] + lmod$coefficients[4]*sinb[,2]
  return(pr_value)
}

fit_pred(2009)
