library(tidyverse)
library(lubridate)
library(geosphere)
library(housingData)
library(circular)

wk <- 27
## data preprocessing
#######################


## load and subset data
load('cdm_data.RData')
data <- as_tibble(data)
sub_data <- data %>%
  filter(planting == 'sentinel',
         month(report_date) < 8,
         year(report_date) > 2007,
         !is.na(lat), !is.na(long)) %>%
  select(lat, long, report_date, state, county, locID) %>%
  mutate(year = year(report_date)) %>%
  arrange(report_date)

## identify first reports

first_reports <- sub_data %>% 
  group_by(year) %>%
  top_n(1, desc(report_date)) %>%
  top_n(1, desc(lat)) %>%
  ungroup()

## counts of reports per year (not many after subsetting...)

n_reports <- sub_data %>%
  group_by(year) %>%
  count()

## match first reports with county centroids

county_centroids <- geoCounty %>%
  select(lon, lat) %>%
  as.matrix()

first_report_centroids <- Reduce(rbind, lapply(2008:2016, function(yr){
  occurrence <- filter(first_reports, year == yr) %>%
    select(long, lat)
  
  county_ix <- distGeo(occurrence, 
                       county_centroids) %>%
    which.min()
  
  out <- geoCounty[county_ix, ] %>%
    select(lon, lat, county, state) %>%
    mutate(year = yr)
  
  return(out)
}))


## check for accuracy (should match)

cbind(first_reports$county,
      as.character(first_report_centroids$county))


## append centroids to report information

first_reports <- select(first_report_centroids, lon, lat, year) %>%
  rename(cc.lat.first = lat,
         cc.long.first = lon) %>%
  merge(first_reports, by = 'year')

north_reports <- filter(geoCounty, state == "MI", rMapCounty == "bay") %>%
  select(lon, lat) %>%
  rename(cc.lat.north = lat,
         cc.long.north = lon) %>%
  merge(2008:2016) %>%
  rename(year = y)

## merge initial occurrence centroids with report data by year

sub_data <- select(first_reports, cc.long.first, cc.lat.first, year) %>%
  merge(sub_data, by = 'year') %>%
  as_tibble()

sub_data <- select(north_reports, cc.long.north, cc.lat.north, year) %>%
  merge(sub_data, by = 'year') %>%
  as_tibble()


## calculate bearing and distance between initial occurrence 
## centroid and report location

post_data <- cbind(sub_data, 
                   bearing.first = bearing(as.matrix(select(sub_data, cc.long.first, cc.lat.first)),
                                           as.matrix(select(sub_data, long, lat))),
                   distance.first = distGeo(as.matrix(select(sub_data, cc.long.first, cc.lat.first)),
                                            as.matrix(select(sub_data, long, lat)))/1000,
                   bearing.north = bearing(as.matrix(select(sub_data, cc.long.north, cc.lat.north)),
                                           as.matrix(select(sub_data, long, lat))),
                   distance.north = distGeo(as.matrix(select(sub_data, cc.long.north, cc.lat.north)),
                                            as.matrix(select(sub_data, long, lat)))/1000) %>%
  as_tibble()

## import wind data

wind08_1 <- read.csv('Raw_Zedx_data/2008_data/2008_FL_Hillsborough_south_wxdata.csv') %>%
  mutate(date = ymd(DATE),
         direction = (WNDDIR + 180)*pi/180)
wind09_1 <- read.csv('Raw_Zedx_data/2009_data/2009_FL_Miami-Dade_wxdata.csv') %>%
  mutate(date = ymd(DATE),
         direction = (WNDDIR + 180)*pi/180)
wind10_1 <- read.csv('Raw_Zedx_data/2010_data/2010_FL_Miami-Dade_wxdata.csv') %>%
  mutate(date = ymd(DATE),
         direction = (WNDDIR + 180)*pi/180)
wind11_1 <- read.csv('Raw_Zedx_data/2011_data/2011_FL_Miami-Dade_wxdata.csv') %>%
  mutate(date = ymd(DATE),
         direction = (WNDDIR + 180)*pi/180)
# wind12_1 <- read.csv('Raw_Zedx_data/2012_data/2012_FL_Miami-Dade_wxdata.csv')
# wind13_1 <- read.csv('Raw_Zedx_data/2013_data/2013_NC_Miami-Dade_wxdata.csv')
wind14_1 <- read.csv('Raw_Zedx_data/2014_data/2014_FL_Miami-Dade_wxdata.csv') %>%
  mutate(date = ymd(DATE),
         direction = (WNDDIR + 180)*pi/180)
# wind15_1 <- read.csv('Raw_Zedx_data/2015_data/2015_SC_Miami-Dade_wxdata.csv')
# wind16_1 <- read.csv('Raw_Zedx_data/2016_data/2016_FL_Miami-Dade_wxdata.csv')


wind08MI <- read.csv('Raw_Zedx_data/2008_data/2008_MI_Bay_1_wxdata.csv') %>%
  mutate(date = ymd(DATE),
         direction = (WNDDIR + 180)*pi/180)
wind09MI <- read.csv('Raw_Zedx_data/2009_data/2009_MI_Bay_wxdata.csv') %>%
  mutate(date = ymd(DATE),
         direction = (WNDDIR + 180)*pi/180)
wind10MI <- read.csv('Raw_Zedx_data/2010_data/2010_MI_Bay_wxdata.csv') %>%
  mutate(date = ymd(DATE),
         direction = (WNDDIR + 180)*pi/180)
wind11MI <- read.csv('Raw_Zedx_data/2011_data/2011_MI_Bay_wxdata.csv') %>%
  mutate(date = ymd(DATE),
         direction = (WNDDIR + 180)*pi/180)
# wind12MI <- read.csv('Raw_Zedx_data/2012_data/2012_MI_Bay_wxdata.csv') %>%
#   mutate(date = ymd(DATE),
#          direction = (WNDDIR + 180)*pi/180)
wind13MI <- read.csv('Raw_Zedx_data/2013_data/2013_MI_Bay_wxdata.csv') %>%
  mutate(date = ymd(DATE),
         direction = (WNDDIR + 180)*pi/180)
# wind14MI <- read.csv('Raw_Zedx_data/2014_data/2014_MI_Bay_wxdata.csv') %>%
#   mutate(date = ymd(DATE),
#          direction = (WNDDIR + 180)*pi/180)
wind15MI <- read.csv('Raw_Zedx_data/2015_data/2015_MI_Bay_wxdata.csv') %>%
  mutate(date = ymd(DATE),
         direction = (WNDDIR + 180)*pi/180)
wind16MI <- read.csv('Raw_Zedx_data/2016_data/2016_MI_Bay_wxdata.csv') %>%
  mutate(date = ymd(DATE),
         direction = (WNDDIR + 180)*pi/180)

windMI <- rbind(wind08MI, 
                wind09MI, 
                wind10MI, 
                wind11MI, 
                wind13MI,
                wind15MI,
                wind16MI) %>%
  filter(month(date) < 8,
         week(date) >= wk) %>%
  select(date, direction) %>%
  group_by(year = year(date))

wind08NY <- read.csv('Raw_Zedx_data/2008_data/2008_NY_Niagara_wxdata.csv') %>%
  mutate(date = ymd(DATE),
         direction = (WNDDIR + 180)*pi/180)
wind09NY <- read.csv('Raw_Zedx_data/2009_data/2009_NY_Niagara_wxdata.csv') %>%
  mutate(date = ymd(DATE),
         direction = (WNDDIR + 180)*pi/180)
wind10NY <- read.csv('Raw_Zedx_data/2010_data/2010_NY_Niagara_wxdata.csv') %>%
  mutate(date = ymd(DATE),
         direction = (WNDDIR + 180)*pi/180)
wind11NY <- read.csv('Raw_Zedx_data/2011_data/2011_NY_Niagara_wxdata.csv') %>%
  mutate(date = ymd(DATE),
         direction = (WNDDIR + 180)*pi/180)
# wind12NY <- read.csv('Raw_Zedx_data/2012_data/2012_NY_Niagara_wxdata.csv') %>%
#   mutate(date = ymd(DATE),
#          direction = (WNDDIR + 180)*pi/180)
# wind13NY <- read.csv('Raw_Zedx_data/2013_data/2013_NY_Niagara_wxdata.csv') %>%
#   mutate(date = ymd(DATE),
#          direction = (WNDDIR + 180)*pi/180)
wind14NY <- read.csv('Raw_Zedx_data/2014_data/2014_NY_Niagara_wxdata.csv') %>%
  mutate(date = ymd(DATE),
         direction = (WNDDIR + 180)*pi/180)
wind15NY <- read.csv('Raw_Zedx_data/2015_data/2015_NY_Niagara_wxdata.csv') %>%
  mutate(date = ymd(DATE),
         direction = (WNDDIR + 180)*pi/180)
# wind16NY <- read.csv('Raw_Zedx_data/2016_data/2016_NY_Niagara_wxdata.csv') %>%
#   mutate(date = ymd(DATE),
#          direction = (WNDDIR + 180)*pi/180)


windNY <- rbind(wind08NY, 
                wind09NY, 
                wind10NY, 
                wind11NY, 
                wind14NY,
                wind15NY) %>%
  filter(month(date) < 8,
         week(date) >= wk) %>%
  select(date, direction) %>%
  group_by(year = year(date))

## for wind data, filter dates to between april 1 and july 30, and 
## convert directions to radians

wind_data_1 <- rbind(wind08_1, 
                     wind09_1, 
                     wind10_1, 
                     wind11_1, 
                     wind14_1) %>%
  filter(month(date) < 8,
         month(date) > 3) %>%
  select(date, direction) %>%
  group_by(year = year(date))

wind_data_2 <- windMI
# wind_data_2 <- windNY

## align axes between wind data and report data via transformation of bearings

bearing.trans <- function(x){
  (x > 0)*x + (x < 0)*(360 + x)
}


post_data_sub <- filter(post_data, year %in% intersect(unique(wind_data_1$year), unique(wind_data_2$year))) %>%
  mutate(direction.circ1 = as.circular(bearing.trans(bearing.first)*pi/180,
                                       zero = 0,
                                       units = 'radians',
                                       rotation = 'counter'),
         direction.circ2 = as.circular(bearing.trans(bearing.north)*pi/180,
                                       zero = 0,
                                       units = 'radians',
                                       rotation = 'counter')) %>%
  as_tibble()


DATA <- Reduce(rbind, lapply(unique(post_data_sub$year), function(yr){
  
  # disease data
  df1 <- filter(post_data_sub, year == yr)
  
  # wind data at southern source
  df2a <- filter(wind_data_1, year == yr)
  
  # wind data at northern source
  df2b <- filter(wind_data_2, year == yr)
  
  # kernel density estimate for wind directions at southern source
  z <-  df1$direction.circ1
  x <- df2a$direction
  density_est1 <- density.circular(x = x, 
                                   z = z, 
                                   bw = 50)
  
  # kernel density estimate for wind directions at northern source
  z <-  df1$direction.circ2
  x <- df2b$direction
  density_est2 <- density.circular(x = x, 
                                   z = z, 
                                   bw = 50)
  
  # append kernel density estimates to disease data
  out <- df1 %>%
    bind_cols(g.theta.1 = density_est1$y,
              g.theta.2 = density_est2$y) %>%
    rename(theta1 = direction.circ1,
           rho1 = distance.first,
           theta2 = direction.circ2,
           rho2 = distance.north) %>%
    select(report_date, lat, long, rho1, theta1, g.theta.1,
           rho2, theta2, g.theta.2)
  
  return(out)
}))

# replace northern source data by NA before 7/1
DATA <- DATA %>%
  mutate(rho2 = na_if(rho2*(week(report_date) >= wk), 0),
         theta2 = na_if(theta2*(week(report_date) >= wk), 0),
         g.theta.2 = na_if(g.theta.2*(week(report_date) >= wk), 0))

save(DATA, file = 'cdm_multisource_test_data.RData')
rm(list = setdiff(ls(), 'DATA'))
#######################

nFolds <- 50
out <- Reduce(rbind, lapply(1:nFolds, function(fold){
# select year
yr <- 2009

## model fitting
######################
load('cdm_multisource_test_data.RData')

# which years available?
unique(year(DATA$report_date))

# subset data and calculate basis expansions and response variables
sub_data <- filter(DATA, year(report_date) == yr)
n <- nrow(sub_data)
sub_data <- mutate(sub_data, 
                   month = month(report_date),
                   day = day(report_date),
                   time = yday(report_date),
                   response.1 = log(1 + rho1/g.theta.1), # response for S source
                   response.2 = log(1 + rho2/g.theta.2), # response for N source
                   sinb1.1 = sin(as.numeric(theta1)), # sin basis 1 for S source
                   sinb2.1 = sin(as.numeric(theta1) + pi/4), # sin basis 2 for S source
                   sinb1.2 = sin(as.numeric(theta2)), # sin basis 1 for N source
                   sinb2.2 = sin(as.numeric(theta2) + pi/4), # sin basis 2 for N source
                   z.reg.1 = rho1/g.theta.1, # regressor 1 to initialize weights
                   z.reg.2 = rho2/g.theta.2, # regressor 2 to initialize weights
                   z.fit = rho1 < rho2) # is southern source closer?

# split data into training and test sets
training_prop <- 0.8
n1 <- n
n2 <- floor(training_prop*n)
# set.seed(32420)
train_idx <- sample(1:n, n2)
# train_idx <- 1:n2
# train_idx <- sort(sub_data$lat, index.return = T)$ix[1:n2]

## TWO SOURCE MODEL

# initialize weights: P(occurrence due to S source) 
# by regressing 1{S closer} on exp{response.1} and exp{response.2} 

train_data <- sub_data[train_idx, ]

const_idx <- which(is.na(train_data$z.reg.2))
update_idx <- which(!is.na(train_data$z.reg.2))

z_glm <- glm(z.fit ~ z.reg.1 + z.reg.2, 
             data = train_data[update_idx, ], 
             family = binomial())
pr_z <- rep(0, length(train_idx))
pr_z[const_idx] <- 1
pr_z[update_idx] <- fitted(z_glm)


# fit model via EM
err <- 1
iter <- 1
while(err > 0.001) {
  # fit model 1 (S source)
  lmod1 <- lm(response.1 ~ time + sinb1.1 + sinb2.1, 
              data = train_data,
              weights = pr_z)
  
  
  # fit model 2 (N source)
  lmod2 <- lm(response.2 ~ time + sinb1.2 + sinb2.2, 
              data = train_data[update_idx, ],
              weights = (1 - pr_z[update_idx]))
  
  # recompute weights
  mean1 <- residuals(lmod1)/summary(lmod1)$sigma 
  mean2 <- residuals(lmod2)/summary(lmod2)$sigma
  pr_z1 <- pr_z
  pr_z1[update_idx] <- (dnorm(mean1[update_idx])*pr_z[update_idx])/(dnorm(mean1[update_idx])*pr_z[update_idx] + dnorm(mean2)*(1 - pr_z[update_idx]))
  
  # check for convergence
  err <- mean(abs(pr_z - pr_z1))
  
  # update weights
  pr_z <- pr_z1
  
  # increment iteration count
  iter <- iter + 1
}

# inspect estimates
summary(lmod1)
summary(lmod2)
beta1 <- coef(lmod1)
beta2 <- coef(lmod2)
round(pr_z, 6)

## ONE SOURCE MODEL

# fit model
lmod <- lm(response.1 ~ time + sinb1.1 + sinb2.1, 
           data = train_data)

# inspect estimates
summary(lmod)
beta <- coef(lmod)

######################


## predictions
######################

const_idx <- which(is.na(sub_data$z.reg.2))
update_idx <- which(!is.na(sub_data$z.reg.2))

z_glm <- glm(z.fit ~ z.reg.1 + z.reg.2, 
             data = sub_data[update_idx, ], 
             family = binomial())
pr_z <- rep(0, n)
pr_z[const_idx] <- 1
pr_z[update_idx] <- fitted(z_glm)
which_pred <- round(pr_z, 0) + (round(pr_z, 0) == 0)*2

# predictions on held-out data
sub_data <- mutate(sub_data, 
       t.hat.S = (response.1 - beta1[1] - beta1[3]*sinb1.1 - beta1[4]*sinb2.1)/beta1[2],
       t.hat.N = (response.2 - beta2[1] - beta2[3]*sinb1.2 - beta2[4]*sinb2.2)/beta2[2],
       pred.S = beta1[1] + beta1[2]*time + beta1[3]*sinb1.1 + beta1[4]*sinb2.1,
       pred.N = beta2[1] + beta2[2]*time + beta2[3]*sinb1.2 + beta2[4]*sinb2.2,
       pred.onesource = beta[1] + beta[2]*time + beta[3]*sinb1.1 + beta[4]*sinb2.1,
       t.hat.onesource = (response.1 - beta[1] - beta[3]*sinb1.1 - beta[4]*sinb2.1)/beta[2]) %>%
  mutate(t.hat.twosource.min = if_else(!is.na(t.hat.N), pmin(t.hat.S, t.hat.N), t.hat.S),
         which.pred.twosource.min = if_else(is.na(t.hat.N), 1, 2*as.numeric(t.hat.S > t.hat.N) + 1*as.numeric(t.hat.S < t.hat.N))) %>%
  cbind(which.pred.twosource.pr = which_pred) %>%
  mutate(t.hat.twosource.pr = if_else(which.pred.twosource.pr == 1, t.hat.S, t.hat.N))

select(sub_data, which.pred.twosource.pr, which.pred.twosource.min)


# plots for two-source model
layout(matrix(1:4, nrow = 2, byrow = T))
plot(sub_data[train_idx, ]$time, sub_data[train_idx, ]$t.hat.twosource.pr,
     xlab = 'Day of year of actual occurrence', ylab = 'Predicted day of year',
     main = paste('Two-source model fit on training data ', yr),
     col = sub_data$which.pred.twosource.pr,
     xlim = c(100, 220),
     ylim = c(50, 300))
abline(b = 1, a = 0)
abline(b = 1, a = 21, lty = 2, col = 4)
abline(b = 1, a = -21, lty = 2, col = 4)


plot(sub_data[-train_idx, ]$time, sub_data[-train_idx, ]$t.hat.twosource.pr,
     xlab = 'Day of year of actual occurrence', ylab = 'Predicted day of year',
     main = paste('Two-source model fit on test data ', yr),
     col = sub_data$which.pred.twosource.pr[-train_idx],
     xlim = c(100, 220),
     ylim = c(50, 300))
abline(b = 1, a = 0)
abline(b = 1, a = 21, lty = 2, col = 4)
abline(b = 1, a = -21, lty = 2, col = 4)

plot(sub_data[train_idx, ]$time, sub_data[train_idx, ]$t.hat.onesource,
     xlab = 'Day of year of actual occurrence', ylab = 'Predicted day of year',
     main = paste('One-source model fit on training data ', yr),
     xlim = c(100, 220),
     ylim = c(50, 300))
abline(b = 1, a = 0)
abline(b = 1, a = 21, lty = 2, col = 4)
abline(b = 1, a = -21, lty = 2, col = 4)

plot(sub_data[-train_idx, ]$time, sub_data[-train_idx, ]$t.hat.onesource,
     xlab = 'Day of year of actual occurrence', ylab = 'Predicted day of year',
     main = paste('One-source model fit on test data ', yr),
     xlim = c(100, 220),
     ylim = c(50, 300))
abline(b = 1, a = 0)
abline(b = 1, a = 21, lty = 2, col = 4)
abline(b = 1, a = -21, lty = 2, col = 4)



# comparison of predictions
# plot(sub_data$t.hat.twosource.pr[-train_idx], 
#      sub_data$t.hat.onesource[-train_idx],
#      xlab = 'Two-source prediction (day of year)',
#      ylab = 'One-source prediction (day of year)',
#      main = 'Prediction vs. prediction (2008)')
# abline(a = 0, b = 1)

# mean error
transmute(sub_data[-train_idx, ], 
          err.1s = abs(t.hat.onesource - time),
          err.2s = abs(t.hat.twosource.pr - time)) %>%
  apply(2, mean)


}))

apply(out, 2, mean)


