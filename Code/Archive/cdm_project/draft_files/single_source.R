library(tidyverse)
library(lubridate)
library(geosphere)
library(housingData)
library(circular)
load('cdm_data.RData')

## load and subset data
data <- as_tibble(data)
sub_data <- data %>%
  filter(planting == 'sentinel',
         month(report_date) < 8,
         !is.na(lat), !is.na(long)) %>%
  select(lat, long, report_date, state, county, locID) %>%
  mutate(year = year(report_date)) %>%
  arrange(report_date)

## identify first report location

first_report <- sub_data %>% 
  group_by(year) %>%
  top_n(1, desc(report_date)) %>%
  top_n(1, desc(lat)) %>%
  ungroup()

## sample size

n_reports <- sub_data %>%
  group_by(year) %>%
  count()

## match first report with county centroids

county_centroids <- geoCounty %>%
  select(lon, lat) %>%
  as.matrix()

first_report_centroid <- Reduce(rbind, 
                                lapply(unique(sub_data$year), 
                                       function(yr){
  occurrence <- filter(first_report, year == yr) %>%
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

cbind(first_report$county,
      as.character(first_report_centroid$county))

## append centroids to report information

first_report <- select(first_report_centroid, lon, lat, year) %>%
  rename(cc.lat = lat,
         cc.long = lon) %>%
  merge(first_report, by = 'year')

## merge initial occurrence centroid with report data by year

sub_data <- select(first_report, cc.long, cc.lat, year) %>%
  merge(sub_data, by = 'year') %>%
  as_tibble()

## calculate bearing, distance between initial occurrence centroid and report location

post_data <- cbind(sub_data, bearing = bearing(as.matrix(select(sub_data, cc.long, cc.lat)),
                                               as.matrix(select(sub_data, long, lat))),
                   distance = distGeo(as.matrix(select(sub_data, cc.long, cc.lat)),
                                      as.matrix(select(sub_data, long, lat)))/1000) %>%
  as_tibble()

## import wind data

wind08 <- read.csv('Raw_Zedx_data/2008_data/2008_FL_Hillsborough_south_wxdata.csv') %>%
  mutate(date = ymd(DATE),
         direction = (WNDDIR + 180)*pi/180)
wind09 <- read.csv('Raw_Zedx_data/2009_data/2009_FL_Miami-Dade_wxdata.csv') %>%
  mutate(date = ymd(DATE),
         direction = (WNDDIR + 180)*pi/180)
wind10 <- read.csv('Raw_Zedx_data/2010_data/2010_FL_Miami-Dade_wxdata.csv') %>%
  mutate(date = ymd(DATE),
         direction = (WNDDIR + 180)*pi/180)
wind11 <- read.csv('Raw_Zedx_data/2011_data/2011_FL_Miami-Dade_wxdata.csv') %>%
  mutate(date = ymd(DATE),
         direction = (WNDDIR + 180)*pi/180)
# wind12 <- read.csv('Raw_Zedx_data/2012_data/2012_FL_Miami-Dade_wxdata.csv')
# wind13 <- read.csv('Raw_Zedx_data/2013_data/2013_NC_Miami-Dade_wxdata.csv')
wind14 <- read.csv('Raw_Zedx_data/2014_data/2014_FL_Miami-Dade_wxdata.csv') %>%
  mutate(date = ymd(DATE),
         direction = (WNDDIR + 180)*pi/180)
# wind15 <- read.csv('Raw_Zedx_data/2015_data/2015_SC_Miami-Dade_wxdata.csv')
# wind16 <- read.csv('Raw_Zedx_data/2016_data/2016_FL_Miami-Dade_wxdata.csv')

## coerce dates and convert dirction to radians

wind_data <- rbind(wind08, 
                   wind09, 
                   wind10, 
                   wind11, 
                   wind14) %>%
  filter(month(date) < 8,
         month(date) > 3) %>%
  select(date, direction) %>%
  group_by(year = year(date))

## align axes between wind data and report data via transformation of bearings

bearing.trans <- function(x){
  (x > 0)*x + (x < 0)*(360 + x)
}

post_data_sub <- filter(post_data, year %in% unique(wind_data$year)) %>%
  mutate(direction.circ = as.circular(bearing.trans(bearing)*pi/180,
                                      zero = 0,
                                      units = 'radians',
                                      rotation = 'counter')) %>%
  as_tibble()

## semiparametric density estimate of wind direction distribution (g(phi))

DATA <- Reduce(rbind, lapply(unique(wind_data$year), function(yr){
df1 <- filter(post_data_sub, year == yr)
df2 <- filter(wind_data, year == yr)
z <-  df1$direction.circ
x <- df2$direction
density_est <- density.circular(x = x, 
                                z = z, 
                                bw = 50)

out <- df1 %>%
  bind_cols(g.theta = density_est$y) %>%
  rename(theta = direction.circ,
         rho = distance)
}))

## plot of density estimate

length(unique(wind_data$year))
layout(matrix(1:6, nrow = 2, byrow = T))
lapply(unique(wind_data$year), function(yr){
x <- filter(wind_data, year == yr)
theta.grid <- seq(from = 0, to = 2*pi, length = 100)
g.theta <- density.circular(x = x$direction, 
                            z = theta.grid, 
                            bw = 20)
plot(x = theta.grid, 
     y = g.theta$y,
     xlab = expression(theta),
     ylab = expression(paste('g(', theta, ')')),
     type = 'l',
     xaxt = 'n',
     main = yr)
axis(1, 
     at = seq(from = 0, to = 2*pi, length = 5), 
     labels = c('0', expression(pi/2), 
                expression(pi), 
                expression(3*pi/2),
                expression(2*pi)))
})

rm(list = setdiff(ls(), c('DATA', 'wind_data')))

## parameter estimates for each year and time residual plots

set.seed(20520)
layout(matrix(1:6, nrow = 2, byrow = T))

Reduce(rbind, lapply(unique(DATA$year), function(yr){
sub_data <- filter(DATA, year == yr)
n <- nrow(sub_data)
nb <- 1
training_prop <- 1
sub_data <- mutate(sub_data, 
               month = month(report_date),
               day = day(report_date),
               time = yday(report_date),
               response = log(1 + rho/g.theta),
               sinb1 = sin(as.numeric(theta)),
               sinb2 = sin(as.numeric(theta) + pi/4))

n1 <- n
n2 <- floor(training_prop*n)
#train_idx <- sort(sub_data$lat, index.return = T)$ix[1:n2]
train_idx <- sample(1:n, n2)
#train_idx <- 1:n2
lmod <- lm(response ~ time + sinb1 + sinb2, 
           data = sub_data[train_idx, ])
beta <- coef(lmod)
out <- data.frame(year = yr, t(beta))

sub_data <- mutate(sub_data,
                   that = (response - beta[1] - beta[3]*sinb1 - beta[4]*sinb2)/beta[2])

# plot(sub_data[train_idx, ]$time, sub_data[train_idx, ]$that,
#      xlab = 'Day of year of actual occurrence', ylab = 'Predicted day of year',
#      main = paste(yr, ' training data'))
# abline(b = 1, a = 0)

# plot(sub_data[-train_idx, ]$time, sub_data[-train_idx, ]$that,
#      xlab = 'Day of year of actual occurrence', ylab = 'Predicted day of year',
#      main = paste(yr, ' test data'))
# abline(b = 1, a = 0)

theta_grid <- seq(0, 2*pi, length = 100)
pred_df <- expand.grid(theta = theta_grid, 
                       t = c(15, 30, 45))
g.theta.pred <- density.circular(x = filter(wind_data, year == yr)$direction, 
                                 z = pred_df$theta, bw = 20)$y
pred_df <- mutate(pred_df,
       sinb1 = sin(as.numeric(theta)),
       sinb2 = sin(as.numeric(theta) + pi/4))

pred_df <- cbind(pred_df, 
      pred_resp = as.matrix(cbind(1, pred_df[, 2:4])) %*% beta,
      g_theta = g.theta.pred) %>%
  mutate(rho_pred = (exp(pred_resp) - 1)*g_theta) %>%
  select(theta, rho_pred, t) %>%
  mutate(lon = rho_pred*cos(theta),
         lat = rho_pred*sin(theta))

plot(x = pred_df$lat[pred_df$t == 45], 
     y = pred_df$lon[pred_df$t == 45], 
     type = 'l',
     xaxt = 'n',
     yaxt = 'n',
     xlab = '',
     ylab = '', 
     main = paste(yr, 'contours'))
points(x = pred_df$lat[pred_df$t == 30], 
     y = pred_df$lon[pred_df$t == 30], 
     type = 'l')
points(x = pred_df$lat[pred_df$t == 15], 
       y = pred_df$lon[pred_df$t == 15], 
       type = 'l')

return(beta)
}))


