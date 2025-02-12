library(tidyverse)
library(geosphere)
library(housingData)
library(circular)
library(stringr)
library(lubridate)
options(stringsAsFactors = F)

## DATA PREPROCESSING
########################

## define functions for string parsing

char_remove <- function(y){y %>%
    str_replace_all('[[\\D]-[., -]]', '') %>%
    str_split('-', simplify = T) %>%
    apply(1, function(x){max(as.numeric(x), na.rm = T)}) %>%
    na_if(-Inf)}

planting_types <- c('ommercial', 'arden', 'entinel', 'esearch')
planting_fun1 <- function(x){
  str_detect(x, planting_types) %>% 
    which() %>%
    c(0) %>%
    max()}
planting_fun2 <- function(x){
  lapply(1:length(x), function(j){planting_fun1(x[j])}) %>%
    Reduce(f = 'c') %>%
    max()
}
planting_fun3 <- function(x){
  x %>% 
    str_split(' ') %>%
    lapply(planting_fun2) %>%
    Reduce(f = 'c')
}

host_types <- c('ucumber', 'antaloupe', 'umpkin', 'quash', 'ermelon')
host_fun1 <- function(x){
  str_detect(x, host_types) %>% 
    which() %>%
    c(0)}
host_fun2 <- function(x){
  hosts <- lapply(1:length(x), function(j){host_fun1(x[j])}) %>%
    Reduce(f = 'c') 
  if(sum(hosts > 0) > 1){
    out <- 6
  }else{
    out <- max(hosts)
  }
  return(out)
}
host_fun3 <- function(x){
  x %>% 
    str_extract_all('[a-z]+') %>%
    lapply(host_fun2) %>%
    Reduce(f = 'c')
}

## read in raw data files and rename/define variables of interest

raw_data09 <- read.csv('cdm/2009raw.csv')
sub_data09 <- raw_data09 %>%
  transmute(report_date = mdy(Report_date),
            symptom_date = mdy(Symptom_date),
            reporter = Reporter,
            state = State,
            county = County,
            planting = planting_fun3(Planting_type),
            host = host_fun3(Host),
            field_size = NA,
            size_unit = NA,
            leaf_area = char_remove(X.Leaf_area),
            incidence = char_remove(X.Incidence),
            occurrence = (incidence != 0),
            lat = Latitude,
            long = Longitude) # variables of interest

raw_data10 <- read.csv('cdm/2010raw.csv')
sub_data10 <- raw_data10 %>%
  transmute(report_date = mdy(Report_date),
            symptom_date = mdy(Symptom_date),
            reporter = Reporter,
            state = State,
            county = County,
            planting = planting_fun3(Planting_type),
            host = host_fun3(Host),
            field_size = NA,
            size_unit = NA,
            leaf_area = char_remove(X.Leaf_area),
            incidence = char_remove(X.Incidence),
            occurrence = (incidence != 0),
            lat = Latitude,
            long = Longitude) # variables of interest

## merge years

data <- rbind(sub_data09,
              sub_data10)

## assign each unique location an ID, define planting and host factors

locations <- distinct(data, lat, long) %>%
  mutate(locID = 1:(nrow(distinct(data, lat, long))))

data <- merge(data, locations, by = c('lat', 'long')) %>%
  mutate(planting = factor(planting, 
                           labels = c('other', 
                                      'commercial',
                                      'garden',
                                      'sentinel',
                                      'research')),
         host = factor(host, 
                       labels = c('other', 
                                  'cucumber',
                                  'cantaloupe',
                                  'pumpkin',
                                  'squash',
                                  'watermelon',
                                  'multiple')))

## write

save(list = 'data', file = 'cdm_data09-10.RData')
# write.csv(data, file = 'cdm_data09.csv')

## DATA POSTPROCESSING
##########################

## subset data
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

wind09 <- read.csv('Raw_Zedx_data/2009_data/2009_FL_Miami-Dade_wxdata.csv') %>%
  mutate(date = ymd(DATE),
         direction = (WNDDIR + 180)*pi/180)
wind10 <- read.csv('Raw_Zedx_data/2010_data/2010_FL_Miami-Dade_wxdata.csv') %>%
  mutate(date = ymd(DATE),
         direction = (WNDDIR + 180)*pi/180)

## coerce dates and convert dirction to radians

wind_data <- rbind(wind09, 
                   wind10) %>%
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

## ANALYSIS USING SINGLE-SOURCE MODEL
#######################################

## choose year

yr <- 2009

## semiparametric density estimates of wind direction distribution (g(phi))

df1 <- filter(post_data_sub, year == yr)
df2 <- filter(wind_data, year == yr)
z <-  df1$direction.circ
x <- df2$direction
density_est <- density.circular(x = x, 
                                z = z, 
                                bw = 50)

sub_data <- df1 %>%
  bind_cols(g.theta = density_est$y) %>%
  rename(theta = direction.circ,
         rho = distance)

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

## compute variables for input to lm() for parameter estimation
sub_data <- mutate(sub_data, 
                   month = month(report_date),
                   day = day(report_date),
                   time = yday(report_date),
                   response = log(1 + rho/g.theta),
                   sinb1 = sin(as.numeric(theta)),
                   sinb2 = sin(as.numeric(theta) + pi/4))

## split data into training and test sets

set.seed(20520)
n <- nrow(sub_data)
training_prop <- 0.8
n1 <- n
n2 <- floor(training_prop*n)
train_idx <- sample(1:n, n2)

## estimate model parameters

lmod <- lm(response ~ time + sinb1 + sinb2, 
           data = sub_data[train_idx, ])
beta <- coef(lmod)

## compute estimated times of occurrence
sub_data <- mutate(sub_data,
                   that = (response - beta[1] - beta[3]*sinb1 - beta[4]*sinb2)/beta[2])

## plot estimated times versus observed times on training set

plot(sub_data[train_idx, ]$time, sub_data[train_idx, ]$that,
     xlab = 'Day of year of actual occurrence', ylab = 'Predicted day of year',
     main = paste(yr, ' training data'))
abline(b = 1, a = 0)

## plot estimated times versus observed times on test set

plot(sub_data[-train_idx, ]$time, sub_data[-train_idx, ]$that,
     xlab = 'Day of year of actual occurrence', ylab = 'Predicted day of year',
     main = paste(yr, ' test data'))
abline(b = 1, a = 0)

## plot contour shapes

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

