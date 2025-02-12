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
         year(report_date) == 2009) %>%
  select(lat, long, report_date, state, county, locID) %>%
  mutate(year = year(report_date)) %>%
  arrange(report_date)

## identify first report location

first_report <- sub_data %>% 
  top_n(1, desc(report_date)) %>%
  ungroup()

## sample size

n_reports <- sub_data %>%
  group_by(year) %>%
  count()

## match first report with county centroids

county_centroids <- geoCounty %>%
  select(lon, lat) %>%
  as.matrix()

first_report_centroid <- Reduce(rbind, lapply(2009, function(yr){
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

wind09 <- read.csv('Raw_Zedx_data/2009_data/2009_FL_Collier_wxdata.csv')

## coerce dates and convert dirction to radians

wind_data <- mutate(wind09,
                    date = ymd(DATE),
                    direction = (WNDDIR + 180)*pi/180) %>%
  filter(month(date) < 8,
         month(date) > 3) %>%
  select(date, direction) %>%
  group_by(year = year(date))

## align axes between wind data and report data via transformation of bearings

bearing.trans <- function(x){
  (x > 0)*x + (x < 0)*(360 + x)
}

post_data_sub <- filter(post_data, year %in% c(2009, 2010, 2011, 2013, 2015)) %>%
  mutate(direction.circ = as.circular(bearing.trans(bearing)*pi/180,
                                      zero = 0,
                                      units = 'radians',
                                      rotation = 'counter')) %>%
  as_tibble()

## semiparametric density estimate of wind direction distribution (g(phi))

z <- post_data_sub$direction.circ
x <- wind_data$direction
density_est <- density.circular(x = x, z = z, bw = 50)
out <- density_est$y

DATA <- post_data_sub %>%
  bind_cols(g.theta = out) %>%
  rename(theta = direction.circ,
         rho = distance)

## plot of density estimate

theta.grid <- seq(from = 0, to = 2*pi, length = 100)
g.theta <- density.circular(x = x, z = theta.grid, bw = 20)
plot(theta.grid, g.theta$y)

rm(list = setdiff(ls(), 'DATA'))

##

n <- nrow(DATA)
nb <- 1
training_prop <- 0.9
DATA <- mutate(DATA, 
       month = month(report_date),
       day = day(report_date),
       time = yday(report_date),
       response = log(1 + rho/g.theta),
       sinb1 = sin(as.numeric(theta)),
       sinb2 = sin(as.numeric(theta) + pi/4))

n1 <- n
n2 <- floor(training_prop*n)
train_idx <- sort(DATA$lat, index.return = T)$ix[1:n2]
#train_idx <- sample(1:n, n2)
#train_idx <- 1:n2
lmod <- lm(response ~ time + sinb1 + sinb2, data = DATA[train_idx, ])
beta <- coef(lmod)

DATA <- mutate(DATA,
       that = (response - beta[1] - beta[3]*sinb1 - beta[4]*sinb2)/beta[2])

layout(t(1:2))
plot(DATA[train_idx, ]$time, DATA[train_idx, ]$that,
     xlab = 'Day of year of actual occurrence', ylab = 'Predicted day of year',
     main = 'Training data')
abline(b = 1, a = 0)

plot(DATA[-train_idx, ]$time, DATA[-train_idx, ]$that,
     xlab = 'Day of year of actual occurrence', ylab = 'Predicted day of year',
     main = 'Test data')
abline(b = 1, a = 0)

