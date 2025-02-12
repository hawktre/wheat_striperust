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
         year(report_date) > 2007) %>%
  select(lat, long, report_date, state, county, locID) %>%
  mutate(year = year(report_date)) %>%
  arrange(report_date)

## separate by year

annual_data <- lapply(2008:2016,
                      function(yr){
                        filter(sub_data, year(report_date) == yr)
                      })

sapply(2008:2016,
       function(yr){
         nrow(filter(sub_data, year(report_date) == yr))
       })

## identify first reports

first_reports <- sub_data %>% 
  filter(year %in% (2008:2015)) %>%
  group_by(year) %>%
  top_n(1, desc(report_date)) %>%
  ungroup()

## counts of reports per year (not many after subsetting...)

n_reports <- sub_data %>%
  group_by(year) %>%
  count()

## match first reports with county centroids

county_centroids <- geoCounty %>%
  select(lon, lat) %>%
  as.matrix()

first_report_centroids <- Reduce(rbind, lapply(2008:2015, function(yr){
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
  rename(cc.lat = lat,
         cc.long = lon) %>%
  merge(first_reports, by = 'year')

## merge initial occurrence centroid with report data by year

sub_data <- select(first_reports, cc.long, cc.lat, year) %>%
  merge(sub_data, by = 'year') %>%
  as_tibble()

## calculate bearing, distance between initial occurrence centroid and report location

post_data <- cbind(sub_data, bearing = bearing(as.matrix(select(sub_data, cc.long, cc.lat)),
        as.matrix(select(sub_data, long, lat))),
      distance = distGeo(as.matrix(select(sub_data, cc.long, cc.lat)),
        as.matrix(select(sub_data, long, lat)))/1000) %>%
  as_tibble()

## import wind data

select(first_reports, year, state, county)

wind09 <- read.csv('Raw_Zedx_data/2009_data/2009_FL_Collier_wxdata.csv')
wind10 <- read.csv('Raw_Zedx_data/2010_data/2010_FL_Miami-Dade_wxdata.csv')
wind11 <- read.csv('Raw_Zedx_data/2011_data/2011_FL_Miami-Dade_wxdata.csv')
# wind12 <- read.csv('Raw_Zedx_data/2012_data/2012_FL_Marion_wxdata.csv')
wind13 <- read.csv('Raw_Zedx_data/2013_data/2013_NC_Johnston_wxdata.csv')
# wind14 <- read.csv('Raw_Zedx_data/2014_data/2014_FL_Marion_wxdata.csv')
wind15 <- read.csv('Raw_Zedx_data/2015_data/2015_SC_Charleston_wxdata.csv')
# wind16 <- read.csv('Raw_Zedx_data/2016_data/2016_FL_Miami-Dade_wxdata.csv')

## coerce dates and convert dirction to radians

wind_data <- rbind(mutate(wind09,
                          date = ymd(DATE),
                          direction = (WNDDIR + 180)*pi/180),
                   mutate(wind10,
                          date = ymd(DATE),
                          direction = (WNDDIR + 180)*pi/180),
                   mutate(wind11,
                          date = ymd(DATE),
                          direction = (WNDDIR + 180)*pi/180),
                   mutate(wind13,
                          date = ymd(DATE),
                          direction = (WNDDIR + 180)*pi/180),
                   mutate(wind15,
                          date = ymd(DATE),
                          direction = (WNDDIR + 180)*pi/180)) %>%
  filter(month(date) < 8,
         month(date) > 3) %>%
  select(date, direction) %>%
  group_by(year = year(date))

# x <- filter(wind_data, year(date) == 2009)$direction
# mle.vonmises(x)

x <- filter(wind_data, year == 2009)$direction
theta.grid <- seq(from = 0, to = 2*pi, length = 100)
g.theta <- density.circular(x = x, z = theta.grid, bw = 20)
plot(g.theta)
plot(theta.grid, g.theta$y)
hist(x)

# vm.mu <- function(x){
#   est <- mle.vonmises(x)
#   return(est$mu)
# }
# 
# vm.kap <- function(x){
#   est <- mle.vonmises(x)
#   return(est$kappa)
# }
# 
# vm_parameters <- summarize(wind_data,
#           mu = vm.mu(direction),
#           kappa = vm.kap(direction))

bearing.trans <- function(x){
  (x > 0)*x + (x < 0)*(360 + x)
}

test.fn <- function(x, y){
  as.numeric(as.circular(bearing.trans(bearing(x, y))*pi/180,
              zero = 0,
              units = "radians",
              rotation = "counter"))*180/pi
}

post_data_sub <- filter(post_data, year %in% c(2009, 2010, 2011, 2013, 2015)) %>%
  mutate(direction.circ = as.circular(bearing.trans(bearing)*pi/180,
                                      zero = 0,
                                      units = 'radians',
                                      rotation = 'counter')) %>%
  # merge(vm_parameters,
  #       by = 'year')  %>%
  as_tibble()

# g.theta <- sapply(1:nrow(post_data_sub),
#        function(row){
#          dvonmises(post_data_sub$direction.circ[row], 
#                    post_data_sub$mu[row], 
#                    post_data_sub$kappa[row])
#        })

yr <- 2009
#lapply(2009:2012, function(yr){
z <- filter(post_data_sub, year == yr)$direction.circ
x <- filter(wind_data, year == yr)$direction
density_est <- density.circular(x = x, z = z, bw = 50)
out <- density_est$y
#return(out)
#})

# DATA <- bind_cols(post_data_sub, g.theta = g.theta) %>%
#   rename(theta = direction.circ,
#          rho = distance) %>%
#   select(-mu, -kappa)
# 
save(DATA, file = 'postprocessed_cdm_data_everything09.RData')


DATA <- filter(post_data_sub, year == yr) %>%
  bind_cols(g.theta = out) %>%
  rename(theta = direction.circ,
         rho = distance)
