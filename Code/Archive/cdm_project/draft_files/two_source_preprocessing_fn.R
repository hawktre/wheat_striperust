library(tidyverse)
library(lubridate)
library(geosphere)
library(housingData)
library(circular)

preprocess_fn <- function(wk, MI){
# wk <- 26 # week of year at which N source appears
# MI <- 1 # N source location (1 = MI, 0 = NY)

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
  select(lat, long, report_date, state, county, locID, host) %>%
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

if(MI == 1){
  wind_data_2 <- windMI
}else{
  wind_data_2 <- windNY
}

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
           rho2, theta2, g.theta.2, host)
  
  return(out)
}))

# replace northern source data by NA before wk
DATA <- DATA %>%
  mutate(rho2 = na_if(rho2*(week(report_date) >= wk), 0),
         theta2 = na_if(theta2*(week(report_date) >= wk), 0),
         g.theta.2 = na_if(g.theta.2*(week(report_date) >= wk), 0))

return(DATA)
#######################
}