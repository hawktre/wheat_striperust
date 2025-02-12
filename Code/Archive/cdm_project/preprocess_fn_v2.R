library(tidyverse)
library(lubridate)
library(geosphere)
library(housingData)
library(circular)

preprocess_fn <- function(){
  # 1st northern location
  imputed_state1 = rep('NY', 9)
  imputed_county1 = rep('niagara', 9)
  imputed_date1 = ymd(c('2008/01/01',
                       '2009/01/01',
                       '2010/01/01',
                       '2011/01/01',
                       '2012/01/01',
                       '2013/01/01',
                       '2014/01/01',
                       '2015/01/01',
                       '2016/01/01'))
  week(imputed_date1) <- rep(23, 9)
  
  # 2nd northern location
  imputed_state2 = rep('OH', 9)
  imputed_county2 = rep('sandusky', 9)
  imputed_date2 = ymd(c('2008/01/01',
                        '2009/01/01',
                        '2010/01/01',
                        '2011/01/01',
                        '2012/01/01',
                        '2013/01/01',
                        '2014/01/01',
                        '2015/01/01',
                        '2016/01/01'))
  week(imputed_date2) <- rep(23, 9)
  
  # 1st southwestern location
  imputed_state3 = rep('TX', 9)
  imputed_county3 = rep('brazos', 9)
  imputed_date3 = ymd(c('2008/01/01',
                        '2009/01/01',
                        '2010/01/01',
                        '2011/01/01',
                        '2012/01/01',
                        '2013/01/01',
                        '2014/01/01',
                        '2015/01/01',
                        '2016/01/01'))
  week(imputed_date3) <- rep(19, 9)
  
  # 2nd southwestern location
  imputed_state4 = rep('TX', 9)
  imputed_county4 = rep('hidalgo', 9)
  imputed_date4 = ymd(c('2008/01/01',
                        '2009/01/01',
                        '2010/01/01',
                        '2011/01/01',
                        '2012/01/01',
                        '2013/01/01',
                        '2014/01/01',
                        '2015/01/01',
                        '2016/01/01'))
  week(imputed_date4) <- rep(19, 9)
  
  ## load and subset data
  
  load('cdm_data.RData')
  data <- as_tibble(data)
  sub_data <- data %>%
    filter(planting == 'sentinel',
           # month(report_date) < 8,
           # year(report_date) > 2007,
           month(symptom_date) < 8,
           year(symptom_date) > 2007,
           !is.na(lat), !is.na(long),
           lat < 48) %>%
    rename(report.date = report_date,
           symptom.date = symptom_date) %>%
    select(lat, long, report.date, symptom.date, 
           state, county, locID, host) %>%
    mutate(year = year(report.date)) %>%
    arrange(report.date)
  
  ## assign unique ID to each sentinel plot
  
  IDs <- sub_data %>% distinct(lat, long) %>% mutate(ID = row_number())
  
  ## merge IDs and select first symptom date associated with each plot
  
  sub_data <- merge(sub_data, IDs, by = c('lat', 'long')) %>%
    group_by(year, ID) %>%
    # top_n(-1, report.date) %>%
    top_n(-1, symptom.date) %>%
    ungroup()
  
  # ## check that assigned ID's match entered ID's
  # n_distinct(sub_data$ID)
  # n_distinct(sub_data$locID)
  
  ## identify first reports
  
  first_occurrences_s <- sub_data %>% 
    group_by(year) %>%
    # top_n(1, desc(report.date)) %>%
    top_n(1, desc(symptom.date)) %>%
    top_n(1, desc(lat)) %>%
    ungroup() %>%
    arrange(year)
  
  first_occurrences_w <- sub_data %>%
    filter(long < -90) %>%
    group_by(year) %>%
    # top_n(1, desc(report.date)) %>%
    top_n(1, desc(symptom.date)) %>%
    top_n(-1, desc(lat)) %>%
    ungroup() %>%
    arrange(year)
  
  first_occurrences_n <- sub_data %>%
    filter(lat > 40, !str_starts(state, 'ON')) %>%
    group_by(year) %>%
    top_n(5, lat) %>%
    # top_n(1, desc(report.date)) %>%
    top_n(1, desc(symptom.date)) %>%
    top_n(1, lat) %>%
    ungroup() %>%
    arrange(year) # observed reports
  
  ## counts of reports per year (not many after subsetting...)
  
  n_reports <- sub_data %>%
    group_by(year) %>%
    count()
  
  ## match first reports with county centroids
  
  county_centroids <- geoCounty %>%
    select(lon, lat) %>%
    as.matrix()
  
  first_occurrence_centroids_s <- Reduce(rbind, 
                                         lapply(2008:2016, function(yr){
                                           occurrence <- filter(first_occurrences_s, year == yr) %>%
                                             select(long, lat)
                                           
                                           county_ix <- distGeo(occurrence, 
                                                                county_centroids) %>%
                                             which.min()
                                           
                                           out <- geoCounty[county_ix, ] %>%
                                             select(lon, lat, county, state) %>%
                                             mutate(year = yr)
                                           
                                           return(out)
                                         }))
  
  first_occurrence_centroids_w <- Reduce(rbind, 
                                         lapply(unique(first_occurrences_w$year), 
                                                function(yr){
                                                  occurrence <- filter(first_occurrences_w, year == yr) %>%
                                                    select(long, lat)
                                                  
                                                  county_ix <- distGeo(occurrence, 
                                                                       county_centroids) %>%
                                                    which.min()
                                                  
                                                  out <- geoCounty[county_ix, ] %>%
                                                    select(lon, lat, county, state) %>%
                                                    mutate(year = yr)
                                                  
                                                  return(out)
                                                }))
  
  first_occurrence_centroids_n <- Reduce(rbind, 
                                         lapply(unique(first_occurrences_n$year), 
                                                function(yr){
                                                  occurrence <- filter(first_occurrences_n, year == yr) %>%
                                                    select(long, lat)
                                                  
                                                  county_ix <- distGeo(occurrence, 
                                                                       county_centroids) %>%
                                                    which.min()
                                                  
                                                  out <- geoCounty[county_ix, ] %>%
                                                    select(lon, lat, county, state) %>%
                                                    mutate(year = yr)
                                                  
                                                  return(out)
                                                }))
  
  
  ## compare county of occurrence with closest county centroid
  
  cbind(first_occurrences_s$county,
        as.character(first_occurrence_centroids_s$county))
  
  cbind(first_occurrences_w$county,
        as.character(first_occurrence_centroids_w$county))
  
  cbind(first_occurrences_n$county,
        as.character(first_occurrence_centroids_n$county))
  
  
  ## append centroids to report information
  
  first_occurrences_s <- select(first_occurrence_centroids_s, lon, lat, year) %>%
    rename(cc.lat.s = lat,
           cc.long.s = lon) %>%
    merge(first_occurrences_s, by = 'year') %>%
    merge(data.frame(year = 2008:2016), by = 'year', all = T)
  
  first_occurrences_w <- select(first_occurrence_centroids_w, lon, lat, year) %>%
    rename(cc.lat.w = lat,
           cc.long.w = lon) %>%
    merge(first_occurrences_w, by = 'year') %>%
    merge(data.frame(year = 2008:2016), by = 'year', all = T)
  
  first_occurrences_n <- select(first_occurrence_centroids_n, lon, lat, year) %>%
    rename(cc.lat.n = lat,
           cc.long.n = lon) %>%
    merge(first_occurrences_n, by = 'year') %>%
    merge(data.frame(year = 2008:2016), by = 'year', all = T)
  
  imputed_occurrences1 <- Reduce(rbind, lapply(1:9, function(yr){
    out <- filter(geoCounty, 
                  state == imputed_state1[yr], 
                  rMapCounty == imputed_county1[yr]) %>%
      mutate(year = 2007 + yr) %>%
      rename(cc.lat.i1 = lat,
             cc.long.i1 = lon) %>%
      select(year, cc.long.i1, cc.lat.i1) %>%
      mutate(lat = NA,
             long = NA,
             report.date = imputed_date1[yr],
             symptom.date = imputed_date1[yr],
             state = imputed_state1[yr],
             county = imputed_county1[yr],
             locID = NA,
             host = NA,
             ID = NA)
    return(out)})) %>%
    merge(data.frame(year = 2008:2016), by = 'year', all = T)
  
  imputed_occurrences2 <- Reduce(rbind, lapply(1:9, function(yr){
    out <- filter(geoCounty, 
                  state == imputed_state2[yr], 
                  rMapCounty == imputed_county2[yr]) %>%
      mutate(year = 2007 + yr) %>%
      rename(cc.lat.i2 = lat,
             cc.long.i2 = lon) %>%
      select(year, cc.long.i2, cc.lat.i2) %>%
      mutate(lat = NA,
             long = NA,
             report.date = imputed_date2[yr],
             symptom.date = imputed_date2[yr],
             state = imputed_state2[yr],
             county = imputed_county2[yr],
             locID = NA,
             host = NA,
             ID = NA)
    return(out)})) %>%
    merge(data.frame(year = 2008:2016), by = 'year', all = T)
  
  imputed_occurrences3 <- Reduce(rbind, lapply(1:9, function(yr){
    out <- filter(geoCounty, 
                  state == imputed_state3[yr], 
                  rMapCounty == imputed_county3[yr]) %>%
      mutate(year = 2007 + yr) %>%
      rename(cc.lat.i3 = lat,
             cc.long.i3 = lon) %>%
      select(year, cc.long.i3, cc.lat.i3) %>%
      mutate(lat = NA,
             long = NA,
             report.date = imputed_date3[yr],
             symptom.date = imputed_date3[yr],
             state = imputed_state3[yr],
             county = imputed_county3[yr],
             locID = NA,
             host = NA,
             ID = NA)
    return(out)})) %>%
    merge(data.frame(year = 2008:2016), by = 'year', all = T)
  
  imputed_occurrences4 <- Reduce(rbind, lapply(1:9, function(yr){
    out <- filter(geoCounty, 
                  state == imputed_state4[yr], 
                  rMapCounty == imputed_county4[yr]) %>%
      mutate(year = 2007 + yr) %>%
      rename(cc.lat.i4 = lat,
             cc.long.i4 = lon) %>%
      select(year, cc.long.i4, cc.lat.i4) %>%
      mutate(lat = NA,
             long = NA,
             report.date = imputed_date4[yr],
             symptom.date = imputed_date4[yr],
             state = imputed_state4[yr],
             county = imputed_county4[yr],
             locID = NA,
             host = NA,
             ID = NA)
    return(out)})) %>%
    merge(data.frame(year = 2008:2016), by = 'year', all = T)
  
  
  ## merge initial occurrence centroids with report data by year
  
  sub_data <- select(first_occurrences_s, cc.long.s, cc.lat.s, year) %>%
    merge(sub_data, by = 'year', all = T) %>%
    as_tibble()
  sub_data <- select(first_occurrences_n, cc.long.n, cc.lat.n, year) %>%
    merge(sub_data, by = 'year', all = T) %>%
    as_tibble()
  sub_data <- select(first_occurrences_w, cc.long.w, cc.lat.w, year) %>%
    merge(sub_data, by = 'year', all = T) %>%
    as_tibble()
  sub_data <- select(imputed_occurrences1, cc.long.i1, cc.lat.i1, year) %>%
    merge(sub_data, by = 'year', all = T) %>%
    as_tibble()
  sub_data <- select(imputed_occurrences2, cc.long.i2, cc.lat.i2, year) %>%
    merge(sub_data, by = 'year', all = T) %>%
    as_tibble()
  sub_data <- select(imputed_occurrences3, cc.long.i3, cc.lat.i3, year) %>%
    merge(sub_data, by = 'year', all = T) %>%
    as_tibble()
  sub_data <- select(imputed_occurrences4, cc.long.i4, cc.lat.i4, year) %>%
    merge(sub_data, by = 'year', all = T) %>%
    as_tibble()
  
    
  ## calculate bearing and distance between initial occurrence 
  ## centroid and report location
  
  sub_data <- cbind(sub_data, 
                    bearing.n = bearing(as.matrix(select(sub_data, cc.long.n, cc.lat.n)),
                                        as.matrix(select(sub_data, long, lat))),
                    distance.n = distGeo(as.matrix(select(sub_data, cc.long.n, cc.lat.n)),
                                         as.matrix(select(sub_data, long, lat)))/1000,
                    bearing.w = bearing(as.matrix(select(sub_data, cc.long.w, cc.lat.w)),
                                        as.matrix(select(sub_data, long, lat))),
                    distance.w = distGeo(as.matrix(select(sub_data, cc.long.w, cc.lat.w)),
                                         as.matrix(select(sub_data, long, lat)))/1000,
                    bearing.s = bearing(as.matrix(select(sub_data, cc.long.s, cc.lat.s)),
                                        as.matrix(select(sub_data, long, lat))),
                    distance.s = distGeo(as.matrix(select(sub_data, cc.long.s, cc.lat.s)),
                                         as.matrix(select(sub_data, long, lat)))/1000,
                    bearing.i1 = bearing(as.matrix(select(sub_data, cc.long.i1, cc.lat.i1)),
                                        as.matrix(select(sub_data, long, lat))),
                    distance.i1 = distGeo(as.matrix(select(sub_data, cc.long.i1, cc.lat.i1)),
                                         as.matrix(select(sub_data, long, lat)))/1000,
                    bearing.i2 = bearing(as.matrix(select(sub_data, cc.long.i2, cc.lat.i2)),
                                         as.matrix(select(sub_data, long, lat))),
                    distance.i2 = distGeo(as.matrix(select(sub_data, cc.long.i2, cc.lat.i2)),
                                          as.matrix(select(sub_data, long, lat)))/1000,
                    bearing.i3 = bearing(as.matrix(select(sub_data, cc.long.i3, cc.lat.i3)),
                                         as.matrix(select(sub_data, long, lat))),
                    distance.i3 = distGeo(as.matrix(select(sub_data, cc.long.i3, cc.lat.i3)),
                                          as.matrix(select(sub_data, long, lat)))/1000,
                    bearing.i4 = bearing(as.matrix(select(sub_data, cc.long.i4, cc.lat.i4)),
                                         as.matrix(select(sub_data, long, lat))),
                    distance.i4 = distGeo(as.matrix(select(sub_data, cc.long.i4, cc.lat.i4)),
                                          as.matrix(select(sub_data, long, lat)))/1000) %>%
    as_tibble()
  
  ## transform bearings to angles
  
  bearing.trans <- function(x){
    (x > 0)*x + (x < 0)*(360 + x)
  }
  
  sub_data <- sub_data %>%
    mutate(theta.s = as.circular(bearing.trans(bearing.s)*pi/180,
                                 zero = 0,
                                 units = 'radians',
                                 rotation = 'counter'),
           theta.w = as.circular(bearing.trans(bearing.w)*pi/180,
                                 zero = 0,
                                 units = 'radians',
                                 rotation = 'counter'),
           theta.n = as.circular(bearing.trans(bearing.n)*pi/180,
                                 zero = 0,
                                 units = 'radians',
                                 rotation = 'counter'),
           theta.i1 = as.circular(bearing.trans(bearing.i1)*pi/180,
                                 zero = 0,
                                 units = 'radians',
                                 rotation = 'counter'),
           theta.i2 = as.circular(bearing.trans(bearing.i2)*pi/180,
                                  zero = 0,
                                  units = 'radians',
                                  rotation = 'counter'),
           theta.i3 = as.circular(bearing.trans(bearing.i3)*pi/180,
                                  zero = 0,
                                  units = 'radians',
                                  rotation = 'counter'),
           theta.i4 = as.circular(bearing.trans(bearing.i4)*pi/180,
                                  zero = 0,
                                  units = 'radians',
                                  rotation = 'counter')) %>%
    as_tibble()
  
  ## import wind data
  
  wind08_s <- read.csv('wind_data/2008s.csv') %>%
    mutate(date = ymd(DATE),
           direction = (WNDDIR + 180)*pi/180)
  wind09_s <- read.csv('wind_data/2009s.csv') %>%
    mutate(date = ymd(DATE),
           direction = (WNDDIR + 180)*pi/180)
  wind10_s <- read.csv('wind_data/2010s.csv') %>%
    mutate(date = ymd(DATE),
           direction = (WNDDIR + 180)*pi/180)
  
  wind_data_s <- rbind(filter(select(wind08_s, -DOY), 
                              date >= first_occurrences_s$symptom.date[1]),
                       filter(wind09_s, 
                              date >= first_occurrences_s$symptom.date[2]),
                       filter(wind10_s, 
                              date >= first_occurrences_s$symptom.date[3])) %>%
    filter(month(date) < 8) %>%
    select(date, direction) %>%
    group_by(year = year(date))
  
  
  wind08_w <- read.csv('wind_data/2008w.csv') %>%
    mutate(date = ymd(DATE),
           direction = (WNDDIR + 180)*pi/180)
  wind09_w <- read.csv('wind_data/2009w.csv') 
  wind09_w <- wind09_w %>% cbind(TIME = str_split_fixed(wind09_w$DATE, 'T', 2)) %>% 
    select(WNDDIR, WNDSPD, TIME.1, TIME.2) %>%
    mutate(DATE = ymd(TIME.1), IHR = hour(hms(TIME.2))) %>% 
    select(DATE, IHR, WNDDIR, WNDSPD, TIME.2) %>%
    mutate(date = ymd(DATE),
           direction = (WNDDIR + 180)*pi/180)
  wind10_w <- read.csv('wind_data/2010w.csv') 
  wind10_w <- wind10_w %>% cbind(TIME = str_split_fixed(wind10_w$DATE, 'T', 2)) %>% 
    select(WNDDIR, WNDSPD, TIME.1, TIME.2) %>%
    mutate(DATE = ymd(TIME.1), IHR = hour(hms(TIME.2))) %>% 
    select(DATE, IHR, WNDDIR, WNDSPD, TIME.2) %>%
    mutate(date = ymd(DATE),
           direction = (WNDDIR + 180)*pi/180)
  
  wind_data_w <- rbind(filter(select(wind08_w, date, direction), 
                              date >= first_occurrences_w$symptom.date[1]),
                       filter(select(wind09_w, date, direction), 
                              date >= first_occurrences_w$symptom.date[2]),
                       filter(select(wind10_w, date, direction), 
                              date >= first_occurrences_w$symptom.date[3])) %>%
    filter(month(date) < 8) %>%
    group_by(year = year(date))
  
  
  wind08_n <- read.csv('wind_data/2008n.csv') %>%
    mutate(date = ymd(DATE),
           direction = (WNDDIR + 180)*pi/180)
  wind09_n <- read.csv('wind_data/2009n.csv') %>%
    mutate(date = ymd(DATE),
           direction = (WNDDIR + 180)*pi/180)
  wind10_n <- read.csv('wind_data/2010n.csv') %>%
    mutate(date = ymd(DATE),
           direction = (WNDDIR + 180)*pi/180)

  wind_data_n <- rbind(filter(select(wind08_n, date, direction), 
                              date >= first_occurrences_n$symptom.date[1]),
                       filter(select(wind09_n, date, direction), 
                              date >= first_occurrences_n$symptom.date[2]),
                       filter(select(wind10_n, date, direction), 
                              date >= first_occurrences_n$symptom.date[3])) %>%
    filter(month(date) < 8) %>%
    group_by(year = year(date))
  
  
  
  wind08_i1 <- read.csv('wind_data/2008i1.csv') %>%
    mutate(date = ymd(DATE),
           direction = (WNDDIR + 180)*pi/180)
  wind09_i1 <- read.csv('wind_data/2009i1.csv') %>%
    mutate(date = ymd(DATE),
           direction = (WNDDIR + 180)*pi/180)
  wind10_i1 <- read.csv('wind_data/2010i1.csv') %>%
    mutate(date = ymd(DATE),
           direction = (WNDDIR + 180)*pi/180)
  
  wind_data_i1 <- rbind(filter(select(wind08_i1, date, direction), 
                              date >= imputed_occurrences1$symptom.date[1]),
                       filter(select(wind09_i1, date, direction), 
                              date >= imputed_occurrences1$symptom.date[2]),
                       filter(select(wind10_i1, date, direction), 
                              date >= imputed_occurrences1$symptom.date[3])) %>%
    filter(month(date) < 8) %>%
    group_by(year = year(date))
  
  wind08_i2 <- read.csv('wind_data/2008i2.csv') %>%
    mutate(date = ymd(DATE),
           direction = (WNDDIR + 180)*pi/180)
  wind09_i2 <- read.csv('wind_data/2009i2.csv') %>%
    mutate(date = ymd(DATE),
           direction = (WNDDIR + 180)*pi/180)
  wind10_i2 <- read.csv('wind_data/2010i2.csv') %>%
    mutate(date = ymd(DATE),
           direction = (WNDDIR + 180)*pi/180)
  
  wind_data_i2 <- rbind(filter(select(wind08_i2, date, direction), 
                               date >= imputed_occurrences2$symptom.date[1]),
                        filter(select(wind09_i2, date, direction), 
                               date >= imputed_occurrences2$symptom.date[2]),
                        filter(select(wind10_i2, date, direction), 
                               date >= imputed_occurrences2$symptom.date[3])) %>%
    filter(month(date) < 8) %>%
    group_by(year = year(date))
  
  wind08_i3 <- read.csv('wind_data/2008i3.csv') 
  wind08_i3 <- wind08_i3 %>% cbind(TIME = str_split_fixed(wind08_i3$DATE, 'T', 2)) %>% 
    select(WNDDIR, WNDSPD, TIME.1, TIME.2) %>%
    mutate(DATE = ymd(TIME.1), IHR = hour(hms(TIME.2))) %>% 
    select(DATE, IHR, WNDDIR, WNDSPD, TIME.2) %>%
    mutate(date = ymd(DATE),
           direction = (WNDDIR + 180)*pi/180)
  wind09_i3 <- read.csv('wind_data/2009i3.csv')
  wind09_i3 <- wind09_i3 %>% cbind(TIME = str_split_fixed(wind09_i3$DATE, 'T', 2)) %>% 
    select(WNDDIR, WNDSPD, TIME.1, TIME.2) %>%
    mutate(DATE = ymd(TIME.1), IHR = hour(hms(TIME.2))) %>% 
    select(DATE, IHR, WNDDIR, WNDSPD, TIME.2) %>%
    mutate(date = ymd(DATE),
           direction = (WNDDIR + 180)*pi/180)
  wind10_i3 <- read.csv('wind_data/2010i3.csv') 
  wind10_i3 <- wind10_i3 %>% cbind(TIME = str_split_fixed(wind10_i3$DATE, 'T', 2)) %>% 
    select(WNDDIR, WNDSPD, TIME.1, TIME.2) %>%
    mutate(DATE = ymd(TIME.1), IHR = hour(hms(TIME.2))) %>% 
    select(DATE, IHR, WNDDIR, WNDSPD, TIME.2) %>%
    mutate(date = ymd(DATE),
           direction = (WNDDIR + 180)*pi/180)
  
  wind_data_i3 <- rbind(filter(select(wind08_i3, date, direction), 
                               date >= imputed_occurrences3$symptom.date[1]),
                        filter(select(wind09_i3, date, direction), 
                               date >= imputed_occurrences3$symptom.date[2]),
                        filter(select(wind10_i3, date, direction), 
                               date >= imputed_occurrences3$symptom.date[3])) %>%
    filter(month(date) < 8) %>%
    group_by(year = year(date))
  
  wind08_i4 <- read.csv('wind_data/2008i4.csv') 
  wind08_i4 <- wind08_i4 %>% cbind(TIME = str_split_fixed(wind08_i4$DATE, 'T', 2)) %>% 
    select(WNDDIR, WNDSPD, TIME.1, TIME.2) %>%
    mutate(DATE = ymd(TIME.1), IHR = hour(hms(TIME.2))) %>% 
    select(DATE, IHR, WNDDIR, WNDSPD, TIME.2) %>%
    mutate(date = ymd(DATE),
           direction = (WNDDIR + 180)*pi/180)
  wind09_i4 <- read.csv('wind_data/2009i4.csv')
  wind09_i4 <- wind09_i4 %>% cbind(TIME = str_split_fixed(wind09_i4$DATE, 'T', 2)) %>% 
    select(WNDDIR, WNDSPD, TIME.1, TIME.2) %>%
    mutate(DATE = ymd(TIME.1), IHR = hour(hms(TIME.2))) %>% 
    select(DATE, IHR, WNDDIR, WNDSPD, TIME.2) %>%
    mutate(date = ymd(DATE),
           direction = (WNDDIR + 180)*pi/180)
  wind10_i4 <- read.csv('wind_data/2010i4.csv') 
  wind10_i4 <- wind10_i4 %>% cbind(TIME = str_split_fixed(wind10_i4$DATE, 'T', 2)) %>% 
    select(WNDDIR, WNDSPD, TIME.1, TIME.2) %>%
    mutate(DATE = ymd(TIME.1), IHR = hour(hms(TIME.2))) %>% 
    select(DATE, IHR, WNDDIR, WNDSPD, TIME.2) %>%
    mutate(date = ymd(DATE),
           direction = (WNDDIR + 180)*pi/180)
  
  wind_data_i4 <- rbind(filter(select(wind08_i4, date, direction), 
                               date >= imputed_occurrences4$symptom.date[1]),
                        filter(select(wind09_i4, date, direction), 
                               date >= imputed_occurrences4$symptom.date[2]),
                        filter(select(wind10_i4, date, direction), 
                               date >= imputed_occurrences4$symptom.date[3])) %>%
    filter(month(date) < 8) %>%
    group_by(year = year(date))
  
  
  ## determine viable years (wind data exists for at least one source)
  yrs <- unique(c(unique(wind_data_i1$year),
                  unique(wind_data_i2$year),
                  unique(wind_data_i3$year),
                  unique(wind_data_i4$year),
                  unique(wind_data_w$year),
                  unique(wind_data_s$year),
                  unique(wind_data_n$year)))
  
  DATA <- Reduce(rbind, lapply(yrs, function(yr){
    print(yr)
    # disease data
    disease_df <- filter(sub_data, year == yr)
    
    # wind data at southern source
    wind_df_s <- filter(wind_data_s, year == yr)
    wind_df_w <- filter(wind_data_w, year == yr)
    wind_df_n <- filter(wind_data_n, year == yr)
    wind_df_i1 <- filter(wind_data_i1, year == yr)
    wind_df_i2 <- filter(wind_data_i2, year == yr)
    wind_df_i3 <- filter(wind_data_i3, year == yr)
    wind_df_i4 <- filter(wind_data_i4, year == yr)
    
    # kernel density estimates for wind directions
    g_theta_s <- density.circular(x = wind_df_s$direction, 
                                    z = disease_df$theta.s, 
                                    bw = 50)$y
    g_theta_n <- density.circular(x = wind_df_n$direction, 
                                    z = disease_df$theta.n, 
                                    bw = 50)$y
    g_theta_w <- density.circular(x = wind_df_w$direction, 
                                    z = disease_df$theta.w, 
                                    bw = 50)$y
    g_theta_i1 <- density.circular(x = wind_df_i1$direction, 
                                    z = disease_df$theta.i1, 
                                    bw = 50)$y
    g_theta_i2 <- density.circular(x = wind_df_i2$direction, 
                                     z = disease_df$theta.i2, 
                                     bw = 50)$y
    g_theta_i3 <- density.circular(x = wind_df_i3$direction, 
                                   z = disease_df$theta.i3, 
                                   bw = 50)$y
    g_theta_i4 <- density.circular(x = wind_df_i4$direction, 
                                   z = disease_df$theta.i4, 
                                   bw = 50)$y
    
    
    # determine first occurrences
    source_s <- filter(first_occurrences_s, year == yr)
    source_w <- filter(first_occurrences_w, year == yr)
    source_n <- filter(first_occurrences_n, year == yr)
    source_i1 <- filter(imputed_occurrences1, year == yr)
    source_i2 <- filter(imputed_occurrences2, year == yr)
    source_i3 <- filter(imputed_occurrences3, year == yr)
    source_i4 <- filter(imputed_occurrences4, year == yr)
    
    # append kernel density estimates to disease data
    out <- disease_df %>%
      mutate(g.theta.s = g_theta_s,
             g.theta.n = g_theta_n,
             g.theta.w = g_theta_w,
             g.theta.i1 = g_theta_i1,
             g.theta.i2 = g_theta_i2,
             g.theta.i3 = g_theta_i3,
             g.theta.i4 = g_theta_i4) %>%
      rename(rho.n = distance.n,
             rho.w = distance.w,
             rho.s = distance.s,
             rho.i1 = distance.i1,
             rho.i2 = distance.i2,
             rho.i3 = distance.i3,
             rho.i4 = distance.i4) %>%
      select(year,
             report.date,
             symptom.date,
             lat, 
             long, 
             host,
             ID,
             rho.n, rho.s, rho.w, rho.i1, rho.i2, rho.i3, rho.i4,
             theta.n, theta.s, theta.w, theta.i1, theta.i2, theta.i3, theta.i4,
             g.theta.n, g.theta.s, g.theta.w, g.theta.i1, g.theta.i2, g.theta.i3, g.theta.i4) %>%
      mutate(rho.n = na_if(rho.n*(symptom.date >= source_n$symptom.date), 0),
             theta.n = na_if(theta.n*(symptom.date >= source_n$symptom.date), 0),
             g.theta.n = na_if(g.theta.n*(symptom.date >= source_n$symptom.date), 0),
             rho.s = na_if(rho.s*(symptom.date >= source_s$symptom.date), 0),
             theta.s = na_if(theta.s*(symptom.date >= source_s$symptom.date), 0),
             g.theta.s = na_if(g.theta.s*(symptom.date >= source_s$symptom.date), 0),
             rho.w = na_if(rho.w*(symptom.date >= source_w$symptom.date), 0),
             theta.w = na_if(theta.w*(symptom.date >= source_w$symptom.date), 0),
             g.theta.w = na_if(g.theta.w*(symptom.date >= source_w$symptom.date), 0),
             rho.i1 = na_if(rho.i1*(symptom.date >= source_i1$symptom.date), 0),
             theta.i1 = na_if(theta.i1*(symptom.date >= source_i1$symptom.date), 0),
             g.theta.i1 = na_if(g.theta.i1*(symptom.date >= source_i1$symptom.date), 0),
             rho.i2 = na_if(rho.i2*(symptom.date >= source_i2$symptom.date), 0),
             theta.i2 = na_if(theta.i2*(symptom.date >= source_i2$symptom.date), 0),
             g.theta.i2 = na_if(g.theta.i2*(symptom.date >= source_i2$symptom.date), 0),
             rho.i3 = na_if(rho.i3*(symptom.date >= source_i3$symptom.date), 0),
             theta.i3 = na_if(theta.i3*(symptom.date >= source_i3$symptom.date), 0),
             g.theta.i3 = na_if(g.theta.i3*(symptom.date >= source_i3$symptom.date), 0),
             rho.i4 = na_if(rho.i4*(symptom.date >= source_i4$symptom.date), 0),
             theta.i4 = na_if(theta.i4*(symptom.date >= source_i4$symptom.date), 0),
             g.theta.i4 = na_if(g.theta.i4*(symptom.date >= source_i4$symptom.date), 0))
    
    return(out)
  }))
  
  return(list(data = DATA,
              sources.s = first_occurrences_s,
              sources.w = first_occurrences_w,
              sources.n = first_occurrences_n,
              sources.i1 = imputed_occurrences1,
              sources.i2 = imputed_occurrences2,
              sources.i3 = imputed_occurrences3,
              sources.i4 = imputed_occurrences4,
              wind.s = wind_data_s,
              wind.w = wind_data_w,
              wind.n = wind_data_n,
              wind.i1 = wind_data_i1,
              wind.i2 = wind_data_i2,
              wind.i3 = wind_data_i3,
              wind.i4 = wind_data_i4))
  #######################
}
