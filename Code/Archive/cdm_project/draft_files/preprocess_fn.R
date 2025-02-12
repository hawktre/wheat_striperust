library(tidyverse)
library(lubridate)
library(geosphere)
library(housingData)
library(circular)

preprocess_fn <- function(){
  imputed_state = rep('MI', 9)
  imputed_county = rep('bay', 9)
  imputed_date = ymd(c('2008/01/01',
                       '2009/01/01',
                       '2010/01/01',
                       '2011/01/01',
                       '2012/01/01',
                       '2013/01/01',
                       '2014/01/01',
                       '2015/01/01',
                       '2016/01/01'))
  week(imputed_date) <- rep(27, 9)
  
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
  
  imputed_occurrences <- Reduce(rbind, lapply(1:9, function(yr){
    out <- filter(geoCounty, 
                  state == imputed_state[yr], 
                  rMapCounty == imputed_county[yr]) %>%
      mutate(year = 2007 + yr) %>%
      rename(cc.lat.i = lat,
             cc.long.i = lon) %>%
      select(year, cc.long.i, cc.lat.i) %>%
      mutate(lat = NA,
             long = NA,
             report.date = imputed_date[yr],
             symptom.date = imputed_date[yr],
             state = imputed_state[yr],
             county = imputed_county[yr],
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
  sub_data <- select(imputed_occurrences, cc.long.i, cc.lat.i, year) %>%
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
                     bearing.i = bearing(as.matrix(select(sub_data, cc.long.i, cc.lat.i)),
                                         as.matrix(select(sub_data, long, lat))),
                     distance.i = distGeo(as.matrix(select(sub_data, cc.long.i, cc.lat.i)),
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
           theta.i = as.circular(bearing.trans(bearing.i)*pi/180,
                                 zero = 0,
                                 units = 'radians',
                                 rotation = 'counter')) %>%
    as_tibble()
  
  ## import wind data
  
  wind08_s <- read.csv('Raw_Zedx_data/2008_data/2008_FL_Collier_wxdata.csv') %>%
    mutate(date = ymd(DATE),
           direction = (WNDDIR + 180)*pi/180)
  wind09_s <- read.csv('Raw_Zedx_data/2009_data/2009_FL_Miami-Dade_wxdata.csv') %>%
    mutate(date = ymd(DATE),
           direction = (WNDDIR + 180)*pi/180)
  wind10_s <- read.csv('Raw_Zedx_data/2010_data/2010_FL_Alachua_wxdata.csv') %>%
    mutate(date = ymd(DATE),
           direction = (WNDDIR + 180)*pi/180)
  wind11_s <- read.csv('Raw_Zedx_data/2011_data/2011_FL_Miami-Dade_wxdata.csv') %>%
    mutate(date = ymd(DATE),
           direction = (WNDDIR + 180)*pi/180)
  # wind12_s <- read.csv('Raw_Zedx_data/2012_data/2012_FL_Collier_wxdata.csv')
  # wind13_s <- read.csv('Raw_Zedx_data/2013_data/2013_FL_Miami-Dade_wxdata.csv')
  wind14_s <- read.csv('Raw_Zedx_data/2014_data/2014_FL_Miami-Dade_wxdata.csv') %>%
    mutate(date = ymd(DATE),
           direction = (WNDDIR + 180)*pi/180)
  # wind15_s <- read.csv('Raw_Zedx_data/2015_data/2015_SC_Miami-Dade_wxdata.csv')
  # wind16_s <- read.csv('Raw_Zedx_data/2016_data/2016_FL_Miami-Dade_wxdata.csv')
  
  wind_data_s <- rbind(filter(select(wind08_s, -DOY), 
                              date >= first_occurrences_s$symptom.date[1]),
                       filter(wind09_s, 
                              date >= first_occurrences_s$symptom.date[2]),
                       filter(wind10_s, 
                              date >= first_occurrences_s$symptom.date[3]), 
                       filter(wind11_s, 
                              date >= first_occurrences_s$symptom.date[4]), 
                       filter(wind14_s, 
                              date >= first_occurrences_s$symptom.date[7])) %>%
    filter(month(date) < 8) %>%
    select(date, direction) %>%
    group_by(year = year(date))
  
  
  wind08_w <- read.csv('Raw_Zedx_data/2008_data/2008_LA_Vernon_wxdata.csv') %>%
    mutate(date = ymd(DATE),
           direction = (WNDDIR + 180)*pi/180)
  # wind09_w <- read.csv('Raw_Zedx_data/2009_data/2009_OK_Payne_wxdata.csv') %>%
  #   mutate(date = ymd(DATE),
  #          direction = (WNDDIR + 180)*pi/180)
  # wind10_w <- read.csv('Raw_Zedx_data/2010_data/2010_TX_Brazos_wxdata.csv') %>%
  #   mutate(date = ymd(DATE),
  #          direction = (WNDDIR + 180)*pi/180)
  wind14_w <- read.csv('Raw_Zedx_data/2014_data/2014_LA_Franklin_wxdata.csv') %>%
    mutate(date = ymd(DATE),
           direction = (WNDDIR + 180)*pi/180)
  wind15_w <- read.csv('Raw_Zedx_data/2015_data/2015_LA_East_Baton_Rouge_wxdata.csv') %>%
    mutate(date = ymd(DATE),
           direction = (WNDDIR + 180)*pi/180)
  
  wind_data_w <- rbind(filter(select(wind08_w, date, direction), 
                              date >= first_occurrences_w$symptom.date[1]),
                       filter(select(wind14_w, date, direction), 
                              date >= first_occurrences_w$symptom.date[7]),
                       filter(select(wind15_w, date, direction), 
                              date >= first_occurrences_w$symptom.date[8])) %>%
    filter(month(date) < 8) %>%
    group_by(year = year(date))
  
  
  wind08_n <- read.csv('Raw_Zedx_data/2008_data/2008_OH_Sandusky_wxdata.csv') %>%
    mutate(date = ymd(DATE),
           direction = (WNDDIR + 180)*pi/180)
  wind09_n <- read.csv('Raw_Zedx_data/2009_data/2009_OH_Huron_wxdata.csv') %>%
    mutate(date = ymd(DATE),
           direction = (WNDDIR + 180)*pi/180)
  wind10_n <- read.csv('Raw_Zedx_data/2010_data/2010_OH_Wayne_wxdata.csv') %>%
    mutate(date = ymd(DATE),
           direction = (WNDDIR + 180)*pi/180)
  wind11_n <- read.csv('Raw_Zedx_data/2011_data/2011_WI_Dane_wxdata.csv') %>%
    mutate(date = ymd(DATE),
           direction = (WNDDIR + 180)*pi/180)
  # wind12_n <- read.csv('Raw_Zedx_data/2012_data/2012_NY_Suffolk_wxdata.csv') %>%
  #   mutate(date = ymd(DATE),
  #          direction = (WNDDIR + 180)*pi/180)
  wind13_n <- read.csv('Raw_Zedx_data/2013_data/2013_NY_Suffolk_wxdata.csv') %>%
    mutate(date = ymd(DATE),
           direction = (WNDDIR + 180)*pi/180)
  wind15_n <- read.csv('Raw_Zedx_data/2015_data/2015_OH_Sandusky_wxdata.csv') %>%
    mutate(date = ymd(DATE),
           direction = (WNDDIR + 180)*pi/180)
  # wind16_n <- read.csv('Raw_Zedx_data/2016_data/2016_NJ_Hunterdon_wxdata.csv') %>%
  #   mutate(date = ymd(DATE),
  #          direction = (WNDDIR + 180)*pi/180)
  
  wind_data_n <- rbind(filter(select(wind08_n, date, direction), 
                              date >= first_occurrences_n$symptom.date[1]),
                       filter(select(wind09_n, date, direction), 
                              date >= first_occurrences_n$symptom.date[2]),
                       filter(select(wind10_n, date, direction), 
                              date >= first_occurrences_n$symptom.date[3]),
                       filter(select(wind11_n, date, direction), 
                              date >= first_occurrences_n$symptom.date[4]),
                       filter(select(wind13_n, date, direction), 
                              date >= first_occurrences_n$symptom.date[6]),
                       filter(select(wind15_n, date, direction), 
                              date >= first_occurrences_n$symptom.date[8])) %>%
    filter(month(date) < 8) %>%
    group_by(year = year(date))
  
  
  
  wind08_i <- read.csv('Raw_Zedx_data/2008_data/2008_MI_Bay_1_wxdata.csv') %>%
    mutate(date = ymd(DATE),
           direction = (WNDDIR + 180)*pi/180)
  wind09_i <- read.csv('Raw_Zedx_data/2009_data/2009_MI_Bay_wxdata.csv') %>%
    mutate(date = ymd(DATE),
           direction = (WNDDIR + 180)*pi/180)
  wind10_i <- read.csv('Raw_Zedx_data/2010_data/2010_MI_Bay_wxdata.csv') %>%
    mutate(date = ymd(DATE),
           direction = (WNDDIR + 180)*pi/180)
  wind11_i <- read.csv('Raw_Zedx_data/2011_data/2011_MI_Bay_wxdata.csv') %>%
    mutate(date = ymd(DATE),
           direction = (WNDDIR + 180)*pi/180)
  # wind12_i <- read.csv('Raw_Zedx_data/2012_data/2012_MI_Bay_wxdata.csv') %>%
  #   mutate(date = ymd(DATE),
  #          direction = (WNDDIR + 180)*pi/180)
  wind13_i <- read.csv('Raw_Zedx_data/2013_data/2013_MI_Bay_wxdata.csv') %>%
    mutate(date = ymd(DATE),
           direction = (WNDDIR + 180)*pi/180)
  # wind14_i <- read.csv('Raw_Zedx_data/2014_data/2014_MI_Bay_wxdata.csv') %>%
  #   mutate(date = ymd(DATE),
  #          direction = (WNDDIR + 180)*pi/180)
  wind15_i <- read.csv('Raw_Zedx_data/2015_data/2015_MI_Bay_wxdata.csv') %>%
    mutate(date = ymd(DATE),
           direction = (WNDDIR + 180)*pi/180)
  wind16_i <- read.csv('Raw_Zedx_data/2016_data/2016_MI_Bay_wxdata.csv') %>%
    mutate(date = ymd(DATE),
           direction = (WNDDIR + 180)*pi/180)
  
  wind_data_i <- rbind(filter(select(wind08_i, date, direction), 
                              date >= imputed_occurrences$symptom.date[1]),
                       filter(select(wind09_i, date, direction), 
                              date >= imputed_occurrences$symptom.date[2]),
                       filter(select(wind10_i, date, direction), 
                              date >= imputed_occurrences$symptom.date[3]),
                       filter(select(wind11_i, date, direction), 
                              date >= imputed_occurrences$symptom.date[4]),
                       filter(select(wind13_i, date, direction), 
                              date >= imputed_occurrences$symptom.date[6]),
                       filter(select(wind15_i, date, direction), 
                              date >= imputed_occurrences$symptom.date[8]),
                       filter(select(wind16_i, date, direction), 
                              date >= imputed_occurrences$symptom.date[9])) %>%
    filter(month(date) < 8) %>%
    group_by(year = year(date))
  
  ## determine viable years (wind data exists for at least one source)
  yrs <- unique(c(unique(wind_data_i$year),
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
    wind_df_i <- filter(wind_data_i, year == yr)
    
    # kernel density estimates for wind directions
    if(nrow(wind_df_s) > 0){
      g_theta_s <- density.circular(x = wind_df_s$direction, 
                                    z = disease_df$theta.s, 
                                    bw = 50)$y
    }else{
      g_theta_s <- NA
    }
    if(nrow(wind_df_n) > 0){
      g_theta_n <- density.circular(x = wind_df_n$direction, 
                                    z = disease_df$theta.n, 
                                    bw = 50)$y
    }else{
      g_theta_n <- NA
    }
    if(nrow(wind_df_w) > 0){
      g_theta_w <- density.circular(x = wind_df_w$direction, 
                                    z = disease_df$theta.w, 
                                    bw = 50)$y
    }else{
      g_theta_w <- NA
    }
    if(nrow(wind_df_i) > 0){
      g_theta_i <- density.circular(x = wind_df_i$direction, 
                                    z = disease_df$theta.i, 
                                    bw = 50)$y
    }else{
      g_theta_i <- NA
    }
    
    # determine first occurrences
    source_s <- filter(first_occurrences_s, year == yr)
    source_w <- filter(first_occurrences_w, year == yr)
    source_n <- filter(first_occurrences_n, year == yr)
    source_i <- filter(imputed_occurrences, year == yr)
    
    # append kernel density estimates to disease data
    out <- disease_df %>%
      mutate(g.theta.s = g_theta_s,
                g.theta.n = g_theta_n,
                g.theta.w = g_theta_w,
                g.theta.i = g_theta_i) %>%
      rename(rho.n = distance.n,
             rho.w = distance.w,
             rho.s = distance.s,
             rho.i = distance.i) %>%
      select(year,
             report.date,
             symptom.date,
             lat, 
             long, 
             host,
             ID,
             rho.n, rho.s, rho.w, rho.i, 
             theta.n, theta.s, theta.w, theta.i,
             g.theta.n, g.theta.s, g.theta.w, g.theta.i) %>%
      mutate(rho.n = na_if(rho.n*(symptom.date >= source_n$symptom.date), 0),
             theta.n = na_if(theta.n*(symptom.date >= source_n$symptom.date), 0),
             g.theta.n = na_if(g.theta.n*(symptom.date >= source_n$symptom.date), 0),
             rho.s = na_if(rho.s*(symptom.date >= source_s$symptom.date), 0),
             theta.s = na_if(theta.s*(symptom.date >= source_s$symptom.date), 0),
             g.theta.s = na_if(g.theta.s*(symptom.date >= source_s$symptom.date), 0),
             rho.w = na_if(rho.w*(symptom.date >= source_w$symptom.date), 0),
             theta.w = na_if(theta.w*(symptom.date >= source_w$symptom.date), 0),
             g.theta.w = na_if(g.theta.w*(symptom.date >= source_w$symptom.date), 0),
             rho.i = na_if(rho.i*(symptom.date >= source_i$symptom.date), 0),
             theta.i = na_if(theta.i*(symptom.date >= source_i$symptom.date), 0),
             g.theta.i = na_if(g.theta.i*(symptom.date >= source_i$symptom.date), 0))
    
    return(out)
  }))
  
  return(list(data = DATA,
              sources.s = first_occurrences_s,
              sources.w = first_occurrences_w,
              sources.n = first_occurrences_n,
              sources.i = imputed_occurrences,
              wind.s = wind_data_s,
              wind.w = wind_data_w,
              wind.n = wind_data_n,
              wind.i = wind_data_i))
  #######################
}
