library(tidyverse)
library(lubridate)
library(ggmap)
library(circular)
load("cdm_data.RData")

# subset: 10 northernmost occurrences each year
northernmost <- as_tibble(data) %>%
  group_by(year(report_date)) %>%
  top_n(10, lat) %>%
  arrange(symptom_date) %>%
  ungroup() %>%
  select(lat,
         long,
         symptom_date,
         report_date,
         state, 
         county,
         planting)

# spatial extent (for plot)
lonrange <- range(northernmost$long)
latrange <- range(northernmost$lat)

# retrieve map (for plot)
mapbox <- make_bbox(lonrange, latrange)
map <- get_map(location = mapbox,
               maptype = 'toner-background')

# plot occurrences
ggmap(map, darken = c(0.5, 'white')) +
  geom_point(data = northernmost,
            aes(x = long,
                y = lat,
                #group = factor(year),
                color = factor(year(symptom_date))))


# investigate points in NC, N Ontario
filter(northernmost, lat < 40)
filter(northernmost, lat > 48)

# drop point (typo in year of report date)
northernmost <- filter(northernmost, lat > 40, lat < 48)

# repeat plot
lonrange <- range(northernmost$long)
latrange <- range(northernmost$lat) + c(-2, 0)
mapbox <- make_bbox(lonrange, latrange)
map <- get_map(location = mapbox,
               maptype = 'toner-background')
ggmap(map, darken = c(0.5, 'white')) +
  geom_point(data = northernmost,
             aes(x = long,
                 y = lat,
                 #group = factor(year),
                 color = factor(year(symptom_date)))) +
  labs(x = 'Latitude',
       y = 'Longitude') +
  guides(color = guide_legend('Year'))

# # date range
# northernmost %>%
#   group_by(year(symptom_date)) %>%
#   summarize(earliest = min(symptom_date),
#             latest = max(symptom_date),
#             avg = mean(symptom_date))
# 
# ggplot(data = northernmost, 
#        aes(x = week(symptom_date),
#            fill = factor(year(symptom_date)),
#            group = factor(year(symptom_date)))) +
#   geom_histogram(position = 'dodge',
#                  binwidth = 2)
# 
# ggplot(data = northernmost, 
#        aes(x = week(symptom_date))) +
#   geom_histogram(binwidth = 2)

# earlier occurrences
northernmost %>%
  group_by(year(symptom_date)) %>%
  arrange(symptom_date) %>%
  top_n(-3, symptom_date) %>%
  summarize(earliest = min(symptom_date),
            latest = max(symptom_date),
            avg = mean(symptom_date))

first_north <- northernmost %>%
  group_by(year(symptom_date)) %>%
  top_n(-1, symptom_date)

# earliest occurrences are in michigan, usually bay county
ggmap(map, darken = c(0.5, 'white')) +
  geom_point(data = filter(first_north),
             aes(x = long,
                 y = lat,
                 color = factor(week(symptom_date)))) +
  labs(x = 'Latitude',
       y = 'Longitude') +
  guides(color = guide_legend('Week'))

# usually around jul 1st overall
mean(week(first_north$symptom_date))
first_north$symptom_date
week('2008-07-01')

# a little later in MI
first_north$symptom_date[first_north$state == 'MI']

# non-MI first northmost occurrences
filter(first_north, state != 'MI')

# given wind availability, NY another option

# check wind data
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
         month(date) >= 7) %>%
  select(date, direction) %>%
  group_by(year = year(date))

layout(matrix(1:8, nrow = 2, byrow = T))
lapply(unique(windMI$year), function(yr){
  x <- filter(windMI, year == yr)
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

# repeat for NY
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
         month(date) >= 7) %>%
  select(date, direction) %>%
  group_by(year = year(date))

layout(matrix(1:6, nrow = 2, byrow = T))
lapply(unique(windNY$year), function(yr){
  x <- filter(windNY, year == yr)
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
