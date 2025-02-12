library(ggmap)
library(tidyverse)
library(lubridate)
library(geosphere)
library(housingData)
library(circular)

load('cdm_data.RData')
data <- as_tibble(data)
sub_data <- data %>%
  filter(# planting == 'sentinel',
    # month(report_date) < 8,
    # year(report_date) > 2007,
    month(symptom_date) <= 6,
    year(symptom_date) > 2007,
    !is.na(lat), !is.na(long),
    lat > 40,
    lat < 45,
    long > -86) %>%
  rename(report.date = report_date,
         symptom.date = symptom_date) %>%
  select(lat, long, report.date, symptom.date, 
         state, county, locID, host, planting) %>%
  mutate(year = year(report.date)) %>%
  arrange(symptom.date)


lonrange <- range(sub_data$long)
latrange <- range(sub_data$lat)
mapbox <- make_bbox(lonrange, latrange)
map <- get_map(location = mapbox,
               maptype = 'toner-background')

ggmap(map, darken = c(0.3, 'white')) +
  geom_point(data = sub_data,
             aes(x = long, 
                 y = lat,
                 color = factor(week(symptom.date),
                                labels = paste('june week ', 1:4)),
                 shape = planting),
             size = 2) +
  facet_wrap(~year) +
  # scale_color_brewer(palette = 'Spectral') +
  guides(color = guide_legend('symptom date'))
