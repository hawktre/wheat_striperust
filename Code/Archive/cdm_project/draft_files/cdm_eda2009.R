options(stringsAsFactors = F)
library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(ggmap)
library(xtable)

raw_data <- read.csv('cdm_data09.csv')
str(raw_data)

sub_data <- raw_data %>%
  transmute(date = as.Date(Symptom_date, origin = '2009-01-01'),
         state = State,
         county = County,
         planting = Planting,
         host = Host,
         field_size = Field.size..Acres.,
         incidence = Dis_incid,
         severity = Dis_sev,
         occurrence = (Dis_sev != 0),
         lat = Latitude,
         long = Longitude) # variables of interest

locations <- distinct(sub_data, lat, long) %>%
  mutate(locID = 1:(nrow(distinct(sub_data, lat, long))))

data <- merge(sub_data, locations, by = c('lat', 'long'))

n_distinct(data$county)
n_distinct(data$state)
n_distinct(data$locID)
n_distinct(data$host)

tbl <- xtabs(occurrence ~ state + host + planting, data = data)
apply(tbl, c(3, 2), sum)


ggplot(aes(x = yday(date),
           color = host),
       data = data) +
  geom_density(adjust = 0.75) +
  facet_wrap(~host*planting, 
             drop = F, nrow = 5)

ggplot(aes(x = yday(date),
           color = host),
       data = data) +
  geom_density(adjust = 0.75)

plot_df1 <- group_by(data, state, host, planting) %>%
  summarize(statewide_occurrences = sum(occurrence),
            mean_sev = mean(severity),
            mean_inc = mean(incidence))
                    
ggplot(aes(x = statewide_occurrences,
           color = host), 
       data = plot_df1) +
  geom_density(adjust = 0.8) +
  facet_wrap(~host*planting, 
             drop = F, nrow = 5)

ggplot(aes(x = statewide_occurrences,
           color = host), 
       data = plot_df1) +
  geom_density(adjust = 2)


plot_df2 <- group_by(data, locID, host, planting) %>%
  summarize(init_date = min(date)) %>% 
  merge(locations, by = 'locID')

bbox <- make_bbox(locations$long, locations$lat)
map <- get_stamenmap(bbox, zoom = calc_zoom(bbox), maptype = 'toner')
ggmap(map, darken = c(0.4, "white"), extent = 'normal') +
  geom_point(aes(x = long, y = lat, color = yday(init_date)),
             data = plot_df2) +
  facet_wrap(~host*planting, drop = F, nrow = 5) +
  scale_color_gradientn(colors = terrain.colors(100)) +
  coord_fixed(expand = F)


ggmap(map, darken = c(0.4, "white"), extent = 'normal') +
  geom_point(aes(x = long, y = lat, color = yday(init_date)),
             data = plot_df2) +
  facet_wrap(~host) +
  scale_color_gradientn(colors = terrain.colors(100)) 

ggmap(map, darken = c(0.4, "white"), extent = 'normal') +
  geom_point(aes(x = long, y = lat, color = yday(init_date)),
             data = plot_df2) +
  scale_color_gradientn(colors = terrain.colors(100)) 
