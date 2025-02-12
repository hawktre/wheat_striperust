options(stringsAsFactors = F)
library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(ggmap)
library(pander)
load('cdm_data.RData')

data[602, 3]
data[602, 3] <- ymd("2008/08/22") #typo

# variable summary
var_summary_tbl <- group_by(data, year(report_date)) %>%
  summarize(n_counties = n_distinct(county),
            n_states = n_distinct(state),
            n_locations = n_distinct(locID),
            n_occurrences = sum(occurrence, na.rm = T),
            n_reports = length(report_date)) %>%
  rename(year = `year(report_date)`)

pander(var_summary_tbl)


## occurrence tables
tbl <- xtabs(occurrence ~ planting + host + year(report_date), data = data)
tbl

# counts by year and host
apply(tbl, c(3, 2), sum)

# counts by year and planting
apply(tbl, c(3, 1), sum)

# counts by host and planting
apply(tbl, c(2, 1), sum)


## distributions of occurrences over time of year, by host and planting
ggplot(aes(x = yday(report_date),
           color = host),
       data = filter(data, planting != 'other')) +
  geom_density(adjust = 1.5) +
  facet_wrap(~planting, nrow = 2,
             scales = 'free_y') 

## distributions of occurrences over time of year by host
ggplot(aes(x = yday(report_date),
           color = host),
       data = filter(data, planting != 'other')) +
  geom_density(adjust = 1.5)


## map of report dates aggregated over all years, plantings, hosts
locations <- select(data, lat, long, locID) %>% distinct()
plot_df1 <- group_by(data, year(report_date), locID, host) %>%
  summarize(first_report = min(report_date)) %>% 
  merge(locations, by = 'locID') %>%
  filter(long < 0) # filter step removes 3 rows

bbox <- make_bbox(plot_df1$long, plot_df1$lat)
map <- get_stamenmap(bbox, zoom = calc_zoom(bbox), maptype = 'toner')

ggmap(map, darken = c(0.4, "white"), extent = 'normal') +
  geom_point(aes(x = long, y = lat, color = yday(first_report)),
             data = plot_df1) +
  scale_color_gradientn(colors = terrain.colors(100)) +
  theme_bw()

## map of report dates by host, aggregated over all years and plantings
ggmap(map, darken = c(0.4, "white"), extent = 'normal') +
  geom_point(aes(x = long, y = lat, color = yday(first_report)),
             data = plot_df1) +
               # filter(plot_df1, 
               #             host %in% levels(data$host)[c(2:4, 7)])) +
  facet_wrap(~host, drop = T, nrow = 2) +
  scale_color_gradientn(colors = terrain.colors(100)) +
  theme_bw()

## map of report dates by year, aggregated over all hosts and plantings
ggmap(map, darken = c(0.4, "white"), extent = 'normal') +
  geom_point(aes(x = long, y = lat, color = yday(first_report)),
             data = plot_df1) +
  facet_wrap(~`year(report_date)`, drop = T, nrow = 3) +
  scale_color_gradientn(colors = terrain.colors(100)) +
  theme_bw()
