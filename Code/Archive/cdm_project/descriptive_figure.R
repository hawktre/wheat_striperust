library(ggmap)
library(tidyverse)
library(lubridate)
load('cdm_data.RData')

sub_data <- filter(data, 
                   year(symptom_date) %in% 2008:2010, 
                   !is.na(lat),
                   !is.na(long),
                   lat < 47,
                   long < 0,
                   long > -110)


lonrange <- range(sub_data$long)
latrange <- range(sub_data$lat)
mapbox <- make_bbox(lonrange, latrange)
map <- get_map(location = mapbox,
               maptype = 'toner-background')


ggmap(map, darken = c(0.4, 'grey')) +
  ## plot errors
  geom_point(data = filter(sub_data, month(symptom_date) < 8), 
             aes(x = long, 
                 y = lat, 
                 shape = planting,
                 # color = factor(month(symptom_date), 
                 #                labels = month.name[1:7])), 
                 color = month(symptom_date)),
             alpha = 0.9) +
  facet_wrap(~year(symptom_date), nrow = 1) +
  scale_color_gradientn(colors = c('yellow', 'orange', 'red', 'black'),
                        breaks = 1:7,
                        labels = month.name[1:7]) +
#  scale_size_manual(values = 1) +
  # scale_colour_gradient2(low = 'yellow', 
  #                        mid = 'red', 
  #                        high = 'black', 
  #                        midpoint = 4.5,
  #                        limits = c(1, 7)) +
  guides(color = guide_legend('Month'),
         shape = guide_legend('Plot type'),
         size = guide_none()) +
  labs(x = 'Longitude', y = 'Latitude') +
  theme(text = element_text(size = 10))

ggsave(filename = 'fig2.tiff',
       width = 7.5, height = 5, units = 'in', 
       dpi = 300)
