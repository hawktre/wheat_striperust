library(ggmap)
library(tidyverse)
library(gganimate)
library(lubridate)
load('cdm_data.RData')
load('pred_data.RData')

sub_data <- filter(data, 
                   year(symptom_date) %in% 2008:2010, 
                   !is.na(lat),
                   !is.na(long),
                   lat < 47,
                   long < 0,
                   long > -110, 
                   planting == 'sentinel') %>%
  select(lat, long, symptom_date, locID) %>%
  group_by(locID, year(symptom_date)) %>%
  arrange(symptom_date) %>%
  slice_min(symptom_date, n = 1) %>%
  slice_sample(n = 1) %>%
  ungroup() %>%
  mutate(Year = paste("Year: ", year(symptom_date), sep = ''))


lonrange <- range(sub_data$long)
latrange <- range(sub_data$lat)
mapbox <- make_bbox(lonrange, latrange)
map <- get_map(location = mapbox,
               maptype = 'toner-background')

pred_df <- rename(pred_df, symptom_date = date)

p <- ggmap(map, darken = c(0.4, 'grey')) +
  geom_point(data = filter(sub_data, 
                           month(symptom_date) < 8, 
                           year(symptom_date) == 2008), 
             aes(x = long, 
                 y = lat, 
                 color = month(symptom_date), 
                 group = week(symptom_date)),
             alpha = 0.9,
             size = 3) +
  facet_wrap(~Year) +
  scale_color_gradientn(colors = c('yellow', 'orange', 'red', 'black'),
                        breaks = 1:7,
                        labels = month.name[1:7]) +
  guides(color = guide_legend('Month')) +
  labs(x = 'Longitude', y = 'Latitude') +
  theme(text = element_text(size = 14)) +
  theme_bw() +
  geom_path(data = filter(pred_df, year(symptom_date) == 2008),
            aes(x = ats2.src1.lon, y = ats2.src1.lat,
                color = month(symptom_date),
                group = week(symptom_date))) +
  geom_path(data = filter(pred_df, year(symptom_date) == 2008),
            aes(x = ats2.src3.lon, y = ats2.src3.lat,
                color = month(symptom_date),
                group = week(symptom_date)),
                linetype = 6) +
  geom_point(data = source_df[c(1,3), ],
             aes(x = cc.long, y = cc.lat),
             color = 'red',
             shape = 3)


p


anim <- p + transition_states(month(symptom_date), 
                      1, 0.25) + 
  enter_fade() + exit_fade()

dirname <- paste(getwd(), '/animations', sep = '')
animate(anim, nframes = 30,
        device = 'png',
        renderer = file_renderer(dir = dirname, 
                                 prefix = '2008ats_n',
                                 overwrite = T),
        width = 6, height = 4, units = 'in',
        res = 200)

# anim_save(filename = 'aos2008_anim.gif')
