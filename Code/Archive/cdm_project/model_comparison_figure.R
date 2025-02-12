library(ggmap)
source('preprocess_fn_v2.R')
source('fit_twosource_fn.R')
source('plot_fn.R')

## preprocess data
preprocess_out <- preprocess_fn()
DATA <- preprocess_out$data 
out_tbl <- NULL

# select year
yrs <- 2008:2010 # possible years

# source information (not year sensitive)
src_type <- 'observed'
wind_src1 <- preprocess_out$wind.s
wind_src2 <- preprocess_out$wind.n
src1_id <- filter(preprocess_out$sources.s,
                  year == yrs[j])$ID
src2_id <- filter(preprocess_out$sources.n,
               year == yrs[j])$ID
# src2_id <- NULL


# source information (year sensitive)
j <- 1 # which year
source_df <- rbind(rename(filter(preprocess_out$sources.s, year == yrs[j]),
                          cc.long = cc.long.s,
                          cc.lat = cc.lat.s),
                   rename(filter(preprocess_out$sources.n, year == yrs[j]),
                          cc.long = cc.long.n,
                          cc.lat = cc.lat.n))


# subset data to year and calculate basis expansions and response variables
sub_data <- filter(DATA, year(symptom.date) == yrs[j]) %>%
  rename(rho.src1 = rho.s,
         theta.src1 = theta.s,
         g.theta.src1 = g.theta.s,
         rho.src2 = rho.n,
         theta.src2 = theta.n,
         g.theta.src2 = g.theta.n) %>%
  mutate(month = month(symptom.date),
         day = day(symptom.date),
         time = yday(symptom.date),
         iresponse.src1 = log(1 + rho.src1), # isotropic response for first source
         iresponse.src2 = log(1 + rho.src2), # isotropic response for second source
         aresponse.src1 = log(1 + rho.src1/g.theta.src1), # anisotropic response for first source
         aresponse.src2 = log(1 + rho.src2/g.theta.src2), # anisotropic response for second source
         sinb1.src1 = sin(as.numeric(theta.src1)), # sin basis 1 for first source
         sinb2.src1 = sin(as.numeric(theta.src1) + pi/4), # sin basis 2 for first source
         sinb1.src2 = sin(as.numeric(theta.src2)), # sin basis 1 for second source
         sinb2.src2 = sin(as.numeric(theta.src2) + pi/4), # sin basis 2 for second source
         z.reg.src1 = rho.src1/g.theta.src1, # first source regressor to initialize weights
         z.reg.src2 = rho.src2/g.theta.src2, # second source regressor to initialize weights
         z.fit = rho.src1 < rho.src2) # is first source closer than second source?
        
sub_data <- sub_data %>%
  group_by(ID) %>%
  sample_n(1)


plot_out <- plot_fn(sub_data, wind_src1, wind_src2, source_df, src1_id, src2_id, src_type)
out_tbl <- rbind(plot_out$tbl, out_tbl) %>% arrange(year)

pred_df <- plot_out$pred
plot_df <- plot_out$plot


lonrange <- range(sub_data$long)
latrange <- range(sub_data$lat)
mapbox <- make_bbox(lonrange, latrange)
map <- get_map(location = mapbox,
               maptype = 'toner-background')


ggmap(map, darken = c(0.8, 'white')) +
  geom_point(data = plot_df, 
             aes(x = long, 
                 y = lat, 
                 # shape = factor(which.pred,
                 #                labels = c('N', 'S')),
                 shape = factor(which.pred,
                                labels = c('S', 'N')),
                 size = abs(error)), 
             alpha = 0.5) +
  facet_wrap(~method) +
  # scale_size_continuous(breaks = c(500, 1000, 2000),
  #                       limits = c(1, 4000),
  #                       range = c(1, 8)) +
  scale_size_continuous(breaks = c(14, 28, 42),
                        range = c(1, 8)) +
  geom_point(data = mutate(source_df[1, ], method = 'aos'),
             aes(x = cc.long, y = cc.lat),
             color = 'red') +
  geom_point(data = mutate(source_df, method = 'ats'),
             aes(x = cc.long, y = cc.lat),
             color = 'red') +
  geom_point(data = mutate(source_df, method = 'its'),
             aes(x = cc.long, y = cc.lat),
             color = 'red') +
  geom_point(data = mutate(source_df[1, ], method = 'ios'),
             aes(x = cc.long, y = cc.lat),
             color = 'red') +
  geom_path(data = mutate(filter(pred_df), #theta < pi/2 | theta > 3*pi/2),
                          method = 'aos'),
            aes(x = aos.lon,
                y = aos.lat,
                color = factor(month(date)),
                group = week(date))) +
  geom_path(data = mutate(filter(pred_df, date > source_df[2, ]$symptom.date), 
                          method = 'ats'),
            aes(x = ats.src2.lon,
                y = ats.src2.lat,
                color = factor(month(date)),
                group = week(date))) +
  geom_path(data = mutate(filter(pred_df), # theta < pi/2 | theta > 3*pi/2),
                          method = 'ats'),
            aes(x = ats.src1.lon,
                y = ats.src1.lat,
                color = factor(month(date)),
                group = week(date))) +
  geom_path(data = mutate(filter(pred_df, date > source_df[2, ]$symptom.date),
                          method = 'its'),
            aes(x = its.src2.lon,
                y = its.src2.lat,
                color = factor(month(date)),
                group = week(date))) +
  geom_path(data = mutate(pred_df, method = 'its'),
            aes(x = its.src1.lon,
                y = its.src1.lat,
                color = factor(month(date)),
                group = week(date))) +
  geom_path(data = mutate(pred_df, method = 'ios'),
            aes(x = ios.lon,
                y = ios.lat,
                color = factor(month(date)),
                group = week(date))) +
  guides(shape = guide_legend('Prediction Source'),
         color = guide_legend('Month'),
         size = guide_legend('Abs. Error (days)')) +
  scale_color_brewer(palette = 'Reds') +
  labs(title = paste(yrs[j], "Time Residuals"),
       x = 'Longitude',
       y = 'Latitude')

# error investigation
filter(pred_df, date > source_df[2, ]$symptom.date) %>%
  select(its.src2.lon, its.src2.lat, date) %>% 
  mutate(color = factor(month(date)), group = week(date), method = 'its')


select(out_tbl, contains('ts.t.rmse'), n, year, src.loc)

mutate(out_tbl, src.type = rep(c('imputed', 'observed'), 3))
