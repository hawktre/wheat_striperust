library(ggmap)
source('preprocess_fn_v2.R')
source('fit_twosource_fn.R')
source('plot_fn.R')
load('cdm_mapbox.RData')
map <- get_map(location = mapbox,
               maptype = 'toner-background')


## preprocess data
preprocess_out <- preprocess_fn()
DATA <- preprocess_out$data 

# select year
yrs <- 2008:2010 # possible years
j <- 1

if(j == 1){
  # source information (not year sensitive)
  src1_type <- 'observed'
  src2_type <- 'imputed'
  src3_type <- 'observed'
  
  wind_src1 <- preprocess_out$wind.s
  wind_src2 <- preprocess_out$wind.i3
  wind_src3 <- preprocess_out$wind.n
  src1_id <- filter(preprocess_out$sources.s,
                    year == yrs[j])$ID
  src2_id <- filter(preprocess_out$sources.i3,
                    year == yrs[j])$ID
  # src2_id <- NULL
  src3_id <- filter(preprocess_out$sources.n,
                    year == yrs[j])$ID
  # src3_id <- NULL
  
  
  # source information (year sensitive)
  source_df <- rbind(rename(filter(preprocess_out$sources.s, year == yrs[j]),
                            cc.long = cc.long.s,
                            cc.lat = cc.lat.s),
                     rename(filter(preprocess_out$sources.i3, year == yrs[j]),
                            cc.long = cc.long.i3,
                            cc.lat = cc.lat.i3),
                     rename(filter(preprocess_out$sources.n, year == yrs[j]),
                            cc.long = cc.long.n,
                            cc.lat = cc.lat.n))
  
  
  # subset data to year and calculate basis expansions and response variables
  sub_data <- filter(DATA, year(symptom.date) == yrs[j]) %>%
    rename(rho.src1 = rho.s,
           theta.src1 = theta.s,
           g.theta.src1 = g.theta.s,
           rho.src2 = rho.i3,
           theta.src2 = theta.i3,
           g.theta.src2 = g.theta.i3,
           rho.src3 = rho.n,
           theta.src3 = theta.n,
           g.theta.src3 = g.theta.n) %>%
    mutate(month = month(symptom.date),
           day = day(symptom.date),
           time = yday(symptom.date),
           iresponse.src1 = log(1 + rho.src1), # isotropic response for first source
           iresponse.src2 = log(1 + rho.src2), # isotropic response for second source
           iresponse.src3 = log(1 + rho.src3), # isotropic response for second source
           aresponse.src1 = log(1 + rho.src1/g.theta.src1), # anisotropic response for first source
           aresponse.src2 = log(1 + rho.src2/g.theta.src2), # anisotropic response for second source
           aresponse.src3 = log(1 + rho.src3/g.theta.src3), # anisotropic response for second source
           sinb1.src1 = sin(as.numeric(theta.src1)), # sin basis 1 for first source
           sinb2.src1 = sin(as.numeric(theta.src1) + pi/4), # sin basis 2 for first source
           sinb1.src2 = sin(as.numeric(theta.src2)), # sin basis 1 for second source
           sinb2.src2 = sin(as.numeric(theta.src2) + pi/4), # sin basis 2 for second source
           sinb1.src3 = sin(as.numeric(theta.src3)), # sin basis 1 for second source
           sinb2.src3 = sin(as.numeric(theta.src3) + pi/4), # sin basis 2 for second source
           z.reg.src1 = rho.src1/g.theta.src1, # first source regressor to initialize weights
           z.reg.src2 = rho.src2/g.theta.src2, # second source regressor to initialize weights
           z.reg.src3 = rho.src3/g.theta.src3, # second source regressor to initialize weights
           z.fit.12 = rho.src1 < rho.src2, # is first source closer than second source?
           z.fit.13 = rho.src1 < rho.src3) # is first source closer than second source?
  
  sub_data <- sub_data %>%
    group_by(ID) %>%
    sample_n(1)
  
  plot_out <- plot_fn(sub_data, wind_src1, wind_src2, wind_src3, source_df, src1_id, NULL, src3_id)
}

if(j == 2){
  # source information (not year sensitive)
  src1_type <- 'observed'
  src2_type <- 'imputed'
  src3_type <- 'observed'
  
  wind_src1 <- preprocess_out$wind.s
  wind_src2 <- preprocess_out$wind.i4
  wind_src3 <- preprocess_out$wind.n
  src1_id <- filter(preprocess_out$sources.s,
                    year == yrs[j])$ID
  src2_id <- filter(preprocess_out$sources.i4,
                    year == yrs[j])$ID
  # src2_id <- NULL
  src3_id <- filter(preprocess_out$sources.n,
                    year == yrs[j])$ID
  # src3_id <- NULL
  
  
  # source information (year sensitive)
  source_df <- rbind(rename(filter(preprocess_out$sources.s, year == yrs[j]),
                            cc.long = cc.long.s,
                            cc.lat = cc.lat.s),
                     rename(filter(preprocess_out$sources.i4, year == yrs[j]),
                            cc.long = cc.long.i4,
                            cc.lat = cc.lat.i4),
                     rename(filter(preprocess_out$sources.n, year == yrs[j]),
                            cc.long = cc.long.n,
                            cc.lat = cc.lat.n))
  
  
  # subset data to year and calculate basis expansions and response variables
  sub_data <- filter(DATA, year(symptom.date) == yrs[j]) %>%
    rename(rho.src1 = rho.s,
           theta.src1 = theta.s,
           g.theta.src1 = g.theta.s,
           rho.src2 = rho.i4,
           theta.src2 = theta.i4,
           g.theta.src2 = g.theta.i4,
           rho.src3 = rho.n,
           theta.src3 = theta.n,
           g.theta.src3 = g.theta.n) %>%
    mutate(month = month(symptom.date),
           day = day(symptom.date),
           time = yday(symptom.date),
           iresponse.src1 = log(1 + rho.src1), # isotropic response for first source
           iresponse.src2 = log(1 + rho.src2), # isotropic response for second source
           iresponse.src3 = log(1 + rho.src3), # isotropic response for second source
           aresponse.src1 = log(1 + rho.src1/g.theta.src1), # anisotropic response for first source
           aresponse.src2 = log(1 + rho.src2/g.theta.src2), # anisotropic response for second source
           aresponse.src3 = log(1 + rho.src3/g.theta.src3), # anisotropic response for second source
           sinb1.src1 = sin(as.numeric(theta.src1)), # sin basis 1 for first source
           sinb2.src1 = sin(as.numeric(theta.src1) + pi/4), # sin basis 2 for first source
           sinb1.src2 = sin(as.numeric(theta.src2)), # sin basis 1 for second source
           sinb2.src2 = sin(as.numeric(theta.src2) + pi/4), # sin basis 2 for second source
           sinb1.src3 = sin(as.numeric(theta.src3)), # sin basis 1 for second source
           sinb2.src3 = sin(as.numeric(theta.src3) + pi/4), # sin basis 2 for second source
           z.reg.src1 = rho.src1/g.theta.src1, # first source regressor to initialize weights
           z.reg.src2 = rho.src2/g.theta.src2, # second source regressor to initialize weights
           z.reg.src3 = rho.src3/g.theta.src3, # second source regressor to initialize weights
           z.fit.12 = rho.src1 < rho.src2, # is first source closer than second source?
           z.fit.13 = rho.src1 < rho.src3) # is first source closer than second source?
  
  sub_data <- sub_data %>%
    group_by(ID) %>%
    sample_n(1)
  
  plot_out <- plot_fn(sub_data, wind_src1, wind_src2, wind_src3, source_df, src1_id, NULL, src3_id)
  
}

if(j == 3){
  # source information (not year sensitive)
  src1_type <- 'observed'
  src2_type <- 'imputed'
  src3_type <- 'imputed'
  
  wind_src1 <- preprocess_out$wind.s
  wind_src2 <- preprocess_out$wind.i3
  wind_src3 <- preprocess_out$wind.i1
  src1_id <- filter(preprocess_out$sources.s,
                    year == yrs[j])$ID
  src2_id <- filter(preprocess_out$sources.i3,
                    year == yrs[j])$ID
  # src2_id <- NULL
  src3_id <- filter(preprocess_out$sources.i1,
                    year == yrs[j])$ID
  # src3_id <- NULL
  
  
  # source information (year sensitive)
  source_df <- rbind(rename(filter(preprocess_out$sources.s, year == yrs[j]),
                            cc.long = cc.long.s,
                            cc.lat = cc.lat.s),
                     rename(filter(preprocess_out$sources.i3, year == yrs[j]),
                            cc.long = cc.long.i3,
                            cc.lat = cc.lat.i3),
                     rename(filter(preprocess_out$sources.i1, year == yrs[j]),
                            cc.long = cc.long.i1,
                            cc.lat = cc.lat.i1))
  
  
  # subset data to year and calculate basis expansions and response variables
  sub_data <- filter(DATA, year(symptom.date) == yrs[j]) %>%
    rename(rho.src1 = rho.s,
           theta.src1 = theta.s,
           g.theta.src1 = g.theta.s,
           rho.src2 = rho.i3,
           theta.src2 = theta.i3,
           g.theta.src2 = g.theta.i3,
           rho.src3 = rho.i1,
           theta.src3 = theta.i1,
           g.theta.src3 = g.theta.i1) %>%
    mutate(month = month(symptom.date),
           day = day(symptom.date),
           time = yday(symptom.date),
           iresponse.src1 = log(1 + rho.src1), # isotropic response for first source
           iresponse.src2 = log(1 + rho.src2), # isotropic response for second source
           iresponse.src3 = log(1 + rho.src3), # isotropic response for second source
           aresponse.src1 = log(1 + rho.src1/g.theta.src1), # anisotropic response for first source
           aresponse.src2 = log(1 + rho.src2/g.theta.src2), # anisotropic response for second source
           aresponse.src3 = log(1 + rho.src3/g.theta.src3), # anisotropic response for second source
           sinb1.src1 = sin(as.numeric(theta.src1)), # sin basis 1 for first source
           sinb2.src1 = sin(as.numeric(theta.src1) + pi/4), # sin basis 2 for first source
           sinb1.src2 = sin(as.numeric(theta.src2)), # sin basis 1 for second source
           sinb2.src2 = sin(as.numeric(theta.src2) + pi/4), # sin basis 2 for second source
           sinb1.src3 = sin(as.numeric(theta.src3)), # sin basis 1 for second source
           sinb2.src3 = sin(as.numeric(theta.src3) + pi/4), # sin basis 2 for second source
           z.reg.src1 = rho.src1/g.theta.src1, # first source regressor to initialize weights
           z.reg.src2 = rho.src2/g.theta.src2, # second source regressor to initialize weights
           z.reg.src3 = rho.src3/g.theta.src3, # second source regressor to initialize weights
           z.fit.12 = rho.src1 < rho.src2, # is first source closer than second source?
           z.fit.13 = rho.src1 < rho.src3) # is first source closer than second source?
  
  sub_data <- sub_data %>%
    group_by(ID) %>%
    sample_n(1)
  
  plot_out <- plot_fn(sub_data, wind_src1, wind_src2, wind_src3, source_df, src1_id, NULL, NULL)
  
  
}



pred_df <- plot_out$pred
plot_df <- plot_out$plot %>%
  mutate(which.pred = if_else(str_detect(method, 'TS-SW'), 
                              which.pred.sw, which.pred.n))
plot_df[str_detect(plot_df$method, 'OS'), ]$which.pred <- 1
plot_df <- filter(plot_df, !str_detect(method, '-N'))


# filter(plot_df, str_detect(method, 'TS-SW'), which.pred == 2)

# if(j == 3){
#   idx <- as.logical(1 - (plot_df$which.pred == 2)*str_detect(plot_df$method, 'TS-SW'))
#   plot_df <- plot_df[idx, ]
# }

rmses <- group_by(plot_df, method) %>%
  summarize(rmse = sqrt(mean(error^2, na.rm = T)))

methodlabels <- paste(rmses$method, ' (RMSE = ', round(rmses$rmse, 2), ' days)', sep = '')


plot_df <- plot_df %>%
  cbind(method.fac = factor(plot_df$method, labels = methodlabels))

methodfactor <- unique(plot_df$method.fac)

# lonrange <- range(sub_data$long)
# latrange <- range(sub_data$lat)
# mapbox <- make_bbox(lonrange, latrange)

ggmap(map, darken = c(0.45, 'grey')) +
  ## plot errors
  geom_point(data = plot_df, 
             aes(x = long, 
                 y = lat, 
                 shape = factor(which.pred,
                                labels = c('Florida', 
                                           'Alternate')),
                 size = abs(error),
                 color = month(symptom.date))) +
  scale_size_binned(breaks = c(7, 14, 21),
                    range = c(1, 5),
                    labels = c('0-1 weeks',
                               '1-2 weeks',
                               '>2 weeks'),
                    guide = 'legend') +
  ## break into panels
  facet_wrap(~method.fac, nrow = 2, dir = 'h') +
  ## add source points
  geom_point(data = mutate(source_df[1, ], method.fac = methodfactor[1]),
             aes(x = cc.long, y = cc.lat),
             color = 'red',
             shape = 3) +
  geom_point(data = mutate(source_df[1, ], method.fac = methodfactor[3]),
             aes(x = cc.long, y = cc.lat),
             color = 'red',
             shape = 3) +
  # geom_point(data = mutate(source_df[c(1, 3), ], method.fac = methodfactor[2]),
  #            aes(x = cc.long, y = cc.lat),
  #            color = 'red',
  #            shape = 3) +
  # geom_point(data = mutate(source_df[c(1, 3), ], method.fac = methodfactor[5]),
  #            aes(x = cc.long, y = cc.lat),
  #            color = 'red',
  #            shape = 3) +
  geom_point(data = mutate(source_df[c(1, 2), ], method.fac = methodfactor[2]),
             aes(x = cc.long, y = cc.lat),
             color = 'red',
             shape = 3) +
  geom_point(data = mutate(source_df[c(1, 2), ], method.fac = methodfactor[4]),
             aes(x = cc.long, y = cc.lat),
             color = 'red',
             shape = 3) +
  ## add contours
  geom_path(data = mutate(filter(pred_df), #theta < pi/2 | theta > 3*pi/2),
                          method.fac = methodfactor[1]),
            aes(x = aos.lon,
                y = aos.lat,
                color = month(date),
                group = week(date))) +
  geom_path(data = mutate(filter(pred_df), #theta < pi/2 | theta > 3*pi/2),
                          method.fac = methodfactor[3]),
            aes(x = ios.lon,
                y = ios.lat,
                color = month(date),
                group = week(date))) +
  geom_path(data = mutate(filter(pred_df, date > source_df[2, ]$symptom.date),
                          method.fac = methodfactor[2]),
            aes(x = ats1.src2.lon,
                y = ats1.src2.lat,
                color = month(date),
                group = week(date)),
            linetype = 6) +
  geom_path(data = mutate(filter(pred_df,
                                 date > source_df[1, ]$symptom.date), #theta < pi/2 | theta > 3*pi/2),
                          method.fac = methodfactor[2]),
            aes(x = ats1.src1.lon,
                y = ats1.src1.lat,
                color = month(date),
                group = week(date))) +
  geom_path(data = mutate(filter(pred_df, date > source_df[2, ]$symptom.date),
                          method.fac = methodfactor[4]),
            aes(x = its1.src2.lon,
                y = its1.src2.lat,
                color = month(date),
                group = week(date)),
            linetype = 6) +
  geom_path(data = mutate(filter(pred_df, date > source_df[1, ]$symptom.date), 
                          method.fac = methodfactor[4]),
            aes(x = its1.src1.lon,
                y = its1.src1.lat,
                color = month(date),
                group = week(date)))  +
  # geom_path(data = mutate(filter(pred_df, date > source_df[3, ]$symptom.date), 
  #                         method.fac = methodfactor[2]),
  #           aes(x = ats2.src3.lon,
  #               y = ats2.src3.lat,
  #               color = month(date),
  #               group = week(date)),
  #           linetype = 6) +
  # geom_path(data = mutate(filter(pred_df, date > source_df[1, ]$symptom.date), 
  #                         method.fac = methodfactor[2]),
  #           aes(x = ats2.src1.lon,
  #               y = ats2.src1.lat,
  #               color = month(date),
  #               group = week(date))) +
  # geom_path(data = mutate(filter(pred_df, date > source_df[3, ]$symptom.date),
  #                         method.fac = methodfactor[5]),
  #           aes(x = its2.src3.lon,
  #               y = its2.src3.lat,
  #               color = month(date),
  #               group = week(date)),
  #           linetype = 6) +
  # geom_path(data = mutate(filter(pred_df, date > source_df[1, ]$symptom.date), 
  #                         method.fac = methodfactor[5]),
  #           aes(x = its2.src1.lon,
  #               y = its2.src1.lat,
  #               color = month(date),
  #               group = week(date))) +
  ## adjust legends
  guides(shape = guide_legend('Prediction Source',
                              override.aes = list(size = 3)),
         linetype = guide_legend('Prediction Source'),
         color = guide_legend('Month'),
         size = guide_legend('Prediction Error')) +
  scale_color_gradientn(colors = c('yellow', 'orange', 'red', 'black'),
                        breaks = 1:7,
                        labels = month.name[1:7]) +
  labs(x = 'Longitude',
       y = 'Latitude') +
  theme(text = element_text(size = 10))

## save as 800x800
ggsave(filename = 'fig3_slideversion.tiff',
       width = 6,
       height = 5,
       units = 'in',
       dpi = 300)

