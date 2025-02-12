source('preprocess_fn.R')
source('fit_twosource_fn.R')

seed <- 42820
yrs <- 2008:2016 # possible years
j <- 1 # which year

## preprocess dat
preprocess_out <- preprocess_fn()
DATA <- preprocess_out$data 


# subset data to year and calculate basis expansions and response variables
sub_data <- filter(DATA, year(symptom.date) == yrs[j]) %>%
  mutate(month = month(symptom.date),
         day = day(symptom.date),
         time = yday(symptom.date),
         iresponse.s = log(1 + rho.s), # isotropic response for S source
         iresponse.n = log(1 + rho.n), # isotropic response for N source
         iresponse.w = log(1 + rho.w), # isotropic response for W source
         iresponse.i = log(1 + rho.i), # isotropic response for I(mputed) source
         aresponse.s = log(1 + rho.s/g.theta.s), # anisotropic response for S source
         aresponse.n = log(1 + rho.n/g.theta.n), # anisotropic response for N source
         aresponse.w = log(1 + rho.w/g.theta.w), # anisotropic response for W source
         aresponse.i = log(1 + rho.i/g.theta.i), # anisotropic response for I(mputed) source
         sinb1.s = sin(as.numeric(theta.s)), # sin basis 1 for S source
         sinb2.s = sin(as.numeric(theta.s) + pi/4), # sin basis 2 for S source
         sinb1.n = sin(as.numeric(theta.n)), # sin basis 1 for N source
         sinb2.n = sin(as.numeric(theta.n) + pi/4), # sin basis 2 for N source
         sinb1.w = sin(as.numeric(theta.w)), # sin basis 1 for W source
         sinb2.w = sin(as.numeric(theta.w) + pi/4), # sin basis 2 for W source
         sinb1.i = sin(as.numeric(theta.i)), # sin basis 1 for Imputed source
         sinb2.i = sin(as.numeric(theta.i) + pi/4), # sin basis 2 for Imputed source
         z.reg.n = rho.n/g.theta.n, # regressor N to initialize weights
         z.reg.s = rho.s/g.theta.s, # regressor S to initialize weights
         z.reg.w = rho.w/g.theta.w, # regressor W to initialize weights
         z.reg.i = rho.i/g.theta.i, # regressor Imputed to initialize weights
         z.fit.sn = rho.s < rho.n, # is S source closer than N source?
         z.fit.sw = rho.s < rho.w, # is S source closer than W source?
         z.fit.si = rho.s < rho.i) # is S source closer than Imputed source?

sub_data <- sub_data %>%
  group_by(ID) %>%
  sample_n(1)

# determine source ID's
s_id <- filter(preprocess_out$sources.s,
               year == yrs[j])$ID
n_id <- filter(preprocess_out$sources.n,
               year == yrs[j])$ID
w_id <- filter(preprocess_out$sources.w,
               year == yrs[j])$ID

# fit IOS model
fit_ios <- lm(iresponse.s ~ time, 
              data = filter(sub_data, ID != s_id))
summary(fit_ios)

# fit AOS model
fit_aos <- lm(aresponse.s ~ time + sinb1.s + sinb2.s, 
              data = filter(sub_data, ID != s_id))
summary(fit_aos)

# fit ITS model
source1_its <- filter(sub_data, !(ID %in% c(s_id, w_id))) %>%
  rename(rho = rho.s, 
         theta = theta.s,
         g.theta = g.theta.s,
         sinb1 = sinb1.s, 
         sinb2 = sinb2.s, 
         response = iresponse.s)
source2_its <- filter(sub_data, !(ID %in% c(s_id, w_id))) %>%
  rename(rho = rho.w, 
         theta = theta.w,
         g.theta = g.theta.w,
         sinb1 = sinb1.w, 
         sinb2 = sinb2.w, 
         response = iresponse.w)
initdata_its <- filter(sub_data, !(ID %in% c(s_id, w_id))) %>%
  rename(z.fit = z.fit.sw,
         z.reg.1 = rho.s,
         z.reg.2 = rho.w)

fit_its_out <- fit_ts(source1_its, 
                      source2_its, 
                      initdata_its, 
                      isotropy = T)
fit_its_s <- fit_its_out$fit1
fit_its_w <- fit_its_out$fit2
pr_z_its <- fit_its_out$prz

summary(fit_its_s)
summary(fit_its_w)

# fit ATS model
source1_ats <- filter(sub_data, !(ID %in% c(s_id, w_id))) %>%
  rename(rho = rho.s, 
         theta = theta.s,
         g.theta = g.theta.s,
         sinb1 = sinb1.s, 
         sinb2 = sinb2.s, 
         response = aresponse.s)
source2_ats <- filter(sub_data, !(ID %in% c(s_id, w_id))) %>%
  rename(rho = rho.w, 
         theta = theta.w,
         g.theta = g.theta.w,
         sinb1 = sinb1.w, 
         sinb2 = sinb2.w, 
         response = aresponse.w)
initdata_ats <- filter(sub_data, !(ID %in% c(s_id, w_id))) %>%
  rename(z.fit = z.fit.sw,
         z.reg.1 = z.reg.s,
         z.reg.2 = z.reg.w)

fit_ats_out <- fit_ts(source1_ats, 
                      source2_ats, 
                      initdata_ats, 
                      isotropy = F)
fit_ats_s <- fit_ats_out$fit1
fit_ats_w <- fit_ats_out$fit2
pr_z_ats <- fit_ats_out$prz

summary(fit_ats_s)
summary(fit_ats_w)

# extract estimates from each model 
betahat_ios <- coef(fit_ios)
betahat_aos <- coef(fit_aos)
betahat_ats1 <- coef(fit_ats_s)
betahat_ats2 <- coef(fit_ats_w)
betahat_its1 <- coef(fit_its_s)
betahat_its2 <- coef(fit_its_w)

# predictions of time of occurrence according to each model
sub_data <- ungroup(sub_data) %>%
  mutate(          t.hat.ats.s = (aresponse.s - betahat_ats1[1] - betahat_ats1[3]*sinb1.s - betahat_ats1[4]*sinb2.s)/betahat_ats1[2],
                   t.hat.ats.w = (aresponse.w - betahat_ats2[1] - betahat_ats2[3]*sinb1.w - betahat_ats2[4]*sinb2.w)/betahat_ats2[2],
                   t.hat.aos = (aresponse.s - betahat_aos[1] - betahat_aos[3]*sinb1.s - betahat_aos[4]*sinb2.s)/betahat_aos[2],
                   t.hat.ios = (iresponse.s - betahat_ios[1])/betahat_ios[2],
                   t.hat.its.s = (iresponse.s - betahat_its1[1])/betahat_its1[2],
                   t.hat.its.w = (iresponse.w - betahat_its2[1])/betahat_its2[2],
                   rho.s.hat.ios = (exp(predict(fit_ios, sub_data)) - 1),
                   rho.s.hat.aos = (exp(predict(fit_aos, sub_data)) - 1)*g.theta.s,
                   rho.s.hat.its = (exp(predict(fit_its_s, sub_data)) - 1),
                   rho.w.hat.its = (exp(predict(fit_its_w, sub_data)) - 1),
                   rho.s.hat.ats = (exp(predict(fit_ats_s, 
                                                rename(sub_data,
                                                       sinb1 = sinb1.s,
                                                       sinb2 = sinb2.s))) - 1)*g.theta.s,
                   rho.w.hat.ats = (exp(predict(fit_ats_w, 
                                                rename(sub_data,
                                                       sinb1 = sinb1.w,
                                                       sinb2 = sinb2.w))) - 1)*g.theta.w)

# round p_i's and choose predictions
which_pred_ats <- data.frame(which.pred.ats = round(pr_z_ats, 0) + (round(pr_z_ats, 0) == 0)*2,
                             ID = source1_ats$ID)
which_pred_ats2 <- data.frame(which.pred.ats.t = replace_na(sub_data$t.hat.ats.s < sub_data$t.hat.ats.w, T),
                              which.pred.ats.rho = replace_na((sub_data$rho.s.hat.ats - sub_data$rho.s) < (sub_data$rho.w.hat.ats - sub_data$rho.w), T),
                              ID = sub_data$ID)
which_pred_its <- data.frame(which.pred.its = round(pr_z_its, 0) + (round(pr_z_its, 0) == 0)*2,
                             ID = source1_its$ID)
which_pred_its2 <- data.frame(which.pred.its.t = replace_na(sub_data$t.hat.its.s < sub_data$t.hat.its.w, T),
                              which.pred.its.rho = replace_na((sub_data$rho.s.hat.its - sub_data$rho.s) < (sub_data$rho.w.hat.its - sub_data$rho.w), T),
                              ID = sub_data$ID)


sub_data_pred <- sub_data %>%
  as_tibble() %>%
  merge(which_pred_ats2, by = "ID", all = F) %>%
  merge(which_pred_its2, by = "ID", all = F) %>%
  # merge(which_pred_ats, by = "ID", all = T) %>%
  # merge(which_pred_its, by = "ID", all = T) %>%
  # mutate(which.pred.ats = replace_na(which.pred.ats, 1),
  #        which.pred.its = replace_na(which.pred.its, 1)) %>%
  mutate(t.resid.aos = t.hat.aos - time,
         t.resid.ios = t.hat.ios - time,
         t.resid.ats = if_else(which.pred.ats.t == T,
                               t.hat.ats.s,
                               t.hat.ats.w) - time,
         t.resid.its = if_else(which.pred.its.t == T,
                               t.hat.its.s,
                               t.hat.its.w) - time,
         rho.resid.its = if_else(which.pred.its.rho == T,
                                 rho.s.hat.its - rho.s,
                                 rho.w.hat.its - rho.w),
         rho.resid.ats = if_else(which.pred.ats.rho == T,
                                 rho.s.hat.ats - rho.s,
                                 rho.w.hat.ats - rho.w),
         # t.resid.ats = if_else(which.pred.ats == 1,
         #                       t.hat.ats.s,
         #                       t.hat.ats.n) - time,
         # t.resid.its = if_else(which.pred.its == 1,
         #                       t.hat.its.s,
         #                       t.hat.its.n) - time,
         # rho.resid.its = if_else(which.pred.its == 1,
         #                         rho.s.hat.its - rho.s,
         #                         rho.n.hat.its - rho.n),
         # rho.resid.ats = if_else(which.pred.ats == 1,
         #                         rho.s.hat.ats - rho.s,
         #                         rho.n.hat.ats - rho.n),
         rho.resid.ios = rho.s.hat.ios - rho.s,
         rho.resid.aos = rho.s.hat.aos - rho.s) %>%
  select(symptom.date, 
         lat, 
         long, 
         ID,
         t.resid.ats,
         t.resid.its,
         t.resid.aos,
         t.resid.ios,
         rho.resid.ats,
         rho.resid.its,
         rho.resid.ios,
         rho.resid.aos,
         # which.pred.ats,
         # which.pred.its,
         which.pred.ats.t,
         which.pred.its.t,
         which.pred.its.rho,
         which.pred.ats.rho,
         rho.s,
         theta.s,
         rho.w,
         theta.w,
         host)

# replace source predictions by NA
names(sub_data_pred[, 5:12])
names(sub_data_pred[, c(5:6, 9:10)])
sub_data_pred[sub_data_pred$ID == s_id, 5:12] <- NA
sub_data_pred[sub_data_pred$ID == w_id, c(5:6, 9:10)] <- NA

############

sub_data_pred %>%
  summarize(ats.t.rmse = sqrt(mean(t.resid.ats^2, na.rm = T)),
            its.t.rmse = sqrt(mean(t.resid.its^2, na.rm = T)),
            aos.t.rmse = sqrt(mean(t.resid.aos^2, na.rm = T)),
            ios.t.rmse = sqrt(mean(t.resid.ios^2, na.rm = T)),
            ats.rho.rmse = sqrt(mean(rho.resid.ats^2, na.rm = T)),
            its.rho.rmse = sqrt(mean(rho.resid.its^2, na.rm = T)),
            aos.rho.rmse = sqrt(mean(rho.resid.aos^2, na.rm = T)),
            ios.rho.rmse = sqrt(mean(rho.resid.ios^2, na.rm = T)),
            n = n())

library(ggmap)
lonrange <- range(sub_data_pred$long)
latrange <- range(sub_data_pred$lat)
mapbox <- make_bbox(lonrange, latrange)
map <- get_map(location = mapbox,
               maptype = 'toner-background')


# residuals 

plot_df <- sub_data_pred %>%
  select(c(1:4, 5:8, 13:14)) %>%
  rename(aos = t.resid.aos,
         ats = t.resid.ats,
         ios = t.resid.ios,
         its = t.resid.its) %>%
  # select(c(1:4, 9:12, 15:16)) %>%
  # rename(aos = rho.resid.aos,
  #        ats = rho.resid.ats,
  #        ios = rho.resid.ios,
  #        its = rho.resid.its) %>%
  # select(c(1:4, 9:12, 13:14)) %>%
  gather('method', 'error', 5:8) %>%
  group_by(method, ID) %>%
  summarize(error = mean(error),
            symptom.date = unique(symptom.date),
            which.pred = if_else(str_detect(method, 'ats'),
                                 unique(which.pred.ats.t),
                                 unique(which.pred.its.t)),
            # which.pred = if_else(str_detect(method, 'ats'),
            #                      unique(which.pred.ats.rho),
            #                      unique(which.pred.its.rho)),
            # which.pred = if_else(str_detect(method, 'ats'),
            #                      unique(which.pred.ats),
            #                      unique(which.pred.its)),
            lat = unique(lat),
            long = unique(long)) %>%
  mutate(which.pred = if_else(str_detect(method, 'os'), T, which.pred))
# mutate(which.pred = if_else(str_detect(method, 'os'), 1, which.pred))

# sources

source_df <- rbind(rename(filter(preprocess_out$sources.s, year == yrs[j]),
                          cc.long = cc.long.s,
                          cc.lat = cc.lat.s),
                   rename(filter(preprocess_out$sources.w, year == yrs[j]),
                          cc.long = cc.long.w,
                          cc.lat = cc.lat.w))

# contours

theta_grid <- seq(0, 2*pi, length = 100)
week_grid <- seq(week(ymd('2009/03/01')), 
                 week(ymd('2009/08/01')),
                 by = 3)
date_grid <- rep(ymd('2009/04/15'), length(week_grid))
week(date_grid) <- week_grid
pred_df <- expand.grid(theta = theta_grid, 
                       date = date_grid) %>%
  mutate(time = yday(date),
         sinb1.s = sin(theta),
         sinb2.s = sin(theta + pi/4),
         sinb1 = sin(theta),
         sinb2 = sin(theta + pi/4))



pred_df <- pred_df %>%
  mutate(pred.ats.s = predict(fit_ats_s, newdata = pred_df),
         pred.ats.w = predict(fit_ats_w, newdata = pred_df),
         pred.its.s = predict(fit_its_s, newdata = pred_df),
         pred.its.w = predict(fit_its_w, newdata = pred_df),
         pred.aos = predict(fit_aos, newdata = pred_df),
         pred.ios = predict(fit_ios, newdata = pred_df),
         g.theta.s = density.circular(x = filter(preprocess_out$wind.s, year == yrs[j])$direction,
                                      z = pred_df$theta, 
                                      bw = 50)$y,
         g.theta.w = density.circular(x = filter(preprocess_out$wind.w, year == yrs[j])$direction,
                                      z = pred_df$theta, 
                                      bw = 50)$y) %>%
  mutate(rho.ats.s = (exp(pred.ats.s) - 1)*g.theta.s,
         rho.ats.w = (exp(pred.ats.w) - 1)*g.theta.w,
         rho.its.s = (exp(pred.its.s) - 1),
         rho.its.w = (exp(pred.its.w) - 1),
         rho.aos = (exp(pred.aos) - 1)*g.theta.s,
         rho.ios = exp(pred.ios) - 1) %>%
  select(theta, date, time, 
         rho.ats.s, rho.ats.w, 
         rho.its.s, rho.its.w, 
         rho.aos, 
         rho.ios)

pred_df <- pred_df %>%
  cbind(ats.s = destPoint(p = source_df[1, 2:3],
                          b = pred_df$theta*360/(2*pi),
                          d = 1000*pred_df$rho.ats.s),
        ats.w = destPoint(p = source_df[2, 2:3],
                          b = pred_df$theta*360/(2*pi),
                          d = 1000*pred_df$rho.ats.w),
        its.s = destPoint(p = source_df[1, 2:3],
                          b = pred_df$theta*360/(2*pi),
                          d = 1000*pred_df$rho.its.s),
        its.w = destPoint(p = source_df[2, 2:3],
                          b = pred_df$theta*360/(2*pi),
                          d = 1000*pred_df$rho.its.w),
        aos = destPoint(p = source_df[1, 2:3],
                        b = pred_df$theta*360/(2*pi),
                        d = 1000*pred_df$rho.aos),
        ios = destPoint(p = source_df[1, 2:3],
                        b = pred_df$theta*360/(2*pi),
                        d = 1000*pred_df$rho.ios))



ggmap(map, darken = c(0.8, 'white')) +
  geom_point(data = plot_df, 
             aes(x = long, 
                 y = lat, 
                 shape = factor(which.pred,
                                labels = c('W', 'S')),
                 # shape = factor(which.pred,
                 #                labels = c('S', 'N')),
                 size = abs(error)), 
             alpha = 0.5) +
  facet_wrap(~method) +
  # scale_size_continuous(breaks = c(500, 1000, 2000),
  #                       limits = c(1, 4000),
  #                       range = c(1, 8)) +
  scale_size_continuous(breaks = c(14, 28, 59),
                        limits = c(1, 250),
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
  geom_path(data = mutate(filter(pred_df, theta < pi/2 | theta > 3*pi/2),
                          method = 'aos'),
            aes(x = aos.lon,
                y = aos.lat,
                color = month(date),
                group = week(date))) +
  geom_path(data = mutate(pred_df, method = 'ats'),
            aes(x = ats.w.lon,
                y = ats.w.lat,
                color = month(date),
                group = week(date))) +
  geom_path(data = mutate(pred_df, method = 'ats'),
            aes(x = ats.s.lon,
                y = ats.s.lat,
                color = month(date),
                group = week(date))) +
  geom_path(data = mutate(pred_df, method = 'its'),
            aes(x = its.w.lon,
                y = its.w.lat,
                color = month(date),
                group = week(date))) +
  geom_path(data = mutate(pred_df, method = 'its'),
            aes(x = its.s.lon,
                y = its.s.lat,
                color = month(date),
                group = week(date))) +
  geom_path(data = mutate(pred_df, method = 'ios'),
            aes(x = ios.lon,
                y = ios.lat,
                color = month(date),
                group = week(date))) +
  guides(shape = guide_legend('Prediction Source'),
         color = guide_colorbar('Month'),
         size = guide_legend('Abs. Error')) +
  labs(title = paste(yrs[j], "Time Residuals"),
       x = 'Longitude',
       y = 'Latitude')

