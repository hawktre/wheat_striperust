plot_fn <- function(sub_data, 
                    wind_src1, wind_src2, wind_src3, 
                    source_df, src1_id, src2_id, src3_id){

# fit IOS model
fit_ios <- lm(iresponse.src1 ~ time, 
              data = filter(sub_data, ID != src1_id))
summary(fit_ios)

# fit AOS model
fit_aos <- lm(aresponse.src1 ~ time + sinb1.src1 + sinb2.src1, 
              data = filter(sub_data, ID != src1_id))
summary(fit_aos)

# fit first ITS model
source1_its1 <- filter(sub_data, !(ID %in% c(src1_id, src2_id))) %>%
  rename(rho = rho.src1, 
         theta = theta.src1,
         g.theta = g.theta.src1,
         sinb1 = sinb1.src1, 
         sinb2 = sinb2.src1, 
         response = iresponse.src1)
source2_its1 <- filter(sub_data, !(ID %in% c(src1_id, src2_id))) %>%
  rename(rho = rho.src2, 
         theta = theta.src2,
         g.theta = g.theta.src2,
         sinb1 = sinb1.src2, 
         sinb2 = sinb2.src2, 
         response = iresponse.src2)
initdata_its1 <- filter(sub_data, !(ID %in% c(src1_id, src2_id))) %>%
  rename(z.reg.1 = rho.src1,
         z.reg.2 = rho.src2,
         z.fit = z.fit.12)

fit_its1_out <- fit_ts(source1_its1, 
                      source2_its1, 
                      initdata_its1, 
                      isotropy = T)
fit_its1_src1 <- fit_its1_out$fit1
fit_its1_src2 <- fit_its1_out$fit2
pr_z_its1 <- fit_its1_out$prz

# fit second ITS model
source1_its2 <- filter(sub_data, !(ID %in% c(src1_id, src3_id))) %>%
  rename(rho = rho.src1, 
         theta = theta.src1,
         g.theta = g.theta.src1,
         sinb1 = sinb1.src1, 
         sinb2 = sinb2.src1, 
         response = iresponse.src1)
source2_its2 <- filter(sub_data, !(ID %in% c(src1_id, src3_id))) %>%
  rename(rho = rho.src3, 
         theta = theta.src3,
         g.theta = g.theta.src3,
         sinb1 = sinb1.src3, 
         sinb2 = sinb2.src3, 
         response = iresponse.src3)
initdata_its2 <- filter(sub_data, !(ID %in% c(src1_id, src3_id))) %>%
  rename(z.reg.1 = rho.src1,
         z.reg.2 = rho.src3,
         z.fit = z.fit.13)

fit_its2_out <- fit_ts(source1_its2, 
                       source2_its2, 
                       initdata_its2, 
                       isotropy = T)
fit_its2_src1 <- fit_its2_out$fit1
fit_its2_src2 <- fit_its2_out$fit2
pr_z_its2 <- fit_its2_out$prz


# fit first ATS model
source1_ats1 <- filter(sub_data, !(ID %in% c(src1_id, src2_id))) %>%
  rename(rho = rho.src1, 
         theta = theta.src1,
         g.theta = g.theta.src1,
         sinb1 = sinb1.src1, 
         sinb2 = sinb2.src1, 
         response = aresponse.src1)
source2_ats1 <- filter(sub_data, !(ID %in% c(src1_id, src2_id))) %>%
  rename(rho = rho.src2, 
         theta = theta.src2,
         g.theta = g.theta.src2,
         sinb1 = sinb1.src2, 
         sinb2 = sinb2.src2, 
         response = aresponse.src2)
initdata_ats1 <- filter(sub_data, !(ID %in% c(src1_id, src2_id))) %>%
  rename(z.reg.1 = z.reg.src1,
         z.reg.2 = z.reg.src2,
         z.fit = z.fit.12)

fit_ats1_out <- fit_ts(source1_ats1, 
                      source2_ats1, 
                      initdata_ats1, 
                      isotropy = F)
fit_ats1_src1 <- fit_ats1_out$fit1
fit_ats1_src2 <- fit_ats1_out$fit2
pr_z_ats1 <- fit_ats1_out$prz

# fit second ATS model
source1_ats2 <- filter(sub_data, !(ID %in% c(src1_id, src3_id))) %>%
  rename(rho = rho.src1, 
         theta = theta.src1,
         g.theta = g.theta.src1,
         sinb1 = sinb1.src1, 
         sinb2 = sinb2.src1, 
         response = aresponse.src1)
source2_ats2 <- filter(sub_data, !(ID %in% c(src1_id, src3_id))) %>%
  rename(rho = rho.src3, 
         theta = theta.src3,
         g.theta = g.theta.src3,
         sinb1 = sinb1.src3, 
         sinb2 = sinb2.src3, 
         response = aresponse.src3)
initdata_ats2 <- filter(sub_data, !(ID %in% c(src1_id, src3_id))) %>%
  rename(z.reg.1 = z.reg.src1,
         z.reg.2 = z.reg.src3,
         z.fit = z.fit.13)

fit_ats2_out <- fit_ts(source1_ats2, 
                       source2_ats2, 
                       initdata_ats2, 
                       isotropy = F)
fit_ats2_src1 <- fit_ats2_out$fit1
fit_ats2_src2 <- fit_ats2_out$fit2
pr_z_ats2 <- fit_ats2_out$prz


# extract estimates from each model 
betahat_ios <- coef(fit_ios)
betahat_aos <- coef(fit_aos)
betahat_ats11 <- coef(fit_ats1_src1)
betahat_ats12 <- coef(fit_ats1_src2)
betahat_its11 <- coef(fit_its1_src1)
betahat_its12 <- coef(fit_its1_src2)
betahat_ats21 <- coef(fit_ats2_src1)
betahat_ats22 <- coef(fit_ats2_src2)
betahat_its21 <- coef(fit_its2_src1)
betahat_its22 <- coef(fit_its2_src2)

# predictions of time of occurrence according to each model
sub_data <- ungroup(sub_data) %>%
  mutate(t.hat.ats1.src1 = (aresponse.src1 - betahat_ats11[1] - betahat_ats11[3]*sinb1.src1 - betahat_ats11[4]*sinb2.src1)/betahat_ats11[2],
         t.hat.ats1.src2 = (aresponse.src2 - betahat_ats12[1] - betahat_ats12[3]*sinb1.src2 - betahat_ats12[4]*sinb2.src2)/betahat_ats12[2],
         t.hat.ats2.src1 = (aresponse.src1 - betahat_ats21[1] - betahat_ats21[3]*sinb1.src1 - betahat_ats21[4]*sinb2.src1)/betahat_ats21[2],
         t.hat.ats2.src3 = (aresponse.src3 - betahat_ats22[1] - betahat_ats22[3]*sinb1.src3 - betahat_ats22[4]*sinb2.src3)/betahat_ats22[2],
         t.hat.aos = (aresponse.src1 - betahat_aos[1] - betahat_aos[3]*sinb1.src1 - betahat_aos[4]*sinb2.src1)/betahat_aos[2],
         t.hat.ios = (iresponse.src1 - betahat_ios[1])/betahat_ios[2],
         t.hat.its1.src1 = (iresponse.src1 - betahat_its11[1])/betahat_its11[2],
         t.hat.its1.src2 = (iresponse.src2 - betahat_its12[1])/betahat_its12[2],
         t.hat.its2.src1 = (iresponse.src1 - betahat_its21[1])/betahat_its21[2],
         t.hat.its2.src3 = (iresponse.src3 - betahat_its22[1])/betahat_its22[2],
         rho.src1.hat.ios = (exp(predict(fit_ios, sub_data)) - 1),
         rho.src1.hat.aos = (exp(predict(fit_aos, sub_data)) - 1)*g.theta.src1,
         rho.src1.hat.its1 = (exp(predict(fit_its1_src1, sub_data)) - 1),
         rho.src2.hat.its1 = (exp(predict(fit_its1_src2, sub_data)) - 1),
         rho.src1.hat.ats1 = (exp(predict(fit_ats1_src1,
                                         rename(sub_data,
                                                sinb1 = sinb1.src1,
                                                sinb2 = sinb2.src1))) - 1)*g.theta.src1,
         rho.src2.hat.ats1 = (exp(predict(fit_ats1_src2,
                                         rename(sub_data,
                                                sinb1 = sinb1.src2,
                                                sinb2 = sinb2.src2))) - 1)*g.theta.src2,
         rho.src1.hat.its2 = (exp(predict(fit_its2_src1, sub_data)) - 1),
         rho.src3.hat.its2 = (exp(predict(fit_its2_src2, sub_data)) - 1),
         rho.src1.hat.ats2 = (exp(predict(fit_ats2_src1,
                                          rename(sub_data,
                                                 sinb1 = sinb1.src1,
                                                 sinb2 = sinb2.src1))) - 1)*g.theta.src1,
         rho.src3.hat.ats2 = (exp(predict(fit_ats2_src2,
                                          rename(sub_data,
                                                 sinb1 = sinb1.src3,
                                                 sinb2 = sinb2.src3))) - 1)*g.theta.src3)

# round p_i's and choose predictions
which_pred_ats1 <- data.frame(which.pred.ats1 = round(pr_z_ats1, 0) + (round(pr_z_ats1, 0) == 0)*2,
                             ID = source1_ats1$ID)
# which_pred_ats2 <- data.frame(which.pred.ats.t = replace_na(sub_data$t.hat.ats.src1 < sub_data$t.hat.ats.src2, T),
#                               which.pred.ats.rho = replace_na((sub_data$rho.src1.hat.ats - sub_data$rho.src1) < (sub_data$rho.src2.hat.ats - sub_data$rho.src2), T),
#                               ID = sub_data$ID)
which_pred_its1 <- data.frame(which.pred.its1 = round(pr_z_its1, 0) + (round(pr_z_its1, 0) == 0)*2,
                             ID = source1_its1$ID)
# which_pred_its2 <- data.frame(which.pred.its.t = replace_na(sub_data$t.hat.its.src1 < sub_data$t.hat.its.src2, T),
#                               which.pred.its.rho = replace_na((sub_data$rho.src1.hat.its - sub_data$rho.src1) < (sub_data$rho.src2.hat.its - sub_data$rho.src2), T),
#                               ID = sub_data$ID)
which_pred_ats2 <- data.frame(which.pred.ats2 = round(pr_z_ats2, 0) + (round(pr_z_ats2, 0) == 0)*2,
                              ID = source1_ats2$ID)
which_pred_its2 <- data.frame(which.pred.its2 = round(pr_z_its2, 0) + (round(pr_z_its2, 0) == 0)*2,
                              ID = source1_its2$ID)


sub_data_pred <- sub_data %>%
  as_tibble() %>%
  # merge(which_pred_ats2, by = "ID", all = T) %>%
  # merge(which_pred_its2, by = "ID", all = T) %>%
  merge(which_pred_ats1, by = "ID", all = T) %>%
  merge(which_pred_its1, by = "ID", all = T) %>%
  merge(which_pred_ats2, by = "ID", all = T) %>%
  merge(which_pred_its2, by = "ID", all = T) %>%
  mutate(which.pred.ats1 = replace_na(which.pred.ats1, 1),
         which.pred.its1 = replace_na(which.pred.its1, 1),
         which.pred.ats2 = replace_na(which.pred.ats2, 1),
         which.pred.its2 = replace_na(which.pred.its2, 1)) %>%
  mutate(t.resid.aos = t.hat.aos - time,
         t.resid.ios = t.hat.ios - time,
         # t.resid.ats = if_else(which.pred.ats.t == T,
         #                       t.hat.ats.src1,
         #                       t.hat.ats.src2) - time,
         # t.resid.its = if_else(which.pred.its.t == T,
         #                       t.hat.its.src1,
         #                       t.hat.its.src2) - time,
         # rho.resid.its = if_else(which.pred.its.rho == T,
         #                         rho.src1.hat.its - rho.src1,
         #                         rho.src2.hat.its - rho.src2),
         # rho.resid.ats = if_else(which.pred.ats.rho == T,
         #                         rho.src1.hat.ats - rho.src1,
         #                         rho.src2.hat.ats - rho.src2),
         t.resid.ats1 = if_else(which.pred.ats1 == 1,
                               t.hat.ats1.src1,
                               t.hat.ats1.src2) - time,
         t.resid.its1 = if_else(which.pred.its1 == 1,
                               t.hat.its1.src1,
                               t.hat.its1.src2) - time,
         rho.resid.its1 = if_else(which.pred.its1 == 1,
                                 rho.src1.hat.its1 - rho.src1,
                                 rho.src2.hat.its1 - rho.src2),
         rho.resid.ats1 = if_else(which.pred.ats1 == 1,
                                 rho.src1.hat.ats1 - rho.src1,
                                 rho.src2.hat.ats1 - rho.src2),
         t.resid.ats2 = if_else(which.pred.ats2 == 1,
                                t.hat.ats2.src1,
                                t.hat.ats2.src3) - time,
         t.resid.its2 = if_else(which.pred.its2 == 1,
                                t.hat.its2.src1,
                                t.hat.its2.src3) - time,
         rho.resid.its2 = if_else(which.pred.its2 == 1,
                                  rho.src1.hat.its2 - rho.src1,
                                  rho.src3.hat.its2 - rho.src3),
         rho.resid.ats2 = if_else(which.pred.ats2 == 1,
                                  rho.src1.hat.ats2 - rho.src1,
                                  rho.src3.hat.ats2 - rho.src3),
         rho.resid.ios = rho.src1.hat.ios - rho.src1,
         rho.resid.aos = rho.src1.hat.aos - rho.src1) %>%
  select(symptom.date, 
         lat, 
         long, 
         ID,
         t.resid.ats1,
         t.resid.its1,
         t.resid.ats2,
         t.resid.its2,
         t.resid.aos,
         t.resid.ios,
         rho.resid.ats1,
         rho.resid.its1,
         rho.resid.ats2,
         rho.resid.its2,
         rho.resid.ios,
         rho.resid.aos,
         which.pred.ats1,
         which.pred.its1,
         which.pred.ats2,
         which.pred.its2,
         # which.pred.ats.t,
         # which.pred.its.t,
         # which.pred.its.rho,
         # which.pred.ats.rho,
         rho.src1,
         theta.src1,
         rho.src2,
         theta.src2,
         rho.src3,
         theta.src3,
         host)

# replace source predictions by NA
names(sub_data_pred[, 5:16])
names(sub_data_pred[, c(5:6, 11:12)])
names(sub_data_pred[, c(7:8, 13:14)])
sub_data_pred[sub_data_pred$ID == src1_id, 5:16] <- NA
sub_data_pred[sub_data_pred$ID == src2_id, c(5:6, 11:12)] <- NA
sub_data_pred[sub_data_pred$ID == src3_id, c(7:8, 13:14)] <- NA

############

# out_tbl <- sub_data_pred %>%
#   summarize(ats1.t.rmse = sqrt(mean(t.resid.ats1^2, na.rm = T)),
#             its1.t.rmse = sqrt(mean(t.resid.its1^2, na.rm = T)),
#             ats2.t.rmse = sqrt(mean(t.resid.ats2^2, na.rm = T)),
#             its2.t.rmse = sqrt(mean(t.resid.its2^2, na.rm = T)),
#             aos.t.rmse = sqrt(mean(t.resid.aos^2, na.rm = T)),
#             ios.t.rmse = sqrt(mean(t.resid.ios^2, na.rm = T)),
#             # ats.rho.rmse = sqrt(mean(rho.resid.ats^2, na.rm = T)),
#             # its.rho.rmse = sqrt(mean(rho.resid.its^2, na.rm = T)),
#             # aos.rho.rmse = sqrt(mean(rho.resid.aos^2, na.rm = T)),
#             # ios.rho.rmse = sqrt(mean(rho.resid.ios^2, na.rm = T)),
#             n = n()) %>%
#   mutate(year = yrs[j], 
#          src2.state = source_df[2, ]$state, 
#          src2.county = source_df[2, ]$county,
#          src2.date = source_df[2, ]$symptom.date,
#          src2.type = src_type)

# residuals 
plot_df <- sub_data_pred %>%
  select(c(1:4, 5:10, 17:20)) %>%
  rename(AOS = t.resid.aos,
         `ATS-SW` = t.resid.ats1,
         `ATS-N` = t.resid.ats2,
         IOS = t.resid.ios,
         `ITS-SW` = t.resid.its1,
         `ITS-N` = t.resid.its2) %>%
  # select(c(1:4, 9:12, 15:16)) %>%
  # rename(aos = rho.resid.aos,
  #        ats = rho.resid.ats,
  #        ios = rho.resid.ios,
  #        its = rho.resid.its) %>%
  # select(c(1:4, 9:12, 13:14)) %>%
  gather('method', 'error', 5:10) %>%
  group_by(method, ID) %>%
  summarize(error = mean(error, na.rm = T),
            symptom.date = unique(symptom.date),
            # which.pred = if_else(str_detect(method, 'ats'),
            #                      unique(which.pred.ats.t),
            #                      unique(which.pred.its.t)),
            # which.pred = if_else(str_detect(method, 'ats'),
            #                      unique(which.pred.ats.rho),
            #                      unique(which.pred.its.rho)),
            which.pred.sw = if_else(str_detect(method, 'ATS-SW'),
                                 unique(which.pred.ats1),
                                 unique(which.pred.its1)),
            which.pred.n = if_else(str_detect(method, 'ATS-N'),
                                    unique(which.pred.ats2),
                                    unique(which.pred.its2)),
            lat = unique(lat),
            long = unique(long)) %>%
  # mutate(which.pred = if_else(str_detect(method, 'os'), T, which.pred))
  mutate(which.pred.sw = if_else(str_detect(method, 'os'), 1, which.pred.sw),
         which.pred.n = if_else(str_detect(method, 'os'), 1, which.pred.n))


# contours
theta_grid <- seq(0, 2*pi, length = 100)
week_grid <- seq(week(ymd(paste(yrs[j], '/03/01', sep = ''))), 
                 week(ymd(paste(yrs[j], '/08/01', sep = ''))),
                 by = 3)
date_grid <- rep(ymd(paste(yrs[j], '/04/15', sep = '')), length(week_grid))
week(date_grid) <- week_grid
pred_df <- expand.grid(theta = theta_grid, 
                       date = date_grid) %>%
  mutate(time = yday(date),
         sinb1.src1 = sin(theta),
         sinb2.src1 = sin(theta + pi/4),
         sinb1 = sin(theta),
         sinb2 = sin(theta + pi/4))



pred_df <- pred_df %>%
  mutate(pred.ats1.src1 = predict(fit_ats1_src1, newdata = pred_df),
         pred.ats1.src2 = predict(fit_ats1_src2, newdata = pred_df),
         pred.its1.src1 = predict(fit_its1_src1, newdata = pred_df),
         pred.its1.src2 = predict(fit_its1_src2, newdata = pred_df),
         pred.ats2.src1 = predict(fit_ats2_src1, newdata = pred_df),
         pred.ats2.src3 = predict(fit_ats2_src2, newdata = pred_df),
         pred.its2.src1 = predict(fit_its2_src1, newdata = pred_df),
         pred.its2.src3 = predict(fit_its2_src2, newdata = pred_df),
         pred.aos = predict(fit_aos, newdata = pred_df),
         pred.ios = predict(fit_ios, newdata = pred_df),
         g.theta.src1 = density.circular(x = filter(wind_src1, year == yrs[j])$direction,
                                         z = pred_df$theta, 
                                         bw = 50)$y,
         g.theta.src2 = density.circular(x = filter(wind_src2, year == yrs[j])$direction,
                                         z = pred_df$theta, 
                                         bw = 50)$y,
         g.theta.src3 = density.circular(x = filter(wind_src3, year == yrs[j])$direction,
                                         z = pred_df$theta, 
                                         bw = 50)$y) %>%
  mutate(rho.ats1.src1 = (exp(pred.ats1.src1) - 1)*g.theta.src1,
         rho.ats1.src2 = (exp(pred.ats1.src2) - 1)*g.theta.src2,
         rho.its1.src1 = (exp(pred.its1.src1) - 1),
         rho.its1.src2 = (exp(pred.its1.src2) - 1),
         rho.ats2.src1 = (exp(pred.ats2.src1) - 1)*g.theta.src1,
         rho.ats2.src3 = (exp(pred.ats2.src3) - 1)*g.theta.src3,
         rho.its2.src1 = (exp(pred.its2.src1) - 1),
         rho.its2.src3 = (exp(pred.its2.src3) - 1),
         rho.aos = (exp(pred.aos) - 1)*g.theta.src1,
         rho.ios = exp(pred.ios) - 1) %>%
  select(theta, date, time, 
         rho.ats1.src1, rho.ats1.src2, 
         rho.its1.src1, rho.its1.src2, 
         rho.ats2.src1, rho.ats2.src3, 
         rho.its2.src1, rho.its2.src3, 
         rho.aos, 
         rho.ios)

pred_df <- pred_df %>%
  cbind(ats1.src1 = destPoint(p = source_df[1, 2:3],
                             b = pred_df$theta*360/(2*pi),
                             d = 1000*pred_df$rho.ats1.src1),
        ats1.src2 = destPoint(p = source_df[2, 2:3],
                             b = pred_df$theta*360/(2*pi),
                             d = 1000*pred_df$rho.ats1.src2),
        its1.src1 = destPoint(p = source_df[1, 2:3],
                             b = pred_df$theta*360/(2*pi),
                             d = 1000*pred_df$rho.its1.src1),
        its1.src2 = destPoint(p = source_df[2, 2:3],
                             b = pred_df$theta*360/(2*pi),
                             d = 1000*pred_df$rho.its1.src2),
        ats2.src1 = destPoint(p = source_df[1, 2:3],
                              b = pred_df$theta*360/(2*pi),
                              d = 1000*pred_df$rho.ats2.src1),
        ats2.src3 = destPoint(p = source_df[3, 2:3],
                              b = pred_df$theta*360/(2*pi),
                              d = 1000*pred_df$rho.ats2.src3),
        its2.src1 = destPoint(p = source_df[1, 2:3],
                              b = pred_df$theta*360/(2*pi),
                              d = 1000*pred_df$rho.its2.src1),
        its2.src3 = destPoint(p = source_df[3, 2:3],
                              b = pred_df$theta*360/(2*pi),
                              d = 1000*pred_df$rho.its2.src3),
        aos = destPoint(p = source_df[1, 2:3],
                        b = pred_df$theta*360/(2*pi),
                        d = 1000*pred_df$rho.aos),
        ios = destPoint(p = source_df[1, 2:3],
                        b = pred_df$theta*360/(2*pi),
                        d = 1000*pred_df$rho.ios))

return(list(pred = pred_df, 
            plot = plot_df))
}
