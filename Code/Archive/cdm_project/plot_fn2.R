plot_fn2 <- function(sub_data, 
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
    mutate(t.hat.ats1 = if_else(which.pred.ats1 == 1,
                                  t.hat.ats1.src1,
                                  t.hat.ats1.src2),
           t.hat.its1 = if_else(which.pred.its1 == 1,
                                  t.hat.its1.src1,
                                  t.hat.its1.src2),
           rho.hat.its1 = if_else(which.pred.its1 == 1,
                                    rho.src1.hat.its1,
                                    rho.src2.hat.its1),
           rho.its1 = if_else(which.pred.its1 == 1,
                                  rho.src1,
                                  rho.src2),
           rho.hat.ats1 = if_else(which.pred.ats1 == 1,
                                    rho.src1.hat.ats1,
                                    rho.src2.hat.ats1),
           rho.ats1 = if_else(which.pred.ats1 == 1,
                                  rho.src1,
                                  rho.src2),
           t.hat.ats2 = if_else(which.pred.ats2 == 1,
                                  t.hat.ats2.src1,
                                  t.hat.ats2.src3),
           t.hat.its2 = if_else(which.pred.its2 == 1,
                                  t.hat.its2.src1,
                                  t.hat.its2.src3),
           rho.hat.its2 = if_else(which.pred.its2 == 1,
                                    rho.src1.hat.its2,
                                    rho.src3.hat.its2),
           rho.its2 = if_else(which.pred.its2 == 1,
                                  rho.src1,
                                  rho.src3),
           rho.hat.ats2 = if_else(which.pred.ats2 == 1,
                                    rho.src1.hat.ats2,
                                    rho.src3.hat.ats2),
           rho.ats2 = if_else(which.pred.ats2 == 1,
                                  rho.src1,
                                  rho.src3),
           rho.hat.ios = rho.src1.hat.ios,
           rho.ios = rho.src1,
           rho.hat.aos = rho.src1.hat.aos,
           rho.aos = rho.src1) %>%
    select(symptom.date,
           time,
           lat, 
           long, 
           ID,
           t.hat.ats1,
           t.hat.its1,
           t.hat.ats2,
           t.hat.its2,
           t.hat.aos,
           t.hat.ios,
           rho.hat.ats1,
           rho.ats1,
           rho.hat.its1,
           rho.its1,
           rho.hat.ats2,
           rho.ats2,
           rho.hat.its2,
           rho.its2,
           rho.hat.ios,
           rho.ios,
           rho.hat.aos,
           rho.aos,
           which.pred.ats1,
           which.pred.its1,
           which.pred.ats2,
           which.pred.its2,
           host)
  
  # replace source predictions by NA
  names(sub_data_pred[, 5:16])
  names(sub_data_pred[, c(5:6, 11:12)])
  names(sub_data_pred[, c(7:8, 13:14)])
  sub_data_pred[sub_data_pred$ID == src1_id, c(6:11, 12, 14, 16, 18, 20, 22)] <- NA
  sub_data_pred[sub_data_pred$ID == src2_id, c(6:7, 12, 14)] <- NA
  sub_data_pred[sub_data_pred$ID == src3_id, c(8:9, 16, 18)] <- NA
  
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
  out1 <- sub_data_pred %>%
    select(c(1:11, 24:28)) %>%
    rename(AOS = t.hat.aos,
           `ATS-SW` = t.hat.ats1,
           `ATS-N` = t.hat.ats2,
           IOS = t.hat.ios,
           `ITS-SW` = t.hat.its1,
           `ITS-N` = t.hat.its2,
           which.ats.n = which.pred.ats2,
           which.its.n = which.pred.its2,
           which.ats.sw = which.pred.ats1,
           which.its.sw = which.pred.its1) %>%
    gather('method', 'pred', 6:11)
  
  out2 <- sub_data_pred %>%
    select(c(1:5, 12:28)) %>% 
    rename(AOS = rho.hat.aos,
           `ATS-SW` = rho.hat.ats1,
           `ATS-N` = rho.hat.ats2,
           IOS = rho.hat.ios,
           `ITS-SW` = rho.hat.its1,
           `ITS-N` = rho.hat.its2,
           rho.ats.n = rho.ats2,
           rho.ats.sw = rho.ats1,
           rho.its.n = rho.its2,
           rho.its.sw = rho.its1,
           which.ats.n = which.pred.ats2,
           which.its.n = which.pred.its2,
           which.ats.sw = which.pred.ats1,
           which.its.sw = which.pred.its1) %>%
    gather('method', 'pred', c(6, 8, 10, 14, 16))
  
  which0 <- out1 %>%
    filter(method == 'IOS') %>%
    select(ID, method) %>%
    mutate(which.pred = 1)
  
  which1 <- out1 %>%
    filter(method == 'AOS') %>%
    select(ID, method) %>%
    mutate(which.pred = 1)
  
  which2 <- out1 %>%
    filter(method == 'ATS-N') %>%
    select(ID, method, which.ats.n) %>%
    rename(which.pred = which.ats.n)
    
  which3 <- out1 %>%
    filter(method == 'ATS-SW') %>%
    select(ID, method, which.ats.sw) %>%
    rename(which.pred = which.ats.sw)
  
  which4 <- out1 %>%
    filter(method == 'ITS-N') %>%
    select(ID, method, which.its.n) %>%
    rename(which.pred = which.its.n)
  
  which5 <- out1 %>%
    filter(method == 'ITS-SW') %>%
    select(ID, method, which.its.sw) %>%
    rename(which.pred = which.its.sw)
  
  which6 <- out2 %>%
    filter(method == 'IOS') %>%
    select(ID, method, rho.ios) %>%
    mutate(which.pred = 1,
           rho = rho.ios)
  
  which7 <- out2 %>%
    filter(method == 'AOS') %>%
    select(ID, method, rho.aos) %>%
    mutate(which.pred = 1,
           rho = rho.aos)
  
  which8 <- out2 %>%
    filter(method == 'ATS-N') %>%
    select(ID, method, which.ats.n, rho.ats.n) %>%
    rename(which.pred = which.ats.n,
           rho = rho.ats.n)
  
  which9 <- out2 %>%
    filter(method == 'ATS-SW') %>%
    select(ID, method, which.ats.sw, rho.ats.sw) %>%
    rename(which.pred = which.ats.sw,
           rho = rho.ats.sw)
  
  which10 <- out2 %>%
    filter(method == 'ITS-N') %>%
    select(ID, method, which.its.n, rho.its.n) %>%
    rename(which.pred = which.its.n,
           rho = rho.its.n)
  
  which11 <- out2 %>%
    filter(method == 'ITS-SW') %>%
    select(ID, method, which.its.sw, rho.its.sw) %>%
    rename(which.pred = which.its.sw,
           rho = rho.its.sw)
  
  out1 <- bind_rows(which0, which1, which2, which3, which4, which5) %>%
    merge(out1, by = c('ID', 'method')) %>%
    select(ID, method, which.pred, time, pred, host)
  
  out2 <- bind_rows(which6, which7, which8, which9, which10, which11) %>%
    merge(out2, by = c('ID', 'method')) %>%
    select(ID, method, which.pred, rho, pred, host)
  
  
  return(list(time = out1,
              dist = out2))
}
