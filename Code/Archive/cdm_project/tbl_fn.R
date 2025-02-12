tbl_fn <- function(sub_data, 
                    wind_src1, wind_src2, 
                    source_df, src1_id, src2_id,
                   src1_type, src2_type){
  
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
           z.fit = z.fit)
  
  fit_its1_out <- fit_ts(source1_its1, 
                         source2_its1, 
                         initdata_its1, 
                         isotropy = T)
  fit_its1_src1 <- fit_its1_out$fit1
  fit_its1_src2 <- fit_its1_out$fit2
  pr_z_its1 <- fit_its1_out$prz
  

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
           z.fit = z.fit)
  
  fit_ats1_out <- fit_ts(source1_ats1, 
                         source2_ats1, 
                         initdata_ats1, 
                         isotropy = F)
  fit_ats1_src1 <- fit_ats1_out$fit1
  fit_ats1_src2 <- fit_ats1_out$fit2
  pr_z_ats1 <- fit_ats1_out$prz
  

  # extract estimates from each model 
  betahat_ios <- coef(fit_ios)
  betahat_aos <- coef(fit_aos)
  betahat_ats11 <- coef(fit_ats1_src1)
  betahat_ats12 <- coef(fit_ats1_src2)
  betahat_its11 <- coef(fit_its1_src1)
  betahat_its12 <- coef(fit_its1_src2)

  # predictions of time of occurrence according to each model
  sub_data <- ungroup(sub_data) %>%
    mutate(t.hat.ats1.src1 = (aresponse.src1 - betahat_ats11[1] - betahat_ats11[3]*sinb1.src1 - betahat_ats11[4]*sinb2.src1)/betahat_ats11[2],
           t.hat.ats1.src2 = (aresponse.src2 - betahat_ats12[1] - betahat_ats12[3]*sinb1.src2 - betahat_ats12[4]*sinb2.src2)/betahat_ats12[2],
           t.hat.aos = (aresponse.src1 - betahat_aos[1] - betahat_aos[3]*sinb1.src1 - betahat_aos[4]*sinb2.src1)/betahat_aos[2],
           t.hat.ios = (iresponse.src1 - betahat_ios[1])/betahat_ios[2],
           t.hat.its1.src1 = (iresponse.src1 - betahat_its11[1])/betahat_its11[2],
           t.hat.its1.src2 = (iresponse.src2 - betahat_its12[1])/betahat_its12[2],
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
                                                   sinb2 = sinb2.src2))) - 1)*g.theta.src2)
           
  # round p_i's and choose predictions
  which_pred_ats1 <- data.frame(which.pred.ats1 = round(pr_z_ats1, 0) + (round(pr_z_ats1, 0) == 0)*2,
                                ID = source1_ats1$ID)
  which_pred_its1 <- data.frame(which.pred.its1 = round(pr_z_its1, 0) + (round(pr_z_its1, 0) == 0)*2,
                                ID = source1_its1$ID)
  
  sub_data_pred <- sub_data %>%
    as_tibble() %>%
    merge(which_pred_ats1, by = "ID", all = T) %>%
    merge(which_pred_its1, by = "ID", all = T) %>%
    mutate(which.pred.ats1 = replace_na(which.pred.ats1, 1),
           which.pred.its1 = replace_na(which.pred.its1, 1)) %>%
    mutate(t.resid.aos = t.hat.aos - time,
           t.resid.ios = t.hat.ios - time,
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
           rho.resid.ios = rho.src1.hat.ios - rho.src1,
           rho.resid.aos = rho.src1.hat.aos - rho.src1) %>%
    select(symptom.date, 
           lat, 
           long, 
           ID,
           t.resid.ats1,
           t.resid.its1,
           t.resid.aos,
           t.resid.ios,
           rho.resid.ats1,
           rho.resid.its1,
           rho.resid.ios,
           rho.resid.aos,
           which.pred.ats1,
           which.pred.its1,
           rho.src1,
           theta.src1,
           rho.src2,
           theta.src2,
           host)
  
  # replace source predictions by NA
  names(sub_data_pred[, 5:12])
  names(sub_data_pred[, c(5:6, 9:10)])
  sub_data_pred[sub_data_pred$ID == src1_id, 5:12] <- NA
  sub_data_pred[sub_data_pred$ID == src2_id, c(5:6, 9:10)] <- NA
  
  ############
  
  out_tbl <- sub_data_pred %>%
    summarize(ats1.t.rmse = sqrt(mean(t.resid.ats1^2, na.rm = T)),
              its1.t.rmse = sqrt(mean(t.resid.its1^2, na.rm = T)),
              aos.t.rmse = sqrt(mean(t.resid.aos^2, na.rm = T)),
              ios.t.rmse = sqrt(mean(t.resid.ios^2, na.rm = T)),
              ats.rho.rmse = sqrt(mean(rho.resid.ats1^2, na.rm = T)),
              its.rho.rmse = sqrt(mean(rho.resid.its1^2, na.rm = T)),
              aos.rho.rmse = sqrt(mean(rho.resid.aos^2, na.rm = T)),
              ios.rho.rmse = sqrt(mean(rho.resid.ios^2, na.rm = T)),
              n = n()) %>%
    mutate(year = unique(source_df$year),
           src2.state = source_df[2, ]$state,
           src2.county = source_df[2, ]$county,
           src2.date = source_df[2, ]$symptom.date,
           src2.type = src2_type,
           src1.state = source_df[1, ]$state,
           src1.county = source_df[1, ]$county,
           src1.date = source_df[1, ]$symptom.date,
           src1.type = src1_type)
  
  return(out_tbl)
  
}
