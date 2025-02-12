estimate_fn <- function(sub_data, 
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
  
  # add CIs
  ci_ios <- confint(fit_ios, level = 0.95)
  betahat_ios <- c(betahat_ios[1], ci_ios[1, ], 
                   betahat_ios[2], ci_ios[2, ],
                   rep(NA, 6))
  ci_its11 <- confint(fit_its1_src1, level = 0.95)
  betahat_its11 <- c(betahat_its11[1], ci_its11[1, ], 
                   betahat_its11[2], ci_its11[2, ],
                   rep(NA, 6))
  ci_its12 <- confint(fit_its1_src2, level = 0.95)
  betahat_its12 <- c(betahat_its12[1], ci_its12[1, ], 
                     betahat_its12[2], ci_its12[2, ],
                     rep(NA, 6))
  ci_its21 <- confint(fit_its2_src1, level = 0.95)
  betahat_its21 <- c(betahat_its21[1], ci_its21[1, ], 
                     betahat_its21[2], ci_its21[2, ],
                     rep(NA, 6))
  ci_its22 <- confint(fit_its2_src2, level = 0.95)
  betahat_its22 <- c(betahat_its22[1], ci_its22[1, ], 
                     betahat_its22[2], ci_its22[2, ],
                     rep(NA, 6))
  ci_aos <- confint(fit_aos, level = 0.95)
  betahat_aos <- c(betahat_aos[1], ci_aos[1, ], 
                   betahat_aos[2], ci_aos[2, ],
                   betahat_aos[3], ci_aos[3, ],
                   betahat_aos[4], ci_aos[4, ])
  ci_ats11 <- confint(fit_ats1_src1, level = 0.95)
  betahat_ats11 <- c(betahat_ats11[1], ci_ats11[1, ], 
                   betahat_ats11[2], ci_ats11[2, ],
                   betahat_ats11[3], ci_ats11[3, ],
                   betahat_ats11[4], ci_ats11[4, ])
  ci_ats12 <- confint(fit_ats1_src2, level = 0.95)
  betahat_ats12 <- c(betahat_ats12[1], ci_ats12[1, ], 
                     betahat_ats12[2], ci_ats12[2, ],
                     betahat_ats12[3], ci_ats12[3, ],
                     betahat_ats12[4], ci_ats12[4, ])
  ci_ats21 <- confint(fit_ats2_src1, level = 0.95)
  betahat_ats21 <- c(betahat_ats21[1], ci_ats21[1, ], 
                     betahat_ats21[2], ci_ats21[2, ],
                     betahat_ats21[3], ci_ats21[3, ],
                     betahat_ats21[4], ci_ats21[4, ])
  ci_ats22 <- confint(fit_ats2_src2, level = 0.95)
  betahat_ats22 <- c(betahat_ats22[1], ci_ats22[1, ], 
                     betahat_ats22[2], ci_ats22[2, ],
                     betahat_ats22[3], ci_ats22[3, ],
                     betahat_ats22[4], ci_ats22[4, ])
  
  names(betahat_ios) <- names(betahat_its11) <- names(betahat_its12) <- names(betahat_its21) <- names(betahat_its22) <- names(betahat_ats11)
  
  out <- rbind(ios = betahat_ios, 
               its.sw.s = betahat_its11,
               its.sw.sw = betahat_its12,
               its.n.s = betahat_its21,
               its.n.n = betahat_its22,
               aos = betahat_aos,
               ats.sw.s = betahat_ats11,
               ats.sw.sw = betahat_ats12,
               ats.n.s = betahat_ats21,
               ats.n.n = betahat_ats22)
  
  return(out)
  }