source('two_source_preprocessing_fn.R')
source('fit_ats_fn.R')
seed <- 41420

yrs <- 2008:2011
wks <- c(25, 26, 29, 29)
locs <- c(0, 1, 0, 1)


set.seed(seed)
yr_out <- lapply(1:4, function(j){
  # preprocess data for selected n source location and time
  DATA <- preprocess_fn(wks[j], locs[j])
  
  # subset data to year and calculate basis expansions and response variables
  sub_data <- filter(DATA, year(report_date) == yrs[j]) %>%
    mutate(month = month(report_date),
           day = day(report_date),
           time = yday(report_date),
           response.0 = log(1 + rho1), # response for isotropic model
           response.1 = log(1 + rho1/g.theta.1), # response for S source
           response.2 = log(1 + rho2/g.theta.2), # response for N source
           sinb1.1 = sin(as.numeric(theta1)), # sin basis 1 for S source
           sinb2.1 = sin(as.numeric(theta1) + pi/4), # sin basis 2 for S source
           sinb1.2 = sin(as.numeric(theta2)), # sin basis 1 for N source
           sinb2.2 = sin(as.numeric(theta2) + pi/4), # sin basis 2 for N source
           z.reg.1 = rho1/g.theta.1, # regressor 1 to initialize weights
           z.reg.2 = rho2/g.theta.2, # regressor 2 to initialize weights
           z.fit = rho1 < rho2) # is southern source closer?
  
  
  ## ----splitting-----------------------------------------------------------
  # split data into training and test sets
  train_data <- sub_data
  
  ## ----fitting-------------------------------------------------------------
  # fit IOS model
  fit_ios <- lm(response.0 ~ time, 
                data = train_data)
  
  # fit AOS model
  fit_aos <- lm(response.1 ~ time + sinb1.1 + sinb2.1, 
                data = train_data)
  
  # fit ATS model
  fit_ats_out <- fit_ats(train_data)
  fit_ats1 <- fit_ats_out$ats1
  fit_ats2 <- fit_ats_out$ats2
  pr_z <- fit_ats_out$prz
  
  # compute combined adjusted R2 for each model
  adj.r.ios <- summary(fit_ios)$adj.r
  adj.r.aos <- summary(fit_aos)$adj.r
  resids <- train_data %>%
    mutate(pred1 = predict(fit_ats1, train_data),
           pred2 = predict(fit_ats2, train_data)) %>%
    transmute(resid1 = response.1 - pred1,
              resid2 = response.2 - pred2,
              resp1 = response.1,
              resp2 = response.2) %>%
    transmute(resid = if_else(round(pr_z, 0) == 1, resid1, resid2),
              resp = if_else(round(pr_z, 0) == 1, resp1, resp2))
  sse <- sum(resids$resid^2)
  sst <- sum((resids$resp - mean(resids$resp))^2)
  adj.r.ats <- 1 - (sse/sst)*((n_train - 1)/(n_train - 8 - 1))
  adj.rs <- c(ios = adj.r.ios,
              aos = adj.r.aos,
              ats = adj.r.ats,
              year = yrs[j],
              rep.num = rep)
  
  # extract estimates from each model 
  betahat_ios <- coef(fit_ios)
  betahat_aos <- coef(fit_aos)
  betahat_ats1 <- coef(fit_ats1)
  betahat_ats2 <- coef(fit_ats2)
  estimates <- c(ios = coef(fit_ios), 
                 aos = coef(fit_aos),
                 ats1 = coef(fit_ats1),
                 ats2 = coef(fit_ats2),
                 year = yrs[j],
                 rep.num = rep)
  
  estimate_tbl <- rbind(data.frame(summary(fit_ios)$coef,
                                    parameter = paste('ios.', 
                                                      names(betahat_ios), 
                                                      sep = ''),
                                    year = yrs[j]),
                        data.frame(summary(fit_aos)$coef,
                                   parameter = paste('aos.', 
                                                     names(betahat_aos), 
                                                     sep = ''),
                                   year = yrs[j]),
                        data.frame(summary(fit_ats1)$coef,
                                   parameter = paste('ats1.', 
                                                     names(betahat_ats1), 
                                                     sep = ''),
                                   year = yrs[j]),
                        data.frame(summary(fit_ats2)$coef,
                                   parameter = paste('ats2.', 
                                                     names(betahat_ats2), 
                                                     sep = ''),
                                   year = yrs[j]))
  
  ## ----predicitons---------------------------------------------------------
  # predictions of time of occurrence according to each model
  sub_data <- mutate(sub_data, 
                     t.hat.ats1 = (response.1 - betahat_ats1[1] - betahat_ats1[3]*sinb1.1 - betahat_ats1[4]*sinb2.1)/betahat_ats1[2],
                     t.hat.ats2 = (response.2 - betahat_ats2[1] - betahat_ats2[3]*sinb1.2 - betahat_ats2[4]*sinb2.2)/betahat_ats2[2],
                     t.hat.aos = (response.1 - betahat_aos[1] - betahat_aos[3]*sinb1.1 - betahat_aos[4]*sinb2.1)/betahat_aos[2],
                     t.hat.ios = (response.0 - betahat_ios[1])/betahat_ios[2]) 
  
  # refit ats model to estimate p_i's for test data
  pr_z_full <- pr_z
  
  # round p_i's and choose time predictions
  which_pred_ats <- round(pr_z_full, 0) + (round(pr_z_full, 0) == 0)*2
  pred_df <- cbind(sub_data,
                   which.pred.ats = which_pred_ats) %>%
    mutate(t.hat.ats = if_else(which.pred.ats == 1, 
                               t.hat.ats1, 
                               t.hat.ats2)) %>%
    select(time, 
           t.hat.ios, 
           t.hat.aos, 
           t.hat.ats, 
           which.pred.ats) %>%
    rename(ios = t.hat.ios,
           aos = t.hat.aos,
           ats = t.hat.ats) %>%
    mutate(year = yrs[j])
  
  out <- list(preds = pred_df,
              ests = estimates,
              est_tbl = estimate_tbl,
              rsqs = adj.rs)
  
  return(out)
})

out <- list(rsqs = Reduce(rbind, lapply(1:4, function(yr){yr_out[[yr]]$rsqs})),
            ests = Reduce(rbind, lapply(1:4, function(yr){yr_out[[yr]]$ests})),
            preds = Reduce(rbind, lapply(1:4, function(yr){yr_out[[yr]]$preds})),
            est_tbl = Reduce(rbind, lapply(1:4, function(yr){yr_out[[yr]]$est_tbl})))

out$est_tbl[, c(6, 5, 1:4)] %>%
  xtable(digits = 3)
out$rsqs[, c(4, 1:3)] %>%
  xtable(digits = 3)
tbl <- out$preds %>%
  as_tibble() %>%
  mutate(sq.resid.ios = (time - ios)^2,
         sq.resid.aos = (time - aos)^2,
         sq.resid.ats = (time - ats)^2,
         resid.ios = (time - ios),
         resid.aos = (time - aos),
         resid.ats = (time - ats)) %>%
  group_by(year) %>%
  summarize(rmse.ios = sqrt(mean(sq.resid.ios)),
            rmse.aos = sqrt(mean(sq.resid.aos)),
            rmse.ats = sqrt(mean(sq.resid.ats)),
            mres.ios = mean(resid.ios),
            mres.aos = mean(resid.aos),
            mres.ats = mean(resid.ats)) %>%
  as.data.frame()

tbl <- tbl[, c(1, 2, 5, 3, 6, 4, 7)]
xtable(tbl, digits = 2)
