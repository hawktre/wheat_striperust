source('two_source_preprocessing_fn.R')
source('fit_ats_fn.R')
source('fit_its_fn.R')
seed <- 42820

yrs <- 2008:2011
wks <- c(25, 26, 29, 29)
locs <- c(0, 1, 0, 1)
set.seed(seed)

j <- 4
nReps <- 100
out2011 <- Reduce(rbind, 
                  lapply(1:nReps, 
                         function(split){
                           # preprocess data for selected n source location and time
                           DATA <- preprocess_fn(wks[j], locs[j])
                           
                           # subset data to year and calculate basis expansions and response variables
                           sub_data <- filter(DATA, year(report_date) == yrs[j]) %>%
                             mutate(month = month(report_date),
                                    day = day(report_date),
                                    time = yday(report_date),
                                    iresponse.1 = log(1 + rho1), # isotropic response for S source
                                    iresponse.2 = log(1 + rho2), # isotropic response for N source
                                    aresponse.1 = log(1 + rho1/g.theta.1), # anisotropic response for S source
                                    aresponse.2 = log(1 + rho2/g.theta.2), # anisotropic response for N source
                                    sinb1.1 = sin(as.numeric(theta1)), # sin basis 1 for S source
                                    sinb2.1 = sin(as.numeric(theta1) + pi/4), # sin basis 2 for S source
                                    sinb1.2 = sin(as.numeric(theta2)), # sin basis 1 for N source
                                    sinb2.2 = sin(as.numeric(theta2) + pi/4), # sin basis 2 for N source
                                    z.reg.1 = rho1/g.theta.1, # regressor 1 to initialize weights
                                    z.reg.2 = rho2/g.theta.2, # regressor 2 to initialize weights
                                    z.fit = rho1 < rho2) # is southern source closer?
                           
                           # split data into training and test sets
                           training_prop <- 0.8
                           n <- nrow(sub_data)
                           n_train <- floor(training_prop*n)
                           train_idx <- sample(1:n, n_train)
                           
                           train_data <- sub_data[train_idx, ]
                           test_data <- sub_data[-train_idx, ]
                           
                           # fit IOS model
                           fit_ios <- lm(iresponse.1 ~ time, 
                                         data = train_data)
                           
                           # fit AOS model
                           fit_aos <- lm(aresponse.1 ~ time + sinb1.1 + sinb2.1, 
                                         data = train_data)
                           
                           # fit ITS model
                           fit_its_out <- fit_its(train_data)
                           fit_its1 <- fit_its_out$its1
                           fit_its2 <- fit_its_out$its2
                           pr_z_its <- fit_its_out$prz
                           
                           # fit ATS model
                           fit_ats_out <- fit_ats(train_data)
                           fit_ats1 <- fit_ats_out$ats1
                           fit_ats2 <- fit_ats_out$ats2
                           pr_z_ats <- fit_ats_out$prz
                           
                           # compute combined adjusted R2 for each model
                           adj.r.ios <- summary(fit_ios)$adj.r
                           adj.r.aos <- summary(fit_aos)$adj.r
                           resids_ats <- train_data %>%
                             mutate(pred1 = predict(fit_ats1, train_data),
                                    pred2 = predict(fit_ats2, train_data)) %>%
                             transmute(resid1 = aresponse.1 - pred1,
                                       resid2 = aresponse.2 - pred2,
                                       resp1 = aresponse.1,
                                       resp2 = aresponse.2) %>%
                             transmute(resid = if_else(round(pr_z_ats, 0) == 1, resid1, resid2),
                                       resp = if_else(round(pr_z_ats, 0) == 1, resp1, resp2))
                           sse_ats <- sum(resids_ats$resid^2)
                           sst_ats <- sum((resids_ats$resp - mean(resids_ats$resp))^2)
                           adj.r.ats <- 1 - (sse_ats/sst_ats)*((n_train - 1)/(n_train - 8 - 1))
                           resids_its <- train_data %>%
                             mutate(pred1 = predict(fit_its1, train_data),
                                    pred2 = predict(fit_its2, train_data)) %>%
                             transmute(resid1 = iresponse.1 - pred1,
                                       resid2 = iresponse.2 - pred2,
                                       resp1 = iresponse.1,
                                       resp2 = iresponse.2) %>%
                             transmute(resid = if_else(round(pr_z_its, 0) == 1, resid1, resid2),
                                       resp = if_else(round(pr_z_its, 0) == 1, resp1, resp2))
                           sse_its <- sum(resids_its$resid^2)
                           sst_its <- sum((resids_its$resp - mean(resids_its$resp))^2)
                           adj.r.its <- 1 - (sse_its/sst_its)*((n_train - 1)/(n_train - 8 - 1))
                           
                           adj.rs <- c(ios = adj.r.ios,
                                       aos = adj.r.aos,
                                       ats = adj.r.ats,
                                       its = adj.r.its,
                                       year = yrs[j])
                           
                           # extract estimates from each model 
                           betahat_ios <- coef(fit_ios)
                           betahat_aos <- coef(fit_aos)
                           betahat_ats1 <- coef(fit_ats1)
                           betahat_ats2 <- coef(fit_ats2)
                           betahat_its1 <- coef(fit_its1)
                           betahat_its2 <- coef(fit_its2)
                           
                           estimates <- c(ios = coef(fit_ios), 
                                          aos = coef(fit_aos),
                                          ats1 = coef(fit_ats1),
                                          ats2 = coef(fit_ats2),
                                          its1 = coef(fit_its1),
                                          its2 = coef(fit_its2),
                                          year = yrs[j])
                           
                           # predictions of time of occurrence according to each model
                           sub_data <- mutate(sub_data, 
                                              t.hat.ats1 = (aresponse.1 - betahat_ats1[1] - betahat_ats1[3]*sinb1.1 - betahat_ats1[4]*sinb2.1)/betahat_ats1[2],
                                              t.hat.ats2 = (aresponse.2 - betahat_ats2[1] - betahat_ats2[3]*sinb1.2 - betahat_ats2[4]*sinb2.2)/betahat_ats2[2],
                                              t.hat.aos = (aresponse.1 - betahat_aos[1] - betahat_aos[3]*sinb1.1 - betahat_aos[4]*sinb2.1)/betahat_aos[2],
                                              t.hat.ios = (iresponse.1 - betahat_ios[1])/betahat_ios[2],
                                              t.hat.its1 = (iresponse.1 - betahat_its1[1])/betahat_its1[2],
                                              t.hat.its2 = (iresponse.2 - betahat_its2[1])/betahat_its2[2],
                                              rho1.hat.ios = (exp(predict(fit_ios, sub_data)) - 1),
                                              rho1.hat.aos = (exp(predict(fit_aos, sub_data)) - 1)*g.theta.1,
                                              rho1.hat.its = (exp(predict(fit_its1, sub_data)) - 1),
                                              rho2.hat.its = (exp(predict(fit_its2, sub_data)) - 1),
                                              rho1.hat.ats = (exp(predict(fit_ats1, sub_data)) - 1)*g.theta.1,
                                              rho2.hat.ats = (exp(predict(fit_ats2, sub_data)) - 1)*g.theta.2)
                           
                           # refit ats and its models to estimate p_i's for test data
                           fit_ats_full <- fit_ats(sub_data)
                           pr_z_ats_full <- fit_ats_full$prz
                           fit_its_full <- fit_its(sub_data)
                           pr_z_its_full <- fit_its_full$prz
                           
                           # round p_i's and choose predictions
                           which_pred_ats <- round(pr_z_ats_full, 0) + (round(pr_z_ats_full, 0) == 0)*2
                           which_pred_its <- round(pr_z_its_full, 0) + (round(pr_z_its_full, 0) == 0)*2
                           
                           
                           pred_df <- cbind(sub_data,
                                            which.pred.ats = which_pred_ats,
                                            which.pred.its = which_pred_its) %>%
                             mutate(t.resid.ats = if_else(which.pred.ats == 1, 
                                                          t.hat.ats1, 
                                                          t.hat.ats2) - time,
                                    t.resid.its = if_else(which.pred.its == 1, 
                                                          t.hat.its1, 
                                                          t.hat.its2) - time,
                                    t.resid.aos = t.hat.aos - time,
                                    t.resid.ios = t.hat.ios - time,
                                    rho.resid.its = if_else(which.pred.its == 1,
                                                            rho1.hat.its - rho1,
                                                            rho2.hat.its - rho2),
                                    rho.resid.ats = if_else(which.pred.ats == 1,
                                                            rho1.hat.ats - rho1,
                                                            rho2.hat.ats - rho2),
                                    rho.resid.ios = rho1.hat.ios - rho1,
                                    rho.resid.aos = rho1.hat.aos - rho1) %>%
                             select(report_date, 
                                    lat, 
                                    long, 
                                    rho1,
                                    theta1,
                                    rho2,
                                    theta2,
                                    host,
                                    t.resid.ats,
                                    t.resid.its,
                                    t.resid.aos,
                                    t.resid.ios,
                                    rho.resid.ats,
                                    rho.resid.its,
                                    rho.resid.ios,
                                    rho.resid.aos,
                                    which.pred.ats,
                                    which.pred.its)
                           
                           pred_df$subset <- 0
                           pred_df$subset[-train_idx] <- 'test'
                           pred_df$subset[train_idx] <- 'train'
                           
                           out <- pred_df %>%
                             filter(subset == 'test') 
                           
                           print(split)
                           return(out)
                         }
                  )
)

# save(out2008, file = 'pred_assessment_2008.RData')
# save(out2009, file = 'pred_assessment_2009.RData')
# save(out2010, file = 'pred_assessment_2010.RData')
save(out2011, file = 'pred_assessment_2011.RData')

# load('pred_assessment_2008.RData')
# load('pred_assessment_2009.RData')
# load('pred_assessment_2010.RData')
# load('pred_assessment_2011.RData')

# bin by month and octants centered on cardinal directions and examine errors (rmse within bin)

octants <- seq(0, 2*pi, length = 9) - 2*pi/16
# out2008 %>%
out2009 %>%
# out2010 %>%
# out2011 %>%
  mutate(direction = cut(theta1 - 2*pi/16, 
                         breaks = octants,
                         labels = c('N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW'))) %>%
  group_by(direction, month(report_date)) %>%
  summarize(ats.t.rmse = sqrt(mean(t.resid.ats^2)),
          its.t.rmse = sqrt(mean(t.resid.its^2)),
          aos.t.rmse = sqrt(mean(t.resid.aos^2)),
          ios.t.rmse = sqrt(mean(t.resid.ios^2)),
          ats.rho.rmse = sqrt(mean(rho.resid.ats^2)),
          its.rho.rmse = sqrt(mean(rho.resid.its^2)),
          aos.rho.rmse = sqrt(mean(rho.resid.aos^2)),
          ios.rho.rmse = sqrt(mean(rho.resid.ios^2)),
          n = n()) %>%
  # select(1:2, 11)
  # select(c(1:2, 3:6)) %>%
  select(c(1:2, 7:10)) %>%
  gather('method', 'rmse', 3:6) %>%
  mutate(sources = factor(str_detect(method, 'ts'), 
                          labels = c('one', 'two')),
         isotropy = factor(str_starts(method, 'i'), 
                           labels = c('anisotropic', 'isotropic'))) %>%
  ggplot(aes(x = `month(report_date)`,
             y = rmse,
             color = isotropy,
             linetype = sources,
             shape = sources)) +
  geom_path() + 
  geom_point() +
  # ylab('time rmse (days)') + scale_y_log10() +
  ylab('distance rmse (km)') + ylim(c(0, 1000)) +
  facet_wrap(~direction)
