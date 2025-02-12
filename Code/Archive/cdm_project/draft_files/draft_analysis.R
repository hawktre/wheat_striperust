source('two_source_preprocessing_fn.R')
source('fit_ats_fn.R')
seed <- 41420

yrs <- 2008:2011
wks <- c(25, 26, 29, 29)
locs <- c(0, 1, 0, 1)

nReps <- 20

rep_out <- lapply(1:nReps, function(rep){
  set.seed(seed + rep)
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
    training_prop <- 0.8
    n <- nrow(sub_data)
    n_train <- floor(training_prop*n)
    train_idx <- sample(1:n, n_train)
    
    train_data <- sub_data[train_idx, ]
    test_data <- sub_data[-train_idx, ]
    
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
    
    ## ----predicitons---------------------------------------------------------
    # predictions of time of occurrence according to each model
    sub_data <- mutate(sub_data, 
                       t.hat.ats1 = (response.1 - betahat_ats1[1] - betahat_ats1[3]*sinb1.1 - betahat_ats1[4]*sinb2.1)/betahat_ats1[2],
                       t.hat.ats2 = (response.2 - betahat_ats2[1] - betahat_ats2[3]*sinb1.2 - betahat_ats2[4]*sinb2.2)/betahat_ats2[2],
                       t.hat.aos = (response.1 - betahat_aos[1] - betahat_aos[3]*sinb1.1 - betahat_aos[4]*sinb2.1)/betahat_aos[2],
                       t.hat.ios = (response.0 - betahat_ios[1])/betahat_ios[2]) 
    
    # refit ats model to estimate p_i's for test data
    fit_ats_out_full <- fit_ats(sub_data)
    pr_z_full <- fit_ats_out_full$prz
    
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
      mutate(year = yrs[j],
             rep.num = rep)
    pred_df$subset <- 'test data'
    pred_df$subset[train_idx] <- 'training data'
    
    out <- list(preds = pred_df,
                ests = estimates,
                rsqs = adj.rs)
    
    return(out)
  })
  out <- list(rsqs = Reduce(rbind, lapply(1:4, function(yr){yr_out[[yr]]$rsqs})),
              ests = Reduce(rbind, lapply(1:4, function(yr){yr_out[[yr]]$ests})),
              preds = Reduce(rbind, lapply(1:4, function(yr){yr_out[[yr]]$preds})))
  return(out)
})

out <- list(rsqs = Reduce(rbind, lapply(1:nReps, function(rep){rep_out[[rep]]$rsqs})),
            ests = Reduce(rbind, lapply(1:nReps, function(rep){rep_out[[rep]]$ests})),
            preds = Reduce(rbind, lapply(1:nReps, function(rep){rep_out[[rep]]$preds})))

save(out, file = 'draft_analysis_holdout.RData')

##########
library(xtable)
# average adjusted R^2 for each model over 20 training/test splits
tbl1 <- out$rsqs %>%
  as_tibble() %>%
  group_by(year) %>%
  summarize_all(c('mean', 'sd')) %>%
  select(-starts_with('rep')) %>%
  as.data.frame()

colnames(tbl1) <- c('Year',
                    'IOS mean',
                    'AOS mean',
                    'ATS mean',
                    'IOS sd',
                    'AOS sd',
                    'ATS sd')

tbl1 <- tbl1[, c(1, 2, 5, 3, 6, 4, 7)]
xtable(tbl1, digits = 3)

# average estimates and variability for each model across 20 training/test splits
tbl2a <- out$ests %>%
  as_tibble() %>%
  group_by(year) %>%
  summarize_all(c('mean', 'sd')) %>%
  select(starts_with("ios"), year) %>%
  as.data.frame()

colnames(tbl2a) <- c('Intercept mean',
                     'Time mean',
                     'Intercept sd',
                     'Time sd',
                     'Year')

tbl2a <- tbl2a[, c(5, 1, 3, 2, 4)]

tbl2b <- out$ests %>%
  as_tibble() %>%
  group_by(year) %>%
  summarize_all(c('mean', 'sd')) %>%
  select(starts_with("aos"), year) %>%
  as.data.frame()

colnames(tbl2b) <- c('Intercept mean',
                     'Time mean',
                     'Basis 1 mean',
                     'Basis 2 mean',
                     'Intercept sd',
                     'Time sd',
                     'Basis 1 sd',
                     'Basis 2 sd',
                     'Year')

tbl2b <- tbl2b[, c(9, 1, 5, 2, 6, 3, 7, 4, 8)]

tbl2c <- out$ests %>%
  as_tibble() %>%
  group_by(year) %>%
  summarize_all(c('mean', 'sd')) %>%
  select(starts_with("ats"), year) %>%
  as.data.frame()

colnames(tbl2c) <- c('Intercept mean',
                     'Time mean',
                     'Basis 1 mean',
                     'Basis 2 mean',
                     'N Intercept mean',
                     'N Time mean',
                     'N Basis 1 mean',
                     'N Basis 2 mean',
                     'Intercept sd',
                     'Time sd',
                     'Basis 1 sd',
                     'Basis 2 sd',
                     'N Intercept sd',
                     'N Time sd',
                     'N Basis 1 sd',
                     'N Basis 2 sd',
                     'Year')

tbl2c <- tbl2c[, c(17, 1, 9, 2, 10, 3, 11, 4, 12, 5, 13, 6, 14, 7, 15, 8, 16)]


tbl2a <- cbind(tbl2a, matrix(rep(NA, 4*12), nrow = 4))
tbl2b <- cbind(tbl2b, matrix(rep(NA, 4*8), nrow = 4))
colnames(tbl2a) <- colnames(tbl2b) <- colnames(tbl2c)

tbl2 <- rbind(IOS = tbl2a, AOS = tbl2b, ATS = tbl2c)
xtable(tbl2, digits = 3)

# average error on time predictions and variability across 20 training/test splits
tbl3 <- out$preds %>%
  as_tibble() %>%
  mutate(sq.resid.ios = (time - ios)^2,
         sq.resid.aos = (time - aos)^2,
         sq.resid.ats = (time - ats)^2,
         resid.ios = (time - ios),
         resid.aos = (time - aos),
         resid.ats = (time - ats)) %>%
  group_by(year, rep.num, subset) %>%
  summarize(rmse.ios = sqrt(mean(sq.resid.ios)),
            rmse.aos = sqrt(mean(sq.resid.aos)),
            rmse.ats = sqrt(mean(sq.resid.ats)),
            mres.ios = mean(resid.ios),
            mres.aos = mean(resid.aos),
            mres.ats = mean(resid.ats)) %>%
  ungroup() %>%
  group_by(subset, year) %>%
  summarize(m.rmse.ios = mean(rmse.ios),
            m.rmse.aos = mean(rmse.aos),
            m.rmse.ats = mean(rmse.ats),
            sd.rmse.ios = sd(rmse.ios),
            sd.rmse.aos = sd(rmse.aos),
            sd.rmse.ats = sd(rmse.ats),
            m.mres.ios = mean(mres.ios),
            m.mres.aos = mean(mres.aos),
            m.mres.ats = mean(mres.ats)) %>%
  as.data.frame()

colnames(tbl3) <- c('Data partition',
                    'Year',
                    'IOS mean RMSE',
                    'AOS mean RMSE',
                    'ATS mean RMSE',
                    'IOS sd RMSE',
                    'AOS sd RMSE',
                    'ATS sd RMSE',
                    'IOS mean error',
                    'AOS mean error',
                    'ATS mean error')

tbl3 <- tbl3[, c(1:2, 3, 6, 9, 4, 7, 10, 5, 8, 11)]

xtable(tbl3)
