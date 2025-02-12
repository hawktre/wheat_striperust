## ----sensitivity, echo = F, warning = F, message = F, cache = T----------
source('two_source_preprocessing_fn.R') # function to preprocess data given second source location and time
source('fit_ats_fn.R') # function to fit ATS model
seed <- 41420 # rng seed
nreps <- 20 # number of repititions 

out <- Reduce(rbind, lapply(1:nreps, function(rep){
  out_inner1 <- Reduce(rbind, lapply(0:1, function(loc){
    out_inner2 <- Reduce(rbind, lapply(25:30, function(w){
      ## ----preprocess, warning=F, message=F------------------------------------
      set.seed(seed + rep)
      DATA <- preprocess_fn(w, loc)
      
      out_inner3 <- Reduce(rbind, lapply(2008:2011, function(yr){
        
        print(c(loc, w, yr, rep))
        ## ----postprocess---------------------------------------------------------
        # subset data to year and calculate basis expansions and response variables
        sub_data <- filter(DATA, year(report_date) == yr) %>%
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
        
        # split data into training and test sets
        n <- nrow(sub_data)
        n_train <- floor(0.8*n)
        train_idx <- sample(1:n, n_train)
        
        train_data <- sub_data[train_idx, ]
        test_data <- sub_data[-train_idx, ]
        
        fit <- try(fit_ats(train_data), silent = T)
        
        if(is.null(attr(fit, 'class'))){ 
        
          fit_ats1 <- fit$ats1
          fit_ats2 <- fit$ats2
          pr_z <- fit$prz
          
        resids <- train_data %>%
          mutate(pred1 = predict(fit_ats1, train_data),
                 pred2 = predict(fit_ats2, train_data)) %>%
          transmute(resid1 = response.1 - pred1,
                    resid2 = response.2 - pred2,
                    resp1 = response.1,
                    resp2 = response.2) %>%
          transmute(resid = if_else(round(pr_z, 0) == 1, resid1, resid2),
                    resp = if_else(round(pr_z, 0) == 1, resp1, resp2))
        
        # compute adjusted R2
        sse <- sum(resids$resid^2)
        sst <- sum((resids$resp - mean(resids$resp))^2)
        rsq.ats <- 1 - (sse/sst)*((n_train - 1)/(n_train - 8 - 1))
        
        
        ## ----predicitons1--------------------------------------------------------
        # extract estimates
        betahat_ats1 <- coef(fit_ats1)
        betahat_ats2 <- coef(fit_ats2)
        
        # refit and use p_i's
        fit_full <- fit_ats(sub_data)
        pr_z_full <- fit_full$prz
        
        # predictions of time of occurrence according to each model
        pred_df <- transmute(sub_data, 
                             t.hat.ats1 = (response.1 - betahat_ats1[1] - betahat_ats1[3]*sinb1.1 - betahat_ats1[4]*sinb2.1)/betahat_ats1[2],
                             t.hat.ats2 = (response.2 - betahat_ats2[1] - betahat_ats2[3]*sinb1.2 - betahat_ats2[4]*sinb2.2)/betahat_ats2[2],
                             time = time)
        
        which_pred_ats_refit <- round(pr_z_full, 0) + (round(pr_z_full, 0) == 0)*2
        pred_df <- cbind(pred_df,
                         which.pred.ats.refit = which_pred_ats_refit) %>%
          mutate(t.hat.ats.refit = if_else(which.pred.ats.refit == 1, 
                                           t.hat.ats1, 
                                           t.hat.ats2)) %>%
          rename(t.hat.ats = t.hat.ats.refit,
                 which.pred.ats = which.pred.ats.refit) %>%
          select(time, 
                 t.hat.ats) %>%
          rename(ats = t.hat.ats) 
        
        pred_df$subset <- 'test'
        pred_df$subset[train_idx] <- 'training'
        
        out <- pred_df %>%
          mutate(sq.resid.ats = (time - ats)^2) %>%
          group_by(subset) %>%
          summarize(ats = sqrt(mean(sq.resid.ats))) %>%
          spread(subset, ats) %>%
          mutate(rsq = rsq.ats,
                 week = w,
                 MI = loc,
                 rep.num = rep,
                 year = yr) 
        }else{
          out <- data.frame(week = w,
                            MI = loc,
                            rep.num = rep,
                            year = yr,
                            rsq = NA,
                            test = NA,
                            training = NA)
        }
        
        return(out)
      }))
    }))
    return(out_inner2)
  }))
  return(out_inner1)
}))


save(out, file = 'second_source_sensitivity.RData')

#####################
load('second_source_sensitivity.RData')

summary <- out %>%
  group_by(year, week, MI) %>%
  summarize(#medRmseTest = median(test, na.rm = T),
            # medRmseTrain = median(training, na.rm = T),
            # medAdjR2 = median(rsq, na.rm = T),
            meanRmseTest = mean(test, na.rm = T),
            meanRmseTrain = mean(training, na.rm = T),
            meanAdjR2 = mean(rsq, na.rm = T),
            sdRmseTest = sd(test, na.rm = T),
            sdRmseTrain = sd(training, na.rm = T),
            sdAdjR2 = sd(rsq, na.rm = T),
            numNA = sum(is.na(test))) %>%
  arrange(year, MI, week) %>%
  ungroup()

# 2008: week 25, michigan
tbl2008 <- filter(summary, year == 2008) %>%
  as.data.frame() 

arrange(tbl2008, desc(`Mean Adj. R^2`)) # wk 29, MI
arrange(tbl2008, desc(`SD Adj. R^2`)) # wk 29, MI
arrange(tbl2008, (`Mean Test RMSE`)) # wk 29, MI
arrange(tbl2008, (`Mean Train RMSE`)) # wk 30, MI
arrange(tbl2008, (`SD Test RMSE`)) # wk 29, MI
arrange(tbl2008, (`SD Train RMSE`)) # wk 30, MI

colnames(tbl2008) <- c('Year',
                       'Week',
                       'Location',
                       'Mean Test RMSE',
                       'Mean Train RMSE',
                       'Mean Adj. R^2',
                       'SD Test RMSE',
                       'SD Train RMSE',
                       'SD Adj. R^2',
                       'Fit failures')

tbl2008 <- tbl2008[, c(1:3, 4, 7, 5, 8, 6, 9, 10)]
xtable(tbl2008)

# 2009: week 26, michigan
tbl2009 <- filter(summary, year == 2009) %>%
  as.data.frame()

colnames(tbl2009) <- c('Year',
                       'Week',
                       'Location',
                       'Mean Test RMSE',
                       'Mean Train RMSE',
                       'Mean Adj. R^2',
                       'SD Test RMSE',
                       'SD Train RMSE',
                       'SD Adj. R^2',
                       'Fit failures')

top_n(tbl2009, 1, (`Mean Adj. R^2`))[, 1:3] # wk 26, NY
top_n(tbl2009, 1, desc(`SD Adj. R^2`))[, 1:3] # wk 29, MI
top_n(tbl2009, 1, desc(`Mean Test RMSE`))[, 1:3] # wk 26, MI
top_n(tbl2009, 1, desc(`Mean Train RMSE`))[, 1:3] # wk 26, MI
top_n(tbl2009, 1, desc(`SD Test RMSE`))[, 1:3] # wk 27, MI
top_n(tbl2009, 1, desc(`SD Train RMSE`))[, 1:3] # wk 27, MI

tbl2009 <- tbl2009[, c(1:3, 4, 7, 5, 8, 6, 9, 10)]
xtable(tbl2009)


# 2010: week 29, new york
tbl2010 <- filter(summary, year == 2010) %>%
  as.data.frame()

colnames(tbl2010) <- c('Year',
                       'Week',
                       'Location',
                       'Mean Test RMSE',
                       'Mean Train RMSE',
                       'Mean Adj. R^2',
                       'SD Test RMSE',
                       'SD Train RMSE',
                       'SD Adj. R^2',
                       'Fit failures')

top_n(tbl2010, 1, (`Mean Adj. R^2`))[, 1:3] # wk 30, NY
top_n(tbl2010, 1, desc(`SD Adj. R^2`))[, 1:3] # wk 30, MI
top_n(tbl2010, 1, desc(`Mean Test RMSE`))[, 1:3] # wk 29, NY
top_n(tbl2010, 1, desc(`Mean Train RMSE`))[, 1:3] # wk 30, MI
top_n(tbl2010, 1, desc(`SD Test RMSE`))[, 1:3] # wk 29, NY
top_n(tbl2010, 1, desc(`SD Train RMSE`))[, 1:3] # wk 30, MI

tbl2010 <- tbl2010[, c(1:3, 4, 7, 5, 8, 6, 9, 10)]
xtable(tbl2010)

# 2011: week 29, michigan
tbl2011 <- filter(summary, year == 2011) %>%
  as.data.frame()

colnames(tbl2011) <- c('Year',
                       'Week',
                       'Location',
                       'Mean Test RMSE',
                       'Mean Train RMSE',
                       'Mean Adj. R^2',
                       'SD Test RMSE',
                       'SD Train RMSE',
                       'SD Adj. R^2',
                       'Fit failures')
top_n(tbl2011, 1, (`Mean Adj. R^2`))[, 1:3] # wk 26, MI
top_n(tbl2011, 1, desc(`SD Adj. R^2`))[, 1:3] # wk 25, MI
top_n(tbl2011, 1, desc(`Mean Test RMSE`))[, 1:3] # wk 30, MI
top_n(tbl2011, 1, desc(`Mean Train RMSE`))[, 1:3] # wk 30, NY
top_n(tbl2011, 1, desc(`SD Test RMSE`))[, 1:3] # wk 30, MI
top_n(tbl2011, 1, desc(`SD Train RMSE`))[, 1:3] # wk 30, MI

tbl2011 <- tbl2011[, c(1:3, 4, 7, 5, 8, 6, 9, 10)]
xtable(tbl2011)
