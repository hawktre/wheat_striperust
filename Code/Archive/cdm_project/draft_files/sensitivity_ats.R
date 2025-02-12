source('two_source_preprocessing_fn.R')
seed <- 40820
nreps <- 2

out <- Reduce(rbind, lapply(1:nreps, function(rep){
  out_inner <- Reduce(rbind, lapply(0:1, function(loc){
    Reduce(rbind, lapply(25:30, function(w){
      ## ----preprocess, warning=F, message=F------------------------------------
      set.seed(seed + rep)
      wk <- w # week of year at which N source appears
      MI <- loc # N source location (1 = MI, 0 = NY)
      DATA <- preprocess_fn(wk, MI)
      rm(list = setdiff(ls(), c("DATA", "MI", "wk")))
      
      ## ----postprocess---------------------------------------------------------
      # select a year
      yr <- 2009
      
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
      
      ## ----fitting2------------------------------------------------------------
      # fit model
      fit_ios <- lm(response.0 ~ time, 
                    data = train_data)
      
      ## ----fitting3------------------------------------------------------------
      # fit model
      fit_aos <- lm(response.1 ~ time + sinb1.1 + sinb2.1, 
                    data = train_data)
      
      ## ----fitting4.init-------------------------------------------------------
      # initial value storage
      pr_z <- rep(0, length(train_idx))
      
      # determine which observations are before N source occurrence
      const_idx <- which(is.na(train_data$z.reg.2))
      update_idx <- which(!is.na(train_data$z.reg.2))
      
      # constrain weights for observations preceding N source
      pr_z[const_idx] <- 1
      
      # regress 1{S closer} on exp{response.1} and exp{response.2} for other observations
      z_glm <- glm(z.fit ~ z.reg.1 + z.reg.2, 
                   data = train_data[update_idx, ], 
                   family = binomial())
      
      
      # assign fitted probabilities for other observations
      pr_z[update_idx] <- fitted(z_glm)
      
      ## ----fitting4.em---------------------------------------------------------
      # storage for iteration change in weights
      pr_z_diff <- 1
      
      # iteration count
      iter <- 0
      
      # repeat until convergence
      while(pr_z_diff > 0.001) {
        # step 1: fit model 1 (S source)
        fit_ats1 <- lm(response.1 ~ time + sinb1.1 + sinb2.1, 
                       data = train_data,
                       weights = pr_z)
        
        
        # step 2: fit model 2 (N source)
        fit_ats2 <- lm(response.2 ~ time + sinb1.2 + sinb2.2, 
                       data = train_data[update_idx, ],
                       weights = (1 - pr_z[update_idx]))
        
        # step 3: E step for observations occurring after N source
        pr_z1 <- pr_z
        mean1 <- residuals(fit_ats1)/summary(fit_ats1)$sigma 
        mean2 <- residuals(fit_ats2)/summary(fit_ats2)$sigma
        pr_z1[update_idx] <- (dnorm(mean1[update_idx])*pr_z[update_idx])/(dnorm(mean1[update_idx])*pr_z[update_idx] + dnorm(mean2)*(1 - pr_z[update_idx]))
        
        # check for convergence
        pr_z_diff <- mean(abs(pr_z - pr_z1))
        
        # update weights
        pr_z <- pr_z1
        
        # increment iteration count
        iter <- iter + 1
      }
      
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
      betahat_ios <- coef(fit_ios)
      betahat_aos <- coef(fit_aos)
      betahat_ats1 <- coef(fit_ats1)
      betahat_ats2 <- coef(fit_ats2)
      
      # refit and use p_i's
      const_idx_full <- which(is.na(sub_data$z.reg.2))
      update_idx_full <- which(!is.na(sub_data$z.reg.2))
      z_glm_full <- glm(z.fit ~ z.reg.1 + z.reg.2, 
                        data = sub_data[update_idx_full, ], 
                        family = binomial())
      pr_z_full <- rep(0, n)
      pr_z_full[const_idx_full] <- 1
      pr_z_full[update_idx_full] <- fitted(z_glm_full)
      pr_z0_full <- pr_z_full
      pr_z_diff_full <- 1
      iter_full <- 0
      while(pr_z_diff_full > 0.001) {
        # step 1: fit model 1 (S source)
        fit_ats1_full <- lm(response.1 ~ time + sinb1.1 + sinb2.1, 
                            data = sub_data,
                            weights = pr_z0_full)
        
        
        # step 2: fit model 2 (N source)
        fit_ats2_full <- lm(response.2 ~ time + sinb1.2 + sinb2.2, 
                            data = sub_data[update_idx_full, ],
                            weights = (1 - pr_z0_full[update_idx_full]))
        
        # step 3: E step for observations occurring after N source
        pr_z1_full <- pr_z0_full
        mean1 <- residuals(fit_ats1_full)/summary(fit_ats1_full)$sigma 
        mean2 <- residuals(fit_ats2_full)/summary(fit_ats2_full)$sigma
        pr_z1_full[update_idx_full] <- (dnorm(mean1[update_idx_full])*pr_z0_full[update_idx_full])/(dnorm(mean1[update_idx_full])*pr_z0_full[update_idx_full] + dnorm(mean2)*(1 - pr_z0_full[update_idx_full]))
        
        # check for convergence
        pr_z_diff_full <- mean(abs(pr_z0_full - pr_z1_full))
        
        # update weights
        pr_z0_full <- pr_z1_full
        
        # increment iteration count
        iter_full <- iter_full + 1
      }
      
      # predictions of time of occurrence according to each model
      pred_df <- transmute(sub_data, 
                           t.hat.ats1 = (response.1 - betahat_ats1[1] - betahat_ats1[3]*sinb1.1 - betahat_ats1[4]*sinb2.1)/betahat_ats1[2],
                           t.hat.ats2 = (response.2 - betahat_ats2[1] - betahat_ats2[3]*sinb1.2 - betahat_ats2[4]*sinb2.2)/betahat_ats2[2],
                           t.hat.aos = (response.1 - betahat_aos[1] - betahat_aos[3]*sinb1.1 - betahat_aos[4]*sinb2.1)/betahat_aos[2],
                           t.hat.ios = (response.0 - betahat_ios[1])/betahat_ios[2],
                           time = time)
      
      which_pred_ats_refit <- round(pr_z0_full, 0) + (round(pr_z0_full, 0) == 0)*2
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
               week = wk,
               loc = MI,
               rep.num = rep) 
      
      return(out)
    }))
  }))
}))

out %>%
  group_by(week, loc) %>%
  summarize(meanRmseTest = mean(test),
            meanRmseTrain = mean(training),
            meanAdjR2 = mean(rsq)) %>%
  arrange(loc, week)
