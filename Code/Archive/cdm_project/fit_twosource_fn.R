fit_ts <- function(source1, source2, initdata, isotropy){
  
  if(isotropy == T){
    model_formula <- formula('response ~ time')
  }else{
    model_formula <- formula('response ~ time + sinb1 + sinb2')
  }
  
  ## ----fitting4.init-------------------------------------------------------
  # initial value storage
  pr_z0 <- rep(0, nrow(source1))
  
  # determine which observations are before second source occurrence
  const_idx <- which(is.na(source2$rho))
  update_idx <- which(!is.na(source2$rho))
  
  # constrain weights for observations preceding N source
  pr_z0[const_idx] <- 1
  
  # regress 1{S closer} on exp{response.1} and exp{response.2} for other observations
  z_glm <- glm(z.fit ~ z.reg.1 + z.reg.2, 
               data = initdata[update_idx, ], 
               family = binomial())
  
  
  # assign fitted probabilities for other observations
  pr_z0[update_idx] <- fitted(z_glm)
  
  ## ----fitting4.em---------------------------------------------------------
  # storage for iteration change in weights
  pr_z_diff <- 1
  
  # iteration count
  iter <- 0
  
  # initial value
  pr_zk <- pr_z0
  
  # repeat until convergence
  while(pr_z_diff > 0.001) {
    # step 1: fit model 1
    fit_1 <- lm(model_formula, 
                   data = source1,
                   weights = pr_zk)
    
    
    # step 2: fit model 2 (N source)
    fit_2 <- lm(model_formula, 
                   data = source2[update_idx, ],
                   weights = (1 - pr_zk[update_idx]))
    
    # step 3: E step for observations occurring after N source
    pr_zkplus1 <- pr_zk
    mean1 <- residuals(fit_1)/summary(fit_1)$sigma 
    mean2 <- residuals(fit_2)/summary(fit_2)$sigma
    pr_zkplus1[update_idx] <- (dnorm(mean1[update_idx])*pr_zk[update_idx])/(dnorm(mean1[update_idx])*pr_zk[update_idx] + dnorm(mean2)*(1 - pr_zk[update_idx]))
    
    # check for convergence
    pr_z_diff <- mean(abs(pr_zk - pr_zkplus1), na.rm = F)
    
    # update weights
    pr_zk <- pr_zkplus1
    
    # increment iteration count
    iter <- iter + 1
  }
  
  out <- list(fit1 = fit_1,
              fit2 = fit_2,
              prz = pr_zk)
  
  return(out)
}
