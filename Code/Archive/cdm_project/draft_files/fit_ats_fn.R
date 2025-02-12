fit_ats <- function(train_data){

## ----fitting4.init-------------------------------------------------------
# initial value storage
pr_z0 <- rep(0, nrow(train_data))

# determine which observations are before N source occurrence
const_idx <- which(is.na(train_data$z.reg.2))
update_idx <- which(!is.na(train_data$z.reg.2))

# constrain weights for observations preceding N source
pr_z0[const_idx] <- 1

# regress 1{S closer} on exp{response.1} and exp{response.2} for other observations
z_glm <- glm(z.fit ~ z.reg.1 + z.reg.2, 
             data = train_data[update_idx, ], 
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
  # step 1: fit model 1 (S source)
  fit_ats1 <- lm(aresponse.1 ~ time + sinb1.1 + sinb2.1, 
                 data = train_data,
                 weights = pr_zk)
  
  
  # step 2: fit model 2 (N source)
  fit_ats2 <- lm(aresponse.2 ~ time + sinb1.2 + sinb2.2, 
                 data = train_data[update_idx, ],
                 weights = (1 - pr_zk[update_idx]))
  
  # step 3: E step for observations occurring after N source
  pr_zkplus1 <- pr_zk
  mean1 <- residuals(fit_ats1)/summary(fit_ats1)$sigma 
  mean2 <- residuals(fit_ats2)/summary(fit_ats2)$sigma
  pr_zkplus1[update_idx] <- (dnorm(mean1[update_idx])*pr_zk[update_idx])/(dnorm(mean1[update_idx])*pr_zk[update_idx] + dnorm(mean2)*(1 - pr_zk[update_idx]))
  
  # check for convergence
  pr_z_diff <- mean(abs(pr_zk - pr_zkplus1), na.rm = F)

  # update weights
  pr_zk <- pr_zkplus1
  
  # increment iteration count
  iter <- iter + 1
}

out <- list(ats1 = fit_ats1,
            ats2 = fit_ats2,
            prz = pr_zk)

return(out)
}
