source('preprocess_fn_v2.R')
source('fit_twosource_fn.R')

## preprocess data
preprocess_out <- preprocess_fn()
DATA <- preprocess_out$data 

# select year
yrs <- 2008:2010 # possible years
j <- 2


# source information (not year sensitive)
src1_type <- 'observed'
src2_type <- 'imputed'
src3_type <- 'observed'

wind_src1 <- preprocess_out$wind.s
wind_src2 <- preprocess_out$wind.i4
wind_src3 <- preprocess_out$wind.n
src1_id <- filter(preprocess_out$sources.s,
                  year == yrs[j])$ID
src2_id <- filter(preprocess_out$sources.i4,
                  year == yrs[j])$ID
# src2_id <- NULL
src3_id <- filter(preprocess_out$sources.n,
                  year == yrs[j])$ID
# src3_id <- NULL


# source information (year sensitive)
source_df <- rbind(rename(filter(preprocess_out$sources.s, year == yrs[j]),
                          cc.long = cc.long.s,
                          cc.lat = cc.lat.s),
                   rename(filter(preprocess_out$sources.i4, year == yrs[j]),
                          cc.long = cc.long.i4,
                          cc.lat = cc.lat.i4),
                   rename(filter(preprocess_out$sources.n, year == yrs[j]),
                          cc.long = cc.long.n,
                          cc.lat = cc.lat.n))


# subset data to year and calculate basis expansions and response variables
sub_data <- filter(DATA, year(symptom.date) == yrs[j]) %>%
  rename(rho.src1 = rho.s,
         theta.src1 = theta.s,
         g.theta.src1 = g.theta.s,
         rho.src2 = rho.i4,
         theta.src2 = theta.i4,
         g.theta.src2 = g.theta.i4,
         rho.src3 = rho.n,
         theta.src3 = theta.n,
         g.theta.src3 = g.theta.n) %>%
  mutate(month = month(symptom.date),
         day = day(symptom.date),
         time = yday(symptom.date),
         iresponse.src1 = log(1 + rho.src1), # isotropic response for first source
         iresponse.src2 = log(1 + rho.src2), # isotropic response for second source
         iresponse.src3 = log(1 + rho.src3), # isotropic response for second source
         aresponse.src1 = log(1 + rho.src1/g.theta.src1), # anisotropic response for first source
         aresponse.src2 = log(1 + rho.src2/g.theta.src2), # anisotropic response for second source
         aresponse.src3 = log(1 + rho.src3/g.theta.src3), # anisotropic response for second source
         sinb1.src1 = sin(as.numeric(theta.src1)), # sin basis 1 for first source
         sinb2.src1 = sin(as.numeric(theta.src1) + pi/4), # sin basis 2 for first source
         sinb1.src2 = sin(as.numeric(theta.src2)), # sin basis 1 for second source
         sinb2.src2 = sin(as.numeric(theta.src2) + pi/4), # sin basis 2 for second source
         sinb1.src3 = sin(as.numeric(theta.src3)), # sin basis 1 for second source
         sinb2.src3 = sin(as.numeric(theta.src3) + pi/4), # sin basis 2 for second source
         z.reg.src1 = rho.src1/g.theta.src1, # first source regressor to initialize weights
         z.reg.src2 = rho.src2/g.theta.src2, # second source regressor to initialize weights
         z.reg.src3 = rho.src3/g.theta.src3, # second source regressor to initialize weights
         z.fit.12 = rho.src1 < rho.src2, # is first source closer than second source?
         z.fit.13 = rho.src1 < rho.src3) # is first source closer than second source?

sub_data <- sub_data %>%
  group_by(ID) %>%
  sample_n(1)



## source-specific datasets

source1_data <- filter(sub_data, !(ID %in% c(src1_id, src2_id, src3_id))) %>%
  rename(rho = rho.src1, 
         theta = theta.src1,
         g.theta = g.theta.src1,
         sinb1 = sinb1.src1, 
         sinb2 = sinb2.src1, 
         response = aresponse.src1) %>%
  select(rho, theta, g.theta, sinb1, sinb2, response, time)
source2_data <- filter(sub_data, !(ID %in% c(src1_id, src2_id, src3_id))) %>%
  rename(rho = rho.src2, 
         theta = theta.src2,
         g.theta = g.theta.src2,
         sinb1 = sinb1.src2, 
         sinb2 = sinb2.src2, 
         response = aresponse.src2) %>%
  select(rho, theta, g.theta, sinb1, sinb2, response, time)
source3_data <- filter(sub_data, !(ID %in% c(src1_id, src2_id, src3_id))) %>%
  rename(rho = rho.src3, 
         theta = theta.src3,
         g.theta = g.theta.src3,
         sinb1 = sinb1.src3, 
         sinb2 = sinb2.src3, 
         response = aresponse.src3) %>%
  select(rho, theta, g.theta, sinb1, sinb2, response, time)

save(list = c('source1_data', 'source2_data', 'source3_data', 'sub_data', 'source_df', 'wind_src1', 'wind_src2', 'wind_src3'), file = '3source_data.RData')

## model fitting
load('3source_data.RData')

model_formula <- formula('response ~ time + sinb1 + sinb2')

pr_z0 <- matrix(c(4/10, 3/10, 3/10),   
               nrow = nrow(source1_data),
               ncol = 3, 
               byrow = T)

const2_idx <- which(is.na(source2_data$rho))
update2_idx <- which(!is.na(source2_data$rho))
pr_z0[const2_idx, 2] <- 0

const3_idx <- which(is.na(source3_data$rho))
update3_idx <- which(!is.na(source3_data$rho))
pr_z0[const3_idx, 3] <- 0

pr_z0 <- pr_z0/apply(pr_z0, 1, sum)
const1_idx <- which(pr_z0[, 1] == 1)
update1_idx <- which(pr_z0[, 1] != 1)

pr_z <- pr_z0
pr_z_diff <- 1
iter <- 1

while(pr_z_diff > exp(-6)){
  pr_z_old <- pr_z
  
  # step 1: fit model 1
  fit_1 <- lm(model_formula, 
              data = source1_data,
              weights = pr_z[, 1])
  
  
  # step 2: fit model 2 (N source)
  fit_2 <- lm(model_formula, 
              data = source2_data,
              weights = pr_z[, 2])
  
  # step 2: fit model 3 (N source)
  fit_3 <- lm(model_formula, 
              data = source3_data,
              weights = pr_z[, 3])
  
  
  mean1 <- fit_1$residuals/summary(fit_1)$sigma^2
  mean2 <- fit_2$residuals/summary(fit_2)$sigma^2
  mean3 <- fit_3$residuals/summary(fit_3)$sigma^2
  
  pr_z[, 1] <- dnorm(mean1)*pr_z[, 1]
  pr_z[update2_idx, 2] <- dnorm(mean2)*pr_z[update2_idx, 2]
  pr_z[update3_idx, 3] <- dnorm(mean3)*pr_z[update3_idx, 3]
  
  
  pr_z[const1_idx, 1] <- 1
  pr_z[const2_idx, 2] <- 0
  pr_z[const3_idx, 3] <- 0
  
  pr_z <- pr_z/apply(pr_z, 1, sum)
  pr_z_diff <- mean(abs(pr_z - pr_z_old))

  iter <- iter + 1
  print(iter)
  }
