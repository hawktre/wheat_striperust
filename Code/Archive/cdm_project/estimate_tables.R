source('preprocess_fn_v2.R')
source('fit_twosource_fn.R')
source('estimate_fn.R')
library(xtable)

## preprocess data
preprocess_out <- preprocess_fn()
DATA <- preprocess_out$data 

# select year
yrs <- 2008:2010 # possible years
j <- 2

est_tbls <- lapply(1:3, function(j){

if(j == 1){
  # source information (not year sensitive)
  src1_type <- 'observed'
  src2_type <- 'imputed'
  src3_type <- 'observed'
  
  wind_src1 <- preprocess_out$wind.s
  wind_src2 <- preprocess_out$wind.i3
  wind_src3 <- preprocess_out$wind.n
  src1_id <- filter(preprocess_out$sources.s,
                    year == yrs[j])$ID
  src2_id <- filter(preprocess_out$sources.i3,
                    year == yrs[j])$ID
  # src2_id <- NULL
  src3_id <- filter(preprocess_out$sources.n,
                    year == yrs[j])$ID
  # src3_id <- NULL
  
  
  # source information (year sensitive)
  source_df <- rbind(rename(filter(preprocess_out$sources.s, year == yrs[j]),
                            cc.long = cc.long.s,
                            cc.lat = cc.lat.s),
                     rename(filter(preprocess_out$sources.i3, year == yrs[j]),
                            cc.long = cc.long.i3,
                            cc.lat = cc.lat.i3),
                     rename(filter(preprocess_out$sources.n, year == yrs[j]),
                            cc.long = cc.long.n,
                            cc.lat = cc.lat.n))
  
  
  # subset data to year and calculate basis expansions and response variables
  sub_data <- filter(DATA, year(symptom.date) == yrs[j]) %>%
    rename(rho.src1 = rho.s,
           theta.src1 = theta.s,
           g.theta.src1 = g.theta.s,
           rho.src2 = rho.i3,
           theta.src2 = theta.i3,
           g.theta.src2 = g.theta.i3,
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
  
  estimate_out <- estimate_fn(sub_data, wind_src1, wind_src2, wind_src3, source_df, src1_id, NULL, src3_id)
}

if(j == 2){
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
  
  estimate_out <- estimate_fn(sub_data, wind_src1, wind_src2, wind_src3, source_df, src1_id, NULL, src3_id)
  
}

if(j == 3){
  # source information (not year sensitive)
  src1_type <- 'observed'
  src2_type <- 'imputed'
  src3_type <- 'imputed'
  
  wind_src1 <- preprocess_out$wind.s
  wind_src2 <- preprocess_out$wind.i3
  wind_src3 <- preprocess_out$wind.i1
  src1_id <- filter(preprocess_out$sources.s,
                    year == yrs[j])$ID
  src2_id <- filter(preprocess_out$sources.i3,
                    year == yrs[j])$ID
  # src2_id <- NULL
  src3_id <- filter(preprocess_out$sources.i1,
                    year == yrs[j])$ID
  # src3_id <- NULL
  
  
  # source information (year sensitive)
  source_df <- rbind(rename(filter(preprocess_out$sources.s, year == yrs[j]),
                            cc.long = cc.long.s,
                            cc.lat = cc.lat.s),
                     rename(filter(preprocess_out$sources.i3, year == yrs[j]),
                            cc.long = cc.long.i3,
                            cc.lat = cc.lat.i3),
                     rename(filter(preprocess_out$sources.i1, year == yrs[j]),
                            cc.long = cc.long.i1,
                            cc.lat = cc.lat.i1))
  
  
  # subset data to year and calculate basis expansions and response variables
  sub_data <- filter(DATA, year(symptom.date) == yrs[j]) %>%
    rename(rho.src1 = rho.s,
           theta.src1 = theta.s,
           g.theta.src1 = g.theta.s,
           rho.src2 = rho.i3,
           theta.src2 = theta.i3,
           g.theta.src2 = g.theta.i3,
           rho.src3 = rho.i1,
           theta.src3 = theta.i1,
           g.theta.src3 = g.theta.i1) %>%
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
  
  estimate_out <- estimate_fn(sub_data, wind_src1, wind_src2, wind_src3, source_df, src1_id, NULL, NULL)
  
}

return(estimate_out)
})

est_tbls[[3]] %>%
  xtable(digits = 3)
