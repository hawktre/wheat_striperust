source('preprocess_fn_v2.R')
source('fit_twosource_fn.R')
source('tbl_fn.R')

## preprocess data
preprocess_out <- preprocess_fn()
DATA <- preprocess_out$data 

# select year
yrs <- 2008:2010 # possible years
j <- 3

sw_tbl <- Reduce(rbind, lapply(1:3, function(j){
# source information 
src1_type <- 'observed'
wind_src1 <- preprocess_out$wind.s
src1_id <- filter(preprocess_out$sources.s,
                  year == yrs[j])$ID


out <- Reduce(rbind, lapply(1:3, function(loc){
# subset data to year and calculate basis expansions and response variables
if(loc == 1){
sub_data <- filter(DATA, year(symptom.date) == yrs[j]) %>%
  rename(rho.src1 = rho.s,
         theta.src1 = theta.s,
         g.theta.src1 = g.theta.s,
         rho.src2 = rho.w,
         theta.src2 = theta.w,
         g.theta.src2 = g.theta.w)

wind_src2 <- preprocess_out$wind.w

src2_id <- filter(preprocess_out$sources.w,
                    year == yrs[j])$ID

src2_type <- 'observed'

source_df <- rbind(rename(filter(preprocess_out$sources.s, year == yrs[j]),
                          cc.long = cc.long.s,
                          cc.lat = cc.lat.s),
                   rename(filter(preprocess_out$sources.w, year == yrs[j]),
                          cc.long = cc.long.w,
                          cc.lat = cc.lat.w))
}

if(loc == 2){
  sub_data <- filter(DATA, year(symptom.date) == yrs[j]) %>%
    rename(rho.src1 = rho.s,
           theta.src1 = theta.s,
           g.theta.src1 = g.theta.s,
           rho.src2 = rho.i3,
           theta.src2 = theta.i3,
           g.theta.src2 = g.theta.i3) 
  
  wind_src2 <- preprocess_out$wind.i3
  
  src2_id <- NULL
  
  src2_type <- 'imputed'
  
  source_df <- rbind(rename(filter(preprocess_out$sources.s, year == yrs[j]),
                            cc.long = cc.long.s,
                            cc.lat = cc.lat.s),
                     rename(filter(preprocess_out$sources.i3, year == yrs[j]),
                            cc.long = cc.long.i3,
                            cc.lat = cc.lat.i3))
}

if(loc == 3){
  sub_data <- filter(DATA, year(symptom.date) == yrs[j]) %>%
    rename(rho.src1 = rho.s,
           theta.src1 = theta.s,
           g.theta.src1 = g.theta.s,
           rho.src2 = rho.i4,
           theta.src2 = theta.i4,
           g.theta.src2 = g.theta.i4) 
  
  wind_src2 <- preprocess_out$wind.i4
  
  src2_id <- NULL
  
  src2_type <- 'imputed'
  
  source_df <- rbind(rename(filter(preprocess_out$sources.s, year == yrs[j]),
                            cc.long = cc.long.s,
                            cc.lat = cc.lat.s),
                     rename(filter(preprocess_out$sources.i4, year == yrs[j]),
                            cc.long = cc.long.i4,
                            cc.lat = cc.lat.i4))
}

sub_data <- sub_data %>%
  mutate(month = month(symptom.date),
         day = day(symptom.date),
         time = yday(symptom.date),
         iresponse.src1 = log(1 + rho.src1), # isotropic response for first source
         iresponse.src2 = log(1 + rho.src2), # isotropic response for second source
         aresponse.src1 = log(1 + rho.src1/g.theta.src1), # anisotropic response for first source
         aresponse.src2 = log(1 + rho.src2/g.theta.src2), # anisotropic response for second source
         sinb1.src1 = sin(as.numeric(theta.src1)), # sin basis 1 for first source
         sinb2.src1 = sin(as.numeric(theta.src1) + pi/4), # sin basis 2 for first source
         sinb1.src2 = sin(as.numeric(theta.src2)), # sin basis 1 for second source
         sinb2.src2 = sin(as.numeric(theta.src2) + pi/4), # sin basis 2 for second source
         z.reg.src1 = rho.src1/g.theta.src1, # first source regressor to initialize weights
         z.reg.src2 = rho.src2/g.theta.src2, # second source regressor to initialize weights
         z.fit = rho.src1 < rho.src2) # is first source closer than second source?

sub_data <- sub_data %>%
  group_by(ID) %>%
  sample_n(1)

out <- tbl_fn(sub_data, 
  wind_src1, wind_src2, 
  source_df, src1_id, src2_id,
  src1_type, src2_type)

return(out)}))

return(out)}))

sw_tbl

n_tbl <- Reduce(rbind, lapply(1:3, function(j){
  # source information 
  src1_type <- 'observed'
  wind_src1 <- preprocess_out$wind.s
  src1_id <- filter(preprocess_out$sources.s,
                    year == yrs[j])$ID
  
  
  out <- Reduce(rbind, lapply(1:3, function(loc){
    # subset data to year and calculate basis expansions and response variables
    if(loc == 1){
      sub_data <- filter(DATA, year(symptom.date) == yrs[j]) %>%
        rename(rho.src1 = rho.s,
               theta.src1 = theta.s,
               g.theta.src1 = g.theta.s,
               rho.src2 = rho.n,
               theta.src2 = theta.n,
               g.theta.src2 = g.theta.n)
      
      wind_src2 <- preprocess_out$wind.n
      
      src2_id <- filter(preprocess_out$sources.n,
                        year == yrs[j])$ID
      
      src2_type <- 'observed'
      
      source_df <- rbind(rename(filter(preprocess_out$sources.s, year == yrs[j]),
                                cc.long = cc.long.s,
                                cc.lat = cc.lat.s),
                         rename(filter(preprocess_out$sources.n, year == yrs[j]),
                                cc.long = cc.long.n,
                                cc.lat = cc.lat.n))
    }
    
    if(loc == 2){
      sub_data <- filter(DATA, year(symptom.date) == yrs[j]) %>%
        rename(rho.src1 = rho.s,
               theta.src1 = theta.s,
               g.theta.src1 = g.theta.s,
               rho.src2 = rho.i2,
               theta.src2 = theta.i2,
               g.theta.src2 = g.theta.i2) 
      
      wind_src2 <- preprocess_out$wind.i2
      
      src2_id <- NULL
      
      src2_type <- 'imputed'
      
      source_df <- rbind(rename(filter(preprocess_out$sources.s, year == yrs[j]),
                                cc.long = cc.long.s,
                                cc.lat = cc.lat.s),
                         rename(filter(preprocess_out$sources.i2, year == yrs[j]),
                                cc.long = cc.long.i2,
                                cc.lat = cc.lat.i2))
    }
    
    if(loc == 3){
      sub_data <- filter(DATA, year(symptom.date) == yrs[j]) %>%
        rename(rho.src1 = rho.s,
               theta.src1 = theta.s,
               g.theta.src1 = g.theta.s,
               rho.src2 = rho.i1,
               theta.src2 = theta.i1,
               g.theta.src2 = g.theta.i1) 
      
      wind_src2 <- preprocess_out$wind.i1
      
      src2_id <- NULL
      
      src2_type <- 'imputed'
      
      source_df <- rbind(rename(filter(preprocess_out$sources.s, year == yrs[j]),
                                cc.long = cc.long.s,
                                cc.lat = cc.lat.s),
                         rename(filter(preprocess_out$sources.i1, year == yrs[j]),
                                cc.long = cc.long.i1,
                                cc.lat = cc.lat.i1))
    }
    
    sub_data <- sub_data %>%
      mutate(month = month(symptom.date),
             day = day(symptom.date),
             time = yday(symptom.date),
             iresponse.src1 = log(1 + rho.src1), # isotropic response for first source
             iresponse.src2 = log(1 + rho.src2), # isotropic response for second source
             aresponse.src1 = log(1 + rho.src1/g.theta.src1), # anisotropic response for first source
             aresponse.src2 = log(1 + rho.src2/g.theta.src2), # anisotropic response for second source
             sinb1.src1 = sin(as.numeric(theta.src1)), # sin basis 1 for first source
             sinb2.src1 = sin(as.numeric(theta.src1) + pi/4), # sin basis 2 for first source
             sinb1.src2 = sin(as.numeric(theta.src2)), # sin basis 1 for second source
             sinb2.src2 = sin(as.numeric(theta.src2) + pi/4), # sin basis 2 for second source
             z.reg.src1 = rho.src1/g.theta.src1, # first source regressor to initialize weights
             z.reg.src2 = rho.src2/g.theta.src2, # second source regressor to initialize weights
             z.fit = rho.src1 < rho.src2) # is first source closer than second source?
    
    sub_data <- sub_data %>%
      group_by(ID) %>%
      sample_n(1)
    
    out <- tbl_fn(sub_data, 
                  wind_src1, wind_src2, 
                  source_df, src1_id, src2_id,
                  src1_type, src2_type)
    
    return(out)}))
  
  return(out)}))

write.csv(file='draft_rmse_tbl.csv', rbind(sw_tbl, n_tbl))



out <- rbind(sw_tbl, n_tbl) %>%
  arrange(year) %>%
  rename(ATS.time = ats1.t.rmse,
         ITS.time = its1.t.rmse,
         AOS.time = aos.t.rmse,
         IOS.time = ios.t.rmse,
         ATS.distance = ats.rho.rmse,
         ITS.distance = its.rho.rmse,
         AOS.distance = aos.rho.rmse,
         IOS.distance = ios.rho.rmse) %>%
  mutate(date1 = as.character(src1.date),
         date2 = as.character(src2.date)) %>%
  select(year, n, src1.county, src1.state, date1,
         src2.county, src2.state, date2, src2.type,
         IOS.time, IOS.distance,
         AOS.time, AOS.distance, 
         ITS.time, ITS.distance, 
         ATS.time, ATS.distance)

write.csv(file='rmse_tbl.csv', out)
library(xtable)
xtable(out)

out %>% 
  select(year, n, src1.county, src1.state, date1,
         IOS.time, IOS.distance,
         AOS.time, AOS.distance) %>%
  distinct() %>%
  xtable()

out %>% 
  select(year, n, src2.county, src2.state, date2, src2.type,
         ITS.time, ITS.distance,
         ATS.time, ATS.distance) %>%
  xtable()
