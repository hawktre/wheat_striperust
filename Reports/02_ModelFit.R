## ---------------------------
##
## Script name: 02_ModelFit.R
##
## Purpose of script: Model fitting for forward model
##
## Author: Trent VanHawkins
##
## Date Created: 2025-04-22
##
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------

## ---------------------------

## view outputs in non-scientific notation

options(scipen = 6, digits = 4) 

## ---------------------------

## load up the packages we will need:  (uncomment as required)

library(tidyverse)
library(here)
source(here("Code/01a_DataFormat_Fun.R"))
source(here("Code/01b_GradDescent_fun.R"))


# Read in the data --------------------------------------------------------
stripe <- readRDS(here("DataProcessed/experimental/stripe_clean.rds"))
wind <- readRDS(here("DataProcessed/wind/wind_clean.rds"))


# Create lagged data frame ------------------------------------------------
stripe.lag <- stripe %>% 
  arrange(plotID_new, plant_num, visit) %>% 
  group_by(plant_id) %>%
  mutate(intensity = intensity/100,
         intensity_prev = lag(intensity),
         date_prev = lag(date)) %>% 
  ungroup() %>% 
  filter(visit != "visit1")


# Create data object for modeling -----------------------------------------
## Create a df of survey periods
survey_periods <- stripe.lag %>% 
  select(plotID_new, visit, date, date_prev) %>% 
  distinct()

mod_dat <- list()

for (i in 1:nrow(survey_periods)){

  #subset the data to the current plot-visit combo
  cur.plot <- survey_periods$plotID_new[i]
  cur.visit <- survey_periods$visit[i]
  
  plot_dat <- stripe.lag %>% 
    filter(plotID_new == cur.plot,
           visit == cur.visit)
  
  #Subset y_prev
  y_prev <-  plot_dat$intensity_prev
  
  #Subset y_cur
  y_cur <- plot_dat$intensity
  
  #distance and direction matrices
  dist_dir <- get_dist_dir(cur.plot, plot_dat)
  
  #wind matrix
  wind_mat <- get_wind_mat(first_day = survey_periods$date_prev[i],
                           last_day = survey_periods$date[i],
                           dir.mat = dist_dir$dir,
                           wind = wind)
  #Define values of kappa to try
  kappa_try <- seq(1,3,0.5)
  
  #Initialize theta
  theta_init <- initialize_theta(y_cur = y_cur, y_prev = y_prev, wind_mat, dist_dir$dist, 
                                 d_0 = 0.01, kappa_try)
  
  mod_dat[[cur.plot]][[cur.visit]] <- list("y_prev" = y_prev,
                                       "y_cur" = y_cur,
                                       "dist_mat" = dist_dir$dist,
                                       "dir_mat" = dist_dir$dir,
                                       "wind_mat" = wind_mat,
                                       "theta_init" = theta_init)
}

# Fit the model ------------------------------------------------------

fit_beta_model(y_current = mod_dat$A1$visit5$y_cur,
               y_prev = mod_dat$A1$visit5$y_prev,
               wind_matrix = mod_dat$A1$visit5$wind_mat,
               dist_matrix = mod_dat$A1$visit5$dist_mat,
               param_init = mod_dat$A1$visit5$theta_init[[1]])
