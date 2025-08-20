## ---------------------------
##
## Script name: 04a_SimFunc.R
##
## Purpose of script: Compile all necessary functions for simulation.
##
## Author: Trent VanHawkins
##
## Date Created: 2025-07-22
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
  
# Function to simulate the data -------------------------------------------
disease_sim <- function(pars, mu, pi) {

  phi <- pars[["phi"]]
  
  # Convert to shape parameters
  alpha <- mu * phi
  beta <- phi * (1 - mu)
  
  # Minimum perceptible value
  min_detectable <- 0.0001
  
  # Draw from beta and clip
  sim_dat <- map2_dbl(alpha, beta, ~max(min_detectable, rbeta(1, shape1 = .x, shape2 = .y)))
  
  # Apply zero inflation
  zeros <- rbinom(length(sim_dat), size = 1, prob = 1 - pi)
  
  return(sim_dat * zeros)
}
  

# Wrapper Function for the simulation -------------------------------------
single_sim <- function(sim_id, dat, forward_mod, kappa_try, output_dir = here("DataProcessed/results/simulation")) {

  base_seed <- 404
  set.seed(base_seed + sim_id)
  tryCatch({
    # Set up indices
    blocks <- dimnames(dat$intensity)[["block"]]
    treats <- dimnames(dat$intensity)[["treat"]]
    visits <- dimnames(dat$intensity)[["visit"]]
    configs <- dimnames(dat$groups)[["config"]]
    
    # 1. Simulate new intensity values
    intensity_sim <- dat$intensity
    for (blk in blocks) {
      for (trt in treats) {
        for (vst in visits[-1]) {
          fit <- forward_mod[block == blk & treat == as.numeric(trt) & visit == as.numeric(vst)]
          pars <- fit[["theta"]][[1]]
          pi <- fit[["pi"]]
          fitted <- fit[["fitted"]][[1]]
          intensity_sim[, blk, trt, vst] <- disease_sim(pars, fitted, pi)
        }
      }
    }
    dat$intensity <- intensity_sim
    
    # 3. Fit Forward Model
    combos_forward <- expand.grid(block = blocks, treat = treats, visit = visits[-1], stringsAsFactors = FALSE)
    forward <- pmap(combos_forward, ~forward_fit(..1, ..2, ..3, dat, dat$dist, kappa_try)) %>% rbindlist()
    
    # 2. Fit backward model
    combos_backward <- expand.grid(config = configs, blk = blocks, trt = treats, vst = visits[-1], stringsAsFactors = FALSE)
    backward <- pmap(combos_backward, ~backward_fit(..1, ..2, ..3, ..4, dat, forward)) %>% rbindlist()
    
    # 4. Source prediction for treat == 1
    backward_t1 <- backward[treat == 1]
    sources_predicted <- pmap(backward_t1, ~source_pred(config = ..1,
                                                        blk = ..2,
                                                        trt = ..3,
                                                        vst = ..4,
                                                        p_mat = ..9,
                                                        mod_dat = dat)) %>% 
      rbindlist()
    
    results_merge <- left_join(backward, forward, by = c("block", "treat", "visit"), suffix = c(".backward", ".forward")) %>% 
      dplyr::select(-c(final_theta, p_mat)) %>% 
      left_join(sources_predicted, by = c("config", "block", "treat", "visit")) %>% 
      mutate(sim = sim_id) %>% 
      dplyr::select(sim, everything())
    # 5. Return results
    return(results_merge)
    
  }, error = function(e) {
    err_file <- file.path(output_dir, paste0("sim_error_", sim_id, ".txt"))
    writeLines(c(
      "Simulation failed:",
      conditionMessage(e),
      capture.output(traceback(max.lines = 10))
    ), con = err_file)
    return(NULL)
  })
}


