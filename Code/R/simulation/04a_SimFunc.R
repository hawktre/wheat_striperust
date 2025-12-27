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
disease_sim <- function(pars, mu, alpha) {

  phi <- exp(pars[["phi"]])
  
  # Convert to shape parameters
  a <- mu * phi
  b <- phi * (1 - mu)
  
  # Minimum and maximum perceptible values
min_detectable <- 0.0001
max_detectable <- 1 - 0.0001  # or 0.9999

# Draw from beta and clip both bounds
sim_dat <- map2_dbl(a, b, ~{
  val <- rbeta(1, shape1 = .x, shape2 = .y)
  pmin(pmax(val, min_detectable), max_detectable)
})
  
  # Apply zero inflation
  zeros <- rbinom(length(sim_dat), size = 1, prob = 1 - alpha)
  
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
          alpha <- fit[["alpha"]]
          fitted <- fit[["fitted"]][[1]]
          intensity_sim[, blk, trt, vst] <- disease_sim(pars, fitted, alpha)
        }
      }
    }
    dat$intensity <- intensity_sim
    
    
    # 3. Fit Forward Model
    combos_forward <- expand.grid(block = blocks, treat = treats, visit = visits[-1], stringsAsFactors = FALSE)
    forward <- pmap(combos_forward, ~forward_fit(..1, ..2, ..3, dat, kappa_try)) |> rbindlist()
    
    
    # 2. Fit backward model
    combos_backward <- expand.grid(config = configs, block = blocks, treat = treats, visit = visits[-1], stringsAsFactors = FALSE) |> 
      filter(!(config == "64" & treat == "4")) |> 
      left_join(forward |> select(block, treat, visit, theta), by = c("block", "treat", "visit"))
    backward <- lapply(seq_len(nrow(combos_backward)), function(row) {
      combo <- combos_backward[row,]
      backward_fit(config = combo$config,
        blk = combo$block,
        trt = combo$treat,
        vst = combo$visit,
        inits = combo$theta[[1]],
        mod_dat = mod_dat,
        tol = 1e-4,
        max_iter = 200)}) |> rbindlist()

    # 4. Source prediction
    sources_predicted <- backward|> 
  select(config, block, treat, visit, n_src, p_mat) |> 
  pmap(~source_pred(config = ..1,
                    blk = ..2,
                    trt = ..3,
                    vst = ..4,
                    n_src = ..5,
                    p_mat = ..6,
                    mod_dat = mod_dat)) |> 
  rbindlist()
    
    results_merge <- left_join(backward, forward |> mutate(visit = as.numeric(visit)), by = c("block", "treat", "visit"), suffix = c(".backward", ".forward")) %>% 
      left_join(sources_predicted, by = c("config", "block", "treat", "visit")) |> 
      mutate(sim = sim_id) |> 
      dplyr::select(sim, everything())
    # 5. Return results
    return(results_merge)
    
  }, error = function(e) {
    # rich, reproducible log
    err_file <- file.path(output_dir, sprintf("sim_error_%s.txt", sim_id))
    cat(
      "Simulation failed:\n",
      "sim_id: ", sim_id, "\n",
      "message: ", conditionMessage(e), "\n",
      "class: ", paste(class(e), collapse = ", "), "\n",
      "call: ", deparse(conditionCall(e)), "\n",
      file = err_file, append = FALSE, sep = ""
    )
    
    # capture the call stack for later interactive debugging
    dump_name <- file.path(output_dir, sprintf("sim_dump_%s", sim_id))
    dump.frames(dump_name, to.file = TRUE)
    
    cat("\nDumped frames to: ", dump_name, ".rda\n", file = err_file, append = TRUE)
    
    # if you’re in an interactive session, stop so you SEE the error
    stop(e)
  })
}


