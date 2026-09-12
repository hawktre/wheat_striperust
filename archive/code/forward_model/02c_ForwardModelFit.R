## ---------------------------
##
## Script name: 04_ModelFit_Hurdle_Logistic.R
##
## Purpose of script: Fit model with logistic autoinfection term
##                    Updated to process and save all grid search profiles.
##
## Author: Trent VanHawkins
##
## Date Created: 2025-05-16
##
## ---------------------------
##
## Notes:
##
##
## ---------------------------

## view outputs in non-scientific notation
options(scipen = 6, digits = 4)

## load up the packages we will need:
library(here)
library(tidyverse)
library(data.table)
library(MASS)

# Read in the data --------------------------------------------------------
mod_dat <- readRDS(here("DataProcessed/experimental/mod_dat_arrays.rds"))
source(here("Code/R/forward_model/02b_ForwardModelFun.R"))


# Set up indices ----------------------------------------------------------
blocks <- dimnames(mod_dat$intensity)[["block"]]
treats <- dimnames(mod_dat$intensity)[["treat"]]
visits <- dimnames(mod_dat$intensity)[["visit"]]
kappa_try <- exp(seq(log(0.5), log(2), length.out = 5))

combos <- expand.grid(
  block = blocks,
  treat = treats,
  visit = visits[2:length(visits)],
  stringsAsFactors = FALSE
)

# Fit the model -----------------------------------------------------------

start <- Sys.time()
# forward_fit now returns ALL valid kappa combinations per profile row
free_fits <- pmap(
  combos,
  ~ forward_fit(..1, ..2, ..3, mod_dat, kappa_try)
) %>%
  rbindlist()
end <- Sys.time()
runtime <- difftime(end, start, units = "mins")
message("Runtime = ", round(runtime, 2), " minutes")




# Filter down to the best profile per combination for plotting ------------
message("Isolating optimal fits for down-stream diagnostics...")
the_best_fits <- free_fits |> 
  filter(converged == TRUE) |> 
  group_by(block, treat, visit) |> 
  # Rounding creates "ties" for practically identical likelihoods
  mutate(loglik_toll = round(neg_loglik, 0)) |> 
  # Slices by the best likelihood tier first, breaks ties with the best gradient
  slice_min(order_by = tibble(loglik_toll, grad_norm), n = 1) |> 
  ungroup() |> 
  dplyr::select(-loglik_toll) |> 
  as.data.table()# Clean up the temporary column

# Save the full execution results containing ALL grid fits
saveRDS(the_best_fits, here("DataProcessed/results/forward_model/forward_fits.rds"))

# Get Fitted values & residuals safely using explicitly named columns -----
the_best_fits$fitted <- pmap(
  list(the_best_fits$block, the_best_fits$treat, the_best_fits$visit, the_best_fits$theta, the_best_fits$alpha),
  ~ get_fitted(
    blk = ..1,
    trt = ..2,
    vst = ..3,
    par = ..4,
    alpha = ..5,
    mod_dat = mod_dat
  )
)

the_best_fits$resid <- pmap(
  list(the_best_fits$block, the_best_fits$treat, the_best_fits$visit, the_best_fits$theta, the_best_fits$alpha, the_best_fits$fitted),
  ~ compute_deviance_resid(
    blk = ..1,
    trt = ..2,
    vst = ..3,
    par = ..4,
    alpha = ..5,
    fitted = ..6,
    data = mod_dat
  )
)


# Faceted Diagnostic Plots -------------------------------------------------
message("Generating block-specific faceted diagnostic plots...")

# Unnest arrays into a tidy dataframe for plotting
plot_dat <- pmap_dfr(
  list(the_best_fits$block, the_best_fits$treat, the_best_fits$visit, the_best_fits$fitted, the_best_fits$resid),
  function(blk, trt, vst, fit_vals, res_vals) {
    obs_vals <- mod_dat$intensity[, blk, trt, vst]
    tibble(
      observed = obs_vals,
      fitted = fit_vals,
      residual = res_vals,
      block = blk,
      treat = trt,
      visit = factor(vst, levels = visits[-1]) # Ensure chronological order
    )
  }
)

# Loop through each block and save/print unique plots
for (b in blocks) {
  block_subset <- plot_dat %>% filter(block == b)
  
  # 1. Fitted vs Observed Plot Faceted by Visit & Treatment
  p_fit_obs <- ggplot(block_subset, aes(x = fitted, y = observed)) +
    geom_point(alpha = 0.5, color = "steelblue", size = 1) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
    facet_grid(visit ~ treat, labeller = label_both) +
    theme_bw() +
    labs(
      title = paste("Block", b, "- Fitted vs. Observed Intensity"),
      x = "Fitted Values",
      y = "Observed Intensity"
    )
  
  # 2. Residuals vs Fitted Plot Faceted by Visit & Treatment
  p_res_fit <- ggplot(block_subset, aes(x = fitted, y = residual)) +
    geom_point(alpha = 0.5, color = "darkorange", size = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    facet_grid(visit ~ treat, labeller = label_both) +
    theme_bw() +
    labs(
      title = paste("Block", b, "- Deviance Residuals vs. Fitted"),
      x = "Fitted Values",
      y = "Deviance Residuals"
    )
  
  # Print plots to the active device screen
  print(p_fit_obs)
  print(p_res_fit)
}


# Parameter Estimate Visualization -----------------------------------------
message("Generating parameter estimation timelines...")

# Unpack parameter data tracking matrices using optimal dataset subsets
param_dat <- the_best_fits[converged == TRUE, {
  pars <- theta[[1]]
  list(
    parameter = names(pars),
    value = as.numeric(pars)
  )
}, by = .(block, treat, visit)]

# Convert visit to a factor for chronological plotting order
param_dat[, visit := factor(visit, levels = visits[-1])]

# 2. Plot parameter timelines faceted by parameter type and treatment
p_params <- ggplot(param_dat, aes(x = visit, y = value, group = block, color = block)) +
  geom_line(alpha = 0.7, linewidth = 0.8) +
  geom_point(size = 1.5) +
  facet_grid(parameter ~ treat, scales = "free_y", labeller = label_both) +
  theme_bw() +
  scale_color_brewer(palette = "Set1") +
  labs(
    title = "Estimated Structural Parameters Across Visits",
    subtitle = "Faceted by Parameter Type and Experimental Treatment (Optimal Profile)",
    x = "Visit Tracking Number",
    y = "Estimated Natural Scale Value",
    color = "Experimental\nBlock"
  ) +
  theme(
    strip.text.y = element_text(angle = 0), # Keep parameter labels horizontal
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

# Display the plot
print(p_params)


