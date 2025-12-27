## ---------------------------
##
## Script name: 
##
## Purpose of script:
##
## Author: Trent VanHawkins
##
## Date Created: 2025-08-05
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


# Read in the data and simulations ----------------------------------------
forward <- readRDS(here("DataProcessed/results/forward_model/forward_fits.rds"))
backward <- readRDS(here("DataProcessed/results/backward_model/backward_fits_shared.rds"))
sims <- readRDS(here("DataProcessed/results/simulation/batch_results/simulation_batch3.rds"))
mod_dat <- readRDS(here("DataProcessed/experimental/mod_dat_arrays.rds"))

# Forward model parameters ------------------------------------------------
forward_sims <- sims %>% 
  select(sim, block, treat, visit, iters, converged.forward, init_kappa, theta, pi) %>% 
  distinct()

not_converged_forward <- forward_sims %>% 
  filter(!converged.forward) %>% 
  group_by(block, treat, visit) %>% 
  summarise(N = n())
  

forward_sims_long <- forward_sims %>%
  filter(converged.forward == T) %>% 
  # Convert each named theta vector into a tibble (with name = param)
  mutate(theta = map(theta, ~ enframe(.x, name = "param", value = "value"))) %>%
  # Unnest into long format
  unnest(theta)

true_params <- forward %>% 
  mutate(theta = map(theta, ~ enframe(.x, name = "param", value = "value"))) %>%
  # Unnest into long format
  unnest(theta) %>% 
  select(block, treat, visit, param, value) 
  
forward_sims_and_truth <- left_join(forward_sims_long, true_params, by = c("block", "treat", "visit", "param"), suffix = c(".sim", ".true"))

bias_df <- forward_sims_and_truth %>% 
  select(sim, block, treat, visit, param, value.sim, value.true) %>% 
  mutate(bias = value.sim - value.true)

gamma_bias <- bias_df %>% 
  filter(param == "gamma", abs(bias) > 100)

bias_summary <- bias_df %>% 
  group_by(block, treat, visit, param) %>% 
  summarise(mean_bias = mean(bias),
            med_bias = median(bias), 
            sd_bias = sd(bias)) %>% 
  arrange(mean_bias)

bias_gamma <- bias_summary %>% filter(param == "gamma")
bias_ngamma <- bias_summary %>% filter(param != "gamma")
# Make a boxplot of parameter distributions -------------------------------
plot_forward_sims <- function(blk){
  forward_sims_and_truth %>%
  mutate(treat_lab = paste0(treat," Source(s)"),
         visit_lab = paste0("Visit ", visit)) %>% 
  filter(block == blk) %>%
  ggplot(aes(x = param, y = value.sim)) +
  geom_boxplot(outlier.alpha = 0.3) +
  geom_point(aes(y = value.true, color = param), 
             shape = 18, size = 3, position = position_dodge(width = 0.75)) +
  facet_grid(visit_lab ~ treat_lab, scales = "free_y") +
  labs(y = "Estimated value", x = "Parameter", color = "True Value",
       title = paste0("Forward Model Simulated Parameter Estimates (Block ", blk, ")")) +
  theme_bw()}

for (letter in LETTERS[1:4]) {
  cur.plt <- plot_forward_sims(letter)
  ggsave(plot = cur.plt, filename = paste0("forward_sim_params_block", letter,".png"), 
         width = 10, height = 8, path = here("DataProcessed/results/simulation/figures/"))
}

forward_sims_and_truth %>% 
  filter(param == "gamma", treat == 1) %>%
  ggplot(aes(x = value.sim))+
  geom_histogram()+
  geom_vline(aes(xintercept = value.true, color = "True Value"), linetype = "dashed")+
  facet_grid(block~visit, scales = "free_x") +
  labs(y = "Count", x = "Estimated Gamma", color = "", 
       title = "Distribution of Estimated Gamma (Srouce(s) = 1)")+
  theme_bw()
  
  

# Make boxplots of dist_acc -----------------------------------------------
sims %>%
  filter(treat == 1) %>%
  mutate(config = factor(config, levels = c("4", "8h", "8v", "16"))) %>% 
  ggplot(aes(x = as.character(visit), y = dist_acc)) +
  geom_jitter(width = 0.2, alpha = 0.2, size = 0.5, color = "steelblue") +
  stat_summary(fun = median, geom = "point", aes(color = "median"), size = 2) +
  facet_grid(config ~ block) +
  labs(x = "Visit", y = "Distance-Weighted Accuracy", color = "") +
  theme_bw()

ggsave(filename = "backward_sims_acc.png", 
       path = here("DataProcessed/results/simulation/figures/"),
       width = 10, height = 10, )

check <- vector("numeric", length = 1000)
for (i in 1:1000) {
  check[i]<-sum(disease_sim(forward$theta[[1]], forward$fitted[[1]]) < 1e-8)
}
sum(check > 0)
lgamma(1e-17)
