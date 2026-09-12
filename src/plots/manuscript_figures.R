library(here)
library(tidyverse)
library(data.table)
library(sf)


# Read in the data -------------------------------------------------------

## Experiment Data
stripe <- readRDS(here("data/processed/experimental/stripe_clean.rds")) %>%
  st_as_sf(coords = c("east", "north"))

## Inoculations
inocs <- readRDS(here("data/processed/experimental/inoc_sp.rds"))

## Grids
design <- readRDS(here("data/processed/experimental/grids_sp.rds"))

## Model data
mod_dat <- readRDS(here("data/processed/experimental/mod_dat_arrays.rds"))
# Read in results ---------------------------------------------------------
## Forward Model
forward <- readRDS(here(
  "output/source_detection/forward_fits.rds"
))

## Backward Model
backward <- readRDS(here(
  "archive/output/backward_model/archived/backward_fits_full.rds"
))


### Separate into single source and multi-source
single <- backward |> filter(n_src.x == 1) |> select(-n_src.y) |> rename(n_src = n_src.x)
multi <- backward |> filter(n_src.x > 1) |> select(-n_src.y) |> rename(n_src = n_src.x)


# Raw Data Plot ----------------------------------------------------------
for(blk in LETTERS[1:4]){
  inoc_tmp <- inocs |> filter (block == blk)
  
  stripe |> 
    filter(block == blk) |> 
    ggplot()+
    geom_sf(aes(color = intensity), size = 2)+
    geom_sf(data = inoc_tmp , shape = 23, aes(fill = "Inoculation Point"), size = 3)+
    labs(x = "East (m)", y = "North (m)", color = "Intensity", fill = "") +
    # Use the labeller argument to append the text to the 'treat' variable
    facet_grid(
      treat ~ date, 
      labeller = labeller(treat = function(x) paste(x, "inoculation(s)"))
    )+
    theme_bw()+
    scale_color_viridis_c()+
    theme(legend.position = "bottom")

  ggsave(paste0("DiseaseIntensity_",blk,".png"), path = here("manuscript/Figures/raw_intensity/"), width = 10, height = 6, units = "in")
}


# Design plot ------------------------------------------------------------

grids <- rbind(
  design$`4`$grid,
  design$`8h`$grid,
  design$`8v`$grid,
  design$`16`$grid,
  design$`64`$grid
) |> 
  mutate(
    name = factor(
      name,
      levels = c(
        "K = 4",
        "K = 8 (Horizontal)",
        "K = 8 (Vertical)",
        "K = 16",
        "K = 64"
      )
    )
  )

for(blk in LETTERS[1:4]){
  inoc_tmp <- inocs |> filter (block == blk)
  
  grids  |> 
  ggplot() +
  geom_sf(fill = "transparent") +
  geom_sf_text(aes(label = grid_id), size = 2.5) +
  geom_sf( data = inoc_tmp, shape = 23, aes(fill = "Inoculation Point"), size = 3) +
  labs(x = "East (m)", y = "North (m)", fill = "") +
  facet_grid(
    treat ~ name,
    labeller = labeller(
      treat = function(x) paste0(x, " inoculation(s)")
    )
  ) +
  theme_classic() +
  theme(legend.position = "bottom")

  ggsave(paste0("BackwardDesign_",blk,".png"), path = here("manuscript/Figures/backward_design/"), width = 8, height = 8, units = "in")
}


# Distance Error Plots ---------------------------------------------------
## Single Source
backward |> 
  mutate(visit = as.numeric(visit),
treat = as.numeric(treat),
config = factor(config, levels = c("4", "8h", "8v", "16", "64"), 
labels = c("K = 4", "K = 8 Horizontal", "K = 8 Vertical", "K = 16", "K = 64"))) |> 
  filter(treat == 1) |> 
  ggplot(aes(x = visit, y = mean_error))+
  geom_point(aes(color = config)) + 
  geom_smooth(aes(color = config), se = F, method = "lm")+
  geom_hline(yintercept = 0, linetype = "dotted")+
  facet_grid(config~block, labeller = labeller(
      block = function(x) paste0("Block ", x)
    ) ) +
  scale_color_viridis_d()+
  theme_bw() +
  theme(legend.position = "bottom")+
  labs(color = "Components",
x = "Week", y = "Error (m)")

## Single Source loess across blocks
backward |> 
  mutate(visit = as.numeric(visit),
treat = as.numeric(treat),
config = factor(config, levels = c("4", "8h", "8v", "16", "64"), 
labels = c("K = 4", "K = 8 Horizontal", "K = 8 Vertical", "K = 16", "K = 64"))) |> 
  filter(treat == 1) |> 
  ggplot(aes(x = visit, y = mean_error))+
  geom_point(aes(color = block), alpha = 0.75) + 
  geom_smooth(color = "white", se = F)+
  facet_wrap(~config) +
  scale_color_viridis_d()+
  theme_dark() +
  theme(legend.position = "bottom")+
  labs(color = "Block",
x = "Week", y = "Error (m)")

## Single source averaged
backward |> 
  group_by(treat, visit, config) |>
  summarise(
    mean_error = mean(mean_error, na.rm = TRUE),
    visit = first(as.numeric(visit)),
    treat = first(as.numeric(treat)),
    config = first(factor(config, levels = c("4", "8h", "8v", "16", "64"), 
                           labels = c("K = 4", "K = 8 Horizontal", "K = 8 Vertical", "K = 16", "K = 64")))
  )|> 
  ungroup() |> 
  filter(treat == 1) |> 
  ggplot(aes(x = visit, y = mean_error))+
  geom_point(aes(color = config)) + 
  geom_smooth(aes(color = config), se = F)+
  geom_hline(yintercept = 0, linetype = "dotted")+
  facet_wrap(~config ) +
  scale_color_viridis_d()+
  theme_bw() +
  theme(legend.position = "bottom")+
  labs(color = "Components",
x = "Week", y = "Error (m)")
