library(tidyverse)
library(here)
library(sf)
library(wesanderson)

# Experimental Data
#Load in the Stripe Rust Data
stripe_dat <- readRDS(here("DataProcessed/experimental/stripe_clean.rds"))
inocs <- readRDS(here("DataProcessed/experimental/inoc_sp.rds"))

intensity_cols <- wes_palette("Zissou1", 7, "continuous")

stripe_map <- stripe_dat %>% 
  select(block, inoculum_total, visit, date, north, east, intensity) %>% 
  mutate(visit = parse_number(visit)) %>% 
  rename(treat = inoculum_total) %>% 
  st_as_sf(coords = c("east", "north")) |> 
  filter(treat == 1, block == "A")

stripe_map |> 
  ggplot()+
  geom_sf(aes(color = intensity), size = 2)+
  geom_sf(data = inocs |> filter(inoculum_total == 1, block == "A"), aes(fill = "Inoculation Point"), shape = 23, size = 3)+
  labs(x = "East (m)", y = "North (m)", color = "Intensity (%)", fill = "") +
  facet_wrap(~date, nrow = 1)+
  theme_bw()+
    scale_color_gradientn(colors = intensity_cols, transform = "sqrt")+
    theme(legend.position = "bottom")

ggsave("DiseaseIntensity_A.png", path = here("Reports/CAS_poster/"), width = 10, height = 3, units = "in")

#Grid Plots
backward <- readRDS(here("DataProcessed/results/backward_model/backward_fits.rds"))
clusters <- readRDS(here("DataProcessed/experimental/clusters.rds"))

grid <- bind_rows(clusters$stripe_4$grid,
                  clusters$stripe_8h$grid,
                  clusters$stripe_8v$grid,
                  clusters$stripe_16$grid)


grid %>% 
  ggplot()+
  geom_sf(fill = "transparent")+
  geom_sf_text(aes(label = grid_id))+
  facet_wrap(~grid_type, nrow = 1)+
  labs(x = "E", y = "N")+
  theme_classic()

# Backward Result Plots
nodes <- 