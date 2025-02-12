library(CircStats)
library(tidyverse)
set.seed(031221)

b_01 <- 4.85; b_11 <- 0.02; b_12 <- 0; b_13 <- 0   
b_02 <- 3.75; b_21 <- 0.03; b_22 <- 0; b_23 <- 0

mu1 <- pi/4; mu2 <- (5*pi)/4  #  location parameter μ is in [0, 2 π)
kappa1  <- 2; kappa2  <- 2# concentration parameter κ > 0

x1 <- 0; y1 <- 0; x2 <- 5000; y2 <- 5000
sources <- data.frame( "x" = c(x1,x2), "y" = c(y1,y2))
N <- 100  # 50, 100, 200
prop1 <- 0.7 * N  # proportion 1
prop2 <- N - prop1  # proportion 2

# time_grid1 <- runif(prop1, min = 100, max = 150)  # daily time grids
# time_grid2 <- runif(prop2, min = 100, max = 150)

time_grid1 <- seq(from = 100, to = 200, length = prop1) %>% sample() # daily time grids
time_grid2 <- seq(from = 100, to = 200, length = prop2) %>% sample()

# angle_grid1 <- runif(prop1, min = 0, max = 2*pi)
# angle_grid2 <- runif(prop2, min = 0, max = 2*pi)

angle_grid1 <- seq(from = 0, to = 2*pi, length = prop1) %>% sample()
angle_grid2 <- seq(from = 0, to = 2*pi, length = prop2) %>% sample()

g_phi1 <- function(phi){dvm(phi, mu1, kappa1)}
g_phi2 <- function(phi){dvm(phi, mu2, kappa2)}
# g_phi1 <- function(phi){return(1)}
# g_phi2 <- function(phi){return(1)}


y_generated1 <- rep(0, length(prop1))  # wrt source 1
for (i in 1:prop1) {
  y_generated1[i] <- b_01 + (b_11 * time_grid1[i]) + (b_12 * sin(angle_grid1[i])) + 
    (b_13 * (sin(angle_grid1[i] + pi/4))) + rnorm(1, 0, 0.5)  # μ = 0, var =  1 or 0.5
}

y_generated2 <- rep(0, length(prop2))  # wrt source 2
for (i in 1:prop2) {
  y_generated2[i] <- b_02 + (b_21 * time_grid2[i]) + (b_22 * sin(angle_grid2[i])) + 
    (b_23 * (sin(angle_grid2[i] + pi/4))) + rnorm(1, 0, 0.5)  # μ = 0, var =  1 or 0.5
}


y_generated1
y_generated2

sim_data1 <- tibble(response = y_generated1,
                    t = time_grid1,
                    rho = angle_grid1) %>%
  mutate(gphi = g_phi1(rho))%>%
  mutate(r = (exp(response) - 1)*gphi) %>%
  mutate(y = r*sin(rho) + y1,
         x = r*cos(rho) + x1) %>%
  mutate(source = 'source1') %>%
  select(t, x, y, source)

sim_data2 <- tibble(response = y_generated2,
                    t = time_grid2,
                    rho = angle_grid2) %>%
  mutate(gphi = g_phi2(rho)) %>%
  mutate(r = (exp(response) - 1)*gphi) %>%
  mutate(y = r*sin(rho) + y2,
         x = r*cos(rho) + x2) %>%
  mutate(source = 'source2') %>%
  select(t, x, y, source)



sim_data <- bind_rows(sim_data1, sim_data2) %>%
  mutate(time = cut(t, 3, labels = paste('time', 1:3)))

sim_data %>%
  ggplot(aes(x = x, y = y)) +
  geom_point(aes(color = time, shape = source)) +
  geom_point(data = sources, color = 'red') +
  scale_color_brewer(type = 'seq', palette = 2) +
  theme_bw()
