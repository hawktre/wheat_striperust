library(circular)
library(tidyverse)

mu <- circular(pi/2)
nu <- circular(pi)
delta <- 4
kappa <- 2

f_fn <- function(x){
  dvonmises(circular(x), mu, delta)
}

g_fn <- function(x, g0){
  dvonmises(circular(x), nu, kappa)*g0
}


curve(f_fn, from = 0, to = 2*pi)
curve(g_fn(x, 100), from = 0, to = 2*pi)

K_fn <- function(rho, phi, g_0, b){
  f_fn(phi)*(b - 1)*(b - 2)/(g_fn(phi, g_0)^2)*((1 + rho/g_fn(phi, g_0))^(-b))
}

rho_vals <- seq(1, 10, length = 50)
phi_vals <- seq(0, 2*pi, length = 50)


z_df <- mutate(expand.grid(rho_vals, phi_vals),
            z = K_fn(Var1, Var2, 10, 4))


ggplot(z_df, aes(x = Var1, y = Var2, z = scale(log(z)))) + 
  geom_contour(binwidth = 0.25, aes(color = stat(level))) + 
  scale_y_continuous(breaks = seq(0, 2*pi, length = 5), 
                     limits = c(0, 2*pi),
                     labels = expression(0, pi/2, pi, 3*pi/2, 2*pi)
                     ) +
  xlim(c(1, 6)) + 
  scale_color_continuous(limits = c(-3, 3)) +
  labs(x = '', y = '') +
  guides(color = guide_colorbar('log(K(x, y))')) +
  coord_polar(theta = "y")

p <- 2
theta <- rvonmises(2, circular(mu), delta) %>% as.numeric()
a <- matrix(rexp(3), nrow = p, ncol = 3, byrow = T)
b <- rexp(3)
c <- rexp(3) 
gamma <- a %*% (cos(theta + c) + sin(theta + b) + 2)

curve(a %*% (cos(x + c) + sin(x + b) + 2), 0, 2*pi)


