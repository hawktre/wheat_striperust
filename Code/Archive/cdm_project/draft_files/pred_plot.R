fit_out <- fit_pred(2009, training_prop = 0.64)
plot(fit_out$true_value, fit_out$pr_value)
abline(a = 0, b = 1)

theta_grid <- seq(0, 2*pi, length = 100)

dates <- ymd(c('2009/05/01', 
               '2009/06/01',
               '2009/07/01',
               '2009/08/01'))

pred_df <- expand.grid(theta = theta_grid, 
                       date = dates)


response <- sapply(1:nrow(pred_df), function(row){
  pred_res(lmod = fit_out$lmod,
         date = pred_df$date[row],
         theta = pred_df$theta[row])
})

# g.theta.pred <- sapply(1:nrow(pred_df),
#                        function(row){
#                          dvonmises(x = circular(pred_df$theta[row],
#                                                type = 'angles',
#                                                units = 'radians'), 
#                                    mu = vm_parameters[1, 2], 
#                                    kappa = as.numeric(vm_parameters[1, 3]))
#                        })

yr <- 2009
x <- filter(wind_data, year == yr)$direction
g.theta.pred <- density.circular(x = x, z = pred_df$theta, bw = 20)$y

# pred_df <- bind_cols(pred_df,
#           pred = response) %>%
#   mutate(g.theta = dvonmises(theta, 
#                              as.circular(vm_parameters[1, 2]), 
#                              as.numeric(vm_parameters[1, 3]))) %>%
#   mutate(rho = (exp(pred) - 1)*g.theta)

pred_df <- bind_cols(pred_df,
                     pred = response, 
                     g.theta = g.theta.pred) %>%
  mutate(rho = (exp(pred) - 1)*g.theta)

pred_latlong <- sapply(1:nrow(pred_df), function(row){
  destPoint(p = first_report_centroid[1, 1:2],
          b = pred_df$theta[row]*360/(2*pi), 
          d = 1000*pred_df$rho[row])
}) %>%
  t() 

colnames(pred_latlong) <- c('pred_lon', 'pred_lat')

pred_df <- bind_cols(pred_df, as.data.frame(pred_latlong))

data09 <- filter(DATA, year == 2009) %>%
  select(lat, long, report_date)

library(ggmap)
lonrange <- range(c(data09$long, pred_df$pred_lon))
latrange <- range(c(data09$lat, pred_df$pred_lat))

mapbox <- make_bbox(lonrange, latrange)
map <- get_map(location = mapbox,
               maptype = 'toner-background')

ggmap(map, darken = c(0.8, 'white')) +
  geom_path(data = pred_df,
            aes(x = pred_lon,
                y = pred_lat,
                group = factor(month(date)),
                color = factor(month(date)))) +
  geom_point(data = data09,
             aes(x = long,
                 y = lat,
                 color = factor(month(report_date)))) +
  scale_color_brewer(type = 'seq',
                     palette = 'Reds') +
  labs(x = 'Longitude',
       y = 'Latitude') +
  theme_bw() +
  guides(color = guide_legend('Month'))

            