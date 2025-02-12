library(tidyverse)
load('3source_data.RData')

model_formula <- formula('response ~ time + sinb1 + sinb2')

pr_z0 <- matrix(c(4/10, 3/10, 3/10),   
                nrow = nrow(source1_data),
                ncol = 3, 
                byrow = T)

# if rho is NA then the source hasn't appeared by the time of that observation
const2_idx <- which(is.na(source2_data$rho))
update2_idx <- which(!is.na(source2_data$rho))
pr_z0[const2_idx, 2] <- 0

const3_idx <- which(is.na(source3_data$rho))
update3_idx <- which(!is.na(source3_data$rho))
pr_z0[const3_idx, 3] <- 0

pr_z0 <- pr_z0/apply(pr_z0, 1, sum)
const1_idx <- which(pr_z0[, 1] == 1)
update1_idx <- which(pr_z0[, 1] != 1)

pr_z <- pr_z0
pr_z_diff <- 1
iter <- 1

while(pr_z_diff > 0.1){
  pr_z_old <- pr_z
  
  # step 1: fit model 1
  fit_1 <- lm(model_formula, 
              data = source1_data,
              weights = pr_z[, 1])
  
  
  # step 2: fit model 2 (N source)
  fit_2 <- lm(model_formula, 
              data = source2_data,
              weights = pr_z[, 2])
  
  # step 2: fit model 3 (N source)
  fit_3 <- lm(model_formula, 
              data = source3_data,
              weights = pr_z[, 3])
  
  
  mean1 <- fit_1$residuals/summary(fit_1)$sigma
  mean2 <- fit_2$residuals/summary(fit_2)$sigma
  mean3 <- fit_3$residuals/summary(fit_3)$sigma
  
  pr_z[, 1] <- dnorm(mean1)*pr_z[, 1]
  pr_z[update2_idx, 2] <- dnorm(mean2)*pr_z[update2_idx, 2]
  pr_z[update3_idx, 3] <- dnorm(mean3)*pr_z[update3_idx, 3]
  
  
  pr_z[const1_idx, 1] <- 1
  pr_z[const2_idx, 2] <- 0
  pr_z[const3_idx, 3] <- 0
  
  pr_z <- pr_z/apply(pr_z, 1, sum)
  pr_z_diff <- mean(abs(pr_z - pr_z_old))
  
  iter <- iter + 1
  print(iter)
}

summary(fit_1)
summary(fit_2)
summary(fit_3)
apply(round(pr_z), 2, sum)

# calculate time predictions
betahat1 <- coef(fit_1)
pred1 <- source1_data %>%
  ungroup() %>%
  mutate(tpred1 = (response - betahat1[1] - betahat1[3]*sinb1 - betahat1[4]*sinb2)/betahat1[2]) %>%
  select(ID, tpred1)
  
betahat2 <- coef(fit_2)
pred2 <- source2_data %>%
  ungroup() %>%
  mutate(tpred2 = (response - betahat2[1] - betahat2[3]*sinb1 - betahat2[4]*sinb2)/betahat2[2]) %>%
  select(ID, tpred2)

betahat3 <- coef(fit_3)
pred3 <- source3_data %>%
  ungroup() %>%
  mutate(tpred3 = (response - betahat3[1] - betahat3[3]*sinb1 - betahat3[4]*sinb2)/betahat3[2]) %>%
  select(ID, tpred3)

preds <- merge(pred1, pred2, by = 'ID') %>%
  merge(pred3, by = 'ID') %>%
  bind_cols(data.frame(p.hat = pr_z)) 

preds <- preds %>%
  select(ID) %>%
  bind_cols(tpred = apply(as.matrix(preds[, 2:4])*round(pr_z), 1, function(x){sum(x, na.rm = T)})) %>%
  merge(select(source1_data, ID, time), by = 'ID') %>%
  mutate(tresid = time - tpred)

sqrt(mean(preds$tresid^2))

plot_data1 <- sub_data %>%
  select(ID, lat, long, symptom.date) %>%
  merge(preds, by = 'ID') %>%
  mutate(which.pred = apply(pr_z, 1, which.max))


theta_grid <- seq(0, 2*pi, length = 100)
week_grid <- seq(week(ymd(paste(yrs[j], '/03/01', sep = ''))), 
                 week(ymd(paste(yrs[j], '/08/01', sep = ''))),
                 by = 3)
date_grid <- rep(ymd(paste(yrs[j], '/04/15', sep = '')), length(week_grid))
week(date_grid) <- week_grid
pred_df <- expand.grid(theta = theta_grid, 
                       date = date_grid) %>%
  mutate(time = yday(date),
         sinb1.src1 = sin(theta),
         sinb2.src1 = sin(theta + pi/4),
         sinb1 = sin(theta),
         sinb2 = sin(theta + pi/4))

pred_df <- pred_df %>%
  mutate(pred.src1 = predict(fit_1, newdata = pred_df),
         pred.src2 = predict(fit_2, newdata = pred_df),
         pred.src3 = predict(fit_3, newdata = pred_df),
         g.theta.src1 = density.circular(x = filter(wind_src1, year == yrs[j])$direction,
                                         z = pred_df$theta, 
                                         bw = 50)$y,
         g.theta.src2 = density.circular(x = filter(wind_src2, year == yrs[j])$direction,
                                         z = pred_df$theta, 
                                         bw = 50)$y,
         g.theta.src3 = density.circular(x = filter(wind_src3, year == yrs[j])$direction,
                                         z = pred_df$theta, 
                                         bw = 50)$y) %>%
  mutate(rho.src1 = (exp(pred.src1) - 1)*g.theta.src1,
         rho.src2 = (exp(pred.src2) - 1)*g.theta.src2,
         rho.src3 = (exp(pred.src3) - 1)*g.theta.src3) %>%
  select(theta, date, time, 
         rho.src1, rho.src2, rho.src3)

pred_df <- pred_df %>%
  cbind(src1 = destPoint(p = source_df[1, 2:3],
                              b = pred_df$theta*360/(2*pi),
                              d = 1000*pred_df$rho.src1),
        src2 = destPoint(p = source_df[2, 2:3],
                              b = pred_df$theta*360/(2*pi),
                              d = 1000*pred_df$rho.src2),
        src3 = destPoint(p = source_df[3, 2:3],
                              b = pred_df$theta*360/(2*pi),
                              d = 1000*pred_df$rho.src3))

library(ggmap)
load('cdm_mapbox.RData')
map <- get_map(location = mapbox,
               maptype = 'toner-background')

plot_data1 <- plot_data1 %>%
  mutate(which.pred = apply(pr_z, 1, which.max))

ggmap(map, darken = c(0.45, 'grey')) +
  ## plot errors
  geom_point(data = plot_data1, 
             aes(x = long, 
                 y = lat, 
                 shape = factor(which.pred,
                                labels = c('Florida', 
                                           'Alternate SW',
                                           'Alternate N')),
                 size = abs(tresid),
                 color = month(symptom.date)))  +
  geom_path(data = pred_df,
            aes(x = src1.lon,
                y = src1.lat,
                color = month(date),
                group = week(date)),
            linetype = 1) +
  geom_path(data = pred_df,
            aes(x = src2.lon,
                y = src2.lat,
                color = month(date),
                group = week(date)),
            linetype = 2) +
  geom_path(data = pred_df,
            aes(x = src3.lon,
                y = src3.lat,
                color = month(date),
                group = week(date)),
            linetype = 3) +
  geom_point(data = source_df,
             aes(x = cc.long, y = cc.lat),
             color = 'red',
             shape = 4,
             size = 3) +
  ## adjust legends
  guides(shape = guide_legend('Prediction Source',
                              override.aes = list(size = 3)),
         linetype = guide_legend('Prediction Source'),
         color = guide_legend('Month'),
         size = guide_legend('Prediction Error')) +
  scale_color_gradientn(colors = c('yellow', 'orange', 'red', 'black'),
                        breaks = 1:7,
                        labels = month.name[1:7]) +
  labs(x = 'Longitude',
       y = 'Latitude') +
  theme(text = element_text(size = 10)) +
  scale_size_binned(breaks = c(7, 14, 21, 28),
                    range = c(1, 5),
                    labels = c('0-1 weeks',
                               '1-2 weeks',
                               '2-3 weeks',
                               '>3 weeks'),
                    guide = 'legend')

ggsave(filename = 'sfig1.tiff',
       width = 7.5, height = 5, units = 'in',
       dpi = 300)
