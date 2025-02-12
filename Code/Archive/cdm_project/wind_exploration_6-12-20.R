library(tidyverse)
library(lubridate)

# 1st southwestern location
imputed_state3 = rep('TX', 9)
imputed_county3 = rep('brazos', 9)
imputed_date3 = ymd(c('2008/01/01',
                      '2009/01/01',
                      '2010/01/01',
                      '2011/01/01',
                      '2012/01/01',
                      '2013/01/01',
                      '2014/01/01',
                      '2015/01/01',
                      '2016/01/01'))
week(imputed_date3) <- rep(19, 9)

wind08_i3 <- read.csv('wind_data/2008i3.csv') 
wind08_i3 <- wind08_i3 %>% cbind(TIME = str_split_fixed(wind08_i3$DATE, 'T', 2)) %>% 
  select(WNDDIR, WNDSPD, TIME.1, TIME.2) %>%
  mutate(DATE = ymd(TIME.1), IHR = hour(hms(TIME.2))) %>% 
  select(DATE, IHR, WNDDIR, WNDSPD, TIME.2) %>%
  mutate(date = ymd(DATE),
         direction = (WNDDIR + 180)*pi/180)
wind09_i3 <- read.csv('wind_data/2009i3.csv')
wind09_i3 <- wind09_i3 %>% cbind(TIME = str_split_fixed(wind09_i3$DATE, 'T', 2)) %>% 
  select(WNDDIR, WNDSPD, TIME.1, TIME.2) %>%
  mutate(DATE = ymd(TIME.1), IHR = hour(hms(TIME.2))) %>% 
  select(DATE, IHR, WNDDIR, WNDSPD, TIME.2) %>%
  mutate(date = ymd(DATE),
         direction = (WNDDIR + 180)*pi/180)
wind10_i3 <- read.csv('wind_data/2010i3.csv') 
wind10_i3 <- wind10_i3 %>% cbind(TIME = str_split_fixed(wind10_i3$DATE, 'T', 2)) %>% 
  select(WNDDIR, WNDSPD, TIME.1, TIME.2) %>%
  mutate(DATE = ymd(TIME.1), IHR = hour(hms(TIME.2))) %>% 
  select(DATE, IHR, WNDDIR, WNDSPD, TIME.2) %>%
  mutate(date = ymd(DATE),
         direction = (WNDDIR + 180)*pi/180)

wind_data_i3 <- rbind(filter(select(wind08_i3, date, direction), 
                             date >= imputed_date3[1]),
                      filter(select(wind09_i3, date, direction), 
                             date >= imputed_date3[2]),
                      filter(select(wind10_i3, date, direction), 
                             date >= imputed_date3[3])) %>%
  filter(month(date) < 8) %>%
  group_by(year = year(date)) %>%
  as_tibble()

wind_data_i3

wind08_i3 %>% head()

ggplot(wind08_i3,
       aes(x = WNDDIR,
           group = IHR,
           color = IHR)) +
  geom_freqpoly(bins = 60) +
  geom_vline(xintercept = 360, color = 2)

ggplot(wind08_i3,
       aes(x = WNDSPD)) +
  geom_freqpoly(bins = 60) 

n_distinct(wind08_i3$WNDSPD)
n_distinct(wind08_i3$WNDDIR)

unique(wind08_i3$WNDSPD)
filter(wind08_i3, WNDSPD > 100)

wind08_i3 %>%
  filter(WNDDIR > 360) %>%
  ggplot(aes(x = IHR)) +
  geom_freqpoly(bins = 24)

unique(filter(wind08_i3, WNDDIR > 360)$WNDSPD)

wind08_i3 %>%
  filter(WNDDIR > 360) %>%
  ggplot(aes(x = WNDSPD)) +
  geom_freqpoly(bins = 24)

wind08_i3 %>%
  filter(WNDDIR < 360,
         month(DATE) %in% 4:7) %>%
  ggplot(aes(x = IHR, 
             group = month(DATE), 
             color = month(DATE))) +
  geom_freqpoly(bins = 24)


wind08_i3 %>%
  filter(WNDDIR < 360,
         month(DATE) %in% 4:7) %>%
  group_by(day = yday(DATE)) %>%
  summarize(n_hrs = n_distinct(IHR)) %>%
  ggplot(aes(x = day, y = n_hrs)) +
  geom_path()

wind08_i3 %>%
  filter(WNDDIR < 360,
         month(DATE) %in% 4:7) %>%
  group_by(day = yday(DATE)) %>%
  summarize(n_hrs = n_distinct(IHR)) %>%
  ggplot(aes(x = n_hrs)) +
  geom_freqpoly(bins = 24)

wind08_i3 %>%
  filter(WNDDIR < 360,
         month(DATE) %in% 4:7) %>%
  group_by(day = yday(DATE)) %>%
  summarize(n_hrs = n_distinct(IHR)) %>%
  summarize(mean(n_hrs < 16))

wind08_s <- read.csv('wind_data/2008s.csv') %>%
  mutate(date = ymd(DATE),
         direction = (WNDDIR + 180)*pi/180)

ggplot(wind08_s,
       aes(x = WNDDIR,
           group = IHR,
           color = IHR)) +
  geom_freqpoly(bins = 30) +
  geom_vline(xintercept = 360, color = 2)

ggplot(wind08_s,
       aes(x = WNDDIR)) +
  geom_freqpoly(bins = 30) +
  geom_vline(xintercept = 360, color = 2)

n_distinct(wind08_s$WNDDIR)

wind08_s %>%
  filter(month(date) %in% 4:7) %>%
  group_by(day = yday(date)) %>%
  summarize(n_hrs = n_distinct(IHR)) %>%
  ggplot(aes(x = day, y = n_hrs)) +
  geom_path()

filter(wind08_s, WNDSPD < 0.5)
range(wind08_s$WNDSPD)
unique(wind08_s$WNDSPD)
hist(wind08_s$WNDSPD)


nrow(wind08_s)
3000/nrow(filter(wind08_s, IHR < 15))

remove <- wind08_s %>%
  mutate(rownum = 1:nrow(wind08_s)) %>%
  filter(IHR < 15) %>%
  sample_n(3000) %>%
  select(rownum)


part <- wind08_s[-remove$rownum, ] %>%
  filter(month(date) %in% 4:7) %>%
  ggplot(aes(x = WNDDIR)) +
  geom_freqpoly()

full <- wind08_s %>%
  filter(month(date) %in% 4:7) %>%
  ggplot(aes(x = WNDDIR)) +
  geom_freqpoly()

grid.arrange(part, full)

wind08_s[-remove$rownum, ] %>%
  filter(month(date) %in% 4:7) %>%
  group_by(yday(date)) %>%
  summarize(n_hrs = n_distinct(IHR)) %>%
  summarize(mean(n_hrs < 23))
