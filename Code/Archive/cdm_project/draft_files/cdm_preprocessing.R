options(stringsAsFactors = F)
library(dplyr)
library(tidyr)
library(stringr)
library(lubridate)

char_remove <- function(y){y %>%
    str_replace_all('[[\\D]-[., -]]', '') %>%
    str_split('-', simplify = T) %>%
    apply(1, function(x){max(as.numeric(x), na.rm = T)}) %>%
    na_if(-Inf)}

planting_types <- c('ommercial', 'arden', 'entinel', 'esearch')
planting_fun1 <- function(x){
  str_detect(x, planting_types) %>% 
    which() %>%
    c(0) %>%
    max()}
planting_fun2 <- function(x){
  lapply(1:length(x), function(j){planting_fun1(x[j])}) %>%
    Reduce(f = 'c') %>%
    max()
}
planting_fun3 <- function(x){
  x %>% 
    str_split(' ') %>%
    lapply(planting_fun2) %>%
    Reduce(f = 'c')
}

host_types <- c('ucumber', 'antaloupe', 'umpkin', 'quash', 'ermelon')
host_fun1 <- function(x){
  str_detect(x, host_types) %>% 
    which() %>%
    c(0)}
host_fun2 <- function(x){
  hosts <- lapply(1:length(x), function(j){host_fun1(x[j])}) %>%
    Reduce(f = 'c') 
  if(sum(hosts > 0) > 1){
    out <- 6
  }else{
    out <- max(hosts)
  }
  return(out)
}
host_fun3 <- function(x){
  x %>% 
    str_extract_all('[a-z]+') %>%
    lapply(host_fun2) %>%
    Reduce(f = 'c')
}

raw_data08 <- read.csv('cdm_reports/2008raw.csv')
sub_data08 <- raw_data08 %>%
  transmute(report_date = mdy(Date.of.Report),
            symptom_date = mdy(Date.of.1st.symptoms),
            reporter = Reporter,
            state = State,
            county = County,
            planting = planting_fun3(Type.of.planting),
            host = host_fun3(Host),
            field_size = NA,
            size_unit = NA,
            incidence = char_remove(X..infected),
            leaf_area = NA,
            occurrence = (incidence != 0),
            lat = Latitude,
            long = Longitude) # variables of interest

raw_data09 <- read.csv('cdm_reports/2009raw.csv')
sub_data09 <- raw_data09 %>%
  transmute(report_date = mdy(Report_date),
            symptom_date = mdy(Symptom_date),
            reporter = Reporter,
            state = State,
            county = County,
            planting = planting_fun3(Planting_type),
            host = host_fun3(Host),
            field_size = NA,
            size_unit = NA,
            leaf_area = char_remove(X.Leaf_area),
            incidence = char_remove(X.Incidence),
            occurrence = (incidence != 0),
            lat = Latitude,
            long = Longitude) # variables of interest


raw_data10 <- read.csv('cdm_reports/2010raw.csv')
str(raw_data10)
sub_data10 <- raw_data10 %>%
  transmute(report_date = mdy(Report_date),
            symptom_date = mdy(Symptom_date),
            reporter = Reporter,
            state = State,
            county = County,
            planting = planting_fun3(Planting_type),
            host = host_fun3(Host),
            field_size = NA,
            size_unit = NA,
            leaf_area = char_remove(X.Leaf_area),
            incidence = char_remove(X.Incidence),
            occurrence = (incidence != 0),
            lat = Latitude,
            long = Longitude) # variables of interest


raw_data11 <- read.csv('cdm_reports/2011raw.csv')
sub_data11 <- raw_data11 %>%
  transmute(report_date = mdy(Report_date),
            symptom_date = mdy(Symptom_date),
            reporter = Reporter,
            state = State,
            county = County,
            planting = planting_fun3(Planting_type),
            host = host_fun3(Host),
            field_size = NA,
            size_unit = NA,
            leaf_area = char_remove(X.Leaf_area),
            incidence = char_remove(X.Incidence),
            occurrence = (incidence != 0),
            lat = Latitude,
            long = Longitude) # variables of interest

raw_data12 <- read.csv('cdm_reports/2012raw.csv')
sub_data12 <- raw_data12 %>%
  transmute(report_date = mdy(Report_date),
            symptom_date = mdy(Symptom_date),
            reporter = Reporter,
            state = State,
            county = County,
            planting = planting_fun3(Planting_type),
            host = host_fun3(Host),
            field_size = NA,
            size_unit = NA,
            leaf_area = char_remove(X.Leaf_area),
            incidence = char_remove(X.Incidence),
            occurrence = (incidence != 0),
            lat = Latitude,
            long = Longitude) # variables of interest

raw_data13 <- read.csv('cdm_reports/2013raw.csv')
sub_data13 <- raw_data13 %>%
  transmute(report_date = mdy(Report_date),
            symptom_date = mdy(Symptom_date),
            reporter = Reporter,
            state = State,
            county = County,
            planting = planting_fun3(Planting_type),
            host = host_fun3(Host),
            field_size = NA,
            size_unit = NA,
            leaf_area = char_remove(X.Leaf_area),
            incidence = char_remove(X.Incidence),
            occurrence = (incidence != 0),
            lat = Latitude,
            long = Longitude) # variables of interest


raw_data14 <- read.csv('cdm_reports/2014raw.csv')
sub_data14 <- raw_data14 %>%
  transmute(report_date = mdy(Report_date),
            symptom_date = mdy(Symptom_date),
            reporter = Reporter,
            state = State,
            county = County,
            planting = planting_fun3(Planting.Type),
            host = host_fun3(Host),
            field_size = Field.Size,
            size_unit = Field.Unit,
            leaf_area = char_remove(Leaf.Area),
            incidence = char_remove(Field.Incidence),
            occurrence = (incidence != 0),
            lat = Latitude,
            long = Longitude) # variables of interest

raw_data15 <- read.csv('cdm_reports/2015raw.csv')
sub_data15 <- raw_data15 %>%
  transmute(report_date = mdy(Report_date),
            symptom_date = mdy(Symptom_date),
            reporter = Reporter,
            state = State,
            county = County,
            planting = planting_fun3(Planting.Type),
            host = host_fun3(Host),
            field_size = Field.Size,
            size_unit = Field.Unit,
            leaf_area = char_remove(Leaf.Area),
            incidence = char_remove(Field.Incidence),
            occurrence = (incidence != 0),
            lat = Latitude,
            long = Longitude) # variables of interest

raw_data16 <- read.csv('cdm_reports/2016raw.csv')
sub_data16 <- raw_data16 %>%
  transmute(report_date = mdy(Report_date),
            symptom_date = mdy(Symptom_date),
            reporter = Reporter,
            state = State,
            county = County,
            planting = planting_fun3(Planting.Type),
            host = host_fun3(Host),
            field_size = Field.Size,
            size_unit = Field.Unit,
            leaf_area = char_remove(Leaf.Area),
            incidence = char_remove(Field.Incidence),
            occurrence = (incidence != 0),
            lat = Latitude,
            long = Longitude) # variables of interest



data <- rbind(sub_data08,
              sub_data09,
              sub_data10,
              sub_data11,
              sub_data12,
              sub_data13,
              sub_data14,
              sub_data15,
              sub_data16)

locations <- distinct(data, lat, long) %>%
  mutate(locID = 1:(nrow(distinct(data, lat, long))))

data <- merge(data, locations, by = c('lat', 'long')) %>%
  mutate(planting = factor(planting, 
                               labels = c('other', 
                                          'commercial',
                                          'garden',
                                          'sentinel',
                                          'research')),
         host = factor(host, 
                           labels = c('other', 
                                      'cucumber',
                                      'cantaloupe',
                                      'pumpkin',
                                      'squash',
                                      'watermelon',
                                      'multiple')))


save(list = 'data', file = 'cdm_data.RData')
write.csv(data, file = 'cdm_data.csv')
