# filter calls that are too short and transform time variables
# Julia Vrtilek, May 2025

# load packages
library(tidyverse)

# set working directory
setwd("/Users/jkvrtilek/Desktop/OSU/PhD/GitHub/pairs")

# combine fundamental frequency measures and spectro_analysis measures
spec <- readRDS("spectro_analysis_2025-05-06.RDS") %>% 
  mutate(ID = paste(sound.files, caller, sep = "_"))

ff <- readRDS("fundfreq_measures_2025-05-08.RDS") %>% 
  mutate(ID = paste(sound.files, caller, sep = "_")) %>% 
  select(ID, maxslope:segments)

d <- inner_join(spec, ff, by = "ID") %>% 
  mutate(date = substring(date, 1, 10)) %>% 
  select(!ID)

# convert seconds to milliseconds
d$duration <- d$duration * 1000
d$time.median <- d$time.median * 1000
d$time.Q25 <- d$time.Q25 * 1000
d$time.Q75 <- d$time.Q75 * 1000
d$time.IQR <- d$time.IQR * 1000

# filter duration, peak frequency, and time variables to remove sounds that are not contact calls
d2 <- d %>% 
  filter(duration > 3) %>% 
  filter(duration < 50) %>% 
  filter(peakf > 10) %>% 
  filter(time.Q25 > 0) %>% 
  filter(time.median > 0) %>% 
  filter(time.Q75 > 0) %>% 
  filter(time.IQR > 0) %>% 
  filter(!is.na(meanslope)) %>% 
  filter(!is.nan(meanslope))

# remove sounds that were manually designated not bat calls
x <- list.files("/Users/jkvrtilek/Desktop/OSU/PhD/GitHub/pairs/not_batcalls/")

not.batcalls <- as.data.frame(x) %>% 
  mutate(x2 = substr(x,1,nchar(x)-5)) %>% 
  separate(x2, into = c("date","numbers"), sep = "_") %>% 
  separate(numbers, into = c("wav","selec"), sep = "-") %>% 
  mutate(filter = paste(date, wav, selec, sep = "_"))

d3 <- d2 %>% 
  filter(!sound.files %in% not.batcalls$filter)

saveRDS(d3, "vocal_data_2024-pairs.RDS")
