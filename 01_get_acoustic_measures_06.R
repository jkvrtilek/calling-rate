# filter calls that are too short and transform time variables
# Julia Vrtilek, May 2025

# load packages
library(tidyverse)
library(viridis)

# set working directory
setwd("/Users/jkvrtilek/Desktop/OSU/PhD/GitHub/calling-rate")

# combine spectro_analysis measures and fundamental frequency measures ----
spec <- readRDS("spectro_analysis_2025-11-16.RDS") %>% 
  mutate(ID = paste(sound.files, caller, sep = "_"))

ff <- readRDS("fundfreq_measures_2025-05-08.RDS") %>% 
  mutate(ID = paste(sound.files, caller, sep = "_")) %>% 
  select(ID, maxslope:segments)

d <- inner_join(spec, ff, by = "ID") %>% 
  mutate(date = substring(date, 1, 10)) %>% 
  select(!ID)

# convert seconds to milliseconds ----
d$duration <- d$duration * 1000
d$time.median <- d$time.median * 1000
d$time.Q25 <- d$time.Q25 * 1000
d$time.Q75 <- d$time.Q75 * 1000
d$time.IQR <- d$time.IQR * 1000

# filter duration, peak frequency, and time variables to remove sounds that are not contact calls ----
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

# remove sounds that were manually designated not bat calls ----
x <- list.files("/Users/jkvrtilek/Desktop/OSU/PhD/GitHub/calling-rate/not_batcalls/")

not.batcalls <- as.data.frame(x) %>% 
  mutate(x2 = substr(x,1,nchar(x)-5)) %>% 
  separate(x2, into = c("date","numbers"), sep = "_") %>% 
  separate(numbers, into = c("wav","selec"), sep = "-") %>% 
  mutate(filter = paste(date, wav, selec, sep = "_"))

d3 <- d2 %>% 
  filter(!sound.files %in% not.batcalls$filter) %>% 
  mutate(undir.pair = if_else(caller < receiver,
                              paste(caller, receiver, sep = "-"),
                              paste(receiver, caller, sep = "-")))

callstats <- d3 %>% 
  group_by(caller) %>% 
  summarize(n = n()) %>% 
  arrange(n)
mean <- mean(callstats$n)

saveRDS(d3, paste("vocal_data_2024-pairs_", Sys.Date(), ".RDS", sep = ""))

# get table of session, batA calls, batB calls ----
d4 <- d3 %>% 
  group_by(date, undir.pair, caller) %>% 
  summarize(numcalls = n()) %>% 
  mutate(batAB = case_when(caller == substr(undir.pair,1,5) ~ "batA",
                           caller == substr(undir.pair,7,11) ~ "batB")) %>% 
  select(undir.pair, batAB, numcalls, date) %>% 
  pivot_wider(names_from = batAB, values_from = numcalls)  %>% 
  mutate(ID = paste(date, undir.pair, sep = "_")) %>% 
  mutate(date = as.Date(date))

empty_sessions <- read_csv("trial_times.csv") %>% 
  mutate(undir.pair = str_to_lower(if_else(batA < batB,
                              paste(batA, batB, sep = "-"),
                              paste(batB, batA, sep = "-")))) %>% 
  mutate(ID = paste(date, undir.pair, sep = "_")) %>% 
  select(ID, undir.pair, date)

d5 <- full_join(d4, empty_sessions) %>% 
  replace(is.na(.), 0) %>% 
  arrange(undir.pair)

write.csv(d5, "calls_per_session.csv", row.names = F)
