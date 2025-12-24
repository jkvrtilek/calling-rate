# filter calls that are too short and transform time variables
# Julia Vrtilek, May 2025

# load packages
library(tidyverse)
library(viridis)

# set working directory
setwd("/Users/jkvrtilek/Desktop/OSU/PhD/GitHub/calling-rate")

# combine fundamental frequency measures and spectro_analysis measures
spec <- readRDS("spectro_analysis_2025-11-16.RDS") %>% 
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

saveRDS(d3, "vocal_data_2024-pairs_2025-12-08.RDS")

# get table of session, batA calls, batB calls ----
d4 <- d3 %>% 
  group_by(date, undir.pair, caller) %>% 
  summarize(numcalls = n()) %>% 
  mutate(batAB = case_when(caller == substr(undir.pair,1,5) ~ "batA",
                           caller == substr(undir.pair,7,11) ~ "batB")) %>% 
  select(date, undir.pair, batAB, numcalls) %>% 
  pivot_wider(names_from = batAB, values_from = numcalls) %>% 
  replace(is.na(.), 0)

write.csv(d4, "calls_per_session.csv", row.names = F)

# make heatmap of calls per directed dyad
temp1 <- d3 %>% 
  group_by(caller, receiver) %>% 
  summarize(numcalls = n()) %>% 
  mutate(pair = paste(caller, receiver))

temp2 <- expand_grid(caller = sort(unique(d3$receiver)),
                    receiver = sort(unique(d3$receiver))) %>% 
  mutate(pair = paste(caller, receiver))

hm <- left_join(temp2, temp1, by = "pair") %>% 
  mutate(caller = caller.x) %>% 
  mutate(receiver = receiver.x) %>% 
  select(caller, receiver, numcalls) %>% 
  mutate(n.calls = case_when(!is.na(numcalls)~numcalls,
                             is.na(numcalls) & caller != receiver~0))

ggplot(hm, aes(caller, receiver, fill = n.calls)) +
  geom_tile() +
  scale_fill_viridis(option="C", limits = c(1,800), na.value = "white") +
  theme_bw()
