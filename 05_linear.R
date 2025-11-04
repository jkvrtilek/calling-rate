# analysis for 2024 paired caller-receiver recordings
# all pairs of 8 bats were recorded with 2 synced mics
# does amount of grooming/foodsharing given predict number of calls?
# Julia Vrtilek and Gerry Carter

library(tidyverse)

setwd("/Users/jkvrtilek/Desktop/OSU/PhD/GitHub/calling-rate")

# load and wrangle vocal data
batcalls <- readRDS("vocal_data_2024-pairs.RDS") %>% 
  group_by(caller) %>% 
  summarize(n.calls = n()) %>% 
  ungroup()

# get allogrooming and allofeeding data
bat.donations <- 
  read.csv('OSU_2024_social_data.csv') %>% 
  mutate(Date = as.Date(Date, format = "%m/%d/%Y")) %>% 
  mutate(caller = tolower(Actor),
         receiver = tolower(Receiver)) %>% 
  filter(caller %in% batcalls$caller) %>% 
  group_by(caller, Behavior) %>% 
  summarize(rate = sum(rate)) %>% 
  ungroup()

# combine
d <- left_join(batcalls, bat.donations, by = "caller")

# aggression ----
d.agg <- d %>% 
  filter(Behavior == "Aggression")

m.agg <- lm(n.calls ~ scale(rate), data = d.agg)
summary(m.agg)

# grooming ----
d.groom <- d %>% 
  filter(Behavior == "Grooming")

m.groom <- lm(n.calls ~ scale(rate), data = d.groom)
summary(m.groom)

# foodsharing ----
d.food <- d %>% 
  filter(Behavior == "Mouthlicking")

m.food <- lm(n.calls ~ scale(rate), data = d.food)
summary(m.food)
