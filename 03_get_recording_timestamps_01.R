# Get timestamps for all WAV files from 2024 pairs that contain at least one usable call
# ONLY WORKS on computer that wrote original files! Acer laptop B405 running Windows.

library(tidyverse)
library(lubridate)
library(clock)

# get list of post-filter usable calls
usable_calls <- readRDS("C:/Users/jkvrt/Downloads/vocal_data_2024-pairs_transformed.RDS") %>%
  select(sound.files, caller) %>%
  mutate(ID = paste(sound.files, caller, sep = "_"))

# get list of all selections made by energy_detector; filter with usable calls
sels <- readRDS("C:/Users/jkvrt/Downloads/selection_tables_05-04-25.RDS") %>%
  filter(selection_length > 0.003) %>%
  separate(source, into = c("date","pair","caller","csv"), sep = "/") %>%
  mutate(ID = paste(sound.files, selec, caller, sep = "_")) %>%
  filter(ID %in% usable_calls$ID)

d <- sels %>%
  mutate(filepath = paste("C:/Users/jkvrt/Documents/2024_pair_recordings",
                          date, pair, caller, sound.files, sep = "/")) %>%
  mutate(start.time = file.info(filepath)$ctime) %>% 
  mutate(end.time = file.info(filepath)$mtime)

options(digits.secs = 6)

saveRDS(d, "C:/Users/jkvrt/Downloads/2024_recording_timestamps.RDS")
