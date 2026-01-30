# Extract start time of RECORDING
# Add start time of SELECTION
# Arrange by time, group by session
# Check all the cases where caller changes (lag/lead) - will get NAs between session
# Get lag time between pairs
# Make intervals plot from white-winged antiphonal calling paper

library(tidyverse)
library(warbleR)
library(clock)
library(data.table)

setwd("/Users/jkvrtilek/Desktop/OSU/PhD/GitHub/calling-rate")

# read list of usable calls with creation and modification timestamps 
d <- readRDS("2024_recording_timestamps_2025-12-08.RDS") %>% 
  mutate(filepath = paste("/Users/jkvrtilek/Desktop/OSU/PhD/Ch3/2024_pair_recordings",
                          date, pair, caller, sound.files, sep = "/")) %>% 
  mutate(length.wav = NA) %>% 
  rename(creation = start.time, modification = end.time)

# get length of WAV files
for(i in 1:nrow(d)) {
  x <- read_sound_file(d$filepath[i])
  d$length.wav[i] <- length(x@left)/500000
  print(i)
}

# check whether (creation time) matches (modification time - length of recording)
test <- d %>% 
  mutate(calc.start = modification - length.wav) %>% 
  mutate(diff = calc.start - creation)

# (creation) consistently LATER than (modification - recording length)
# this makes sense because mic was triggered to create the file when call began
# use (modification - recording length) for best results:
# onset of call = modification time of WAV file - length of WAV file + time past beginning of recording

# get lag time between each pair of consecutive calls within session ----
# filter by whether the caller changed between those calls
# MAXIMUM LAG FILTER to make graph interpretable = 10 s
d2 <- d %>% 
  mutate(session = paste(pair,date,sep="_")) %>% 
  mutate(call.start = modification - length.wav + start) %>% 
  mutate(duration = end - start) %>% 
  mutate(call.end = call.start + duration) %>%
  dplyr::select(sound.files:start,date:caller,call.start,call.end,session) %>% 
  arrange(call.start) %>% 
  group_by(session) %>% 
  mutate(lag.time = call.start - lag(call.end)) %>% 
  mutate(lagsec = as.numeric(lag.time)) %>% 
  mutate(caller_change = lag(caller) != caller) %>% 
  filter(caller_change == T) %>% 
  filter(lag.time < 5) %>% # to make plot readable - excludes 123 calls
  filter(lag.time > .0018) # length of time it would take call to travel between bats - excludes 3 calls

p <- ggplot(d2, aes(x = lagsec)) + 
  geom_histogram(binwidth = 0.1, col = "white") +
  scale_x_continuous(breaks = seq(0, 5, by=1)) +
  geom_vline(xintercept = 0.5, color = "red") +
  xlab("interval (s)") +
  ylab("number of calls") +
  theme_bw()
p

# save plot
ggsave(
  "results/call_intervals.jpg",
  plot = p,
  scale = 1,
  width = 6.5,
  height = 2.3,
  units = c("in", "cm", "mm", "px"),
  dpi = 300)


# probability from 2008 white-winged vampire bat paper ----
# number of antiphonal calls/total number of calls
n.antiph <- d2 %>% 
  filter(lag.time < 0.5) %>% 
  nrow()

n.total <- nrow(d)

diaemus.comp <- n.antiph/n.total

# get 95% CI
binom.test(n.antiph, n.total)


# probability that bat GETS response to call ----
# number of antiphonal calls bat received/number of calls bat made
n.received.anti <- d2 %>% 
  filter(lag.time < 0.5) %>% 
  separate(pair, into = c("batA","batB"), remove = F) %>% 
  mutate(focal.bat = case_when(caller == batA ~ batB,
                               caller == batB ~ batA)) %>% 
  group_by(focal.bat) %>% 
  summarize(n.r = n())

n.made.all <- d %>% 
  mutate(focal.bat = caller) %>% 
  group_by(focal.bat) %>% 
  summarize(n.made = n())

pget <- n.made.all %>% 
  left_join(n.received.anti, by = "focal.bat")

prob.get <- pget %>% 
  add_row(focal.bat = "OVERALL", n.r = sum(pget$n.r, na.rm = T), n.made = sum(pget$n.made)) %>% 
  mutate(focal.bat = fct_rev(factor(focal.bat,
                                    levels = c("alana","beast","creep","kazul",
                                               "lurch","mochi","yikes","OVERALL")))) %>% 
  mutate(prob = n.r/n.made) %>% 
  replace(is.na(.), 0) %>% 
  mutate(low = NA) %>% 
  mutate(high = NA) %>% 
  mutate(panel = "getting response to a call")

for (m in 1:nrow(prob.get)) {
  prob.get$low[m] <- binom.test(prob.get$n.r[m], prob.get$n.made[m])$conf.int[[1]]
  prob.get$high[m] <- binom.test(prob.get$n.r[m], prob.get$n.made[m])$conf.int[[2]]
}


# probability that bat GIVES response to call ----
# number of antiphonal calls bat made/total number of calls bat received
n.made.anti <- d2 %>% 
  filter(lag.time < 0.5) %>% 
  mutate(focal.bat = caller) %>% 
  group_by(focal.bat) %>% 
  summarize(n.m = n())

n.received.all <- d %>% 
  separate(pair, into = c("batA","batB"), remove = F) %>% 
  mutate(focal.bat = case_when(caller == batA ~ batB,
                               caller == batB ~ batA)) %>% 
  group_by(focal.bat) %>% 
  summarize(n.received = n())

pgive <- n.received.all %>% 
  left_join(n.made.anti, by = "focal.bat")

prob.give <- pgive %>% 
  add_row(focal.bat = "OVERALL", n.m = sum(pgive$n.m, na.rm = T), n.received = sum(pgive$n.received)) %>%
  mutate(focal.bat = fct_rev(factor(focal.bat,
                                    levels = c("alana","beast","creep","kazul",
                                               "lurch","mochi","quark","yikes","OVERALL")))) %>% 
  mutate(prob = n.m/n.received) %>% 
  replace(is.na(.), 0) %>% 
  mutate(low = NA) %>% 
  mutate(high = NA) %>% 
  mutate(panel = "giving response to a call")

for (m in 1:nrow(prob.give)) {
  prob.give$low[m] <- binom.test(prob.give$n.m[m], prob.give$n.received[m])$conf.int[[1]]
  prob.give$high[m] <- binom.test(prob.give$n.m[m], prob.give$n.received[m])$conf.int[[2]]
}


# probabilities plot ----
pp <- ggplot() +
  geom_errorbarh(data = prob.give, aes(y=focal.bat, xmin=low, xmax=high, height = 0.3),
                 size=1, color = "gray") +
  geom_errorbarh(data = prob.get, aes(y=focal.bat, xmin=low, xmax=high, height = 0.3),
                 size=1, color = "gray") +
  geom_point(data = prob.give, aes(x=prob, y=focal.bat, size=n.m)) +
  geom_point(data = prob.get, aes(x=prob, y=focal.bat, size=n.made)) +
  facet_grid(~panel) +
  ylab("focal bat") +
  xlab("probability") +
  theme_bw()
pp

# save plot
ggsave(
  "results/getgive_probability.jpg",
  plot = pp,
  scale = 1,
  width = 9,
  height = 4,
  units = c("in", "cm", "mm", "px"),
  dpi = 300)

# correlation between giving response and number of calls made
n.made.all %>% 
  left_join(prob.give, by = "focal.bat") %>% 
  ggplot(aes(x = n.made, y = prob)) +
  geom_point() +
  geom_smooth(method = "lm") +
  ylab("probability of giving response to a call") +
  xlab("number of calls made") +
  theme_bw()


# simulation ----

# get start time of each session
session.time <- read_csv("trial_times.csv") %>% 
  separate(date, into = c("y","m","d"), sep = "-", remove = F) %>% 
  rename(creation = start.time) %>% 
  mutate(time = as.POSIXct(paste(date, creation))) %>% 
  # get session ID
  mutate(month = as.integer(m),
         year = as.integer(y) - 2000) %>% 
  mutate(batA = tolower(batA),
         batB = tolower(batB)) %>% 
  mutate(date = paste(month, d, year, sep = "-"),
         pair = paste(batA,batB,sep="-"),
         session = paste(pair,date,sep="_")) %>% 
  # number trials
  arrange(time) %>% 
  mutate(trial = row_number()) %>% 
  select(session, time, trial)

d3 <- d %>%
  mutate(session = paste(pair,date,sep="_")) %>% 
  mutate(call.start = modification - length.wav + start) %>% 
  mutate(duration = end - start) %>% 
  dplyr::select(sound.files:start,date:caller,call.start,duration,session) %>% 
  left_join(session.time, by = "session") %>% 
  # get time between start of recording session and start of call
  mutate(time.from.start = call.start - time) %>%
  # get time call ended
  mutate(call.offset = time.from.start + duration/60) %>% # time.from.start is in mins, duration is in secs
  # prepare to split data
  # designate "left" and "right" bat for each session
  separate(pair, into = c("left","right"), remove = F) %>% 
  # get whether bat that made each call is "left" or "right"
  mutate(LR = case_when(caller == left ~ "left",
                        caller == right ~ "right",
                        .default = NA)) %>% 
  select(session, trial, LR, caller, time.from.start, call.offset) %>% 
  group_by(trial) %>% 
  arrange(time.from.start, .by_group = T) %>% 
  ungroup()

# split data into "left" and "right" bats
left <- d3 %>%
  filter(LR == "left") %>% 
  mutate(rand.trial = NA)
right <- d3 %>%
  filter(LR == "right")
  
# prep for loop
sims <- 5000

results <- setNames(data.frame(matrix(ncol = 1, nrow = sims)), "antiphon")

for (j in 1:sims) {
  
  set.seed(123+j)
  
  # assign random "trial" number to each real trial
  r.trial <- right %>% 
    select(trial) %>% 
    distinct() %>% 
    mutate(rand.trial = sample(1:84, length(unique(right$trial)), replace = F))
  
  # assign random trial # to "right" bat calls
  r.sim <- right %>% 
    left_join(r.trial, by = "trial")
  
  # keep real trial # for "left" bat calls
  left$rand.trial <- left$trial
  
  # combine such that real left-bat trials are paired with a random right-bat trial
  sim <- rbind(r.sim,left)
  
  # get all trials where a bat was randomly paired with itself
  same.bat <- sim %>% 
    select(rand.trial,LR,caller) %>% 
    distinct() %>% 
    pivot_wider(names_from = LR, values_from = caller) %>% 
    filter(left == right)
  
  sim.antiph <- sim %>% 
    filter(!(rand.trial %in% same.bat$rand.trial)) %>% # remove paired-with-self
    group_by(rand.trial) %>% 
    arrange(time.from.start, .by_group = T) %>% 
    mutate(caller_change = lag(caller) != caller) %>% 
    filter(caller_change == T) %>% 
    mutate(lag.time = time.from.start - lag(call.offset)) %>% 
    mutate(lagsec = as.numeric(lag.time)) %>% 
    filter(lagsec < 0.5) %>% 
    filter(lagsec > 0.0018) # average length of call - overlapping calls don't count
  
  results$antiphon[j] <- nrow(sim.antiph)
  
  print(j)
}

range95 <- round(quantile(results$antiphon, probs= c(0.025, 0.975), na.rm=T),3)

obs <- n.antiph/n.total
exp.low <- range95[1]/n.total
exp.high <- range95[2]/n.total
