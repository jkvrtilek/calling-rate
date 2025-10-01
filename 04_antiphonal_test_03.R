# Extract start time of RECORDING (end time - length of wav)
# Then add start time of SELECTION
# Arrange by start time, group by session
# Check all the cases where caller changes (lag/lead) - should get NAs between session
# Get lag time between pairs
# Make intervals plot from white-winged antiphonal calling paper

library(tidyverse)
library(warbleR)
library(clock)
library(data.table)

setwd("/Users/jkvrtilek/Desktop/OSU/PhD/GitHub/pairs")

# read list of usable calls with creation and modification timestamps 
d <- readRDS("2024_recording_timestamps.RDS") %>% 
  mutate(filepath = paste("/Users/jkvrtilek/Desktop/OSU/PhD/Ch3/2024_pair_recordings",
                          date, pair, caller, sound.files, sep = "/")) %>% 
  mutate(length.wav = NA)

# get length of WAV files
for(i in 1:nrow(d)) {
  x <- read_sound_file(d$filepath[i])
  d$length.wav[i] <- length(x@left)/500000
  print(i)
}

# check whether creation time matches modification time - length of recording
test <- d %>% 
  mutate(calc.start = end.time - length.wav) %>% 
  mutate(diff = calc.start - start.time)

# creation time consistently LATER than mtime - recording length
# this makes sense because mic was triggered to create file after call began
# use mtime - recording length for best results

# get lag time between each pair of consecutive calls within session ----
# filter by whether the caller changed between those calls
# MAXIMUM LAG FILTER to make graph interpretable = 10 s
d2 <- d %>% 
  mutate(session = paste(pair,date,sep="_")) %>% 
  mutate(rec.start = end.time - length.wav) %>% 
  mutate(call.start = rec.start + start) %>% 
  mutate(duration = end - start) %>% 
  mutate(call.end = call.start + duration) %>% 
  dplyr::select(sound.files:start,date:caller,call.start,call.end,session) %>% 
  arrange(call.start) %>% 
  group_by(session) %>% 
  mutate(lag.time = call.start - lag(call.end)) %>% 
  mutate(lagsec = as.numeric(lag.time)) %>% 
  mutate(caller_change = lag(caller) != caller) %>% 
  filter(caller_change == T) %>% 
  filter(lag.time < 10) %>% # to make plot readable - excludes 80 calls
  filter(lag.time > .0018) # length of time it would take call to travel between bats - does not exclude any calls

p <- ggplot(d2, aes(x = lagsec)) + 
  geom_histogram(binwidth = 0.25) +
  scale_x_continuous(breaks = seq(0, 10, by=1)) +
  geom_vline(xintercept = 0.5, color = "red") +
  xlab("interval (s)") +
  ylab("number of calls") +
  theme_bw()
p


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
                                    levels = c("alana","beast","creep","kazul","lurch","mochi","yikes","OVERALL")))) %>% 
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
                                    levels = c("alana","beast","creep","kazul","lurch","mochi","quark","yikes","OVERALL")))) %>% 
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
  geom_point(data = prob.give, aes(x=prob, y=focal.bat), size=3) +
  geom_point(data = prob.get, aes(x=prob, y=focal.bat), size=3) +
  geom_errorbarh(data = prob.give, aes(y=focal.bat, xmin=low, xmax=high), size=1) +
  geom_errorbarh(data = prob.get, aes(y=focal.bat, xmin=low, xmax=high), size=1) +
  facet_grid(~panel) +
  ylab("focal bat") +
  xlab("probability") +
  theme_bw()
pp

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

d3 <- d %>%
  mutate(session = paste(pair,date,sep="_")) %>% 
  mutate(use.start = end.time - length.wav) %>% 
  mutate(call.start = use.start + start) %>% 
  mutate(duration = end - start) %>% 
  dplyr::select(sound.files:start,date:caller,call.start,duration,session)

sims <- 5000

results <- setNames(data.frame(matrix(ncol = 2, nrow = sims)),
                    c("antiphon","overlap"))

# get earliest and latest call per session
d4 <- d3 %>% 
  group_by(session) %>% 
  mutate(earliest = min(call.start)) %>% 
  mutate(latest = max(call.start)) %>% 
  mutate(seclength = latest - earliest) %>% # number of seconds between onset of first and last calls in session
  mutate(add.sec = 0) %>% 
  ungroup()

# make list of one df per session
bysession <- d4 %>% 
  group_split(session)

for (j in 1:sims) {
  
  for (l in 1:length(bysession)) {
    
    bysession[[l]]$add.sec <- runif(nrow(bysession[[l]]), min = 0, max = bysession[[l]]$seclength[1])
    
  }
    
    # simulate data within blocks
    sim.d <- bysession %>% 
      rbindlist() %>% 
      mutate(call.start.sim = earliest + add.sec) %>% 
      mutate(call.end.sim = call.start.sim + duration) %>% 
      mutate(caller_change = lag(caller) != caller) %>% 
      filter(caller_change == F) %>% 
      arrange(call.start.sim) %>%
      mutate(overlap = lag(call.end.sim) > call.start.sim)
  
  sim.antiph <- sim.d %>% 
    arrange(call.start.sim) %>% 
    group_by(session) %>% 
    mutate(lag.time = call.start.sim - lag(call.end.sim)) %>% 
    mutate(lagsec = as.numeric(lag.time)) %>% 
    mutate(caller_change = lag(caller) != caller) %>% 
    filter(caller_change == T) %>% 
    filter(lagsec < 0.5) %>% 
    filter(lagsec > 0.0018)
  
  results$antiphon[j] <- nrow(sim.antiph)
  
  results$overlap[j] <- sim.d %>% 
    filter(overlap == T) %>% 
    nrow()
  
  print(j)
}

range95 <- round(quantile(results$antiphon, probs= c(0.025, 0.975), na.rm=T),3)

obs <- n.antiph/n.total
exp.low <- range95[1]/n.total
exp.high <- range95[2]/n.total
