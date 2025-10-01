# analysis for 2024 paired caller-receiver recordings
# do different bats call different amounts? 
# does calling rate depend on recipient?
# all pairs of 8 bats were recorded with 2 synced mics
# Julia Vrtilek and Gerry Carter

# load packages

library(devtools)
# install.packages("cmdstanr", repos = c('https://stan-dev.r-universe.dev', getOption("repos")))
# install_github('ctross/STRAND@phosphorescent_desert_buttons')
library(STRAND)

# for data wranging
library(tidyverse)
library(igraph)
library(broom)

# for permutation tests
library(vegan)
library(asnipe)

# for Bayesian models
library(brms)
library(performance)
library(tidybayes)

setwd("/Users/jkvrtilek/Desktop/OSU/PhD/GitHub/calling-rate")
# setwd("~/Dropbox (Personal)/Dropbox/_working/_ACTIVE/__students/Julia Vrtilek/2025/bat pairs experiment")


# load and wrangle vocal data ----
batcalls <- readRDS("vocal_data_2024-pairs.RDS") %>% 
  group_by(caller, receiver) %>% 
  summarize(n.calls = n()) %>% 
  ungroup() %>% 
  # quark never called; adding this row makes names appear alphabetically
  add_row(caller = "quark", 
          receiver = "yikes", n.calls = 0) %>% 
  arrange(caller)

# get count of calls as matrix
# rows are callers, cols are receivers
mcalls <- 
  batcalls %>% 
  graph_from_data_frame() %>% 
  as_adjacency_matrix(attr= 'n.calls', sparse=F)

# get bats used in caller-receiver study 
bats.used <- rownames(mcalls)
bats.used


# load and wrangle social data ----

# get allogrooming and allofeeding data
bat.donations <- 
  read.csv('OSU_2024_social_data.csv') %>% 
  mutate(Date = as.Date(Date, format = "%m/%d/%Y")) %>% 
  mutate(Actor = tolower(Actor),
         Receiver = tolower(Receiver)) %>% 
  filter(Actor %in% bats.used) %>% 
  filter(Receiver %in% bats.used) 

# make foodsharing matrix
mf <-
  bat.donations %>%
  filter(Behavior == "Mouthlicking") %>% 
  group_by(Actor, Receiver) %>%
  summarize(rate= sum(rate, na.rm=T)) %>%
  ungroup() %>% 
  graph_from_data_frame() %>% 
  as_adjacency_matrix(attr= 'rate', sparse=F)

# make grooming matrix
mg <-
  bat.donations %>%
  filter(Behavior == "Grooming") %>%
  group_by(Actor, Receiver) %>%
  summarize(rate= sum(rate, na.rm=T)) %>%
  ungroup() %>% 
  graph_from_data_frame() %>% 
  as_adjacency_matrix(attr= 'rate', sparse=F)

# make affiliation matrix
ma <- mf + mg

# make kinship matrix
raw.kin <- read.table('KING.txt')
colnames(raw.kin) <- tolower(colnames(raw.kin))
rownames(raw.kin) <- tolower(rownames(raw.kin))

# get kinship as genetic similarity
km <- 
  raw.kin %>% 
  select(all_of(bats.used)) %>% 
  filter(row.names(raw.kin) %in% bats.used) %>% 
  as.matrix()
# note that all these kinship values are negative


# for MCMC 
nchains <- 6
warmup <- 1500
chain_length <- 3000



# BRMS: fit GLMMs that assume actor and receiver have independent effects in calling -----------

matrix_to_df <- function(m){
  data.frame(row=rownames(m)[row(m)], col=colnames(m)[col(m)], value=c(m), stringsAsFactors = F)
}

# get calls
n.calls2 <- 
  mcalls %>% 
  matrix_to_df() %>% 
  filter(row!=col) %>% 
  group_by(row, col) %>% 
  summarize(n.calls= sum(value, na.rm=T)) %>% 
  rename(caller=row, receiver = col) %>% 
  ungroup()

# groom
groom2 <- 
  mg %>% 
  matrix_to_df() %>% 
  filter(row!=col) %>% 
  group_by(row, col) %>% 
  summarize(groom= sum(value, na.rm=T)) %>% 
  rename(caller=row, receiver = col) %>% 
  ungroup()

# feed
feed2 <- 
  mf %>% 
  matrix_to_df() %>% 
  filter(row!=col) %>% 
  group_by(row, col) %>% 
  summarize(feed= sum(value, na.rm=T)) %>% 
  rename(caller=row, receiver = col) %>% 
  ungroup()

# kinship
kinship2 <- 
  km %>% 
  matrix_to_df() %>% 
  filter(row!=col) %>% 
  group_by(row, col) %>% 
  summarize(kinship= sum(value, na.rm=T)) %>% 
  rename(caller=row, receiver = col) %>% 
  ungroup()

# combine
d2 <- 
  n.calls2 %>% 
  full_join(groom2) %>% 
  full_join(feed2) %>% 
  full_join(kinship2) 

# plot variables
d2 %>% 
  pivot_longer(cols= n.calls:kinship) %>% 
  ggplot(aes(x=value, color=name))+
  facet_wrap(~name, scales= "free")+
  geom_histogram(fill="light blue", color="black")

# plot correlations
d2 %>% 
  pivot_longer(cols= groom:kinship) %>% 
  ggplot(aes(x=value, y= n.calls, color=name))+
  facet_wrap(~name, scales= "free", nrow=3)+
  geom_point(size=3, alpha=0.8)+
  geom_smooth(method= "lm")

# fit normal GLMMs 
fit1 <-
  brm(n.calls ~ 
        scale(groom) +
        (1|caller)+
        (1|receiver),
      data = d2, 
      family = negbinomial(),
      seed = 123,
      cores = nchains,
      chains = nchains,
      iter = chain_length,
      warmup = warmup)
fit2 <-
  brm(n.calls ~ 
        scale(feed) +
        (1|caller)+
        (1|receiver),
      data = d2, 
      family = negbinomial(),
      seed = 123,
      cores = nchains,
      chains = nchains,
      iter = chain_length,
      warmup = warmup)
fit3 <-
  brm(n.calls ~ 
        scale(kinship) +
        (1|caller)+
        (1|receiver),
      data = d2, 
      family = negbinomial(),
      seed = 123,
      cores = nchains,
      chains = nchains,
      iter = chain_length,
      warmup = warmup)

# I bet these divergent transitions are due to the outlier case
pp_check(fit1, ndraws=100)+xlim(0, 1000)
pp_check(fit2, ndraws=100)+xlim(0, 1000)
pp_check(fit3, ndraws=100)+xlim(0, 1000)

# get posteriors
t1 <- 
  fit1 %>% 
  spread_draws(b_scalegroom) %>% 
  mutate(model = "Allogrooming") %>% 
  pivot_longer(b_scalegroom, names_to = 'term', values_to= 'coeff') 

t2 <- 
  fit2 %>% 
  spread_draws(b_scalefeed) %>% 
  mutate(model = "Allofeeding") %>% 
  pivot_longer(b_scalefeed, names_to = 'term', values_to= 'coeff') 

t3 <- 
  fit3 %>% 
  spread_draws(b_scalekinship) %>% 
  mutate(model = "Kinship") %>% 
  pivot_longer(b_scalekinship, names_to = 'term', values_to= 'coeff') 

# plot models
(models.plot <- 
    rbind(t1, t2, t3) %>% 
    mutate(model= fct_rev(model)) %>%  
    ggplot(aes(y = model, x = coeff, fill=model)) +
    stat_halfeye(.width = c(0.95), linewidth= 5, size=5)+
    geom_vline(xintercept = 0)+
    ylab("")+
    xlab("estimated effect on calling rate (coefficient)")+
    theme_bw()+
    theme(legend.position= 'none', 
          axis.text.y = element_text(size = 12, vjust = 0.5, hjust = 0.5),
          axis.text.x = element_text(size = 12),
          axis.title.x = element_text(size = 12)
          ))

# save plot
ggsave(
  "results/brms_models_plot.pdf",
  plot = models.plot,
  scale = 1,
  width = 5.5,
  height = 6.5,
  units = c("in", "cm", "mm", "px"),
  dpi = 300)

# get model results
t1 <- summary(fit1)$fixed %>% rownames_to_column("estimate") %>% mutate(model= "allogrooming given") %>% relocate(model)
t2 <- summary(fit2)$fixed %>% rownames_to_column("estimate") %>% mutate(model= "allofeeding given") %>% relocate(model)
t3 <- summary(fit3)$fixed %>% rownames_to_column("estimate") %>% mutate(model= "kinship") %>% relocate(model)

# save results
rbind(t1,t2,t3) %>% 
  write.csv("results/BRMS_model_results.csv")



# STRAND: fit social relations model ---------
# correlated effects of actor and receiver

### create function to save model diagnostics
get_diagnostics <- function(model = model.here, 
                            label= "label here"){
  
  model$fit$summary() %>% 
    select(rhat, ess_bulk, ess_tail) %>% 
    summarize(min.rhat= min(rhat, na.rm=T), 
              max.rhat= max(rhat, na.rm=T),
              min.ess_bulk= min(ess_bulk, na.rm=T), 
              max.ess_bulk= max(ess_bulk, na.rm=T),
              min.ess_tail= min(ess_tail, na.rm=T),
              max.ess_tail= max(ess_tail, na.rm=T)) %>% 
    mutate(model = label) %>% 
    relocate(model)
}

### function to plot posteriors from STRAND
# note: block parameters 1 = intercept
plot_posteriors <- function(results= results, 
                            title= "title here", 
                            effect= "effect here"){
  
  # get samples for focal and target effects SD
  t1 <- 
    results$samples$srm_model_samples$focal_target_sd %>% 
    as_tibble() %>% 
    rename(focal_effects_sd= V1, target_effects_sd= V2) %>% 
    pivot_longer(focal_effects_sd: target_effects_sd,
                 names_to = "label", 
                 values_to = "value") 
  
  # get samples for dyadic effect SD
  t2 <- 
    results$samples$srm_model_samples$dyadic_sd %>% 
    as_tibble() %>% 
    rename(value= V1) %>% 
    mutate(label= "dyadic_effect_sd") 
  
  
  # get samples for dyadic effect coeff
  t3 <- 
    results$samples$srm_model_samples$dyadic_coeffs %>% 
    as_tibble() %>% 
    rename(value= V1) %>% 
    mutate(label= "dyadic_effect_coeff") 
  
  
  # get samples for focal-target effects rho (generalized reciprocity)
  t4 <- 
    results$samples$srm_model_samples$focal_target_L %>% 
    as_tibble() %>% 
    rename(value= `2.1`) %>% 
    select(value) %>% 
    mutate(label= "focal-target_effects_rho")
  
  # get samples for dyadic effects rho (dyadic reciprocity)
  t5 <- 
    results$samples$srm_model_samples$dyadic_L %>% 
    as_tibble() %>% 
    rename(value= `2.1`) %>% 
    select(value) %>% 
    mutate(label= "dyadic_effects_rho") 
  
  # combine results
  (tplot <- 
      rbind(t1,t2,t3,t4,t5) %>% 
      mutate(label= case_when(
        label== 'focal_effects_sd' ~ "2. caller identity (SD)",
        label=='target_effects_sd' ~ "3. receiver identity (SD)",
        label=='dyadic_effect_sd' ~ "4. pair identity (SD)",
        label=='dyadic_effect_coeff' ~ paste0("1. ",effect),
        label=='dyadic_effects_rho' ~ "5. within-pair correlation",
        label=='focal-target_effects_rho' ~ "6. generalized correlation")) %>% 
      mutate(label= fct_rev(label)) %>% 
      ggplot(aes(y = label, x = value, fill=label)) +
      stat_halfeye(.width = c(0.95), linewidth= 5, size=5)+
      geom_vline(xintercept = 0)+
      ylab("")+
      xlab("model coefficient")+
      ggtitle(title)+
      theme_bw()+
      theme(legend.position= 'none', 
            axis.text.y = element_text(size = 10, vjust = 0.5, hjust = 0.5)))
  
  return(tplot)
  
}

# make the STRAND data structure ----
chars <- read.csv('campus_bat_chars.csv', stringsAsFactors = F)
chars$Bat.name <- tolower(chars$Bat.name)

age <- 
  chars %>% 
  filter(Bat.name %in% bats.used) %>% 
  select(Age) %>% 
  mutate(Age= standardize(Age))

rownames(age) <- bats.used

# scale and store in list
dyad = list(Kinship = scale(km),
            Feeding = scale(mf),
            Grooming = scale(mg))

### We need to figure out whether to scale the matrix or the vector

# combine into df
dat = make_strand_data(
  outcome = list(calls=mcalls),
  individual_covariates = age,
  dyadic_covariates = dyad,
  outcome_mode = "poisson",
  link_mode = "log",
  check_standardization = F) 


# model effect of allogrooming given -----
sfit1 <- 
  fit_social_relations_model(
    data=dat,
    focal_regression = ~ 1,
    target_regression = ~ 1,
    dyad_regression = ~ Grooming,
    mode="mcmc",
    stan_mcmc_parameters = list(
      chains = nchains,
      iter_warmup = chain_length,
      iter_sampling = warmup))

# summarize results
res1 <-  summarize_strand_results(sfit1)
res1$summary

# save diagnostics
diagnostics1 <- get_diagnostics(sfit1, "allogrooming given")

# save results
strand1 <- 
  res1$summary %>% 
  as_tibble() %>% 
  mutate(Model = "allogrooming given")

# plot results
(plot1a <- strand_caterpillar_plot(res1, normalized=FALSE,  only_slopes=F))
(plot1b <- strand_caterpillar_plot(res1, 
                                   submodels = "Dyadic effects", 
                                   normalized=T,  
                                   only_slopes=T))
(plot1c <- plot_posteriors(res1, "allogrooming given", "allogrooming given"))

# save plot
ggsave(
  "results/allogrooming.pdf",
  plot = plot1c,
  scale = 1,
  width = 5.5,
  height = 6.5,
  units = c("in", "cm", "mm", "px"),
  dpi = 300)

# Normalized effect is
# df$Diff = df$HI - df$LI
# if (normalized == TRUE) {
#  df$Median = df$Median/df$Diff
#  df$LI = df$LI/df$Diff
#  df$HI = df$HI/df$Diff

# there is no clear effect of allogrooming on calling rate after controlling for reciprocity in calling


# model effect of feeding given -----
sfit2 <- 
  fit_social_relations_model(
    data=dat,
    focal_regression = ~ 1,
    target_regression = ~ 1,
    dyad_regression = ~ Feeding,
    mode="mcmc",
    stan_mcmc_parameters = list(
      chains = nchains,
      iter_warmup = chain_length,
      iter_sampling = warmup))

# summarize results
res2 <-  summarize_strand_results(sfit2)
res2$summary

# save diagnostics
diagnostics2 <- get_diagnostics(sfit2, "allofeeding given")

# save results
strand2 <- 
  res2$summary %>% 
  as_tibble() %>% 
  mutate(Model = "allofeeding given")

# plot results
(plot2a <- strand_caterpillar_plot(res2, normalized=FALSE,  only_slopes=F))
(plot2b <- strand_caterpillar_plot(res2, 
                                   submodels = "Dyadic effects", 
                                   normalized=T,  
                                   only_slopes=T)) 
(plot2c <- plot_posteriors(res2, "allofeeding given", "allofeeding given"))

# save plot
ggsave(
  "results/allofeeding.pdf",
  plot = plot2c,
  scale = 1,
  width = 5.5,
  height = 6.5,
  units = c("in", "cm", "mm", "px"),
  dpi = 300)


# model effect of kinship ----
sfit3 <- 
  fit_social_relations_model(
    data=dat,
    focal_regression = ~ 1,
    target_regression = ~ 1,
    dyad_regression = ~ Kinship,
    mode="mcmc",
    stan_mcmc_parameters = list(
      chains = nchains,
      iter_warmup = chain_length,
      iter_sampling = warmup))

# summarize results
res3 <-  summarize_strand_results(sfit3)

# save diagnostics
diagnostics3 <- get_diagnostics(sfit3, "kinship")

# save results
strand3 <- 
  res3$summary %>% 
  as_tibble() %>% 
  mutate(Model = "kinship")

# plot results
(plot3a <- strand_caterpillar_plot(res3, normalized=FALSE,  only_slopes=F))
(plot3b <- strand_caterpillar_plot(res3, 
                                   submodels = "Dyadic effects", 
                                   normalized=T,  
                                   only_slopes=T)) 
(plot3c <- plot_posteriors(res3, "genetic relatedness", "relatedness"))

# save plot
ggsave(
  "results/kinship.pdf",
  plot = plot3c,
  scale = 1,
  width = 5.5,
  height = 6.5,
  units = c("in", "cm", "mm", "px"),
  dpi = 300)

# plot dyadic effect coeff for allogrooming, allofeeding, kinship ----
# get samples for dyadic effect coeff
tf <- 
  res2$samples$srm_model_samples$dyadic_coeffs %>% 
  as_tibble() %>% 
  rename(value= V1) %>% 
  mutate(label= "Allofeeding")

tg <- 
  res1$samples$srm_model_samples$dyadic_coeffs %>% 
  as_tibble() %>% 
  rename(value= V1) %>% 
  mutate(label= "Allogrooming")

tk <- 
  res3$samples$srm_model_samples$dyadic_coeffs %>% 
  as_tibble() %>% 
  rename(value= V1) %>% 
  mutate(label= "Kinship")

(dyad_coeff_plot <- 
    rbind(tf,tg,tk) %>% 
    mutate(label= fct_rev(label)) %>% 
    ggplot(aes(y = label, x = value, fill=label)) +
    stat_halfeye(.width = c(0.95), linewidth= 5, size=5)+
    geom_vline(xintercept = 0)+
    ylab("")+
    xlab("model coefficient")+
    theme_bw()+
    theme(legend.position= 'none', 
          axis.text.y = element_text(size = 12, vjust = 0.5, hjust = 0.5),
          axis.text.x = element_text(size = 12),
          axis.title.x = element_text(size = 12)
    ))

# save plot
ggsave(
  "results/summary.pdf",
  plot = dyad_coeff_plot,
  scale = 1,
  width = 5.5,
  height = 6.5,
  units = c("in", "cm", "mm", "px"),
  dpi = 300)

# wrangle data to measure effects for allogrooming and allofeeding RECEIVED 
# scale and store in list
# transpose the matrix to flip given and received
dyad = list(Feeding.received = scale(t(mf)),
            Grooming.received = scale(t(mg)))

# combine into df
dat = make_strand_data(
  outcome = list(calls=mcalls),
  individual_covariates = age,
  dyadic_covariates = dyad,
  outcome_mode = "poisson",
  link_mode = "log",
  check_standardization = T) 


# model: effect of allogrooming received -----
sfit4 <- 
  fit_social_relations_model(
    data=dat,
    focal_regression = ~ 1,
    target_regression = ~ 1,
    dyad_regression = ~ Grooming.received,
    mode="mcmc",
    stan_mcmc_parameters = list(
      chains = nchains,
      iter_warmup = chain_length,
      iter_sampling = warmup))

# summarize results
res4 <-  summarize_strand_results(sfit4)

# save diagnostics
diagnostics4 <- get_diagnostics(sfit4, "allogrooming received")

# save results
strand4 <- 
  res4$summary %>% 
  as_tibble() %>% 
  mutate(Model = "allogrooming received")

# plot results
(plot4a <- strand_caterpillar_plot(res4, normalized=FALSE,  only_slopes=F))
(plot4b <- strand_caterpillar_plot(res4, 
                                   submodels = "Dyadic effects", 
                                   normalized=T,  
                                   only_slopes=T)) 
(plot4c <- plot_posteriors(res4, "allogrooming received", "allogrooming received"))


# model effect of feeding received -----
sfit5 <- 
  fit_social_relations_model(
    data=dat,
    focal_regression = ~ 1,
    target_regression = ~ 1,
    dyad_regression = ~ Feeding.received,
    mode="mcmc",
    stan_mcmc_parameters = list(
      chains = nchains,
      iter_warmup = chain_length,
      iter_sampling = warmup))

# summarize results
res5 <-  summarize_strand_results(sfit5)

# get diagnostics
diagnostics5 <- get_diagnostics(sfit5, "allofeeding received")

# save results
strand5 <- 
  res5$summary %>% 
  as_tibble() %>% 
  mutate(Model = "allofeeding received")

# plot results
(plot5a <- strand_caterpillar_plot(res5, normalized=FALSE,  only_slopes=F))
(plot5b <- strand_caterpillar_plot(res5, 
                                   submodels = "Dyadic effects", 
                                   normalized=T,  
                                   only_slopes=T)) 
(plot5c <- plot_posteriors(res5, "allofeeding received", "allofeeding received"))


# plot all
library(patchwork)
plot1c+ plot2c+ plot3c+ plot4c+plot5c

### save all STRAND results
rbind(strand1, strand2, strand3, strand4, strand5) %>% 
  relocate(Model) %>% 
  write.csv(file= "results/STRAND_model_results.csv", row.names= FALSE)

### save all STRAND diagnostics
rbind(diagnostics1, diagnostics2, diagnostics3, diagnostics4, diagnostics5) %>% 
  write.csv(file= "results/STRAND_model_diagnostics.csv", row.names= FALSE)

### save workspace
timestamp <- substr(gsub(x=gsub(":","",Sys.time()), pattern=" ", replacement="_"), start=1, stop=15)
save.image(file= paste("model_fit_workspace_", timestamp, ".Rdata", sep=""))

