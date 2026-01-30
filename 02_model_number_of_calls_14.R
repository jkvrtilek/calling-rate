# analysis for 2024 paired caller-receiver recordings
# does contact calling rate depend on social relationship?
# all pairs of 8 bats were recorded with 2 synced mics

# Julia Vrtilek and Gerry Carter

# setup ----

# install packages (need to do this once)
# install.packages("cmdstanr", repos = c('https://stan-dev.r-universe.dev', getOption("repos")))
# library(cmdstanr)
# install_cmdstan()
# library(devtools)
# devtools::install_github('ctross/STRAND@beauty_in_the_dissonance')

# set working directory
setwd("/Users/jkvrtilek/Desktop/OSU/PhD/GitHub/calling-rate")

# load packages for data wrangling
library(tidyverse)
library(igraph)
library(broom)

# load packages for Bayesian models
library(cmdstanr)
library(STRAND)
library(brms)
library(performance)
library(tidybayes)

# packages for plotting
library(patchwork)
library(ggdist)


# load and wrangle vocal data ----
batcalls <- readRDS("vocal_data_2024-pairs_2025-12-08.RDS") %>% 
  group_by(caller, receiver) %>% 
  summarize(n.calls = n()) %>% 
  ungroup() %>% 
  # add quark who never called
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


# load and wrangle social data ----
# get allogrooming and allofeeding data
bat.donations <- 
  read.csv('OSU_2024_social_data.csv') %>% 
  mutate(Date = as.Date(Date, format = "%m/%d/%Y")) %>% 
  mutate(Actor = tolower(Actor),
         Receiver = tolower(Receiver)) %>% 
  filter(Actor %in% bats.used) %>% 
  filter(Receiver %in% bats.used) 

# make feeding matrix
fm <-
  bat.donations %>%
  filter(Behavior == "Mouthlicking") %>% 
  group_by(Actor, Receiver) %>%
  summarize(rate= sum(rate, na.rm=T)) %>%
  ungroup() %>% 
  mutate(rate= scale(rate)) %>% 
  graph_from_data_frame() %>% 
  as_adjacency_matrix(attr= 'rate', sparse=F)

# make grooming matrix
gm <-
  bat.donations %>%
  filter(Behavior == "Grooming") %>%
  group_by(Actor, Receiver) %>%
  summarize(rate= sum(rate, na.rm=T)) %>%
  ungroup() %>% 
  mutate(rate= scale(rate)) %>% 
  graph_from_data_frame() %>% 
  as_adjacency_matrix(attr= 'rate', sparse=F)

# make kinship matrix
raw.kin <- read.table('KING.txt') 
colnames(raw.kin) <- tolower(colnames(raw.kin))
rownames(raw.kin) <- tolower(rownames(raw.kin))

km <- 
  raw.kin %>% 
  select(all_of(bats.used)) %>% 
  filter(row.names(raw.kin) %in% bats.used) %>% 
  scale() %>% 
  as.matrix()


# set number and length of MCMC chains----- 
nchains <- 4
warmup <- 2000
chain_length <- 2000


# BRMS: fit GLMMs that assume actor and receiver have independent effects in calling -----------

# create function to convert matrix to dataframe
matrix_to_df <- function(m){
  data.frame(row=rownames(m)[row(m)], col=colnames(m)[col(m)], value=c(m), stringsAsFactors = F)
}

# get call counts
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
  gm %>% 
  matrix_to_df() %>% 
  filter(row!=col) %>% 
  group_by(row, col) %>% 
  summarize(groom= sum(value, na.rm=T)) %>% 
  rename(caller=row, receiver = col) %>% 
  ungroup()

# feed
feed2 <- 
  fm %>% 
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

# plot correlations
d2 %>% 
  pivot_longer(cols= groom:kinship) %>% 
  ggplot(aes(x=value, y= n.calls, color=name))+
  facet_wrap(~name, scales= "free", nrow=3)+
  geom_point(size=3, alpha=0.8)+
  geom_smooth(method= "lm")+
  xlab("standardized value")

## fit allogrooming GLMM ----
fit1 <-
  brm(n.calls ~ 
        groom +
        (1|caller)+
        (1|receiver)+
        (1|mm(caller,receiver)),
      data = d2, 
      family = negbinomial(),
      seed = 123,
      cores = nchains,
      chains = nchains,
      iter = chain_length + warmup,
      warmup = warmup,
      control = list(adapt_delta = 0.95))

# compare to model without dyadic random effects
fit1b <-
  brm(n.calls ~ 
        groom +
        (1|caller)+
        (1|receiver),
      data = d2, 
      family = negbinomial(),
      seed = 123,
      cores = nchains,
      chains = nchains,
      iter = chain_length + warmup,
      warmup = warmup,
      control = list(adapt_delta = 0.95))

# leave one out cross-validation
loo(fit1, fit1b)
# models are basically equivalent

## fit allofeeding GLMM ----
fit2 <-
  brm(n.calls ~ 
        feed +
        (1|caller)+
        (1|receiver)+
        (1|mm(caller,receiver)),
      data = d2, 
      family = negbinomial(),
      seed = 123,
      cores = nchains,
      chains = nchains,
      iter = chain_length + warmup,
      warmup = warmup,
      control = list(adapt_delta = 0.95))

# compare to model without dyadic random effects
fit2b <-
  brm(n.calls ~ 
        feed +
        (1|caller)+
        (1|receiver),
      data = d2, 
      family = negbinomial(),
      seed = 123,
      cores = nchains,
      chains = nchains,
      iter = chain_length + warmup,
      warmup = warmup,
      control = list(adapt_delta = 0.95))

# leave-one-out CV
loo(fit2, fit2b)
# models are equivalent

## fit kinship GLMM ----
fit3 <-
  brm(n.calls ~ 
        kinship +
        (1|caller)+
        (1|receiver)+
        (1|mm(caller,receiver)),
      data = d2, 
      family = negbinomial(),
      seed = 123,
      cores = nchains,
      chains = nchains,
      iter = chain_length+ warmup,
      warmup = warmup,
      control = list(adapt_delta = 0.95))

# compare to model without dyadic random effects
fit3b <-
  brm(n.calls ~ 
        kinship +
        (1|caller)+
        (1|receiver),
      data = d2, 
      family = negbinomial(),
      seed = 123,
      cores = nchains,
      chains = nchains,
      iter = chain_length + warmup,
      warmup = warmup,
      control = list(adapt_delta = 0.95))

# leave-one-out CV
loo(fit3, fit3b)
# models are equivalent

# posterior predictive checks
pp_check(fit1, ndraws=100)+xlim(0, 1000)
pp_check(fit2, ndraws=100)+xlim(0, 1000)
pp_check(fit3, ndraws=100)+xlim(0, 1000)
pp_check(fit1b, ndraws=100)+xlim(0, 1000)
pp_check(fit2b, ndraws=100)+xlim(0, 1000)
pp_check(fit3b, ndraws=100)+xlim(0, 1000)

# get posteriors from simpler model with independent caller and receiver effects
b1 <- 
  fit1b %>% 
  spread_draws(b_groom) %>% 
  mutate(model = "Allogrooming") %>% 
  pivot_longer(b_groom, names_to = 'term', values_to= 'coeff') 

b2 <- 
  fit2b %>% 
  spread_draws(b_feed) %>% 
  mutate(model = "Allofeeding") %>% 
  pivot_longer(b_feed, names_to = 'term', values_to= 'coeff') 

b3 <- 
  fit3b %>% 
  spread_draws(b_kinship) %>% 
  mutate(model = "Kinship") %>% 
  pivot_longer(b_kinship, names_to = 'term', values_to= 'coeff') 

## plot brms models-------
(brms.plot <- 
    rbind(b1, b2, b3) %>% 
    mutate(label= "independent fixed effects on calling rate") %>% 
    mutate(model= fct_rev(model)) %>%  
    ggplot(aes(y = model, x = coeff, fill=model)) +
    facet_wrap(~label)+
    stat_halfeye(.width = c(0.8, 0.95))+
    geom_vline(xintercept = 0, linetype="dashed", color= "grey")+
    ylab("")+
    xlab("estimated effect on call count (scaled coefficient)")+
    theme_bw()+
    coord_cartesian(xlim = c(-2,2))+
    theme(legend.position= 'none', 
          axis.text.y = element_text(size = 10, vjust = 0.5, hjust = 0.5),
          axis.text.x = element_text(size = 10),
          axis.title.x = element_text(size = 10)
    ))

# save plot
ggsave(
  "results/brms_models_plot.pdf",
  plot = brms.plot,
  scale = 1,
  width = 5.5,
  height = 6.5,
  units = c("in", "cm", "mm", "px"),
  dpi = 300)


## get and save brms model results -----
brm1 <- summary(fit1b)$fixed %>% rownames_to_column("estimate") %>% mutate(model= "Allogrooming") %>% relocate(model)
brm2 <- summary(fit2b)$fixed %>% rownames_to_column("estimate") %>% mutate(model= "Allofeeding") %>% relocate(model)
brm3 <- summary(fit3b)$fixed %>% rownames_to_column("estimate") %>% mutate(model= "Kinship") %>% relocate(model)
(brm_results <- 
    rbind(brm1,brm2,brm3) %>% 
    filter(estimate!= "Intercept") %>% 
    mutate(Type= "Caller-Receiver Model") %>% 
    select(Type, Predictor= model, Estimate, Low_95= `l-95% CI`, High_95= `u-95% CI`, Rhat, Bulk_ESS, Tail_ESS)) 

brm_results %>% write.csv("results/brms_model_results.csv")

# STRAND: fit social relations model ---------

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

### create function to plot posterior distributions from STRAND
# note: block parameters 1 = intercept
plot_posteriors <- function(results= results, 
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
      mutate(type = case_when(
        label=='dyadic_effects_rho' ~ "correlations between random effects (calling and receiving)",
        label=='focal-target_effects_rho' ~ "correlations between random effects (calling and receiving)",
        label=='dyadic_effect_coeff' ~ "adjusted fixed effect on calling rate (scaled coefficient)",
        TRUE ~ "random effects on calling rate (standard deviations)")) %>% 
      mutate(type= fct_relevel(type, "random effects on calling rate (standard deviations)", 
                               "correlations between random effects (calling and receiving)")) %>% 
      mutate(label= case_when(
        label== 'focal_effects_sd' ~ "Caller",
        label=='target_effects_sd' ~ "Receiver",
        label=='dyadic_effect_sd' ~ "Dyad",
        label=='dyadic_effect_coeff' ~ effect,
        label=='dyadic_effects_rho' ~ "Dyad-level",
        label=='focal-target_effects_rho' ~ "Group-level")) %>% 
      mutate(label= fct_relevel(label, "Caller", "Receiver", "Dyad")) %>% 
      mutate(label= fct_rev(label)) %>% 
      ggplot(aes(y = label, x = value, fill=label)) +
      facet_wrap(~type, scales= "free", ncol=1)+
      stat_halfeye(.width = c(0.8, 0.95))+
      ylab("")+
      xlab("estimate in social relations model")+
      geom_vline(xintercept = 0, linetype="dashed", color= "grey")+
      theme_bw()+
      coord_cartesian(xlim = c(-1,5))+
      theme(legend.position= 'none', 
          axis.text.y = element_text(size = 10, vjust = 0.5, hjust = 0.5),
          axis.text.x = element_text(size = 10),
          axis.title.x = element_text(size = 10)))
  
  return(tplot)
  
}

## make the STRAND data structure ----
dyad = list(Kinship = km,
            Feeding = fm,
            Grooming = gm)

# combine into df
dat = make_strand_data(
  outcome = list(calls=mcalls),
  dyadic_covariates = dyad,
  outcome_mode = "poisson",
  link_mode = "log",
  check_data_organization = TRUE,
  check_standardization = TRUE) 

## model effect of allogrooming given on call count -----
sfit1 <- 
  fit_social_relations_model(
    data=dat,
    focal_regression = ~ 1,
    target_regression = ~ 1,
    dyad_regression = ~ Grooming,
    mode="mcmc",
    mcmc_parameters = list(
      seed = 123,
      chains = nchains,
      parallel_chains = nchains,
      iter_warmup = chain_length,
      iter_sampling = warmup, 
      cores= nchains))

# save diagnostics
(diagnostics1 <- get_diagnostics(sfit1, "allogrooming"))

# get results
res1 <-  summarize_strand_results(sfit1)

# save results
(strand1 <- 
  res1$summary %>% 
  as_tibble() %>% 
  filter(Variable != "error sd") %>% 
  mutate(Model = "allogrooming"))

### plot allogrooming results-----
(strand_caterpillar_plot(res1, normalized=F))
(s.plot1 <- plot_posteriors(res1, effect= "Allogrooming"))

# save plot
ggsave(
  "results/STRAND_model_allogrooming.pdf",
  plot = s.plot1,
  scale = 1,
  width = 4.5,
  height = 6,
  units = c("in", "cm", "mm", "px"),
  dpi = 300)

## model effect of allofeeding given on call count -----
sfit2 <- 
  fit_social_relations_model(
    data=dat,
    focal_regression = ~ 1,
    target_regression = ~ 1,
    dyad_regression = ~ Feeding,
    mode="mcmc",
    mcmc_parameters = list(
      seed = 123,
      chains = nchains,
      parallel_chains = nchains,
      iter_warmup = chain_length,
      iter_sampling = warmup, 
      cores= nchains))

# save diagnostics
(diagnostics2 <- get_diagnostics(sfit2, "allofeeding"))

# get results
res2 <-  summarize_strand_results(sfit2)

# save results
(strand2 <- 
    res2$summary %>% 
    as_tibble() %>% 
    filter(Variable != "error sd") %>% 
    mutate(Model = "allofeeding"))

### plot allofeeding results--------
(strand_caterpillar_plot(res2, normalized=F))
res2$summary
(s.plot2 <- plot_posteriors(res2, "Allofeeding"))

# save plot
ggsave(
  "results/STRAND_model_allofeeding.pdf",
  plot = s.plot2,
  scale = 1,
  width = 4.5,
  height = 6,
  units = c("in", "cm", "mm", "px"),
  dpi = 300)


## model effect of kinship on call count ----
sfit3 <- 
  fit_social_relations_model(
    data=dat,
    focal_regression = ~ 1,
    target_regression = ~ 1,
    dyad_regression = ~ Kinship,
    mode="mcmc",
    mcmc_parameters = list(
      seed = 123,
      chains = nchains,
      parallel_chains = nchains,
      iter_warmup = chain_length,
      iter_sampling = warmup, 
      cores= nchains))

# save diagnostics
(diagnostics3 <- get_diagnostics(sfit3, "kinship"))

# get results
res3 <-  summarize_strand_results(sfit3)

# save results
(strand3 <- 
    res3$summary %>% 
    as_tibble() %>% 
    filter(Variable != "error sd") %>% 
    mutate(Model = "kinship"))

### plot kinship results--------
strand_caterpillar_plot(res3, normalized=F)
(s.plot3 <- plot_posteriors(res3,"Kinship"))

# save plot
ggsave(
  "results/STRAND_model_kinship.pdf",
  plot = s.plot3,
  scale = 1,
  width = 5.5,
  height = 6.5,
  units = c("in", "cm", "mm", "px"),
  dpi = 300)


## get samples for dyadic effect coeff for plotting ----
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

## plot STRAND models ----
(strands.plot <- 
    rbind(tf,tg,tk) %>% 
    mutate(label= fct_rev(label)) %>% 
    ggplot(aes(y = label, x = value, fill=label)) +
    stat_halfeye(.width = c(0.8, 0.95))+
    geom_vline(xintercept = 0, linetype="dashed", color= "grey")+
    ylab("")+
   xlab("estimated effect on call count (model coefficient)")+
    theme_bw()+
    coord_cartesian(xlim = c(-2,2))+
    theme(legend.position= 'none', 
          axis.text.y = element_text(size = 10, vjust = 0.5, hjust = 0.5),
          axis.text.x = element_text(size = 10),
          axis.title.x = element_text(size = 10)
    ))

# save plot
ggsave(
  "results/STRAND_model_summary.pdf",
  plot = strands.plot,
  scale = 1,
  width = 5.5,
  height = 6.5,
  units = c("in", "cm", "mm", "px"),
  dpi = 300)


### combine plots -----

# create main results plot for paper
(t1 <- 
   plot_posteriors(res1, "Allogrooming"))
(fig1 <- brms.plot + t1 + plot_layout(ncol=2, widths= c(1,1.2)))

# save plot
ggsave(
  "results/Fig1.pdf",
  plot = fig1,
  scale = 1,
  width = 9.3,
  height = 4,
  units = c("in", "cm", "mm", "px"),
  dpi = 300)

# compare vocal reciprocity with allogrooming reciprocity-------
# model effect of kinship on allogrooming

# make grooming matrix
gm2 <-
  bat.donations %>%
  filter(Behavior == "Grooming") %>%
  group_by(Actor, Receiver) %>%
  summarize(rate= sum(rate, na.rm=T)) %>%
  ungroup() %>% 
  graph_from_data_frame() %>% 
  as_adjacency_matrix(attr= 'rate', sparse=F)

# get covariates
dyad2 = list(Kinship = km)

# combine into df
dat2 = make_strand_data(
  outcome = list(calls=gm2),
  dyadic_covariates = dyad2,
  outcome_mode = "poisson",
  link_mode = "log",
  check_data_organization = TRUE,
  check_standardization = TRUE) 

sfit4 <- 
  fit_social_relations_model(
    data=dat2,
    focal_regression = ~ 1,
    target_regression = ~ 1,
    dyad_regression = ~ Kinship,
    mode="mcmc",
    mcmc_parameters = list(
      seed = 123,
      chains = nchains,
      parallel_chains = nchains,
      iter_warmup = chain_length,
      iter_sampling = warmup, 
      cores= nchains))

# save diagnostics
(diagnostics4 <- get_diagnostics(sfit4, "kinship"))

# get results
res4 <-  summarize_strand_results(sfit4)

# save results
(strand4 <- 
    res4$summary %>% 
    as_tibble() %>% 
    filter(Variable != "error sd"))

# plot results
plot_posteriors2 <- function(results= results, 
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
      mutate(type = case_when(
        label=='dyadic_effects_rho' ~ "correlations between random effects (grooming given and received)",
        label=='focal-target_effects_rho' ~ "correlations between random effects (grooming given and received)",
        label=='dyadic_effect_coeff' ~ "adjusted fixed effect on allogrooming rate (scaled coefficient)",
        TRUE ~ "random effects on allogrooming rate (standard deviations)")) %>% 
      mutate(type= fct_relevel(type, "random effects on allogrooming rate (standard deviations)", "correlations between random effects (grooming given and received)")) %>% 
      mutate(label= case_when(
        label== 'focal_effects_sd' ~ "Groomer",
        label=='target_effects_sd' ~ "Receiver",
        label=='dyadic_effect_sd' ~ "Dyad",
        label=='dyadic_effect_coeff' ~ effect,
        label=='dyadic_effects_rho' ~ "Dyad-level",
        label=='focal-target_effects_rho' ~ "Group-level")) %>% 
      mutate(label= fct_relevel(label, "Groomer", "Receiver", "Dyad")) %>% 
      mutate(label= fct_rev(label)) %>% 
      ggplot(aes(y = label, x = value, fill=label)) +
      facet_wrap(~type, scales= "free", ncol=1)+
      stat_halfeye(.width = c(0.8, 0.95))+
      ylab("")+
      xlab("estimate in social relations model")+
      geom_vline(xintercept = 0, linetype="dashed", color= "grey")+
      coord_cartesian(xlim = c(-1,5))+
      theme_bw()+
      theme(legend.position= 'none', 
            axis.text.y = element_text(size = 10, vjust = 0.5, hjust = 0.5),
            axis.text.x = element_text(size = 10),
            axis.title.x = element_text(size = 10)))
  
  return(tplot)
  
}

strand_caterpillar_plot(res4, normalized=F)
(s.plot4 <- plot_posteriors2(res4,"Kinship"))

(figS3 <- s.plot3+ s.plot4)

# save plot
ggsave(
  "results/FigS3.pdf",
  plot = fig2 ,
  scale = 1,
  width = 10,
  height = 7,
  units = c("in", "cm", "mm", "px"),
  dpi = 300)

# get samples for saving STRAND models credible intervals
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

# save STRAND summary results-----------

# save STRAND results
rbind(strand1, strand2, strand3) %>% 
  relocate(Model) %>% 
  write.csv(file= "results/STRAND_model_HPDI.csv", row.names= FALSE)

# save all STRAND diagnostics
rbind(diagnostics1, diagnostics2, diagnostics3) %>% 
  write.csv(file= "results/STRAND_model_diagnostics.csv", row.names= FALSE)

# get predictors with 95% credible intervals ----

# get brm predictors
brm_results2 <- 
  brm_results %>% 
  select(Type, Predictor, Estimate, Low_95, High_95)
  
# save all predictors
rbind(tf, tg, tk) %>% 
  group_by(label) %>% 
    summarize(
      Estimate = mean (value),
      Low_95  = quantile(value, probs = 0.025), # Calculate 2.5% percentile of mpg
      High_95  = quantile(value, probs = 0.975)  # Calculate 97.5th percentile
    ) %>% 
  rename(Predictor= label) %>% 
  mutate(Type= "Social Relations Model") %>% 
  rbind(brm_results2) %>% 
  select(Type, Predictor, Estimate, Low_95, High_95) %>% 
  mutate(IRR = exp(Estimate), 
         IRR_low = exp(Low_95), 
         IRR_high = exp(High_95)) %>% 
  arrange(Type, Predictor) %>% 
  mutate(across(where(is.numeric), round, digits = 2)) %>% 
  write.csv(file= "results/predictors_results.csv", row.names= FALSE)

### save workspace
timestamp <- substr(gsub(x=gsub(":","",Sys.time()), pattern=" ", replacement="_"), start=1, stop=15)
save.image(file= paste("model_fit_workspace_", timestamp, ".Rdata", sep=""))
