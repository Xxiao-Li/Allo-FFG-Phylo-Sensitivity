
# this script is to used to run models in HPC

rm(list=ls())

library(here)
library(tidyverse)  
library(ape)
library(brms)
library(rstan)
library(cmdstanr)

#----- 1. load dataset ------
spe_tax_phylo_BM_FFG_sensi_subCAS <- readRDS(here("data/6_combined/spe_tax_phylo_BM_FFG_sensi_subCAS.rds"))

phylo_tree <- readRDS(here("data/6_combined/phylo_tree.rds"))


#----- 2. prepare data for brms model ------
dat <- spe_tax_phylo_BM_FFG_sensi_subCAS %>%
  transmute(
    log_sensi = log_sensi,
    log_rel_sensi = log_rel_sensi,   
    log_mass  = log_mass,
    log_mass_c  = log_mass_c,
    group = factor(group),
    FG = factor(FG5),
    chem_group = factor(chem_group),
    duration  =  duration,
    CAS       =  as.character(CAS),
    latin     =  latin
  ) %>%
  filter(duration == "24")    %>%
  mutate(phylo = latin) 

splist <- unique(dat$latin)
tree_sub <- keep.tip(phylo_tree, intersect(phylo_tree$tip.label, splist))
tree_A <- ape::vcv.phylo(tree_sub, corr = TRUE)
dat$phylo <- factor(dat$phylo, levels = rownames(tree_A))

#----- 3. run brms models ------
## run brms model for a specific chemical group
CG <- unique(dat$chem_group)
chem_target <- CG[1] # change index to run model for different chemical groups
print(chem_target)

dat_sub <- dat %>%
  filter(chem_group == chem_target) %>% 
  mutate(FG = relevel(factor(FG), ref = "Cyanobacteria")) 
ref_FG <- levels(dat_sub$FG)[1]

pri <- c(
  prior(normal(0, 1),   class = "Intercept"),
  prior(normal(0, 0.5), class = "b"),
  prior(student_t(3, 0, 0.5), class = "sd"),
  prior(student_t(3, 0, 0.5), class = "sigma")
)

pri2 <- c( # for M4
  prior(normal(0, 1), class = "Intercept"),
  prior(normal(0, 0.5), class = "b"),
  prior(student_t(3, 0, 0.5), class = "sigma"))

thr <- threading(threads = 15)    # number of threads for a chain
Sys.setenv(STAN_NUM_THREADS = 15) 

chain <- 4
cor <- 4
ite <-10000
warm <- 5000

out_dir <- here("model_output", "bayesian_models")
dir.create(out_dir, showWarnings = FALSE)

# model M1: sensitivity ~ mass + (1 + CAS) + (1 | phylogeny)
fit_M1 <- brm(
  log_rel_sensi ~ log_mass_c + (1 | CAS) + (1 | gr(phylo, cov = tree_A)),
  data = dat_sub, data2 = list(tree_A = tree_A),
  family = gaussian(), prior = pri,
  backend = "cmdstanr", chains = chain, cores = cor, 
  threads = thr, seed = 1007,
  iter = ite, warmup = warm, 
  control = list(adapt_delta = 0.98, max_treedepth = 15))

fit_file_name <- paste0(out_dir, "/fit_M1_", chem_target, ".rds")
saveRDS(fit_M1, here(fit_file_name))


# model M2: sensitivity ~ mass + FG + (1 + CAS) + (1 | phylogeny)
fit_M2 <- brm(
  log_rel_sensi ~ log_mass_c + FG + (1 | CAS) + (1 | gr(phylo, cov = tree_A)),
  data = dat_sub, data2 = list(tree_A = tree_A),
  family = gaussian(), prior = pri,
  backend = "cmdstanr", chains = chain, cores = cor,
  threads = thr, seed = 1007,
  iter = ite, warmup = warm, 
  control = list(adapt_delta = 0.98, max_treedepth = 15))

fit_file_name <- paste0(out_dir, "/fit_M2_", chem_target, ".rds")
saveRDS(fit_M2, here(fit_file_name))


# model M3: sensitivity ~ mass * FG + (1 + CAS) + (1 | phylogeny)
fit_M3 <- brm(
  log_rel_sensi ~ log_mass_c * FG + (1 | CAS) + (1 | gr(phylo, cov = tree_A)),
  data = dat_sub, data2 = list(tree_A = tree_A),
  family = gaussian(), prior = pri,
  threads = thr, seed = 1007,
  iter = ite, warmup = warm, 
  control = list(adapt_delta = 0.98, max_treedepth = 15))

fit_file_name <- paste0(out_dir, "/fit_M3_", chem_target, ".rds")
saveRDS(fit_M3, here(fit_file_name))


# model M4: sensitivity ~ mass * FG + (1 + CAS)
fit_M4 <- brm(
  log_rel_sensi ~ log_mass_c * FG + (1 | CAS),
  data = dat_sub,
  family = gaussian(), prior = pri2,
  threads = thr, seed = 1007,
  iter = ite, warmup = warm, 
  control = list(adapt_delta = 0.98, max_treedepth = 15))

fit_file_name <- paste0(out_dir, "/fit_M4_", chem_target, ".rds")
saveRDS(fit_M4, here(fit_file_name))

