#----------------------------------------------------------------------------
# Test H3 after randomly removing 5%, 10%, or 15% of sensitivity records.
# Main-analysis H3 fits provide the 0% baseline.
#----------------------------------------------------------------------------

rm(list = ls())

library(here)
library(tidyverse)
library(ape)
library(brms)
library(cmdstanr)

Sys.unsetenv("STANCFLAGS")

#---- 0. Read array task ID ----------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
job_id <- as.integer(args[1])

#---- 1. Load data and trees ---------------------------------------------------
dat_list <- list(
  algae = readRDS(here("data/6_combined_meta_data/algae_meta_data.rds")),
  invertebrate = readRDS(here("data/6_combined_meta_data/invertebrate_meta_data.rds")),
  fish = readRDS(here("data/6_combined_meta_data/fish_meta_data.rds"))
)

tree_list <- list(
  algae = readRDS(here("data/6_combined_meta_data/algae_genus_tree.rds")),
  invertebrate = readRDS(here("data/6_combined_meta_data/invertebrate_genus_tree.rds")),
  fish = readRDS(here("data/6_combined_meta_data/fish_genus_tree.rds"))
)

#---- 2. Define jobs -----------------------------------------------------------
chem_groups <- c(
  "Industrial_Chemical",
  "P_Arthropod_control",
  "P_Microbial_control",
  "P_Mixed_target",
  "P_Plant_control",
  "Pharmaceutical"
)

drop_rates <- c(0.05, 0.10, 0.15)
n_repeats <- 5

available_groups <- map_dfr(names(dat_list), function(tax) {
  tibble(
    tax = tax,
    chem_group = intersect(
      chem_groups, unique(as.character(dat_list[[tax]]$chem_group))
    )
  )
}) %>%
  mutate(group_id = row_number())

job_table <- available_groups %>%
  crossing(drop_rate = drop_rates, repeat_id = seq_len(n_repeats))

stopifnot(length(args) == 1L, !is.na(job_id),
          job_id >= 1L, job_id <= nrow(job_table))

job <- job_table[job_id, ]

tax <- job$tax
cg <- job$chem_group
drop_rate_i <- job$drop_rate
repeat_id_i <- job$repeat_id
group_id_i <- job$group_id
seed_i <- 2026L + job_id
exclusion_seed <- 2026L + 100L * group_id_i + repeat_id_i

message("Job: ", job_id, " / ", nrow(job_table))
message("Taxonomic group: ", tax)
message("Chemical group: ", cg)
message("Drop rate: ", drop_rate_i, "; repeat: ", repeat_id_i)

#---- 3. Prepare data and phylogenetic covariance ------------------------------
genus_mat <- ape::vcv(tree_list[[tax]], corr = TRUE)

df_full <- dat_list[[tax]] %>%
  filter(chem_group == cg) %>%
  drop_na(log_rel_sensi, log_mass, CAS, genus, FFG)

shared_genera <- intersect(unique(as.character(df_full$genus)), rownames(genus_mat))

df_full <- df_full %>%
  filter(genus %in% shared_genera) %>%
  mutate(
    CAS = droplevels(factor(CAS)),
    genus_phylo = factor(genus, levels = shared_genera),
    FFG = droplevels(factor(FFG))
  )

# Keep the full-data means and SDs across exclusion rates.
scaling <- tibble(
  taxon = tax, chem_group = cg, n = nrow(df_full),
  mean_sensi = mean(df_full$log_rel_sensi),
  sd_sensi = sd(df_full$log_rel_sensi),
  mean_mass = mean(df_full$log_mass),
  sd_mass = sd(df_full$log_mass)
)

stopifnot(scaling$sd_sensi > 0, scaling$sd_mass > 0)

#---- 4. Randomly exclude records ---------------------------------------------
# Use the same random ordering across drop rates within each repeat.
set.seed(exclusion_seed)

df_subset <- df_full %>%
  group_by(CAS) %>%
  mutate(drop_order = sample.int(n())) %>%
  filter(drop_order > floor(n() * drop_rate_i)) %>%
  ungroup() %>%
  select(-drop_order) %>%
  mutate(
    log_rel_sensi_z = (log_rel_sensi - scaling$mean_sensi) / scaling$sd_sensi,
    log_mass_z = (log_mass - scaling$mean_mass) / scaling$sd_mass
  ) %>%
  droplevels()

used_genera <- intersect(
  unique(as.character(df_subset$genus_phylo)), shared_genera
)

df_subset$genus_phylo <- factor(df_subset$genus_phylo, levels = used_genera)
genus_mat_subset <- genus_mat[used_genera, used_genera, drop = FALSE]

#---- 5. Define formula --------------------------------------------------------
formula_i <- bf(
  log_rel_sensi_z ~ 0 + FFG + FFG:log_mass_z +
    (1 | CAS) +
    (1 | gr(genus_phylo, cov = genus_mat))
)

#---- 6. Define priors ---------------------------------------------------------
priors_i <- c(
  prior(exponential(1), class = "sd", group = "CAS"),
  prior(normal(0, 1.5), class = "sd", group = "genus_phylo"),
  prior(student_t(3, 0, 0.5), class = "sigma"),
  prior(normal(0, 0.5), class = "b")
)

#---- 7. Fit model -------------------------------------------------------------
fit_start <- Sys.time()
fit_model <- brm(
  formula = formula_i,
  data = df_subset,
  data2 = list(genus_mat = genus_mat_subset),
  family = gaussian(),
  prior = priors_i,
  chains = 4,
  cores = 4,
  iter = 8000,
  warmup = 4000,
  backend = "cmdstanr",
  metric = "dense_e",
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  output_dir = NULL,
  seed = seed_i,
  refresh = 100
)

#---- 8. Save model fit --------------------------------------------------------
attr(fit_model, "standardisation") <- scaling
attr(fit_model, "sampling_settings") <- tibble(
  job_id = job_id, seed = seed_i,
  chains = 4L, iter = 8000L, warmup = 4000L, thin = 1L,
  backend = "cmdstanr", metric = "dense_e",
  adapt_delta = 0.99, max_treedepth = 15L,
  fit_hours = as.numeric(difftime(Sys.time(), fit_start, units = "hours"))
)
attr(fit_model, "exclusion") <- tibble(
  drop_rate = drop_rate_i, repeat_id = repeat_id_i, seed = exclusion_seed,
  n_full = nrow(df_full), n_retained = nrow(df_subset),
  n_dropped = nrow(df_full) - nrow(df_subset)
)
attr(fit_model, "session_info") <- capture.output(sessionInfo())

out_dir <- here("data/8_bayesian_model_fits/random_exclusion_test")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fit_name <- paste(
  "H3", tax, cg, paste0("drop", round(drop_rate_i * 100)),
  paste0("rep", repeat_id_i), sep = "_"
)

fit_path <- file.path(out_dir, paste0(fit_name, ".rds"))
saveRDS(fit_model, fit_path)
message("Saved fit to: ", fit_path)
