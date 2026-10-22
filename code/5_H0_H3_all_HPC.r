# Run H0-H3 models for algae, invertebrates and fish on HPC

rm(list = ls())

library(here)
library(tidyverse)
library(ape)
library(brms)
library(cmdstanr)

Sys.unsetenv("STANCFLAGS")

#---- 0. Define array jobs ---------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

chem_groups <- c(
  "Industrial_Chemical",
  "P_Arthropod_control",
  "P_Microbial_control",
  "P_Mixed_target",
  "P_Plant_control",
  "Pharmaceutical"
)

job_table <- expand_grid(
  taxon = c("algae", "invertebrate", "fish"),
  chem_group = chem_groups,
  model_type = c("H0", "H1", "H2", "H3")
)

root <- here("data/8_bayesian_model_fits")

#---- 1. Load data -----------------------------------------------------------
job_id <- as.integer(args[1])
stopifnot(length(args) == 1L, !is.na(job_id),
          job_id >= 1L, job_id <= nrow(job_table))

taxon_i <- job_table$taxon[job_id]
cg <- job_table$chem_group[job_id]
model_type_i <- job_table$model_type[job_id]
seed_i <- 2026L + job_id
model_name <- paste0(model_type_i, "_", taxon_i, "_", cg)

message("Taxon: ", taxon_i, " | Chemical group: ", cg,
        " | Model: ", model_type_i, " | Array task: ", job_id)

meta_df <- readRDS(here("data/6_combined_meta_data", paste0(taxon_i, "_meta_data.rds")))
tree <- readRDS(here("data/6_combined_meta_data", paste0(taxon_i, "_genus_tree.rds")))

#---- 2. Prepare data and standardise variables ------------------------------
genus_mat <- ape::vcv(tree, corr = TRUE)

df_subset <- meta_df %>%
  filter(chem_group == cg) %>%
  drop_na(log_rel_sensi, log_mass, CAS, genus, FFG)

shared_genera <- intersect(unique(as.character(df_subset$genus)), rownames(genus_mat))
genus_mat_subset <- genus_mat[shared_genera, shared_genera, drop = FALSE]

df_subset <- df_subset %>%
  filter(genus %in% shared_genera) %>%
  mutate(
    CAS = droplevels(factor(CAS)),
    genus_phylo = factor(genus, levels = shared_genera),
    FFG = droplevels(factor(FFG))
  )

scaling <- tibble(
  taxon = taxon_i, chem_group = cg, n = nrow(df_subset),
  mean_sensi = mean(df_subset$log_rel_sensi),
  sd_sensi = sd(df_subset$log_rel_sensi),
  mean_mass = mean(df_subset$log_mass),
  sd_mass = sd(df_subset$log_mass)
)

stopifnot(scaling$sd_sensi > 0, scaling$sd_mass > 0)

# Standardise within each taxon and chemical-use group
df_subset <- df_subset %>%
  mutate(
    log_rel_sensi_z = (log_rel_sensi - scaling$mean_sensi) / scaling$sd_sensi,
    log_mass_z = (log_mass - scaling$mean_mass) / scaling$sd_mass
  )

#---- 3. Define formula ------------------------------------------------------
if (model_type_i == "H0") {

  formula_i <- bf(
    log_rel_sensi_z ~ 0 +
      (1 | CAS) +
      (1 | gr(genus_phylo, cov = genus_mat))
  )

} else if (model_type_i == "H1") {

  formula_i <- bf(
    log_rel_sensi_z ~ 0 + log_mass_z +
      (1 | CAS) +
      (1 | gr(genus_phylo, cov = genus_mat))
  )

} else if (model_type_i == "H2") {

  formula_i <- bf(
    log_rel_sensi_z ~ 0 + FFG + log_mass_z +
      (1 | CAS) +
      (1 | gr(genus_phylo, cov = genus_mat))
  )

} else if (model_type_i == "H3") {

  formula_i <- bf(
    log_rel_sensi_z ~ 0 + FFG + FFG:log_mass_z +
      (1 | CAS) +
      (1 | gr(genus_phylo, cov = genus_mat))
  )
}

#---- 4. Define priors -------------------------------------------------------
priors_i <- c(
  prior(exponential(1), class = "sd", group = "CAS"),
  prior(normal(0, 1.5), class = "sd", group = "genus_phylo"),
  prior(student_t(3, 0, 0.5), class = "sigma")
)

if (model_type_i != "H0") {
  priors_i <- c(priors_i, prior(normal(0, 0.5), class = "b"))
}

#---- 5. Fit and save model --------------------------------------------------
fit_dir <- file.path(root, paste0(taxon_i, "_fits"))
dir.create(fit_dir, recursive = TRUE, showWarnings = FALSE)

message("Observations: ", nrow(df_subset), " | CAS: ", nlevels(df_subset$CAS),
        " | Genera: ", nlevels(df_subset$genus_phylo))

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

# Save scaling and settings with the fit
attr(fit_model, "standardisation") <- scaling
attr(fit_model, "sampling_settings") <- tibble(
  job_id = job_id, seed = seed_i,
  chains = 4L, iter = 8000L, warmup = 4000L, thin = 1L,
  backend = "cmdstanr", metric = "dense_e",
  adapt_delta = 0.99, max_treedepth = 15L,
  fit_hours = as.numeric(difftime(Sys.time(), fit_start, units = "hours"))
)
attr(fit_model, "session_info") <- capture.output(sessionInfo())

fit_path <- file.path(fit_dir, paste0(model_name, ".rds"))
saveRDS(fit_model, fit_path)
message("Saved fit to: ", fit_path)
