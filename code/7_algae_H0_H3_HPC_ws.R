
# this script is to run the algae H0-H3 models on HPC with solubility included as an additional fixed effect.

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

#---- 1. Load data -------------------------------------------------------------

algae_df <- readRDS(here("data/6_combined_meta_data/algae_meta_data.rds"))
algae_tree <- readRDS(here("data/6_combined_meta_data/algae_genus_tree.rds"))

#---- 2. Define jobs -----------------------------------------------------------

chem_groups_algae <- c(
  "Industrial_Chemical",
  "P_Arthropod_control",
  "P_Microbial_control",
  "P_Mixed_target",
  "P_Plant_control",
  "Pharmaceutical"
)

model_types <- c(
  "H0",
  "H1",
  "H2",
  "H3"
)

job_table <- expand_grid(
  chem_group = chem_groups_algae,
  model_type = model_types
) %>%
  mutate(
    job_id = row_number(),
    tax = "algae",
    family = "gaussian",
    adapt_delta = 0.95,
    max_treedepth = 15,
    iter = 16000L,
    warmup = 6000L,
    chains = if_else(chem_group == "P_Arthropod_control", 6L, 4L),
    cores = if_else(chem_group == "P_Arthropod_control", 6L, 4L)
  )

job <- job_table %>%
  filter(job_id == !!job_id)

if (nrow(job) != 1) {
  stop("job_id does not match exactly one row in job_table.")
}

tax <- job$tax
cg <- job$chem_group
model_type_i <- job$model_type
family_name_i <- job$family
adapt_delta_i <- job$adapt_delta
max_treedepth_i <- job$max_treedepth
iter_i <- job$iter
warmup_i <- job$warmup
chains_i <- job$chains
cores_i <- job$cores

#---- 3. Prepare data ----------------------------------------------------------

prepare_model_data <- function(df) {
  df %>%
    filter(
      !is.na(log_rel_sensi),
      !is.na(log_mass),
      !is.na(CAS),
      !is.na(genus),
      !is.na(chem_group),
      !is.na(FFG),
      !is.na(log_molar_solubility)
    ) %>%
    mutate(
      CAS = as.character(CAS),
      CAS = str_remove(CAS, "^0+"),
      CAS = factor(CAS),
      genus_phylo = factor(genus),
      chem_group = factor(chem_group),
      FFG = factor(FFG)
    )
}

algae_df <- prepare_model_data(algae_df)

#---- 4. Phylogenetic covariance matrix ----------------------------------------

make_genus_cov <- function(tree) {
  ape::vcv(tree, corr = TRUE)
}

genus_mat_algae <- make_genus_cov(algae_tree)

match_data_and_genus_mat <- function(df, genus_mat) {
  
  shared_genera <- intersect(
    as.character(unique(df$genus_phylo)),
    rownames(genus_mat)
  )
  
  df2 <- df %>%
    filter(genus_phylo %in% shared_genera) %>%
    mutate(
      genus_phylo = factor(genus_phylo, levels = shared_genera),
      CAS = droplevels(factor(CAS)),
      FFG = droplevels(factor(FFG)),
      chem_group = droplevels(factor(chem_group))
    )
  
  mat2 <- genus_mat[shared_genera, shared_genera, drop = FALSE]
  
  list(
    data = df2,
    genus_mat = mat2
  )
}

#---- 5. Select subset ---------------------------------------------------------

df_subset_raw <- algae_df %>%
  filter(chem_group == cg) %>%
  droplevels()

matched <- match_data_and_genus_mat(
  df = df_subset_raw,
  genus_mat = genus_mat_algae
)

df_subset <- matched$data %>%
  mutate(
    log_rel_sensi_c =
      log_rel_sensi - mean(log_rel_sensi, na.rm = TRUE),
    
    log_mass_c =
      log_mass - mean(log_mass, na.rm = TRUE),
    
    log_molar_solubility_c =
      log_molar_solubility - mean(log_molar_solubility, na.rm = TRUE)
  ) %>%
  droplevels()

genus_mat_subset <- matched$genus_mat

#---- 5b. Basic model-data summaries ------------------------------------------

n_obs <- nrow(df_subset)
n_CAS <- n_distinct(df_subset$CAS)
n_genus <- n_distinct(df_subset$genus_phylo)
n_FFG <- n_distinct(df_subset$FFG)

mean_y <- mean(df_subset$log_rel_sensi_c, na.rm = TRUE)
mean_x <- mean(df_subset$log_mass_c, na.rm = TRUE)
mean_sol <- mean(df_subset$log_molar_solubility_c, na.rm = TRUE)

sd_y <- sd(df_subset$log_rel_sensi_c, na.rm = TRUE)
sd_x <- sd(df_subset$log_mass_c, na.rm = TRUE)
sd_sol <- sd(df_subset$log_molar_solubility_c, na.rm = TRUE)

raw_mean_y <- mean(df_subset$log_rel_sensi, na.rm = TRUE)
raw_mean_x <- mean(df_subset$log_mass, na.rm = TRUE)
raw_mean_sol <- mean(df_subset$log_molar_solubility, na.rm = TRUE)

message("Taxon: ", tax)
message("Chemical group: ", cg)
message("Model type: ", model_type_i)
message("n_obs: ", n_obs)
message("n_CAS: ", n_CAS)
message("n_genus: ", n_genus)
message("n_FFG: ", n_FFG)
message("raw mean log_rel_sensi: ", raw_mean_y)
message("raw mean log_mass: ", raw_mean_x)
message("raw mean log_molar_solubility: ", raw_mean_sol)
message("sd log_rel_sensi_c: ", sd_y)
message("sd log_mass_c: ", sd_x)
message("sd log_molar_solubility_c: ", sd_sol)

if (n_obs == 0) {
  stop("No observations available after filtering.")
}

if (is.na(sd_y) || sd_y == 0) {
  stop("Response SD is zero or NA.")
}

if (is.na(sd_sol) || sd_sol == 0) {
  stop("Solubility SD is zero or NA in this subset.")
}

#---- 6. Define formula --------------------------------------------------------

if (model_type_i == "H0") {
  
  formula_i <- bf(
    log_rel_sensi_c ~ 0 + log_molar_solubility_c +
      (1 | CAS) +
      (1 | gr(genus_phylo, cov = genus_mat))
  )
  
} else if (model_type_i == "H1") {
  
  formula_i <- bf(
    log_rel_sensi_c ~ 0 + log_mass_c + log_molar_solubility_c +
      (1 | CAS) +
      (1 | gr(genus_phylo, cov = genus_mat))
  )
  
} else if (model_type_i == "H2") {
  
  formula_i <- bf(
    log_rel_sensi_c ~ 0 + FFG + log_mass_c + log_molar_solubility_c +
      (1 | CAS) +
      (1 | gr(genus_phylo, cov = genus_mat))
  )
  
} else if (model_type_i == "H3") {
  
  formula_i <- bf(
    log_rel_sensi_c ~ 0 + FFG + FFG:log_mass_c + log_molar_solubility_c +
      (1 | CAS) +
      (1 | gr(genus_phylo, cov = genus_mat))
  )
  
} else {
  stop("Unknown model_type: ", model_type_i)
}

#---- 7. Define priors ---------------------------------------------------------

cas_prior_rate_i <- signif(1 / sd_y, 6)
cas_prior_mean_i <- 1 / cas_prior_rate_i

priors_i <- c(
  prior_string(
    paste0("exponential(", cas_prior_rate_i, ")"),
    class = "sd",
    group = "CAS"
  ),
  prior(exponential(2), class = "sd", group = "genus_phylo"),
  prior(student_t(3, 0, 0.3), class = "sigma"),
  prior(normal(0, 0.3), class = "b")
)

family_i <- gaussian()

#---- 8. Fit model -------------------------------------------------------------

fit_model <- brm(
  formula = formula_i,
  data = df_subset,
  data2 = list(genus_mat = genus_mat_subset),
  family = family_i,
  prior = priors_i,
  chains = chains_i,
  cores = cores_i,
  iter = iter_i,
  warmup = warmup_i,
  backend = "cmdstanr",
  control = list(
    adapt_delta = adapt_delta_i,
    max_treedepth = max_treedepth_i
  ),
  seed = 123 + job_id
)

#---- 9. Save model fit --------------------------------------------------------

base_dir <- here("results/algae_H0_H3_with_solubility")
out_dir <- file.path(base_dir, "fits")
job_dir <- file.path(base_dir, "job_info")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(job_dir, recursive = TRUE, showWarnings = FALSE)

safe_cg <- gsub("[^A-Za-z0-9_]+", "_", cg)

fit_name <- paste0(
  model_type_i,
  "_algae_",
  safe_cg,
  "_with_solubility"
)

fit_path <- file.path(out_dir, paste0(fit_name, ".rds"))

saveRDS(fit_model, fit_path)

message("Saved fit to: ", fit_path)

#---- 10. Save job information -------------------------------------------------

job_info <- list(
  job = job,
  formula = formula_i,
  priors = priors_i,
  n_obs = n_obs,
  n_CAS = n_CAS,
  n_genus = n_genus,
  n_FFG = n_FFG,
  mean_y = mean_y,
  mean_x = mean_x,
  mean_sol = mean_sol,
  sd_y = sd_y,
  sd_x = sd_x,
  sd_sol = sd_sol,
  raw_mean_y = raw_mean_y,
  raw_mean_x = raw_mean_x,
  raw_mean_sol = raw_mean_sol,
  cas_prior_rate = cas_prior_rate_i,
  cas_prior_mean = cas_prior_mean_i
)

job_info_path <- file.path(job_dir, paste0(fit_name, "_job_info.rds"))

saveRDS(job_info, job_info_path)

message("Saved job info to: ", job_info_path)

