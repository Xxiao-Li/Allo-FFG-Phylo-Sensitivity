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

invertebrate_df <- readRDS(here("data/6_combined_meta_data/invertebrate_meta_data.rds"))
invertebrate_tree <- readRDS(here("data/6_combined_meta_data/invertebrate_genus_tree.rds"))

#---- 2. Define jobs -----------------------------------------------------------

chem_groups_invertebrate <- c(
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
  chem_group = chem_groups_invertebrate,
  model_type = model_types
) %>%
  mutate(
    job_id = row_number(),
    tax = "invertebrate",
    prior_type = "gaussian_CAS_scaled_phylo2_fixed03",
    family = "gaussian",
    adapt_delta = 0.95,
    max_treedepth = 15L,
    iter = 16000L,
    warmup = 6000L,
    chains = 4L,
    cores = 4L
  )


job <- job_table %>%
  filter(job_id == !!job_id)

tax <- job$tax
cg <- job$chem_group
model_type_i <- job$model_type
prior_type_i <- job$prior_type
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
      !is.na(FFG)
    ) %>%
    mutate(
      CAS = factor(CAS),
      genus_phylo = factor(genus),
      chem_group = factor(chem_group),
      FFG = factor(FFG)
    )
}

invertebrate_df <- prepare_model_data(invertebrate_df)

#---- 4. Phylogenetic covariance matrix ----------------------------------------

make_genus_cov <- function(tree) {
  ape::vcv(tree, corr = TRUE)
}

genus_mat_invertebrate <- make_genus_cov(invertebrate_tree)

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

df_subset_raw <- invertebrate_df %>%
  filter(chem_group == cg) %>%
  droplevels()

matched <- match_data_and_genus_mat(
  df = df_subset_raw,
  genus_mat = genus_mat_invertebrate
)

df_subset <- matched$data %>%
  mutate(
    log_rel_sensi_c = log_rel_sensi - mean(log_rel_sensi, na.rm = TRUE),
    log_mass_c = log_mass - mean(log_mass, na.rm = TRUE)
  ) %>%
  droplevels()

genus_mat_subset <- matched$genus_mat

n_obs <- nrow(df_subset)
n_CAS <- n_distinct(df_subset$CAS)
n_genus <- n_distinct(df_subset$genus_phylo)
n_FFG <- n_distinct(df_subset$FFG)

mean_y <- mean(df_subset$log_rel_sensi_c, na.rm = TRUE)
mean_x <- mean(df_subset$log_mass_c, na.rm = TRUE)
sd_y <- sd(df_subset$log_rel_sensi_c, na.rm = TRUE)
sd_x <- sd(df_subset$log_mass_c, na.rm = TRUE)


#---- 6. Define formula --------------------------------------------------------

if (model_type_i == "H0") {
  
  formula_i <- bf(
    log_rel_sensi_c ~ 0 +
      (1 | CAS) +
      (1 | gr(genus_phylo, cov = genus_mat))
  )
  
} else if (model_type_i == "H1") {
  
  formula_i <- bf(
    log_rel_sensi_c ~ 0 + log_mass_c +
      (1 | CAS) +
      (1 | gr(genus_phylo, cov = genus_mat))
  )
  
} else if (model_type_i == "H2") {
  
  formula_i <- bf(
    log_rel_sensi_c ~ 0 + FFG + log_mass_c +
      (1 | CAS) +
      (1 | gr(genus_phylo, cov = genus_mat))
  )
  
} else if (model_type_i == "H3") {
  
  formula_i <- bf(
    log_rel_sensi_c ~ 0 + FFG + FFG:log_mass_c +
      (1 | CAS) +
      (1 | gr(genus_phylo, cov = genus_mat))
  )
  
} else {
  stop("Unknown model_type: ", model_type_i)
}


#---- 7. Define family ---------------------------------------------------------

family_i <- gaussian()

#---- 8. Define priors ---------------------------------------------------------

cas_prior_rate_i <- signif(1 / sd_y, 6)
cas_prior_mean_i <- 1 / cas_prior_rate_i

priors_i <- c(
  prior_string(
    paste0("exponential(", cas_prior_rate_i, ")"),
    class = "sd",
    group = "CAS"
  ),
  prior(normal(0, 2), class = "sd", group = "genus_phylo"),
  prior(student_t(3, 0, 0.3), class = "sigma")
)

if (model_type_i != "H0") {
  priors_i <- c(
    priors_i,
    prior(normal(0, 0.3), class = "b")
  )
}

#---- 9. Fit model -------------------------------------------------------------

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

#---- 10. Save model fit -------------------------------------------------------

base_dir <- here("results/invertebrate_final_H0_H3")
out_dir <- file.path(base_dir, "fits")
job_dir <- file.path(base_dir, "job_info")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(job_dir, recursive = TRUE, showWarnings = FALSE)

safe_cg <- gsub("[^A-Za-z0-9_]+", "_", cg)

fit_name <- paste0(
  model_type_i,
  "_invertebrate_",
  safe_cg
)

fit_path <- file.path(out_dir, paste0(fit_name, ".rds"))

saveRDS(fit_model, fit_path)

message("Saved fit to: ", fit_path)

