# This script is running on HPC with array jobs.
# This script calculates Pagel's lambda after randomly removing 0%, 5%, 10%,
# or 15% of predicted sensitivity values.

rm(list = ls())

library(here)
library(tidyverse)
library(ape)
library(phytools)

set.seed(2026)

#---- 0. Load data and basic settings ------------------------------------------
algae_df <- readRDS(here("data/6_combined_meta_data/algae_meta_data.rds"))
invertebrate_df <- readRDS(here("data/6_combined_meta_data/invertebrate_meta_data.rds"))
fish_df <- readRDS(here("data/6_combined_meta_data/fish_meta_data.rds"))

algae_tree <- readRDS(here("data/6_combined_meta_data/algae_genus_tree.rds"))
invertebrate_tree <- readRDS(here("data/6_combined_meta_data/invertebrate_genus_tree.rds"))
fish_tree <- readRDS(here("data/6_combined_meta_data/fish_genus_tree.rds"))

dat_list <- list(
  algae = algae_df,
  invertebrate = invertebrate_df,
  fish = fish_df
)

tree_list <- list(
  algae = algae_tree,
  invertebrate = invertebrate_tree,
  fish = fish_tree
)

chem_groups <- sort(unique(c(algae_df$chem_group,
                             invertebrate_df$chem_group,
                             fish_df$chem_group)))
tax_groups <- c("algae", "invertebrate", "fish")

drop_rates <- c(0, 0.05, 0.10, 0.15)
n_repeats <- 20

#---- 1. Define SLURM array task ------------------------------------------------
task_grid <- tidyr::expand_grid(
  chem_group = chem_groups,
  tax_group = tax_groups
)

task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))

if (is.na(task_id) || !task_id %in% seq_len(nrow(task_grid))) {
  stop("Invalid SLURM_ARRAY_TASK_ID")
}

this_chem_group <- task_grid$chem_group[task_id]
this_tax_group <- task_grid$tax_group[task_id]

message("Running SLURM task: ", task_id, " / ", nrow(task_grid))
message("Chemical-use group: ", this_chem_group)
message("Taxonomic group: ", this_tax_group)

dat_this <- dat_list[[this_tax_group]]
tree_this <- tree_list[[this_tax_group]]

#---- 2. Output directory -------------------------------------------------------
fits_dir <- here("data/7_phylogenetic_signal_by_tax_group", "lambda_random_exclusion")
dir.create(fits_dir, recursive = TRUE, showWarnings = FALSE)

#---- 3. Functions --------------------------------------------------------------

lambda_one_trait <- function(vec_named, tree, min_n = 4) {
  
  vec_named <- vec_named[!is.na(vec_named)]
  keep <- intersect(names(vec_named), tree$tip.label)
  
  if (length(keep) < min_n) {
    return(tibble(lambda = NA_real_, logL = NA_real_, p = NA_real_,
                  n_genus_used = length(keep), status = "too_few_genera"))
  }
  
  tree_k <- keep.tip(tree, keep)
  x_k <- vec_named[tree_k$tip.label]
  
  fit <- try(
    phylosig(tree_k, x_k, method = "lambda", test = TRUE),
    silent = TRUE
  )
  
  if (inherits(fit, "try-error")) {
    return(tibble(lambda = NA_real_, logL = NA_real_, p = NA_real_,
                  n_genus_used = length(keep), status = "phylosig_failed"))
  }
  
  tibble(
    lambda = unname(fit$lambda),
    logL = unname(fit$logL),
    p = unname(fit$P),
    n_genus_used = length(keep),
    status = "ok"
  )
}

compute_lambda_one_run <- function(dat, tree, chemgroup, taxgroup, drop_rate, repeat_id) {
  
  set.seed(2026 + repeat_id)
  
  dat_sub <- dat %>%
    filter(chem_group == chemgroup) %>%
    group_by(CAS) %>%
    mutate(drop_order = sample.int(n())) %>%
    filter(drop_order > floor(n() * drop_rate)) %>%
    ungroup() %>%
    select(-drop_order)
  
  genus_cas_summary <- dat_sub %>%
    group_by(genus, CAS) %>%
    summarise(val = mean(log_rel_sensi, na.rm = TRUE), .groups = "drop")
  
  sample_size <- dat_sub %>%
    group_by(CAS) %>%
    summarise(n_records = n(), n_species = n_distinct(latin),
              n_genera = n_distinct(genus), .groups = "drop")
  
  wide_trait <- genus_cas_summary %>%
    pivot_wider(names_from = CAS, values_from = val) %>%
    rename(phylo_name = genus)
  
  cas_cols <- setdiff(colnames(wide_trait), "phylo_name")
  
  map_dfr(cas_cols, function(cas) {
    vec <- setNames(wide_trait[[cas]], wide_trait$phylo_name)
    lambda_one_trait(vec, tree) %>% mutate(CAS = cas, .before = 1)
  }) %>%
    left_join(sample_size, by = "CAS") %>%
    mutate(drop_rate = drop_rate, repeat_id = repeat_id,
           chem_group = chemgroup, tax_group = taxgroup, .before = 1)
}

#---- 4. Run deletion sensitivity analysis -------------------------------------

run_grid <- bind_rows(
  tibble(drop_rate = 0, repeat_id = 0),
  expand_grid(drop_rate = drop_rates[drop_rates > 0], repeat_id = seq_len(n_repeats))
)

lambda_results <- pmap_dfr(run_grid, function(drop_rate, repeat_id) {
  message("Drop rate: ", drop_rate, "; repeat: ", repeat_id)
  compute_lambda_one_run(dat_this, tree_this, this_chem_group,
                         this_tax_group, drop_rate, repeat_id)
}) %>%
  group_by(drop_rate, repeat_id) %>%
  mutate(p_adj_BH = p.adjust(p, method = "BH")) %>%
  ungroup()

lambda_summary <- lambda_results %>%
  group_by(drop_rate, repeat_id) %>%
  summarise(median_lambda = median(lambda, na.rm = TRUE),
            proportion_significant = mean(p_adj_BH < 0.05, na.rm = TRUE),
            .groups = "drop")

#---- 5. Save results -----------------------------------------------------------
out_file <- file.path(fits_dir, paste0("lumbda_random_exclusion_test_", this_chem_group,
                                       "_", this_tax_group, ".rds"))

saveRDS(list(results = lambda_results, summary = lambda_summary), out_file)


