#----------------------------------------------------------------------------
# This script is running on HPC with array jobs.
# This script will calculate the phylogenetic signal (Pagel’s λ) in species 
# relative sensitivity estimated within each taxonomic group for each CAS.
#----------------------------------------------------------------------------

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

full_df <- readRDS(here("data/6_combined_meta_data/full_meta_data.rds"))

algae_tree <- readRDS(here("data/6_combined_meta_data/algae_genus_tree.rds"))
invertebrate_tree <- readRDS(here("data/6_combined_meta_data/invertebrate_genus_tree.rds"))
fish_tree <- readRDS(here("data/6_combined_meta_data/fish_genus_tree.rds"))

full_tree <- readRDS(here("data/6_combined_meta_data/full_genus_tree.rds"))

dat_list <- list(
  algae = algae_df,
  invertebrate = invertebrate_df,
  fish = fish_df,
  all_taxa = full_df
)

tree_list <- list(
  algae = algae_tree,
  invertebrate = invertebrate_tree,
  fish = fish_tree,
  all_taxa = full_tree
)

chem_groups <- sort(unique(full_df$chem_group))

tax_groups <- c("algae", "invertebrate", "fish", "all_taxa")

#---- 1. Define SLURM array task ------------------------------------------------
# This script analyses 24 groups:
# 6 chemical-use groups × 3 separate taxonomic groups = 18
# 6 chemical-use groups × all taxa pooled = 6
# Total = 24

task_grid <- tidyr::expand_grid(
  chem_group = chem_groups,
  tax_group = tax_groups
)

task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))

this_chem_group <- task_grid$chem_group[task_id]
this_tax_group <- task_grid$tax_group[task_id]

message("Running SLURM task: ", task_id, " / ", nrow(task_grid))
message("Chemical-use group: ", this_chem_group)
message("Taxonomic group: ", this_tax_group)

dat_this <- dat_list[[this_tax_group]]
tree_this <- tree_list[[this_tax_group]]

#---- 2. Output directory -------------------------------------------------------

fits_dir <- here("results/phylogenetic_signal_by_tax_group")
dir.create(fits_dir, recursive = TRUE, showWarnings = FALSE)

#---- 3. Build genus × CAS table ------------------------------------------------
# For each chemical-use × taxonomic group:
#   1. average log_rel_sensi at genus × CAS level
#   2. retain sample-size information for each CAS
#
# The phylogeny is genus-level, so each genus contributes one mean sensitivity value per CAS.

build_wide_nomask <- function(dat, chemgroup) {
  
  dat_sub <- dat %>%
    filter(chem_group == chemgroup) %>%
  
  genus_cas_summary <- dat_sub %>%
    group_by(genus, CAS) %>%
    summarise(
      val = mean(log_rel_sensi, na.rm = TRUE),
      .groups = "drop"
    )
  
  wide_trait <- genus_cas_summary %>%
    select(genus, CAS, val) %>%
    pivot_wider(
      names_from = CAS,
      values_from = val
    ) %>%
    rename(phylo_name = genus) %>%
    distinct()
  
  rownames(wide_trait) <- wide_trait$phylo_name
  
  sample_size <- dat_sub %>%
    group_by(CAS) %>%
    summarise(
      n_records_raw = n(),
      n_species_raw = n_distinct(latin),
      n_genera_data = n_distinct(genus),
      .groups = "drop"
    )
  
  list(
    wide_trait = wide_trait,
    sample_size = sample_size
  )
}

lambda_one_trait <- function(vec_named, tree, min_n = 4) {
  
  vec_named <- vec_named[!is.na(vec_named)]
  
  keep <- intersect(names(vec_named), tree$tip.label)
  
  tree_k <- drop.tip(tree, setdiff(tree$tip.label, keep))
  x_k <- vec_named[tree_k$tip.label]
  
  fit <- try(
    phylosig(tree_k, x_k, method = "lambda", test = TRUE),
    silent = TRUE
  )
  
  if (inherits(fit, "try-error")) {
    return(tibble(
      lambda = NA_real_,
      logL = NA_real_,
      p = NA_real_,
      n_genus_used = length(keep),
      status = "phylosig_failed"
    ))
  }
  
  tibble(
    lambda = unname(fit$lambda),
    logL = unname(fit$logL),
    p = unname(fit$P),
    n_genus_used = length(keep),
    status = "ok"
  )
}

compute_lambda_one_group <- function(dat, tree, chemgroup, taxgroup, fits_dir) {

  message("Chemical-use group: ", chemgroup)
  message("Taxonomic group: ", taxgroup)
  
  built <- build_wide_nomask(dat = dat, chemgroup = chemgroup)
  
  wide <- built$wide_trait
  sample_size <- built$sample_size
  
  cas_cols <- setdiff(colnames(wide), "phylo_name")
  
  df_out <- purrr::map_dfr(cas_cols, function(cas) {
    
    vec <- setNames(wide[[cas]], wide$phylo_name)
    
    out <- lambda_one_trait(vec_named = vec, tree = tree)
    
    out %>%
      mutate(
        chem_group = chemgroup,
        tax_group = taxgroup,
        CAS = cas,
        trait_col = "log_rel_sensi",
        phylo_col = "genus",
        .before = 1
      )
  }) %>%
    left_join(sample_size, by = "CAS") %>%
    mutate(
      p_adj_BH_group = p.adjust(p, method = "BH")
    ) %>%
    relocate(
      n_records_raw,
      n_species_raw,
      n_genera_data,
      n_genus_used,
      .after = CAS
    ) %>%
    arrange(desc(lambda))
  
  out_file <- file.path(
    fits_dir,
    paste0(
      "phytools_lambda_by_CAS_",
      chemgroup, "_",
      taxgroup,
      ".rds"
    )
  )
  
  saveRDS(df_out, out_file)
  
  message("Saved result to: ", out_file)
  
  df_out
}

lambda_results <- compute_lambda_one_group(
  dat = dat_this,
  tree = tree_this,
  chemgroup = this_chem_group,
  taxgroup = this_tax_group,
  fits_dir = fits_dir
)

message("Finished task.")
