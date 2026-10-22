# Phylogenetic imputation of masked genus-CAS sensitivity values.
# Scenarios: 30% and 50% masking.

rm(list = ls())

library(here)
library(tidyverse)
library(ape)
library(Rphylopars)

set.seed(2026)

invertebrate_df <- readRDS(
  here("data/6_combined_meta_data/invertebrate_meta_data.rds")
)

inv_tree <- readRDS(
  here("data/6_combined_meta_data/invertebrate_genus_tree.rds")
)

chem_groups <- c(
  "Industrial_Chemical", "Pharmaceutical", "P_Arthropod_control",
  "P_Plant_control", "P_Microbial_control", "P_Mixed_target"
)

mask_fracs <- c(0.30, 0.50)

fits_dir <- here("data/9_phylopars_pred")
dir.create(fits_dir, recursive = TRUE, showWarnings = FALSE)


build_and_mask_wide <- function(dat, chemgroup, mask_frac, tree) {
  
  # Build genus-level wide matrix.
  wide <- dat %>%
    filter(chem_group == chemgroup) %>%
    group_by(genus, CAS) %>%
    summarise(
      val = mean(log_sensi, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_wider(names_from = CAS, values_from = val) %>%
    rename(species = genus) %>%
    distinct() %>%
    as.data.frame()
  
  # Keep only genera present in the tree.
  common_genera <- intersect(wide$species, tree$tip.label)
  
  wide <- wide %>%
    filter(species %in% common_genera)
  
  rownames(wide) <- wide$species
  
  trait_mv <- wide
  cas_cols <- setdiff(names(trait_mv), "species")
  
  # Mask a fraction of available genus-CAS values per CAS.
  masked_pairs <- map_dfr(cas_cols, function(cn) {
    
    genus_nonmiss <- trait_mv$species[!is.na(trait_mv[[cn]])]
    
    n_mask <- max(1, floor(mask_frac * length(genus_nonmiss)))
    n_mask <- min(n_mask, length(genus_nonmiss))
    
    tibble(
      species = sample(genus_nonmiss, n_mask),
      CAS = cn
    )
  })
  
  if (nrow(masked_pairs) > 0) {
    for (j in seq_len(nrow(masked_pairs))) {
      trait_mv[
        trait_mv$species == masked_pairs$species[j],
        masked_pairs$CAS[j]
      ] <- NA
    }
  }
  
  tree_pruned <- keep.tip(tree, trait_mv$species)
  
  list(
    wide_full = wide,
    trait_mv = trait_mv,
    masked_pairs = masked_pairs,
    tree_pruned = tree_pruned
  )
}


fit_phylopars_fast <- function(trait_mv, tree) {
  
  phylopars(
    trait_data = trait_mv,
    tree = tree,
    model = "lambda",
    pheno_error = FALSE,
    phylo_correlated = FALSE,
    pheno_correlated = FALSE
  )
}


fit_and_save_all <- function(mask_frac) {
  
  all_results <- list()
  
  for (g in chem_groups) {
    
    message("Fitting: ", g, " | mask ", mask_frac)
    
    wm <- build_and_mask_wide(
      dat = invertebrate_df,
      chemgroup = g,
      mask_frac = mask_frac,
      tree = inv_tree
    )
    
    fit <- fit_phylopars_fast(
      trait_mv = wm$trait_mv,
      tree = wm$tree_pruned
    )
    
    res <- list(
      chem_group = g,
      mask_frac = mask_frac,
      fit = fit,
      wide_full = wm$wide_full,
      trait_mv = wm$trait_mv,
      masked_pairs = wm$masked_pairs,
      tree_pruned = wm$tree_pruned,
      lambda = suppressWarnings(as.numeric(fit$model$lambda))
    )
    
    all_results[[g]] <- res
    
    saveRDS(
      res,
      file = file.path(
        fits_dir,
        paste0(
          "phylopars_fit_", g, "_mask",
          round(mask_frac * 100), ".rds"
        )
      )
    )
  }
  
  saveRDS(
    all_results,
    file = file.path(
      fits_dir,
      paste0(
        "phylopars_fits_all_mask",
        round(mask_frac * 100), ".rds"
      )
    )
  )
  
  invisible(all_results)
}


all_results_list <- lapply(mask_fracs, fit_and_save_all)

names(all_results_list) <- paste0(
  "mask",
  round(mask_fracs * 100)
)

saveRDS(
  all_results_list,
  file = file.path(fits_dir, "phylopars_fits_all_masks.rds")
)

