#----------------------------------------------------------------------------
# Calculate Pagel's lambda before and after removing lineage means,
# and within major invertebrate lineages.
#----------------------------------------------------------------------------

rm(list = ls())

library(here)
library(tidyverse)
library(ape)
library(phytools)

#---- 0. Load data and settings -----------------------------------------------
inv_df <- readRDS(here("data/6_combined_meta_data/invertebrate_meta_data.rds"))
inv_tree <- readRDS(here("data/6_combined_meta_data/invertebrate_genus_tree.rds"))

chem_groups <- c(
  "Industrial_Chemical",
  "P_Arthropod_control",
  "P_Microbial_control",
  "P_Mixed_target",
  "P_Plant_control",
  "Pharmaceutical"
)

task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
cg <- chem_groups[task_id]
message("Chemical group: ", cg)

#---- 1. Assign lineages ------------------------------------------------------
crustacean_classes <- c(
  "branchiopoda", "cephalocarida", "copepoda", "hexanauplia", "ichthyostraca",
  "malacostraca", "maxillopoda", "ostracoda", "remipedia", "thecostraca"
)

inv_df <- inv_df %>%
  filter(genus %in% inv_tree$tip.label) %>%
  mutate(
    tax_class = str_to_lower(`class`),
    tax_phylum = str_to_lower(phylum),
    major_lineage = case_when(
      tax_class == "insecta" ~ "Insects",
      tax_class %in% crustacean_classes ~ "Crustaceans",
      tax_phylum == "arthropoda" ~ str_to_title(tax_class),
      TRUE ~ str_to_title(tax_phylum)
    )
  )

lineage_counts <- inv_df %>%
  distinct(major_lineage, genus) %>%
  count(major_lineage, name = "n_genera")

testable_lineages <- lineage_counts %>%
  filter(n_genera >= 10) %>%
  pull(major_lineage)

inv_raw <- inv_df %>%
  filter(chem_group == cg) %>%
  group_by(major_lineage, genus, CAS) %>%
  summarise(val = mean(log_rel_sensi), .groups = "drop")

inv_centred <- inv_raw %>%
  group_by(major_lineage, CAS) %>%
  mutate(val = val - mean(val)) %>%
  ungroup()

#---- 2. Calculate lambda -----------------------------------------------------
calc_lambda <- function(dat) {
  split(dat, as.character(dat$CAS), drop = TRUE) %>%
    imap_dfr(function(x, cas) {
      if (nrow(x) < 4) return(NULL)
      
      tree <- keep.tip(inv_tree, x$genus)
      trait <- setNames(x$val, x$genus)
      
      fit <- try(
        phylosig(tree, trait[tree$tip.label], method = "lambda", test = TRUE),
        silent = TRUE
      )
      if (inherits(fit, "try-error")) return(NULL)
      
      tibble(CAS = cas, lambda = unname(fit$lambda),
             p = unname(fit$P), n_genera = nrow(x))
    })
}

results <- bind_rows(
  calc_lambda(inv_raw) %>%
    mutate(analysis = "raw", major_lineage = "All"),
  calc_lambda(inv_centred) %>%
    mutate(analysis = "centred", major_lineage = "All"),
  map_dfr(testable_lineages, function(lineage) {
    calc_lambda(filter(inv_raw, major_lineage == lineage)) %>%
      mutate(analysis = "within_lineage", major_lineage = lineage)
  })
) %>%
  group_by(analysis, major_lineage) %>%
  mutate(p_adj_BH = p.adjust(p, method = "BH")) %>%
  ungroup() %>%
  mutate(chem_group = cg, .before = 1)

#---- 3. Save results ---------------------------------------------------------
out_dir <- here("data/7_phylogenetic_signal_by_tax_group/lambda_invertebrate_lineages")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

saveRDS(results, file.path(
  out_dir, paste0("lambda_invertebrate_lineages_", cg, ".rds")
))

