# Combine predicted sensitivity, chemical-use groups, 
# species taxonomy, phylogeny, body mass and functional groups 
# into a single dataset for downstream analysis.

rm(list = ls())

library(here)
library(tidyverse)

source(here("code/functions.r"))
options(dplyr.summarise.inform = FALSE)


out_dir <- here("data/6_combined_meta_data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ----- 1. Prepare chemical sensitivity data -----------------------------------
# 1.1 Read predicted sensitivity
predicted_sensitivity <- read_csv(
  here("data/2_predicted_sensitivity/predictions.csv"),
  show_col_types = FALSE) %>%
  mutate(latin = str_squish(str_to_lower(species)),
         CAS = as.character(CAS))

prediction_summary <- predicted_sensitivity %>%
  summarise(n_rows = n(),
            n_species = n_distinct(latin),
            n_CAS = n_distinct(CAS),
            n_duration = n_distinct(duration))

print(prediction_summary)


# 1.2 Add chemical-use groups
chem_groups_source <- read_csv(
  here("data/6_combined_meta_data/chem_groups.csv"),
  show_col_types = FALSE
) %>%
  mutate(CAS = str_squish(as.character(CAS)),
         broad_group = str_squish(use_group),
         p_sub_group = str_squish(Pesticide_sub_group)) %>%
  select(CAS, broad_group, p_sub_group) %>%
  filter(!is.na(broad_group), broad_group != "Metal")

chem_groups_source %>%
  count(broad_group, sort = TRUE) %>%
  print()


# 1.3 Map PPDB pesticide types to final chemical groups
chem_groups_source <- chem_groups_source %>%
  mutate(
    p_sub_group = str_replace_all(p_sub_group, "Insectiicde", "Insecticide"),
    p_target_group = map_chr(p_sub_group, map_subgroup_to_target,
                             tag_to_group = tag_to_group)
  )

keep_pesticide <- chem_groups_source %>%
  filter(broad_group == "Pesticide") %>%
  filter(p_target_group %in% c("arthropod-control", "microbial-control",
                               "plant-control", "mixed-target")) %>%
  pull(CAS)

chem_groups_source <- chem_groups_source %>%
  filter(broad_group != "Pesticide" | CAS %in% keep_pesticide) %>%
  mutate(
    final_group = case_when(
      broad_group == "Industrial Chemical" ~ "Industrial_Chemical",
      broad_group == "Pharmaceutical" ~ "Pharmaceutical",
      p_target_group == "arthropod-control" ~ "P_Arthropod_control",
      p_target_group == "microbial-control" ~ "P_Microbial_control",
      p_target_group == "plant-control" ~ "P_Plant_control",
      p_target_group == "mixed-target" ~ "P_Mixed_target",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(final_group))

chemical_group_summary <- chem_groups_source %>%
  summarise(n_CAS = n_distinct(CAS), .by = final_group) %>%
  arrange(final_group)

print(chemical_group_summary)


# 1.4 Keep selected chemicals and add chemical-use groups
sensitivity_subCAS <- predicted_sensitivity %>%
  inner_join(
    chem_groups_source %>% select(CAS, chem_group = final_group),
    by = "CAS"
  )

selected_prediction_summary <- sensitivity_subCAS %>%
  summarise(n_rows = n(),
            n_species = n_distinct(latin),
            n_CAS = n_distinct(CAS),
            n_chem_groups = n_distinct(chem_group))

print(selected_prediction_summary)


# ----- 2. Read species taxonomy, phylogeny, body mass and final FFG -----------
spe_tax_clean_tree_mass_FFG <- readRDS(
  file.path(out_dir, "spe_tax_clean_tree_mass_FFG.rds")
) %>%
  mutate(latin = str_squish(str_to_lower(latin)),
         genus = str_squish(str_to_lower(genus)),
         group = str_squish(str_to_lower(group)))

trait_summary <- spe_tax_clean_tree_mass_FFG %>%
  group_by(group) %>%
  summarise(n_species = n_distinct(latin),
            n_genus = n_distinct(genus),
            n_FFG = n_distinct(FFG),
            .groups = "drop")

print(trait_summary)
unique(spe_tax_clean_tree_mass_FFG$FFG)


# -----3. Combine predicted sensitivity with species traits --------------------
spe_tax_clean_tree_mass_FFG_sens <- sensitivity_subCAS %>%
  semi_join(spe_tax_clean_tree_mass_FFG %>% distinct(latin), by = "latin") %>%
  left_join(spe_tax_clean_tree_mass_FFG, by = "latin") %>%
  filter(
    duration == case_when(
      group == "algae" ~ 72,
      group %in% c("invertebrate", "fish") ~ 96,
      TRUE ~ NA_real_
    )
  ) %>%
  mutate(log_sensi = predicted,
         sensi = 10^log_sensi) %>%
  select(-predicted)

combined_summary <- spe_tax_clean_tree_mass_FFG_sens %>%
  summarise(n_rows = n(),
            n_species = n_distinct(latin),
            n_genus = n_distinct(genus),
            n_FFG = n_distinct(FFG),
            n_CAS = n_distinct(CAS),
            n_chem_groups = n_distinct(chem_group))

print(combined_summary)


# ----- 4. Calculate relative sensitivity --------------------------------------
ref_sensitivity <- spe_tax_clean_tree_mass_FFG_sens %>%
  filter(latin == "daphnia magna") %>%
  distinct(CAS, sensi) %>%
  group_by(CAS) %>%
  summarise(ref_sensi = mean(sensi, na.rm = TRUE), .groups = "drop")

spe_tax_clean_tree_mass_FFG_sens <- spe_tax_clean_tree_mass_FFG_sens %>%
  left_join(ref_sensitivity, by = "CAS") %>%
  mutate(rel_sensi = sensi / ref_sensi,
         log_rel_sensi = log10(rel_sensi))

spe_tax_clean_tree_mass_FFG_sens <- spe_tax_clean_tree_mass_FFG_sens %>%
  filter(!is.na(log_mass), !is.na(log_rel_sensi), is.finite(log_rel_sensi))

relative_sensitivity_summary <- spe_tax_clean_tree_mass_FFG_sens %>%
  summarise(n_rows = n(),
            n_species = n_distinct(latin),
            n_CAS = n_distinct(CAS))

print(relative_sensitivity_summary)


# ----- 5. Prepare full meta data for downstream analysis ----------------------
full_meta_data <- spe_tax_clean_tree_mass_FFG_sens %>%
  select(chem_group, CAS, latin, genus, family, order, class, phylum, kingdom,
         group, log_mass, BMsource, BM_match_level, FFG, FFG_source,
         ref_sensi, log_sensi, rel_sensi, log_rel_sensi, duration) %>%
  rename(tax_group = group) %>%
  distinct()


# ----- 6. Check alignment between phylogeny tree and meta data ----------------
full_tree <- readRDS(file.path(out_dir, "full_genus_tree.rds"))

full_tree_check <- tibble(
  item = c(
    "n_genus_in_dataset",
    "n_tips_in_full_tree",
    "genus_in_tree_not_dataset",
    "genus_in_dataset_not_tree"
  ),
  value = c(
    n_distinct(full_meta_data$genus),
    length(full_tree$tip.label),
    length(setdiff(full_tree$tip.label, unique(full_meta_data$genus))),
    length(setdiff(unique(full_meta_data$genus), full_tree$tip.label))
  )
)

print(full_tree_check)


# ----- 7. Save full metadata and taxon-specific metadata for downstream analysis ----
saveRDS(full_meta_data, file.path(out_dir, "full_meta_data.rds"))

save_taxon_dataset <- function(dat, taxon, tree_file, output_file) {
  
  dat_taxon <- dat %>%
    filter(tax_group == taxon)
  
  genus_tree <- readRDS(file.path(out_dir, tree_file))
  
  check_tbl <- tibble(
    group = taxon,
    n_rows = nrow(dat_taxon),
    n_CAS = n_distinct(dat_taxon$CAS),
    n_species = n_distinct(dat_taxon$latin),
    n_genus = n_distinct(dat_taxon$genus),
    n_FFG = n_distinct(dat_taxon$FFG),
    n_tree_tips = length(genus_tree$tip.label),
    genus_in_tree_not_dataset = length(setdiff(genus_tree$tip.label,
                                               unique(dat_taxon$genus))),
    genus_in_dataset_not_tree = length(setdiff(unique(dat_taxon$genus),
                                               genus_tree$tip.label))
  )
  
  saveRDS(dat_taxon, file.path(out_dir, output_file))
  
  check_tbl
}

algae_check <- save_taxon_dataset(
  full_meta_data, "algae", "algae_genus_tree.rds", "algae_meta_data.rds"
)

invertebrate_check <- save_taxon_dataset(
  full_meta_data, "invertebrate",
  "invertebrate_genus_tree.rds", "invertebrate_meta_data.rds"
)

fish_check <- save_taxon_dataset(
  full_meta_data, "fish", "fish_genus_tree.rds", "fish_meta_data.rds"
)

taxon_output_summary <- bind_rows(algae_check, invertebrate_check, fish_check)
print(taxon_output_summary)


# 8. Dataset summary
# Number of species, genera, functional groups, chemicals and species-chemical pairs
full_meta_data_summary <- full_meta_data %>%
  summarise(
    n_species = n_distinct(latin),
    n_genus = n_distinct(genus),
    n_FFG = n_distinct(FFG),
    n_chemicals = n_distinct(CAS),
    n_species_chemical_pairs = n_distinct(latin, CAS)
  )

print(full_meta_data_summary)


# Summary by broad taxonomic group
taxon_meta_data_summary <- full_meta_data %>%
  group_by(tax_group) %>%
  summarise(
    n_rows = n(),
    n_species = n_distinct(latin),
    n_genus = n_distinct(genus),
    n_FFG = n_distinct(FFG),
    n_chemicals = n_distinct(CAS),
    n_species_chemical_pairs = n_distinct(latin, CAS),
    .groups = "drop"
  )

print(taxon_meta_data_summary)


# Summary by taxonomic group and chemical-use group
taxon_chemical_group_summary <- full_meta_data %>%
  group_by(tax_group, chem_group) %>%
  summarise(
    n_species = n_distinct(latin),
    n_genus = n_distinct(genus),
    n_FFG = n_distinct(FFG),
    n_chemicals = n_distinct(CAS),
    n_species_chemical_pairs = n_distinct(latin, CAS),
    .groups = "drop"
  )

print(taxon_chemical_group_summary)

