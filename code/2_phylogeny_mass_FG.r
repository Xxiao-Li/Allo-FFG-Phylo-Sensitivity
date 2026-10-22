# Build phylogeny tree, collect species body mass, and assign species 
# functional groups. 
# The final output includes a cleaned dataset with taxonomy, body mass, and FFG, 
# as well as pruned phylogeny trees for each group and a combined tree.


rm(list = ls())

library(here)
library(tidyverse)
library(rotl)
library(ape)

options(dplyr.summarise.inform = FALSE)

out_dir <- here("data/6_combined_meta_data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ----- 1. Taxonomy and name corrections ---------------------------------------
spe_tax_clean <- readRDS(
  here("data/1_observed_sensitivity/combined_observed_sensitivity/combined_observed_toxicity.rds")
) %>%
  select(latin, kingdom, phylum, `class`, order, family, genus, group) %>%
  mutate(
    latin = str_squish(str_to_lower(latin)),
    across(c(kingdom, phylum, `class`, order, family, genus), ~ str_squish(str_to_lower(.x))),
    group = str_squish(str_to_lower(group))
  ) %>%
  filter(!is.na(latin), latin != "", !is.na(genus), genus != "", !is.na(group)) %>%
  distinct()

print(spe_tax_clean %>%
        group_by(group) %>%
        summarise(n_species = n_distinct(latin), n_genus = n_distinct(genus)))

cat("Unique species across groups:", n_distinct(spe_tax_clean$latin), "\n",
    "Unique genera across groups:", n_distinct(spe_tax_clean$genus), "\n")


# ----- 2. Match genera to Open Tree and build trees ---------------------------
groups <- c("algae", "invertebrate", "fish")

matches <- setNames(lapply(groups, function(g) {
  genera <- sort(unique(spe_tax_clean$genus[spe_tax_clean$group == g]))
  tnrs_match_names(
    names = genera,
    context_name = if (g == "algae") "All life" else "Animals",
    do_approximate_matching = TRUE
  )
}), groups)

tx_full <- bind_rows(matches) %>%
  filter(!is.na(ott_id)) %>%
  distinct(ott_id, .keep_all = TRUE)

make_tree <- function(tx) {
  ids <- unique(as.integer(tx$ott_id[!is.na(tx$ott_id)]))
  ids <- unname(ids[is_in_tree(ids)])
  tree <- tol_induced_subtree(ott_ids = ids)
  tip_id <- as.numeric(sub(".*_ott", "", tree$tip.label))
  tree$tip.label <- str_to_lower(tx$search_string[match(tip_id, tx$ott_id)])
  compute.brlen(tree, method = "Grafen")
}

# raw_trees <- lapply(matches, make_tree)
# full_tree_raw <- make_tree(tx_full)
# # save phylogeny trees before body mass and FFG
# saveRDS(raw_trees$algae, "data/3_phylogeny/algae_genus_tree_raw.rds")
# saveRDS(raw_trees$invertebrate, "data/3_phylogeny/invertebrate_genus_tree_raw.rds")
# saveRDS(raw_trees$fish, "data/3_phylogeny/fish_genus_tree_raw.rds")
# saveRDS(full_tree_raw, "data/3_phylogeny/full_genus_tree_raw.rds")
# directly read the saved tress to save time
raw_trees <- list(
  algae = readRDS("data/3_phylogeny/algae_genus_tree_raw.rds"),
  invertebrate = readRDS("data/3_phylogeny/invertebrate_genus_tree_raw.rds"),
  fish = readRDS("data/3_phylogeny/fish_genus_tree_raw.rds")
)
full_tree_raw <- readRDS("data/3_phylogeny/full_genus_tree_raw.rds")

print(tibble(
  group = names(raw_trees),
  tips = map_int(raw_trees, ~ length(.x$tip.label))
))

spe_tax_clean_tree <- spe_tax_clean %>%
  filter(genus %in% full_tree_raw$tip.label)

print(spe_tax_clean_tree %>%
        summarise(n_species = n_distinct(latin), n_genus = n_distinct(genus)))

print(spe_tax_clean_tree %>%
        group_by(group) %>%
        summarise(n_species = n_distinct(latin), n_genus = n_distinct(genus)))


# ----- 3. Body mass: species match, then genus mean ---------------------------
sp_BM <- readRDS(here("data/4_body_mass/BM_species_summary.rds")) %>%
  mutate(latin = str_squish(str_to_lower(latin)),
         genus = str_squish(str_to_lower(genus))) %>%
  rename(BMsource = sources)

genus_avg <- sp_BM %>%
  group_by(genus) %>%
  summarise(
    genus_mass = mean(mean_mass_g, na.rm = TRUE),
    genus_source = paste(unique(na.omit(BMsource)), collapse = "; "),
    .groups = "drop"
  ) %>%
  mutate(genus_mass = if_else(is.nan(genus_mass), NA_real_, genus_mass))

spe_mass <- spe_tax_clean_tree %>%
  distinct(latin, genus) %>%
  left_join(sp_BM %>% select(latin, mean_mass_g, BMsource), by = "latin") %>%
  left_join(genus_avg, by = "genus") %>%
  mutate(
    BM_match_level = case_when(
      !is.na(mean_mass_g) ~ "species",
      !is.na(genus_mass) ~ "genus",
      TRUE ~ "unmatched"
    ),
    BMsource = if_else(BM_match_level == "genus", genus_source, BMsource),
    mean_mass_g = coalesce(mean_mass_g, genus_mass)
  ) %>%
  select(-genus_mass, -genus_source) %>%
  distinct(latin, .keep_all = TRUE)

print(spe_mass %>% count(BM_match_level))

mass_fixes <- c(
  "nitzschia closterium" = 1.58e-09,
  "neanthes arenaceodentata" = 0.0147,
  "lysmata seticaudata" = 0.77,
  "odontesthes argentinensis" = 92.5,
  "oryzias melastigma" = 1.73,
  "oryzias dancena" = 1.73,
  "oryzias sinensis" = 1.73,
  "pagrus major" = 104.9
)

spe_mass <- spe_mass %>%
  mutate(
    mean_mass_g = if_else(latin %in% names(mass_fixes), unname(mass_fixes[latin]), mean_mass_g),
    log_mass = log10(mean_mass_g)
  ) %>%
  filter(!str_detect(latin, "ulva"), is.finite(log_mass))

spe_tax_clean_tree_mass <- spe_tax_clean_tree %>%
  left_join(
    spe_mass %>% select(latin, mean_mass_g, BMsource, BM_match_level, log_mass),
    by = "latin"
  ) %>%
  filter(!is.na(log_mass))

print(spe_tax_clean_tree_mass %>%
        summarise(
          n_species = n_distinct(latin),
          n_genus = n_distinct(genus)
        ))

print(spe_tax_clean_tree_mass %>%
        group_by(group) %>%
        summarise(
          n_species = n_distinct(latin),
          n_genus = n_distinct(genus),
          n_BM_species = sum(BM_match_level == "species"),
          n_BM_genus = sum(BM_match_level == "genus")
        ))


# ----- 4. Final functional feeding groups -------------------------------------
# 4.1 Read FFG files
algae_FG <- read_csv(here("data/5_FFG/algae_FG.csv"), show_col_types = FALSE) %>%
  mutate(
    latin = str_squish(str_to_lower(latin)),
    FG = str_squish(FG),
    algae.tax.group = str_squish(algae.tax.group),
    FG_source = "Alina"
  )

invert_fish_FG <- read_csv(here("data/5_FFG/invert_vert_FG.csv"), show_col_types = FALSE) %>%
  mutate(
    latin = str_squish(str_to_lower(latin)),
    FG = str_squish(FG),
    FG_source = str_squish(FG_source),
    algae.tax.group = NA_character_)

FFG_source <- bind_rows(algae_FG, invert_fish_FG)

# Check species without any FFG record
missing_FG_species <- spe_tax_clean_tree_mass %>%
  anti_join(FFG_source %>% distinct(latin), by = "latin") %>%
  select(group, latin, genus, family, `class`, order, phylum) %>%
  distinct() %>%
  arrange(group, latin)

print(missing_FG_species %>% count(group))

# 4.2 Join and convert directly to final FFG
algae_groups <- c(
  "Blue-green" = "Cyanobacteria",
  "Green" = "Green algae",
  "Charophytes" = "Green algae",
  "Euglenoids" = "Green algae",
  "Diatoms" = "Diatom",
  "Chrysophytes" = "Other algae",
  "Cryptomonads" = "Other algae",
  "Dinoflagellates" = "Other algae",
  "Haptophytes" = "Other algae",
  "Stramenopiles" = "Other algae"
)

spe_tax_clean_tree_mass_FFG <- spe_tax_clean_tree_mass %>%
  left_join(FFG_source, by = "latin") %>%
  filter(!is.na(FG), !FG %in% c("Unknown", "Parasitic")) %>%
  mutate(
    FFG = if_else(group == "algae", dplyr::recode(algae.tax.group, !!!algae_groups), FG),
    FFG_source = FG_source) %>%
  select(-FG, -FG_source, -algae.tax.group) %>%
  separate_rows(FFG, sep = "\\s*,\\s*") %>%
  mutate(
    FFG = str_trim(FFG),
    FFG = str_replace_all(FFG, "\\s*/\\s*", "_"),
    FFG = str_squish(FFG),
    FFG = if_else(
      group == "fish" & FFG %in% c("Detritivore_Herbivore", "Herbivore", "Detritivore"),
      "Detritivore_Herbivore", FFG
    ),
    FFG = case_when(
      group == "invertebrate" ~ paste0("I_", FFG),
      group == "fish" ~ paste0("F_", FFG),
      TRUE ~ FFG
    ),
    FFG = if_else(phylum == "arthropoda", paste0(FFG, "_A"), FFG)
  ) %>%
  filter(!is.na(FFG), FFG != "") %>%
  distinct()

print(spe_tax_clean_tree_mass_FFG %>%
        distinct(group, latin, genus, FFG) %>%
        group_by(group, FFG) %>%
        summarise(n_species = n_distinct(latin), n_genus = n_distinct(genus)) %>%
        arrange(group, FFG))

print(spe_tax_clean_tree_mass_FFG %>%
        group_by(group) %>%
        summarise(
          n_species = n_distinct(latin),
          n_genus = n_distinct(genus),
          n_FFG = n_distinct(FFG)
        ))

print(spe_tax_clean_tree_mass_FFG %>%
        summarise(
          n_species = n_distinct(latin),
          n_genus = n_distinct(genus),
          n_FFG = n_distinct(FFG)
        ))

unique(spe_tax_clean_tree_mass_FFG$FFG)
spe_tax_clean_tree_mass_FFG %>%
  distinct(FFG, latin) %>%
  count(FFG, name = "n_species", sort = TRUE) %>%
  print(n = Inf)

# ----- 5. Prune trees to final dataset ----------------------------------------
final_genus <- spe_tax_clean_tree_mass_FFG %>%
  distinct(group, genus)

prune <- function(tree, tips) {
  keep.tip(tree, intersect(tree$tip.label, tips))
}

final_trees <- setNames(lapply(groups, function(g) {
  prune(raw_trees[[g]], final_genus$genus[final_genus$group == g])
}), groups)

full_tree <- prune(full_tree_raw, final_genus$genus)


# 6. Save all----
saveRDS(
  spe_tax_clean_tree_mass_FFG,
  file.path(out_dir, "spe_tax_clean_tree_mass_FFG.rds")
)

saveRDS(final_trees$algae, file.path(out_dir, "algae_genus_tree.rds"))
saveRDS(final_trees$invertebrate, file.path(out_dir, "invertebrate_genus_tree.rds"))
saveRDS(final_trees$fish, file.path(out_dir, "fish_genus_tree.rds"))
saveRDS(full_tree, file.path(out_dir, "full_genus_tree.rds"))
