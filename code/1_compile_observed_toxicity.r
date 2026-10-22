# Compile and compare observed sensitivity datasets

rm(list = ls())

library(here)
library(tidyverse)
library(readxl)

source(here("code/functions.r"))

data_dir <- here("data/1_observed_sensitivity")

# ----- 1. ADORE ---------------------------------------------------------------
adore_dir <- file.path(data_dir, "ADORE")

adore <- read_csv(file.path(adore_dir, "ecotox_mortality_processed.csv"),
                  show_col_types = FALSE) %>%
  transmute(
    latin = clean_latin(tax_gs), CAS = as.character(test_cas),
    duration = result_obs_duration_mean,
    conc_log10mgL = result_conc1_mean_log,
    tax_group, tax_class, tax_order, tax_family, tax_genus
  )

taxonomy_adore <- read_csv(file.path(adore_dir, "taxonomy_adore.csv"),
                           show_col_types = FALSE) %>%
  select(-any_of("species"))

adore <- adore %>%
  left_join(taxonomy_adore, by = c("latin" = "original_name")) %>%
  mutate(
    phylum = if_else(`class` %in% c("Malacostraca", "Maxillopoda"),
                     "Arthropoda", phylum),
    CAS_norm = normalize_cas(CAS)
  ) %>%
  standardise_names() %>%
  filter(!is.na(CAS_norm), CAS_norm != "",
         !is.na(latin), latin != "",
         is.finite(conc_log10mgL), !is.na(duration))

adore_species <- adore %>%
  distinct(latin, genus, family, `class`, order, phylum, kingdom, tax_group)

cat("\nADORE species composition:\n")
print(adore_species %>% count(tax_group, sort = TRUE))
cat("\nADORE summary:\n")
print(dataset_summary(adore))
cat("\nADORE species per chemical:\n")
print(species_per_chemical(adore))


# ----- 2. Liang and Sinclair --------------------------------------------------
inv_dir <- file.path(data_dir, "Liang_Sinclair")
inv_files <- list.files(inv_dir, pattern = "\\.xlsx?$", full.names = TRUE)

inv_tox_merge <- map_dfr(inv_files, ~ read_excel(.x, sheet = "merge")) %>%
  transmute(
    latin = clean_latin(Species.Scientific.Name),
    CAS = as.character(CAS.Number),
    conc_log10mgL = log10(Conc.1.Mean..Standardized.),
    duration = 24 * Observed.Duration.Mean..Days.
  )

invert_extra_taxonomy <- read_csv(
  file.path(inv_dir, "invert_extra_taxonomy.csv"),
  show_col_types = FALSE
) %>% select(-any_of("species"))

inv_tox_merge <- inv_tox_merge %>%
  left_join(invert_extra_taxonomy, by = "latin") %>%
  mutate(tax_group = "invertebrate", CAS_norm = normalize_cas(CAS)) %>%
  standardise_names() %>%
  filter(!is.na(CAS_norm), CAS_norm != "",
         !is.na(latin), latin != "",
         is.finite(conc_log10mgL), !is.na(duration)) %>%
  distinct()

inv_species <- inv_tox_merge %>%
  distinct(latin, genus, family, `class`, order, phylum, kingdom)

cat("\nLiang–Sinclair species composition:\n")
print(inv_species %>% count(phylum, sort = TRUE))
cat("\nLiang–Sinclair summary:\n")
print(dataset_summary(inv_tox_merge))
cat("\nLiang–Sinclair species per chemical:\n")
print(species_per_chemical(inv_tox_merge))

# ----- 3. Cross-check ------------------------------------------
datasets <- list(ADORE = adore, Liang_Sinclair = inv_tox_merge)

coverage_summary <- imap_dfr(
  datasets, ~ dataset_summary(.x) %>% mutate(dataset = .y, .before = 1)
)
cat("\nCoverage by dataset:\n")
print(coverage_summary)

# Returns one row per chemical, species, or chemical–species pair,
# with membership in each of the datasets.
crosscheck <- function(columns) {
  imap_dfr(datasets, function(x, dataset_name) {
    x %>%
      filter(if_all(all_of(columns), ~ !is.na(.x) & .x != "")) %>%
      distinct(across(all_of(columns))) %>%
      mutate(dataset = dataset_name)
  }) %>%
    mutate(present = TRUE) %>%
    pivot_wider(names_from = dataset, values_from = present,
                values_fill = FALSE) %>%
    mutate(
      status = case_when(
        ADORE & Liang_Sinclair ~ "ADORE + Liang–Sinclair only",
        ADORE ~ "ADORE only",
        Liang_Sinclair ~ "Liang–Sinclair only"
      )
    )
}

chemical_crosscheck <- crosscheck("CAS_norm")
species_crosscheck <- crosscheck("latin")
pair_crosscheck <- crosscheck(c("CAS_norm", "latin"))

chemical_summary <- chemical_crosscheck %>%
  count(status, name = "n_chemicals")
species_summary <- species_crosscheck %>%
  count(status, name = "n_species")
pair_summary <- pair_crosscheck %>%
  count(status, name = "n_pairs")

cat("\nChemical overlap:\n")
print(chemical_summary)
cat("\nSpecies overlap:\n")
print(species_summary)
cat("\nChemical–species pair overlap:\n")
print(pair_summary)


# 3.1 Species overlap by taxonomic group
species_taxonomy <- bind_rows(
  species_groups(adore) %>% mutate(priority = 1L),
  species_groups(inv_tox_merge) %>% mutate(priority = 2L)
) %>%
  arrange(priority) %>%
  distinct(latin, .keep_all = TRUE) %>%
  select(latin, taxon_group)

species_crosscheck <- species_crosscheck %>%
  left_join(species_taxonomy, by = "latin")

species_group_summary <- species_crosscheck %>%
  count(status, taxon_group, name = "n_species") %>%
  pivot_wider(names_from = taxon_group, values_from = n_species,
              values_fill = 0)

crustacean_summary <- species_crosscheck %>%
  filter(taxon_group %in% c("Crustacean", "Noncrustacean")) %>%
  count(status, taxon_group, name = "n_species") %>%
  pivot_wider(names_from = taxon_group, values_from = n_species,
              values_fill = 0)

cat("\nSpecies overlap by taxonomic group:\n")
print(species_group_summary)
cat("\nInvertebrate overlap:\n")
print(crustacean_summary)

# ----- 4. Combine and save all three datasets --------------------------------
out_dir <- here("data/1_observed_sensitivity/combined_observed_sensitivity")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Keep all observations and their source-specific information.
combined_observed_sensitivity <- bind_rows(
  adore %>% mutate(dataset = "ADORE"),
  inv_tox_merge %>% mutate(dataset = "Liang–Sinclair")
) %>%
  add_group() %>%
  select(dataset, CAS, CAS_norm, latin_raw, latin, conc_log10mgL, duration,
         group, tax_group, kingdom, phylum, `class`, order, family, genus)

print(
  combined_observed_sensitivity %>%
    summarise(
      observations = n(),
      chemicals = n_distinct(CAS_norm, na.rm = TRUE),
      species = n_distinct(latin, na.rm = TRUE),
      genera = n_distinct(genus, na.rm = TRUE),
      chemical_species_pairs = nrow(
        distinct(
          filter(
            combined_observed_sensitivity,
            !is.na(CAS_norm), !is.na(latin), latin != ""
          ),
          CAS_norm, latin
        )
      )
    )
)

print(combined_observed_sensitivity %>%
        count(dataset, group, name = "observations"))

saveRDS(combined_observed_sensitivity, 
        file = file.path(out_dir, "combined_observed_toxicity.rds"))

# Machine learning input: CAS, latin, log10(mg/L), duration in hours.
observed_sensitivity <- combined_observed_sensitivity %>%
  transmute(
    CAS = CAS_norm,
    latin,
    conc = conc_log10mgL,
    duration
  )

write_csv(
  observed_sensitivity,
  file = file.path(out_dir, "combined_observed_toxicity.csv")
)

