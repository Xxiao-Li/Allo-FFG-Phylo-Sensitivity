clean_latin <- function(x) str_to_lower(str_squish(str_replace_all(x, "_", " ")))

clean_name <- function(x) {
  x <- str_to_lower(str_squish(as.character(x)))
  x <- str_replace_all(x, "\\(.*?\\)|\\[.*?\\]", " ")
  x <- str_squish(str_replace_all(x, "[^a-z -]", " "))
  parts <- str_match(x, "^([a-z]+) ([a-z][a-z-]+)(?: |$)")
  valid <- !is.na(parts[, 2]) & !is.na(parts[, 3]) &
    !parts[, 3] %in% c("sp", "spp", "cf", "aff", "var", "subsp", "ssp", "forma")
  if_else(valid, paste(parts[, 2], parts[, 3]), NA_character_)
}

name_fixes <- c(
  "lithodes antarcticus" = "lithodes antarctica",
  "viviparus bengalensis" = "filopaludina bengalensis",
  "acanthomysis sculpta" = "holmesimysis sculpta",
  "acartia clausi" = "acartia clausii",
  "acartia lillgeborgi" = "acartia lilljeborgii",
  "acheilognathus moriokae" = "acheilognathus melanogaster",
  "achnanthes minutissima" = "achnanthidium minutissimum",
  "nuria danrica" = "esomus danrica",
  "oreochromis mossambicus x urolepis" = "oreochromis mossambicus",
  "oreochromis mossambicus x niloticus" = "oreochromis mossambicus",
  "enteromorpha intestinalis" = "ulva intestinalis"
)

standardise_names <- function(df) {
  df %>% mutate(latin_raw = latin, latin = recode(latin, !!!name_fixes),
                latin = clean_name(latin), genus = str_extract(latin, "^[a-z]+"))
}

normalize_cas <- function(x) {
  digits <- str_remove_all(as.character(x), "[^0-9]")
  if_else(!is.na(digits) & str_detect(digits, "^\\d{4,}$"),
          str_replace(digits, "^(\\d+)(\\d{2})(\\d)$", "\\1-\\2-\\3"),
          NA_character_)
}

dataset_summary <- function(x) {
  tibble(
    chemicals = n_distinct(x$CAS_norm, na.rm = TRUE),
    species = n_distinct(x$latin[!is.na(x$latin) & x$latin != ""]),
    pairs = x %>% filter(!is.na(CAS_norm), !is.na(latin), latin != "") %>%
      distinct(CAS_norm, latin) %>% nrow()
  )
}

species_per_chemical <- function(x) {
  x %>% filter(!is.na(CAS_norm), !is.na(latin), latin != "") %>%
    distinct(CAS_norm, latin) %>%
    count(CAS_norm, name = "n_species", sort = TRUE)
}

add_group <- function(x) {
  x %>% mutate(
    group = case_when(
      str_to_lower(coalesce(phylum, "")) %in% c(
        "chlorophyta", "charophyta", "rhodophyta", "ochrophyta", "haptophyta",
        "cryptophyta", "cyanobacteria", "myzozoa", "euglenozoa") |
        str_detect(str_to_lower(coalesce(tax_group, "")), "alga|phytoplankton") ~ "algae",
      str_to_lower(coalesce(`class`, "")) %in% c(
        "actinopterygii", "cephalaspidomorphi", "chondrichthyes",
        "petromyzonti", "elasmobranchii", "myxini", "hyperoartia", "sarcopterygii") |
        str_detect(str_to_lower(coalesce(tax_group, "")), "fish") ~ "fish",
      str_to_lower(coalesce(kingdom, "")) == "animalia" |
        str_detect(str_to_lower(coalesce(tax_group, "")),
                   "invert|crusta|insect|mollusc|worm") ~ "invertebrate",
      TRUE ~ NA_character_
    )
  )
}

species_groups <- function(x) {
  add_group(x) %>%
    filter(!is.na(latin), latin != "") %>%
    mutate(
      class_lower = str_to_lower(coalesce(as.character(`class`), "")),
      taxon_group = case_when(
        group == "algae" ~ "Algae",
        group == "fish" ~ "Fish",
        group == "invertebrate" ~
          if_else(str_detect(str_to_lower(coalesce(tax_group, "")), "crusta") |
                    class_lower %in% c("malacostraca", "maxillopoda", "branchiopoda",
                                       "copepoda", "ostracoda"),
                  "Crustacean", "Noncrustacean"),
        TRUE ~ "Other or unclassified"
      )
    ) %>%
    distinct(latin, taxon_group)
}


tag_to_group <- c(
  # Arthropod control
  "Insecticide" = "arthropod-control",
  "Acaricide" = "arthropod-control",
  "Miticide" = "arthropod-control",
  "Larvicide" = "arthropod-control",
  "Termiticide" = "arthropod-control",
  "Insect Growth Regulator" = "arthropod-control",
  
  # Plant control
  "Herbicide" = "plant-control",
  "Plant Growth Regulator" = "plant-control",
  "Plant Growth Regulator precursor" = "plant-control",
  "Desiccant" = "plant-control",
  "Defoliant" = "plant-control",
  "Graminicide" = "plant-control",
  
  # Microorganism control
  "Fungicide" = "microbial-control",
  "Microbiocide" = "microbial-control",
  "Algicide" = "microbial-control",
  "bacteria" = "microbial-control",
  "bacteria, fungi" = "microbial-control",
  
  # Nematode control
  "Nematicide" = "nematode-control",
  
  # Mollusc control
  "Molluscicide" = "mollusc-control",
  
  # Fish / vertebrate control
  "Rodenticide" = "fish-vertebrate-control",
  "Avicide" = "fish-vertebrate-control",
  "Piscicide" = "fish-vertebrate-control",
  "Lampricide" = "fish-vertebrate-control",
  
  # Veterinary antiparasitic / drug
  "Veterinary substance" = "veterinary-drug",
  "Crop antibiotic" = "veterinary-drug",
  
  # Repellent / chemical signals
  "Repellent" = "repellent-signal",
  "Semiochemical" = "repellent-signal",
  "Attractant" = "repellent-signal",
  
  # Fumigant / soil sterilant
  "Fumigant" = "fumigant-soil-sterilant",
  "Soil sterilant" = "fumigant-soil-sterilant",
  
  # Auxiliary groups
  "Metabolite" = "auxiliary",
  "Adjuvant" = "auxiliary"
)

map_subgroup_to_target <- function(x, tag_to_group) {
  
  if (is.na(x) || str_squish(x) == "") return(NA_character_)
  
  tags <- str_split(x, pattern = ";", simplify = FALSE)[[1]] %>%
    str_squish()
  
  tags <- tags[tags != "Other substance"]
  
  if (length(tags) == 0) return(NA_character_)
  if (length(tags) == 1 && identical(tags, "Mixed")) return("mixed-target")
  
  groups_all <- unique(na.omit(unname(tag_to_group[tags])))
  main_groups <- setdiff(groups_all, "auxiliary")
  has_auxiliary <- "auxiliary" %in% groups_all
  
  if (length(main_groups) == 0 && !has_auxiliary) return(NA_character_)
  if (length(main_groups) == 0 && has_auxiliary) return("auxiliary")
  if (length(main_groups) == 1) return(main_groups[1])
  if (length(main_groups) > 1) return("mixed-target")
  
  return(NA_character_)
}