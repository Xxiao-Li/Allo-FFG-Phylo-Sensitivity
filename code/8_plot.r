#==============================================================================
# 0_plot.R
# - Figure 2: Full tree + six chemical-use groups
# - Figure 3: Phylogenetic signal by taxonomic group and chemical-use group
# - Figure 4: Leave-one-out model comparison
# - Figure 5: Standardised FFG-specific slopes
# - Figure 6: Pairwise FFG-specific slope comparison
# - Figure 7: Pairwise FFG-specific intercept comparison
# - Figure 8: Model explanatory power and random-effect variance partitioning
# - Figures S3-S16
#==============================================================================

rm(list = ls())

library(here)
library(tidyverse)
library(ape)
library(ggtree)
library(patchwork)
library(scales)
library(ggbeeswarm)
library(paletteer)
library(car)
library(emmeans)
library(brms)
library(ggdist)
library(glue)
library(paletteer)

#=======1. Load data and general settings=======================================
algae_df <- readRDS(here("data/6_combined_meta_data/algae_meta_data.rds"))
invertebrate_df <- readRDS(here("data/6_combined_meta_data/invertebrate_meta_data.rds"))
fish_df <- readRDS(here("data/6_combined_meta_data/fish_meta_data.rds"))

full_df <- readRDS(here("data/6_combined_meta_data/full_meta_data.rds"))
full_tree <- readRDS(here("data/6_combined_meta_data/full_genus_tree.rds"))

dir.create(here("figures"), showWarnings = FALSE, recursive = TRUE)

chem_gr_order <- c(
  "Industrial_Chemical",
  "Pharmaceutical",
  "P_Arthropod_control",
  "P_Plant_control",
  "P_Microbial_control",
  "P_Mixed_target"
)

chem_gr_label <- c(
  Industrial_Chemical   = "Industrial chemical",
  Pharmaceutical        = "Pharmaceutical",
  P_Arthropod_control   = "P: Arthropods",
  P_Plant_control       = "P: Plants",
  P_Microbial_control   = "P: Microbes",
  P_Mixed_target        = "P: Mixed"
)
chem_gr_label_2row <- c(
  Industrial_Chemical   = "Industrial\nchemical",
  Pharmaceutical        = "Pharmaceutical",
  P_Arthropod_control   = "P: Arthropods",
  P_Plant_control       = "P: Plants",
  P_Microbial_control   = "P: Microbes",
  P_Mixed_target        = "P: Mixed"
)
chem_gr_label_order_2row <- unname(chem_gr_label_2row[chem_gr_order])

tax_group_label <- c("Algae", "Invertebrate", "Fish")

tax_group_order <- c("Algae", "Invertebrate", "Fish")

tax_cols <- c(
  Algae        = "#1b9e77",
  Invertebrate = "#d95f02",
  Fish         = "#7570b3"
)

ffg_label_vec <- c(
  "Greenalgae"             = "Green algae",
  "Diatom"                  = "Diatom",
  "Otheralgae"             = "Other algae",
  
  "I_Filterer_A"            = "Arthropod filterer",
  "I_Predator_A"            = "Arthropod predator",
  "I_Shredder_A"            = "Arthropod shredder",
  "I_Gatherer_A"            = "Arthropod gatherer",
  "I_Herbivore_A"           = "Arthropod herbivore",
  
  "I_Gatherer"              = "Non-arthropod gatherer",
  "I_Filterer"              = "Non-arthropod filterer",
  "I_Predator"              = "Non-arthropod predator",
  "I_Herbivore"             = "Non-arthropod herbivore",
  
  "F_Invertivore"           = "Invertivore",
  "F_Piscivore"             = "Piscivore",
  "F_Omnivore"              = "Omnivore",
  "F_Detritivore_Herbivore" = "Detritivore/herbivore",
  "F_Planktivore"           = "Planktivore"
)
ffg_order_algae <- c("Diatom", "Green algae", "Other algae")

ffg_order_invertebrate <- c(
  "Arthropod filterer",
  "Arthropod gatherer",
  "Arthropod herbivore",
  "Arthropod predator",
  "Arthropod shredder",
  "Non-arthropod filterer",
  "Non-arthropod gatherer",
  "Non-arthropod herbivore",
  "Non-arthropod predator"
)

ffg_order_fish <- c(
  "Detritivore/herbivore",
  "Invertivore",
  "Omnivore",
  "Piscivore",
  "Planktivore"
)

ffg_order_plot <- rev(unique(c(
  ffg_order_algae,
  ffg_order_invertebrate,
  ffg_order_fish
)))

pick_mode <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_character_)
  names(sort(table(x), decreasing = TRUE))[1]
}

standardise_tax_group <- function(x) {
  x_clean <- as.character(x) %>%
    stringr::str_to_lower() %>%
    stringr::str_trim() %>%
    stringr::str_replace_all("[[:space:]]+", " ")
  
  dplyr::case_when(
    x_clean %in% c("algae") ~ "Algae",
    x_clean %in% c("invertebrate", "invertebrates") ~ "Invertebrate",
    x_clean %in% c("fish", "vertebrate", "vertebrates") ~ "Fish",
    x_clean %in% c("all taxa", "all_taxa", "all-taxa", "alltaxa", "full") ~ "All taxa",
    TRUE ~ NA_character_
  )
}

plot_theme <- function(fs) {
  theme_bw(base_size = fs) +
    theme(
      panel.grid.minor = element_blank(),
      axis.title = element_text(color = "black", size = fs),
      axis.text = element_text(color = "black", size = fs),
      axis.ticks = element_line(color = "black"),
      legend.title = element_text(color = "black", size = fs),
      legend.text = element_text(color = "black", size = fs),
      strip.background = element_rect(fill = "white", color = "black"),
      strip.text = element_text(size = fs, face = "bold")
    )
}


#=======Figure 2: Full tree + six chemical-use groups===========================
x_break_map <- list(
  Industrial_Chemical = c(-1, 0, 1, 2),
  Pharmaceutical = c(-1, 0, 1, 2),
  P_Arthropod_control = c(-1, 0, 1, 2),
  P_Plant_control = c(-1, 0, 1, 2),
  P_Microbial_control = c(-1, 0, 1, 2),
  P_Mixed_target = c(-1, 0, 1, 2)
)

fs <- 18
fig2_dat <- full_df %>%
  mutate(
    tax_group = case_when(
      tax_group == "algae" ~ "Algae",
      tax_group == "invertebrate" ~ "Invertebrate",
      tax_group == "fish" ~ "Fish",
      TRUE ~ NA_character_),
    chem_group = factor(chem_group, levels = chem_gr_order),
    tax_group  = factor(tax_group, levels = tax_group_order),
    genus = str_to_lower(str_trim(as.character(genus))) ) %>%
  droplevels()

# Genus-level taxonomic group
sp_traits <- fig2_dat %>%
  filter(!is.na(tax_group)) %>%
  group_by(genus) %>%
  summarise( tax_group = pick_mode(as.character(tax_group)),
             .groups = "drop") %>%
  mutate(tax_group = factor(tax_group, levels = tax_group_order))

# Make sure colours match taxonomic group names
tax_cols <- tax_cols[tax_group_order]
names(tax_cols) <- tax_group_order

# OTU groups for branch colouring
otus <- split(as.character(sp_traits$genus),
              as.character(sp_traits$tax_group))

otus <- otus[lengths(otus) > 0]

tree_g <- groupOTU(full_tree, otus, group_name = "tax_group")

p_tree <- ggtree(tree_g, aes(colour = tax_group), linewidth = 0.85) +
  theme_tree() +
  scale_colour_manual(values = tax_cols, breaks = tax_group_order,
    name = "Taxonomic group", drop = FALSE,
    na.value = "lightgrey") +
  theme(plot.margin = margin(2, 3, 2, 2),
        legend.position = "bottom")

p_tree

# Extract tree tip y-positions for alignment with right panels
tip_y <- p_tree$data %>%
  filter(isTip) %>%
  transmute(genus = label, y = y)
y_rng <- range(tip_y$y, na.rm = TRUE)

# Summarise genus-level sensitivity
summarise_species_sensitivity <- function(d, cg) {
  genus_summary <- d %>%
    filter(chem_group == cg) %>%
    group_by(genus) %>%
    summarise(
      med  = median(log_rel_sensi, na.rm = TRUE),
      q25  = quantile(log_rel_sensi, 0.25, na.rm = TRUE, names = FALSE),
      q75  = quantile(log_rel_sensi, 0.75, na.rm = TRUE, names = FALSE),
      lo95 = quantile(log_rel_sensi, 0.025, na.rm = TRUE, names = FALSE),
      hi95 = quantile(log_rel_sensi, 0.975, na.rm = TRUE, names = FALSE),
      n_cas = n_distinct(CAS),
      .groups = "drop")
  
  tibble(genus = full_tree$tip.label) %>%
    mutate(genus = str_to_lower(str_trim(as.character(genus)))) %>%
    left_join(sp_traits, by = "genus") %>%
    left_join(genus_summary, by = "genus") %>%
    left_join(tip_y, by = "genus") %>%
    mutate(tax_group = factor(tax_group, levels = tax_group_order)) %>%
    filter(!is.na(y))
}

make_right_panel <- function(cg) {
  
  df <- summarise_species_sensitivity(fig2_dat, cg)
  
  xb <- x_break_map[[cg]]
  if (is.null(xb)) xb <- c(-1, 0, 1)
  
  xlim_i <- range(xb, na.rm = TRUE)
  
  grp_q <- df %>%
    group_by(tax_group) %>%
    summarise(
      med_g  = median(med, na.rm = TRUE),
      q25_g  = quantile(med, 0.25, na.rm = TRUE, names = FALSE),
      q75_g  = quantile(med, 0.75, na.rm = TRUE, names = FALSE),
      q025_g = quantile(med, 0.025, na.rm = TRUE, names = FALSE),
      q975_g = quantile(med, 0.975, na.rm = TRUE, names = FALSE),
      y_min  = min(y, na.rm = TRUE),
      y_max  = max(y, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      y_range = y_max - y_min,
      
      # keep label clearly inside each taxonomic block
      y_lab = case_when(
        tax_group == "Algae" ~ y_min + 0.40 * y_range,
        TRUE ~ y_min + 0.50 * y_range
      ),
      
      x_lab = xlim_i[1] + 0.12 * diff(xlim_i),
      lab = sprintf("%.2f [%.2f, %.2f]", med_g, q25_g, q75_g)
    )
  
  ggplot(df, aes(y = y)) +
    
    # 95% range across genus-level median sensitivity
    geom_rect(data = grp_q, 
      aes(xmin = q025_g, xmax = q975_g,
          ymin = y_min, ymax = y_max,
          fill = tax_group), 
      inherit.aes = FALSE, alpha = 0.3, color = NA) +
    # IQR across genus-level median sensitivity
    geom_rect(data = grp_q, 
      aes(xmin = q25_g, xmax = q75_g,
          ymin = y_min, ymax = y_max,
          fill = tax_group), 
      inherit.aes = FALSE, alpha = 0.7, color = NA) +
    # Median [IQR] label
    geom_text(data = grp_q,
      aes(x = x_lab, y = y_lab, label = lab, color = tax_group),
      inherit.aes = FALSE, angle = 90, size = 5, fontface = "bold",
      hjust = 0.3, vjust = 0.1, lineheight = 0.85) +
    # Genus-level IQR
    geom_segment(
      aes(x = q25, xend = q75, yend = y),
      linewidth = 0.25, color = "grey70", alpha = 0.85) +
    # Genus-level median
    geom_point(aes(x = med),
      size = 0.45, color = "grey65", alpha = 0.85) +
   # geom_vline(xintercept = 0, linewidth = 0.25, color = "grey60", 
    #           linetype = "dashed") +
    scale_fill_manual(values = tax_cols, breaks = tax_group_order,
                      guide = "none") +
    scale_color_manual(values = tax_cols, breaks = tax_group_order, 
                       guide = "none" ) +
    scale_y_continuous(limits = y_rng, expand = expansion(mult = c(0.01, 0.01))) +
    scale_x_continuous(breaks = xb, labels = as.character(xb), 
                       expand = expansion(mult = c(0.04, 0.04))) +
    coord_cartesian(xlim = xlim_i, clip = "off") +
    labs(title = chem_gr_label[[cg]], 
         x = NULL, y = NULL) +
    theme_minimal(base_size = fs) +
    theme(
      plot.title = element_text(
        size = fs - 2, face = "bold", hjust = 0.5, margin = margin(b = 4) ),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(2, 1, 8, 1), 
      axis.line.x = element_line(color = "gray70", linewidth = 0.3),
      
      axis.title = element_text(size = fs + 2, color = "black"),
      axis.text = element_text(size = fs, color = "black"),
      
      # move x-axis tick labels lower
      axis.text.x = element_text(size = fs , color = "black", margin = margin(t = -10)),
      
      legend.title = element_text(size = fs, color = "black"),
      legend.text = element_text(size = fs, color = "black")
    )
    
}

right_panels <- map(chem_gr_order, make_right_panel)

right_panels <- map2(right_panels, seq_along(right_panels), \(p, i) {
  if (i == 3) {
    p +
      labs(x = expression(log[10] * "Relative sensitivity")) +
      theme(axis.title.x = element_text(size = fs, margin = margin(t = 16)))
  } else {
    p +
      labs(x = NULL) +
      theme(axis.title.x = element_blank())
  }
})

right_with_xlab <- wrap_plots(right_panels, nrow = 1)

fig2 <- p_tree + right_with_xlab +
  plot_layout(widths = c(1, 5.2), guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(size = fs),
    legend.text = element_text(size = fs),
    legend.key.height = unit(4, "mm"),
    legend.key.width = unit(10, "mm"),
    plot.margin = margin(4, 12, 4, 4)
  )

fig2

ggsave(filename = here("figures", "Figure2.jpeg"),
       plot = fig2, width = 16, height = 10, dpi = 600)

ggsave(filename = here("figures", "Figure2.pdf"),
       plot = fig2, width = 16, height = 10, dpi = 600)


# Summaries for revising Fig. 2 text
fig2_summary <- map_dfr(chem_gr_order, \(cg) {
  
  summarise_species_sensitivity(fig2_dat, cg) %>%
    filter(!is.na(med), !is.na(tax_group)) %>%
    group_by(tax_group) %>%
    summarise(
      n_genus = n_distinct(genus),
      median = median(med, na.rm = TRUE),
      q25    = quantile(med, 0.25,  na.rm = TRUE, names = FALSE),
      q75    = quantile(med, 0.75,  na.rm = TRUE, names = FALSE),
      q025   = quantile(med, 0.025, na.rm = TRUE, names = FALSE),
      q975   = quantile(med, 0.975, na.rm = TRUE, names = FALSE),
      .groups = "drop"
    ) %>%
    mutate(
      chem_group = cg,
      chem_group_label = chem_gr_label[[cg]],
      .before = 1
    )
})

fig2_summary



#=======Figure 3: Phylogenetic signal within each taxonomic group and chemical-use group=======
fits_dir <- here("data/7_phylogenetic_signal_by_tax_group")
tax_group_order_fig3 <- c("Algae", "Invertebrate", "Fish")

plot_df <- list.files(fits_dir,
                        pattern = "^phytools_lambda_by_CAS_.*\\.rds$",
                        full.names = TRUE) %>%
  map_dfr(readRDS) %>%
  mutate(
    tax_group_raw = tax_group,
    tax_group = standardise_tax_group(tax_group),
    tax_group = factor(tax_group, levels = tax_group_order_fig3)) %>%
  filter(tax_group_raw %in% c("algae", "invertebrate", "fish"))

m_lm <- lm(lambda ~ tax_group * chem_group, data = plot_df)
anova_lm <- car::Anova(m_lm, type = 3)
emm_chem_by_tax <- emmeans(m_lm, ~ chem_group | tax_group)

chem_letters_by_tax <- emmeans:::cld.emmGrid(
  emm_chem_by_tax,
  adjust = "BH",
  Letters = letters,
  sort = TRUE,
  reversed = TRUE
) %>%
  as.data.frame() %>%
  mutate(
    .group = str_replace_all(.group, " ", ""),
    chem_group = factor(chem_group, levels = chem_gr_order),
    tax_group = factor(
      tax_group,
      levels = c("Algae", "Invertebrate", "Fish"))
  )

letter_pos <- chem_letters_by_tax %>%
  mutate(x_pos = 0.86)

pal_raw <- as.character(paletteer::paletteer_d("ggthemes::excel_Vapor_Trail"))

pal <- setNames(pal_raw[seq_along(chem_gr_order)], chem_gr_order)

fs <- 18

fig3 <- ggplot(plot_df, aes(x = lambda, y = chem_group, 
                            color = tax_group)) +
  geom_boxplot(width = 0.55, alpha = 0.20, 
               outlier.shape = NA, linewidth = 0.9) +
  ggbeeswarm::geom_quasirandom(size = 1.6, alpha = 0.60, width = 0.18) +
  geom_text(
    data = letter_pos,
    aes(x = x_pos, y = chem_group, label = .group),
    inherit.aes = FALSE, hjust = 0, vjust = 0.5, size = 7, color = "black") +
  facet_wrap(~ tax_group, nrow = 1, drop = TRUE) +
  scale_x_continuous(
    breaks = scales::breaks_extended(n = 3),
    limits = c(0, 1),
    expand = expansion(mult = c(0.02, 0.12))) +
  scale_y_discrete(
    labels = chem_gr_label[chem_gr_order],
    drop = FALSE) +
  scale_colour_manual(
    values = tax_cols[tax_group_order_fig3],
    breaks = tax_group_order_fig3,
    drop = FALSE) +
  labs(x = "Phylogenetic signal\n(Pagel's \u03BB)", y = NULL) +
  theme_bw(base_size = fs) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "none",
    axis.title = element_text(color = "black", size = fs),
    axis.text = element_text(color = "black", size = fs),
    axis.ticks = element_line(color = "black"),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = fs, face = "bold"),
    panel.spacing.x = unit(1.2, "cm")
  )

fig3

ggsave(filename = here("figures", "Figure3.jpeg"),
       plot = fig3, width = 14, height = 8, dpi = 600)

ggsave(filename = here("figures", "Figure3.pdf"),
       plot = fig3, width = 14, height = 8, dpi = 600)

# Reported summaries in the text
print(anova_lm)

# entire range
round(range(plot_df$lambda, na.rm = TRUE),2) 


#median and IQR for each scenario
plot_df %>%
  group_by(tax_group) %>%
  summarise(
    median_lambda = round(median(lambda, na.rm = TRUE), 4),
    q25_lambda = round(quantile(lambda, 0.25, na.rm = TRUE, names = FALSE), 4),
    q75_lambda = round(quantile(lambda, 0.75, na.rm = TRUE, names = FALSE), 4),
    .groups = "drop"
  ) %>%
  arrange(tax_group)

signal_full_summary <- plot_df %>%
  group_by(tax_group, chem_group) %>%
  summarise(
    median_lambda = round(median(lambda, na.rm = TRUE), 4),
    q25_lambda = round(quantile(lambda, 0.25, na.rm = TRUE, names = FALSE), 4),
    q75_lambda = round(quantile(lambda, 0.75, na.rm = TRUE, names = FALSE), 4),
    .groups = "drop"
  ) %>%
  arrange(tax_group, chem_group)

print(signal_full_summary)

# Post-hoc comparisons: chemical groups within taxonomic groups
emm_chem_by_tax <- emmeans(m_lm, ~ chem_group | tax_group)
pairs(emm_chem_by_tax, adjust = "BH")

#=======Figure 4: Leave-one-out model comparison===============================
taxon_order <- c("algae", "invertebrate", "fish")
taxon_label_vec <- c(algae = "Algae", invertebrate = "Invertebrate", fish = "Fish")
chem_gr_label_vec <- chem_gr_label

fig4_data <- read_csv(here("data/8_bayesian_model_fits/bayesian_model_inference/LOO_model_comparisons.csv")) %>%
  mutate(
    taxon = factor(taxon, levels = taxon_order),
    chem_gr = factor(chem_gr, levels = chem_gr_order),
    hypothesis = factor(hypothesis, levels = c("H0", "H1", "H2", "H3")),
    model_id = factor(model_id, levels = c("M0", "M1", "M2", "M3")),
    taxon_label_plot = factor(
      taxon_label_vec[as.character(taxon)],
      levels = unname(taxon_label_vec[taxon_order])),
    chem_gr_label_plot = factor(
      chem_gr_label_vec[as.character(chem_gr)],
      levels = unname(chem_gr_label_vec[chem_gr_order])
    )
  )

fs <- 18

fig4 <- ggplot(fig4_data, aes(elpd_diff, model_id, colour = taxon_label_plot)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, colour = "grey70") +
  geom_errorbar(aes(xmin = elpd_diff - se_diff, xmax = elpd_diff + se_diff),
                orientation = "y", width = 0.25, linewidth = 1.5) +
  geom_point(size = 5) +
  facet_grid(taxon_label_plot ~ chem_gr_label_plot, drop = FALSE, scales = "free_x",
             labeller = labeller(taxon_label_plot = label_wrap_gen(12),
                                 chem_gr_label_plot = label_wrap_gen(16))) +
  scale_colour_manual(
    values = tax_cols[unname(taxon_label_vec[taxon_order])],
    guide = "none") +
  scale_x_continuous(breaks = scales::breaks_extended(n = 3)) +
  scale_y_discrete(
    labels = c(
      M0 = "M0",
      M1 = "M1",
      M2 = "M2",
      M3 = "M3")
  ) +
  labs(x = expression(bold(Delta * ELPD~"relative to the best-supported model")),
       y = "Candidate model") +
  theme_bw(base_size = fs) +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = fs, face = "bold", colour = "black"),
    axis.text = element_text(size = fs - 2, colour = "black"),
    legend.text = element_text(size = fs, colour = "black"),
    strip.text = element_text(size = fs, face = "bold", lineheight = 0.9),
    legend.position = "bottom",
    panel.spacing.x = unit(0.6, "cm")
  )

fig4

ggsave(here("figures/Figure4.jpeg"),
       fig4, width = 15, height = 9, dpi = 600)
ggsave(here("figures/Figure4.pdf"),
       fig4, width = 15, height = 9, dpi = 600)

# summary in the text
loo_summary <- fig4_data %>%
  group_by(taxon, chem_gr) %>%
  summarise(
    best_model = model_id[which.max(elpd_diff)],
    best_hypothesis = hypothesis[which.max(elpd_diff)],
    elpd_diff_best = max(elpd_diff),
    se_diff_best = se_diff[which.max(elpd_diff)],
    .groups = "drop"
  ) %>%
  mutate(
    taxon_label = taxon_label_vec[as.character(taxon)],
    chem_gr_label = chem_gr_label_vec[as.character(chem_gr)],
    .before = 1
  )
loo_summary

#=======Figure 5: Standardised FFG-specific slopes=============================
rope_cols <- c(
  "Inside ±0.10"   = "#1b9e77",
  "Inside ±0.20"   = "#66a61e",
  "Overlaps ±0.20" = "grey45",
  "Outside ±0.20"  = "#d95f02"
)

fig_slope_data <- read.csv(here(
  "data", "8_bayesian_model_fits", "bayesian_model_inference",
  "standardised_FFG_slope_summary.csv")) %>%
  mutate(
    taxon_label_plot = factor(standardise_tax_group(taxon), levels = tax_group_order),
    chem_gr = factor(chem_gr, levels = chem_gr_order),
    chem_gr_label_plot = factor(chem_gr_label_2row[as.character(chem_gr)],
                                levels = chem_gr_label_order_2row),
    FFG_plot = forcats::fct_rev(forcats::fct_inorder(FFG)),
    hdi_rope_status = case_when(
      hdi_inside_0.10  ~ "Inside ±0.10",
      hdi_inside_0.20  ~ "Inside ±0.20",
      hdi_outside_0.20 ~ "Outside ±0.20",
      TRUE ~ "Overlaps ±0.20"),
    hdi_rope_status = factor(
      hdi_rope_status,
      levels = c("Inside ±0.10", "Inside ±0.20",
                 "Overlaps ±0.20", "Outside ±0.20")
    )
  )
fig_slope_data <- fig_slope_data %>%
  mutate(
    FFG_plot = dplyr::recode(FFG, !!!ffg_label_vec, .default = FFG),
    FFG_plot = factor(FFG_plot, levels = ffg_order_plot)
  )

fs <- 20

fig5 <- ggplot(fig_slope_data, aes(beta_std_median, FFG_plot)) +
  annotate("rect", xmin = -0.20, xmax = 0.20, ymin = -Inf, ymax = Inf,
           fill = "grey60", alpha = 0.22) +
  # annotate("rect", xmin = -0.10, xmax = 0.10, ymin = -Inf, ymax = Inf,
  #          fill = "grey45", alpha = 0.22) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, colour = "grey35") +
  geom_errorbar(aes(xmin = hdi_low, xmax = hdi_high, colour = hdi_rope_status),
                orientation = "y", width = 0.28, linewidth = 1.4) +
  geom_point(aes(colour = hdi_rope_status), size = 4) +
  facet_grid(taxon_label_plot ~ chem_gr_label_plot, scales = "free", space = "free_y") +
  scale_colour_manual(values = rope_cols, name = "95% HDI vs ROPE") +
  scale_x_continuous(breaks = scales::breaks_pretty(n = 3),
                     labels = scales::label_number(accuracy = 0.1)) +
  labs(
    x = expression("Standardised FFG-specific slope (" * beta[std] * ")"),
    y = NULL
  ) +
  theme_bw(base_size = fs) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.title = element_text(size = fs, face = "bold", colour = "black"),
    axis.text.x = element_text(size = fs, colour = "black"),
    axis.text.y = element_text(size = fs, colour = "black"),
    strip.text = element_text(size = fs-2, face = "bold", colour = "black"),
    strip.background = element_rect(fill = "grey90", colour = "grey40"),
    legend.position = "bottom",
    legend.title = element_text(size = fs, face = "bold"),
    legend.text = element_text(size = fs),
    panel.spacing.x = unit(0.35, "cm"),
    panel.spacing.y = unit(0.35, "cm")
  )

fig5

ggsave(here("figures", "Figure5.jpeg"),
       fig5, width = 16, height = 12, dpi = 600)

ggsave(here("figures", "Figure5.pdf"),
       fig5, width = 16, height = 12, dpi = 600)

# summaries used in text in result 
slope_rope_summary <- fig_slope_data %>%
  group_by(taxon_label_plot, chem_gr_label_plot, hdi_rope_status) %>%
  summarise(
    n = n(),
    .groups = "drop_last"
  ) %>%
  mutate(
    prop = n / sum(n)
  ) %>%
  ungroup()

write_csv(slope_rope_summary, here("data", "8_bayesian_model_fits", "bayesian_model_inference",
                                "slope_hdi_rope_summary.csv"))

# Overall number of slopes outside practical-effect threshold
fig_slope_data %>%
  summarise(
    n_slopes = n(),
    n_inside_010 = sum(hdi_inside_0.10, na.rm = TRUE),
    n_inside_020 = sum(hdi_inside_0.20, na.rm = TRUE),
    n_outside_020 = sum(hdi_outside_0.20, na.rm = TRUE),
    prop_outside_020 = n_outside_020 / n_slopes
  )

# Slopes with 95% HDI entirely outside ±0.20
fig_slope_data %>%
  filter(hdi_outside_0.20) %>%
  arrange(taxon_label_plot, chem_gr, desc(abs(beta_std_median))) %>%
  select(taxon_label_plot, chem_gr_label_plot, FFG_plot, beta_std_median, 
         hdi_low, hdi_high, hdi_rope_status)

#=======Figure 6: Pairwise FFG-specific slope comparison heatmap================
fig_slope_heatmap_data <- read.csv(here("data", "8_bayesian_model_fits", "bayesian_model_inference",
                                        "pairwise_FFG_slope_contrast_summary.csv")) %>%
  mutate(
    taxon_label_plot = factor(standardise_tax_group(taxon), levels = tax_group_order),
    chem_gr = factor(chem_gr, levels = chem_gr_order),
    chem_gr_label_plot = factor(chem_gr_label_2row[as.character(chem_gr)],
                                levels = unname(chem_gr_label_2row[chem_gr_order])),
    FFG_1_plot = dplyr::recode(str_trim(as.character(FFG_1)), 
                               !!!ffg_label_vec, .default = FFG_1),
    FFG_2_plot = dplyr::recode(str_trim(as.character(FFG_2)), 
                               !!!ffg_label_vec, .default = FFG_2),
    contrast_label = paste(FFG_1_plot, "\u2212", FFG_2_plot),
    hdi_rope_star = case_when(hdi_outside_0.50 ~ "**", 
                              hdi_outside_0.20 ~ "*", 
                              TRUE ~ "")
  )

fs <- 20

fig6 <- ggplot(fig_slope_heatmap_data,
               aes(x = chem_gr_label_plot, y = contrast_label, 
                   fill = slope_diff_median)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  geom_text(aes(label = hdi_rope_star), hjust = 0.5, vjust = 0.7, 
            size = 9, colour = "black") +
  facet_grid(taxon_label_plot ~ ., scales = "free_y", 
             space = "free_y", switch = "y") +
  scale_fill_gradient2(midpoint = 0, 
                       name = expression(atop("Median", Delta * beta[std]))) +
  labs(x = NULL, y = NULL,
       caption = "* 95% HDI entirely outside ±0.20;  ** 95% HDI entirely outside ±0.50") +
  theme_bw(base_size = fs) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 40, hjust = 1, colour = "black", size = fs - 1),
    axis.text.y = element_text(colour = "black", size = fs - 2),
    axis.title = element_text(face = "bold", colour = "black", size = fs),
    strip.placement = "outside",
    strip.text.y.left = element_text(face = "bold", angle = 90, size = fs, colour = "black"),
    strip.background = element_rect(fill = "grey85", colour = "grey30"),
    plot.caption = element_text(hjust = 1, size = fs, colour = "black"),
    legend.position = "right",
    legend.title = element_text(size = fs, face = "bold"),
    legend.text = element_text(size = fs - 1)
  )

fig6

ggsave(here("figures", "Figure6.jpeg"), 
       fig6, width = 14, height = 15, dpi = 600)

ggsave(here("figures", "Figure6.pdf"), 
       fig6, width = 14, height = 15, dpi = 600)

#=======Figure 7: Pairwise FFG-specific intercept comparison heatmap============
fig_intercept_heatmap_data <- read.csv(here("data", "8_bayesian_model_fits", "bayesian_model_inference",
                                            "pairwise_FFG_intercept_contrast_summary.csv")) %>%
  mutate(
    taxon_label_plot = factor(standardise_tax_group(taxon), 
                              levels = tax_group_order),
    chem_gr = factor(chem_gr, levels = chem_gr_order),
    chem_gr_label_plot = factor(chem_gr_label_2row[as.character(chem_gr)],
                                levels = unname(chem_gr_label_2row[chem_gr_order])),
    FFG_1_plot = dplyr::recode(str_trim(as.character(FFG_1)), 
                               !!!ffg_label_vec, .default = FFG_1),
    FFG_2_plot = dplyr::recode(str_trim(as.character(FFG_2)), 
                               !!!ffg_label_vec, .default = FFG_2),
    contrast_label = paste(FFG_1_plot, "\u2212", FFG_2_plot),
    hdi_rope_star = case_when(
      hdi_outside_0.80 ~ "***",
      hdi_outside_0.50 ~ "**",
      hdi_outside_0.20 ~ "*",
      TRUE ~ ""
    )
  )

fs <- 20

fig7 <- ggplot(fig_intercept_heatmap_data,
               aes(x = chem_gr_label_plot, y = contrast_label, fill = d_baseline_median)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  geom_text(aes(label = hdi_rope_star), hjust = 0.5, vjust = 0.7, 
            size = 9, colour = "black") +
  facet_grid(taxon_label_plot ~ ., scales = "free_y", 
             space = "free_y", switch = "y") +
  scale_fill_gradient2(midpoint = 0, name = expression(atop("Median", italic(d)))) +
  labs(x = NULL, y = NULL,
       caption = "* 95% HDI entirely outside ±0.20;  ** 95% HDI entirely outside ±0.50;  *** 95% HDI entirely outside ±0.80") +
  theme_bw(base_size = fs) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 40, hjust = 1, colour = "black", size = fs - 1),
    axis.text.y = element_text(colour = "black", size = fs - 2),
    axis.title = element_text(face = "bold", colour = "black", size = fs),
    strip.placement = "outside",
    strip.text.y.left = element_text(face = "bold", angle = 90, size = fs, colour = "black"),
    strip.background = element_rect(fill = "grey85", colour = "grey30"),
    plot.caption = element_text(hjust = 0.8, size = fs, colour = "black"),
    legend.position = "right",
    legend.title = element_text(size = fs, face = "bold"),
    legend.text = element_text(size = fs - 1)
  )

fig7

ggsave(here("figures", "Figure7.jpeg"), 
       fig7, width = 14, height = 15, dpi = 600)

ggsave(here("figures", "Figure7.pdf"), 
       fig7, width = 14, height = 15, dpi = 600)


#=======Figure 8: Model explanatory power and random-effect variance partitioning ============
r2_marginal_conditional_best <- read_csv(here(
  "data/8_bayesian_model_fits/bayesian_model_inference/r2_marginal_conditional_best_model.csv"))

random_variance_partition_best <- read_csv(here(
  "data/8_bayesian_model_fits/bayesian_model_inference/random_variance_partition_best_model.csv"))

chem_levels <- rev(chem_gr_order)
chem_labels <- rev(unname(chem_gr_label[chem_gr_order]))

fs <- 18
theme_fig8 <- plot_theme(fs) +
  theme(
    panel.grid.major.y = element_blank(),
    legend.position = "bottom",
    legend.box = "horizontal",
    axis.title = element_text(size = fs),
    axis.text = element_text(size = fs),
    strip.background = element_rect(fill = "grey90", colour = "grey40"),
    strip.text = element_text(size = fs, face = "bold"),
    plot.title = element_text(size = fs, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = fs)
  )

dodge_fig8 <- position_dodge(width = 0.55)

r2_plot_data <- r2_marginal_conditional_best %>%
  mutate(
    taxon_label_plot = factor(
      standardise_tax_group(taxon),
      levels = tax_group_order),
    chem_gr_label_plot = factor(
      chem_gr,
      levels = chem_levels,
      labels = chem_labels),
    r2_type = case_when(
      str_detect(r2_type, regex("marginal", ignore_case = TRUE)) ~ "Marginal R²",
      str_detect(r2_type, regex("conditional", ignore_case = TRUE)) ~ "Conditional R²",
      TRUE ~ as.character(r2_type)),
    r2_type = factor(r2_type, levels = c("Marginal R²", "Conditional R²"))
  )


p_r2 <- ggplot(r2_plot_data, aes(Estimate, chem_gr_label_plot,
                                 colour = taxon_label_plot, shape = r2_type)) +
  geom_errorbar(aes(xmin = Q2.5, xmax = Q97.5), orientation = "y",
                height = 0.45, linewidth = 1, position = dodge_fig8) +
  geom_point(size = 5, position = dodge_fig8) +
  facet_grid(taxon_label_plot ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_colour_manual(values = tax_cols, guide = "none") +
  scale_shape_manual(values = c("Marginal R²" = 16, "Conditional R²" = 17),
                     name = "Bayesian R²") +
  scale_x_continuous(breaks = seq(0, 1, 0.25), limits = c(0, 1),
                     expand = expansion(mult = c(0.02, 0.03))) +
  labs(x = expression(R^2), y = NULL,
       title = "(a) Marginal and conditional Bayesian R²") +
  theme_fig8 +
  theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 90))

variance_plot_data <- random_variance_partition_best %>%
  filter(!str_detect(group, regex("residual", ignore_case = TRUE))) %>%
  mutate(
    taxon_label_plot = factor(
      standardise_tax_group(taxon),
      levels = tax_group_order),
    chem_gr_label_plot = factor(
      chem_gr,
      levels = chem_levels,
      labels = chem_labels),
    variance_source = case_when(
      group == "CAS" ~ "Chemical identity",
      group == "genus_phylo" ~ "Phylogeny",
      TRUE ~ as.character(group)),
    variance_source = factor(
      variance_source,
      levels = c("Chemical identity", "Phylogeny")),
    prop_estimate = variance_proportion,
    prop_lower = pmax(variance_proportion_lower, 0),
    prop_upper = pmin(variance_proportion_upper, 1)
  )

p_var <- ggplot(variance_plot_data, aes(prop_estimate, chem_gr_label_plot,
                                        colour = taxon_label_plot, shape = variance_source)) +
  geom_errorbar(aes(xmin = prop_lower, xmax = prop_upper), orientation = "y",
                height = 0.45, linewidth = 1, position = dodge_fig8) +
  geom_point(size = 5, position = dodge_fig8) +
  facet_grid(taxon_label_plot ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_colour_manual(values = tax_cols, guide = "none") +
  scale_shape_manual(values = c("Chemical identity" = 16, "Phylogeny" = 17),
                     name = "Variance source") +
  scale_x_continuous(breaks = seq(0, 1, 0.25), labels = percent_format(accuracy = 1),
                     limits = c(0, 1), expand = expansion(mult = c(0.02, 0.03))) +
  labs(x = "Proportion of random-effect variance", y = NULL,
       title = "(b) Random-effect variance partitioning") +
  theme_fig8 +
  theme(
    strip.placement = "outside",
    strip.text.y.left = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

fig8 <- p_r2 | plot_spacer() | p_var

fig8 <- fig8 +
  plot_layout(widths = c(1, 0.08, 1), guides = "collect") &
  theme(legend.position = "bottom", 
        legend.box = "horizontal", 
        legend.spacing.x = unit(5, "cm"),
        legend.box.spacing = unit(0.8, "cm"))

fig8

ggsave(here("figures", "Figure8.jpeg"),
       fig8, width = 14, height = 10, dpi = 600)

ggsave(here("figures", "Figure8.pdf"),
       fig8, width = 14, height = 10, dpi = 600)

#=======Figure S3: Posterior fitted relationships ==============================
base_dir <- here("data", "8_bayesian_model_fits")

fit_dirs <- c(
  algae        = file.path(base_dir, "algae_fits"),
  invertebrate = file.path(base_dir, "invertebrate_fits"),
  fish         = file.path(base_dir, "fish_fits")
)

taxon_label_vec <- c(
  algae        = "Algae",
  invertebrate = "Invertebrate",
  fish         = "Fish"
)

chem_gr_label_vec <- c(
  Industrial_Chemical   = "Industrial chemical",
  Pharmaceutical        = "Pharmaceutical",
  P_Arthropod_control   = "P: Arthropods",
  P_Plant_control       = "P: Plants",
  P_Microbial_control   = "P: Microbes",
  P_Mixed_target        = "P: Mixed"
)

model_files <- expand_grid(
  taxon = taxon_order,
  chem_gr = chem_gr_order,
  hypothesis = "H3"
) %>%
  mutate(
    file = file.path(
      fit_dirs[taxon],
      glue::glue("{hypothesis}_{taxon}_{chem_gr}.rds")
    ),
    taxon_label = unname(taxon_label_vec[taxon]),
    chem_gr_label = unname(chem_gr_label_vec[chem_gr])
  )

best_model_fits <- model_files %>%
  filter(file.exists(file)) %>%
  mutate(fit = map(file, readRDS))

ffg_cols <- c(
  setNames(c("#a1d99b", "#41ab5d", "#006d2c"),
           c("Diatom", "Green algae", "Other algae")),
  setNames(colorRampPalette(c("#feedde", "#fd8d3c", "#a63603"))(9),
           c("I_Filterer", "I_Filterer_A", "I_Gatherer", "I_Gatherer_A",
             "I_Herbivore", "I_Herbivore_A", "I_Predator", "I_Predator_A",
             "I_Shredder_A")),
  setNames(colorRampPalette(c("#efedf5", "#9e9ac8", "#54278f"))(5),
           c("F_Detritivore_Herbivore", "F_Invertivore", "F_Omnivore",
             "F_Piscivore", "F_Planktivore"))
)

ffg_cols <- setNames(
  ffg_cols,
  dplyr::recode(names(ffg_cols), !!!ffg_label_vec, .default = names(ffg_cols))
)

make_h3_pred_grid <- function(fit, n_mass = 100) {
  
  dat <- fit$data %>%
    filter(!is.na(log_mass_c), !is.na(FFG)) %>%
    mutate(FFG = as.character(FFG))
  
  ffg_levels <- if (is.factor(fit$data$FFG)) {
    levels(fit$data$FFG)
  } else {
    sort(unique(dat$FFG))
  }
  
  dat %>%
    group_by(FFG) %>%
    summarise(
      log_mass_c = list(seq(min(log_mass_c), max(log_mass_c),
                            length.out = n_mass)),
      .groups = "drop"
    ) %>%
    unnest(log_mass_c) %>%
    mutate(FFG = factor(FFG, levels = ffg_levels))
}

predict_h3 <- function(fit, taxon, chem_gr, taxon_label, chem_gr_label, file,
                       n_mass = 100) {
  
  newdata <- make_h3_pred_grid(fit, n_mass)
  
  fitted(
    fit,
    newdata = newdata,
    re_formula = NA,
    summary = TRUE,
    probs = c(0.025, 0.975)
  ) %>%
    as.data.frame() %>%
    bind_cols(newdata, .) %>%
    transmute(
      taxon, chem_gr, taxon_label, chem_gr_label, file,
      FFG,
      log_mass_c,
      pred_mean  = Estimate,
      pred_lower = Q2.5,
      pred_upper = Q97.5
    )
}

figS3_data <- best_model_fits %>%
  select(fit, taxon, chem_gr, taxon_label, chem_gr_label, file) %>%
  pmap_dfr(predict_h3) %>%
  mutate(
    taxon_label_plot = factor(
      taxon_label,
      levels = unname(taxon_label_vec[taxon_order])
    ),
    chem_gr_label_plot = factor(
      chem_gr_label,
      levels = unname(chem_gr_label_vec[chem_gr_order])
    ),
    FFG_plot = dplyr::recode(
      as.character(FFG),
      !!!ffg_label_vec,
      .default = as.character(FFG)
    ),
    FFG_plot = factor(FFG_plot, levels = names(ffg_cols))
  )

figS3_data %>%
  distinct(chem_gr, chem_gr_label_plot) %>%
  arrange(chem_gr)

fs <- 18

figS3 <- ggplot(
  figS3_data,
  aes(log_mass_c, pred_mean,
      colour = FFG_plot, fill = FFG_plot, group = FFG_plot)
) +
  geom_hline(yintercept = 0, linetype = "dashed",
             linewidth = 0.3, colour = "grey80") +
  geom_ribbon(aes(ymin = pred_lower, ymax = pred_upper),
              alpha = 0.08, linewidth = 0, colour = NA) +
  geom_line(linewidth = 1) +
  facet_grid(
    taxon_label_plot ~ chem_gr_label_plot,
    drop = FALSE,
    labeller = labeller(
      taxon_label_plot = label_wrap_gen(12),
      chem_gr_label_plot = label_wrap_gen(16)
    )
  ) +
  scale_colour_manual(values = ffg_cols, name = NULL, drop = TRUE) +
  scale_fill_manual(values = ffg_cols, name = NULL, drop = TRUE) +
  guides(
    colour = guide_legend(ncol = 5),
    fill = guide_legend(ncol = 5)
  ) +
  labs(
    x = expression("Centred " * log[10] * " body mass"),
    y = expression("Centred predicted " * log[10] * " relative sensitivity")
  ) +
  plot_theme(fs) +
  theme(
    panel.grid.major = element_line(colour = "grey80", linewidth = 0.2),
    axis.title = element_text(size = fs, face = "bold", colour = "black"),
    axis.text = element_text(size = fs - 2, colour = "black"),
    strip.text = element_text(size = fs, face = "bold", lineheight = 0.9),
    legend.position = "bottom",
    legend.text = element_text(size = fs, colour = "black"),
    panel.spacing.x = unit(0.6, "cm")
  )

figS3

ggsave(here("figures", "FigureS3.jpeg"),
       figS3, width = 15, height = 9, dpi = 600)


#=======Figure S4-S6: Posterior distributions of H3 slopes and intercepts=======
base_dir <- here("data", "8_bayesian_model_fits")

fit_dirs <- c(
  algae        = file.path(base_dir, "algae_fits"),
  invertebrate = file.path(base_dir, "invertebrate_fits"),
  fish         = file.path(base_dir, "fish_fits")
)

taxon_label_vec <- c(
  algae        = "Algae",
  invertebrate = "Invertebrate",
  fish         = "Fish"
)

chem_gr_label_vec <- c(
  Industrial_Chemical   = "Industrial chemical",
  Pharmaceutical        = "Pharmaceutical",
  P_Arthropod_control   = "P: Arthropods",
  P_Plant_control       = "P: Plants",
  P_Microbial_control   = "P: Microbes",
  P_Mixed_target        = "P: Mixed"
)

h3_model_files <- expand_grid(
  taxon = taxon_order,
  chem_gr = chem_gr_order,
  hypothesis = "H3"
) %>%
  mutate(
    file = file.path(
      fit_dirs[taxon],
      glue::glue("{hypothesis}_{taxon}_{chem_gr}.rds")
    ),
    taxon_label = unname(taxon_label_vec[taxon]),
    chem_gr_label = unname(chem_gr_label_vec[chem_gr]),
    file_exists = file.exists(file)
  )

h3_model_files %>%
  count(taxon, chem_gr, file_exists)

extract_h3_slopes_intercepts <- function(file, taxon, chem_gr, taxon_label, chem_gr_label,
                                         hypothesis, file_exists) {
  
  fit <- readRDS(file)
  dat <- fit$data

  ffg_levels <- if (is.factor(dat$FFG)) levels(dat$FFG) else sort(unique(dat$FFG))
  
  newdata_0 <- tibble(log_mass_c = 0, FFG = ffg_levels)
  newdata_1 <- tibble(log_mass_c = 1, FFG = ffg_levels)
  
  if (is.factor(dat$FFG)) {
    newdata_0 <- mutate(newdata_0, FFG = factor(FFG, levels = levels(dat$FFG)))
    newdata_1 <- mutate(newdata_1, FFG = factor(FFG, levels = levels(dat$FFG)))
  }
  
  pred_0 <- as.matrix(posterior_epred(fit, newdata = newdata_0,
                                      re_formula = NA, summary = FALSE))
  pred_1 <- as.matrix(posterior_epred(fit, newdata = newdata_1,
                                      re_formula = NA, summary = FALSE))
  
  colnames(pred_0) <- colnames(pred_1) <- paste0("V", seq_len(ncol(pred_0)))
  
  bind_rows(
    as_tibble(pred_1 - pred_0) %>% mutate(.draw = row_number(), parameter_type = "Slope"),
    as_tibble(pred_0)          %>% mutate(.draw = row_number(), parameter_type = "Intercept")
  ) %>%
    pivot_longer(starts_with("V"), names_to = ".row", values_to = "value") %>%
    mutate(.row = as.integer(str_remove(.row, "V"))) %>%
    left_join(
      newdata_0 %>% mutate(.row = row_number()) %>% select(.row, FFG),
      by = ".row"
    ) %>%
    transmute(
      taxon, chem_gr, taxon_label, chem_gr_label, hypothesis, file,
      .draw, FFG, parameter_type, value
    )
}

posterior_slopes_intercepts <- h3_model_files %>%
  pmap_dfr(extract_h3_slopes_intercepts) %>%
  mutate(
    taxon_label_plot = factor(
      taxon_label,
      levels = unname(taxon_label_vec[taxon_order])
    ),
    chem_gr_label_plot = factor(
      chem_gr_label,
      levels = unname(chem_gr_label_vec[chem_gr_order])
    ),
    FFG_plot = dplyr::recode(
      as.character(FFG),
      !!!ffg_label_vec,
      .default = as.character(FFG)
    ),
    parameter_type = factor(parameter_type, levels = c("Slope", "Intercept"))
  )

posterior_slopes_intercepts_summary <- posterior_slopes_intercepts %>%
  group_by(taxon, chem_gr, taxon_label_plot, chem_gr_label_plot,
           FFG, FFG_plot, parameter_type) %>%
  summarise(
    estimate = median(value, na.rm = TRUE),
    lower = quantile(value, 0.025, na.rm = TRUE),
    upper = quantile(value, 0.975, na.rm = TRUE),
    prob_positive = mean(value > 0, na.rm = TRUE),
    .groups = "drop"
  )

get_ffg_levels <- function(taxon_i, dat) {
  order_i <- switch(
    taxon_i,
    algae = ffg_order_algae,
    invertebrate = ffg_order_invertebrate,
    fish = ffg_order_fish,
    rev(unique(as.character(dat$FFG_plot)))
  )
  rev(order_i[order_i %in% unique(dat$FFG_plot)])
}

plot_one_parameter_panel <- function(dat, parameter_name,
                                     show_y_axis = TRUE,
                                     show_xlab = FALSE,
                                     fs = 12) {
  
  dat_sub <- dat %>%
    filter(parameter_type == parameter_name)
  
  ggplot(dat_sub, aes(value, FFG_plot)) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      linewidth = 0.4,
      colour = "grey60"
    ) +
    ggdist::stat_slab(
      fill = "grey80",
      colour = NA,
      alpha = 0.9,
      scale = 0.85
    ) +
    ggdist::stat_pointinterval(
      point_interval = ggdist::median_qi,
      .width = 0.95,
      colour = "black",
      point_size = 2.4,
      linewidth = 0.8
    ) +
    scale_x_continuous(
      breaks = scales::breaks_pretty(n = 3),
      labels = scales::label_number(accuracy = 0.1)
    ) +
    labs(
      x = if (show_xlab) parameter_name else NULL,
      y = NULL
    ) +
    theme_bw(base_size = fs) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey88", linewidth = 0.25),
      axis.title.x = element_text(size = fs, face = "bold", colour = "black"),
      axis.text.x = element_text(size = fs - 2, colour = "black"),
      axis.text.y = element_text(size = fs, colour = "black"),
      legend.position = "none",
      plot.margin = margin(3, 3, 3, 3)
    ) +
    {
      if (!show_y_axis) {
        theme(
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()
        )
      }
    }
}

plot_one_chem_panel <- function(dat, chem_title, panel_letter,
                                show_xlab = FALSE, fs = 12) {
  
  title_plot <- ggplot() +
    annotate(
      "text",
      x = 0.5, y = 0.5,
      label = paste0("(", panel_letter, ") ", chem_title),
      fontface = "bold",
      size = fs / 3,
      hjust = 0.5,
      vjust = 0.5
    ) +
    theme_void() +
    theme(
      plot.margin = margin(t = 2, r = 0, b = 8, l = 0)
    )
  
  p_slope <- plot_one_parameter_panel(
    dat = dat,
    parameter_name = "Slope",
    show_y_axis = TRUE,
    show_xlab = show_xlab,
    fs = fs
  )
  
  p_intercept <- plot_one_parameter_panel(
    dat = dat,
    parameter_name = "Intercept",
    show_y_axis = FALSE,
    show_xlab = show_xlab,
    fs = fs
  )
  
  title_plot / (p_slope | p_intercept) +
    plot_layout(heights = c(0.22, 1))
}

plot_one_taxon_slopes_intercepts <- function(taxon_i, fs = 12) {
  
  taxon_dat <- posterior_slopes_intercepts %>%
    filter(taxon == taxon_i)
  
  map2(chem_gr_order, letters[seq_along(chem_gr_order)], \(chem_i, panel_letter) {
    
    dat_i <- taxon_dat %>%
      filter(chem_gr == chem_i) %>%
      mutate(FFG_plot = factor(
        FFG_plot,
        levels = get_ffg_levels(taxon_i, taxon_dat)
      ))
    
      plot_one_chem_panel(
        dat = dat_i,
        chem_title = chem_gr_label_vec[[chem_i]],
        panel_letter = panel_letter,
        show_xlab = chem_i %in% chem_gr_order[5:6],
        fs = fs
      )
  }) %>%
    wrap_plots(nrow = 3, ncol = 2)
}

fs <- 20

figS4 <- plot_one_taxon_slopes_intercepts("algae", fs = fs)
figS5 <- plot_one_taxon_slopes_intercepts("invertebrate", fs = fs+2)
figS6 <- plot_one_taxon_slopes_intercepts("fish", fs = fs+2)

figS4
figS5
figS6

ggsave(here("figures", "FigureS4.jpeg"),
       figS4, width = 13, height = 10, dpi = 600)

ggsave(here("figures", "FigureS5.jpeg"),
       figS5, width = 17, height = 20, dpi = 600)

ggsave(here("figures", "FigureS6.jpeg"),
       figS6, width = 15, height = 16, dpi = 600)

#=======Figure S7: Phylogenetic signal across all taxa ========================
fits_dir <- here("data", "7_phylogenetic_signal_by_tax_group")

figS7_data <- list.files(fits_dir, "^phytools_lambda_by_CAS_.*\\.rds$", full.names = TRUE) %>%
  map_dfr(readRDS) %>%
  mutate(
    tax_group = standardise_tax_group(tax_group),
    chem_group = factor(chem_group, levels = chem_gr_order)
  ) %>%
  filter(tax_group == "All taxa", !is.na(chem_group), !is.na(lambda))

figS7_letters <- lm(lambda ~ chem_group, data = figS7_data) %>%
  emmeans(~ chem_group) %>%
  emmeans:::cld.emmGrid(adjust = "BH", Letters = letters, sort = TRUE, reversed = TRUE) %>%
  as.data.frame() %>%
  mutate(
    .group = str_replace_all(.group, " ", ""),
    chem_group = factor(chem_group, levels = chem_gr_order),
    x_pos = 1.03
  )

fs <- 18

figS7 <- ggplot(figS7_data, aes(lambda, chem_group)) +
  geom_boxplot(width = 0.55, alpha = 0.20, 
               outlier.shape = NA, linewidth = 0.9) +
  ggbeeswarm::geom_quasirandom(size = 1.6, alpha = 0.60, 
                               width = 0.18) +
  geom_text(data = figS7_letters, aes(x_pos, chem_group, label = .group),
            inherit.aes = FALSE, hjust = 0, size = 7, colour = "black") +
  facet_wrap(~ tax_group, nrow = 1) +
  scale_x_continuous(breaks = c(0, 0.5, 1), limits = c(0, 1.15),
                     expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_discrete(limits = chem_gr_order, 
                   labels = chem_gr_label[chem_gr_order], drop = FALSE) +
  labs(x = "Phylogenetic signal\n(Pagel's \u03BB)", y = NULL) +
  coord_cartesian(clip = "off") +
  plot_theme(fs) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text = element_text(size = fs, face = "bold"),
    plot.margin = margin(5.5, 35, 5.5, 5.5)
  )

figS7

ggsave(here("figures", "FigureS7.jpeg"), 
       figS7, width = 7.5, height = 5.5, dpi = 600)



#=======Figure S8-S10: PhyloPars genus-level prediction =======================
fits_dir <- here("data", "9_phylopars_pred")
observed_sensi <- readRDS(here("data", "6_combined_meta_data", "combined_observed_sensitivity.rds"))

mask_fracs <- c(0.30, 0.50, 0.70)
fig_names  <- c("FigureS8", "FigureS9", "FigureS10")

chem_title <- setNames(
  paste0("(", letters[seq_along(chem_gr_order)], ") ", chem_gr_label[chem_gr_order]),
  chem_gr_order
)

get_genus <- function(x) as.character(x) %>% str_to_lower() %>% str_squish() %>% word(1)
safe_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) return(NA_real_)
  unname(cor(x[ok], y[ok], method = "pearson"))
}
fmt <- function(x) ifelse(is.na(x), "NA", formatC(x, format = "f", digits = 2))

obs_long <- observed_sensi %>%
  filter(is.finite(conc)) %>%
  transmute(species = get_genus(latin), CAS = as.character(CAS_norm), observed = conc) %>%
  filter(!is.na(species), species != "", !is.na(CAS)) %>%
  group_by(species, CAS) %>%
  summarise(observed = mean(observed, na.rm = TRUE), .groups = "drop")

read_results <- function(mask_frac) {
  tag <- paste0("mask", round(mask_frac * 100))
  readRDS(file.path(fits_dir, paste0("phylopars_fits_all_", tag, ".rds")))
}

make_eval_df <- function(res) {
  phylo_long <- as.data.frame(res$fit$anc_recon[res$fit$tree$tip.label, , drop = FALSE]) %>%
    rownames_to_column("species") %>%
    pivot_longer(-species, names_to = "CAS", values_to = "phylo_pred") %>%
    mutate(species = get_genus(species), CAS = as.character(CAS))
  
  ml_long <- res$wide_full %>%
    mutate(species = get_genus(species)) %>%
    pivot_longer(-species, names_to = "CAS", values_to = "ML_pred") %>%
    mutate(CAS = as.character(CAS))
  
  masked_key <- res$masked_pairs %>%
    transmute(species = get_genus(species), CAS = as.character(CAS)) %>%
    distinct()
  
  phylo_long %>%
    semi_join(masked_key, by = c("species", "CAS")) %>%
    left_join(ml_long, by = c("species", "CAS")) %>%
    left_join(obs_long, by = c("species", "CAS")) %>%
    filter(is.finite(phylo_pred), is.finite(ML_pred))
}

make_phylopars_panel <- function(res, chem_group, fs = 18) {
  df <- make_eval_df(res)
  
  r_text <- paste0(
    "Phylo vs ML = ", fmt(safe_cor(df$phylo_pred, df$ML_pred)), "\n",
    "Phylo vs obs = ", fmt(safe_cor(df$phylo_pred, df$observed)), "\n",
    "ML vs obs = ", fmt(safe_cor(df$ML_pred, df$observed))
  )
  
  ggplot(df, aes(ML_pred, phylo_pred)) +
    geom_point(alpha = 0.4, size = 1.5, colour = "grey30") +
    geom_point(data = filter(df, is.finite(observed)), aes(observed, phylo_pred),
               inherit.aes = FALSE, size = 1.8, alpha = 0.7, colour = "blue", shape = 1) +
    geom_abline(slope = 1, intercept = 0, linewidth = 1.1, linetype = "dashed", colour = "orange") +
    annotate("text", x = -Inf, y = Inf, label = r_text,
             hjust = -0.05, vjust = 1.2, size = fs / 4, fontface = "bold") +
    labs(title = chem_title[[chem_group]],
         x = expression("ML- or observed sensitivity"),
         y = expression("Phylogeny-predicted sensitivity")) +
    theme_bw(base_size = fs) +
    theme(panel.grid = element_blank(),
          legend.position = "none",
          plot.title = element_text(face = "bold", hjust = 0.5))
}

for (i in seq_along(mask_fracs)) {
  res_list <- read_results(mask_fracs[i])
  
  fig_i <- chem_gr_order[chem_gr_order %in% names(res_list)] %>%
    map(~ make_phylopars_panel(res_list[[.x]], .x)) %>%
    wrap_plots(nrow = 2, ncol = 3)
  
  ggsave(here("figures", paste0(fig_names[i], ".jpeg")),
         fig_i, width = 15, height = 11, dpi = 600)
}


#=======Figure S11datascores of phylogenetic predictions =======================
fits_dir <- here("data", "9_phylopars_pred")
observed_sensi <- readRDS(here("data", "6_combined_meta_data", "combined_observed_sensitivity.rds"))

mask_fracs <- c(0.30, 0.50, 0.70)

get_genus <- function(x) as.character(x) %>% str_to_lower() %>% str_squish() %>% word(1)

obs_long <- observed_sensi %>%
  filter(is.finite(conc)) %>%
  transmute(species = get_genus(latin), CAS = as.character(CAS_norm), observed = conc) %>%
  filter(!is.na(species), species != "", !is.na(CAS)) %>%
  group_by(species, CAS) %>%
  summarise(observed = mean(observed, na.rm = TRUE), .groups = "drop")

read_results <- function(mask_frac) {
  tag <- paste0("mask", round(mask_frac * 100))
  readRDS(file.path(fits_dir, paste0("phylopars_fits_all_", tag, ".rds")))
}

make_z_df <- function(res, chem_group, mask_frac) {
  phylo_mu <- as.data.frame(res$fit$anc_recon[res$fit$tree$tip.label, , drop = FALSE]) %>%
    rownames_to_column("species") %>%
    pivot_longer(-species, names_to = "CAS", values_to = "phylo_pred") %>%
    mutate(species = get_genus(species), CAS = as.character(CAS))
  
  phylo_se <- as.data.frame(res$fit$anc_var[res$fit$tree$tip.label, , drop = FALSE]) %>%
    rownames_to_column("species") %>%
    pivot_longer(-species, names_to = "CAS", values_to = "phylo_var") %>%
    mutate(species = get_genus(species), CAS = as.character(CAS),
           phylo_se = sqrt(pmax(phylo_var, 0))) %>%
    select(species, CAS, phylo_se)
  
  ml_long <- res$wide_full %>%
    mutate(species = get_genus(species)) %>%
    pivot_longer(-species, names_to = "CAS", values_to = "ML") %>%
    mutate(CAS = as.character(CAS))
  
  masked_key <- res$masked_pairs %>%
    transmute(species = get_genus(species), CAS = as.character(CAS)) %>%
    distinct()
  
  base_df <- phylo_mu %>%
    semi_join(masked_key, by = c("species", "CAS")) %>%
    left_join(phylo_se, by = c("species", "CAS")) %>%
    left_join(ml_long, by = c("species", "CAS")) %>%
    left_join(obs_long, by = c("species", "CAS")) %>%
    filter(is.finite(phylo_pred), is.finite(phylo_se), phylo_se > 0)
  
  bind_rows(
    transmute(base_df, chem_group, mask_frac, truth_type = "ML",       z = (ML - phylo_pred) / phylo_se),
    transmute(base_df, chem_group, mask_frac, truth_type = "Observed", z = (observed - phylo_pred) / phylo_se)
  ) %>%
    filter(is.finite(z))
}

z_all <- map_dfr(mask_fracs, \(mf) {
  res_list <- read_results(mf)
  map_dfr(chem_gr_order[chem_gr_order %in% names(res_list)],
          ~ make_z_df(res_list[[.x]], chem_group = .x, mask_frac = mf))
}) %>%
  mutate(
    chem_group = factor(chem_group, levels = chem_gr_order),
    truth_type = factor(truth_type, levels = c("ML", "Observed")),
    mask_frac = factor(mask_frac, levels = mask_fracs,
                       labels = paste0("Mask ", c(30, 50, 70), "%"))
  ) %>%
  filter(abs(z) <= 6)

z_sum <- z_all %>%
  group_by(chem_group, truth_type, mask_frac) %>%
  summarise(
    n = n(),
    z_mean = mean(z),
    z_sd = sd(z),
    cover80 = mean(abs(z) <= 1.2816),
    cover95 = mean(abs(z) <= 1.96),
    .groups = "drop"
  )

fs <- 16

figS11 <- ggplot(z_all, aes(z, chem_group)) +
  ggridges::geom_density_ridges(scale = 1.05, rel_min_height = 0.01,
                                alpha = 0.85, linewidth = 0.25) +
  geom_vline(xintercept = 0, linewidth = 0.4) +
  facet_grid(truth_type ~ mask_frac, switch = "y") +
  scale_y_discrete(labels = chem_gr_label[chem_gr_order]) +
  scale_x_continuous(limits = c(-7, 7), breaks = seq(-6, 6, 2)) +
  labs(x = "Standardised prediction error", y = NULL) +
  theme_bw(base_size = fs) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white"),
    strip.placement = "outside",
    axis.text = element_text(size = fs, colour = "black"),
    axis.title = element_text(size = fs, colour = "black"),
    strip.text = element_text(size = fs + 1, colour = "black", face = "bold")
  )

figS11

ggsave(here("figures", "FigureS11.jpeg"), 
       figS11, width = 14, height = 8, dpi = 600)






