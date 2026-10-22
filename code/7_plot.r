# To plot figures in the main text

rm(list = ls())

library(here)
library(tidyverse)
library(ggtree)
library(patchwork)
library(scales)
library(emmeans)

#---- 1. Load data and general settings ---------------------------------------
full_df <- readRDS(here("data/6_combined_meta_data/full_meta_data.rds"))
full_tree <- readRDS(here("data/6_combined_meta_data/full_genus_tree.rds"))

inference_dir <- here("data/8_bayesian_model_fits/bayesian_model_inference")
phylo_dir <- here("data/7_phylogenetic_signal_by_tax_group")
figure_dir <- here("figures")

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

full_tree$tip.label <- str_to_lower(str_trim(full_tree$tip.label))

chem_gr_order <- c(
  "Industrial_Chemical",
  "Pharmaceutical",
  "P_Arthropod_control",
  "P_Plant_control",
  "P_Microbial_control",
  "P_Mixed_target")

chem_gr_label <- c(
  Industrial_Chemical   = "Industrial chemical",
  Pharmaceutical        = "Pharmaceutical",
  P_Arthropod_control   = "P: Arthropods",
  P_Plant_control       = "P: Plants",
  P_Microbial_control   = "P: Microbes",
  P_Mixed_target        = "P: Mixed")

taxon_order <- c("algae", "invertebrate", "fish")
tax_group_order <- c("Algae", "Invertebrate", "Fish")
taxon_label_vec <- setNames(tax_group_order, taxon_order)

tax_cols <- c(
  Algae        = "#1b9e77",
  Invertebrate = "#d95f02",
  Fish         = "#7570b3")

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
  "F_Planktivore"           = "Planktivore")

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
  "Non-arthropod predator")

ffg_order_fish <- c(
  "Detritivore/herbivore",
  "Invertivore",
  "Omnivore",
  "Piscivore",
  "Planktivore")

ffg_order_plot <- rev(unique(c(
  ffg_order_algae,
  ffg_order_invertebrate,
  ffg_order_fish)))

pick_mode <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_character_)
  names(sort(table(x), decreasing = TRUE))[1]
}

standardise_tax_group <- function(x) {
  x <- str_to_lower(str_trim(as.character(x)))
  dplyr::recode(x, algae = "Algae", invertebrate = "Invertebrate", fish = "Fish")
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
      strip.text = element_text(size = fs, face = "bold"))
}

#---- Figure 2: Full tree + six chemical-use groups ---------------------------
fs <- 18
fig2_dat <- full_df %>%
  mutate(
    tax_group = standardise_tax_group(tax_group),
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
summarise_genus_sensitivity <- function(d, cg) {
  genus_summary <- d %>%
    filter(chem_group == cg) %>%
    group_by(genus) %>%
    summarise(
      med  = median(log_rel_sensi, na.rm = TRUE),
      q25  = quantile(log_rel_sensi, 0.25, na.rm = TRUE, names = FALSE),
      q75  = quantile(log_rel_sensi, 0.75, na.rm = TRUE, names = FALSE),
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
  
  df <- summarise_genus_sensitivity(fig2_dat, cg)
  
  xb <- c(-1, 0, 1, 2)
  xlim_i <- range(xb)
  
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
      .groups = "drop") %>%
    mutate(
      y_range = y_max - y_min,
      
      # keep label clearly inside each taxonomic block
      y_lab = case_when(
        tax_group == "Algae" ~ y_min + 0.40 * y_range,
        TRUE ~ y_min + 0.50 * y_range),
      
      x_lab = xlim_i[1] + 0.12 * diff(xlim_i),
      lab = sprintf("%.2f [%.2f, %.2f]", med_g, q25_g, q75_g))
  
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
      legend.text = element_text(size = fs, color = "black"))
  
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
    plot.margin = margin(4, 12, 4, 4))

fig2

ggsave(filename = file.path(figure_dir, "Figure_2.jpeg"),
       plot = fig2, width = 16, height = 10, dpi = 600)

ggsave(filename = file.path(figure_dir, "Figure_2.pdf"),
       plot = fig2, width = 16, height = 10, dpi = 600)

# Summaries for revising Fig. 2 text
fig2_summary <- map_dfr(chem_gr_order, \(cg) {
  
  summarise_genus_sensitivity(fig2_dat, cg) %>%
    filter(!is.na(med), !is.na(tax_group)) %>%
    group_by(tax_group) %>%
    summarise(
      n_genus = n_distinct(genus),
      median = median(med, na.rm = TRUE),
      q25    = quantile(med, 0.25,  na.rm = TRUE, names = FALSE),
      q75    = quantile(med, 0.75,  na.rm = TRUE, names = FALSE),
      q025   = quantile(med, 0.025, na.rm = TRUE, names = FALSE),
      q975   = quantile(med, 0.975, na.rm = TRUE, names = FALSE),
      .groups = "drop") %>%
    mutate(
      chem_group = cg,
      chem_group_label = chem_gr_label[[cg]],
      .before = 1)
})

fig2_summary



#---- Figure 3: Phylogenetic signal within each taxonomic group and chemical-use group -----
plot_df <- list.files(phylo_dir,
                      pattern = "^phytools_lambda_by_CAS_.*\\.rds$",
                      full.names = TRUE) %>%
  map_dfr(readRDS) %>%
  mutate(
    tax_group = factor(standardise_tax_group(tax_group), levels = tax_group_order),
    chem_group = factor(chem_group, levels = chem_gr_order)) %>%
  filter(!is.na(tax_group), !is.na(chem_group), is.finite(lambda))

m_lm <- lm(lambda ~ tax_group * chem_group, data = plot_df)
anova_lm <- car::Anova(m_lm, type = 3)
print(anova_lm)

emm_chem_by_tax <- emmeans(m_lm, ~ chem_group | tax_group)
pairs(emm_chem_by_tax, adjust = "BH")

chem_letters_by_tax <- multcomp::cld(
  emm_chem_by_tax,
  adjust = "BH",
  Letters = letters,
  sort = TRUE,
  reversed = TRUE) %>%
  as.data.frame() %>%
  mutate(
    .group = str_replace_all(.group, " ", ""),
    chem_group = factor(chem_group, levels = chem_gr_order),
    tax_group = factor(
      tax_group,
      levels = c("Algae", "Invertebrate", "Fish")))

letter_pos <- chem_letters_by_tax %>%
  mutate(x_pos = 0.86)

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
    labels = chem_gr_label,
    drop = FALSE) +
  scale_colour_manual(
    values = tax_cols[tax_group_order],
    breaks = tax_group_order,
    drop = FALSE) +
  labs(x = "Phylogenetic signal\n(Pagel's \u03BB)", y = NULL) +
  plot_theme(fs) +
  theme(
    legend.position = "none",
    panel.spacing.x = unit(1.2, "cm"))

fig3

ggsave(filename = file.path(figure_dir, "Figure_3.jpeg"),
       plot = fig3, width = 14, height = 8, dpi = 600)

ggsave(filename = file.path(figure_dir, "Figure_3.pdf"),
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
    .groups = "drop") %>%
  arrange(tax_group)

signal_full_summary <- plot_df %>%
  group_by(tax_group, chem_group) %>%
  summarise(
    median_lambda = round(median(lambda, na.rm = TRUE), 4),
    q25_lambda = round(quantile(lambda, 0.25, na.rm = TRUE, names = FALSE), 4),
    q75_lambda = round(quantile(lambda, 0.75, na.rm = TRUE, names = FALSE), 4),
    .groups = "drop") %>%
  arrange(tax_group, chem_group)

print(signal_full_summary)

# Post-hoc comparisons: chemical groups within taxonomic groups
pairs(emm_chem_by_tax, adjust = "BH")



#---- Figure 4: Leave-one-out model comparison --------------------------------
fig4_data <- read_csv(
  file.path(inference_dir, "LOO_model_comparisons.csv"),
  show_col_types = FALSE) %>%
  mutate(
    taxon = factor(taxon, levels = taxon_order),
    chem_gr = factor(chem_gr, levels = chem_gr_order),
    hypothesis = factor(hypothesis, levels = c("H0", "H1", "H2", "H3")),
    model_id = factor(model_id, levels = c("M0", "M1", "M2", "M3")),
    taxon_label_plot = factor(
      taxon_label_vec[as.character(taxon)],
      levels = unname(taxon_label_vec[taxon_order])),
    chem_gr_label_plot = fct_reorder(chem_gr_label, as.integer(chem_gr))
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
  labs(x = expression(bold(Delta * ELPD~"relative to the best-supported model")),
       y = "Candidate model") +
  plot_theme(fs) +
  theme(
    axis.title = element_text(size = fs, face = "bold", colour = "black"),
    axis.text = element_text(size = fs - 2, colour = "black"),
    legend.text = element_text(size = fs, colour = "black"),
    strip.text = element_text(size = fs, face = "bold", lineheight = 0.9),
    legend.position = "bottom",
    panel.spacing.x = unit(0.6, "cm"))

fig4

ggsave(file.path(figure_dir, "Figure_4.jpeg"),
       fig4, width = 15, height = 9, dpi = 600)
ggsave(file.path(figure_dir, "Figure_4.pdf"),
       fig4, width = 15, height = 9, dpi = 600)

# summary in the text
loo_summary <- fig4_data %>%
  filter(selected) %>%
  select(taxon_label_plot, chem_gr_label_plot,
         model_id, hypothesis, elpd_loo, se_elpd_loo)

loo_summary



#---- Figure 5: Standardised FFG-specific slopes ------------------------------
chem_gr_label_2row <- setNames(str_replace(chem_gr_label, " ", "\n"), names(chem_gr_label))

rope_cols <- c("Inside ±0.10" = "#1b9e77", "Inside ±0.20" = "#66a61e",
               "Overlaps ±0.20" = "grey45", "Outside ±0.20" = "#d95f02")
fs <- 20

fig_slope_data <- read_csv(file.path(inference_dir, "standardised_FFG_slope_summary.csv"),
                           show_col_types = FALSE) %>%
  mutate(
    taxon_label_plot = factor(standardise_tax_group(taxon), levels = tax_group_order),
    chem_gr = factor(chem_gr, levels = chem_gr_order),
    chem_gr_label_plot = factor(chem_gr, levels = chem_gr_order,
                                labels = unname(chem_gr_label_2row[chem_gr_order])),
    hdi_rope_status = factor(case_when(
      hdi_inside_0.10 ~ "Inside ±0.10", hdi_inside_0.20 ~ "Inside ±0.20",
      hdi_outside_0.20 ~ "Outside ±0.20", TRUE ~ "Overlaps ±0.20"
    ), levels = names(rope_cols)),
    FFG_plot = dplyr::recode(FFG, !!!ffg_label_vec, .default = FFG),
    FFG_plot = factor(FFG_plot, levels = unique(c(ffg_order_plot, FFG_plot)))
  )

fig5 <- ggplot(fig_slope_data, aes(beta_std_median, FFG_plot)) +
  annotate("rect", xmin = -0.20, xmax = 0.20, ymin = -Inf, ymax = Inf,
           fill = "grey60", alpha = 0.22) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, colour = "grey35") +
  geom_errorbar(aes(xmin = hdi_low, xmax = hdi_high, colour = hdi_rope_status),
                orientation = "y", width = 0.28, linewidth = 1.4) +
  geom_point(aes(colour = hdi_rope_status), size = 4) +
  facet_grid(taxon_label_plot ~ chem_gr_label_plot, scales = "free", space = "free_y") +
  scale_colour_manual(values = rope_cols, name = "95% HDI vs ROPE") +
  scale_x_continuous(breaks = scales::breaks_pretty(n = 3),
                     labels = scales::label_number(accuracy = 0.1)) +
  labs(x = expression("Standardised body-mass slope (" * beta[std] * ")"), y = NULL) +
  plot_theme(fs) +
  theme(
    panel.grid.major.y = element_blank(), axis.title = element_text(face = "bold"),
    strip.text = element_text(size = fs - 2, face = "bold", colour = "black"),
    strip.background = element_rect(fill = "grey90", colour = "grey40"),
    legend.position = "bottom", legend.title = element_text(face = "bold"),
    panel.spacing.x = unit(0.35, "cm"), panel.spacing.y = unit(0.35, "cm")
  )
fig5

ggsave(file.path(figure_dir, "Figure_5.jpeg"), fig5, width = 16, height = 12, dpi = 600)
ggsave(file.path(figure_dir, "Figure_5.pdf"), fig5, width = 16, height = 12, dpi = 600)

# Summaries used in the text
write_csv(fig_slope_data %>%
            select(taxon, chem_gr, FFG, beta_std_median, hdi_low, hdi_high,
                   hdi_inside_0.10, hdi_inside_0.20, hdi_outside_0.20),
          file.path(inference_dir, "Figure5_slope_summary.csv"))

#---- Figure 6: Pairwise FFG-specific slope comparison heatmap ----------------
heatmap_theme <- plot_theme(fs) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 40, hjust = 1, size = fs - 1),
    axis.text.y = element_text(size = fs - 2), axis.title = element_text(face = "bold"),
    strip.placement = "outside", strip.text.y.left = element_text(face = "bold", angle = 90),
    strip.background = element_rect(fill = "grey85", colour = "grey30"),
    plot.caption = element_text(hjust = 1, size = fs, colour = "black"),
    legend.position = "right", legend.title = element_text(face = "bold"),
    legend.text = element_text(size = fs - 2)
  )

make_contrast_heatmap <- function(data, value, legend_title, caption) {
  plot_data <- data %>%
    mutate(
      taxon_label_plot = factor(standardise_tax_group(taxon), levels = tax_group_order),
      chem_gr = factor(chem_gr, levels = chem_gr_order),
      FFG_1_plot = dplyr::recode(str_trim(as.character(FFG_1)), !!!ffg_label_vec, .default = FFG_1),
      FFG_2_plot = dplyr::recode(str_trim(as.character(FFG_2)), !!!ffg_label_vec, .default = FFG_2),
      contrast_label = paste(FFG_1_plot, "\u2212", FFG_2_plot),
      hdi_rope_star = coalesce(as.character(effect_marker), "")
    )
  
  ggplot(plot_data, aes(chem_gr, contrast_label, fill = .data[[value]])) +
    geom_tile(colour = "white", linewidth = 0.25) +
    geom_text(aes(label = hdi_rope_star), hjust = 0.5, vjust = 0.7, 
              size = 9, colour = "black") +
    facet_grid(taxon_label_plot ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_x_discrete(labels = chem_gr_label_2row) +
    scale_fill_gradient2(midpoint = 0, name = legend_title) +
    labs(x = NULL, y = NULL, caption = caption) + heatmap_theme
}

fig_slope_heatmap_data <- read_csv(
  file.path(inference_dir, "pairwise_FFG_slope_contrast_summary.csv"), show_col_types = FALSE)
fig6 <- make_contrast_heatmap(
  fig_slope_heatmap_data, "slope_diff_median", expression(atop("Median", Delta * beta[std])),
  "* 95% HDI entirely outside ±0.20;  ** 95% HDI entirely outside ±0.50")
fig6

ggsave(file.path(figure_dir, "Figure_6.jpeg"), fig6, width = 14, height = 15, dpi = 600)
ggsave(file.path(figure_dir, "Figure_6.pdf"), fig6, width = 14, height = 15, dpi = 600)


#---- Figure 7: Pairwise FFG-specific intercept comparison heatmap ------------
fig_intercept_heatmap_data <- read_csv(
  file.path(inference_dir, "pairwise_FFG_intercept_contrast_summary.csv"), show_col_types = FALSE)
fig7 <- make_contrast_heatmap(
  fig_intercept_heatmap_data, "d_baseline_median", expression(atop("Median", italic(d))),
  paste("* 95% HDI entirely outside ±0.20;  ** 95% HDI entirely outside ±0.50;",
        " *** 95% HDI entirely outside ±0.80"))
fig7

ggsave(file.path(figure_dir, "Figure_7.jpeg"), fig7, width = 15, height = 15, dpi = 600)
ggsave(file.path(figure_dir, "Figure_7.pdf"), fig7, width = 15, height = 15, dpi = 600)


#---- Figure 8: Model explanatory power and random-effect variance partitioning ----
r2_marginal_conditional_best <- read_csv(file.path(inference_dir, "r2_marginal_conditional_best_model.csv"), show_col_types = FALSE)

random_variance_partition_best <- read_csv(file.path(inference_dir, "random_variance_partition_best_model.csv"), show_col_types = FALSE)

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
    plot.title = element_text(size = fs, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = fs))

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
    r2_type = dplyr::recode(r2_type,
                            "Marginal R2" = "Marginal R²", "Conditional R2" = "Conditional R²"),
    r2_type = factor(r2_type, levels = c("Marginal R²", "Conditional R²")))

p_r2 <- ggplot(r2_plot_data, aes(Estimate, chem_gr_label_plot,
                                 colour = taxon_label_plot, shape = r2_type)) +
  geom_errorbar(aes(xmin = Q2.5, xmax = Q97.5), orientation = "y",
                width = 0.45, linewidth = 1, position = dodge_fig8) +
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
  filter(group %in% c("CAS", "genus_phylo")) %>%
  mutate(
    taxon_label_plot = factor(
      standardise_tax_group(taxon),
      levels = tax_group_order),
    chem_gr_label_plot = factor(
      chem_gr,
      levels = chem_levels,
      labels = chem_labels),
    variance_source = dplyr::recode(group,
                                    CAS = "Chemical identity", genus_phylo = "Phylogeny"),
    variance_source = factor(
      variance_source,
      levels = c("Chemical identity", "Phylogeny")),
    prop_estimate = variance_proportion,
    prop_lower = pmax(variance_proportion_lower, 0),
    prop_upper = pmin(variance_proportion_upper, 1))

p_var <- ggplot(variance_plot_data, aes(prop_estimate, chem_gr_label_plot,
                                        colour = taxon_label_plot, shape = variance_source)) +
  geom_errorbar(aes(xmin = prop_lower, xmax = prop_upper), orientation = "y",
                width = 0.45, linewidth = 1, position = dodge_fig8) +
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
    axis.ticks.y = element_blank())

fig8 <- p_r2 | plot_spacer() | p_var

fig8 <- fig8 +
  plot_layout(widths = c(1, 0.08, 1), guides = "collect") &
  theme(legend.position = "bottom",
        legend.box = "horizontal",
        legend.spacing.x = unit(5, "cm"),
        legend.box.spacing = unit(0.8, "cm"))

fig8

ggsave(file.path(figure_dir, "Figure_8.jpeg"),
       fig8, width = 14, height = 10, dpi = 600)

ggsave(file.path(figure_dir, "Figure_8.pdf"),
       fig8, width = 14, height = 10, dpi = 600)



























