# Figures in the appendix

rm(list = ls())

library(here)
library(tidyverse)
library(brms)
library(patchwork)

#---- 1. Paths and general settings ------------------------------------------
model_dir <- here("data/8_bayesian_model_fits")
lambda_dir <- here("data/7_phylogenetic_signal_by_tax_group/lambda_random_exclusion")
exclusion_dir <- file.path(model_dir, "random_exclusion_test")
figure_dir <- here("figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

taxon_order <- c("algae", "invertebrate", "fish")
taxon_labels <- c(algae = "Algae", invertebrate = "Invertebrate", fish = "Fish")
tax_order <- unname(taxon_labels[taxon_order])
tax_cols <- c(Algae = "#1b9e77", Invertebrate = "#d95f02", Fish = "#7570b3")

chem_gr_order <- c(
  "Industrial_Chemical", "Pharmaceutical", "P_Arthropod_control",
  "P_Plant_control", "P_Microbial_control", "P_Mixed_target"
)

chem_gr_labels <- c(
  Industrial_Chemical = "Industrial chemical", Pharmaceutical = "Pharmaceutical",
  P_Arthropod_control = "P: Arthropods", P_Plant_control = "P: Plants",
  P_Microbial_control = "P: Microbes", P_Mixed_target = "P: Mixed"
)

models <- expand_grid(chem_gr = chem_gr_order, taxon = taxon_order) %>%
  mutate(file = file.path(model_dir, paste0(taxon, "_fits"),
                          paste0("H3_", taxon, "_", chem_gr, ".rds")))

fs <- 18
comparison_theme <- theme_bw(base_size = fs) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "top",
    axis.title = element_text(colour = "black"),
    axis.text = element_text(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text = element_text(face = "bold"),
    panel.spacing.x = grid::unit(1.2, "cm")
  )

#---- Figure S3: Posterior predictive checks ----------------------------------
set.seed(123)
ppc <- vector("list", nrow(models))

for (i in seq_len(nrow(models))) {
  message("PPC: ", models$taxon[i], " / ", models$chem_gr[i])
  
  fit <- readRDS(models$file[i])
  y <- fit$data$log_rel_sensi_z
  yrep <- posterior_predict(fit, ndraws = 30, re_formula = NULL)
  
  values <- c(list(Observed = y),
              set_names(asplit(yrep, 1), paste0("rep", seq_len(nrow(yrep)))))
  
  # Keep density curves to reduce memory use.
  ppc[[i]] <- map_dfr(values, function(x) {
    d <- density(x)
    tibble(x = d$x, density = d$y)
  }, .id = "draw") %>%
    mutate(taxon = models$taxon[i], chemical = models$chem_gr[i])
  
  rm(fit, y, yrep, values)
  gc()
}

ppc_df <- bind_rows(ppc) %>%
  mutate(
    taxon = factor(taxon, levels = taxon_order, labels = tax_order),
    chemical = factor(chemical, levels = chem_gr_order,
                      labels = unname(chem_gr_labels[chem_gr_order]))
  )

fig_S3 <- ggplot(ppc_df, aes(x, density, group = draw)) +
  geom_line(data = filter(ppc_df, draw != "Observed"),
            colour = "grey70", linewidth = 0.3, alpha = 0.6) +
  geom_line(data = filter(ppc_df, draw == "Observed"),
            colour = "black", linewidth = 0.6) +
  facet_grid(chemical ~ taxon, scales = "free", switch = "y") +
  labs(x = "Standardised log10 relative sensitivity", y = "Density") +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.x = element_text(face = "bold", size = 11),
    strip.text.y.left = element_text(angle = 90, size = 10),
    panel.spacing = grid::unit(0.6, "lines")
  )

fig_S3

ggsave(file.path(figure_dir, "Figure_S3.jpeg"), fig_S3, width = 9, height = 12, dpi = 600)

#---- Figure S4: Phylogenetic signal after random exclusion --------------------
lambda_results <- list.files(
  lambda_dir, pattern = "^l[au]mbda_random_exclusion_test_.*\\.rds$", full.names = TRUE) %>%
  map(readRDS) %>%
  map_dfr("results") %>%
  mutate(tax_group = factor(dplyr::recode(tax_group, !!!taxon_labels), levels = tax_order))

lambda_baseline <- lambda_results %>%
  filter(drop_rate == 0, status == "ok", is.finite(lambda)) %>%
  distinct(tax_group, chem_group, CAS, lambda_full = lambda)

lambda_comparison <- lambda_results %>%
  filter(drop_rate > 0, status == "ok", is.finite(lambda)) %>%
  inner_join(lambda_baseline, by = c("tax_group", "chem_group", "CAS")) %>%
  mutate(
    delta_lambda = lambda - lambda_full,
    exclusion = factor(drop_rate, levels = c(0.05, 0.10, 0.15),
                       labels = c("5% excluded", "10% excluded", "15% excluded"))
  )

lambda_comparison %>%
  group_by(tax_group, exclusion) %>%
  summarise(
    n_chemicals = n_distinct(CAS),
    median_change = median(delta_lambda),
    median_abs_change = median(abs(delta_lambda)),
    .groups = "drop"
  ) %>%
  print(n = Inf)

lambda_medians <- lambda_results %>%
  filter(status == "ok", is.finite(lambda)) %>%
  group_by(tax_group, chem_group, drop_rate, repeat_id) %>%
  summarise(median_lambda = median(lambda), .groups = "drop") %>%
  mutate(
    chem_group = factor(chem_group, levels = chem_gr_order),
    exclusion = factor(drop_rate, levels = c(0, 0.05, 0.10, 0.15),
                       labels = c("0", "5", "10", "15"))
  )

fig_S4 <- ggplot(lambda_medians, aes(exclusion, median_lambda, colour = tax_group)) +
  geom_boxplot(aes(group = interaction(exclusion, tax_group)),
               position = position_dodge(width = 0.7), width = 0.5,
               outlier.shape = NA, linewidth = 0.7, fill = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.7),
             size = 1.8, alpha = 0.6, show.legend = FALSE) +
  facet_wrap(~ chem_group, ncol = 3, labeller = labeller(chem_group = chem_gr_labels)) +
  scale_colour_manual(values = tax_cols, breaks = tax_order, drop = FALSE) +
  scale_y_continuous(breaks = scales::breaks_pretty(n = 6)) +
  labs(x = "Predicted sensitivity data excluded (%)",
       y = expression("Median phylogenetic signal" * "(Pagel's " * lambda * ")"),
       colour = NULL) +
  theme_bw(base_size = fs) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey92", linewidth = 0.4),
    margins = margin(6, 6, 6, 6),
    legend.position = "top",
    legend.key = element_blank(),
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black"),
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text = element_text(face = "bold"),
    panel.spacing = grid::unit(0.3, "cm")
  )

fig_S4

ggsave(file.path(figure_dir, "Figure_S4.jpg"), fig_S4, width = 12, height = 6, dpi = 600)

#---- Figure S5: FFG-specific allometric relationships -------------------------
ffg_label_vec <- c(
  Greenalgae = "Green algae", Diatom = "Diatom", Otheralgae = "Other algae",
  I_Filterer_A = "Arthropod filterer", I_Gatherer_A = "Arthropod gatherer",
  I_Herbivore_A = "Arthropod herbivore", I_Predator_A = "Arthropod predator",
  I_Shredder_A = "Arthropod shredder", I_Filterer = "Non-arthropod filterer",
  I_Gatherer = "Non-arthropod gatherer", I_Herbivore = "Non-arthropod herbivore",
  I_Predator = "Non-arthropod predator", F_Invertivore = "Invertivore",
  F_Piscivore = "Piscivore", F_Omnivore = "Omnivore",
  F_Detritivore_Herbivore = "Detritivore/Herbivore", F_Planktivore = "Planktivore"
)
ffg_order <- list(
  algae = c("Diatom", "Green algae", "Other algae"),
  invertebrate = c(
    "Arthropod filterer", "Arthropod gatherer", "Arthropod herbivore",
    "Arthropod predator", "Arthropod shredder", "Non-arthropod filterer",
    "Non-arthropod gatherer", "Non-arthropod herbivore", "Non-arthropod predator"
  ),
  fish = c("Detritivore/Herbivore", "Invertivore", "Omnivore", "Piscivore", "Planktivore")
)

posterior_slopes_intercepts <- vector("list", nrow(models))
mass_ranges <- vector("list", nrow(models))

for (i in seq_len(nrow(models))) {
  message("Allometry: ", models$taxon[i], " / ", models$chem_gr[i])
  
  fit <- readRDS(models$file[i])
  ffg <- unique(as.character(fit$data$FFG))
  
  mass_ranges[[i]] <- fit$data %>%
    group_by(FFG) %>%
    summarise(
      mass_min = min(log_mass_z),
      mass_max = max(log_mass_z),
      .groups = "drop"
    ) %>%
    mutate(taxon = models$taxon[i], chem_gr = models$chem_gr[i])
  
  # Fixed-effect intercepts at mean log10 body mass.
  newdata <- tibble(FFG = ffg, log_mass_z = 0)
  intercept <- posterior_linpred(fit, newdata = newdata, re_formula = NA)
  
  # Change in the fixed-effect mean per one-SD increase in log10 body mass.
  newdata$log_mass_z <- 1
  slope <- posterior_linpred(fit, newdata = newdata, re_formula = NA) - intercept
  
  colnames(intercept) <- ffg
  colnames(slope) <- ffg
  
  posterior_slopes_intercepts[[i]] <- map_dfr(
    list(Intercept = intercept, Slope = slope),
    function(x) {
      as_tibble(x) %>%
        mutate(.draw = row_number()) %>%
        pivot_longer(-.draw, names_to = "FFG", values_to = "value")
    },
    .id = "parameter_type"
  ) %>%
    mutate(taxon = models$taxon[i], chem_gr = models$chem_gr[i])
  
  rm(fit, intercept, slope)
  gc()
}

posterior_slopes_intercepts <- bind_rows(posterior_slopes_intercepts)
mass_ranges <- bind_rows(mass_ranges)

allometry_draws <- posterior_slopes_intercepts %>%
  select(taxon, chem_gr, FFG, .draw, parameter_type, value) %>%
  pivot_wider(names_from = parameter_type, values_from = value) %>%
  left_join(mass_ranges, by = c("taxon", "chem_gr", "FFG"))

allometry_summary <- allometry_draws %>%
  group_by(taxon, chem_gr, FFG) %>%
  group_modify(function(dat, key) {
    x <- seq(first(dat$mass_min), first(dat$mass_max), length.out = 60)
    mu <- outer(dat$Slope, x) + dat$Intercept
    tibble(
      log_mass_z = x, estimate = apply(mu, 2, median),
      lower = apply(mu, 2, quantile, 0.025), upper = apply(mu, 2, quantile, 0.975)
    )
  }) %>%
  ungroup() %>%
  mutate(
    taxon_label_plot = factor(taxon, levels = taxon_order, labels = tax_order),
    chem_gr_label_plot = factor(chem_gr, levels = chem_gr_order,
                                labels = unname(chem_gr_labels[chem_gr_order])),
    FFG_plot = dplyr::recode(FFG, !!!ffg_label_vec)
  )

ffg_cols <- c(
  setNames(c("#A1D99B", "#41AB5D", "#006D2C"), ffg_order$algae),
  setNames(grDevices::colorRampPalette(c("#FEE6CE", "#E6550D", "#A63603"))(9),
           ffg_order$invertebrate),
  setNames(c("#EFEDF5", "#BCBDDC", "#9E9AC8", "#756BB1", "#54278F"), ffg_order$fish)
)

chem_labels_S5 <- chem_gr_labels
chem_labels_S5["Industrial_Chemical"] <- "Industrial\nchemical"

fs <- 22

fig_S5 <- ggplot(
  allometry_summary,
  aes(log_mass_z, estimate, group = FFG, colour = FFG_plot, fill = FFG_plot)
) +
  geom_hline(yintercept = 0, linetype = "dotted",
             linewidth = 0.3, colour = "grey65") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.10, colour = NA) +
  geom_line(linewidth = 0.9) +
  facet_grid(
    taxon_label_plot ~ chem_gr,
    labeller = labeller(chem_gr = chem_labels_S5)
  ) +
  scale_colour_manual(values = ffg_cols, breaks = names(ffg_cols),
                      name = NULL, drop = FALSE) +
  scale_fill_manual(values = ffg_cols, breaks = names(ffg_cols),
                    name = NULL, drop = FALSE) +
  scale_x_continuous(breaks = scales::breaks_pretty(n = 4)) +
  scale_y_continuous(breaks = scales::breaks_pretty(n = 4)) +
  labs(
    x = "Standardised log10 body mass",
    y = "Predicted standardised log10 relative sensitivity"
  ) +
  theme_bw(base_size = fs) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey85", linewidth = 0.3),
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black"),
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text = element_text(face = "bold"),
    panel.spacing = grid::unit(0.4, "cm"),
    legend.position = "bottom",
    legend.text = element_text(size = fs - 4),
    legend.key.width = grid::unit(0.8, "cm")
  ) +
  guides(colour = guide_legend(ncol = 5), fill = guide_legend(ncol = 5))

fig_S5

ggsave(file.path(figure_dir, "Figure_S5.jpeg"), fig_S5,
       width = 18, height = 11, dpi = 600)

# #---- Figure S5: Alternative single-panel plot --------------------------------
# slope_hdi <- posterior_slopes_intercepts %>%
#   filter(parameter_type == "Slope") %>%
#   group_by(taxon, chem_gr, FFG) %>%
#   summarise(hdi = list(bayestestR::hdi(value, ci = 0.95)), .groups = "drop") %>%
#   mutate(includes_zero = map_lgl(hdi, ~ .x$CI_low[1] <= 0 & .x$CI_high[1] >= 0)) %>%
#   select(-hdi)
# 
# allometry_lines <- allometry_summary %>%
#   left_join(slope_hdi, by = c("taxon", "chem_gr", "FFG"))
# 
# fs <- 18
# 
# fig_S5_base <- ggplot(
#   allometry_lines,
#   aes(log_mass_z, estimate, linetype = chem_gr,
#       group = interaction(taxon, chem_gr, FFG))
# ) +
#   geom_line(data = filter(allometry_lines, !includes_zero),
#             aes(colour = taxon_label_plot), linewidth = 0.9) +
#   geom_line(data = filter(allometry_lines, includes_zero),
#             aes(colour = "HDI includes zero"), linewidth = 1) +
#   # geom_line(data = filter(allometry_lines, !includes_zero),
#   #           aes(colour = taxon_label_plot), linewidth = 0.6) +
#   scale_colour_manual(
#     values = c(tax_cols, "HDI includes zero" = "grey80"),
#     breaks = c(tax_order), name = "Taxonomic group"
#   ) +
#   scale_linetype_manual(
#     values = setNames(c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash"),
#                       chem_gr_order),
#     breaks = chem_gr_order, labels = unname(chem_gr_labels[chem_gr_order]),
#     name = ""
#   ) +
#   labs(x = "Standardised log10 body mass",
#        y = "Predicted standardised\nlog10 relative sensitivity") +
#   theme_bw(base_size = fs) +
#   theme(
#     panel.grid.minor = element_blank(),
#     axis.text = element_text(colour = "black"),
#     legend.text = element_text(size = fs-4),
#     legend.title = element_text(size = fs-4, face = "bold"),
#     legend.key.width = grid::unit(1.2, "cm"),
#     legend.key.height = grid::unit(0.6, "cm"),
#     legend.background = element_rect(fill = NA, colour = NA)
#   )
# 
# get_legend <- function(p) {
#   g <- ggplotGrob(p + theme(legend.position = "right"))
#   g$grobs[[which(g$layout$name == "guide-box-right")]]
# }
# 
# tax_legend <- get_legend(fig_S9_base + guides(
#   linetype = "none",
#   colour = guide_legend(override.aes = list(linetype = "solid", linewidth = 0.8))
# ))
# 
# chem_legend <- get_legend(fig_S9_base + guides(
#   colour = "none",
#   linetype = guide_legend(override.aes = list(colour = "black", linewidth = 0.8))
# ))
# 
# fig_S5_lines <- fig_S9_base + theme(legend.position = "none") +
#   #  inset_element(tax_legend, left = 0.01, bottom = 0.75, right = 0.34, top = 0.99) +
#   inset_element(chem_legend, left = 0.5, bottom = 0.11, right = 0.99, top = 0.41)
# 
# fig_S5_lines
# 
# ggsave(file.path(figure_dir, "Figure_S5_lines.jpeg"), fig_S5_lines,
#        width = 6, height = 4.5, dpi = 600)


#---- Figure S6: H3 slopes after random exclusion ------------------------------
exclusion_files <- tibble(
  file = list.files(exclusion_dir, pattern = "^H3_.*\\.rds$", full.names = TRUE)
) %>%
  mutate(filename = basename(file)) %>%
  extract(
    filename, into = c("tax", "chem_group", "drop_percent", "repeat_id"),
    regex = "^H3_(algae|invertebrate|fish)_(.+)_drop(5|10|15)_rep([0-9]+)\\.rds$",
    convert = TRUE
  )

baseline_files <- exclusion_files %>%
  distinct(tax, chem_group) %>%
  mutate(file = file.path(model_dir, paste0(tax, "_fits"),
                          paste0("H3_", tax, "_", chem_group, ".rds")))

extract_H3_slopes <- function(file) {
  fit <- readRDS(file)
  
  # Both variables were standardised using full-data means and SDs.
  slopes <- fixef(fit, robust = TRUE) %>%
    as.data.frame() %>%
    rownames_to_column("term") %>%
    filter(str_detect(term, ":log_mass_z$|^log_mass_z:")) %>%
    transmute(term, beta_std = Estimate)
  
  rm(fit)
  gc()
  slopes
}

baseline_slopes <- baseline_files %>%
  mutate(slopes = map(file, extract_H3_slopes)) %>%
  select(tax, chem_group, slopes) %>%
  unnest(slopes) %>%
  rename(beta_std_full = beta_std)

exclusion_slopes <- exclusion_files %>%
  mutate(slopes = map(file, extract_H3_slopes)) %>%
  select(tax, chem_group, drop_percent, repeat_id, slopes) %>%
  unnest(slopes)

slope_comparison <- exclusion_slopes %>%
  inner_join(baseline_slopes, by = c("tax", "chem_group", "term")) %>%
  mutate(
    delta_beta_std = beta_std - beta_std_full,
    same_direction = sign(beta_std) == sign(beta_std_full),
    tax_group = factor(recode(tax, !!!taxon_labels), levels = tax_order),
    exclusion = factor(drop_percent, levels = c(5, 10, 15),
                       labels = c("5% excluded", "10% excluded", "15% excluded"))
  )

slope_comparison %>%
  group_by(tax_group, exclusion) %>%
  summarise(
    n_estimates = n(),
    median_abs_change = median(abs(delta_beta_std)),
    q95_abs_change = quantile(abs(delta_beta_std), 0.95),
    prop_same_direction = mean(same_direction),
    .groups = "drop"
  ) %>%
  print(n = Inf)

slope_limits <- range(slope_comparison$beta_std_full, slope_comparison$beta_std, na.rm = TRUE)

fig_S6 <- ggplot(slope_comparison, aes(beta_std_full, beta_std, colour = tax_group)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed",
              linewidth = 0.8, colour = "black") +
  geom_point(size = 1.6, alpha = 0.6) +
  facet_wrap(~ exclusion, nrow = 1, drop = FALSE) +
  scale_colour_manual(values = tax_cols, breaks = tax_order, drop = FALSE) +
  coord_equal(xlim = slope_limits, ylim = slope_limits) +
  labs(x = expression(beta[std] ~ "(complete dataset)"),
       y = expression(beta[std] ~ "(reduced dataset)"), colour = NULL) +
  comparison_theme +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 3)))

fig_S6

ggsave(file.path(figure_dir, "Figure_S6.jpeg"), fig_S6, width = 12, height = 6, dpi = 600)

#---- Figures S7-S9: Posterior distributions of H3 slopes and intercepts -------
ffg_label_vec <- c(
  Greenalgae = "Green algae", Diatom = "Diatom", Otheralgae = "Other algae",
  I_Filterer_A = "Arthropod filterer", I_Gatherer_A = "Arthropod gatherer",
  I_Herbivore_A = "Arthropod herbivore", I_Predator_A = "Arthropod predator",
  I_Shredder_A = "Arthropod shredder", I_Filterer = "Non-arthropod filterer",
  I_Gatherer = "Non-arthropod gatherer", I_Herbivore = "Non-arthropod herbivore",
  I_Predator = "Non-arthropod predator", F_Invertivore = "Invertivore",
  F_Piscivore = "Piscivore", F_Omnivore = "Omnivore",
  F_Detritivore_Herbivore = "Detritivore/Herbivore", F_Planktivore = "Planktivore"
)

ffg_order <- list(
  algae = c("Diatom", "Green algae", "Other algae"),
  invertebrate = c(
    "Arthropod filterer", "Arthropod gatherer", "Arthropod herbivore",
    "Arthropod predator", "Arthropod shredder", "Non-arthropod filterer",
    "Non-arthropod gatherer", "Non-arthropod herbivore", "Non-arthropod predator"
  ),
  fish = c("Detritivore/Herbivore", "Invertivore", "Omnivore", "Piscivore", "Planktivore")
)

extract_h3_slopes_intercepts <- function(file) {
  fit <- readRDS(file)
  ffg_levels <- levels(fit$data$FFG)
  newdata_0 <- tibble(log_mass_z = 0, FFG = factor(ffg_levels, levels = ffg_levels))
  newdata_1 <- mutate(newdata_0, log_mass_z = 1)
  
  pred_0 <- posterior_epred(fit, newdata = newdata_0, re_formula = NA)
  pred_1 <- posterior_epred(fit, newdata = newdata_1, re_formula = NA)
  colnames(pred_0) <- colnames(pred_1) <- ffg_levels
  
  draws <- bind_rows(
    as_tibble(pred_1 - pred_0) %>% mutate(.draw = row_number(), parameter_type = "Slope"),
    as_tibble(pred_0) %>% mutate(.draw = row_number(), parameter_type = "Intercept")
  ) %>%
    pivot_longer(-c(.draw, parameter_type), names_to = "FFG", values_to = "value")
  
  mass_range <- fit$data %>%
    group_by(FFG) %>%
    summarise(mass_min = min(log_mass_z), mass_max = max(log_mass_z), .groups = "drop") %>%
    mutate(FFG = as.character(FFG))
  
  rm(fit)
  gc()
  list(draws = draws, mass_range = mass_range)
}

h3_results <- models %>%
  mutate(result = map(file, extract_h3_slopes_intercepts))

mass_ranges <- h3_results %>%
  transmute(taxon, chem_gr, ranges = map(result, "mass_range")) %>%
  unnest(ranges)

posterior_slopes_intercepts <- h3_results %>%
  transmute(taxon, chem_gr, draws = map(result, "draws")) %>%
  unnest(draws) %>%
  mutate(
    taxon_label_plot = factor(taxon, levels = taxon_order, labels = tax_order),
    chem_gr_label_plot = factor(chem_gr, levels = chem_gr_order,
                                labels = unname(chem_gr_labels[chem_gr_order])),
    FFG_plot = dplyr::recode(FFG, !!!ffg_label_vec),
    parameter_type = factor(parameter_type, levels = c("Slope", "Intercept"))
  )

rm(h3_results)
gc()

posterior_slopes_intercepts_summary <- posterior_slopes_intercepts %>%
  group_by(taxon, chem_gr, taxon_label_plot, chem_gr_label_plot, FFG, FFG_plot, parameter_type) %>%
  summarise(
    estimate = median(value), lower = quantile(value, 0.025),
    upper = quantile(value, 0.975), prob_positive = mean(value > 0), .groups = "drop"
  )

plot_parameter <- function(dat, parameter_name, show_y_axis, show_xlab, fs) {
  axis_labels <- c(
    Slope = "Standardised\nbody-mass slope",
    Intercept = "FFG baseline\n(response SD units)"
  )
  
  p <- ggplot(filter(dat, parameter_type == parameter_name), aes(value, FFG_plot)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, colour = "grey60") +
    ggdist::stat_slab(fill = "grey80", colour = NA, alpha = 0.9, scale = 0.85) +
    ggdist::stat_pointinterval(
      point_interval = ggdist::median_qi, .width = 0.95,
      colour = "black", point_size = 2.4, linewidth = 0.8
    ) +
    scale_x_continuous(breaks = scales::breaks_pretty(n = 3),
                       labels = scales::label_number(accuracy = 0.1)) +
    labs(x = if (show_xlab) axis_labels[[parameter_name]] else NULL, y = NULL) +
    theme_bw(base_size = fs) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey88", linewidth = 0.25),
      axis.title.x = element_text(face = "bold", colour = "black", size = fs - 2),
      axis.text.x = element_text(size = fs - 4, colour = "black"),
      axis.text.y = element_text(colour = "black", size = fs - 2),
      legend.position = "none", plot.margin = margin(6, 6, 6, 6)
    )
  
  if (!show_y_axis) {
    p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  }
  
  p
}
plot_taxon_posteriors <- function(taxon_i, fs = 20) {
  taxon_dat <- posterior_slopes_intercepts %>%
    filter(taxon == taxon_i) %>%
    mutate(FFG_plot = factor(FFG_plot, levels = rev(ffg_order[[taxon_i]])))
  
  map(seq_along(chem_gr_order), function(i) {
    dat <- filter(taxon_dat, chem_gr == chem_gr_order[i])
    
    title_plot <- ggplot() +
      annotate("text", x = 0.5, y = 0.5,
               label = paste0("(", letters[i], ") ", chem_gr_labels[[chem_gr_order[i]]]),
               fontface = "bold", size = fs / 3) +
      theme_void() + theme(plot.margin = margin(2, 0, 8, 0))
    
    p_slope <- plot_parameter(dat, "Slope", TRUE, i > 4, fs)
    p_intercept <- plot_parameter(dat, "Intercept", FALSE, i > 4, fs)
    
    title_plot / (p_slope | p_intercept) + plot_layout(heights = c(0.22, 1))
  }) %>%
    wrap_plots(nrow = 3, ncol = 2)
}

fig_S7 <- plot_taxon_posteriors("algae", fs = 18)
fig_S8 <- plot_taxon_posteriors("invertebrate", fs = 20)
fig_S9 <- plot_taxon_posteriors("fish", fs = 20)

fig_S7
fig_S8
fig_S9

ggsave(file.path(figure_dir, "Figure_S7.jpeg"), fig_S7, width = 13, height = 10, dpi = 600)
ggsave(file.path(figure_dir, "Figure_S8.jpeg"), fig_S8, width = 17, height = 20, dpi = 600)
ggsave(file.path(figure_dir, "Figure_S9.jpeg"), fig_S9, width = 15, height = 16, dpi = 600)

#---- Figure S10: Phylogenetic signal after removing lineage means --------------
lineage_dir <- here("data/7_phylogenetic_signal_by_tax_group/lambda_invertebrate_lineages")

chem_labels <- c(
  Industrial_Chemical = "Industrial chemicals",
  P_Arthropod_control = "Arthropod control",
  P_Microbial_control = "Microbial control",
  P_Mixed_target = "Mixed target",
  P_Plant_control = "Plant control",
  Pharmaceutical = "Pharmaceuticals"
)

lineage_results <- map_dfr(
  file.path(lineage_dir,
            paste0("lambda_invertebrate_lineages_", names(chem_labels), ".rds")),
  readRDS
) %>%
  mutate(chem_label = factor(chem_labels[chem_group], levels = unname(chem_labels)))

centred_comparison <- lineage_results %>%
  filter(analysis %in% c("raw", "centred")) %>%
  select(chem_label, CAS, analysis, lambda) %>%
  pivot_wider(names_from = analysis, values_from = lambda) %>%
  drop_na(raw, centred)

lineage_theme <- theme_bw(base_size = fs) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    axis.text = element_text(size = fs, colour = "black"),
    axis.title = element_text(size = fs, colour = "black"),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(size = fs, face = "bold"),
    panel.spacing.x = unit(0.6, "cm")
  )

lambda_lim <- max(1, max(centred_comparison$raw, centred_comparison$centred))

fig_S10 <- ggplot(centred_comparison, aes(raw, centred)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_point(size = 2, alpha = 0.65, colour = "#d95f02") +
  facet_wrap(~ chem_label, ncol = 3, drop = FALSE) +
  coord_equal(xlim = c(0, lambda_lim), ylim = c(0, lambda_lim)) +
  labs(x = expression(lambda ~ "(before removing lineage means)"),
       y = expression(lambda ~ "(after removing lineage means)")) +
  lineage_theme

fig_S10

ggsave(file.path(figure_dir, "Figure_S10.jpeg"), fig_S10,
       width = 11, height = 8, dpi = 600)

#---- Figure S11: Phylogenetic signal within invertebrate lineages --------------
lineage_cols <- c(
  Crustaceans = "#1b9e77", Insects = "#7570b3",
  Annelida = "#66a61e", Mollusca = "#e7298a"
)

within_lineage <- lineage_results %>%
  filter(analysis == "within_lineage", major_lineage %in% names(lineage_cols)) %>%
  mutate(major_lineage = factor(major_lineage, levels = names(lineage_cols)))

fig_S11 <- ggplot(within_lineage, aes(lambda, major_lineage, colour = major_lineage)) +
  geom_boxplot(orientation = "y", width = 0.55, outlier.shape = NA) +
  ggbeeswarm::geom_quasirandom(
    groupOnX = FALSE, size = 1.6, alpha = 0.6, width = 0.18
  ) +
  facet_wrap(~ chem_label, ncol = 3, drop = FALSE) +
  scale_x_continuous(breaks = scales::breaks_extended(n = 3)) +
  scale_y_discrete(drop = FALSE) +
  scale_colour_manual(values = lineage_cols, drop = FALSE) +
  labs(x = expression("Phylogenetic signal (Pagel's" ~ lambda * ")"), y = NULL) +
  lineage_theme

fig_S11

ggsave(file.path(figure_dir, "Figure_S11.jpeg"), fig_S11,
       width = 11, height = 8, dpi = 600)

#------Figure S12-S13: PhyloPars genus-level prediction =======================
fits_dir <- here("data", "9_phylopars_pred")
observed_sensi <- readRDS(here("data", "1_observed_sensitivity", "combined_observed_sensitivity", "combined_observed_toxicity.rds"))

mask_fracs <- c(0.30, 0.50)
fig_names  <- c("Figure_S12", "Figure_S13")

chem_title <- setNames(
  paste0("(", letters[seq_along(chem_gr_order)], ") ", chem_gr_labels[chem_gr_order]),
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
  transmute(species = get_genus(latin), CAS = as.character(CAS_norm), observed = conc_log10mgL) %>%
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




#---- Figure S14: Standardised errors of phylogenetic predictions --------------
fits_dir <- here("data/9_phylopars_pred")
mask_fracs <- c(0.30, 0.50)

get_genus <- function(x) str_to_lower(as.character(x)) %>% str_squish() %>% word(1)

matrix_long <- function(x, value) {
  as.data.frame(x) %>%
    rownames_to_column("species") %>%
    pivot_longer(-species, names_to = "CAS", values_to = value) %>%
    mutate(species = get_genus(species))
}

obs_long <- readRDS(here("data", "1_observed_sensitivity", "combined_observed_sensitivity", "combined_observed_toxicity.rds"))%>%
  transmute(species = get_genus(latin), CAS = as.character(CAS_norm), observed = conc_log10mgL) %>%
  group_by(species, CAS) %>%
  summarise(observed = mean(observed), .groups = "drop")

z_all <- map_dfr(mask_fracs, function(mf) {
  results <- readRDS(file.path(
    fits_dir, paste0("phylopars_fits_all_mask", round(mf * 100), ".rds")
  ))
  
  map_dfr(intersect(chem_gr_order, names(results)), function(cg) {
    res <- results[[cg]]
    tips <- res$fit$tree$tip.label
    
    masked <- res$masked_pairs %>%
      transmute(species = get_genus(species), CAS = as.character(CAS))
    
    ml <- res$wide_full %>%
      mutate(species = get_genus(species)) %>%
      pivot_longer(-species, names_to = "CAS", values_to = "ML")
    
    matrix_long(res$fit$anc_recon[tips, , drop = FALSE], "prediction") %>%
      semi_join(masked, by = c("species", "CAS")) %>%
      left_join(
        matrix_long(res$fit$anc_var[tips, , drop = FALSE], "variance"),
        by = c("species", "CAS")
      ) %>%
      left_join(ml, by = c("species", "CAS")) %>%
      left_join(obs_long, by = c("species", "CAS")) %>%
      filter(is.finite(prediction), is.finite(variance), variance > 0) %>%
      pivot_longer(c(ML, observed), names_to = "truth_type", values_to = "truth") %>%
      mutate(z = (truth - prediction) / sqrt(variance),
             chem_group = cg, mask_frac = mf) %>%
      filter(is.finite(z))
  })
}) %>%
  mutate(
    chem_group = factor(chem_group, levels = chem_gr_order),
    truth_type = factor(truth_type, levels = c("ML", "observed"),
                        labels = c("ML", "Observed")),
    mask_frac = factor(mask_frac, levels = mask_fracs,
                       labels = paste0("Mask ", mask_fracs * 100, "%"))
  )

fig_S14 <- ggplot(z_all, aes(z, chem_group)) +
  ggridges::geom_density_ridges(scale = 1.05, rel_min_height = 0.01,
                                alpha = 0.85, linewidth = 0.25) +
  geom_vline(xintercept = 0, linewidth = 0.4) +
  facet_grid(truth_type ~ mask_frac, switch = "y") +
  scale_y_discrete(labels = chem_gr_labels, drop = FALSE) +
  scale_x_continuous(breaks = seq(-6, 6, 2)) +
  coord_cartesian(xlim = c(-7, 7)) +
  labs(x = "Standardised prediction error", y = NULL) +
  theme_bw(base_size = fs) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white"),
    strip.placement = "outside",
    axis.text = element_text(size = fs, colour = "black"),
    axis.title = element_text(size = fs, colour = "black"),
    strip.text = element_text(size = fs, face = "bold", colour = "black")
  )

fig_S14

ggsave(file.path(figure_dir, "Figure_S14.jpeg"), fig_S14,
       width = 10, height = 8, dpi = 600)
