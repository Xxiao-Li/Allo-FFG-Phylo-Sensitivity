# To analyse outputs of all bayesian models in 5_H0_H3_all_HPC.R

rm(list = ls())

library(here)
library(tidyverse)
library(brms)
library(loo)
library(posterior)
library(bayestestR)
library(glue)

set.seed(2026)

#---- 1. Load model fits files and basic settings ------------------------------
base_dir <- here("data", "8_bayesian_model_fits")
output_dir <- file.path(base_dir, "bayesian_model_inference")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

fit_dirs <- c(algae = file.path(base_dir, "algae_fits"),
              invertebrate = file.path(base_dir, "invertebrate_fits"),
              fish = file.path(base_dir, "fish_fits"))

taxon_order <- c("algae", "invertebrate", "fish")
taxon_label <- c(algae = "Algae", invertebrate = "Invertebrate", fish = "Fish")

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
  P_Arthropod_control   = "P: Arthropod",
  P_Plant_control       = "P: Plant",
  P_Microbial_control   = "P: Microbial",
  P_Mixed_target        = "P: Mixed target")

hypothesis_order <- c("H0", "H1", "H2", "H3")

hypothesis_info <- tibble(
  hypothesis = c("H0", "H1", "H2", "H3"),
  model_id = c("M0", "M1", "M2", "M3"),
  hypothesis_label = c(
    "Baseline model",
    "Global allometry",
    "Group baselines",
    "Group allometry"),
  formula_structure = c(
    "Sensitivity ~ 0 + REchem + REphylo",
    "Sensitivity ~ 0 + Mass + REchem + REphylo",
    "Sensitivity ~ 0 + FFG + Mass + REchem + REphylo",
    "Sensitivity ~ 0 + FFG + FFG:Mass + REchem + REphylo"))

model_files <- expand_grid(
  taxon = taxon_order,
  chem_gr = chem_gr_order,
  hypothesis = hypothesis_order) %>%
  mutate(file = file.path(fit_dirs[taxon], 
                          glue("{hypothesis}_{taxon}_{chem_gr}.rds"))) %>%
  left_join(hypothesis_info, by = "hypothesis") %>%
  mutate(taxon_label = unname(taxon_label[taxon]), 
         chem_gr_label = unname(chem_gr_label[chem_gr]))
print(model_files)

# Read one model at a time
read_and_apply <- function(file, fun, ...) {
  fit <- readRDS(file)
  out <- fun(fit, ...)
  rm(fit)
  gc()
  out
}


#---- 2. Summary model settings and priors -------------------------------------
# priors, formula, chains, iterations, warmup, adapt_delta, treedepth
extract_model_settings <- function(fit) {
  settings <- attr(fit, "sampling_settings")
  tibble(
    formula = paste(deparse(formula(fit), width.cutoff = 500), collapse = " "),
    family = family(fit)$family,
    link = family(fit)$link,
    nobs = nobs(fit),
    chains = settings$chains,
    iter = settings$iter,
    warmup = settings$warmup,
    thin = settings$thin,
    metric = settings$metric,
    adapt_delta = settings$adapt_delta,
    max_treedepth = settings$max_treedepth
  )
}

model_settings <- model_files %>%
  mutate(settings = map(file, ~ read_and_apply(.x, extract_model_settings))) %>%
  select(-file) %>%
  unnest(settings)

glimpse(model_settings)

write_csv(model_settings, file.path(output_dir, "model_settings_summary.csv"))

model_priors <- model_files %>%
  mutate(priors = map(file, ~ read_and_apply(.x, function(fit) {
    as_tibble(prior_summary(fit))
  }))) %>%
  select(taxon, chem_gr, taxon_label, chem_gr_label, 
         hypothesis, model_id, hypothesis_label, priors ) %>%
  unnest(priors)

glimpse(model_priors)

write_csv(model_priors, file.path(output_dir, "model_priors_summary.csv"))



#---- 3. Summary model diagnostics for all models ------------------------------
extract_diagnostics_one_model <- function(fit) {
  diag_df <- posterior::summarise_draws(
    posterior::as_draws_array(fit), "rhat", "ess_bulk", "ess_tail"
  )
  settings <- attr(fit, "sampling_settings")
  sampler <- nuts_params(fit) %>%
    arrange(Chain, Iteration) %>%
    group_by(Chain) %>%
    summarise(
      n_divergent = sum(Value[Parameter == "divergent__"]),
      max_treedepth = max(Value[Parameter == "treedepth__"]),
      n_treedepth_hits = sum(Value[Parameter == "treedepth__"] >= settings$max_treedepth),
      bfmi = mean(diff(Value[Parameter == "energy__"])^2) /
        var(Value[Parameter == "energy__"]),
      .groups = "drop"
    )
  tibble(
    max_rhat = max(diag_df$rhat, na.rm = TRUE),
    min_bulk_ess = min(diag_df$ess_bulk, na.rm = TRUE),
    min_tail_ess = min(diag_df$ess_tail, na.rm = TRUE),
    n_undefined_rhat = sum(!is.finite(diag_df$rhat)),
    n_undefined_ess = sum(!is.finite(diag_df$ess_bulk) | !is.finite(diag_df$ess_tail)),
    n_divergent = sum(sampler$n_divergent),
    max_treedepth = max(sampler$max_treedepth),
    n_treedepth_hits = sum(sampler$n_treedepth_hits),
    min_bfmi = min(sampler$bfmi)
  )
}

model_diagnostics <- model_files %>%
  mutate(diagnostics = map(file, ~ read_and_apply(.x, extract_diagnostics_one_model))) %>%
  unnest(diagnostics) %>%
  select(taxon, chem_gr, taxon_label, chem_gr_label,
         hypothesis, model_id, hypothesis_label,
         max_rhat, min_bulk_ess, min_tail_ess,
         n_undefined_rhat, n_undefined_ess, n_divergent, max_treedepth,
         n_treedepth_hits, min_bfmi)

glimpse(model_diagnostics)

write_csv(model_diagnostics, file.path(output_dir, "model_diagnostics_summary.csv"))



#---- 4. LOO comparison -------------------------------------------------------
loo_dir <- file.path(output_dir, "loo")
dir.create(loo_dir, recursive = TRUE, showWarnings = FALSE)

run_loo_from_file <- function(file, ndraws_loo = 2000) {
  loo_file <- file.path(
    loo_dir,
    paste0(tools::file_path_sans_ext(basename(file)), "_loo_", ndraws_loo, ".rds")
  )
  
  # Reuse saved LOO if the model file has not changed.
  if (file.exists(loo_file) &&
      file.info(loo_file)$mtime >= file.info(file)$mtime) {
    message("Loading saved LOO: ", basename(file))
    return(readRDS(loo_file))
  }
  
  message("Calculating LOO: ", basename(file))
  fit <- readRDS(file)
  
  loo_obj <- loo(
    fit, ndraws = ndraws_loo, cores = 1,
    pointwise = TRUE, moment_match = FALSE, save_psis = FALSE
  )
  
  saveRDS(loo_obj, loo_file)
  
  rm(fit)
  gc()
  loo_obj
}

compare_loo_one_group <- function(df, ndraws_loo = 2000) {
  stopifnot(setequal(df$hypothesis, hypothesis_order))
  
  message("Running LOO for: ", df$taxon[1], " / ", df$chem_gr[1])
  
  loo_list <- map(df$file, run_loo_from_file, ndraws_loo = ndraws_loo)
  names(loo_list) <- df$hypothesis
  
  loo_cmp <- loo_compare(loo_list) %>%
    as.data.frame() %>%
    rownames_to_column("hypothesis") %>%
    select(hypothesis, elpd_diff, se_diff)
  
  df %>%
    mutate(
      elpd_loo = map_dbl(loo_list, ~ .x$estimates["elpd_loo", "Estimate"]),
      se_elpd_loo = map_dbl(loo_list, ~ .x$estimates["elpd_loo", "SE"]),
      p_loo = map_dbl(loo_list, ~ .x$estimates["p_loo", "Estimate"]),
      se_p_loo = map_dbl(loo_list, ~ .x$estimates["p_loo", "SE"]),
      max_pareto_k = map_dbl(loo_list,
                             ~ max(.x$diagnostics$pareto_k, na.rm = TRUE)),
      n_bad_pareto_k = map_int(loo_list,
                               ~ sum(.x$diagnostics$pareto_k > 0.7, na.rm = TRUE))
    ) %>%
    left_join(loo_cmp, by = "hypothesis") %>%
    arrange(desc(elpd_loo)) %>%
    mutate(
      rank = row_number(),
      selected = rank == 1,
      selected_hypothesis = hypothesis[1]
    )
}

loo_results <- model_files %>%
  group_by(taxon, chem_gr) %>%
  group_split() %>%
  map_dfr(compare_loo_one_group, ndraws_loo = 2000)

loo_results <- loo_results %>%
  select(taxon, chem_gr, taxon_label, chem_gr_label,
         hypothesis, model_id, hypothesis_label,
         elpd_loo, se_elpd_loo, elpd_diff, se_diff, p_loo, se_p_loo,
         max_pareto_k, n_bad_pareto_k,
         rank, selected, selected_hypothesis) %>%
  arrange(taxon, chem_gr, rank)

write_csv(loo_results, file.path(output_dir, "LOO_model_comparisons.csv"))



#---- 5. Extract posterior estimates for H3 ------------------------------------
best_model_files <- model_files %>% filter(hypothesis == "H3")

posterior_fixed_effects_best <- best_model_files %>%
  mutate(fixed_effects = map(file, ~ read_and_apply(.x, function(fit) {
    fixef(fit, probs = c(0.025, 0.975)) %>%
      as.data.frame() %>% rownames_to_column("parameter") %>% as_tibble()
  }))) %>%
  select(taxon, chem_gr, taxon_label, chem_gr_label,
         hypothesis, model_id, hypothesis_label, fixed_effects) %>%
  unnest(fixed_effects)

write_csv(posterior_fixed_effects_best,
          file.path(output_dir, "posterior_fixed_effect_estimates_best_model.csv"))

#---- 6. Standardised FFG-specific slopes -------------------------------------
summarise_posterior_effect <- function(x, thresholds = c(0.20, 0.50), ci = 0.95) {
  hdi <- bayestestR::hdi(x, ci = ci)
  low <- hdi$CI_low[1]
  high <- hdi$CI_high[1]
  width <- high - low
  out <- tibble(
    median = median(x, na.rm = TRUE),
    q025 = quantile(x, 0.025, na.rm = TRUE), q975 = quantile(x, 0.975, na.rm = TRUE),
    hdi_low = low, hdi_high = high, hdi_width = width,
    pd = max(mean(x > 0, na.rm = TRUE), mean(x < 0, na.rm = TRUE))
  )
  
  for (thr in thresholds) {
    thr_name <- sprintf("%.2f", thr)
    overlap <- max(0, min(high, thr) - max(low, -thr))
    out[[paste0("prob_abs_gt_", thr_name)]] <- mean(abs(x) > thr, na.rm = TRUE)
    out[[paste0("hdi_rope_prop_", thr_name)]] <-
      if (width == 0) as.numeric(abs(low) <= thr) else overlap / width
    out[[paste0("hdi_inside_", thr_name)]] <- low >= -thr & high <= thr
    out[[paste0("hdi_outside_", thr_name)]] <- high < -thr | low > thr
  }
  out
}

extract_standardised_slopes <- function(fit) {
  draws <- posterior::as_draws_df(fit, variable = "^b_", regex = TRUE) %>% as_tibble()
  slope_cols <- names(draws) %>% str_subset("^b_FFG.*:log_mass_z$")
  scaling <- attr(fit, "standardisation")
  
  draws %>%
    select(.draw, all_of(slope_cols)) %>%
    pivot_longer(cols = all_of(slope_cols), names_to = "parameter", values_to = "beta_std") %>%
    mutate(
      parameter = str_remove(parameter, "^b_"), FFG = str_remove(parameter, ":log_mass_z$"),
      FFG = str_remove(FFG, "^FFG"), beta_raw = beta_std * scaling$sd_sensi / scaling$sd_mass
    ) %>%
    select(FFG, parameter, .draw, beta_raw, beta_std)
}

standardised_slope_draws <- best_model_files %>%
  mutate(slopes = map(file, ~ read_and_apply(.x, extract_standardised_slopes))) %>%
  select(taxon, chem_gr, taxon_label, chem_gr_label, slopes) %>%
  unnest(slopes)

standardised_slope_summary <- standardised_slope_draws %>%
  group_by(taxon, chem_gr, taxon_label, chem_gr_label, FFG, parameter) %>%
  summarise(effect = list(summarise_posterior_effect(beta_std, thresholds = c(0.10, 0.20, 0.50))),
            .groups = "drop") %>%
  unnest(effect) %>%
  rename(beta_std_median = median, beta_std_q025 = q025, beta_std_q975 = q975) %>%
  mutate(effect_marker = case_when(hdi_outside_0.50 ~ "**", hdi_outside_0.20 ~ "*", TRUE ~ ""))

write_csv(standardised_slope_summary, file.path(output_dir, "standardised_FFG_slope_summary.csv"))

standardised_slope_general_summary <- standardised_slope_summary %>%
  summarise(
    n_slopes = n(), median_abs_beta_std = median(abs(beta_std_median), na.rm = TRUE),
    n_hdi_inside_0.20 = sum(hdi_inside_0.20, na.rm = TRUE),
    n_hdi_outside_0.20 = sum(hdi_outside_0.20, na.rm = TRUE),
    n_hdi_outside_0.50 = sum(hdi_outside_0.50, na.rm = TRUE),
    prop_hdi_inside_0.20 = n_hdi_inside_0.20 / n_slopes,
    prop_hdi_outside_0.20 = n_hdi_outside_0.20 / n_slopes,
    prop_hdi_outside_0.50 = n_hdi_outside_0.50 / n_slopes
  )

write_csv(standardised_slope_general_summary,
          file.path(output_dir, "general_standardised_FFG_slope_summary.csv"))

#---- 7. Pairwise FFG slope and baseline contrasts ----------------------------
extract_pairwise_slope_contrasts <- function(df) {
  slope_wide <- df %>%
    select(.draw, FFG, beta_std) %>% pivot_wider(names_from = FFG, values_from = beta_std)
  ffg_levels <- setdiff(names(slope_wide), ".draw")
  if (length(ffg_levels) < 2) return(tibble())
  
  combn(ffg_levels, 2, simplify = FALSE) %>% map_dfr(\(pair) {
    FFG_1 <- pair[1]
    FFG_2 <- pair[2]
    slope_wide %>% transmute(
      FFG_1 = FFG_1, FFG_2 = FFG_2, .draw = .draw,
      beta_std_1 = .data[[FFG_1]], beta_std_2 = .data[[FFG_2]],
      slope_diff = beta_std_1 - beta_std_2
    )
  })
}

pairwise_slope_draws <- standardised_slope_draws %>%
  group_by(taxon, chem_gr, taxon_label, chem_gr_label) %>%
  group_modify(~ extract_pairwise_slope_contrasts(.x)) %>% ungroup()

pairwise_slope_summary <- pairwise_slope_draws %>%
  group_by(taxon, chem_gr, taxon_label, chem_gr_label, FFG_1, FFG_2) %>%
  summarise(effect = list(summarise_posterior_effect(slope_diff, thresholds = c(0.20, 0.50))),
            .groups = "drop") %>%
  unnest(effect) %>%
  rename(slope_diff_median = median, slope_diff_q025 = q025, slope_diff_q975 = q975) %>%
  mutate(effect_marker = case_when(hdi_outside_0.50 ~ "**", hdi_outside_0.20 ~ "*", TRUE ~ ""))

write_csv(pairwise_slope_summary, file.path(output_dir, "pairwise_FFG_slope_contrast_summary.csv"))

extract_pairwise_intercept_contrasts <- function(fit) {
  draws <- posterior::as_draws_df(fit, variable = "^b_FFG|^sigma$", regex = TRUE) %>% as_tibble()
  intercept_cols <- names(draws) %>%
    str_subset("^b_FFG") %>% str_subset(":log_mass_z$", negate = TRUE)
  if (length(intercept_cols) < 2) return(tibble())
  ffg_levels <- str_remove(intercept_cols, "^b_FFG")
  sd_sens <- attr(fit, "standardisation")$sd_sensi
  
  combn(seq_along(intercept_cols), 2, simplify = FALSE) %>% map_dfr(\(pair) {
    alpha_1 <- draws[[intercept_cols[pair[1]]]]
    alpha_2 <- draws[[intercept_cols[pair[2]]]]
    tibble(
      FFG_1 = ffg_levels[pair[1]], FFG_2 = ffg_levels[pair[2]], .draw = draws$.draw,
      alpha_1 = alpha_1, alpha_2 = alpha_2, contrast_std = alpha_1 - alpha_2,
      contrast_raw = contrast_std * sd_sens, sigma = draws$sigma, d_baseline = contrast_std / sigma
    )
  })
}

pairwise_intercept_draws <- best_model_files %>%
  mutate(contrasts = map(file, ~ read_and_apply(.x, extract_pairwise_intercept_contrasts))) %>%
  select(taxon, chem_gr, taxon_label, chem_gr_label, contrasts) %>%
  unnest(contrasts)

pairwise_intercept_summary <- pairwise_intercept_draws %>%
  group_by(taxon, chem_gr, taxon_label, chem_gr_label, FFG_1, FFG_2) %>%
  summarise(
    contrast_raw_median = median(contrast_raw, na.rm = TRUE),
    contrast_raw_q025 = quantile(contrast_raw, 0.025, na.rm = TRUE),
    contrast_raw_q975 = quantile(contrast_raw, 0.975, na.rm = TRUE),
    effect = list(summarise_posterior_effect(d_baseline, thresholds = c(0.20, 0.50, 0.80))),
    .groups = "drop"
  ) %>%
  unnest(effect) %>%
  rename(d_baseline_median = median, d_baseline_q025 = q025, d_baseline_q975 = q975) %>%
  mutate(effect_marker = case_when(hdi_outside_0.80 ~ "***", hdi_outside_0.50 ~ "**",
                                   hdi_outside_0.20 ~ "*", TRUE ~ ""))

write_csv(pairwise_intercept_summary, file.path(output_dir, "pairwise_FFG_intercept_contrast_summary.csv"))



#---- 8. Marginal/conditional R2 and random variation partitioning -------------
n_draws_r2  <- 2000
chunk_size  <- 2000
response_var <- "log_rel_sensi_z"

calc_r2_chunked <- function(fit, response_var, ndraws = 1000, chunk_size = 1000, re_formula = NA) {
  
  dat <- fit$data
  y <- dat[[response_var]]
  n <- length(y)
  
  n_available <- posterior::ndraws(
    posterior::as_draws_df(fit, variable = "sigma")
  )
  draw_ids <- sample(seq_len(n_available), size = min(ndraws, n_available))
  d <- length(draw_ids)
  
  sum_mu <- sum_mu2 <- sum_res <- sum_res2 <- numeric(d)
  chunks <- split(seq_len(n), ceiling(seq_len(n) / chunk_size))
  
  for (idx in chunks) {
    mu_i <- posterior_epred(fit, newdata = dat[idx, , drop = FALSE],
                            re_formula = re_formula, draw_ids = draw_ids) %>%
      as.matrix()
    res_i <- sweep(mu_i, 2, y[idx], "-")
    sum_mu   <- sum_mu   + rowSums(mu_i)
    sum_mu2  <- sum_mu2  + rowSums(mu_i^2)
    sum_res  <- sum_res  + rowSums(res_i)
    sum_res2 <- sum_res2 + rowSums(res_i^2)
    rm(mu_i, res_i)
    gc()
  }
  
  var_mu  <- (sum_mu2  - sum_mu^2  / n) / (n - 1)
  var_res <- (sum_res2 - sum_res^2 / n) / (n - 1)
  r2_draws <- var_mu / (var_mu + var_res)
  
  tibble(
    Estimate  = mean(r2_draws, na.rm = TRUE),
    Est.Error = sd(r2_draws, na.rm = TRUE),
    Q2.5      = quantile(r2_draws, 0.025, na.rm = TRUE),
    Q97.5     = quantile(r2_draws, 0.975, na.rm = TRUE)
  )
}

calc_r2_from_fit <- function(fit, response_var, ndraws, chunk_size) {
  bind_rows(
    calc_r2_chunked(fit, response_var, ndraws, chunk_size, re_formula = NA)   %>%
      mutate(r2_type = "Marginal R2", .before = 1),
    calc_r2_chunked(fit, response_var, ndraws, chunk_size, re_formula = NULL) %>%
      mutate(r2_type = "Conditional R2", .before = 1)
  )
}

# Keep the same R2 calculation as the previous script
r2_marginal_conditional_best <- best_model_files %>%
  mutate(r2 = map(file, ~ read_and_apply(.x, calc_r2_from_fit,
                                         response_var = response_var, ndraws = n_draws_r2, chunk_size = chunk_size))) %>%
  select(taxon, chem_gr, taxon_label, chem_gr_label,
         hypothesis, model_id, hypothesis_label, formula_structure, r2) %>%
  unnest(r2)

glimpse(r2_marginal_conditional_best)
write_csv(r2_marginal_conditional_best,
          file.path(output_dir, paste0("r2_marginal_conditional_best_model.csv")))

extract_random_variance_partition <- function(fit) {
  
  draws <- posterior::as_draws_df(fit, variable = "^b_|^sigma$|^sd_", regex = TRUE) %>%
    as_tibble() %>%
    mutate(.draw_id = row_number())
  
  sd_cols <- names(draws)[str_detect(names(draws), "^sd_")]
  
  draws %>%
    select(.draw_id, all_of(sd_cols)) %>%
    pivot_longer(
      cols = all_of(sd_cols),
      names_to = "component",
      values_to = "sd_value"
    ) %>%
    mutate(
      group = component %>%
        str_remove("^sd_") %>%
        str_remove("__.*$"),
      variance = sd_value^2
    ) %>%
    group_by(.draw_id, group) %>%
    summarise(variance = sum(variance, na.rm = TRUE), .groups = "drop") %>%
    group_by(.draw_id) %>%
    mutate(
      total_random_variance = sum(variance, na.rm = TRUE),
      prop_draw = variance / total_random_variance
    ) %>%
    ungroup() %>%
    filter(is.finite(prop_draw)) %>%
    group_by(group) %>%
    summarise(
      variance_estimate = median(variance, na.rm = TRUE),
      variance_lower = quantile(variance, 0.025, na.rm = TRUE),
      variance_upper = quantile(variance, 0.975, na.rm = TRUE),
      variance_proportion = median(prop_draw, na.rm = TRUE),
      variance_proportion_lower = quantile(prop_draw, 0.025, na.rm = TRUE),
      variance_proportion_upper = quantile(prop_draw, 0.975, na.rm = TRUE),
      variance_proportion_sd = sd(prop_draw, na.rm = TRUE),
      n_prop_draws = n(),
      n_unique_prop = n_distinct(prop_draw),
      .groups = "drop"
    )
}

# Proportions are within random-effect variance only (CAS + phylogeny)
random_variance_partition_best <- best_model_files %>%
  mutate(variance_components = map(file,
                                   ~ read_and_apply(.x, extract_random_variance_partition))) %>%
  select(
    taxon, chem_gr, taxon_label, chem_gr_label,
    hypothesis, model_id, hypothesis_label, formula_structure,
    variance_components
  ) %>%
  unnest(variance_components)

glimpse(random_variance_partition_best)

write_csv(random_variance_partition_best, 
          file.path(output_dir, "random_variance_partition_best_model.csv"))

