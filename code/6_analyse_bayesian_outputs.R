
rm(list=ls())

library(here)
library(tidyverse)
library(brms)
library(loo)
library(tidybayes)
library(posterior)
library(bayestestR)
library(performance)
library(glue)
library(fs)
library(patchwork)
library(scales)
library(ggdist)


`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

set.seed(2026)

#=======1. Load model fits files and basic settings========================
base_dir <- here("data", "8_bayesian_model_fits")
output_dir <- here("data", "8_bayesian_model_fits", "bayesian_model_inference")

dir_create(output_dir)

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
  P_Arthropod_control   = "P: Arthropod control",
  P_Plant_control       = "P: Plant control",
  P_Microbial_control   = "P: Microbial control",
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
    "Sensitivity ~ REchem + REphylo",
    "Sensitivity ~ Mass + REchem + REphylo",
    "Sensitivity ~ Mass + FFG + REchem + REphylo",
    "Sensitivity ~ Mass * FFG + REchem + REphylo"))

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

model_fits <- model_files %>%
  mutate(fit = map(file, readRDS))

#=======2. Summary model settings and priors===============================
# priors, formula, chains, iterations, warmup, adapt_delta, treedepth
get_control_value <- function(fit, value_name, default_value = NA_real_) {
  value <- tryCatch(
    fit$stan_args$control[[value_name]],
    error = function(e) NULL
  )
  if (is.null(value)) default_value else value
}

model_settings <- model_fits %>%
  mutate(
    formula = map_chr(fit,
      ~ deparse(formula(.x), width.cutoff = 500) %>% paste(collapse = " ")),
    family = map_chr(fit, ~ family(.x)$family),
    link = map_chr(fit, ~ family(.x)$link),
    nobs = map_int(fit, nobs),
    chains = map_int(fit, ~ .x$fit@sim$chains),
    iter = map_int(fit, ~ .x$fit@sim$iter),
    warmup = map_int(fit, ~ .x$fit@sim$warmup),
    thin = map_int(fit, ~ .x$fit@sim$thin),
    adapt_delta = map_dbl(fit,
      ~ get_control_value(.x, "adapt_delta", default_value = NA_real_)),
    max_treedepth = map_int(fit,
      ~ as.integer(get_control_value(.x, "max_treedepth", default_value = NA_integer_)))
  ) %>%
  select(
    taxon, chem_gr, taxon_label, chem_gr_label,
    hypothesis, model_id, hypothesis_label,
    formula_structure,
    formula, family, link, nobs,
    chains, iter, warmup, thin,
    adapt_delta, max_treedepth)

glimpse(model_settings)

write_csv(model_settings, file.path(output_dir, "model_settings_summary.csv"))

model_priors <- model_fits %>%
  mutate(priors = map(fit, ~ prior_summary(.x) %>%
        as_tibble())) %>%
  select(taxon, chem_gr, taxon_label, chem_gr_label, 
         hypothesis, model_id, hypothesis_label, priors ) %>%
  unnest(priors)

glimpse(model_priors)

write_csv(model_priors, file.path(output_dir, "model_priors_summary.csv"))

#=======3. Summary model diagnostics for all models======
extract_diagnostics_one_model <- function(fit) {
  summ <- summary(fit)
  diag_tables <- list(fixed = summ$fixed, random = summ$random, 
                      spec_pars = summ$spec_pars)
  diag_df <- diag_tables %>%
    compact() %>%
    map(as.data.frame) %>%
    bind_rows(.id = "parameter_group")
  rhat_col <- intersect(c("Rhat", "Rhat."), names(diag_df))[1]
  bulk_col <- intersect(c("Bulk_ESS", "Bulk.ESS"), names(diag_df))[1]
  tail_col <- intersect(c("Tail_ESS", "Tail.ESS"), names(diag_df))[1]
  nuts <- nuts_params(fit)
  tibble(
    max_rhat = if (!is.na(rhat_col)) max(diag_df[[rhat_col]], na.rm = TRUE) else NA_real_,
    min_bulk_ess = if (!is.na(bulk_col)) min(diag_df[[bulk_col]], na.rm = TRUE) else NA_real_,
    min_tail_ess = if (!is.na(tail_col)) min(diag_df[[tail_col]], na.rm = TRUE) else NA_real_,
    n_divergent = sum(
      nuts$Parameter == "divergent__" & nuts$Value == 1, na.rm = TRUE),
    max_treedepth = max(
      nuts$Value[nuts$Parameter == "treedepth__"], na.rm = TRUE ))
}

model_diagnostics <- model_fits %>%
  mutate(diagnostics = map(fit, extract_diagnostics_one_model)) %>%
  unnest(diagnostics) %>%
  select(taxon, chem_gr, taxon_label, chem_gr_label,
         hypothesis, model_id, hypothesis_label,
         max_rhat, min_bulk_ess, min_tail_ess,
         n_divergent, max_treedepth)

glimpse(model_diagnostics)

write_csv(model_diagnostics, file.path(output_dir, "model_diagnostics_summary.csv"))

#=======4. LOO analysis per taxonomic group and chemical-use group==============
mem.maxVSize(Inf)

# Function to run LOO on one model file at a time
run_loo_from_file <- function(file, ndraws_loo = 1000) {
  fit <- readRDS(file)
  loo_obj <- tryCatch(
    {
      loo(fit, moment_match = FALSE, ndraws = ndraws_loo,
        save_psis = FALSE)
    },
    error = function(e) {
      message("LOO failed for file: ", file)
      message("Reason: ", e$message)
      return(NULL)
    }
  )
  
  rm(fit)
  gc()
  loo_obj
}

# Compare H1-H3 for one taxon × chemical-use group
compare_loo_one_group <- function(df, ndraws_loo = 1000) {
  
  taxon_i <- unique(df$taxon)
  chem_i <- unique(df$chem_gr)
  
  message("Running LOO for: ", taxon_i, " / ", chem_i)
  
  # Check same sample size without keeping models in memory
  nobs_vec <- map_int(df$file, ~ {
    fit_tmp <- readRDS(.x)
    n <- nobs(fit_tmp)
    rm(fit_tmp)
    gc()
    n
  })
  
  # Run LOO one model at a time
  loo_list <- map(df$file, run_loo_from_file, ndraws_loo = ndraws_loo)
  names(loo_list) <- as.character(df$hypothesis)
  
  # If any LOO failed, skip this group
  if (any(map_lgl(loo_list, is.null))) {
    warning(glue("Skipping {taxon_i} / {chem_i} because at least one LOO failed."))
    return(tibble())
  }
  
  loo_cmp <- loo_compare(loo_list) %>%
    as.data.frame() %>%
    rownames_to_column("hypothesis") %>%
    as_tibble() %>%
    select(hypothesis, elpd_diff, se_diff)
  
  loo_summary <- df %>%
    mutate(
      nobs = nobs_vec,
      loo_obj = loo_list,
      elpd_loo = map_dbl(
        loo_obj,
        ~ .x$estimates["elpd_loo", "Estimate"]),
      se_elpd_loo = map_dbl(
        loo_obj,
        ~ .x$estimates["elpd_loo", "SE"]),
      p_loo = map_dbl(
        loo_obj,
        ~ .x$estimates["p_loo", "Estimate"]),
      se_p_loo = purrr::map_dbl(
        loo_obj,
        ~ .x$estimates["p_loo", "SE"]),
      max_pareto_k = map_dbl(
        loo_obj,
        ~ max(.x$diagnostics$pareto_k, na.rm = TRUE)),
      n_bad_pareto_k = map_int(
        loo_obj,
        ~ sum(.x$diagnostics$pareto_k > 0.7, na.rm = TRUE))
    ) %>%
    select(-loo_obj) %>%
    left_join(loo_cmp, by = "hypothesis") %>%
    arrange(desc(elpd_loo)) %>%
    mutate(
      rank = row_number(),
      selected = rank == 1,
      selected_hypothesis = hypothesis[selected][1],
      selected_model_id = model_id[selected][1],
      selected_hypothesis_label = hypothesis_label[selected][1]
    )
  rm(loo_list)
  gc()
  
  loo_summary
}

# Run LOO comparison
loo_results <- model_files %>%
  group_by(taxon, chem_gr) %>%
  group_split() %>%
  map_dfr(compare_loo_one_group, ndraws_loo = 2000) %>%
  select(taxon, chem_gr, taxon_label, chem_gr_label, 
         hypothesis, model_id, hypothesis_label,
         elpd_loo, se_elpd_loo, elpd_diff, se_diff, p_loo, se_p_loo,
         max_pareto_k, n_bad_pareto_k,
         rank, selected, selected_hypothesis) %>%
  arrange(taxon, chem_gr, rank)

write_csv(loo_results, file.path(output_dir, "LOO_model_comparisons.csv"))

#=======5. Extract posterior estimates for the best-supported model=============
# use LOO-selected best model per taxon and chemical group.
best_model_files <- loo_results %>%
  filter(selected) %>%
  distinct(taxon, chem_gr, hypothesis, .keep_all = TRUE) # all H3

best_model_fits <- model_fits %>%
  inner_join(best_model_files %>%
               select(taxon, chem_gr, hypothesis, model_id),
             by = c("taxon", "chem_gr", "hypothesis", "model_id"))

# Fixed-effect posterior summaries for best models only
posterior_fixed_effects_best <- best_model_fits %>%
  mutate(fixed_effects = map(fit, ~ fixef(.x, probs = c(0.025, 0.975)) %>%
                          as.data.frame() %>%
                          rownames_to_column("parameter") %>%
                          as_tibble())) %>%
  select(taxon, chem_gr, taxon_label, chem_gr_label,
         hypothesis, model_id, hypothesis_label, fixed_effects) %>%
  unnest(fixed_effects)

write_csv(posterior_fixed_effects_best,
  file.path(output_dir, "posterior_fixed_effect_estimates_best_model.csv"))

#=======6. Standardised FFG-specific slopes ====================================
calc_hdi_rope <- function(x, rope = 0.20, ci = 0.95) {
  
  x <- x[is.finite(x)]
  
  if (length(x) == 0) {
    return(tibble(
      hdi_low = NA_real_,
      hdi_high = NA_real_,
      hdi_width = NA_real_,
      hdi_rope_prop = NA_real_
    ))
  }
  
  hdi_x <- bayestestR::hdi(x, ci = ci)
  
  hdi_low <- hdi_x$CI_low[1]
  hdi_high <- hdi_x$CI_high[1]
  hdi_width <- hdi_high - hdi_low
  
  overlap_width <- max(
    0,
    min(hdi_high, rope) - max(hdi_low, -rope)
  )
  
  tibble(
    hdi_low = hdi_low,
    hdi_high = hdi_high,
    hdi_width = hdi_width,
    hdi_rope_prop = overlap_width / hdi_width
  )
}

summarise_posterior_effect <- function(x, thresholds = c(0.20, 0.50), ci = 0.95) {
  
  x <- x[is.finite(x)]
  
  hdi_main <- calc_hdi_rope(x, rope = thresholds[1], ci = ci)
  
  out <- tibble(
    median = median(x, na.rm = TRUE),
    q025 = quantile(x, 0.025, na.rm = TRUE),
    q975 = quantile(x, 0.975, na.rm = TRUE),
    hdi_low = hdi_main$hdi_low,
    hdi_high = hdi_main$hdi_high,
    hdi_width = hdi_main$hdi_width,
    pd = max(
      mean(x > 0, na.rm = TRUE),
      mean(x < 0, na.rm = TRUE)
    )
  )
  
  for (thr in thresholds) {
    
    thr_name <- sprintf("%.2f", thr)
    rope_res <- calc_hdi_rope(x, rope = thr, ci = ci)
    
    out[[paste0("prob_abs_gt_", thr_name)]] <- mean(abs(x) > thr, na.rm = TRUE)
    out[[paste0("hdi_rope_prop_", thr_name)]] <- rope_res$hdi_rope_prop
    out[[paste0("hdi_inside_", thr_name)]] <- rope_res$hdi_rope_prop == 1
    out[[paste0("hdi_outside_", thr_name)]] <- rope_res$hdi_rope_prop == 0
  }
  
  out
}

extract_standardised_slopes <- function(fit, taxon, chem_gr, taxon_label, chem_gr_label) {
  
  dat <- fit$data
  
  sd_mass <- sd(dat$log_mass_c, na.rm = TRUE)
  sd_sens <- sd(dat$log_rel_sensi_c, na.rm = TRUE)
  
  if (!is.finite(sd_mass) || !is.finite(sd_sens) || sd_mass == 0 || sd_sens == 0) {
    warning(glue::glue("Skipping {taxon} / {chem_gr}: invalid SD."))
    return(tibble())
  }
  
  draws <- posterior::as_draws_df(fit) %>%
    as_tibble() %>%
    mutate(.draw = row_number())
  
  slope_cols <- names(draws) %>%
    str_subset("^b_FFG.*:log_mass_c$")
  
  if (length(slope_cols) == 0) {
    warning(glue::glue("Skipping {taxon} / {chem_gr}: no FFG-specific slopes."))
    return(tibble())
  }
  
  draws %>%
    select(.draw, all_of(slope_cols)) %>%
    pivot_longer(
      cols = all_of(slope_cols),
      names_to = "parameter",
      values_to = "beta_raw"
    ) %>%
    mutate(
      parameter = str_remove(parameter, "^b_"),
      FFG = parameter %>%
        str_remove(":log_mass_c$") %>%
        str_remove("^FFG"),
      beta_std = beta_raw * sd_mass / sd_sens,
      taxon = taxon,
      chem_gr = chem_gr,
      taxon_label = taxon_label,
      chem_gr_label = chem_gr_label
    ) %>%
    select(
      taxon, chem_gr, taxon_label, chem_gr_label,
      FFG, parameter, .draw, beta_raw, beta_std
    )
}

standardised_slope_draws <- best_model_fits %>%
  select(fit, taxon, chem_gr, taxon_label, chem_gr_label) %>%
  pmap_dfr(extract_standardised_slopes)

standardised_slope_summary <- standardised_slope_draws %>%
  group_by(taxon, chem_gr, taxon_label, chem_gr_label, FFG, parameter) %>%
  summarise(
    effect = list(summarise_posterior_effect(
      beta_std,
      thresholds = c(0.10, 0.20, 0.50)
    )),
    .groups = "drop"
  ) %>%
  unnest(effect) %>%
  rename(
    beta_std_median = median,
    beta_std_q025 = q025,
    beta_std_q975 = q975
  ) %>%
  mutate(
    effect_marker = case_when(
      hdi_outside_0.50 ~ "**",
      hdi_outside_0.20 ~ "*",
      TRUE ~ ""
    )
  )

write_csv(standardised_slope_summary,
  file.path(output_dir, "standardised_FFG_slope_summary.csv"))

standardised_slope_general_summary <- standardised_slope_summary %>%
  summarise(
    n_slopes = n(),
    median_abs_beta_std = median(abs(beta_std_median), na.rm = TRUE),
    
    n_hdi_inside_0.20 = sum(hdi_inside_0.20, na.rm = TRUE),
    n_hdi_outside_0.20 = sum(hdi_outside_0.20, na.rm = TRUE),
    n_hdi_outside_0.50 = sum(hdi_outside_0.50, na.rm = TRUE),
    
    prop_hdi_inside_0.20 = n_hdi_inside_0.20 / n_slopes,
    prop_hdi_outside_0.20 = n_hdi_outside_0.20 / n_slopes,
    prop_hdi_outside_0.50 = n_hdi_outside_0.50 / n_slopes
  )

write_csv(standardised_slope_general_summary,
  file.path(output_dir, "general_standardised_FFG_slope_summary.csv"))

#=======7. Pairwise FFG slope and baseline contrasts ===========================
# check pairwise difference in standardised slopes between FFGs.
extract_pairwise_slope_contrasts <- function(df) {
  
  slope_wide <- df %>%
    select(.draw, FFG, beta_std) %>%
    pivot_wider(
      names_from = FFG,
      values_from = beta_std
    )
  
  ffg_levels <- setdiff(names(slope_wide), ".draw")
  
  if (length(ffg_levels) < 2) {
    return(tibble())
  }
  
  combn(ffg_levels, 2, simplify = FALSE) %>%
    map_dfr(\(pair) {
      
      FFG_1 <- pair[1]
      FFG_2 <- pair[2]
      
      slope_wide %>%
        transmute(
          FFG_1 = FFG_1,
          FFG_2 = FFG_2,
          .draw = .draw,
          beta_std_1 = .data[[FFG_1]],
          beta_std_2 = .data[[FFG_2]],
          slope_diff = beta_std_1 - beta_std_2
        )
    })
}

pairwise_slope_draws <- standardised_slope_draws %>%
  group_by(taxon, chem_gr, taxon_label, chem_gr_label) %>%
  group_modify(~ extract_pairwise_slope_contrasts(.x)) %>%
  ungroup()

pairwise_slope_summary <- pairwise_slope_draws %>%
  group_by(taxon, chem_gr, taxon_label, chem_gr_label, FFG_1, FFG_2) %>%
  summarise(
    effect = list(summarise_posterior_effect(
      slope_diff,
      thresholds = c(0.20, 0.50)
    )),
    .groups = "drop"
  ) %>%
  unnest(effect) %>%
  rename(
    slope_diff_median = median,
    slope_diff_q025 = q025,
    slope_diff_q975 = q975
  ) %>%
  mutate(
    effect_marker = case_when(
      hdi_outside_0.50 ~ "**",
      hdi_outside_0.20 ~ "*",
      TRUE ~ ""
    )
  )

write_csv(pairwise_slope_summary,
  file.path(output_dir, "pairwise_FFG_slope_contrast_summary.csv"))



extract_pairwise_intercept_contrasts <- function(fit, taxon, chem_gr,
                                                 taxon_label, chem_gr_label) {
  
  draws <- posterior::as_draws_df(fit) %>%
    as_tibble() %>%
    mutate(.draw = row_number())
  
  intercept_cols <- names(draws) %>%
    str_subset("^b_FFG") %>%
    str_subset(":log_mass_c$", negate = TRUE)
  
  if (length(intercept_cols) < 2) {
    warning(glue::glue("Skipping {taxon} / {chem_gr}: fewer than two FFG intercepts."))
    return(tibble())
  }
  
  if (!"sigma" %in% names(draws)) {
    warning(glue::glue("Skipping {taxon} / {chem_gr}: sigma not found."))
    return(tibble())
  }
  
  intercept_wide <- draws %>%
    select(.draw, sigma, all_of(intercept_cols)) %>%
    rename_with(
      ~ .x %>%
        str_remove("^b_") %>%
        str_remove("^FFG"),
      all_of(intercept_cols)
    )
  
  ffg_levels <- setdiff(names(intercept_wide), c(".draw", "sigma"))
  
  combn(ffg_levels, 2, simplify = FALSE) %>%
    map_dfr(\(pair) {
      
      FFG_1 <- pair[1]
      FFG_2 <- pair[2]
      
      intercept_wide %>%
        transmute(
          taxon = taxon,
          chem_gr = chem_gr,
          taxon_label = taxon_label,
          chem_gr_label = chem_gr_label,
          FFG_1 = FFG_1,
          FFG_2 = FFG_2,
          .draw = .draw,
          alpha_1 = .data[[FFG_1]],
          alpha_2 = .data[[FFG_2]],
          contrast_raw = alpha_1 - alpha_2,
          sigma = sigma,
          d_baseline = contrast_raw / sigma
        )
    })
}

pairwise_intercept_draws <- best_model_fits %>%
  select(fit, taxon, chem_gr, taxon_label, chem_gr_label) %>%
  pmap_dfr(extract_pairwise_intercept_contrasts)

pairwise_intercept_summary <- pairwise_intercept_draws %>%
  group_by(taxon, chem_gr, taxon_label, chem_gr_label, FFG_1, FFG_2) %>%
  summarise(
    contrast_raw_median = median(contrast_raw, na.rm = TRUE),
    contrast_raw_q025 = quantile(contrast_raw, 0.025, na.rm = TRUE),
    contrast_raw_q975 = quantile(contrast_raw, 0.975, na.rm = TRUE),
    
    effect = list(summarise_posterior_effect(
      d_baseline,
      thresholds = c(0.20, 0.50, 0.80)
    )),
    .groups = "drop"
  ) %>%
  unnest(effect) %>%
  rename(
    d_baseline_median = median,
    d_baseline_q025 = q025,
    d_baseline_q975 = q975
  ) %>%
  mutate(
    effect_marker = case_when(
      hdi_outside_0.80 ~ "***",
      hdi_outside_0.50 ~ "**",
      hdi_outside_0.20 ~ "*",
      TRUE ~ ""
    )
  )

write_csv(pairwise_intercept_summary,
  file.path(output_dir, "pairwise_FFG_intercept_contrast_summary.csv"))




#=======8. Marginal/conditional R2 and random variation partitioning============
n_draws_r2  <- 2000
chunk_size  <- 2000
response_var <- "log_rel_sensi_c"

calc_r2_chunked <- function(fit, response_var, ndraws = 1000, chunk_size = 1000, re_formula = NA) {
  
  dat <- fit$data
  y <- dat[[response_var]]
  n <- length(y)
  
  n_available <- posterior::ndraws(fit)
  draw_ids <- sample(seq_len(n_available), size = min(ndraws, n_available))
  d <- length(draw_ids)
  
  sum_mu <- sum_mu2 <- sum_res <- sum_res2 <- numeric(d)
  chunks <- split(seq_len(n), ceiling(seq_len(n) / chunk_size))
  
  for (idx in chunks) {
    mu_i <- posterior_epred(fit, newdata = dat[idx, , drop = FALSE],
      re_formula = re_formula, summary = FALSE, draw_ids = draw_ids) %>%
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

r2_marginal_conditional_best <- best_model_fits %>%
  mutate(
    r2 = map(fit, 
             \(x) calc_r2_from_fit(fit = x, response_var = response_var, 
                                   ndraws = n_draws_r2, chunk_size = chunk_size)
             )) %>%
  select(taxon, chem_gr, taxon_label, chem_gr_label,
         hypothesis, model_id, hypothesis_label, formula_structure, r2) %>%
  unnest(r2)

glimpse(r2_marginal_conditional_best)
write_csv(r2_marginal_conditional_best,
          file.path(output_dir, paste0("r2_marginal_conditional_best_model.csv")))


extract_random_variance_partition <- function(fit) {
  
  draws <- posterior::as_draws_df(fit) %>%
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

random_variance_partition_best <- best_model_fits %>%
  mutate(variance_components = map(fit, extract_random_variance_partition)) %>%
  select(
    taxon, chem_gr, taxon_label, chem_gr_label,
    hypothesis, model_id, hypothesis_label, formula_structure,
    variance_components
  ) %>%
  unnest(variance_components)

glimpse(random_variance_partition_best)

write_csv(random_variance_partition_best, 
          file.path(output_dir, "random_variance_partition_best_model.csv"))


