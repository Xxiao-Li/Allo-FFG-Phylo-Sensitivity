#!/usr/bin/env Rscript

# ============================================================
# brms pooled models per taxonomic group (algae / invertebrate / vertebrate)
# One SLURM array task = ONE model (MODEL_ID) for ONE taxon (TAXON)
#
# Task mapping (default):
#   taxa = 3 (algae, invertebrate, vertebrate)
#   models = 8 (M0.1 ... M3.2)
#   TASK_ID = SLURM_ARRAY_TASK_ID (or TASK_ID env)
#
#   taxon_index = ((TASK_ID - 1) %/% 8) + 1
#   model_id    = ((TASK_ID - 1) %%  8) + 1
#
# Models:
#   1 M0.1: y ~ 1 + chem_group + (1 | CAS) + (1 | order/family)
#   2 M0.2: + (1 | gr(genus, cov = A_genus))
#   3 M1.1: + z_log_mass
#   4 M1.2: + chem_group * z_log_mass
#   5 M2.1: + z_log_mass + FFG
#   6 M2.2: + chem_group*z_log_mass + chem_group*FFG
#   7 M3.1: + z_log_mass*FFG
#   8 M3.2: + chem_group*FFG*z_log_mass
#
# Inputs (as provided):
#   trees:
#     here("data/6_combined/algae_genus_tree.rds")
#     here("data/6_combined/invertebrate_genus_tree.rds")
#     here("data/6_combined/vertebrate_genus_tree.rds")
#   data:
#     here("data/6_combined/alage_mass_FFG_sens.rds")   # (typo in filename kept)
#     here("data/6_combined/invertebrate_mass_FFG_sens.rds")
#     here("data/6_combined/vertebrate_mass_FFG_sens.rds")
#
# Output:
#   - saves ONLY the fitted brms model object (fit) as .rds
# ============================================================

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(ape)
  library(brms)
  library(cmdstanr)
})

# --------------------------
# SLURM/task mapping
# --------------------------
TASK_ID <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", Sys.getenv("TASK_ID", "1")))

taxa <- c("algae", "invertebrate", "vertebrate")
N_MODELS <- 8L

taxon_index <- ((TASK_ID - 1L) %/% N_MODELS) + 1L
MODEL_ID    <- ((TASK_ID - 1L) %%  N_MODELS) + 1L
stopifnot(taxon_index >= 1L, taxon_index <= length(taxa))
stopifnot(MODEL_ID %in% 1:8)

TAXON <- taxa[taxon_index]

# --------------------------
# knobs (slurm env optional)
# --------------------------
CHAINS <- as.integer(Sys.getenv("CHAINS", "4"))
CORES  <- as.integer(Sys.getenv("CORES",  "4"))
THREADS <- as.integer(Sys.getenv("THREADS","8"))
GRAINSIZE <- as.integer(Sys.getenv("GRAINSIZE","0"))

ITER   <- as.integer(Sys.getenv("ITER",   "10000"))
WARMUP <- as.integer(Sys.getenv("WARMUP", "5000"))
ADAPT_DELTA <- as.numeric(Sys.getenv("ADAPT_DELTA", "0.99"))
MAX_TREEDEPTH <- as.integer(Sys.getenv("MAX_TREEDEPTH", "18"))
SEED <- as.integer(Sys.getenv("SEED", "1007"))

# response column
RESP <- Sys.getenv("RESP", "log_rel_sensi")

Sys.setenv(
  STAN_NUM_THREADS = THREADS,
  OMP_NUM_THREADS  = 1,
  MKL_NUM_THREADS  = 1,
  OPENBLAS_NUM_THREADS = 1
)

thr  <- if (GRAINSIZE > 0) brms::threading(threads = THREADS, grainsize = GRAINSIZE) else brms::threading(threads = THREADS)
ctrl <- list(adapt_delta = ADAPT_DELTA, max_treedepth = MAX_TREEDEPTH)

# --------------------------
# input files (fixed)
# --------------------------
tree_files <- list(
  algae         = here("data/6_combined/algae_genus_tree.rds"),
  invertebrate  = here("data/6_combined/invertebrate_genus_tree.rds"),
  vertebrate    = here("data/6_combined/vertebrate_genus_tree.rds")
)

data_files <- list(
  algae         = here("data/6_combined/alage_mass_FFG_sens.rds"),   
  invertebrate  = here("data/6_combined/invertebrate_mass_FFG_sens.rds"),
  vertebrate    = here("data/6_combined/vertebrate_mass_FFG_sens.rds")
)

# --------------------------
# output directory
# --------------------------
out_dir <- here("model_output", "bayesian_models", "pooled_M1_M8")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --------------------------
# load inputs
# --------------------------
tree0 <- readRDS(tree_files[[TAXON]])

df_path <- data_files[[TAXON]]

meta0 <- readRDS(df_path)

# --------------------------
# prepare modelling data
# --------------------------
dat <- meta0 %>%
  transmute(
    y          = .data[[RESP]],
    chem_group = factor(chem_group),
    CAS        = factor(as.character(CAS)),
    genus      = factor(as.character(genus)),
    order      = factor(as.character(order)),
    family     = factor(as.character(family)),
    FFG        = factor(as.character(FG5)),
    log_mass   = as.numeric(log_mass)
  ) %>%
  droplevels()

# global y-centering within taxon dataset (keeps chem_group interpretable)
y_mean <- mean(dat$y, na.rm = TRUE)
dat <- dat %>% mutate(y = y - y_mean)

# z-score log_mass within taxon dataset
dat <- dat %>% mutate(z_log_mass = as.numeric(scale(log_mass)))

# --------------------------
# build A_genus from pruned tree (Grafen brlen), then cov2cor
# (Only needed for models 2..8)
# --------------------------
A_genus <- NULL
if (MODEL_ID != 1) {
  g_all <- sort(unique(as.character(dat$genus)))
  keep <- intersect(g_all, tree0$tip.label)

  tree_sub <- ape::keep.tip(tree0, keep)
  tree_sub <- ape::multi2di(tree_sub)
  tree_sub <- ape::compute.brlen(tree_sub, method = "Grafen")
  
  V_genus <- ape::vcv(tree_sub)
  A_genus <- cov2cor(V_genus)
  
  # align genus factor levels to A_genus
  dat <- dat %>%
    filter(as.character(genus) %in% rownames(A_genus)) %>%
    mutate(genus = factor(as.character(genus), levels = rownames(A_genus))) %>%
    droplevels()
  
  stopifnot(identical(levels(dat$genus), rownames(A_genus)))
}

# --------------------------
# model formulas (always use order/family)
# --------------------------
RE_no_phylo <- "(1 | CAS) + (1 | order/family)"
RE_phylo    <- "(1 | CAS) + (1 | order/family) + (1 | gr(genus, cov = A_genus))"

form_str <- switch(
  as.character(MODEL_ID),
  "1" = paste0("y ~ 1 + chem_group + ", RE_no_phylo),                       # M0.1
  "2" = paste0("y ~ 1 + chem_group + ", RE_phylo),                          # M0.2
  "3" = paste0("y ~ 1 + chem_group + z_log_mass + ", RE_phylo),             # M1.1
  "4" = paste0("y ~ 1 + chem_group * z_log_mass + ", RE_phylo),             # M1.2
  "5" = paste0("y ~ 1 + chem_group + z_log_mass + FFG + ", RE_phylo),       # M2.1
  "6" = paste0("y ~ 1 + chem_group * z_log_mass + chem_group * FFG + ", RE_phylo), # M2.2
  "7" = paste0("y ~ 1 + chem_group + z_log_mass * FFG + ", RE_phylo),       # M3.1
  "8" = paste0("y ~ 1 + chem_group * FFG * z_log_mass + ", RE_phylo)        # M3.2
)

form <- bf(as.formula(form_str))

model_lab <- sprintf("M%02d", MODEL_ID)
cat("\n============================\n")
cat("TASK_ID   :", TASK_ID, "\n")
cat("TAXON     :", TAXON, "\n")
cat("MODEL_ID  :", MODEL_ID, " (", model_lab, ")\n", sep = "")
cat("FORMULA   :", form_str, "\n")
cat("n         :", nrow(dat), "\n")
cat("levels    : cg=", nlevels(dat$chem_group),
    " FFG=", nlevels(dat$FFG),
    " genus=", nlevels(dat$genus),
    " CAS=", nlevels(dat$CAS), "\n", sep = "")
cat("============================\n\n")


# --------------------------
# priors (all RE SDs = 0.5 * sdy)
# --------------------------
sdy <- sd(dat$y, na.rm = TRUE)
if (!is.finite(sdy) || sdy <= 0) stop("Non-positive SD for centered y.")

# fixed RE SD scale (half-normal implied by normal(0, ...) on sd class)
SD_SCALE <- 0.5 * sdy

pri <- c(
  # Intercept is grand mean because y is globally centered
  set_prior("normal(0, 0.5)", class = "Intercept"),
  
  # all fixed effects (chem_group, FFG, z_log_mass, interactions)
  set_prior(sprintf("normal(0, %.6f)", sdy), class = "b"),
  
  # residual SD: mean ~ sdy under gamma(2, 2/sdy)
  set_prior(sprintf("gamma(2, %.6f)", 2 / sdy), class = "sigma"),
  
  # random-effect SDs (all 0.5*sdy)
  set_prior(sprintf("normal(0, %.6f)", SD_SCALE), class = "sd", group = "CAS"),
  set_prior(sprintf("normal(0, %.6f)", SD_SCALE), class = "sd", group = "order"),
  set_prior(sprintf("normal(0, %.6f)", SD_SCALE), class = "sd", group = "order:family")
)

# add phylogenetic genus SD prior only when phylogeny term is present (models 2..8)
if (MODEL_ID != 1) {
  pri <- c(pri, set_prior(sprintf("normal(0, %.6f)", SD_SCALE), class = "sd", group = "genus"))
}

print(pri)


# --------------------------
# fit + save (save ONLY fit object)
# --------------------------
fit_tag <- paste0(
  TAXON, "_", model_lab,
  "_n", nrow(dat),
  "_cg", nlevels(dat$chem_group),
  "_FFG", nlevels(dat$FFG),
  "_g", nlevels(dat$genus),
  "_CAS", nlevels(dat$CAS)
)

fit_file <- file.path(out_dir, paste0(fit_tag, ".rds"))

fit <- brm(
  form,
  data    = dat,
  data2   = if (MODEL_ID == 1) NULL else list(A_genus = A_genus),
  family  = gaussian(),
  prior   = pri,
  backend = "cmdstanr",
  chains  = CHAINS, cores = CORES,
  threads = thr,
  iter    = ITER, warmup = WARMUP,
  control = ctrl,
  seed    = SEED,
  refresh = 100
)

print(summary(fit))
saveRDS(fit, fit_file)
message("Saved fit: ", fit_file)
