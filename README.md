# Allo-FFG-Phylo-Sensitivity

This repository contains the code and data for the paper **“Phylogenetic, allometric, and functional group structure in chemical sensitivity across aquatic taxa”.**

Chemical risk assessment across diverse communities remains challenging because most chemical-species pairs are empirically untested. Identifying phylogenetic, allometric, and functional group structure in chemical sensitivity could support cross-species extrapolation and community risk assessment, yet how these factors explain sensitivity across diverse taxa and chemicals remains poorly understood. We examined how phylogenetic relatedness, body mass, and functional group shape relative sensitivity across 580 aquatic species and 1103 chemicals spanning six chemical-use groups. We quantified phylogenetic signal in species relative sensitivity and fitted Bayesian phylogenetic hierarchical models to test global allometry, functional group-specific baseline sensitivity, and functional group-specific allometry.

## Folder structure

The `code` folder contains all R scripts.

The `data` folder contains raw dataset, processed data and model outputs. Large Bayesian fits are stored separately from the repository. Obtain the required datasets and saved fits from the associated data archive and preserve their expected paths. **Data archive:**

The `figures` folder contains the generated main and supplementary figures.

## Workflow

1.  **Compile observed toxicity data.** Harmonise toxicity records, species names, and CAS identifiers. Script: `code/1_compile_observed_toxicity.r`

2.  **Predict species sensitivity.** Generate predictions of species chemical sensitivity using the developed machine learning model in <https://github.com/rivm-syso/adore_pairwise_learning/>. Save the output as `data/2_predicted_sensitivity/predictions.csv`.

3.  **Prepare traits and analysis datasets.** Combine taxonomy, phylogenies, body mass, functional groups, chemical-use groups, and predictions. Calculate relative sensitivity and save datasets for algae, invertebrates, and fish.\
    Scripts: `code/2_phylogeny_mass_FG.r`; `code/3_combine_data.r`

4.  **Estimate phylogenetic signal.** Calculate chemical-specific Pagel’s λ for each taxonomic group. Script: `code/4_phylogenetic_signal_HPC.r`

5.  **Fit Bayesian models.** Compare H0–H3 using standardised sensitivity and body mass, with chemical-identity and phylogenetic random effects. Script: `code/5_H0_H3_all_HPC.r`

6.  **Summarise model outputs.** Check convergence, compare models, and summarise slopes, pairwise contrasts, Bayesian R², and random-effect variance. Script: `code/6_analyse_H0_H3_fits.r`

7.  **Test robustness to predicted estimates.** Repeat the λ and H3 analyses after excluding 5%, 10%, or 15% of records within each chemical. Scripts: `code/8_appendix_phylogenetic_signal_random_exclusion_test_HPC.r`; `code/8_appendix_H3_random_exclusion_HPC.r`

8.  **Analyse invertebrate lineages.** Compare λ before and after removing lineage means and estimate λ within major lineages. Script: `code/8_appendix_lambda_invertebrate_lineages.r`

9.  **Run phylogenetic prediction tests.** Use PhyloPars to predict invertebrate genus–chemical values after masking 30% or 50% of available values. Script: `code/8_appendix_phylopars_prediction.r`

10. **Generate figures.** Produce the main and supplementary figures from the saved datasets and analysis outputs. Scripts: `code/7_plot.r`; `code/9_appendix_plot.r`
