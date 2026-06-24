# Allo-FFG-Phylo-Sensitivity

This repository contains the analysis workflow for the paper **“Phylogenetic, allometric, and functional structure of species chemical sensitivity across aquatic taxa.”**

The workflow tests how phylogenetic relatedness, body mass, functional feeding groups (FFGs), and chemical-use groups structure variation in aquatic species chemical sensitivity. It integrates species sensitivity data with body mass, FFG classifications, chemical-use information, and genus-level phylogenies. The main analyses estimate phylogenetic signal, fit Bayesian hierarchical models, test the influence of water solubility, and evaluate phylogeny-based prediction using PhyloPars.

## Folder structure

The `code` folder contains all R scripts.

The `data` folder contains raw dataset, processed data and model outputs. Large Bayesian model-fit files are not stored directly in this repository and should be downloaded from the associated data archive and placed in the expected folders.

The `figures` folder contains the generated main and supplementary figures.

## Workflow

The analysis workflow consists of the following steps.

1.  Compile and harmonise observed species sensitivity data.\
    Script: `code/1_compile_observed_sensitivity_data.rmd`

2.  Predict species sensitivity using a machine-learning model.\
    The machine-learning model was developed in <https://doi.org/10.1021/acs.est.5c01289>, and the associated code is available at <https://github.com/rivm-syso/adore_pairwise_learning/>.

3.  Prepare species metadata for downstream analyses.\
    This step combines taxonomy, genus-level phylogenies, body mass, functional feeding groups, chemical-use groups, and predicted sensitivity data into taxon-specific datasets for algae, invertebrates, and fish.\
    Scripts: `code/2_prepare_phylogeny_mass_FFG.rmd`; `code/3_prepare_data_for_brms.rmd`

4.  Estimate phylogenetic signal.\
    Pagel’s λ is estimated for each CAS identity, taxonomic group, and chemical-use group.\
    Script: `code/4_phylogenetic_signal_HPC.R`

5.  Fit Bayesian hierarchical models.\
    Models are fitted separately for algae, invertebrates, and fish. The candidate models compare increasing biological structure, from a baseline model including chemical identity and phylogeny to models including body mass, FFG identity, and FFG-specific allometric slopes.\
    Scripts: `code/5_algae_H0_H3_HPC.R`; `code/5_invertebrate_H0_H3_HPC.R`; `code/5_fish_H0_H3_HPC.R`

6.  Summarise Bayesian model outputs.\
    This step includes model comparison, posterior slope and intercept summaries, pairwise FFG contrasts, Bayesian R², and random-effect variance partitioning.\
    Script: `code/6_analyse_bayesian_outputs.R`

7.  Repeat the Bayesian modelling workflow with chemical water solubility.\
    Chemical water solubility is included as an additional predictor, and the corresponding model outputs are summarised.\
    Scripts: `code/7_algae_H0_H3_HPC_ws.R`; `code/7_invertebrate_H0_H3_HPC_ws.R`; `code/7_fish_H0_H3_HPC_ws.R`; `code/7_analyse_bayesian_outputs_ws.R`

8.  Run PhyloPars prediction tests.\
    These tests evaluate genus-level phylogenetic prediction of species sensitivity.\
    Script: `code/8_phylopars_prediction_HPC.R`

9.  Generate the main and supplementary figures.\
    Script: `code/0_plot.R`
