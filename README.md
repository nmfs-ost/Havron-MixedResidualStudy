
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->
<!-- badges: end -->

Residual analysis for mixed effects models. Code to implement this
projects can be found on
[Github](https://github.com/Cole-Monnahan-NOAA/mixed_resids)

Manuscript:

- [Draft](https://docs.google.com/document/d/19Y39GqVRAmoIEegxgzyf6HXVWYgWGyr5x8cyHe3qHK8/edit)
- [Tables](https:://andrea-havron-github-io/tables/TMB-validation-tables.pdf)
- [Figures](https:://andrea-havron-github-io/articles/TMB-validation-figure.html)

Project Repo Structure:

- [code](https://github.com/Cole-Monnahan-NOAA/mixed_resids/tree/main/code):
  Demo code and test snippets

- [R](https://github.com/Cole-Monnahan-NOAA/mixed_resids/tree/main/R): R
  scripts used to run models

  - [make_plots](https://github.com/Cole-Monnahan-NOAA/mixed_resids/blob/main/R/make_plots.R):
    Initial code used to make plots. Final code now moved to docs/

  - [model_fns](https://github.com/Cole-Monnahan-NOAA/mixed_resids/blob/main/R/model_fns.R):
    Core functions used to run model iterations

  - [resid_fns](https://github.com/Cole-Monnahan-NOAA/mixed_resids/blob/main/R/resid_fns.R):
    Functions used to calculate residuals

  - [run_analysis](https://github.com/Cole-Monnahan-NOAA/mixed_resids/blob/main/R/run_analysis.R):
    Script used to run analysis in parallel

  - [run_sample_sizes](https://github.com/Cole-Monnahan-NOAA/mixed_resids/blob/main/R/run_sample_sizes.R):
    Script used to run models for increasing sample sizes

  - [sim_data](https://github.com/Cole-Monnahan-NOAA/mixed_resids/blob/main/R/sim_data.R):
    Functions used to simulate data

  - [startup](https://github.com/Cole-Monnahan-NOAA/mixed_resids/blob/main/R/startup.R):
    Prepares the R environment for the analysis

  - [utils](https://github.com/Cole-Monnahan-NOAA/mixed_resids/blob/main/R/utils.R):
    Helper functions

- [results](https://github.com/Cole-Monnahan-NOAA/mixed_resids/tree/main/results):
  Saved results as .RMD files

- [src](https://github.com/Cole-Monnahan-NOAA/mixed_resids/tree/main/src):
  TMB models

This is an ongoing analysis, please check back later for more
information.
