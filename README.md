# Bayesian Hierarchical Spatial Difference-in-Difference Models

This repository contains code to reproduce the simulation experiment tables presented in the Appendix of our paper, *Bayesian Hierarchical Spatial Difference-in-Difference Models*. The R Markdown file **reproducibility-guide-spatial-DiD-simulation-code.Rmd** includes all code used to obtain the posterior samples, credible intervals, and model performance metrics reported in these tables.

## Getting Started

1. Download the ZIP file and open the folder titled `spatialDiD-simulation-code-main`
2. Open the R Project file titled `spatialDiD-simulation-code-tutorial.Rproj`
3. Open the R Markdown file titled `reproducibility-guide-spatial-DiD-simulation-code.Rmd`
4. Knit the file to HTML — the code will begin to run automatically

## Repository Structure

The following folders are required to run the code and must remain in the same directory as the R Project file:

- `dgp-cls-homoskedastic` — simulated data for the traditional DiD model
- `dgp-us-homoskedastic` — simulated data for the uniform spatial DiD model
- `dgp-svc-homoskedastic` — simulated data for the spatially varying coefficient DiD model

Each folder contains:

- `data/` — the simulated dataset for that model type, along with the adjacency matrix for California counties
- `src/` — JAGS model files used to fit each of the three models to the simulated data, along with the original scripts used to generate and prepare the data. 
- `test/` — original R scripts for fitting each model (`model_fit/`) and evaluating model performance using the WAIC (`model_choice/`). 

The data generation and preparation scripts included in the `src/` folder, and the model fit and model choice scripts included in the `test/` folder are not needed to run the reproducibility guide but are provided for users who wish to view these materials or individually run each model on each dataset and compare model performance.

## Requirements

This project requires the following software and packages:

- **JAGS:** must be installed on your machine before running the code
- **R packages:** `tidyverse`, `rjags`, `coda`, `kableExtra`, `loo`

## Just Want to View the Output?

To view the reproduced tables without re-running the script, simply open the HTML file included in the repository in your browser.

---

*Questions, comments, or concerns? Feel free to reach out to me (Valentina) at valentina214@g.ucla.edu*!

