# spatialDiD-simulation-code

The code included in this repository is designed to reproduce the tables initially generated for the simulation experiment provided in the Appendix of our paper, *Bayesian Hierarchical Spatial Difference-in-Difference Models*. For transparency, the code used to obtain the posterior samples, credible intervals, and model performance metrics presented in these tables is also included the R Markdown file titled **reproducibility-guide-spatial-DiD-simulation-code.Rmd**. To create the document, simply:

1. Download the ZIP file and open the folder titled 'spatialDiD-simulation-code-main.'
2. Open the `R Project` file titled 'spatialDiD-simulation-code-tutorial.Rproj.'
3. Once inside the R project, open the R Markdown file titled 'reproducibility-guide-spatial-DiD-simulation-code.Rmd.'
4. Knit the file to HTML and the code should begin to run.

Please note that the Bayesian modeling in this study was done using JAGS, and thus requires installation of the JAGS software on the computer of the user. Several R packages are also required, including `tidyverse`, `rjags`, `coda`, `kableExtra` and `loo`.

If you wish to simply view the code and reproduced tables without re-running the entirescript, simply open the HTML file included in the repository in your chosen browser.

Thank you! Please feel free to reach out to me (Valentina) at valentina214@g.ucla.edu if you have any questions, comments, or concerns! 
