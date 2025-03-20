# spaCRPower: Simulation and Analysis of spatial phenotype analysis of CRISPR-Cas9 screens

This is an R package for simulating pooled genetic screens and single cell phenotyping. The question
is for a given library size and a cell-level classifier, on average how many cells of different genotypes
and how many wells are needed to accurately classify a genotype as having the phenotype?

<img width="916" alt="image" src="https://github.com/user-attachments/assets/c80c2baf-9c31-402f-ad16-bff10172c2df" />

The analysis a Bayesian multiple-linear-regression with a horseshoe prior to estimate geneotypes built on the `brms` and `stan` ecosystem.


## Quick Setup

    install.packages("remotes"
    remotes::install_github("maomlab/spaCRPower")

Then look in the vignettes.


<img width="1147" alt="image" src="https://github.com/user-attachments/assets/226cf378-267c-432f-b1b1-4a37062a9209" />
