# spaCRPower: Simulation and Analysis of spatial phenotype analysis of CRISPR-Cas9 screens

This is an R package for simulating pooled genetic screens and single cell phenotyping. The question
is for a given library size and a cell-level classifier, on average how many cells of different genotypes
and how many wells are needed to accurately classify a genotype as having the phenotype?

<img width="1182" alt="image" src="https://github.com/user-attachments/assets/af14ee26-6afa-4aed-a2ba-64692cf205af" />


## Quick Setup

    install.packages("remotes"
    remotes::install_github("maomlab/spaCRPower")

To begin, simulate a screen

    data <- spaCRPower::simulate_screen(
      n_genes_per_library = 200,
      gene_abundance_factor_mu = 1,
      gene_abundance_factor_var = .01,
      gene_hit_rate = 0.05,
      n_wells_per_screen = 80,
      well_abundance_factor_mu = 1,
      well_abundance_factor_var = .01,
      n_genes_per_well = 10,
      n_cells_per_well_mu = 850,
      n_cells_per_well_var = 1000,
      class_pos_mu = 0.99,
      class_pos_var = 0.0001,
      class_neg_mu = 0.01,
      class_neg_var = 0.0001)

investigate properties of the simulated screen

    library(patchwork)
    plot1 <- plot_genes_per_well(data)
    plot2 <- plot_wells_per_gene(data)
    plot3 <- plot_positivity_rate_by_well(data)
    plot4 <- plot_well_gene_heatmap(data)
    (plot1 + plot2) / (plot3 + plot4)

<img width="1147" alt="image" src="https://github.com/user-attachments/assets/226cf378-267c-432f-b1b1-4a37062a9209" />

predict the which genotypes are hits:

    model_data <- prepare_model_data(data)
    model <- compile_model(model_data)
    model_fit <- fit_model(model_data, model)
    model_estimate <- gather_model_estimate(model_fit)
    model_evaluation <- evaluate_model_fit(data, model_estimate)

Or sweep over parameters:

    scores <- spaCRPower::scan_parameters(
      n_genes_per_library = 200,
      gene_abundance_factor_mu = 1,
      gene_abundance_factor_var = .01,
      gene_hit_rate = 0.05,
      n_wells_per_screen = seq(5, 500, length.out = 100),
      well_abundance_factor_mu = 1,
      well_abundance_factor_var = .01,
      n_genes_per_well = c(1, 3, 10, 30),
      n_cells_per_well_mu = 850,
      n_cells_per_well_var = 1000,
      class_pos_mu = 0.99,
      class_pos_var = 0.0001,
      class_neg_mu = 0.01,
      class_neg_var = 0.0001,
      progress_file = "parameter_scan.tsv",
      verbose = TRUE,
      refresh = 0)

![parameter_scan_plot](https://github.com/user-attachments/assets/dd13bbb5-e4f6-4981-bfba-47bbe6c30373)



