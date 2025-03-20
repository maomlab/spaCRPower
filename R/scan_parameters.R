

#' Scan and evaluate simulated screens
#'
#' For each combination of given parameters, generate a simulated dataset
#' fit a model to it and estimate the predictive accuracy.
#'
#' @param n_genes_per_library integer number of genes in the library
#' @param gene_abundance_factor_mu numeric
#' @param gene_abundance_factor_var numeric
#' @param gene_hit_rate numeric
#' @param n_wells_per_screen integer
#' @param well_abundance_factor_mu numeric
#' @param well_abundance_factor_var numeric
#' @param n_genes_per_well integer
#' @param n_cells_per_well_mu numeric
#' @param n_cells_per_well_var numeric
#' @param class_pos_mu numeric
#' @param class_pos_var numeric
#' @param class_neg_mu numeric
#' @param class_neg_var numeric
#' @param verbose = FALSE
#'
#'
#' @export
scan_parameters <- function(
    n_genes_per_library,
    gene_abundance_factor_mu,
    gene_abundance_factor_var,
    gene_hit_rate,
    n_wells_per_screen,
    well_abundance_factor_mu,
    well_abundance_factor_var,
    n_genes_per_well,
    n_cells_per_well_mu,
    n_cells_per_well_var,
    class_pos_mu,
    class_pos_var,
    class_neg_mu,
    class_neg_var,
    verbose = FALSE,
    ...) {

  parameter_sets <- tidyr::expand_grid(
    n_genes_per_library = n_genes_per_library,
    gene_abundance_factor_mu = gene_abundance_factor_mu,
    gene_abundance_factor_var = gene_abundance_factor_var,
    gene_hit_rate = gene_hit_rate,
    n_wells_per_screen = n_wells_per_screen,
    well_abundance_factor_mu = well_abundance_factor_mu,
    well_abundance_factor_var = well_abundance_factor_var,
    n_genes_per_well = n_genes_per_well,
    n_cells_per_well_mu = n_cells_per_well_mu,
    n_cells_per_well_var = n_cells_per_well_var,
    class_pos_mu = class_pos_mu,
    class_pos_var = class_pos_var,
    class_neg_mu = class_neg_mu,
    class_neg_var = class_neg_var) |>
    dplyr::mutate(
      param_index = dplyr::row_number())

  if (verbose) {
    cat(
      "Generating and estimating ", nrow(parameter_sets), " ",
      "different datasets\n", sep = "")
  }

  model <- NULL

  parameter_sets |>
    dplyr::rowwise() |>
    dplyr::do({
      params <- .
      if (verbose) {
        cat(
          "params: \n",
          "  n_genes_per_library = ", params$n_genes_per_library[1], "\n",
          "  gene_abundance_factor_mu = ", params$gene_abundance_factor_mu[1], "\n",
          "  gene_abundance_factor_var = ", params$gene_abundance_factor_var[1], "\n",
          "  gene_hit_rate = ", params$gene_hit_rate[1], "\n",
          "  n_wells_per_screen = ", params$n_wells_per_screen[1], "\n",
          "  well_abundance_factor_mu = ", params$well_abundance_factor_mu[1], "\n",
          "  well_abundance_factor_var = ", params$well_abundance_factor_var[1], "\n",
          "  n_genes_per_well = ", params$n_genes_per_well[1], "\n",
          "  n_cells_per_well_mu = ", params$n_cells_per_well_mu[1], "\n",
          "  n_cells_per_well_var = ", params$n_cells_per_well_var[1], "\n",
          "  class_pos_mu = ", params$class_pos_mu[1], "\n",
          "  class_pos_var = ", params$class_pos_var[1], "\n",
          "  class_neg_mu = ", params$class_neg_mu[1], "\n",
          "  class_neg_var = ", params$class_neg_var[1], "\n", sep = "")
      }
      data <- simulate_screen(
        n_genes_per_library = params$n_genes_per_library[1],
        gene_abundance_factor_mu = params$gene_abundance_factor_mu[1],
        gene_abundance_factor_var = params$gene_abundance_factor_var[1],
        gene_hit_rate = params$gene_hit_rate[1],
        n_wells_per_screen = params$n_wells_per_screen[1],
        well_abundance_factor_mu = params$well_abundance_factor_mu[1],
        well_abundance_factor_var = params$well_abundance_factor_var[1],
        n_genes_per_well = params$n_genes_per_well[1],
        n_cells_per_well_mu = params$n_cells_per_well_mu[1],
        n_cells_per_well_var = params$n_cells_per_well_var[1],
        class_pos_mu = params$class_pos_mu[1],
        class_pos_var = params$class_pos_var[1],
        class_neg_mu = params$class_neg_mu[1],
        class_neg_var = params$class_neg_var[1])

      model_data <- prepare_model_data(data)

      if (is.null(model)) {
        if (verbose) {
          cat("Compiling model for the first time time...\n")
        }
        model <<- compile_model(model_data)
      }

      model_fit <- fit_model(model_data, model, ...)
      model_estimate <- gather_model_estimate(model_fit)
      model_evaluation <- evaluate_model_fit(data, model_estimate)

      if (verbose) {
        cat(
          "model_ap: ", model_evaluation$model_ap[1], " ",
          "model_auroc: ", model_evaluation$model_auroc[1],
          "\n", sep = "")
      }

      model_evaluation |>
        dplyr::mutate(
          parameters = list(params),
          data = list(data),
          model = list(model_fit))
    }) |>
    dplyr::ungroup()
}
