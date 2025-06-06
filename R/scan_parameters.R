

#' Scan and evaluate simulated screens
#'
#' For each combination of given parameters, generate a simulated dataset
#' fit a model to it and estimate the predictive accuracy.
#'
#' @param n_genes_per_library integer number of genes in the library
#' @param gene_abundance_alpha numeric
#' @param gene_hit_rate numeric
#' @param n_wells_per_screen integer
#' @param well_abundance_factor_mu numeric
#' @param well_abundance_factor_var numeric
#' @param imaging_n_cells_per_well_lambda numeric,
#' @param imaging_n_cells_per_well_nu numeric,
#' @param class_pos_mu numeric
#' @param class_pos_var numeric
#' @param class_neg_mu numeric
#' @param class_neg_var numeric
#' @param sequencing_n_cells_per_well_lambda numeric
#' @param sequencing_n_cells_per_well_nu numeric
#' @param pcr_factor_mu numeric
#' @param pcr_factor_var numeric
#' @param n_reads_total numeric
#' @param progress_file numeric
#' @param verbose bool
#'
#' @examples
#' \dontrun{
#'
#' scores <- spaCRPower::scan_parameters(
#'     n_genes_in_library = 200,
#'     gene_abundance_alpha = 1000,
#'     gene_hit_rate = 0.05,
#'     n_wells_per_screen = 20,
#'     well_abundance_factor_mu = 10,
#'     well_abundance_factor_var = 2,
#'     imaging_n_cells_per_well_lambda = 1000,
#'     imaging_n_cells_per_well_nu = 1,
#'     class_pos_mu = 0.99,
#'     class_pos_var = .0001,
#'     class_neg_mu = 0.01,
#'     class_neg_var = .0001,
#'     sequencing_n_cells_per_well_lambda = 1000,
#'     sequencing_n_cells_per_well_nu = 1,
#'     pcr_factor_mu = 0.1,
#'     pcr_factor_var = 0.01,
#'     n_reads_total = 1000,
#'     progress_file = "parameter_scan.tsv",
#'     verbose = TRUE,
#'     refresh = 0)
#' }
#'
#' @export
scan_parameters <- function(
    n_genes_in_library,
    gene_abundance_alpha,
    gene_hit_rate,
    n_wells_per_screen,
    well_abundance_factor_mu,
    well_abundance_factor_var,
    imaging_n_cells_per_well_lambda,
    imaging_n_cells_per_well_nu,
    class_pos_mu,
    class_pos_var,
    class_neg_mu,
    class_neg_var,
    sequencing_n_cells_per_well_lambda,
    sequencing_n_cells_per_well_nu,
    pcr_factor_mu,
    pcr_factor_var,
    n_reads_total,
    progress_file = NULL,
    verbose = FALSE,
    ...) {

  parameter_sets <- tidyr::expand_grid(
    n_genes_per_library = n_genes_per_library,
    gene_abundance_alpha = gene_abundance_alpha,
    gene_hit_rate = gene_hit_rate,
    n_wells_per_screen = n_wells_per_screen,
    well_abundance_factor_mu = well_abundance_factor_mu,
    well_abundance_factor_var = well_abundance_factor_var,
    imaging_n_cells_per_well_lambda = imaging_n_cells_per_well_lambda,
    imaging_n_cells_per_well_nu = imaging_n_cells_per_well_nu,
    class_pos_mu = class_pos_mu,
    class_pos_var = class_pos_var,
    class_neg_mu = class_neg_mu,
    class_neg_var = class_neg_var,
    sequencing_n_cells_per_well_lambda = sequencing_n_cells_per_well_lambda,
    sequencing_n_cells_per_well_nu = sequencing_n_cells_per_well_nu,
    pcr_factor_mu = pcr_factor_mu,
    pcr_factor_var = pcr_factor_var,
    n_reads_total = n_reads_total) |>
    dplyr::mutate(
      param_index = dplyr::row_number())

  if (verbose) {
    cat(
      "Generating and estimating ", nrow(parameter_sets), " ",
      "different datasets\n", sep = "")

    if (!is.null(progress_file)){
      if (!file.exists(progress_file)) {
        cat("Writing progress to ", progress_file, "\n", sep = "")
      } else {
        cat(
          "Warning: requesting writing to progress file that alread exists: ",
          progress_file, "\n", sep = "")
      }
    }
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
          "  gene_abundance_alpha = ", params$gene_abundance_alpha[1], "\n",
          "  gene_hit_rate = ", params$gene_hit_rate[1], "\n",
          "  n_wells_per_screen = ", params$n_wells_per_screen[1], "\n",
          "  well_abundance_factor_mu = ", params$well_abundance_factor_mu[1], "\n",
          "  well_abundance_factor_var = ", params$well_abundance_factor_var[1], "\n",
          "  imaging_n_cells_per_well_lambda = ", params$imaging_n_cells_per_well_lambda[1], "\n",
          "  imaging_n_cells_per_well_nu = ", params$imaging_n_cells_per_well_nu[1], "\n",
          "  class_pos_mu = ", params$class_pos_mu[1], "\n",
          "  class_pos_var = ", params$class_pos_var[1], "\n",
          "  class_neg_mu = ", params$class_neg_mu[1], "\n",
          "  class_neg_var = ", params$class_neg_var[1], "\n",
          "  sequencing_n_cells_per_well_lambda = ", params$sequencing_n_cells_per_well_lambda[1], "\n",
          "  sequencing_n_cells_per_well_nu = ", params$sequencing_n_cells_per_well_nu[1], "\n",
          "  pcr_factor_mu = ", params$pcr_factor_mu[1], "\n",
          "  pcr_factor_var = ", params$pcr_factor_var[1], "\n",
          "  n_reads_total = ", params$n_reads_total[1], "\n",
          sep = "")
      }
      data <- simulate_screen(
        n_genes_per_library = params$n_genes_per_library[1],
        gene_abundance_alpha = params$gene_abundance_alpha[1],
        gene_hit_rate = params$gene_hit_rate[1],
        n_wells_per_screen = params$n_wells_per_screen[1],
        well_abundance_factor_mu = params$well_abundance_factor_mu[1],
        well_abundance_factor_var = params$well_abundance_factor_var[1],
        imaging_n_cells_per_well_lambda = params$imaging_n_cells_per_well_lambda[1],
        imaging_n_cells_per_well_nu = params$imaging_n_cells_per_well_nu[1],
        class_pos_mu = params$class_pos_mu[1],
        class_pos_var = params$class_pos_var[1],
        class_neg_mu = params$class_neg_mu[1],
        class_neg_var = params$class_neg_var[1],
        sequencing_n_cells_per_well_lambda = params$sequencing_n_cells_per_well_lambda[1],
        sequencing_n_cells_per_well_nu = params$sequencing_n_cells_per_well_nu[1],
        pcr_factor_mu = params$pcr_factor_mu[1],
        pcr_factor_var = params$pcr_factor_var[1],
        n_reads_total = params$n_reads_total[1])

      model_data <- prepare_model_data(data)

      if (is.null(model)) {
        if (verbose) {
          cat("Compiling model for the first time...\n")
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

      model_evaluation <- data.frame(params) |>
        dplyr::bind_cols(model_evaluation)

      if (!is.null(progress_file)) {
        if (!file.exists(progress_file)) {
          model_evaluation |>
            readr::write_tsv(file = progress_file)
        } else {
          model_evaluation |>
            readr::write_tsv(
              file = progress_file,
              append = TRUE)
        }
      }

      model_evaluation |>
        dplyr::mutate(
          data = list(data),
          model = list(model_fit))
    }) |>
    dplyr::ungroup()
}
