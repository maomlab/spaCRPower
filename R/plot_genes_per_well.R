
#' Plot genes per well histogram
#'
#' From the given dataset, plot the distribution of genes per well
#' as a histogram
#'
#' @param data data.frame with columns \[`well`, `gene`,
#'   `n_cells_per_gene_per_well`\] and each row is a well x gene pair
#'
#' @returns `ggplot2::ggplot` object
#'
#' @examples
#' \dontrun{
#'    data <- simulate_data(...)
#'    plot <- plot_genes_per_well(data)
#'    ggplot2::ggsave(
#'      filename = "genes_per_well.png",
#'      plot = plot,
#'      width = 5,
#'      height = 4)
#' }
#'
#' @export
plot_genes_per_well <- function(data) {
  plot_data <- data |>
    dplyr::group_by(well) |>
    dplyr::summarize(
      n_genes = sum(gene_in_well > 0),
      .groups = "drop")

  plot_data_mean <- plot_data |>
    dplyr::summarize(
      mean_n_genes = mean(n_genes, na.rm =TRUE),
      .groups = "drop")

  max_height <- plot_data |>
    dplyr::count(n_genes) |>
    purrr::pluck("n") |>
    max()

  ggplot2::ggplot() +
    ggplot2::theme_bw() +
    ggplot2::geom_bar(
      data = plot_data,
      mapping = ggplot2::aes(
        x = n_genes)) +
    ggplot2::geom_vline(
      data = plot_data_mean,
      mapping = ggplot2::aes(
        xintercept = mean_n_genes)) +
    ggrepel::geom_text_repel(
      data = plot_data_mean,
      mapping = ggplot2::aes(
        x = mean_n_genes,
        label = signif(mean_n_genes, digits = 2),
        y = 1.1 * max_height),
      direction = "x",
      nudge_x = 1) +
    ggplot2::scale_y_continuous("Well Count") +
    ggplot2::scale_x_continuous("Genes per Well Count")
}
