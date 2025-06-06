
#' Plot well x gene Reads heatmap
#'
#' Given a dataset plot heatmap where each box corresponds to a well
#' (row) by gene (column) pair. The values are the `log10(reads + 1)` on
#' a scale from blue to yellow.
#'
#' @param data data.frame with columns \[`gene`, `well`,
#'   `n_reads_per_gene_per_well`\]
#'
#' @returns `ggplot2::ggplot` object
#'
#' @examples
#' \dontrun{
#'    data <- simulate_screen(...)
#'    plot <- plot_well_gene_reads_heatmap(data)
#'    ggplot2::ggsave(
#'      filename = "well_gene_reads_heatmap.png",
#'      plot = plot,
#'      width = 5,
#'      height = 4)
#' }
#'
#' @export
plot_well_gene_reads_heatmap <- function(data) {

  data |>
    dplyr::select(gene, well, n_reads_per_gene_per_well) |>
    tidyr::pivot_wider(
      id_cols = well,
      names_from = "gene",
      values_from = "n_reads_per_gene_per_well")

  ggplot2::ggplot(data = data) +
    ggplot2::theme_bw() +
    ggplot2::ggtitle("Read Count for Gene by Well") +
    ggplot2::geom_tile(
      mapping = ggplot2::aes(
        x = gene,
        y = well,
        fill = log10(n_reads_per_gene_per_well + 1))) +
    viridis::scale_fill_viridis(
      "Log10(Read Count + 1)",
      option = "cividis") +
    ggplot2::scale_x_continuous(
      "Gene ID",
      expand = c(0, 0)) +
    ggplot2::scale_y_continuous(
      "Well ID",
      expand = c(0, 0)) +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid = ggplot2::element_blank())
}
