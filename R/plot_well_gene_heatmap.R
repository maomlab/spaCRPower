
#' Plot well x gene heatmap
#'
#' Given a dataset plot heatmap where each box corresponds to a well
#' (row) by gene (column) pair. The values are the `log10(count + 1)` on
#' a scale from blue to yellow.
#'
#' @param data data.frame with columns \[`gene`, `well`, `count`\]
#'
#' @returns `ggplot2::ggplot` object
#'
#' @examples
#' \dontrun{
#'    data <- simulate_data(...)
#'    plot <- plot_well_gene_heatmap(data)
#'    ggplot2::ggsave(
#'      filename = "well_gene_heatmap.png",
#'      plot = plot,
#'      width = 5,
#'      height = 4)
#' }
#'
#' @export
plot_well_gene_heatmap <- function(data) {

  data |>
    dplyr::select(gene, well, count) |>
    tidyr::pivot_wider(
      id_cols = well,
      names_from = "gene",
      values_from = "count")

  ggplot2::ggplot(data = data) +
    ggplot2::theme_bw() +
    ggplot2::ggtitle("Cell Count for Gene by Well")
    ggplot2::geom_tile(
      mapping = ggplot2::aes(
        x = gene,
        y = well,
        fill = log10( count + 1))) +
    viridis::scale_fill_viridis(
      "Log10(Cell Count + 1)",
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
